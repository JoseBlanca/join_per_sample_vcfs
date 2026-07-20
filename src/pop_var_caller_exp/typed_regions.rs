//! The `type-regions` subcommand: run step 3's walk over a reference and
//! write the typed-region partition to a file. This module owns the `Args`
//! struct, the `run_typed_regions` driver, and its `#[non_exhaustive]`
//! error enum.
//!
//! **Complete.** The knob surface and its
//! [parsers](crate::pop_var_caller_exp::cli::parsers); the error enum; the
//! output writer ([`write_header`] + [`write_row`], carrying the T4 coordinate
//! conversion); the fallible setup ([`prepare_walk_inputs`]); and the driver
//! ([`run_typed_regions`]), which streams the partition to `<output>.tmp` and
//! joins the reference verification **before** renaming it into place.
//! See `doc/devel/ng/impl_plan/typed_regions_cli.md`.

use std::fs::File;
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use clap::Args;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::WindowedRefSeq;
use crate::ng::reference_info::{
    ReferenceInfoCache, ReferenceInfoError, VerificationHandle,
    read_reference_verifying_or_creating_fai,
};
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_FLANK_BP, DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY,
    DEFAULT_MIN_SCORE, MAX_MOTIF_LEN, MinCopies, RejectionCounts, SsrSegmentCriteria,
};
use crate::ng::region_typing::{
    DEFAULT_MAX_STR_LEN, DEFAULT_WINDOW_BP, GenomeRegions, RegionKind, TypedRegion,
    TypedRegionConfig, TypedRegionCounts, TypedRegionError, TypedRegionIterator,
};
use crate::ng::tandem_repeat::{
    DEFAULT_MATCH_REWARD, DEFAULT_MIN_COPIES, DEFAULT_MISMATCH_PENALTY, PeriodRange,
    PeriodRangeError, RepeatInterval, ScanParams,
};
use crate::ng::types::{Bp, ContigId, GenomeRegion};
use crate::pop_var_caller::common::DEFAULT_BUFFERED_IO_CAPACITY;
use crate::regions::{BedError, ContigBounds};

/// `type-regions` arguments — the authoritative knob list; `run_typed_regions`
/// (Milestone E) translates it into a `TypedRegionConfig`. Every knob defaults
/// to the library `pub const` that defines it, so the CLI's `Default` tracks the
/// short-read settings ng ships (spec §2.1, §2.3).
///
/// `--min-copies` is a table rather than a scalar, so it carries its own
/// [`value_parser`](crate::pop_var_caller_exp::cli::parsers::parse_min_copies)
/// (a `MinCopies` cannot derive a clap parser on its own).
#[derive(Debug, Args, Clone)]
pub struct TypedRegionsArgs {
    /// Reference FASTA. Its `.fai` is read, or created if absent (spec T1).
    #[arg(long)]
    pub reference: PathBuf,

    /// Output partition file (a plain-text BED-shaped TSV).
    #[arg(long)]
    pub output: PathBuf,

    /// BED of regions to emit; omit for the whole genome. A narrow BED still
    /// scans the whole contig, so it costs scan time, not memory or
    /// correctness (spec T10).
    #[arg(long)]
    pub regions: Option<PathBuf>,

    /// Narrowest STR period classified. One range detects and classifies
    /// (spec §2.2); the default of 1 types homopolymers as period-1 loci.
    ///
    /// Bounded at the CLI to 1..=6 — see `max_period` for why that is a usage
    /// error rather than a propagated one.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "Advanced"
    )]
    pub min_period: u8,

    /// Widest STR period classified — the microsatellite ceiling.
    ///
    /// **Bounded at the CLI, unlike the `--max-str-len`/`--flank-bp` pair.** That
    /// pair is propagated because the walk refuses it with a real *error* (T3);
    /// this one is not, because the walk guards the period ceiling with a release
    /// `assert!` — deliberately, since the knob is swept and sweeps run in
    /// `--release`. An assert is right for a library invariant and wrong for a
    /// flag: `--max-period 7` would be a backtrace, and spec §6 says a typo
    /// deserves a message. So clap rejects it as a usage error before `run`.
    #[arg(
        long,
        default_value_t = DEFAULT_MAX_PERIOD,
        value_parser = clap::value_parser!(u8).range(1..=MAX_MOTIF_LEN as i64),
        help_heading = "Advanced"
    )]
    pub max_period: u8,

    /// Satellite cap and scan margin (bp), one field: a tract longer than this
    /// is a `Satellite`, not an STR. Must be `>= --flank-bp` (the walk enforces
    /// it, spec T3).
    #[arg(long, default_value_t = DEFAULT_MAX_STR_LEN, help_heading = "Advanced")]
    pub max_str_len: u64,

    /// Window core (bp) — a memory knob only; it must not change the output.
    #[arg(long, default_value_t = DEFAULT_WINDOW_BP, help_heading = "Advanced")]
    pub window_bp: u64,

    /// Flank (bp) a locus needs each side, and the bundle radius (spec §2.4).
    ///
    /// Bounded below at 1 for the same reason as the period range: the walk
    /// guards `flank_bp >= 1` with a release `assert!` (at 0 every locus fails
    /// its own flank test), and that guard is reached on the **first window of
    /// any reference**, so a typo would be an immediate backtrace (spec §6).
    #[arg(
        long,
        default_value_t = DEFAULT_FLANK_BP,
        value_parser = clap::value_parser!(u64).range(1..),
        help_heading = "Advanced"
    )]
    pub flank_bp: u64,

    /// Purity floor in `[0, 1]`: a tract matching less than this fraction of a
    /// perfect motif tiling is not classified.
    ///
    /// Range-checked here because clap has no float range parser and the walk
    /// guards this with a release `assert!`. `nan` is the case that matters
    /// most — it parses happily and would otherwise pass *every* tract.
    #[arg(
        long,
        default_value_t = DEFAULT_MIN_PURITY,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_purity,
        help_heading = "Advanced"
    )]
    pub min_purity: f32,

    /// Score floor on the scanner's alignment score (`0` = a no-op for scanner
    /// output).
    #[arg(long, default_value_t = DEFAULT_MIN_SCORE, help_heading = "Advanced")]
    pub min_score: i32,

    /// Per-period copy-number floor: exactly six comma-separated values, one per
    /// period 1..=6. Below its period's floor a tract is not classified as an
    /// STR. Any other count is a hard parse error.
    #[arg(
        long,
        value_parser = crate::pop_var_caller_exp::cli::parsers::parse_min_copies,
        default_value = "6,4,4,3,3,3",
        help_heading = "Advanced"
    )]
    pub min_copies: MinCopies,

    /// Scanner match reward.
    #[arg(long, default_value_t = DEFAULT_MATCH_REWARD, help_heading = "Advanced")]
    pub scan_match_reward: i32,

    /// Scanner mismatch penalty.
    #[arg(long, default_value_t = DEFAULT_MISMATCH_PENALTY, help_heading = "Advanced")]
    pub scan_mismatch_penalty: i32,

    /// Scanner minimum copies — a permissive detection floor. The stricter,
    /// per-period classification floor is --min-copies.
    #[arg(long, default_value_t = DEFAULT_MIN_COPIES, help_heading = "Advanced")]
    pub scan_min_copies: u32,
}

/// Errors from the `type-regions` subcommand. The house pattern: a
/// `#[non_exhaustive]` `thiserror` enum, `#[from]` where a single source type
/// identifies the failure. `main_exp` walks the source chain to render it (spec
/// §6). **Not** a `--max-str-len`/`--flank-bp` variant — the walk's
/// [`TypedRegionError`] already carries both numbers (spec T3).
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum TypedRegionsCliError {
    /// Reading the reference, verifying it against its `.fai`, or writing a
    /// missing sibling `.fai` failed (spec T1).
    #[error("could not read the reference")]
    Reference(#[from] ReferenceInfoError),

    /// A contig is longer than 4 Gb, so its length does not fit the `u32` a
    /// `RegionSet` needs (spec T2). ng fails rather than folds.
    #[error("contig '{contig}' is {len} bp, too long for a 32-bit region set (> 4 Gb)")]
    ContigTooLong { contig: String, len: u64 },

    /// Parsing the `--regions` BED failed.
    #[error("could not read the --regions BED")]
    Bed(#[from] BedError),

    /// `--min-period` exceeds `--max-period`, so the range would scan nothing.
    ///
    /// The *bounds* of each flag are clap's to reject (see
    /// [`TypedRegionsArgs::max_period`]); this is the cross-flag rule, which clap
    /// cannot express, so it surfaces as an error rather than the panic
    /// `PeriodRange::new`'s caller would otherwise take (spec §6).
    #[error("--min-period and --max-period do not form a range")]
    PeriodRange(#[from] PeriodRangeError),

    /// The walk's fallible setup or streaming failed — including the
    /// `--max-str-len`/`--flank-bp` guard, which the walk carries (spec T3).
    #[error("the typed-region walk failed")]
    Walk(#[from] TypedRegionError),

    /// Writing the output file, or renaming the temp file into place, failed.
    #[error("could not write the output file")]
    Output(#[from] std::io::Error),
}

// ---------------------------------------------------------------------
// The output format (spec §3)
// ---------------------------------------------------------------------

/// The file's columns, in order — and the order the row formatter must follow.
///
/// The first four are always filled, so the file is a valid BED to any tool that
/// reads the first three; the rest are filled only by the kinds that have them
/// (spec §3.1). This is also what our own `--regions` accepts back, because
/// `RegionSet::from_bed_path` ignores `#` comments and every column past the
/// third (spec §3.2).
const COLUMNS: [&str; 9] = [
    "chrom", "start", "end", "kind", "motif", "period", "copies", "purity", "members",
];

/// Render a [`MinCopies`] table in `--min-copies` syntax — the six per-period
/// floors, comma-separated — so the header round-trips straight back into a
/// re-run of the command that produced it.
fn format_min_copies(min_copies: &MinCopies) -> String {
    (1..=MAX_MOTIF_LEN)
        .map(|period| {
            let period = u8::try_from(period).expect("MAX_MOTIF_LEN fits in a u8");
            min_copies.for_period(period).to_string()
        })
        .collect::<Vec<_>>()
        .join(",")
}

/// Write the `##` provenance block, then the `#`-prefixed column header.
///
/// **Every resolved config value appears** (spec §3.4) — not just the flags the
/// user typed, or a run at defaults would record nothing and be incomparable to
/// one that set them explicitly.
///
/// The block is written from an **exhaustive destructuring** of
/// [`TypedRegionConfig`], [`SsrSegmentCriteria`] and [`ScanParams`] (no `..`), so
/// a knob added to any of the three cannot be silently forgotten here. Precisely:
/// a *new* field fails this function's compile outright, and a field bound but
/// never written is an `unused_variables` warning — which the pre-commit check
/// promotes to an error (`clippy … -D warnings`, `scripts/precommit-check.sh`).
/// The device stops at the leaves: a new field on `MinCopies` or `PeriodRange`
/// would not trip it.
///
/// **No `date` and no `reference_md5`** (spec §3.4, §6). The walk is a pure
/// function of reference + regions + config, so the same three inputs must give a
/// byte-identical file: a timestamp would break that outright, and the digest is
/// only known when the background FASTA verification joins — *after* the walk,
/// whereas this header is written before it, so recording it would force the walk
/// to block on the very pass that runs beside it (spec T1).
///
/// `window_bp` is recorded but is **not a comparison key**: it is a memory knob
/// that must not change the output, so two files differing only in `window_bp`
/// should be identical below the header.
pub fn write_header<W: io::Write>(
    out: &mut W,
    reference: &Path,
    config: &TypedRegionConfig,
) -> io::Result<()> {
    // Exhaustive on purpose — see the doc above. Do not replace with `..`.
    let TypedRegionConfig {
        scan,
        max_str_len,
        window_bp,
        criteria,
    } = config;
    let SsrSegmentCriteria {
        periods,
        min_purity,
        min_score,
        flank_bp,
        min_copies,
    } = criteria;
    let ScanParams {
        match_reward,
        mismatch_penalty,
        min_copies: scan_min_copies,
    } = scan;

    writeln!(out, "## tool: pop_var_caller_exp type-regions")?;
    writeln!(out, "## version: {}", env!("CARGO_PKG_VERSION"))?;
    writeln!(out, "## reference: {}", reference.display())?;
    writeln!(out, "## min_period: {}", periods.min())?;
    writeln!(out, "## max_period: {}", periods.max())?;
    writeln!(out, "## max_str_len: {}", max_str_len.get())?;
    writeln!(out, "## window_bp: {}", window_bp.get())?;
    writeln!(out, "## flank_bp: {flank_bp}")?;
    writeln!(out, "## min_purity: {min_purity}")?;
    writeln!(out, "## min_score: {min_score}")?;
    writeln!(out, "## min_copies: {}", format_min_copies(min_copies))?;
    writeln!(out, "## scan_match_reward: {match_reward}")?;
    writeln!(out, "## scan_mismatch_penalty: {mismatch_penalty}")?;
    writeln!(out, "## scan_min_copies: {scan_min_copies}")?;
    writeln!(out, "#{}", COLUMNS.join("\t"))?;
    Ok(())
}

/// The `kind` column's vocabulary (spec §3.1). `awk '$4=="ssr_locus"'` is the
/// whole query language, so these strings are part of the format.
fn kind_label(kind: &RegionKind) -> &'static str {
    match kind {
        RegionKind::SsrSegment(_) => "ssr_locus",
        RegionKind::SsrBundle { .. } => "ssr_bundle",
        RegionKind::Generic => "generic",
        RegionKind::Satellite => "satellite",
    }
}

/// **T4 — the one place the row's two coordinate systems meet.**
///
/// A [`GenomeRegion`] hull is 1-based **inclusive**; BED is 0-based
/// **half-open**. So the start loses one and the end does not move: 1-based
/// `[1, 10]` is ten bases, and so is 0-based `[0, 10)`.
///
/// A bundle's member `RepeatInterval`s are **not** converted — they are already
/// 0-based half-open *and* already re-based to contig coordinates by the walk, so
/// they are written through unchanged (see [`members_json`]). That is the trap:
/// one row, two coordinate systems, in opposite directions. The conversion lives
/// here, once, or it gets written twice and one of them is wrong — and a wrong
/// one is a *wrong region*, not a panic.
fn hull_to_bed(region: &GenomeRegion) -> (u64, u64) {
    let start = region.start.get();
    let end = region.end.get();
    // Both halves of the invariant, because both failures are silently-invalid
    // BED rather than a crash: a zero start would underflow, and an inverted
    // region would emit `99 50`. `GenomeRegion` has public fields and no
    // validating constructor, so nothing upstream enforces this for us.
    assert!(
        start >= 1 && end >= start,
        "a GenomeRegion is 1-based inclusive and non-empty; got [{start}, {end}] (spec §4)"
    );
    (start - 1, end)
}

/// A bundle's members as a **JSON array** — the one genuinely nested field, and
/// the only cell that carries JSON (spec §3.1). Fixed key order and no
/// incidental whitespace, so the cell serialises deterministically for §6's
/// byte-identity.
///
/// The intervals are written through **unshifted**: they are already 0-based
/// half-open contig coordinates (T4, and see [`hull_to_bed`] for the hull, which
/// is not).
fn members_json(tracts: &[RepeatInterval]) -> String {
    let mut out = String::from("[");
    for (i, tract) in tracts.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&format!(
            r#"{{"start":{},"end":{},"period":{}}}"#,
            tract.start, tract.end, tract.period
        ));
    }
    out.push(']');
    out
}

/// Resolve a [`ContigId`] to its name through the run's one contig table.
///
/// Every row is named from this table rather than from
/// [`SsrSegment::chrom`](crate::ng::region_typing::segment_criteria::SsrSegment::chrom):
/// `Generic`, `Satellite` and `SsrBundle` carry no segment, so the table is
/// needed anyway (spec §4), and naming every row from one source is what stops
/// the two ever disagreeing.
///
/// The lookup cannot miss: T8 requires the *same* `ContigList` to feed both
/// `WindowedRefSeq` and `GenomeRegions`, so every `ContigId` the walk emits
/// indexes this table. A miss is an internal wiring bug, not user input.
fn contig_name(contigs: &ContigList, contig: ContigId) -> &str {
    contigs
        .entries
        .get(contig.get() as usize)
        .map(|entry| entry.name.as_str())
        .expect("the contig table that named the walk's regions must contain every ContigId it emits (spec T8)")
}

/// Write one [`TypedRegion`] as a BED row: the four always-filled columns, then
/// the typed columns with `.` for absent (spec §3.1, §4).
///
/// `copies` is **derived and fractional** (T5): the number of repeats is stored
/// nowhere — it is `tract_len / period`, and a tract is rarely a whole number of
/// motifs, which is part of what purity measures. trf-mod's own BED calls the
/// fractional value `copyNum`, so a fractional column is the field's convention.
///
/// Floats render via `{:?}`, which is the shortest **round-tripping** form *and*
/// keeps the decimal point on whole values — `10.0`, not `10`, matching spec
/// §3.1's worked row. Plain `{}` (the catalog writer's choice) would drop it and
/// make the column's type depend on its value.
pub fn write_row<W: io::Write>(
    out: &mut W,
    region: &TypedRegion,
    contigs: &ContigList,
) -> io::Result<()> {
    let (start, end) = hull_to_bed(&region.region);
    let chrom = contig_name(contigs, region.region.contig);
    let kind = kind_label(&region.kind);

    match &region.kind {
        RegionKind::SsrSegment(segment) => {
            // T6: `Motif` has no `Display`, only a `Debug` that prints quotes.
            // `motif()` returns by value (it is `Copy`), so bind it before
            // borrowing its bytes.
            let motif = segment.motif();
            let motif_str = std::str::from_utf8(motif.as_bytes())
                .expect("motif bytes are ASCII by construction");
            let period = segment.period();
            let copies = segment.tract_len() as f64 / period as f64;
            writeln!(
                out,
                "{chrom}\t{start}\t{end}\t{kind}\t{motif_str}\t{period}\t{copies:?}\t{:?}\t.",
                segment.purity_fraction()
            )
        }
        RegionKind::SsrBundle { tracts } => writeln!(
            out,
            "{chrom}\t{start}\t{end}\t{kind}\t.\t.\t.\t.\t{}",
            members_json(tracts)
        ),
        RegionKind::Generic | RegionKind::Satellite => {
            writeln!(out, "{chrom}\t{start}\t{end}\t{kind}\t.\t.\t.\t.\t.")
        }
    }
}

// ---------------------------------------------------------------------
// The fallible setup (spec T1, T2, T8; arch §4)
// ---------------------------------------------------------------------

/// The reference-derived inputs the walk needs, assembled once.
///
/// **Holds one [`ContigList`], which is the point** (spec T8). The walk's
/// `over_regions` is fallible because a [`GenomeRegions`] was validated against
/// *a* contig table and nothing ties that one to the reference's; if they
/// disagree the failure is `UnknownContig(ContigId)`, which renders as an index,
/// not a name. Feeding this one table to both `GenomeRegions` (here) and
/// `WindowedRefSeq` (in the driver) makes that error unreachable — but only if
/// the caller uses *this* list rather than reading the reference again.
pub struct WalkInputs {
    /// The reference's contig table — the one table, used twice.
    pub contigs: ContigList,
    /// What to emit. A BED chooses what is emitted, never what is scanned
    /// (spec T10).
    pub regions: GenomeRegions,
    /// The background FASTA verification, when a `.fai` was read in the
    /// foreground. **`join()` it before the output is published** (spec T1, §6) —
    /// it is what catches a stale `.fai`, and joining after the rename would
    /// publish a partition built against a wrong contig table.
    pub verify: Option<VerificationHandle>,
}

/// Narrow the reference's contig lengths to the `u32` a [`RegionSet`] needs,
/// **failing rather than folding** (spec T2).
///
/// `ContigEntry::length` is `u64`; `ContigBounds::length` is `u32`. The only
/// precedent in the tree casts (`e.length as u32`) and it is a test. B2 deleted
/// exactly that shape from ng and recorded the rule: a silent narrowing would
/// make a >4 Gb contig report a wrong length, and a region set built on it would
/// clamp spans to nonsense. A contig that large is unrepresentable to `RegionSet`
/// whatever we do, so failing is the honest answer, not a limitation we chose.
fn contig_bounds(contigs: &ContigList) -> Result<Vec<ContigBounds<'_>>, TypedRegionsCliError> {
    contigs
        .entries
        .iter()
        .map(|entry| {
            let length =
                u32::try_from(entry.length).map_err(|_| TypedRegionsCliError::ContigTooLong {
                    contig: entry.name.clone(),
                    len: entry.length,
                })?;
            Ok(ContigBounds {
                name: &entry.name,
                length,
            })
        })
        .collect()
}

/// Read the reference and resolve what to emit — the walk's fallible setup.
///
/// The reference goes through `reference_info`'s batteries-included entry point
/// (spec T1), which has two cases:
///
/// - ***`.fai` present*** → the index is read in the foreground (cheap, so the
///   walk can start at once) and the FASTA is verified against it on a
///   background thread, returning `Some(handle)`. A stale `.fai` is caught
///   there rather than by producing a wrong partition.
/// - ***`.fai` absent*** → the FASTA is read once, the sibling `.fai` is written
///   so the next run takes the fast path, and `None` comes back. A `.fai`-write
///   failure is fatal.
///
/// The `#[must_use]` handle is returned for the caller to join **before it
/// publishes output**, never after.
///
/// **Known wart, on the error paths only.** If the region set fails to build
/// (a bad BED, an over-long contig) after a `.fai` was read, the pending handle
/// is dropped un-joined and `VerificationHandle`'s `Drop` prints its "dropped
/// without join()" warning beside the real error. Nothing is published on those
/// paths, so abandoning the check is legitimate — the warning is misleading, not
/// a correctness problem. Joining instead would stall a *failing* run on a whole
/// FASTA read, which is worse. A clean fix needs a way to abandon a handle
/// deliberately (a `VerificationHandle::abandon`), which lives in
/// `reference_info` rather than here.
pub fn prepare_walk_inputs(
    args: &TypedRegionsArgs,
    cache: &Arc<ReferenceInfoCache>,
) -> Result<WalkInputs, TypedRegionsCliError> {
    let (info, verify) = read_reference_verifying_or_creating_fai(cache, args.reference.clone())?;
    let contigs = info.contig_list();

    // `bounds` borrows `contigs`; the region set owns its data, so the one table
    // is free to travel on to `WindowedRefSeq` (T8). The block is for legibility
    // — under NLL the borrow would end at `bounds`'s last use anyway.
    let regions = {
        let bounds = contig_bounds(&contigs)?;
        match args.regions.as_deref() {
            Some(bed) => GenomeRegions::from_bed_path(bed, &bounds)?,
            None => GenomeRegions::whole_contigs(&bounds),
        }
    };

    Ok(WalkInputs {
        contigs,
        regions,
        verify,
    })
}

/// Translate the args into the walk's config.
///
/// The only fallible part is the period *range*: each flag's bounds are clap's
/// to reject, but the cross-flag `min <= max` rule is not something clap can
/// express, so it surfaces here as an error rather than the panic
/// `PeriodRange::new`'s caller would otherwise take (spec §6).
fn walk_config(args: &TypedRegionsArgs) -> Result<TypedRegionConfig, TypedRegionsCliError> {
    Ok(TypedRegionConfig {
        scan: ScanParams {
            match_reward: args.scan_match_reward,
            mismatch_penalty: args.scan_mismatch_penalty,
            min_copies: args.scan_min_copies,
        },
        max_str_len: Bp(args.max_str_len),
        window_bp: Bp(args.window_bp),
        criteria: SsrSegmentCriteria {
            periods: PeriodRange::new(args.min_period, args.max_period)?,
            min_purity: args.min_purity,
            min_score: args.min_score,
            flank_bp: args.flank_bp,
            min_copies: args.min_copies,
        },
    })
}

/// `<output>.tmp` — the name the partition is built under before it is published.
fn tmp_path_for(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_os_string();
    name.push(".tmp");
    PathBuf::from(name)
}

/// Report what the walk tallied, to stderr (spec §6, T9).
///
/// **All five rejection counters, labelled, with none of them hard-coded away.**
/// Four are structurally zero today, but *which* four moves whenever an upstream
/// stage changes, so the summary states them flatly rather than explaining any
/// one of them — a hidden zero is worse than an explained one, and prose about
/// "no impure tracts in this genome" would be a wrong answer wearing a
/// measurement's clothes. The destructuring is exhaustive so a counter added
/// later cannot go unreported.
fn report_counts(counts: &TypedRegionCounts) {
    // Exhaustive on BOTH structs (no `..`), so a counter added to either cannot
    // go unreported — the same device `write_header` uses for the config.
    let TypedRegionCounts {
        spans,
        ssr_loci,
        ssr_bundles,
        ssr_bundle_bp,
        generic,
        satellites,
        satellite_bp,
        repeat_bp_with_no_locus,
        rejected_by_reason,
    } = counts;
    let RejectionCounts {
        copy_floor,
        purity,
        compound,
        no_clean_trim,
        flank_clamped,
    } = rejected_by_reason;

    eprintln!(
        "counts: spans={spans} ssr_loci={ssr_loci} ssr_bundles={ssr_bundles} \
         ssr_bundle_bp={ssr_bundle_bp} generic={generic} satellites={satellites} \
         satellite_bp={satellite_bp} repeat_bp_with_no_locus={repeat_bp_with_no_locus}"
    );
    // Stated as what it is: these do NOT partition `repeat_bp_with_no_locus`
    // (overlapping rejected repeats are each charged, and bases with no locus for
    // a reason that is not a rejection — bundled, capped, out of scope — are in
    // the total but under no reason here).
    eprintln!(
        "rejected_bp (does not partition repeat_bp_with_no_locus): copy_floor={copy_floor} \
         purity={purity} compound={compound} no_clean_trim={no_clean_trim} \
         flank_clamped={flank_clamped}"
    );
}

/// Run the walk and write the partition.
///
/// **Streams** — the walk holds one window and three coordinates, and collecting
/// a whole-genome partition would put back exactly the memory it exists to avoid
/// (spec §6). It needs no sort either: the walk emits in genomic order.
///
/// **Writes atomically, and the verification joins into that same barrier.** The
/// rows go to `<output>.tmp`; only a clean `VerificationHandle::join` renames it
/// into place. Both halves matter and for the same reason — *a wrong file here is
/// silently valid*: a partition cut off halfway through a chromosome is
/// indistinguishable from a complete partition of a smaller genome, and one built
/// against a stale `.fai` is indistinguishable from a correct one. So the join
/// goes **before** the rename (spec T1, §6), never after.
pub fn run_typed_regions(args: &TypedRegionsArgs) -> Result<(), TypedRegionsCliError> {
    // The config depends on nothing but the args, so it is built FIRST: a bad
    // flag pair should not cost a whole reference read (and, until
    // `reference_info` grows a way to abandon a handle deliberately, should not
    // trip the dropped-without-join warning either — see `prepare_walk_inputs`).
    let config = walk_config(args)?;

    let cache = Arc::new(ReferenceInfoCache::new());
    let WalkInputs {
        contigs,
        regions,
        verify,
    } = prepare_walk_inputs(args, &cache)?;

    eprintln!(
        "type-regions: reference={} regions={} output={}",
        args.reference.display(),
        args.regions
            .as_ref()
            .map_or_else(|| "<whole genome>".to_string(), |b| b.display().to_string()),
        args.output.display(),
    );

    // The reference and the row namer take the same contig table (T8) — the walk
    // needs it by value, so it is cloned from the one `contig_list()` D1 read.
    // A clone of one value cannot disagree with itself, which is the property T8
    // is protecting.
    let reference = WindowedRefSeq::new(args.reference.clone(), contigs.clone());
    // Fallible setup: the contig cross-check and the `--max-str-len`/`--flank-bp`
    // pair both surface here as errors, before any work (spec T3).
    let mut walk = TypedRegionIterator::over_regions(reference, regions, config.clone())?;

    let tmp_path = tmp_path_for(&args.output);
    let file = File::create(&tmp_path)?;
    let mut out = BufWriter::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file);

    write_header(&mut out, &args.reference, &config)?;
    for region in walk.by_ref() {
        write_row(&mut out, &region?, &contigs)?;
    }

    // The commit point. The background FASTA check overlapped the whole walk, so
    // this costs almost nothing — and a stale `.fai` aborts here with the
    // partition still under its temp name, never published.
    if let Some(handle) = verify {
        handle.join()?;
    }

    let file = out.into_inner().map_err(|e| e.into_error())?;
    file.sync_all()?;
    std::fs::rename(&tmp_path, &args.output)?;

    report_counts(walk.counts());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};
    use clap::Parser;

    /// The subcommand parses through the top-level CLI, and every knob resolves
    /// to its short-read §2.3 default. Mirrors `ssr-catalog`'s
    /// `parses_through_the_top_level_cli`.
    #[test]
    fn parses_through_the_top_level_cli() {
        let cli = Cli::try_parse_from([
            "pop_var_caller_exp",
            "type-regions",
            "--reference",
            "r.fa",
            "--output",
            "regions.tsv",
        ])
        .unwrap();
        let PopVarCallerExpCommand::TypeRegions(args) = cli.cmd;
        assert_eq!(args.reference, PathBuf::from("r.fa"));
        assert_eq!(args.output, PathBuf::from("regions.tsv"));
        assert_eq!(args.regions, None);

        // The short-read §2.3 defaults resolve through the CLI.
        assert_eq!(args.min_period, 1, "period-1 homopolymers classified");
        assert_eq!(args.max_period, 6);
        assert_eq!(args.max_str_len, 100, "short-read satellite cap");
        assert_eq!(args.window_bp, 100_000);
        assert_eq!(args.flank_bp, 30, "short-read flank");
        assert_eq!(args.min_purity, 0.8);
        assert_eq!(args.min_score, 0);
        assert_eq!(args.scan_match_reward, 2);
        assert_eq!(args.scan_mismatch_penalty, 7);
        assert_eq!(args.scan_min_copies, 2);

        // The `--min-copies` table resolves through its own value_parser, to the
        // short-read floors — and to the SAME table the library ships, so a floor
        // sweep (spec §10) that moves one without the other fails here rather than
        // silently applying two different defaults.
        let floors: Vec<u32> = (1..=6).map(|p| args.min_copies.for_period(p)).collect();
        assert_eq!(floors, vec![6, 4, 4, 3, 3, 3], "the short-read floors");
        assert_eq!(
            args.min_copies,
            MinCopies::default(),
            "the CLI default must track MinCopies::default()"
        );
    }

    // ---- the header block (spec §3.4) ----------------------------------

    fn header_of(config: &TypedRegionConfig) -> String {
        let mut buf = Vec::new();
        write_header(&mut buf, Path::new("/refs/ref.fa"), config).expect("writing to a Vec");
        String::from_utf8(buf).expect("the header is UTF-8")
    }

    /// The header is **byte-identical across two calls** on the same config —
    /// the determinism §6 makes the regression anchor. A `date` field (which
    /// `ssr-catalog`'s header carries) would break this outright.
    #[test]
    fn the_header_is_byte_identical_across_calls() {
        let config = TypedRegionConfig::default();
        assert_eq!(header_of(&config), header_of(&config));
    }

    /// **Every knob appears, resolved** — a file at defaults must record the same
    /// keys as one that set them explicitly, or two files cannot be compared
    /// (spec §3.4). `window_bp` is among them: recorded for reproducibility even
    /// though it is a memory knob, not a comparison key.
    #[test]
    fn the_header_carries_every_resolved_knob() {
        let header = header_of(&TypedRegionConfig::default());
        for key in [
            "tool",
            "version",
            "reference",
            "min_period",
            "max_period",
            "max_str_len",
            "window_bp",
            "flank_bp",
            "min_purity",
            "min_score",
            "min_copies",
            "scan_match_reward",
            "scan_mismatch_penalty",
            "scan_min_copies",
        ] {
            assert!(
                header.contains(&format!("## {key}: ")),
                "the header must record '{key}':\n{header}"
            );
        }
    }

    /// The knobs are recorded at their **resolved values**, not their names only.
    ///
    /// Compares **whole lines**, not substrings: `contains("## max_period: 6")`
    /// is satisfied by `## max_period: 60`, so a drifted floor (`3` → `30`) would
    /// leave this green while the header — whose whole job is making two sweep
    /// files comparable — no longer says what it claims.
    #[test]
    fn the_header_records_the_short_read_defaults() {
        let header = header_of(&TypedRegionConfig::default());
        let lines: Vec<&str> = header.lines().collect();
        for line in [
            "## min_period: 1",
            "## max_period: 6",
            "## max_str_len: 100",
            "## window_bp: 100000",
            "## flank_bp: 30",
            "## min_purity: 0.8",
            "## min_score: 0",
            "## min_copies: 6,4,4,3,3,3",
            "## scan_match_reward: 2",
            "## scan_mismatch_penalty: 7",
            "## scan_min_copies: 2",
            "## reference: /refs/ref.fa",
        ] {
            assert!(
                lines.contains(&line),
                "expected the exact line '{line}' in:\n{header}"
            );
        }
    }

    /// **No `date`, no `reference_md5`** (spec §3.4, §6) — the two fields that
    /// would make the file non-reproducible, or force the walk to block on the
    /// background verification pass.
    #[test]
    fn the_header_carries_no_timestamp_and_no_digest() {
        let header = header_of(&TypedRegionConfig::default());
        assert!(!header.contains("date"), "a timestamp breaks determinism");
        assert!(
            !header.contains("md5"),
            "the digest lands only when the background verification joins, after the walk"
        );
    }

    /// The `min_copies` line is in **`--min-copies` syntax**, so the header feeds
    /// straight back into a re-run of the command that produced the file.
    #[test]
    fn the_header_min_copies_round_trips_through_the_flag_parser() {
        use crate::pop_var_caller_exp::cli::parsers::parse_min_copies;
        let config = TypedRegionConfig::default();
        let header = header_of(&config);
        let recorded = header
            .lines()
            .find_map(|l| l.strip_prefix("## min_copies: "))
            .expect("the header records min_copies");
        assert_eq!(
            parse_min_copies(recorded).expect("the recorded table re-parses"),
            config.criteria.min_copies,
            "the header must round-trip into the flag that produced it"
        );
    }

    /// The column header is `#`-prefixed and names the nine columns in order, so
    /// the file is a valid BED to anything reading the first three (spec §3.1).
    #[test]
    fn the_column_header_names_the_columns_in_order() {
        let header = header_of(&TypedRegionConfig::default());
        let columns = header.lines().last().expect("a last line");
        assert_eq!(
            columns,
            "#chrom\tstart\tend\tkind\tmotif\tperiod\tcopies\tpurity\tmembers"
        );
    }

    /// A non-default config is recorded at *its* values, not the defaults — the
    /// header is a function of the resolved config, not a constant.
    #[test]
    fn the_header_follows_a_non_default_config() {
        use crate::ng::tandem_repeat::PeriodRange;
        use crate::ng::types::Bp;
        let config = TypedRegionConfig {
            max_str_len: Bp(250),
            criteria: SsrSegmentCriteria {
                periods: PeriodRange::new(2, 5).expect("2..=5 is valid"),
                min_copies: MinCopies::new([9, 8, 7, 6, 5, 4], 3),
                flank_bp: 44,
                ..SsrSegmentCriteria::default()
            },
            ..TypedRegionConfig::default()
        };
        let header = header_of(&config);
        let lines: Vec<&str> = header.lines().collect();
        for line in [
            "## min_period: 2",
            "## max_period: 5",
            "## max_str_len: 250",
            "## flank_bp: 44",
            "## min_copies: 9,8,7,6,5,4",
        ] {
            assert!(
                lines.contains(&line),
                "expected the exact line '{line}' in:\n{header}"
            );
        }
    }

    // ---- the row formatter and the T4 conversion (spec §3.1, §4, T4/T5/T6) ----

    use crate::fasta::ContigEntry;
    use crate::ng::region_typing::segment_criteria::{Motif, SsrSegment};
    use crate::ng::types::Position;

    /// One row, parsed back into its fields. The **oracle's other half**: the
    /// round-trip is only a check if the parse is independent of the format
    /// function, so this reads the columns positionally rather than reusing any
    /// writer helper.
    #[derive(Debug, PartialEq)]
    struct ParsedRow {
        chrom: String,
        start: u64,
        end: u64,
        kind: String,
        motif: Option<String>,
        period: Option<usize>,
        copies: Option<f64>,
        purity: Option<f32>,
        /// `(start, end, period)` per member.
        members: Option<Vec<(u64, u64, u8)>>,
    }

    fn cell(raw: &str) -> Option<&str> {
        if raw == "." { None } else { Some(raw) }
    }

    /// Parse the `members` JSON cell. Deliberately a small hand parser: the cell
    /// is a fixed shape, and depending on a JSON library here would test the
    /// library rather than our serialisation.
    fn parse_members(raw: &str) -> Vec<(u64, u64, u8)> {
        let inner = raw
            .strip_prefix('[')
            .and_then(|r| r.strip_suffix(']'))
            .expect("the members cell is a JSON array");
        if inner.is_empty() {
            return Vec::new();
        }
        inner
            .split("},{")
            .map(|object| {
                let object = object.trim_start_matches('{').trim_end_matches('}');
                let mut start = None;
                let mut end = None;
                let mut period = None;
                for field in object.split(',') {
                    let (key, value) = field.split_once(':').expect("a JSON key:value pair");
                    let value = value.trim();
                    match key.trim().trim_matches('"') {
                        "start" => start = Some(value.parse().expect("start is a number")),
                        "end" => end = Some(value.parse().expect("end is a number")),
                        "period" => period = Some(value.parse().expect("period is a number")),
                        other => panic!("unexpected key in the members cell: {other}"),
                    }
                }
                (
                    start.expect("start"),
                    end.expect("end"),
                    period.expect("period"),
                )
            })
            .collect()
    }

    fn parse_row(line: &str) -> ParsedRow {
        let columns: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            columns.len(),
            COLUMNS.len(),
            "every row has one cell per column: {line}"
        );
        ParsedRow {
            chrom: columns[0].to_string(),
            start: columns[1].parse().expect("start is a number"),
            end: columns[2].parse().expect("end is a number"),
            kind: columns[3].to_string(),
            motif: cell(columns[4]).map(str::to_string),
            period: cell(columns[5]).map(|c| c.parse().expect("period is a number")),
            copies: cell(columns[6]).map(|c| c.parse().expect("copies is a number")),
            purity: cell(columns[7]).map(|c| c.parse().expect("purity is a number")),
            members: cell(columns[8]).map(parse_members),
        }
    }

    fn row_of(region: &TypedRegion, contigs: &ContigList) -> ParsedRow {
        let mut buf = Vec::new();
        write_row(&mut buf, region, contigs).expect("writing to a Vec");
        let text = String::from_utf8(buf).expect("the row is UTF-8");
        let line = text.strip_suffix('\n').expect("a row ends in a newline");
        assert!(!line.contains('\n'), "one region is exactly one row");
        parse_row(line)
    }

    fn test_contigs() -> ContigList {
        ContigList {
            entries: vec![
                ContigEntry {
                    name: "chr1".to_string(),
                    length: 100_000,
                    md5: None,
                },
                ContigEntry {
                    name: "chr2".to_string(),
                    length: 50_000,
                    md5: None,
                },
            ],
        }
    }

    /// A hull at 1-based inclusive `[start, end]` on `chr1`.
    fn hull(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    fn tract(start: u64, end: u64, period: u8) -> RepeatInterval {
        RepeatInterval {
            start,
            end,
            period,
            score: 100,
        }
    }

    /// **THE ROUND-TRIP ORACLE (plan C2, spec §7).** Format each kind, parse the
    /// row back, and assert it equals the input field for field — including a
    /// bundle, without which the member coordinates are never exercised.
    ///
    /// This is what pins **T4's two directions in one row**: the hull's start
    /// loses one (1-based inclusive → 0-based half-open) while the members, which
    /// are already 0-based half-open contig coordinates, do not move. An
    /// off-by-one in either is a *wrong region*, not a panic — so it is checked
    /// against the input rather than against the writer's own idea of itself.
    #[test]
    fn every_kind_round_trips_through_a_row() {
        let contigs = test_contigs();

        // --- Generic: the region IS the whole claim; everything else absent ---
        let generic = TypedRegion {
            region: hull(1, 999),
            kind: RegionKind::Generic,
        };
        let row = row_of(&generic, &contigs);
        assert_eq!(row.chrom, "chr1");
        assert_eq!(row.start, 0, "1-based inclusive 1 is 0-based 0 (T4)");
        assert_eq!(row.end, 999, "the end does not move (T4)");
        assert_eq!(row.kind, "generic");
        assert_eq!(
            (row.motif, row.period, row.copies, row.purity, row.members),
            (None, None, None, None, None),
            "Generic fills nothing but the region"
        );

        // --- Satellite: likewise, and the extent IS the claim ---
        let satellite = TypedRegion {
            region: hull(2041, 4240),
            kind: RegionKind::Satellite,
        };
        let row = row_of(&satellite, &contigs);
        assert_eq!((row.start, row.end), (2040, 4240));
        assert_eq!(row.kind, "satellite");
        assert_eq!(row.members, None);

        // --- SsrSegment: the region IS the tract; motif/period/copies/purity ---
        let segment = SsrSegment::new(
            "chr1".into(),
            1000, // 1-based inclusive
            1029,
            Motif::new(b"CAG").expect("CAG is a valid motif"),
            1.0,
        )
        .expect("a valid segment");
        let locus = TypedRegion {
            region: hull(1000, 1029),
            kind: RegionKind::SsrSegment(segment),
        };
        let row = row_of(&locus, &contigs);
        assert_eq!((row.start, row.end), (999, 1029), "the hull shifts (T4)");
        assert_eq!(row.kind, "ssr_locus");
        assert_eq!(row.motif.as_deref(), Some("CAG"), "T6: bytes, not Debug");
        assert_eq!(row.period, Some(3));
        // T5: derived, and fractional in general — 30 bp of CAG is 10.0 copies.
        assert_eq!(row.copies, Some(10.0));
        assert_eq!(row.purity, Some(1.0));
        assert_eq!(row.members, None);
        // The span and `copies` come from two different sources — the hull and
        // the segment's own `tract_len` — and nothing in `write_row` ties them.
        // They agree only because the walk builds the hull from the segment; if
        // a locus were ever clipped, the row would claim a span its copy count
        // contradicts, with no panic.
        let RegionKind::SsrSegment(ref emitted) = locus.kind else {
            unreachable!("the fixture is a locus")
        };
        assert_eq!(
            row.end - row.start,
            emitted.tract_len(),
            "the row's span must be the segment's tract, or `copies` describes a \
             different region than the coordinates do"
        );

        // --- SsrBundle: the hull, plus members that must NOT shift (T4) ---
        // Members are 0-based half-open contig coordinates already.
        let members = vec![tract(1990, 2010, 2), tract(2020, 2040, 5)];
        let bundle = TypedRegion {
            // The hull in 1-based inclusive terms: [1991, 2040].
            region: hull(1991, 2040),
            kind: RegionKind::SsrBundle {
                tracts: members.clone().into_boxed_slice(),
            },
        };
        let row = row_of(&bundle, &contigs);
        assert_eq!(
            (row.start, row.end),
            (1990, 2040),
            "the hull shifts by one at the start (T4)"
        );
        assert_eq!(row.kind, "ssr_bundle");
        assert_eq!(
            row.members,
            Some(vec![(1990, 2010, 2), (2020, 2040, 5)]),
            "members are already 0-based half-open — they must come back UNSHIFTED (T4)"
        );
        assert_eq!(
            (row.motif, row.period, row.copies, row.purity),
            (None, None, None, None),
            "a bundle fills only its members"
        );
    }

    /// **The emitted rows are byte-for-byte the spec's worked example** (§3.1).
    ///
    /// The round-trip oracle parses cells back into numbers, so it cannot see how
    /// a value is *rendered* — `"10"` and `"10.0"` both parse to `10.0`. This
    /// pins the rendering itself against the format the spec publishes, which is
    /// what a consumer's `awk`/pandas reads.
    #[test]
    fn the_rows_render_exactly_as_the_spec_example() {
        let contigs = test_contigs();
        let line_of = |region: &TypedRegion| {
            let mut buf = Vec::new();
            write_row(&mut buf, region, &contigs).expect("writing to a Vec");
            String::from_utf8(buf).expect("UTF-8")
        };

        let segment = SsrSegment::new(
            "chr1".into(),
            1000,
            1029,
            Motif::new(b"CAG").expect("valid"),
            1.0,
        )
        .expect("a valid segment");
        assert_eq!(
            line_of(&TypedRegion {
                region: hull(1000, 1029),
                kind: RegionKind::SsrSegment(segment),
            }),
            "chr1\t999\t1029\tssr_locus\tCAG\t3\t10.0\t1.0\t.\n",
            "a whole copy count keeps its decimal point"
        );

        assert_eq!(
            line_of(&TypedRegion {
                region: hull(1, 999),
                kind: RegionKind::Generic,
            }),
            "chr1\t0\t999\tgeneric\t.\t.\t.\t.\t.\n"
        );

        assert_eq!(
            line_of(&TypedRegion {
                region: hull(2041, 4240),
                kind: RegionKind::Satellite,
            }),
            "chr1\t2040\t4240\tsatellite\t.\t.\t.\t.\t.\n"
        );

        assert_eq!(
            line_of(&TypedRegion {
                region: hull(1991, 2040),
                kind: RegionKind::SsrBundle {
                    tracts: vec![tract(1990, 2010, 2), tract(2020, 2040, 5)].into_boxed_slice(),
                },
            }),
            "chr1\t1990\t2040\tssr_bundle\t.\t.\t.\t.\t\
             [{\"start\":1990,\"end\":2010,\"period\":2},\
             {\"start\":2020,\"end\":2040,\"period\":5}]\n"
        );
    }

    /// **The two directions are not the same direction** — the sharpest form of
    /// T4. A bundle whose hull and whose first member describe the *same* first
    /// base must render them as the *same* BED number: the hull via `-1`, the
    /// member via no shift. If either conversion were applied to the other, these
    /// two cells would differ by one.
    #[test]
    fn the_hull_and_its_first_member_agree_on_the_first_base() {
        let contigs = test_contigs();
        // 0-based half-open member [1990, 2010) == 1-based inclusive [1991, 2010].
        let first = tract(1990, 2010, 2);
        let region = TypedRegion {
            region: hull(1991, 2040),
            kind: RegionKind::SsrBundle {
                tracts: vec![first, tract(2020, 2040, 5)].into_boxed_slice(),
            },
        };
        let row = row_of(&region, &contigs);
        let members = row.members.expect("a bundle has members");
        assert_eq!(
            row.start, members[0].0,
            "the hull starts at its first member; a shift applied to the wrong one \
             makes these differ by exactly 1 — silently, and it is a wrong region"
        );
    }

    /// A single-base region survives the conversion: 1-based inclusive `[n, n]`
    /// is 0-based half-open `[n-1, n)`, one base wide, not zero or two.
    #[test]
    fn a_single_base_region_keeps_its_width() {
        let contigs = test_contigs();
        let region = TypedRegion {
            region: hull(1, 1),
            kind: RegionKind::Generic,
        };
        let row = row_of(&region, &contigs);
        assert_eq!((row.start, row.end), (0, 1));
        assert_eq!(row.end - row.start, 1, "one base wide");
    }

    /// The row is named from the run's one contig table, so a region on the
    /// second contig is not silently labelled with the first one's name (T8).
    #[test]
    fn rows_are_named_from_the_shared_contig_table() {
        let contigs = test_contigs();
        let region = TypedRegion {
            region: GenomeRegion {
                contig: ContigId(1),
                start: Position(10),
                end: Position(20),
            },
            kind: RegionKind::Generic,
        };
        assert_eq!(row_of(&region, &contigs).chrom, "chr2");
    }

    /// The `members` cell is deterministic — fixed key order, no incidental
    /// whitespace — which is what §6's byte-identity rests on for a bundle.
    #[test]
    fn the_members_cell_serialises_deterministically() {
        let tracts = [tract(1990, 2010, 2), tract(2020, 2040, 5)];
        assert_eq!(
            members_json(&tracts),
            r#"[{"start":1990,"end":2010,"period":2},{"start":2020,"end":2040,"period":5}]"#
        );
        assert_eq!(members_json(&tracts), members_json(&tracts));
    }

    /// `copies` is fractional in general (T5): a tract that is not a whole number
    /// of motifs reports the fraction, the way trf-mod's `copyNum` does.
    #[test]
    fn copies_is_fractional_when_the_tract_is_not_whole_motifs() {
        let contigs = test_contigs();
        // 1-based inclusive [1000, 1028] = 29 bp of a period-3 motif.
        let segment = SsrSegment::new(
            "chr1".into(),
            1000,
            1028,
            Motif::new(b"CAG").expect("valid"),
            0.95,
        )
        .expect("a valid segment");
        let row = row_of(
            &TypedRegion {
                region: hull(1000, 1028),
                kind: RegionKind::SsrSegment(segment),
            },
            &contigs,
        );
        let copies = row.copies.expect("a locus reports copies");
        assert!(
            (copies - 29.0 / 3.0).abs() < 1e-12,
            "29 bp of a 3 bp motif is {} copies, got {copies}",
            29.0 / 3.0
        );
        assert!(
            copies.fract() > 0.0,
            "the column is fractional, not truncated: {copies}"
        );
    }

    // ---- the fallible setup (spec T1, T2, T8) ---------------------------

    /// A two-contig FASTA on disk, one line per contig, with **no** `.fai`.
    /// Returns the temp dir (kept alive by the caller) and the FASTA path.
    fn write_test_fasta() -> (tempfile::TempDir, PathBuf) {
        use std::io::Write as _;
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta = dir.path().join("ref.fa");
        let mut file = std::fs::File::create(&fasta).expect("create the FASTA");
        writeln!(file, ">chr1").expect("write");
        writeln!(file, "{}", "ACGTACGTAC".repeat(10)).expect("write"); // 100 bp
        writeln!(file, ">chr2").expect("write");
        writeln!(file, "{}", "TTGGCCAATT".repeat(5)).expect("write"); // 50 bp
        (dir, fasta)
    }

    fn args_for(reference: &Path, regions: Option<PathBuf>) -> TypedRegionsArgs {
        args_for_output(reference, regions, Path::new("out.tsv"))
    }

    /// Args built by **parsing them**, so every test runs against the same
    /// defaults resolution a user gets, not a hand-filled struct.
    fn args_for_output(
        reference: &Path,
        regions: Option<PathBuf>,
        output: &Path,
    ) -> TypedRegionsArgs {
        let mut argv = vec![
            "pop_var_caller_exp".to_string(),
            "type-regions".to_string(),
            "--reference".to_string(),
            reference.display().to_string(),
            "--output".to_string(),
            output.display().to_string(),
        ];
        if let Some(bed) = regions {
            argv.push("--regions".to_string());
            argv.push(bed.display().to_string());
        }
        let cli = Cli::try_parse_from(argv).expect("the args parse");
        let PopVarCallerExpCommand::TypeRegions(args) = cli.cmd;
        args
    }

    /// **T2 — the contig-length narrowing fails rather than folds.** A contig
    /// past `u32::MAX` is unrepresentable to a `RegionSet` whatever we do, so it
    /// is an error naming the contig, not an `as` cast that silently reports a
    /// wrong length.
    #[test]
    fn a_contig_longer_than_four_gb_is_rejected_by_name() {
        let contigs = ContigList {
            entries: vec![
                ContigEntry {
                    name: "small".to_string(),
                    length: 1_000,
                    md5: None,
                },
                ContigEntry {
                    name: "huge".to_string(),
                    length: u64::from(u32::MAX) + 1,
                    md5: None,
                },
            ],
        };
        let err = contig_bounds(&contigs).expect_err("a >4 Gb contig must be refused");
        match err {
            TypedRegionsCliError::ContigTooLong { contig, len } => {
                assert_eq!(contig, "huge", "the error names the offending contig");
                assert_eq!(len, u64::from(u32::MAX) + 1, "and its real length");
            }
            other => panic!("expected ContigTooLong, got {other:?}"),
        }
    }

    /// The boundary is exact: `u32::MAX` fits, one more does not.
    #[test]
    fn a_contig_of_exactly_u32_max_still_fits() {
        let contigs = ContigList {
            entries: vec![ContigEntry {
                name: "edge".to_string(),
                length: u64::from(u32::MAX),
                md5: None,
            }],
        };
        let bounds = contig_bounds(&contigs).expect("u32::MAX is representable");
        assert_eq!(bounds[0].length, u32::MAX);
    }

    /// The whole-genome path builds one whole-contig span per contig.
    ///
    /// Asserts the **extents**, not merely which contigs appear: membership
    /// alone would pass if `whole_contigs` produced 1-base spans, spans of the
    /// wrong length, or spans in the wrong order. The fixture knows the answer
    /// exactly (chr1 = 100 bp, chr2 = 50 bp), so it is asserted.
    #[test]
    fn the_whole_genome_path_covers_each_contig_end_to_end() {
        let (_dir, fasta) = write_test_fasta();
        let cache = Arc::new(ReferenceInfoCache::new());
        let inputs =
            prepare_walk_inputs(&args_for(&fasta, None), &cache).expect("the setup succeeds");

        assert_eq!(
            inputs.contigs.entries.len(),
            2,
            "both contigs are in the table"
        );
        let spans: Vec<(u32, u64, u64)> = inputs
            .regions
            .iter()
            .map(|r| (r.contig.get(), r.start.get(), r.end.get()))
            .collect();
        assert_eq!(
            spans,
            vec![(0, 1, 100), (1, 1, 50)],
            "one whole-contig span each, 1-based inclusive, in contig order"
        );
    }

    /// The `--regions` path builds, and resolves BED names against **this**
    /// contig table (T8) — a span named `chr2` comes back on the contig the
    /// table calls `chr2`, not on index 0.
    #[test]
    fn a_bed_builds_regions_resolved_against_the_same_contig_table() {
        use std::io::Write as _;
        let (dir, fasta) = write_test_fasta();
        let bed = dir.path().join("want.bed");
        let mut file = std::fs::File::create(&bed).expect("create the BED");
        writeln!(file, "chr2\t10\t40").expect("write");
        drop(file);

        let cache = Arc::new(ReferenceInfoCache::new());
        let inputs =
            prepare_walk_inputs(&args_for(&fasta, Some(bed)), &cache).expect("the setup succeeds");

        let regions: Vec<_> = inputs.regions.iter().collect();
        assert_eq!(regions.len(), 1, "one requested span: {regions:?}");
        let contig_id = regions[0].contig.get() as usize;
        assert_eq!(
            inputs.contigs.entries[contig_id].name, "chr2",
            "the BED's chrom resolved through the SAME table the walk will use (T8)"
        );
    }

    /// **T8 — one table, used twice.** Every `ContigId` the region set carries
    /// indexes the contig table the caller is handed, so the list that resolved
    /// the regions is the same one that will open the reference. If they could
    /// differ, the walk's failure would be `UnknownContig(ContigId)` — an index,
    /// not a name.
    /// **Tied by content, not by count.** A bare `id < entries.len()` bounds
    /// check is true *by construction* — `whole_contigs` derives each id from
    /// the slice index of whatever table it was handed — so it could never fail
    /// and would still pass against a permuted table or a different reference's
    /// table of the same size. Asserting each span's extent against *that entry's
    /// length* is what a mismatched table actually breaks (chr1 100 vs chr2 50).
    ///
    /// The cross-module half of T8 — that E1 hands *this* list to
    /// `WindowedRefSeq` — is not observable from inside `prepare_walk_inputs`;
    /// this guards the in-function half.
    #[test]
    fn each_region_matches_the_returned_tables_entry() {
        let (_dir, fasta) = write_test_fasta();
        let cache = Arc::new(ReferenceInfoCache::new());
        let inputs =
            prepare_walk_inputs(&args_for(&fasta, None), &cache).expect("the setup succeeds");

        let regions: Vec<_> = inputs.regions.iter().collect();
        assert_eq!(
            regions.len(),
            inputs.contigs.entries.len(),
            "one whole-contig span per table entry"
        );
        for region in &regions {
            let entry = inputs
                .contigs
                .entries
                .get(region.contig.get() as usize)
                .expect("every ContigId indexes the returned table");
            assert_eq!(
                region.end.get(),
                entry.length,
                "the span on '{}' must run to that entry's own length — a permuted \
                 or foreign table breaks here, where a bounds check would not",
                entry.name
            );
        }
    }

    /// A BED naming a contig the reference does not have is a `Bed` error — the
    /// one error arm of the setup that is neither the reference nor the
    /// narrowing. `RegionSet` rejects it up front, which is what lets the walk
    /// itself have so few failure modes.
    ///
    /// **Known noise on this path** (see `prepare_walk_inputs`): the background
    /// verification handle is dropped un-joined here, so `VerificationHandle`'s
    /// `Drop` prints its "dropped without join()" warning alongside the real
    /// error. Nothing was published, so the abandonment is legitimate — the
    /// warning is misleading rather than wrong.
    #[test]
    fn a_bed_naming_an_unknown_contig_is_a_bed_error() {
        use std::io::Write as _;
        let (dir, fasta) = write_test_fasta();
        let bed = dir.path().join("unknown.bed");
        let mut file = std::fs::File::create(&bed).expect("create the BED");
        writeln!(file, "chrZ\t10\t40").expect("write");
        drop(file);

        let cache = Arc::new(ReferenceInfoCache::new());
        // Matched rather than `expect_err`, which would need `Debug` on
        // `WalkInputs` — and that would mean deriving it through a live thread
        // handle for a test message.
        match prepare_walk_inputs(&args_for(&fasta, Some(bed)), &cache) {
            Err(TypedRegionsCliError::Bed(_)) => {}
            Err(other) => panic!("expected a Bed error, got {other:?}"),
            Ok(_) => panic!("a BED naming an absent contig must be refused"),
        }
    }

    /// **T1, the `.fai`-absent case** — the reference is scanned once and the
    /// sibling `.fai` is **written**, so the next run takes the fast path. No
    /// handle comes back: there is nothing left to verify.
    #[test]
    fn an_absent_fai_is_created_beside_the_reference() {
        let (_dir, fasta) = write_test_fasta();
        // Ask the module for the convention rather than re-deriving it: a
        // hand-built `with_extension` only happens to agree while the fixture is
        // named `.fa`, and would silently check an unrelated path if renamed.
        let fai = crate::ng::reference_info::sibling_fai_path(&fasta);
        assert!(!fai.exists(), "the fixture starts without a .fai");

        let cache = Arc::new(ReferenceInfoCache::new());
        let inputs =
            prepare_walk_inputs(&args_for(&fasta, None), &cache).expect("the setup succeeds");

        assert!(fai.exists(), "the sibling .fai must be written: {fai:?}");
        assert!(
            inputs.verify.is_none(),
            "a FASTA scanned in the foreground is already verified — nothing to join"
        );
    }

    /// **T1, the `.fai`-present case** — the index is read in the foreground and
    /// the FASTA is verified on a background thread, so a handle comes back for
    /// the caller to join before publishing (§6).
    #[test]
    fn a_present_fai_returns_a_verification_handle() {
        let (_dir, fasta) = write_test_fasta();
        let cache = Arc::new(ReferenceInfoCache::new());
        // First run writes the .fai.
        prepare_walk_inputs(&args_for(&fasta, None), &cache).expect("the first run succeeds");

        // Second run finds it, and verifies in the background.
        let inputs =
            prepare_walk_inputs(&args_for(&fasta, None), &cache).expect("the second run succeeds");
        let handle = inputs
            .verify
            .expect("a .fai read in the foreground leaves a verification pending");
        handle
            .join()
            .expect("the reference matches its own fresh .fai");
    }

    // ---- E2: the file, end to end (spec §7, arch §9) --------------------

    /// A two-contig reference with deliberate structure: aperiodic filler, a
    /// clean `(CAG)*10` locus, two tracts closer than the flank radius (a
    /// **bundle**, without which the member coordinates are never exercised),
    /// and a 300 bp array past the satellite cap. Verified to produce all four
    /// kinds.
    fn e2e_reference() -> (tempfile::TempDir, PathBuf) {
        use std::io::Write as _;
        // 16 bp, aperiodic, and free of any homopolymer run >= 6 — so the filler
        // itself is never classified at `--min-period 1`.
        const FILLER: &str = "ACGTTGCAAGCTTGCA";

        let mut ctg1 = String::new();
        ctg1.push_str(&FILLER.repeat(13));
        ctg1.push_str(&"CAG".repeat(10)); // a locus
        ctg1.push_str(&FILLER.repeat(13));
        ctg1.push_str(&"AT".repeat(10)); // ┐ 16 bp apart, inside the 30 bp
        ctg1.push_str(FILLER); // │ flank radius, so they
        ctg1.push_str(&"GACA".repeat(6)); // ┘ bundle rather than classify
        ctg1.push_str(&FILLER.repeat(13));
        ctg1.push_str(&"AT".repeat(150)); // 300 bp: past the 100 bp cap
        ctg1.push_str(&FILLER.repeat(13));
        let ctg2 = FILLER.repeat(20);

        let dir = tempfile::tempdir().expect("tempdir");
        let fasta = dir.path().join("ref.fa");
        let mut file = std::fs::File::create(&fasta).expect("create the FASTA");
        writeln!(file, ">ctg1\n{ctg1}\n>ctg2\n{ctg2}").expect("write");
        (dir, fasta)
    }

    /// The rows of a written partition, header lines skipped.
    fn read_rows(path: &Path) -> Vec<ParsedRow> {
        std::fs::read_to_string(path)
            .expect("the output exists")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(parse_row)
            .collect()
    }

    /// Drive the walk directly, for the round-trip's other side.
    fn walk_directly(args: &TypedRegionsArgs) -> (ContigList, Vec<TypedRegion>) {
        let cache = Arc::new(ReferenceInfoCache::new());
        let inputs = prepare_walk_inputs(args, &cache).expect("setup");
        let config = walk_config(args).expect("config");
        let reference = WindowedRefSeq::new(args.reference.clone(), inputs.contigs.clone());
        let walk = TypedRegionIterator::over_regions(reference, inputs.regions, config)
            .expect("the walk starts");
        let regions = walk.map(|r| r.expect("no read fails")).collect();
        if let Some(handle) = inputs.verify {
            handle.join().expect("the fixture's .fai is fresh");
        }
        (inputs.contigs, regions)
    }

    /// **THE ROUND-TRIP, on a real reference through the whole stack** (plan E2,
    /// spec §7): the rows written to the file parse back equal, field for field,
    /// to what the iterator itself emitted — including a bundle's members.
    ///
    /// C2's oracle pinned the formatter against hand-derived expectations; this
    /// pins the *file* against the *walk*, so a dropped row, a header line
    /// miscounted as a row, or a truncated write is caught too.
    #[test]
    fn the_written_file_round_trips_against_the_iterators_own_output() {
        let (dir, fasta) = e2e_reference();
        let output = dir.path().join("out.tsv");
        let args = args_for_output(&fasta, None, &output);

        run_typed_regions(&args).expect("the run succeeds");
        let rows = read_rows(&output);
        let (contigs, regions) = walk_directly(&args);

        assert_eq!(rows.len(), regions.len(), "one row per emitted region");
        assert!(!regions.is_empty(), "the fixture produces regions");
        for (row, region) in rows.iter().zip(&regions) {
            assert_eq!(row, &row_of(region, &contigs), "row vs the walk's region");
            // T4's conversion re-derived INLINE, not by calling `hull_to_bed` —
            // asserting against the helper the row was built with would be
            // `f(x) == f(x)`, green through any off-by-one edit inside it.
            assert_eq!(
                (row.start, row.end),
                (region.region.start.get() - 1, region.region.end.get()),
                "the hull's start loses one, its end does not move"
            );
        }

        // The fixture must exercise every kind, or the round-trip is weaker than
        // it looks — a bundle above all (spec §7: "a fixture with a bundle is
        // required, or the member coordinates are never exercised").
        let kinds: std::collections::BTreeSet<&str> =
            rows.iter().map(|r| r.kind.as_str()).collect();
        for kind in ["generic", "ssr_locus", "ssr_bundle", "satellite"] {
            assert!(
                kinds.contains(kind),
                "the fixture must produce a {kind} row"
            );
        }
        let bundle = rows
            .iter()
            .find(|r| r.kind == "ssr_bundle")
            .expect("a bundle row");
        assert!(
            bundle.members.as_ref().is_some_and(|m| m.len() >= 2),
            "a bundle carries >= 2 members: {bundle:?}"
        );
    }

    /// **The partition invariant, in the file** (spec §7): the rows' spans are
    /// contiguous, non-overlapping, and reconstruct each requested contig end to
    /// end. This is the test that catches a dropped `generic` row, a header line
    /// miscounted as a row, or an off-by-one at T4.
    #[test]
    fn the_files_spans_reconstruct_the_requested_regions() {
        let (dir, fasta) = e2e_reference();
        let output = dir.path().join("out.tsv");
        let args = args_for_output(&fasta, None, &output);
        run_typed_regions(&args).expect("the run succeeds");

        let (contigs, _) = walk_directly(&args);
        let rows = read_rows(&output);

        for entry in &contigs.entries {
            let spans: Vec<(u64, u64)> = rows
                .iter()
                .filter(|r| r.chrom == entry.name)
                .map(|r| (r.start, r.end))
                .collect();
            assert!(!spans.is_empty(), "contig {} has rows", entry.name);

            let mut expected_start = 0;
            for (start, end) in &spans {
                assert_eq!(
                    *start, expected_start,
                    "contig {}: a gap or overlap at {start}",
                    entry.name
                );
                assert!(end > start, "contig {}: empty span at {start}", entry.name);
                expected_start = *end;
            }
            assert_eq!(
                expected_start, entry.length,
                "contig {} must be covered end to end",
                entry.name
            );
        }
    }

    /// **Determinism** (spec §6): the same reference, regions and config give a
    /// byte-identical file. The header carries no timestamp and no digest
    /// precisely so this holds.
    #[test]
    fn two_runs_are_byte_identical() {
        let (dir, fasta) = e2e_reference();
        let first = dir.path().join("first.tsv");
        let second = dir.path().join("second.tsv");

        run_typed_regions(&args_for_output(&fasta, None, &first)).expect("run 1");
        run_typed_regions(&args_for_output(&fasta, None, &second)).expect("run 2");

        assert_eq!(
            std::fs::read(&first).expect("first"),
            std::fs::read(&second).expect("second"),
            "the walk is a pure function of reference + regions + config"
        );
    }

    /// **`--regions` read-back** (spec §3.2): the file is a BED our own
    /// `--regions` accepts, so a selection cut from it walks exactly those spans.
    /// Otherwise §3.2's claim is an argument rather than a property.
    ///
    /// **Feeding the file back WHOLE would prove nothing.** Its rows abut by
    /// construction (that is the partition invariant), and `RegionSet` coalesces
    /// abutting spans — so the region set would collapse to one whole-contig
    /// span per contig, i.e. exactly the whole-genome path, and the test would
    /// reduce to determinism plus "the BED parsed". So the selection here is
    /// **non-adjacent** rows (`grep`-ing out the interesting kinds, which is the
    /// workflow §3.2 exists for), which is the only shape that exercises the
    /// multi-span path.
    #[test]
    fn a_selection_cut_from_the_output_walks_exactly_those_spans() {
        use std::io::Write as _;
        let (dir, fasta) = e2e_reference();
        let first = dir.path().join("first.tsv");
        run_typed_regions(&args_for_output(&fasta, None, &first)).expect("the first run");

        // The `grep -P '\tsatellite\t'` workflow: pick the non-generic findings,
        // which the fixture places far apart, so they survive coalescing.
        let wanted: Vec<(String, u64, u64)> = read_rows(&first)
            .into_iter()
            .filter(|r| r.kind == "ssr_locus" || r.kind == "satellite")
            .map(|r| (r.chrom, r.start, r.end))
            .collect();
        assert!(
            wanted.len() >= 2,
            "the selection must be more than one span, or coalescing hides the point"
        );

        let bed = dir.path().join("selection.bed");
        let mut file = std::fs::File::create(&bed).expect("create the BED");
        for (chrom, start, end) in &wanted {
            writeln!(file, "{chrom}\t{start}\t{end}").expect("write");
        }
        drop(file);

        let second = dir.path().join("second.tsv");
        run_typed_regions(&args_for_output(&fasta, Some(bed), &second))
            .expect("a selection cut from our own output is a valid --regions BED");

        let got: Vec<(String, u64, u64)> = read_rows(&second)
            .into_iter()
            .map(|r| (r.chrom, r.start, r.end))
            .collect();
        assert_eq!(
            got, wanted,
            "the emitted spans are exactly the requested ones — each finding \
             returned whole, and nothing outside the selection"
        );
    }

    /// **The knobs whose library guard is a release `assert!` are bounded at the
    /// CLI, so a typo is a usage error and never a backtrace** (spec §6).
    ///
    /// These bounds are the only thing standing between a user and a panic —
    /// without them `--max-period 7` trips `classify`'s period-ceiling assert,
    /// `--flank-bp 0` trips its flank assert on the *first window of any
    /// reference*, and `--min-purity nan` trips its purity assert (and a NaN
    /// floor would pass every tract rather than reject impure ones). Deleting
    /// any of the three `value_parser`s must fail here.
    #[test]
    fn knobs_guarded_by_library_asserts_are_usage_errors_at_the_cli() {
        // The `=` form throughout, so clap reads every case as a VALUE — bare
        // `--min-purity -0.5` is rejected too, but as an unknown argument (`-0`
        // parses as a flag), which would be testing clap's tokeniser rather than
        // these bounds.
        for bad in [
            "--max-period=7",    // above MAX_MOTIF_LEN
            "--min-period=0",    // period 0 is meaningless
            "--flank-bp=0",      // every locus would fail its own flank test
            "--min-purity=1.5",  // outside [0, 1]
            "--min-purity=-0.5", //
            "--min-purity=nan",  // parses as f32, and would pass every tract
            "--min-purity=inf",  //
        ] {
            let err = Cli::try_parse_from([
                "pop_var_caller_exp",
                "type-regions",
                "--reference",
                "r.fa",
                "--output",
                "o.tsv",
                bad,
            ])
            .expect_err("'{bad}' must be refused at parse time");
            assert_eq!(
                err.kind(),
                clap::error::ErrorKind::ValueValidation,
                "'{bad}' must be a usage failure, not a run panic: {err}"
            );
        }
    }

    /// The cross-flag period rule clap cannot express: `--min-period` above
    /// `--max-period` is an error from the driver, not the panic
    /// `PeriodRange::new`'s caller would otherwise take.
    #[test]
    fn a_min_period_above_max_period_is_an_error_not_a_panic() {
        let (dir, fasta) = e2e_reference();
        let output = dir.path().join("out.tsv");
        let mut args = args_for_output(&fasta, None, &output);
        args.min_period = 5;
        args.max_period = 2;

        match run_typed_regions(&args) {
            Err(TypedRegionsCliError::PeriodRange(PeriodRangeError::MinExceedsMax {
                min,
                max,
            })) => assert_eq!((min, max), (5, 2), "it names both numbers"),
            Err(other) => panic!("expected a PeriodRange error, got {other:?}"),
            Ok(()) => panic!("an empty period range must be refused"),
        }
        assert!(!output.exists(), "nothing is published");
    }

    /// **T3 — the flag pair is a CLI error, not a panic.** `--max-str-len` below
    /// `--flank-bp` is refused by the walk before any work, and the CLI just
    /// propagates it: a second check here would be a second place for the rule to
    /// drift, and the walk's message already carries both numbers.
    #[test]
    fn a_margin_narrower_than_the_flank_is_an_error_not_a_panic() {
        let (dir, fasta) = e2e_reference();
        let output = dir.path().join("out.tsv");
        let mut args = args_for_output(&fasta, None, &output);
        args.max_str_len = 10;
        args.flank_bp = 30;

        match run_typed_regions(&args) {
            Err(TypedRegionsCliError::Walk(TypedRegionError::MarginNarrowerThanFlank {
                max_str_len,
                flank_bp,
            })) => {
                assert_eq!((max_str_len, flank_bp), (10, 30), "it names both numbers");
            }
            Err(other) => panic!("expected the walk's flank-pair error, got {other:?}"),
            Ok(()) => panic!("a margin narrower than the flank must be refused"),
        }
        assert!(
            !output.exists(),
            "nothing is published when the walk refuses to start"
        );
    }

    /// **T1 — a stale `.fai` aborts AT THE JOIN BARRIER, with nothing
    /// published.** The verification runs beside the walk and is joined at the
    /// commit point, so a reference that no longer matches its index fails with
    /// the partition still under its temp name. Joining *after* the rename would
    /// publish a partition built against a wrong contig table — silently wrong,
    /// which is the whole reason the barrier sits where it does.
    ///
    /// **The staleness is chosen so the walk cannot catch it first.** The
    /// reference is *extended*, so every span the stale table describes still
    /// reads cleanly and the walk runs to completion — only the background check,
    /// which recomputes the whole index, can tell. (A *shrunken* reference aborts
    /// too, but through the walk's own read: see the companion test, which exists
    /// so this one is not credited with catching that case.)
    #[test]
    fn a_stale_fai_aborts_at_the_join_with_nothing_published() {
        use std::io::Write as _;
        let (dir, fasta) = e2e_reference();

        // A first run creates the sibling `.fai`.
        let warm = dir.path().join("warm.tsv");
        run_typed_regions(&args_for_output(&fasta, None, &warm)).expect("the first run");

        // Append to the last contig: the stale table still describes readable
        // spans, so the walk succeeds — but the file no longer matches the index.
        let mut file = std::fs::OpenOptions::new()
            .append(true)
            .open(&fasta)
            .expect("open the FASTA for append");
        writeln!(file, "ACGTACGTAC").expect("write");
        drop(file);

        let output = dir.path().join("stale.tsv");
        let err = run_typed_regions(&args_for_output(&fasta, None, &output))
            .expect_err("a stale .fai must abort the run");
        assert!(
            matches!(err, TypedRegionsCliError::Reference(_)),
            "the barrier, not the walk, must catch this: {err:?}"
        );
        assert!(
            !output.exists(),
            "the partition must NOT be published — the join precedes the rename"
        );
    }

    /// A reference that has *shrunk* under its index also aborts before
    /// publishing, but through the walk's own read (`UnexpectedEof`) rather than
    /// the verification barrier.
    ///
    /// Recorded as its own case so the barrier test above cannot be credited
    /// with catching it: the two failures travel different paths, and only the
    /// extended-reference case actually exercises the join.
    #[test]
    fn a_truncated_reference_also_aborts_before_publishing() {
        use std::io::Write as _;
        let (dir, fasta) = e2e_reference();
        let warm = dir.path().join("warm.tsv");
        run_typed_regions(&args_for_output(&fasta, None, &warm)).expect("the first run");

        let mut file = std::fs::File::create(&fasta).expect("truncate the FASTA");
        writeln!(file, ">ctg1\nACGTACGTACGTACGTACGT\n>ctg2\nACGTACGTACGT").expect("write");
        drop(file);

        let output = dir.path().join("short.tsv");
        let err = run_typed_regions(&args_for_output(&fasta, None, &output))
            .expect_err("a reference shorter than its index must abort");
        assert!(
            matches!(err, TypedRegionsCliError::Walk(_)),
            "the walk's own read catches this one: {err:?}"
        );
        assert!(!output.exists(), "and still nothing is published");
    }

    /// `--min-copies` is accepted when it carries exactly six values.
    #[test]
    fn min_copies_accepts_exactly_six_values() {
        let cli = Cli::try_parse_from([
            "pop_var_caller_exp",
            "type-regions",
            "--reference",
            "r.fa",
            "--output",
            "o.tsv",
            "--min-copies",
            "9,5,4,3,3,3",
        ])
        .expect("six values parse");
        let PopVarCallerExpCommand::TypeRegions(args) = cli.cmd;
        assert_eq!(args.min_copies.for_period(1), 9);
        assert_eq!(args.min_copies.for_period(2), 5);
    }

    /// **A wrong count is a clap usage failure**, not a `run` error — it is
    /// rejected during parsing, so `run_typed_regions` is never reached
    /// (spec §2.1).
    #[test]
    fn min_copies_with_a_wrong_count_is_a_cli_usage_error() {
        for bad in ["6,4,4,3,3", "6,4,4,3,3,3,3", "6"] {
            let err = Cli::try_parse_from([
                "pop_var_caller_exp",
                "type-regions",
                "--reference",
                "r.fa",
                "--output",
                "o.tsv",
                "--min-copies",
                bad,
            ])
            .expect_err("a wrong count must be rejected at parse time");
            assert_eq!(
                err.kind(),
                clap::error::ErrorKind::ValueValidation,
                "'{bad}' must fail as a value-validation usage error: {err}"
            );
        }
    }
}
