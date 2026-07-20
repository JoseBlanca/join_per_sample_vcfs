//! The `type-regions` subcommand: run step 3's walk over a reference and
//! write the typed-region partition to a file. This module owns the `Args`
//! struct, the `run_typed_regions` driver, and its `#[non_exhaustive]`
//! error enum.
//!
//! **In progress (Milestone B complete).** The full knob surface — including
//! the `--min-copies` table and its
//! [parser](crate::pop_var_caller_exp::cli::parsers::parse_min_copies) — and the
//! real error enum are here; `run_typed_regions` is wired in Milestone E.
//! See `doc/devel/ng/impl_plan/typed_regions_cli.md`.

use std::io;
use std::path::{Path, PathBuf};

use clap::Args;
use thiserror::Error;

use crate::fasta::ContigList;
use crate::ng::reference_info::ReferenceInfoError;
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_FLANK_BP, DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY,
    DEFAULT_MIN_SCORE, MAX_MOTIF_LEN, MinCopies, SsrSegmentCriteria,
};
use crate::ng::region_typing::{
    DEFAULT_MAX_STR_LEN, DEFAULT_WINDOW_BP, RegionKind, TypedRegion, TypedRegionConfig,
    TypedRegionError,
};
use crate::ng::tandem_repeat::{
    DEFAULT_MATCH_REWARD, DEFAULT_MIN_COPIES, DEFAULT_MISMATCH_PENALTY, RepeatInterval, ScanParams,
};
use crate::ng::types::{ContigId, GenomeRegion};
use crate::regions::BedError;

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
    #[arg(long, default_value_t = DEFAULT_MIN_PERIOD, help_heading = "Advanced")]
    pub min_period: u8,

    /// Widest STR period classified — the microsatellite ceiling.
    #[arg(long, default_value_t = DEFAULT_MAX_PERIOD, help_heading = "Advanced")]
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
    #[arg(long, default_value_t = DEFAULT_FLANK_BP, help_heading = "Advanced")]
    pub flank_bp: u64,

    /// Purity floor in `[0, 1]`: a tract matching less than this fraction of a
    /// perfect motif tiling is not classified.
    #[arg(long, default_value_t = DEFAULT_MIN_PURITY, help_heading = "Advanced")]
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

/// Run the walk and write the partition. Wired in Milestone E1; Milestones B–D
/// build its inputs (the args surface, the output writer, the reference setup).
pub fn run_typed_regions(_args: &TypedRegionsArgs) -> Result<(), TypedRegionsCliError> {
    todo!("run_typed_regions is implemented in Milestone E1")
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
