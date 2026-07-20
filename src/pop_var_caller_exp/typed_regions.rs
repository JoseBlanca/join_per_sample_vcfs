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

use crate::ng::reference_info::ReferenceInfoError;
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_FLANK_BP, DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY,
    DEFAULT_MIN_SCORE, MAX_MOTIF_LEN, MinCopies, SsrSegmentCriteria,
};
use crate::ng::region_typing::{
    DEFAULT_MAX_STR_LEN, DEFAULT_WINDOW_BP, TypedRegionConfig, TypedRegionError,
};
use crate::ng::tandem_repeat::{
    DEFAULT_MATCH_REWARD, DEFAULT_MIN_COPIES, DEFAULT_MISMATCH_PENALTY, ScanParams,
};
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
