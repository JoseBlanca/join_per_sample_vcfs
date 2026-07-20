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

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::ng::reference_info::ReferenceInfoError;
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_FLANK_BP, DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY,
    DEFAULT_MIN_SCORE, MinCopies,
};
use crate::ng::region_typing::{DEFAULT_MAX_STR_LEN, DEFAULT_WINDOW_BP, TypedRegionError};
use crate::ng::tandem_repeat::{
    DEFAULT_MATCH_REWARD, DEFAULT_MIN_COPIES, DEFAULT_MISMATCH_PENALTY,
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
