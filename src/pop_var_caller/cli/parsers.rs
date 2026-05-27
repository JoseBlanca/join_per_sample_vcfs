//! `clap` `value_parser` helpers for the cohort CLI subcommands.
//!
//! Each function is the CLI-side gatekeeper for one knob, mirroring
//! the range constants enforced engine-side in
//! [`crate::var_calling`]'s `Config::new` / `Config::validate`
//! constructors (Task 1). Two-layer validation — defence in depth:
//!
//! * **CLI side (this module):** fast, human-readable error rendered
//!   by clap before any `.psp` is opened. Mirrors
//!   [`crate::pop_var_caller::cli::parse_mismatch_fraction`].
//! * **Engine side:** typed `ConfigError` variants produced by each
//!   stage's `Config::new` / `Config::validate`. Library users who
//!   skip the CLI still get the same guard.
//!
//! Each parser returns `Result<T, String>` because clap expects that
//! shape from `value_parser`. The error string surfaces verbatim to
//! the user.
//!
//! Convention: function names match the CLI flag they back —
//! `parse_ploidy` parses `--ploidy`, `parse_max_gq_phred` parses
//! `--max-gq-phred`, etc.

use crate::psp::writer::{MAX_BLOCK_TARGET_BYTES, MIN_BLOCK_TARGET_BYTES};
use crate::var_calling::contamination_estimation::{
    C_S_INIT_RANGE_MAX, MIN_COHORT_MINOR_FRACTION_RANGE_MAX_EXCLUSIVE,
    MIN_MAJOR_FRACTION_RANGE_MIN_EXCLUSIVE, STABILITY_TOLERANCE_RANGE_MAX,
};
use crate::var_calling::dust_filter::{MAX_DUST_WINDOW, SD_WLEN};
use crate::var_calling::per_group_merger::{
    MAX_ALLELES_PER_VAR_CAP, MAX_BITMASK_ALLELES, MAX_PLOIDY,
};
use crate::var_calling::posterior_engine::{
    CONVERGENCE_THRESHOLD_RANGE_MAX, GQ_PHRED_RANGE_MAX, GQ_PHRED_RANGE_MIN_EXCLUSIVE,
    MAX_ITERATIONS_RANGE_MAX, PSEUDOCOUNT_RANGE_MAX,
};

// ---------------------------------------------------------------------
// Generic helpers — call sites pass a closure for the range predicate.
// ---------------------------------------------------------------------

/// Parse an `f64`, reject non-finite values, then apply `check`.
/// `name` is the flag name shown in the error, `range` is a
/// human-readable range description (e.g. `"(0.0, 0.1]"`).
fn parse_f64_with(
    s: &str,
    name: &str,
    check: impl Fn(f64) -> bool,
    range: &str,
) -> Result<f64, String> {
    let v: f64 = s
        .parse()
        .map_err(|e| format!("{name}: not a number ({e})"))?;
    if !v.is_finite() {
        return Err(format!("{name} must be finite, got `{s}`"));
    }
    if !check(v) {
        return Err(format!("{name} must be in {range}, got `{s}`"));
    }
    Ok(v)
}

/// Parse a `u32`, then check `range.contains(value)`.
fn parse_u32_in(s: &str, name: &str, lo: u32, hi: u32) -> Result<u32, String> {
    let v: u32 = s
        .parse()
        .map_err(|e| format!("{name}: not a non-negative integer ({e})"))?;
    if !(lo..=hi).contains(&v) {
        return Err(format!("{name} must be in {lo}..={hi}, got `{s}`"));
    }
    Ok(v)
}

/// Parse a `u32`, then check `v >= lo`. Use this when the only
/// constraint is a lower bound; pairing `parse_u32_in` with `u32::MAX`
/// as the upper bound advertises a policy that isn't there.
fn parse_u32_min(s: &str, name: &str, lo: u32) -> Result<u32, String> {
    let v: u32 = s
        .parse()
        .map_err(|e| format!("{name}: not a non-negative integer ({e})"))?;
    if v < lo {
        return Err(format!("{name} must be >= {lo}, got `{s}`"));
    }
    Ok(v)
}

/// Parse a `u32` with no range constraint — use sparingly, when every
/// non-negative integer is genuinely admissible (e.g. DUST threshold,
/// where large values just mean "don't mask anything").
fn parse_u32_any(s: &str, name: &str) -> Result<u32, String> {
    s.parse()
        .map_err(|e| format!("{name}: not a non-negative integer ({e})"))
}

/// Parse a `u8`, then check `lo..=hi`.
fn parse_u8_in(s: &str, name: &str, lo: u8, hi: u8) -> Result<u8, String> {
    let v: u8 = s
        .parse()
        .map_err(|e| format!("{name}: not an integer in 0..=255 ({e})"))?;
    if !(lo..=hi).contains(&v) {
        return Err(format!("{name} must be in {lo}..={hi}, got `{s}`"));
    }
    Ok(v)
}

/// Parse a `usize`, then check `lo..=hi`.
fn parse_usize_in(s: &str, name: &str, lo: usize, hi: usize) -> Result<usize, String> {
    let v: usize = s
        .parse()
        .map_err(|e| format!("{name}: not a non-negative integer ({e})"))?;
    if !(lo..=hi).contains(&v) {
        return Err(format!("{name} must be in {lo}..={hi}, got `{s}`"));
    }
    Ok(v)
}

// ---------------------------------------------------------------------
// Per-stage parsers — ranges match the engine-side validation table
// in `doc/devel/implementation_plans/pop_var_caller_cohort_cli.md`.
// ---------------------------------------------------------------------

// ---- Stage 5 (per-group merger) ---------------------------------

/// `--ploidy`: `1..=MAX_PLOIDY`.
pub fn parse_ploidy(s: &str) -> Result<u8, String> {
    parse_u8_in(s, "ploidy", 1, MAX_PLOIDY)
}

/// `--max-alleles-per-var`: `2..=MAX_ALLELES_PER_VAR_CAP`.
pub fn parse_max_alleles(s: &str) -> Result<usize, String> {
    parse_usize_in(s, "max-alleles-per-var", 2, MAX_ALLELES_PER_VAR_CAP)
}

/// `--max-alleles-lh-calc`: `2..=MAX_BITMASK_ALLELES`. Absolute
/// ceiling on the unified allele set; backstops the likelihood
/// routine's fixed-width bitmask + stack scratch.
pub fn parse_max_alleles_lh_calc(s: &str) -> Result<usize, String> {
    parse_usize_in(s, "max-alleles-lh-calc", 2, MAX_BITMASK_ALLELES)
}

// ---- Stage 4 (grouper) ------------------------------------------

/// `--var-group-max-span`: `>= 1`. No engine-side upper bound; a
/// value of `0` makes every non-empty group exceed the cap
/// (degenerate).
pub fn parse_var_group_max_span(s: &str) -> Result<u32, String> {
    parse_u32_min(s, "var-group-max-span", 1)
}

// ---- Stage 3 (DUST filter) --------------------------------------

/// `--complexity-window`: `SD_WLEN..=MAX_DUST_WINDOW` (`3..=65535`).
pub fn parse_dust_window(s: &str) -> Result<u32, String> {
    parse_u32_in(s, "complexity-window", SD_WLEN, MAX_DUST_WINDOW)
}

/// `--complexity-threshold`: any `u32`. The DUST scorer treats large
/// thresholds as "don't mask anything"; no engine-side upper bound
/// applies.
pub fn parse_dust_threshold(s: &str) -> Result<u32, String> {
    parse_u32_any(s, "complexity-threshold")
}

// ---- Stage 6 (posterior engine) ---------------------------------

/// `--inbreeding-coefficient`: finite, `[0.0, 1.0]`.
pub fn parse_inbreeding_coefficient(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "inbreeding-coefficient",
        |v| (0.0..=1.0).contains(&v),
        "[0.0, 1.0]",
    )
}

/// `--em-convergence-threshold`: finite, `(0.0, 0.1]`.
pub fn parse_em_convergence_threshold(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "em-convergence-threshold",
        |v| 0.0 < v && v <= CONVERGENCE_THRESHOLD_RANGE_MAX,
        "(0.0, 0.1]",
    )
}

/// `--em-max-iterations`: `1..=MAX_ITERATIONS_RANGE_MAX`.
pub fn parse_em_max_iterations(s: &str) -> Result<u32, String> {
    parse_u32_in(s, "em-max-iterations", 1, MAX_ITERATIONS_RANGE_MAX)
}

/// `--max-gq-phred`: finite, `(10.0, 200.0]`.
pub fn parse_max_gq_phred(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "max-gq-phred",
        |v| GQ_PHRED_RANGE_MIN_EXCLUSIVE < v && v <= GQ_PHRED_RANGE_MAX,
        "(10.0, 200.0]",
    )
}

/// `--min-qual`: finite, `[0.0, 1000.0]`. `0.0` disables the filter.
pub fn parse_min_qual_phred(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "min-qual",
        |v| (0.0..=1000.0).contains(&v),
        "[0.0, 1000.0]",
    )
}

/// `--min-alt-obs-per-sample`: `0..=100`. `0` and `1` are no-ops.
pub fn parse_min_alt_obs_per_sample(s: &str) -> Result<u32, String> {
    parse_u32_in(s, "min-alt-obs-per-sample", 0, 100)
}

/// `--min-mapq-diff-t`: any finite or `-inf` float; finite values
/// must lie in `[-50.0, 50.0]` to catch obvious unit confusion (the
/// effective range under any realistic MAPQ profile is `[-20, 20]`).
/// `-inf` disables the filter without needing `--no-mapq-diff-filter`.
pub fn parse_min_mapq_diff_t(s: &str) -> Result<f32, String> {
    let value: f32 = s
        .parse()
        .map_err(|_| format!("--min-mapq-diff-t: not a valid float: {s:?}"))?;
    if value.is_nan() {
        return Err("--min-mapq-diff-t: NaN is not a valid threshold".into());
    }
    if value.is_finite() && !(-50.0..=50.0).contains(&value) {
        return Err(format!(
            "--min-mapq-diff-t: {value} is outside [-50.0, 50.0] (use -inf to disable)"
        ));
    }
    Ok(value)
}

/// Shared body for every Dirichlet-pseudocount parser. Range is the
/// same across the posterior engine and the contamination side-pass
/// (`PSEUDOCOUNT_RANGE_MAX = 1000.0` in both modules); the only thing
/// that varies between flags is the *name* embedded in the error.
fn parse_pseudocount_with_name(s: &str, flag: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        flag,
        |v| 0.0 < v && v <= PSEUDOCOUNT_RANGE_MAX,
        "(0.0, 1000.0]",
    )
}

/// `--ref-pseudocount`: finite, `(0.0, 1000.0]`.
pub fn parse_ref_pseudocount(s: &str) -> Result<f64, String> {
    parse_pseudocount_with_name(s, "ref-pseudocount")
}

/// `--snp-alt-pseudocount`: finite, `(0.0, 1000.0]`.
pub fn parse_snp_alt_pseudocount(s: &str) -> Result<f64, String> {
    parse_pseudocount_with_name(s, "snp-alt-pseudocount")
}

/// `--indel-alt-pseudocount`: finite, `(0.0, 1000.0]`.
pub fn parse_indel_alt_pseudocount(s: &str) -> Result<f64, String> {
    parse_pseudocount_with_name(s, "indel-alt-pseudocount")
}

/// `--compound-alt-pseudocount`: finite, `(0.0, 1000.0]`.
pub fn parse_compound_alt_pseudocount(s: &str) -> Result<f64, String> {
    parse_pseudocount_with_name(s, "compound-alt-pseudocount")
}

// ---- Side-pass (contamination estimator) ------------------------

/// `--stability-tolerance`: finite, `(0.0, 0.1]`.
pub fn parse_stability_tolerance(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "stability-tolerance",
        |v| 0.0 < v && v <= STABILITY_TOLERANCE_RANGE_MAX,
        "(0.0, 0.1]",
    )
}

/// `--min-major-fraction`: finite, `(0.5, 1.0]`.
pub fn parse_min_major_fraction(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "min-major-fraction",
        |v| MIN_MAJOR_FRACTION_RANGE_MIN_EXCLUSIVE < v && v <= 1.0,
        "(0.5, 1.0]",
    )
}

/// `--min-cohort-minor-fraction`: finite, `[0.0, 0.5)`.
pub fn parse_min_cohort_minor_fraction(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "min-cohort-minor-fraction",
        |v| (0.0..MIN_COHORT_MINOR_FRACTION_RANGE_MAX_EXCLUSIVE).contains(&v),
        "[0.0, 0.5)",
    )
}

/// `--c-s-init`: finite, `[0.0, 0.5]`.
pub fn parse_c_s_init(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "c-s-init",
        |v| (0.0..=C_S_INIT_RANGE_MAX).contains(&v),
        "[0.0, 0.5]",
    )
}

/// `--q-b-init-per-class`: finite, `[0.0, 1.0]`.
pub fn parse_q_b_init_per_class(s: &str) -> Result<f64, String> {
    parse_f64_with(
        s,
        "q-b-init-per-class",
        |v| (0.0..=1.0).contains(&v),
        "[0.0, 1.0]",
    )
}

/// `--block-size`: `>= 1`. The side-pass needs a positive heartbeat
/// to make progress.
pub fn parse_block_size(s: &str) -> Result<u32, String> {
    parse_u32_min(s, "block-size", 1)
}

/// `--min-depth`: `>= 1`.
pub fn parse_min_depth(s: &str) -> Result<u32, String> {
    parse_u32_min(s, "min-depth", 1)
}

/// `--min-batch-size-for-contamination`: `>= 2`. Singletons are
/// unidentifiable per [`crate::var_calling::contamination_estimation`].
///
/// Flag name mirrors the engine-side field name
/// [`crate::var_calling::contamination_estimation::ContaminationEstimationConfig::min_batch_size_for_contamination`]
/// and the artefact parameter-map key — same concept, same name in
/// all three places. Mi14 from the 2026-05-19 cohort CLI review.
pub fn parse_min_batch_size_for_contamination(s: &str) -> Result<u32, String> {
    parse_u32_min(s, "min-batch-size-for-contamination", 2)
}

/// `--min-cohort-minor-count`: `>= 1`. A zero floor admits
/// every site on the count axis and degenerates the informative-site
/// filter to its fraction axis alone — defensible only if explicitly
/// chosen, which we don't expose. Default is 2 per
/// [`crate::var_calling::contamination_estimation::DEFAULT_MIN_COHORT_MINOR_COUNT`].
pub fn parse_min_cohort_minor_count(s: &str) -> Result<u32, String> {
    parse_u32_min(s, "min-cohort-minor-count", 1)
}

/// `--stability-blocks`: `>= 1`.
pub fn parse_stability_blocks(s: &str) -> Result<u32, String> {
    parse_u32_min(s, "stability-blocks", 1)
}

// ---- Stage 1 (PSP writer) ---------------------------------------

/// `--block-target-bytes`: `MIN_BLOCK_TARGET_BYTES..=MAX_BLOCK_TARGET_BYTES`
/// (16 KiB..=16 MiB). PSP writer auto-flushes when an open block's
/// projected uncompressed payload reaches this value.
///
/// The range matches the regime covered by the 2026-05-27 sweep on
/// real per-sample PSPs. Below the floor zstd loses meaningful
/// compression context (the sweep showed 64 KiB at +75% disk already
/// past the knee); above the ceiling sits the legacy hardcoded
/// 16 MiB default that was retired for cohort-step memory blow-up.
pub fn parse_block_target_bytes(s: &str) -> Result<usize, String> {
    parse_usize_in(
        s,
        "block-target-bytes",
        MIN_BLOCK_TARGET_BYTES,
        MAX_BLOCK_TARGET_BYTES,
    )
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- f64-shaped parsers ----------------------------------------

    #[test]
    fn em_convergence_threshold_accepts_default_and_max() {
        assert!((parse_em_convergence_threshold("1e-4").unwrap() - 1e-4).abs() < 1e-12);
        assert_eq!(
            parse_em_convergence_threshold("0.1").unwrap(),
            CONVERGENCE_THRESHOLD_RANGE_MAX
        );
    }

    #[test]
    fn em_convergence_threshold_rejects_zero_negative_nan_inf_above_max() {
        assert!(parse_em_convergence_threshold("0").is_err());
        assert!(parse_em_convergence_threshold("-0.001").is_err());
        assert!(parse_em_convergence_threshold("nan").is_err());
        assert!(parse_em_convergence_threshold("inf").is_err());
        assert!(parse_em_convergence_threshold("0.5").is_err());
        assert!(parse_em_convergence_threshold("hello").is_err());
    }

    #[test]
    fn inbreeding_coefficient_boundaries() {
        parse_inbreeding_coefficient("0.0").unwrap();
        parse_inbreeding_coefficient("0.5").unwrap();
        parse_inbreeding_coefficient("1.0").unwrap();
        assert!(parse_inbreeding_coefficient("-0.001").is_err());
        assert!(parse_inbreeding_coefficient("1.001").is_err());
        assert!(parse_inbreeding_coefficient("nan").is_err());
    }

    #[test]
    fn max_gq_phred_boundaries() {
        parse_max_gq_phred("99").unwrap();
        parse_max_gq_phred("200").unwrap();
        assert!(parse_max_gq_phred("10").is_err());
        assert!(parse_max_gq_phred("200.001").is_err());
        assert!(parse_max_gq_phred("nan").is_err());
    }

    #[test]
    fn min_qual_phred_boundaries() {
        parse_min_qual_phred("0").unwrap();
        parse_min_qual_phred("30").unwrap();
        parse_min_qual_phred("1000").unwrap();
        assert!(parse_min_qual_phred("-0.001").is_err());
        assert!(parse_min_qual_phred("1000.001").is_err());
        assert!(parse_min_qual_phred("nan").is_err());
        assert!(parse_min_qual_phred("inf").is_err());
    }

    #[test]
    fn min_alt_obs_per_sample_boundaries() {
        parse_min_alt_obs_per_sample("0").unwrap();
        parse_min_alt_obs_per_sample("1").unwrap();
        parse_min_alt_obs_per_sample("2").unwrap();
        parse_min_alt_obs_per_sample("100").unwrap();
        assert!(parse_min_alt_obs_per_sample("101").is_err());
        assert!(parse_min_alt_obs_per_sample("-1").is_err());
        assert!(parse_min_alt_obs_per_sample("abc").is_err());
    }

    #[test]
    fn pseudocount_boundaries() {
        // Every per-flag entry delegates to the same body, so the
        // range check only needs one boundary sweep. We probe each
        // entry once so a future regression that drops the wiring
        // for any one of them still surfaces.
        parse_ref_pseudocount("0.01").unwrap();
        parse_snp_alt_pseudocount("10.0").unwrap();
        parse_indel_alt_pseudocount("1000.0").unwrap();
        parse_compound_alt_pseudocount("0.001").unwrap();
        assert!(parse_ref_pseudocount("0").is_err());
        assert!(parse_snp_alt_pseudocount("-1").is_err());
        assert!(parse_indel_alt_pseudocount("1000.1").is_err());
        assert!(parse_compound_alt_pseudocount("inf").is_err());
    }

    #[test]
    fn pseudocount_error_carries_flag_name() {
        // The whole point of the per-flag split (Mi1): the rendered
        // error string identifies *which* knob took the bad value.
        let err = parse_snp_alt_pseudocount("2000").unwrap_err();
        assert!(
            err.contains("snp-alt-pseudocount"),
            "error should name the flag, got: {err}"
        );
    }

    #[test]
    fn stability_tolerance_boundaries() {
        parse_stability_tolerance("1e-3").unwrap();
        parse_stability_tolerance("0.1").unwrap();
        assert!(parse_stability_tolerance("0").is_err());
        assert!(parse_stability_tolerance("0.5").is_err());
        assert!(parse_stability_tolerance("nan").is_err());
    }

    #[test]
    fn min_major_fraction_boundaries() {
        parse_min_major_fraction("0.95").unwrap();
        parse_min_major_fraction("1.0").unwrap();
        assert!(parse_min_major_fraction("0.5").is_err());
        assert!(parse_min_major_fraction("1.001").is_err());
    }

    #[test]
    fn min_cohort_minor_fraction_boundaries() {
        parse_min_cohort_minor_fraction("0.0").unwrap();
        parse_min_cohort_minor_fraction("0.005").unwrap();
        assert!(parse_min_cohort_minor_fraction("0.5").is_err()); // exclusive upper
        assert!(parse_min_cohort_minor_fraction("-0.001").is_err());
    }

    #[test]
    fn c_s_init_boundaries() {
        parse_c_s_init("0.0").unwrap();
        parse_c_s_init("0.02").unwrap();
        parse_c_s_init("0.5").unwrap();
        assert!(parse_c_s_init("0.501").is_err());
        assert!(parse_c_s_init("-0.001").is_err());
    }

    #[test]
    fn q_b_init_per_class_boundaries() {
        parse_q_b_init_per_class("0.0").unwrap();
        parse_q_b_init_per_class("0.333").unwrap();
        parse_q_b_init_per_class("1.0").unwrap();
        assert!(parse_q_b_init_per_class("1.001").is_err());
    }

    // ---- u32 / u8 / usize parsers ----------------------------------

    #[test]
    fn ploidy_boundaries() {
        parse_ploidy("1").unwrap();
        parse_ploidy("2").unwrap();
        parse_ploidy("8").unwrap();
        assert!(parse_ploidy("0").is_err());
        assert!(parse_ploidy("9").is_err());
        assert!(parse_ploidy("-1").is_err());
        assert!(parse_ploidy("abc").is_err());
    }

    #[test]
    fn max_alleles_boundaries() {
        parse_max_alleles("2").unwrap();
        parse_max_alleles("6").unwrap();
        parse_max_alleles("16").unwrap();
        assert!(parse_max_alleles("1").is_err());
        assert!(parse_max_alleles("17").is_err());
    }

    #[test]
    fn em_max_iterations_boundaries() {
        parse_em_max_iterations("1").unwrap();
        parse_em_max_iterations("50").unwrap();
        parse_em_max_iterations("500").unwrap();
        assert!(parse_em_max_iterations("0").is_err());
        assert!(parse_em_max_iterations("501").is_err());
    }

    #[test]
    fn var_group_max_span_rejects_zero() {
        parse_var_group_max_span("1").unwrap();
        parse_var_group_max_span("10000").unwrap();
        assert!(parse_var_group_max_span("0").is_err());
    }

    #[test]
    fn dust_window_boundaries() {
        parse_dust_window("3").unwrap();
        parse_dust_window("64").unwrap();
        parse_dust_window("65535").unwrap();
        assert!(parse_dust_window("2").is_err());
        assert!(parse_dust_window("65536").is_err());
    }

    #[test]
    fn dust_threshold_admits_any_u32() {
        parse_dust_threshold("0").unwrap();
        parse_dust_threshold("20").unwrap();
        parse_dust_threshold("4294967295").unwrap(); // u32::MAX
        assert!(parse_dust_threshold("-1").is_err());
    }

    #[test]
    fn block_size_rejects_zero() {
        parse_block_size("1").unwrap();
        parse_block_size("1000").unwrap();
        assert!(parse_block_size("0").is_err());
    }

    #[test]
    fn min_depth_rejects_zero() {
        parse_min_depth("1").unwrap();
        parse_min_depth("10").unwrap();
        assert!(parse_min_depth("0").is_err());
    }

    #[test]
    fn min_batch_size_for_contamination_rejects_below_two() {
        parse_min_batch_size_for_contamination("2").unwrap();
        parse_min_batch_size_for_contamination("5").unwrap();
        assert!(parse_min_batch_size_for_contamination("0").is_err());
        assert!(parse_min_batch_size_for_contamination("1").is_err());
    }

    #[test]
    fn stability_blocks_rejects_zero() {
        parse_stability_blocks("1").unwrap();
        parse_stability_blocks("3").unwrap();
        assert!(parse_stability_blocks("0").is_err());
    }

    #[test]
    fn block_target_bytes_accepts_range_and_rejects_outside() {
        // Both endpoints accepted.
        assert_eq!(
            parse_block_target_bytes(&MIN_BLOCK_TARGET_BYTES.to_string()).unwrap(),
            MIN_BLOCK_TARGET_BYTES
        );
        assert_eq!(
            parse_block_target_bytes(&MAX_BLOCK_TARGET_BYTES.to_string()).unwrap(),
            MAX_BLOCK_TARGET_BYTES
        );
        // A representative point in the middle (the new default).
        parse_block_target_bytes("1048576").unwrap(); // 1 MiB
        // One below the floor and one above the ceiling — both errors,
        // and the error mentions the flag name (so the user can tell
        // which knob misbehaved).
        let below = parse_block_target_bytes(&(MIN_BLOCK_TARGET_BYTES - 1).to_string());
        assert!(below.is_err());
        assert!(below.unwrap_err().contains("block-target-bytes"));
        assert!(parse_block_target_bytes(&(MAX_BLOCK_TARGET_BYTES + 1).to_string()).is_err());
        // Garbage input.
        assert!(parse_block_target_bytes("0").is_err());
        assert!(parse_block_target_bytes("abc").is_err());
        assert!(parse_block_target_bytes("-1").is_err());
    }
}
