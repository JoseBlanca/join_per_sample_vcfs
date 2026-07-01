//! Per-sample summary statistics carried in the `.psp` metadata section.
//!
//! The Stage-1 pileup computes two per-sample summaries the hidden-paralog
//! filter needs (architecture
//! `doc/devel/architecture/hidden_paralog_psp_integration.md`):
//!
//! - a **coverage-by-GC histogram** — a raw 2-D count matrix of tiled
//!   windows keyed by `(GC fraction, covered-bases mean depth)`, from
//!   which var-calling fits the depth∼GC curve and single-copy scale;
//! - **observed heterozygosity** — four counts (confident het / hom-alt /
//!   ambiguous / total variant sites) from a rough per-site binomial-LR
//!   genotype call.
//!
//! The TOML document model + serialisation lives here; the accumulators
//! that *produce* the summaries from the Stage-1 stream live in submodules
//! ([`coverage`], [`het`]). The model fit that *consumes* the histogram
//! lives downstream in var-calling.

pub mod coverage;
pub mod het;

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Schema version of the sample-summary TOML document. A flat counter
/// (no major/minor split): every bump is breaking, and a reader accepts
/// only `1..=SAMPLE_SUMMARY_VERSION`.
///
/// - v1: `coverage_by_gc` + `heterozygosity`.
/// - v2: adds `coverage_by_gc.callable_positions` (the het-rate denominator,
///   P1 of the hidden-paralog filter). A required field (no
///   `#[serde(default)]`), so a stale on-disk v1 document is rejected at the
///   TOML layer as a missing-field [`SampleSummaryError::ParseToml`] *before*
///   the version guard in [`SampleSummary::validate`] runs — the guard
///   catches only future/zeroed versions, not v1 on disk. Pre-alpha, no
///   backwards-compatibility promise: regenerate old summary sections by
///   re-running `pileup`.
pub const SAMPLE_SUMMARY_VERSION: u16 = 2;

/// Default tile / GC covariate window in bp (architecture Premise 3). The
/// `pileup --gc-window-bp` flag overrides it.
pub const DEFAULT_GC_WINDOW_BP: u32 = 500;

/// Default number of GC bins (≈ 2 % GC resolution). Tuning, not
/// architecture (Premise 2) — recorded in the `.psp` so it can change
/// without breaking readers.
pub const DEFAULT_GC_BINS: u32 = 50;

/// Default depth-bin width in covered-bases-mean-depth units.
pub const DEFAULT_DEPTH_BIN_WIDTH: f64 = 0.5;

/// Default number of regular depth bins (so the regular range is
/// `0..DEFAULT_DEPTH_BINS * DEFAULT_DEPTH_BIN_WIDTH` = `0..100×`, with an
/// overflow bin above — generous for high-copy paralogs at typical depth).
pub const DEFAULT_DEPTH_BINS: u32 = 200;

/// Default minimum total fragment depth for a site to enter the het call.
pub const DEFAULT_HET_MIN_DEPTH: u32 = 4;

/// Default per-read error rate `ε` for the het hom-model.
pub const DEFAULT_HET_ERROR_RATE: f64 = 0.02;

/// Default het confidence margin `M` (nats) = `ln 10`, i.e. 10:1 odds.
pub const DEFAULT_HET_LR_MARGIN: f64 = std::f64::consts::LN_10;

/// The per-sample summary document stored (TOML, then zstd-framed) in the
/// `.psp` metadata section. Serialises with kebab-case keys; `version`
/// precedes the two tables so the TOML is well-formed (top-level scalars
/// before tables).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[non_exhaustive]
pub struct SampleSummary {
    /// Document schema version; see [`SAMPLE_SUMMARY_VERSION`].
    pub version: u16,
    pub coverage_by_gc: CoverageByGcHistogram,
    pub heterozygosity: HetCounts,
}

/// Raw coverage-by-GC histogram: a row-major `[gc_bin][depth_bin]` count
/// matrix over tiled windows, plus the bin schemes needed to interpret a
/// cell. The curve / single-copy-scale fit is downstream (these are
/// sufficient statistics, not a fitted model).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[non_exhaustive]
pub struct CoverageByGcHistogram {
    /// Tile width in bp — the GC covariate window, equal to the analysis
    /// window (architecture Premise 3).
    pub window_bp: u32,
    /// Number of GC bins, uniform over the closed unit interval `[0, 1]`
    /// (a tile of GC fraction `g` lands in bin `min(floor(g * gc_bins),
    /// gc_bins - 1)`).
    pub gc_bins: u32,
    /// Width of one depth bin in covered-bases-mean-depth units.
    pub depth_bin_width: f64,
    /// Number of regular depth bins. One overflow bin (for depth at or
    /// above `depth_bins * depth_bin_width`) follows them, so each GC row
    /// holds `depth_bins + 1` cells.
    pub depth_bins: u32,
    /// Tiles folded into the histogram (its support).
    pub n_tiles: u64,
    /// Tiles skipped because they had no GC-defined (non-`N`) covered
    /// positions, so no `(GC, depth)` sample could be formed.
    pub n_skipped_tiles: u64,
    /// Grand total of GC-defined (non-`N`) covered positions summed across
    /// every folded tile — the sample's callable-position count. It is the
    /// denominator of the observed-heterozygosity **rate**
    /// `Hobs = n_het_sites / callable_positions` the hidden-paralog filter
    /// consumes (spec §3; *not* `n_het/(n_het+n_hom_alt)`, which tracks
    /// reference divergence and inverts). The binned `counts` matrix cannot
    /// recover it — a tile's per-position covered count is lost when the
    /// tile collapses to one `(GC, mean depth)` cell — so it is kept as its
    /// own running total. Every folded tile contributes at least one covered
    /// position, so `callable_positions >= n_tiles`.
    pub callable_positions: u64,
    /// Row-major `[gc_bin][depth_bin]` counts. Length is exactly
    /// `gc_bins * (depth_bins + 1)`.
    pub counts: Vec<u32>,
}

/// Observed-heterozygosity counts over variant sites, from the rough
/// per-site **binomial het-vs-hom likelihood-ratio** genotype call
/// (architecture Premise 1b): each variant site is confident-het
/// (`logLR > +margin`), confident-hom-alt (`logLR < −margin`), or
/// ambiguous (`|logLR| ≤ margin`). The observed-het **rate** the paralog
/// filter consumes is `Hobs = n_het_sites / callable_positions` (the
/// denominator lives on [`CoverageByGcHistogram::callable_positions`]); the
/// confident ratio `n_het / (n_het + n_hom_alt)` is a *reference-divergence*
/// proxy, **not** `Hobs` — it is dominated by the hom-alt count and inverts
/// (spec §3). `n_ambiguous_sites` is the per-sample uncertainty weight
/// (large at low coverage). Counts, not a logLR histogram, because the
/// classification margin is a *settled* threshold (contrast the
/// still-calibrating coverage curve — Premise 2).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[non_exhaustive]
pub struct HetCounts {
    /// Variant sites called confident heterozygous (`logLR > +margin`).
    pub n_het_sites: u64,
    /// Variant sites called confident homozygous-alt (`logLR < −margin`).
    pub n_hom_alt_sites: u64,
    /// Variant sites where het vs hom-alt is unresolved
    /// (`|logLR| ≤ margin`) — the uncertainty signal.
    pub n_ambiguous_sites: u64,
    /// Total variant sites; must equal
    /// `n_het_sites + n_hom_alt_sites + n_ambiguous_sites`.
    pub n_variant_sites: u64,
    /// Minimum total depth for a site to be considered a candidate
    /// (recorded for reproducibility of the rough genotype call).
    pub min_depth: u32,
    /// Per-read error rate `ε` in the hom-alt model `Binomial(n, 1 − ε)`.
    /// In the open interval `(0, 1)` (the LR's `ln ε` / `ln(1 − ε)` are
    /// `−∞` at the endpoints).
    pub error_rate: f64,
    /// Confidence margin `M` (nats) splitting het / hom-alt / ambiguous.
    /// `>= 0`.
    pub lr_margin: f64,
}

/// Failure modes for building / parsing a [`SampleSummary`].
///
/// `#[non_exhaustive]` (matching the crate's on-disk artefact convention,
/// e.g. `ContaminationArtefactError`) so future variants land additively.
/// The `toml` source types are exposed directly, also matching that
/// convention — this crate is not a published library, so the dependency
/// is not part of an external API contract.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SampleSummaryError {
    /// `toml` could not serialise the document (should not happen for our
    /// schema; the surface exists so the failure is never swallowed).
    #[error("failed to serialise sample summary to TOML")]
    SerializeToml {
        #[source]
        source: toml::ser::Error,
    },
    /// The metadata-section bytes are not valid sample-summary TOML.
    #[error("failed to parse sample summary TOML")]
    ParseToml {
        #[source]
        source: toml::de::Error,
    },
    /// The metadata-section bytes are not valid UTF-8, so they cannot be
    /// the TOML document. Carries the underlying [`std::str::Utf8Error`]
    /// rather than flattening it into a string.
    #[error("metadata section is not valid UTF-8")]
    NotUtf8 {
        #[source]
        source: std::str::Utf8Error,
    },
    /// The document's version is not one this reader recognises. Every
    /// version bump is breaking (the version is a flat counter, not
    /// major/minor), so a reader accepts only `1..=SAMPLE_SUMMARY_VERSION`
    /// — both `0` (never written; a zeroed/truncated blob) and any future
    /// version are refused.
    #[error("unsupported sample-summary version {got} (this reader supports 1..={supported})")]
    UnsupportedVersion { got: u16, supported: u16 },
    /// `counts.len()` disagrees with `gc_bins * (depth_bins + 1)`.
    #[error(
        "coverage histogram counts length {got} != gc_bins {gc_bins} * (depth_bins {depth_bins} + 1) = {expected}"
    )]
    HistogramShapeMismatch {
        got: usize,
        gc_bins: u32,
        depth_bins: u32,
        expected: usize,
    },
    /// A field violated a value invariant (zero where positive required,
    /// non-finite float, error rate outside `(0, 1)`, negative margin, or
    /// the het class counts not summing to `n_variant_sites`).
    #[error("invalid sample-summary field {field}: {reason}")]
    InvalidField { field: &'static str, reason: String },
}

impl SampleSummary {
    /// Serialise to the TOML bytes that fill the `.psp` metadata section.
    /// Validates first, so a malformed summary never reaches disk.
    pub fn to_toml_bytes(&self) -> Result<Vec<u8>, SampleSummaryError> {
        self.validate()?;
        let body = toml::to_string_pretty(self)
            .map_err(|source| SampleSummaryError::SerializeToml { source })?;
        Ok(body.into_bytes())
    }

    /// Parse and validate a summary from metadata-section bytes.
    pub fn from_toml_bytes(bytes: &[u8]) -> Result<Self, SampleSummaryError> {
        let text =
            std::str::from_utf8(bytes).map_err(|source| SampleSummaryError::NotUtf8 { source })?;
        let summary: SampleSummary =
            toml::from_str(text).map_err(|source| SampleSummaryError::ParseToml { source })?;
        summary.validate()?;
        Ok(summary)
    }

    /// Check every structural and value invariant. Run on both the
    /// serialise and parse paths so a bad summary is rejected at the
    /// producer and never trusted at the consumer.
    pub fn validate(&self) -> Result<(), SampleSummaryError> {
        // The version is a flat counter (no major/minor); a reader accepts
        // only versions it recognises. `0` is never written — it is the
        // value a zeroed or truncated blob deserialises to — so it is
        // refused alongside any future version.
        if self.version == 0 || self.version > SAMPLE_SUMMARY_VERSION {
            return Err(SampleSummaryError::UnsupportedVersion {
                got: self.version,
                supported: SAMPLE_SUMMARY_VERSION,
            });
        }
        self.coverage_by_gc.validate()?;
        self.heterozygosity.validate()?;
        Ok(())
    }
}

impl CoverageByGcHistogram {
    fn validate(&self) -> Result<(), SampleSummaryError> {
        let bad = |field, reason: String| SampleSummaryError::InvalidField { field, reason };
        if self.window_bp == 0 {
            return Err(bad("coverage-by-gc.window-bp", "must be >= 1".to_string()));
        }
        if self.gc_bins == 0 {
            return Err(bad("coverage-by-gc.gc-bins", "must be >= 1".to_string()));
        }
        if self.depth_bins == 0 {
            return Err(bad("coverage-by-gc.depth-bins", "must be >= 1".to_string()));
        }
        if !(self.depth_bin_width.is_finite() && self.depth_bin_width > 0.0) {
            return Err(bad(
                "coverage-by-gc.depth-bin-width",
                format!("must be finite and > 0, got {}", self.depth_bin_width),
            ));
        }
        // `depth_bins + 1` (overflow bin) per GC row. Use u64/usize to
        // avoid overflow on hostile bin counts before the length check.
        let expected = (self.gc_bins as usize)
            .checked_mul(self.depth_bins as usize + 1)
            .ok_or_else(|| {
                bad(
                    "coverage-by-gc.counts",
                    "gc_bins * (depth_bins + 1) overflows usize".to_string(),
                )
            })?;
        if self.counts.len() != expected {
            return Err(SampleSummaryError::HistogramShapeMismatch {
                got: self.counts.len(),
                gc_bins: self.gc_bins,
                depth_bins: self.depth_bins,
                expected,
            });
        }
        // Every folded (non-skipped) tile carries at least one GC-defined
        // covered position, and `callable_positions` sums exactly those
        // positions, so it can never be smaller than the folded-tile count.
        // A document violating this is corrupt (or was built by a producer
        // that miscounts) and would yield a nonsensical het rate.
        if self.callable_positions < self.n_tiles {
            return Err(bad(
                "coverage-by-gc.callable-positions",
                format!(
                    "{} < n-tiles {} (each folded tile has >= 1 covered position)",
                    self.callable_positions, self.n_tiles,
                ),
            ));
        }
        Ok(())
    }
}

impl HetCounts {
    fn validate(&self) -> Result<(), SampleSummaryError> {
        let bad = |field, reason: String| SampleSummaryError::InvalidField { field, reason };
        // The three class counts must sum to the recorded total. Checked
        // with `checked_add` so a hostile/corrupt overflowing triple is
        // caught rather than wrapping to a value that happens to match.
        let sum = self
            .n_het_sites
            .checked_add(self.n_hom_alt_sites)
            .and_then(|s| s.checked_add(self.n_ambiguous_sites));
        if sum != Some(self.n_variant_sites) {
            return Err(bad(
                "heterozygosity.n-variant-sites",
                format!(
                    "{} != n-het-sites {} + n-hom-alt-sites {} + n-ambiguous-sites {}",
                    self.n_variant_sites,
                    self.n_het_sites,
                    self.n_hom_alt_sites,
                    self.n_ambiguous_sites,
                ),
            ));
        }
        if !(self.error_rate.is_finite() && self.error_rate > 0.0 && self.error_rate < 1.0) {
            return Err(bad(
                "heterozygosity.error-rate",
                format!("must be finite in (0, 1), got {}", self.error_rate),
            ));
        }
        if !(self.lr_margin.is_finite() && self.lr_margin >= 0.0) {
            return Err(bad(
                "heterozygosity.lr-margin",
                format!("must be finite and >= 0, got {}", self.lr_margin),
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> SampleSummary {
        SampleSummary {
            version: SAMPLE_SUMMARY_VERSION,
            coverage_by_gc: CoverageByGcHistogram {
                window_bp: 500,
                gc_bins: 2,
                depth_bin_width: 0.5,
                depth_bins: 3,
                n_tiles: 42,
                n_skipped_tiles: 1,
                callable_positions: 4200,
                // gc_bins (2) * (depth_bins + 1 = 4) = 8 cells.
                counts: vec![0, 1, 2, 3, 4, 5, 6, 7],
            },
            heterozygosity: HetCounts {
                n_het_sites: 100,
                n_hom_alt_sites: 120,
                n_ambiguous_sites: 30,
                n_variant_sites: 250,
                min_depth: 4,
                error_rate: 0.02,
                lr_margin: std::f64::consts::LN_10,
            },
        }
    }

    #[test]
    fn round_trips_through_toml_bytes() {
        let s = sample();
        let bytes = s.to_toml_bytes().expect("serialise");
        let back = SampleSummary::from_toml_bytes(&bytes).expect("parse");
        assert_eq!(s, back);
    }

    #[test]
    fn wire_keys_are_kebab_case() {
        let body = String::from_utf8(sample().to_toml_bytes().unwrap()).unwrap();
        for key in [
            "[coverage-by-gc]",
            "window-bp",
            "gc-bins",
            "depth-bin-width",
            "n-skipped-tiles",
            "callable-positions",
            "[heterozygosity]",
            "n-het-sites",
            "n-hom-alt-sites",
            "n-ambiguous-sites",
            "n-variant-sites",
            "error-rate",
            "lr-margin",
        ] {
            assert!(body.contains(key), "wire key {key:?} absent from:\n{body}");
        }
        // snake_case must not leak.
        for forbidden in [
            "window_bp",
            "gc_bins",
            "callable_positions",
            "n_het_sites",
            "n_hom_alt_sites",
            "lr_margin",
        ] {
            assert!(!body.contains(forbidden), "snake_case {forbidden:?} leaked");
        }
    }

    #[test]
    fn rejects_histogram_shape_mismatch() {
        let mut s = sample();
        s.coverage_by_gc.counts.push(99); // now 9, expected 8
        let err = s.to_toml_bytes().expect_err("shape mismatch must fail");
        assert!(
            matches!(
                err,
                SampleSummaryError::HistogramShapeMismatch {
                    got: 9,
                    expected: 8,
                    ..
                }
            ),
            "got {err:?}"
        );
    }

    #[test]
    fn rejects_error_rate_out_of_range() {
        for bad_eps in [0.0, 1.0, -0.1, f64::NAN] {
            let mut s = sample();
            s.heterozygosity.error_rate = bad_eps;
            let err = s.to_toml_bytes().expect_err("bad error-rate must fail");
            assert!(
                matches!(err, SampleSummaryError::InvalidField { field, .. } if field == "heterozygosity.error-rate"),
                "error_rate {bad_eps} -> {err:?}"
            );
        }
    }

    #[test]
    fn rejects_count_sum_mismatch() {
        let mut s = sample();
        // n_variant no longer equals het + hom_alt + ambiguous.
        s.heterozygosity.n_variant_sites += 1;
        let err = s.to_toml_bytes().expect_err("sum mismatch must fail");
        assert!(
            matches!(err, SampleSummaryError::InvalidField { field, .. } if field == "heterozygosity.n-variant-sites"),
            "got {err:?}"
        );
    }

    #[test]
    fn rejects_future_version_on_parse() {
        let mut s = sample();
        // Serialise a valid doc, then bump the version in the text and
        // re-parse — simulates a newer producer.
        let bytes = s.to_toml_bytes().unwrap();
        s.version = SAMPLE_SUMMARY_VERSION + 1;
        let bumped = String::from_utf8(bytes).unwrap().replace(
            &format!("version = {SAMPLE_SUMMARY_VERSION}"),
            &format!("version = {}", SAMPLE_SUMMARY_VERSION + 1),
        );
        let err = SampleSummary::from_toml_bytes(bumped.as_bytes()).expect_err("future version");
        assert!(
            matches!(err, SampleSummaryError::UnsupportedVersion { .. }),
            "got {err:?}"
        );
    }

    #[test]
    fn rejects_non_finite_depth_bin_width() {
        let mut s = sample();
        s.coverage_by_gc.depth_bin_width = 0.0;
        let err = s.validate().expect_err("zero width must fail");
        assert!(
            matches!(err, SampleSummaryError::InvalidField { field, .. } if field == "coverage-by-gc.depth-bin-width"),
            "got {err:?}"
        );
    }

    /// `version = 0` is the zeroed/truncated-blob default, never written —
    /// it must be refused, not silently trusted.
    #[test]
    fn validate_rejects_version_zero() {
        let mut s = sample();
        s.version = 0;
        let err = s.validate().expect_err("version 0 must be rejected");
        assert!(
            matches!(err, SampleSummaryError::UnsupportedVersion { got: 0, .. }),
            "got {err:?}"
        );
    }

    /// Each positivity guard fires at `0`. `gc_bins = 0` is the sharp case:
    /// it also makes the expected `counts` length `0`, so without the
    /// dedicated guard an empty histogram would pass the length check.
    #[test]
    fn validate_rejects_zero_bin_dimensions() {
        for (mutate, field) in [
            (
                (|s: &mut SampleSummary| s.coverage_by_gc.window_bp = 0) as fn(&mut SampleSummary),
                "coverage-by-gc.window-bp",
            ),
            (
                |s: &mut SampleSummary| {
                    s.coverage_by_gc.gc_bins = 0;
                    s.coverage_by_gc.counts.clear();
                },
                "coverage-by-gc.gc-bins",
            ),
            (
                |s: &mut SampleSummary| s.coverage_by_gc.depth_bins = 0,
                "coverage-by-gc.depth-bins",
            ),
        ] {
            let mut s = sample();
            mutate(&mut s);
            let err = s.validate().expect_err("zero dimension must fail");
            assert!(
                matches!(err, SampleSummaryError::InvalidField { field: f, .. } if f == field),
                "expected InvalidField {field:?}, got {err:?}"
            );
        }
    }

    /// The `callable_positions >= n_tiles` invariant is enforced: a document
    /// claiming fewer callable positions than folded tiles (impossible — each
    /// folded tile carries >= 1 covered position) is rejected. Without this
    /// test a future refactor that drops or inverts the check would pass CI.
    #[test]
    fn validate_rejects_callable_positions_below_n_tiles() {
        let mut s = sample(); // n_tiles: 42, callable_positions: 4200
        s.coverage_by_gc.callable_positions = 41; // < n_tiles 42
        let err = s.validate().expect_err("callable < n_tiles must fail");
        assert!(
            matches!(err, SampleSummaryError::InvalidField { field, .. }
                     if field == "coverage-by-gc.callable-positions"),
            "got {err:?}"
        );
    }

    /// The accept side of the boundary: `callable_positions == n_tiles` (one
    /// covered position per tile) is the tight legitimate edge and must pass.
    /// Pins the `>=` against an off-by-one flip to `>`.
    #[test]
    fn validate_accepts_callable_positions_equal_to_n_tiles() {
        let mut s = sample();
        s.coverage_by_gc.n_tiles = 8;
        s.coverage_by_gc.callable_positions = 8; // exactly one covered pos/tile
        s.validate()
            .expect("callable == n_tiles is the valid boundary");
    }

    /// Non-UTF-8 metadata bytes surface as `NotUtf8` (carrying the cause),
    /// not a panic or a mislabelled field error.
    #[test]
    fn from_toml_bytes_rejects_invalid_utf8() {
        let err = SampleSummary::from_toml_bytes(&[0xff, 0xfe, 0x00]).expect_err("invalid utf-8");
        assert!(
            matches!(err, SampleSummaryError::NotUtf8 { .. }),
            "got {err:?}"
        );
    }

    /// Valid UTF-8 that is not TOML surfaces as `ParseToml`.
    #[test]
    fn from_toml_bytes_rejects_malformed_toml() {
        let err =
            SampleSummary::from_toml_bytes(b"this is not = = toml").expect_err("malformed toml");
        assert!(
            matches!(err, SampleSummaryError::ParseToml { .. }),
            "got {err:?}"
        );
    }

    /// `validate` fires on the *parse* path, not just on serialise: a
    /// parseable, valid-version document that violates a value invariant is
    /// rejected by `from_toml_bytes`.
    #[test]
    fn from_toml_bytes_rejects_value_invariant_violation() {
        let body = String::from_utf8(sample().to_toml_bytes().unwrap())
            .unwrap()
            // Break the count identity (250 = 100 + 120 + 30 -> 999).
            .replace("n-variant-sites = 250", "n-variant-sites = 999");
        let err = SampleSummary::from_toml_bytes(body.as_bytes())
            .expect_err("value invariant must be enforced on parse");
        assert!(
            matches!(err, SampleSummaryError::InvalidField { field, .. } if field == "heterozygosity.n-variant-sites"),
            "got {err:?}"
        );
    }

    proptest::proptest! {
        /// Round-trip over random valid documents, including boundary values
        /// (`gc_bins`/`depth_bins` at 1, `u32::MAX` counts, the count
        /// identity, error-rate near the open-interval edges, `min_depth =
        /// 0`). Catches float-format / large-integer drift and kebab-case
        /// key collisions a single hand-built example misses.
        #[test]
        fn round_trips_arbitrary_valid_summary(
            gc_bins in 1u32..6,
            depth_bins in 1u32..6,
            depth_bin_width in 0.001f64..1e6,
            error_rate in 0.000_1f64..0.999_9,
            lr_margin in 0.0f64..50.0,
            n_het in 0u64..1000,
            n_hom in 0u64..1000,
            n_amb in 0u64..1000,
            window_bp in 1u32..100_000,
        ) {
            let cells = (gc_bins as usize) * (depth_bins as usize + 1);
            let s = SampleSummary {
                version: SAMPLE_SUMMARY_VERSION,
                coverage_by_gc: CoverageByGcHistogram {
                    window_bp,
                    gc_bins,
                    depth_bin_width,
                    depth_bins,
                    n_tiles: n_het + 7,
                    n_skipped_tiles: 3,
                    callable_positions: n_het + 7,
                    counts: vec![u32::MAX; cells],
                },
                heterozygosity: HetCounts {
                    n_het_sites: n_het,
                    n_hom_alt_sites: n_hom,
                    n_ambiguous_sites: n_amb,
                    n_variant_sites: n_het + n_hom + n_amb,
                    min_depth: 0,
                    error_rate,
                    lr_margin,
                },
            };
            let bytes = s.to_toml_bytes().expect("serialise");
            let back = SampleSummary::from_toml_bytes(&bytes).expect("parse");
            proptest::prop_assert_eq!(back, s);
        }
    }
}
