//! Per-sample summary statistics carried in the `.psp` metadata section.
//!
//! The Stage-1 pileup computes two per-sample summaries the hidden-paralog
//! filter needs (architecture
//! `doc/devel/architecture/hidden_paralog_psp_integration.md`):
//!
//! - a **coverage-by-GC histogram** — a raw 2-D count matrix of tiled
//!   windows keyed by `(GC fraction, covered-bases mean depth)`, from
//!   which var-calling fits the depth∼GC curve and single-copy scale;
//! - **observed heterozygosity** — two counts (`n_het_sites`,
//!   `n_variant_sites`) from a rough per-site genotype call.
//!
//! The TOML document model + serialisation lives here; the accumulators
//! that *produce* the summaries from the Stage-1 stream live in submodules
//! ([`coverage`]; het in a later step). The model fit that *consumes* the
//! histogram lives downstream in var-calling.

pub mod coverage;

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Schema version of the sample-summary TOML document. A flat counter
/// (no major/minor split): every bump is breaking, and a reader accepts
/// only `1..=SAMPLE_SUMMARY_VERSION`.
pub const SAMPLE_SUMMARY_VERSION: u16 = 1;

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
    /// Row-major `[gc_bin][depth_bin]` counts. Length is exactly
    /// `gc_bins * (depth_bins + 1)`.
    pub counts: Vec<u32>,
}

/// Observed-heterozygosity counts: among sites where the sample carries a
/// real minor allele (`n_variant_sites`), how many were called
/// heterozygous (`n_het_sites`). `Hobs = n_het_sites / n_variant_sites`
/// is formed downstream; storing the counts preserves the support so a
/// low-evidence estimate can be down-weighted.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[non_exhaustive]
pub struct HetCounts {
    /// Sites called heterozygous (minor-allele VAF in the het band).
    pub n_het_sites: u64,
    /// Sites with a real minor allele (the denominator of `Hobs`).
    pub n_variant_sites: u64,
    /// Minimum total depth for a site to be considered (recorded for
    /// reproducibility of the rough genotype call).
    pub min_depth: u32,
    /// Lower edge of the het VAF band (inclusive).
    pub het_vaf_lo: f64,
    /// Upper edge of the het VAF band (inclusive).
    pub het_vaf_hi: f64,
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
    /// non-finite float, out-of-order VAF band, het exceeding variant
    /// sites).
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
        Ok(())
    }
}

impl HetCounts {
    fn validate(&self) -> Result<(), SampleSummaryError> {
        let bad = |field, reason: String| SampleSummaryError::InvalidField { field, reason };
        if self.n_het_sites > self.n_variant_sites {
            return Err(bad(
                "heterozygosity.n-het-sites",
                format!(
                    "{} exceeds n-variant-sites {}",
                    self.n_het_sites, self.n_variant_sites
                ),
            ));
        }
        for (field, v) in [
            ("heterozygosity.het-vaf-lo", self.het_vaf_lo),
            ("heterozygosity.het-vaf-hi", self.het_vaf_hi),
        ] {
            if !v.is_finite() || !(0.0..=1.0).contains(&v) {
                return Err(bad(field, format!("must be finite in [0, 1], got {v}")));
            }
        }
        if self.het_vaf_lo > self.het_vaf_hi {
            return Err(bad(
                "heterozygosity.het-vaf-lo",
                format!("{} exceeds het-vaf-hi {}", self.het_vaf_lo, self.het_vaf_hi),
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
                // gc_bins (2) * (depth_bins + 1 = 4) = 8 cells.
                counts: vec![0, 1, 2, 3, 4, 5, 6, 7],
            },
            heterozygosity: HetCounts {
                n_het_sites: 100,
                n_variant_sites: 250,
                min_depth: 4,
                het_vaf_lo: 0.3,
                het_vaf_hi: 0.7,
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
            "[heterozygosity]",
            "n-het-sites",
            "n-variant-sites",
            "het-vaf-lo",
        ] {
            assert!(body.contains(key), "wire key {key:?} absent from:\n{body}");
        }
        // snake_case must not leak.
        for forbidden in ["window_bp", "gc_bins", "n_het_sites"] {
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
    fn rejects_inverted_vaf_band() {
        let mut s = sample();
        s.heterozygosity.het_vaf_lo = 0.8;
        s.heterozygosity.het_vaf_hi = 0.6;
        let err = s.to_toml_bytes().expect_err("inverted band must fail");
        assert!(
            matches!(err, SampleSummaryError::InvalidField { field, .. } if field == "heterozygosity.het-vaf-lo"),
            "got {err:?}"
        );
    }

    #[test]
    fn rejects_het_exceeding_variant_sites() {
        let mut s = sample();
        s.heterozygosity.n_het_sites = s.heterozygosity.n_variant_sites + 1;
        let err = s.to_toml_bytes().expect_err("het > variant must fail");
        assert!(
            matches!(err, SampleSummaryError::InvalidField { field, .. } if field == "heterozygosity.n-het-sites"),
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
            // Make n-het-sites (100) exceed n-variant-sites (250 -> 50).
            .replace("n-variant-sites = 250", "n-variant-sites = 50");
        let err = SampleSummary::from_toml_bytes(body.as_bytes())
            .expect_err("value invariant must be enforced on parse");
        assert!(
            matches!(err, SampleSummaryError::InvalidField { field, .. } if field == "heterozygosity.n-het-sites"),
            "got {err:?}"
        );
    }

    proptest::proptest! {
        /// Round-trip over random valid documents, including boundary values
        /// (`gc_bins`/`depth_bins` at 1, `u32::MAX` counts, `lo == hi`,
        /// `min_depth = 0`). Catches float-format / large-integer drift and
        /// kebab-case key collisions a single hand-built example misses.
        #[test]
        fn round_trips_arbitrary_valid_summary(
            gc_bins in 1u32..6,
            depth_bins in 1u32..6,
            depth_bin_width in 0.001f64..1e6,
            lo in 0.0f64..=1.0,
            band in 0.0f64..=1.0,
            n_het in 0u64..1000,
            extra in 0u64..1000,
            window_bp in 1u32..100_000,
        ) {
            let het_vaf_lo = lo;
            let het_vaf_hi = (lo + band).min(1.0);
            let cells = (gc_bins as usize) * (depth_bins as usize + 1);
            let s = SampleSummary {
                version: SAMPLE_SUMMARY_VERSION,
                coverage_by_gc: CoverageByGcHistogram {
                    window_bp,
                    gc_bins,
                    depth_bin_width,
                    depth_bins,
                    n_tiles: u64::from(n_het) + 7,
                    n_skipped_tiles: 3,
                    counts: vec![u32::MAX; cells],
                },
                heterozygosity: HetCounts {
                    n_het_sites: n_het,
                    n_variant_sites: n_het + extra,
                    min_depth: 0,
                    het_vaf_lo,
                    het_vaf_hi,
                },
            };
            let bytes = s.to_toml_bytes().expect("serialise");
            let back = SampleSummary::from_toml_bytes(&bytes).expect("parse");
            proptest::prop_assert_eq!(back, s);
        }
    }
}
