//! TRF-mod invocation + BED parsing (architecture
//! [`ssr_catalog.md`](../../../doc/devel/architecture/ssr_catalog.md) §3).
//!
//! **This increment: the parsed-record type only.** [`TrfRecord`] is the
//! reduced BED row that [`super::postprocess`] consumes. The locate / spawn /
//! parse functions (`locate_trf_mod`, `run_on_contig`, `version`,
//! `parse_bed_line`) land in the next increment, now that the `trf-mod` binary
//! is available in the dev container.
//!
//! trf-mod emits a tab-separated BED with 10 columns:
//! `ctg  start  end  period  copyNum  fracMatch  fracGap  score  entropy
//! pattern`. We keep `start`/`end` (0-based half-open), the authoritative
//! `period`, `fracMatch` (a sanity field — purity is recomputed from the tract,
//! §4), `score` (an early accept-gate), and the consensus `pattern` (kept for
//! sanity; the catalog motif comes from the reference tract, §3). `copyNum` /
//! `fracGap` / `entropy` are parsed-and-dropped.

/// One parsed TRF-mod BED row, reduced to what the catalog needs.
/// Coordinates are **0-based half-open** (`[start, end)`), as trf-mod emits.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TrfRecord {
    /// Tract start (0-based, inclusive).
    pub(crate) start: u32,
    /// Tract end (0-based, exclusive).
    pub(crate) end: u32,
    /// Authoritative repeat period (motif length in bases). May differ from
    /// `pattern.len()` (architecture §3), so the catalog motif is taken as the
    /// first `period` bases of the tract, not from `pattern`.
    pub(crate) period: u16,
    /// TRF percent-match fraction in `[0, 1]` — a sanity field only; the
    /// catalog `purity_fraction` is recomputed from the trimmed tract (§4).
    pub(crate) frac_match: f32,
    /// TRF alignment score — used as an early accept-gate.
    pub(crate) score: i32,
    /// TRF consensus pattern (sanity only; not the catalog motif).
    pub(crate) pattern: Box<[u8]>,
}

#[cfg(test)]
impl TrfRecord {
    /// Test constructor — builds a record from the fields `postprocess`
    /// reads. (`frac_match` is unused by the pipeline, so it is fixed to 1.0.)
    pub(crate) fn for_test(start: u32, end: u32, period: u16, score: i32, pattern: &[u8]) -> Self {
        Self {
            start,
            end,
            period,
            frac_match: 1.0,
            score,
            pattern: pattern.to_vec().into_boxed_slice(),
        }
    }
}
