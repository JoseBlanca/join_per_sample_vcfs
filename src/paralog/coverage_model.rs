//! The per-sample single-copy coverage model (Q2).
//!
//! Answers "at a window of this GC content, what read depth would *one
//! copy* produce in this sample?" by fitting three quantities from the
//! sample's stored [`CoverageByGcHistogram`]:
//!
//! - `single_copy_scale` — the sample's typical one-copy depth level, taken
//!   as the **global mode** of the covered-window depth distribution (arch
//!   Premise 1: the mode marks the one-copy peak without assuming
//!   single-copy windows are an outright majority, as the prototype's median
//!   does). A `mode/median ∉ [0.5, 1.5]` sanity guard rejects the fit when
//!   the mode landed on a wrong copy-number peak.
//! - `gc_bias_curve` — a per-GC multiplier, the weighted median of the
//!   relative depth over the single-copy band, gap-filled and smoothed.
//! - `single_copy_depth_sd` (σ₀) — the robust spread (`1.4826·MAD`) of a
//!   one-copy window's relative depth, fit from the data (≈ 0.28 on
//!   tomato2), feeding the Normal coverage term of the H1/H2 likelihood.
//!
//! Their product `single_copy_scale · gc_bias_curve(GC)` is the expected
//! single-copy depth; dividing an observed window depth by it gives the
//! **relative copy number** (`1.0 = one copy`). The `(GC, depth)` at a
//! locus is measured fresh from the `.psp` body at scoring time — the
//! histogram is only the training data.
//!
//! Faithful in mechanism to the prototype
//! `benchmarks/tomato2/src/build_gc_normalization.py`, adapted to fit from
//! the binned 2-D histogram rather than the raw per-window table, and to
//! anchor on the mode rather than the median (arch Premise 1, spec §4).

use thiserror::Error;

use crate::sample_summary::CoverageByGcHistogram;

/// Default minimum tile count for a GC bin's curve value to be trusted;
/// thinner bins are gap-filled from neighbours (prototype `>= 50`).
pub const DEFAULT_MIN_BIN_COUNT: u64 = 50;

/// Default lower edge of the single-copy relative-depth band. Windows below
/// it are deletions / low-mappability edges, excluded from the GC-curve and
/// σ₀ fits (prototype `rel > 0.4`).
pub const DEFAULT_SINGLE_COPY_LO: f64 = 0.4;

/// Default upper edge of the single-copy relative-depth band. Windows above
/// it are duplications, excluded from the fits (prototype `rel < 1.6`).
pub const DEFAULT_SINGLE_COPY_HI: f64 = 1.6;

/// Default half-not-quite: the GC-curve smoothing window is this many bins
/// wide (a centred median smooth; prototype smooths over `i-2..=i+2`).
pub const DEFAULT_SMOOTH_WINDOW: usize = 5;

/// Default `mode/median` acceptance bounds. A healthy sample sits at
/// `≈ 0.83–0.98` (measured on tomato2); outside `[0.5, 1.5]` the mode landed
/// on a wrong copy-number peak (e.g. a 2× peak → ratio ≈ 2), so the fit is
/// rejected rather than silently mis-scaling (arch Premise 1, F-9 guard).
pub const DEFAULT_MODE_MEDIAN_RATIO_LO: f64 = 0.5;
/// Upper `mode/median` acceptance bound. See [`DEFAULT_MODE_MEDIAN_RATIO_LO`].
pub const DEFAULT_MODE_MEDIAN_RATIO_HI: f64 = 1.5;

/// Default cap on the fraction of covered positions allowed in the depth
/// histogram's **overflow** bin before the fit is rejected as out-of-range.
///
/// The regular depth bins span `0 .. depth_bins·depth_bin_width` (the default
/// scheme = `0..100×`); everything above lands in a single overflow column.
/// The overflow column legitimately holds a small tail — the high-copy
/// (`≥ 2×`) duplication windows the filter is meant to catch — so a few
/// percent is expected. When the fraction is large, the sample's own
/// *single-copy* depth has reached the top of the range: the single-copy peak
/// itself overflows, the regular bins hold only noise, and `mode_depth`
/// anchors `single_copy_scale` on that noise (a single ~100× human sample fit
/// its scale to a 145-position stray bin at 76× while 99.9% of positions
/// overflowed — every ordinary variant then read as ~1.3 copies and was
/// flagged as a hidden paralog). Past this point the coverage signal cannot
/// separate single-copy from a collapsed paralog anyway — a 2× paralog and the
/// single-copy peak both overflow — so rejecting the fit (→ the sample is
/// carried absent, its coverage inert) is the honest outcome, not a
/// workaround. `0.20` = reject once more than a fifth of covered positions
/// exceed the histogram's depth range.
pub const DEFAULT_MAX_OVERFLOW_FRACTION: f64 = 0.20;

/// Acceptance bounds on the `mode/median` depth ratio: the fit is rejected
/// unless `mode/median ∈ [lo, hi]`. A well-formed pair has `0 < lo <= hi`;
/// [`new`] enforces that at the construction boundary, and [`DEFAULT`] is the
/// tomato2-measured `[0.5, 1.5]` (a healthy sample sits at ≈ 0.83–0.98, so a
/// ratio near 2 means the mode landed on a 2× peak).
///
/// [`new`]: ModeMedianRatioBounds::new
/// [`DEFAULT`]: ModeMedianRatioBounds::DEFAULT
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModeMedianRatioBounds {
    /// Lower acceptance bound (`> 0`).
    pub lo: f64,
    /// Upper acceptance bound (`>= lo`).
    pub hi: f64,
}

impl ModeMedianRatioBounds {
    /// The default `[0.5, 1.5]` acceptance window (arch Premise 1, F-9 guard).
    pub const DEFAULT: Self = Self {
        lo: DEFAULT_MODE_MEDIAN_RATIO_LO,
        hi: DEFAULT_MODE_MEDIAN_RATIO_HI,
    };

    /// Build bounds, returning `None` for a degenerate window (`lo <= 0`,
    /// `lo > hi`, or a non-finite endpoint). Use for bounds derived from
    /// input; the trusted [`DEFAULT`] is a compile-time literal.
    ///
    /// [`DEFAULT`]: ModeMedianRatioBounds::DEFAULT
    pub fn new(lo: f64, hi: f64) -> Option<Self> {
        (lo.is_finite() && hi.is_finite() && lo > 0.0 && lo <= hi).then_some(Self { lo, hi })
    }

    /// Whether `ratio` is within `[lo, hi]` (inclusive).
    pub fn contains(&self, ratio: f64) -> bool {
        (self.lo..=self.hi).contains(&ratio)
    }
}

/// Tuning knobs for [`SingleCopyCoverageModel::fit`]. [`Default`] reproduces
/// the prototype's choices.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoverageFitConfig {
    /// Minimum tile count for a GC bin's curve value to be trusted.
    pub min_bin_count: u64,
    /// Lower edge of the single-copy relative-depth band `(lo, hi)`.
    pub single_copy_lo: f64,
    /// Upper edge of the single-copy relative-depth band `(lo, hi)`.
    pub single_copy_hi: f64,
    /// Width (in GC bins) of the centred median smooth applied to the curve.
    pub smooth_window: usize,
    /// If set, use this σ₀ instead of fitting it from the data (the
    /// per-sample MAD fit is the default; a cohort-wide constant override
    /// stays available — arch Premise 1).
    pub single_copy_depth_sd_override: Option<f64>,
    /// `mode/median` must fall within these bounds or the fit is rejected.
    pub mode_median_ratio_bounds: ModeMedianRatioBounds,
    /// Reject the fit when more than this fraction of covered positions land in
    /// the depth histogram's overflow bin (the sample's single-copy depth has
    /// exceeded the histogram range). See [`DEFAULT_MAX_OVERFLOW_FRACTION`].
    pub max_overflow_fraction: f64,
}

impl Default for CoverageFitConfig {
    fn default() -> Self {
        Self {
            min_bin_count: DEFAULT_MIN_BIN_COUNT,
            single_copy_lo: DEFAULT_SINGLE_COPY_LO,
            single_copy_hi: DEFAULT_SINGLE_COPY_HI,
            smooth_window: DEFAULT_SMOOTH_WINDOW,
            single_copy_depth_sd_override: None,
            mode_median_ratio_bounds: ModeMedianRatioBounds::DEFAULT,
            max_overflow_fraction: DEFAULT_MAX_OVERFLOW_FRACTION,
        }
    }
}

/// Why fitting a [`SingleCopyCoverageModel`] failed. Each variant marks a
/// sample whose depth distribution cannot be anchored, raised rather than
/// silently mis-scaling (arch Premise 1).
#[derive(Debug, Error, PartialEq)]
#[non_exhaustive]
pub enum CoverageModelError {
    /// The histogram has no folded tiles (or none in a regular depth bin),
    /// so there is no depth distribution to find a mode in.
    #[error("coverage histogram has no usable tiles (degenerate sample)")]
    NoTiles,
    /// The densest depth bin is the bottom bin, so the fitted single-copy
    /// scale is below one full depth bin — the sample is essentially
    /// uncovered and has no resolvable single-copy peak. Rejected rather
    /// than anchoring on a sub-bin scale that would inflate every window's
    /// relative copy number.
    #[error(
        "depth mode is the bottom bin: single-copy scale {scale:.3} is below one depth bin ({depth_bin_width:.3}) — sample too shallow to anchor"
    )]
    DepthModeAtBottomBin { scale: f64, depth_bin_width: f64 },
    /// `mode/median` fell outside the acceptance bounds — the mode landed on
    /// a wrong copy-number peak (a value near `2` means a 2× peak).
    #[error(
        "mode/median {ratio:.3} outside [{lo:.3}, {hi:.3}]: the depth mode landed on a wrong copy-number peak"
    )]
    ModeMedianRatioOutOfBounds { ratio: f64, lo: f64, hi: f64 },
    /// Too large a fraction of covered positions landed in the depth
    /// histogram's overflow bin: the sample's single-copy depth has reached
    /// the top of the histogram range, so the single-copy peak is unresolvable
    /// and the coverage signal can no longer separate single-copy from a
    /// collapsed paralog. Rejected rather than anchoring the scale on the noise
    /// left in the near-empty regular bins.
    #[error(
        "depth range exceeded: {overflow_fraction:.3} of covered positions are in the overflow bin (limit {limit:.3}) — single-copy depth is above the histogram range"
    )]
    DepthRangeExceeded { overflow_fraction: f64, limit: f64 },
    /// No GC bin had enough single-copy-band tiles to seed the GC curve, so
    /// it cannot be gap-filled from any anchor.
    #[error("no GC bin has enough single-copy tiles to fit the GC-bias curve")]
    NoGcCurveAnchor,
    /// A config knob was invalid (empty single-copy band, `smooth_window`
    /// `< 1`, or non-finite band edges).
    #[error("invalid coverage-fit config: {reason}")]
    InvalidConfig { reason: &'static str },
}

/// What one copy's read depth looks like in a single sample, as a function
/// of GC content. `relative_copy_number == 1.0` means one copy.
///
/// Fit with [`fit`]; query with [`relative_copy_number`]. Holds only the
/// three fitted quantities plus the GC binning needed to interpolate the
/// curve — the training histogram is not retained.
///
/// [`fit`]: SingleCopyCoverageModel::fit
/// [`relative_copy_number`]: SingleCopyCoverageModel::relative_copy_number
#[derive(Debug, Clone, PartialEq)]
pub struct SingleCopyCoverageModel {
    /// Analysis window width in bp, inherited from the histogram.
    window_bp: u32,
    /// The sample's typical one-copy depth level (the depth-distribution
    /// mode). Always `> 0`.
    single_copy_scale: f64,
    /// Number of GC bins the curve is indexed by (= the histogram's).
    gc_bins: u32,
    /// Per-GC-bin coverage multiplier, gap-filled + smoothed. Index `i`
    /// corresponds to GC-bin centre `(i + 0.5) / gc_bins`.
    gc_bias_curve: Vec<f64>,
    /// Robust SD of a one-copy window's relative depth (σ₀). Always `> 0`.
    single_copy_depth_sd: f64,
}

impl SingleCopyCoverageModel {
    /// Fit the model from a sample's stored coverage-by-GC histogram.
    ///
    /// # Errors
    ///
    /// [`CoverageModelError`] when the sample is degenerate (no tiles, a
    /// bottom-bin depth mode, a mode/median that fails the sanity guard, or
    /// no GC bin dense enough to seed the curve) or the config is invalid.
    pub fn fit(
        hist: &CoverageByGcHistogram,
        cfg: &CoverageFitConfig,
    ) -> Result<Self, CoverageModelError> {
        validate_config(cfg)?;
        let layout = HistogramLayout::new(hist);

        // 0. Depth-range guard: if the single-copy peak has overflowed the
        //    histogram range, the regular bins hold only noise and the mode
        //    anchor is meaningless. Reject before fitting anything.
        let depth_marginal = layout.depth_marginal(hist);
        let regular_total: u64 = depth_marginal.iter().sum();
        let overflow_total = layout.overflow_total(hist);
        let covered_total = regular_total + overflow_total;
        if covered_total > 0 {
            let overflow_fraction = overflow_total as f64 / covered_total as f64;
            if overflow_fraction > cfg.max_overflow_fraction {
                return Err(CoverageModelError::DepthRangeExceeded {
                    overflow_fraction,
                    limit: cfg.max_overflow_fraction,
                });
            }
        }

        // 1. Single-copy scale = mode of the marginal depth distribution
        //    (sub-bin-refined), guarded against a wrong-peak landing.
        let single_copy_scale = mode_depth(&depth_marginal, hist.depth_bin_width)?;
        let median =
            median_depth(&layout, hist, &depth_marginal).ok_or(CoverageModelError::NoTiles)?;
        let ratio = single_copy_scale / median;
        let bounds = cfg.mode_median_ratio_bounds;
        if !bounds.contains(ratio) {
            return Err(CoverageModelError::ModeMedianRatioOutOfBounds {
                ratio,
                lo: bounds.lo,
                hi: bounds.hi,
            });
        }

        // 2. GC-bias curve = per-GC-bin weighted median of rel over the
        //    single-copy band, gap-filled from neighbours, then smoothed.
        let gc_bias_curve = fit_gc_curve(hist, &layout, single_copy_scale, cfg)?;

        // 3. σ₀ = the override, or the robust SD of the GC-corrected
        //    relative copy number over single-copy-band tiles.
        let single_copy_depth_sd = match cfg.single_copy_depth_sd_override {
            Some(sd) => sd,
            None => fit_sigma0(hist, &layout, single_copy_scale, &gc_bias_curve, cfg),
        };

        Ok(Self {
            window_bp: hist.window_bp,
            single_copy_scale,
            gc_bins: hist.gc_bins,
            gc_bias_curve,
            single_copy_depth_sd,
        })
    }

    /// Expected depth of a single-copy window at this GC content:
    /// `single_copy_scale · gc_bias_curve(gc_fraction)`. `gc_fraction` is
    /// clamped to `[0, 1]`.
    pub fn expected_single_copy_depth(&self, gc_fraction: f64) -> f64 {
        self.single_copy_scale * self.gc_multiplier(gc_fraction)
    }

    /// Observed window depth expressed in copies (`1.0 = one copy`):
    /// `observed_depth / expected_single_copy_depth(gc_fraction)`.
    pub fn relative_copy_number(&self, gc_fraction: f64, observed_depth: f64) -> f64 {
        observed_depth / self.expected_single_copy_depth(gc_fraction)
    }

    /// The fitted (or overridden) σ₀ — the one-copy relative-depth SD.
    pub fn single_copy_depth_sd(&self) -> f64 {
        self.single_copy_depth_sd
    }

    /// The sample's single-copy depth level (the depth-distribution mode).
    pub fn single_copy_scale(&self) -> f64 {
        self.single_copy_scale
    }

    /// The analysis window width in bp.
    pub fn window_bp(&self) -> u32 {
        self.window_bp
    }

    /// The GC-bias multiplier at `gc_fraction`, linearly interpolated
    /// between GC-bin centres and clamped at the ends.
    fn gc_multiplier(&self, gc_fraction: f64) -> f64 {
        let gc_bins = self.gc_bias_curve.len();
        debug_assert!(gc_bins >= 1, "curve has at least one bin after a fit");
        let g = gc_fraction.clamp(0.0, 1.0);
        // GC-bin centre for index i is (i + 0.5) / gc_bins. Position g on
        // that centre axis, then linearly interpolate between the two
        // straddling centres (clamped to the end values beyond the extremes).
        let x = g * gc_bins as f64 - 0.5;
        if x <= 0.0 {
            return self.gc_bias_curve[0];
        }
        if x >= (gc_bins - 1) as f64 {
            return self.gc_bias_curve[gc_bins - 1];
        }
        let lo = x.floor() as usize;
        let frac = x - lo as f64;
        self.gc_bias_curve[lo] * (1.0 - frac) + self.gc_bias_curve[lo + 1] * frac
    }
}

/// Static geometry of a [`CoverageByGcHistogram`]: bin counts and the
/// row-stride, so cell lookups and bin-centre maths don't re-derive them.
struct HistogramLayout {
    gc_bins: usize,
    depth_bins: usize,
    /// Cells per GC row = regular depth bins + one overflow bin.
    row_stride: usize,
}

impl HistogramLayout {
    fn new(hist: &CoverageByGcHistogram) -> Self {
        let depth_bins = hist.depth_bins as usize;
        Self {
            gc_bins: hist.gc_bins as usize,
            depth_bins,
            row_stride: depth_bins + 1,
        }
    }

    /// Count in cell `(gc_bin, depth_bin)`; `depth_bin == depth_bins` is the
    /// overflow column.
    fn cell(&self, hist: &CoverageByGcHistogram, gc_bin: usize, depth_bin: usize) -> u64 {
        u64::from(hist.counts[gc_bin * self.row_stride + depth_bin])
    }

    /// Total count in the overflow column across all GC bins (positions whose
    /// mean depth exceeded the histogram's regular range).
    fn overflow_total(&self, hist: &CoverageByGcHistogram) -> u64 {
        (0..self.gc_bins)
            .map(|gc_bin| self.cell(hist, gc_bin, self.depth_bins))
            .sum()
    }

    /// Marginal depth histogram over the **regular** bins (overflow
    /// excluded): `out[d] = Σ_gc counts[gc][d]`.
    fn depth_marginal(&self, hist: &CoverageByGcHistogram) -> Vec<u64> {
        let mut out = vec![0u64; self.depth_bins];
        for gc_bin in 0..self.gc_bins {
            for (depth_bin, slot) in out.iter_mut().enumerate() {
                *slot += self.cell(hist, gc_bin, depth_bin);
            }
        }
        out
    }
}

/// Representative depth of regular depth bin `d`: the bin centre.
fn depth_bin_center(d: usize, depth_bin_width: f64) -> f64 {
    (d as f64 + 0.5) * depth_bin_width
}

/// Validate the fit config's shape-level invariants up front.
fn validate_config(cfg: &CoverageFitConfig) -> Result<(), CoverageModelError> {
    let bad = |reason| Err(CoverageModelError::InvalidConfig { reason });
    if !(cfg.single_copy_lo.is_finite() && cfg.single_copy_hi.is_finite()) {
        return bad("single-copy band edges must be finite");
    }
    if cfg.single_copy_lo >= cfg.single_copy_hi {
        return bad("single-copy band is empty (lo >= hi)");
    }
    if cfg.single_copy_lo < 0.0 {
        return bad("single-copy band lower edge is negative");
    }
    if cfg.smooth_window < 1 {
        return bad("smooth_window must be >= 1");
    }
    // A `ModeMedianRatioBounds` built via `new` is already well-formed, but
    // the pub fields let a literal set a degenerate window, so re-check here.
    let bounds = cfg.mode_median_ratio_bounds;
    if !(bounds.lo.is_finite()
        && bounds.hi.is_finite()
        && bounds.lo > 0.0
        && bounds.lo <= bounds.hi)
    {
        return bad("mode/median bounds must be finite with 0 < lo <= hi");
    }
    Ok(())
}

/// The single-copy scale: the mode of the marginal depth distribution,
/// refined to sub-bin precision by parabolic interpolation around the
/// densest regular bin. Rejects an all-empty distribution ([`NoTiles`]) or a
/// mode at the bottom bin ([`DepthModeAtBottomBin`]).
///
/// On a tie between two equally-dense bins, `max_by_key` keeps the *last*
/// (higher-depth) bin — a deterministic, intentional convention (a real
/// single-copy peak biases high against the low-mappability shoulder), so a
/// future refactor must preserve it.
///
/// [`NoTiles`]: CoverageModelError::NoTiles
/// [`DepthModeAtBottomBin`]: CoverageModelError::DepthModeAtBottomBin
fn mode_depth(depth_marginal: &[u64], depth_bin_width: f64) -> Result<f64, CoverageModelError> {
    let peak = (0..depth_marginal.len())
        .filter(|&i| depth_marginal[i] > 0)
        .max_by_key(|&i| depth_marginal[i])
        .ok_or(CoverageModelError::NoTiles)?;

    // A mode in the bottom depth bin means the single-copy scale would be
    // below one full bin — an essentially-uncovered sample. Reject it rather
    // than anchor on a sub-bin scale (which would read every window as
    // high-copy). The parabolic branch below can only *raise* the scale, so
    // this is the sole sub-one-bin case.
    if peak == 0 {
        return Err(CoverageModelError::DepthModeAtBottomBin {
            scale: 0.5 * depth_bin_width,
            depth_bin_width,
        });
    }

    // Parabolic sub-bin refinement using the peak and its two neighbours:
    // the vertex offset of the parabola through (−1, c₋), (0, c₀), (+1, c₊).
    // Only when the triple is concave (a real interior peak); otherwise the
    // bin centre. Offset is clamped to ±0.5 so it never leaves the bin.
    let offset = if peak + 1 < depth_marginal.len() {
        let cm = depth_marginal[peak - 1] as f64;
        let c0 = depth_marginal[peak] as f64;
        let cp = depth_marginal[peak + 1] as f64;
        let denom = cm - 2.0 * c0 + cp;
        if denom < 0.0 {
            (0.5 * (cm - cp) / denom).clamp(-0.5, 0.5)
        } else {
            0.0
        }
    } else {
        0.0
    };

    // peak >= 1 and offset >= -0.5, so this is always >= 1.0 * width > 0.
    Ok((peak as f64 + 0.5 + offset) * depth_bin_width)
}

/// The median depth over all folded tiles, including the overflow column
/// (assigned its lower edge, a rank-safe underestimate since the median
/// lands in a regular bin unless > 50 % of tiles overflow). Takes the
/// already-computed regular-bin `marginal` so the matrix is not swept twice.
/// `None` only if there are no tiles at all.
fn median_depth(
    layout: &HistogramLayout,
    hist: &CoverageByGcHistogram,
    marginal: &[u64],
) -> Option<f64> {
    let mut bins: Vec<(f64, u64)> = Vec::with_capacity(layout.depth_bins + 1);
    for (d, &count) in marginal.iter().enumerate() {
        if count > 0 {
            bins.push((depth_bin_center(d, hist.depth_bin_width), count));
        }
    }
    // Overflow column: total across GC rows, placed at its lower edge.
    let overflow: u64 = (0..layout.gc_bins)
        .map(|gc_bin| layout.cell(hist, gc_bin, layout.depth_bins))
        .sum();
    if overflow > 0 {
        bins.push((layout.depth_bins as f64 * hist.depth_bin_width, overflow));
    }
    weighted_median(&mut bins)
}

/// The weighted median of `(value, weight)` pairs: the smallest value whose
/// cumulative weight reaches half the total. `None` if the total weight is
/// zero. Sorts `pairs` in place by value.
fn weighted_median(pairs: &mut [(f64, u64)]) -> Option<f64> {
    let total: u128 = pairs.iter().map(|&(_, w)| u128::from(w)).sum();
    if total == 0 {
        return None;
    }
    pairs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let half = total.div_ceil(2); // first value with cumulative >= ceil(total/2)
    let mut acc: u128 = 0;
    for &(value, weight) in pairs.iter() {
        acc += u128::from(weight);
        if acc >= half {
            return Some(value);
        }
    }
    // Total > 0 guarantees the loop returns; unreachable in practice.
    pairs.last().map(|&(v, _)| v)
}

/// Fit the per-GC-bin bias curve: for each GC bin, the weighted median of
/// `rel = depth / scale` over the single-copy band; sparse bins gap-filled
/// by linear interpolation from anchors, then a centred median smooth.
fn fit_gc_curve(
    hist: &CoverageByGcHistogram,
    layout: &HistogramLayout,
    scale: f64,
    cfg: &CoverageFitConfig,
) -> Result<Vec<f64>, CoverageModelError> {
    // Raw per-GC-bin median rel; NaN where the band count is too thin.
    let mut raw = vec![f64::NAN; layout.gc_bins];
    let mut scratch: Vec<(f64, u64)> = Vec::with_capacity(layout.depth_bins);
    for (gc_bin, slot) in raw.iter_mut().enumerate() {
        scratch.clear();
        let mut band_count: u64 = 0;
        for depth_bin in 0..layout.depth_bins {
            let count = layout.cell(hist, gc_bin, depth_bin);
            if count == 0 {
                continue;
            }
            let rel = depth_bin_center(depth_bin, hist.depth_bin_width) / scale;
            if rel > cfg.single_copy_lo && rel < cfg.single_copy_hi {
                scratch.push((rel, count));
                band_count += count;
            }
        }
        if band_count >= cfg.min_bin_count
            && let Some(m) = weighted_median(&mut scratch)
        {
            *slot = m;
        }
    }

    gap_fill(&mut raw).ok_or(CoverageModelError::NoGcCurveAnchor)?;
    Ok(smooth_median(&raw, cfg.smooth_window))
}

/// Fill `NaN` entries in place by linear interpolation between the nearest
/// valid anchors, extending the end anchors flat. Returns `None` (leaving
/// `curve` untouched) if there is no valid anchor at all.
fn gap_fill(curve: &mut [f64]) -> Option<()> {
    let anchors: Vec<usize> = (0..curve.len()).filter(|&i| curve[i].is_finite()).collect();
    let first = *anchors.first()?;
    // PANIC-FREE: `first()?` returned above, so `anchors` is non-empty.
    let last = *anchors.last().unwrap();
    // Flat-extend before the first / after the last anchor.
    for i in 0..first {
        curve[i] = curve[first];
    }
    for i in (last + 1)..curve.len() {
        curve[i] = curve[last];
    }
    // Linearly interpolate each interior gap between consecutive anchors.
    for pair in anchors.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if b == a + 1 {
            continue;
        }
        let (va, vb) = (curve[a], curve[b]);
        for (offset, slot) in curve[(a + 1)..b].iter_mut().enumerate() {
            let frac = (offset + 1) as f64 / (b - a) as f64;
            *slot = va * (1.0 - frac) + vb * frac;
        }
    }
    Some(())
}

/// A centred median smooth with the given window width: each output is the
/// median of the `window`-wide neighbourhood clamped to the ends. The two
/// central values are averaged when the neighbourhood is even-length (which
/// the clamped ends always are for an odd `window` — matching the
/// prototype's `np.median`). A `window` of 1 is the identity.
fn smooth_median(curve: &[f64], window: usize) -> Vec<f64> {
    let half = window / 2;
    let n = curve.len();
    let mut scratch: Vec<f64> = Vec::with_capacity(window);
    (0..n)
        .map(|i| {
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(n);
            scratch.clear();
            scratch.extend_from_slice(&curve[lo..hi]);
            scratch.sort_unstable_by(f64::total_cmp);
            let mid = scratch.len() / 2;
            if scratch.len().is_multiple_of(2) {
                0.5 * (scratch[mid - 1] + scratch[mid])
            } else {
                scratch[mid]
            }
        })
        .collect()
}

/// σ₀ = `1.4826 · MAD` of the GC-corrected relative copy number over
/// single-copy-band tiles.
///
/// Tiles are selected by the **un-corrected** `rel` band (matching the
/// prototype, which filters `rel` then measures spread on the corrected
/// value), and the MAD is taken over the corrected `rcn = rel / mult`.
///
/// The result is floored at the **sampling** spread — the relative SD a
/// one-copy window's depth would have from read-sampling alone, `√scale /
/// scale = 1/√scale` (a count with mean `scale` reads has Poisson SD
/// `√scale`). A one-copy window cannot be *tighter* than this, so it is the
/// honest lower bound when the band's MAD is uninformative: a band that
/// collapses to a single dominant depth bin has weighted MAD 0, and the true
/// spread is then unknowable from the histogram.
///
/// The earlier floor used the histogram's own *resolution* (`0.5 ·
/// depth_bin_width / scale`) instead. That shrinks as `1/scale`, so at deep
/// coverage — where the band MAD does collapse — it drove σ₀ to ~1e-3
/// (≈ `0.25/100` at 100× with the 0.5-wide bins), some 30× under the true
/// ~`1/√scale` ≈ 0.1 single-copy spread. The `Normal(1, σ₀)` coverage term
/// then became absurdly overconfident: an ordinary variant window scored
/// dozens of σ from single-copy and was flagged as a hidden paralog, so the
/// filter deleted almost the whole callset of any single high-depth sample
/// whose fit slipped past the mode/median guard (GIAB HG004 @≈100×: 725 → 106
/// calls, recall 0.98 → 0.02). The sampling floor is depth-appropriate at every
/// tier and keeps σ₀ strictly positive so the Normal coverage term never spikes.
fn fit_sigma0(
    hist: &CoverageByGcHistogram,
    layout: &HistogramLayout,
    scale: f64,
    curve: &[f64],
    cfg: &CoverageFitConfig,
) -> f64 {
    debug_assert_eq!(curve.len(), layout.gc_bins, "curve is fit per GC bin");
    // Read-sampling (Poisson) relative-depth SD: the tightest a one-copy
    // window can plausibly be, hence the honest σ₀ floor. `scale > 0` is an
    // invariant of the fit (mode-anchored, bottom-bin rejected).
    let sampling_floor = 1.0 / scale.sqrt();

    // Collect (relative_copy_number, weight) over single-copy-band tiles.
    // `curve` is index-aligned with the GC bins (asserted above).
    let mut vals: Vec<(f64, u64)> = Vec::new();
    for (gc_bin, &mult) in curve.iter().enumerate() {
        if mult <= 0.0 {
            continue;
        }
        for depth_bin in 0..layout.depth_bins {
            let count = layout.cell(hist, gc_bin, depth_bin);
            if count == 0 {
                continue;
            }
            let rel = depth_bin_center(depth_bin, hist.depth_bin_width) / scale;
            if rel > cfg.single_copy_lo && rel < cfg.single_copy_hi {
                vals.push((rel / mult, count));
            }
        }
    }
    let Some(center) = weighted_median(&mut vals) else {
        return sampling_floor;
    };
    let mut devs: Vec<(f64, u64)> = vals.iter().map(|&(v, w)| ((v - center).abs(), w)).collect();
    let mad = weighted_median(&mut devs).unwrap_or(0.0);
    (1.4826 * mad).max(sampling_floor)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a histogram by binning `(gc_fraction, depth, count)` tiles with
    /// the same cell arithmetic the accumulator uses. `callable_positions`
    /// is set to the tile total (`>= n_positions`, satisfying the invariant).
    fn hist(
        gc_bins: u32,
        depth_bins: u32,
        width: f64,
        window_bp: u32,
        tiles: &[(f64, f64, u32)],
    ) -> CoverageByGcHistogram {
        let row_stride = depth_bins as usize + 1;
        let mut counts = vec![0u32; gc_bins as usize * row_stride];
        let mut n = 0u64;
        for &(gc, depth, count) in tiles {
            let gb = ((gc * gc_bins as f64) as usize).min(gc_bins as usize - 1);
            let db = ((depth / width) as usize).min(depth_bins as usize);
            counts[gb * row_stride + db] += count;
            n += u64::from(count);
        }
        CoverageByGcHistogram {
            window_bp,
            gc_bins,
            depth_bin_width: width,
            depth_bins,
            n_positions: n,
            n_skipped_tiles: 0,
            callable_positions: n,
            counts,
        }
    }

    fn approx(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() <= eps
    }

    /// A GC-flat single-copy sample: the scale anchors on the mode (bin
    /// centre 10.5), `relative_copy_number(scale) == 1.0` by construction,
    /// and a window at twice the scale reads as two copies.
    #[test]
    fn fits_scale_and_relative_copy_number() {
        // Symmetric single-copy spread around depth-bin 10 in both GC bins.
        let tiles = [
            (0.25, 9.5, 100),
            (0.25, 10.5, 300),
            (0.25, 11.5, 100),
            (0.75, 9.5, 100),
            (0.75, 10.5, 300),
            (0.75, 11.5, 100),
        ];
        let h = hist(2, 40, 1.0, 500, &tiles);
        let cfg = CoverageFitConfig {
            single_copy_depth_sd_override: Some(0.28),
            ..Default::default()
        };
        let model = SingleCopyCoverageModel::fit(&h, &cfg).expect("fit");

        assert!(
            approx(model.single_copy_scale(), 10.5, 1e-9),
            "{}",
            model.single_copy_scale()
        );
        // The anchor: one scale's worth of depth is exactly one copy.
        assert!(approx(
            model.relative_copy_number(0.25, model.single_copy_scale()),
            1.0,
            1e-9
        ));
        // Twice the expected single-copy depth → two copies.
        let two_copy_depth = 2.0 * model.expected_single_copy_depth(0.25);
        assert!(approx(
            model.relative_copy_number(0.25, two_copy_depth),
            2.0,
            1e-9
        ));
        assert_eq!(model.single_copy_depth_sd(), 0.28);
        assert_eq!(model.window_bp(), 500);
    }

    /// The fit is invariant to the histogram's **count scale**. The sliding
    /// window (M3) folds one sample per *covered position* rather than per
    /// *tile* — ~`window_bp`× more counts of the same shape (and correlated
    /// between neighbours). Because the fit reads only the density (the
    /// arg-max depth mode and the weighted-MAD σ₀), a `W`×-denser histogram of
    /// the same shape produces the identical `single_copy_scale` and σ₀. This
    /// is why the sliding-window switch needs no fit-code change: the absolute
    /// σ₀ shrinks only because the *distribution* is smoother, not because the
    /// denser counts bias the estimator.
    #[test]
    fn fit_is_invariant_to_count_scale() {
        // A single-copy peak at depth-bin 10 with a symmetric spread, in two GC
        // bins so the GC curve has anchors.
        // Counts high enough per GC bin (sum 100 > min_bin_count) that even the
        // "tile scale" is a valid fit; the shape is what matters.
        let shape = [
            (0.25, 8.5, 10u32),
            (0.25, 9.5, 20),
            (0.25, 10.5, 40),
            (0.25, 11.5, 20),
            (0.25, 12.5, 10),
            (0.75, 8.5, 10),
            (0.75, 9.5, 20),
            (0.75, 10.5, 40),
            (0.75, 11.5, 20),
            (0.75, 12.5, 10),
        ];
        // "Tile scale" vs "position scale": the same shape, 500× more counts.
        let dense: Vec<(f64, f64, u32)> =
            shape.iter().map(|&(gc, d, c)| (gc, d, c * 500)).collect();
        let cfg = CoverageFitConfig::default(); // MAD fit, no σ₀ override.

        let sparse =
            SingleCopyCoverageModel::fit(&hist(2, 40, 1.0, 500, &shape), &cfg).expect("sparse fit");
        let denser =
            SingleCopyCoverageModel::fit(&hist(2, 40, 1.0, 500, &dense), &cfg).expect("dense fit");

        assert!(
            approx(sparse.single_copy_scale(), denser.single_copy_scale(), 1e-9),
            "mode invariant to count scale: {} vs {}",
            sparse.single_copy_scale(),
            denser.single_copy_scale(),
        );
        assert!(
            approx(
                sparse.single_copy_depth_sd(),
                denser.single_copy_depth_sd(),
                1e-9
            ),
            "σ₀ invariant to count scale: {} vs {}",
            sparse.single_copy_depth_sd(),
            denser.single_copy_depth_sd(),
        );
        // The MAD fit produced a real, positive σ₀ (not the override path).
        assert!(sparse.single_copy_depth_sd() > 0.0);
    }

    /// A GC-dependent single-copy shift (low-GC windows read low, high-GC
    /// high) is absorbed by `gc_bias_curve`: both extremes normalise to one
    /// copy, and the expected single-copy depth rises with GC.
    #[test]
    fn gc_bias_curve_flattens_gc_shift() {
        let tiles = [
            (0.125, 9.5, 150),  // GC bin 0: reads low
            (0.375, 10.5, 150), // GC bin 1: neutral
            (0.625, 10.5, 150), // GC bin 2: neutral
            (0.875, 11.5, 150), // GC bin 3: reads high
        ];
        let h = hist(4, 30, 1.0, 500, &tiles);
        // Isolate the curve fit from the smoother (window 1 = identity).
        let cfg = CoverageFitConfig {
            smooth_window: 1,
            single_copy_depth_sd_override: Some(0.28),
            ..Default::default()
        };
        let model = SingleCopyCoverageModel::fit(&h, &cfg).expect("fit");

        assert!(approx(model.relative_copy_number(0.125, 9.5), 1.0, 1e-9));
        assert!(approx(model.relative_copy_number(0.375, 10.5), 1.0, 1e-9));
        assert!(approx(model.relative_copy_number(0.875, 11.5), 1.0, 1e-9));
        // The curve is doing real work: low-GC expects less depth than high-GC.
        assert!(model.expected_single_copy_depth(0.125) < model.expected_single_copy_depth(0.875));
    }

    /// σ₀ is `max(1.4826·MAD, 1/√scale)` — the robust MAD spread, floored at
    /// the read-sampling (Poisson) relative SD. This shallow five-bin spread
    /// (scale 10.5) has a MAD-based σ₀ of `1.4826·0.095 ≈ 0.141`, but a
    /// one-copy window at depth 10.5 cannot be tighter than `1/√10.5 ≈ 0.309`,
    /// so the sampling floor governs.
    #[test]
    fn sigma0_fit_from_spread() {
        let tiles = [
            (0.5, 8.5, 100),
            (0.5, 9.5, 100),
            (0.5, 10.5, 140), // the mode
            (0.5, 11.5, 100),
            (0.5, 12.5, 100),
        ];
        let h = hist(1, 30, 1.0, 500, &tiles);
        let model = SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default()).expect("fit");
        let mad_based = 1.4826 * (1.0 / 10.5);
        let sampling_floor = 1.0 / 10.5_f64.sqrt();
        assert!(
            sampling_floor > mad_based,
            "the floor governs at this depth"
        );
        assert!(
            approx(model.single_copy_depth_sd(), sampling_floor, 1e-6),
            "sigma0 = {}",
            model.single_copy_depth_sd()
        );
    }

    /// Regression for the high-depth σ₀ collapse: when the single-copy band
    /// collapses to one depth bin the MAD is 0, and the honest σ₀ is the
    /// sampling floor `1/√scale`, NOT the histogram-resolution floor
    /// `0.5·depth_bin_width/scale` (which shrinks as `1/scale` and drove σ₀ to
    /// ~1e-3 at deep coverage, making the H1 coverage term absurdly
    /// overconfident).
    #[test]
    fn sigma0_floors_at_sampling_spread_when_band_collapses() {
        // All mass in one depth bin at ~64× (mode == median, no overflow).
        let h = hist(
            2,
            200,
            0.5,
            500,
            &[(0.25, 64.25, 1000), (0.75, 64.25, 1000)],
        );
        let model = SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default()).expect("fit");
        let sampling_floor = 1.0 / 64.25_f64.sqrt(); // ≈ 0.1247
        let old_resolution_floor = 0.5 * 0.5 / 64.25; // ≈ 0.0039
        assert!(
            approx(model.single_copy_depth_sd(), sampling_floor, 1e-6),
            "sigma0 = {} (want sampling floor {sampling_floor})",
            model.single_copy_depth_sd()
        );
        assert!(
            model.single_copy_depth_sd() > 10.0 * old_resolution_floor,
            "sampling floor must dwarf the old resolution floor"
        );
    }

    /// The mode/median guard rejects a duplication-rich sample whose densest
    /// depth bin is the 2× peak while the median stays at single-copy.
    #[test]
    fn rejects_wrong_copy_number_peak() {
        let mut tiles = vec![];
        // Single-copy spread across bins 8..=12 (the median lands here).
        for d in 8..=12 {
            tiles.push((0.25, d as f64 + 0.5, 60));
        }
        // A taller, concentrated 2× peak at depth 20 → the mode.
        tiles.push((0.25, 20.5, 100));
        let h = hist(2, 30, 1.0, 500, &tiles);
        let err = SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default())
            .expect_err("wrong-peak sample must be rejected");
        assert!(
            matches!(err, CoverageModelError::ModeMedianRatioOutOfBounds { .. }),
            "got {err:?}"
        );
    }

    /// The depth-range guard rejects a sample whose single-copy depth has
    /// exceeded the histogram range: most covered positions pile into the
    /// overflow bin, leaving only noise in the regular bins. (The ≈100× human
    /// GIAB samples hit exactly this — 99.9% overflow of the 0..100× range.)
    #[test]
    fn rejects_when_depth_range_exceeded() {
        // Range 0..10× (10 bins × 1.0). A little single-copy mass in-range,
        // most positions at ~50× → the overflow bin (fraction 0.8 > 0.2).
        let h = hist(2, 10, 1.0, 500, &[(0.5, 5.5, 100), (0.5, 50.0, 400)]);
        let err = SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default())
            .expect_err("out-of-range sample must be rejected");
        match err {
            CoverageModelError::DepthRangeExceeded {
                overflow_fraction,
                limit,
            } => {
                assert!(
                    approx(overflow_fraction, 0.8, 1e-9),
                    "got {overflow_fraction}"
                );
                assert!(approx(limit, DEFAULT_MAX_OVERFLOW_FRACTION, 1e-9));
            }
            other => panic!("got {other:?}"),
        }
    }

    /// A small overflow tail (the high-copy duplication windows the filter is
    /// meant to catch) is under the limit and does NOT trip the range guard.
    #[test]
    fn small_overflow_tail_is_accepted() {
        // ~9% overflow (10 of 110) around a single-copy peak at 5×.
        let tiles = [
            (0.25, 4.5, 200),
            (0.25, 5.5, 300),
            (0.25, 6.5, 200),
            (0.75, 4.5, 200),
            (0.75, 5.5, 300),
            (0.75, 6.5, 200),
            (0.25, 50.0, 70), // overflow
            (0.75, 50.0, 70), // overflow
        ];
        let h = hist(2, 10, 1.0, 500, &tiles);
        let err_is_range = matches!(
            SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default()),
            Err(CoverageModelError::DepthRangeExceeded { .. })
        );
        assert!(
            !err_is_range,
            "a small overflow tail must not trip the guard"
        );
    }

    /// An empty (all-`N` / no-tile) histogram has no depth distribution to
    /// anchor, so the fit fails with `NoTiles` rather than dividing by zero.
    #[test]
    fn rejects_empty_histogram() {
        let h = hist(2, 30, 1.0, 500, &[]);
        let err = SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default())
            .expect_err("empty histogram must be rejected");
        assert_eq!(err, CoverageModelError::NoTiles);
    }

    /// The σ₀ override bypasses the MAD fit exactly.
    #[test]
    fn sigma0_override_is_used() {
        let tiles = [(0.5, 10.5, 500)];
        let h = hist(1, 30, 1.0, 500, &tiles);
        let cfg = CoverageFitConfig {
            single_copy_depth_sd_override: Some(0.5),
            ..Default::default()
        };
        let model = SingleCopyCoverageModel::fit(&h, &cfg).expect("fit");
        assert_eq!(model.single_copy_depth_sd(), 0.5);
    }

    /// An empty single-copy band (`lo >= hi`) is a config error, caught
    /// before any fitting.
    #[test]
    fn rejects_invalid_config() {
        let h = hist(1, 30, 1.0, 500, &[(0.5, 10.5, 500)]);
        let cfg = CoverageFitConfig {
            single_copy_lo: 1.6,
            single_copy_hi: 0.4,
            ..Default::default()
        };
        let err = SingleCopyCoverageModel::fit(&h, &cfg).expect_err("empty band must fail");
        assert!(
            matches!(err, CoverageModelError::InvalidConfig { .. }),
            "got {err:?}"
        );
    }

    /// `weighted_median` returns the first value whose cumulative weight
    /// reaches half the total, and `None` on empty/zero-weight input.
    #[test]
    fn weighted_median_basics() {
        assert_eq!(weighted_median(&mut []), None);
        assert_eq!(weighted_median(&mut [(5.0, 0)]), None);
        // Weights 1,1,8: half = 5, crossed at the third value.
        assert_eq!(
            weighted_median(&mut [(1.0, 1), (2.0, 1), (3.0, 8)]),
            Some(3.0)
        );
        // Symmetric 1,1: half = 1, crossed at the first value.
        assert_eq!(weighted_median(&mut [(1.0, 1), (2.0, 1)]), Some(1.0));
    }

    /// `gap_fill` linearly interpolates interior gaps and flat-extends the
    /// ends; a curve with no valid anchor returns `None`.
    #[test]
    fn gap_fill_interpolates_and_extends() {
        let mut c = [f64::NAN, 1.0, f64::NAN, f64::NAN, 4.0, f64::NAN];
        gap_fill(&mut c).expect("has anchors");
        assert_eq!(c, [1.0, 1.0, 2.0, 3.0, 4.0, 4.0]);

        let mut all_nan = [f64::NAN, f64::NAN];
        assert!(gap_fill(&mut all_nan).is_none());
    }

    /// A single anchor fills the whole curve flat (the `first == last` path:
    /// no interior interpolation, both ends flat-extended).
    #[test]
    fn gap_fill_single_anchor_fills_flat() {
        let mut c = [f64::NAN, f64::NAN, 2.5, f64::NAN];
        gap_fill(&mut c).expect("one anchor");
        assert_eq!(c, [2.5, 2.5, 2.5, 2.5]);
    }

    /// A bottom-bin-only sample (essentially uncovered) is rejected rather
    /// than anchored on a sub-one-bin scale — the reachable
    /// `DepthModeAtBottomBin` guard (the old dead `ZeroScale` never fired).
    #[test]
    fn rejects_bottom_depth_bin_only_sample() {
        // All tiles in depth bin 0 (depth 0.5·width), so the mode is peak 0.
        let h = hist(1, 30, 2.0, 500, &[(0.5, 0.5, 500)]);
        let err = SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default())
            .expect_err("bottom-bin sample must be rejected");
        assert!(
            matches!(err, CoverageModelError::DepthModeAtBottomBin { .. }),
            "got {err:?}"
        );
    }

    /// When the single-copy band collapses to one depth bin the MAD is 0, so
    /// σ₀ falls back to the read-sampling (Poisson) floor `1/√scale` — the
    /// honest single-copy spread, not the depth-shrinking histogram resolution
    /// (see [`sigma0_floors_at_sampling_spread_when_band_collapses`] for the
    /// high-depth regime where the two diverge catastrophically).
    #[test]
    fn sigma0_floors_at_sampling_spread_when_band_collapses_low_depth() {
        let h = hist(1, 30, 1.0, 500, &[(0.5, 10.5, 500)]);
        let model = SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default()).expect("fit");
        // scale = 10.5 → sampling floor = 1/√10.5 ≈ 0.309.
        assert!(
            approx(model.single_copy_depth_sd(), 1.0 / 10.5_f64.sqrt(), 1e-12),
            "sigma0 = {}",
            model.single_copy_depth_sd()
        );
    }

    /// Every GC bin below `min_bin_count` leaves the curve with no anchor, so
    /// the fit fails with `NoGcCurveAnchor` (exercising the path at the `fit`
    /// level, not just `gap_fill` in isolation).
    #[test]
    fn rejects_when_all_gc_bins_thin() {
        // 10 in-band tiles per GC bin, below the default min_bin_count 50.
        let h = hist(2, 30, 1.0, 500, &[(0.25, 10.5, 10), (0.75, 10.5, 10)]);
        let err = SingleCopyCoverageModel::fit(&h, &CoverageFitConfig::default())
            .expect_err("all-thin curve must fail");
        assert_eq!(err, CoverageModelError::NoGcCurveAnchor);
    }

    /// Each `InvalidConfig` sub-check fires: non-finite band edge, zero
    /// smooth window, and degenerate mode/median bounds (both `lo == 0` and
    /// `lo > hi`, constructed as literals to bypass `new`'s validation).
    #[test]
    fn rejects_all_invalid_config_shapes() {
        let h = hist(1, 30, 1.0, 500, &[(0.5, 10.5, 500)]);
        let cfgs = [
            CoverageFitConfig {
                single_copy_lo: f64::NAN,
                ..Default::default()
            },
            CoverageFitConfig {
                smooth_window: 0,
                ..Default::default()
            },
            CoverageFitConfig {
                mode_median_ratio_bounds: ModeMedianRatioBounds { lo: 0.0, hi: 1.5 },
                ..Default::default()
            },
            CoverageFitConfig {
                mode_median_ratio_bounds: ModeMedianRatioBounds { lo: 1.5, hi: 0.5 },
                ..Default::default()
            },
        ];
        for cfg in cfgs {
            assert!(
                matches!(
                    SingleCopyCoverageModel::fit(&h, &cfg),
                    Err(CoverageModelError::InvalidConfig { .. })
                ),
                "config {cfg:?} should be rejected"
            );
        }
    }

    /// The parabolic refinement shifts the scale off the bin centre toward
    /// the heavier neighbour on an asymmetric peak (the happy-path test uses
    /// a symmetric triple where the offset is exactly 0).
    #[test]
    fn mode_depth_parabolic_refines_asymmetric_peak() {
        // cm=100, c0=300, cp=200: denom = -300, offset = +1/6.
        let tiles = [(0.5, 9.5, 100), (0.5, 10.5, 300), (0.5, 11.5, 200)];
        let h = hist(1, 40, 1.0, 500, &tiles);
        let cfg = CoverageFitConfig {
            single_copy_depth_sd_override: Some(0.28),
            ..Default::default()
        };
        let m = SingleCopyCoverageModel::fit(&h, &cfg).expect("fit");
        assert!(
            approx(m.single_copy_scale(), 10.5 + 1.0 / 6.0, 1e-9),
            "scale = {}",
            m.single_copy_scale()
        );
    }

    /// On a tie between two equally-dense depth bins the mode snaps to the
    /// higher-depth bin (the documented `max_by_key` last-max convention).
    #[test]
    fn mode_depth_breaks_ties_toward_higher_bin() {
        let tiles = [(0.5, 9.5, 100), (0.5, 11.5, 100)];
        let h = hist(1, 30, 1.0, 500, &tiles);
        let cfg = CoverageFitConfig {
            single_copy_depth_sd_override: Some(0.28),
            ..Default::default()
        };
        let m = SingleCopyCoverageModel::fit(&h, &cfg).expect("fit");
        assert!(
            approx(m.single_copy_scale(), 11.5, 1e-9),
            "scale = {}",
            m.single_copy_scale()
        );
    }

    /// `gc_multiplier` clamps beyond the GC-bin-centre axis rather than
    /// extrapolating — GC below the first / above the last centre reuses the
    /// end bins (checked through `expected_single_copy_depth`).
    #[test]
    fn gc_multiplier_clamps_beyond_bin_centres() {
        let tiles = [
            (0.125, 9.5, 150),
            (0.375, 10.5, 150),
            (0.625, 10.5, 150),
            (0.875, 11.5, 150),
        ];
        let h = hist(4, 30, 1.0, 500, &tiles);
        let cfg = CoverageFitConfig {
            smooth_window: 1,
            single_copy_depth_sd_override: Some(0.28),
            ..Default::default()
        };
        let m = SingleCopyCoverageModel::fit(&h, &cfg).expect("fit");
        assert_eq!(
            m.expected_single_copy_depth(-1.0),
            m.expected_single_copy_depth(0.0)
        );
        assert_eq!(
            m.expected_single_copy_depth(2.0),
            m.expected_single_copy_depth(1.0)
        );
    }

    /// `smooth_median` averages the two central values on an even-length
    /// neighbourhood (an upper-median would return 3.0 here, not 2.0).
    #[test]
    fn smooth_median_averages_even_neighbourhood() {
        assert_eq!(smooth_median(&[1.0, 3.0], 2), vec![2.0, 2.0]);
    }

    /// `weighted_median` takes the lower of the two middles on an even total
    /// (the `div_ceil(2)` crossing) — a pinned convention the curve and σ₀
    /// medians depend on.
    #[test]
    fn weighted_median_lower_middle_on_even_total() {
        // total 4, half = ceil(4/2) = 2: cumulative reaches 2 at value 2.0.
        assert_eq!(
            weighted_median(&mut [(1.0, 1), (2.0, 1), (3.0, 1), (4.0, 1)]),
            Some(2.0)
        );
    }

    /// `ModeMedianRatioBounds::new` accepts a well-formed window (incl. the
    /// `lo == hi` edge) and rejects degenerate ones.
    #[test]
    fn mode_median_bounds_new_validates() {
        assert!(ModeMedianRatioBounds::new(0.5, 1.5).is_some());
        assert!(ModeMedianRatioBounds::new(0.5, 0.5).is_some()); // lo == hi ok
        assert!(ModeMedianRatioBounds::new(1.5, 0.5).is_none()); // lo > hi
        assert!(ModeMedianRatioBounds::new(0.0, 1.5).is_none()); // lo not > 0
        assert!(ModeMedianRatioBounds::new(f64::NAN, 1.5).is_none()); // non-finite
        assert_eq!(
            ModeMedianRatioBounds::DEFAULT,
            ModeMedianRatioBounds { lo: 0.5, hi: 1.5 }
        );
    }
}
