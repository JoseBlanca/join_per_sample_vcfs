# Code review — hidden-paralog filter Q2 (`SingleCopyCoverageModel::fit`)

## 1. Scope
- **Reviewed:** the new `src/paralog/coverage_model.rs` (the per-sample
  single-copy coverage-model fit) + a re-export line in `src/paralog/mod.rs`.
- **Against:** branch `tomato2-paralog-filter`, uncommitted diff.

## 2. Method
Two parallel review sub-agents per `ai/skills/rust-code-review`: (a) reliability
+ refactor_safety with a numerical-correctness emphasis; (b) naming + errors +
idiomatic + defaults + smells + module_structure. Verification (dev container):
`cargo test --lib paralog` 24 passed (was 13; +11 challenge/boundary tests);
`cargo clippy --lib --tests -D warnings` clean; `cargo fmt --check` clean.

## 3. Verdict
**Sound after fixes.** The reviewers confirmed the core numerics (parabolic
vertex + concavity guard, bin-centre convention, `div_ceil` lower-median
crossing, overflow-at-lower-edge rank safety, gap-fill indices). Three
result-affecting defects and one API-contract defect were found and fixed; no
issue survived.

## 4. Findings (all fixed)

- **B/M — σ₀ silently floored to `1e-3` on band collapse** (`fit_sigma0`).
  When the single-copy band collapses to one depth bin (MAD = 0), σ₀ was
  clamped to a hardcoded `1e-3` — a Normal ~1000× tighter than the fitted
  ~0.28, making a 1%-off single-copy window score as astronomically improbable
  (silent, sample-wide mis-calibration). → **Fixed:** floor at the histogram's
  own resolution, `0.5·depth_bin_width/single_copy_scale` (≈ 0.04 at production
  binning) — the honest "unknowable below one bin" σ₀, still strictly positive.
- **M — `ZeroScale` error variant was unreachable/dead.** The scale is always
  `≥ 1.0·width > 0` (the negative parabolic offset is gated on `peak > 0`), so
  the advertised "bottom-bin sample rejected" contract was false — such a
  sample was silently accepted with `scale = 0.5·width`, inflating every
  window's copy number. → **Fixed:** replaced with a reachable
  `DepthModeAtBottomBin { scale, depth_bin_width }` gated on `peak == 0`.
- **M (design) — `mode_median_ratio_bounds: (f64, f64)` primitive obsession.**
  A bare tuple with an ordered-pair invariant enforced only far away in
  `validate_config`, diverging from Q1's own `GridSpec` precedent. → **Fixed:**
  introduced `ModeMedianRatioBounds { lo, hi }` with a checked `new()`, a
  `DEFAULT` const, and `contains()`.
- **Mi — σ₀ band applied to the corrected `rcn`, not `rel` as the prototype
  does.** → **Fixed:** band on the un-corrected `rel`, take the MAD on the
  corrected `rcn` (matches `build_gc_normalization.py`).
- **Mi — `smooth_median` upper-median on even-length neighbourhoods** (differs
  from the prototype's `np.median` at the clamped curve ends). → **Fixed:**
  average the two central values for even lengths.
- **Mi — mode tie-break (`max_by_key` last-max) undocumented/unpinned.** →
  **Fixed:** documented the higher-bin convention + a pinning test.
- **Efficiency — `depth_marginal` swept twice** (fit + `median_depth`). →
  **Fixed:** `median_depth` takes the precomputed marginal.
- **Nits — `SIGMA0_FLOOR` promoted to a resolution-derived value; defensive
  `curve[gc_bin.min(len-1)]` simplified to a `debug_assert_eq!` + direct index;
  `anchors.last().unwrap()` given a `PANIC-FREE` comment.** All fixed.
- **Skipped w/ rationale:** the `(f64, u64)` value/weight pairs → newtype (Low
  confidence; `weighted_median` is a genuinely generic utility with the meaning
  documented on the fn).

## 5. Tests added (11)
Bottom-bin rejection, σ₀ resolution floor on collapse, `NoGcCurveAnchor` at the
fit level, all `InvalidConfig` sub-checks, asymmetric parabolic refinement,
mode tie-break, `gc_multiplier` clamp, even-length `smooth_median` averaging,
`weighted_median` even-total lower-middle, single-anchor `gap_fill`, and
`ModeMedianRatioBounds::new` validation.

## 6. What's good
Names say what values are (`HistogramLayout`, `depth_marginal`, `mode_depth`);
`#[non_exhaustive]` error enum matches the crate convention; defaults live in
documented `pub const`s; the module is a clean pure-statistics peer depending on
`sample_summary` only (arch Premise 0).

Audit trail: `tmp/review_2026-07-01_paralog-q2/{reliability,naming_etc}.md`.
