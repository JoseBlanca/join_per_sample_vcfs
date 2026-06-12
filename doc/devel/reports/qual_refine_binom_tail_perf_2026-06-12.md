# QUAL-refine binomial tail: exact incomplete-beta above the cap

**Date:** 2026-06-12
**Area:** `src/vcf/qual_refine.rs` (variant-QUAL refinement)
**Builds on:** `1e5212e` (2026-06-11, *"bound QUAL-refine binomial tail to
avoid O(depth) hang"*), which first bounded the cost with a normal
approximation above an `n` cap.

---

## 1. Background

`refine_qual` runs on **every emitted record** and calls `tail_phred` three
times (allele-balance + two strand/position-bias penalties), each computing a
two-sided binomial-tail p-value `binom_two_sided_p(k, n, p)`. The original code
enumerated **every** outcome `0..=n`, where `n` is the **cohort-summed read
depth** — O(cohort depth) per record, and a multi-minute hang on the
`num_obs = u32::MAX` encoder-overflow tests (they need `> i32::MAX` to trip the
encoder's depth guard, so they couldn't avoid the loop).

`1e5212e` bounded this: keep the exact discrete sum for `n ≤ 2000` (byte-
identical QUAL for all validated single-sample depths — the depth sweep ran at
≤301x) and switch to a continuity-corrected **normal approximation** above the
cap.

## 2. Why the normal approximation wasn't enough

A scipy-calibrated harness (`tmp/qual_binom_calib.py`) graded both methods
against `scipy.stats.binomtest` (which the exact discrete sum matches to 0.0000
phred). The symmetric normal approximation is accurate near `p = 0.5` but
**30–97 phred wrong** at the extreme `p` values `refine_qual` actually produces
just above the cap:

- balance `p` is clamped to `[1e-6, 0.999]`, bias `p` to `[0.01, 0.99]`;
- a skewed binomial's two-sided Sterne tail is **not symmetric**, and the phred
  (log-tail) metric magnifies the deep-tail divergence;
- e.g. `n` just over 2000 with `p = 0.01` (a low-frequency variant in a large
  cohort, or strong REF strand bias) lands ~35 phred off in the decision band.

The variance threshold for ≤0.5-phred agreement is `np(1−p) ≈ 2500` —
effectively unreachable. The shipped `1e5212e` test only checked `p ∈ {0.3,
0.5}`, where the normal approximation happens to be fine.

## 3. The change: exact tail via the regularized incomplete beta

Kept `1e5212e`'s structure (exact discrete sum for `n ≤ EXACT_TAIL_MAX_N`,
preserving byte-identical QUAL) but **replaced the normal approximation above
the cap with the exact binomial CDF in closed form**:

- `P(X ≤ k) = I_{1−p}(n−k, k+1)` via a self-contained continued-fraction
  `reg_incomplete_beta` (Numerical Recipes `betai`/`betacf`, reusing the
  existing Lanczos `ln_gamma`; no new dependency — matches the file's
  "self-contained math, no external deps" ethos);
- the two-sided Sterne set ("outcomes no more likely than `k`") is the near tail
  through `k` plus a far tail beyond the opposite-flank index of equal pmf; the
  binomial is unimodal, so that far cutoff is one **binary search** over
  `ln_binom_pmf`.

Result: above the cap the tail is now **exact at every `p`** (not approximate)
and still **O(log n)** — independent of read depth. `binom_two_sided_p_normal`
and `erfc` are removed.

The cap is now a **byte-identity boundary, not an accuracy compromise**: below
it the exact discrete sum is bit-for-bit-identical to the pre-`1e5212e` QUAL;
above it the beta method is exact where `1e5212e` was approximate.

## 4. Validation

The discrete sum is the reference (== `scipy.stats.binomtest` to 0.0000 phred).
Tests added:

- **`beta_tail_matches_enumeration_across_grid`** — `binom_two_sided_p_beta`
  vs the enumeration across `p ∈ {1e-6 … 0.999}`, `n` up to 4000 (incl. `>`
  cap), `k` over both tails. Decision band (penalty ≤ 100 phred): agree to
  `< 1e-3` phred. Deeper: both saturate past 100 (the only divergence is
  ~4.4 phred at penalty ≈186, where QUAL clamps to 0 regardless).
- **`production_dispatch_matches_path_by_cap`** — the entry point is bit-for-bit
  the discrete sum at/below the cap and the beta method above it (locks in the
  byte-identity boundary).
- **`reg_incomplete_beta_matches_scipy_golden`** / `…_boundaries_and_symmetry`
  / **`binom_cdf_matches_direct_summation`** — anchor the `betainc` primitive
  (pinned scipy values, `I₀=0`/`I₁=1`/`Iₓ(a,b)=1−I₁₋ₓ(b,a)`, CDF vs direct sum).
- **`beta_invariants_at_huge_depth`** — at `n = 2,000,000`: p = 1 at the mode,
  monotone away, always finite.
- **`large_n_is_bounded_and_sane`** (kept from `1e5212e`) — `u32::MAX` depth
  returns instantly, pinned at the clamp floor.
- **`proptest_beta_matches_enumeration`** — random `(n ≤ 3000, p, k)` fuzz.

The differential test caught a real edge-case bug during development (far-tail
range started at `mode+1` and a `mode ≥ n` guard wrongly zeroed the opposite
tail at `n=1, p=0.5, k=0`); fixed by searching `[mode, n]` and special-casing
`k == mode → p = 1`.

## 5. Result

- Two overflow tests run in 0.00 s (were >11 min before `1e5212e`).
- Full `cargo test --lib`: **all pass, 0 skipped**.
- `cargo clippy --all-targets -- -D warnings` / `doc` / `fmt` clean.
- Behavioural anchors (`balance_penalty_grows_with_depth_at_fixed_low_vaf`,
  `balanced_call_pays_no_penalty`) and `1e5212e`'s lint tidies unchanged.

## 6. Notes

- **Numeric output:** QUAL is byte-identical to the pre-fix code for `n ≤ 2000`
  (all validated depths). Above the cap it is now *exact* where `1e5212e` was
  approximate, so any change vs `1e5212e` there is a correction toward the true
  binomial tail.
- **Separate discussion (unchanged):** `refine_qual` uses *cohort-summed* depth
  for an allele-balance test; whether that should be per-sample is a modelling
  question, not addressed here.
- **Calibration harness:** `tmp/qual_binom_calib.py` (scipy) is the throwaway
  audit tool; the durable verification is the test suite above.
</content>
