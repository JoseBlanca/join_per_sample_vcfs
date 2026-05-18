# Performance Review: posterior_engine

**Date:** 2026-05-18
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** `src/var_calling/posterior_engine.rs` (4000 lines) + the three
submodules `posterior_engine/{backends.rs, interp.rs, shape.rs}` —
Stage 6 of the multi-stage SNP caller. Branch `review/vcf-writer` at
`9cf0010`.
**Verdict:** **Run experiments** — multiple plausible candidates, one
algorithmic restructuring with measurable upside, and a clear
measurement plan. No new sampling profile was collected this round;
findings build on the May 17 / May 18 profile + bench evidence quoted
verbatim from the in-tree implementation reports.
**Hot-path evidence:** existing benchmark + profile artefacts:

```
| Backend                              | µs/record   |
| ExactMath                            | 24.41       |
| InterpUnivariateMath (scalar interp) | 22.33 (-8.5%)|
| InterpUnivariateSimdMath (wide native)| 18.15 (-25.6%)|
```

(Fixture: 64 samples × 10 000 records biallelic SNP, `c_s=3%`,
`q_b=[0.6,0.3,0.1]`; via [examples/profile_posterior_engine.rs](../../../examples/profile_posterior_engine.rs).
Source: [posterior_engine_simd_analysis_2026-05-18.md](../implementations/posterior_engine_simd_analysis_2026-05-18.md).)

The pre-SIMD May-17 flamegraph attributed ~40 % self-time to `ln`/`exp`
and ~30 % to `e_step`. After the wide-native swap the relative
distribution has not been re-flamegraphed, so several findings below
mark themselves *Speculative* on confidence grounds.

---

## 1. Scope and constraints

**What was reviewed:** the posterior_engine module on the current
branch (commit `9cf0010`). Four in-scope files:

- [src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs) (4000 lines; EM driver, both `e_step` flavours, mixture pre-pass, log-sum-exp helpers, summarisation)
- [src/var_calling/posterior_engine/backends.rs](../../../src/var_calling/posterior_engine/backends.rs) (217 lines; `MathBackend` trait + three concrete backends)
- [src/var_calling/posterior_engine/interp.rs](../../../src/var_calling/posterior_engine/interp.rs) (237 lines; IEEE-decomposition scalar ln/exp tables)
- [src/var_calling/posterior_engine/shape.rs](../../../src/var_calling/posterior_engine/shape.rs) (281 lines; per-`(ploidy, n_alleles)` shape cache)

**Performance intent and targets:** per-record EM (3–5 typical
iterations + a final E-step) producing posteriors, allele frequencies,
GQ/QUAL summaries. Production = x86_64 server (max priority); dev =
Apple Silicon M5 (aarch64/NEON, testing only). Cohorts: 10⁵–10⁶ records
× tens to thousands of samples. Accuracy budget: 1e-4 on posterior
parity; current backend has 50× margin (max ~2e-6).

**Hot-path evidence available:** existing bench + profile artefacts
quoted above. No fresh samply or DHAT was run for this review; the
per-category sub-agents and the cross-codebase research worked from the
quoted evidence.

**Deliberately out of scope:**
- Rayon-over-records parallelism (deferred follow-up, separate plan).
- Cohort CLI wiring.
- Legacy gVCF-merger code under `src/gvcf_parser.rs` etc.
- BED-region filter, Tabix `.tbi` index.

**Categories dispatched** (per the perf-review skill):
- `methodology` — always.
- `allocations` — every per-record `Vec::new` / `vec!` / `collect` is a candidate.
- `data_layout` — `log_likelihoods` is sample-major; SIMD path gathers it strided.
- `hot_loops` — the inner-loop content is the core of the review.

Two categories skipped: `concurrency` (single-threaded today;
rayon-over-records is its own follow-up plan) and `io_and_syscalls`
(pure compute, no I/O after the upstream merger).

Audit trail at [tmp/perf_review_2026-05-18_posterior_engine/](../../../tmp/perf_review_2026-05-18_posterior_engine/) (four
per-category finding files, retained as evidence).

**External research grounded against local source:** the GATK
(`gatk/`), freebayes (`freebayes/`), and bcftools (`bcftools/`) source
trees are vendored at the project root. The following were verified
in-tree (not from web blog posts):

- GATK's `JacobianLogTable` exists at
  [gatk/src/main/java/org/broadinstitute/hellbender/utils/MathUtils.java:406-423](../../../gatk/src/main/java/org/broadinstitute/hellbender/utils/MathUtils.java#L406-L423)
  — 80 001 entries over `[0, 8]` step `1e-4` tabulating
  `log10(1 + 10^-x)`. Consumed by `approximateLog10SumLog10` and used as
  the default path (the exact form is only for unit tests).
- bcftools' `pl2p` table at
  [bcftools/mcall.c:59-64](../../../bcftools/mcall.c#L59-L64) — Phred→prob,
  not log; no log-sum-exp table. EM is in
  [bcftools/em.c](../../../bcftools/em.c) (`ITER_MAX = 50`, `ITER_TRY = 10`,
  `EPS = 1e-5`); no warm-start across sites.
- freebayes caches `lgamma(n + 1)` and per-genotype `permutationsln` —
  the same shape our [shape.rs](../../../src/var_calling/posterior_engine/shape.rs) `log_multinomial_coeffs` already
  covers.
- **No SIMD in any of the three.** `wide::f64x4::{ln,exp}` puts our
  default backend ahead of all three reference implementations.
- None implements SQUAREM, Anderson acceleration, or warm-start.

These ground-truths reshape the finding priority: log-sum-exp tabulation
(GATK precedent, in-tree) is the cheapest big-impact change; SQUAREM /
warm-start (paper-only, no caller precedent) drop to Speculative.

## 2. Verdict

**Run experiments.** Three Hot-path findings have high confidence, name
the exact site, and propose a falsifiable bench. Two of them
(`log_indep` cross-batch hoist, homogeneous-fixation short-circuit) are
algorithmic — bit-identical when triggered. One (softplus_neg table for
`log_sum_exp_2_x4`) introduces a small approximation that the existing
parity harness already gates.

Eight Likely findings — most are mechanical wins (per-engine scratch
buffer, `#[inline(always)]`, `#[cold]`) blocked only on a bench
comparison. The single biggest open question is whether the
`log_likelihoods` transpose-on-entry (data_layout L8/L9) actually pays
at `n_genotypes = 3`; the recommendation is to run it but with a low
merge threshold (≥ 5 % wall time).

Several Speculative items (SQUAREM, warm-start across records, adaptive
genotype-set pruning à la Octopus, AVX-512 `f64x8` upgrade, PGO) are
listed for the team to triage after the Hot-path candidates have been
measured and applied. They don't get a measurement plan against the
current bench because the dominant inner-loop cost shifts as the
Hot-path findings land.

## 3. Measurement plan

In the order they unblock other findings:

1. **Fresh samply profile of the contam-on workload, post-SIMD.**
   `samply record -- target-container/release/examples/profile_posterior_engine`
   with `POSTERIOR_BACKEND=simd`. Today's profile evidence (May 17) is
   *pre*-SIMD swap; relative self-times after the swap have not been
   re-measured. The May 18 SIMD analysis used a different driver
   (`profile_posterior_engine` directly) but didn't capture a
   flamegraph. Without this, several Likely findings cannot be
   promoted to Hot-path. Threshold: a flamegraph where each of the top
   5 self-time frames is named and >= 5 % of inclusive time.
   ([feedback_samply_no_permission_prompt.md] memory: run on host,
   blanket permission granted.)

2. **DHAT pass on the same fixture** with `--bench` and a small (~100
   record) drain. Today's allocations findings are pattern-match only;
   DHAT will name the per-record allocation sites and rank them so we
   know whether the per-engine `RecordScratch` lift (Likely L5–L7) is
   worth the API churn. Threshold to merge those findings: ≥ 5 fewer
   mallocs per record on the DHAT-after pass.

3. **In-bench self-validation assertion** in `bench_posterior_drain`:
   `debug_assert!((2..=10).contains(&r.diagnostics.iterations))` —
   guards every following measurement against a silent regression to
   `DidNotConverge` or to the trivial-record path. One-line change;
   no perf cost in release builds. See methodology finding M1.

4. **Triple back-to-back bench runs** in the same shell session for
   every code-level finding (per the existing `feedback_bench_cpu_governor`
   memory: the i7-1260P's powersave governor swings ~2× across thermal
   states). Cross-session comparisons are unreliable.

5. **Parity-budget regression** for any approximation finding
   (Hot-path H3 softplus_neg, Likely L1 branchless multiply, etc.):
   the existing `tests_math_backend_accuracy` 12-fixture harness must
   stay within 1e-4 max posterior error. No new test required.

6. **`cargo asm`** for the two Likely codegen findings (`#[inline(always)]`
   on the `_x4` helpers, `#[cold]` on `classify_nonfinite`). Each only
   merges if the asm change is visible and the bench moves ≥ 0.5 %.

## 4. Build / toolchain configuration

Baseline is solid; no code-level finding is blocked on build config.
Per the methodology sub-agent:

- `[profile.release]`: `lto = "fat"`, `codegen-units = 1`,
  `panic = "abort"`, `debug = "line-tables-only"`, `opt-level = 3`
  (default).
- `[profile.bench]` inherits `release` with `debug = true` for
  flamegraph symbol resolution.
- Toolchain pinned at `rust-toolchain.toml` channel `1.95` — cross-commit
  baselines stay comparable.
- `.cargo/config.toml`: `target-cpu = x86-64-v3` (Linux) / `apple-m1`
  (macOS aarch64). Explicit microarchitectures, not `native` —
  binaries portable within their deployment class.
- `mimalloc` wired behind an opt-in `alloc-mimalloc` feature; activated
  in `benches/var_calling_perf.rs`.

Two Speculative build-config findings worth scheduling but not blocking
on:

- **PGO** — release profile already has all the prerequisites (LTO fat,
  cg-units 1). Train against the contam-on
  [examples/profile_posterior_engine.rs](../../../examples/profile_posterior_engine.rs)
  drain; rebench. Threshold: ≥ 5 % single-thread win at no parity cost.
- **`target-cpu = x86-64-v4`** — would let the autovectoriser emit
  AVX-512 in scalar reduction loops (`m_step_p_hat`, `summarise_posteriors`).
  Verify production CPU model first; if uniformly Skylake-X+, raise the
  floor (one-line change in `.cargo/config.toml`). Won't affect the
  explicit `f64x4` SIMD path — `wide::f64x4` doesn't widen to AVX-512;
  the parallel candidate there is an explicit `f64x8` backend (see S5
  below).

`group.sample_size(10)` on `bench_posterior_engine` is below criterion's
default of 100. Intentional given the per-sample budget (~1.5 s); cross-commit
comparisons off this group need the revert-experiment discipline of the
methodology checklist.

## 5. Code-level findings

Severity bands: **Hot-path** (high-confidence + named site + matching
profile) → **Likely** (mechanism is real, measurement plan defined) →
**Speculative** (plausible, lower-confidence) → **Note** (recorded,
no action).

### Hot-path

#### H1 — [posterior_engine.rs:1691-1699](../../../src/var_calling/posterior_engine.rs#L1691-L1699) — Hidden SIMD slowdown: `log_sum_exp_2_x4` is invoked on the homozygous branch even when `log_f_v` is all `-∞` (the **default** config)

- **Confidence:** High.
- **Mechanism:** `PosteriorEngineConfig::with_project_defaults()` sets
  `fixation_index_default = 0.0` and `fixation_index_overrides = None`,
  so `log_f_per_sample = safe_ln(0.0) = -∞` for every sample on the
  default path. The scalar `log_sum_exp_2` ([posterior_engine.rs:1722-1731](../../../src/var_calling/posterior_engine.rs#L1722-L1731))
  short-circuits on `a == NEG_INFINITY`. The SIMD sibling
  `log_sum_exp_2_x4` ([posterior_engine.rs:2053-2059](../../../src/var_calling/posterior_engine.rs#L2053-L2059))
  does not. Every iteration × every homozygous genotype × every 4-sample
  batch dispatches three transcendentals whose closed-form answer on the
  default path is just `term1 = log_one_minus_f_v + log_indep_v`.
- **Hot-path evidence:** the May 17 flamegraph put `ln`/`exp` at ~40 % self-time. Even after the SIMD swap, every wasted SIMD `exp`/`ln` here is on the dominant path.
- **Measurement plan:** add the short-circuit; rebench the contam-on
  fixture triple-back-to-back; expect ≥ 5 % wall-clock reduction.
  Parity test: `tests_math_backend_accuracy` must pass unchanged (this
  is *more* accurate than the current code, not less).
- **Complexity cost:** add a record-static `homogeneous_f_zero: bool`
  to `EmContext` and hoist the branch above the genotype loop. ~10
  lines. No `unsafe`, no new dep. **Subsumed by H4 below in a
  stronger form** — implement H4 instead of H1 as a standalone if both
  are accepted.

#### H2 — [posterior_engine.rs:1679-1689](../../../src/var_calling/posterior_engine.rs#L1679-L1689) — `log_indep` recomputed per `(batch_idx, g_idx)` in `e_step_simd`, identical across all batches in one EM iteration

- **Confidence:** High. The cross-batch invariance was already noted in
  the May 18 SIMD analysis report (§4.4 "Lane = sample → most
  per-genotype work is scalar splat"), but the fix has not landed.
- **Mechanism:** `log_indep = log_multinomial_coeffs[g_idx] + Σ k ·
  log_p_effective[a]`. None of those terms touches `s0` or the lane
  index; for the 64-sample bench, that's 16 redundant evaluations per
  (iter, genotype). Hoist into a `Vec<f64>` of length `n_genotypes`
  computed once per E-step call (between the `log_p_effective` build
  and the batch loop), then splat from there.
- **Measurement plan:** add a per-engine scratch slot for
  `log_indep_per_g`; rebench. Threshold: ≥ 2 % wall-clock or a clear
  flamegraph delta on the genotype-inner block.
- **Complexity cost:** one extra `Vec<f64>` of length `n_genotypes`
  in `EmScratch` (no allocation in steady state once L5 lands). One
  extra pre-pass over genotypes per E-step call. Bit-identical.

#### H3 — [posterior_engine.rs:2053-2059](../../../src/var_calling/posterior_engine.rs#L2053-L2059) — `log_sum_exp_2_x4` does 3 transcendentals per call where `m + softplus_neg(|a-b|)` needs 1

- **Confidence:** High (mechanism); Medium that the win survives H4.
- **Mechanism:** identity `log_sum_exp(a, b) = max(a, b) +
  ln(1 + exp(-|a - b|))`. Today: `exp(a-m) + exp(b-m)` then `ln(sum)` —
  3 transcendentals. With a tabulated `softplus_neg(x) = ln(1 + exp(-x))`
  over `x ∈ [0, 30]` (linear interp, 1024 bins → 8 kB, ≤ 5e-6 absolute
  error; above 30 return 0), this collapses to one table read + interp.
- **Local-source precedent (verified):**
  [gatk/src/main/java/org/broadinstitute/hellbender/utils/MathUtils.java:406-423](../../../gatk/src/main/java/org/broadinstitute/hellbender/utils/MathUtils.java#L406-L423)
  tabulates `log10(1 + 10^-x)` over `[0, 8]` step `1e-4` as the default
  path. Same idea in natural log space.
- **Measurement plan:**
  1. Build the table in [interp.rs](../../../src/var_calling/posterior_engine/interp.rs);
     add a `MathBackend::softplus_neg_x4` method (or inline it).
  2. Sweep ULP: `a, b ∈ [-30, 30]`, assert absolute error ≤ 5e-6 vs
     the `wide`-native baseline; assert end-to-end posterior parity
     ≤ 1e-4 via `tests_math_backend_accuracy`.
  3. Triple-back-to-back contam-on bench; threshold ≥ 5 % wall-clock.
- **Complexity cost:** new 8 kB `LazyLock` static + one `MathBackend`
  method (defaulted to scalar fallback for `ExactMath`). The
  approximation tax is real (~5e-6 vs the current ~2e-6 from
  `wide`-native ln/exp), so document the budget headroom explicitly.
- **Order with H4:** apply H4 first. If H4 alone removes ≥ 80 % of
  the `log_sum_exp_2_x4` calls (the default-config `f_s = 0` case),
  this finding drops to Speculative — the remaining calls are
  override-set workloads, which are rare.

#### H4 — Posterior engine — homogeneous-fixation hoist (covers H1 + a larger structural saving)

- **Site:** [posterior_engine.rs:1672-1689 + 1691-1699](../../../src/var_calling/posterior_engine.rs#L1672-L1699)
  (the `for g_idx in 0..n_genotypes { … }` inner loop body of
  `e_step_simd`) plus the scalar `e_step` body at
  [:1568-1578](../../../src/var_calling/posterior_engine.rs#L1568-L1578).
- **Confidence:** High (mechanism). Subsumes H1.
- **Mechanism:** When `fixation_index_overrides = None` (default), every
  sample shares the same `f_s = fixation_index_default`. Then:
  - `log_f` and `log_one_minus_f` are scalar constants (not n_samples
    values).
  - `log_prior[g]` depends only on `(log_f, log_one_minus_f,
    log_indep[g], log_p_effective[homozygous_allele_for[g]])` — *all
    sample-invariant within an EM iteration*.
  - Therefore the entire per-batch `log_prior_v` block reduces to one
    `f64x4::splat(log_prior[g])` after a single scalar precompute per
    `(record, iter, g)`.

  This removes `n_homozygous_g × n_batches × n_iters` `log_sum_exp_2_x4`
  calls per record (for 64 samples × 2 homozygous gts × 3–5 iters ≈ 96–160
  calls per record → 1–2 scalar calls per iter).
- **Hot-path evidence:** the May 17 flamegraph quoted `e_step` at
  ~30 % self-time and `ln`/`exp` at ~40 %. The transcendentals removed
  here are the same ones that flamegraph attributed cost to. Post-SIMD
  numbers haven't been re-flamegraphed (see measurement plan §1).
- **Measurement plan:** add a `homogeneous_fixation: bool` to
  `EmContext` (computed once at `run_em_for_record` entry); when true,
  compute a scratch `log_prior_per_g: Vec<f64>` of length `n_genotypes`
  once per EM iteration; the batch loop just splats from it. Run the
  contam-on bench triple-back-to-back. Threshold: ≥ 10 % wall-clock
  reduction on the default-config path. Parity: bit-identical when
  triggered (purely algebraic refactor).
- **Complexity cost:** ~30 lines. One new `bool` on `EmContext`, one
  scratch slot, one branch. Bit-identical when triggered; the
  heterogeneous path remains as-is.
- **Rationale for ranking above H1:** H1 is a special case of H4 (just
  the `f_s = 0` short-circuit); H4 covers the broader case (any constant
  `f_s` across samples, including non-zero overrides applied uniformly,
  e.g. when the CLI ships a single global `--inbreeding` flag).

### Likely

#### L1 — [posterior_engine.rs:1869-1873](../../../src/var_calling/posterior_engine.rs#L1869-L1873) — `if k != 0` branch in `accumulate_expected_counts` general path blocks autovectorisation on tri-/tetra-allelic records

- **Confidence:** Medium. Biallelic-diploid hot path is specialised
  elsewhere; this affects records with `n_alleles ≥ 3`.
- **Mechanism:** `weight * (k as f64)` is `0.0` when `k == 0`, and
  `acc += 0.0` is a no-op under IEEE 754. The branch prevents the inner
  loop from being a straight-line vector store.
- **Measurement plan:** `cargo asm` on the function (confirm FMA
  emitted, no conditional jump); synthetic triallelic-diploid bench;
  threshold 2 %+ reduction on that subset.
- **Complexity cost:** one-line change.

#### L2 — [backends.rs:127-137](../../../src/var_calling/posterior_engine/backends.rs#L127-L137) — `InterpUnivariateMath::HAS_LANE_4 = false` leaves a one-line speedup on the table

- **Confidence:** Medium.
- **Mechanism:** the scalar interp backend gets the `f64::ln`/`f64::exp`-via-table
  fast path but skips the SIMD code path entirely. The default
  `MathBackend::ln_x4`/`exp_x4` fallback is "call scalar four times" —
  which would still route through the (already-fast) `interp::ln_approx`
  / `exp_approx`. Cost: one `const HAS_LANE_4: bool = true;`.
  This isn't a "real" SIMD path — it's letting `e_step_simd` run with
  the four scalar `ln_approx` calls per lane — but it removes the
  scalar-vs-SIMD branch divergence and pairs nicely with H4.
- **Measurement plan:** flip the flag; rerun all `posterior_engine`
  unit tests + the contam-on bench. Threshold: ≥ 1 % wall-clock on
  `interp_univariate` (non-SIMD-default) sub-bench, or no regression.
- **Complexity cost:** one line. No new tests required (the existing
  `tests_math_backend_accuracy` already covers it).

#### L3 — [posterior_engine.rs:2037, 2052, 2067](../../../src/var_calling/posterior_engine.rs#L2037) — `#[inline]` on `log_sum_exp_*_x4` should be `#[inline(always)]`

- **Confidence:** Medium.
- **Mechanism:** the helpers take a `wide::f64x4` and call through a
  generic `MathBackend` method. With `lto = "fat"` and `codegen-units
  = 1` the optimiser usually inlines, but `#[inline]` is a hint.
  `#[inline(always)]` guarantees the cross-impl path and avoids spilling
  the SIMD vector to memory if the optimiser changes its mind.
- **Measurement plan:** `cargo asm` on `e_step_simd` before/after;
  confirm no `call …log_sum_exp_2_x4…` remains. Rebench; threshold ≥ 0.5 % wall-clock.
- **Complexity cost:** attribute swap, no body change.

#### L4 — [posterior_engine.rs:826-831, 1007-1014, 1069-1074, 1594-1599, 1604-1611, 1720-1726, 1735-1741, 1789-1794](../../../src/var_calling/posterior_engine.rs#L826) — `#[cold]` missing on `classify_nonfinite` / non-finite-error build sites

- **Confidence:** Medium.
- **Mechanism:** Err arms include `RecordLocus` construction +
  `classify_nonfinite`; today they sit inline in the hot path, competing
  for icache. `#[cold]` on `classify_nonfinite` and a single shared
  cold helper (e.g. `non_finite_posterior_err(...)`) moves the error
  build out of the hot block.
- **Measurement plan:** `cargo asm` before/after; bench threshold
  0.5 %+.
- **Complexity cost:** one attribute + (optionally) one small helper
  function.

#### L5 — [posterior_engine.rs:1318, 1281-1295, 1319-1321](../../../src/var_calling/posterior_engine.rs#L1318) — `EmScratch::new` + 8 satellite `Vec`s allocated per record; promote to a per-engine `RecordScratch`

- **Confidence:** High (pattern); Medium (measurable on top of 18.15 µs).
- **Mechanism:** `EmScratch::new` builds 5 `Vec`s; the surrounding
  scope builds 8 more (`log_f_per_sample`, `log_one_minus_f_per_sample`,
  `allele_classes`, `pseudocounts`, `compound_mask`, `p_hat`,
  `f_hat_compound`, `posteriors`). All sized by `n_alleles` / `n_samples`
  / `n_genotypes`, all rebuilt at the same shape across records. Move
  to a per-`PosteriorEngine` `RecordScratch` with `resize_to(...)` at
  function entry; the EM loop already overwrites every cell so
  `Vec::resize` keeps capacity in the steady state.
- **Measurement plan:** add `posterior_engine_em_scratch_reuse`
  Criterion sub-bench rolling all 13 allocations into one fix; back-to-back
  before/after. DHAT confirms ≥ 10 fewer mallocs/record. Threshold
  3 %+ wall-clock improvement to merge.
- **Complexity cost:** medium. Adds a `&mut RecordScratch` parameter
  threaded through `run_em_for_record`. `p_hat` / `f_hat_compound` /
  `posteriors` are currently `move`d into `PosteriorRecord`; lifting
  them to scratch swaps "1 malloc + zero-init" for "1 memcpy out". The
  bench decides whether memcpy ≤ malloc+init at the dominant shape;
  if not, leave those three on the per-record path and lift only the
  internal scratches.

#### L6 — [posterior_engine.rs:1803-1819, 1878-1894](../../../src/var_calling/posterior_engine.rs#L1803-L1894) — `m_step_p_hat` and `m_step_f_hat_compound` `.collect()` a fresh `Vec<f64>` per EM iteration

- **Confidence:** High (pattern). With 3–5 iters × 10⁵–10⁶ records that's
  6 × 10⁵ to 1 × 10⁷ `Vec` allocations per cohort.
- **Mechanism:** the M-steps return owned `Vec<f64>` and the caller
  drops/replaces. Rewrite to write into `&mut [f64]` from the scratch.
  The `max_abs_diff` check needs both buffers, so double-buffer with two
  scratch slots and rotate.
- **Measurement plan:** rolled into L5's bench.
- **Complexity cost:** double-buffer rotation needs a clear comment;
  the M-step signatures change.

#### L7 — [posterior_engine.rs:761-768, 898-919](../../../src/var_calling/posterior_engine.rs#L919) — Mixture pre-pass allocates `mixture_log_likelihoods` and (SIMD path) `fallback_log_likelihoods.to_vec()` per contam-on record

- **Confidence:** High.
- **Mechanism:** `run_em_for_record` already owns the input
  `log_likelihoods: Vec<f64>`. Pass it as `&mut Vec<f64>` to the mixture
  function, which overwrites in place. Removes the malloc + memcpy.
  Lifting `n_obs` / `mean_err` / `c_s_all` / `mean_err_all` into the same
  per-engine scratch removes 3–4 more mallocs per contam-on record.
- **Measurement plan:** contam-on sub-bench only. Threshold 3 %+.
- **Complexity cost:** signature change; the in-source comment at
  [posterior_engine.rs:763-766](../../../src/var_calling/posterior_engine.rs#L763)
  documents that a similar hoist (`n_obs` / `mean_err`) was already
  worth ~20 % of contam-on self-time per a prior samply run — so the
  rule of thumb has been demonstrated once on this exact code.

#### L8 — [posterior_engine.rs:1703-1708](../../../src/var_calling/posterior_engine.rs#L1703-L1708) — Strided gather of `log_likelihoods` in `e_step_simd` batch loop

- **Confidence:** Medium. Flagged but unmeasured by the May 18 SIMD
  analysis report.
- **Mechanism:** 4 scalar reads with stride `n_genotypes * 8` bytes
  (24 B diploid biallelic; 80 B compound 10-genotype). A one-shot
  transpose to batch-of-4-samples-major at record entry (length
  `(n_samples/4) * n_genotypes * 4`) makes each lane load a contiguous
  32 B `vmovupd`. Cost: O(n_samples * n_genotypes) reads/writes per
  record, paid once vs. the gather firing `n_genotypes` times per EM
  iteration plus a final pass.
- **Measurement plan:** behind a const/feature flag, materialise the
  transposed buffer; bench. Threshold ≥ 5 % per-record wall-time
  reduction. If < 2 %, the transpose memory traffic doesn't pay back
  at `n_g = 3` — revert.
- **Complexity cost:** medium — the scalar `e_step` tail and the
  M-step / summariser still want sample-major; cleanest packaging
  builds the transpose only when `M::HAS_LANE_4` and stores it in the
  per-engine scratch (depends on L5 plumbing).

#### L9 — [posterior_engine.rs:1743](../../../src/var_calling/posterior_engine.rs#L1743) — Symmetric scatter of `posteriors` write-back in `e_step_simd`

- **Confidence:** Medium.
- **Mechanism:** mirror of L8 on the write side. Each lane writes to a
  different cache line. Same fix: write batch-major into a scratch,
  transpose to sample-major once after `run_em_loop` returns (the public
  `PosteriorRecord.posteriors` is a sample-major contract — keep the
  output shape).
- **Measurement plan:** combine with L8 in the same fix; merge only
  if combined wall-time ≥ 5 % reduction. Risk: more API surface change
  than L8 alone (3 readers of `posteriors`: `m_step_p_hat`,
  `summarise_posteriors`, the public record field).

#### L10 — [benches/var_calling_perf.rs:686](../../../benches/var_calling_perf.rs#L686) — Bench has no self-validation assertion on `diagnostics.iterations`

- **Confidence:** Medium.
- **Mechanism:** a refactor that silently flips the engine into the
  trivial-record path (`n_alleles < 2`) or the `DidNotConverge` short
  return would still produce 10 000 results with plausible-looking
  wall time. The fixture is engineered for 3–5 iterations; assert that.
- **Measurement plan:** add the `debug_assert!`; rerun once; passes.
- **Complexity cost:** one line.

#### L11 — [benches/var_calling_perf.rs:665](../../../benches/var_calling_perf.rs#L665) — `iter_batched` setup `merged.to_vec()` contaminates `--profile-time` flamegraphs

- **Confidence:** Medium. The project already works around this via
  [examples/profile_posterior_engine.rs](../../../examples/profile_posterior_engine.rs);
  the constraint isn't documented in the bench file itself.
- **Mechanism:** Criterion times only the body, but a sampling
  profiler records the setup too. `MergedRecord` `Clone` invokes
  several inner `Vec` clones; per-iteration setup performs
  O(records × inner-vecs) heap ops.
- **Measurement plan:** N/A — documentation finding. Future profile
  evidence sourced from the bench (vs. the example driver) must include
  a note that the bench's flamegraph is unreliable.
- **Complexity cost:** add a 6-line comment at the bench fn entry. Or
  refactor the bench to `iter_batched_ref` if the engine API can grow a
  `&[MergedRecord]` consumer (larger change; defer).

### Speculative

#### S1 — [posterior_engine.rs:1448, 2068-2084](../../../src/var_calling/posterior_engine.rs#L2068-L2084) — `log_sum_exp_slice_x4` two-pass → one-pass online softmax (Milakov–Gimelshein 2018)

- **Confidence:** Low for biallelic-diploid (`n_g = 3`, slice is short
  → two-pass is 6 reads vs one-pass 3 reads, win marginal). Higher for
  tetra-allelic (`n_g = 15`).
- **Measurement plan:** bench tetra-allelic only. File as a
  tetra-only specialisation if so.

#### S2 — [posterior_engine.rs:1136](../../../src/var_calling/posterior_engine.rs#L1136) — `log_post_unnorm_lane: Vec<wide::f64x4>` 32-B alignment

- Likely matters only after L8/L9 land and the gather/scatter is no
  longer the dominant memory-traffic cost. `perf stat -e
  L1-dcache-load-misses` decides.

#### S3 — [posterior_engine.rs:1731-1745](../../../src/var_calling/posterior_engine.rs#L1731-L1745) — `posteriors` scatter write per (g, lane)

- Subsumed by L9 (which does it on a batch boundary). Skip standalone.

#### S4 — SQUAREM acceleration of the EM map

- **Source:** Varadhan & Roland JSS 2021 ([arxiv 1810.11163](https://arxiv.org/abs/1810.11163)). No production caller uses it (verified absent in gatk, freebayes, bcftools, [octopus population_model.cpp](https://github.com/luntergroup/octopus/blob/master/src/core/models/genotype/population_model.cpp)). Would compress 3–5 EM steps into 1–2 outer SQUAREM steps (= 3–6 EM evals) — net wash on the typical path. Real win is on
  `DidNotConverge` cases (50-iter cap).
- **Complexity:** ~100 lines; needs a log-likelihood-sum pass per
  outer step (we don't track total LL today, just `max |Δp̂|`).
- **Defer** until L5–L9 + H1–H4 land and we know whether convergence is
  still hitting the cap.

#### S5 — AVX-512 `f64x8` backend with `multiversion` runtime dispatch

- **Source:** the `wide` 0.7 crate doesn't widen `f64x4` to AVX-512; the
  AVX-512 lane width would come from `wide::f64x8` (0.8+) or `pulp` /
  `std::simd`. Production = x86_64 server, possibly Sapphire Rapids /
  Zen 4-5.
- **Complexity:** new backend, multiversion gate, CI exercising both
  paths. Apple Silicon stays at 2-lane NEON.
- **Expected:** 30-60 % on AVX-512 hosts (twice the lane count). The
  reduction across genotypes stays scalar.
- **Defer** behind the production-CPU verification needed for the
  `x86-64-v4` build-config change anyway.

#### S6 — Warm-start across records sharing the same `(ploidy, n_alleles, compound_mask)` shape

- Initial `p̂` and `f̂_C` from the previous record's converged values
  when records are neighbouring. EM converges to the same fixed point
  from any positive starting interior point — correctness-neutral if
  the convergence threshold is honoured. No production-caller precedent.
- **Complexity:** thread a previous-solution `Option` through the
  Stage 6 driver; gate by `(chrom_id, distance)` or by batch slot.
  Conflicts with rayon work-stealing — keep within a thread's chunk.

#### S7 — Adaptive genotype-set pruning (Octopus's pattern)

- After 1–2 EM iterations, drop genotypes with posterior < ε from
  subsequent samples. Source: [octopus population_model.cpp:508-551](https://github.com/luntergroup/octopus/blob/master/src/core/models/genotype/population_model.cpp).
- **Expected:** zero on biallelic-diploid (n_g = 3); 2–3× reduction
  in inner-loop trips on tri-/tetra-allelic after iter 2.
- **Complexity:** correctness story — a dropped genotype can't
  recover. Octopus puts it behind a "computationally feasible" gate;
  we'd want a hard floor (e.g. ε = 1e-15) so the kept-set is provably
  ≥ the converged set. ~80 lines.

#### S8 — Cross-record tensor-style batching (same-shape micro-batches)

- **Site:** architectural — would touch `run_em_for_record` and the
  upstream driver, not a single file:line.
- **Idea:** group consecutive records that share `(ploidy, n_alleles,
  compound_mask)` into a batch of N (typical N = 4–256), materialise a
  3-D tensor `log_likelihoods: [batch_idx][sample][genotype]` plus
  per-record state `p_hat: [batch_idx][allele]`, and run the EM loop
  in lockstep across the batch. The hot inner loop becomes a clean
  tensor expression (NumPy-`einsum`-shaped), which gives LLVM a much
  better autovectorisation surface than the current hand-written
  `wide::f64x4` + trait-dispatch shape. Diploid biallelic is the
  dominant shape (~80–90 % of records), so a single batch type covers
  most of the work.
- **Confidence:** Low. The mechanism is real but the headline gain is
  bounded.
- **What this does NOT do:** widen SIMD lanes. We are already at
  `f64x4`; batching N records still runs at 4 lanes — the lanes are
  now (record × sample) cells instead of (sample) cells, but the lane
  count is the same. The direct lane-width lever is S5 (AVX-512
  `f64x8`), not this finding. The NumPy analogy ("batching makes
  per-call overhead vanish") doesn't fully transfer to Rust because
  monomorphised inlined dispatch already has ~zero per-call overhead.
- **What this DOES do:** (a) lets the compiler autovec a clean tensor
  loop more aggressively than today's hand-written SIMD; (b) amortises
  per-record setup (`log_f` / `log_one_minus_f` / pseudocounts /
  shape lookup) — though L5's `RecordScratch` already removes most of
  that cost; (c) reorganises memory access into long stride-1 sweeps
  that are friendlier to L2 / hardware prefetchers if the batch fits.
- **Three costs to weigh:**
  1. **Shape fragmentation.** Records arrive in genome order with shape
     mixing. Same-shape batching needs either a sliding window that
     flushes on shape change (clean, but cuts effective batch size when
     tri-allelic interleave with biallelic) or an out-of-order sort
     (mangles output order; needs sort-back). The dominant biallelic
     stretch is fine; sparse non-biallelic interleaves are the
     pathological case.
  2. **EM iteration mismatch.** Per-record convergence is 3–50 iters;
     lockstep batches stall on the long-tail records. Either keep going
     (~5–10 % wasted FMAs) or mask writes per record (more bookkeeping).
  3. **Transpose memcpy in / out.** For a 1024-record batch ×
     1000 samples × 3 genotypes that's ~24 MB copied per batch. At
     ~50 GB/s memory bandwidth that's ~0.5 ms vs ~18 ms batch EM →
     ~3 % overhead. The batch buffer itself is reused (one allocation
     per engine); the recurring cost is the memcpy, not malloc.
- **Expected:** **5–15 %** wall-clock on the dominant biallelic shape
  if the autovec win is real. Same order as H2 / H4 already in the
  Hot-path band, which are far cheaper to land. The bigger lever for
  cross-record parallelism is **rayon-over-records** (deferred,
  separate plan): 16–128× near-linear core scaling on the production
  server, against the same per-thread SIMD code. Tensor batching only
  matters after rayon to push per-thread efficiency further.
- **Measurement plan (cheap test first):**
  1. **Micro-batch of 4 same-shape records** (so all 4 × `n_samples`
     samples are tiled in L1). One synthetic Criterion fixture; build
     a 3-D tensor by hand, hand-write the inner-loop kernel without
     the engine state machine, compare against the current per-record
     loop on the same input. Threshold: ≥ 5 % wall-clock improvement
     on `cargo bench -- micro_batch_4_records`. If not, the megabatch
     won't be faster either — close the finding here.
  2. If the micro-batch wins: prototype a sliding-window batcher
     keyed on `(ploidy, n_alleles)` in `examples/profile_posterior_engine`
     (not the engine itself yet) and rerun the 64 × 10 000 contam-on
     bench. Threshold: ≥ 10 % wall-clock and same parity-budget pass.
  3. Only at step 3 commit to the engine API change.
- **Complexity cost:** large. The current `PosteriorEngine` is a
  per-record state machine; batching reshapes the public iterator
  surface, the upstream driver contract (Stage 5 → Stage 6 streaming),
  and the convergence / non-finite-error reporting (which record's
  posterior failed?). The bookkeeping for EM-iteration mismatch and
  shape flushing is non-trivial. Estimated 2–4 weeks of work for a 5–15 %
  win, vs ~1–3 days each for H2/H3/H4/L5–L9 totalling a comparable
  range — strongly defer.
- **Order:** evaluate **after** Hot-path findings + rayon-over-records
  + AVX-512 `f64x8` (S5) have landed and we have a flamegraph of the
  remaining hotspots. If at that point the per-thread per-record EM is
  still ≥ 30 % of CPU time and shows clear autovec headroom (uneven
  SIMD utilisation per `cargo asm`), revisit. Otherwise close.

#### S9 — PGO (build-config; see §4 above)

#### S10 — `target-cpu = x86-64-v4` floor (build-config; see §4 above)

#### S11 — Streaming-workload bench variant (cache-cold)

- Single-fixture criterion bench keeps the full input set in L2/L3
  across iterations; production thrashes both inputs and the shape
  cache. Add a streaming variant driven by the live merger. If the
  cache-cold variant is ≥ 15 % slower than the cache-warm baseline,
  current bench numbers under-estimate production wall time.

### Note

- **`#[inline]` const-dispatch on `M::HAS_LANE_4`** at
  [posterior_engine.rs:1457-1477, 1341-1361](../../../src/var_calling/posterior_engine.rs#L1457-L1477)
  is already correct — `M::HAS_LANE_4` is a `const` after
  monomorphisation, DCE-removes the unused branch.
- **Static interp tables (4 kB total)** at
  [interp.rs:51-68](../../../src/var_calling/posterior_engine/interp.rs#L51-L68)
  fit comfortably in L1. No defect.
- **`ExactMath` / `InterpUnivariateMath` / `InterpUnivariateSimdMath`**
  in [backends.rs](../../../src/var_calling/posterior_engine/backends.rs)
  are zero-sized; trait already takes/returns `wide::f64x4` directly
  (the May 18 SIMD analysis report's "trait surface refactor" is
  already done — no action).
- **Dead `_x4` helpers from May 17** are also gone — `interp.rs` is
  237 lines, ends at the test module. No cleanup needed.
- **Out-of-cache shape rebuild** at [shape.rs:72-74](../../../src/var_calling/posterior_engine/shape.rs#L72-L74)
  allocates fresh on every call — production max_alleles = 6 keeps it
  inside the cache; only revisit if a workload pushes past
  `n_alleles ≤ 16` or `ploidy ≤ 8`.

## 6. Out-of-scope observations

- **`thread_local!` `SHAPE_CACHE` interaction with rayon-over-records**
  ([shape.rs:62-65](../../../src/var_calling/posterior_engine/shape.rs#L62-L65)).
  Each worker thread will pay a one-time fill cost. Concurrency category
  for the deferred rayon plan — not actionable in this review.
- **`PosteriorRecord::compound_frequencies: Vec<Option<f64>>`** at
  [posterior_engine.rs:1366-1376](../../../src/var_calling/posterior_engine.rs#L1366-L1376):
  per-record alloc that can't live in scratch (the field outlives the
  next record). If DHAT names it in the top 10 sites, a parallel
  encoding (`Vec<f64>` + `compound_mask`) saves the malloc; otherwise
  leave it.

## 7. What's already good

Three specific patterns to call out:

- **`GenotypeShape` per-`(ploidy, n_alleles)` cache** ([shape.rs:42-60](../../../src/var_calling/posterior_engine/shape.rs#L42-L60))
  carries `log_multinomial_coeffs`, `nonzero_pairs`, `homozygous_allele_for`
  — exactly the precomputations freebayes does at higher per-call cost
  (no shape-keyed cache there). The thread-local slot-array is a tiny,
  branch-free win over a `HashMap`.
- **`MathBackend` trait + monomorphisation** ([backends.rs:53-89](../../../src/var_calling/posterior_engine/backends.rs#L53-L89))
  with `HAS_LANE_4` as a `const` lets the engine pick the SIMD vs scalar
  body at compile time. The unused branch DCEs. The `wide::f64x4`-native
  ln/exp swap landed in two lines without touching the trait surface.
- **Static interp tables in `LazyLock`** ([interp.rs:51-68](../../../src/var_calling/posterior_engine/interp.rs#L51-L68))
  — 257 × 8 B × 2 = 4 kB sits in L1; the IEEE-decomposition trick keeps
  the tables small and the fallback for subnormal / non-finite inputs
  is explicit. The 256-bin tuning (down from the initial 1024) on the
  Stage 5 measurement is a textbook example of profile-driven sizing.

---

### Author response convention

Address each finding by its identifier (e.g. "H1", "L2") with one of:
`applied in <commit>` / `experiment shows no gain — closing` /
`disputed because …` / `deferred to <issue>` / `won't fix because …`.
The "experiment shows no gain" path is expected and welcome — that is
what the measurement plan is for.
