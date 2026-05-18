# Post-H4 sampling profile — posterior engine

**Date:** 2026-05-18.
**Branch:** `perf/posterior-engine-h4-fixation-hoist` @ `013b49f` (one
commit ahead of `main` @ `0170655`).
**Workload:** `POSTERIOR_BACKEND=simd
./target-container/release/examples/profile_posterior_engine` — 30
drains × 10 000-record contam-on biallelic-SNP fixture (64 samples,
default `PosteriorEngineConfig::with_project_defaults()`, three-allele
representative contamination `c_s = 3 %, q_b = [0.6, 0.3, 0.1]`).
**Host:** i7-1260P, `kernel.perf_event_paranoid = 1`, single CPU pinned
with `taskset -c 0` to suppress cross-core scheduler chatter.
**Sampler:** `perf record -F 1999 --call-graph=dwarf -- <binary>`,
reported with `perf report --stdio --no-children --no-inline -F
overhead,symbol -g none`.

This is the post-H4 profile the perf review's measurement plan §3
called for. It supersedes the May-17 flamegraph (pre-SIMD-swap,
pre-H4) for ranking the remaining Likely findings.

## Top self-time frames (single-core, 14.78 G cycles, 6 409 samples)

```
   41.25 %  e_step_simd
   19.37 %  profile_posterior_engine::main         [driver noise]
    6.18 %  _int_malloc
    5.48 %  __memmove_avx_unaligned_erms            [driver noise]
    4.17 %  malloc
    3.12 %  __ieee754_log_fma
    2.30 %  summarise_posteriors
    2.28 %  cfree@GLIBC_2.2.5
    2.18 %  __log10_finite                         [summarise_posteriors's Phred log10]
    1.29 %  _int_free_create_chunk
    1.13 %  _int_free_chunk
    1.04 %  bridge_producer_consumer::helper       [rayon idle workers]
    0.97 %  unlink_chunk.isra.0
    0.90 %  m_step_p_hat
    0.82 %  malloc_consolidate
    0.64 %  __libc_calloc
    0.56 %  unify_alleles                          [per_group_merger one-shot setup]
    0.40 %  _raw_spin_lock                         [kernel]
```

Driver-noise frames (`profile_posterior_engine::main`,
`__memmove_avx_unaligned_erms`) account for ~25 % of total cycles —
they're the per-drain `merged.to_vec()` clone the example does to
preserve fixture state across drains. Not engine cost. The remaining
~75 % is what the production cohort pipeline will actually pay.

## Cross-checked against the cross-core profile

The first profile (unpinned, also captured today at
[tmp/perf_h4.data](../../../tmp/perf_h4.data)) showed `~37 %` of
cycles in `native_queued_spin_lock_slowpath` /
`crossbeam_deque::Stealer::steal` /
`rayon::iter::plumbing::bridge_producer_consumer::helper`. Pinning to
one core dropped that block to ~1 %. The unpinned numbers were
kernel-side scheduler chatter from rayon worker threads waking
periodically on a futex — **not** engine work. The pinned profile is
the trustworthy signal.

## What this confirms and what it refutes from the May-18 perf review

Cross-referenced against
[perf_posterior_engine_2026-05-18.md](../reviews/perf_posterior_engine_2026-05-18.md).

### Confirms

- **H4 worked end-to-end.** The 18.15 → 9.77 µs/record drop (~−46 %
  in this profile run; the back-to-back triple-trial averages were
  −14 % to −32 %) is consistent with the −10–15 % bench threshold
  the review predicted, plus the bonus of running on a warmer thermal
  state than the prior baseline.
- **H3 (softplus_neg table) demotes to Speculative,** exactly as the
  review forecast ("Order *after* H4: H4 may demote H3 to Speculative
  by removing most of the call sites"). Transcendental cost is now
  `__ieee754_log_fma 3.12 % + __log10_finite 2.18 % = ~5.3 %`, and
  of that, only the `log_fma` slice is in the EM hot path
  (`log10_finite` is `summarise_posteriors`'s Phred conversion,
  outside the EM loop). A softplus_neg table replacing the few
  remaining `log_sum_exp_2_x4` calls (in the *heterogeneous*-fixation
  branch only, since H4 already short-circuits the homogeneous case)
  would shave at most ~2–3 % wall-clock on the default-config
  workload — not worth the 8 kB + parity-budget churn.
- **L4 (`#[cold]` on `classify_nonfinite`) — drop.** The symbol does
  not appear in the top 30. Branch predictor already handles it; the
  layout is fine.
- **L2 (`InterpUnivariateMath::HAS_LANE_4 = true`) — defer.** Only
  affects the non-default scalar interp backend; doesn't move the
  default-config profile.
- **L10 (in-bench `debug_assert!(diagnostics.iterations …)`) — keep.**
  One-line discipline irrespective of the profile.

### Refutes / promotes

- **L5–L7 (per-engine `RecordScratch`) → Hot-path.** The review
  ranked these as Likely pending a DHAT pass. The sampling profile is
  enough to promote: aggregate allocator self-time
  (`_int_malloc + malloc + cfree + _int_free_chunk + _int_free_create_chunk +
  malloc_consolidate + __libc_calloc + unlink_chunk`) is **~16 % of
  total cycles** — second only to `e_step_simd` itself, and the only
  remaining single category with a clear, contained fix. Expected
  win: 5–10 % wall-clock, plus the engine becomes friendlier to
  future rayon-over-records work because the steady-state malloc
  count drops to ~zero.

### Still pending / unmoved

- **L3 (`#[inline(always)]` on `log_sum_exp_*_x4`)** — profile doesn't
  separate this from the inlined parent `e_step_simd`; `cargo asm`
  remains the right check, but not a priority while
  RecordScratch is on the table.
- **L8/L9 (transpose `log_likelihoods` to batch-of-4 layout)** —
  `e_step_simd` is 41.25 % of cycles, which is where the transpose
  would help, but the profile doesn't decompose what fraction of
  that is the strided gather vs the SIMD compute vs the splat / lane
  load. A flat-symbol view can't tell us. Defer until RecordScratch
  lands; if `e_step_simd` is still ≥ 30 % afterwards, a finer-grained
  profile (per-line, or `perf annotate`) is the next measurement
  step.

## Where the cycles go (engine-attributable, normalised)

Removing the driver noise (`main` + `memmove` ≈ 25 %) and the
one-shot fixture setup (`unify_alleles` 0.56 %) leaves ~74 % of
cycles attributable to the engine work proper. Within that envelope:

| Category | Share of engine cycles (approx.) |
|---|---|
| `e_step_simd` (incl. its SIMD ln/exp via `wide::f64x4`) | 56 % |
| glibc allocator (malloc / free family) | 22 % |
| Transcendentals outside `_simd` (`__ieee754_log_fma`, `__log10_finite`) | 7 % |
| `summarise_posteriors` (Phred GQ / QUAL) | 3 % |
| `m_step_p_hat` + other M-steps | < 2 % |
| Misc + rayon idle | residual |

The headline: post-H4, **the engine is allocator-bound after EM
compute, and transcendental cost is a tail term.** Pre-SIMD, the May-17
flamegraph put ln/exp at ~40 %; post-`wide`-native + H4 that's
shrunk to ~7 % and the allocator has taken its place at ~22 % of
engine cycles.

## Wall-time corroboration (triple back-to-back, perf-off)

From the H2+H4 commit message (`013b49f`), measured by
`profile_posterior_engine` without perf instrumentation:

```
trial 1: 14.59 → 9.89  µs/record  (−32.2 %)   [BEFORE = main; AFTER = H4]
trial 2: 17.30 → 12.69 µs/record  (−26.6 %)
trial 3: 14.84 → 12.80 µs/record  (−13.7 %)
```

The cross-trial spread reflects the powersave-governor caveat
([feedback_bench_cpu_governor](../../../home/jose/.claude/projects/-home-jose-devel-join-per-sample-vcfs/memory/feedback_bench_cpu_governor.md)).
Within-trial pairs are the trustworthy signal.

## Profile artefacts

Saved for re-analysis:

- `tmp/perf_h4.data` — unpinned (rayon-noisy) record on the
  cpu_atom / E-core event slot.
- `tmp/perf_h4_pinned.data` — single-core pinned (`taskset -c 0`),
  cpu_core / P-core event slot. **This is the trustworthy one.**
- `tmp/profile_posterior_engine_after_h4` — release binary of the
  H4 commit (`013b49f`), used for both `perf record` runs and for
  the wall-time triples above.
- `tmp/profile_posterior_engine_before_h4` — release binary of
  `main` (`0170655`), used as the wall-time baseline.

## Recommended next move

Per the review's update protocol, **promote L5–L7
(`RecordScratch`) from Likely to Hot-path** and apply next.
Expected win: 5–10 % wall-clock on the default-config contam-on
workload. After it lands, re-profile and decide whether L8/L9
(transpose) is still worth pursuing.
