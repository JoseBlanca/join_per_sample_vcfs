# Stage 6 — posterior engine performance (review-application wave 2)

**Date:** 2026-05-18 (afternoon session).
**Branch:** all work landed on `main` after fast-forward merges.
**Tip:** `9594533`.

Picks up where the morning's review report
([perf_posterior_engine_2026-05-18.md](../reviews/perf_posterior_engine_2026-05-18.md))
left off. That report ran four parallel per-category sub-agents +
deep web + local-source research and verdicted "Run experiments" with
4 Hot-path, 11 Likely, 11 Speculative findings. This session applied
the top-of-list Hot-path items, re-profiled, and used the new profile
to reshape the priority list for next session.

## Cumulative result

| Snapshot | µs/record (contam-on, 64×10 000) | Source |
|---|---|---|
| May-17 (pre-SIMD swap) | ~49.7 | [perf 2026-05-17](./posterior_engine_perf_2026-05-17.md) |
| May-18 morning (post `wide::f64x4` swap) | 24.41 / 22.33 / 18.15 (Exact / Interp / SIMD) | [simd analysis](./posterior_engine_simd_analysis_2026-05-18.md) |
| **May-18 afternoon, post H4** | 14.59→9.89 µs (best trial), 17.30→12.69 (warm) | `013b49f` |
| **May-18 afternoon, post RecordScratch** | 9.06–9.33 µs (perf-on, single-core pinned) | `9594533` |

**End-to-end this session:** ~12–18 µs/record before the wave → ~9–10 µs after.
Subject to powersave-governor noise — the cycle-distribution
comparison below is the trustworthy signal.

## What landed

Three commits, all on `main`:

### 1. `013b49f` — H2 + H4: log_indep cache + homogeneous-fixation hoist

**H2** (review §5 Hot-path) — cache `log_indep[g] = log_multinomial_coeffs[g]
+ Σ k·log_p_effective[a]` in `EmScratch::log_indep_per_g` once per
E-step call. Both `e_step` (scalar) and `e_step_simd` (batch + tail)
now read from the cache instead of rebuilding the dot product per
(sample / batch, genotype).

**H4** (review §5 Hot-path) — when `fixation_index_overrides = None`
and every sample shares the same `f_s` (default config:
`fixation_index_default = 0.0`), the entire per-genotype log-prior is
sample-invariant. Precompute it into `EmScratch::log_prior_per_g` once
per EM iteration; the hot inner loop splats one scalar per genotype
instead of running `log_sum_exp_2_x4` (3 SIMD transcendentals) per
(batch, homozygous genotype). H1 (the missing SIMD `-∞` short-circuit)
is subsumed by H4.

**Measurement** — three interleaved BEFORE / AFTER trials in one
shell session:

| Trial | BEFORE (main) | AFTER (H4) | Δ |
|---|---|---|---|
| 1 | 14.59 µs | 9.89 µs | −32.2 % |
| 2 | 17.30 µs | 12.69 µs | −26.6 % |
| 3 | 14.84 µs | 12.80 µs | −13.7 % |

Conservative bound: −13.7 %. Clean win.

### 2. `666a432` — Post-H4 sampling profile (docs)

`POSTERIOR_BACKEND=simd taskset -c 0 perf record -F 1999
--call-graph=dwarf` of the H4 binary. Single-core pinning was
essential: an unpinned run had ~37 % cycles in
`native_queued_spin_lock_slowpath` /
`rayon::iter::plumbing::bridge_producer_consumer::helper` — kernel
scheduler chatter from idle rayon worker threads in the fixture
builder, not engine work. Pinned dropped that to ~1 %.

Pinned-profile findings (against the review's Likely list):

| Review finding | Pre-profile band | Post-profile decision |
|---|---|---|
| H3 (softplus_neg table) | Hot-path | **Demote to Speculative** — `__ieee754_log_fma` only 3.12 % of cycles after H4 |
| L5–L7 (per-engine RecordScratch) | Likely | **Promote to Hot-path** — allocator ~16 % of engine cycles |
| L4 (`#[cold]` on `classify_nonfinite`) | Likely | **Drop** — symbol invisible (< 0.3 %) |
| L2 (`InterpUnivariateMath::HAS_LANE_4 = true`) | Likely | **Defer** — non-default backend only |
| L3 (`#[inline(always)]` on `_x4` helpers) | Likely | **Keep, low priority** — needs `cargo asm` check |
| L10 (in-bench `debug_assert` on iters) | Likely | **Keep** — one-line discipline |

Saved as `posterior_engine_post_h4_profile_2026-05-18.md`.

### 3. `9594533` — RecordScratch lift (L5–L7)

Lifted **21 per-record allocations** into a `RecordScratch` field on
`PosteriorEngine`. Replaces `EmScratch::new(...)` + 8 satellite
`Vec::collect`/`vec![…; n]` + 5 inner allocations + 5 output Vecs
with a single `scratch.resize_to(...)` at function entry.
Steady-state records of the same shape allocate ~zero.

Structural changes worth noting beyond the cycle savings:

- **`EmContext<'_>` is now Copy and slice-free** (just scalars + `&GenotypeShape`).
  The slice data (`compound_mask`, `pseudocounts`,
  `log_f_per_sample`, `log_one_minus_f_per_sample`) moved into
  scratch; EM reads them through the same `&mut RecordScratch` it
  already receives. Sidesteps the borrow conflict between holding
  `&EmContext` and `&mut scratch` simultaneously.
- **M-step uses double-buffering**: `m_step_p_hat` writes into
  `scratch.p_hat_next`; `run_em_loop` `mem::swap`s the two buffers
  after the convergence test.
- **Mixture pre-pass returns `Result<(), _>`** and fills
  `scratch.mixture_log_likelihoods`. `run_em_for_record` then
  `mem::swap`s that buffer into the local `log_likelihoods` owner so
  EM reads it via a borrow disjoint from `&mut scratch`.
- **Output Vecs are `.clone()`d from scratch** into `PosteriorRecord`.
  5 mallocs + memcpys per record — the minimum compatible with the
  owned `PosteriorRecord` contract.

**Measurement** — cycle distribution comparison (single-core pinned):

| Symbol | H4 baseline | post-RecordScratch | Δ |
|---|---|---|---|
| `_int_malloc` | 6.18 % | 2.81 % | **−3.37 pp** |
| `_int_free_create_chunk` | 1.29 % | 0.57 % | **−0.72 pp** |
| `cfree@GLIBC` | 2.28 % | 1.79 % | −0.49 pp |
| `malloc` | 4.17 % | 4.19 % | ~0 |
| `_int_free_chunk` | 1.13 % | 1.43 % | +0.30 pp |
| `unlink_chunk.isra.0` | 0.97 % | 1.66 % | +0.69 pp |
| `__libc_calloc` | 0.64 % | — (<0.3 %) | dropped below threshold |
| **Allocator total** | **~17.5 %** | **~13.3 %** | **~−4.2 pp** |
| `__memmove_avx_unaligned_erms` | 5.48 % | 6.16 % | +0.68 pp (output clones) |
| `e_step_simd` | 41.25 % | 42.20 % | ~unchanged (relative-share) |

**Wall-time was inconclusive** on this dev laptop —
`feedback_bench_cpu_governor`'s thermal-state caveat fully active.
Triple back-to-back ranged −17 % to +12 % across trials. The cycle
distribution above is the trustworthy comparison: ~4 pp of cycles
freed from the allocator, ~0.7 pp refunded to memmove for the
output clones, net ~3.5 pp.

Honest disclosure: the review predicted "5–10 % wall-clock". The
actual was smaller — partly because the allocator was already in the
glibc tcmalloc fast path (small-bin reuse is cheap), partly because
the 5 mandatory output clones refund some of the savings.

**Strategic value beyond the cycle count:**

- Per-engine ownership of all per-record state is the prerequisite
  for **rayon-over-records** — each worker thread gets its own
  `RecordScratch` with no shared mutable state. This was the
  single-biggest motivation for landing it.
- `EmContext<'_>` (Copy, scalar-only) makes future per-record-metadata
  additions cheap — no new function-arg plumbing.

## Post-RecordScratch profile + hardware counters

Re-profiled after `9594533` lands. Pinned single-core `perf record`
captured 6063 samples / 13.9 G cycles.

```
   42.20 %  e_step_simd
   20.91 %  profile_posterior_engine::main         [driver noise]
    6.16 %  __memmove_avx_unaligned_erms            [driver noise + 5 output clones]
    4.19 %  malloc
    3.06 %  __ieee754_log_fma
    2.97 %  summarise_posteriors
    2.81 %  _int_malloc
    2.63 %  __log10_finite                         [summarise_posteriors's Phred log10]
    1.79 %  cfree@GLIBC
    1.66 %  unlink_chunk.isra.0
    1.43 %  _int_free_chunk
    1.25 %  m_step_p_hat
    0.86 %  bridge_producer_consumer::helper       [rayon idle, residual]
    0.84 %  malloc_consolidate
    0.65 %  unify_alleles                          [per_group_merger one-shot]
```

**Hardware-counter analysis** (perf stat, same single-core pinned
workload):

```
       2,996254941 seconds time elapsed

    13.865.927.036  cpu_core/cycles/
    41.172.115.378  cpu_core/instructions/           #    2,97  insn per cycle
     5.388.257.983  cpu_core/branches/
         5.747.219  cpu_core/branch-misses/          #    0,11% of all branches
    12.252.137.697  cpu_core/L1-dcache-loads/
        41.088.117  cpu_core/L1-dcache-load-misses/  #    0,34% of all L1-dcache accesses
        26.311.232  cpu_core/LLC-loads/
        23.401.146  cpu_core/LLC-load-misses/        #   88,94% of all LL-cache accesses
         2.007.549  cpu_core/dTLB-load-misses/
```

The numbers tell a sharp story:

- **IPC = 2.97** — near-saturating superscalar execution. The
  i7-1260P P-core architectural max is ~5; sustained ~3–4 is
  typical for tight numeric code. The CPU is *not* stalling on
  memory.
- **L1-dcache-load-miss rate = 0.34 %** — the engine's hot working
  set fits in L1. The 64 sample × 3 genotype `log_likelihoods` (1.5
  KB), the `posteriors` output (1.5 KB), the scratch buffers (few
  KB total), and the interp tables (4 KB) all coexist within the
  48 KB L1d.
- **Branch misses = 0.11 %** — predictor is happy; control flow
  is straight-line at the cycle level.
- **LLC pressure essentially absent** — 26 M LLC loads vs 12 B
  L1 loads is 0.21 %.

## What this changes for the priority list

The May-18 review's Likely band ranked **L8/L9 (transpose
`log_likelihoods` to batch-of-4 layout)** as a candidate to attack
`e_step_simd`'s 41 % cycle share. The whole pitch was "strided
gather hits multiple cache lines per lane". But:

- The L1 hit rate is 99.66 %.
- LLC traffic is negligible.
- IPC is at 2.97 (near-saturated).

**The transpose would re-arrange data that's already arriving in
time.** It would not move IPC up — execution units are already
busy. Strong empirical refutation of L8/L9's predicate; should drop
from Likely to **Note** in the perf review.

What's actually left as a single-thread lever:

| Lever | What it attacks | Expected | Cost |
|---|---|---|---|
| Rayon-over-records | Parallelism (not single-thread) | 16–128× on a server | Small (RecordScratch already prep'd this) |
| AVX-512 `f64x8` backend (S5) | Wider SIMD = fewer instructions | ~2× on capable hosts | New backend, multiversion crate, CI matrix |
| H3 softplus_neg table | `log_sum_exp_2_x4` transcendentals | < 2 % (only heterogeneous-fixation path) | 8 KB table, parity-budget check |
| SQUAREM acceleration | EM iteration count | Only helps cap-bound records | ~100 LOC, log-likelihood-sum pass needed |
| L3 `#[inline(always)]` on `_x4` helpers | Codegen | < 1 % if anything | Attribute swap, `cargo asm` check |

Of these, only **rayon-over-records** is a clear order-of-magnitude
move; the rest are at most 2× per-thread and add real complexity.

## Updated review-finding status

Reflects what we learned this session against
[perf_posterior_engine_2026-05-18.md](../reviews/perf_posterior_engine_2026-05-18.md):

| Finding | Pre-session band | Status |
|---|---|---|
| **H1** (SIMD `-INF` short-circuit) | Hot-path | **Applied** (subsumed by H4) — `013b49f` |
| **H2** (cross-batch `log_indep` reuse) | Hot-path | **Applied** — `013b49f` |
| **H3** (softplus_neg table) | Hot-path | **Demoted to Speculative** — transcendental cost now ~3 % |
| **H4** (homogeneous-fixation hoist) | Hot-path | **Applied** — `013b49f`, the headline win |
| **L1** (autovec branch in `accumulate_expected_counts` general path) | Likely | unchanged — applies to ≥3-allele records only |
| **L2** (`InterpUnivariateMath::HAS_LANE_4 = true`) | Likely | **Deferred** — non-default backend |
| **L3** (`#[inline(always)]` on `_x4`) | Likely | unchanged — cargo asm check pending |
| **L4** (`#[cold]` on `classify_nonfinite`) | Likely | **Dropped** — symbol invisible in profile |
| **L5–L7** (per-engine `RecordScratch`) | Likely | **Applied** — `9594533`, ~4 pp allocator share saved |
| **L8** (transpose `log_likelihoods`) | Likely | **Demoted to Note** — L1 hit 99.66 %, no room for cache-layout wins |
| **L9** (transpose `posteriors` scatter) | Likely | **Demoted to Note** — same reason as L8 |
| **L10** (in-bench `debug_assert` on iters) | Likely | unchanged — one-line discipline, still worth doing |
| **L11** (`iter_batched` profile-contamination doc) | Likely | unchanged |
| **S4** (SQUAREM) | Speculative | unchanged |
| **S5** (AVX-512 `f64x8` backend) | Speculative | **Promoted to Likely** — IPC=2.97 says wider lanes are the only single-thread lever left |
| **S6** (warm-start across records) | Speculative | unchanged |
| **S7** (Octopus-style genotype pruning) | Speculative | unchanged |
| **S8** (cross-record tensor-style batching) | Speculative | unchanged — see HW counters: cache layout already maximal |
| **S9** (PGO) | Speculative | unchanged |
| **S10** (`x86-64-v4` floor) | Speculative | unchanged |
| **S11** (streaming-workload bench) | Speculative | unchanged |

## Suggested next moves

In rough priority order, deferring per the user's call:

1. **Rayon-over-records.** Order-of-magnitude lever (16–128× on a
   server). The current per-thread cost (9 µs) stays put; throughput
   scales near-linearly with core count. `RecordScratch` is already
   shaped for per-thread ownership. *Deferred this session per the
   user's call — separate plan.*
2. **AVX-512 `f64x8` backend (was S5, now Likely).** ~2× per-thread
   on capable hosts. Pairs naturally with rayon — same code,
   different lane width.
3. **L10 + L3 mechanical cluster.** Tiny PRs; land alongside any
   future change.
4. **Re-evaluate after rayon lands.** Profile shape may shift once
   the throughput is in a different regime.

## Artefacts retained

In `tmp/` (gitignored):

- `perf_h4_pinned.data` — H4 baseline single-core profile.
- `perf_recordscratch.data` / `perf_post_rs_v2.data` —
  post-RecordScratch profiles (two captures, consistent).
- `samply_post_rs.json.gz` — samply profile of the
  post-RecordScratch binary (185 KB gzipped JSON;
  https://profiler.firefox.com).
- `profile_posterior_engine_before_h4` / `..._after_h4` /
  `..._after_recordscratch` — release binaries from each major
  commit, kept for any future bisect / comparison run.

## Files touched this session

Three commits' worth, all on `main`:

- [src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs) — H2 + H4 + RecordScratch refactor (~960 line diff across the two perf commits).
- [src/var_calling/posterior_engine/shape.rs](../../../src/var_calling/posterior_engine/shape.rs) — `#[allow(dead_code)]` on `GenotypeShape::n_genotypes` (test-only read after RecordScratch removed the runtime use).
- [doc/devel/reports/implementations/posterior_engine_post_h4_profile_2026-05-18.md](posterior_engine_post_h4_profile_2026-05-18.md) — pinned-profile capture between H4 and RecordScratch.
- This file — wave-2 session report.
- [PROJECT_STATUS.md](../../../PROJECT_STATUS.md) — Open-items and last-completed-task block.

## Validation matrix

For each commit, inside `./scripts/dev.sh`:

| Gate | `013b49f` (H2+H4) | `9594533` (RecordScratch) |
|---|---|---|
| `cargo fmt --all -- --check` | clean | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean | clean |
| `cargo test --lib` (passing) | 774 | 774 |
| `cargo test --all-targets --all-features` integration | 113 | 113 |
| 12-fixture parity battery | passes (50× margin) | passes (50× margin) |
| Proptests (simplex / permutation / sum-to-one) | pass at any (ploidy, n_alleles) | pass at any (ploidy, n_alleles) |
