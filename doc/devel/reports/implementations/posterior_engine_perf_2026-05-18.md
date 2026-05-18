# Stage 6 — posterior engine performance (samply-guided cycle)

Date: 2026-05-18.

Picks up where the 2026-05-17 SIMD analysis
([posterior_engine_simd_analysis_2026-05-18.md](./posterior_engine_simd_analysis_2026-05-18.md))
left off — that work landed the SIMD backend as the default and
revealed that `wide`'s native `f64x4::ln` / `f64x4::exp` were the
right primitives. This session uses the samply profile of that new
default state to guide three concrete optimisations.

The posterior engine source is at
[src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs).

## Cumulative result (this session)

Triple back-to-back inside one shell session on biallelic-contam-on,
64 samples × 10 000 records × 30 runs (the `examples/profile_posterior_engine`
driver):

| Backend | µs/record (contam-on) | vs `ExactMath` | vs `InterpUnivariateMath` |
|---|---|---|---|
| `ExactMath` | ~24.8 | — | +10 % |
| `InterpUnivariateMath` (scalar interp) | 22.4 | −10 % | — |
| **`InterpUnivariateSimdMath` (default)** | **~15.0** | **−40 %** | **−33 %** |

Versus the *yesterday-end* SIMD baseline (the post-default-flip state
of commit `c0f05bd`, which already had the wide-native lane-of-4
primitives), this session shaved another **~6 µs/record** (~30 %) off
through allocator + algorithmic wins, all single-threaded.

Accuracy: the `interp_univariate_simd_clears_approximate_parity_budget`
harness passes unchanged after every commit. 50× margin under the 1e-4
parity budget holds.

## Step-by-step

Three commits landed; one substantial refactor was tried and
abandoned.

| Commit | Step | Outcome |
|---|---|---|
| `5058709` | (carry-over from yesterday) drop dead `interp::*_x4` + tighten `MathBackend::{ln_x4, exp_x4}` to `f64x4` signature | code clarity; perf-neutral within noise (16-byte binary delta) |
| samply record | **profile** the new SIMD default | found 4 actionable hot regions: malloc 20 %, scalar mixture pre-pass 12 %, `m_step_p_hat` inner loop 7 %, scattered scalar work |
| `04bdd9e` | (1) hoist `compute_mixture_log_likelihoods` per-sample allocs | **−6.2 %** |
| `(abandoned)` | (2) engine-owned `RecordScratch` arena reused across records | wash-to-mild-regression (−2.8 % POST, 14-pair A/B); reverted |
| `b7dc557` | (3) SIMD mixture pre-pass (`compute_mixture_log_likelihoods_simd`) | **−7.2 %** |
| `6917b08` | (4) biallelic-diploid `m_step_p_hat` fast path | **−6.5 %** |
| samply re-record | **profile** the new state | confirmed SIMD math intrinsics now dominate at 68 %; remaining hot lines all < 1 % each |

Stacked the three wins compound: end-of-day is ~30 % faster than
the start-of-day SIMD default on the contam-on bench.

## Empirical findings

Four surprises reshape what we'd guess for the next cycle.

### 1. The arena refactor was a wash, not the −20 % we expected

The samply profile showed ~20 % of total in libc malloc/memcpy, so an
engine-owned arena reusing every per-record `Vec` looked like an easy
20 %. It wasn't.

A complete 12-field `RecordScratch` struct passed through
`run_em_for_record` (plus a top-level destructure to satisfy the
borrow checker against the live `log_likelihoods` slice) measured
**~2.8 % slower** in POST over 14 A/B pairs. Probable causes:

- Glibc's malloc tcache makes small (16–256 byte) allocations
  effectively free for the size class our per-record `Vec`s fall into.
- Most of the samply "malloc" cost was the *initialization* memcpy
  from `vec![0.0; n]`, not malloc/free traffic.
- The added complexity (12-arg `compute_mixture`, destructure-at-top,
  more `&mut Vec` plumbing through `run_em_loop` / `e_step` /
  `e_step_simd`) probably mildly pessimised register allocation.

The narrower hoist (the `n_obs` / `mean_err` per-sample-loop allocs
in `compute_mixture`, commit `04bdd9e`) *did* pay — those were 128
allocs per record, two orders of magnitude more volume than the rest.
**Lesson: malloc volume matters more than total bytes; the long tail
of "one Vec per record" sites doesn't move the needle.**

The arena diff (~350 lines) is discarded. The `m_step_p_hat` /
`m_step_f_hat_compound` half of it (write-into-buffer instead of
`.collect`) went with it, since on its own it was perf-neutral.

### 2. samply was right about the mixture pre-pass

samply predicted ~12 % combined in `compute_mixture_log_likelihoods`
(3.8 % at the inner `math.ln(mix)` + ~8 % in scalar `interp::*_approx`
called from there). The SIMD rewrite measured **−7.2 %**; the gap
(samply's 12 % minus measured 7 %) is the new pre-pass overhead
(two extra Vecs per record for `c_s_all` / `mean_err_all`) plus the
`fallback_log_likelihoods.to_vec()` initial copy. Both are candidates
for a future hoist if the rayon work doesn't make them irrelevant
first.

### 3. The biallelic-diploid fast path captured exactly samply's range

samply predicted ~5-7 % for `m_step_p_hat`'s inner triple loop. The
hand-coded biallelic-diploid path (no inner loop, no branch on
`if k != 0`, just two FMAs per sample) measured **−6.5 %** — square
in the predicted band. POST won 14 of 16 A/B pairs.

The general path is unchanged for polyploid / multiallelic shapes,
so this is a pure dominant-case specialization.

### 4. The post-cycle profile shape says single-threaded is mostly done

After three wins, the new samply profile shows:

- **SIMD math intrinsics (avx + fma + avx2): 68 %** — `wide::f64x4::ln`
  and `f64x4::exp` polynomial work (~30 FMAs per call). Productive
  cycles; can only shrink via algorithmic change (3D `log_mix` LUT,
  fewer-term polynomial via sleef) or via multi-core.
- **libc malloc/memcpy: 21 %** — virtually all per-record output
  Vecs (posteriors, p_hat, best_genotype, gq_phred, compound_freqs)
  plus the new SIMD-mixture pre-pass Vecs. Hard to attack without
  changing the `PosteriorRecord` ownership model.
- **Everything else: < 1 % per line, < 11 % total.**

The micro-optimisation well is drained. The next big lever is
**rayon over records** — near-linear core scaling, multiplicative on
top of what we have.

## Decisions made

1. **Three commits ship; arena commit abandoned.** Each landed commit
   measured a real, A/B-confirmed win (POST wins ≥ 11 of 16 pairs in
   each case). The arena work was reverted intact rather than being
   committed as "code clarity" — the complexity-vs-perf trade-off
   wasn't justified for a neutral result.
2. **All four originally-samply-identified opportunities are now
   acted upon.** Wins (1, 3, 4) shipped; (2) was empirically
   rejected. We don't have a fifth obvious target.
3. **Next session: rayon over records.** The post-cycle samply
   profile makes the case clean — 68 % of cost is genuine SIMD math
   that can only be parallelised at the record granularity, not
   sped up further per-thread.

## Infrastructure delivered

- `compute_mixture_log_likelihoods_simd<M: MathBackend>` in
  [posterior_engine.rs:843-1054](../../../src/var_calling/posterior_engine.rs#L843-L1054)
  — lane-of-4 SIMD pre-pass with per-lane skip masks for the
  `c_s <= 0` and `n_a == 0` cases, defensive `mix <= 0` lane-level
  error reporting, and a scalar tail for `n_samples % 4`. Dispatched
  in `run_em_for_record` when `M::HAS_LANE_4`.
- `accumulate_expected_counts` helper in
  [posterior_engine.rs:1818-1888](../../../src/var_calling/posterior_engine.rs#L1818-L1888)
  factoring out `m_step_p_hat`'s inner triple loop with a
  biallelic-diploid fast path (no branch, two FMAs per sample).
- Two narrower allocations hoisted in
  [compute_mixture_log_likelihoods](../../../src/var_calling/posterior_engine.rs#L765-L770)
  (`n_obs`, `mean_err`) so the EM hot path no longer pays
  `2 × n_samples` mallocs per record.
- samply tooling pipeline: `addr2line` + a Python aggregator script
  reusing the syms.json sidecar for libc-side resolution. Saved at
  `/tmp/aggregate_v2.py` for the next cycle. Two profiles archived
  at `/tmp/profile_simd.json.gz` (start-of-session) and
  `/tmp/profile2.json.gz` (post-three-wins).

## Validation

All run inside the dev container:

```sh
./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings
./scripts/dev.sh cargo test --release --lib var_calling::posterior_engine
./scripts/dev.sh cargo test --test posterior_engine_integration
./scripts/dev.sh cargo build --release --example profile_posterior_engine
```

85 lib tests + 7 integration tests pass. `cargo fmt --check` shows
pre-existing diffs in regions this session didn't touch — out of
scope for these commits.

A/B benches all done back-to-back inside one shell session
(per `feedback_bench_cpu_governor` — CPU swings ~10-15 % between
thermal states, so all comparisons are inside one warmed-up window).

## Branch state

`perf/posterior-samply` is ready to merge to main. Commits:

| SHA | Commit |
|---|---|
| `04bdd9e` | `posterior_engine: hoist compute_mixture per-sample allocations` |
| `b7dc557` | `posterior_engine: SIMD mixture pre-pass for HAS_LANE_4 backends` |
| `6917b08` | `posterior_engine: biallelic-diploid m_step_p_hat fast path` |

(The branch's first commit `5058709` was already on main from
yesterday's work via the `perf/simd-analysis` → main merge.)

## Out-of-scope follow-ups (in priority order)

1. **Rayon over records** — the next session's headline. Near-linear
   core scaling on a 4–8-core machine. The per-record EM is already a
   pure function over `&MergedRecord` plus engine state; the work is
   wiring `par_iter` and threading thread-local scratch arenas.
   Expected: 4–8× on multi-core production hardware.
2. **3-D `log_mix` table for the SIMD mixture pre-pass** — still
   deferred from the 2026-05-17 plan. Per-`(p_own, p_contam, c_s)`
   tuple, the `ln(mix)` lookup could be a 3-D linear-interp gather
   instead of a polynomial. Only worth prototyping if the mixture
   pre-pass is still the bottleneck after rayon. Likely irrelevant
   then.
3. **Faster SIMD `ln` / `exp`** — `sleef` or hand-tuned shorter
   polynomials. Trading ~10× accuracy margin (we have 50×) for fewer
   FMA cycles. Moderate effort, ~10-20 % per-thread gain estimate.
4. **`PosteriorRecord` ownership rework** — the remaining 21 %
   malloc is dominated by output Vecs (posteriors, p_hat,
   best_genotype, gq_phred, compound_frequencies). A
   `&mut PosteriorRecord` output parameter pattern would let the
   caller reuse buffers across records. Invasive (touches the public
   API), low expected payoff (~5 %).
5. **Scattered scalar work**: per-lane finiteness check + scatter in
   `e_step_simd` (0.6 %), `summarise_posteriors` GQ `log10` (0.5 %).
   Each < 1 % individually; not worth chasing in isolation.

## File touch list

- [`src/var_calling/posterior_engine.rs`](../../../src/var_calling/posterior_engine.rs)
  — hoist (in `compute_mixture_log_likelihoods`); add
  `compute_mixture_log_likelihoods_simd` + dispatch in
  `run_em_for_record`; extract `accumulate_expected_counts` with
  biallelic-diploid fast path.
- This report.

## Reproducing the benches

```sh
./scripts/dev.sh cargo build --release --example profile_posterior_engine

# Triple back-to-back to control for CPU-state drift:
POSTERIOR_BACKEND=exact  ./target-container/release/examples/profile_posterior_engine 2>&1 | tail -1
POSTERIOR_BACKEND=interp ./target-container/release/examples/profile_posterior_engine 2>&1 | tail -1
POSTERIOR_BACKEND=simd   ./target-container/release/examples/profile_posterior_engine 2>&1 | tail -1
```

Expected on the 12th-gen Intel Core i7-1260P used this session
(thermal state varies; rerun the triple back-to-back if numbers
swing):

```
exact:  ~24 µs/record
interp: ~22 µs/record
simd:   ~15 µs/record
```

## Reproducing the samply profile

```sh
POSTERIOR_BACKEND=simd samply record \
    --save-only --no-open --unstable-presymbolicate \
    -r 2000 -o /tmp/profile.json.gz \
    ./target-container/release/examples/profile_posterior_engine

gunzip -fc /tmp/profile.json.gz > /tmp/profile.json
# Then aggregate hot source lines via the script saved at /tmp/aggregate_v2.py.
```

The host runs rootless podman (samply must be invoked outside the
container; `perf_event_paranoid=2` is enough for user-space sampling).
