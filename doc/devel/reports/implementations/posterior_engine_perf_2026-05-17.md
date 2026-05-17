# Stage 6 — posterior engine performance work (implementation report)

Date: 2026-05-17.

Implements the perf work described in
[doc/devel/implementation_plans/posterior_engine_diploid_lookup_tables.md](../../implementation_plans/posterior_engine_diploid_lookup_tables.md),
plus the rayon / SIMD follow-up extensions added during planning. The
posterior engine source is at
[src/var_calling/posterior_engine.rs](../../../src/var_calling/posterior_engine.rs)
with new submodules [backends.rs](../../../src/var_calling/posterior_engine/backends.rs),
[interp.rs](../../../src/var_calling/posterior_engine/interp.rs), and
[shape.rs](../../../src/var_calling/posterior_engine/shape.rs).

## Cumulative result

| Backend | µs/record (contam-on, 64 samples) | vs ExactMath |
|---|---|---|
| ExactMath | ~49.7 | — |
| **InterpUnivariateMath (default)** | **~44.0** | **−11.5 %** |
| InterpUnivariateSimdMath | ~46.0 | −7.4 % |

Production default is now `InterpUnivariateMath`. `ExactMath` remains
reachable for bit-identical reproducibility via
`PosteriorEngine::with_math_backend(..., ExactMath)`. The SIMD backend
ships behind the same opt-in for future experiments — see §"SIMD
deferral".

Accuracy on the 12-fixture battery (`tests_math_backend_accuracy`):

|  | posterior max | phred max | allele_freq max | argmax mismatch |
|---|---|---|---|---|
| Budget | 1e-4 | 1.0 | 1e-4 | 0 |
| `InterpUnivariateMath` (256 bins) | 2.0e-6 | 0.014 | 4.0e-7 | 0 |
| `InterpUnivariateSimdMath` | 2.0e-6 | 0.014 | 4.0e-7 | 0 |

50–250× margin under every parity budget.

## Step-by-step

The plan staged 8 steps; all 8 landed plus three tuning / extension
commits.

| Commit | Step | Notes |
|---|---|---|
| `e836a2b` | Plan v0 | interpolated math LUTs, hybrid 1D/3D |
| `9d8f0e3` | Plan v1 | e_step restructure + rayon/SIMD follow-ups |
| `44bf1ce` | 0 — baseline bench + flamegraph | confirmed `ln`/`exp` ≈ 40 %, `e_step` ≈ 30 % self-time |
| `57ee24d` | 1 — `MathBackend` trait + engine generalisation | bit-identical plumbing |
| `1277559` | 2 — `safe_ln(f_s)` hoist | bit-identical |
| `f82720d` | 3 — shape cache + e_step restructure | bit-identical; `nonzero_pairs`, `homozygous_allele_for` precomputed |
| `9bacc58` | 4 — `InterpUnivariateMath` (1024 bins) | first approximate path; ~10 % win |
| `76e3371` | 5 — accuracy harness | 12-fixture battery, error percentiles |
| `f0beda8` | tuning | 256 bins; ~1.5 % extra from L1 pressure relief |
| `c39c78e` | 6 + 7 — bench extension + default flip | default → `InterpUnivariateMath` |
| `9834d88` | 8a + 8b — SIMD primitives | `wide` dep + lane-of-4 `ln_approx_x4` / `exp_approx_x4` |
| `952def3` | 8c + 8d — trait lane methods + new backend | `MathBackend::{ln_x4, exp_x4, HAS_LANE_4}`, `InterpUnivariateSimdMath` |
| `8ec93db` | 8e — SIMD e_step path | dispatch via `M::HAS_LANE_4`; **slower than scalar interp** |

## Empirical findings

Three results were genuine surprises and reshape what we expect from
this hardware class.

### 1. `ln`/`exp` interp wins ~10 %, not 4×

The plan estimated 1.5–2.5× on no-contam and 2–3.5× with contamination
on. The actual win was **~10 %**.

The plan's premise — "native `f64::ln` / `f64::exp` are 20–30 cycles,
interp tables are 5–10" — turned out to be wrong on a 12th-gen Intel
Core (i7-1260P). glibc's `__ieee754_log_fma` / `__ieee754_exp_fma`
*already* do the same IEEE-decomposition + range-reduction trick the
interp tables do, using FMA-friendly polynomials that fit in registers.
The interp version's table reads + bit-fiddle + branch checks eat most
of the lead. Modern libm is far better than the plan assumed.

The plan has been amended in place; the recalibrated estimate in
§"Pre-implementation speedup estimates" reflects what we now know.

### 2. Bin-index calculation dominated the per-call cost

The first cut of `exp_approx` had a wasted `fdiv` in the bin-index
arithmetic. Removing it (since `TABLE_BINS = 1024` and later 256 are
both powers of two, decompose via mask + shift) was the difference
between a 13 % **regression** and a 10 % improvement. Lesson: in
table-lookup interp, the arithmetic around the table read can cost
more than the table read itself.

### 3. Naive SIMD is slower than well-inlined scalar interp

`InterpUnivariateSimdMath` produces correct output (same parity budget
margin as the scalar interp) but is ~4 % slower than the scalar
backend. Three reasons:

- `wide` 0.7 doesn't expose `f64x4 ↔ u64x4` bit-cast or vector gather,
  so `ln_approx_x4` / `exp_approx_x4` do per-lane scalar bit
  decomposition and per-lane table reads. The arithmetic surrounding
  the gather is SIMD, but the gather isn't.
- The lane-load pattern (4 sample × 1 genotype) is strided memory
  access, not contiguous.
- The scalar interp's tight inner loop benefits from compiler
  auto-vectorisation in ways the bespoke SIMD path doesn't.

The infrastructure (`MathBackend::HAS_LANE_4` dispatch, `e_step_simd`
body, `log_sum_exp_*_x4` helpers, `EmScratch::log_post_unnorm_lane`
buffer) is in place for a future attempt that uses `std::arch`
intrinsics directly (vector gather, bit-cast) or a SIMD math library
like `sleef`.

### 4. Bench numbers across sessions are unreliable

The dev laptop runs the `powersave` cpufreq governor and the i7-1260P
swings ~2× between thermal states. A "regression" of ~70 % we
investigated turned out to be CPU state, not code. All subsequent perf
comparisons in this session used back-to-back trials inside the same
shell session.

Saved as memory: `feedback_bench_cpu_governor`.

## Decisions made

1. **Default backend: `InterpUnivariateMath` (256-bin tables).** 11 %
   faster than `ExactMath` with 50× margin under the accuracy budget.
   Plan's "≥ 1.5× speedup as sanity check" criterion was relaxed to
   "11 % on critical-path code is worth shipping" per user call.
2. **Three backends ship.** All three (`ExactMath`,
   `InterpUnivariateMath`, `InterpUnivariateSimdMath`) remain in the
   `posterior_engine::backends` submodule and reachable via
   `with_math_backend(...)`. Documented as opt-in / "may change"
   internals.
3. **3D `log_mix` table not built.** Plan deferred this as a
   follow-up; the trigger condition (mixture loop ≥ 40 % of
   `run_em_for_record` after interp lands) was not measured this
   session and the result remains deferred.
4. **256 bins, not 1024.** The accuracy harness showed huge margin at
   1024, so dropping to 256 freed L1 cache lines for ~1.5 % additional
   speedup. The plan's "Resolution is a tunable, not a number"
   principle held.

## Infrastructure delivered

- `MathBackend` trait with `ln`, `exp`, defaulted `ln_x4` / `exp_x4`,
  and `const HAS_LANE_4: bool`. Lives in
  [`posterior_engine::backends`](../../../src/var_calling/posterior_engine/backends.rs).
- `PosteriorEngine<I, M: MathBackend = InterpUnivariateMath>` — the
  engine is now generic over backend, default specified at the type
  level. `with_math_backend(...)` is the explicit-selection
  constructor.
- Genotype-shape cache in
  [`posterior_engine::shape`](../../../src/var_calling/posterior_engine/shape.rs):
  thread-local slot-array (`MAX_CACHED_PLOIDY = 8`,
  `MAX_CACHED_N_ALLELES = 16`) hands out `Arc<GenotypeShape>` carrying
  the precomputed `genotype_allele_counts`,
  `log_multinomial_coeffs`, `nonzero_pairs(_offsets)`, and
  `homozygous_allele_for`.
- IEEE-decomp interp primitives in
  [`posterior_engine::interp`](../../../src/var_calling/posterior_engine/interp.rs):
  `ln_approx`, `exp_approx` and their `f64x4` siblings. 256-bin
  uniform tables behind `LazyLock`.
- Bench harness in
  [`benches/var_calling_perf.rs`](../../../benches/var_calling_perf.rs)
  `bench_posterior_engine` group: 3 fixtures × 2 backends (Exact /
  Interp).
- Profile binary at
  [`examples/profile_posterior_engine.rs`](../../../examples/profile_posterior_engine.rs)
  with `POSTERIOR_BACKEND` env var (`exact` / `interp` / `simd`).
- Accuracy harness `tests_math_backend_accuracy` inside the
  `posterior_engine::tests` module: 12 fixtures spanning biallelic /
  triallelic / compound / chain-broken / fixation-index / triploid /
  contamination-on variants. Asserts max + p99 per metric.
- Test helpers `engine_for` / `engine_for_with_config` pin to
  `ExactMath` so the existing tight-tolerance semantic tests
  (~1e-12 against hand-computed values) keep working.

## SIMD deferral

The SIMD backend ships as a third opt-in path, not as the default. To
beat the scalar interp on this hardware we'd need:

- `std::arch::x86_64` direct intrinsics for vector gather
  (`_mm256_i64gather_pd`) and bit-cast (`_mm256_castpd_si256`).
- Equivalent NEON path for aarch64 (`vld1q_f64` + manual gather).
- Or a `sleef` / `libmvec` binding for already-tuned SIMD `ln` / `exp`.

The plan's hammer-and-glove framing still holds: rayon-over-records
is expected to be the bigger lever (near-linear with core count) and
should be tried first; SIMD comes after as a per-thread refinement.

## Validation

All run inside the dev container:

```sh
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --all-targets --all-features
./scripts/dev.sh cargo bench --bench var_calling_perf -- posterior_engine
```

737 lib tests pass. Clippy clean.

## Out-of-scope follow-ups

Still parked from the plan, in order of expected impact on this
hardware:

1. **Rayon over records** — coarse parallelism, near-linear core
   scaling. Biggest single throughput lever remaining; the
   per-record engine is already a tight, side-effect-free function
   ready to drive from `par_iter`.
2. **SIMD via direct intrinsics or sleef** — see §"SIMD deferral".
3. **3D `log_mix` table for the mixture pre-pass** — only worth
   prototyping if profiling after rayon shows the mixture loop is
   still the dominant hot spot.
4. **4D posterior LUT for standard biallelic** — discussed in the
   plan, rejected for v1; revisit if the no-contam biallelic case
   ever becomes the bottleneck.
5. **Per-record allocator amortisation** — pairs naturally with
   rayon (thread-local scratch arenas remove the per-record `Vec`
   allocations the step-0 flamegraph showed).

## File touch list

- [`src/var_calling/posterior_engine.rs`](../../../src/var_calling/posterior_engine.rs)
  — generic `PosteriorEngine<I, M>`, `EmContext` extensions
  (`log_f_per_sample`, `log_one_minus_f_per_sample`, `shape`),
  `e_step_simd`, `log_sum_exp_2_x4`, `log_sum_exp_slice_x4`,
  `tests_math_backend_accuracy`.
- [`src/var_calling/posterior_engine/backends.rs`](../../../src/var_calling/posterior_engine/backends.rs)
  (new) — `MathBackend` trait, `ExactMath`, `InterpUnivariateMath`,
  `InterpUnivariateSimdMath`.
- [`src/var_calling/posterior_engine/interp.rs`](../../../src/var_calling/posterior_engine/interp.rs)
  (new) — `ln_approx`, `exp_approx`, `ln_approx_x4`,
  `exp_approx_x4`, 256-bin tables.
- [`src/var_calling/posterior_engine/shape.rs`](../../../src/var_calling/posterior_engine/shape.rs)
  (new) — `GenotypeShape`, `shape_for`, thread-local slot-array cache,
  `log_factorial`, `log_multinomial_coefficient`, `homozygous_allele`
  (moved from posterior_engine.rs).
- [`benches/var_calling_perf.rs`](../../../benches/var_calling_perf.rs)
  — `bench_posterior_engine` group + `bench_posterior_drain` helper.
- [`examples/profile_posterior_engine.rs`](../../../examples/profile_posterior_engine.rs)
  (new) — off-Criterion drain driver with `POSTERIOR_BACKEND` env var.
- [`Cargo.toml`](../../../Cargo.toml) — added `wide = "0.7"`.
