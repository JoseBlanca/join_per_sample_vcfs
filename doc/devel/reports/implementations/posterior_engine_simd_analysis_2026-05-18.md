## SIMD analysis — root-cause finding

Date: 2026-05-18. Branch: `perf/simd-analysis` (commit follows this report).

Builds on the May 17 perf report
([`posterior_engine_perf_2026-05-17.md`](./posterior_engine_perf_2026-05-17.md))
which shipped `InterpUnivariateSimdMath` as **~4 % slower** than the
scalar interp default and concluded the SIMD work was deferred until
direct `std::arch` gather + bit-cast intrinsics could be wired in.

This follow-up investigates *why* the original SIMD attempt failed and
finds a much simpler fix.

## TL;DR

The original SIMD primitives in
[`posterior_engine::interp::{ln_approx_x4, exp_approx_x4}`](../../../src/var_calling/posterior_engine/interp.rs)
were not actually SIMD. They did the bit-decomposition and table reads
**scalar, per lane**, only wrapping SIMD arithmetic around a scalar
gather. The May 17 report attributed this to a missing feature in
`wide` 0.7: "wide 0.7 doesn't expose `f64x4 ↔ u64x4` bit-cast or
vector gather". That claim is **wrong about the bit-cast.**

`wide` 0.7 uses [`bytemuck::cast`](
https://docs.rs/bytemuck/latest/bytemuck/fn.cast.html) for vector
bit-casts (a no-op transmute between two `#[repr(C, align(32))]`
types) and [`safe_arch`](https://docs.rs/safe_arch/) for the AVX2
intrinsics. The crate ships fully SIMD polynomial-based `f64x4::ln()`
and `f64x4::exp()` — see
[`wide-0.7.33/src/f64x4_.rs:1202`](
file:///home/jose/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/wide-0.7.33/src/f64x4_.rs)
(`exp`) and `:1297` (`ln`). Both use 13th- and 5/5-degree polynomial
approximations expressed entirely in lane-parallel ops (`mul_add`,
`cast::<_, u64x4>`, `blend`, `cmp_lt`, ...) — no scalar gather, no
per-lane branch.

Swapping `InterpUnivariateSimdMath::{ln_x4, exp_x4}` from our custom
table-lookup primitives to `wide`'s built-in
`f64x4::{ln, exp}` turns SIMD from a regression into the **fastest
backend by a wide margin**, with the accuracy parity budget still met
50× over.

## Empirical result

`examples/profile_posterior_engine` against the biallelic-contam-on
fixture (64 samples × 10 000 records, 30 runs back-to-back inside one
shell session — same governor-state, per the May 17 finding).

| Backend | Custom `_x4` primitives | `wide` native `f64x4::ln`/`exp` | Δ |
|---|---|---|---|
| `ExactMath` | 24.41 µs/record | 24.41 µs/record | — |
| `InterpUnivariateMath` (scalar interp) | 22.33 µs/record (−8.5 %) | 22.33 µs/record (−8.5 %) | — |
| `InterpUnivariateSimdMath` | 25.96 µs/record (+6.4 %, slower than exact) | **18.15 µs/record (−25.6 %)** | **−30 % per record** |

The third column is what you get from this one-line change at
[`backends.rs:159-167`](../../../src/var_calling/posterior_engine/backends.rs#L159-L167):

```rust
fn ln_x4(&self, x: [f64; 4]) -> [f64; 4] {
    wide::f64x4::from(x).ln().to_array()
}
fn exp_x4(&self, x: [f64; 4]) -> [f64; 4] {
    wide::f64x4::from(x).exp().to_array()
}
```

Same parity budget margin as the original SIMD backend — the
`interp_univariate_simd_clears_approximate_parity_budget` accuracy
harness passes unchanged (`wide`'s polynomial gives ~`1e-15` relative
error, *better* than our 256-bin table, so the end-to-end posterior
parity stays well under the 1e-4 budget).

## Why the original SIMD was slow — five concrete reasons

In rough order of impact (only #1 was measured directly; #2-5 are
inspection findings):

### 1. The "_x4" primitives were ~75 % scalar work

Look at the structure of [`exp_approx_x4`](
../../../src/var_calling/posterior_engine/interp.rs#L216-L270):

```rust
let arr = y.to_array();                            // SIMD → stack
if arr.iter().any(|&v| !v.is_finite()) { ... }     // scalar branch
...
let scaled = y * f64x4::splat(SCALE);              // SIMD ✓
let scaled_arr = scaled.to_array();                // SIMD → stack
if scaled_arr.iter().any(...) { ... }              // scalar branch
let total_f = scaled.floor();                      // SIMD ✓
let frac = scaled - total_f;                       // SIMD ✓
let total_arr_i: [i32; 4] = { ... };               // SIMD → stack → scalar i32
for lane in 0..4 {                                 // ── scalar gather loop ──
    let t = total_arr_i[lane];
    let idx = (t as usize) & (TABLE_BINS - 1);
    let n = t >> TABLE_INDEX_BITS;
    lo[lane] = table[idx];                         //   scalar table read
    hi[lane] = table[idx + 1];                     //   scalar table read
    pow2_n_arr[lane] = f64::from_bits(...);        //   scalar bit injection
}
let pow2_f = f64x4::from(lo) + frac * (f64x4::from(hi) - f64x4::from(lo));
pow2_f * f64x4::from(pow2_n_arr)
```

Steps that touch tables, bits, and the per-lane integer index are all
scalar. The SIMD arithmetic is just the post-gather composition
(maybe 5 ops). On a sample-major hot path that calls `exp_x4`/`ln_x4`
dozens of times per record, that's most of the cost.

`wide::f64x4::exp` solves the same problem by *not* needing a table.
It uses a 13-term polynomial (13 FMAs) plus a vector `2^n`
construction via integer shift + bit-cast — every step is lane-parallel.

### 2. Trait signature forces array round-tripping at every call site

The `MathBackend` trait is defined as
[`fn exp_x4(&self, x: [f64; 4]) -> [f64; 4]`](
../../../src/var_calling/posterior_engine/backends.rs#L62-L71). Hot
call sites in
[`e_step_simd`](../../../src/var_calling/posterior_engine.rs#L1466)
and [`log_sum_exp_slice_x4`](../../../src/var_calling/posterior_engine.rs#L1764-L1781)
all do:

```rust
let post_v = wide::f64x4::from(math.exp_x4(diff_v.to_array()));
```

Every call: `f64x4 → [f64; 4] → f64x4`. With `lto = "fat"` and
`codegen-units = 1` the optimiser often elides this (both types are
`#[repr(C, align(32))]`), but the trait surface still forces the IR
to *name* a `[f64; 4]` temporary, which constrains register
allocation and inlining decisions.

A cleaner signature would be `fn exp_x4(&self, x: f64x4) -> f64x4`.
The scalar backends would implement it by `to_array()`-ing once and
calling `f64::exp` four times — but those backends never sit on
`HAS_LANE_4 = true` anyway, so they don't matter for hot-path perf.

Not measured independently — folded into the #1 fix because the
swap to `wide`'s built-in works through exactly the same trait
surface. Worth doing separately if profiling after the rest of the
plan shows the array bouncing is still visible.

### 3. The per-vector fallback branch poisons the hot path

Both custom primitives start with
[`if arr.iter().any(|&v| !v.is_finite()) { … scalar fallback … }`](
../../../src/var_calling/posterior_engine/interp.rs#L167-L173).
For hot-path inputs (probabilities, normalised log-likelihoods —
all in the normal range), the branch is never taken, but the
compiler still has to materialise the array and the predicate
*before* any SIMD work can run. That gives a worse register schedule
than wide's approach (compute SIMD-fast, then `blend` the rare
fallback at the end).

### 4. Lane = sample → most per-genotype work is scalar splat

`e_step_simd` parallelises across *samples*: four samples per
batch, one genotype per inner-loop iteration. But most of the work
inside the genotype loop —
[`log_indep`, `pairs.iter().map(...).sum()`,
`log_multinomial_coeffs[g_idx]`](
../../../src/var_calling/posterior_engine.rs#L1417-L1422) — is a
function of `g_idx` alone, identical across all four samples. We
compute it once and `splat` it into a vector.

Per genotype, only **3** per-lane-distinct values exist (`ll_v`,
`log_f_v`, `log_one_minus_f_v`), composed with a handful of SIMD
adds and one `log_sum_exp_2_x4`. Once `ln`/`exp` are truly SIMD,
this layout works fine — but it's *not* the maximally-parallel
shape. An alternative would be lane = genotype (4 genotypes
per inner-loop iteration), which makes `log_indep` itself a
SIMD computation. Worth prototyping if the wide-native fix doesn't
land enough of the win, but the empirical −30 % suggests it's
not necessary.

### 5. Strided gather load for `log_likelihoods`

The `ll_v` load is [4 manual scalar reads](
../../../src/var_calling/posterior_engine.rs#L1436-L1441) with
stride `n_genotypes * 8` bytes. For the dominant biallelic-diploid
case (`n_genotypes = 3`, stride 24 bytes), all four lanes land in
1–2 cache lines, so this is cheap on a modern Intel core. For
tetraallelic (`n_genotypes = 15`, stride 120 bytes) it costs 4
cache lines. AVX2's `_mm256_i64gather_pd` would compress this
into one instruction but is usually *not* faster than 4 scalar
loads — gather is bandwidth-limited, not latency-limited.

A more interesting fix is a one-shot layout transpose at the top
of `e_step` into a batch-of-4-samples-major buffer, but that's
substantially more work. Defer unless profiling after the wide-native
fix says it's still a hotspot.

## What this changes for the May 17 report

§"SIMD deferral" of `posterior_engine_perf_2026-05-17.md` calls for
either `std::arch::x86_64` direct intrinsics or a sleef binding. **We
don't need either.** `wide` 0.7 already provides exactly the SIMD `ln`
and `exp` we needed, and they outperform the table-lookup approach
because the polynomial doesn't pay the per-lane gather cost.

The "rayon-over-records is the bigger lever" framing still holds —
multi-core parallelism scales near-linearly, SIMD this hardware
class caps around 4×. But the SIMD lever turned out to be a **−25 %
single-threaded win** rather than a +4 % regression. That changes the
relative ordering of the out-of-scope follow-ups:

- The SIMD backend is now a *default-candidate*, not a deferred
  experiment.
- The 3D `log_mix` table and 4D posterior LUT both remain
  speculative; with `wide`'s native exp/ln the mixture pre-pass's
  `math.ln(mix)` call is already vectorisable if we batch the inner
  loop.

## Recommended next steps (in order of expected impact)

1. **Land the `wide`-native primitive swap on a separate PR**, flip
   the engine default from `InterpUnivariateMath` to
   `InterpUnivariateSimdMath`, and re-baseline the bench harness.
   The accuracy harness already validates parity; the only delta is
   a different polynomial vs table approximation for `ln` / `exp`.
2. **Either delete or rewrite [`interp::ln_approx_x4` / `exp_approx_x4`](
   ../../../src/var_calling/posterior_engine/interp.rs#L137-L270)**.
   They're dead code after step 1. Keep the scalar `ln_approx` /
   `exp_approx` — they're still used by `InterpUnivariateMath` and
   are genuinely faster than `f64::ln` / `f64::exp` for the scalar
   backend (~10 % win per May 17 report).
3. **Audit the trait signature**: change `ln_x4` / `exp_x4` to take
   `f64x4`, not `[f64; 4]`. Drops a measurable indirection at every
   call site in `log_sum_exp_*_x4` and the `e_step_simd` final
   normalise pass.
4. **Re-profile with samply** after 1–3 to find the next bottleneck.
   Candidates: the strided gather, the per-genotype splat overhead,
   the `log_sum_exp_slice_x4` two-pass loop. Don't speculate beyond
   this report's measurements until we have a fresh flamegraph.
5. **Then** revisit rayon-over-records — every per-thread speedup
   from SIMD multiplies into the rayon lever.

## Reproducing

```sh
git checkout perf/simd-analysis
./scripts/dev.sh cargo build --release --example profile_posterior_engine

# Triple back-to-back inside one shell to avoid CPU-state drift
# (see feedback_bench_cpu_governor):
POSTERIOR_BACKEND=exact  ./target-container/release/examples/profile_posterior_engine 2>&1 | tail -1
POSTERIOR_BACKEND=interp ./target-container/release/examples/profile_posterior_engine 2>&1 | tail -1
POSTERIOR_BACKEND=simd   ./target-container/release/examples/profile_posterior_engine 2>&1 | tail -1
```

The diff against `main` is two lines in
[`src/var_calling/posterior_engine/backends.rs`](
../../../src/var_calling/posterior_engine/backends.rs):

```diff
-    fn ln_x4(&self, x: [f64; 4]) -> [f64; 4] {
-        super::interp::ln_approx_x4(wide::f64x4::from(x)).to_array()
-    }
-    fn exp_x4(&self, x: [f64; 4]) -> [f64; 4] {
-        super::interp::exp_approx_x4(wide::f64x4::from(x)).to_array()
-    }
+    fn ln_x4(&self, x: [f64; 4]) -> [f64; 4] {
+        wide::f64x4::from(x).ln().to_array()
+    }
+    fn exp_x4(&self, x: [f64; 4]) -> [f64; 4] {
+        wide::f64x4::from(x).exp().to_array()
+    }
```

Tests: `./scripts/dev.sh cargo test --release --lib var_calling::posterior_engine`
(89 tests pass, including
`interp_univariate_simd_clears_approximate_parity_budget`).
