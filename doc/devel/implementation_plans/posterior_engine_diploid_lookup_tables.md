# Posterior engine: interpolated lookup tables (hybrid)

Proposal date: 2026-05-17. Revised same day (twice).

> **Status: draft for discussion.** Open questions are flagged inline
> (`**Open question**`) and collected at the bottom under
> "Decisions to confirm before implementation".

> **Revision note.** Earlier drafts of this plan proposed
> *structural* specialisation of the diploid case (precomputed
> discrete-input tables, restructured inner loops). After
> discussion we pivoted to a different strategy: replace expensive
> `ln`/`exp` calls with **univariate interpolated lookup tables**
> over continuous arguments. A multivariate 3D `log_mix` table
> was considered for v1 and is now deferred to a follow-up
> revisited only if the univariate-only approach falls short. The
> rationale is in §"Why interpolated tables, not structural
> specialisation".

## Domain intent

The Stage 6 EM loop in
[src/var_calling/posterior_engine.rs](../../src/var_calling/posterior_engine.rs)
spends most of its CPU time in a handful of `ln` and `exp` calls
inside three nested hot paths:

- the mixture-likelihood pre-pass at
  [posterior_engine.rs:765-799](../../src/var_calling/posterior_engine.rs#L765-L799)
  — one `mix.ln()` per `(sample, genotype, allele)` cell;
- the E-step at
  [posterior_engine.rs:1217-1244](../../src/var_calling/posterior_engine.rs#L1217-L1244)
  — `safe_ln`, `log_sum_exp_2`, `log_sum_exp_slice`, plus one
  normalising `exp` per cell;
- `log_p_effective` rebuild per E-step iteration — one `safe_ln`
  per allele.

Native `f64::ln` and `f64::exp` cost ~20–30 cycles each on x86-64.
A linear-interpolation lookup table over a small grid is ~5–10
cycles per call and carries a *controlled* relative error
(`~1e-5` with a 1024-entry table, `~1e-7` with 4096 entries) —
well inside the 1 Phred / `1e-4`-posterior approximate-parity
target we've set.

The hybrid this plan proposes:

1. **Univariate interpolated `ln_approx` / `exp_approx`** — drop-in
   replacements for `f64::ln` / `f64::exp` at every hot call site.
   This is the bulk of the speedup and is independent of ploidy.
2. **Two preparatory hoists** (a `safe_ln(f_s)` move and a tiny
   `(ploidy, n_alleles)`-keyed genotype-shape cache) — small,
   independent wins worth doing alongside.

The diploid-only specialisation from the previous draft is
dropped. The interpolation strategy works equally well at any
ploidy, so there's no reason to restrict it. A 3D `log_mix`
table that bakes the `(1-c_s)*p_own + c_s*q_b` mixture
arithmetic into a single trilinear lookup was also considered;
the bench-vs-arithmetic trade-off is unclear a priori (8 cache
reads vs ~3 multiplies and a 1D `ln` lookup), so we **start with
univariate-only** and revisit the 3D table only if the bench
after item 1 lands shows the mixture loop is still the
bottleneck. See §"Out-of-scope follow-ups / 3D `log_mix` table".

## Why now

Same as the previous draft: Stage 6 is the bottleneck on the
critical path, the EM inner loop is the hottest function, and
the contamination side-pass adds a second `ln`-heavy hot loop.
Native `ln`/`exp` show up as the top self-time entries in any
profile of `run_em_for_record`.

## Why interpolated tables, not structural specialisation

The previous draft proposed precomputing diploid-specific
genotype tables and restructuring the E-step inner loop. After
discussion:

- **Bigger reachable speedup.** The structural rewrite collapses
  branches and saves a few multiplies; the per-call `ln` / `exp`
  cost is unchanged. Interp tables attack the cycles those calls
  actually spend.
- **Ploidy-agnostic.** No diploid/polyploid dispatch needed.
- **Smaller surface change.** A few drop-in helper functions; no
  parallel `e_step_diploid`, no `run_em_diploid`, no shape cache
  redesign.
- **Composes with structural work later.** If a future pass wants
  to specialise diploid as well, the interp tables stay in place
  underneath.
- **What's lost.** The structural rewrite gave bit-identical or
  near-bit-identical output. Interp tables introduce a real (but
  bounded) approximation. We've already decided the 1 Phred /
  `1e-4` tolerance is acceptable, so this isn't a regression — but
  it is a meaningful change in mindset.

## Why not a multidimensional posterior LUT

The full-ambition idea — table the per-sample posterior as a
function of `(eps_REF, eps_ALT, depth, balance, c_s, q_b, p_hat,
F)` and read out the answer with multidim linear interp —
doesn't work in practice. Sketch of the analysis:

- **7–8 axes** for the biallelic-with-contamination case. Even at
  10 bins per axis the table is ~100 MB; at 20 bins it's
  multi-GB.
- **Linear interp cost grows as `2^D`.** A 7D lookup reads 128
  corner cells. If the table is bigger than L2, each lookup is
  cache-miss-bound and likely slower than native `ln`/`exp`.
- **`p_hat` is iterated by EM.** It's not a static input; the
  table would still be queried inside the loop.
- **Per-`n_alleles` tables** — dimensionality grows with allele
  count.

A 4D table for the standard biallelic case (~4 MB, fits in L2)
*is* feasible and could be revisited as a follow-up, but it ducks
the contamination case — which is where most of the speedup is
wanted. The hybrid plan below captures most of the realistic win
without those constraints.

## Approximate parity

Unchanged from the previous draft. Targets:

- **Posteriors**: `1e-4` absolute everywhere.
- **`gq_phred` / `qual_phred`**: `1.0` Phred absolute.
- **`allele_frequencies` / `compound_frequencies`**: `1e-4`
  absolute.
- **`best_genotype`**: must match when the top posterior margin is
  ≥ `0.05`; skipped on borderline fixtures.
- **`diagnostics.iterations` / `final_max_delta_p`**: not
  compared.

The posterior is orientative — produced by an approximate model
— so a few-ULP approximation error per `ln`/`exp` call, summed
across the EM trajectory, is comfortably below the model's own
uncertainty.

## Out of scope

- **Structural diploid specialisation.** Dropped per the
  rationale above.
- **Multidim posterior LUT** (4D and up). Listed as a
  follow-up.
- **Algorithmic changes.** No change to the EM update rules, the
  HWE-with-`F` prior, the convergence criterion, the error
  taxonomy, or the `PosteriorRecord` shape.
- **`MathBackend` becoming a public API surface.** The pluggable
  backend (§"Pluggable math backend") is an internal mechanism
  for A/B testing and for runtime-selectable behaviour, not a
  user-facing knob initially.
- **Parallelism.** Per-record EM stays single-threaded.

## The interpolated tables

### A. Univariate `ln_approx` and `exp_approx`

Two new helpers in a new submodule `posterior_engine::interp`:

```rust
// src/var_calling/posterior_engine/interp.rs

/// Approximate `f64::ln` via linear interpolation over a uniform
/// grid in `[X_MIN, X_MAX]`. Inputs outside the range fall back
/// to the native implementation.
pub(crate) fn ln_approx(x: f64) -> f64 { … }

/// Approximate `f64::exp` analogously over `[Y_MIN, Y_MAX]`.
pub(crate) fn exp_approx(y: f64) -> f64 { … }
```

**Domain selection.** For `ln_approx`, the realistic input range
is determined by where it's called from:

- `safe_ln(p)` for `p ∈ (0, 1]` — needs the small-`x` tail.
- `mix.ln()` for `mix ∈ (0, 1]` — same.
- `log_sum_exp_*` internals reduce to `ln(1 + ε)` for small `ε`,
  but the engine evaluates the natural form with `m + (...).ln()`
  where the argument is in `[1, k]` for `k = n_genotypes` —
  needs the modest-`x` regime.

Cleanest implementation: exploit the IEEE 754 decomposition.
`f64::ln(x) = e * LN_2 + ln(m)` where `m ∈ [1, 2)` is the
mantissa and `e` the (unbiased) exponent. We only need to table
`ln(m)` over `[1, 2)`. The same trick works for `exp_approx`:
`exp(y) = 2^(y * LOG2_E)`; split into integer + fractional
parts, look up `2^f` for `f ∈ [0, 1)`.

**Why uniform spacing inside the mantissa sub-table is fine.**
The theoretical optimum for linear interpolation is to space
input points so that *outputs* are equally spaced — equivalently,
spacing inversely proportional to `sqrt(|f''(x)|)`. Lazy uniform
spacing over a wide `x` range would be wasteful because curvature
varies wildly. But the IEEE decomposition itself acts as the
adaptive transform: inside `m ∈ [1, 2)`, `f''(m) = -1/m²` ranges
only from `-1` to `-0.25` (a 4× variation), and for `2^f` over
`f ∈ [0, 1)` curvature varies by 2×. Uniform spacing on the
sub-domain already gives near-uniform output error, with no
extra branch in the lookup. Non-uniform spacing would add a
search-or-branch cost that swamps the per-call interp time.

**Resolution is a tunable, not a number in this plan.** The
right number of entries is whatever the smallest table is that
keeps the **error distribution** (not just the max) below the
parity budget on real fixtures. The process is: start at e.g.
1024 entries per sub-table, run the accuracy harness on the full
fixture battery + the real-data spot check, and double until the
max absolute output error on posteriors is comfortably below
`1e-4` and the Phred error is below `1.0`. With 1024 entries we
expect ~`1e-6` relative error per `ln` call; 4096 gives ~`1e-7`.
Tables stay small enough to fit in L1 in either case (8 KB and
32 KB respectively).

**Out-of-range handling.** For inputs outside the tabled domain
(very small `x`, very large negative `y`), fall back to native
`f64::ln` / `f64::exp`. These cases are rare but must remain
correct.

**Edge cases.**

- `ln_approx(0.0)` → `f64::NEG_INFINITY` (matches `safe_ln`'s
  contract).
- `ln_approx(x < 0)` → `f64::NEG_INFINITY` (matches `safe_ln`).
  In practice this is unreachable in the engine because the
  upstream code guards probabilities to be non-negative.
- `exp_approx(f64::NEG_INFINITY)` → `0.0`.
- `exp_approx(f64::INFINITY)` → `f64::INFINITY`.
- `NaN` inputs propagate through native fallback.

### B. (Deferred to follow-up) — 3D `log_mix_approx`

A 3D table indexed by `(p_own, c_s, q_b)` that returns
`ln((1-c_s)*p_own + c_s*q_b)` directly was considered for v1 and
is **not pursued**. The reasoning is captured in §"Out-of-scope
follow-ups / 3D `log_mix` table", together with the criteria for
revisiting it if `InterpUnivariateMath` doesn't deliver enough
speedup.

### C. (Not pursued in v1) — additional ≤3D tables

Other intermediates considered and rejected:

- `log_prior(k, p, F)` — the IBD prior contribution. After
  `safe_ln(f_s)` is hoisted (item D below) and `log_p_effective`
  is cached per E-step iteration, the prior cost is one
  `log_sum_exp_2` per homozygous genotype. Replacing the inner
  `exp`/`ln` of `log_sum_exp_2` via the univariate approx
  already captures the win; a dedicated table doesn't add much.
- `log_factorial(n)` — `n ∈ [0, ploidy]` is small; the loop is
  unrolled by the compiler. Not worth a table.

These can be revisited as follow-ups if profiling identifies them
later.

## Preparatory hoists (kept from previous draft)

### D. Hoist `safe_ln(f_s)` and `safe_ln(1 - f_s)` to `EmContext`

Today these are computed per sample inside `e_step` at
[posterior_engine.rs:1214-1215](../../src/var_calling/posterior_engine.rs#L1214-L1215),
every iteration. They're functions of `fixation_indices`, which
is record-static. Move to record setup:

```rust
struct EmContext<'a> {
    // … existing fields …
    log_f_per_sample: &'a [f64],            // safe_ln(f_s)
    log_one_minus_f_per_sample: &'a [f64],  // safe_ln(1 - f_s)
}
```

Bit-identical to today's output. Trivial to validate.

### E. `(ploidy, n_alleles)`-keyed genotype-shape cache + e_step restructure

The step-0 flamegraph shows `e_step` itself at ~30 % of self-time
on top of the ~40 % spent in `ln`/`exp`. Item A targets the
math; this item targets the EM inner-loop arithmetic by
**extending `GenotypeShape`** with two precomputed fields that
let `e_step` skip three avoidable costs per cell:

1. **The `if k == 0` skip-branch** in the per-allele sum.
2. **The `homozygous_allele(gt_counts)` linear scan** done per
   cell to pick the IBD prior branch.
3. **The cast-and-multiply by `k as f64`** when `k = 1` (a
   wasted multiply by 1.0).

`run_em_for_record` today rebuilds three pure-function-of-`(ploidy,
n_alleles)` artefacts on every record:

- `genotype_order(ploidy, n_alleles)` ([per_group_merger.rs:355](../../src/var_calling/per_group_merger.rs#L355));
- the flat `genotype_allele_counts` table;
- `log_multinomial_coeffs`.

We add two more, both also pure functions of `(ploidy, n_alleles)`:

- `nonzero_pairs` — flat list of `(allele_idx, count)` pairs with
  `count > 0`, concatenated per genotype; for diploid each
  genotype has 1 or 2 entries.
- `homozygous_allele_for[g]` — `Some(a)` if genotype `g` is
  homozygous for allele `a`, `None` otherwise.

Adjacent biallelic-SNP records share the same `(2, 2)` shape;
caching is a free win. Implementation — a fixed-size slot array
indexed by `(ploidy, n_alleles)`:

```rust
struct GenotypeShape {
    n_genotypes: usize,
    genotype_allele_counts: Vec<u32>,
    log_multinomial_coeffs: Vec<f64>,
    /// Flat list of `(allele_idx, count)` pairs with count > 0,
    /// concatenated per genotype. For diploid: 1 entry (homozygous)
    /// or 2 entries (heterozygous) per genotype.
    nonzero_pairs: Vec<(u8, u32)>,
    /// `(start, len)` into `nonzero_pairs`, per genotype.
    nonzero_pairs_offsets: Vec<(u32, u8)>,
    /// Per genotype: `Some(a)` if homozygous for allele `a`,
    /// `None` otherwise. Replaces the per-cell
    /// `homozygous_allele(gt_counts)` scan in `e_step`.
    homozygous_allele_for: Vec<Option<u8>>,
}

const MAX_CACHED_PLOIDY: usize = 8;       // covers diploid through octoploid
const MAX_CACHED_N_ALLELES: usize = 16;   // generous headroom above Stage 5's max_alleles
const CACHE_SLOTS: usize = (MAX_CACHED_PLOIDY + 1) * (MAX_CACHED_N_ALLELES + 1);

thread_local! {
    static SHAPE_CACHE: RefCell<[Option<Arc<GenotypeShape>>; CACHE_SLOTS]>
        = RefCell::new(std::array::from_fn(|_| None));
}

fn shape_for(ploidy: u8, n_alleles: usize) -> Arc<GenotypeShape> {
    let p = ploidy as usize;
    if p > MAX_CACHED_PLOIDY || n_alleles > MAX_CACHED_N_ALLELES {
        // Out-of-bounds shapes are rare; build without caching to
        // avoid widening the slot array. No correctness risk —
        // GenotypeShape construction is pure.
        return Arc::new(build_genotype_shape(ploidy, n_alleles));
    }
    let idx = p * (MAX_CACHED_N_ALLELES + 1) + n_alleles;
    SHAPE_CACHE.with(|cell| {
        let mut slots = cell.borrow_mut();
        if let Some(shape) = &slots[idx] {
            return Arc::clone(shape);
        }
        let shape = Arc::new(build_genotype_shape(ploidy, n_alleles));
        slots[idx] = Some(Arc::clone(&shape));
        shape
    })
}
```

No eviction, no LRU bookkeeping; the slot array is small enough
(~144 `Option<Arc<...>>` entries per thread) that all realistic
shapes fit. Out-of-bounds shapes fall back to the uncached path —
no correctness risk because `GenotypeShape` construction is
pure. Bit-identical to today's output.

If a future profile points the bottleneck at the cache lookup
itself (it won't — the lookup is once per record, dwarfed by
the EM loop), revisit with a HashMap+LRU. The slot array is
the simpler shape and matches the (heavily concentrated)
real-world `(ploidy, n_alleles)` distribution.

**The matching `e_step` rewrite.** With the cache populated, the
hot inner loop is rewritten to consume `nonzero_pairs` and
`homozygous_allele_for` directly:

```rust
for g in 0..shape.n_genotypes {
    let (start, len) = shape.nonzero_pairs_offsets[g];
    let pairs = &shape.nonzero_pairs[start as usize..start as usize + len as usize];

    // log_indep — sum is at most `ploidy` terms (1 or 2 for diploid),
    // and the `k == 0` branch is gone because zero entries aren't listed.
    let log_indep = shape.log_multinomial_coeffs[g]
        + pairs
            .iter()
            .map(|&(a, k)| f64::from(k) * scratch.log_p_effective[a as usize])
            .sum::<f64>();

    // log_prior — the IBD/non-IBD branch reads from a precomputed
    // flag instead of re-scanning the count row.
    let log_prior = match shape.homozygous_allele_for[g] {
        Some(a) => log_sum_exp_2(
            log_one_minus_f + log_indep,
            log_f + scratch.log_p_effective[a as usize],
        ),
        None => log_one_minus_f + log_indep,
    };

    scratch.log_post_unnorm[g] = ll_row[g] + log_prior;
}
```

**Numerical parity.** Bit-identical to today's output. Skipping
`k == 0` entries doesn't change float arithmetic
(`0.0 * log_p + y == y` exactly under IEEE 754). The
`homozygous_allele_for` lookup substitutes for the scan but
returns the same value. No re-association of the inner sum.

**Ploidy-agnostic.** Works at any ploidy — polyploid genotypes
simply have longer `nonzero_pairs` entries. The deeper
diploid-only specialisation (a `match` over a
`DiploidGenotypeCounts` enum that hard-codes `k_a ∈ {0, 1, 2}`)
remains available as a follow-up if profiling after this lands
still points at `e_step`.

## Pluggable math backend

To compare different interp strategies (and to keep the exact
implementation around for parity testing and as a safety
fallback), introduce a `MathBackend` trait monomorphised into the
engine:

```rust
// src/var_calling/posterior_engine/backends.rs (re-exported as
// `posterior_engine::backends`)

/// Math backend trait. Public surface because it appears in
/// `PosteriorEngine`'s type parameter bound.
pub trait MathBackend: Sync {
    fn ln(&self, x: f64) -> f64;
    fn exp(&self, x: f64) -> f64;
}

/// Bit-identical baseline. Use this when reproducibility against
/// the unoptimised engine matters.
pub struct ExactMath;
impl MathBackend for ExactMath {
    fn ln(&self, x: f64) -> f64 { x.ln() }
    fn exp(&self, x: f64) -> f64 { x.exp() }
}

/// Interpolated `ln` / `exp` via static lookup tables.
///
/// Public for benchmarking and integration testing. The internals
/// (table layout, resolution, the precise set of backends in this
/// module) may change without a major version bump; for stable
/// output across versions, construct the engine via
/// `PosteriorEngine::new()` / `::with_config()` and don't name
/// the backend explicitly.
pub struct InterpUnivariateMath { /* static tables */ }
impl MathBackend for InterpUnivariateMath { … }
```

Both impl types live in a `pub mod backends` submodule. The
segregation keeps `PosteriorEngine`'s top-level public surface
focused on the engine itself, while the trait + types are
reachable from benches and integration tests as
`posterior_engine::backends::{MathBackend, ExactMath,
InterpUnivariateMath}`. The stability note on
`InterpUnivariateMath` signals that callers who want stable
output across versions should use the default constructors.

The trait surface intentionally stops at `ln` and `exp`. The
mixture loop calls `backend.ln(mix)` directly with the
arithmetically computed `mix`. If a future revisit adds a 3D
`log_mix` table (see follow-ups), the trait gains a defaulted
`log_mix(p_own, c_s, q_b)` method and the mixture call site
switches to that — a strictly additive change.

`PosteriorEngine` becomes generic over the backend:

```rust
pub struct PosteriorEngine<I, M: MathBackend = ExactMath> { … }

impl<I> PosteriorEngine<I, ExactMath> {
    pub fn new(upstream: I) -> Self;            // default = exact
    pub fn with_config(upstream: I, config: PosteriorEngineConfig) -> Self;
}

impl<I, M: MathBackend> PosteriorEngine<I, M> {
    pub fn with_math_backend(upstream: I, config: PosteriorEngineConfig, math: M) -> Self;
}
```

Monomorphisation gives zero-overhead static dispatch; the trait
methods inline. `Default` keeps the existing call sites
compiling.

**Open question:** should the production default flip to an
interp backend once the bench is run, or stay on `ExactMath` and
require opt-in? Proposal: flip after the accuracy harness passes
on a representative real-data sample, with `ExactMath` available
as an escape hatch.

**Open question:** is `MathBackend: Sync` enough, or do we need
`Send` for future parallelism? `Sync` is the read-only-data
contract; the tables are static so `Sync` suffices today.

## Test strategy

### Accuracy harness

A new test module `tests_math_backend_accuracy` runs every
backend over a fixture battery and asserts the approximate-parity
targets against `ExactMath`:

```rust
fn run_accuracy_battery<M: MathBackend>(math: M) -> AccuracyReport { … }

struct AccuracyReport {
    // Per-cell error distribution across all (fixture, sample,
    // genotype) cells — not just the max. With non-uniform output
    // error from uniform-input tables, the max can come from one
    // pathological corner; the percentiles tell us whether that
    // corner is representative of real data.
    posterior_error: ErrorDistribution,
    phred_error: ErrorDistribution,
    allele_freq_error: ErrorDistribution,
    fixtures_total: usize,
    fixtures_with_argmax_mismatch_above_margin: usize,  // must be 0
    fixtures_borderline_skipped: usize,
}

struct ErrorDistribution {
    max: f64,
    p99: f64,
    p95: f64,
    p50: f64,
    mean: f64,
    /// Where the worst case occurred — useful for diagnosing which
    /// fixture / region of input space dominates the error.
    argmax_fixture: String,
}
```

The fixture battery covers:

- biallelic SNP, single sample, strong-REF and strong-ALT;
- biallelic SNP, 100 samples, mixed evidence;
- triallelic SNP;
- compound allele (chain-anchored);
- chain-broken sample on a compound;
- record with non-zero `fixation_index_default`;
- record with per-sample fixation index overrides;
- **triploid** (exercises the polyploid path — the interp
  backend must work here unchanged);
- **contamination on**, all samples `c_s > 0`;
- **contamination on**, mixed `c_s = 0` and `c_s > 0`;
- **contamination on**, triallelic site;
- **contamination on**, compound allele.

For each backend, the test asserts on the **distribution**:

- `posterior_error.max ≤ 1e-4` **and** `posterior_error.p99 ≤ 5e-5`
  — the p99 cap catches a backend that meets the max but has a
  fat tail;
- `phred_error.max ≤ 1.0` **and** `phred_error.p99 ≤ 0.5`;
- `allele_freq_error.max ≤ 1e-4`;
- `fixtures_with_argmax_mismatch_above_margin == 0`.

The full report (all percentiles, all error families, plus
`argmax_fixture` for each) is printed at INFO level and pasted
into the PR description. This is the artefact that drives the
resolution-tuning loop: if `posterior_error.p99` is well below
the budget but `posterior_error.max` is close to it, the worst
case is a corner — look at `argmax_fixture` and decide whether
that corner is representative; if it's a synthetic pathology,
keep the resolution; if it's a real-data pattern, bump the
resolution.

### Production-data accuracy spot check

A separate `--ignored` test (manual run) takes a small real-data
sample (~1 000 records from a representative cohort) and runs
both `ExactMath` and the candidate backends, reporting the same
metrics. This catches distribution-shift issues the synthetic
fixtures miss (e.g. real `c_s` distributions, real `q_b` shapes).

**Fixture location.** While this is experimental work, the real
data lives in the project-local `tmp/` (git-ignored per
[CLAUDE.md §"Scratch space"](../../CLAUDE.md)). The
`--ignored` test reads from a hard-coded `tmp/...` path and
skips with a clear message if the file isn't present. Once the
experiment converges on a backend choice and the spot check
becomes something we want to be reproducible across machines, we
decide separately how to handle it (committed `.psp`-style
frozen fixture, downloader script, etc.) — out of scope for this
plan.

### Existing tests

All ~50 existing tests in
[posterior_engine.rs](../../src/var_calling/posterior_engine.rs#L1538)
continue to use `ExactMath` (via the `Default` parameter) and
must pass unchanged. The interp backends are tested via the
accuracy harness, not by re-running every existing test.

### Cache test

A test for the genotype-shape cache (item E) asserts repeated
lookups within bounds return `Arc`-equal pointers (cache hit),
and that out-of-bounds shapes (e.g. `n_alleles >
MAX_CACHED_N_ALLELES`) return a freshly built shape each time
(uncached fallback). The slot-array layout means there's no LRU
behaviour to test.

## Bench strategy

The bench harness is built **first** (step 0 of the
implementation order), before any engine changes, so we have a
real baseline to measure subsequent steps against — and so we
catch issues in the bench infrastructure itself early. A new
bench file `benches/posterior_engine_perf.rs` (or a new group
inside the existing `benches/var_calling_perf.rs`) runs the same
engine code across every available backend. Initial shape at
step 0:

```rust
fn bench_posterior_engine(c: &mut Criterion) {
    for (name, records) in fixtures() {
        let mut group = c.benchmark_group(name);
        group.bench_function("baseline", |b| { run_today(b, &records) });
    }
}
```

`run_today` calls `run_em_for_record` as it exists pre-changes
(no `MathBackend` trait, no caching). After step 1 lands, this
becomes `bench_function("exact", run::<ExactMath>)`; after step 6,
a `bench_function("interp_univariate", run::<InterpUnivariateMath>)`
line is added. A third `interp_log_mix_3d` line can be added
later if the follow-up 3D table is implemented; the bench
harness shape doesn't change.

**One-off flamegraph at step 0.** Run `samply` (or
`cargo flamegraph`) against the heaviest fixture
(biallelic-contam-on, 64 samples, 10 000 records) and capture
the resulting profile in the step-0 PR description. Confirm
`ln` / `exp` / mixture-loop code shows up at the top of self-time;
if it doesn't, the plan's premise is wrong and we revise
before building the trait. Subsequent steps don't require a
new flamegraph unless a bench result surprises us.

Fixtures (each generating a stream of `MergedRecord`s of the
appropriate shape):

- **biallelic SNP, 64 samples, 10 000 records, no contamination**
  — the dominant production case.
- **biallelic SNP, 64 samples, 10 000 records, contamination on**
  — exercises the mixture pre-pass.
- **triallelic SNP, 64 samples, 1 000 records, contamination on**
  — wider mixture loop.
- **6-allele site, 64 samples, 1 000 records** — wide EM, no
  contamination.
- **biallelic SNP, 64 samples, 10 000 records, `ploidy = 3`** —
  polyploid control.

Each fixture is run against every backend; Criterion's
`benchmark_group` gives a side-by-side comparison.

**Reporting.** A short script (or just a copy-pasted criterion
table in the PR description) summarises:

- absolute wall-time per (fixture, backend) — with the step-0
  baseline as the anchor row;
- relative speedup of each interp backend vs `ExactMath`
  (which itself should match the step-0 baseline within noise);
- relative speedup of each step in the implementation order
  vs the step-0 baseline (useful to confirm steps 2 and 3 don't
  regress);
- correlated with the accuracy report.

The plan does not commit to a specific speedup target; we ship
`InterpUnivariateMath` if it clears the accuracy budget on the
real-data spot check, regardless of headline speedup. The
mixture-loop bench numbers are also the **trigger condition** for
the 3D `log_mix` follow-up: if the bench shows the mixture loop
is still the dominant hot spot after `InterpUnivariateMath`
lands (say, ≥ 40 % of `run_em_for_record` time in the
biallelic-contam-on fixture — informed by the step-0 flamegraph
breakdown), the 3D table becomes worth prototyping. Otherwise we
leave it deferred.

**Pre-implementation speedup estimates, informed by the step-0
flamegraph** (`ln`+`exp` ≈ 40 %, `e_step` arithmetic ≈ 30 % of
self-time on the contam-on fixture):

- **Item A alone (interp `ln`/`exp`)**: attacks ~40 %.
  Conservative estimate ~1.4× on contam-on (math-heavy);
  ~1.3× on no-contam (less math).
- **Item E alone (shape cache + e_step restructure)**: attacks
  the ~30 % `e_step` arithmetic via removed branches, removed
  scans, and zero-skipping. Estimate ~1.15× on its own.
- **A + E combined** (the v1 target): plausibly **~1.6-1.8× on
  contam-on, ~1.4-1.6× on no-contam**. These are guesses; the
  bench from step 6 is the source of truth.

## Implementation order

Step-by-step so each commit is reviewable on its own:

0. **Baseline bench + hot-path sanity check.** Before touching
   any engine code, build the bench harness (`bench_posterior_engine`)
   with the fixture battery and run it against the **current**
   `run_em_for_record` (no `MathBackend` trait, no changes — just
   measure today's code). Also capture one flamegraph (via
   `samply` per [CLAUDE.md §"Profiling"](../../CLAUDE.md), or
   `cargo flamegraph` if `perf_event_paranoid` permits) to
   confirm `ln`/`exp` actually dominate the EM loop and the
   mixture pre-pass. Purposes:

   - Anchor a baseline number per fixture (wall-time / record).
   - Validate that the bench infrastructure (fixtures,
     Criterion config, sample sizes) is stable before we depend
     on it for decisions in later steps.
   - Falsify the hot-path assumption cheaply — if the flamegraph
     points somewhere other than `ln`/`exp` (e.g. allocations,
     genotype enumeration), pause and revise the plan before
     building the trait.

   Deliverables: the bench file committed, baseline numbers
   recorded in the PR description, the flamegraph attached or
   summarised. The bench file initially measures the engine as
   a single function call (no backend dispatch yet); step 1
   converts it to use the trait.

1. **`MathBackend` trait + `ExactMath` impl + engine
   generalisation** (no behaviour change; `Default = ExactMath`).
   Land independently. All existing tests pass. Re-run the
   step-0 bench — numbers should be identical to baseline (this
   confirms the trait introduction doesn't accidentally regress
   via failed inlining).
2. **Hoist `safe_ln(f_s)` to `EmContext`** (item D). Bit-identical;
   ploidy-independent. Land independently. Re-run the bench —
   should be neutral-to-slightly-positive; catches accidental
   regression.
3. **Genotype-shape cache + e_step restructure** (item E). The
   cache adds `nonzero_pairs` and `homozygous_allele_for` to
   `GenotypeShape`; the e_step inner loop is rewritten to
   consume them (removed `k == 0` branch, removed
   `homozygous_allele` scan). Bit-identical. Land
   independently. Re-run the bench — expect a small positive on
   wide-allele fixtures where setup cost matters.
4. **`InterpUnivariateMath` impl** (item A) — the `ln_approx` /
   `exp_approx` tables, plumbed through the trait. No production
   wiring yet — the backend exists and is unit-tested for
   per-call accuracy, but the engine still defaults to
   `ExactMath`.
5. **Accuracy harness** (the `tests_math_backend_accuracy`
   module). Runs `InterpUnivariateMath` against the fixture
   battery; gates on the approximate-parity targets.
6. **Extend the bench harness** to add the
   `interp_univariate` line alongside `exact`. Reports
   side-by-side speedups. (The harness already exists from
   step 0; this is a small additive change.)
7. **Review the bench + accuracy results and decide on the
   default backend.** This is a judgement call that depends on
   the concrete numbers from steps 5 and 6, not on a pre-set
   threshold. Two outcomes:

   - **Flip the default to `InterpUnivariateMath`** if the
     speedup is meaningful (subjective, but as a sanity check:
     ≥ 1.5× wall-time reduction on the biallelic-contam-on
     bench) *and* the accuracy harness clears its budget with
     margin to spare on both the synthetic battery and the
     real-data spot check. `ExactMath` remains as an opt-in
     escape hatch (via `with_math_backend(ExactMath)`).
   - **Keep `ExactMath` as the default** if the speedup is
     marginal, or the accuracy harness is close to its limits
     (max error approaching the budget), or the real-data spot
     check surprised us. `InterpUnivariateMath` stays available
     as an explicit opt-in for callers that want to trade
     precision for speed.

   Either way, no code outside the engine module needs to
   change — the trait + default makes this a single-line edit.
   The decision can be recorded in the PR description for step
   7 along with the bench/accuracy tables that motivated it.

Step 4 is the only one that introduces approximation. Everything
before is bit-identical or trivially equivalent. The 3D
`log_mix` table (the previously-listed step 7) is **not part of
this plan**; if the bench numbers from step 6 meet the trigger
condition for it (see §"Bench strategy / Reporting"), it gets
its own follow-up plan.

## API shape

- `pub mod backends` (re-exported from `posterior_engine`) —
  contains the trait and both impls:
  - `pub trait MathBackend` — public because it appears in the
    engine's type parameter bound.
  - `pub struct ExactMath` (zero-sized) — bit-identical baseline.
  - `pub struct InterpUnivariateMath` — interpolated, with a
    docstring noting that its internals (table layout,
    resolution) may change without a major version bump.
- `pub fn PosteriorEngine::with_math_backend(...)` — escape
  hatch for explicit backend selection.
- The existing constructors (`::new`, `::with_config`) keep
  their signatures and use the defaulted backend.

The internals of the backends (table resolutions, layouts) are
private to the `backends` submodule. Resolutions are `const`
in v1; if a future need arises to make them runtime-tunable we
can add a builder.

## Validation

Inside the dev container (`./scripts/dev.sh`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo bench --bench var_calling_perf -- posterior_engine`
  (informational — captures the speed × accuracy trade-off table
  to paste into the PR).

## Assumptions / silent choices

- **Universal interp domains via float decomposition.**
  Tabling `ln(m)` for `m ∈ [1, 2)` and `2^f` for `f ∈ [0, 1)`
  covers every IEEE 754 finite input. The alternative — tabling
  `ln(x)` over a fixed `x` range and falling back outside —
  works but has worse uniform-accuracy properties. We choose the
  decomposition approach.
- **Static (`'static`) tables.** Computed by a `const fn` or by
  a `OnceLock`-initialised global; no per-engine allocation. The
  table is shared by all `PosteriorEngine<_, InterpUnivariateMath>`
  instances within a process.
- **Trait dispatch is monomorphised.** All hot calls go through
  the trait method, which the compiler inlines. The trait object
  form (`Box<dyn MathBackend>`) is not used in production; if a
  test ever wants dynamic dispatch (e.g. iterating over backends
  in a single test), it can do so locally.
- **No `unsafe`.** All table indexing is bounds-checked; the hot
  loop bounds-check is one branch per call, well-predicted, and
  cheap relative to the lookup itself.
- **Table resolutions are `const`** but their *values* are
  decided by the accuracy harness, not by this plan. Starting
  points: 1024 entries for the `ln` mantissa sub-table, 1024 for
  the `exp` fractional sub-table. The doubling/halving loop (see
  §"The interpolated tables / A / Resolution is a tunable")
  converges to the smallest table that meets
  `posterior_error.max ≤ 1e-4` and `posterior_error.p99 ≤ 5e-5`
  on the full fixture battery including the real-data spot
  check. Changing resolution later is a one-constant edit.
- **Uniform spacing inside the sub-tables, not over raw `x`.** The
  IEEE 754 decomposition is itself the adaptive transform —
  uniform spacing on `m ∈ [1, 2)` gives near-uniform output error
  on `ln(x)` for any `x`. Non-uniform spacing inside the
  sub-table would add a per-call branch/search cost that swamps
  the interp itself. For the 3D `log_mix` table, if uniform
  doesn't meet the budget, the escalation is axis transform →
  piecewise-uniform → drop the table (see §"The interpolated
  tables / B").

## Risks

- **Accuracy regression on borderline records.** Linear interp
  has its largest error where the second derivative is largest —
  for `ln(x)`, that's near `x = 0`. `safe_ln(p)` for very small
  `p` could see noticeable error. Mitigation: the accuracy
  harness includes "near-degenerate" fixtures (samples with
  essentially zero posterior on some genotype) and we measure
  the tail of the error distribution, not just the mean.
- **Out-of-domain fallback masking the problem.** If `ln_approx`
  silently falls back to `f64::ln` for some inputs, a bug in the
  table could be invisible because the fallback covers it.
  Mitigation: the accuracy harness asserts that the in-domain
  path is exercised at least once per fixture (counted via a
  thread-local flag flipped at table-hit time, asserted post-run).
- **Cache miss on unusual `(ploidy, n_alleles)` workloads.** The
  slot array caps caching at `ploidy ≤ 8` and `n_alleles ≤ 16`;
  shapes outside that fall back to per-record construction
  (uncached but correct). Real workloads almost never exceed
  these bounds. Mitigation if a future workload does: bump the
  bounds (one-constant edit) or switch the cache to a HashMap +
  LRU (re-evaluate the trade-off then).
- **Generic-explosion compile time.** Adding a type parameter to
  the engine means every site that names `PosteriorEngine`
  monomorphises per backend. We have one production caller, so
  this is fine in practice; flag if compile time worsens
  visibly.

## Decisions to confirm before implementation

These are the inline `**Open question**`s consolidated:

1. **Table resolution — starting points and tuning loop.**
   The plan commits to a *process* (double until the error
   distribution comfortably clears the budget), not to specific
   numbers. Are the proposed starting points (1024 for the `ln`
   and `exp` mantissa/fraction tables) reasonable, or do we want
   to start higher / lower? Also: are the per-cell error caps
   (`p99 ≤ 5e-5` on posteriors, `p99 ≤ 0.5` on Phred) the right
   secondary gates, or just headline `max` thresholds?
2. **Trigger condition for the 3D `log_mix` follow-up.** The
   plan defers the 3D table and revisits it only if the bench
   numbers say the mixture loop is still the dominant hot spot
   after `InterpUnivariateMath` lands. Proposed trigger: the
   mixture pre-pass is ≥ 40 % of `run_em_for_record` wall-time
   on the biallelic-contam-on bench fixture. Is 40 % the right
   threshold, or should it be higher / lower?
3. **Production default backend.** Decided post-measurement
   (see implementation order step 7), not in advance. The
   criteria — meaningful speedup *and* comfortable accuracy
   margin — are deliberately judgement-based: we pick the
   default that matches what the bench and accuracy harness
   actually show. The plan's only commitment is that whichever
   backend isn't the default remains available via
   `with_math_backend(...)`.
4. **Real-data accuracy fixture.** *Decided for the
   experimentation phase:* fixture lives in `tmp/` (git-ignored
   scratch). The `--ignored` test reads from a hard-coded path
   and skips cleanly if absent. Revisit only if the spot check
   needs to become reproducible across machines (committed
   `.psp` fixture, downloader script, etc.) — out of scope for
   this plan.
5. **`MathBackend` visibility.** *Decided:* trait + both impls
   live in a `pub mod backends` submodule. The trait must be
   `pub` (it appears in the engine's type-parameter bound); the
   impls are `pub` because the bench harness lives in a separate
   crate and `pub(crate)` would block it from constructing them.
   The submodule segregation plus a stability docstring on
   `InterpUnivariateMath` signals that the set of backends and
   their internals are optimisation details that may change.
6. **Cache cap for genotype shapes.** *Decided:* fixed slot
   array indexed by `(ploidy, n_alleles)` with bounds `ploidy ≤
   8` and `n_alleles ≤ 16`. No LRU, no eviction; out-of-bounds
   shapes fall back to per-record construction. Rationale: the
   cache is touched once per record (far from the hot path), the
   `(ploidy, n_alleles)` distribution is heavily concentrated,
   and the slot array is simpler with no LRU bookkeeping. If a
   future workload pushes past these bounds or if profiling
   identifies the cache lookup itself as a bottleneck (it
   shouldn't), revisit with a HashMap + LRU.

## Out-of-scope follow-ups

- **3D `log_mix` table for the mixture pre-pass.** Indexes
  `(p_own, c_s, q_b)` and returns `ln((1-c_s)*p_own + c_s*q_b)`
  via trilinear interpolation, folding the mixture arithmetic
  into the lookup itself. Considered for v1 and deferred because
  the cost trade-off (8 cache reads + trilinear interp vs ~3
  multiplies + one `ln_approx` call) isn't obviously favourable
  on modern CPUs, and the 3D non-uniform-curvature problem (the
  function has a near-singularity as `mix → 0`, plus real-data
  clustering on each axis) makes accuracy tuning fiddly.
  **Trigger to revisit:** the bench from step 6 of the
  implementation order shows the mixture pre-pass is still ≥ 40 %
  of `run_em_for_record` wall-time on the biallelic-contam-on
  fixture (i.e. the univariate `ln_approx` alone didn't shift the
  bottleneck away from this loop). If triggered, sketch:
  - Starting resolution `16 × 8 × 16` (16 KB, L1-resident).
  - Trilinear interpolation; 8 corner reads per lookup.
  - Escalation path if uniform spacing fails the accuracy
    budget: (a) axis transform — table in `ln(c_s + δ)` instead
    of `c_s`; (b) piecewise-uniform `c_s` axis with finer
    resolution near 0; (c) abandon.
  - Adds a defaulted `log_mix(p_own, c_s, q_b)` method on
    `MathBackend` (default: arithmetic + `self.ln(...)`); a new
    `InterpLogMix3DMath` impl overrides it. Mixture call site
    switches from `self.ln(mix)` to `self.log_mix(p_own, c_s,
    q_b)`. Strictly additive — no rewrite of v1 code.
- **4D posterior LUT for the standard biallelic case.** Discussed
  in §"Why not a multidimensional posterior LUT" and rejected
  for v1 because it ducks the contamination case. Revisit if the
  no-contam biallelic case is the bottleneck after this lands.
- **Diploid structural specialisation.** The previous draft's
  proposal. Compose with this work cleanly; revisit if the
  bench shows further headroom.
- **Adaptive grid for `ln_approx` near `x = 0`.** If the tail
  error becomes a problem, switch to a piecewise table with
  finer resolution near zero.
- **Per-batch / per-chromosome warm caches for `log_mix`** when
  `c_s` and `q_b` are constant within a batch. Collapses the
  deferred 3D lookup to 1D in `p_own`; potentially worth a
  specialised path if the 3D follow-up is triggered.
- **Rayon over records — coarse-grained parallelism (the "hammer").**
  `PosteriorEngine` is single-threaded today; each `next()` call
  drives one record through `run_em_for_record`. Records are
  independent (no cross-record state during the EM), so the
  natural parallel shape is to buffer `K` upstream `MergedRecord`s
  per `next()` refill and process them with `rayon::par_iter`.
  Output order can be preserved by tagging records with the
  upstream index and reassembling, the same pattern Stage 5
  ([per_group_merger.rs:484](../../src/var_calling/per_group_merger.rs#L484))
  already uses. Expected speedup: near-linear with core count up
  to memory-bandwidth saturation on the `MergedRecord` clones
  and the per-record `Vec` allocations. Likely the **biggest
  single throughput improvement available** to Stage 6, and
  also the easiest to land — most of the engine code is unchanged,
  the surface change is in the iterator driver. **Cleanly
  multiplicative with the interp tables + e_step restructure**:
  the per-thread cost stays exactly what those optimisations
  brought it down to, and rayon multiplies the total throughput
  by core count. (Note: the cohort-level `p_hat` aggregation
  across samples is *inside* one record's EM, so it doesn't
  cross record boundaries — `p_hat` doesn't need synchronisation
  in this design.)
- **SIMD over samples within a record — fine-grained
  parallelism (the "glove").** Inside one record, the E-step
  does identical work for each of the `n_samples` samples. The
  per-cell `exp` at the normalisation step, the per-sample
  `log_sum_exp_slice`, and the elementwise `ll[s, g] +
  log_prior[s, g]` all SIMD-vectorise across samples in lanes
  of 4 (AVX2) or 8 (AVX-512). `log_indep[g]` is a scalar that
  broadcasts perfectly to every lane. The per-sample data
  (`log_f[s]`, `log_one_minus_f[s]`, `ll[s, g]`) needs to be
  laid out genotype-major (contiguous in the sample dimension);
  one per-record transpose of the engine buffers (a few KB,
  L1-resident) sets that up. **Composes multiplicatively with
  the interp tables**: a vectorised `exp_approx_x4` (single
  vector gather + lane-parallel interpolation) is the natural
  endgame for item A's interp `exp`. Expected on-thread
  speedup ~2-3× on top of items A + E. Surface change is large
  — buffer layout, `wide`/`std::simd` dependency, AVX2/NEON
  multiversioning — but the math contract is unchanged.
  **Composes with rayon-over-records as well**: rayon multiplies
  throughput by core count; SIMD multiplies per-thread
  throughput by lane count. They're orthogonal — the "hammer"
  parallelises tasks, the "glove" tightens each task.
- **Per-record allocator amortisation.** The engine allocates
  `posteriors`, `p_hat`, `f_hat_compound`, and several `Vec`s
  per `run_em_for_record` call. With rayon-over-records, a
  per-thread scratch arena (or a `thread_local` `Vec` pool)
  removes the malloc/free traffic the step-0 flamegraph showed
  in the bench-loop allocator share. Independent of the math
  story, but the natural pairing with the rayon follow-up.

## File touch list

- `src/var_calling/posterior_engine.rs` — main engine change.
  Adds `MathBackend` trait usage at the `safe_ln` / `ln` / `exp`
  call sites; adds the `safe_ln(f_s)` hoist and the genotype-shape
  cache; adds the generic type parameter on `PosteriorEngine`.
- `src/var_calling/posterior_engine/backends.rs` (new, re-exported
  as `pub mod backends` from `posterior_engine`) — `MathBackend`
  trait, `ExactMath`, `InterpUnivariateMath`. Per-backend unit
  tests for per-call accuracy (`ln_approx(x) ≈ x.ln()` to within
  table tolerance).
- `src/var_calling/posterior_engine/interp.rs` (new, private) —
  the 1D linear-interpolation primitive used by
  `InterpUnivariateMath`.
- `benches/var_calling_perf.rs` (or new
  `benches/posterior_engine_perf.rs`) — bench harness. Lands at
  step 0 with the baseline ("today's code") run; grows to
  side-by-side backend comparison after step 6.

No other files are touched. `genotype_order` in
[per_group_merger.rs](../../src/var_calling/per_group_merger.rs#L355)
is reused as-is by the genotype-shape cache.
