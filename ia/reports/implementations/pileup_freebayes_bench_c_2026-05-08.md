# Pileup walker — freebayes-style alternative benchmark, Scope C

**Date:** 2026-05-08
**Investigator:** Claude (research-only session)
**Plan:** [`ia/feature_implementation_plans/pileup_freebayes_style_benchmark.md`](../../feature_implementation_plans/pileup_freebayes_style_benchmark.md)
**Motivating review:** [`ia/reviews/pileup_freebayes_comparison_2026-05-08.md`](../../reviews/pileup_freebayes_comparison_2026-05-08.md)
§3 `A3` and §6.

## Verdict

**Do not proceed to Scope A.** The microbenchmark refutes the
plan's going-in prediction ("tie within noise") in the strongest
possible direction: our open-record bookkeeping data structure
wins decisively at every cell above the smallest, and the gap
widens with both coverage depth and footprint length. At the
realistic-Illumina mixed cell (10% short-footprint events,
coverage=30) the production-shape impl is **2.05× faster** than a
freebayes-shape flat-pool impl. At deletion-heavy cells the gap
reaches **17×**.

Per the plan's exit criteria, this is the "OursStyle wins by ≥ 15%
at high coverage / deletion-heavy" outcome — the strongest possible
defence of the design. The freebayes review's `A3` and §6 arguments
were architectural; this scope adds empirical backing.

## What was built

[`benches/freebayes_bookkeeping.rs`](../../../benches/freebayes_bookkeeping.rs)
— a criterion microbenchmark of two open-record bookkeeping data
structures driven by the same synthetic event stream:

- **`OursStyle`** — `BTreeMap<u32, OurSlot>` keyed by anchor.
  Events whose anchor falls inside an existing record's footprint
  merge into that record; events at fresh anchors open new slots.
  Ageing uses an early-break drain that relies on the
  disjoint-footprint invariant (see file's inline comment for the
  proof sketch).
- **`FreebayesStyle`** — `Vec<PoolEntry>`, one entry per event.
  No merging on insert. Ageing is `pool.retain(|e| e.anchor +
  e.ref_span > walker_pos)` — full linear sweep per step. Per-step
  also does a linear filter to count alleles overlapping the
  walker (freebayes' per-position genotyping query, kept under
  `black_box` so the optimiser cannot elide it).

Both impls process the same pre-built event stream. **No real
walker, no real reads, no fold work, no scalar arithmetic, no
reference fetch.** Per-event work is identical (`u32` increments
on a fixed allele-id buckets array of size 2). The only thing
being measured is the data-structure access pattern.

A `validate_parity()` runs at every `cargo bench` startup,
asserting both impls drain the same total event count for each
distribution × coverage combination at moderate coverages.

### Bench parameters

- `SPAN = 50_000` — number of walker steps. Matches
  [`benches/pileup_walker_scaling.rs`](../../../benches/pileup_walker_scaling.rs).
- Coverages: 10, 30, 100, 500, 1000.
- Distributions:
  - **`pure_snp`** — every event has `ref_span = 1`. Records open
    and close in adjacent steps; active set bounded by one step's
    merge.
  - **`mixed`** — 90% `span = 1`, 10% `span ∈ [2, 10]` uniform.
    Mimics typical Illumina BAQ-adjusted output.
  - **`deletion_heavy`** — 50% `span = 1`, 50% `span ∈ [10, 100]`
    uniform. Pathological regime that exercises the linear-pool
    cost in `FreebayesStyle`.
- Criterion config: `sample_size = 10`, `measurement_time = 3 s`.

### Discipline

Both runs taken inside the dev container at the same machine
session, ~10 minutes apart. Per the noise-floor rule from
[`pileup_lazy_cigar_2026-05-07.md`](pileup_lazy_cigar_2026-05-07.md)
("A note on noise and baselines"), deltas under 15% on a
single-run pair are noise. Two-run agreement on the headline
ratios (mixed: 2.05× vs 1.97× across runs; deletion_heavy at
cov=100: 16.8× vs 14.6×) confirms the gaps are stable.

## Results

Median wall time per criterion cell. Each (impl, distribution,
coverage) cell measures the full 50 000-step run on a freshly-built
event stream.

### `pure_snp` — all events `ref_span = 1`

| Coverage | Ours run1 | Ours run2 | FB run1 | FB run2 | FB / Ours (mean) |
|---:|---:|---:|---:|---:|---:|
|  10 |   1.92 ms |   1.99 ms |   1.96 ms |   1.93 ms | **0.99×** (tied) |
|  30 |   3.88 ms |   3.69 ms |   4.65 ms |   4.26 ms | **1.18×** |
| 100 |  10.73 ms |  10.38 ms |  13.67 ms |  14.49 ms | **1.33×** |
| 500 |  46.41 ms |  44.45 ms |  69.16 ms |  63.66 ms | **1.46×** |
| 1000|  91.63 ms |  88.06 ms | 129.01 ms | 138.81 ms | **1.49×** |

`OursStyle` ties `FreebayesStyle` at the smallest cell (one event
per step on average for either impl, both essentially measuring
allocator + per-step overhead) but pulls ahead linearly with
coverage. The widening gap is the per-position overlap query in
`FreebayesStyle` — `pool.iter().filter(...)` over `coverage`
entries every step — that has no analogue in `OursStyle`'s
"the slot is the answer" lookup.

### `mixed` — 90% `span=1`, 10% `span ∈ [2, 10]`

| Coverage | Ours run1 | Ours run2 | FB run1 | FB run2 | FB / Ours (mean) |
|---:|---:|---:|---:|---:|---:|
|   10 |   2.06 ms |   2.08 ms |   3.81 ms |   3.87 ms | **1.86×** |
|   30 |   3.89 ms |   4.30 ms |   8.32 ms |   8.46 ms | **2.05×** |
|  100 |  12.18 ms |  13.02 ms |  25.54 ms |  26.57 ms | **2.07×** |
|  500 |  59.48 ms |  61.31 ms | 134.61 ms | 131.53 ms | **2.20×** |
| 1000 | 114.07 ms | 117.50 ms | 256.75 ms | 243.18 ms | **2.16×** |

The realistic-Illumina regime. A stable **~2×** gap across the
full coverage sweep. Even at 10% short-footprint events, the
`FreebayesStyle` pool grows enough that the per-step linear sweep
+ overlap query become dominant.

### `deletion_heavy` — 50% `span=1`, 50% `span ∈ [10, 100]`

| Coverage | Ours run1 | Ours run2 | FB run1 | FB run2 | FB / Ours (mean) |
|---:|---:|---:|---:|---:|---:|
|   10 |   2.04 ms |   1.73 ms |  22.10 ms |  18.67 ms | **10.85×** |
|   30 |   4.46 ms |   4.64 ms |  59.20 ms |  55.72 ms | **12.62×** |
|  100 |  12.49 ms |  12.30 ms | 209.25 ms | 179.12 ms | **15.66×** |
|  500 |  55.09 ms |  62.72 ms | 997.40 ms | 1065.20 ms | **17.51×** |
| 1000 | 111.80 ms | 115.64 ms | 1899.10 ms | 1971.70 ms | **17.02×** |

The pathological regime, designed precisely to stress
`FreebayesStyle`'s linear-pool weakness. The result is exactly
that: each long-footprint event keeps an entry alive for 10-100
walker steps, the pool grows ~50× larger than the open-records
set, and `pool.retain` + per-position filter cost dominate. The
ratio plateaus around 17× because at high coverage both impls are
also paying for per-event push/insert work, which scales the
same way on both sides.

`OursStyle`'s per-step time stays sub-linear in pool size because
merge-on-overlap absorbs every Match event during a deletion's
footprint into the deletion's existing slot — the BTreeMap
doesn't grow during the long deletion's life.

## Interpretation

### What the gap actually measures

Three forces compound to widen the gap:

1. **`FreebayesStyle` does an `O(N)` linear sweep every step.**
   `pool.retain` walks the entire pool. At deletion-heavy
   coverage=1000 the pool size hovers around 50 × 1000 = 50 000
   entries, so each of the 50 000 walker steps does a 50 000-entry
   sweep — `O(N²)` overall.

2. **`FreebayesStyle` does a second `O(N)` linear filter every step**
   for the per-position genotyping query. Same `N`, doubles the
   per-step cost.

3. **`OursStyle`'s merge-on-overlap collapses many events into
   one slot.** During a 50-base deletion's lifetime, `coverage` Match
   events at each of the 50 walker positions on top of it all merge
   into the single deletion slot — the BTreeMap stays at one entry
   for the whole deletion. The pool, by contrast, gets one entry per
   event regardless.

The two `O(N)` sweeps and the merging together explain the
~17× plateau in deletion-heavy: it's not that `OursStyle` is
faster per-event; it's that `OursStyle` does much less per-event
work because most events fold into existing structure.

### Why this is honest about the design choice

The architectural argument in
[`pileup_freebayes_comparison_2026-05-08.md` §6](../../reviews/pileup_freebayes_comparison_2026-05-08.md)
was: "freebayes' lazy-pool model is correct for freebayes' output
contract; ours is correct for ours." It sat next to a per-record
close boundary and said "this is necessary for what we emit, and
the cost of maintaining it is small." Scope C now produces the
empirical claim: **the cost isn't just small — the alternative is
much more expensive.** The choice is performance-positive, not
performance-neutral.

### Caveats

This is a microbenchmark. It measures pure data-structure cost
in isolation. Specifically:

- **No fold work.** The real walker spends ~367 ns per fold
  iteration on `apply_events_to_ref`, scalar arithmetic, and
  per-allele bookkeeping
  ([`pileup_fold_cache_2026-05-07.md`](pileup_fold_cache_2026-05-07.md)).
  A freebayes-style production walker would carry the same fold
  cost; the bench's data-structure gap would be *diluted* in the
  full pipeline by the fold work both impls share.
- **No real reads.** The synthetic event stream gives every event
  `anchor = walker_pos`, matching the cursor's emission shape but
  abstracting away the cursor walk itself.
- **The freebayes-shape impl is a clean Rust port.** Freebayes'
  C++ uses `vector<Allele*>` (one indirection per access),
  `Allele*` nullification then `erase(remove(NULL))` (two passes
  vs. our one `retain`), and `map<long, deque<RegisteredAlignment>>`
  (deque indirection). All would amplify the gap further on
  freebayes' actual code, not narrow it.
- **The 2× mixed gap is the most decision-relevant number.** The
  17× deletion-heavy gap is a stress test that real BAQ-adjusted
  Illumina won't produce; the 2× mixed gap is the realistic
  steady-state expectation for the bookkeeping component.

So the report's claim is bounded:

> **In a microbenchmark stripped of all per-fold work and fold-time
> allocations, the production walker's BTreeMap-of-disjoint-records
> data structure outperforms a freebayes-style flat-allele-pool
> data structure by 1.4-2.2× at typical Illumina depths and by an
> order of magnitude on long-footprint stress tests.**

The gap when full fold work is layered on (Scope B) would be
smaller in proportion (because fold cost is shared) but would
still strictly favour `OursStyle` — the bookkeeping component of
total time only ever shrinks; it never goes negative.

## What the data says about the next scope

The plan's exit criteria for Scope C:

| Outcome | Action |
|---|---|
| Tie within noise across all cells | Don't proceed; close issue with paragraph in review. |
| FreebayesStyle wins ≥ 15% at deep / deletion-heavy | Proceed to A. |
| OursStyle wins ≥ 15% at deep / deletion-heavy | Don't proceed; strongest defence of design. |

The third row matches the result. **Do not proceed to Scope A.**

If a future signal does motivate revisiting (e.g. a real-data
profile shows the BTreeMap traversal is hot in production), the
right next step would still not be Scope A's degraded-output
walker — it would be a more focused experiment on whether
*different* open-record data structures (e.g. `Vec<OpenSlot>` with
tombstones, or `IndexMap`, or a slab allocator backing the BTreeMap)
beat the current `BTreeMap<u32, OpenSlot>`. Scope A and B as
specified in the plan only measure the freebayes-shape vs.
ours-shape distinction; they cannot inform a within-shape choice.

## Artefacts left in the tree

- [`benches/freebayes_bookkeeping.rs`](../../../benches/freebayes_bookkeeping.rs)
  — kept as a regression guard. The bench takes ~5 minutes total
  to run and exercises a corner of the design space the
  production bench (`pileup_walker_scaling.rs`) does not. If the
  open-record data structure is changed in a future commit, this
  bench will catch a regression of the merge-on-overlap property
  by showing `deletion_heavy` cells climbing toward
  `FreebayesStyle` wall times.

- `Cargo.toml` — one new `[[bench]]` stanza.

No production code touched. No public API change.

## Follow-ups

None forced by this report. Two optional ones:

1. **Add the C result to the freebayes review's `A3`.** A
   one-paragraph update citing this report's median gap on the
   `mixed` and `deletion_heavy` distributions, so the
   architecture argument is paired with the empirical one.
   (See companion edit accompanying this report.)

2. **If/when long-read input lands**, re-run this bench with a
   distribution that includes `span ∈ [1000, 5000]` events
   (currently outside scope because `MAX_RECORD_SPAN = 5000`).
   The gap would widen further at that regime, since the
   freebayes-shape impl's per-step `O(N²)` sweep cost grows
   quadratically with span. Recording the prediction here so a
   future engineer can verify.

## Reproducing

```sh
./scripts/dev.sh cargo bench --bench freebayes_bookkeeping
```

Expected wall time: ~5 minutes for one full pass. Two passes
recommended; the discipline rule from
[`pileup_lazy_cigar_2026-05-07.md`](pileup_lazy_cigar_2026-05-07.md)
applies — deltas under 15% on a single-run pair are noise.
