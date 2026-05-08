# Pileup walker — freebayes-style alternative benchmark

Implementation plan for the empirical follow-up to
[`ia/reviews/pileup_freebayes_comparison_2026-05-08.md`](../reviews/pileup_freebayes_comparison_2026-05-08.md)
§3 `A3` and §6: confirm or refute the architectural choice of the
record-closure model over freebayes' lazy-allele-pool model with
benchmark numbers.

This is an **experimental / research plan**, not a feature plan.
The primary deliverable is a yes/no answer to "does our walker's
per-record close boundary cost meaningfully more than freebayes'
per-step pool sweep?" Code may or may not survive — most of it
is throwaway. The user-facing system remains unchanged
regardless of the outcome.

## Status

- **Scope C (open-record bookkeeping microbenchmark):** open. Run first.
- **Scope A (degraded-output walker):** conditional on C.
- **Scope B (full-output walker):** conditional on A.

## Motivation

The freebayes review's §6 argument for keeping the record-closure
model is currently architecture-only: "our output contract requires
per-record close boundaries; freebayes' position-snapshot model does
not." That argument is correct as far as it goes, but it does not
address the secondary question some readers will ask:

> **Even if we have to pay the per-record close-boundary cost, how big
> is it compared to freebayes' per-step pool sweep?**

If the answer is "tiny" the architecture argument stands alone. If
the answer is "actually, freebayes' approach is faster at our depths"
the architecture argument becomes a trade-off rather than a free
choice, and the report needs an honest paragraph about it.

The point of this plan is to spend a bounded amount of effort
producing the number, by escalating only as far as the cheap
experiments justify.

## What we already know

From the existing perf history, before any of the work in this plan
runs:

- **Per-position cost on our walker is 9-16 µs** depending on op
  count and read length, measured on
  [`benches/pileup_walker_scaling.rs`](../../benches/pileup_walker_scaling.rs)
  with 50 kb window and 30× coverage. From
  [`pileup_lazy_cigar_2026-05-07.md`](../reports/implementations/pileup_lazy_cigar_2026-05-07.md).
- **Per-fold budget is ~367 ns**, of which `events_overlapping`
  takes ~150-200 ns. The rest is `apply_events_to_ref` allocation
  + scalar arithmetic + allele-bucket lookup. From
  [`pileup_fold_cache_2026-05-07.md`](../reports/implementations/pileup_fold_cache_2026-05-07.md).
- **Open-record bookkeeping is not the dominant cost.** The
  fold-cache study found 0% redundant re-folds on bench fixtures —
  each `(record, read)` pair is folded exactly once on real-shape
  workloads. The BTreeMap of open records is small and accessed a
  small constant number of times per step.
- **Run-to-run noise floor is 15-25%** on the multi-op bench
  fixtures (same bench file documents this). Any delta below 15%
  is noise.

The headline implication: **most of our per-position time is
per-event work that a freebayes-style walker would also have to
do.** The bookkeeping overhead — which is what differs between the
two designs — is a thin slice. Predicted outcome of every scope in
this plan, going in: **tie within noise on Scope C; widening
margin in our favour on B as the per-fold work grows.** Findings
that contradict either prediction are the interesting ones to act
on.

## Common ground across all scopes

Every scope uses the same bench infrastructure, language, and
machine.

- **Bench framework:** criterion, same as
  [`benches/pileup_walker_scaling.rs`](../../benches/pileup_walker_scaling.rs).
- **Build profile:** release, inside the dev container
  (`./scripts/dev.sh cargo bench --bench …`).
- **Discipline:** two runs at every codepoint before drawing
  conclusions; deltas under 15% on a single-run pair are noise.
  This rule is recorded in the lazy-CIGAR report's "A note on noise
  and baselines" section and applies to every measurement here.
- **Reporting:** each scope produces a short report under
  [`ia/reports/implementations/`](../reports/implementations/)
  named `pileup_freebayes_bench_<scope>_<date>.md`, with the same
  tables-and-narrative shape as the lazy-CIGAR report. The report
  is what survives even if the code does not.
- **Code disposition:** Scope-C code lives in `benches/` alongside
  the existing harness. Scope-A / Scope-B walker code lives behind
  a private module gated by `cfg(feature = "freebayes_bench")` so
  it does not affect the production crate's surface area or
  compile time. Whether to delete it post-measurement is a Scope-end
  decision; default is delete unless we kept finding new questions
  to throw at it.

## Scope C — open-record bookkeeping microbenchmark

### What it measures

The pure data-structure cost difference between:

- **`BTreeMap`-keyed open-record set** (our model): close-boundary
  lookup, targeted insert/update at anchor, range-erase on close.
- **Flat `Vec<Allele>` pool with per-step sweep** (freebayes' model):
  linear filter for "alleles overlapping `walker_pos`",
  linear filter for "alleles whose footprint is past `walker_pos`",
  push-back insert.

No fold work. No scalars. No phase chains. No real reads. No
walker. The two implementations are driven by the same synthetic
event stream and produce the same per-step counts (records-closed,
records-emitted), with the counts cross-checked for parity.

### What it builds

A new file `benches/freebayes_bookkeeping.rs` with:

1. **Synthetic event generator.** Deterministic — seeded RNG —
   produces a stream of `(walker_pos, anchor, ref_span, allele_id,
   read_id)` tuples. Parameters:
   - `span` (number of walker positions traversed),
   - `coverage` (events per position),
   - `footprint_dist` — distribution of `ref_span` (1 for SNP-like,
     wider for deletion-like, mixed for realistic).
   - `compound_rate` — fraction of positions where multiple events
     fall on the same anchor (forces REF widening on our side / a
     larger pool entry on theirs).

2. **`OursStyle` impl.**
   ```rust
   struct OursStyle {
       open: BTreeMap<u32 /* anchor */, OpenSlot>,
   }
   ```
   `OpenSlot` carries a small `SmallVec<AlleleEntry>` (1-4 entries
   typical) and the current `ref_span`. Two operations:
   - `step(walker_pos, events)` — for each event: locate or create
     the open slot at `event.anchor`, possibly widen its
     `ref_span`, insert/update the allele entry. Then range-erase
     all slots with `anchor + ref_span ≤ walker_pos`, counting
     them as closed.
   - `finalize()` — drain everything still open.

3. **`FreebayesStyle` impl.**
   ```rust
   struct FreebayesStyle {
       pool: Vec<AlleleEntry>,
   }
   ```
   Two operations:
   - `step(walker_pos, events)` — for each event: push to `pool`.
     Then do a `pool.retain(|a| a.anchor + a.ref_span > walker_pos)`
     sweep, counting drops as closed. Then linear-filter
     `pool.iter().filter(|a| a.anchor <= walker_pos && walker_pos
     < a.anchor + a.ref_span)` to count "alleles overlapping
     walker_pos" (the freebayes per-position genotyping query).
   - `finalize()` — same.

4. **Parity oracle.** Each `step` call, both impls produce the
   same counts (records closed, records overlapping). A debug
   assert in the bench fixture catches divergence; bench code
   itself runs in release with the assert compiled out.

5. **Criterion harness.** Bench groups parameterised over
   `coverage ∈ {10, 30, 100, 500, 1000}` and over three
   `footprint_dist`s: `pure_snp` (span=1), `mixed` (90% span=1, 10%
   span=2-10 uniform), `deletion_heavy` (50% span=1, 50% span=10-100
   uniform). `span = 50_000` to match the existing bench. Output:
   per-step time and total time per `(impl, coverage, dist)` cell.

### What success / failure looks like

Three outcomes:

- **Tie within noise across all cells** (predicted). Conclude C is
  consistent with the existing fold-cache finding: bookkeeping is
  not a meaningful cost driver at our depths. **Do not proceed
  to A.** Update the freebayes review's `A3` with a one-paragraph
  summary citing the C report. Done.

- **`FreebayesStyle` wins at high coverage / deletion-heavy cells
  by ≥ 15%.** This is the interesting outcome. It would mean the
  `BTreeMap` overhead is non-trivial when the open-record set is
  large (deletion-heavy regions hold many records open
  simultaneously) and the per-step lookup cost dominates.
  **Proceed to A.**

- **`OursStyle` wins at high coverage / deletion-heavy cells by
  ≥ 15%.** The opposite of the above — `Vec::retain` linear sweep
  is paying its way at depth. Update the review with the result;
  **do not proceed.** This would be the strongest possible
  defence of the design as is.

In all three outcomes the C report is the deliverable; the bench
file may or may not be kept depending on whether we want
regression-guard coverage on the bookkeeping cost.

### Cost estimate

- ~1 day to build the harness and the two impls.
- ~1 day to run, double-check, write the report.
- **Total: 2 days.**

## Scope A — degraded-output walker

### What it measures

Same comparison as Scope C, but with real `PreparedRead` input
flowing through the actual CIGAR cursor and active-set machinery
on both sides. The output is degraded: per-position depth counts
only — no scalars, no phase chains, no REF widening. This isolates
the cost of the *walker shape* (close-boundary vs. lazy-aging)
from the per-fold scalar work, while still exercising the real
read-decomposition path.

### What it builds

A new private module `src/per_sample_caller/pileup/freebayes_bench.rs`
gated behind `#[cfg(feature = "freebayes_bench")]`, plus a
companion bench in `benches/freebayes_walker_degraded.rs`.

The freebayes-style walker:

```rust
pub(super) struct FreebayesStyleWalker<R: RefBaseFetcher> {
    active: ActiveSet,         // reused unchanged
    allele_pool: Vec<PoolEntry>,
    ref_fetcher: R,
    /* + the existing per-step plumbing minus folded_reads, open_records */
}

struct PoolEntry {
    anchor: u32,
    ref_span: u32,
    read_id: u32,
    /* enough state to count "this entry was observed" — no scalars */
}
```

Steps per walker position:
1. Admit reads as today.
2. For each active read, query the cursor for events at `walker_pos`.
3. Push each event into `allele_pool` as a `PoolEntry`.
4. Sweep: `allele_pool.retain(|e| e.anchor + e.ref_span > walker_pos)`.
5. Emit a per-position depth count: how many entries in
   `allele_pool` overlap `walker_pos`.
6. Advance.

The output is **just per-position depth** — a single u32 per
position. The production walker, for the comparison, runs in a
mode that emits the same per-position depth (drop the scalar
accumulation — a `cfg(feature = "freebayes_bench")` branch in
`open_record::process_position` that bypasses the fold and just
reports `affected_records.iter().map(record_depth).sum()`).

### What success / failure looks like

- **Tie within noise** (predicted again, but with a wider margin
  than C because both impls now share the cursor walk and fold
  bypass). **Do not proceed to B.**
- **`FreebayesStyle` wins by ≥ 15%.** The walker shape itself is
  costing us. **Proceed to B** to confirm the gap survives full
  output parity.
- **`OursStyle` wins by ≥ 15%.** Confirms the close-boundary
  bookkeeping is paying its way once real read decomposition
  is in the loop. Update the review; **do not proceed.**

### Cost estimate

- ~2 days to build the freebayes-style walker (the active-set,
  cursor, ref-fetch plumbing all reuse our code; the new code is
  the pool + sweep + degraded emit).
- ~1 day to wire the production walker's degraded-output mode.
- ~1 day to run and report.
- **Total: 4 days.**

### Honesty caveat

A degraded-output comparison is informative but not conclusive: it
measures walker-shape cost without measuring scalar / chain / REF
widening cost. A freebayes-style walker that produced full output
would also have to do that work, just structured differently. The
gap A measures is therefore an **upper bound** on the walker-shape
gap in B — if A is a tie, B will be a tie or better for
`OursStyle` (because B adds work to `FreebayesStyle` that
`OursStyle` already amortises through `folded_reads`). If A is a
freebayes-style win, B's outcome depends on how the missing
machinery distributes.

## Scope B — full-output walker

### What it measures

Honest end-to-end comparison: a freebayes-style walker that produces
**bit-identical `PileupRecord` output** to the production walker, on
the same `PreparedRead` stream. This means adding back, on the
freebayes side, all the machinery our spec requires and freebayes
itself does not have:

- per-allele scalar accumulation (the five scalars per allele);
- REF-span widening when overlapping events compound, with allele
  literals rewritten against the wider span;
- phase chain id allocation and propagation matching our model;
- close-boundary emission (because the output unit is a closed
  `PileupRecord`, not a per-position genotype call — freebayes' own
  output shape is not what we are measuring against).

The freebayes-shape part is the *internal data structure choice*
(flat allele pool with lazy aging, no `BTreeMap` of open records).
The output, the input, and the per-allele math are identical to
the production walker.

### What it builds

Continuation of the `freebayes_bench` feature module from Scope A,
extended to:

1. Carry per-allele scalars on `PoolEntry` and update them per
   contributing read.
2. Detect and execute REF widening when two pool entries' footprints
   overlap on the same anchor (or migrate the entries to share an
   anchor with widened span — the equivalent of our
   `open_record::extend_ref_span`).
3. Allocate / propagate phase chain ids using the existing
   `SlotAllocator` (reused unchanged).
4. Group pool entries that share `(anchor, ref_span)` into emit
   buckets when the close boundary is reached, and produce
   `PileupRecord` values byte-identical to the production walker's
   for the same input.

The bench emits both walkers' `PileupRecord` streams, asserts
parity (same chrom/pos/REF/alts/scalars/chain_slots in the same
order), and times the production of the stream end-to-end.

### What success / failure looks like

- **Tie within noise** (most likely). Confirms the design choice
  is performance-neutral and rests on output-contract fit, exactly
  as the review's §6 already argues. The win is converting "§6 is
  an architecture argument" into "§6 is an architecture argument
  with empirical backing." Update the review.

- **`FreebayesStyle` wins by ≥ 15%.** Genuinely surprising. Would
  warrant a follow-up pass on the production walker to import the
  flat-pool data structure (without giving up our output shape —
  the win is in the data structure, not the output contract).

- **`OursStyle` wins by ≥ 15%.** Confirms our design choice
  decisively. Stronger version of the C/A outcome.

### Cost estimate

- Scope-A work as prerequisite.
- ~3-4 days to add scalar accumulation, REF widening, chain
  propagation, and emit-bucket logic to the freebayes-style walker.
- ~1-2 days to drive both walkers on the existing `pileup_walker_scaling`
  fixtures plus a deletion-heavy fixture (the one S2 already added
  to test for; we need the deletion-heavy regime in B specifically
  because that is where the open-record set is largest and where
  the `BTreeMap` overhead would, if anywhere, become visible).
- ~1 day to run and report.
- **Total: 5-7 days on top of Scope A → 9-11 days end-to-end.**

### Honesty caveat

Even at scope B, this comparison still has one structural unfairness:
the production walker has been through several perf passes (lazy
CIGAR, binary-search cursor, S5 caps), while the freebayes-style
walker would be a fresh implementation. To make B a fair test, the
freebayes-style walker must use the same `CigarCursor`, the same
mate-overlap resolver, the same scalar arithmetic — only the
open-state data structure differs. This is the discipline the
plan commits to.

A 15% B-side win on a fresh implementation that hasn't been perf-passed
should be treated as "investigate further," not "rewrite."

## Suggested execution order

1. **Run Scope C** — 2 days, low risk. Definitive answer on the
   pure data-structure question. Likely outcome closes the issue.
2. **If C says proceed, run Scope A** — 4 days. Likely outcome
   closes the issue.
3. **If A says proceed, run Scope B** — 5-7 days. Definitive
   end-to-end answer.

Each scope produces its own report and a clean stop point. The
plan is bounded in cost: worst case is 9-11 days of work spread
across 2-3 weeks (allowing for runs / write-ups / interrupts);
best case is 2 days and a paragraph in the review.

## Validation across all scopes

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib pileup` — must remain green at every commit.
  Scope-A and Scope-B walker code must not affect the production
  walker's tests; the gating feature flag is what enforces that.
- Parity oracle on every bench: each impl's output cross-checked
  against the other (Scope C: same closed/overlapping counts;
  Scope A: same per-position depth; Scope B: bit-identical
  `PileupRecord` stream).

## Out of scope

- **Long-read input.** `MAX_RECORD_SPAN = 5000` and the CRAM decoder
  limits both stay as-is; the existing perf history already shows
  the production walker is not the long-read bottleneck. The
  freebayes-style walker would inherit the same constraint via the
  shared `PreparedRead` input.
- **Calibration / accuracy comparison.** This plan measures
  performance, not output quality. Both walkers, by construction,
  produce identical output (Scope B parity oracle); accuracy work
  is its own track.
- **Adopting the freebayes-style walker as production.** Even if
  Scope B shows it wins, the decision to switch is downstream of
  this plan — it would need its own design pass on chain semantics,
  REF widening, and the `MAX_RECORD_SPAN` story under a different
  data structure. This plan only produces the empirical input to
  that decision.

## Reports to be produced

| Scope | Report path | Trigger |
|---|---|---|
| C | `ia/reports/implementations/pileup_freebayes_bench_c_<date>.md` | Always |
| A | `ia/reports/implementations/pileup_freebayes_bench_a_<date>.md` | Only if C says proceed |
| B | `ia/reports/implementations/pileup_freebayes_bench_b_<date>.md` | Only if A says proceed |

Each report follows the structure of
[`pileup_lazy_cigar_2026-05-07.md`](../reports/implementations/pileup_lazy_cigar_2026-05-07.md):
methodology, results table, what the data says, follow-ups. The
freebayes review's `A3` and §6 get a one-paragraph update at the
end of each scope citing that scope's report.
