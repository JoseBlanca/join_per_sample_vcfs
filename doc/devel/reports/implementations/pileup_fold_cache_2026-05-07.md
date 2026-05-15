# Pileup walker — fold-cache investigation (negative result)

**Date:** 2026-05-07
**Investigator:** Claude (research-only session)
**Motivation:** post-implementation hypothesis from the lazy-CIGAR
work (`ia/reports/implementations/pileup_lazy_cigar_2026-05-07.md`)
that the open-record fold redundantly recomputes
`events_overlapping → apply_events_to_ref → find_or_create_allele_index`
across walker steps when neither the record's footprint nor the
contributor's relevant state has changed.

**Verdict:** **Don't ship a cache.** On the multi-op L=5000 bench
fixture (and every other fixture in
`benches/pileup_walker_scaling.rs`), **0 of ~1.5 M fold iterations
are re-folds**. The hypothesised cache has zero hit rate on the
workloads the bench measures. The cache mechanism *would* help on
deletion-heavy workloads (we measured 80 % redundancy on a
synthetic deletion fixture), but those workloads aren't what the
walker is being tuned for and trip `MAX_RECORD_SPAN` at the
bench's coverage anyway.

## Methodology

The fold loop in [`open_record::process_position`][1] iterates over
`(affected record, contributor)` pairs, computes the contributor's
new haplotype/scalar contribution against the record's footprint,
and either replaces an existing entry in `folded_reads` (a re-fold)
or installs a fresh one. A re-fold is *redundant* when the new
`(allele_index, FiveScalars)` pair is byte-equal to the prior one
— the subtract-then-add round-trip cancels out, but the cursor
walk and `apply_events_to_ref` both still ran.

[1]: ../../../src/per_sample_caller/pileup/open_record.rs

I added a temporary `fold_metrics` module to `open_record.rs`
(test-only atomics, gated behind `cfg(test)`) plus an `#[ignore]`d
harness in `pileup/tests.rs` that ran the same fixtures
`benches/pileup_walker_scaling.rs` uses, plus a deletion variant.
The harness counted, per fold iteration:

- `iterations` — every (affected record, contributor) pass.
- `nonempty_window` — those whose `events_overlapping` returned a
  non-empty window (the early-return short-circuit *did* fire if
  this is below `iterations`).
- `refolds` — those that found a prior `FoldedReadState`.
- `redundant` — re-folds where the new (allele seq, contribution)
  equals the prior.
- `overlapping_events_returned` — sum of events emitted by
  `events_overlapping` in the fold (proxy for cursor/apply work).

The instrumentation has since been reverted now that the result is
recorded here. Reproducing requires re-adding the counters; the
table below is the answer they gave on the run that motivated the
revert.

## Results

`span = 50_000`, `coverage = 30` (matching the bench).

| Fixture                          | walker_ms | iters     | non-empty (%) | re-folds (%)   | redundant (% of re-folds) | redundant (% of non-empty) |
|----------------------------------|----------:|----------:|--------------:|---------------:|--------------------------:|---------------------------:|
| single L=150                     |       512 | 1 500 000 |        100 %  |          0 %   |                       n/a |                       0 %  |
| single L=5000                    |       482 | 1 500 000 |        100 %  |          0 %   |                       n/a |                       0 %  |
| multi  L=150  (5 ops)            |       538 | 1 500 000 |        100 %  |          0 %   |                       n/a |                       0 %  |
| multi  L=1500 (59 ops)           |       537 | 1 500 000 |        100 %  |          0 %   |                       n/a |                       0 %  |
| **multi L=5000 (199 ops)**       |   **551** | 1 500 000 |        100 %  |        **0 %** |                       n/a |                       0 %  |
| del L=150 cyc=50 del=2 cov=4     |        98 |   194 618 |        100 %  |          8.2 % |                  100.0 %  |                     8.2 %  |
| del L=150 cyc=20 del=8 cov=4     |        75 |   146 630 |        100 %  |         80.6 % |                   99.8 %  |                    80.5 %  |

(`coverage = 4` on the deletion rows so chained deletions don't
widen one record past `MAX_RECORD_SPAN`. The bench keeps `cov = 30`
which is why it only ships insertion-multi-op fixtures, not
deletion-multi-op.)

## Why every bench fold is fresh

The closure rule is `pos + ref_span ≤ walker_pos`. A read that has a
`Match` at walker_pos `k` creates or updates the record anchored at
`k` with `ref_span = 1`; that record closes the moment the walker
advances to `k + 1`. So a pure-Match read folds into 150 different
records (one per position), each exactly once.

Insertions don't widen records (footprint = 1, just the anchor
base). The multi-op bench fixture uses `[M(50), I(2)] × cycles`,
so its records remain 1-base-wide too. Same outcome: every fold
is fresh.

A re-fold *requires* a record's footprint to span multiple ref
bases, so the same read can be active at two consecutive walker
steps inside that footprint. Only deletions widen footprints in
practice — and the bench deliberately avoids them because at
30× coverage the chained-deletion fixture trips
`MAX_RECORD_SPAN = 5000`.

## Why "100 % of refolds are redundant" on the deletion fixture

For a given `(record, read)` pair within an open record's lifetime,
the inputs to the fold are immutable except for three things:

1. `record.ref_seq` (changes only on widening).
2. `bq_zero_in_window` (set when the read becomes a match-only
   mate-overlap loser).
3. `bq_override_at_walker_pos` (per-step BQ rewrite for S7
   match-only winners).

The deletion fixture has no mate overlap and no widening *after*
the first fold (each record opens with its full deletion footprint
on the first event), so all three are constant across re-folds.
Inputs unchanged ⇒ outputs unchanged ⇒ every re-fold is redundant.

## Why caching is still the wrong move

Two reasons:

1. **The bench fixture has 0 % re-folds.** The cache would never
   fire on the workload the user wants to optimise. Adding any
   cache-validity check (even a fast one — version counter, hash,
   or `prev == new`) costs cycles per iteration that aren't
   amortised because no skip ever happens.

2. **Real-data refold rates fall between the deletion fixtures and
   zero.** Illumina BAQ-adjusted reads have ~0.1 % indel-base rate;
   ~10 % of those are deletions; deletion footprints span ~2-10
   ref bases. So on real Illumina data we'd expect a re-fold rate
   in the low single digits — the `cyc=50 del=2 cov=4` row's 8.2 %
   is an upper bound for typical deletion density. At those rates,
   even a free cache would cap the win at ~5 %, well below the
   ≥ 15 % threshold the user set.

   For long-read deletion-heavy workloads (the `cyc=20 del=8` row,
   80 % re-fold), the cache *would* be a clear win — but long-read
   support isn't in scope yet (`MAX_RECORD_SPAN` ceiling not
   raised; CRAM decoder limits not addressed) and the bench
   fixture cannot represent it without lifting the cap. Revisit
   if/when long-read input lands.

## What the data says about the secondary "are events_overlapping
+ apply_events_to_ref the bottleneck?" question

The harness doesn't directly time those functions, but the call
counts give an upper bound. On multi-op L=5000:

- 1 500 000 fold iterations in 551 ms ⇒ ~367 ns/iteration walker-side.
- 1 529 700 events emitted ⇒ 1.02 events per fold on average.

A single `events_overlapping` call on a 199-op CIGAR via the
binary-search cursor is `partition_point` (~100 ns) plus a short
forward walk (~50 ns) plus the per-event push. With one event
emitted per call, that's roughly 150-200 ns — half the per-fold
budget. The other half lives in `apply_events_to_ref` (a small
allocation and a copy of ~100 bytes), the `subtract/add`
arithmetic, and the per-iteration overhead of the fold loop.

So *if* one were to optimise these calls further, the headroom is
plausible but bounded by the per-fold budget; we'd be fighting for
single-digit-percent wins against a ~15-25 % run-to-run noise floor.
Out of scope for this investigation; flagged as a follow-up if
future profiling shows the fold dominating real-data wall time.

## Bench compile fix (incidental)

`benches/pileup_walker_scaling.rs` was stale: `run` gained a
`&WalkerConfig` parameter (S5 — per-column depth caps) and the
bench wasn't updated. Fixed in this branch by importing
`WalkerConfig` and passing `&WalkerConfig::default()` to `run`.
The fix is part of the working tree changes left behind by this
investigation.

## Artefacts left in the tree

- `benches/pileup_walker_scaling.rs` — bench-compile fix only:
  `run` now requires `&WalkerConfig` (S5), and the bench wasn't
  updated to pass one.

The investigation's instrumentation (the `fold_metrics` module
in `open_record.rs` and the `fold_cache_metrics` harness in
`pileup/tests.rs`) was reverted on conclusion. Re-adding it from
this report is straightforward if/when the question is reopened
(e.g. when long-read input lands and `MAX_RECORD_SPAN` is raised).
