# Pileup walker — lazy CIGAR migration

**Date:** 2026-05-07
**Implementer:** Claude (rust-feature-implementation skill)
**Plan:** [`ia/feature_implementation_plans/pileup_lazy_cigar.md`](../../feature_implementation_plans/pileup_lazy_cigar.md)
**Motivating review:** [`ia/reviews/pileup_samtools_comparison_2026-05-07.md`](../../reviews/pileup_samtools_comparison_2026-05-07.md)

## What was built

The walker no longer materialises a `Vec<ReadEvent>` per active
read on admission and no longer clones it once per walker step.
It now drives a `CigarCursor` (in
[`src/per_sample_caller/pileup/cigar_cursor.rs`](../../../src/per_sample_caller/pileup/cigar_cursor.rs))
that holds a small per-CIGAR-op offset table (~50 entries even
for long reads) and produces events on demand via two
intentionally-distinct query shapes:

- `events_at(walker_pos, read)` — events anchored exactly at
  `walker_pos`. Walker uses this once per step.
- `events_overlapping(lo, hi, read)` — events whose footprint
  intersects `[lo, hi)`. Open-record fold uses this to assemble
  haplotype strings for `apply_events_to_ref`.

`ReadContribution` lost its `full_window_events: Vec<ReadEvent>`
field; it now carries `read_id` and a `bq_zero_in_window: bool`
flag set by the mate-overlap resolver for match-only losers. The
fold pulls the live read by id from `&ActiveSet` and queries the
cursor.

`decompose` is preserved behind `#[cfg(test)]` as the parity
oracle the cursor's tests assert against. The `event_cursor`
residual-events invariant on expiry is gone — single-fold per
`(record, read)` was already enforced by
`OpenPileupRecord.folded_reads`.

## Validation

Run inside the dev container.

| Command | Result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | pileup scope clean (9 pre-existing diagnostics in `variant_grouping.rs` / `decompression_pool.rs` remain — out of scope) |
| `cargo test --lib pileup` | 60 passed, 0 failed (was 55; +5 cursor tests) |
| `cargo test --lib` | 150 passed, 0 failed |
| `cargo bench --bench pileup_walker_scaling` | see below |

## Performance

Same fixture as the baseline measurement
([`benches/pileup_walker_scaling.rs`](../../../benches/pileup_walker_scaling.rs)):
50 kb window, 30× coverage, pure-Match reads. Both runs use the
release profile inside the container.

| read length | eager (baseline 2026-05-07) | cursor (this commit) | speedup |
|---|---|---|---|
|   150 |   1.42 s |   480 ms |   **2.96×** |
|   500 |   3.02 s |   483 ms |   **6.25×** |
|  1500 |   8.70 s |   457 ms |   **19×**   |
|  5000 |  70.9 s  |   399 ms |   **178×**  |

Per-position time on the cursor walker is essentially flat at
~8–10 µs across the whole L sweep, against the eager walker's
28 → 1418 µs ramp. L=5000 is now slightly faster than L=150 in
absolute terms because at fixed (window, coverage) longer reads
mean fewer admit/expire cycles.

Plan targets:

- L=150 within ±10% of baseline — **3× faster than baseline**
  (no regression at typical Illumina; production runs benefit
  too).
- L=1500 ≥3× faster — **19× faster**.
- L=5000 ≥10× faster — **178× faster**.

All clear by an order of magnitude or more.

## Multi-op baseline (added later same day)

Recorded after the early-break landed (commit `401fe94`), as a
forward-looking guide for future cursor changes (the
binary-search-on-offsets idea queued in the plan). Same span /
coverage as the main bench; CIGAR is `[M(50), I(2)] × cycles +
M(remainder)`, so each read carries many small ops scaling with
`L`.

| L | CIGAR ops | time | per-position |
|---|---|---|---|
|  150 |   5 | 436 ms |  8.7 µs |
|  500 |  19 | 551 ms | 11.0 µs |
| 1500 |  59 | 602 ms | 12.0 µs |
| 5000 | 199 | 803 ms | 16.1 µs |

For comparison, the single-op fixture reported 470 / 452 / 456 /
400 ms on the same run — per-position cost stays ~9 µs across L
because there's only one op to walk per query. The multi-op
fixture climbs from 8.7 µs to 16.1 µs as ops grow 5 → 199, which
is the cost the early-break already trims and the
binary-search-on-offsets would trim further.

The fixture uses **insertions** rather than deletions: an early
draft with `[M(50), D(2)]` cycles tripped the walker's
`MAX_RECORD_SPAN = 5000` cap because high-coverage overlapping
deletions chain-widen the same open record. Insertions have
footprint span 1 (just the anchor base), so they never widen
records. The walker safety cap correctly fired on the first
draft — the walker did the right thing; the fixture was the
wrong workload.

## Auto-selecting binary search (added later same day)

Followed up the multi-op-baseline measurement with a binary-search
implementation of the cursor's op walk. The cursor now carries a
`mode: CursorMode` field set once at construction:
`CursorMode::auto_select(cigar.len())` returns `BinarySearch` when
the CIGAR has more than `BINARY_SEARCH_OP_THRESHOLD = 16` ops and
`Linear` otherwise. Public methods dispatch internally — callers
see no API change.

Two private impls per public method:

- **Linear** (`events_at_linear`, `events_overlapping_linear`):
  walk every op left-to-right with the early break. Per-op match
  arms inlined manually rather than delegated to a helper — the
  bench showed `#[inline]` / `#[inline(always)]` on a shared
  helper still left a ~15 % function-call tax on hot Linear
  workloads. Manual inlining duplicates ~80 lines but recovers
  full baseline performance.
- **Binary search** (`events_at_binary`,
  `events_overlapping_binary`): use `offsets.partition_point` to
  find the first relevant op and walk forward from there. The
  shared helpers `emit_event_for_op_*` (annotated `#[inline]`)
  back these — per-call setup dominates either way, so the
  helper-call cost is invisible.

The CIGAR walk's structure makes a `max_deletion_len` field
unnecessary: a deletion at op `j` has footprint ending at exactly
`offsets[j+1].ref_pos`, so the lower bound `last op with
ref_pos ≤ lo` is safe even for deletions anchored before `lo` —
their footprints all end at or before `lo` by construction.

Tests parameterise the parity oracle over both modes, so each
implementation is checked against `decompose`. An additional
test (`auto_select_picks_mode_by_op_count`) pins the threshold
boundary at 16 ops.

### Measured impact

Comparing the auto-select cursor against a fresh baseline run of
commit `5886d0c` (linear + early-break, no binary search). Same
benchmark fixtures, same machine, same span / coverage as the
prior tables.

| Workload | Mode picked | Baseline | Cursor | Δ |
|---|---|---|---|---|
| single L=150  | Linear | 453 ms | 449 ms | ~0 |
| single L=500  | Linear | 448 ms | 463 ms | +3 % |
| single L=1500 | Linear | 427 ms | 458 ms | +7 % |
| single L=5000 | Linear | 473 ms | 431 ms | -9 % |
| multi L=150  (5 ops)   | Linear        | 537 ms | 497 ms | **-7 %** |
| multi L=500  (19 ops)  | BinarySearch  | 460 ms | 512 ms | +11 % |
| **multi L=1500** (59 ops)  | BinarySearch  | 665 ms | 481 ms | **-28 %** |
| **multi L=5000** (199 ops) | BinarySearch  | 948 ms | 547 ms | **-42 %** |

The big wins are exactly where binary search was supposed to
help: long, indel-rich CIGARs at L=1500 and L=5000 both pick the
binary path and gain 28–42 %. Linear-mode workloads sit within
single-digit percent of the baseline (one is +7, one is -9 — both
inside the run-to-run noise band of ~15 % we measured on
unchanged code).

The single +11 % at multi L=500 looks like the only remaining
"regression"; it's against a baseline (460 ms) that landed at the
fast tail of its distribution while the cursor's number (512 ms)
is mid-distribution. Across the broader dataset multi L=500 is
within noise of baseline either way.

### A note on noise and baselines

The earlier sections of this report compared against the
single-point numbers from commit `5886d0c`'s original run. A
fresh re-run of that same code under current system load showed
~15–25 % run-to-run variance on the multi-op fixture (e.g. multi
L=150 was 436 ms originally and 537 ms on the fresh re-run —
same code). Several "regressions" we chased down during the
binary-search work turned out to be artefacts of that fast-tail
baseline. The measurements above use the fresh-baseline run as
the comparison point, which is the honest one.

For future cursor work, the discipline is: take **two runs at
each codepoint** before drawing conclusions, and treat any delta
under ~15 % on a single-run pair as noise.

## Trade-offs and follow-ups

- **Public API unchanged.** `PreparedRead`, `PileupRecord`,
  `AlleleObservation`, `RefBaseFetcher`, `run`, and `RunSummary`
  all keep their signatures. The lazy machinery is private to
  the `pileup` module.
- **Memory footprint dropped sharply.** Per-active-read state
  went from `O(span)` (the eager event vec) to
  `O(num_cigar_ops)` (the cursor's offset table). The
  per-walker-step contribution clones disappeared entirely.
- **`MAX_RECORD_SPAN = 5000` ceiling not raised.** Long-read
  *input* support (CRAM decoder limits, BAQ scaling, span
  ceiling) remains out of scope; this work removes the walker
  itself as a bottleneck.
- **Multi-op (insertion-heavy) benchmark fixture added** in a
  follow-up commit, then used to validate the auto-select
  binary-search work above. The original "indel-heavy" stretch
  goal is partially covered (insertions exercise op count;
  deletions still don't because high-coverage deletion fixtures
  trip `MAX_RECORD_SPAN`). A deletion-heavy variant remains a
  follow-up if/when long-read CRAM input lands.
- **`decompose` retained as test-only oracle.** Could be deleted
  in a follow-up if the cursor's own unit tests grow to cover
  every case `decompose`'s tests do, but parity vs. a
  reference impl is cheap insurance.
