# Stage 2 — streaming fold+compact producer (2026-05-31)

Implementation report for Stage 2 of
[cohort_produce_streaming_columnar.md](../../../doc/devel/implementation_plans/cohort_produce_streaming_columnar.md):
replace the batch "grow the span to the variant target, then compact"
loader with a **streaming fold+compact** producer that holds only the
open variant group, fixing the cohort caller's memory footprint.

## Plan (as executed)

Sub-staged to keep each step gated:

- **2.0** — flip `BlockColumnReader` / `ColumnSpanReader` to *own* their
  `PspReader` (the producer needs `Vec<ColumnSpanReader>` without a
  self-referential borrow).
- **2.1** — `SpanColumnSource` trait + `StreamingBlockLoader`: read the
  cohort forward to a shared watermark, fold, and compact *closed*
  variant groups each round, keeping only the open-group straddle.
  Validated by equivalence tests against the batch loader.
- **2.2a** — rewire `BlockIterator` onto the sources + the streaming
  loader; `produce_block` via `fill_block`. The end-to-end memory fix.
- **2.2b** — delete the now-dead batch safe-cut/carryover machinery.

## Assumptions / decisions

- **Watermark = min over sources of `peek_next_span`** (one block per
  round). The per-block memory floor is "one decoded block × N" — held
  in the readers regardless — so aligning reads to block edges is as
  lean as a smaller fixed step and avoids redundant re-folds.
- **`load_chunk_from_iters` is retained** as the streaming loader's
  independent batch oracle for the equivalence tests, even though it is
  off the production path. The user's "delete old code" applies to the
  *producer* complexity (safe-cut / carryover / NoSafeGap), which is
  gone; an independent oracle earns its keep.
- **Block cut = the open group's start** (or watermark+1 when all
  closed) — a clean boundary by construction, so the VCF is independent
  of where blocks fall (byte-identity).
- `target_variants` is `.max(1)` at the producer so a `0` (disabled)
  target still cuts at the first closed group rather than emitting empty.

## Changes made

- **psp** ([reader.rs](../../../src/psp/reader.rs)): `BlockColumnReader`
  now owns its `PspReader` (`into_column_blocks` / `into_reader`).
- **[loader.rs](../../../src/var_calling/cohort_block/loader.rs)**:
  `SpanColumnSource` trait; `StreamingBlockLoader` (`fill_block` →
  kept-variant count, `compact_closed_prefix`, `find_block_cut`).
  Fixed a latent `is_kept` stale-`true` bug in the shared fold
  (`rebuild_position_union_and_is_kept`) exposed by scratch reuse.
- **[column_span_reader.rs](../../../src/var_calling/cohort_block/column_span_reader.rs)**:
  own the reader; `detached`/`reset`/`block_index`; `SpanColumnSource`
  impl.
- **[driver.rs](../../../src/var_calling/cohort_block/driver.rs)**:
  `BlockIterator` holds `Vec<ColumnSpanReader>` + `StreamingBlockLoader`;
  `enter_current_interval` resets the sources per covered interval;
  `next_block` loops over intervals; `produce_block` calls `fill_block`.
  New `ChunkDriverError::StreamRead`.
- **2.2b deletions**: `chunk_boundaries.rs` (the safe-cut module),
  `ChunkDriverError::FinaliseChunkBoundaries`, the
  `run_finalise_chunk_boundaries` test helper, `MAX_CHUNK_SPAN_GROWTH`.

## Tests added

- `loader.rs`: a `Vec<PileupRecord>`-backed `SpanColumnSource` fixture +
  equivalence tests asserting the streaming loader's kept records (across
  all emitted blocks) equal the batch loader's, for **aligned**,
  **cross-sample misaligned**, and **clamped** regions, at several
  variant targets.

## Validation results

- `cargo test --lib`: **1052 passed, 0 failed, 1 ignored**;
  `cohort_cli_integration`: 20/20. `cargo clippy --lib`: clean (one
  pre-existing, unrelated `needless_range_loop` in `dust_filter.rs`).
- **End-to-end byte-identical** (header-stripped md5) to the prior
  driver at N=8 and N=26, serial and at 8 threads.
- **Memory (real 26-sample tomato cohort, 8 threads, peak RSS):**

  | N | before (batch) | after (streaming) | main |
  |---|---:|---:|---:|
  | 8 | 765 MB | 397 MB | 154 MB |
  | 26 | 3963 MB | 550 MB | 397 MB |

  The N=26 gap to `main` closes from ~9× to **1.4×**.
- **Wall (N=26, 8 threads):** 19.3 s → 15.2 s (bonus: no re-fold on
  span growth, fewer allocations). Still 1.76× `main` — the *producer*
  is now the serial floor (single producer thread; serial DUST), which
  Stage 3 + future produce parallelism address. Memory was the goal and
  is met.

## Tradeoffs and follow-ups

- The remaining wall gap to `main` is the single-threaded producer.
  **Stage 3** (DUST-ahead queue) removes the ~10 s serial DUST from its
  critical path; beyond that, parallel producers over independent
  covered intervals are the lever (future, noted in the plan).
- The branch's VCF still differs from `main` by the known, pre-existing
  filter-order divergence — orthogonal to this work.
- `load_chunk_from_iters` + its tests remain as the batch oracle.
