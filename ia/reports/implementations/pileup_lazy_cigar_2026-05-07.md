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
- **Indel-heavy benchmark variant deferred.** The plan's commit 3
  stretch goal of an indel-heavy fixture was not added — the
  pure-Match numbers were already a clean win across the
  decision-relevant range. Worth adding alongside any future
  long-read work to catch regressions in the indel path
  specifically.
- **`decompose` retained as test-only oracle.** Could be deleted
  in a follow-up if the cursor's own unit tests grow to cover
  every case `decompose`'s tests do, but parity vs. a
  reference impl is cheap insurance.
