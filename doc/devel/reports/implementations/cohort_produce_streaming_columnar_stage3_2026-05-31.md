# Stage 3 — DUST-ahead queue (2026-05-31)

Implementation report for Stage 3 of
[cohort_produce_streaming_columnar.md](../../../doc/devel/implementation_plans/cohort_produce_streaming_columnar.md):
take the serial sdust DUST mask off the producer's critical path. After
Stage 2 (the memory fix) the single-threaded producer was the remaining
wall floor, and its largest serial piece was the inline DUST pass
(~10 s at N=26, constant in N). A background thread now precomputes the
DUST masks for the covered intervals and feeds them to the producer
through a bounded queue; the producer slices its block's span out of the
interval mask instead of DUSTing inline.

## Plan (as executed)

1. **Types + helpers** — `DustChromPlan` (the covered-interval plan
   shared by producer and DUST thread), `IntervalDustMask` (queue
   message), and the free functions `build_dust_plans` /
   `dust_mask_for_interval` / `spawn_dust_ahead` / `slice_interval_mask`.
2. **Producer rewire** — `BlockIterator` plans the covered intervals
   once at construction, spawns the DUST-ahead thread (gated on
   `!no_complexity_filter`), and in `produce_block` slices the current
   interval mask rather than calling `sdust_mask_for_span`.
3. **Lifecycle** — the producer's `ChromCursor` loses its
   `ManualEvictChromRefFetcher`; the DUST thread owns one per chromosome.
   `Drop` joins the thread after releasing the receiver.
4. **Verify** — tests, clippy, byte-identity + perf gate.

## Domain intent

DUST depends only on the reference sequence + the region — nothing about
the PSP data or the fold — and the covered intervals are known up front
from the in-memory block indices. So the masks precompute trivially off
the producer thread. The contract is **byte-identical VCF** to the Stage 2
branch: the mask fed to `partition_window` for a block's span
`[range.start, safe_end)` must equal the prior inline
`sdust_mask_for_span` over that span.

## Why slicing is byte-identical

`sdust_mask_for_span`'s per-position verdict is documented byte-identical
to a whole-contig scan at every position inside the requested span (for
any span past its `>= window`-unmasked reset barrier; the `MIN_DUST_HALO`
+ barrier rule is unchanged). So:

- Precomputing the mask over a whole covered interval and **slicing** it
  to a block's sub-span gives the same per-position verdict as the inline
  per-block call — `slice_interval_mask` clips each run to the span,
  re-emitting a run that straddles a block boundary into both blocks.
- The DUST thread computes each interval in `DUST_AHEAD_SUBSPAN` (1 Mb)
  sub-spans and **coalesces** masked runs at the boundaries, which
  reproduces a single whole-interval scan exactly while bounding the
  resident reference buffer to ~one sub-span + halo (not a whole
  chromosome arm). `evict_before` keeps the buffer compact.

## Assumptions / decisions

- **Lockstep by shared plan, not recomputation.** The producer and the
  DUST thread iterate the *same* `Vec<DustChromPlan>` (built once via
  `build_dust_plans`, the same skip rules — drop zero-length and
  uncovered contigs — the producer's chromosome advance applies). The Nth
  `recv()` is therefore the Nth covered interval the producer enters; a
  `debug_assert` checks the `(chrom_id, interval)` match each enter.
- **Per-interval granularity, sliced by the producer.** Block cut points
  are dynamic (they depend on the fold), so the thread cannot precompute
  per-block. It DUSTs the *covered intervals* and the producer slices —
  exactly as the plan prescribes.
- **Bounded queue.** `sync_channel(DUST_AHEAD_QUEUE_CAP = 4)` — the thread
  can run at most a few interval masks ahead; masks are cheap and the
  resident reference is one sub-span, so memory stays flat.
- **`no_complexity_filter`** spawns no thread; `dust_rx` is `None` and the
  producer leaves the mask empty (matches the inline `dust_fetcher: None`
  path).
- **Error propagation.** A fetch/open failure in the thread is queued as a
  `ChunkDriverError` and surfaced by the producer when it reaches that
  interval; an unexpected early disconnect surfaces as the new
  `DustAheadGone`.

## Changes made

All in [driver.rs](../../../src/var_calling/cohort_block/driver.rs):

- New `DustChromPlan` / `IntervalDustMask` types; `DUST_AHEAD_QUEUE_CAP`,
  `DUST_AHEAD_SUBSPAN` consts.
- `build_dust_plans` (covered-interval plans from block indices),
  `dust_mask_for_interval` (sub-span DUST + coalesce, sub-span size
  parameterised for testability), `spawn_dust_ahead` (the background
  thread + bounded channel), `slice_interval_mask` (per-block slice with a
  forward sweep cursor).
- `ChromCursor`: dropped `dust_fetcher` and the now-unused `chrom_name` /
  `chrom_length` fields.
- `BlockIterator`: dropped `chromosomes` / `next_chrom_idx`; added
  `chrom_plans`, `next_plan_idx`, `dust_rx`, `dust_handle`,
  `cur_interval_mask`, `mask_sweep`. `new` plans intervals + spawns the
  thread; `advance_chrom` walks the plans; `enter_current_interval` (now
  fallible) `recv`s the interval mask; `produce_block` slices it.
- `Drop for BlockIterator` releases the receiver then joins the thread.
- New `ChunkDriverError::DustAheadGone { chrom_id }`.

## Tests added

In `driver.rs` (`cargo test --lib`):

- `slice_interval_mask_clips_and_advances_across_blocks` — consecutive
  block spans each get their runs; the sweep advances monotonically.
- `slice_interval_mask_reemits_run_straddling_block_boundary` — a run
  across a block cut is clipped into both blocks (the byte-identity case).
- `slice_interval_mask_empty_when_span_misses_every_run`.
- `dust_mask_for_interval_matches_whole_interval_scan` — on a temp FASTA
  (high-complexity flanks around a poly-A tract that crosses the
  boundaries), the sub-span-chunked + coalesced mask equals a single
  whole-interval `sdust_mask_for_span`, for sub-span sizes 7/13/41/whole.

## Validation results

- `cargo test --lib`: **1056 passed**, 0 failed, 1 ignored (1052 + 4 new);
  `--tests` integration suites all green (cohort_cli 20/20, +5 others).
  `cargo clippy --lib`: clean (one pre-existing, unrelated
  `needless_range_loop` in `dust_filter.rs`). The `psp_writer_perf`
  criterion bench panics under `--all-targets` — pre-existing, untouched,
  unrelated to this change.
- **Byte-identical** (header-stripped: drop `^##`, then md5) to the Stage 2
  branch baselines, **serial and 8 threads**, on the real tomato cohort:

  | N  | md5 (serial = 8 threads) | baseline |
  |----|--------------------------|----------|
  | 8  | `8f117ac7…` | `8f117ac7fef9f83d240644b1d0def3c8` ✓ |
  | 26 | `a3164c38…` | `a3164c386c1b21cf94c6adc5cf4c60d3` ✓ |

- **Wall + RSS (real tomato cohort, 8 threads, `perf_main_vs_branch.py`):**

  | N  | main wall | branch wall | main RSS | branch RSS |
  |----|----------:|------------:|---------:|-----------:|
  | 8  | 6.0 s | 9.1 s | 147 MB | 83 MB |
  | 26 | 9.1 s | 9.8 s | 395 MB | 294 MB |

  At N=26/8 threads the branch is **9.8 s vs main 9.1 s (1.06×)** — down
  from Stage 2's 1.76× (15.2 vs 8.7 s): the wall gap essentially closed,
  DUST having left the producer's critical path. Peak RSS dropped further
  (Stage 2 550 MB → 294 MB, now **below** main's 395 MB) because the
  producer no longer holds a reference fetcher and the DUST thread's
  sub-span buffer is bounded. The branch-vs-main md5 difference is the
  known, pre-existing filter-order divergence — orthogonal to this work;
  the byte-identity gate above is against the branch's own output.

## Tradeoffs and follow-ups

- **Low N still produce-bound.** At N=8 the branch (9.1 s) trails main
  (6.0 s): produce dominates and the single producer thread is the floor.
  Closing that needs **parallel producers over independent covered
  intervals** — the next lever, noted in the plan and out of scope here.
- The DUST thread shuts down on receiver drop (`Drop` joins it); worst
  case it finishes the in-flight interval before its next `send` fails —
  bounded, sub-second.
- `load_chunk_from_iters` remains the streaming loader's batch oracle
  (unchanged from Stage 2).
