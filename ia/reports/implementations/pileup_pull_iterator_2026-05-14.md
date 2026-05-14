# Pileup walker: push-channel → pull-iterator (implementation report)

Implementation date: 2026-05-14. Plan:
[ia/feature_implementation_plans/pileup_pull_iterator.md](../../feature_implementation_plans/pileup_pull_iterator.md).

## Plan

Replace the pileup walker's push-based output (`SyncSender<PileupRecord>`,
caller-spawned consumer thread, bounded channel for backpressure) with a
pull-shaped `Iterator<Item = Result<PileupRecord, WalkerError>>`. The
emitted record stream is identical: same ordering, same lifecycle-mark
attachment rule, same error variants minus the now-unreachable
`ChannelClosed`.

Same-day-conversation rationale: the bounded channel (capacity 64) sat
between walker and writer; per-record writer work is cheap (varint
accumulation), only the per-block work is heavy (zstd L9 on 16 MiB
blocks). During the heavy step the walker filled 64 slots and stalled —
overlap was on the order of 0.01 % of a block. No multi-sample driver
exists yet, so cross-sample parallelism will come from rayon-over-samples
when it's built, leaving within-sample threading redundant.

## Assumptions / silent choices

- **Iterator error semantics: terminate on first error.** Once `next()`
  returns `Err`, all subsequent calls return `None`. Matches the
  push-shape's behaviour where `run` returned `Err` and stopped emitting
  at the same moment. Pinned by the new test
  `walker_iterator_returns_none_after_yielding_error`
  ([src/per_sample_caller/pileup/tests.rs:1037-1059](../../../src/per_sample_caller/pileup/tests.rs#L1037-L1059)).
- **Blanket `&T: RefSeqFetcher` forwarding impl** added to keep existing
  call sites (`run(reads, &fetcher, &config)`) working unchanged. The
  trait only takes `&self` so the forwarding is zero-cost.
- **`pending` is `VecDeque<PileupRecord>`** rather than `Vec` with a
  cursor. Per-tick batches are small; `pop_front` is O(1) and preserves
  the tick-internal ordering the lifecycle-mark stamp rule depends on.
- **`stamp_lifecycle_marks` helper deleted.** With the inversion, the
  walker now stamps marks directly via
  `out.get_mut(batch_start_idx)` inside `close_aged_records_into` and
  `flush_chromosome_into` ([src/per_sample_caller/pileup/walker.rs:497-526](../../../src/per_sample_caller/pileup/walker.rs#L497-L526)
  and lines 528-559). The helper had no remaining callers.

## Changes made

### `src/per_sample_caller/pileup/walker.rs`

- Replaced the `pub fn run<I, F>(reads, ref_fetcher, tx, config) ->
  Result<RunSummary, WalkerError>` push-API with:
  - `pub struct PileupWalker<I, F>` holding the peekable read iterator,
    the ref fetcher (by value, with a blanket impl letting callers pass
    `&F`), the existing `WalkerState`, a `VecDeque<PileupRecord>`
    pending buffer, and a `done` flag.
  - `PileupWalker::new(reads, ref_fetcher, config)` constructor.
  - `PileupWalker::summary(&self) -> RunSummary`.
  - `impl Iterator for PileupWalker { type Item = Result<PileupRecord, WalkerError>; … }`.
  - A new private `fill_pending(&mut self)` method that hosts the outer
    walker loop from the old `run`. It loops until at least one record
    is in `pending`, or until end-of-input (in which case it flushes
    the final chromosome and sets `done = true`).
  - `pub fn run<R, F>(reads, ref_fetcher, config) -> PileupWalker<R::IntoIter, F>`
    as a thin `IntoIterator` constructor that matches the previous
    entry-point name.
- Renamed `close_aged_records(tx)` → `close_aged_records_into(out: &mut VecDeque<PileupRecord>)`
  and `flush_chromosome(tx)` → `flush_chromosome_into(out: &mut VecDeque<PileupRecord>)`.
  Both now push finalised records onto `out` in emission order and
  stamp lifecycle marks onto the first record of this tick's batch via
  `out.get_mut(batch_start)`.

### `src/per_sample_caller/pileup/mod.rs`

- Re-exported `PileupWalker` alongside `run` and `RunSummary`.
- Removed the `DEFAULT_OUTPUT_CHANNEL_CAPACITY` constant (no longer
  applicable).
- Added a blanket `impl<T: RefSeqFetcher + ?Sized> RefSeqFetcher for &T`
  forwarding impl so callers can pass either an owned fetcher or `&F`.

### `src/per_sample_caller/pileup/errors.rs`

- Deleted the `WalkerError::ChannelClosed { context: String }` variant.
  No production path can produce it after the rewrite.

### `src/per_sample_caller/pileup/open_record.rs`

- Deleted the `pub(super) fn stamp_lifecycle_marks(...)` helper. The
  walker now stamps marks inline; there are no remaining callers.

### `src/per_sample_caller/pileup/tests.rs`

- Rewrote `drive_walker_with_config` to drive the iterator directly
  ([src/per_sample_caller/pileup/tests.rs:113-138](../../../src/per_sample_caller/pileup/tests.rs#L113-L138)).
  No more `mpsc::sync_channel`, no collector thread, no `drop(tx)` dance.
- Replaced the four channel-allocating error-variant tests with a small
  `first_walker_error` helper that consumes the iterator until it
  surfaces an error.
- Deleted `run_returns_channel_closed_when_receiver_dropped_mid_stream`
  (variant gone).
- Added `walker_iterator_returns_none_after_yielding_error` pinning the
  terminal-on-first-error contract.
- Updated `out_of_order_input_is_a_hard_error` and
  `chromosome_id_regression_is_a_hard_error` to use `first_walker_error`.

### `examples/dhat_pileup.rs`

- Direct iteration: a single `for item in run(reads, &ref_fetcher,
  &WalkerConfig::default())` loop counts emitted records. No channel,
  no collector thread.

### `benches/pileup_walker_scaling.rs`

- `setup_walker` removed (no per-iter channel + thread setup needed).
- `drive_walker(reads, ref_fetcher)` now iterates the walker to
  exhaustion inside the timed body. The previous explanatory comment
  about the per-iter channel + thread spawn (round-2 L2 in
  `ia/reviews/perf_pileup_2026-05-12.md`) is replaced with a note that
  the rewrite removed that overhead.

## Tests added/updated

- **Updated** `drive_walker_with_config` — single helper now drives the
  iterator; all 60+ behaviour tests that depend on it are exercised
  unchanged.
- **Updated** `out_of_order_input_is_a_hard_error`,
  `chromosome_id_regression_is_a_hard_error`,
  `zero_ref_span_input_is_a_hard_error`,
  `open_record_widening_past_max_record_span_errors`, three
  `admit_rejects_*` tests, and
  `fasta_fetch_failure_propagates_as_walker_error_fasta` — all use the
  new `first_walker_error` helper.
- **Deleted** `run_returns_channel_closed_when_receiver_dropped_mid_stream`
  (variant no longer exists).
- **Added** `walker_iterator_returns_none_after_yielding_error` — pins
  the new terminal-on-first-error contract: after `next()` returns
  `Err`, every subsequent call returns `None`.

## Validation results

All commands run inside the dev container (`./scripts/dev.sh`):

- `cargo fmt` — applied (one wrap on the imports block in `walker.rs`
  and one on the `flush_chromosome_into` signature; both formatting-only).
- `cargo clippy --all-targets --all-features -- -D warnings` —
  **clean**, no warnings.
- `cargo test --tests --lib --examples --all-features` — **all green**:
  - 474 unit tests (lib) — pass.
  - 25 + 26 + 17 + 8 + 12 = 88 integration tests across `tests/` — pass.
  - 0 doc tests and 0 example "tests" (examples have no `#[test]`).
- `cargo check --benches --all-features` — clean.

`cargo test --all-targets` itself trips an unrelated pre-existing
panic in `benches/gvcf_perf.rs:10` (`Benchmark file not found:
"/home/jose/analyses/g2psol/source_data/TS.vcf.gz"` — a missing
local fixture, not a regression from this change). Verified by
re-running the focused command above.

## Tradeoffs and follow-ups

- **Channel overhead removed.** The bench's timed body no longer pays
  the per-iter `mpsc::sync_channel` allocation + `thread::spawn` +
  `drop(tx) + join()` teardown. Numbers will shift slightly downward
  on the next bench run; this is the desired direction and was one of
  the design motivations.
- **Out-of-scope follow-ups** (unchanged from the plan):
  - A tiny prefetch adapter (1–2-slot buffer) if profiling later shows
    jitter matters. Do not pre-build.
  - Stage 1 multi-sample driver (rayon over BAMs). The natural place
    for cross-sample parallelism after this rewrite.
- **API surface** is strictly narrower now: one fewer public constant
  (`DEFAULT_OUTPUT_CHANNEL_CAPACITY`), one fewer error variant
  (`ChannelClosed`), one fewer free function used externally
  (`stamp_lifecycle_marks` was `pub(super)` so no external impact).
  The new public item is `PileupWalker`.
