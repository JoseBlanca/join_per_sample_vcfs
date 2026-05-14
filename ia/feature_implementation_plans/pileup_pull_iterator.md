# Pileup walker: push-channel → pull-iterator

Proposal date: 2026-05-14.

## Domain intent

Replace the pileup walker's current push-based output (`SyncSender<PileupRecord>`
on a bounded channel, with the consumer running on a separate thread) with a
pull-based shape: the walker becomes an `Iterator<Item = Result<PileupRecord,
WalkerError>>` that the consumer (the `.psp` writer, today; Stage 1 callers
in general) drives directly by calling `next()`.

The behaviour the walker produces — exactly the same record stream, in exactly
the same order, with the same lifecycle-mark attachment rules and the same
error variants (minus the now-unreachable `ChannelClosed`) — does not change.
Only the **threading boundary and the API surface** change.

## Why now

See the discussion in the user-facing thread that led to this plan; the
short version:

- The bounded channel (capacity 64) sits between walker and writer; the
  expensive writer work is per-block zstd L9 on 16 MiB blocks (~530 k records
  per block on SNP-heavy data). Per-record work on the writer is cheap
  (varint accumulation). During a block flush the walker fills 64 slots and
  stalls — overlap window is on the order of 0.01 % of the block. The
  threading boundary is not buying material parallelism for this workload.
- No multi-sample driver exists yet. When one is added, the natural shape is
  rayon-over-samples, which provides compute/IO overlap at the sample level.
  Within-sample threading then becomes redundant.
- Pull-shape gives easier tests (`.collect()` in unit tests), simpler error
  propagation (no `RecvError`-vs-producer-error reconciliation), and one
  fewer thread per sample at cohort scale.

## Out of scope

- Building the multi-sample driver. That's a separate piece.
- Changing the psp writer's record-side API (`write_record(&PileupRecord)`)
  or its compression / block format.
- Touching the walker's internal logic (closure rule, mate overlap, per-column
  depth caps, slot lifecycle). All of that stays bit-for-bit identical.
- Performance optimization beyond what the structural change does on its own.

## API shape

```rust
// src/per_sample_caller/pileup/walker.rs

pub struct PileupWalker<I, F>
where
    I: Iterator<Item = PreparedRead>,
    F: RefSeqFetcher,
{ /* fields private */ }

impl<I, F> PileupWalker<I, F> { /* … */
    pub fn new(reads: I, ref_fetcher: F, config: &WalkerConfig) -> Self;
    /// Snapshot of the cumulative counters. Safe to call mid-stream; the
    /// final summary is the value observed after `next()` has returned
    /// `None`.
    pub fn summary(&self) -> RunSummary;
}

impl<I, F> Iterator for PileupWalker<I, F>
where
    I: Iterator<Item = PreparedRead>,
    F: RefSeqFetcher,
{
    type Item = Result<PileupRecord, WalkerError>;
}

/// Thin constructor that mirrors the previous `run` entry point's
/// `IntoIterator` ergonomics.
pub fn run<R, F>(
    reads: R,
    ref_fetcher: F,
    config: &WalkerConfig,
) -> PileupWalker<R::IntoIter, F>
where
    R: IntoIterator<Item = PreparedRead>,
    F: RefSeqFetcher;
```

To preserve `&F` calling style without forcing every caller to a blanket
move, add a blanket impl in `mod.rs`:

```rust
impl<T: RefSeqFetcher + ?Sized> RefSeqFetcher for &T {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32)
        -> Result<Vec<u8>, std::io::Error>
    {
        (**self).fetch(chrom_id, start_1based, length)
    }
}
```

Callers may then construct either `PileupWalker<_, MockFasta>` (moved) or
`PileupWalker<_, &MockFasta>` (borrowed); both are valid.

## Internal restructure

`WalkerState` is preserved verbatim — same fields, same per-step logic.
The only changes are at the seam:

- `close_aged_records(&mut self, tx: &SyncSender<PileupRecord>)` →
  `close_aged_records_into(&mut self, out: &mut VecDeque<PileupRecord>)`.
  Drains aged records into a finalised batch and stamps lifecycle marks
  onto the **first record of this tick's batch**, indexed via
  `out.get_mut(batch_start_idx)`. Behaviour is identical to today's stamp
  on `records[0]`.
- `flush_chromosome(&mut self, tx)` → `flush_chromosome_into(&mut self, out)`
  with the same shape.
- The outer loop that today lives inside `run` moves into a new
  `PileupWalker::fill_pending(&mut self) -> Result<(), WalkerError>` method.
  It loops until at least one record sits in `pending`, or until end-of-input
  (in which case it flushes the final chromosome and sets `done = true`).

`Iterator::next()`:

1. If `pending` non-empty, pop front and return `Ok(record)`.
2. If `done`, return `None`.
3. Otherwise call `fill_pending()`. On `Ok(())`, pop front (or return `None`
   if end-of-input drained empty). On `Err(e)`, set `done = true` and
   return `Some(Err(e))` once — subsequent calls return `None`.

The pending buffer is `VecDeque<PileupRecord>`: a single walker tick may
emit 0–many records (e.g. a wide deletion at an earlier anchor unblocks
several narrower records simultaneously), and callers consume them one at
a time across `next()` invocations. Order within a tick is preserved
(`pop_front` matches the order `close_aged_records_into` pushed), so the
"lifecycle marks attach to the first record of this tick's batch" invariant
is unchanged.

## Error-variant cleanup

`WalkerError::ChannelClosed` becomes unreachable. Delete the variant and
the regression test that pins it
(`run_returns_channel_closed_when_receiver_dropped_mid_stream`). All other
variants and their tests are preserved.

`DEFAULT_OUTPUT_CHANNEL_CAPACITY` (and its doc comment) is deleted from
`pileup/mod.rs` — no longer applicable. No production caller references it
(only the examples and tests construct channels, and they're rewritten too).

## Test strategy

Existing test coverage is preserved — only the harness changes. The
`drive_walker_with_config` helper in [src/per_sample_caller/pileup/tests.rs](../../src/per_sample_caller/pileup/tests.rs)
becomes:

```rust
pub fn drive_walker_with_config(
    reads: Vec<PreparedRead>,
    ref_fetcher: MockFasta,
    config: &WalkerConfig,
) -> (Vec<PileupRecord>, RunSummary) {
    let mut walker = run(reads, &ref_fetcher, config);
    let records: Vec<_> = (&mut walker)
        .map(|r| r.expect("walker yielded error"))
        .collect();
    let summary = walker.summary();
    (records, summary)
}
```

No `mpsc::sync_channel`, no `thread::spawn`, no `drop(tx)` dance. Tests
calling `super::run` directly (the error-variant regressions) are rewritten
to iterate and extract the error with `.find_map(|r| r.err())` or
equivalent.

New test:

- `next_after_error_returns_none` — pin that the iterator stops yielding
  once it has reported an error, even if `next()` is called again.

The `ChannelClosed`-specific test is deleted.

## Examples and benches

Both rewrite analogously:

- [examples/dhat_pileup.rs](../../examples/dhat_pileup.rs): replace the
  `(tx, rx)` + collector thread with a direct `for record in run(...)` loop
  that counts records.
- [benches/pileup_walker_scaling.rs](../../benches/pileup_walker_scaling.rs):
  `setup_walker` no longer pre-allocates a channel; `drive_walker`
  iterates and counts. The benches retain their current parametrisation
  (read length, multi-op, etc.); only the harness changes.

A side-effect of the bench rewrite: the previous timed body included the
`drop(tx) + collector.join()` teardown as walker-attributable work
(comment at [benches/pileup_walker_scaling.rs:188-197](../../benches/pileup_walker_scaling.rs#L188-L197)).
The new body removes that overhead, so the bench numbers will shift
slightly downward. This is expected — it's the channel cost going away
— and is itself a positive signal for the rewrite. Re-baseline the
bench numbers post-rewrite.

## Validation

Inside the dev container (`./scripts/dev.sh`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo build --examples`
- `cargo build --benches`

Bench wall-time comparison is informational (not required for correctness)
and tracked separately; this plan does not block on a bench number target.

## Assumptions / silent choices

- **Blanket `&T: RefSeqFetcher` impl.** The current trait already takes
  `&self`, so this is a zero-cost forwarding impl that lets callers continue
  passing `&fetcher` without forcing a move. Decided over making the iterator
  generic with a lifetime parameter `'f` over `&'f F` because the former is
  strictly more flexible (callers can move *or* borrow) at no behavioural cost.
- **Iterator error semantics: terminate on first error.** Once `next()` returns
  `Err(e)`, all subsequent calls return `None`. This matches the previous
  shape where `run` returned `Err` and stopped emitting. Resumable error
  handling was never a feature of the old API.
- **`pending` is `VecDeque<PileupRecord>`** rather than `Vec` with a cursor.
  Per-tick batches are small (typically 1–few records); `VecDeque::pop_front`
  is O(1) and the order-preservation it provides is what the lifecycle-mark
  stamping rule relies on.
- **`flush_chromosome_into` is called at end-of-input even when there are
  no open records.** Cheap (no-op path) and preserves the existing
  `flush_chromosome(tx)` shape at end-of-loop in today's `run`.

## Risks

- **Borrow-checker friction** between `Peekable::peek(&mut self.reads)` and
  field accesses on `self.state` / `self.pending`. The original `run` works
  because `reads_iter` and `state` are separate locals; inside the
  iterator they are sibling fields, so split borrows must work via
  direct field access. The expected fix if it fights us is to copy peeked
  fields into a local before the call that needs `&mut self.state`, as the
  chromosome-boundary branch in the plan already does.
- **Bench numbers may drop slightly** due to channel overhead going away.
  This is desired (and was identified as one of the reasons for the
  rewrite), not a regression to chase. Confirm the drop is in the
  expected direction; investigate only if numbers move the other way.

## Out-of-scope follow-ups

- **Tiny optional prefetch adapter** (e.g. a 1–2-slot buffer fed by a
  helper thread) if profiling later shows jitter matters. Don't pre-build
  it.
- **Stage 1 multi-sample driver** (rayon over BAMs). The natural place for
  cross-sample parallelism after this rewrite.

## File touch list

- `src/per_sample_caller/pileup/walker.rs` — main change.
- `src/per_sample_caller/pileup/mod.rs` — blanket impl, delete
  `DEFAULT_OUTPUT_CHANNEL_CAPACITY`, re-export `PileupWalker`.
- `src/per_sample_caller/pileup/errors.rs` — delete `ChannelClosed` variant.
- `src/per_sample_caller/pileup/tests.rs` — rewrite `drive_walker_with_config`
  and the error-variant regressions; delete the `ChannelClosed` test.
- `examples/dhat_pileup.rs` — direct iteration.
- `benches/pileup_walker_scaling.rs` — direct iteration in setup/drive.
