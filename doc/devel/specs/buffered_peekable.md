# `BufferedPeekable` — shared peekable adapter for fallible iterators

**Status:** Draft, 2026-04-27. Specification of a small shared
utility used by Stage 1 (per-sample caller) and Stage 4 (variant
grouping) of the calling pipeline. Lives in `src/` and is referenced
from [per_sample_caller.md](per_sample_caller.md) and (eventually)
the cohort-side spec.

## Purpose

A streaming wrapper around any `Iterator<Item = Result<T, E>>` that
adds two things at once:

1. **A `peek()` operation** in the transposed shape
   `Result<Option<&T>, E>`, so callers see the head of the stream
   without consuming it and can propagate errors by value.
2. **A configurable in-memory buffer** of already-pulled items, so
   the inner iterator runs in batched bursts (`buffer_size` items
   per refill) rather than one record per outer call. This keeps the
   inner iterator's hot loop hot and amortises any per-call overhead.

Both behaviours are needed today by the project's fallible parsers
(VCF lines parsed into `Variant`, CRAM records decoded by noodles).
Standard library `Peekable<I>` does not give either: with
`Item = Result<T, E>` it returns `Option<&Result<T, E>>`, leaving
errors trapped inside the buffer, and it never batch-pulls.

## Why it exists in the project

Two concrete callers, today and tomorrow:

- **Stage 1, per-sample caller.** Wraps each input CRAM's record
  iterator in `BufferedPeekable` and merges them by `(ref_id, pos)`
  with peek-and-scan ([per_sample_caller.md:153](per_sample_caller.md#L153)).
  CRAM record iteration is already block-buffered inside noodles, so
  Stage 1 will use a small `buffer_size` (1 is sufficient — the
  buffer is there for peek semantics, not batching).

- **Stage 4, variant grouping.** Already does peek-and-scan with a
  hand-rolled `peek_variant()` and an internal `vars_buffer` of 100
  variants ([gvcf_parser.rs:569](../../src/gvcf_parser.rs#L569),
  [variant_grouping.rs:172](../../src/variant_grouping.rs#L172)).
  The hand-rolled version predates this utility and should be
  retrofitted — see §"Follow-up: retrofit `VarIterator`" below.

Two callers is enough to justify extracting the abstraction; this is
not speculative reuse.

## API

```rust
pub const DEFAULT_BUFFER_SIZE: usize = 100;

pub struct BufferedPeekable<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    inner: I,
    buffer: VecDeque<Result<T, E>>,
    buffer_size: usize,
}

impl<I, T, E> BufferedPeekable<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    /// Wrap `inner` with the default buffer size.
    pub fn new(inner: I) -> Self;

    /// Wrap `inner` with an explicit `buffer_size`. Must be ≥ 1.
    /// The internal `VecDeque` is pre-allocated to this capacity so
    /// it never reallocates.
    pub fn with_buffer_size(inner: I, buffer_size: usize) -> Self;

    /// Look at the next item without consuming it.
    ///
    /// - `Ok(Some(&t))` — head of the stream is item `t`; further
    ///   peek calls return the same reference until `next()` is
    ///   called.
    /// - `Ok(None)` — the inner iterator is exhausted; subsequent
    ///   calls keep returning `Ok(None)`.
    /// - `Err(e)` — the inner iterator yielded an error at this
    ///   point in the stream. The error is *consumed by this peek*
    ///   (see §"Peek consumes errors" below); the next peek looks
    ///   at whatever follows.
    pub fn peek(&mut self) -> Result<Option<&T>, E>;
}

impl<I, T, E> Iterator for BufferedPeekable<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    type Item = Result<T, E>;
    fn next(&mut self) -> Option<Self::Item>;
}
```

Two constructors:

- `new(inner)` is the common case — `DEFAULT_BUFFER_SIZE = 100`,
  matches the existing `VarIterator::fill_buffer` target, fine for
  any caller that wants both peek and batched parsing.
- `with_buffer_size(inner, n)` is for callers who want explicit
  control. Stage 1 calls `with_buffer_size(decoder, 1)` because the
  CRAM decoder underneath already does its own batching; a fat outer
  buffer would just delay records without saving work.

`buffer_size` of 0 is rejected at construction (degenerate — there's
no head to peek). 1 is allowed and gives pure peek semantics with no
batching.

## Behaviour

### Internal invariant

The buffer is the **prefix of the inner iterator's output that has
been pulled but not yet consumed by the caller**, in the same order.
It contains both successful and failed pulls verbatim:

```
caller has consumed k items, inner has produced k + b items
                              └────────── b items in buffer ─────────┘
buffer = [inner.output[k], inner.output[k+1], ..., inner.output[k+b-1]]
```

This invariant is what makes the unified `VecDeque<Result<T, E>>`
work: errors and items keep their natural order, and there is no
separate "pending error" slot to keep in sync.

### Refill

A refill runs whenever `peek` or `next` is called and the buffer is
empty. It pulls from the inner iterator until the buffer holds
`buffer_size` items or the inner iterator returns `None`:

```
loop {
    if buffer.len() >= buffer_size { break; }
    match inner.next() {
        Some(item)  => buffer.push_back(item),  // Ok or Err, both go in
        None        => break,                   // exhausted
    }
}
```

Errors do **not** stop the refill loop. The refill happily pushes
`Err` items into the buffer alongside `Ok` items; the order is
preserved and the caller sees them at their correct point in the
stream. The only thing that stops a refill is the buffer reaching
`buffer_size` or the inner iterator returning `None`.

### Peek consumes errors

`peek()` returns `Result<Option<&T>, E>`. When the front of the
buffer is `Err(_)`, peek **takes** the error from the buffer and
returns it by value:

```rust
pub fn peek(&mut self) -> Result<Option<&T>, E> {
    if self.buffer.is_empty() { self.refill(); }
    if matches!(self.buffer.front(), Some(Err(_))) {
        match self.buffer.pop_front() {
            Some(Err(e)) => return Err(e),
            _ => unreachable!(),
        }
    }
    Ok(self.buffer.front().map(|r| match r {
        Ok(t) => t,
        Err(_) => unreachable!("error case handled above"),
    }))
}
```

This breaks the usual "peek does not advance" contract for the
specific case of errors. The justification: an error is not a record
in the stream the caller is iterating; it is an out-of-band event
that has to be handled once and then is gone. Surfacing the same
error on every subsequent peek would either require `Clone` on `E`
(which `std::io::Error` does not satisfy) or borrowing the error
(which prevents callers from propagating it).

Concretely, two peek calls with no `next()` in between can return
different things only when the first one popped an error. The second
peek then sees whatever was behind that error in the stream — the
next record, the next error, or `None`.

### Next

`next()` always behaves as a regular `Iterator::next` — it returns
`Option<Result<T, E>>`. If the buffer's front is the same item that
peek would have shown, `next` returns it (consuming it). If peek had
just popped an error, `next` returns whatever comes next.

Calling `next` after the inner iterator has returned `None` keeps
returning `None` (we do not implement `FusedIterator` explicitly,
but the behaviour is fused as a side effect of "empty inner returns
None").

### Error retry semantics

After an error has surfaced (via peek-consumes or via `next`), the
buffer is in whatever state it was minus that error. The next
`peek`/`next` will refill if needed by calling `inner.next()` again.
This matches Rust's `Iterator` convention that an iterator may be
called after returning `Err` and will produce the next thing the
underlying source gives.

If the underlying error is persistent (e.g. corrupt file), the
caller will see the same kind of error repeatedly; it is the
caller's job to stop iterating when an error is fatal.

## Comparison with `std::iter::Peekable`

| | `std::iter::Peekable<I>` | `BufferedPeekable<I, T, E>` |
|---|---|---|
| `peek()` shape | `Option<&T>` | `Result<Option<&T>, E>` |
| Works for `Item = Result<T, E>` | yes, but `peek()` returns `Option<&Result<T, E>>` | yes, transposed |
| Buffer size | 1 (head only) | configurable, default 100 |
| Errors | trapped inside the head slot | consumed once on peek, by value |

`BufferedPeekable` is strictly a superset for fallible iterators.
For infallible iterators std's `Peekable` is fine and we should not
re-roll it.

## Tests

Behavioural tests, all on a synthetic inner iterator that produces a
known sequence:

- Empty inner: `peek` and `next` both return `Ok(None)` / `None`
  forever.
- All `Ok`: peek returns the same head until `next` consumes; ordered
  delivery.
- Single error: `[Ok, Ok, Err, Ok]` — verify order is preserved and
  peek-consumes-error works.
- Error first: `[Err, Ok]` — verify error surfaces on first peek and
  the buffer recovers to deliver `Ok` next.
- Multiple errors: `[Err, Err, Ok]` — verify each error surfaces in
  turn.
- `buffer_size = 1`: minimal buffering, peek-only behaviour.
- `buffer_size = N` larger than the entire inner stream: refill
  drains inner exactly once, then `next` delivers everything until
  exhaustion.
- `with_buffer_size(_, 0)` is rejected (panic or `Result` on
  construction — pick at implementation time).
- After exhaustion, further calls return `None` / `Ok(None)`.

No criterion benchmarks for v1 — the adapter is a thin wrapper and
its overhead is dwarfed by inner-iterator work in every realistic
caller.

## File layout

Single file, single module:

```
src/buffered_peekable.rs    — struct, impls, tests
```

Re-exported from `lib.rs` so callers in the same crate use it as
`crate::buffered_peekable::BufferedPeekable`.

## Follow-up: retrofit `VarIterator`

`VarIterator` in `src/gvcf_parser.rs` predates this adapter and
implements its own buffer + peek. The retrofit is small but it is a
real change to `VarIterator`'s public API and so is **out of scope
for the initial `BufferedPeekable` introduction** — we land the
adapter, use it from Stage 1 first, and come back to `VarIterator`
when we touch that area.

The shape of the retrofit, when we do it:

1. Remove `vars_buffer: VecDeque<Variant>` and
   `peek_variant(&mut self) -> VcfResult<Option<&Variant>>` from
   `VarIterator`.
2. Shrink `VarIterator::Iterator::next` to "read one line, parse,
   return one `VcfResult<Variant>`" — no batching, no buffer.
3. At every call site that constructs a `VarIterator` and uses
   `peek_variant` (mainly inside `variant_grouping.rs`), wrap the
   iterator in `BufferedPeekable::new(...)` and use `.peek()` and
   `.next()` from the adapter instead.
4. Verify that the new behaviour matches the old: same external
   semantics, errors surface in source order rather than eagerly
   (this is a small *fix*, not a regression — the existing
   `fill_buffer` propagates errors before delivering records that
   came chronologically earlier; the adapter delivers them in
   order).
5. Run the existing `variant_grouping.rs` and integration tests; no
   external behaviour change is expected for callers.

This follow-up is tracked here so that whoever next opens
`gvcf_parser.rs` for any reason has a one-stop reference for the
intended end state.
