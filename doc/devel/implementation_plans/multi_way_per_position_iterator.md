# Multi-way per-position iterator (cohort-side merge over `.psp` files)

Proposal date: 2026-05-14.

## Domain intent

The first cohort-side component of the multi-sample SNP caller (see
[ia/specs/calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)).
It sits between Stage 2 (the per-sample `.psp` reader) and Stage 3
(the DUST low-complexity filter), and its job is exactly one thing:

```
sample_A.psp ─► PspReader ─┐
sample_B.psp ─► PspReader ─┼─► PerPositionMerger ─► (chrom_id, pos, per-sample pileups)
sample_C.psp ─► PspReader ─┘
```

Given N per-sample iterators that each yield `PileupRecord`s in
genomic order, the merger yields one item per *genomic position any
sample covers*, with an N-slot vector telling each downstream stage
which samples had a record there. Stages 3 (DUST), 4 (grouping), 5
(per-group), and 6 (posterior) all consume this stream — none of
them re-opens `.psp` files individually.

Pipeline-overview diagram: [calling_pipeline_architecture.md:33-48](../specs/calling_pipeline_architecture.md#L33-L48).
Stage-3 placement: [calling_pipeline_architecture.md:883-895](../specs/calling_pipeline_architecture.md#L883-L895).
Grouping precedent (analogous merge over the old gVCF path): [src/variant_grouping.rs](../../src/variant_grouping.rs).

## Why now

Stage 1 is feature-complete as a kit of parts (CRAM input → BAQ →
walker → `pileup_to_psp` seam → `PspWriter`/`PspReader`). Every
downstream stage (DUST, grouping, per-group, posterior) consumes the
multi-way merge, so building it unblocks all of them at once. It is
also the place where the chromosome-id-space agreement across `.psp`
files becomes a hard precondition — pinning that here, once,
simplifies everything that follows.

## Algorithm: linear-scan k-way merge

For each output position the merger does an O(N) scan over peeked
heads to find the min `(chrom_id, pos)`, pulls successors only from
the readers that were tied at that min, and yields a `PerPositionPileups`.

**Why linear scan and not a heap.** The textbook "heap wins at large
N" rule of thumb assumes each output position only advances a few
readers — the heap's `O(log N)` per push/pop pays off because most
readers are skipped. In the cohort access pattern this caller targets
(WGS data, `.psp` records emitted at every covered reference
position, N ≈ thousands), every reader is tied at nearly every output
position. With `k ≈ N` per position the totals are:

| | Per output position | Genome total (P positions, M ≈ N·P records) |
|---|---|---|
| Linear scan | ~2N (find min + take ties) | `O(N·P + M)` ≈ `2·N·P` |
| Min-heap | ~2k log N (pop + push tied) | `O(M log N)` ≈ `N·P·log N` |

At N=2500 the heap does ~5× more work than the linear scan. The
crossover the rule of thumb predicts simply does not happen in this
data shape — the deciding variable is sparsity (`k/N`), not N.

The heap *would* win if a future input distribution made the access
pattern sparse: exome / targeted-panel `.psp` files, a variant-only
`.psp` schema, or cohorts where samples cover disjoint contig sets.
None of those is in scope today (the format spec emits a record per
covered reference position — [calling_pipeline_architecture.md:529-684](../specs/calling_pipeline_architecture.md#L529-L684)),
so the merger ships with linear scan only; see "Out-of-scope
follow-ups" for the criteria that would trigger building the heap.

## Out of scope

- The heap variant. Not built; revisit only on a data-shape change.
- The DUST filter (Stage 3). The merger emits unfiltered items;
  DUST sits on top of the iterator as a separate adaptor.
- Parallel per-reader prefetch / threaded readahead. The pull-iterator
  pattern already decided against per-sample threading
  ([pileup_pull_iterator.md:30-31](pileup_pull_iterator.md#L30-L31)),
  and cross-sample parallelism happens at the rayon-over-groups level
  in Stage 5, not here.
- A "self-bundling" constructor that takes file paths and manages
  `PspReader` ownership for the caller. The merger is generic on
  `Iterator<Item = Result<PileupRecord, PspReadError>>`; the caller
  is responsible for opening readers and handing iterators in. A
  thin helper for the path-based workflow can be added later if its
  shape becomes obvious (see "Out-of-scope follow-ups").

## API shape

New module: `src/cohort/` (first occupant of the cohort-side
namespace). The Stage 3/4/5/6 code lands here over time.

```rust
// src/cohort/mod.rs
pub mod per_position_merger;

// src/cohort/per_position_merger.rs
use crate::per_sample_caller::pileup::PileupRecord;
use crate::per_sample_caller::psp::header::ParsedChromosome;
use crate::per_sample_caller::psp::PspReadError;

pub struct PerPositionMerger<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{ /* fields private */ }

impl<I> PerPositionMerger<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    /// Construct a merger over `readers`. `sample_names` and
    /// `chromosomes` are taken from the caller-validated headers
    /// (see "Header validation" below for what the caller must
    /// check before calling this).
    ///
    /// Returns an error only on per-reader prefetch failure (the
    /// first `next()` on any reader). Empty `readers` is allowed
    /// and yields an immediately-exhausted iterator.
    pub fn new(
        readers: Vec<I>,
        sample_names: Vec<String>,
        chromosomes: Vec<ParsedChromosome>,
    ) -> Result<Self, MergerError>;

    pub fn sample_names(&self) -> &[String];
    pub fn chromosomes(&self) -> &[ParsedChromosome];
    pub fn n_samples(&self) -> usize;
}

impl<I> Iterator for PerPositionMerger<I>
where
    I: Iterator<Item = Result<PileupRecord, PspReadError>>,
{
    type Item = Result<PerPositionPileups, MergerError>;
}

/// One emitted item: a single `(chrom_id, pos)` with a per-sample
/// slot for every sample passed to `PerPositionMerger::new`.
/// Samples without a record at this position carry `None`.
#[derive(Debug, Clone, PartialEq)]
pub struct PerPositionPileups {
    pub chrom_id: u32,
    pub pos: u32,
    /// Indexed by sample order (same as `sample_names()`).
    /// `per_sample.len() == n_samples()`.
    pub per_sample: Vec<Option<PileupRecord>>,
}

#[derive(thiserror::Error, Debug)]
pub enum MergerError {
    #[error("reader {sample_idx} ({sample_name}) failed: {source}")]
    Reader {
        sample_idx: usize,
        sample_name: String,
        #[source]
        source: PspReadError,
    },
    #[error("reader {sample_idx} ({sample_name}) emitted out-of-order record: {chrom_id}:{pos}")]
    OutOfOrder {
        sample_idx: usize,
        sample_name: String,
        chrom_id: u32,
        pos: u32,
    },
}
```

The generic-iterator shape matches the pull-iterator precedent
([src/per_sample_caller/pileup/walker.rs](../../src/per_sample_caller/pileup/walker.rs))
and makes the merger trivially testable with synthetic
`std::iter`-based fixtures — no on-disk `.psp` required.

Typical caller wiring:

```rust
let mut readers: Vec<PspReader<BufReader<File>>> = open_inputs(&paths)?;
let chromosomes = check_chromosome_agreement(&readers)?; // see below
let sample_names = readers.iter().map(|r| r.header().sample.clone()).collect();
let iters: Vec<_> = readers.iter_mut().map(|r| r.records()).collect();
let merger = PerPositionMerger::new(iters, sample_names, chromosomes)?;
for pileups in merger {
    let pileups = pileups?;
    // hand to DUST / grouping / ...
}
```

The merger does not own the `PspReader`s — they live on the caller's
stack alongside the iterators. This keeps lifetimes plain (no
self-referential bundling) at the cost of asking the caller to hold
two parallel `Vec`s for the duration of the merge.

## Algorithm details

State kept by the merger:

```rust
struct PerPositionMerger<I> {
    readers: Vec<I>,
    /// Peeked next record per reader. `None` once that reader has
    /// produced `None`. Same length as `readers`.
    heads: Vec<Option<PileupRecord>>,
    sample_names: Vec<String>,
    chromosomes: Vec<ParsedChromosome>,
    /// Last `(chrom_id, pos)` yielded, for per-reader monotonicity
    /// checks. `None` before the first yield.
    last_emitted: Option<(u32, u32)>,
    /// `true` once we have surfaced an error or exhausted every
    /// reader. Latches `next()` to `None` thereafter.
    done: bool,
}
```

`PerPositionMerger::new` prefetches one record from every reader to
populate `heads`. A reader that yields an immediate error short-
circuits construction with `MergerError::Reader`.

`Iterator::next()`:

1. If `done`, return `None`.
2. Find the minimum `(chrom_id, pos)` across non-`None` entries of
   `heads` with a single linear scan. If every entry is `None`, set
   `done = true` and return `None`.
3. Validate strict monotonicity: if `last_emitted` is `Some(lo)` and
   the new min ≤ lo, surface `MergerError::OutOfOrder` from the
   offending reader. (One reader's `RecordsIter` already guarantees
   per-reader monotonicity internally; this check defends against
   pathological mocks and against drift if a future `RecordsIter`
   variant relaxes that property.)
4. Build a `PerPositionPileups` of length `n_samples()`. For each
   reader index `i`:
   - If `heads[i]` matches the min `(chrom_id, pos)`, move it into
     `pileups.per_sample[i]` (so the slot is `Some(record)`) and
     refill `heads[i]` from `readers[i].next()`.
   - Otherwise leave the slot as `None`.
5. Update `last_emitted` and return `Ok(pileups)`.

A reader error during refill (step 4) surfaces as
`MergerError::Reader { sample_idx, sample_name, source }`, with
`done` latched so subsequent `next()` calls return `None`. This
matches the walker's terminate-on-first-error contract.

Memory: `O(N)` peeked records + one in-flight `PerPositionPileups`.
No buffering beyond a single record per reader.

## Header validation

The merger's correctness depends on all input `.psp` files
agreeing on the chromosome id space — `chrom_id == 7` must mean
the same chromosome for every sample. This is a *caller-side*
precondition. The merger only takes the post-validation
`chromosomes` vector and assumes it is the common ground.

A standalone helper for the path-based opener lives alongside the
merger:

```rust
/// Verify every reader agrees on chromosome ordering / lengths /
/// MD5s with `readers[0]`, and return the agreed-upon list.
///
/// Returns `Err(MergerError::ChromosomeMismatch { … })` on the
/// first divergence.
pub fn check_chromosome_agreement<R: Read + Seek>(
    readers: &[PspReader<R>],
) -> Result<Vec<ParsedChromosome>, MergerError>;
```

Agreement is checked field-by-field on `ParsedChromosome`
(`name`, `length`, `md5`). Reference-string mismatches that
*don't* result in chromosome divergence (e.g. two `.psp` files
naming the same reference build differently) are not fatal here;
the chromosome list is the operative contract.

Add a `ChromosomeMismatch` variant to `MergerError`:

```rust
#[error("sample {sample_idx} ({sample_name}) chromosome {chrom_id} disagrees with sample 0: {detail}")]
ChromosomeMismatch {
    sample_idx: usize,
    sample_name: String,
    chrom_id: u32,
    detail: String,
},
```

## Test strategy

All tests are unit-level in `cohort/per_position_merger.rs`'s
`#[cfg(test)]` module. They use synthetic iterators
(`std::iter`-based `Result<PileupRecord, PspReadError>` streams),
so no on-disk `.psp` round-trip is required — that's covered by
the writer/reader contract tests.

Required cases:

- **Empty cohort.** `readers = vec![]` yields an immediately
  exhausted iterator. `n_samples() == 0`.
- **Single reader.** Identity pass-through: each `PileupRecord`
  becomes a `PerPositionPileups` of length 1 with
  `per_sample[0] = Some(r)`.
- **Two readers, fully overlapping positions.** Every item has
  both slots `Some`.
- **Two readers, disjoint positions.** Every item has exactly
  one slot `Some`; ordering follows the merged sequence of all
  positions.
- **Two readers, partial overlap.** Mixed items. Covers the
  "WGS-like" expected shape.
- **Multi-chromosome.** Two readers, both spanning two chromosomes.
  Confirms chrom_id is the major sort key and pos is the minor one.
- **Three readers, k=2 tied.** Confirms the "advance only tied
  readers" invariant — the third reader's head is not consumed
  when it's not at the current min.
- **Reader error mid-stream.** A synthetic iterator that yields
  `Err(PspReadError::…)` after some records. Confirms
  `MergerError::Reader { sample_idx, sample_name, … }` surfaces
  and that subsequent `next()` calls return `None`.
- **Reader error on prefetch.** First record from one reader is
  `Err`; `new()` returns `MergerError::Reader`.
- **Out-of-order detection.** A synthetic iterator that yields
  a backwards `(chrom_id, pos)` is caught with
  `MergerError::OutOfOrder`.
- **Emission order across `next()`.** Items are yielded in strictly
  increasing `(chrom_id, pos)`. (Pinned with an explicit assertion
  loop in one of the multi-reader tests.)

`check_chromosome_agreement` has its own tests, using
in-memory `PspReader<Cursor<Vec<u8>>>` constructed with the
existing `writer_header(…)` helper at varying chromosome
configurations:

- Identical chromosome lists across two readers → `Ok`.
- Differing length on one chromosome → `Err(ChromosomeMismatch)`
  with the right `sample_idx` and `chrom_id`.
- Differing MD5 → `Err(ChromosomeMismatch)`.
- Differing name (renamed contig) → `Err(ChromosomeMismatch)`.
- Differing number of chromosomes → `Err(ChromosomeMismatch)`.

## Validation

Inside the dev container (`./scripts/dev.sh`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo build --examples`
- `cargo build --benches`

No new benchmark is added in this plan. Benchmarking is justified
once Stages 3-5 exist and can drive realistic load through the
merger; benching it in isolation would over-fit a synthetic
generator.

## Assumptions / silent choices

- **`Vec<Option<PileupRecord>>` per emitted item**, not a sparse
  representation. WGS-style coverage means most slots are `Some` at
  most positions, so the `Option` discriminator overhead is small
  and downstream stages benefit from O(1) "what did sample j see?"
  lookups by index. A sparse encoding (`SmallVec<[(idx, rec); 8]>`)
  would win for exome/panel data — pair it with a heap merger if
  the access pattern ever shifts that way.
- **Merger is generic on the iterator type**, not coupled to
  `RecordsIter<'_, R>` specifically. Mirrors the
  [pull-iterator precedent](pileup_pull_iterator.md) and makes
  testing trivial. Concrete type for production callers will be
  `RecordsIter<'_, BufReader<File>>` once the path-opening helper
  lands; nothing in the merger needs to know that.
- **Sample-name uniqueness is not checked.** Two `.psp` files with
  the same sample name is unusual (replicates, mislabelled
  inputs) but not the merger's policy call. Stage 5/6 may want
  to reject that — they have the context to do so.
- **Per-reader monotonicity is re-checked at the merge layer**, even
  though `RecordsIter` already enforces it within one block. Cheap
  (one tuple comparison per emitted item) and defends against drift
  and pathological mocks.
- **Errors latch.** Once any reader fails or the iterator emits an
  error, all subsequent `next()` calls return `None`. Same shape as
  the walker's pull-iterator contract.

## Risks

- **Lifetime ergonomics at the call site.** Holding `Vec<PspReader>`
  and `Vec<RecordsIter<'_, _>>` in parallel works (the iterators
  borrow from the readers) but the caller must keep both alive for
  the merger's lifetime. If this turns out to be awkward in
  practice — e.g. when the path-based helper lands and wants to
  return a single bundled value — revisit with a `self_cell` or
  manual `Pin` wrapper. Don't pre-build that abstraction.
- **Pileups clone cost in tests.** `PerPositionPileups: Clone` is only
  needed for `PartialEq`-driven assertions in tests; production
  paths move the value out. Confirm production calls don't end up
  cloning by accident.

## Out-of-scope follow-ups

- **Heap-based merger** (`PerPositionHeapMerger`) behind the same
  public surface. Build only if a future input pattern makes
  `k/N` small in practice: exome / panel datasets, a variant-only
  `.psp` schema, cohorts with disjoint contig coverage. The trigger
  is data shape, not N alone — the table at "Algorithm" sets the
  threshold quantitatively.
- **Path-based opener helper.** `open_psp_files(paths) ->
  Result<PerPositionMerger<…>, _>` that owns the readers itself.
  Wait until a real caller (the Stage 1 CLI subcommand, or the
  cohort entry point) has shipped and the ergonomic shape is
  obvious. Likely needs `self_cell` or similar to bundle readers
  and iterators in a single owned value.
- **Parallel per-reader prefetch.** A small per-reader prefetch
  thread or async readahead. Only worth it if profiling shows the
  merger is IO-bound on `PspReader::next` rather than on
  downstream stages. Defer until profiles exist.
- **Sparse `PerPositionPileups` encoding.** `SmallVec<[(sample_idx, record); K]>`
  instead of `Vec<Option<…>>`. Paired with the heap merger if the
  cohort access pattern ever becomes sparse.

## File touch list

- `src/cohort/mod.rs` — new module root, declares
  `pub mod per_position_merger;`.
- `src/cohort/per_position_merger.rs` — new file: `PerPositionMerger`,
  `PerPositionPileups`, `MergerError`, `check_chromosome_agreement`,
  full `#[cfg(test)]` module.
- `src/lib.rs` — `pub mod cohort;` added alongside the existing
  module declarations.
