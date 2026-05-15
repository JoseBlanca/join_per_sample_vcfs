# Multi-way per-position iterator — implementation report

Date: 2026-05-15.
Plan: [multi_way_per_position_iterator.md](../../implementation_plans/multi_way_per_position_iterator.md).
Slice: first cohort-side stage of the multi-sample SNP caller. Takes
N per-sample `.psp` iterators (each in genomic order) and emits one
`PerPositionPileups` item per `(chrom_id, pos)` any sample covers,
with N slots saying which samples had a record there.

## Plan

Followed the plan as written. Linear-scan k-way merge keyed on
`(chrom_id, pos)`; per-reader head buffer; tied-readers-only refill;
strict-monotonicity check across the emit stream; latch-on-error
contract matching the walker. New module
[src/cohort/](../../../../src/cohort/) (first occupant of the
cohort-side namespace), wired into [src/lib.rs](../../../../src/lib.rs#L19).

## Assumptions

Silent choices the plan left open:

- **Defensive `sample_names.len() == readers.len()` check at
  construction.** Surfaces as `MergerError::SampleCountMismatch`
  rather than panicking. The plan framed the two vectors as "taken
  from caller-validated headers", which would have allowed skipping
  the check; the user asked for defensiveness, so it's a hard
  precondition now. `chromosomes` is metadata only and is *not*
  length-checked (there's nothing to check it against — chromosome
  count is independent of sample count).
- **`MergerError` is not `Clone` or `PartialEq`.** Its `Reader`
  variant carries a `PspReadError`, which transitively carries
  `std::io::Error` — neither of those is `Clone`/`PartialEq`-able.
  Tests assert on variants via `matches!` and field destructuring
  rather than equality.
- **`MergerError::Reader` boxes its inner `PspReadError`.** Clippy's
  `result_large_err` lint flagged the variant at ≥128 bytes
  inline. `Box<PspReadError>` brings the variant down to ~40 bytes
  and keeps the merger's per-emit `Result` cheap.
- **`ChromosomeMismatch.detail` is a human-readable string**, not
  structured `got`/`expected` fields. Downstream code matches on
  the variant tag (`ChromosomeMismatch { .. }`); the detail is for
  display. Structured per-field divergence wasn't requested and
  would have ballooned the variant.
- **`check_chromosome_agreement(&[])` returns `Ok(vec![])`**, mirror
  of the merger's empty-cohort behaviour.
- **Out-of-order "offending reader" is the first head matching the
  regressing key.** When multiple readers tie at the same
  monotonicity-violating position, the lowest sample index is
  blamed. The plan said "the offending reader" without disambiguating;
  in production this only fires on pathological mocks, so the
  tiebreak is informational only.

## Changes made

- [src/cohort/mod.rs](../../../../src/cohort/mod.rs) — new module
  root; declares `pub mod per_position_merger;` and a one-paragraph
  doc explaining the cohort-side stages.
- [src/cohort/per_position_merger.rs](../../../../src/cohort/per_position_merger.rs)
  — the merger plus types:
  - `PerPositionPileups { chrom_id, pos, per_sample: Vec<Option<PileupRecord>> }`
    — emitted item.
  - `PerPositionMerger<I>` — generic on
    `Iterator<Item = Result<PileupRecord, PspReadError>>`. Owns
    readers + peeked heads + last-emitted key + done latch.
  - `MergerError` — `SampleCountMismatch` / `Reader` / `OutOfOrder` /
    `ChromosomeMismatch`. `thiserror::Error + Debug`,
    `#[non_exhaustive]`.
  - `check_chromosome_agreement<R: Read + Seek>(&[PspReader<R>])`
    — verifies count, name, length, MD5 agreement against
    `readers[0]`; returns a clone of the agreed list.
  - Custom `Debug` for `PerPositionMerger` (mirrors the
    `PspReader<R>` pattern at
    [reader.rs:73-81](../../../../src/per_sample_caller/psp/reader.rs#L73-L81)
    so the impl doesn't require `I: Debug`).
- [src/lib.rs](../../../../src/lib.rs#L19) — `pub mod cohort;` added
  alphabetically alongside the existing module declarations.

## Tests added/updated

All 18 tests live in [src/cohort/per_position_merger.rs](../../../../src/cohort/per_position_merger.rs)'s
`#[cfg(test)] mod tests` block (project precedent). Synthetic
`std::vec::IntoIter<Result<PileupRecord, PspReadError>>` streams
drive the merger — no on-disk `.psp` round-trip — and
`PspReader<Cursor<Vec<u8>>>` built from the existing `writer_header`
fixture drives `check_chromosome_agreement`.

Merger:
- `empty_cohort_yields_immediately` — empty `readers` → 0 samples,
  immediately exhausted.
- `single_reader_identity_pass_through` — each record becomes a
  one-slot `PerPositionPileups`.
- `two_readers_fully_overlapping` — both slots `Some` at every
  position.
- `two_readers_disjoint_positions` — exactly one slot `Some` per
  emit; merged-sequence ordering.
- `two_readers_partial_overlap` — WGS-like shape; checked slot
  presence position-by-position.
- `multi_chromosome_chrom_is_major_sort_key` — chrom 0 fully drained
  before any chrom 1 emit.
- `three_readers_k2_tied_does_not_consume_unrelated_head` — the
  "advance only tied readers" invariant: third reader's head is
  preserved when not at the current min.
- `reader_error_mid_stream_surfaces_and_latches_done` — three
  successful emits, then a refill error on reader B; subsequent
  `next()` calls return `None`. Note: refill errors abort the
  in-progress emission by design.
- `reader_error_on_prefetch_aborts_construction` — first `next()`
  on a reader yields `Err` → `PerPositionMerger::new` returns
  `MergerError::Reader`.
- `out_of_order_record_is_detected` — synthetic reader regresses
  `(0,5) → (0,3)`; merger surfaces `MergerError::OutOfOrder`.
- `emission_order_is_strictly_increasing` — windowed assertion
  across all emits.
- `sample_count_mismatch_rejected_at_construction` — defensive
  length check.

`check_chromosome_agreement`:
- `chromosome_agreement_empty_slice_is_ok` — `&[]` → `Ok(vec![])`.
- `chromosome_agreement_identical_lists_passes` — two readers built
  from the same `writer_header(2)` → `Ok`.
- `chromosome_agreement_differing_length_fails` — mutate
  `chromosomes[0].length` → `ChromosomeMismatch` with `detail`
  containing "length".
- `chromosome_agreement_differing_md5_fails` — mutate MD5.
- `chromosome_agreement_differing_name_fails` — mutate name.
- `chromosome_agreement_differing_count_fails` — push an extra
  chromosome; `detail` contains "count".

## Validation results

All run inside the dev container (`./scripts/dev.sh`):

- `cargo clippy --all-targets --all-features -- -D warnings` —
  clean (after boxing `PspReadError` for `result_large_err`, using
  `.iter().flatten()` for the head scan, and looping over
  `per_sample.iter_mut().enumerate()` instead of a raw range).
- `cargo test --all-features --lib` — **516 passed, 0 failed**
  (includes 18 new cohort tests).
- `cargo test --all-features --tests` — **all 12 integration tests
  pass** in the affected suites; full set is green.
- `cargo build --examples` — clean.
- `cargo build --benches` — clean.
- `cargo fmt --check` — the three files in this change (`src/lib.rs`,
  `src/cohort/mod.rs`, `src/cohort/per_position_merger.rs`) are
  format-clean. Pre-existing drift in
  `src/per_sample_caller/pileup/walker.rs` and
  `src/per_sample_caller/psp/reader.rs` (last touched in commits
  `b0a1c13` / `a805f70`) is unrelated to this change and was left
  alone; fix in a separate cleanup pass.

Not run: the `gvcf_perf` bench's actual execution requires a fixture
file outside the repo (`/home/jose/analyses/g2psol/source_data/TS.vcf.gz`)
and panics at startup without it. This is pre-existing environmental
state, not a regression.

## Tradeoffs and follow-ups

Deferred per the plan's "Out-of-scope follow-ups" section:

- **Heap-based merger.** Linear scan only; the heap variant is
  worth building only if a future data shape makes `k/N` small
  (exome / panel, variant-only `.psp` schema, disjoint contig
  coverage).
- **Path-based opener helper.** `open_psp_files(paths) ->
  Result<PerPositionMerger<…>, _>` that owns the readers itself.
  Wait until a real caller (Stage 1 CLI subcommand, cohort entry
  point) drives the ergonomic shape.
- **Parallel per-reader prefetch.** Only worth it if profiling
  shows the merger is IO-bound on `PspReader::next`. Defer until
  profiles exist.
- **Sparse `PerPositionPileups` encoding.** `SmallVec<[(idx,
  record); K]>` instead of `Vec<Option<...>>`. Pair with the heap
  merger if access patterns ever shift sparse.
- **No benchmark added.** The plan defers benchmarking until
  Stages 3-5 can drive realistic load through the merger;
  benching in isolation would over-fit a synthetic generator.

Risks called out in the plan that survived implementation:

- **Lifetime ergonomics at the call site.** The caller must hold
  `Vec<PspReader>` and `Vec<RecordsIter<'_, _>>` in parallel for
  the merger's lifetime. The path-based opener helper (above)
  would hide this; it's deliberately not built yet.
- **`PerPositionPileups: Clone` cost in tests only.** Production
  paths move the value out — confirmed by inspection (the only
  `.clone()` on emitted items lives in test-side equality
  assertions).
