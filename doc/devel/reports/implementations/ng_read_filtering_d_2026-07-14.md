# Implementation report: ng read filtering — Milestone D (the ReadFilter iterator)

**Date:** 2026-07-14
**Feature:** ng step 1 — read filtering (final milestone; **Step 1 complete**)
**Plan:** [read_filtering.md](../../ng/impl_plan/read_filtering.md) (Milestone D, steps D1–D3)
**Spec / arch:** [spec §2.6/§5](../../ng/spec/read_filtering.md), [arch §3](../../ng/arch/read_filtering.md)

## 1. Plan

Wire the cascade + a record source into the end-to-end driver: `ReadFilter::new`
(fail-fast contig validation), the `Iterator<Item = MappedRead>` pipeline, the
running `ReadFilterCounts`, and the fixture-driven drop-count tests. Carries the
deferred enforcement items.

## 2. Decisions

- **`header()` added to `RecordSource`** so `new` can validate contigs (a small
  trait extension; the arch's illustrative sketch omitted it, but `new`'s
  fail-fast contract needs it). All impls (BAM/CRAM/fake) provide it.
- **Fatal-error model = latch + `finish`.** The item is a bare `MappedRead`
  (spec §5), so a fatal condition (`read_next`/`decode`/#8-fetch error) is latched
  in `fatal_error`, `next()` returns `None` and stays stopped, and
  `finish() -> Result<ReadFilterCounts, ReadFilterError>` surfaces it. A `Drop`
  guard debug-asserts the error was observed (catches the "iterate in a `for`
  loop and drop" misuse in dev/test); `#[must_use]` on `finish`.
- **`ReadFilterCounts::record_drop`** is the exhaustive-`match` enforcement of the
  `DropReason` ↔ counts 1:1 mapping (carried from A/B).

## 3. Changes made

- **`src/ng/read/filtering.rs`** — `ReadFilterError` (Source/Decode/Reference),
  `ReadFilter<S, R>` + `new` + `counts` + `finish` + the `Iterator`/`Drop` impls +
  a `fail()` helper; `ReadFilterCounts::record_drop`; `header()` on the trait +
  impls. Removed the now-obsolete `#[cfg_attr(not(test), allow(dead_code))]` on the
  `verdict_*` functions (now called by `next()`).
- **`src/ng/read/mod.rs`** — re-export the step-1 public surface.

## 4. Tests added (9)

- `read_filter_bam_fixture_matches_hand_counted_drops` /
  `read_filter_cram_fixture_matches_hand_counted_drops` — the **port anchor**: a
  small known BAM / CRAM with two kept reads plus one read per mapped drop reason,
  asserting the exact `ReadFilterCounts`. The bad-CIGAR read is *also* high-mismatch,
  so the exact counts **discriminate the #9-before-#8 order** (it must land in
  `bad_cigar`, not `high_mismatch_fraction`).
- `read_filter_charges_an_unmapped_read_end_to_end` — the #5 counter (kept out of
  the BAM/CRAM fixture because a realistic unmapped read has MAPQ 0 → #2, and a
  fake MAPQ doesn't survive a CRAM round-trip).
- `read_filter_new_rejects_a_contig_missing_from_the_reference` — fail-fast `new`.
- `read_filter_source_read_error_is_fatal` / `_decode_error_is_fatal` /
  `_reference_error_mid_stream_is_fatal` — each `ReadFilterError` variant, latched
  and surfaced by `finish` (and "stays stopped").
- `read_filter_over_an_empty_source_yields_nothing_and_zero_counts`;
  `read_filter_counts_is_a_running_tally_before_exhaustion`.

## 5. Validation

Dev container: `cargo fmt -- --check` (ng clean), `cargo clippy --lib` (clean),
`cargo test --lib -- ng::read::filtering` → **35 tests pass**.

## 6. Tradeoffs and follow-ups

- **Fatal-error observability (from the review Blocker):** because the surface is
  `Item = MappedRead` (spec §5), a bare `for`-loop that never calls `finish` still
  loses a latched fatal error silently in *release* builds. Strongly mitigated
  (Drop debug-assert + `#[must_use]` + docs); the residual is inherent to the spec's
  design and is surfaced to the owner. The real pipeline consumer must call
  `finish`.
- **Deferred:** the `record_source.rs` submodule split (the noodles adapters +
  traits, now ~264 lines) — a pure-move refactor, a tracked follow-up.
