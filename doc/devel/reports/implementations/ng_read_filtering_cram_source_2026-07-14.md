# Implementation report: ng read filtering — CRAM RecordSource

**Date:** 2026-07-14
**Feature:** ng step 1 — read filtering (Milestone C follow-on)
**Context:** The owner asked for CRAM support after Milestone C landed BAM-only (the CRAM `RecordSource` had been deferred as a sibling). This adds it.

## 1. Plan

Add `CramRecordSource<R>`, a `RecordSource` sibling to `BamRecordSource`, so read
filtering accepts CRAM input. The seam (`RawRecord`/`RecordSource` traits,
`NoodlesRawRecord`) is unchanged — CRAM is purely a second source impl.

## 2. Assumptions / decisions

- **Container-buffered, not per-record-reuse.** CRAM decodes at *container*
  granularity (a container holds slices whose records decode together against the
  reference), so `CramRecordSource` buffers one container's owned `RecordBuf`s and
  yields them one per `read_next`, decoding the next container when the buffer
  drains. Unlike `BamRecordSource` it does not refill one record in place; the
  container buffer's allocations are what get reused across containers.
- **Replicates noodles' internal `Records::read_container_records`.** noodles'
  own whole-file `records()` iterator borrows the reader mutably, so it can't be
  held in a struct (lending-iterator). The source therefore drives
  `read_container` + `slice.decode_blocks()` + `slice.records(repo, …)` +
  `RecordBuf::try_from_alignment_record` directly, mirroring noodles' logic. It
  passes its own `fasta::Repository` (the reader's is private; only the slice
  decoder needs it) — so a `cram::io::Reader` built without a repository works.
- **EOF is latched (`exhausted`).** noodles' `read_container` is not idempotent
  past the CRAM EOF marker (a second call reads past the end and errors), so an
  `exhausted` flag makes `read_next` return `Ok(false)` idempotently after the
  first EOF — matching the BAM source.

## 3. Changes made

- **`src/ng/read/filtering.rs`** — `CramRecordSource<R>` (fields `reader`,
  `header`, `reference_sequence_repository`, `source_file_index`, the reused
  `container` buffer, `buffered_records`, `exhausted`), `new`, the
  `refill_from_next_container` helper (returning a named `ContainerRefill` enum),
  and the `RecordSource` impl. Updated the `BamRecordSource` doc to point at the
  new sibling. Added `noodles_cram`/`noodles_fasta` imports.

## 4. Tests added (2)

- `cram_record_source_reads_flag_mapq_and_decodes_through_a_real_cram` — a CRAM
  round-trip via the repo's `build_cram`/`build_fasta` helpers; a record carries
  non-reference bases (`C`, `G`) so CRAM's substitution-feature decode path is
  exercised, and the decoded sequence is asserted **byte-for-byte**.
- `cram_record_source_matches_noodles_records_across_multiple_containers` —
  10,241 records (> the 10,240/container threshold → ≥ 2 containers), asserting
  the source yields **all** records **and** an exact match against noodles' own
  `reader.records()` iterator, in order. This exercises the drain-and-refill loop
  (the whole reason the CRAM source differs from BAM) and guards against noodles
  container-decode semantic drift. Also asserts post-EOF idempotence.

## 5. Validation

Dev container: `cargo fmt -- --check` (ng clean), `cargo clippy --lib` (clean),
`cargo test --lib -- ng::read::filtering` → **26 tests pass**.

## 6. Tradeoffs and follow-ups

- The container-decode logic mirrors a noodles-internal state machine; the
  equivalence-vs-`records()` test is the guard against a noodles bump changing it
  silently.
- Deferred to Milestone D: a fatal read-error test on the sources; the
  `record_source` submodule split (both sources + traits now justify it).
