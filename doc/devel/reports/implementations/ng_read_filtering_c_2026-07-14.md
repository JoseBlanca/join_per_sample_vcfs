# Implementation report: ng read filtering — Milestone C (the record-source seam)

**Date:** 2026-07-14
**Feature:** ng step 1 — read filtering
**Plan:** [read_filtering.md](../../ng/impl_plan/read_filtering.md) (Milestone C, steps C1–C2)
**Spec / arch:** [spec §2.5/§5/§7](../../ng/spec/read_filtering.md), [arch §3](../../ng/arch/read_filtering.md)

## 1. Plan

The input-edge seam that lets the flag/MAPQ cascade run before decode: the
`RawRecord` / `RecordSource` traits (C1) + a test fake, and the ng-owned noodles
adapters (C2). C1 and C2 land together under Checkpoint C.

## 2. Assumptions / deviations (recorded, not silent)

- **`decode` is fallible** — `RawRecord::decode(&self) -> io::Result<MappedRead>`,
  not the spec's illustrative infallible `MappedRead`. The reused
  `record_buf_to_mapped_read` is fallible; a decode failure (a corrupt record:
  unmapped flag clear yet no position) surfaces as a fatal `Err`, consistent with
  the #8-fetch and `read_next` error model (spec §7), not a panic. This is
  adapting to the reused API's real shape.
- **BAM only; CRAM deferred** — `BamRecordSource` fits the
  one-record-into-a-reused-buffer shape (`read_record_buf`). CRAM does not (noodles
  decodes CRAM at container granularity into owned `Vec<RecordBuf>` and consults a
  reference at decode time), so a CRAM `RecordSource` is a documented sibling for a
  later step. Surfaced at Checkpoint C; recorded in the plan (C2) and the code.
- **Visibility lift** — `record_buf_to_mapped_read` widened `pub(super)`→`pub(crate)`
  (commented at the definition) so the ng adapter can reuse the decode path — the
  ng → existing-code dependency spec §7 (decision a) calls for.

## 3. Changes made

- **`src/ng/read/filtering.rs`** — `RawRecord` (`flag`/`mapq`/`decode`) and
  `RecordSource` (`read_next` into a reused buffer, `Default` buffer bound);
  `NoodlesRawRecord` (wraps a `RecordBuf` + `source_file_index`; `flag`/`mapq` cheap
  field reads, `decode` reuses `record_buf_to_mapped_read`, `0xFF`→`MapQual(0)`);
  `BamRecordSource<R>` (`read_next` via `read_record_buf`, true buffer reuse,
  `Ok(0)`→`Ok(false)` EOF; unfiltered).
- **`src/bam/alignment_input.rs`** — `record_buf_to_mapped_read` visibility widened
  to `pub(crate)` (+ rationale comment).

## 4. Tests added (6)

- `fake_source_drives_the_seam` — a `FakeRecord`/`FakeSource` drive read_next →
  pre-decode → decode → post-decode with no BAM.
- `noodles_raw_record_reads_flag_mapq_and_decodes` — hand-built `RecordBuf` → flag,
  mapq, decoded fields, `source_file_index`.
- `noodles_raw_record_maps_unavailable_mapq_to_zero` — `0xFF` → `MapQual(0)`.
- `noodles_raw_record_decode_errors_on_a_record_with_no_position` — decode Err is
  `io::ErrorKind::InvalidData`.
- `bam_record_source_reads_flag_mapq_and_decodes_through_a_real_bam` — in-memory BAM
  round-trip: two records, pre-decode reads + decode + EOF.
- `bam_record_source_reuses_the_buffer_without_leaking_a_prior_record` — a 40-base
  then a 10-base record through the same buffer; no stale tail.

## 5. Validation

Dev container: `cargo fmt -- --check` (ng clean), `cargo clippy --lib` (clean),
`cargo test --lib -- ng::read::filtering` → **24 tests pass**.

## 6. Tradeoffs and follow-ups

- **CRAM `RecordSource`** — a tracked follow-up (a container-buffer source).
- **Deferred to Milestone D:** a fatal-read-error test (robust error injection fits
  the iterator/fixture layer); splitting the noodles adapters into a
  `record_source` submodule (when CRAM/D lands); the `DropReason`↔counts exhaustive
  match (tally site).
