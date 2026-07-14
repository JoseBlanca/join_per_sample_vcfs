# Fix Application Report: ng read filtering — Milestone C

**Date:** 2026-07-14
**Source review:** [ng_read_filtering_c_2026-07-14.md](ng_read_filtering_c_2026-07-14.md)
**Source state reviewed against:** working tree on `main`, HEAD `0d309fc`
**Execution mode:** non-interactive (CRAM deferral surfaced to owner at Checkpoint C)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 2 · Minors: 6 · Nits: 4

### Outcome totals
- Applied: 8 · Deferred: 3 · Disputed: 0

### Validation summary
- `cargo fmt -- --check` (ng) → clean
- `cargo clippy --lib` → clean
- `cargo test --lib -- ng::read::filtering` → 24 tests pass
- Performance check → Skipped (no `benches/` path touched; `read_next` is allocation-free)

### Unresolved high-priority findings
- **M2** deferred to Milestone D (fatal read-error test) — the `?` propagation was verified correct; robust error injection fits the iterator/fixture layer.

## 2. Findings table

| ID | Severity | Title | Final status | Files |
|---|---|---|---|---|
| M1 | Major | buffer-reuse no-leak untested | Applied | `read/filtering.rs` |
| M2 | Major | fatal read-error path untested | Deferred (→ D) | None |
| Mi1 | Minor | `NoodlesRawRecord` pub fields | Applied | `read/filtering.rs` |
| Mi2 | Minor | gratuitous `Clone` | Applied | `read/filtering.rs` |
| Mi3 | Minor | decode-error test under-asserts | Applied | `read/filtering.rs` |
| Mi4 | Minor | `Default` bound undocumented | Applied | `read/filtering.rs` |
| Mi5 | Minor | CRAM deferral undocumented + mis-cited | Applied | `read/filtering.rs`, plan, PROJECT_STATUS |
| Mi6 | Minor | adapters share file with cascade | Deferred (→ CRAM/D) | None |
| Nit-a | Nit | `Ok(false)` buffer state | Applied | `read/filtering.rs` |
| Nit-b | Nit | `FakeSource.next` → `next_index` | Applied | `read/filtering.rs` |
| Nit-c | Nit | redundant `drop(writer)` | Applied | `read/filtering.rs` |
| Nit-d | Nit | test-fixture duplication | Deferred (3rd site) | None |

## 3. Questions asked and answers
None (the CRAM deferral was recorded as a scope decision and surfaced at the Checkpoint-C summary rather than blocking).

## 4. Per-finding log (abbreviated)

- **M1 — Applied.** `bam_record_source_reuses_the_buffer_without_leaking_a_prior_record`: a 40-base record then a 10-base record through the same reused buffer; decoding the second asserts `seq.len() == 10` (no stale tail).
- **M2 — Deferred to D.** The `read_next` `Err` is a one-line `?` propagation, verified correct against the noodles source. A robust malformed-stream test needs error injection at a specific BGZF offset (or a failing reader), which fits the Milestone-D iterator + fixture, where the run-level error boundary lives. Recorded as a D carry.
- **Mi1 — Applied.** `record` and `source_file_index` → `pub(crate)` (source-managed internal state; `NoodlesRawRecord` stays `pub` because it is `BamRecordSource`'s associated `Record` type).
- **Mi2 — Applied.** Dropped `Clone` (kept `Debug, Default`).
- **Mi3 — Applied.** Decode-error test asserts `io::ErrorKind::InvalidData`.
- **Mi4 — Applied.** Documented the `Default` bound (the D iterator seeds one reusable buffer).
- **Mi5 — Applied.** Corrected the code's CRAM rationale (the mis-cited spec §2.5), and recorded the BAM-only scope + CRAM-deferral in the plan (C2, D3) and PROJECT_STATUS; surfaced to the owner at the Checkpoint-C summary.
- **Mi6 — Deferred.** Split the noodles adapters into a `record_source` submodule when the CRAM source or the D driver lands.
- **Nit-a/b/c — Applied.** `Ok(false)` buffer-state doc; `FakeSource.next`→`next_index`; scoped the writer in a block (removing the redundant `drop`).
- **Nit-d — Deferred.** Extract a shared BAM test-fixture module at the 3rd site.

## 5. Deferred findings to carry forward
- M2 (fatal read-error test), Mi6 (submodule split), Nit-d (fixture dedup), plus the carried `DropReason`↔counts enforcement → Milestone D / CRAM follow-up.

## 6–8. Disputed / failed / blocked
None.

## 9. Performance check
Skipped — no `benches/` path touched; `BamRecordSource::read_next` reuses one `RecordBuf` (no per-read allocation).

## 10–11. Commands
- `./scripts/dev.sh cargo fmt -- --check` → clean
- `./scripts/dev.sh cargo clippy --lib` → clean
- `./scripts/dev.sh cargo test --lib -- ng::read::filtering` → 24 pass

## 12. Notes
- The CRAM deferral is the one scope decision of this milestone; it does not change the seam's design (a CRAM `RecordSource` is an additive sibling). Surfaced to the owner at Checkpoint C.
