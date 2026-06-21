# SSR `ssr-call` reading layer — Phase 0 + Phase 1 (per-sample cursor)

**Date:** 2026-06-21 · **Branch:** `ssr-cohort` · **Skill:** rust-feature-implementation

Implements the first phases of the Stage-2 (`ssr-call`) reading & merge layer.
Spec [ssr_cohort_mark2.md §4.1](../../specs/ssr_cohort_mark2.md); architecture
[ssr_call_reading.md](../../architecture/ssr_call_reading.md); plan
[ssr_call_reading.md](../../implementation_plans/ssr_call_reading.md).

## Plan

Build the reading layer in incremental, independently-tested phases. This report
covers **Phase 0** (scaffolding) and **Phase 1** (the per-sample cursor, split into a
psp enabler + the cursor itself). Phases 2–5 (merger, driver, two-pass, prefetch pool)
are still open.

## Assumptions / decisions recorded

- **`LocusId.chrom_id` is cohort-global**, assigned by the merger from catalog
  chromosome order (Phase 2 builds the per-file→global map); the catalog's **0-based
  half-open** frame is canonical. The cursor converts per-file records into it.
- **Coordinate inversion:** container records are 1-based (`[start,end)`); the Stage-1
  writer added `+1` to both bounds, so the cursor subtracts `1` from both — verified
  against [driver.rs:210-211](../../../src/ssr/pileup/driver.rs#L210-L211).
- **`observed → seq_counts`:** the container/Stage-1 field name `observed` is kept; the
  cursor's decode adapter renames into `SampleEvidence.seq_counts` (cursor-local).
- **`evidence_at` returns `Result<Option<…>>`, not the arch doc's bare `Option`** — a
  lazy block refill can fail, so decode errors must propagate; `Option` stays for
  Present/Absent, the monotonic-order violation stays a `debug_assert`/panic.
- **Same-catalog md5 check deferred to the merger** (Phase 2), which holds the catalog
  and validates all inputs before building cursors — the cursor stays single-sample.
- **Owning iterator added rather than refactoring `RecordsIter`** — the borrowing
  iterator can't be stored beside the reader it reads (self-referential). A separate
  owning type leaves the production SNP path byte-for-byte untouched (lowest risk in a
  byte-identity-critical file).

## Changes

**Phase 0 — scaffolding** (commit `a0babd2`):
- [src/ssr/cohort/types.rs](../../../src/ssr/cohort/types.rs) — `LocusId` (global chrom
  id + 0-based half-open, lexicographic `Ord`), `SsrQc`, `SampleEvidence`, sparse-SoA
  `CohortLocus` (`present: Vec<u32>` ∥ `samples`) + builder.
- [src/pop_var_caller/ssr_call.rs](../../../src/pop_var_caller/ssr_call.rs) —
  `SsrCallArgs`/`SsrCallCliError`/`run_ssr_call` stub, wired into the CLI + dispatch.

**Phase 1a — owning typed iterator** (commit `30f30ec`):
- [src/psp/reader.rs](../../../src/psp/reader.rs) — `OwnedRecordsIter<R, S>` +
  `PspReader::into_records_of::<S>()`; mirrors `RecordsIter`'s sequential state machine
  but owns the reader, decoding one block at a time.

**Phase 1b — the cursor** (commit `4bfeb0d`):
- [src/ssr/cohort/reader.rs](../../../src/ssr/cohort/reader.rs) —
  `SampleEvidenceCursor<R>` (`held` + `last_query` guard; `evidence_at`; `advance`;
  `adapt`) + `SsrCohortReadError`.

## Tests added

- **types (4):** `LocusId` ordering precedence; `CohortLocus` empty/build/parallel
  vectors/ascending-index guard.
- **CLI (3):** parse + required-arg + stub-ok.
- **owning iterator (2):** owning vs borrowing identical across a multi-block fixture;
  SSR-schema-over-SNP-file mismatch → one error then `None`.
- **cursor (9):** coordinate conversion + QC carry-through; Absent before / between /
  after stored loci; exhaustion → permanent `None`; empty-observed locus present;
  multi-block crossing; chrom remap; skip panic; rewind panic.

## Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings`
  clean.
- `cargo test --lib` → **1147 passed, 0 failed, 2 ignored** (+18 over the pre-Phase-0
  baseline of 1129).
- `ssr-call --help` renders.

## Tradeoffs / follow-ups

- **Phase 1 cursor is synchronous inline decode** (no pool); the shared decode-priority
  pool + prefetched futures are Phase 5, profiling-gated (Q-R4↔Q-R6).
- **Next — Phase 2:** the catalog-driven k-way merger (builds the per-file→global chrom
  map, the same-catalog md5 check, and emits one `CohortLocus` at a time).
- Then Phase 3 (driver + worker stub + seq writer), Phase 4 (two-pass re-read), Phase 5
  (prefetch pool).
