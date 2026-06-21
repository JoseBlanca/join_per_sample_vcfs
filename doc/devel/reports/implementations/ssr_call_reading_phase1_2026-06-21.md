# SSR `ssr-call` reading layer — Phases 0–3 (scaffolding → cursor → merger → driver)

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

**Phase 2 — the merger** (commit `412536d`):
- [src/ssr/cohort/merge.rs](../../../src/ssr/cohort/merge.rs) — `CohortMerger<R, C>`,
  a catalog-driven k-way merge as `Iterator<Item = Result<(u64, CohortLocus)>>`
  (sparse-omit, monotonic seq). `from_parts` owns the same-catalog md5 check + the
  chromosome-id reconciliation (cohort-global ids from the shared table; per-file
  remaps by name). `SsrMergeError`.

**Phase 3 — the driver** (commit `8915e5e`):
- [src/ssr/cohort/driver.rs](../../../src/ssr/cohort/driver.rs) — `run` (open +
  merge + write), `write_dump` (generic, in-memory-testable), `format_locus` (pure);
  `SsrCallConfig`, `SsrCallError`. `CohortMerger::open` (file-backed) + `chrom_names`.
  `run_ssr_call` now drives the cohort. **Single-threaded by design** (see Tradeoffs).
  Shared test fixtures lifted to
  [src/ssr/cohort/test_support.rs](../../../src/ssr/cohort/test_support.rs).

## Tests added

- **types (4):** `LocusId` ordering precedence; `CohortLocus` empty/build/parallel
  vectors/ascending-index guard.
- **CLI (3):** parse + required-arg + stub-ok.
- **owning iterator (2):** owning vs borrowing identical across a multi-block fixture;
  SSR-schema-over-SNP-file mismatch → one error then `None`.
- **cursor (9):** coordinate conversion + QC carry-through; Absent before / between /
  after stored loci; exhaustion → permanent `None`; empty-observed locus present;
  multi-block crossing; chrom remap; skip panic; rewind panic.
- **merger (7):** present-sample gather + catalog-order + sparse-omit + dense seqs;
  catalog-frame + evidence carry-through; cross-chromosome ordering; no-inputs; md5
  mismatch; missing md5 param; unknown catalog chrom.
- **driver (5):** `format_locus` render + unknown-chrom fallback; `write_dump` header +
  one row/locus in catalog order; `run` over real temp files; `run` errors on missing
  inputs.

## Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings`
  clean.
- `cargo test --lib` → **1159 passed, 0 failed, 2 ignored** (+30 over the pre-Phase-0
  baseline of 1129).
- `ssr-call --help` renders; `ssr-call` runs end-to-end (temp-file test + binary smoke).

## Tradeoffs / follow-ups

- **Phase 3 driver is single-threaded** — the producer/queue/worker-pool/writer topology
  is deliberately *not* built. It exists to overlap expensive EM work with decode; there
  is no EM yet (genotyping doc), and decode parallelism is the separate Phase-5 prefetch
  pool, so a thread pool around a formatting stub would be unverifiable complexity.
  `--threads` / `--queue-depth` are accepted but reserved.
- **Output is a TSV dump**, not a VCF — a placeholder + reading-layer inspection tool
  until the EM + VCF land (the genotyping doc).
- **Cursor decode is synchronous inline** (no pool); the shared decode-priority pool +
  prefetched futures are Phase 5, profiling-gated (Q-R4↔Q-R6).
- **Remaining:** Phase 4 (two-pass re-read — the pre-pass consumes the merge stream,
  then genotyping re-reads it), Phase 5 (prefetch pool), and the genotyping EM + VCF
  (separate plan).
