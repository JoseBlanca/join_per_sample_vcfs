# Code Review: ng read filtering — Milestone C (the record-source seam)
**Date:** 2026-07-14
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Milestone C diff — the `RawRecord`/`RecordSource` seam + noodles adapters in `src/ng/read/filtering.rs`, and a `pub(super)`→`pub(crate)` lift in `src/bam/alignment_input.rs`
**Status:** Approve-with-changes

---

### 1. Scope
- In-scope: the seam traits, `NoodlesRawRecord`, `BamRecordSource`, and their tests; the one-line visibility change.
- Reviewed against: working tree on `main`, HEAD `0d309fc`.
- Categories dispatched: reliability, errors, naming, idiomatic, refactor_safety, module_structure, smells, extras (8 sub-agents). Verified the noodles contracts against the registry source.

### 2. Verdict
Approve-with-changes. **0 Blocker, 2 Major** (both test-coverage — the seam logic was verified correct against the noodles source), plus Minors/Nits.

### 3. Execution status
- fmt/clippy `--lib` clean; `cargo test --lib -- ng::read::filtering` → 24 pass (after fixes).

### 4. Open questions
1. CRAM deferral (extras Mi) — resolved: recorded as a scope decision (plan C2 + code + PROJECT_STATUS), surfaced to the owner at Checkpoint C.

### 5. Top 3 priorities
1. **M1** — buffer-reuse no-leak invariant untested.
2. **M2** — `BamRecordSource` fatal read-error path untested.
3. **CRAM deferral** documentation + mis-citation.

### 6. Findings

**Major**
- **M1: no buffer-reuse no-leak test.** The real-BAM test reused the buffer across two same-length records and never decoded the second, so a stale-tail leak would go uncaught. (Reliability verified noodles fully overwrites the record, so no defect today — a regression guard against a noodles bump.) **Applied:** `bam_record_source_reuses_the_buffer_without_leaking_a_prior_record` (40-base then 10-base, decode the second). *(reliability)*
- **M2: `BamRecordSource` fatal read-error path untested.** The `Err`-is-fatal contract (`?` propagation) has no malformed-stream test. **Deferred to Milestone D** — robust error injection (a failing reader / truncated stream at the right BGZF offset) fits the iterator/fixture layer, where the error boundary lives; the `?` propagation itself was verified correct. *(reliability)*

**Minor**
- **Mi1: `NoodlesRawRecord` `pub` fields over-expose source-managed state** (`source_file_index` is re-stamped every read; `record` is refilled). **Applied** — fields → `pub(crate)`. *(idiomatic, refactor_safety, module_structure)*
- **Mi2: gratuitous `Clone` on `NoodlesRawRecord`** (deep-copies the `RecordBuf` the buffer exists to avoid copying). **Applied** — dropped `Clone`. *(idiomatic)*
- **Mi3: decode-error test under-asserts** (`is_err()` only). **Applied** — asserts `io::ErrorKind::InvalidData`. *(extras)*
- **Mi4: `Default` bound rationale undocumented.** **Applied** — doc note (the D iterator seeds one reusable buffer). *(refactor_safety)*
- **Mi5: CRAM deferral undocumented upstream + mis-cited spec §2.5.** **Applied** — corrected the code rationale, recorded the scope decision in the plan (C2, D3) and PROJECT_STATUS. *(extras)*
- **Mi6: adapters share the file with the pure cascade.** **Deferred** — split the noodles adapters into a `record_source` submodule when the CRAM source or the D driver lands. *(module_structure)*

**Nits**
- `read_next` `Ok(false)` buffer state unspecified — **Applied** (doc note).
- `FakeSource.next` bare participle — **Applied** (`next_index`).
- redundant `drop(writer)` — **Applied** (scoped the writer in a block).
- test-fixture duplication vs `segment_reader` helpers — **Deferred** (extract at the 3rd BAM-fixture site).

### 7. Out of scope observations
Mid-stream read error lacks the failing file's identity — deferred to D, where `ReadFilter` owns the file and can add context.

### 8. Missing tests added now
- `bam_record_source_reuses_the_buffer_without_leaking_a_prior_record` (M1).
- decode-error kind assertion (Mi3).

### 9. What's good
- The noodles `read_record_buf` contract was verified against the registry source (fully overwrites the buffer, `Ok(0)` = EOF, corruption → `Err`) — the reuse and EOF claims hold.
- The visibility lift is the minimal `pub(crate)` (not `pub`), documented at the definition, with the ng→bam direction and no back-reference.
- The in-memory-BAM round-trip test exercises the real noodles writer *and* reader, so an EOF-convention change would fail it.

### 10. Commands to re-verify
- `./scripts/dev.sh cargo test --lib -- ng::read::filtering`
- `./scripts/dev.sh cargo clippy --lib`
