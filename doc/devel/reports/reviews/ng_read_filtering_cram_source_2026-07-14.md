# Code Review + Fixes: ng read filtering — CRAM RecordSource
**Date:** 2026-07-14
**Reviewer:** rust-code-review skill (orchestrator, focused 4-category pass)
**Scope:** `CramRecordSource<R>` + its tests in `src/ng/read/filtering.rs`
**Status:** Approve-with-changes → fixes applied

---

## Scope & dispatch
Focused review of the CRAM source (a parallel sibling to the already-reviewed `BamRecordSource`). Categories: reliability, refactor_safety, extras, idiomatic (the container-decode replication + reference-threading being the real risks). Reviewed against the working tree at HEAD `a93d52d`.

## Verdict
Approve-with-changes. **1 Blocker, 2 Major, several Minor/Nits** — all test-coverage or clarity; the container-decode logic was verified faithful to noodles' own `Records::read_container_records` line-by-line.

## Findings and resolution

**Blocker**
- **B1 (reliability): multi-container spanning untested.** The initial CRAM test had 2 records, which fit in ONE container (noodles packs 10,240/container), so the drain-and-refill loop — the whole reason `CramRecordSource` exists — never ran; a broken refill would silently drop every record past container 1. **Applied:** `cram_record_source_matches_noodles_records_across_multiple_containers` (10,241 records → ≥ 2 containers), asserting full-count + exact equality to noodles' own iterator.

**Major**
- **M1 (refactor_safety): re-implemented noodles-internal logic is brittle.** `refill_from_next_container` copies noodles' private container-decode state machine; a noodles bump could change semantics without a compile error. **Applied:** the multi-container test doubles as an **equivalence guard** vs noodles' own `reader.records()` — a divergence now fails a test.
- **M2 (extras): reference-diff decode + sequence content untested.** The old test asserted only `seq.len()`, and all-'A'-vs-all-'A' emits no CRAM substitution features. **Applied:** a record now carries non-reference bases (`C`,`G`) and the decoded sequence is asserted byte-for-byte.

**Minor / Nits (all applied)**
- Inverted `Ok(true)=EOF` bool on `refill_from_next_container` → replaced with a named `ContainerRefill { Decoded, EndOfInput }` enum (idiomatic).
- Double-collect + `flatten` → a single `for`-loop with `push` (idiomatic; one buffer, one pass).
- Field `buffered` (bare participle) → `buffered_records` (naming).
- `new` repo-precedence undocumented → documented (the passed repository is authoritative; the reader's baked-in copy is unused).
- `BamRecordSource` doc still said "CRAM … for a later step" → updated to point at `CramRecordSource`.
- Post-EOF idempotence: surfaced a **real robustness gap** — noodles' `read_container` errors past the CRAM EOF marker, so a second `read_next` after EOF would `Err`. **Fixed** with an `exhausted` latch so EOF is idempotently `Ok(false)` (matching the BAM source); asserted in the test.

**Deferred**
- Empty-container / multi-slice defensive branch (unreachable via the noodles writer) — left in with the equivalence test covering the real paths.
- Fatal read-error test + the `record_source` submodule split → Milestone D.

## Validation
`cargo fmt -- --check` (ng) clean; `cargo clippy --lib` clean; `cargo test --lib -- ng::read::filtering` → **26 pass**.
