# Code Review: ng read filtering — Milestone B (the two-phase cascade)
**Date:** 2026-07-14
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Milestone B working-tree diff — `verdict_pre_decode` / `verdict_post_decode` + tests in `src/ng/read/filtering.rs`
**Status:** Approve-with-changes

---

### 1. Scope
- Reviewed: the two-phase filter cascade (the decision logic, pure) and its unit tests.
- Reviewed against: working tree on `main`, HEAD `e3a26e6`.
- In-scope file: `src/ng/read/filtering.rs` (`verdict_pre_decode`, `verdict_post_decode`, the `// ----- verdict_*` test blocks).
- Out of scope: the Milestone-A types; `crate::bam::alignment_input` reuse targets (`classify_pre_decode`, `read_exceeds_mismatch_fraction`, `cigar_is_bad`, `cigar_ref_span`, `MappedRead`, `FLAG_*`) — verified against, not reviewed; `ref_seq.rs`; the `ReadFilter` iterator (Milestone D).
- Categories dispatched: reliability, errors, naming, idiomatic, refactor_safety, smells, extras (7 sub-agents). Skipped: defaults (reviewed in A), module_structure (single-file), unsafe_concurrency, tooling.

### 2. Verdict
Approve-with-changes. **1 Major** (an error-model design decision, resolved with the owner), convergent Minors/Nits.

### 3. Execution status
- `cargo fmt -- --check` (ng), `cargo clippy --lib`, `cargo test --lib -- ng::read::filtering` — clean/passing. The two `verdict_*` fns carry `#[cfg_attr(not(test), allow(dead_code))]` until Milestone D wires them.

### 4. Open questions and assumptions
1. **M1 error model** (below) — resolved by the owner: keep the fatal model.

### 5. Top 3 priorities
1. **M1** — `OutOfBounds` on #8's fetch is fatal (aborts the run) vs production's skip-and-keep. *Owner decision: keep fatal.*
2. **Casts** — `read.pos as u32` (u64→u32) silently truncates on the (impossible) >4.29 Gbp position.
3. **Test coverage** — the `OutOfBounds` fatal path, and the full pre-decode cascade order, were unpinned.

### 6. Findings

**Major**
- **M1: src/ng/read/filtering.rs — #8 `OutOfBounds` fetch is fatal.** A read whose reference window runs past a contig's 3′ end makes #8's `fetch_raw_into` return `RefSeqError::OutOfBounds`, which is propagated as fatal and aborts the whole run; production (`read_processor.rs`) skips #8 and keeps such a read. *Resolution (owner, 2026-07-14): keep the fatal model* — a validly-aligned read cannot cover positions the contig lacks, so `OutOfBounds` signals a malformed record and is treated like corrupt input. Applied: a test (`high_mismatch_fetch_past_contig_end_is_fatal`) and a doc note making the decision explicit. *(Categories: reliability, extras)*

**Minor**
- **Mi1: `read.pos as u32` / `read.ref_id as u32`.** Unchecked narrowing; `read.pos` is `u64`. Silent truncation would fetch the wrong window and mis-verdict the read. **Applied:** `u32::try_from(...).expect(...)` with a `PANIC-FREE` note (fails loudly on the impossible, consistent with the fatal model). *(Categories: reliability, errors, idiomatic, smells)*
- **Mi2: `DropReason`↔`ReadFilterCounts` still unenforced.** No exhaustive `match` yet. **Deferred to Milestone D** (the tally site). *(Categories: refactor_safety)*
- **Mi3: `verdict_pre_decode`/`verdict_post_decode` are noun-named.** naming.md prefers verbs. **Disputed/keep** — the arch doc §3 mandates these exact names. *(Categories: naming)*
- **Mi4: stale module doc header** ("milestone A lands types only"). **Applied** — header now reflects A+B. *(Categories: extras)*
- **Mi5: pre-decode order only partially pinned.** **Applied** — added `pre_decode_charges_filters_in_full_cascade_order` pinning the full total order. *(Categories: refactor_safety)*
- **Mi6: `mismatch_bq_floor` co-dependent field.** **Disputed/keep** — deliberate flat mirror of production `AlignmentMergedReaderConfig` (spec §4), same as Milestone A. *(Categories: idiomatic)*

**Nits**
- Redundant `max_read_mismatch_fraction: None` in the too-short test (**Applied** — removed).
- Poly-A reference literal repeated 7× (**Applied** — `poly_a_ref(n)` helper).
- Test helpers `pre`/`post`/`mapped` slightly terse (kept — section-commented, clear).
- Test configs use `..default()` (kept — acceptable for tests).

### 7. Out of scope observations
The reused-buffer no-stale-bytes assertion belongs to the `ReadFilter` pass (Milestone D, plan D3). A differential/golden test tying `verdict_pre_decode` to production `classify_pre_decode` is best done via the D fixture (real records), not a hand-built `RecordBuf` here.

### 8. Missing tests added now
- `pre_decode_charges_filters_in_full_cascade_order` — full #1→#6 total order.
- `high_mismatch_fetch_past_contig_end_is_fatal` — the `OutOfBounds` fatal path.

### 9. What's good
- The #7→#9→#8 ordering choice is pinned by two attribution tests, each with a companion assertion proving the losing filter *would* have fired — so the ordering intent is regression-protected.
- `FLAG_*`/`DEFAULT_*` and the predicates are used by name (imported), so a production change flows through rather than silently drifting from a re-hardcoded copy.
- #8 allocates nothing per read (reused `ref_buf`) and a read dropped at #7/#9 or with #8 disabled pays no reference work.
- The extras reviewer independently confirmed the left-alignment tally-invariance claim is mathematically sound (legal shifts only re-pair equal bases).

### 10. Commands to re-verify
- `./scripts/dev.sh cargo test --lib -- ng::read::filtering`
- `./scripts/dev.sh cargo clippy --lib`
