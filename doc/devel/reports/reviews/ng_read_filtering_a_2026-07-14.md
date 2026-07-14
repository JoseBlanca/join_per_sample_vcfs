# Code Review: ng read filtering — Milestone A (types + scaffold)
**Date:** 2026-07-14
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Milestone A working-tree diff — `src/ng/read/` (new), `src/ng/types.rs` additions, `src/ng/mod.rs` wiring
**Status:** Approve-with-changes

---

### 1. Scope
- Reviewed: uncommitted Milestone-A diff (pure type foundations, no filtering logic).
- Reviewed against: working tree on `main`, HEAD `8ab92c7`.
- In-scope files: `src/ng/types.rs` (additions), `src/ng/read/mod.rs`, `src/ng/read/filtering.rs`, `src/ng/mod.rs`.
- Out of scope: `ref_seq.rs`/`WindowedRefSeq` (committed `8ab92c7`); unrelated var_calling/paralog/examples tree changes; the `crate::bam::alignment_input` reuse target (verified, not reviewed).
- Categories dispatched: reliability, naming, defaults, idiomatic, errors, refactor_safety, module_structure, smells (8 parallel sub-agents). `unsafe_concurrency`/`tooling`/`extras` skipped (no unsafe/concurrency, no `Cargo.toml` change, no parser/hot-path/public-crate surface in a types-only diff).

### 2. Verdict
Approve-with-changes. No Blocker, no Major. Convergent Minors/Nits only.

### 3. Execution status
- `cargo fmt` (ng files), `cargo clippy --lib`, `cargo test --lib -- ng::read::filtering ng::types` — all clean/passing (see impl report). No "Needs verification" findings.

### 4. Open questions and assumptions
None requiring a user decision. The flat-config-vs-nested question (Mi3) is resolved by spec §4 (deliberate production mirror).

### 5. Top 3 priorities
1. **Mi1** — `counts_default_is_all_zero` is tautological + covers 3/10 fields.
2. **Mi2** — `.expect()` in `Default` lacks the `// PANIC-FREE:` audit marker and its invariant lives cross-module.
3. **Mi4** — `DropReason`↔`ReadFilterCounts` 1:1 mapping is doc-only, unenforced (forward-looking; mapping site is Milestone B).

### 6. Findings

**Minor**

- **Mi1: src/ng/read/filtering.rs — `counts_default_is_all_zero` weak.** Opens with `assert_eq!(default, default)` (tautological) and checks only `kept`/`duplicate`/`bad_cigar`. Replace with an explicit all-zero struct literal so every counter is pinned and a new field forces the test to update. *(Categories: reliability, smells)*
- **Mi2: src/ng/read/filtering.rs — `.expect()` in `Default`.** Justified in prose but not marked `// PANIC-FREE:` (not greppable), and the named invariant depends on `DEFAULT_MAX_READ_MISMATCH_FRACTION` in another module. Add the marker; note the cross-module dependency. *(Categories: errors, defaults, reliability)*
- **Mi3: src/ng/read/filtering.rs — co-dependent config fields.** `max_read_mismatch_fraction: Option<_>` + `mismatch_bq_floor` encode one filter; the floor's doc says it is only meaningful when the `Option` is `Some` (invariant-in-a-comment). *Resolution:* deliberate mirror of production `AlignmentMergedReaderConfig` (spec §4); documented on the struct. Keep as-is. *(Categories: idiomatic, smells)*
- **Mi4: src/ng/read/filtering.rs — `DropReason`↔`ReadFilterCounts` unenforced.** Holds today (9↔9, in order) but no exhaustive `match`/destructure guards it yet. Defer to Milestone B (the tally site); B acceptance criterion. *(Categories: reliability, refactor_safety, smells)*

**Nits**
- Per-field default values not announced in field docs (only struct-level).
- `get()` duplicated on `pub`-field newtypes — kept (convention: `ng_step_interfaces.md` §1 mandates uniform `.get()`).
- `ReadFilterConfig` is `Copy` at ~20–24 bytes — acceptable (value-like).
- Optional `TryFrom<f32>` on `MismatchFraction` alongside `try_new`.
- Missing ±INFINITY case in the `MismatchFraction` rejection test.

### 7. Out of scope observations
Unrelated uncommitted tree changes (var_calling/paralog/examples) left untouched. Pre-existing crate-wide gate reds (documented in `PROJECT_STATUS.md`) are not from this diff.

### 8. Missing tests to add now
- `counts_default_is_all_zero` → assert against an explicit all-zero literal (all 10 fields).
- `mismatch_fraction_rejects_out_of_range` → add `f32::INFINITY` / `f32::NEG_INFINITY`.

### 9. What's good
- `ReadFilterConfig::default()` sources every threshold from the reused `DEFAULT_*` constants — no magic numbers, and the port anchor is a real test.
- Newtype convention applied cleanly: unconstrained types `pub` + `.get()`; `MismatchFraction` hides its field behind a checked `try_new` returning a typed error.
- Names match the arch vocabulary exactly (naming: No findings), so nothing renames when the rest of the ng vocabulary lands.
- `#[non_exhaustive]` + thiserror `DomainError` matches the established `RefSeqError` pattern.

### 10. Commands to re-verify
- `./scripts/dev.sh cargo test --lib -- ng::read::filtering ng::types`
- `./scripts/dev.sh cargo clippy --lib`
