# Fix Application Report: ng read filtering — Milestone A

**Date:** 2026-07-14
**Source review:** [ng_read_filtering_a_2026-07-14.md](ng_read_filtering_a_2026-07-14.md)
**Source state reviewed against:** working tree on `main`, HEAD `8ab92c7`
**Execution mode:** non-interactive (plan-driven)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 4 · Nits: 5

### Outcome totals
- Applied: 3 · Applied with adaptation: 0 · Already fixed: 0 · Deferred: 1 · Disputed: 2 · Failed validation: 0 · Blocked: 0 · Superseded: 0

(The 5 Nits: 1 Applied — ±INFINITY test; 2 Disputed/keep — `.get()` convention, per-field default docs; 2 no-change — `Copy` size, optional `TryFrom`.)

### Validation summary
- `cargo fmt -- --check` (ng files) → clean
- `cargo clippy --lib` → clean (no warnings on in-scope files)
- `cargo test --lib -- ng::read::filtering ng::types` → 6 in-scope tests pass
- Performance check → Skipped (types-only, no `benches/` path touched)

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `counts_default_is_all_zero` weak | Apply | Applied | `read/filtering.rs` | Pass |
| Mi2 | Minor | `.expect()` marker + cross-module invariant | Apply | Applied | `read/filtering.rs` | Pass |
| Mi3 | Minor | co-dependent config fields | Dispute | Disputed | None | N/A |
| Mi4 | Minor | `DropReason`↔counts unenforced | Defer | Deferred | None | N/A |
| Nit-inf | Nit | missing ±INFINITY test | Apply | Applied | `types.rs` | Pass |
| Nit-doc | Nit | ng/mod.rs header stale (cross-cat) | Apply | Applied | `mod.rs` | Pass |
| Nit-get | Nit | `.get()` on pub-field newtypes | Dispute | Disputed | None | N/A |
| Nit-* | Nit | per-field default docs / Copy size / TryFrom | Defer | Deferred | None | N/A |

## 3. Questions asked and answers
None.

## 4. Per-finding log

### Mi1 — `counts_default_is_all_zero` weak
- **Final status:** Applied. Replaced the tautological `assert_eq!(default, default)` + 3-field spot check with an explicit all-zero struct literal (no `..`) pinning all 10 counters. Doubles as a refactor-safety guard (new field → test must update).

### Mi2 — `.expect()` in `Default`
- **Final status:** Applied. Added the `// PANIC-FREE:` marker and documented that the invariant rests on `DEFAULT_MAX_READ_MISMATCH_FRACTION` (0.10) and is exercised by `default_config_reproduces_the_production_filter_policy`. Moving the range check into the production constant was rejected — out of this step's scope (would edit `crate::bam::alignment_input`).

### Mi3 — co-dependent config fields
- **Final status:** Disputed. The flat `max_read_mismatch_fraction` + `mismatch_bq_floor` layout is the specified type (spec §4, arch §2.2) and a deliberate mirror of production `AlignmentMergedReaderConfig`; the struct doc records the mirror. Nesting would diverge from the spec and the port target — a design change, not a fix.

### Mi4 — `DropReason`↔`ReadFilterCounts` unenforced
- **Final status:** Deferred to Milestone B. The mapping site (tally/`match`) lands with the cascade; A3 is deliberately no-logic. Recorded as a B acceptance criterion (exhaustive `match`/destructure, mirroring production `FilterCounts::record_drop`/`merge`).

### Nit-inf — missing ±INFINITY test
- **Final status:** Applied. Added `f32::INFINITY` / `f32::NEG_INFINITY` rejection assertions to `mismatch_fraction_rejects_out_of_range`.

### Nit-doc — stale ng/mod.rs header
- **Final status:** Applied. Header now lists `types`, `ref_seq`, and the step-1 `read` module.

### Nit-get / per-field docs / Copy / TryFrom
- **Final status:** Disputed (`.get()` convention: mandated by `ng_step_interfaces.md` §1) / Deferred (cosmetic; no material value this milestone).

## 5. Deferred findings to carry forward
- Mi4 — `DropReason`↔`ReadFilterCounts` enforcement → Milestone B.

## 6. Disputed findings to return to reviewer
- Mi3 — flat config layout is spec-mandated (production mirror).
- Nit-get — `.get()` on pub-field newtypes is the mandated uniform-accessor convention.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no Apply touched perf-sensitive code (types-only, off any `benches/` path).

## 10. Commands run
- `./scripts/dev.sh cargo fmt -- ...` / `-- --check`
- `./scripts/dev.sh cargo clippy --lib`
- `./scripts/dev.sh cargo test --lib -- ng::read::filtering ng::types`

## 11. Command results
- fmt → clean (ng files); clippy `--lib` → clean; test → 6 in-scope pass.

## 12. Notes
- Foundations-C `WindowedRefSeq` was committed first (`8ab92c7`) to give this milestone a clean, read-filtering-only diff (per user direction). Unrelated var_calling/paralog/examples changes were left untouched in the tree.
