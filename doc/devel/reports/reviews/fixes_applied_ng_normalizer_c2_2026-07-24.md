# Fix Application Report: ng_normalizer_c2_2026-07-24.md

**Date:** 2026-07-24
**Source review:** `doc/devel/reports/reviews/ng_normalizer_c2_2026-07-24.md`
**Source state reviewed against:** branch `ng-locus-evidence`, uncommitted C2 diff
**Execution mode:** non-interactive (plan-driven)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 5 · Nits: ~3

### Outcome totals
- Applied: 6 · Disputed (won't-fix Nits): 2 · Failed validation: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 2328 passed / 4 ignored; `left_align_structured` 18 → 19
- Performance check → skipped (not a bench path; 1c is 1a in a bounded loop)

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status |
|---|---|---|---|---|
| Mi1 | Minor | `to_fixpoint` misleading `to_*` name | Apply | Applied (→ `drive_to_fixpoint`) |
| Mi2 | Minor | cap has no regression guard | Apply | Applied (compile-time assert) |
| Mi3 | Minor | no already-leftmost 1c test | Apply | Applied |
| Mi4 | Minor | panic message not reproducible | Apply | Applied |
| Mi5 | Minor | agreement test ad-hoc tuple | Apply | Applied (→ `Case`) |
| N1 | Nit | stub name is a verb phrase | Apply | Applied (→ `NonConvergingNormalizer`) |
| N2 | Nit | `FIXPOINT_MAX_ITERATIONS` pub | Dispute | Won't fix (1b symmetry) |
| N3 | Nit | iteration/pass vocab drift | Dispute | Won't fix (distinct concepts) |

## 3. Questions asked and answers
None.

## 4. Per-finding log (condensed)

- **Mi1 Applied:** `to_fixpoint` → `drive_to_fixpoint` (the `to_*` idiom implies a by-value conversion; this mutates in place and can panic).
- **Mi2 Applied:** `const _: () = assert!(FIXPOINT_MAX_ITERATIONS >= 2, …)` — a compile-time floor so a regression lowering the cap below what a one-pass fixpoint needs fails the build, not silently in production.
- **Mi3 Applied:** `an_already_leftmost_input_reaches_the_fixpoint_without_shifting`.
- **Mi4 Applied:** the fixpoint panic now appends the alignment/read/reference (matching the property oracle's diagnostic), so a non-convergence is reproducible from the failure line.
- **Mi5 Applied:** the agreement test uses the shared `Case` struct instead of a positional 4-tuple with two adjacent `&[u8]` fields.
- **N1 Applied:** `NeverStabilizes` → `NonConvergingNormalizer` (types are nouns).
- **N2/N3 Won't fix:** `FIXPOINT_MAX_ITERATIONS` stays `pub` for symmetry with 1b's `pub MAX_PASSES`; the iteration/pass vocabulary difference is correct (a 1c iteration is one whole 1a application; a 1b pass is one shift sweep).

## 5–8.
Deferred: none. Disputed to return to reviewer: none. Failed validation: none. Blocked: none.

## 9. Performance check
Skipped — not a bench path.

## 10–12. Commands / notes
- `cargo fmt`, `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --lib`.
- Note: the review confirmed the fixpoint/panic logic correct and both `#[should_panic]` tests genuine (they cannot pass for the wrong reason).
