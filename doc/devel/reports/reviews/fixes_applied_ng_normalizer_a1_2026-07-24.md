# Fix Application Report: ng_normalizer_a1_2026-07-24.md

**Date:** 2026-07-24
**Source review:** `doc/devel/reports/reviews/ng_normalizer_a1_2026-07-24.md`
**Source state reviewed against:** branch `ng-locus-evidence`, uncommitted A1 diff
**Execution mode:** non-interactive (plan-driven)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 0 · Minors: 2 · Nits: 2

### Outcome totals
- Applied: 1 · Applied with adaptation: 0 · Already fixed: 0 · Deferred: 0 · Disputed: 1 · Failed validation: 0 · Blocked: 0 · Superseded: 0 · Awaiting answer: 0
- Nits: 2 won't-fix (convention).

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 2264 passed / 4 ignored
- `cargo doc` / `cargo audit` → not run (no doc/dep surface change; deferred to milestone host verify)
- Performance check → skipped (no perf-sensitive code; trait declaration + test-only anchors)

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| Mi1 | Minor | `LeadingDeletionStripper` guard/empty branches untested | Apply | Applied | `src/ng/alignment/mod.rs` | Pass |
| Mi2 | Minor | Trait lacks `# Examples` doc test | Dispute | Disputed | None | N/A |
| Nit-idiom | Nit | Swappable `&[u8]` params | Dispute | Disputed | None | N/A |
| Nit-xcat | Nit | Unchecked `+=` in test anchor | Dispute | Disputed | None | N/A |

## 3. Questions asked and answers
None.

## 4. Per-finding log

### Mi1 — guard/empty branches untested
- **Final status:** Applied. Added two `#[test]`s: `leading_deletion_stripper_leaves_a_non_leading_deletion_alignment_untouched` and `leading_deletion_stripper_leaves_an_empty_cigar_untouched`, both driven through `normalize_via`. Verified they pass; verified they fail if the guard is widened (the non-deletion test) / removed (the empty test) — the point of the anchor.
- **Files changed:** `src/ng/alignment/mod.rs`.

### Mi2 — missing `# Examples` doc test
- **Final status:** Disputed (recorded exemption). The two sibling aligner traits carry no `# Examples` doc test; documenting via prose + `#[cfg(test)]` anchors is the module's established pattern, and A1 ships no implementor to write an example against without restating anchor code. Consistent-with-siblings, not a design change. The reviewer explicitly offered this exemption.

### Nit-idiom — swappable `&[u8]` params
- **Final status:** Disputed. Mirrors `MarginalAligner::marginal_probability` and the spec §7 bare-`&[u8]` mandate; a local change would diverge from the trait family.

### Nit-xcat — unchecked `+=` in anchor
- **Final status:** Disputed. `#[cfg(test)]` with controlled small inputs; no shipping path.

## 5–8.
Deferred: none. Disputed to return to reviewer: none (all adjudicated). Failed validation: none. Blocked: none.

## 9. Performance check
Skipped — no Apply touched perf-sensitive code.

## 10–12. Commands / notes
- Commands: `cargo fmt`, `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --lib`.
- Note: impl report for A1 is folded into the Milestone A report (per prior ng practice, per-milestone impl reports).
