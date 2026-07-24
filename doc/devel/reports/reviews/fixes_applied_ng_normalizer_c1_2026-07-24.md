# Fix Application Report: ng_normalizer_c1_2026-07-24.md

**Date:** 2026-07-24
**Source review:** `doc/devel/reports/reviews/ng_normalizer_c1_2026-07-24.md`
**Source state reviewed against:** branch `ng-locus-evidence`, uncommitted C1 diff
**Execution mode:** non-interactive (plan-driven)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 3 · Minors: 3 · Nits: ~6

### Outcome totals
- Applied: 8 · Applied as documented-difference: 1 (M3) · Disputed (won't-fix Nits): 3 · Failed validation: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 2322 passed / 4 ignored; `left_align_repeated` 12 → 17
- Performance check → skipped (not on a bench path)

### Unresolved high-priority findings
- **M3 (complex-indel trimming)** — resolved as a *documented difference* + a pinned test, **flagged for the owner at Checkpoint C**: whether 1b should also trim (for full 1a parity) is a scope decision. Not a defect that blocks the commit; the code is correct and honest as written.

## 2. Findings table

| ID | Severity | Title | Decision | Final status |
|---|---|---|---|---|
| M1 | Major | report not `#[must_use]` | Apply | Applied |
| M2 | Major | canonicalize skipped on no-shift path | Apply | Applied |
| M3 | Major | no D/I parsimony trim | Document + flag | Applied (doc + test); owner decision pending |
| Mi1 | Minor | `is_indel` non-exhaustive | Apply | Applied |
| Mi2 | Minor | `with_max_passes(0)` silent | Apply | Applied |
| Mi3 | Minor | helper duplication under-documented | Apply | Applied |
| N* | Nit | `# Panics`, `UNREACHABLE:`, `result`→`shiftable` | Apply | Applied |
| N* | Nit | `as usize`, `nonzero` var, `=/X` untested | Dispute | Won't fix (recorded) |

## 3. Questions asked and answers
None (M3 flagged for the Checkpoint C pause rather than asked mid-step).

## 4. Per-finding log (condensed)

- **M1 Applied:** `#[must_use = "an ExhaustedCap result is not known to be leftmost; handle it or discard it explicitly"]` on `ConvergenceReport`. The trait `normalize`'s `let _ =` is now the documented single sanctioned discard.
- **M2 Applied:** `shift_pass` builds the shifted cigar unconditionally, canonicalizes, and returns `moved = canonical != *cigar`. Existing pass-count semantics preserved (already-leftmost-canonical → passes 1; k-shift → passes k+1); a non-canonical-but-leftmost input now canonicalizes in one extra pass. Verified by `an_already_leftmost_but_non_canonical_input_is_consolidated` and the new non-canonical agreement case (1b == 1a → `[I3, M2]`).
- **M3 Applied (documented difference):** corrected the `canonicalize` doc and module doc to state 1b does not trim an overlapping D/I; pinned by `a_complex_indel_overlap_is_not_trimmed_and_may_differ_from_1a_in_spelling` (1b `[D2,I1,M3]` vs 1a `[D1,M4]`, both leftmost). Reasoning: full trim reimplements `normalize_alleles` trim+shift interleaving (large, risky, beyond C1's shift+cap+report scope); on real reads complex indels are rare, so the residual 1a/1b difference is negligible there; the difference is itself something the D1 screen legitimately surfaces. **Flagged for owner at Checkpoint C.**
- **Mi1 Applied:** `is_indel` rewritten as an exhaustive match.
- **Mi2 Applied:** `with_max_passes` doc note + a pinning test.
- **Mi3 Applied:** a shared debt comment on the duplicated CIGAR helpers pointing at their frozen-production twins.
- **Nits Applied:** `# Panics` heading; `UNREACHABLE:` casing; `result` → `shiftable`.
- **Nits Won't fix:** `reference_offset as usize` (ng one-width convention, consistent with A2/B1; the `debug_assert` bounds it); `nonzero` local (mirrors `build_cigar`, debt note points there); `SeqMatch`/`SeqMismatch` untested (unsettled producer op-subset; current fixtures/producers use `Match`).

## 5. Deferred findings to carry forward
- **M3** — owner decision on full parsimony trimming for 1b (Checkpoint C).

## 6–8.
Disputed to return to reviewer: none. Failed validation: none. Blocked: none.

## 9. Performance check
Skipped — not on a bench path.

## 10–12. Commands / notes
- `cargo fmt`, `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --lib`.
- Note: the 50k-case fuzz that surfaced M2/M3 was the reviewer's, run then reverted; it confirmed no read-consumption or convergence corruption. `the_output_always_consumes_exactly_the_read` pins that invariant against future regression.
