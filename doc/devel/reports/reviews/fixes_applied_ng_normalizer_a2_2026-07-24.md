# Fix Application Report: ng_normalizer_a2_2026-07-24.md

**Date:** 2026-07-24
**Source review:** `doc/devel/reports/reviews/ng_normalizer_a2_2026-07-24.md`
**Source state reviewed against:** branch `ng-locus-evidence`, uncommitted A2 diff
**Execution mode:** non-interactive (plan-driven)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1 · Majors: 4 · Minors: 2 · Nits: 3 (+ challenge-pass gaps)

### Outcome totals
- Applied: 7 · Applied with adaptation: 0 · Already fixed: 0 · Deferred: 0 · Disputed (won't-fix Nits): 3 · Failed validation: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 2292 passed / 4 ignored; `leftmost_property` 13 → 20 tests
- Performance check → skipped (test-only `#[cfg(test)]` module, no perf-sensitive path)

### Unresolved high-priority findings
None. The oracle's logic was confirmed correct by two independent traces; all fixes hardened its test suite (the ng "tests that can't fail" lesson) or improved its diagnostics/API.

## 2. Findings table

| ID | Severity | Title | Decision | Final status |
|---|---|---|---|---|
| B1 | Blocker | `reference_offset` untested (shift-invariant fixture) | Apply | Applied |
| M1 | Major | soft-clip read cursor untested | Apply | Applied |
| M2 | Major | multi-base insertion far-edge untested | Apply | Applied |
| M3 | Major | indel-arm `left_neighbour` reset untested (false-positive risk) | Apply | Applied |
| M4 | Major | match-required clause not independently pinned | Apply | Applied |
| Mi1 | Minor | opaque OOB panic on malformed alignment | Apply | Applied |
| Mi2 | Minor | `kind: &'static str` stringly-typed | Apply | Applied |
| N1 | Nit | `as usize` not width-change-protected | Dispute | Won't fix (ng convention) |
| N2 | Nit | two near-identical indel arms | Dispute | Won't fix (below threshold) |
| N3 | Nit | `fn alignment` factory is a noun | Dispute | Won't fix (test idiom) |

## 3. Questions asked and answers
None.

## 4. Per-finding log (condensed)

- **B1 Applied:** replaced `the_reference_offset_is_respected` with `a_shiftable_deletion_is_found_using_the_reference_offset_not_the_stretch_start` — `CCCCGAAAAT` @ offset 4; hand-verified that reading from `reference[0..]` lands the crossed column on `C` and flips the verdict.
- **M1 Applied:** `a_soft_clip_advances_the_read_cursor_so_the_crossed_read_base_is_correct` — `XXGCAAT`/`GCAAAT`; dropping the clip advance lands the crossed read base on `C`. Kept the clip-not-aligned second assertion.
- **M2 Applied:** `a_rightmost_period_two_insertion_is_detected_as_shiftable` — `Insertion(2)`; a `read[read_position]` mutation (dropping `+ n - 1`) reads `C` and misses it.
- **M3 Applied:** `an_indel_abutting_a_preceding_indel_is_not_treated_as_slidable` — two sub-cases (`D D`, `I D`) after a mismatch M; each false-positives if its arm's `left_neighbour_is_aligned = false` reset is dropped.
- **M4 Applied:** `a_deletion_whose_far_edge_matches_but_left_column_mismatches_is_leftmost` — `TCA`/`TA`; far edge matches, crossed column mismatches, so only clause 1 keeps it leftmost.
- **Mi1 Applied:** `consumed_read`/`consumed_reference` helpers + two `debug_assert!`s + a documented well-formedness precondition on `find_shiftable_indel`.
- **Mi2 Applied:** `IndelKind { Insertion, Deletion }` enum with `label()`; `ShiftableIndel.kind: IndelKind`; panic uses `.label()`; test assertions use the variants.
- **Extra challenge-pass tests applied:** `an_empty_cigar_is_leftmost`, `a_trailing_insertion_at_the_read_end_is_evaluated_in_bounds`, `a_skip_advances_the_reference_cursor`, `a_second_deletion_uses_the_first_deletions_reference_advance`.
- **N1/N2/N3 Won't fix:** recorded above; all Nits the reviewers themselves deprioritized, and each would diverge from an established ng convention/idiom.

## 5–8.
Deferred: none. Disputed to return to reviewer: none. Failed validation: none. Blocked: none.

## 9. Performance check
Skipped — `#[cfg(test)]` module.

## 10–12. Commands / notes
- `cargo fmt`, `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --lib`.
- Note: every discriminating fixture was hand-traced to confirm the named mutation actually flips the verdict — the fixes are not just "more tests" but tests that provably fail under the mutation they target.
