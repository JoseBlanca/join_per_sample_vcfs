# Fix Application Report: ng_normalizer_b1_2026-07-24.md

**Date:** 2026-07-24
**Source review:** `doc/devel/reports/reviews/ng_normalizer_b1_2026-07-24.md`
**Source state reviewed against:** branch `ng-locus-evidence`, uncommitted B1 diff
**Execution mode:** non-interactive (plan-driven)
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0 · Majors: 2 · Minors: 6 · Nits: 2

### Outcome totals
- Applied: 8 · Applied with adaptation: 0 · Disputed (won't-fix Nit): 1 · Applied (Nit): 1 · Deferred: 0 · Failed validation: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --lib` → 0, 2305 passed / 4 ignored; `left_align_structured` 7 → 13
- Performance check → skipped (a two-line wrapper delegating to an already-benched production function; no new hot path)

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status |
|---|---|---|---|---|
| M1 | Major | `reference_offset` boundary/panic untested | Apply | Applied |
| M2 | Major | no empty-cigar coverage | Apply | Applied |
| Mi1 | Minor | precondition unasserted | Apply | Applied |
| Mi2 | Minor | doc bundles two preconditions | Apply | Applied |
| Mi3 | Minor | already-leftmost not byte-checked | Apply | Applied |
| Mi4 | Minor | no leading-deletion-at-nonzero-offset | Apply | Applied |
| Mi5 | Minor | cross-block propagation indirect | Apply | Applied |
| Mi6 | Minor | back-reference into pileup::walker | Apply | Applied (debt note) |
| N1 | Nit | `as usize` lossy cast | Dispute | Won't fix (ng convention) |
| N2 | Nit | duplicate `Case` struct | Apply | Applied (hoisted) |

## 3. Questions asked and answers
None.

## 4. Per-finding log (condensed)

- **M1/M2/Mi3/Mi4/Mi5 Applied:** six new tests — `an_offset_past_the_reference_panics` (`#[should_panic(expected = "reference_offset")]`; the panic is real in release via the slice's own bounds check, so this is not the false-confidence-on-a-debug-assert case the project warns against), `an_offset_equal_to_the_reference_length_does_not_panic`, `an_empty_cigar_is_left_untouched`, `an_already_leftmost_deletion_comes_back_byte_identical` (idempotence), `a_deletion_rolling_to_the_read_start_keeps_a_nonzero_reference_offset` (the offset-invariance crux), `two_deletions_across_a_block_merge_and_left_align` (cross-block + merge, also added to the byte-parity list).
- **Mi1/Mi2 Applied:** `debug_assert!(offset <= reference.len(), …)` + the `normalize` doc reworded to separate the graceful cigar/read guard from the panicking `reference_offset` bound.
- **Mi6 Applied:** an explicit debt comment on the `use crate::pileup::walker::indel_norm::left_align_indels` import, naming the behavioural back-reference into the pipeline stage (mirrors the `CigarOp` note; deliberate under reuse-over-rewrite; production frozen).
- **N2 Applied:** hoisted the duplicated local `Case` struct to a single shared one.
- **N1 Won't fix:** `reference_offset as usize` kept for consistency with A2 (`leftmost_property.rs`) and ng's one-width `as usize` convention; lossless on the 64-bit target, and the new `debug_assert` bounds it.
- **Name verdict:** `StructuredLeftAligner` kept (both name-reviewers concur: standard genomics vocabulary, family-consistent with `RepeatedLeftAligner`).

## 5–8.
Deferred: none. Disputed to return to reviewer: none. Failed validation: none. Blocked: none.

## 9. Performance check
Skipped — wrapper only, no new perf-sensitive path.

## 10–12. Commands / notes
- `cargo fmt`, `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test --lib`.
- Note: the property oracle grading the port passed on every fixture — the A2 oracle and production's `left_align_indels` agree, so there is no oracle-vs-production divergence to record for these inputs.
