# Fix Application Report: pileup_2026-05-06.md

**Date:** 2026-05-06
**Source review:** `ia/reviews/pileup_2026-05-06.md`
**Source state reviewed against:** `b5e1cba` (main; review committed)
**Execution mode:** interactive
**Overall status:** In progress

---

## 1. Executive summary

### Review totals
- Blockers: 2 (B1, B2)
- Majors: 5 (M1, M2, M3, M4, M5)
- Minors: 6 (Mi1, Mi2, Mi3, Mi4, Mi5, Mi6)
- Nits: ~5 small items

### Outcome totals (running)
- Applied: 0
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` — pending
- `cargo clippy --all-targets --all-features -- -D warnings` — pre-existing failures in unrelated files; pileup-scoped rerun pending
- `cargo test --all-targets --all-features` — pileup tests pass; full suite uses unrelated bench fixture
- `cargo doc --no-deps` — not required for this scope
- `cargo audit` — not required for this scope

### Unresolved high-priority findings
*(populated as findings progress)*

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | REF re-folded once per walker step inside open record | Apply | Applied | No | `src/per_sample_caller/pileup/{open_record,walker,tests}.rs` | Pass | No |
| B2 | Blocker | Indel mate-overlap doesn't collapse to one observation | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| M1 | Major | `find_overlapping` early break unsound across heterogeneous spans | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| M2 | Major | Mate-overlap tie-break uses `alignment_start`, not `is_first_mate` | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| M3 | Major | `record_widen_events` counter conflates open and widen | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| M4 | Major | Lifecycle marks attach only to first drained record | Ask | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
| M5 | Major | `windows(2)` only inspects adjacent contributor pairs | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| Mi1 | Minor | `checked_sub.unwrap_or(0)` masks precondition violation | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| Mi2 | Minor | `walker_pos == 0` guard is dead code | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| Mi3 | Minor | `flush_chromosome` accepts unused `_fasta` | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| Mi4 | Minor | `from_slot_counters` awkward shape | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| Mi5 | Minor | `FiveScalars::zero` duplicates `Default::default` | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| Mi6 | Minor | `widen` always appends `extra_bases` to every allele | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |
| Nits | Nit | Clippy errors + small style items | Apply | _pending_ | No | _pending_ | _pending_ | _pending_ |

## 3. Questions asked and answers

*(populated as Asks are issued)*

## 4. Per-finding log

### B1 — REF re-folded once per walker step inside open record

- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Defect confirmed by a regression test added test-first. The review's suggested fix (split "affected = newly opened/widened" + "fold every active read once") was viable but more invasive; an equivalent narrower fix tracks per-record `(read_id → previous contribution)` and subtracts the prior contribution before adding the new one on every fold, so each `(record, read)` pair contributes exactly one net observation across the record's lifetime. The narrower shape preserves the existing `process_position` flow (events at walker_pos drive find_overlapping → fold) while removing the multiplicative inflation.
- **Implementation summary:**
  - `OpenPileupRecord` gained an `AHashMap<u32, FoldedReadState>` field tracking each contributing read's last folded `(allele_index, contribution)`.
  - `process_position` now subtracts a read's prior contribution from its old bucket before adding to the new bucket (sub-/add-helpers added).
  - `ReadContribution` gained a `read_id` field plumbed from the active set in `walker.rs::process_position`.
  - Old `find_or_create_allele(&mut OpenPileupRecord) -> &mut OpenAllele` was replaced by `find_or_create_allele_index(&mut OpenPileupRecord) -> usize`, freeing the borrow so callers can keep mutating the record.
- **Review suggestion used verbatim?:** No.
- **Adaptation:** Subtract-old/add-new instead of "split affected and fold once". Equivalent observable behavior; smaller diff.
- **Verification performed:** Added regression test `deletion_record_does_not_double_count_ref_reads` that fails before fix (REF.num_obs=4) and passes after (REF.num_obs=1, fwd=1). Ran the full pileup suite (50 tests pass).
- **Files changed:**
  - `src/per_sample_caller/pileup/open_record.rs`
  - `src/per_sample_caller/pileup/walker.rs`
  - `src/per_sample_caller/pileup/tests.rs`
- **Tests added or modified:** Added `deletion_record_does_not_double_count_ref_reads` (tests.rs); updated `find_or_create_allele_returns_same_bucket_on_match` to use the renamed `find_or_create_allele_index`.
- **Validation:**
  - `cargo test --lib pileup` → exit 0, 50 passed
  - `cargo fmt --check` → exit 0 (after one auto-format pass)
- **User input:** None.
- **Follow-up:** Phantom-slot residual on a read that switches buckets across compound-event widens (a slot may persist on the bucket the read left, since other reads sharing the slot might still be there). Tests do not currently exercise this; flagged in §11.
- **Residual risk:** Low for the SNP / single-deletion / pure-REF mixes the existing tests cover. Compound events on a single read crossing walker steps are not exercised; a future test should verify the bucket-switch path.

## 5. Deferred findings to carry forward

*(populated at end)*

## 6. Disputed findings to return to reviewer

*(populated at end)*

## 7. Failed-validation findings

*(populated at end)*

## 8. Blocked-by-context-mismatch findings

*(populated at end)*

## 9. Commands run

*(running log)*

## 10. Command results

*(running log)*

## 11. Notes

- Pre-existing clippy errors in `variant_grouping.rs` and `decompression_pool.rs` are out of scope per the source review's §7 ("Out of Scope Observations"). Pileup-scoped clippy will be checked module-locally.
- Tests are run inside the project's dev container via `./scripts/dev.sh` (see `CLAUDE.md`).
