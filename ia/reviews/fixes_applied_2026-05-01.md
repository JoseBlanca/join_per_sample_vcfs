# Fix Application Report: per_sample_caller_cram_input_2026-04-29.md

**Date:** 2026-05-01
**Source review:** `ia/reviews/per_sample_caller_cram_input_2026-04-29.md`
**Source state reviewed against:** commit b2729fe (branch main)
**Execution mode:** interactive
**Overall status:** In progress

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 5
- Minors: 6
- Nits: 6

### Outcome totals
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
- `cargo fmt --check` → not run
- `cargo clippy --all-targets --all-features -- -D warnings` → not run
- `cargo test --all-targets --all-features` → not run
- `cargo doc --no-deps` → not run
- `cargo audit` → not run

### Unresolved high-priority findings
- (to be filled in)

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | silent `as u32` truncation on position and length values | Apply | Applied | No | `cram_input.rs`, `errors.rs`, `record_specs.rs`, `cram_files.rs` | Pass | No |
| M2 | Major | `MultipleSampleNames` reused for within-file SM disagreement | Apply | Applied | No | `cram_input.rs`, `errors.rs` | Pass | No |
| M3 | Major | `decode_md5_hex` silently tolerates malformed `M5` | Apply | TBD | No | TBD | TBD | TBD |
| M4 | Major | `peek_head_keys` reports `MalformedRecord` for legitimate kept-unmapped reads | Apply | TBD | No | TBD | TBD | TBD |
| M5 | Major | `OutOfOrderRead` carries packed `(ref_id, pos)` keys | Apply | TBD | No | TBD | TBD | TBD |
| Mi1 | Minor | empty-input check returns `CramInputError::Io` | Apply | TBD | No | TBD | TBD | TBD |
| Mi2 | Minor | defensive `.unwrap_or_default()` on `canonical_path` masks an invariant | Apply | TBD | No | TBD | TBD | TBD |
| Mi3 | Minor | iterator behaviour after `Some(Err)` is unspecified | Apply | TBD | No | TBD | TBD | TBD |
| Mi4 | Minor | `MalformedRecord` errors with empty `qname` context | Apply | TBD | No | TBD | TBD | TBD |
| Mi5 | Minor | `min_read_length` filter pays for full record decode before rejecting | Apply | TBD | No | TBD | TBD | TBD |
| Mi6 | Minor | MAPQ-missing collapses to 0 — needs doc and regression test | Apply | TBD | No | TBD | TBD | TBD |
| N1 | Nit | `OpenCram::path_for_errors` rename | Apply | TBD | No | TBD | TBD | TBD |
| N2 | Nit | `window: VecDeque` → `Vec` | Apply | TBD | No | TBD | TBD | TBD |
| N3 | Nit | WindowKey clone | Apply | TBD | No | TBD | TBD | TBD |
| N4 | Nit | `encode_order_key` removal | Apply (subsumed by M5) | TBD | No | TBD | TBD | TBD |
| N5 | Nit | Implementation file size — won't fix | Won't fix (per review) | Deferred | No | None | N/A | No |
| N6 | Nit | `Container::default()` aesthetics | No concrete change | Deferred | No | None | N/A | No |

## 3. Questions asked and answers

None. The review's Open Questions 1-4 were resolved at review time by Jose; all remaining decisions are in-line in the findings.

## 4. Per-finding log

### M1 — silent `as u32` truncation on position and length values
- **Severity:** Major
- **Initial decision:** Apply (decision: widen to u64)
- **Final status:** Applied
- **Reasoning:** Bug confirmed at the five call sites named in the review. Widening to `u64` removes every truncation path; `usize → u64` is lossless on every supported target.
- **Implementation summary:** Widened `MappedRead::pos`, `MappedRead::mate_pos`, `ContigEntry::length`, `WindowKey::pos`, `HeadKey::pos`, `PerFileOrder`, `window_anchor`, `DuplicateReadAcrossFiles.pos`, `RecordSpec::pos`/`mate_pos`, `ContigSpec::length`, and `HeaderOverrides::length_overrides` to `u64`. Replaced the truncating casts in `extract_header`, `validate_fasta_agreement`, `head_key`, and `record_buf_to_mapped_read` with `as u64`. In test fixtures, `usize::try_from` guards the (only meaningful on 32-bit) lossy step.
- **Review suggestion used verbatim?:** No
- **Adaptation:** `OutOfOrderRead`'s `prev_pos`/`this_pos` are already `u64` and untouched here — they are restructured by M5. `encode_order_key`'s parameter widened to `u64` so the function still compiles before M5 deletes it.
- **Verification performed:** Targeted regression test `a15_position_above_u32_max_round_trips_through_mapped_read` exercises a position above `u32::MAX` end-to-end through `from_open_crams`. Existing 22 tests still green.
- **Files changed:** `src/per_sample_caller/cram_input.rs`, `src/per_sample_caller/errors.rs`, `src/per_sample_caller/record_specs.rs`, `src/per_sample_caller/cram_files.rs`
- **Tests added or modified:** `a15_position_above_u32_max_round_trips_through_mapped_read`
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 23 passed
  - `./scripts/dev.sh cargo test --lib` → 0, 84 passed
  - `cargo fmt --check` → 0, clean (after `cargo fmt`)
  - `cargo clippy --all-targets --all-features` → no new warnings in `per_sample_caller` (pre-existing warnings in out-of-scope `decompression_pool.rs` etc. unchanged)
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M2 — `MultipleSampleNames` reused for within-file SM disagreement
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Bug confirmed at the call site: when a single CRAM has two `@RG`s with different `SM` values, the variant `MultipleSampleNames` was emitted with `path_a == path_b`, producing a self-referential error message and conflating two distinct failure modes.
- **Implementation summary:** Added `CramInputError::MultipleSampleNamesInFile { path, rg_a, sm_a, rg_b, sm_b }` and rewrote `extract_single_sample_name` to track the *(rg_id, sm)* pair as it iterates and emit the new variant on first disagreement. Cross-file `MultipleSampleNames` is unchanged and still raised by the cross-file checks in `CramMergedReader::new`.
- **Review suggestion used verbatim?:** Yes (variant shape + emission point match the review).
- **Adaptation:** None.
- **Verification performed:** Added `b4_sample_tag_handling` sub-case (4) that builds a single CRAM with two disagreeing `@RG`s and asserts the new variant + the two `(rg_id, sm)` pairs match the expected order.
- **Files changed:** `src/per_sample_caller/cram_input.rs`, `src/per_sample_caller/errors.rs`
- **Tests added or modified:** `b4_sample_tag_handling` (new sub-case asserting `MultipleSampleNamesInFile`)
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 23 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None


## 5. Deferred findings to carry forward
(to be filled in)

## 6. Disputed findings to return to reviewer
(to be filled in)

## 7. Failed-validation findings
(to be filled in)

## 8. Blocked-by-context-mismatch findings
(to be filled in)

## 9. Commands run
(to be filled in)

## 10. Command results
(to be filled in)

## 11. Notes
- Baseline before fixes: `./scripts/dev.sh cargo test --lib per_sample_caller` → 22 passed.
- Each finding is implemented as a separate commit so the per-finding diff is reviewable in isolation.
