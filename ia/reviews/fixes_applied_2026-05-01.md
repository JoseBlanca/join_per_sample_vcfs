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
| M3 | Major | `decode_md5_hex` silently tolerates malformed `M5` | Apply | Applied | No | `cram_input.rs`, `errors.rs` | Pass | No |
| M4 | Major | `peek_head_keys` reports `MalformedRecord` for legitimate kept-unmapped reads | Apply | Applied | No | `cram_input.rs` | Pass | No |
| M5 | Major | `OutOfOrderRead` carries packed `(ref_id, pos)` keys | Apply | Applied | No | `cram_input.rs`, `errors.rs` | Pass | No |
| Mi1 | Minor | empty-input check returns `CramInputError::Io` | Apply | Applied | No | `cram_input.rs`, `errors.rs` | Pass | No |
| Mi2 | Minor | defensive `.unwrap_or_default()` on `canonical_path` masks an invariant | Apply | Applied | No | `cram_input.rs` | Pass | No |
| Mi3 | Minor | iterator behaviour after `Some(Err)` is unspecified | Apply | Applied | No | `cram_input.rs` | Pass | No |
| Mi4 | Minor | `MalformedRecord` errors with empty `qname` context | Apply | Applied | No | `cram_input.rs`, `errors.rs` | Pass | No |
| Mi5 | Minor | `min_read_length` filter pays for full record decode before rejecting | Apply | Applied | No | `cram_input.rs` | Pass | No |
| Mi6 | Minor | MAPQ-missing collapses to 0 — needs doc and regression test | Apply | Applied | No | `cram_input.rs` | Pass | No |
| N1 | Nit | `OpenCram::path_for_errors` rename | Apply | TBD | No | TBD | TBD | TBD |
| N2 | Nit | `window: VecDeque` → `Vec` | Apply | TBD | No | TBD | TBD | TBD |
| N3 | Nit | WindowKey clone | Apply | TBD | No | TBD | TBD | TBD |
| N4 | Nit | `encode_order_key` removal | Apply (subsumed by M5) | Superseded | No | `cram_input.rs` | (covered by M5) | No |
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

### Mi6 — MAPQ-missing collapses to 0 — needs doc and regression test
- **Severity:** Minor
- **Initial decision:** Apply (decision: missing MAPQ → MAPQ 0 → reject under any non-zero minimum)
- **Final status:** Applied
- **Reasoning:** The behaviour is correct — `unwrap_or(0)` matches the project's coverage-target rationale — but it was undocumented and not pinned. A future refactor swapping `unwrap_or(0)` for `unwrap_or(u8::MAX)` would silently let unknown-MAPQ reads through.
- **Implementation summary:** Extended the doc comment on `DEFAULT_MIN_MAPQ` with the missing-MAPQ rule and rationale (low coverage targets, bcftools/freebayes convention). Added `a7b_missing_mapq_is_treated_as_zero_and_filtered_by_default` exercising the `mapping_quality()==None` path: builds a record with `spec.mapq=0` (the `record_spec` guard maps this to no `MappingQuality` set) and asserts it lands in the `low_mapq` bucket under default config.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** New test exercises the `None` mapping-quality path; existing `a7_min_mapq_filter` continues to pass.
- **Files changed:** `src/per_sample_caller/cram_input.rs`
- **Tests added or modified:** `a7b_missing_mapq_is_treated_as_zero_and_filtered_by_default`
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 29 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### Mi5 — `min_read_length` filter pays for full record decode before rejecting
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `record_buf_to_mapped_read` clones qname, walks the CIGAR vector, uppercases SEQ, and copies QUAL — all work that is thrown away when the post-decode min-read-length check rejects the record. SEQ length is available on the `RecordBuf` directly, so the check can run earlier.
- **Implementation summary:** Moved the `min_read_length` check from after `record_buf_to_mapped_read` to immediately after the order/dup checks consume the head. Reads `record.sequence().as_ref().len()` directly. Bucket increment and `continue` semantics unchanged.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** `a11_min_read_length_drops_short_and_empty` continues to pass — same input set (long, short, empty), same expected counts (`too_short == 2`).
- **Files changed:** `src/per_sample_caller/cram_input.rs`
- **Tests added or modified:** None — `a11` already exercises the filter and the bucket count.
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 28 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### Mi4 — `MalformedRecord` errors with empty `qname` context
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** When the qname is unknown (peek failed before the record could be inspected), the variant emitted `qname=''` which doesn't help localise the failure.
- **Implementation summary:** Changed `MalformedRecord.qname` from `String` to `Option<String>`. Updated the `#[error]` format string to omit the `(qname='…')` clause when `None` via a small `match` expression in the format args. Updated all five call sites: the two consumed-on-`fail` paths in `Iterator::next` use `Some(qname)` from the head; the three call sites where the record cannot be inspected (`refill_heads` peek error, `peek_head_keys` peek error, and the defensive `head_key()==None` arm) pass `None` (or `Some(name)` when the record name is recoverable).
- **Review suggestion used verbatim?:** Yes — picked the `Option<String>` shape.
- **Adaptation:** None; format string interpolation uses thiserror's expression-form to render the qname clause conditionally.
- **Verification performed:** Added `p2_malformed_record_message_omits_empty_qname`, asserting the rendered message contains `qname='R1'` when set and contains no `qname` substring when the field is `None`.
- **Files changed:** `src/per_sample_caller/cram_input.rs`, `src/per_sample_caller/errors.rs`
- **Tests added or modified:** `p2_malformed_record_message_omits_empty_qname`
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 28 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### Mi3 — iterator behaviour after `Some(Err)` is unspecified
- **Severity:** Minor
- **Initial decision:** Apply (decision: fuse on first error)
- **Final status:** Applied
- **Reasoning:** When `OutOfOrderRead` or `DuplicateReadAcrossFiles` was returned the offending head was not consumed; a caller using `for r in reader { … }` without breaking on error would loop the same error.
- **Implementation summary:** Added a `fused: bool` field; the constructors initialise it to `false`. Added a private `fail(self, e)` helper that flips `fused` to `true` and returns `Some(Err(e))`. Every `return Some(Err(...))` site in `Iterator::next` now returns `self.fail(...)`. `next()` early-returns `None` once `fused`. Implemented `std::iter::FusedIterator` for combinator correctness. Added a doc comment on the `impl Iterator` block.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** Added `a4b_iterator_returns_none_after_first_error`, mirroring `a4`'s fixture: drains the first ok record, then asserts the error variant, then asserts every subsequent `next()` returns `None`.
- **Files changed:** `src/per_sample_caller/cram_input.rs`
- **Tests added or modified:** `a4b_iterator_returns_none_after_first_error`
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 27 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### Mi2 — defensive `.unwrap_or_default()` on `canonical_path` masks an invariant
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `canonical_path` is `Some` whenever `canonical_contigs` or `canonical_sample` is `Some`, by construction. The two error paths previously substituted `PathBuf::new()` if the invariant ever broke, hiding a real bug. The existing call sites at the end of the function already use `.expect`.
- **Implementation summary:** Replaced both `.unwrap_or_default()` calls with `.expect("canonical_path is set when canonical_contigs is Some")` / `.expect("canonical_path is set when canonical_sample is Some")`.
- **Review suggestion used verbatim?:** Yes (chose the minimal `.expect` route, not the larger refactor).
- **Adaptation:** None.
- **Verification performed:** All 26 tests still pass; the `.expect` strings document the invariant in code.
- **Files changed:** `src/per_sample_caller/cram_input.rs`
- **Tests added or modified:** None — the change is a defensive assertion on an existing invariant that none of the current tests exercise.
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 26 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### Mi1 — empty-input check returns `CramInputError::Io`
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** A programming error (empty inputs) was being squeezed through the `Io` variant with a synthesized `io::Error`, which (a) emitted "I/O error on '': …" and (b) was indistinguishable from a real I/O failure for callers that match on the variant.
- **Implementation summary:** Added `CramInputError::NoInputs`. The empty-input guard in `CramMergedReader::new` now returns it directly.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** Added `a16_empty_inputs_returns_no_inputs`, which calls `new(&[], …)` and asserts the variant via a `match`.
- **Files changed:** `src/per_sample_caller/cram_input.rs`, `src/per_sample_caller/errors.rs`
- **Tests added or modified:** `a16_empty_inputs_returns_no_inputs`
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 26 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M5 — `OutOfOrderRead` carries packed `(ref_id, pos)` keys
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Bug confirmed: `OutOfOrderRead.prev_pos` and `this_pos` carried `(ref_id << 32) | pos` packed values, producing user-facing positions like `4294967396` and silently colliding when the cross-chromosome offset was `< 2^32`.
- **Implementation summary:** Replaced packed fields with structured `prev_ref_id`, `prev_pos`, `this_ref_id`, `this_pos`. The `#[error]` format string now renders `(ref_id=…, pos=…)` for both prev and this. The packing helper `encode_order_key` is gone (this also resolves Nit "encode_order_key removal" — N4).
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** Updated `a4_out_of_order_within_a_single_stream` to assert on each of the four structured fields. The previous assertion only checked the qname.
- **Files changed:** `src/per_sample_caller/cram_input.rs`, `src/per_sample_caller/errors.rs`
- **Tests added or modified:** `a4_out_of_order_within_a_single_stream` (now asserts the four structured fields)
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 25 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### N4 — `encode_order_key` removal
- **Severity:** Nit
- **Initial decision:** Apply (subsumed by M5)
- **Final status:** Superseded (by M5)
- **Reasoning:** `encode_order_key` was the helper that built the packed-u64 OutOfOrderRead key; with M5 restructuring the variant, the function has no remaining callers.
- **Implementation summary:** Removed alongside the M5 fix.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** Tests still green; no remaining references in tree.
- **Files changed:** `src/per_sample_caller/cram_input.rs`
- **Tests added or modified:** None.
- **Validation:** Covered by M5's run.
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M4 — `peek_head_keys` reports `MalformedRecord` for legitimate kept-unmapped reads
- **Severity:** Major
- **Initial decision:** Apply (decision: drop unmapped unconditionally)
- **Final status:** Applied
- **Reasoning:** With `drop_unmapped=false`, an unmapped record (`reference_sequence_id` = `None`, `alignment_start` = `None`) survives `refill_heads`, then trips `peek_head_keys` and is reported as `MalformedRecord`. Removing the toggle eliminates the path: unmapped reads are always dropped at `classify_pre_decode` and never reach `peek_head_keys`.
- **Implementation summary:** Removed `CramMergedReaderConfig::drop_unmapped`. `classify_pre_decode` now drops `FLAG_UNMAPPED` unconditionally at the same hit-rate-ordered position. `peek_head_keys`'s defensive `head_key()==None` arm is kept and re-commented as defence in depth against malformed CRAM bytes (input that claims a read is mapped but omits both fields).
- **Review suggestion used verbatim?:** Yes for the config and `classify_pre_decode` change; chose review's option (a) for the defensive arm in `peek_head_keys`.
- **Adaptation:** None.
- **Verification performed:** Updated `a9_each_flag_drop_one_at_a_time` to drop the `unmapped` row from the cases table (no toggle to flip). Renamed and rewrote `a10_all_optional_flag_drops_disabled_keeps_everything_except_unmapped` so its expectation matches the new contract (5 records, `unmapped == 1`). Added `a14b_truly_unmapped_record_is_filtered_not_errored` to pin the new behaviour.
- **Files changed:** `src/per_sample_caller/cram_input.rs`
- **Tests added or modified:** `a9_each_flag_drop_one_at_a_time` (cases table updated), `a10_all_optional_flag_drops_disabled_keeps_everything_except_unmapped` (renamed + new expectation), `a14b_truly_unmapped_record_is_filtered_not_errored` (new)
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 25 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M3 — `decode_md5_hex` silently tolerates malformed `M5`
- **Severity:** Major
- **Initial decision:** Apply (decision: hard error on malformed M5)
- **Final status:** Applied
- **Reasoning:** Bug confirmed: `md5_from_reference_sequence` collapsed every malformed `M5` to `None`, which the wildcard `ContigEntry::eq` rule then accepted as "absent". A subtly wrong CRAM (truncated, non-hex character, lowercase typo) was silently treated as MD5-clean, contradicting the "errors must not pass silently" design principle.
- **Implementation summary:** Promoted `md5_from_reference_sequence` to fallible (`Result<Option<[u8;16]>, CramInputError>`). Threaded the contig name and CRAM path so the new error variant `CramInputError::MalformedMd5 { path, contig, detail }` can be emitted with usable diagnostics. `extract_header` now propagates the error via `?`. A present-but-malformed `M5` is therefore a hard error; an absent `M5` still returns `Ok(None)` and acts as the documented wildcard.
- **Review suggestion used verbatim?:** No
- **Adaptation:** Detail string distinguishes "wrong length" from "non-hex character" rather than always reporting the byte count, since both inputs are real-world failure modes; the test exercises both branches.
- **Verification performed:** Added `b4b_malformed_md5_rejected` covering three sub-cases: a non-hex 7-byte M5, a 30-char hex M5, and a 32-char string with a non-hex character — all three produce `MalformedMd5`. Existing CRAM-with-good-MD5 tests still pass.
- **Files changed:** `src/per_sample_caller/cram_input.rs`, `src/per_sample_caller/errors.rs`
- **Tests added or modified:** `b4b_malformed_md5_rejected`
- **Validation:**
  - `./scripts/dev.sh cargo test --lib per_sample_caller` → 0, 24 passed
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
