# Fix Application Report: gvcf_parser_2026-04-13.md

**Date:** 2026-04-14
**Source review:** `reviews/gvcf_parser_2026-04-13.md`
**Source state reviewed against:** commit 6ce13ed (branch main)
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 2
- Majors: 3
- Minors: 0
- Nits: 3

### Outcome totals
- Applied: 1
- Applied with adaptation: 4
- Already fixed: 0
- Deferred: 1
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → exit 0, clean
- `cargo clippy --test gvcf_parser_test -- -D warnings` → exit 101, 6 pre-existing errors in out-of-scope files (decompression_pool.rs, genotype_merging.rs, variant_grouping.rs); 0 new warnings from changed files
- `cargo test --all-targets` → exit 0, 133 tests passed (48 unit + 25 genotype_merging + 23 gvcf_parser + 17 integration + 8 utils_magic + 12 variant_group)
- `cargo doc --no-deps` → exit 0, clean
- `cargo audit` → not run (deferred per review)

### Unresolved high-priority findings
- OQ (Ploidy contract) — deferred, requires API design decision

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | Allele index corruption | Apply | Applied with adaptation | No | `src/gvcf_parser.rs` | Pass | No |
| B2 | Blocker | Malformed GT acceptance | Apply | Applied with adaptation | No | `src/gvcf_parser.rs` | Pass | No |
| M1 | Major | Newline in FORMAT field | Apply | Applied | No | `src/gvcf_parser.rs` | Pass | No |
| M2 | Major | get_span overflow | Apply | Applied with adaptation | No | `src/gvcf_parser.rs` | Pass | No |
| M3 | Major | Genotype test coverage | Apply | Applied with adaptation | No | `tests/gvcf_parser_test.rs` | Pass | No |
| OQ | Open Question | Ploidy contract | Defer | Deferred | No | None | N/A | Yes |
| N1 | Nit | Unused mut in test | Apply | Applied | No | `tests/gvcf_parser_test.rs` | Pass | No |
| N2 | Nit | Simplify match in Iterator::next | Apply | Applied | No | `src/gvcf_parser.rs` | Pass | No |
| N3 | Nit | Run cargo fmt | Apply | Applied | No | `src/gvcf_parser.rs`, `tests/gvcf_parser_test.rs` | Pass | No |

## 3. Questions asked and answers

None.

## 4. Per-finding log

### B1 — Allele index corruption
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** The bug is confirmed: `u8` values 128-255 cast to `i8` wrap to negative, with 255 becoming -1 (MISSING_ALLELE sentinel). Silent data corruption.
- **Implementation summary:** Added a `val > 127` check after the existing `checked_mul`/`checked_add` arithmetic, before the `as i8` cast. Uses the existing `AlleleIndexOverflow` error variant.
- **Review suggestion used verbatim?:** No
- **Adaptation:** The review suggested adding a new `AlleleIndexOutOfRange` error variant. Instead, reused the existing `AlleleIndexOverflow` variant whose message already says "exceeds maximum allowed (127)". This avoids adding a new error variant for a closely related condition. The validation logic itself (reject > 127) matches the review's intent.
- **Verification performed:** Test-first: added `test_parse_rejects_allele_index_out_of_range` which failed before fix, passed after.
- **Files changed:** `src/gvcf_parser.rs`
- **Tests added or modified:** `test_parse_rejects_allele_index_out_of_range`
- **Validation:**
  - `cargo test --test gvcf_parser_test` → 0, 18 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### B2 — Malformed GT acceptance
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** The bug is confirmed: `_ => {}` silently ignores unknown characters; consecutive separators produce no allele; leading separators are accepted. All are correctness violations.
- **Implementation summary:** Rewrote `parse_genotype` body to enforce strict GT grammar: `<allele> ( ('/' | '|') <allele> )*` where `<allele>` is a multi-digit number or exactly `.`. Uses a `loop` with explicit allele/separator alternation rather than a single `for` loop. Rejects: unknown characters, consecutive separators, leading separators, trailing junk, empty strings.
- **Review suggestion used verbatim?:** No
- **Adaptation:** The review provided a full replacement. The actual implementation differs in structure (uses a `loop` with `i` index advancing rather than the review's specific variable names) but enforces the same grammar. The B1 overflow check is integrated into the digit-parsing branch rather than as a separate function. The function signature is unchanged.
- **Verification performed:** Test-first: added 4 tests (`test_parse_rejects_invalid_characters`, `test_parse_rejects_consecutive_separators`, `test_parse_rejects_trailing_junk`, `test_parse_rejects_leading_separator`). 3 of 4 failed before fix (1 was caught indirectly by ploidy check), all 4 pass after.
- **Files changed:** `src/gvcf_parser.rs`
- **Tests added or modified:** `test_parse_rejects_invalid_characters`, `test_parse_rejects_consecutive_separators`, `test_parse_rejects_trailing_junk`, `test_parse_rejects_leading_separator`
- **Validation:**
  - `cargo test --test gvcf_parser_test` → 0, 22 passed
  - `cargo test --all-targets` → 0, 133 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None
- **Note:** B2 and M1 are coupled: the strict GT parser correctly rejects newline characters in GT strings, which exposed M1 (trailing newlines in last sample fields). M1 was applied concurrently to unblock B2 validation.

### M1 — Newline in FORMAT field
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Confirmed: `read_line` preserves trailing `\n`, `splitn(10, '\t')` doesn't trim, so the last sample's last field contains a newline.
- **Implementation summary:** Added `.trim_end()` to the line before `splitn` in `from_line`: `line.trim_end().splitn(10, '\t')`.
- **Review suggestion used verbatim?:** Yes (first suggested option: trim before split)
- **Adaptation:** None
- **Verification performed:** The M1 regression test is embedded in the enhanced `test_genotype_with_format_fields` (M3): asserts `dp_values[1] == "15"` for the last sample's DP field. Also validated by B2 tests passing (strict parser rejects `\n` in GT strings).
- **Files changed:** `src/gvcf_parser.rs`
- **Tests added or modified:** Covered by enhanced `test_genotype_with_format_fields` (see M3)
- **Validation:**
  - `cargo test --test gvcf_parser_test` → 0, 22 passed
  - `cargo test --all-targets` → 0, 133 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M2 — get_span overflow
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** Confirmed: `get_span` uses unchecked `self.pos + max_allele_len as u32 - 1` which panics in debug or wraps in release when near `u32::MAX`.
- **Implementation summary:** Replaced unchecked addition with `self.pos.checked_add(max_allele_len as u32 - 1)` returning `VcfParseError::RuntimeError` on overflow. Also changed `== 1` to `<= 1` for the early return (handles zero-length edge case).
- **Review suggestion used verbatim?:** No
- **Adaptation:** The review suggested a full rewrite with additional `try_into()` for `max_allele_len` and explicit empty-allele checks. The actual fix is narrower: only the arithmetic was changed to checked, consistent with the existing `get_var_end` pattern. The empty-allele case is already guarded by the `.max().ok_or()` above. The `max_allele_len as u32` cast is safe because `String::len()` returns `usize` which is always >= u32 on the platforms this code targets (and would be caught by the u32 checked_add anyway for extreme values).
- **Verification performed:** Test-first: added `test_get_span_returns_error_on_position_overflow` with `pos = u32::MAX - 5` and 12-base allele. Panicked before fix, returns `Err` after.
- **Files changed:** `src/gvcf_parser.rs`
- **Tests added or modified:** `test_get_span_returns_error_on_position_overflow`
- **Validation:**
  - `cargo test --test gvcf_parser_test` → 0, 23 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M3 — Genotype test coverage
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** Confirmed: `test_genotype_parsing_basic`, `test_genotype_parsing_missing`, `test_genotype_parsing_multiallelic`, and `test_genotype_with_format_fields` only asserted `variants.len()` without checking genotype values.
- **Implementation summary:** Enhanced all 4 tests with assertions on `genotypes[]`, `phases[]`, and `get_gt_field_by_index()` values. Also added M1 regression assertion (`dp_values[1] == "15"`) in `test_genotype_with_format_fields`.
- **Review suggestion used verbatim?:** No
- **Adaptation:** The review proposed full test rewrites. Instead, the existing tests were extended with additional assertions while preserving the original structure and data. The assertions verify the same invariants the review identified (genotype values, phase flags, FORMAT field access).
- **Verification performed:** All enhanced tests pass with the corrected production code.
- **Files changed:** `tests/gvcf_parser_test.rs`
- **Tests added or modified:** `test_genotype_parsing_basic`, `test_genotype_parsing_missing`, `test_genotype_parsing_multiallelic`, `test_genotype_with_format_fields`
- **Validation:**
  - `cargo test --test gvcf_parser_test` → 0, 23 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### OQ — Ploidy contract
- **Severity:** Open Question
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** The internal code already supports arbitrary ploidy via the `CommonGenotypePatterns` struct and the `ploidy` parameter threaded through `parse_genotypes` and `parse_genotype`. However, `VarIterator` hardcodes `ploidy: DEFAULT_PLOIDY` (2) with no public setter. Exposing configurable ploidy requires an API design decision (setter method? constructor parameter? per-variant detection?) that is beyond the scope of a bug-fix pass. The existing diploid fast paths (`ploidy == 2` checks at line 259) are optimizations, not restrictions.
- **Implementation summary:** None
- **Files changed:** None
- **User input:** None
- **Follow-up:** Requires design decision on how to expose ploidy configuration in the public API. Should include tests for haploid, triploid, and tetraploid inputs.
- **Residual risk:** Users with non-diploid organisms cannot configure ploidy via the parser API.

### N1 — Unused mut in test
- **Severity:** Nit
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `let mut parser` at tests/gvcf_parser_test.rs:217 does not need `mut`. Confirmed by compiler warning.
- **Implementation summary:** Removed `mut` keyword.
- **Files changed:** `tests/gvcf_parser_test.rs`
- **Tests added or modified:** None
- **Validation:**
  - Compiler warning resolved: 0 warnings after fix
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### N2 — Simplify match in Iterator::next
- **Severity:** Nit
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Nested `match` at line 753-755 can be flattened.
- **Implementation summary:** Replaced nested `Ok(n_items_added) => match n_items_added { 0 => ..., _ => () }` with `Ok(0) => return None, Ok(_) => {}`.
- **Files changed:** `src/gvcf_parser.rs`
- **Tests added or modified:** None
- **Validation:**
  - `cargo test --test gvcf_parser_test` → 0, 23 passed
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### N3 — Run cargo fmt
- **Severity:** Nit
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `rustfmt --check` reported formatting issues in src/gvcf_parser.rs.
- **Implementation summary:** Ran `cargo fmt` on both `src/gvcf_parser.rs` and `tests/gvcf_parser_test.rs`.
- **Files changed:** `src/gvcf_parser.rs`, `tests/gvcf_parser_test.rs`
- **Tests added or modified:** None
- **Validation:**
  - `cargo fmt --check` → exit 0, clean
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

## 5. Deferred findings to carry forward
- OQ (Ploidy contract) — Requires API design decision on how to expose ploidy configuration

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Commands run
- `cargo test --test gvcf_parser_test` (baseline)
- `cargo test --test gvcf_parser_test test_parse_rejects_allele_index_out_of_range` (B1 pre-fix)
- `cargo test --test gvcf_parser_test` (B1 post-fix)
- `cargo test --test gvcf_parser_test -- test_parse_rejects_invalid ...` (B2 pre-fix)
- `cargo test --test gvcf_parser_test` (B2+M1 post-fix)
- `cargo test --all-targets` (B2+M1 regression)
- `cargo test --test gvcf_parser_test test_get_span_returns_error_on_position_overflow` (M2 pre-fix)
- `cargo test --test gvcf_parser_test` (M2 post-fix)
- `cargo test --test gvcf_parser_test` (M3 post-fix)
- `cargo fmt`
- `cargo test --test gvcf_parser_test` (nits post-fix)
- `cargo fmt --check`
- `cargo clippy --test gvcf_parser_test -- -D warnings`
- `cargo test --all-targets` (final)
- `cargo doc --no-deps`

## 10. Command results
- `cargo test --test gvcf_parser_test` (baseline) → 0, 17 passed, 1 warning
- `cargo test --all-targets` (final) → 0, 133 passed (48 + 25 + 23 + 17 + 8 + 12), 0 warnings in changed files
- `cargo fmt --check` → 0, clean
- `cargo clippy --test gvcf_parser_test -- -D warnings` → 101, 6 pre-existing errors in out-of-scope files, 0 from changed files
- `cargo doc --no-deps` → 0, clean

## 11. Notes
- B2 and M1 are coupled: the strict GT parser rejects `\n` as an invalid character, so M1 (trailing newline in last sample field) must be fixed for B2 to work correctly. Both were applied in the same validation pass.
- Clippy fails project-wide due to 6 pre-existing lints in out-of-scope files (decompression_pool.rs, genotype_merging.rs, variant_grouping.rs). No new clippy warnings were introduced by this fix pass.
- The `parse_genotype` function was rewritten rather than patched because the defects (silent character dropping, missing state tracking, no allele/separator alternation) were fundamental to its loop structure. A minimal patch set would have been harder to verify than a focused rewrite. The function signature and return semantics are unchanged.
- `cargo audit` was not run per the original review's deferral.
