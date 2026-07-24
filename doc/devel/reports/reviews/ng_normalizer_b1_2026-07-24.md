# Code Review: ng normalizer — B1 (algorithm 1a, `StructuredLeftAligner`)

**Date:** 2026-07-24
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** working-tree — new file `src/ng/alignment/left_align_structured.rs` + `pub mod` registration in `src/ng/alignment/mod.rs`
**Status:** Approve-with-changes

---

### 1. Scope
- Algorithm 1a: `StructuredLeftAligner`, an `AlignmentNormalizer` wrapping production's `left_align_indels` (plan B1; spec §6; arch §5).
- In-scope: [left_align_structured.rs](../../../../src/ng/alignment/left_align_structured.rs), the `pub mod` line.
- Categories dispatched (8, parallel): reliability, errors, naming, idiomatic, refactor_safety, smells, extras, module_structure.

### 2. Verdict
Approve-with-changes.

### 3. Execution status
- `cargo fmt --check` → 0, clean · `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean · `cargo test --lib` → 0, 2305 passed / 4 ignored (post-fix); `left_align_structured` 13/13.

### 4. Open questions and assumptions
1. **Both anchor oracles confirmed genuine, not vacuous** (reliability): byte-parity recomputes `expected` from the raw `(offset, read, reference)` inputs (so a transposition / dropped offset / wrong slice diverges), and the property test's rightmost-placement cases would fail a do-nothing port. Intent (wrapper not naive pass; both behaviours; both verifications) confirmed (extras).
2. **Name `StructuredLeftAligner`**: both name-reviewers say keep — "left-aligner" is standard genomics vocabulary and the whole normalizer sub-family shares it.

### 5. Top priorities
1. Pin the `reference_offset` boundary/panic behaviour and add empty/idempotence/leading-deletion coverage (reliability M1/M2 + Minors). Applied.
2. Assert the documented `reference_offset` precondition + reword the doc (errors/extras/naming, convergent). Applied.

### 6. Findings

**Major**
- **M1: `reference_offset` vs `reference.len()` boundary/panic untested** (reliability). The panic-on-`offset > len` is the intended fail-fast, but nothing pins it (a future clamp would silently reintroduce a wrong variant), and the legal `offset == len` empty-slice boundary is unexercised. **Applied**: `an_offset_past_the_reference_panics` (`#[should_panic]`, real in release via slice bounds) + `an_offset_equal_to_the_reference_length_does_not_panic`.
- **M2: no empty-read/empty-cigar coverage for `StructuredLeftAligner`** (reliability). **Applied**: `an_empty_cigar_is_left_untouched`.

**Minor**
- **Mi1: precondition documented but unasserted** (errors + extras, convergent). **Applied**: `debug_assert!(offset <= reference.len(), …)` matching the module's named-precondition pattern.
- **Mi2: doc bundles two preconditions, attributing production's graceful guard to both** (naming). **Applied**: reworded to separate the graceful cigar/read guard from the panicking `reference_offset` bound.
- **Mi3: already-leftmost checked for leftmost-ness but not byte-identity** (reliability). **Applied**: `an_already_leftmost_deletion_comes_back_byte_identical` (idempotence).
- **Mi4: no leading-deletion-at-nonzero-offset case** (reliability). **Applied**: `a_deletion_rolling_to_the_read_start_keeps_a_nonzero_reference_offset`.
- **Mi5: cross-block propagation asserted only indirectly** (reliability). **Applied**: `two_deletions_across_a_block_merge_and_left_align` + the same fixture added to the byte-parity list.
- **Mi6: behavioural back-reference into `pileup::walker`** (module_structure). **Applied**: recorded as explicit debt in a comment on the import, alongside the existing `CigarOp` note; acceptable (reuse-over-rewrite, production frozen).

**Nits** (recorded)
- refactor_safety: `reference_offset as usize` lossy cast → **won't fix**, consistent with A2's disposition (ng one-width `as usize` convention, lossless on the 64-bit target; the new `debug_assert` bounds it). smells: duplicate `Case` struct across the two table tests → **applied** (hoisted to one shared `Case`).

### 7. Out of scope observations
None.

### 8. Missing tests to add now
All applied (6 new tests); `left_align_structured` 7 → 13.

### 9. What's good
- The wrapper re-derives nothing: byte-parity against `left_align_indels` proves the two load-bearing behaviours (merge + cross-block propagation) are inherited, exactly the plan's "reuse over rewrite" intent.
- The property oracle grading the port passes on every fixture — the A2 oracle and production's left-aligner agree, no divergence to record.

### 10. Commands to re-verify
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::alignment::left_align_structured`
