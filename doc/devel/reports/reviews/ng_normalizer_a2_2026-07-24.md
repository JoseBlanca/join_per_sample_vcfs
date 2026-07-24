# Code Review: ng normalizer — A2 (the leftmost property checker / oracle)

**Date:** 2026-07-24
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** working-tree — new file `src/ng/alignment/leftmost_property.rs` + `#[cfg(test)] mod` registration in `src/ng/alignment/mod.rs`
**Status:** Approve-with-changes

---

### 1. Scope
- New `#[cfg(test)]` module: the leftmost property checker that grades algorithms 1a/1b/1c against the *definition* of leftmost (plan A2; spec §6; arch §Test & bench shape).
- In-scope: [leftmost_property.rs](../../../../src/ng/alignment/leftmost_property.rs) (all), the 2-line mod.rs registration.
- Categories dispatched (8, parallel): reliability, errors, naming, idiomatic, refactor_safety, smells, extras, module_structure.

### 2. Verdict
Approve-with-changes.

### 3. Execution status
- `cargo fmt --check` → 0, clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean.
- `cargo test --lib` → 0, 2292 passed / 4 ignored (post-fix); `leftmost_property` 20/20.

### 4. Open questions and assumptions
1. Oracle *logic* correctness: two reviewers independently traced the shift condition against spec §6 and confirmed **no false positive** (on a genuinely-leftmost alignment) and **no false negative** (a fixed point of one-base slides is leftmost). The findings are all in the test suite guarding it, not the logic.

### 5. Top 3 priorities
1. **B1** — `reference_offset` handling is untested (a shift-invariant fixture lets the offset be dropped with the whole suite still green). Applied.
2. **M-suite** — soft-clip cursor, multi-base insertion far-edge, indel-arm `left_neighbour` reset, and the match-required clause each have a surviving mutant. Applied.
3. **Mi (errors)** — malformed alignment panics opaquely; a named well-formedness `debug_assert` improves the oracle's own diagnostics. Applied.

### 6. Findings

**Blocker**
- **B1: leftmost_property.rs — `the_reference_offset_is_respected` does not pin `reference_offset`.** (reliability, High) Its `TTGAAAAT` homopolymer fixture is shift-invariant, so dropping the offset passes both assertions and the entire 13-test suite — the oracle would silently mis-grade every mid-stretch placement. **Applied**: replaced with `a_shiftable_deletion_is_found_using_the_reference_offset_not_the_stretch_start` (prefix `CCCC`, so dropping the offset flips the verdict).

**Major**
- **M1: soft-clip read-cursor advance untested** (reliability). Same shift-invariance blind spot. **Applied**: `a_soft_clip_advances_the_read_cursor_so_the_crossed_read_base_is_correct` (`GCAAT` fixture where the wrong cursor lands on `C`).
- **M2: multi-base insertion far-edge untested** (reliability). Only `Insertion(1)` non-leading cases existed; the `+ n - 1` collapses at `n=1`. **Applied**: `a_rightmost_period_two_insertion_is_detected_as_shiftable`.
- **M3: indel-arm `left_neighbour_is_aligned = false` resets untested** (reliability). Dropping the deletion-arm reset makes a **false-positive** oracle that fails a correct normalizer. **Applied**: `an_indel_abutting_a_preceding_indel_is_not_treated_as_slidable` (deletion + insertion sub-cases).
- **M4: the "crossed column must match" clause not independently pinned** (extras). The lone mismatch fixture (`CAA`/`CT`) has both clauses false. **Applied**: `a_deletion_whose_far_edge_matches_but_left_column_mismatches_is_leftmost` (`TCA`/`TA`).

**Minor**
- **Mi1: malformed alignment panics opaquely** (errors). **Applied**: named well-formedness `debug_assert` (`consumed_read`/`consumed_reference`) + documented precondition.
- **Mi2: `kind: &'static str` stringly-typed** (idiomatic). **Applied**: `IndelKind` enum + `label()`; the `pub(crate)` witness API the future 1a/1b/1c test modules match on is now exhaustive-checkable.

**Nits** (won't-fix, recorded)
- refactor_safety: `n as usize` not width-change-protected → matches ng's established `as usize` convention; leave.
- smells: two near-identical indel arms (below the 3-occurrence threshold; the meaningful difference is the far-edge source) → leave.
- naming: `fn alignment` test factory is a noun → established test-factory idiom; leave.
- Also applied opportunistically from the challenge pass: empty-cigar, trailing-insertion-at-read-end (far-edge boundary), `Skip`-arm, and chained-deletion reference-advance tests.

### 7. Out of scope observations
None.

### 8. Missing tests to add now
All applied — see Findings. Test count in `leftmost_property`: 13 → 20.

### 9. What's good
- The oracle is derived from the definition, not from `normalize_alleles` — so it is genuinely independent of algorithm 1a, which the module doc states and the plan requires.
- The exhaustive wildcard-free `match` over `CigarOp` means a 10th op is a compile error, not a silent fall-through.

### 10. Commands to re-verify
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::alignment::leftmost_property`
