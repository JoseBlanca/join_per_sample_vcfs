# Code Review: ng normalizer — C1 (algorithm 1b, `RepeatedLeftAligner`)

**Date:** 2026-07-24
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** working-tree — new file `src/ng/alignment/left_align_repeated.rs` + `pub mod` registration in `src/ng/alignment/mod.rs`
**Status:** Approve-with-changes

---

### 1. Scope
- Algorithm 1b: `RepeatedLeftAligner` — freebayes' repeated capped simple passes, reporting exhaustion (plan C1; spec §6/§8; arch §5). An independent left-aligner (not wrapping 1a/production).
- In-scope: [left_align_repeated.rs](../../../../src/ng/alignment/left_align_repeated.rs) (~580 lines), the `pub mod` line.
- Categories dispatched (8, parallel): reliability, errors, naming, idiomatic, refactor_safety, smells, extras, module_structure.

### 2. Verdict
Approve-with-changes.

### 3. Execution status
- `cargo fmt --check` → 0, clean · `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean · `cargo test --lib` → 0, 2322 passed / 4 ignored (post-fix); `left_align_repeated` 12 → 17.
- The reliability reviewer ran a **50,000-case differential fuzz** (1b vs 1a vs the property oracle), temporarily appended then reverted.

### 4. Open questions and assumptions
1. **No corruption or convergence bug** — the fuzz confirmed 1b always converges, always preserves read consumption, always produces leftmost output; the CIGAR surgery and the `1..=max_passes` off-by-one are sound; the `unreachable!` is genuinely unreachable.
2. **Independence confirmed** — the non-test core calls none of `left_align_indels`/`left_align_cigar`/`normalize_alleles`/`StructuredLeftAligner`. 1b is a rival, not a wrapper.
3. **Intent confirmed** — from-description not transliteration; exhaustion reported not swallowed; oracle-verified.

### 5. Top priorities
1. `#[must_use]` on `ConvergenceReport` — enforce the "cannot be dropped" thesis at compile time. Applied.
2. `shift_pass` must canonicalize on the no-shift path (Finding 1). Applied.
3. Correct the "never in spelling" overclaim re: complex-indel trimming (Finding 2). Applied (documented) + flagged at Checkpoint C.

### 6. Findings

**Major**
- **M1: `ConvergenceReport`/`left_align` not `#[must_use]`** (errors + idiomatic, convergent). The module's headline guarantee ("a caller cannot drop it by accident") was documented but not enforced — a bare `left_align(…);` statement silently discards `ExhaustedCap`. **Applied**: `#[must_use = "…"]` on the enum; the trait's `let _ =` remains the one sanctioned discard.
- **M2: `shift_pass` skips `canonicalize` on the no-shift path** (reliability, fuzz-found). An already-leftmost but non-canonically-spelled input (split same-kind gaps, I-before-D, zero-length ops) was returned verbatim while 1a always canonicalizes — scattering cohort buckets. **Applied**: `shift_pass` now always canonicalizes and reports movement by comparing to the input; verified by `an_already_leftmost_but_non_canonical_input_is_consolidated` and a new non-canonical agreement case.
- **M3: no deletion/insertion parsimony trimming** (reliability, fuzz-found). 1b can emit `D2 I1` where 1a trims to `D1`; both leftmost (the oracle grades placement, not parsimony), but different spelling. The `canonicalize` doc's "differ only in placement, never in spelling" was false. **Applied as a documented difference, not a code change**: corrected the doc, added a module-doc note, and pinned the behaviour (`a_complex_indel_overlap_is_not_trimmed_and_may_differ_from_1a_in_spelling`). Full trimming reimplements `normalize_alleles`' trim+shift interleaving — substantial, risky, beyond C1's "shift+cap+report" scope. **Flagged for owner decision at Checkpoint C** (whether 1b should trim for full 1a parity, given the D1 screen).

**Minor**
- **Mi1: `is_indel` used `matches!` with implicit `_`** (refactor_safety) — the `unreachable!` depends on it. **Applied**: exhaustive match.
- **Mi2: `with_max_passes(0)` silently no-ops** (reliability). **Applied**: doc note + `a_zero_cap_does_nothing_and_reports_exhaustion` test.
- **Mi3: helper duplication under-documented** (smells) — `op_len`/`with_len`/`is_indel`/`canonicalize` mirror private `indel_norm` twins. **Applied**: a shared debt note (the duplication itself is justified — production is frozen, a rival can't lean on it).

**Nits** (applied unless noted)
- `# Panics` doc on `left_align` — applied. `UNREACHABLE:` marker casing — applied. `result` → `shiftable` rename — applied. `reference_offset as usize` → **won't fix** (ng convention, consistent with A2/B1; the `debug_assert` bounds it). `nonzero` var → **won't fix** (mirrors `build_cigar`, and the debt note points there). `SeqMatch`/`SeqMismatch` handling untested → documented as low-priority (the affine aligner's emitted op subset is unsettled; current fixtures/producers use `Match`).

### 7. Out of scope observations
None.

### 8. Missing tests to add now
All high-value ones applied (12 → 17): Finding-1 consolidation, the documented complex-indel difference, read-consumption round-trip invariant, cap-exactly-shift-distance boundary, zero cap, non-canonical agreement case.

### 9. What's good
- 1b is a genuinely independent left-aligner — the fuzz + grep confirm it leans on nothing it competes with, which is what makes the property oracle's grading of it meaningful.
- Exhaustion is a first-class returned value with `#[must_use]`, the concrete fix for freebayes' swallowed-flag failure.

### 10. Commands to re-verify
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::alignment::left_align_repeated`
