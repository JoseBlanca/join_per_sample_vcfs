# Code Review: ng normalizer — A1 (`AlignmentNormalizer` trait)

**Date:** 2026-07-24
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** working-tree diff for plan step A1 — the `AlignmentNormalizer` trait + `#[cfg(test)]` compile anchors in `src/ng/alignment/mod.rs`
**Status:** Approve-with-changes

---

### 1. Scope
- Reviewed: working-tree diff (step A1 of [alignment_normalization.md](../../ng/impl_plan/alignment_normalization.md)).
- Against: branch `ng-locus-evidence`, uncommitted A1 changes.
- In-scope: the new `pub trait AlignmentNormalizer: Sized` and the new `A1:` test anchors (`IdentityNormalizer`, `LeadingDeletionStripper`, `normalize_via`, tests) in [src/ng/alignment/mod.rs](../../../../src/ng/alignment/mod.rs).
- Out of scope: everything else in mod.rs (plans 1–2).
- Categories dispatched: reliability, errors, naming, idiomatic, refactor_safety, smells, defaults, extras (8 parallel sub-agents). module_structure/unsafe_concurrency/tooling: N/A (single-file, no unsafe/concurrency, no Cargo changes).

### 2. Verdict
Approve-with-changes.

### 3. Execution status
- `cargo fmt --check` → 0, clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean.
- `cargo test --lib` → 0, 2264 passed / 0 failed / 4 ignored (post-fix).
- Findings labeled "Needs verification": 0.

### 4. Open questions and assumptions
1. Doc-test convention for trait-only surfaces without a shipped implementor (affects Mi2). Resolved below in favour of the module's own precedent.

### 5. Top 3 priorities
1. **Mi1** — add the two challenge tests for `LeadingDeletionStripper`'s guard + empty-cigar branches. Applied.

### 6. Findings

**Minor**

- **Mi1: src/ng/alignment/mod.rs — `LeadingDeletionStripper`'s guard and empty-cigar branches untested.** (reliability, High) The anchor's `if let Some(CigarOp::Deletion(dropped))` guard and the empty-cigar no-op are unexercised; a regression that stripped unconditionally or indexed `cigar[0]` would pass the suite. Fix: two challenge tests (guard-skips-non-deletion; empty-cigar boundary). **Applied.**
- **Mi2: src/ng/alignment/mod.rs — public trait lacks a `# Examples` doc test.** (extras, High) The extras "documentation completeness" item asks every `pub` item for a runnable `# Examples`. **Disputed → recorded exemption (Nit).** The two sibling aligner traits (`BestPathAligner`, `MarginalAligner`) — the module's own precedent — carry no `# Examples` doc test; they document via prose + `#[cfg(test)]` anchors. A1 deliberately ships no implementor, so an example would restate the anchor code. The reviewer offered exactly this exemption. Not applied.

**Nits**
- Idiomatic: `normalize`'s two same-typed `&[u8]` params (`read`, `reference`) are positionally swappable — but this mirrors the sibling `MarginalAligner::marginal_probability` signature and the spec-mandated bare-`&[u8]` shape (§7). Whole-trait-family convention; no local fix. Not applied.
- Cross-category (from reliability): `reference_offset += u64::from(dropped)` is an unchecked add — test-only, controlled inputs. No action.

### 7. Out of scope observations
None.

### 8. Missing tests to add now
- `leading_deletion_stripper_leaves_a_non_leading_deletion_alignment_untouched` — guard skips a non-deletion first op. **Added.**
- `leading_deletion_stripper_leaves_an_empty_cigar_untouched` — empty-cigar boundary, no `.remove(0)`. **Added.**

### 9. What's good
- The `LeadingDeletionStripper` anchor exercises precisely the behaviour the whole-`Alignment` signature exists for (moving `reference_offset`) — the load-bearing design point of A1 has a test, not just prose.
- `normalize_via` drives both anchors through a generic bound, so the `Sized`/static-dispatch supertrait is exercised rather than merely declared (the recurring "test that can actually fail" lesson).

### 10. Commands to re-verify
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib`
