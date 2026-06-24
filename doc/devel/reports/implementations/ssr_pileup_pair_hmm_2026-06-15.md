# SSR Stage-1 — `pair_hmm` slow-path forward scorer

**Date:** 2026-06-15 · **Skill:** rust-feature-implementation · **Branch:** `ssr-architecture`

The stage's central net-new numerical machinery
([plan §5](../../../doc/devel/implementation_plans/ssr_pileup.md), arch §5) — and
its flagged main risk. Self-contained (byte slices + quals + constants), so built
ahead of triage/candidate wiring and tested in isolation.

## 1. Plan

`forward(read, quals, candidate_seq, scratch, model) -> f64`: a 3-state
(Match/Insertion/Deletion) log-space forward giving `ln P(read | candidate_seq)`
summed over alignments, under sequencing error only (no stutter). Pure scorer, no
traceback → a rolling two-row scratch. Dindel per-Q emission; affine gaps.

## 2. Assumptions / deferrals (surfaced)

- **Constant gap-open, not the homopolymer-run-indexed table.** Arch §5.4 adopts
  Dindel's per-run-length gap-open verbatim; I used a single constant
  (`GAP_OPEN_PROB = 2.9e-5`, Dindel's short-run base). Arch §14 lists the
  homopolymer extrapolation as a *calibration* item — the structural forward is
  what's load-bearing; the refinement swaps the constant for a position-indexed
  lookup. Documented in the module.
- **Unbanded full `O(read × hap)` DP.** Banding (`PAIR_HMM_BAND_BP`, arch §5.5) is
  an optimization to add once the slow path is shown to bind; the unbanded forward
  is the correct `W → ∞` limit (correctness-first). Documented.
- **Insertion emission = `ln(1/4)`** (Durbin background `qₐ`); a constant per
  inserted base, so it only shifts scores when candidates differ in insertion
  count. **Match `begin` state = `Match[0][0]`** (no separate begin/end τ);
  leading/trailing indels reached via the I/D boundary cells.
- **Q0 match floored** to the smallest positive `f64` (ε=1 ⇒ `ln 0` would
  annihilate the row); Q0 is pathological, real data is Q≥2.
- **Scope: the `read × haplotype → Qᵣ` scorer only.** No candidate-set wrapper
  (`score_candidates` needs `CandidateAllele` from `candidate_generation`), no
  shared-flank lattice.

## 3. Changes made

- **New** [src/ssr/pileup/pair_hmm.rs](../../../src/ssr/pileup/pair_hmm.rs):
  `forward`, `HmmModel` (precomputed log transitions), `PairHmmScratch`
  (grow-and-keep two rows), `EMISSION_LN`/`INS_EMIT_LN`/`GAP_EXTEND_PROB`
  process-wide lookups (the BAQ `Q2P` *pattern*, not its code), `ln_sum_exp2/3`.
- [src/ssr/pileup/mod.rs](../../../src/ssr/pileup/mod.rs): `pub mod pair_hmm;`.

## 4. Tests added (8, in-module)

Two **hand-computable exact-value** cases (single-base match → `ln(1−ε₄₀) +
ln(1−2·gap_open)`; single-base mismatch → `ln(ε/3) + …`); match ≫ mismatch;
exact-length haplotype beats a 2-bp-longer one; longer exact match is more
negative (more `(1−ε)` factors) yet both ≤ 0; low-Q mismatch less penalized than
high-Q; empty read scores a finite all-deletion path; scratch reuse across a long
then short haplotype stays correct.

## 5. Validation results

Run in the dev container (`./scripts/dev.sh`):
- `cargo fmt -- --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ssr::pileup::pair_hmm` — **8 passed, 0 failed** (1075 others
  filtered out, untouched).

## 6. Tradeoffs and follow-ups

- Homopolymer-indexed gap-open + banding deferred (calibration/optimization, arch
  §14); the constant-gap unbanded forward is correct and is the de-risking target.
- `score_candidates` (the per-read candidate-set wrapper returning the `Qᵣ`
  distribution) lands once `candidate_generation` supplies `CandidateAllele`.
- **Container generalization (task 3) still open.**
