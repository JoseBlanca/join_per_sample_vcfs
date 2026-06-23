# Code Review: ssr-call genotyping+pre-pass — Milestone B
**Date:** 2026-06-23
**Reviewer:** rust-code-review skill (orchestrator, focused inline pass)
**Scope:** Milestone B diff (commit `d324a0c`: B1 rung_ladder, B2 stutter, B3 pair_hmm)
**Status:** Approve-with-changes

---

## 1. Scope
- **Reviewed:** the Milestone B diff on branch `ssr-cohort`.
- **In-scope:** [rung_ladder.rs](../../../../src/ssr/cohort/rung_ladder.rs) (B1 additions),
  [stutter.rs](../../../../src/ssr/cohort/stutter.rs) (B2 additions),
  [pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs) (B3, new),
  [mod.rs](../../../../src/ssr/cohort/mod.rs).
- **Out of scope:** A1 types (reviewed), the reading layer.
- **Categories applied:** reliability, errors, naming, idiomatic, refactor_safety,
  smells, extras (numerical correctness / determinism / stable output). Skipped
  unsafe_concurrency (no `unsafe`/threads), tooling (no manifest change).

## 2. Verdict
**Approve-with-changes.** The three primitives are correct on their tested surface and
the math is sound — the kernel sums to 1, the align closed forms are exact, the
heuristic gate hits every documented outcome. Findings are robustness/clarity Minors
plus one numerical-honesty note and a missing test.

## 3. Execution status
- `cargo fmt --check` → pass.
- `cargo clippy --all-targets --all-features -- -D warnings` → pass (a `type_complexity`
  lint was fixed during implementation with a `SeqCount` alias).
- `cargo test --all-features` → **1206 lib pass**, integration + doctests green.
- Findings "Needs verification": 0.

## 4. Open questions and assumptions
1. `align_subst`'s unequal-length branch (M-class candidate, filed Minor) — is the
   slight super-unit total mass across differing-length `obs` acceptable for v1, given
   the equal-length core is exactly normalized and the EM normalizes per read? Filed as
   Minor with a documentation fix; raise to a model decision only if F2 shows drift.

## 5. Top 3 priorities
1. **Mi1** — `is_clear_peak` computes `length + 1`, which panics on overflow in debug
   for a pathological length; use `checked_add`.
2. **Mi2** — document that `align_subst` is exactly normalized only for equal-length
   `obs`/`v`; the flank-gap branch adds a small unnormalized correction.
3. **Mi3** — `Unresolved(Merged)` is returned for *any* zero-clear-peak locus, including
   structureless noise, not only a 1-apart het; document the conflation (or split later).

## 6. Findings

### Minor

- `src/ssr/cohort/rung_ladder.rs` (`is_clear_peak`) — **[Minor]** `length + 1` can overflow
- **Confidence:** High (behavior); Low (reachability — lengths are tract-bounded)
- **Problem:** `histogram.get(&(length + 1))` adds to a `u16`; for `length == u16::MAX`
  this panics in debug / wraps in release. Tract lengths are far below that, so it is a
  latent guard, not a live bug.
- **Why it matters:** Cheap to make total; a wrap would silently read the wrong neighbour.
- **Suggested fix:** `length.checked_add(1).and_then(|u| histogram.get(&u)).copied().unwrap_or(0)`.

- `src/ssr/cohort/pair_hmm.rs` (`align_subst`) — **[Minor]** "probability" is normalized only for equal length
- **Confidence:** High
- **Problem:** For equal-length `obs`/`v` the substitution product is a proper
  distribution (`Σ_obs = 1`). The unequal-length `banded_forward` adds flank-gap paths
  on top of the diagonal mass, so `Σ_obs align_subst(obs|v)` slightly exceeds 1 across
  differing lengths. The doc comment calls it "the probability `P(obs|v)`" without that
  caveat.
- **Why it matters:** In the `Qᵣ = Σ_Δ S_θ·Σ_v align` sum, differing-length terms carry a
  small extra mass; harmless at the tested scale (EM normalizes per read), but the claim
  should match the math.
- **Suggested fix:** note in the doc comment that the equal-length closed form is exactly
  normalized and the flank-gap branch is an unnormalized boundary-slop correction (proper
  transition normalization is an F2 refinement).

- `src/ssr/cohort/rung_ladder.rs` (`resolve_confident_genotype`) — **[Minor]** `Merged` conflates two cases
- **Confidence:** High
- **Problem:** Zero clear peaks at adequate depth returns `Unresolved(Merged)`, which is
  right for a 1-apart het but also captures a structureless / pure-noise locus. The reason
  enum has no "no structure" variant.
- **Why it matters:** Downstream diagnostics that count `Merged` as "masked het" would
  over-count. Both correctly mean "not seeded", so it is cosmetic for now.
- **Suggested fix:** document the conflation on the `Merged` arm; split into a
  `NoStructure` reason only if D1's diagnostics need to distinguish them.

### Nits
- `src/ssr/cohort/rung_ladder.rs` (`representative_sequence`) — the
  `a_c.cmp(b_c).then_with(|| b_seq.cmp(a_seq))` tie-break is correct but hard to read;
  a comment naming the rule (most-supported, then lexicographically smallest) would help.
- `src/ssr/cohort/pair_hmm.rs` (`banded_forward`) — allocates the full `(m+1)·(n+1)`
  buffer; compute is band-limited but memory is not. Fine for v1 (perf is F1); worth a
  note that true memory-banding is deferred.

## 7. Out of scope observations
- `benches/psp_writer_perf.rs:386` — pre-existing bench-harness panic, unchanged.

## 8. Missing tests to add now
- `reach_variants_pure_contraction_shortens_the_tiling` — **input:** a pure `(CA)×5`
  with `Δ = −2`. **Bug it catches:** a sign error or off-by-one in the contraction path
  (the tests cover `+2` and an un-contractable `−3`, but not a *successful* contraction).
  **Body:** assert one variant equal to `(CA)×3`.

## 9. What's good
- `s_theta` is written as the exact scoring inverse of the simulator's draw and a test
  pins `Σ_Δ S_θ = 1` to 1e-12 — the forward/score contract is machine-checked
  ([stutter.rs](../../../../src/ssr/cohort/stutter.rs)).
- `align_subst` keeps the clean majority on a closed form and only drops to DP for the
  rare length-slop case, with `obs == v` short-circuited
  ([pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs)).
- The interior-vs-flank indel test encodes the identifiability invariant (in-tract
  indels not competed) as an executable assertion.
- The merged-het heuristic (zero clear peaks ⇒ adjacent alleles cancelled each other's
  prominence) is a genuinely clever, test-pinned signal
  ([rung_ladder.rs](../../../../src/ssr/cohort/rung_ladder.rs)).

## 10. Commands to re-verify
- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-features`

### Author response convention
Address each finding by ID (Mi1–Mi3, Nits) with `fixed in <commit>` / `disputed` /
`deferred`.
