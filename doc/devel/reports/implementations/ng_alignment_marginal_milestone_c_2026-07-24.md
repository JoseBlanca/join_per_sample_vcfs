# ng alignment — marginal aligners, Milestone C (algorithm 6 + the synthetic-truth proof)

**Date:** 2026-07-24
**Plan:** [alignment_marginal.md](../../ng/impl_plan/alignment_marginal.md) (plan 2 of 3), steps C1–C2, at **Checkpoint C — plan complete**.
**Design authority:** spec [alignment.md](../../ng/spec/alignment.md) §9, §5, §7, §10.3; arch [alignment.md](../../ng/arch/alignment.md) §5.
**Method:** plan-driven — one implement → review → apply-fixes → commit loop per step. Per-category review outputs in gitignored `tmp/review_2026-07-24_ng-marginal-c{1,2}/`.

---

## 1. What landed

- **C1** (`7829eb5`) — algorithm 6, `SsrWholeReadMarginal` in [ssr_marginal_whole_read.rs](../../../../src/ng/alignment/ssr_marginal_whole_read.rs): a whole-read two-regime forward pair-HMM. Scores the whole read against a full flank-repeat-flank reference in one forward pass, marginalizing the repeat length rather than measuring it. **New code, no port oracle.**
- **C2** (this commit) — the synthetic-truth proof (tests only): each marginal (5 and 6, separately) scores its own generating allele highest.

## 2. The C1 design decision (owner-confirmed at Checkpoint B)

Algorithm 6 was the one genuinely under-specified step: the spec fixes the *what* (whole read vs flank-repeat-flank, one flat rate, fixed synthetic qualities, repeat-aware, returns `LogProb`, `Context` = geometry) but not the exact gap recurrence, and there is no production counterpart. **Surfaced to the owner via AskUserQuestion; the owner chose the two-regime forward** — the sum-reduction analog of the tract delimiter (algorithm 3): stiff flank gaps, soft tract gaps, summed (forward) not maximised (Viterbi). The two regimes are what make it repeat-aware and let a read score highest against its own allele.

**Reuse over rewrite** (arch §5): it reuses production's `TransitionCosts` (algorithm 3's exact two-regime log-space gap model, via its `pub(super)` accessors — the same seam algorithm 4 uses, so **no edit to algorithm 3**), inheriting the provisional gap-open calibration and the documented match→match inconsistency by reuse; and the module's log-space `FlatEmission` (one flat rate + a fixed synthetic quality — which for algorithm 6 genuinely *is* the "reuse the Emission component" the plan asked for, unlike algorithm 5). It runs a **log-space** 3-state (M/I/D) forward with a stable log-sum-exp, so it produces `LogProb` directly. **Unbanded** (full matrix), matching how algorithm 3 first shipped.

**The defining contrast with algorithm 5:** algorithm 5 *forbids* interior gaps (they are the stutter model's to explain); algorithm 6 *permits* soft tract gaps precisely so it can marginalize the length without a prior measurement.

## 3. C2 — each algorithm computes what it claims

Deliberately **not** a 5-vs-6 bake-off — that scores a measured repeat against a whole read and spans the genotyping (spec §10.3). C2 proves each marginal *in isolation*: a read from a k-unit allele (with two real substitution errors) scores higher than every neighbour (k±1, k±2), by a **clear margin** (not merely argmax), **stably across the flat error rate** (two ε) and **across allele length** (two truths). Algorithm 5's far neighbours (k±2 = ±6 bases, beyond the flank slop) are unreachable → −∞ by design; algorithm 6's soft tract gaps make them reachable but dearer.

## 4. Verification — no oracle, so pinned by reasoning + anchors

Algorithm 6 has no byte-parity oracle. Its correctness rests on:
- the recurrence being **hand-traced correct cell-by-cell** against algorithm 3's Viterbi (the review did this);
- two **absolute anchor tests** (added at C1): a single-cell case pinning `match_ln + ln_match_to_match`, and a two-path case where a `max` (a Viterbi masquerading as a forward) would score ~0.49 nats too low — the one property separating a marginal from a best path;
- a direct **tract/flank junction** test (`is_tract_column` extracted to a testable free fn);
- the C2 synthetic-truth proof above.

## 5. Review outcomes

- **C1** (3 sub-agents / 8 checklists): recurrence confirmed correct. Review's **Blocker** was the recurring "tests that cannot fail" — all tests were relational, so a monotonicity-preserving regression (notably `f64::max` in the forward's place) passed them all, silently corrupting downstream likelihood ratios. Applied: the two absolute anchors, the junction test, an empty-input test, `is_tract_column` extraction, saturating `usize::try_from` (sibling-consistent, 32-bit safety), noun-form scratch fields, and a doc refresh. The sibling-to-sibling `TransitionCosts` reuse is an accepted seam (2 consumers; lift to a peer file if a 3rd appears).
- **C2** (1 test-validity sub-agent): **Major** — both tests *claimed* two substitution errors but the second mutation was a **no-op** (it wrote the base already present), so each exercised only one error. Applied: real second errors, a strict-margin assertion (over the best neighbour, not bare argmax), a hoisted scratch, and a second truth length.

Every finding across both steps was applied. No stop-and-ask arose in C (the C1 Blocker was test-coverage, not a design flaw); the one genuine design fork was resolved by the owner at Checkpoint B before C1.

## 6. Validation (container, via `./scripts/dev.sh`)

- `cargo fmt --check` → clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` → clean.
- `cargo test --lib` → **2260 passed, 0 failed**, 4 ignored.

The two pre-existing, unrelated red project-wide commands (`cargo test --all-targets` bench panic; `cargo doc --no-deps` intra-doc links) are the documented Standing items; untouched.

## 7. Plan complete

Both marginals — algorithm 5 (a byte-parity port of `align_subst`) and algorithm 6 (the whole-read two-regime forward) — run on synthetic truth and each scores its own generating allele highest. **Plan 2 of 3 (`alignment_marginal.md`) is complete.** The 5-versus-6 comparison is handed to the genotyping, which owns the harness that can run it (spec §10.3). Plan 3 ([`alignment_normalization.md`](../../ng/impl_plan/alignment_normalization.md)) — the normalizer interface and algorithms 1a/1b/1c — remains, independent of this work.

**Checkpoint C — hard pause. Plan 2 complete.**
