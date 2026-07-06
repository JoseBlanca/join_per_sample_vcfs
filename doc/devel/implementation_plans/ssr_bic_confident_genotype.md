# SSR pre-pass D1 — BIC confident-genotype gate (implementation plan)

*Build order for spec [ssr_bic_confident_genotype.md](../specs/ssr_bic_confident_genotype.md)
and arch [../architecture/ssr_bic_confident_genotype.md](../architecture/ssr_bic_confident_genotype.md).
Branch `ssr-interruptions` (stacked on Phase-1). Each step: implement (types first,
per `ai/skills/rust-feature-implementation`) → code-review (`ai/skills/rust-code-review`)
→ apply fixes (`ai/skills/apply-code-review-fixes`) → commit → **pause for sign-off**.*

---

## How the order was chosen

- **Types before behaviour.** `GateParams` and the (small) `UnresolvedReason` reshape
  land first, compiling against the old body, so the interface is reviewable before the
  algorithm moves.
- **Gate before attribution before wiring.** The BIC gate (D1a) is self-contained and
  unit-testable on hand-built samples. The attribution swap (D1b) is what makes the
  gate's same-length verdict *matter* for `ε`; it is meaningless without D1a and is
  tested by the `ε`-recovery move. Wiring + determinism (D1c) come last, once both
  behaviours exist.
- **Simulator-gated, not byte-identity-gated.** The win is `ε` de-inflation on injected
  same-length hets and no regression on `ssr_tomato1`; each behavioural step carries its
  own sim assertion (spec §6).
- **Measure the win explicitly (D1d).** A dedicated validation step re-runs `ssr-call`
  on the tomato `.ssr.psp` and reports `ε` and confident-set same-length-het recovery,
  honestly separating D1's pre-pass win from the out-of-scope F_IS downstream tail.

Preconditions already in place (do **not** rebuild): Phase-1 set-valued rungs
(`RungSeq.samples` per-sequence tally, `seqs_at`), `read_likelihood`/`read_given_genotype`
(C2), `nearest_called_by_sequence` (sequence attribution), `LikelihoodScratch`,
`FALLBACK_SHAPE`, `FixedPointAccum`.

---

## The steps

### D1a. Types — `GateParams` + gate signature + enum reshape.  ☐ types
Add `GateParams` (arch §2.2) with `dev_default()` (exposed dev values; `het_admission_cost`
conservative, F2-pinned). Change `resolve_confident_genotype`'s signature to take
`&GateParams` + `&mut LikelihoodScratch` (arch §2.1), *keeping the old heuristic body*
so it compiles and the existing tests pass. Remove `UnresolvedReason::DosageInconsistent`
(model subsumes it — spec §3); adjust its `match`/test sites. **No behaviour change yet.**
*Tests:* existing `rung_ladder.rs` tests updated to the new signature still green;
`DosageInconsistent` removal compiles clean. *Commit.* **Pause.**

### D1b. The BIC gate body.  ☐ arch ☐ plan
Replace the heuristic body with the 1-vs-2-allele BIC test (arch §2.3): top-M candidate
selection from the sample's own sequences, the shared `Qᵣ` score matrix, `argmax_single`
/ `argmax_pair` (fixed-order sum, bytes/index tie-break), the `het_admission_cost·ln(n)`
admission, and the sequence-keyed recurrence guard (`all_recurrent`, arch §2.4). Keep the
depth skip. *Tests (arch §6):* homozygote → 1 peak; length-separated het → 2 peaks;
**same-length het → 2 same-length peaks, distinct bytes** (the new capability); 1-apart
balanced → 1 peak ("contributes nothing"); hom+heavy-stutter → 1 peak (dosage subsumed);
non-recurrent same-length minority → `NonRecurrent`; thin → `Thin`. *Review → apply →
commit.* **Pause.**

### D1c. Sequence-aware attribution + wiring + determinism.  ☐ arch ☐ plan
Swap `accumulate_locus`'s read attribution from `nearest_parent` to
`nearest_called_by_sequence` (arch §3), threading `seed.eps` + the scratch. Thread
`GateParams` + a per-thread `LikelihoodScratch` through `run_prepass_stats`/`run_prepass`
(arch §4); the integer `reduce` is unchanged. Driver/D2-seam call site passes
`GateParams::dev_default()`. *Tests:* on an injected same-length-het simulator with known
`ε`, recovered `ε` matches truth where the pre-swap heuristic over-estimates it; the
existing checkpoint-2 recovery + `separated_hets_contribute_two_length_bins` tests still
pass (length-separated effect-neutral); **`--threads 1` vs `K` byte-identical** params
(extend the existing `prepass_is_byte_identical_across_thread_counts`). *Review → apply →
commit.* **Pause.**

### D1d. Validate on `ssr_tomato1`.  ☐ plan
Re-run **only** `ssr-call` on the existing per-sample `.ssr.psp` (Stage 1 unaffected;
`DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks`, in-container). Report,
against the Phase-1 baseline:
1. **frozen `ε`** — expect a *fall* on this interruption-rich cohort (the direct win);
2. **confident-set same-length-het placement** — the pre-pass now seeds same-length hets
   (P1.5 `seq_concordance.py`, `statSTR --use-length` off);
3. **no regression** — length-genotype concordance (96.5 %) and SNP e2e tests hold.
Write `doc/devel/reports/ssr_bic_confident_genotype_validation_2026-07-…md`, stating
plainly which part of the 46 % downstream het tail D1 addresses (pre-pass ε/θ) and which
it does not (the F_IS genotyping-prior tail, roadmap E1/E2 — spec §6). *Commit report.*
**Pause.**

---

## Exit criteria (honest, per spec §6)

- **Committed:** injected-`ε` recovery correct on the same-length-het simulator (the
  heuristic over-estimates, the BIC gate recovers truth); `ssr_tomato1` frozen `ε`
  de-inflates; **no regression** on common-loci length concordance or SNP e2e tests;
  `--threads 1..K` byte-identical pre-pass params.
- **Reported, not gated on D1 alone:** any rise in downstream same-length-het
  concordance — welcome, but partly gated by the out-of-scope F_IS prior; D1 is not
  failed if the full 46 % tail does not close.
- **Untouched, verified:** `em.rs`, `candidate_set.rs`, `G₀`, `inbreeding.rs`,
  `vcf_out.rs`, the pre-pass reduce, the sufficient-statistic structs, Stage 1 (spec §8).

## What this plan is *not*

- **Not** the burn-in loop (roadmap D2) — D1 scores the gate with a **coded seed**;
  co-evolving it with the pre-pass's own fitted params is D2. The `GateParams` call site
  is D2's seam (arch §4).
- **Not** a change to the read likelihood, the candidate set, the genotype prior, or the
  genotyping EM — the gate reuses `Qᵣ` and the sample's own sequences (spec §8).
- **Not** F2 calibration — dev defaults ship; `het_admission_cost`, seed `ε`/level/`λ`,
  recurrence `k`, min-depth are pinned on the simulator in F2.
