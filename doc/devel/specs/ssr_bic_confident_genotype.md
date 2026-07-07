# SSR pre-pass D1 — model-based confident-genotype resolution (the BIC 1-vs-2-allele test)

*Status: proposal, 2026-07-06, branch `ssr-interruptions` (stacked on the Phase-1
interrupted-repeat recall). Companion to [`ssr_cohort_mark2.md`](ssr_cohort_mark2.md) §4.3/§4.4 (the settled
intent) and to [`ssr_interrupted_repeat_recall.md`](ssr_interrupted_repeat_recall.md)
§5.5 (why Phase 1 unlocks this but does not deliver it). Code-layout companion:
[`architecture/ssr_bic_confident_genotype.md`](../architecture/ssr_bic_confident_genotype.md);
build order:
[`implementation_plans/ssr_bic_confident_genotype.md`](../implementation_plans/ssr_bic_confident_genotype.md).
Where this and Mark-2 disagree on intent, Mark-2 wins; this doc only *sharpens the
gate* Mark-2 §4.3 already settled at the "one-allele-vs-two-allele, not a generic
goodness-of-fit" level (decision Q-P7).*

---

## 0. Glossary (the terms this doc leans on)

- **Pre-pass** — the pass that estimates the frozen chemistry parameters (per-base
  error `ε`, stutter shape, stutter level) from the cohort's *confident genotypes*,
  before per-locus genotyping runs. Lives in
  [`prepass.rs`](../../../src/ssr/cohort/prepass.rs); its gate lives in
  [`rung_ladder.rs`](../../../src/ssr/cohort/rung_ladder.rs).
- **Confident genotype** — a (sample, locus) whose genotype the pre-pass trusts
  enough to *label every read* against and pool into the chemistry estimate: a
  homozygote (one allele) or a well-resolved heterozygote (two alleles). Mark-2
  calls this the **CG-seed**.
- **`ε` (per-base error)** — the within-tract substitution rate: given a read sits
  at a called allele's length, the fraction of tract bases that differ from that
  allele. Measured off confident genotypes' faithful reads.
- **Stutter** — a PCR replication slip that changes the tract by a whole number of
  repeat units (a *length* change). `Qᵣ` (below) models it.
- **`Qᵣ(obs | cand)`** — the read likelihood: the probability an observed tract
  `obs` arose from a candidate allele `cand` via stutter (length change) plus base
  error (composition change). The production model is HipSTR-style
  ([`read_likelihood`](../../../src/ssr/cohort/likelihood.rs), roadmap C2). This is
  the likelihood the BIC test reuses — **we do not reinvent a scoring function**.
- **BIC** — the Bayesian Information Criterion, a model-selection rule that rewards
  fit (log-likelihood) and penalises free parameters (`k·ln n`). Here the two models
  are "one allele" vs "two alleles"; the second allele must earn its extra parameter.
- **Same-length het** — a heterozygote whose two alleles have the *same tract
  length* but differ in *composition* (e.g. a pure `(CA)₉` and an interrupted
  `(CA)₄CG(CA)₄`). The case a length histogram structurally cannot see.

---

## 1. What is wrong today, and the evidence for it

The pre-pass picks its confident genotypes with a **length-histogram heuristic**
([`resolve_confident_genotype`](../../../src/ssr/cohort/rung_ladder.rs)): it looks
for clear local maxima in a sample's *length* histogram (one repeat-unit length
standing prominently above its ±1 neighbours), keeps the top `ploidy` of them, and
guards the result with separation (peaks ≥ 2 units apart), dosage balance, cohort
recurrence, and a depth floor. One clear length peak ⇒ homozygote for that length's
dominant sequence.

That gate is blind to composition, so it **structurally mislabels a same-length
heterozygote as a homozygote**. A sample that is genuinely `pure-18 / interrupted-18`
shows a *single* length peak; the heuristic calls it homozygous for the majority
composition (`pure-18`), and the minority allele's on-length reads — same length,
different bases — are then scored as the majority allele's *base-error halo*. Each
differing base is counted as a substitution, so `ε` is **inflated**, concentrated on
interruption-rich cohorts. The pre-pass code already flags exactly this in its
`BIAS NOTE` ([`prepass.rs:223`](../../../src/ssr/cohort/prepass.rs)).

This is not hypothetical. The Phase-1 validation
([`ssr_interrupted_repeat_p1_validation_2026-07-06.md`](../reports/ssr_interrupted_repeat_p1_validation_2026-07-06.md))
measured the tail on the tomato cohort with the sequence-genotype metric
(`statSTR --use-length` off): of the HipSTR same-length-het cells now visible,
Phase 1 recovered **54 %** as same-length hets (up from 26 % pre-Phase-1) — leaving a
**46 % undercall tail** still called homozygous. Phase 1 was the prerequisite (it
made the caller *able* to nominate two same-length sequences at all); this gate is
the change that turns that ability into cleaner chemistry.

**Why the heuristic cannot be patched in place.** The reason it needs a separation
guard and a dosage-balance guard at all is that a length histogram cannot tell an
allele from its own stutter skirt: a −1 stutter satellite and a real −1 allele look
identical. The fix is not another histogram rule — it is to score reads under a model
that *predicts the stutter skirt* (`Qᵣ`), so "is this a second allele or just
stutter?" becomes a likelihood question, and "two alleles at one length" becomes
askable. That is the Mark-2 §4.3 model-based test, decision Q-P7.

---

## 2. The model form — one allele vs two, scored with `Qᵣ`  *(settled at the §4.3 level; this is the concrete form)*

For one (sample, locus) with adequate depth, the gate compares two hypotheses over
the sample's reads, both scored with the **existing** read likelihood `Qᵣ` and the
existing genotype-likelihood mixer
([`read_given_genotype`](../../../src/ssr/cohort/likelihood.rs)):

- **One-allele model `M₁` (homozygote).** Pick the single allele `A` that best
  explains all the reads. Its genotype is `(A, A)`, so `P(read | M₁) =
  Qᵣ(read | A)` (mixed with the uniform outlier floor `λ/D`, exactly as the EM does).
  Free parameters: **one** allele identity, `k₁ = 1`.
- **Two-allele model `M₂` (heterozygote).** Pick the pair `(A, B)` that best explains
  all the reads. `P(read | M₂) = ½·[Qᵣ(read | A) + Qᵣ(read | B)] + λ/D` (the diploid
  mixing proportion is fixed at ½ — it is **not** a free parameter). Free parameters:
  **two** allele identities, `k₂ = 2`.

The candidate alleles `A`, `B` are drawn from the **sample's own observed distinct
sequences** (the tract sequences in its `seq_counts`, already byte-sorted by the
Stage-1 contract), not the full cohort candidate set — the gate asks "which one or
two of the sequences this sample actually produced best explain its reads." Phase 1
is what makes two *same-length* sequences available to be picked here.

Both models are scored by their **maximised total log-likelihood** over the reads:

```text
ln L̂₁ = max_A   Σ_reads  count · ln[ (1−λ)·Qᵣ(read | A) + λ/D ]
ln L̂₂ = max_(A,B) Σ_reads count · ln[ (1−λ)·½·(Qᵣ(read|A)+Qᵣ(read|B)) + λ/D ]
```

(`count` is the read's multiplicity in `seq_counts`; `D` is the number of distinct
sequences at the locus.) Because `M₁ ⊂ M₂` (set `B = A`), `ln L̂₂ ≥ ln L̂₁` always;
the question is only whether the gain is worth the extra allele.

### 2.1 The BIC decision and the purity tuning

Admit the second allele — call the sample a **confident het** — iff the
log-likelihood gain beats a complexity penalty:

```text
2·(ln L̂₂ − ln L̂₁)  >  penalty · ln(n)          (n = total reads scored)
    ⇔   ln L̂₂ − ln L̂₁  >  het_admission_cost · ln(n)
```

with `het_admission_cost = penalty / 2` exposed as a single tunable coefficient.
Standard BIC for one extra parameter would put `penalty = k₂ − k₁ = 1` (so
`het_admission_cost = ½`). **We deliberately run it higher** — the penalty is
**tuned conservatively for purity** (Mark-2 §4.3, decision Q-P7): a hidden het
sneaked into the confident set poisons `ε`/`θ` for the whole cohort, whereas a real
homozygote wrongly discarded merely costs a little data. So the dev default sets
`het_admission_cost` well above ½ (a provisional dev value pinned on the simulator in
F2, §6 below), and the calibration target is explicitly asymmetric: **let almost no
false het through, even at the cost of discarding real homozygotes.**

The comparison is a pure function of the locus's reads and the seed params (§4); it
is computed entirely within one locus, so it never crosses a thread boundary (§5).

### 2.2 The same-length het — the case this delivers

This is the case the length heuristic cannot see and the reason D1 exists. Take a
`pure-18 / interrupted-18` sample:

- Under `M₁ = (pure-18, pure-18)`: the `pure-18` reads score as faithful (`Δ = 0`,
  no mismatch); the length skirts score as `pure-18`'s stutter; but the
  `interrupted-18` reads score as `Qᵣ(interrupted-18 | pure-18)` — same length, so
  `Δ = 0`, but several substitutions, each costing a factor ≈ `(ε/3)/(1−ε)`. They
  are explained only weakly (the base-error tail), so they drag `ln L̂₁` down.
- Under `M₂ = (pure-18, interrupted-18)`: the `interrupted-18` reads now score as
  *faithful* to their own allele — a large per-read likelihood jump (§5.5 of the
  interruption spec puts it at ≳ Phred 20 per read even at an inflated `ε ≈ 0.02`).

So `ln L̂₂ − ln L̂₁` is large and the second allele is admitted. Crucially this
signal comes from *composition* (`Δ = 0` for both alleles), so it is **robust to a
rough seed stutter level** — the very failure mode that makes 1-apart hets delicate
(§2.3) does not touch the same-length case. This is why D1 both *can* and *safely
does* recover same-length hets.

### 2.3 The 1-apart merged het — the case this must still reject  *(design corrected in D1b — see below)*

There are **two** distinct "1-apart" situations, and they are handled by two different
mechanisms. Getting this wrong is easy, so it is spelled out.

**(i) A homozygote with heavy −1 stutter** (`A`, with a −1 satellite from stutter, no
real second allele). Under `M₁ = (A, A)` the satellite reads are *already explained* as
`A`'s down-stutter, so `M₂ = (A, A−1)` earns little. Crucially, `M₂` also **pays the
diploid mixing penalty on every faithful `A` read**: a het scores each majority read as
`½·(Qᵣ(read|A) + Qᵣ(read|A−1)) ≈ ½·Qᵣ(read|A)`, i.e. ≈ `ln ½ ≈ −0.69` nats *per faithful
read*. With most reads faithful, that penalty dominates the small satellite gain, so the
BIC test **keeps one allele** → a confident homozygote. This is what "a hom+heavy-stutter
contributes as a hom, not a het" means, and it is the mirror of §2.4 — the mixing penalty
is the same force that stops an error halo from inventing an allele.

**(ii) A genuine balanced `(A, A−1)` heterozygote** (~half the reads at each). Here the
`A−1` reads are far more than stutter predicts, so the two-allele gain *overwhelms* the
mixing penalty and the **BIC test correctly admits two alleles** — it *is* a real het.
But a 1-apart het is **not a clean chemistry seed**: `A`'s down-stutter lands on `A−1`
and `A−1`'s up-stutter lands on `A`, so the two alleles' skirts are entangled and cannot
be cleanly labelled (Mark-2 §4.3 requires a seed het's alleles to be **≥ 2 units apart**).
So D1 keeps an explicit **length-separation clean-seed guard**: when BIC admits a
length-separated pair closer than `separation_min` (dev 2), the sample is left
`Unresolved(Merged)` rather than seeded — we neither seed it as a het (skirts entangled)
nor collapse it to a homozygote (that would bury a real allele's reads in `A`'s stutter
and re-inflate the level). This guard is **scoped to length-separated pairs only**: a
**same-length** het (§2.2) has a length gap of 0 and bypasses it entirely, because its
two alleles are told apart by *composition*, not length, so there is no skirt overlap.

> **Design correction (D1b, 2026-07-07).** An earlier draft of this section claimed the
> separation guard was fully *subsumed* by "the model + the penalty" (that a 1-apart het
> would simply fail the BIC penalty). That is wrong for case (ii): the BIC test *admits*
> a genuine balanced 1-apart het (it really is two alleles), so a separation guard is
> still needed — not to decide "is it a het" but to decide "is it a *clean seed*." The
> guard survives, scoped to exclude same-length pairs. Only the **dosage-balance**
> heuristic is fully subsumed (case (i) is handled by the mixing penalty). §3's table row
> is corrected to match.

### 2.4 A low-frequency error vs a real same-length allele — the mirror bias, and why the model does not fall into it

This is the sharpest risk in the whole design, and its exact mirror image of the bug
we are fixing. The **heuristic** errs one way: it labels a same-length minority
allele's reads as *base error* and so **inflates `ε`**. A **naïve** likelihood test
would err the other way: label a cluster of `ε`-error reads as *a second allele*, and
so **inflate the allele count and deflate `ε`** (every error read reclassified as some
allele's faithful read — and, worse, self-reinforcing, since the deflated `ε` makes the
next such call easier). D1 must land *between* these, and it does so on three
mechanisms, in increasing order of decisiveness.

**(a) The one-allele model already expects the error halo — so error reads do not
force an allele.** `M₁`'s `Qᵣ(B | A)` for a same-length read `B` differing from `A` at
`m` positions is ≈ `(1−level)·(1−ε)^(L−m)·(ε/3)^m`: M₁ *predicts* a thin halo of
near-`A` reads at rate `ε`. A handful of them is exactly what M₁ expects; they earn the
second allele nothing. The question is never "are there any non-`A` reads" (there
always are) but "are there *far more, and far more concentrated,* than `ε` predicts."

**(b) `M₂` can add only *one* allele, so it pays off only for a *spike*, not a
*spread*.** This is the structural crux. Random sequencing error scatters the
non-faithful mass across ~`3L` different single-base neighbours of `A`, each at
frequency ~`ε/3` — a **flat spread**. A real same-length allele concentrates ~half the
reads on **one** composition `B` — a **spike**. Because the two-allele model admits a
single extra allele, it can capture the spike (large `ln L̂₂ − ln L̂₁`) but only
~`1/(3L)` of the spread (a gain far below the penalty). So the gate discriminates on
the *shape* of the excess mismatch mass — concentrated (allele) vs dispersed (error) —
**not on its magnitude alone**, and this shape argument is *largely robust to the exact
seed `ε`*: getting `ε` somewhat wrong rescales both hypotheses' halos together and does
not turn a flat spread into a spike. The candidate search reinforces this: candidates
are the sample's **top-M sequences by read support** (§2, arch §2.2), so a scattered
error variant — low count by construction — never even enters the pair search, while a
~50 %-dosage `B` ranks high.

**(c) Cohort recurrence is the decisive backstop against `ε` deflation.** Even a
per-sample spike can be a *systematic* artefact (a homopolymer miscall, a mapping
edge). The sequence-keyed recurrence guard (§3) is what makes the design safe: a real
interruption allele recurs at the **same composition** across ≥ `k` carriers (§5.2
systematic signal), whereas a sporadic error lands at a **different base each time** and
never recurs. An error-driven "allele" is therefore rejected `NonRecurrent`, so it never
enters the confident set and never touches `ε`. A **private** (single-sample)
same-length variant is *deliberately excluded* from the chemistry seed for exactly this
reason — indistinguishable from an artefact without orthogonal data, so the pre-pass,
which wants purity not recall, must not seed on it (the genotyping EM downstream is a
separate stage and may still call it). This is the same recall/purity trade the
`min_same_length_reads` / distinct-sample admission already makes at the Phase-1
candidate level.

**What remains genuinely ambiguous — and is treated honestly.** A same-length variant
that (i) concentrates like an allele *and* (ii) recurs like an allele *is* seeded as
one. If it is in truth a *systematic chemistry artefact* that reproduces across
samples, no single-cohort caller can distinguish it from a biological allele without
orthogonal evidence (HipSTR has the same limit). We accept this: it is a far smaller
and rarer error than the systematic `ε` inflation the heuristic commits on *every*
interruption-rich locus, and the penalty is tuned (§2.1) so that only strongly-supported,
strongly-recurrent spikes clear it. The residual is surfaced, not hidden — the paralog /
universal-het FP audit (interruption spec §10, already run in P1.5) is the cross-cutting
check that a same-length het is not being called in *every* carrier (the collapsed-CNV
signature).

**The calibration that pins this (F2).** The penalty's operating point is set on a
simulator built for *this* failure: a cohort of true homozygotes with an injected `ε`
(sweeping `ε` up to interruption-rich levels) and **no** extra same-length alleles. The
gate must invent **zero** alleles and recover the injected `ε` **undeflated**; the
largest `het_admission_cost` that still recovers the genuine same-length hets of §2.2
without inventing alleles here is the dev default. This measures the mirror bias
directly, rather than trusting the argument above.

---

## 3. Which guards survive, which the model subsumes

The heuristic's four guards were length-histogram proxies for "is this really a het."
Under D1 the model **subsumes one** (dosage balance), **re-scopes one** (separation, to
length-separated pairs only — D1b correction, §2.3), and **keeps two** (min depth, and
cohort recurrence generalised to sequence):

| heuristic guard | fate under D1 | why |
|---|---|---|
| **separation ≥ 2 units** | **kept, but re-scoped to length-separated pairs only** | Still needed as a *clean-seed* guard: the BIC test *admits* a genuine balanced 1-apart het, but a 1-apart het's skirts are entangled and cannot be cleanly labelled (Mark-2 §4.3), so it is left `Unresolved(Merged)` rather than seeded. **Same-length** hets (gap 0) bypass the guard — composition, not length, separates them (§2.2/§2.3). |
| **dosage balance** (minor ≥ ratio·major) | **dropped (subsumed)** | A homozygote-plus-heavy-stutter no longer needs a height rule: `Qᵣ` *predicts* the stutter satellite (so `M₁` explains it) and `M₂` also pays the diploid mixing penalty on every faithful read, so the second allele earns no net gain. The model is the balance test (§2.3 case (i)). |
| **min depth** | **kept** (a pre-BIC skip) | Too few reads to distinguish one allele from two; skip, do not guess (unchanged `Thin`). |
| **cohort recurrence** | **kept, generalised to sequence** | An allele the BIC test admits must still recur across samples — a per-carrier substitution artefact is sporadic, a real interruption allele recurs (Phase 1's per-sequence sample tally, `RungSeq.samples`, supplies this for same-length alleles; the existing per-length `peak_recurrence` covers length-separated ones). This is the cohort-level false-allele defence the single-sample likelihood cannot provide, and a second purity lever. |

The `UnresolvedReason` enum keeps `Thin` and `NonRecurrent`; `Merged` now means "BIC
admitted a length-adjacent pair too close to seed cleanly" (§2.3 case (ii));
`DosageInconsistent` is **removed** (the model subsumes it). A same-length het that
fails recurrence is `NonRecurrent`, not silently dropped.

---

## 4. What params the gate scores with — the D1/D2 boundary  *(the one genuine scope decision — flag for sign-off)*

`Qᵣ` needs `ε`, a stutter shape, and a stutter level to score. The pre-pass is what
*measures* those — a bootstrap. Mark-2 §4.3/§4.4 resolves this by having the gate
"co-evolve with the parameters in the **burn-in loop**": the gate re-runs as the
parameters settle. **That burn-in loop is roadmap D2, explicitly out of D1's scope**
(the roadmap dependency graph is D1 → D2; the prompt names the burn-in as a piece D1
must not touch).

So D1 scores the gate with a **coded seed parameter set** — the "burn-in start" §4.3
already names: the existing `FALLBACK_SHAPE` (up ½ / down ½ / decay 0.1), a seed `ε`,
a seed stutter level, and a seed outlier `λ` (literature/dev defaults, exposed, no
hidden constants; pinned in F2). The gate is **prior hygiene, low-stakes** — Mark-2
§4.3 is explicit that the EM can overrule a biased prior — and, decisively, the case
D1 targets (the same-length het, §2.2) is *robust to a rough seed* because its signal
is composition, not length. A rough seed is therefore adequate to deliver the ε
de-contamination D1 is for.

The single pass this keeps (heuristic → BIC, both one sweep) matches D1's roadmap
shape. **The data-driven co-evolution — re-run the gate with the pre-pass's own fitted
params, one or more iterations — is the documented D2 upgrade**, noted here so the
boundary is explicit and so D2 has a clean seam to extend (the gate is already a pure
function of `(reads, params)`; D2 only changes *which* params it is handed and *how
many times*). **Sign-off point:** if you want D1 to instead do a two-pass bootstrap
(heuristic-seed → BIC-refine) so the gate scores with data-driven rather than coded
params, that is a small, self-contained change to the pre-pass wiring — but it starts
to blur the D1/D2 line, so I have scoped D1 as single-pass-with-coded-seed and left
the two-pass as the D2 entry point. Confirm which you want.

---

## 5. The determinism story

SSR output must be byte-identical across thread counts. D1 preserves this by
construction:

- **The BIC decision never crosses a thread boundary.** `run_prepass_stats` maps each
  locus to exactly one thread (`par_iter().fold`); the whole gate — the candidate
  search, both models' log-likelihoods, the comparison — runs within `accumulate_locus`
  for that locus. So the float log-likelihoods are computed atomically per (sample,
  locus) and are identical run-to-run and thread-to-thread. **Fixed sorted-order f64
  summation suffices here; `FixedPointAccum` is only required for floats reduced
  *across* threads** (the D2 burn-in `ℓ_pen` barrier — not built here). The
  determinism contract's phrase "BIC log-liks via FixedPointAccum" refers to that
  cross-thread D2 reduce; D1's within-locus decision meets the same guarantee with
  sorted-order summation, and the plan proves it with a `--threads 1 vs K` byte-
  identity test.
- **Every reduction is fixed-order.** Reads are summed in `seq_counts` order (byte-
  sorted, thread-independent). The candidate search iterates the sample's distinct
  sequences in byte order.
- **Every tie breaks on bytes / lower index.** The best-allele and best-pair argmax
  break likelihood ties on the lexicographically smaller sequence (then lower
  candidate index) — never on float-equality drift. The read→allele attribution in
  `accumulate_locus` uses the existing `nearest_called_by_sequence`, whose same-length
  tie-break is already the shared, documented rule.
- **The sufficient statistics stay integer.** The pre-pass accumulators
  (`SlipProfile`, `SampleStutterStats`, `purity_slip`, the allele-spread sums) are
  unchanged integer counts; D1 only changes *which reads land in which bin* (via the
  gate's verdict and sequence-aware attribution, §7), not the reduce.

---

## 6. Success metrics and calibration  *(the regression gate is concordance, not byte-identity)*

Pre-alpha, SSR output byte-identity is **not** the gate; the regression gate is the
`ssr_tomato1` benchmark concordance plus the SNP caller's end-to-end tests. D1's
*win* is measured on two axes:

1. **`ε` de-inflation on interruption-rich cohorts (the direct, attributable win).**
   On the tomato cohort (and a purpose-built interruption-rich simulator), the frozen
   `ε` the pre-pass emits should *fall* once same-length hets stop dumping their
   minority reads into the halo. This is the cleanest measurable D1 effect and the
   primary acceptance signal. A simulator with injected same-length hets and known
   `ε` is the ground-truth check (recover the injected `ε`, where the heuristic gate
   over-estimates it).
2. **Same-length-het recovery in the confident set (the enabling win).** The pre-pass
   should now *place* same-length hets into its confident-genotype set (previously
   impossible), measured with the P1.5 sequence-concordance instrument
   (`tmp/ssr_p15/seq_concordance.py`, `statSTR --use-length` off).

**A caveat this doc states honestly.** The P1.5 report attributes much of the 46 %
downstream het-undercall tail to the tomato cohort's extreme apparent inbreeding
(`F_IS ≈ 0.821`) tilting the *genotyping-time* prior toward homozygosity — a
mechanism in the genotyping EM / outer `F` loop (roadmap E1/E2), **not** the pre-pass
gate. D1 removes the *pre-pass ε/θ contamination* and lets the pre-pass *see*
same-length hets; the portion of the downstream tail driven by `F_IS` is a separate,
out-of-scope lever. So D1's committed target is **ε de-inflation + the confident-set
recovery**, with any downstream het-concordance rise a welcome secondary effect,
not the pass/fail line. Overclaiming a full 46 %→0 recovery from D1 alone would be
wrong, and the plan's exit criteria are written accordingly (§ plan).

**No regression on the common loci** is a hard requirement: the length-genotype
concordance (96.5 %) and the SNP caller's tests must not move adversely.

All numeric knobs — `het_admission_cost`, the seed `ε`/level/`λ`, the min-depth skip,
the recurrence `k` — are exposed dev defaults, pinned on the simulator in **F2**
(Mark-2 §4.3 "penalty + separation = F2 calibration"). D1 ships working dev values,
not final ones.

---

## 7. The second, coupled change — sequence-aware read attribution in the pre-pass

Resolving a same-length het is only half the fix; the reads must then be **split
between the two same-length alleles by composition**, or `ε` stays contaminated. The
pre-pass's read→allele step (`accumulate_locus`) currently uses length-only
attribution (`nearest_parent`), which ties on two same-length peaks and cannot split
them. D1 switches that step to the **already-existing** sequence-aware attribution
(`nearest_called_by_sequence`, [`attribution.rs`](../../../src/ssr/cohort/attribution.rs)),
so an `interrupted-18` read attributes to the `interrupted-18` allele (faithful, no
mismatch) rather than to `pure-18` (all its interruption bases counted as
substitutions). *This* is the mechanical step that de-inflates `ε`. Both changes live
in D1's named touch points (`rung_ladder.rs` + `prepass.rs`); nothing downstream
moves.

---

## 8. Out of scope — the pieces D1 must not touch (a change here is a scope-creep smell)

- **The genotyping EM** ([`em.rs`](../../../src/ssr/cohort/em.rs)) and its per-locus
  θ refit — D1 only changes the *pre-pass* seed of the frozen params, not the EM.
- **The base measure `G₀`** and its decay fit — untouched; D1 does not change the
  genotype prior.
- **The candidate set** ([`candidate_set.rs`](../../../src/ssr/cohort/candidate_set.rs))
  — the gate's hypotheses come from the *sample's own observed sequences*, not from
  `assemble_candidates`; the cohort candidate machinery is unchanged.
- **The burn-in loop** (roadmap D2) — D1 is the gate; the co-evolution loop that
  re-runs it with fitted params is D2 (§4).
- **The read likelihood `Qᵣ`** — reused verbatim; D1 adds no new scoring math.
- **Stage 1** (the `.ssr.psp` producer) — unaffected; re-run only `ssr-call`.

If implementing any step tugs on one of these, **stop and ask** — it means the design
missed something.

---

## 9. Summary — the one thing to hold onto

*The confident-genotype gate stops being a length histogram and becomes a likelihood
comparison: score the reads under the best one allele vs the best two, admit the
second only when it beats a purity-tuned BIC penalty, and split same-length reads by
composition. Same-length hets — invisible to a length histogram, unlocked by Phase 1
— finally enter the confident set correctly, and `ε`/`θ` on interruption-rich cohorts
stop being contaminated by their minority reads.*
