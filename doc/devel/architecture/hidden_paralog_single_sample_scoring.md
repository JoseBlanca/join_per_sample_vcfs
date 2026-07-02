# Making the hidden-paralog score work from one individual up

**Status:** design, 2026-07-02, branch `tomato2-paralog-filter`. Reformulates
where the per-locus paralog score gets its **heterozygosity expectation** and its
**inbreeding term**, so the same statistic degrades gracefully from a single
sample to a large cohort — no cohort-wide allele-frequency quantity, no
minimum-sample gate.

**Supersedes** Premise 3 of
[hidden_paralog_locus_statistic.md](hidden_paralog_locus_statistic.md) (the
`F = 1 − Hobs/Hexp` derivation) and removes the `min_samples` gate. Everything
else in that document (the H1-vs-H2 likelihood ratio, the coverage model, the
empirical-Bayes prior/FDR) stands.

---

## 0. Why (plain English, no formulas yet)

The hidden-paralog filter drops SNP calls that are better explained as **two
reference-collapsed gene copies** piling their reads onto one position than as
**one real variant**. Two footprints give it away, and the second only in the
individuals that carry the extra copy:

- **excess read depth** — two copies' worth of reads land where the reference has
  one;
- **a heterozygous-looking allele balance** — the fixed difference between the two
  collapsed copies reads out as a "heterozygote."

The strongest evidence is the two **together**: *the heterozygous samples are the
over-covered ones.*

As first built, the filter needed a **cohort**. It computed a per-sample
inbreeding coefficient from a cohort-wide *expected heterozygosity* (built from
every locus's allele frequencies), and it refused to score any locus with fewer
than 20 usable samples. **Neither is fundamental.** This reformulation removes
both, so the identical statistic runs on **one individual** — with only the
coverage footprint to go on, and only when the depth is high enough to see it —
and scales smoothly up to a large cohort at full power, with no discontinuity in
between. The single individual gives up power, not correctness.

---

## 1. Glossary

- **Heterozygous (het) / homozygous** — a diploid site with two different alleles
  (het) or two identical ones (hom). "Hom-alt" = homozygous for the
  non-reference allele.
- **Relative copy number** — a window's read depth expressed in copies
  (`1.0 = one copy`), after dividing out the sample's own sequencing depth and GC
  bias. `≈ 2` is the collapsed-paralog tell.
- **σ₀** — how much a *single-copy* window's relative depth scatters around 1.0
  for a given sample; fit from that sample's own coverage histogram. Small at high
  depth, large at low depth. It is what tells the score whether a "2×" is real or
  noise.
- **SFS prior** — *site-frequency-spectrum* prior: the distribution of allele
  frequencies expected under neutral mutation–drift equilibrium
  (`∝ 1/(p(1−p))`). A prior over how common a variant is, used *without* looking
  at cohort data.
- **F (inbreeding coefficient)** — how much *less* heterozygous a population is
  than random mating would predict at the same allele frequencies. Near 0 for an
  outbred population, near 1 for a selfer (tomato). A genotype-level tilt, one
  number for the cohort.
- **H1 / H2** — the two hypotheses the score weighs at a locus: **H1** = a real,
  single-copy variant; **H2** = a hidden (reference-collapsed) paralog.
- **LR (likelihood ratio)** — `log P(data | H2) − log P(data | H1)`. Positive
  favours a hidden paralog.
- **π** — the estimated fraction of *loci* (genome-wide) that are hidden paralogs;
  one number per run, for the empirical-Bayes posterior and the FDR cut.

---

## 2. What is per-individual, what is cohort, what is across-loci

The reformulation rests on one observation: **the evidence separates cleanly by
what it needs.**

- **Per-individual — works at n=1, the load-bearing signal.** The single-copy
  coverage model (each sample's own depth histogram → expected `1.0` and its
  spread `σ₀`), and, per locus, that sample's window depth and allele counts.
  Nothing cohort-wide.
- **Per-locus, *inferred* not supplied — the heterozygosity expectation.** "How
  many hets do we expect *here*" depends on the allele frequency at this locus.
  The score does **not** need that supplied: **H1 marginalises over the unknown
  frequency under the SFS prior** — it integrates over all frequencies, per locus,
  and lets the data pick. No per-sample het baseline, no cohort allele
  frequencies. (§4.)
- **Cohort-level, weak — the inbreeding coefficient F.** A genotype-level tilt on
  the prior, one number for the cohort, and a *secondary* knob (measured on
  tomato2: swapping a wrong `F` for the exact one moved the flagged set < 10 %,
  π unchanged). It is **not identifiable from one individual**. So we take it from
  the caller's single cohort inbreeding coefficient — the same one the SNP caller
  already uses — not a per-individual quantity. (§5.)
- **Across *loci*, not across *individuals* — π.** The fraction of loci that are
  paralogs, for the empirical-Bayes posterior and the FDR cut, estimated from the
  distribution of LRs over loci. It works with one individual and many loci. (No
  change.)

So the score is **per-individual + per-locus**, with one **across-loci**
calibration (π) on top. The only thing that was cohort-*wide* (the expected-het
`Hexp` from allele frequencies) is removed.

---

## 3. Change 1 — drop the `min_samples` gate; the statistic self-gates

The old gate skipped any locus scored on fewer than 20 usable samples. It is
unnecessary, because **the likelihood ratio already collapses to "keep" when the
data can't distinguish the two hypotheses.**

When the coverage can't tell 1× from 2× — `σ₀` large relative to the gap — and the
reads are few (so the allele-balance binomial is weak — see §4), every branch of
H1 and H2 assigns the data nearly the same likelihood. `log P(H2) − log P(H1) → 0`
(the ratio → 1). A locus at the null can't cross a positive FDR cut, so it is
**kept**, with no gate. And `σ₀` is *fit from the same window-depth distribution*
the score uses, so it already encodes the per-window uncertainty at that sample's
depth — the self-gating is correctly calibrated, not a hand-tuned threshold. If
anything, the marginal (Occam) penalty on H2's extra parameters tips an
uninformative locus slightly toward *keep* — the safe direction.

**Consequence to accept:** with the gate gone, many uninformative loci fold LR≈0
into the histogram the empirical-Bayes uses. This mildly *dilutes* π (more loci
look null). That is the **conservative** direction — a slightly stricter cut — and
those loci wouldn't be flagged anyway. Watch it in validation (§7); do not add a
hard gate back.

---

## 4. Change 2 — the het expectation comes from H1's SFS marginalisation, joint with coverage

"Is this het *unexpected*?" needs an expected heterozygosity at the locus, which
needs the locus's allele frequency. The score gets it the right way — by
**marginalising H1 over the unknown frequency** under the SFS prior
(`∫ P(data | p) · SFS(p) dp`), exactly as it already does. No per-sample het
baseline is supplied; the expectation is *inferred* per locus.

**Why this is immune to the paralog hiding itself (the load-bearing detail).** A
collapsed paralog makes every carrier look het, so the locus's *apparent* allele
frequency comes out ≈ 0.5 — "a common variant." If the het expectation alone drove
the verdict, the paralog would explain away its own hets and escape. It cannot,
because **H1's coverage term is `Normal(1, σ₀)` regardless of frequency**: a
paralog's ≈ 2× depth makes H1's total likelihood low *no matter what frequency it
invokes*, so H2 wins. Two things fall out, and they are exactly the operator's
two-expectation intuition:

- the **SFS marginalisation** lets a *genuine* common variant (real high `p`,
  normal coverage) be fully explained by H1 → **not** flagged;
- the **coverage** is what makes a *paralog* (hets riding on 2×) unexplainable by
  H1 → flagged.

An unexpected het only becomes evidence when it rides on unexpected coverage —
which is why the score is a *joint* coverage+allele likelihood, not a het test and
a coverage test bolted together.

**Do not feed H1 the caller's per-locus frequency point estimate.** The caller's
`p̂` is the thing the paralog inflates; although the coverage term would still
catch it, a contaminated `p̂` needlessly weakens H1. We reuse the caller's
*machinery* (an SFS prior + a cohort inbreeding coefficient), not its per-locus
number.

**Optional refinement — calibrate the SFS prior to the cohort.** The prior is
today the theoretical `∝ 1/(p(1−p))` with a `1/2N` floor. We *may* replace it with
the cohort's *observed* frequency spectrum (from the caller's allele frequencies)
when a cohort exists — a better match for a selfing or bottlenecked population. It
is cohort-dependent (so it **defaults to the theoretical prior at n=1**) and, like
`F`, probably a secondary knob under coverage. Wire it as an *optional*
calibration, measure whether it moves the drop set, keep the theoretical prior as
the default. **Deferred; not part of the core change.**

---

## 5. Change 3 — F is the caller's cohort inbreeding coefficient, not a per-individual `Hexp`-derived value

`F` is the genotype-level inbreeding tilt on H1's prior. The old derivation formed
it per sample from `Hobs_s / Hexp`, where `Hexp` was a cohort-wide expected
heterozygosity accumulated from every locus's allele frequencies during the main
pass. **Replace it with the single cohort inbreeding coefficient the SNP caller
already uses** (`--inbreeding-coefficient` / its default), applied to every
sample.

Why this is the right call, not a shortcut:

- **A correct per-individual `F` is not achievable AF-free.** `F` is *by
  definition* a deviation from a random-mating expectation, which needs a
  population reference. From one individual there is no such reference.
- **The obvious AF-free per-individual proxy is a measured trap.** Using the het
  rate *over variant loci* (`n_het / n_variant`) as the baseline **inverts**: its
  denominator is dominated by the hom-alt count, i.e. by distance from the
  reference, not by inbreeding, so a divergent-but-outbred sample looks inbred and
  its *real* hets get over-flagged. Premise 3 measured exactly this on tomato2
  (Spearman **−0.28** with the true cohort `F`; most-outbred samples came out at
  `F ≈ 0.83`). We deliberately avoid it.
- **`F` is a weak knob anyway** (Premise 3: swapping it moved the flagged set
  < 10 %, π unchanged), so a single cohort value costs little accuracy while
  buying consistency with the caller's genotype model and a clean n=1 default.
- **It removes a whole global quantity.** `Hexp`, its main-pass accumulator, and
  the per-sample `obs_het` all disappear; `F` becomes a single number known up
  front, so the two-pass flow gets simpler (§8).

At n=1 there is nothing to estimate: `F` is whatever the operator supplies (or the
default), and coverage carries the verdict.

---

## 6. Graceful degradation, n=1 → many

The point is a **smooth transition**, not a cutoff:

- **n = 1:** only the coverage footprint is available (the single sample is het or
  not; there is no cross-individual het pattern). The filter detects only
  **strong** signals — a clear ≈ 2× at good depth — and its empirical-Bayes
  calibration is weak (the LR distribution is less bimodal, so π is noisier and
  the fixed-prior fallback fires more often). It *works*, but it is
  coverage-dominant with soft calibration.
- **n growing:** the het-rides-on-coverage coupling sharpens (more carriers to
  correlate), the LR distribution separates, π calibrates. Full power at cohort
  scale — the R1/T1 result.

Nothing switches on at a threshold; power rises continuously with n and with depth.

**Policy left open (not a bug):** at very low n we *may* prefer to **annotate**
the paralog score on survivors rather than hard-drop, since the FDR guarantee is
soft there. That is an operator choice, decided in the plan, not a property of the
statistic.

---

## 7. Assumptions to state, and validation without a truth set

**Assumption (make it explicit in the code docs):** the paralog fraction is
**small**. The inbreeding coefficient's absence and the empirical-Bayes π both
tolerate a little self-reference — the calibration includes a few of the very loci
it is estimating — only because realistic populations carry few collapsed-paralog
sites. State this where π and the prior are documented; if a future dataset
violates it, the escape hatch is to re-estimate on the survivors.

**Validation has no truth set** (tomato2 has no benchmark). Judge the change on
**profile coherence**, not accuracy:

- the drop set must still carry the **coverage-excess + het-excess** signature
  (T1's mean-DP and het-fraction split);
- **freebayes** must still emit the dropped loci (the cross-caller
  "leaves-them-in" story);
- **π** must stay sane and the EM converge on the cohort case;
- diff the **old vs new** drop sets to *understand* the shift (the new `F` and the
  gone gate will move it — expected);
- a **high-coverage small-n** check (e.g. one deep tomato sample) to confirm the
  coverage-only path fires on clear paralogs and stays quiet elsewhere.

"Less aggressive" is only bad if it starts *keeping* loci with the clear paralog
profile; the profile is the yardstick.

---

## 8. What is removed / simplified

- `HexpAccumulator` (main-pass accumulation) — **gone**.
- Per-sample `obs_het` in the pre-pass — **gone** (only fed the old `F`).
- The `hexp` threaded through calibrate + write pass, and the `SinkOutput::Spill`
  carrying it — **gone**; `F` is a single up-front number.
- The `min_samples` gate in the calibrate pass — **gone**.
- The pre-pass shrinks to *just* the per-sample coverage model.

The paralog verdict becomes: per-sample coverage model (up front) + H1/H2 LR with
an SFS-prior het expectation and a cohort `F` + across-loci π. Per-individual and
per-locus, one global calibration, no cohort allele frequencies.
