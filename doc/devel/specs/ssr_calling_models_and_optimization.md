# SSR calling — statistical models, the options we can turn on, and the plan to find the best combination

*Status: reference + plan, 2026-07-08 (rev. 2026-07-09 — the `MARG-SFS` genotype-prior
probe folded in as §3.2b / §4 finding 4, null), branch `ssr-bic-emission` (all emission
models unified behind `PVC_SSR_EMIT_MODEL`). This document has two jobs. First, it
is a **single map** of every statistical model and tunable option now in the SSR
caller, so we can reason about them together instead of one report at a time.
Second, it lays out the **plan to get a gold-standard truth set** and use it to
choose the optimal combination — because our current yardsticks cannot settle the
one question that matters most. Companion detail lives in the per-model specs cited
throughout; where any of them disagree with [`ssr_cohort_mark2.md`](ssr_cohort_mark2.md)
on intent, Mark-2 wins.*

---

## 0. The vocabulary this document leans on

Defined once, in plain terms, because the rest of the document reuses them.

- **Microsatellite (SSR / STR)** — a short DNA motif repeated in tandem (e.g.
  `AT AT AT …`). "Calling" it means measuring, per plant, how many copies each
  chromosome carries — its **length** in repeat units.
- **Stutter** — the dominant error: the PCR replication slips and produces a read a
  whole repeat-unit or two off the true length. It is **systematic per sequencing
  chemistry** (the same library prep slips the same way in every plant), not random
  per plant. It is what makes SSRs hard.
- **Segregation** — the signal that tells a real allele from stutter: a real allele
  is *read-dominant in the plants that carry it and absent in the rest*; stutter is
  a *consistent low shoulder present in every plant*. It is a **cohort-level**
  signal — no single plant's reads reveal it, but the cohort's do.
- **Read likelihood `Lr(read | allele)`** — the probability an observed tract arose
  from a candidate allele, via stutter (a length change) plus base error (a
  composition change). The caller's HipSTR-style model
  ([`likelihood.rs`](../../../src/ssr/cohort/likelihood.rs)); every model below
  reuses it rather than inventing a new one.
- **Genotype likelihood `Lg(s, g)`** — the likelihood of plant `s`'s reads under
  genotype `g`, `P(reads_s | G_s = g)` — the standard "genotype likelihood" (the VCF
  `GL` field). Built from `Lr`, once per locus; the code stores its natural log as
  the array `data_ll[s][g]`. This is the single read-level quantity every emission
  model consumes — so a comparison between models changes only the *decision*, never
  the reads.
- **Inbreeding coefficient `F`** — how much more homozygous a population is than
  random mating predicts. Tomato self-pollinates, so it is high here (`F_IS ≈ 0.82`)
  and almost every true genotype is homozygous. Threaded per plant.
- **Marginal likelihood** — the probability of the data under a model after
  *averaging over* the model's unknowns (here: the genotypes and the allele
  frequency) rather than fixing them at their best-fit values. Averaging charges a
  model for its flexibility — an automatic Occam penalty.
- **BIC (Bayesian Information Criterion)** — a model-selection rule that rewards fit
  (log-likelihood) and charges a fixed cost per extra free parameter (`k·ln n`). The
  extra allele must earn its parameter.
- **Site-frequency-spectrum (SFS) prior** — the population-genetics prior on whether
  a site carries a variant and, if so, in how many of the cohort's chromosomes. The
  neutral (Ewens) form is governed by a diversity parameter `θ`; a small `θ` puts
  most of the prior mass on the site being monomorphic and the rest on rare
  (low-carrier-count) alleles. (Its role in emission is §3.3; its tested-and-rejected
  role as a genotype-prior base measure is §3.2.)
- **Silver standard** — a *confidently-labelled subset* of loci we treat as truth
  in the absence of a gold standard, built from the reads (§5). **Not** a truth set.
- **Concordance** — the fraction of individual genotype calls (one plant at one
  locus) on which two callers *agree*, counted over the calls both of them made. It
  measures **agreement between the two callers, not accuracy**: they can agree and
  both be wrong, so agreeing with a reference caller that is itself imperfect — here
  HipSTR, which is biased on this cohort (§5) — is not evidence of being correct.

> **Notation note.** Earlier SSR specs and the source code write the read likelihood
> as **`Qᵣ`** and the genotype log-likelihood array as **`data_ll`**. This document
> uses **`Lr`** (read likelihood) and **`Lg`** (genotype likelihood) instead: `Q`
> collides with the Phred quality score `Q`, and `data_ll` says nothing about *which*
> data. Clarity outranks matching the old symbols. When the SSR work next touches
> those files, the code and the remaining specs should adopt `Lr` / `Lg` too, so we
> end up with one name per concept again.

---

## 1. Where the choices live — the pipeline, and a name for each choice

`ssr-call` (Stage 2) turns per-plant `.ssr.psp` evidence into a cohort VCF in five
steps. The rest is the settled base model (§2); the choices cluster at two steps.
This section names every choice so the rest of the document can refer to each by a
short **handle** (in bold below).

```
  pre-pass   →  candidate   →  genotyping   →  emission           →  VCF
 (freeze        assembly +     (per-locus       (per-locus
  chemistry)     admission      EM)              emit / drop)
     │                            │                  │
 confident-gt                genotype prior:    emission model:
 gate (§3.1)                 plug-in | MARG     heuristic | BIC | freebayes
                             (§3.2)             + clean-hom null (§3.4)
                                                + KEEP no-call   (§3.5)
```

**The choices, by name** (each is the default's alternative unless stated; detail in
the cited section):

- **confident-genotype gate** — in the chemistry pre-pass, which plants get to teach
  the caller its error/stutter parameters (§3.1).
- **genotype prior** — how the genotyping EM builds each plant's prior: **plug-in**
  (default) or **marginalized** (`MARG`) — plus the marginalized prior's **SFS
  base-measure** variant (`PVC_SSR_MARG_SFS`), tested and rejected (§3.2).
- **emission model** — how the caller decides whether to write a locus at all:
  **heuristic** (default), **BIC**, or **freebayes** (§3.3).
- **emission null** — what the "monomorphic" hypothesis uses for its stutter:
  **genotype-driven** (default) or **clean-hom** (§3.4).
- **per-sample no-call** (`KEEP`) — under BIC/freebayes, whether to still no-call
  shaky per-plant genotypes for quality (§3.5).

Two facts the diagram cannot show but the reader needs:

- **The EM always does the genotyping; BIC and freebayes are not alternative
  genotypers.** They run *after* the EM and change only the *emit / drop* decision —
  the genotypes written to the row are the EM's either way. (In freebayes-the-tool
  the Bayesian model replaces an EM outright; here we borrow only its emission test.)
- **BIC and freebayes differ in how they treat the allele frequency** — this is the
  "freebayes-vs-EM" relationship. **BIC reuses the EM's fitted frequency**;
  **freebayes marginalizes the frequency itself** against the SFS prior (§0), so its
  emit decision does not lean on the EM's estimate. Both still emit the EM's
  genotypes, which is why the **genotype prior** and the **emission model** are
  orthogonal knobs whose interaction is the whole story (§4). It also makes the
  freebayes path **statistically incoherent, by design**: the site QUAL and the
  reported genotypes then come from *two different* frequency priors (the SFS for
  the emit call, the EM's HWE/Dirichlet for the genotypes). We emit only when the two
  models agree (§3.3), so it holds where they do, but it is not one unified posterior.
  The cheap route to coherence — put the SFS prior in the genotyper too — is the
  probe in §8; the fully-coherent alternative is the memory-costly freebayes
  genotyper we set aside (§7).

---

## 2. The base model (always on, not a toggle)

These are settled and shared by every option; the detail is in
[`ssr_cohort_mark2.md`](ssr_cohort_mark2.md).

- **The read model `Lr`** — stutter (length) + base error (composition), the
  bake-off winner ("Model A"). Per **sample-group (chemistry)** error `ε`, stutter
  shape, and level; **per-locus** refinement of the shape and level
  ([`refine_theta_locus`](../../../src/ssr/cohort/stutter.rs)), each shrunk toward
  its chemistry prior so a locus departs from the chemistry baseline only when its
  own data justifies it.
- **Candidate assembly + admission** — build the per-locus allele set and gate the
  locus on depth (`lowDepth`), periodicity (`notPeriodic`), and allele count
  (`tooManyAlleles`).
- **The chemistry pre-pass** — freeze `ε`, the stutter shape/level, the base measure
  `G₀`, and per-individual `F` from the cohort's confident genotypes, before
  genotyping runs.

---

## 3. The tunable options (the catalog)

Every option is **off by default and byte-identical when off** — the default VCF is
the settled heuristic caller. Turning one on is a runtime environment variable.

| # | option | env var | default | axis |
|---|---|---|---|---|
| 3.1 | confident-genotype gate | *(pre-pass, design)* | — | chemistry purity |
| 3.2 | genotype prior | `PVC_SSR_MARGINALIZED_PRIOR` | plug-in | genotyping |
| 3.2b | genotype-prior base measure (SFS) | `PVC_SSR_MARG_SFS` | mode-centred `G₀` | genotyping *(tested → null, §4)* |
| 3.3 | emission model | `PVC_SSR_EMIT_MODEL` | `heuristic` | emission |
| 3.4 | emission null | `PVC_SSR_NULL_FROM_HOMS` | genotype-driven | emission |
| 3.5 | per-sample no-call under a model | `PVC_SSR_KEEP_FP_CONTROL` | off | per-sample QC |
| 3.6 | cohort-recurrence rescue *(superseded)* | `PVC_SSR_COHORT_FP_CONTROL` | off | per-sample QC |

### 3.1 Confident-genotype gate (chemistry pre-pass)

*What it does.* When the pre-pass decides which plants are "confident" enough to
teach it the chemistry, it replaces a length-histogram heuristic with a **likelihood
test**: score each plant's reads under the best *one* allele versus the best *two*,
and admit a second allele (call it a heterozygote) only when the fit gain beats a
purity-tuned BIC penalty. This lets it see a **same-length heterozygote** (two
alleles of equal length but different composition — an interrupted repeat vs a pure
one) that a length histogram structurally cannot, and stops those minority reads
from inflating the error rate `ε`. Detail:
[`ssr_bic_confident_genotype.md`](ssr_bic_confident_genotype.md), validated
[report](../reports/ssr_bic_confident_genotype_validation_2026-07-07.md).

*Why it is in this map.* It is a **different** BIC test from the emission one (§3.3):
this one runs in the pre-pass and cleans the *chemistry estimate*; the emission BIC
runs later and decides *which loci to emit*. They share only the idea "one allele vs
two, judged by likelihood." Keeping the distinction explicit avoids confusion.

### 3.2 Genotype prior — plug-in vs marginalized (`PVC_SSR_MARGINALIZED_PRIOR`)

*What it does.* The genotyping loop (expectation-maximisation) needs a prior on each
plant's genotype, built from the cohort allele frequency. Two forms:

- **Plug-in (default)** — estimate one allele frequency, plug it into a Wright
  genotype prior (with `F`), call each plant against it.
- **Marginalized + leave-one-out (opt-in)** — average over the frequency's
  uncertainty (a Dirichlet-multinomial marginal), and give each plant a prior
  computed from *all the other* plants (leave-one-out), starting from a flat first
  step. This is the SNP path's empirical-Bayes fix, ported to SSR.

*Why it matters.* At ~3 reads per plant the plug-in prior falls into a **homozygous-
reference trap**: a real-but-rare allele keeps the frequency low, a low frequency
biases every shallow plant toward hom-ref, and calling them hom-ref keeps the
frequency low. The marginalized prior escapes the trap and **re-genotypes the
lone-carrier recurrent heterozygotes as variant** — which turns out to be exactly
the cells the whole recall/precision tension is about (§4).

*A third variant — the SFS base measure (`PVC_SSR_MARG_SFS`), tested and rejected.*
The marginalized prior integrates each plant's genotype over the frequency `p` drawn
from a Dirichlet whose **base measure** is the mode-centred `G₀` (§4.3). One can ask
whether swapping that base for the **Ewens `θ/k` SFS** — the same neutral spectrum that
gives the freebayes *emission* its edge (§3.3) — produces better *genotypes*. It does
so cheaply: the faithful finite-`k` SFS prior on `p` *is* a symmetric `Dir(θ/k)`, so the
graft is the existing marginalized machinery with a flat sparse `θ/k` base in place of
`G₀` — no joint/cohort-scale structure (reusing the freebayes count-vector marginal
literally per plant would be the memory-costly full genotyper we rejected, so we did not).
The **null verdict** is in §4 (finding 4): at matched emission and matched false-positive
rate the SFS base is indistinguishable from `G₀`. Kept as a documented dead end, not a
recommended toggle. Detail:
[`ssr_marg_sfs_genotype_prior_2026-07-09.md`](../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md).

### 3.3 Emission model — `PVC_SSR_EMIT_MODEL = heuristic | bic | freebayes`

*What it does.* Decides, per locus, whether to write it to the VCF as variable or
drop it as monomorphic. All three consume the **same** genotype likelihoods `Lg`.

- **`heuristic` (default)** — emit iff some plant's called genotype carries a
  whole-unit length variant, after an allele-balance **no-call** step
  (`apply_fp_control`) demotes imbalanced (stutter-shaped) heterozygotes. A stack of
  per-plant rules; leaks false positives and has no clean knob.
- **`bic`** — model selection: is the cohort's reads better explained by
  *monomorphic* (one allele + stutter) or *polymorphic* (an alt segregates)? Emit
  iff `2·(ln_marginal − ln_monomorphic) > n_alt·ln(N) + 2·margin`. Knob
  `PVC_SSR_BIC_MARGIN` (raise it → fewer, surer emissions). Detail:
  [`ssr_bic_emission_2026-07-08.md`](../reports/ssr_bic_emission_2026-07-08.md).
- **`freebayes`** — the same monomorphic-vs-polymorphic question, but Bayesian: the
  posterior `P(monomorphic | data)`, integrating the allele frequency against the
  Ewens SFS prior (§0) and `F`-aware. `QUAL = −10·log10 P(monomorphic)`; emit iff
  `QUAL ≥ PVC_SSR_FREEBAYES_MIN_QUAL` (the precision/recall knob), with
  `PVC_SSR_FREEBAYES_THETA` the SFS `θ` (default 0.01). The small-`θ` prior is what
  makes the test conservative — a locus must show strong cohort segregation to
  overturn the monomorphic default — so the **lone-carrier** false positives (an
  apparent allele in one or two shallow plants; §4) are rejected as rare,
  weakly-evidenced alleles, and more sharply than BIC's flat per-allele penalty
  (freebayes' edge). Detail:
  [`ssr_freebayes_marginal_emission.md`](ssr_freebayes_marginal_emission.md).

The `bic` and `freebayes` paths **skip** `apply_fp_control` by default — the site
model owns the emit decision — which is what §3.5 exists to reconsider.

### 3.4 Emission null — clean-hom (`PVC_SSR_NULL_FROM_HOMS`)

*What it does.* A false positive is a *monomorphic null the model can't explain its
own stutter with*: if the null's stutter rate is too low, the shoulder looks like an
allele. By default the null's stutter comes from the genotype-driven estimate, which
is deflated at exactly the false-positive loci (the shoulder reads got attributed to
the phantom allele). This option instead learns the null's per-locus stutter **from
the clearly-homozygous plants** (`attribute_clean_homs`) — plants whose off-length
reads are unambiguously stutter — so the monomorphic model can fully **absorb the
shoulder** and stops ceding it to a phantom allele. Emission-only; genotypes are
untouched. Knobs `PVC_SSR_NULL_HOM_MINDEPTH` (4), `PVC_SSR_NULL_HOM_FRAC` (0.75).

### 3.5 Per-sample no-call under a model (`PVC_SSR_KEEP_FP_CONTROL`)

*What it does.* Under `bic`/`freebayes`, the **site** is decided by the model, but
the emitted row's per-plant genotypes can still get the allele-balance **no-call**
(`apply_fp_control`) for genotype quality — the model decides the locus, the no-call
cleans the cells. This is the switch that recovers HipSTR concordance (§4), and it is
**not** independent of §3.2 the way we first assumed.

### 3.6 Cohort-recurrence rescue — superseded (`PVC_SSR_COHORT_FP_CONTROL`)

*What it was.* A heuristic that spared an imbalanced heterozygote the no-call when
its minority allele was read-dominant in ≥2 other plants
([`ssr_cohort_recurrence_aware_gates.md`](ssr_cohort_recurrence_aware_gates.md)).
Kept for reference, but **superseded** by the emission models — the model-selection
test is the principled version of the same "does this recur as a real allele?"
question. Not part of the recommended combinations.

---

## 4. What we have learned so far — and the tension our yardsticks expose

Measured on the silver confident core (561 confident-true loci, 8850
confident-monomorphic) and by HipSTR concordance; full numbers in
[`ssr_emission_model_comparison_2026-07-08.md`](../reports/ssr_emission_model_comparison_2026-07-08.md).

1. **`freebayes` ≥ `bic`, modestly.** Same-or-slightly-more recall at every matched
   false-positive rate, and it reaches the ~87% recall ceiling at 0.28% FP where BIC
   needs 0.35%. It is also the cleaner instrument — one real per-locus QUAL, no
   arbitrary margin unit. The toy-prototype's "2–3× better than BIC" shrank to a
   small-but-consistent edge once BIC ran on the real chemistry.
2. **The marginalized prior is what unlocks high recall** — plain BIC and the
   heuristic both plateau at ~81% because they emit the plug-in EM's hom-ref-collapsed
   genotypes; MARG re-genotypes the lone-carrier hets and lifts both models to ~87%.
3. **The central finding — recall and concordance are the *same cells*, and they do
   not compose.** The extra silver recall and the HipSTR concordance are two views of
   one set: the **lone-carrier recurrent heterozygotes**. Emit their genotype and you
   gain silver recall but lose HipSTR agreement (HipSTR, `F`-blind on a selfer,
   no-calls or hom-calls them); no-call them and you gain agreement but lose the
   recall. Concretely: `freebayes + MARG + KEEP` collapses recall to ~54% (the
   no-call removes exactly what MARG recovered), while `freebayes + MARG, no KEEP`
   holds ~87% recall at ~88% concordance. **You cannot bank both at 3 reads/plant.**
4. **The SFS frequency prior does not improve *genotyping* over the marginalized DM
   prior (§3.2b) — a null result.** Swapping the marginalized prior's base measure from
   the mode-centred `G₀` to the flat Ewens `θ/k` SFS leaves genotypes unchanged where it
   matters: on the freebayes frontier (QUAL swept, so a recall gap at matched FP *is* a
   genotyping gap) DM-marg and SFS-marg are indistinguishable at every `MIN_QUAL ≥ 10`
   (recall within 0.7 pt, FP within 0.03%); both beat plug-in, neither beats the other.
   The one difference — the sparser SFS base is ~3× cleaner on FP at the *unfiltered*
   point (it suppresses the lowest-confidence lone-carrier false positives) — is
   **redundant with the `MIN_QUAL` knob**, which prunes the same loci. This confirms
   directly what finding 1 implied: the SFS prior's value lives in **marginalizing the
   frequency in the emit decision**, not in a better genotype frequency model — so the
   memory-costly full freebayes *genotyper* is **definitively not worth building** (the
   cheapest graft of its frequency model is inert). Detail:
   [`ssr_marg_sfs_genotype_prior_2026-07-09.md`](../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md).

The current recommendation (two regimes, from the comparison report):

- **Precision/QC default:** `freebayes`, MARG **off**, KEEP **on**, `MIN_QUAL ≈ 20` —
  matches the heuristic on recall/precision/concordance (96.5%) but from one
  principled test with a clean knob, emitting ~46% more genotyped loci.
- **Recall/discovery regime (opt-in):** `freebayes`, MARG **on**, KEEP **off**,
  `MIN_QUAL ≈ 20` — ~84% recall / 0.19% FP / ~88% concordance, for maximum real-
  variant discovery at the cost of HipSTR agreement.

**But which regime is *right* depends on a fact we do not know:** are those
lone-carrier hets real variants (→ the discovery regime is undercalling) or
systematic stutter (→ the precision regime is correct)? Our yardsticks cannot tell us.

---

## 5. Why our yardsticks cannot finish the job

Both current measures are **proxies**, and both are blind to exactly the contested
cells:

- **The silver standard is read-grounded, so it is partly circular.** It labels a
  locus `true100`/`false100` by whether an allele is read-*dominant* across plants —
  the same segregation signal the callers use. A model that keys on segregation will
  agree with it almost by construction. It is a good *relative* yardstick (it ranks
  the options consistently) but not an *absolute* one, and it cannot adjudicate a
  cell it is itself unsure about.
- **HipSTR concordance is agreement, not accuracy, against a biased reference.**
  HipSTR is `F`-blind and over-calls dinucleotide heterozygotes on a selfer, so
  "disagreeing with HipSTR" is not the same as "wrong" — and the structured
  discordance (ours skews *shorter*) is consistent with HipSTR over-expanding, not
  with us erring. Using it as truth would systematically punish the correct behaviour.

The one question that decides the default — **are the lone-carrier recurrent hets
real?** — is exactly the question on which both proxies are silent or biased. That is
the definition of needing a gold standard.

---

## 6. The gold standard — what, how, and what it resolves

*What we need.* An **orthogonal, higher-accuracy genotype** for a subset of loci and
plants — one produced by a method that does not share the short-read caller's failure
mode — to serve as truth. It does **not** need to be genome-wide; a few hundred loci
chosen to span the motif lengths, the depth range, and — crucially — the **contested
lone-carrier-het cells** is enough to anchor everything.

*The options, best first.*

1. **Long reads (PacBio HiFi) spanning the whole tract.** A single read that covers
   the repeat plus both flanks measures its length directly, with no stutter-vs-allele
   ambiguity — as close to truth as short-read data can be checked against, and the
   modern STR gold standard. Best if any of these accessions (or close relatives) can
   be HiFi-sequenced, even at modest depth on a subset.
2. **Capillary electrophoresis (fragment sizing).** The classical microsatellite gold
   standard: PCR a locus, size the product. Laborious per locus, but decisive, and a
   few dozen loci genotyped this way across the contested set would settle the tension.
3. **A deep re-sequence of a few plants.** Raising a handful of plants to high depth
   turns their now-ambiguous 3-read cells into confident calls, which act as truth for
   the same cells at low depth. Cheapest to obtain from existing pipelines, but it only
   resolves cells those particular plants carry.

*What it resolves.* Three things our proxies cannot:

- **The recall/concordance tie (§4).** On the gold subset we can finally ask: at the
  contested lone-carrier-het cells, is the alt real? That single answer picks the
  default regime — if real, the discovery regime is right and HipSTR concordance is a
  red herring; if stutter, the precision regime is right.
- **The silver standard's circularity (§5).** Calibrate the confident core against
  gold — measure how often `true100`/`false100` are actually right — and we can trust
  (or correct) the silver numbers on the loci gold does not cover.
- **Absolute error rates.** Convert every combination's *relative* silver ranking into
  a real false-positive and false-negative rate and a real genotype accuracy.

---

## 7. The optimization — the combination matrix and how gold picks the winner

*The axes.* The options compose into a matrix, but the evidence already prunes it:

| axis | values to test | pruned by |
|---|---|---|
| emission model | `freebayes`, `bic` | `heuristic` is the thing we are replacing; keep as the reference point only |
| genotype prior | MARG on, MARG off | the **SFS base-measure** variant (§3.2b) is silver-null (§4, finding 4) → dropped from the sweep; gold may confirm it once, but do not carry it in the matrix |
| per-sample no-call | KEEP on, KEEP off | **MARG on + KEEP on is dead** (recall → 54%); do not test it |
| emission null | clean-hom on, off | test on both models |
| knob | `MIN_QUAL` / `margin` sweep | continuous — report the frontier, not points |

So the live combinations are, per emission model: `{MARG off × KEEP on}` (the
precision/QC regime), `{MARG on × KEEP off}` (the discovery regime), each `× clean-hom
{on, off}`, swept over the QUAL/margin knob. A dozen frontiers, not hundreds of points.

*The metrics, now with a truth anchor.* For every combination, on the **gold subset**:
true false-positive rate, true false-negative rate (recall), and genotype accuracy at
the cell level. On the **silver core** (the loci gold does not cover): the same recall
/ FP frontier we have now, but *calibrated* by how well silver matched gold where they
overlap. HipSTR concordance drops from a decision axis to a sanity check.

*The decision procedure.*

1. Build the gold subset (§6), deliberately over-sampling the contested cells.
2. Score every live combination against gold → true FP/FN/accuracy; against calibrated
   silver → the extrapolated frontier.
3. Read off the **Pareto frontier** (no combination beats it on all three metrics).
4. The gold answer to "are the contested cells real" selects **which point on the
   frontier is the default** — the precision regime if they are mostly stutter, the
   discovery regime if they are mostly real. Publish the runner-up as the documented
   opt-in for the other use case.
5. Fold the winning combination's toggles into the caller's defaults (they are already
   implemented and byte-identical when off, so this is a config change, not a rewrite).

*Expected outcome, stated honestly.* The gold standard will most likely confirm that
**both regimes are legitimate and the choice is use-case-dependent** — discovery vs
precision — and its real value is telling us *which to make the default* and *what the
absolute error rates are*, not revealing a single combination that dominates. The
depth floor (§8) means it will not hand us a free precision jump.

---

## 8. What remains after gold, and the honest ceiling

- **Depth is the wall.** At median ~3 reads/plant the recall/precision/concordance
  tension is a genuine information limit, not a modelling gap — the gold standard will
  *confirm* the wall on this panel and quantify it, but the real fix for it is deeper
  data, not a better prior. The clearest single follow-up beyond emission is a **better
  per-locus null stutter estimate** (§3.4 is the first step) and, more fundamentally,
  a deeper cohort.
- **The clean-hom null (§3.4) and the pre-pass gate (§3.1)** interact with the
  emission choice and have not been swept in the full matrix yet — §7 folds them in.
- **A coherent, cheap frequency prior — the `MARG-SFS` probe (DONE, null).** The
  freebayes path decides the site under the SFS prior but reports the EM's genotypes
  under a different one (§1) — a deliberate statistical incoherence. We grafted the SFS
  prior into the genotyper cheaply (the marginalized prior's base measure swapped `G₀ →
  θ/k`, keeping the per-locus EM, *not* the memory-costly joint genotyper of §7 —
  `PVC_SSR_MARG_SFS`, §3.2b) to test whether it improves *genotypes*, not just emission.
  **It does not** (§4, finding 4): at matched emission and FP the SFS base is
  indistinguishable from `G₀`. This resolves the incoherence question the pragmatic way —
  the two priors give the same genotypes, so the incoherence is cosmetic — and **closes
  the freebayes-genotyper question for good**: the cheapest graft of its frequency model
  is inert, so the expensive one is not worth building. The remaining statistical
  incoherence (site vs genotype prior) is therefore a documented non-issue, not a
  to-do.
- **Ploidy.** The genotyping is diploid-only today; a polyploid grating is a separate,
  documented follow-up and changes the genotype-configuration priors above.

---

## 9. Summary — the one thing to hold onto

*We now have a principled, tunable stack — a marginalized cohort prior for
genotyping and a marginal-likelihood emission test (freebayes-style, BIC as the
sibling) — that matches the old heuristic and extends its reach. But at 3 reads per
plant the last recall lives in the lone-carrier heterozygotes, and emitting them
trades directly against HipSTR agreement; the same cells, two views. Our read-grounded
silver standard and the `F`-blind HipSTR reference cannot say whether those cells are
real, and that single fact decides the default. The next move is not another model —
it is a gold standard (long reads first) over a few hundred loci, over-sampling the
contested cells, to convert our relative frontier into absolute error rates and to
tell us which regime to ship.*
