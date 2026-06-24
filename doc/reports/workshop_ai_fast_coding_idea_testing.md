# Testing ideas before committing to them: an AI-assisted bake-off

*A worked example for a workshop on AI in programming and science.*

**Session date:** 2026-06-24 · **Project:** a from-scratch multi-sample variant
caller (the SSR/microsatellite genotyping module) · **Tool:** an AI coding agent
driving the editor, compiler, and test runner under human direction.

---

## 1. The point of this report

In science and in engineering we constantly face forks where the "right" choice
is not obvious from first principles. The literature suggests one thing, our
intuition suggests another, and the honest answer is *we don't know until we
measure*. Traditionally we resolve these forks by **argument** — we pick the
option that sounds best, or the one a respected paper used, because actually
*building and comparing* all the candidates is too expensive to justify.

Fast, AI-assisted coding changes that economics. When implementing a full
candidate solution drops from *weeks* to *part of an afternoon*, it becomes
cheaper to **build all the candidates and let the data decide** than to argue
about them. This session is a concrete instance: we had three credible ways to
compute one core quantity, we built all three, raced them on synthetic data with
known answers, and let the measurements — not seniority or citation count —
choose the winner. Along the way the experiment surfaced two problems that
neither the literature nor our intuition had flagged.

---

## 2. The fork in the road

The microsatellite genotyper scores every sequencing read against every candidate
allele with a **read likelihood** `Qᵣ(obs | candidate)` — the probability that an
observed repeat tract arose from a candidate allele through PCR/sequencing
*stutter* (the length slips characteristic of repeats) plus ordinary base error.
This single function drives the whole genotyping step, so its shape decides
accuracy, calibration, and speed.

There were three credible ways to compute it:

- **Model A — "follow the literature."** Mirror **HipSTR**, the field-standard
  microsatellite caller: an explicit stutter probability that separates *in-frame*
  (whole-repeat-unit) slips from *out-of-frame* (partial) slips, each with its own
  geometric distribution.
- **Model B — "fix what we have."** Keep our current approach (sum over whole-unit
  slips, score the leftover by substitution) and just repair a known weakness in
  how its parameters were estimated.
- **Model C — "our own intuition."** A single alignment with **two gap costs** — a
  cheap one for inserting/deleting a whole repeat unit, a stiff one for a single
  base — applied uniformly. Stutter is folded *into* the alignment rather than
  handled as a separate step. This was the in-house idea we found elegant and
  wanted to believe in.

The conventional move is to pick one — most likely A, because "HipSTR does it that
way" — and move on. Instead we built all three behind a single swappable interface
and measured them.

---

## 3. The experiment

We wrote a small **simulator** that generates cohorts of reads from *known* true
genotypes, so we always have a ground-truth answer to score against. To avoid the
trap of a model winning only on data that matches its own assumptions, we
generated under **three increasingly realistic regimes** and scored every dataset
with every model (a 3×3 matrix):

- **G1 — clean:** only whole-unit stutter (the easy case).
- **G2 — out-of-frame:** add the partial-length slips real data shows.
- **G3 — messy:** out-of-frame slips + sequence interruptions + a heavier
  base-error tail (the case that actually matters).

We measured four things per cell: **genotype concordance** (fraction of samples
called correctly), **allele-length error** (how far off the wrong calls are),
**calibration** (does the model's stated confidence match its real accuracy —
Expected Calibration Error, ECE), and **runtime**.

---

## 4. Results

### 4.1 The 3×3 scoring × generative matrix

16-sample, 3-locus synthetic cohort, depth 60, parameters taken from ground truth
(so the comparison is of the *models*, not of parameter estimation). `scoring ms`
is a debug build — the **ratio** between models is the meaningful figure, not the
absolute number.

| model × generative      | concordance | mean &#124;Δ&#124; | max &#124;Δ&#124; | ECE      | scoring ms |
|-------------------------|------------:|-------------------:|------------------:|---------:|-----------:|
| B classic  / G1 clean   |      1.0000 |             0.0000 |                 0 |   0.0000 |       1740 |
| A HipSTR   / G1 clean   |      1.0000 |             0.0000 |                 0 |   0.0000 |          9 |
| C 2-penalty/ G1 clean   |      0.6667 |             0.3229 |                 2 |   0.2981 |        339 |
| B classic  / G2 oof     |      1.0000 |             0.0000 |                 0 |   0.0000 |       2367 |
| A HipSTR   / G2 oof     |      1.0000 |             0.0000 |                 0 |   0.0000 |         10 |
| C 2-penalty/ G2 oof     |      0.0208 |             5.8750 |                11 |   0.9769 |        682 |
| B classic  / G3 messy   |      1.0000 |             0.0000 |                 0 |   0.0166 |       8522 |
| A HipSTR   / G3 messy   |      1.0000 |             0.0000 |                 0 |   0.0001 |         34 |
| C 2-penalty/ G3 messy   |      0.0417 |             5.5000 |                11 |   0.9382 |       1948 |

### 4.2 Robustness: concordance vs sequencing depth (on the messy G3 case)

| depth | A HipSTR | B classic | C 2-penalty |
|------:|---------:|----------:|------------:|
|    10 |    0.708 |     0.729 |       0.062 |
|    20 |    0.979 |     0.979 |       0.042 |
|    40 |    0.979 |     0.979 |       0.062 |
|    80 |    1.000 |     1.000 |       0.042 |

### 4.3 What the data said

- **Our own intuition (C) was wrong — and we only know because we built it.**
  On paper the two-penalty alignment is elegant. In practice it collapsed: on the
  messy, realistic case it called barely 4% of genotypes correctly, missing by ~5
  repeat units, with near-zero calibration. Folding stutter into the alignment lets
  noisy reads drag the genotype call toward the wrong allele. No amount of
  whiteboard argument would have settled this as decisively as one run did.
- **The literature choice (A) and the incremental fix (B) tied on accuracy** —
  both perfect across clean, out-of-frame, and messy data. The tie-break came from
  the two secondary axes: **A is far better calibrated** (its confidence is
  honest — ECE 0.0001 vs 0.017 on messy data) and **~100× faster per evaluation**.
  That speed matters because this function is ~75% of the genotyper's runtime.

**Decision:** productionize **A (HipSTR)**; keep **B** as a reference baseline;
**delete C**. The choice rests on numbers anyone can reproduce, not on who argued
hardest.

### 4.4 Two discoveries the experiment handed us for free

Empirical testing doesn't just rank the options you brought — it surfaces things
you didn't know to look for:

1. **Model C had a hidden length bias.** The first run had C calling *every*
   sample as the single longest allele. Investigating revealed its scoring
   systematically over-credited longer candidates. A standard normalization fixed
   that specific bug — and *then* C still lost, for the deeper reason above. Two
   distinct lessons from one cheap experiment.
2. **A pre-existing filter is too aggressive.** Running the messy case exposed
   that an upstream "is this locus a clean repeat?" check rejects a locus the
   moment a few out-of-frame reads appear — filtering data that is perfectly
   callable. That was nobody's hypothesis going in; it fell out of looking at why
   every model scored zero on G2/G3 at first. It's now logged as a separate fix.

---

## 5. The cost of doing it this way

The striking part is how little this empirical detour cost.

| What | Amount |
|---|---|
| Candidate models fully implemented & compared | **3** |
| Net code committed (interface + winner + reference + test harness + simulator extensions) | **~1,400 lines** |
| Code written, tested, then **deliberately discarded** (the losing Model C) | **~310 lines + ~10 tests** |
| Automated tests added (net) | **+20** (1,305 → 1,325) |
| Evaluation harness (simulator → score → genotype → compare, with metrics) | **~640 lines, committed and reproducible** |
| Compiler / lint / full-test-suite cycles run | **dozens**, each in seconds |
| Human-guided wall-clock | **a single working session** |

The "wasted" work — building Model C only to delete it — is the whole point. In a
slower world that ~310-line cost is exactly what deters you from testing the idea
at all, so you argue instead. Here it was cheap enough that *finding out our
favourite idea was wrong* was a good trade.

**Conventional estimate for the same scope.** Implementing three statistical
read-likelihood models (one requiring a study of the HipSTR C++ source, one a
novel dynamic-programming alignment), plus a ground-truth simulator with
calibration metrics, plus a cross-generative evaluation and the debugging of two
emergent issues, is realistically **2–4 focused engineer-weeks** for an
experienced bioinformatics developer. Most teams, facing that cost, would skip the
bake-off and just implement the literature default — and would never have learned
that their own elegant idea fails, that the default is also the fastest and
best-calibrated, or that an upstream filter is mis-tuned.

*(These figures are an honest after-the-fact accounting of one session, not a
benchmark. The runtime numbers are debug-build ratios; absolute timings need a
release-build measurement. Validation here is on synthetic data with known truth —
the appropriate tool for choosing a model's *shape*; calibrating absolute
parameter values against real data is a separate, later step.)*

---

## 6. The takeaway for the workshop

The scientific method runs on experiments, but in software we often substitute
argument for experiment because experiments are expensive to build. AI-assisted
coding lowers that cost enough to **restore the experiment as the default**:

- **Test intuitions instead of defending them.** Our in-house idea was appealing
  and wrong; we found out in one run rather than after shipping it.
- **Verify the literature rather than deferring to it.** "HipSTR does it this way"
  turned out to be right *here* — but we confirmed it on our own data and learned
  *why* (calibration + speed), instead of assuming it.
- **Let measurements break ties that opinion can't.** Two models tied on the
  headline metric; the decision came from calibration and runtime, measured.
- **Cheap experiments pay unexpected dividends.** The same run that ranked the
  models also exposed two unrelated bugs.

The deliverable of the session is not just "we chose Model A." It is a **committed,
re-runnable experiment**: anyone can re-execute the bake-off, change the
assumptions, add a fourth model, and watch the table update. That reproducibility —
an experiment you keep, not an argument you had — is what the speed buys you.

---

### Appendix — where the artefacts live

- Implementation report: `doc/devel/reports/implementations/ssr_stutter_scoring_model_bakeoff_2026-06-24.md`
- The swappable model interface + the production winner: `src/ssr/cohort/read_model/`
- The evaluation harness + metrics: `src/ssr/cohort/bakeoff.rs`
- The ground-truth simulator + generative regimes: `src/ssr/cohort/sim.rs`
- The original plan: `doc/devel/implementation_plans/ssr_stutter_scoring_model_bakeoff.md`
