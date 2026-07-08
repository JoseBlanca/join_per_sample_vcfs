# Cohort-recurrence-aware admission gates for the SSR caller

**Status:** draft for discussion — 2026-07-08. No code. Written to frame a
design conversation, not to prescribe an implementation.

## 1. The problem in one paragraph

Our microsatellite caller is more precise than HipSTR on the `ssr_tomato1`
benchmark — it makes about a quarter as many false variant calls — but it drops
some microsatellite variants that are genuinely real. When we restrict attention
to the loci that the reads *confidently* say are polymorphic (see §2), three
decision gates account for nearly all the loss. Each of those gates decides
whether to keep or reject a variant using **one plant's evidence at a time**, and
never asks whether the same allele shows up cleanly in the *other* plants. On
this cohort that matters a lot: tomato self-pollinates and is sequenced shallowly
(a median of three reads per plant per locus), so no single plant is convincing
on its own — but the cohort together often is. This spec proposes teaching those
gates to look across the cohort before they reject a call.

## 2. What the benchmark evidence says

We built two read-grounded "silver standards" — sets of loci we can label
confidently without a truth set, one from our own pileup reads and one from
HipSTR's per-read fields — and scored where each caller loses real variants.
(Method and numbers: the error-signals dashboard and the drop-attribution
diagnostic on this branch.)

Of the **561 loci both silver standards agree are truly polymorphic**, our caller
misses **108 (19%)**. Attributing each miss to the gate that caused it:

| Gate that dropped the real variant | share of the 108 misses | motif skew |
|---|--:|---|
| Allele-balance no-call (the "FP-control" step) | 43% | mostly dinucleotide |
| The genotyping loop settling on homozygous-reference | 33% | almost entirely dinucleotide |
| The periodicity filter rejecting the whole locus | 19% | impure / interrupted repeats |
| A variant length never entering the candidate list | 3% | — |

The false positives we *avoid* — which is where our precision advantage comes
from — are sporadic slippage artifacts that do **not** recur as clean calls. So
the losses and the wins separate on the same axis: **recurrence**.

### 2.1 What the first experiment already showed

Before proposing new code we ran the one change that needed none: the genotyping
loop's cohort-frequency prior (§6.3) already ships behind a switch, so we turned
it on and scored it against the silver standard. The result reordered the plan.

Turning it on made recall **worse**, not better: the confident-real misses rose
from 19% to 43%. It recovered 14 real loci the default missed but newly lost 149
it had been calling. The drop-attribution explains why — the prior lifted 327
loci out of the homozygous-reference trap, but the **allele-balance gate then
no-called 907 more calls**, so the loci the loop rescued died at the very next
gate. Precision barely moved (it was already excellent).

Two conclusions, both now evidence rather than guess:

- **The allele-balance gate is the binding constraint.** It caps recall no
  matter what the loop does upstream, so it is the change to make *first*.
- **The gates are coupled** (this spec's open question, §9.4, confirmed): the
  frequency prior can only help *after* the allele-balance gate stops discarding
  what it recovers. So the prior is not dead — it is blocked, and worth
  re-testing once the gate ahead of it is fixed.

## 3. The core idea, and the trap inside it

**Real versus artifact separates by recurrence — but only recurrence of the
right *shape*.**

A real microsatellite allele, in an inbred cohort, recurs as a **homozygous or
balanced** call: several plants carry it, and in each the reads sit squarely on
that length. A slippage artifact recurs only as a **lopsided minority shoulder**:
a thin tail of reads one repeat unit off the true length.

The trap — and the reason a naive "spare the call if the allele recurs" rule
would backfire — is that **slippage is systematic per sequencing chemistry, not
per plant.** A dinucleotide tract slips the same way in every library prepared
the same way. So "this minority allele appears in several plants" is satisfied by
systematic slippage *just as readily* as by a real variant, and a gate that
spared calls on bare recurrence would re-admit exactly the false positives we
currently reject — throwing away the precision advantage.

Two consequences shape the whole design:

- The cohort signal must be **recurrence of a clean allele** (homozygous or
  balanced), not recurrence of mere presence.
- We already model slippage **per sample-group (chemistry)** — the stutter shape
  and level are fit per group, not per plant. That gives a sharper definition
  than a hand-drawn "shape" rule: a clean recurrence is read support for the
  allele **in excess of what that group's slippage model predicts**, seen across
  independent plants. Reusing the existing per-chemistry model keeps the new
  signal consistent with how the caller already reasons about slippage.

This is the same quantity the silver standard already computes to label a locus
`true100`. So the signal is not a new invention — it is the internal counterpart
of the yardstick we are grading against.

## 4. Why this is cheap architecturally

Stage 1 (`ssr-pileup`) reads one alignment file per plant and writes that plant's
observed tract lengths and quality counts into its `.ssr.psp`. It genotypes
nothing and gates nothing; it cannot see the cohort, because it processes one
plant at a time. **None of the three gates lives there.**

All three gates run in Stage 2 (`ssr-call`), inside the per-locus genotyping
function, and they operate on the `CohortLocus` — a structure that already holds
**every present plant's reads at that locus**, side by side, at the moment each
gate fires. The allele-balance step even has the whole cohort call in scope
already; it simply looks only at the one plant it is deciding.

The practical upshot: **making these gates cohort-aware needs no change to
Stage 1, no re-running the pileup, and no change to the `.ssr.psp` format.** It is
a Stage-2 change to decision logic, on data that is already assembled. That is
what makes this attractive.

## 5. The shared signal (proposed)

Compute, **once per locus**, a per-candidate-allele **cohort-support score**: for
each candidate allele other than the reference, how many plants show that allele
as a clean call — read support above the plant's own chemistry slippage
prediction — and how strong that support is in total.

Stated plainly, the score answers: *"Setting slippage aside, how many independent
plants really carry this length?"*

The precise definition is deliberately left open for §9 to resolve, but the
inputs are all present in the `CohortLocus`: each plant's observed tract lengths
and counts, and the frozen per-group slippage model. Two candidate definitions to
weigh:

- **Count of clean carriers** — the number of plants whose excess-over-slippage
  support for the allele clears a small bar. Simple, robust, matches the silver
  standard's recurrence test.
- **Summed excess evidence** — the total read support for the allele across the
  cohort after subtracting each plant's predicted slippage. Continuous, so it
  degrades gracefully, but harder to threshold.

Computing this once and threading it into every gate — rather than giving each
gate its own ad-hoc recurrence check — keeps the gates consistent with one
another and gives us a single knob to reason about and validate.

## 6. How each gate would use the signal

The three gates need **different kinds of work**, not just different amounts.
Two are new heuristic code that consumes the shared signal (§6.1, §6.2); the
third is already-written machinery that only needs to be switched on and measured
(§6.3). None of the three is an open-ended design problem. The experiment in §2.1
sets the order: the **allele-balance gate (§6.1) goes first**, because it is the
gate everything else bottlenecks on.

### 6.1 Allele-balance no-call — the lead change (heuristic)

Today: a heterozygous call whose two alleles are supported by lopsided read
counts is scaled down and, below a quality floor, converted to a no-call. This is
a pure function of the one plant's read split.

Proposed: **before** demoting an imbalanced heterozygote, check whether its
minority allele is a clean, recurrent call elsewhere in the cohort. If it is,
keep the heterozygote untouched instead of no-calling it.

Concretely, and staying inside what the step can already see (every plant's
called alleles, their deconvolved per-allele read support, and quality — no read
model): in a single pre-pass over the locus, count for each allele how many
plants carry it as a **confident, non-stutter call** — either **homozygous for
it**, or a **balanced heterozygote** (both alleles well-supported). Then, when an
imbalanced heterozygote would be no-called, spare it if its minority allele's
count clears a small bar (a first cut: at least two such plants; balanced means a
minor fraction of at least 0.4).

Counting corroboration **only** from homozygous or balanced carriers is what
keeps this stutter-proof. Systematic slippage only ever appears as a lopsided
minority shoulder, so it can never present as a homozygous-alt or a balanced-het
carrier — a stutter allele's corroboration count stays at zero, and the gate
still rejects it. "The allele appears elsewhere" would re-admit slippage; "the
allele stands on its own feet elsewhere" does not.

This is a local change at a single, well-isolated step that already has the whole
cohort call in scope, needs no read model and no new data, and its effect is
directly measurable against the silver standard. It lands behind a switch
(default off → byte-identical) until the benchmark confirms it, matching how the
frequency prior was introduced.

### 6.2 Periodicity filter — the medium one (heuristic)

Today: the whole locus is rejected when too large a fraction of its pooled reads
sit off the modal-length grid — a guard against non-repeat regions.

Proposed: do not count off-grid support against a locus when that support forms a
**coherent, recurrent alternate length** across plants (a real interrupted or
impure repeat) rather than scattered noise. The gate already pools reads across
the cohort; the change is to weigh *systematic* off-grid signal differently from
*scattered* off-grid signal, using the same cohort-support score.

Still a heuristic, but it touches locus admission rather than a single plant's
call, so it needs more care about what "coherent" means.

### 6.3 The genotyping loop — mostly built already; the question is empirical

The genotyping loop is not plant-by-plant. It is an iterative fit
(expectation-maximisation — it alternates between estimating the cohort's allele
frequencies and re-scoring each plant's genotype under those frequencies), and it
**already** borrows strength across the cohort through the shared frequency it
feeds back as each plant's genotype prior. The dinucleotide collapse to
homozygous-reference is therefore not a gate ignoring the cohort. It is a
**self-limiting-frequency problem**: a real-but-rare allele keeps the estimated
frequency low, a low frequency biases every shallow plant toward
homozygous-reference, and calling them homozygous-reference keeps the frequency
low. The loop settles in the wrong basin.

**We have already solved this once, in the SNP caller.** The SNP path added an
empirical-Bayes, leave-one-out frequency prior with a flat first step. In plain
terms it does three things: it starts the fit from a neutral frequency rather
than the homozygous-reference fixed point; it lets each plant borrow the cohort's
belief in an allele *computed from all the other plants*, so that plant's own
shallow reads cannot drag the estimate down; and a convergence guard forbids the
fit from stopping at the flat starting point. Those three together are exactly
what climbs out of the bad basin.

**And that machinery is already ported to the SSR loop.** It is
`run_pi_em_marginalized` — described in the code as "the SSR analogue of the SNP
engine's `e_step_cohort_loo`" — carrying the same flat first iteration,
leave-one-out per-plant prior, and convergence guard. It is complete and tested,
sitting behind an off-by-default toggle (`PVC_SSR_MARGINALIZED_PRIOR`).

So this gate is **not** an unsolved design problem. The open question is
empirical, and it is one we could not answer before but can now. We benchmarked
this prior once, and the verdict was deferred: on `ssr_tomato1` it emitted about
**30% fewer** of HipSTR's polymorphic loci, and we could not tell whether that
was a loss or a gain. The reason we could not tell is the reason it was
inconclusive: the comparison was against HipSTR, whose polymorphic set we now
know is contaminated with false positives, and the marginalised prior is also
**inbreeding-aware** — it correctly demotes the lone heterozygotes that a selfer
should rarely produce. So "fewer emissions" mixed together *real variants
recovered from the homozygous-reference trap* and *HipSTR false positives
correctly rejected*, with no way to separate them.

The silver standard separates them. So the concrete next step for this gate is
not new design — it is to re-run the existing marginalised prior and score it
against the two numbers the earlier benchmark lacked: does it **recover the
confident-real misses** the plug-in loop drops (the recall we are chasing), and
does the **silver-standard false-positive rate hold** (the precision we are
protecting)?

The one genuine unknown is depth. The trap-escape was designed and validated on
SNP data at ordinary depth; whether it is strong enough at three reads per plant
is exactly what the re-benchmark would tell us. If it is not, the shared
cohort-support signal (§5) offers a lever the SNP path did not need: seed the fit
from the clean-recurrence score, so it begins even further outside the bad basin.
That stays inside the probability model — an initialisation choice, not a
post-hoc override — so the determinism guarantees the genotyping core carries are
preserved.

## 7. How we will know it worked

Every change is gated on the same acceptance test, run on `ssr_tomato1`:

- **Recall:** the count of confident-real misses (silver-standard `true100` loci
  we drop) must fall.
- **Precision:** the silver-standard false-positive rate must **not** rise. This
  is the guardrail — the whole point is to recover real variants *without*
  re-admitting the systematic-slippage false positives, and the precision number
  is where that would show up.

A change that lifts recall while holding precision is a win; one that trades
precision for recall is not, because precision is the advantage we are protecting.

## 8. Scope and non-goals

- **Stage 1 is untouched.** No re-pileup, no `.ssr.psp` format change.
- **The read model and the slippage model are reused, not rewritten.** The new
  signal is built from what they already produce.
- **We are not trying to match HipSTR's recall.** A large share of the loci
  HipSTR calls that we miss are HipSTR's own false positives; matching them would
  cost precision. The target is our own confident-real misses, nothing more.
- **Incremental.** The three gates land and are measured one at a time. The
  allele-balance gate goes first (§2.1 showed it is the bottleneck); the
  frequency prior is re-tested once that gate no longer discards what it
  recovers.

## 9. Open questions for discussion

1. **The exact cohort-support score.** Clean-carrier count or summed
   excess-over-slippage evidence — and what bar counts as "clean" at three reads
   per plant, where the score is itself noisy.
2. **Does the existing marginalised prior recover the trap misses? (§6.3)** This
   is now a measurement, not a design: re-run `PVC_SSR_MARGINALIZED_PRIOR`, score
   it against the silver standard's confident-real misses and false-positive
   rate, and read off whether the ported SNP-path trap-escape is strong enough at
   three reads per plant. Only if it is not do we reach for the seed-from-
   recurrence lever — and even that is an initialisation choice, not new
   machinery. This is arguably the cheapest experiment of the three gates,
   because the code already exists.
3. **Chemistry handling.** The slippage model is per sample-group. Does the
   cohort-support score aggregate across groups, or should recurrence *across*
   groups count for more than recurrence within one group (which a group-wide
   library artifact could fake)?
4. **Gate interaction and ordering.** The gates run in sequence. If two of them
   become cohort-aware, do their signals need to agree, and does one change shift
   what the next one sees?
5. **Circularity.** The score is built from the same reads the gates judge. We
   need to be sure a plant's own shaky evidence does not both create the
   recurrence signal and get spared by it.
