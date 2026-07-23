# ng — the read alignment module

*Status: design spec, 2026-07-23. Defines `src/ng/alignment/` — the algorithms that line a read up
against a reference sequence. This module knows only about alignment: what an alignment is, how to
compute one, and how to score one. It has no knowledge of which caller step calls it or why.
Grounded in the frozen production code (`src/ssr/pileup/alignment.rs`, `src/ssr/cohort/`), which ng
ports and adapts but never edits. **No code yet.** Naming: **STR** in prose, `ssr` in code.
Code-facing companion: [`../arch/alignment.md`](../arch/alignment.md) (the types and trait
signatures).*

---

## 1. What this module is

A variant caller decides what a sequencing read says by lining the read's bases up against a
stretch of reference. That line-up is an **alignment**. This module owns the algorithms that
produce and score them, and nothing else.

It provides two independent things:

- **Alignment algorithms** — given a read and a reference sequence, either return the single most
  probable line-up, or return the total probability of the read summed over every possible line-up
  (§3–§5).
- **Alignment normalization** — given an alignment that already exists, rewrite it into a canonical
  form so that equivalent alignments get an identical spelling (§6).

These are separate tasks with separate interfaces. Computing an alignment and tidying an existing
one have nothing in common beyond the word "alignment".

**Out of scope — deliberately.** This module does not know which step calls it, does not decide
which algorithm a given read deserves, does not fetch reads or reference, does not tally
observations, and does not generate candidate alleles or genotypes. A caller picks an algorithm and
supplies the inputs; interpreting the result is the caller's business. Keeping that boundary strict
is what lets an algorithm be swapped and measured without touching anything else.

**One shape deliberately absent: per-base alignment confidence.** A third thing an alignment can
yield is, for each read base, how sure we are it sits at the right reference position — which is what
base alignment quality (BAQ) computes, with a forward-backward pass, to cap quality scores near
ambiguous indels. It is neither an alignment nor one total probability, so neither interface below
covers it. **Decided 2026-07-23: this module does not provide it**, on production experience that BAQ
is not worth carrying forward. Recorded rather than omitted, because it is a real answer shape and
the decision may be revisited — if it is, it needs a third interface, and it is the same machinery a
posterior over the repeat's boundary would need.

---

## 2. Vocabulary

The reader knows the biology; these are the alignment and statistics terms this spec leans on.

- **Best-path alignment.** Out of all the ways a read could line up against a reference, the single
  most probable one. The result is *one* line-up. (The classic Viterbi / Smith-Waterman answer.)
- **Marginal alignment.** Instead of picking the best line-up, add up the probabilities of *all* of
  them. The result is one number — the total probability the reference produced this read, with the
  line-up uncertainty averaged out rather than guessed. The best-path score is simply the largest
  single term in that sum. Computing the sum is the *forward algorithm*.
- **Affine gap.** The usual way to price an insertion or deletion: one cost to open the gap, a
  smaller cost per base to extend it — so one long gap is cheaper than several short ones.
- **Tract and flanks.** In a microsatellite, the **tract** is the repeat itself (the run of units)
  and the **flanks** are the unique sequence either side.
- **Stutter.** During PCR (and, less often, sequencing) the polymerase gains or loses **whole
  repeat units**. A read from a true 10-unit allele routinely shows 9 or 11. An alignment model that
  prices such a change as an ordinary gap gets microsatellites badly wrong.
- **Emission.** How the algorithm scores one read base against one reference base — the cost of a
  match versus a mismatch. It can be driven by the read's own per-base quality score, or by a single
  flat error rate.
- **Banding.** Restricting the alignment to line-ups that stay within a fixed distance of the
  expected diagonal. Cuts the work from "read length × reference length" to "read length × band
  width", and is safe whenever the read's position is already pinned closely enough.

---

## 3. What distinguishes one algorithm from another

Two properties define an algorithm, and both are intrinsic to it — neither depends on who calls it.

**What it returns.** Either the single best line-up, or the marginal probability over all line-ups.
This is not a tuning choice: the two answer different questions and are computed differently. They
share a recurrence — the marginal replaces "take the best predecessor" with "add up all
predecessors" — but the best-path version also needs a traceback to recover the line-up, and the two
have different performance ceilings (§4). They are separate algorithms.

**What sequence changes it models.** Either plain affine gaps, or a repeat-aware model that (a)
prices gaps differently in the flanks and inside the tract, and (b) treats a whole-unit slip as its
own kind of event rather than as a gap of that many bases.

That gives four combinations, of which three are in the initial set:

| | **affine gaps** | **repeat-aware** (flank/tract split + whole-unit slips) |
|---|---|---|
| **best single alignment** | affine best-path aligner (§4.1) | repeat-aware best-path aligner (§4.2) |
| **marginal over all alignments** | *not in the initial set* (§5.3) | repeat-aware marginal aligner (§5.1) |

Two further properties are **swappable knobs on any of them**, not separate algorithms:

- **Emission** — per-base quality, or a flat error rate. Which scores better is measured, not
  assumed, so both are configurations of the same algorithm.
- **Banding** — whether the matrix is restricted to line-ups near the expected diagonal, and how
  wide that band is. A performance choice that does not change what an algorithm *means*, only what
  it costs — **provided the band is wide enough**. Too narrow and it silently loses long alleles
  rather than merely costing accuracy: the one failure mode here that does not announce itself.

  **The floor is arithmetic; the headroom above it is not.** These are *global* alignments — the read
  must be consumed entirely, and so must the reference stretch. So when the two differ in length, the
  path is **forced** that far off the diagonal before it can reach the end. The band therefore has a
  hard floor of **|read length − reference length|**, computed per read, below which the true path is
  provably lost. Above the floor it still needs **headroom**, because an optimal path can stray
  further than the minimum and come back — an insertion, then a deletion. How much headroom is an
  empirical choice, not a derivation (§9).

  **Do not derive the band from the slip cutoff.** An earlier draft of this spec set the width from
  the largest slip worth scoring, times the period. That is wrong twice over. It is **numerically too
  small**: production's own long-allele test aligns a 52-base read against a 28-base frame, which
  needs 24 cells of deviation where that formula allows 22 — the test would fail. And it is **wrong
  in kind**, because the slip cutoff is a *scoring* parameter. Letting it set the ruler's band would
  make the measurement blind to any length the scoring model finds implausible, which is exactly the
  coupling §4.2 exists to prevent: a 15-unit expansion is data, and the ruler's job is to measure it,
  not to pre-judge it. It would also mishandle the reads that run off their own end mid-repeat, since
  those are the ones differing most in length from the frame.

  **Production is not uniformly banded, and a port must not assume it is.** Its repeat delimiter
  fills its whole matrix — no band at all. Its repeat marginal restricts the *computation* to a band
  but still allocates the full matrix. The consequence for the build order is in §10.3.

---

## 4. Best-path algorithms

Both return the single most probable line-up of a read against a reference sequence, and both need a
traceback to recover it. They differ only in how they price sequence change.

### 4.1 Affine best-path aligner

Standard match / mismatch / insertion / deletion alignment under affine gap costs. One gap model
everywhere along the reference. This is the general-purpose aligner: given a read and a reference
stretch, where does the read go and what does it show.

**It scores with per-base qualities** (decided 2026-07-23) — the read already carries a confidence
for every base, and throwing that away to align it would discard information we paid for. The
consequence is deliberate and worth naming: it **rules out the fast best-path cores** (see the note
below), because those work in fixed costs and cannot consume a per-base quality. Speed is the price
of using the read's own information. A flat-cost configuration of this same algorithm remains
available, and is what the repeat-aware comparison uses as its quality-blind end (§10).

### 4.2 Repeat-aware best-path aligner

The same alignment, but with **two gap regimes along the reference**: a stiff gap in the flanks and
a soft gap inside the tract. The reasoning is that flanks are clean unique sequence and should hold
the tract's edges in place, while length variation inside the tract is expected — it is the repeat
changing size, not an error. With a single stiff gap everywhere, any allele more than about one unit
from the reference collapses back onto the reference; production hit exactly that failure before
splitting the regimes.

Because the flanks stay anchored, the two flank-to-tract boundary columns of the resulting line-up
say where the read's repeat starts and ends — which is how this aligner measures a read's repeat.
Production's tract delimiter is this algorithm
([alignment.rs:171](../../../../src/ssr/pileup/alignment.rs#L171)); the details below are what an
implementation has to get right, and each one is load-bearing.

**Which reference columns count as inside the tract.** A gap is priced at the tract rate when it
inserts beside, or deletes, a **tract** base — that is, when the reference column it touches lies
after the last left-flank base and no later than the last tract base. Every other column, including
both junctions' outer sides, keeps the stiff flank rate. This is what keeps the junctions anchored
while the interior is free to change length.

**The costs.** Production's values, as probabilities (the algorithm itself works in logarithms):

| transition | in the flanks | inside the tract |
|---|---|---|
| open a gap (match → insertion, or match → deletion) | 2.9 × 10⁻⁵ | 1 × 10⁻² |
| extend an open gap | e⁻¹ ≈ 0.368 | e⁻¹ ≈ 0.368 (shared) |
| close a gap (back to match) | 1 − extend | 1 − extend (shared) |

The tract gap-open is roughly 350× softer than the flank one. The flank value is the Dindel model's
base rate for short homopolymer runs; the tract value was set by the investigation that found the
collapse described above. Only the **open** cost switches regime: extension and closing are shared,
because the open is the dominant fixed cost and the extension is a knob nobody has yet needed to
turn. **Both numbers are provisional development calibration**, to be reconciled against the
measured slip rate — treat them as starting points, not as findings.

**A known inconsistency, to reproduce or fix deliberately — not by accident.** Production computes
the match→match probability once, from the *flank* gap-open (`1 − 2 × 2.9e-5`), and never recomputes
it under the tract regime. So inside the tract the three transitions leaving a match sum to about
1.02 rather than 1 (0.999942 + 0.01 + 0.01); in the flanks they sum to exactly 1.

**Be precise about what that skews.** It is *direction*-symmetric — insertion and deletion get the
same tract open — so it does not favour longer over shorter. But normalising would divide all three
by 1.02 **at every departure from a match**, and a path carrying an *L*-base gap passes through about
*L* fewer matches than the diagonal does. So the un-normalised model quietly penalises the
length-changing path by roughly 1.02 per gapped base — about 1.2× at a 10 bp change, 1.8× at 30 bp.
The bias is toward **the reference length**, and it grows with the size of the change. That is the
same class of bias the tract regime was introduced to remove; it is dominated today by the 350× open
ratio, which is why it has never shown up in results, but anyone who narrows that ratio or lengthens
the tract inherits an unstated pull toward reference. An implementer should know it is there:
reproducing it silently and "fixing" it silently both end with a parity test failing for a reason
nobody can find.

**How to price a length change inside the tract — two models, to be compared.** Biologically the
common event inside a tract is gaining or losing a **whole unit**, so pricing a unit slip separately
from an ordinary indel is the more faithful model, and it is what HipSTR does. Production's
delimiter does not — it charges one flat per-base gap anywhere inside the tract. Both are worth
building and measuring; neither is settled.

| | **one penalty** (production today) | **two penalties** (HipSTR-shaped) |
|---|---|---|
| a whole-unit slip | priced as a gap of that many bases | its own, much cheaper price |
| any other length change | the same per-base price | a separate, dearer price |
| does it need the motif? | no — only the tract boundaries | yes — it must recognize a unit |

The **one-penalty** model is content-agnostic by design: it measures whatever the read shows,
verbatim, including repeats that are interrupted or whose length is not a whole number of units.

The **two-penalty** model is closer to the biology. Two things must be right for it to work at all:

- **Gaining and losing must be priced separately.** Stutter is markedly asymmetric — contractions
  are more common than expansions. HipSTR carries separate expansion and contraction probabilities
  and fits them per locus, and the fitted values are contraction-biased. Production's scoring model
  already does this ([hipstr.rs](../../../../src/ssr/cohort/read_model/hipstr.rs) derives its
  expansion and contraction mass from a separate direction split, with a test asserting a one-unit
  contraction outscores a one-unit expansion). A two-penalty delimiter must inherit that asymmetry;
  a symmetric one would be a step backwards from production's scoring. **The parameters, the
  formulas and the two silent conversion traps are set out in full in §5.2** — this aligner and the
  repeat-aware marginal use the same stutter model, so it is written once and shared.
- **Out-of-frame changes must keep a route.** Making unit slips cheap is only half the model. HipSTR
  pairs its in-frame slip model with a **separate out-of-frame model** for changes that are not unit
  multiples. Without that second route, a read whose repeat is interrupted or not a whole number of
  units has nowhere to go, and the aligner will force it onto the nearest tidy length.

Two further facts about the pair decide how the comparison between them has to be run.

**At period 1 the two models do *not* become one.** An earlier draft of this spec said they did, and
excluded homopolymers from the comparison on that basis — which skipped the very loci where the indel
deficit lives. What actually collapses at period 1 is only the arithmetic in-frame/out-of-frame split
(see below). Four differences survive, and all four are live at a homopolymer:

- **direction** — the two-penalty model prices gaining and losing separately; the one-penalty model
  uses the same gap cost for an insertion and a deletion, so it is **structurally unable** to say
  that losing a base is likelier than gaining one;
- **size decay** — an affine gap costs `open + (n−1)·extend`; a geometric costs
  `geom·(1−geom)^(n−1)`. These agree only at one particular `geom`, which nothing sets;
- **placement multiplicity** — the two-penalty model divides by the number of reachable placements
  and the one-penalty model does not, and at period 1 that count is at its largest (a base inserted
  anywhere in a run of *n* gives the same sequence *n*+1 ways);
- **what counts as stutter at all** — see immediately below.

The direction difference is worth dwelling on: the recorded production failure is stutter reads
pulled to reference at homopolymers, real stutter is contraction-biased, and a direction-symmetric
model cannot express that bias. That is a hypothesis, not a finding — but it is why period 1 is the
case this comparison most needs, not the one to skip.

**The in-frame test is arithmetic only — it never asks *what* was inserted.** A length change is
called in-frame when it is a multiple of the period, full stop. So inserting a `C` into a poly-A
tract is classified as one unit of stutter exactly as inserting an `A` is, though only the `A`
actually is one. The composition *is* caught, but downstream and through the wrong channel: the
resized candidate tiles the motif, so the odd base scores a mismatch and pays a **base-error**
price rather than an out-of-frame one. At period 2 and above the same event would route to the
out-of-frame sub-model and get its own, rarer price — same biology, two prices, decided by period.
The effect is not exclusive to period 1 (inserting `GG` into a `CA` repeat is divisible by the
period too), but at period 1 it happens to roughly three of every four single-base insertions
instead of occasionally. **Consequence for the comparison:** at period 1 the two sub-cases must be
scored apart — an indel *of the repeat's own base* (genuine stutter, where direction asymmetry
should favour the two-penalty model) and an indel of *a different base* (not stutter, where the
two-penalty model mis-routes it and the one-penalty model's plain gap may well be the better
answer). Reported together, the two effects cancel in the average.

**What the comparison is testing, and the evidence that already exists.** The risk with a
two-penalty model in a *ruler* is that making tidy in-frame lengths cheaper pulls the measurement
toward the answer the model already expects, quietly rounding interruptions away — the aligner stops
being neutral about what it measures. That is the hypothesis to test, not a reason to skip the test.
There is directly relevant prior evidence: production built a two-penalty pair-HMM for the
**scoring** job and eliminated it, because "even self-affinity-normalized it collapsed on
out-of-frame reads" ([read_model/](../../../../src/ssr/cohort/read_model/) records the comparison and
its outcome). That was a different job, so it does not settle this one — but it sets the bar, and it
names the out-of-frame route as the part most likely to break.

Whichever model wins here, whole-unit slips are **also** modelled downstream, by the genotyping
likelihood that scores an already-measured repeat against a candidate allele (§5.1). Keeping some
separation between the two is what lets base error and stutter stay separately identifiable: if the
ruler and the model of length change share every assumption, a stutter assumption biases the very
measurement meant to test it.

**Emission (per-base quality by default).** Each read base is scored against a reference base using
the read's own quality score. With the error probability ε = 10^(−Q/10) for quality Q, a match scores
1 − ε and a mismatch ε/3 — the error split evenly across the three other bases. Two details carry
over from production: a Q0 base has ε = 1, so its match score is floored to a tiny positive number
rather than zero, which would otherwise annihilate every path through it; and an **inserted** read
base, having no reference base to compare against, scores against a uniform base composition (1/4).
Swapping this component for a flat error rate gives the quality-blind configuration (§3).

**Reading the repeat off the line-up.** Walk the traceback and record, for each reference base
consumed, how many read bases were consumed before it. The read's repeat is the stretch of read
between the count at the first tract base and the count at the first right-flank base. Insertions are
counted as soon as they occur, which is what makes the assignment rule below concrete.

**Determinism.** Two reads of the same molecule must measure the same repeat, and the answer must
never depend on thread count or the order reads arrive in. Two rules settle every tie:

- When several predecessors score equally, prefer **match, then deletion, then insertion** — both
  within the alignment and when choosing which state the alignment ends in.
- A gap landing exactly on a flank/tract junction belongs to the block on its **5′ side**: an
  insertion at the left junction joins the flank, one at the right junction joins the repeat.

Without both rules, equally-scoring line-ups can measure different repeats from identical input.

> **Faster cores — considered and not taken.** Wavefront alignment (WFA), Myers' bit-vector method
> and the difference-recurrence core (ksw2) are far faster than a plain alignment matrix, but each
> computes a best path *only*, in fixed costs. Their specific tricks — a furthest-reaching frontier
> indexed by score, ±1 differences packed into machine words, saturating integer maxima — all depend
> on taking a maximum, and a sum destroys every one of them. None can express per-base quality
> emissions or the two-regime gaps either, so none can implement §4.1.
>
> **Be careful what that closes.** It rules out *those three cores*, not vectorisation as a technique.
> Anti-diagonal and striped vectorisation are reduction-agnostic, and GATK — cited in §5 as the home
> of the whole-read forward — ships a vectorised forward pair-HMM. So a fast *marginal* remains
> possible; it is simply a different kernel (floating-point, log-domain) rather than one of these.
> What is closed is the shortcut, not the speed.

---

## 5. Marginal algorithms

These return the total probability of the read summed over all line-ups, not a line-up. They need no
traceback, and they are banded for performance.

### 5.1 Repeat-aware marginal aligner

The probability that one sequence produced another, summed over the ways it could have, with gaps
allowed only near the ends. It scores with a **single flat error rate**, not per-base qualities
([pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs) is production's).

Two cases, and the common one is degenerate:

- **equal lengths** — there is exactly one way to line the two sequences up, so the sum has a single
  term and the result is a straight base-by-base comparison under the error rate.
- **different lengths** — the difference has to be absorbed by gaps, and those may fall only within
  a couple of bases of either end. The sum then runs over the few ways that can happen: a genuine
  forward pass.

**Gaps in the interior are forbidden, and that restriction is load-bearing.** An indel in the middle
of the repeat is not competed against a gap at the end; it simply scores far below one. That is
deliberate. Length change inside a repeat is what the **stutter** model explains (§5.2), and letting
this algorithm explain it too would make the two indistinguishable — you could no longer tell a
sequencing error from a slip. Confining gaps to the ends is what keeps the base-error rate and the
stutter rate separately estimable. Production allows two bases either side.

**A known normalisation defect, to reproduce or fix deliberately.** Over equal-length inputs this is
a proper distribution. The unequal-length case adds end-gap paths on top of it, so summed over all
differing-length outputs the total mass slightly exceeds one. It is harmless where it is used today
(the genotyping normalises per read) and the equal-length case dominates overwhelmingly — but an
implementer should know it is there rather than reproduce it unawares.

**What this algorithm is not: it does not model stutter.** It does not decide how likely a length
change was. That is a separate job, and the boundary matters:

> **Where the boundary falls.** Production's read-likelihood model computes *P(observed repeat |
> candidate allele)* in two steps: how likely is this length change (stutter), times how well the
> letters match once the candidate has been stretched to the observed length. **Only the second step
> is alignment** — the first is a model of how polymerase slips. Production already keeps the two in
> different files, the stutter half sitting with the genotyping code. This module owns the comparison
> only; the genotyping likelihood composes it. An earlier draft of this spec bundled both and called
> the bundle an alignment algorithm, which put a genetics model inside a module that is supposed to
> know only about alignment.
>
> **A consequence worth knowing before porting.** Because the stutter half always stretches the
> candidate to the observed length, the *unequal-length* branch of this comparison is **never reached
> from that caller** — production exercises it only from a different call site. A test built around
> the read-likelihood model would therefore never touch the forward pass at all. The unequal-length
> case has to be driven deliberately, or it ships untested.

### 5.2 The stutter model — shared, and not an alignment algorithm

How likely each length change is. It is set out here despite not being an alignment algorithm,
because **two consumers share it**: the two-penalty best-path aligner in this module (§4.2), and the
genotyping likelihood outside it (§5.1). One description, so the two cannot drift apart. It is
HipSTR's model.

**The stutter model has two regimes, and that split is its defining structure.** A read's length
change is either a whole number of repeat units or it is not, and the two mean different things:

- **in frame** — the change is a whole number of units. This is slippage, the common event, and its
  size is measured in **units**.
- **out of frame** — it is not a whole number of units. This is a sequencing indel or an
  interruption, not slippage; it is rarer, it gets its **own** parameters, and its size is measured
  in **base pairs**.

Each regime splits again by direction, because stutter is asymmetric: **losing units is more common
than gaining them**.

**Which regime applies is decided by arithmetic alone** — is the length change a multiple of the
period? — and **never by what was actually inserted**. So an insertion that happens to be
period-sized is treated as slippage whether or not its bases are the repeat unit. The composition is
caught downstream, as a base-error mismatch against the re-tiled candidate, rather than as an
out-of-frame event. The mis-routing is worst at period 1, where it catches roughly three of every
four single-base insertions; §4.2 works through what that costs and what the comparison must do
about it.

**The parameters — seven.** HipSTR's field names are in brackets:

| quantity | what it is |
|---|---|
| `equal` [`log_equal_`] | probability the read shows the allele's length unchanged |
| `in_up`, `in_down` [`in_up_`, `in_down_`] | probability of an in-frame expansion / contraction, of any size |
| `in_geom` [`in_geom_`] | how fast in-frame slip size decays, per **unit** |
| `out_up`, `out_down` [`out_up_`, `out_down_`] | the same two, for out-of-frame changes |
| `out_geom` [`out_geom_`] | how fast out-of-frame size decays, per **base pair** |

The unchanged mass is whatever the four direction masses leave:
`equal = 1 − in_up − in_down − out_up − out_down`.

**The distribution.** Write Δ for the read's length change in base pairs and *p* for the repeat's
unit length. Both regimes are a direction mass times a geometric over size:

```
in frame  (Δ divisible by p), n = Δ/p units:
    n = 0   →  equal
    n > 0   →  in_up    · in_geom  · (1 − in_geom )^(n − 1)
    n < 0   →  in_down  · in_geom  · (1 − in_geom )^(|n| − 1)

out of frame  (Δ not divisible by p), e = Δ − Δ/p  (truncated division):
    e > 0   →  out_up   · out_geom · (1 − out_geom)^(e − 1)
    e < 0   →  out_down · out_geom · (1 − out_geom)^(|e| − 1)
```

A geometric with success probability *g* puts mass *g* on the first step and multiplies by (1 − *g*)
for each step after, so a **larger** *g* concentrates the mass on single-unit slips. HipSTR ships
0.95 — nineteen slips in twenty are exactly one unit. Note *g* is a **success** probability, not a
decay; the two are complements, and confusing them is the first of the traps below.

**Why out-of-frame sizes are re-indexed.** The out-of-frame geometric is indexed by
`e = Δ − Δ/p` (truncated division), not by Δ. The reason is **not** double-counting — the two regimes
are disjoint by construction, since a change is out-of-frame precisely when it is *not* a multiple of
the period, so no length can reach both. What the re-indexing does is **compress the ranks**: it maps
the out-of-frame values onto consecutive integers so the geometric's support has no gaps. At period 3
the out-of-frame Δ values 1, 2, 4, 5, 7 map to e = 1, 2, 3, 4, 5. Without it the geometric would be
evaluated at indices that skip the multiples, distorting the distribution.

**A slip can land in more than one place — in frame.** In a pure repeat, adding a unit anywhere gives
the same sequence and there is nothing to choose. In an **interrupted** repeat the placements give
genuinely different sequences, so the model enumerates them and sums with equal weight (one over the
number of placements). If a slip cannot be reached from the allele at all — contracting away more
units than exist — the term is zero. **Production does this only for in-frame slips**; out-of-frame
resizing uses a single placement at the end of the repeat, a documented simplification of the same
class as the out-of-frame mass below.

**Slips past a cutoff score zero**, rather than being explained by an implausibly large slip. Such a
read falls to the model's outlier handling instead. **Mind the units:** production applies one
constant to the *unit* count in the in-frame branch and to the re-indexed *base-pair* count in the
out-of-frame branch. One number, two scales — decide deliberately whether that is intended rather
than inheriting it.

**Where the parameters come from.** HipSTR does **not** rely on its shipped values in normal use: it
fits expansion, contraction and the geometric **per locus by expectation-maximization**, jointly with
the genotype, so the rates adapt to the motif and the repeat's length. Production derives them per call instead,
from the per-locus stutter shape (its direction split and decay) and the per-read stutter level. One
piece is currently a placeholder rather than an estimate: **the out-of-frame mass is a fixed small
fraction of the in-frame mass** (5%). A real out-of-frame estimator — binning out-of-frame reads
separately in the parameter pre-pass — is a known and unfinished follow-up, and any comparison
involving out-of-frame reads inherits that weakness.

**Two conversion traps, both silent.**

1. **A decay is not a geometric success probability.** If a parameter is expressed as the probability
   of *continuing* to the next step (mean size 1/(1 − decay)), it is the complement of HipSTR's
   geometric, whose mean is 1/geom: `geom = 1 − decay`. Getting this backwards inverts the size
   distribution — large slips become common — and nothing crashes.
2. **Clamps carry weight.** The geometric probabilities are held strictly inside (0, 1) and `equal`
   is floored, so a parameter combination that would otherwise drive the unchanged mass negative
   degrades instead of producing a negative probability.

**On the defaults — quote them in matched sets.** HipSTR has two, and mixing them yields a pairing
that exists nowhere:

| | in-frame geom | in up / down | out-of-frame geom | out up / down |
|---|---|---|---|---|
| shipped default | **0.95** | 0.05 / 0.05 | **0.95** | 0.01 / 0.01 |
| the EM's starting point | **0.9** | 0.1 / 0.1 | **0.8** | 0.01 / 0.01 |

The second row is an *initialisation the fitting immediately moves away from*, not a default. (An
earlier draft of this spec paired 0.9 with 0.05/0.05, which is one number from each row.)

Note the shipped row also makes expansion and contraction *equal* — but that symmetry is a starting
point, not a claim: HipSTR's fitted values are contraction-biased. And note that HipSTR keeps the
in-frame and out-of-frame geometrics as **independent** parameters, with genuinely different values
in its EM start. **Production ties them to a single value** — a second undeclared placeholder
alongside the out-of-frame mass, and it should be recorded as one rather than mistaken for a fitted
result.

Finally, stutter depends strongly on library chemistry: PCR-free preparation reduces it several-fold
against PCR-plus, and by different factors per motif length. So these parameters belong to a **sample
group**, not to the cohort — which is the shape the cohort parameter design already uses.

### 5.3 Affine marginal aligner — not in the initial set

The marginal counterpart of §4.1 (the pair-HMM forward that GATK uses to score a read against a
candidate haplotype) is **not built initially**: no current consumer needs it. Production's
general-purpose caller derives its read likelihoods from summed per-base error probabilities over
already-decomposed observations, without an alignment step. Building this algorithm is therefore an
addition to compare against that approach, not a port of it — recorded in §9, not assumed.

---

## 6. Alignment normalization

A separate operation with its own interface: it takes an alignment that already exists and rewrites
it, rather than computing one.

The same insertion or deletion can be written at several equivalent positions when it sits in or
near a repeat — the gap can slide without changing a single base of the result. Left-alignment
shifts every gap as far left as it can go without introducing a mismatch, and merges adjacent gaps
of the same kind. Nothing about the alignment's quality changes; only its spelling. That matters
because equivalent alignments spelled differently look like different variants, and their supporting
reads scatter instead of pooling.

This does no dynamic programming and shares no machinery with §4 and §5. It is in this module
because it operates on alignments, and it is a separate interface because it is a separate task.

**Three algorithms — and neither existing one is provably leftmost.** An earlier draft of this spec
said there was nothing here to compare, on the belief that production does a naive single pass and
freebayes iterates to convergence. Both halves were wrong.

- **1a — one structured pass.** Production's: a single right-to-left traversal that **merges**
  consecutive indels and **carries an indel across alignment blocks** rather than stopping at the
  first. A port of GATK's left-aligner. The failure it is often assumed to have — "shifting one gap
  opens room for the next, and a single pass misses it" — is the case the merge and the propagation
  exist to handle.
- **1b — repeated simple passes.** freebayes': re-run a simple shift pass until nothing moves,
  capped at 20 iterations. It **reports** exhaustion, and its own caller ignores the report — so
  freebayes ships non-convergent results silently too.
- **1c — the structured pass to a fixpoint.** 1a in a loop until nothing moves, with the iteration
  cap as a safety net that **fails loudly** rather than returning a half-normalised result. A thin
  wrapper over 1a, not a third implementation of the shifting itself.

**This comparison has an oracle — the only one in this module that does.** "Leftmost" is a
*definition*, not a preference. For any output, take each indel, shift it one base left, and ask
whether the result still represents the same read against the same reference. If it does, the output
was not leftmost. So these three can be graded **against the definition** rather than against each
other, which no other comparison here can do — every other one scores algorithms against a rival or
against production. The property test is worth more than the bake-off, and it is what would justify
1c winning by construction rather than by preference.

**Run the cheap screen first.** Before asking which is better, ask whether they differ at all: run
them over the same real reads and count the outputs that disagree. If that count is near zero,
normalisation cannot explain any difference in calling, and a whole avenue closes for the price of a
single run.

**What this is and is not a fix for.** Production calls indels worse than freebayes, and the failures
concentrate at homopolymers and short repeats — exactly where indel placement is most ambiguous, and
therefore where two normalizers are likeliest to disagree. That makes normalization a live hypothesis
for part of the gap: differently-spelled equivalents split one variant's supporting reads across
several weak candidates and none survives. The screen above is what tests it.

**What the freebayes comparison does and does not establish.** freebayes has **no stutter model at
all**, so a stutter model is evidently not *necessary* to beat production here — worth knowing. But it
does **not** identify production's defect. A competitor that lacks a mechanism can avoid a failure for
entirely unrelated reasons: different candidate generation, a different prior, different filters. And
the recorded diagnosis is that stutter reads are *mis-assigned to reference* at repeats — a statement
about how production handles them, not about a feature it lacks. So normalization is **a** candidate,
the cheapest to test and the only one this module owns; production's reference-favouring genotype
prior and its filters are the others, and they live elsewhere. *(An earlier draft of this spec said
the freebayes comparison "rules out" stutter as the cause. It does not — absence of a mechanism in
the competitor is not evidence about the incumbent's defect — and that inference should not be
repeated.)*

---

## 7. The interface — chosen at compile time, never dynamically dispatched

These algorithms run **per read across a whole cohort** — millions to billions of calls. So which
algorithm runs is resolved at **compile time**, through generic type parameters that the compiler
specializes, never through a trait object (`Box<dyn …>`): a virtual call plus its indirection on
every read is a cost this path cannot carry. This matches how the existing read-likelihood models
are already swapped — concrete types behind a generic
([`read_model/`](../../../../src/ssr/cohort/read_model/)).

Every algorithm takes a **reusable scratch buffer** for its matrices, owned by the caller and reused
across reads. Allocating a matrix per read is the other cost this path cannot carry.

**Three interfaces, one per shape of transformation.** What separates them is what goes in and what
comes out — nothing about who calls them or why:

| interface | takes | returns |
|---|---|---|
| `AlignmentNormalizer` | an alignment, and the two sequences it relates | that alignment, canonically spelled (§6) |
| `BestPathAligner` | two sequences + the locus context | one alignment (§4) |
| `MarginalAligner` | two sequences + the locus context | one log-probability (§5) |

The locus context is what varies as the caller walks the genome — for a repeat-aware algorithm, where
the flanks end and how the repeat slips. It is empty for the affine algorithms, which need nothing
locus-specific.

The two aligners share a noun deliberately: they are the same recurrence under two reductions — take
the best predecessor, or add up all of them — and naming one of them something other than an aligner
would hide that.

*Signatures below are illustrative. What this section fixes is the contract: the three shapes are
separate interfaces, the marginal is returned in logarithms, the emission model is a component
rather than a variant, dispatch is static, and scratch is caller-owned and reused.*

```rust
/// Line a read up against a reference sequence the single most probable way.
pub trait BestPathAligner {
    /// Reused alignment buffers (matrix, traceback) — allocated once, not per read.
    type Scratch: Default;
    /// What reading the line-up yields — a placement, or a measured repeat span.
    type Output;
    /// What varies per locus — nothing for the affine aligner, the repeat's geometry and
    /// stutter model for the repeat-aware ones.
    type Context;
    fn align(&self, read: &[u8], quality: &[u8], reference: &[u8],
             context: &Self::Context, scratch: &mut Self::Scratch) -> Self::Output;
}

/// The total probability the reference produced this read, summed over every line-up,
/// as a natural logarithm. Negative infinity means no line-up is possible.
pub trait MarginalAligner {
    type Scratch: Default;
    type Context;
    fn marginal_probability(&self, read: &[u8], reference: &[u8],
                            context: &Self::Context, scratch: &mut Self::Scratch) -> LogProb;
}

/// Rewrite an existing alignment into its canonical (left-most) spelling. The read is
/// required: how far a gap may shift is bounded by matching on BOTH sequences.
pub trait AlignmentNormalizer {
    fn normalize(&self, alignment: &mut Alignment, read: &[u8], reference: &[u8]);
}
```

**The marginal is returned as a logarithm, and that is a contract, not an implementation detail.**
Two reasons, and it is worth being exact about which one does the work, because the obvious one is
weaker than it looks.

**A single alignment does not underflow.** A 300-base read at a per-base error of 10⁻³ scores about
0.74; even at an absurd 10⁻¹ it reaches only ~10⁻¹⁴, against a floor near 10⁻³⁰⁸. Underflow would
need thousands of bases. So "a whole read would underflow" is false, and it is not the reason.

**The product across reads does underflow, and it fails invisibly.** A genotype likelihood is a
product over every read at a locus, and that is where the magnitudes collapse. That arithmetic is the
caller's — but handing it a linear probability makes underflow the caller's problem to rediscover at
every call site. **The decisive reason is the failure mode:** in logarithms an impossible line-up
reaches negative infinity, which a caller can see and act on; in linear space it reaches zero, which
is indistinguishable from a value that merely underflowed. One is a fact, the other is a silent loss.

This does **not** force a slow inner loop. Adding probabilities in logarithms is far more expensive
per cell than multiplying them, so an implementation is free to run its matrix in ordinary
probabilities and take a single logarithm at the end (safe at repeat sizes), or to rescale each row
and accumulate the scale factors (the standard trick when the sequences are long enough to
underflow). **The contract fixes the returned value; how the matrix reaches it is the algorithm's
business.** Production's plain-probability result becomes a logarithm at this boundary — a one-line
adapter, and parity against it stays checkable in either space.

**What is fixed at construction, and what travels per call.** The split is *how often the thing
changes*:

- **Fixed when the aligner is built — the emission model.** Per-base quality or a flat error rate.
  It is the experiment's configuration, and it carries the base-error rate, which is a property of
  the sample group rather than of any locus. This is what makes with-quality and without-quality two
  configurations of one algorithm rather than two algorithms.
- **Supplied per call — everything that varies per locus.** For a repeat-aware aligner that is where
  the flanks end, and the stutter model. These change every time the caller moves to the next repeat.

**Why the locus is not constructor state**, which an earlier draft of this spec claimed: holding it
would force a **new aligner per locus** — across millions of loci, and the geometry carries an
allocated motif. It would also be wrong for any implementation whose stutter parameters depend on
the allele being scored, since stutter rate rises with repeat length. Production settled this the
same way: its read-likelihood model is a **stateless value**, and everything that varies arrives
through a per-call context. This design follows that shape rather than inventing another.

A `quality` argument appears on the best-path aligner and not the marginal because production's
marginal scores with a single error rate rather than per-base qualities (§5.1); an implementation
wanting per-base qualities would take them the same way.

---

## 8. Reuse over rewrite — the map to production

| what | existing code | ng reuse |
|---|---|---|
| repeat-aware best-path | `delimit_read`, `ViterbiScratch` ([ssr/pileup/alignment.rs](../../../../src/ssr/pileup/alignment.rs)) | the model for §4.2; its output must additionally report which flank anchored ([`read_preparation_ssr.md`](read_preparation_ssr.md) §3) |
| repeat-aware marginal | `align_subst` ([pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs)) | §5.1 — **this function alone**. **Returns a plain probability; this interface returns its logarithm** (§7) — convert at the boundary |
| the stutter distribution | `HipstrModel` ([ssr/cohort/read_model/](../../../../src/ssr/cohort/read_model/)) | §5.2 — the *genotyping* model, **not an algorithm in this module**. It composes the row above; this module documents its distribution only because the two-penalty best-path aligner (§4.2) needs the same parameters |
| emission models | Dindel per-base-quality table; flat error rate ([alignment.rs](../../../../src/ssr/pileup/alignment.rs), [pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs)) | the swappable emission component (§3) |
| left-alignment | freebayes `LeftAlign.cpp` (vendored, read-only reference) | §6 — the algorithm is standard; implement from the description, not by transliteration |
| affine best-path | — | new (§4.1); optional fast-core swap |
| affine marginal | — | **not built** (§5.3) |

---

## 9. Open questions

- **Is an affine marginal aligner worth building?** (§5.3) Production's general-purpose caller uses
  summed per-base error probabilities over decomposed observations instead. Adding the alignment-based
  marginal is a comparison to run, and the result decides whether the fourth cell of §3's table ever
  gets filled.
- **Does the repeat-aware marginal score the measured repeat, or the whole read?** Production scores
  the already-measured repeat against each candidate allele — cheaper, and the calibration winner in
  production's read-model comparison. Scoring the whole read against a full flank-tract-flank sequence
  in one pass is the alternative (§10, algorithm 6): more principled, and it inherits no delimitation,
  so it does not carry forward a mistake the ruler made.

  **It scores with one flat error rate, not per-base qualities** — following production, which pairs a
  flat rate at scoring with a base-quality gate at read ingestion. That choice is load-bearing beyond
  this module: per-base qualities are per *read*, so keeping them would stop identical-sequence reads
  collapsing into one row with a count, turning the per-locus observation table from a tally into
  per-read records. A flat rate keeps the tally. **A residual cost stays even so:** the whole-read
  design needs more of the read than the repeat, so identical-repeat reads with differing flanks no
  longer collapse as tightly. Build and compare it with fixed synthetic qualities first — that settles
  whether it wins before anyone pays for the storage it implies.
- **One tract penalty or two?** (§4.2) Production charges one flat per-base gap inside the tract;
  the more biologically faithful model prices a whole-unit slip separately from everything else, as
  HipSTR does. Both get built and compared, at **period 1 as well as period 2 and above** — the two
  regimes exercise different parts of the model and neither is skippable (§10.3). The two-penalty
  model must carry the contraction/expansion asymmetry and a working out-of-frame route, or it will
  fail for a reason that is its implementation rather than its premise.
- **Default emission for the repeat-aware best-path aligner.** Production uses per-base quality; a
  flat cost is the comparison. Both ship behind the swappable component; which is the default is
  measured, not assumed.
- **How much headroom above the band's floor?** (§3) The floor — the difference in length between
  read and reference — is forced, and needs no decision. The headroom above it does: it covers a path
  that strays further than the minimum and returns, and nothing derives it. Production's forward
  allows two bases; whether that suits the delimiter is untested. This is the one parameter whose
  wrong value loses long alleles **silently** rather than merely costing accuracy, so it wants a named
  constant carrying its reasoning plus a test at the extremes — never a literal.

---

## 10. The algorithms to build

Eight entries, grouped by the interface each one implements — three normalizers plus five aligners.
All of them get built and then compared; §10.3 says in what order, and why that order is not the list
order.

### 10.1 The list

| # | interface | algorithm | where it comes from |
|---|---|---|---|
| 1a | `AlignmentNormalizer` | one structured pass: merge consecutive indels, carry one across blocks | §6 — production's, a GATK port |
| 1b | `AlignmentNormalizer` | repeated simple passes until nothing moves, capped | §6 — freebayes' |
| 1c | `AlignmentNormalizer` | 1a to a fixpoint, failing loudly if the cap is hit | §6 — a wrapper over 1a |
| 2 | `BestPathAligner`, affine | banded affine matrix, per-base quality | §4.1 — **new code**; production has no general-purpose aligner |
| 3 | `BestPathAligner`, repeat-aware | two gap regimes, flat per-base gap inside the repeat | §4.2 — production's delimiter; **the parity anchor** |
| 4 | `BestPathAligner`, repeat-aware | two penalties: whole-unit slips priced apart from everything else, direction-asymmetric, with an out-of-frame route | §4.2 — HipSTR-shaped |
| 5 | `MarginalAligner`, repeat-aware | one sequence against another, gaps confined to the ends | §5.1 — production's `align_subst`. The **stutter half** of production's read-likelihood model is not part of it: that is the genotyping's (§5.1's boundary note) |
| 6 | `MarginalAligner`, repeat-aware | forward over the whole read against flank-repeat-flank, one flat error rate | §9 — HipSTR/GATK-shaped |

**Not in the set: a forward over a repeat-loop graph.** Encoding the locus as a graph whose repeat
unit has a self-loop, and running a forward pass over it, would yield `P(read | N units)` for
**every** N at once instead of one pass per candidate. It is dropped from the initial set for three
reasons, recorded so it can be picked up if the trigger below fires:

- **It is an optimisation of algorithm 6's regime, not of the incumbent's.** Amortising across
  lengths only pays when each per-candidate score is expensive. Production's current read-likelihood
  model is nearly free — a closed-form stutter term, then algorithm 5 degenerating to a base-by-base
  comparison with no matrix filled at all — so against it the graph pass saves nothing.
- **It does not fit any interface here** (§7). Its input is a locus structure plus a set of unit
  counts, not two sequences; its output is a value per length, not one value. It would need a fourth
  interface, for one algorithm.
- **It has no reference implementation anywhere.** ExpansionHunter uses this graph with a *best-path*
  aligner; a forward over it is an extension, not a port — research risk rather than porting risk.

**Its trigger:** if algorithm 6 wins the marginal comparison, then amortising an expensive per-candidate
forward across lengths becomes worth having, and that is when this earns a description, an interface,
and a place in the list.

### 10.2 Two things that look like algorithms but are compositions

Neither needs its own code, and noticing that removes two entries from the build:

- **The flat-cost repeat aligner** (the GangSTR approach — realign every read against a synthetic
  repeat reference) is **algorithm 2 in its flat-cost configuration, pointed at a
  synthetically-built reference**. The repeat-awareness lives in how the reference is constructed, not
  in the aligner. It comes free once 2 exists, and it is the quality-blind end of the repeat-aware
  comparison.
- **The affine marginal** (§5.3) is **algorithm 2's matrix with the sum reduction instead of the
  maximum** — the same recurrence, the other reduction. Cheap once 2 exists. Whether it is worth
  *using* is a separate question, because its comparator is not another aligner: it is production's
  summed base-quality product, which does no alignment and lives in the caller.

### 10.3 How they are compared, and in what order

**Compared on synthetic reads with known truth**, which is the only setting where the measured repeat
can be scored against the repeat that was actually simulated. Two requirements on the comparison:

- **Measure calibration, not only accuracy.** Production's current marginal was chosen over its
  closest rival *on calibration*, after the two tied on genotype concordance — accuracy alone would
  not have separated them.
- **Compare 3 against 4 at period 1 *and* at period 2 and above — the two regimes test different
  parts of the model.** Period 1 is where the indel deficit lives and where direction asymmetry,
  size decay and placement multiplicity are all live; period 2 and above is the **only** place the
  out-of-frame route is exercised at all. Neither is skippable. At period 1, score indels **of the
  repeat's own base** separately from indels of **a different base**: the first is genuine stutter,
  the second is mis-routed as stutter by the arithmetic in-frame test, and averaged together the two
  effects cancel (§4.2). *(An earlier draft excluded period 1 entirely, on the false premise that the
  two models are identical there.)*
- **The marginal comparison does not fit inside this module.** Algorithms 5 and 6 are not
  head-to-head rivals: 5 scores a measured repeat against a candidate, 6 scores a whole read. The
  comparison that matters — production's current read-likelihood model against a whole-read forward —
  spans this module *and* the genotyping that composes it (§5.1). Its harness therefore lives with
  the genotyping, not here; this module's own tests only establish that each algorithm computes what
  it claims.

**Build order**, which is not the list order:

1. **3, and unbanded.** The parity anchor, and what the STR locus generator is blocked on today.
   Production's delimiter fills its whole matrix (§3), so an **unbanded** port is the only version
   that can be shown byte-identical to it — and that parity is the one hard oracle this module has.
   **Banding follows as its own change**, at the width §3 sets — the per-read length difference plus
   headroom — carrying a test that the output on the parity fixture is unchanged, and a case at the
   long-allele extreme where a too-narrow band would show. A performance commit with its own
   evidence, rather than a behaviour change smuggled into the port.
2. **5** — the marginal that pairs with it, and a small one: it is a single function, most of whose
   work is a base-by-base comparison. Its unequal-length forward must be tested deliberately, because
   the caller it was written for never reaches it (§5.1).
3. **4** — the first genuine comparison, and the one most likely to move indel calling at repeats.
4. **6** — the remaining marginal comparison, built with fixed synthetic qualities so it is settled
   before anyone pays for the storage it would imply (§9).
5. **1a, 1b, 1c** — independent of everything else, and cheap. Build the **property test first**
   (§6): it grades all three against the definition of leftmost, so it is worth more than the
   comparison between them. Then run the differ-at-all screen before anything expensive.
6. **2 last.** The mode that would call it — re-aligning a read whose placement is not trusted — has
   no trigger yet (`read_preparation.md` §4). An aligner built before that decision cannot fire, so
   building it earlier buys nothing.

**How the plans group this.** The build is split into three plans, one per interface, because an
interface plus its implementations is a coherent unit to build and review:

1. [`alignment_best_path.md`](../impl_plan/alignment_best_path.md) — the module skeleton, the shared
   aligner types, then algorithms 3, 4 and (gated) 2. First, because algorithm 3 is the only
   externally-blocking item in the module and the only byte-parity oracle in it.
2. [`alignment_marginal.md`](../impl_plan/alignment_marginal.md) — `LogProb`, the marginal
   interface, algorithms 5 and 6.
3. [`alignment_normalization.md`](../impl_plan/alignment_normalization.md) — the normalizer
   interface and algorithms 1a/1b/1c. Last because nothing depends on it, though it is short and can
   be pulled forward if the indel-placement question becomes urgent.

Grouping this way separates 3 from 5, which the order above pairs. That is harmless: algorithm 4 is
compared against **3**, not against 5, so nothing in plan 1 waits on plan 2.
