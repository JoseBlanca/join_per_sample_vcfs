# SSR interrupted-repeat recall — set-valued rungs (same-length alleles)

*Status: proposal, 2026-07-06. Companion to
[`ssr_cohort_mark2.md`](ssr_cohort_mark2.md) — it resolves that spec's one
open S1 threshold (§9, "the same-length cohort frequency threshold") and
brings the `ssr-call` implementation into line with what §5 and §7 already
mandate. Code-layout companion:
[`architecture/ssr_interrupted_repeat_recall.md`](../architecture/ssr_interrupted_repeat_recall.md);
build order:
[`implementation_plans/ssr_interrupted_repeat_recall.md`](../implementation_plans/ssr_interrupted_repeat_recall.md).*

**Goal: match HipSTR's correctness on interrupted repeats, in two phases.**
Phase 1 recovers the dropped loci by making alleles sequence-keyed (the
structural fix — no new statistics). Phase 2 closes the remaining modelling gap
against HipSTR by giving each allele its own stutter rate — estimated by
covariate-informed hierarchical shrinkage (purity one covariate among length,
period, and motif), so the rate responds to the allele and not just its length.
Phase 1 is a conformance fix; Phase 2 is the refinement Phase 1 makes
measurable. The two are sequenced deliberately (§8): Phase 1 produces the called
impure alleles and the per-allele slip counts that Phase 2 fits on.

## 1. The problem, in one sentence

Our SSR caller silently drops real length-variable loci whose polymorphism is
an **interior interruption** — a base that breaks the repeat, like the `TC` in
`TGTGTG`**`TC`**`TGTGTG` — because the genotyper collapses every tract of the
same length into a single allele, so an interruption that changes *composition*
without changing *length* becomes invisible.

## 2. What an interruption is, and why it hides

A microsatellite (**SSR** — short tandem repeat) is a motif repeated head to
tail: `(TA)n`, `(TTA)n`. An **interruption** is a base inside the tract that
does not fit the motif — `TATATA`**`C`**`TATATA` — so the tract is *impure*.
Interruptions are common and heritable, so a locus can be polymorphic *for the
interruption itself*: some samples carry the pure repeat, others the
interrupted version, at the **same tract length**.

That "same length" is the whole problem. Our Stage-2 genotyper (`ssr-call`)
organises evidence on a **length ladder**: it bins every observed tract by its
length in motif units, and — in the current code — treats one length as one
allele. Two sequences of equal length but different interior therefore land in
the same bin and are genotyped as identical. The interruption polymorphism is
averaged away, the locus looks monomorphic, and the emit rule ("emit only if
variable") drops it with no VCF row at all.

This contradicts the Mark-2 spec, which already says the opposite (§4).

## 3. Evidence — the drop, localized

On the `ssr_tomato1` benchmark (51 tomato samples, our caller vs HipSTR), of
HipSTR's 1,578 length-variable loci, 765 have no overlapping row in our VCF.
About a quarter carry an interior interruption. Six were traced end to end.

**The drop is in Stage-2 (`ssr-call`), not Stage-1 (`ssr-pileup`).** All six
loci are present in the per-sample pileup files with deep coverage; all six
reach the genotyper's emit gate as *admitted, non-variable* and are dropped.
The clearest case, at `SL4.0ch01:4279322` (motif `TTA`), pooled across the
cohort:

| observed tract sequence | length | reads | samples |
|---|---|---|---|
| `TTATTATTATTGTTATTA` (interrupted, `TTG`) | 18 bp | 235 | 50 |
| `TTATTATTATTATTATTA` (interruption resolved to pure `TTA`) | 18 bp | 26 | **17** |
| length-changing stutter variants (13/15/16 bp) | — | 1 each | 1 each |

Seventeen of fifty samples carry the pure-repeat haplotype — a common
polymorphism HipSTR calls without trouble (34/50 non-reference there). Both
real alleles are **18 bp**, so our length ladder cannot tell them apart: the
locus is length-monomorphic and gets dropped. The only length-changing reads
are singletons — sequencing noise, correctly ignored.

(A second, distinct failure mode also drops interruption loci: a genuine
*length* alt exists but is carried by only a handful of samples and is
suppressed at low per-sample depth under the cohort's high apparent inbreeding.
That is a threshold/prior issue, not a structural one; it overlaps the separate
short-period effort and is **out of scope here** — see §11.)

## 4. What the spec already says (and the code does not do)

This proposal is a conformance fix because Mark-2 already specifies the target
behaviour. Two settled passages:

**§5 (candidate assembly).** *"A rung is length-keyed and holds a set of
sequences — every distinct same-length sequence whose cohort frequency clears a
threshold (often one sequence, sometimes several: interruption / substitution
variants). Each such sequence is an independent allele … The below-threshold
same-length cloud is per-base sequencing error — not promoted to candidates."*

**§7 (stutter).** *"Stutter is a per-allele operation, not membership of a
global length lattice … Each candidate — pure or impure — anchors its own
stutter ladder … Impure alleles are treated identically to pure ones."*

So the design of record is: **a candidate allele is a sequence, not a length; a
rung holds a set of them; interruptions are first-class.** The implementation
took a documented v1 shortcut and stopped at one sequence per rung. Four
places carry that shortcut:

1. **Candidate nomination** — [`candidate_set.rs`](../../../src/ssr/cohort/candidate_set.rs),
   `cohort_representative` returns only the single most-supported sequence at a
   nominated length. This is the shortcut that discards the interruption allele.
2. **EM seeding** — [`em_init.rs`](../../../src/ssr/cohort/em_init.rs),
   `candidate_of_length` maps a nominated length to the *first* candidate at
   that length, so a second same-length candidate gets no seed mass.
3. **Per-locus stutter refit** — [`em.rs`](../../../src/ssr/cohort/em.rs),
   `attribute_locus` / `nearest_parent` assign each read to a called allele *by
   length*, so two same-length called alleles are indistinguishable to the
   refit.
4. **Allele-balance FP term** — [`vcf_out.rs`](../../../src/ssr/cohort/vcf_out.rs),
   `allele_balance` short-circuits to "balanced" when the two called alleles have
   equal `genotype_units` (it reads a same-length het as a homozygote) and, behind
   that, attributes reads to alleles *by length* too. So the false-heterozygote
   defence Mark-2 §6 relies on is **silently disabled for exactly the same-length
   het this proposal adds** — the one axis where `G₀` and recurrence give no
   separation, so balance carries the most load. It must instead trigger on allele
   *identity* and split reads by *composition*, consuming the E-step's own
   per-candidate `Qᵣ` (Mark-2 §6 amendment, 2026-07-06). This is the same
   sequence-aware attribution as item 3, applied to the FP term.

The remaining machinery is already in conformance, which is what keeps this fix
small:

- The **read likelihood** already scores same-length sequences differently.
  `Qᵣ(obs | candidate)` in
  [`read_model/hipstr.rs`](../../../src/ssr/cohort/read_model/hipstr.rs) reaches
  its slip variants **from the candidate's own sequence** (per-allele, Mark-2
  §7), and
  for a zero-length-change read it scores the base composition through the
  substitution term `align_subst(obs, cand, ε)`. Two 18 bp candidates that
  differ at the interruption base get different likelihoods already.
- The **EM and genotype enumeration** already index by candidate, not by
  length: genotypes are unordered pairs of candidate indices, π is
  per-candidate. Two candidates that share a length are already legal downstream.
- The **frequency prior** `G₀` gives same-length candidates equal weight (both
  sit at the modal length), which is correct — it is a length-distance prior and
  has nothing to say about the composition axis.

## 5. Phase 1 — recover interrupted loci (sequence-keyed alleles)

Phase 1 is purely structural: make the allele a *sequence*, not a length, so
the machinery that already scores sequences per-allele (§4) actually receives
the same-length candidates it needs. No new statistics, no new parameters that
need fitting. This alone stops the drops.

### 5.1 Rungs become set-valued; nomination emits sequences

The rung ladder already stores the distinct sequences per length
(`seqs_by_length` in
[`rung_ladder.rs`](../../../src/ssr/cohort/rung_ladder.rs)); the fix uses them.
Replace the "one representative per length" step with "every same-length
sequence that clears the admission threshold (§5.2) becomes its own candidate."
Concretely, `cohort_representative(length) -> Option<sequence>` becomes
`cohort_alleles(length) -> Vec<sequence>`, and nomination unions the results as
it already unions across samples. The reference allele stays seeded
unconditionally. `MAX_CANDIDATE_ALLELES` still caps the union; a locus that
exceeds it is no-called as before.

On a clean (pure) locus every length holds exactly one sequence, so this is a
no-op — which is why the 96.7 % genotype concordance on already-called loci is
expected to hold.

Two placement points, settled:
- **The threshold filters *nomination*, not rung storage.** `build_rungs` keeps
  every observed sequence; only candidate promotion applies §5.2. The
  below-threshold error cloud still reaches the likelihood regardless, because
  the read model scores each sample's *raw observed sequences*, not the rung's
  stored set — so filtering the nomination source loses nothing.
- **Same-length candidates are ordered deterministically** by sequence bytes, so
  the ALT column order and the cross-thread byte-identity contract both hold.

### 5.2 The same-length admission threshold (resolving Mark-2 §9)

This is the one genuinely new decision. A same-length sequence must clear a bar
that a real interruption allele passes and a sequencing-error variant does not.
The discriminator is **recurrence**, exactly as it is for admitting a rung
(Mark-2 §5): a real interruption-resolved haplotype recurs across its carriers,
whereas a substitution error is sporadic and lands at a different base each
time. A same-length sequence is promoted to an independent candidate when it
clears **all** of:

- **cohort read count** ≥ `min_same_length_reads` — enough total support to be
  more than noise;
- **sample recurrence** ≥ `min_same_length_samples` — seen in at least this
  many distinct samples (the sporadic-vs-systematic test);
- **within-length fraction** ≥ `min_same_length_fraction` — at least this share
  of the reads at its length, so a single deeply-sequenced sample cannot alone
  manufacture an allele.

All three are new `CandidateCfg` fields with conservative defaults to confirm on
the benchmark (starting proposal: 8 reads, 3 samples, 0.10). The
`4279322` pure allele (26 reads, 17 samples, 26/262 ≈ 0.10 of the 18 bp reads)
clears them comfortably; the singleton noise variants fail on every count.

**Admission is recurrence-based only — no single-sample "height" path (settled).**
Mark-2 §5 also admits a rung by *height* (a local maximum in one sample), which
rescues rare private alleles. That path does **not** translate to same-length
variants: two alleles at one length cannot form a local maximum on the length
ladder, so there is no height signal to key on — recurrence is the only
discriminator available. The accepted consequence: a genuine **private
(one-sample) same-length allele** — e.g. a lone sample homozygous for an
interruption nobody else carries — is **not admitted**, even with deep clean
support. This is a deliberate precision-first limit: from a single sample, a
same-length substitution is indistinguishable from a paralog or a
sample-specific systematic artifact, so we decline to call it rather than risk
the false positive.

Sequences below the bar are **not discarded** — per Mark-2 §5 they re-enter the
likelihood as the substitution-error products of the admitted candidates, so
nothing is lost, it simply is not called.

### 5.3 Seeding and stutter refit become sequence-aware

- **Seeding.** `candidate_of_length` must return the candidate whose *sequence*
  a sample's peak actually matches, not the first candidate at that length. When
  several candidates share a length, pick the one the sample's representative
  (most-supported) sequence at that length equals, falling back to the closest by
  the substitution score the likelihood uses. For a **homozygote** this seeds the
  right same-length allele. For a **same-length heterozygote it cannot seed a het
  at all**: the sample shows a single length peak with one representative
  sequence, so the putative-genotype tally still assigns *both* allele copies to
  the majority composition (a hom seed) and the minority same-length allele starts
  only at its `G₀` pseudocount floor — accepted, because recovering the het is the
  EM's job, not the seed's (below), exactly as the seed deliberately mislabels a
  length-merged het and lets the self-correcting E-step fix it (Mark-2 §4.3).
- **Stutter refit.** `nearest_parent` must attribute a read to the nearest
  *called allele*, breaking a same-length tie by composition (the read scores
  higher against the allele it matches). This keeps the per-locus stutter shape
  `θ_locus` estimated from the reads that genuinely belong to each allele.

The tie-break "by composition" **reuses the read model's `align_subst`** (exact
match, else its substitution score) — the same metric the likelihood uses — so
attribution and scoring never disagree and the choice stays deterministic.

Both are local changes: they replace a length lookup with a sequence match over
the small candidate set.

**Same-length het recovery — the EM, not the seed (a bounded limitation).**
Because the seed cannot express a same-length het, a het carrier starts from a hom
seed and its minority allele from the pseudocount floor; the call is recovered
only if the E-step moves posterior weight onto the het genotype. It can — the read
likelihood is composition-aware (§4), so genotype `(A,B)` outscores `(A,A)` once a
sample carries clear minority-composition reads, and the floor keeps `π_B > 0` so
the correction is never trapped at zero (Mark-2 §4.3). Recovery is strong at
moderate depth but **weakens at the low-depth / high-`F_IS` tail** — a handful of
minority-composition reads against an inbreeding-tilted, homozygosity-favouring
prior with a small cohort `π_B`. This is the **same** suppression as the length-alt
failure mode (§3, scoped out in §11), here applying to same-length hets, which are
**in scope**. Consequence: Phase 1 reliably **emits** an interruption locus (some
carriers are same-length homozygotes), but a fraction of its **het** carriers at
low depth may still be called homozygous — a genotype-concordance effect, not a
drop. It is *measured, not assumed* (§9). If the benchmark shows material het
undercall, the mitigation is a per-sample same-length het seed proposal (split a
length peak whose reads divide across two §5.2-admitted sequences), which converges
with the D1 two-allele test and is a follow-up, not Phase 1.

### 5.4 VCF representation — same-length alleles follow the HipSTR form (settled)

Two same-length alleles are distinguished in the output the way **HipSTR** does
it: by the **REF/ALT sequence** plus the `GT` allele index — *not* by a
repeat-count field. Our VCF already writes actual tract sequences in REF/ALT
(Mark-2 §4.5), so same-length ALTs emit with **no format change**.

Our VCF's own length annotation is **`REPCN`** (repeat copy number per called
allele); we do **not** emit `BPDIFFS`. `REPCN` is the direct analogue of HipSTR's
length fields — `BPDIFFS` (INFO) and `GB` (per-sample FORMAT), both a pure
`len(allele) − len(REF)` (verified in
[`HipSTR/src/seq_stutter_genotyper.cpp`](../../../HipSTR/src/seq_stutter_genotyper.cpp)).
All of them **collapse** for a same-length allele: it has a 0 bp difference from
REF and prints the same copy number as its same-length sibling. This is not a
defect — it is exactly what HipSTR emits: its record for `4279322` carries three
ALTs with `BPDIFFS = -3, 0, 0`, the two `0` entries being same-length,
sequence-distinct alleles it tells apart only by the REF/ALT string + `GT`.
GangSTR, by contrast, keys an allele *by* repeat count and renders each ALT as a
pure motif tiling, so it **cannot represent a same-length variant at all** — which
is why matching HipSTR (the tool this benchmark targets) is the right choice.

**Follow HipSTR: keep `REPCN`, document the ambiguity — do not remove it (settled).**
No length annotation — `REPCN`, `BPDIFFS`, or `GB` — can distinguish two alleles of
the same length; the distinction lives only in the sequence, which REF/ALT + `GT`
already carry unambiguously. HipSTR keeps its ambiguous length fields rather than
drop them, and `REPCN` is correct and useful on the pure-repeat majority, so we do
the same: keep it, and state in the header/spec that it is a length annotation
**ambiguous across same-length alleles**. Removing it would lose a standard,
human- and tool-facing field to "fix" an ambiguity that removal does not resolve
(the sequence is the only carrier either way). No new field is added in v1;
downstream consumers must read identity from `GT` + REF/ALT, not from `REPCN`
(see §9 for the concordance-metric consequence).

### 5.5 A prerequisite, not a payoff — Phase 1 unblocks a later pre-pass fix but does not deliver it

The same length-collapse blindsides the **pre-pass that fits the shared stutter
and error parameters** (Mark-2 §4.3/§4.4), one layer before genotyping. The
pre-pass learns the per-base error `ε` and the stutter shape/level from confident
genotypes, and its gate today is a **length-histogram heuristic**
([`resolve_confident_genotype`](../../../src/ssr/cohort/rung_ladder.rs)): it finds
a clear single peak in a sample's *length* histogram, calls it a homozygote for
the peak's dominant sequence, and treats every other read as either that allele's
stutter (a length slip) or its sequencing-error halo (same length, different base).

A sample that is actually **heterozygous for two same-length alleles** —
pure-18 / interrupted-18 — shows one length peak, so the length-based gate
mislabels it a homozygote for the majority composition. Its second allele's
on-length reads sit at the peak length but differ in sequence, so they are
miscounted as the error halo and **inflate `ε`** (each differing base scored as a
substitution — a bias the pre-pass code already flags in a `BIAS NOTE`). The
stutter **level is *not* pushed up** by this: those same-length reads are counted
as *faithful* (Δ = 0), not as slips, so if anything they mildly dilute the slip
fraction — the contamination is essentially an `ε` inflation, concentrated on
interruption-rich cohorts.

**The fix is a separate change this proposal does not make.** Catching it needs
the pre-pass gate to become sequence-aware — the Mark-2 **model-based
one-allele-vs-two-allele test** (the 1-vs-2-peak BIC resolution test, Mark-2 §4.3 /
roadmap **D1**). D1's milestone shipped, but with the length-histogram heuristic
above as a stand-in; the **BIC form — the part that would fix this — is not yet
built** (the code marks it as the intended layer), and it is **not in Phase 1's
touch points** (§7). Phase 1 is a genuine **prerequisite** for it — the two-allele
test can only propose "two alleles at one length" once the candidate machinery will
nominate two same-length sequences, which is Phase 1's whole point — but **landing
Phase 1 alone changes nothing in the length-based gate**, so the `ε` contamination
**persists until the BIC-form gate lands separately**. (An earlier draft called
cleaner `ε`/`θ` a *free payoff* of Phase 1; that was wrong — Phase 1 *unlocks* the
fix, the later gate change *delivers* it.)

One consequence to keep in view: because the contamination persists, the frozen
`ε` that **Phase 1's own read model consumes** may be mildly inflated on
interruption-rich data. This does **not** block Phase 1 — the same-length
discrimination has margin to spare (a single interruption base costs a read a
factor `(ε/3)/(1−ε)` against the wrong composition, ≳ Phred 20 per read even at an
inflated `ε ≈ 0.02`, §4) — but it makes the BIC-form gate a worthwhile **follow-up**
for parameter calibration on interruption-rich cohorts, not a Phase 1 deliverable.

### 5.6 What Phase 1 leaves — one stutter rate shared across an allele set

Phase 1 makes interrupted alleles callable, but it still scores every allele at a
locus with the **same stutter rate** (Mark-2's per-locus `θ_locus`, adjusted only
by length); pure and impure alleles are separated by composition (`ε`) alone.
That recovers the loci, but it is not as correct as HipSTR, which gives each
haplotype its own stutter behaviour. And the gap is not only about purity: two
alleles at one locus can genuinely differ in how they slip for reasons a single
per-locus number cannot express. Closing that is Phase 2.

## 6. Phase 2 — per-allele stutter (HipSTR-parity)

### 6.1 The goal

Stutter is a property of the **specific allele**, not just of its length. An
interruption stabilises a tract — the polymerase slips less on `(TA)₃C(TA)₃` than
on a clean `(TA)₇` — but purity is not the only driver: motif composition
(AT-rich tracts slip more) and length move the rate too, and two alleles at one
locus can differ for reasons none of those name. HipSTR gives each haplotype its
own stutter behaviour. Phase 2 does the same — estimated in a way that survives
the fact that most alleles are seen in only a handful of reads.

### 6.2 The obstacle — you cannot fit a free stutter model per allele

The tempting move — give every distinct allele its own free stutter parameters —
fails on data. Stutter is learned from the **skirt** of reads around a
confidently-called allele, and a rare allele contributes only a few reads
cohort-wide; there is no signal to fit even a one-parameter rate per sequence,
let alone the three-parameter shape. This is why the Mark-2 pre-pass pools all
confident genotypes of a period before fitting (Mark-2 §4.3). Crucially, **data
starvation is not special to impure alleles** — a rare *pure* allele is just as
thin — so the cure has to be general, not an impurity special case.

### 6.3 The model — a covariate-informed default, per-allele deviation earned from data

The general cure is **hierarchical partial pooling**: each allele has its own
stutter rate, shrunk toward a default by how much data that allele carries. A
data-rich allele keeps its own estimate; a thin one falls back to the default.
The one decision that makes this *correct* rather than merely general is **what
the default is**.

A plain pooled default — "the average rate over all alleles of this period" —
is wrong for exactly the case shrinkage is meant to save. That pool is dominated
by **pure** alleles (they are the majority), so a thin impure allele shrinks
toward a pure-like, **too-high** rate — precisely when its own data cannot pull
it back. The fix is to make the default a **prediction from the known drivers of
stutter**, not a blind average:

> `default_rate(allele) = f(length, period, motif, purity)`

fit as a small regression on the data-rich pooled skirts, with each allele's rate
shrunk toward its own covariate prediction:

> `rate(allele) = shrink( this allele's skirt evidence, toward f(covariates), by its read count )`

Now a thin impure allele **borrows strength from every other impure allele** in
the cohort (through the purity term in `f`), so its default is correctly lower
than the pure-dominated pool; a data-rich allele can still **deviate** from what
its covariates predict, capturing idiosyncrasy the covariates miss. **Purity is
one covariate among several, not a privileged truth** — which matches the fact
that two alleles are not more alike just because both happen to be impure. Its
coefficient is only *identifiable*, though, if it survives the two checks in §6.5
(that purity is separable from length, and which stutter parameter it moves).

This **provisionally** keeps the stutter **shape** (the up/down/decay profile)
pooled — there is rarely enough per-allele data to move a three-parameter shape —
and lets the covariates and shrinkage act on the **level** (the overall slip
rate), where the bulk of the signal is *expected*. Whether that is the right lever
is **not assumed**: §6.5(b) flags that an interruption may act on the slip-size
*decay* rather than the level, and the §8-step-5 measurement decides where the
per-allele term attaches. The read model already computes the level per (sample,
candidate) in its scoring loop, so applying a per-allele level is one extra factor
at an existing site; a pure allele at an average locus lands on today's arithmetic.

### 6.4 Fitting it in the pre-pass — why it must follow Phase 1

`f` is fit where the shared parameters are already fit: the confident-genotype
**pre-pass** (Mark-2 §4.3/§4.4), whose skirts are the raw slip evidence. Two
things make Phase 1 a hard prerequisite:

- The **per-allele slip counts** come from Phase 1's sequence-aware attribution
  (§5.3) — each read attributed to its parent allele *as a sequence*, its length
  offset the slip. Without it there are no per-allele skirts to regress.
- The **purity covariate** can only be fit once impure alleles are being called,
  which is Phase 1's whole point. Before Phase 1 there are zero called impure
  alleles, so `f` could carry only a literature guess for purity.

Guard rails, so nothing misfires on thin data. Shrinkage toward the covariate
default *is* the guard by construction: an allele with little data becomes its
prediction, not a noisy fit. With no impure alleles anywhere the purity term
never activates and the model collapses to Phase 1 identically. And whether a
per-allele deviation term **beyond** the covariates earns its place is itself
measurement-gated (§8): if length + period + motif + purity explain nearly all
the variation, the per-allele residual rarely moves and we do not pay for it —
build the covariate-informed default first, add the per-allele deviation only if
the data show residual variation the covariates miss.

### 6.5 Two model-form questions to settle by measurement, before writing types

The §6.3 form — purity as one covariate, acting on the level — rests on two
assumptions that Phase 1's own data can and must check first (this is what the
§8-step-5 gate is for). Settle both before writing Phase 2 types; each has a clean
test and a defined fallback.

**(a) Is purity separable from length?** Interruptions accumulate with tract
length, so purity and length **co-vary** across a catalog — a plain
`rate ~ length + purity` regression cannot cleanly attribute variance to one or the
other (the coefficients trade off; a purity term can be a length effect in
disguise). The identifying variation is exactly the one Phase 1 newly supplies:
**same-length pure-vs-impure allele pairs** — the interruption polymorphism itself,
two alleles at one length differing only in purity. Test purity **at fixed length**
on those pairs, not on the length-confounded marginal. If impure alleles slip less
than pure alleles *of the same length*, purity is a real, separable effect and
earns its covariate; if not, drop it and let length + period + motif carry the
default. So the same-length contrast is not only what Phase 1 recovers — it is also
the natural experiment that makes the purity coefficient identifiable.

**(b) Does purity move the *level* or the *decay*?** An interruption breaks the
uninterrupted run the polymerase slips along, and slip *magnitude* grows with run
length — so an interruption should suppress **multi-unit** slips more than ±1
slips: it steepens the geometric **decay** (a *shape* effect), rather than merely
scaling the overall rate (a *level* effect). HipSTR carries this in its per-locus
in-frame geometric; §6.3's "level, not shape" lever cannot represent "impure
alleles still slip, but only by ±1." Test: bin each allele's slip counts by slip
**size** (|Δ|) and compare the *distribution* pure-vs-impure, not just the total
slip fraction. If the **decay** differs, the purity term must attach to the
shape/decay (or to both level and decay); if only the total differs, level alone is
right. Resolve this before choosing where the per-allele purity term lives
(`param_estimation.rs` and the read-model `level` vs `decay` site, §7).

Neither question blocks Phase 1; both are Phase 2 design inputs. The measurement
gate already stops Phase 2 if purity moves nothing (§8 step 5) — these two checks
**extend** it from "does purity matter?" to "*how* does purity enter?", so the
model form is a data decision, not a default carried into code.

## 7. Touch points

**Phase 1** (structural — sequence-keyed alleles):

| File | Change |
|---|---|
| `src/ssr/cohort/rung_ladder.rs` | expose the set of admitted sequences per length (`cohort_alleles`); keep `cohort_support` for `G₀` |
| `src/ssr/cohort/candidate_set.rs` | nominate every same-length sequence clearing §5.2; new `CandidateCfg` thresholds |
| `src/ssr/cohort/em_init.rs` | `candidate_of_length` → sequence match, not first-at-length |
| `src/ssr/cohort/em.rs` | `attribute_locus` / `nearest_parent` → attribute reads to the nearest called allele **by sequence**, same-length tie broken by composition (the stutter refit; effect-neutral for same-length alleles but shares the shortcut) |
| `src/ssr/cohort/vcf_out.rs` | **allele-balance made sequence-aware** — trigger the balance check on distinct candidate *indices*, not equal `genotype_units`, and split reads by composition via the E-step's per-candidate `Qᵣ` (Mark-2 §6 amendment); the length-keyed short-circuit currently disables the FP defence for same-length hets. ALT emission needs no change (already candidate sequences, §5.4); confirm multiple same-length ALTs format correctly and that `REPCN` (our only length annotation — no `BPDIFFS`, §5.4) prints without error for same-length ALTs |

`G₀` ([`allele_freq_prior.rs`](../../../src/ssr/cohort/allele_freq_prior.rs))
and the genotype EM need no change in Phase 1 (§4).

**Phase 2** (statistical — per-allele stutter rate):

| File | Change |
|---|---|
| `src/ssr/cohort/param_estimation.rs` | fit the covariate default `f(length, period, motif, purity)` on the pooled skirts (purity kept only if it survives §6.5a); a per-allele purity measure (interruption count / impure-base fraction); the per-allele shrinkage toward `f` of whichever stutter parameter §6.5b selects — the level, the decay, or both |
| `src/ssr/cohort/read_model/` | apply the per-allele term at the existing per-candidate site — to `level`, or the `decay`/shape term if §6.5b selects it |
| `src/ssr/cohort/em.rs` | feed each called allele's slip counts (from §5.3 attribution) into the per-allele rate estimator; shrink toward the covariate default on thin data |

## 8. Implementation plan

Two phases, sequenced because Phase 1 is the measurement instrument for Phase 2
(§6.3). Each phase ships and validates on its own; the caller is correct and
better after Phase 1, and HipSTR-parity after Phase 2. Work each phase
incrementally — settle the types, then the implementation, then validate — and
pause between phases for review.

### Phase 1 — sequence-keyed alleles

1. **Types first.** Turn the per-length representative into a per-length *set*:
   `rung_ladder.rs` gains `cohort_alleles(length) -> Vec<sequence>`;
   `CandidateCfg` gains the three §5.2 thresholds. Settle these signatures
   before touching the callers.
2. **Nomination.** `candidate_set.rs` unions the admitted same-length sequences
   (§5.1–§5.2). Unit-test the new admission with a pure-vs-interrupted fixture.
3. **Seed + attribution.** Make `candidate_of_length` (seed) and `nearest_parent`
   (refit) sequence-aware (§5.3).
4. **Validate** (§9): unit fixture, then the `ssr_tomato1` benchmark
   + dashboard, then the threshold sweep. **Exit criterion:** interrupted-repeat
   drops recover as variable loci **and** the 847 already-common PASS loci hold
   ≈96.7 % genotype concordance.

### Phase 2 — per-allele stutter rate

Only after Phase 1 lands and its benchmark is clean.

5. **Measure first.** With Phase 1 calling impure alleles, read out the observed
   per-allele slip rates (Phase 1's attribution already records them) and settle
   the §6.5 model-form questions before writing any Phase 2 type: (i) **does purity
   lower the rate at *fixed length*** — tested on Phase 1's same-length
   pure-vs-impure pairs, so the purity effect is identified free of the length
   confound (§6.5a); (ii) **does purity move the overall level or the slip-size
   *decay*** — tested by comparing the |Δ| slip-size distribution pure-vs-impure,
   not just the total (§6.5b); (iii) is there per-allele variation left after
   length + period + motif + purity? **Gate:** include a covariate only where it
   moves the data (drop purity if it fails (i)); attach the purity term to the
   level or the decay per (ii); add the per-allele deviation term only if a residual
   survives (iii). If nothing moves the rate, Phase 2 is unnecessary and we stop —
   the honest outcome, not a failure.
6. **Types first.** A per-allele purity measure + the covariate default
   `f(length, period, motif, purity)` and a per-allele rate field in the stutter
   parameters (`param_estimation.rs`), shrinking toward `f`.
7. **Estimate + apply.** Fit `f` on the pooled skirts and the per-allele
   deviation in the pre-pass (§6.4); apply the per-allele rate to `level` in the
   read model (§6.3).
8. **Validate.** Re-run the benchmark. **Exit criterion:** genotype-quality
   calibration on impure loci improves (or matches HipSTR) with no regression on
   pure loci and no loss of the Phase 1 recall gain.

## 9. Validation plan

Applied per phase (§8):

1. **Unit** — a cohort fixture with a same-length interruption polymorphism
   (pure vs interrupted at one length) must yield two candidates and a variable
   call; the pre-fix code yields one candidate and a dropped locus. **The fixture
   must include a genuine same-length *heterozygote* carrier (pure/interrupted at
   one length) at realistic low depth (~2–5 reads), and the test must assert that
   carrier is called *het* — not merely that the locus is emitted.** The seed
   cannot represent a same-length het (§5.3), so this is the case that exercises
   the EM's recovery path; a fixture of same-length *homozygotes* alone would pass
   without ever testing it. Phase 2 adds a fixture where the impure allele's
   known-lower slip rate must be recovered.
2. **Benchmark** — re-run `ssr-call` on `ssr_tomato1` and the ours-vs-HipSTR
   dashboard. Phase 1 success: the interruption drops (the ~23 % of 765 that
   carry an interior interruption) recover as emitted variable loci, **and** the
   847 already-common PASS loci keep ≈96.7 % genotype concordance. **Report
   per-sample *het* concordance on the recovered interruption loci, not only locus
   recovery** — the het-undercall tail (§5.3) shows up there, not in the emit
   count. **Concordance must be scored on the *sequence* genotype, not the length
   genotype:** a same-length interruption het is a HET by sequence but collapses to
   a HOMOZYGOTE by length, so a length-based metric would mis-score a *correct*
   recovery as a regression. Concretely — if comparing through TRtools, read
   compareSTR's `metric-conc-seq` (not `metric-conc-len`) and pass statSTR
   `--use-length` off (`GetLengthGenotypes` / the default `uselength=True` collapse
   same-length alleles; `GetStringGenotypes` / mergeSTR keep them); the ours-vs-HipSTR
   dashboard must likewise key genotypes on `GT` + REF/ALT sequence, not `REPCN`
   (§5.4). Phase 2 success: improved calibration on impure loci, no pure-locus
   regression. Report both numbers each time — recovered loci and concordance
   movement.
3. **Threshold sweep** — the §5.2 defaults (8 reads / 3 samples / 0.10) are a
   starting point; sweep them against recovered-loci vs new-false-positive count
   and pick the knee, the same way the length-alt thresholds were set.

## 10. Risks and how the design bounds them

- **False heterozygotes from substitution error promoted to an allele** (Phase
  1). The §5.2 recurrence bar is the primary guard (sporadic error fails
  recurrence). Population recurrence also drives a non-recurrent candidate's
  frequency toward zero. Because `G₀` cannot penalise a same-length variant (both
  sit at the modal length, §4), the §5.2 threshold carries more of the load here
  than it does for length alts — hence three joined conditions, not one. **The
  allele-balance term (Mark-2 §6) is *not* a working backstop here until §7's
  vcf_out change lands:** as written it is length-keyed and disabled for exactly
  the same-length het (item 4 above), so this proposal must make it sequence-aware
  or the depth-driven-het defence does not apply to same-length alleles at all.
- **Systematic same-length artifacts read as real interruption alleles** (Phase
  1). Recurrence separates a real allele from *random* error, but not from a
  *systematic* same-length variant — a collapsed paralog / CNV copy, a consistent
  mismap, or systematic damage carrying a fixed substitution at the modal length
  recurs across carriers and clears §5.2. Neither `G₀` (no length separation) nor
  sequence-aware allele balance (a paralog sits near a clean 1/ploidy) rejects it,
  so Phase 1 can convert a currently-silent drop into an emitted false het in every
  carrier. This is the sharpest new FP mode; it needs a coverage-excess /
  cohort-heterozygosity guard adapted from the SNP paralog filter, tracked in
  `doc/devel/TODO.txt`. Until that guard exists, gate Phase 1 on a benchmark FP
  audit that flags cohort-universal same-length hets.
- **Same-length het undercall at low depth** (Phase 1). The seed cannot express a
  same-length het, so recovery rests entirely on the EM; at the low-depth /
  high-`F_IS` tail a het carrier can stay called homozygous (§5.3). Bounded by the
  pseudocount floor + the composition-aware likelihood, and *measured* by the §9
  per-sample het-concordance check. The seed-level mitigation (a same-length het
  proposal, converging with the D1 two-allele test) is added only if the benchmark
  shows material undercall — it is not Phase 1. This is a concordance effect, not a
  drop: the locus is still emitted.
- **Over-fitting the per-allele rate on thin data** (Phase 2). Guarded by
  shrinkage toward the covariate default and the §6.4 measurement gate: a thin
  allele becomes its `f(covariates)` prediction, not a noisy fit; the per-allele
  deviation term is added only if a residual survives the covariates.
- **Regression on already-called loci.** On pure loci Phase 1 is a no-op (one
  sequence per length) and Phase 2 lands a pure allele at an average locus on
  today's arithmetic (the covariate default reduces to the current shared rate
  there). The regression gate is the `ssr_tomato1` concordance (must hold
  ≈96.7 %) plus the SNP caller's end-to-end tests (Mark-2's stated regression
  gate — pre-alpha, no byte-identity).
- **Candidate-set blow-up on messy loci** (Phase 1). `MAX_CANDIDATE_ALLELES`
  already caps the union and no-calls the overflow; set-valued rungs raise the
  count only where genuine same-length structure exists, and the §5.2 bar
  filters the error cloud before the cap sees it.

## 11. Out of scope (named, not hidden)

- **Rare *length*-alt suppression** — a second failure mode behind the 765 drops
  (§3): a real *length* alt carried by few samples, suppressed at low per-sample
  depth under a high apparent `F_IS`. This is a prior/threshold question that
  overlaps the separate short-period/dinucleotide effort; it needs no
  candidate-set change and is tracked there. (The **same-length het** analogue of
  this suppression is **in scope** here — the seed cannot express a same-length
  het, so it is documented and validated in §5.3/§9/§10; only the *length*-alt
  version is deferred to that effort.)
- **Tract-boundary disagreements** with HipSTR (our catalog trims to the pure
  core, HipSTR anchors on an impure prefix) — a catalog-coordinate matter, not a
  genotyping one; it shows up as position-shifted calls, not drops.
- **Other HipSTR-vs-ours modelling gaps** unrelated to interruptions — the
  out-of-frame (single-base) stutter term we fold into `ε`, and HipSTR's
  quality-weighted alignment emission — are deliberate Mark-2 simplifications
  tracked in that spec (§6/§9), not part of this interruption work.
