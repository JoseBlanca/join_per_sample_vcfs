# SSR allele model — Mark 2 (empirical candidate alleles)

**Status:** draft, 2026-06-17, branch `pileup-offladder`. Built **from scratch**,
one agreed premise at a time (the way the other SSR docs were grown). It exists
because a code-vs-spec divergence surfaced while wiring off-ladder — the as-built
rung builder reconstructs candidates as **pure tilings** (`motif × L`), which does
not reproduce an **impure** reference tract — and probing it exposed a deeper
problem: a **reference-anchored ladder is reference-biased** (if the population
sits off the reference's whole-unit ladder, the whole real ladder reads as
"off-ladder"). The Mark-2 model removes that bias by making candidate alleles
**empirical** (assembled from observed reads, identity by sequence), with the
reference demoted to a **coordinate frame** only.

This **supersedes the Mark-1 design** in [ssr_pileup.md](ssr_pileup.md) (§2
reference-anchored ladder), [ssr_shared_types.md](ssr_shared_types.md) §2 (the
`OnLadder`/`OffLadder` enum), and [ssr_offladder.md](ssr_offladder.md) (off-ladder
wiring) — those are on hold pending this model settling; the corrections feed back
once it does. Spec §4.2/§4.3/§5.1 are likewise owed an amendment.

This doc is written **as the dialogue settles each step** — sections firm up only
once agreed; §7 is the live agenda of open details.

---

## 1. What a reference repeat *is* (agreed 2026-06-17)

A reference repeat **locus** is fundamentally a **location** — `(chrom, start,
end)` — a region of the reference that **TRF flagged as a tandem repeat**. That
location is its identity: it is how a locus is matched across samples in the
cohort coordinate-merge. The catalog additionally embeds the **reference bases**
there (tract + flanks), but only as the **coordinate frame** Stage 1 uses to
delimit a read's repeat region from its flanks (§4) — it asserts nothing about
which alleles are real.

The **motif** and the **repeat count** are *not* part of a locus's identity:

- At **catalog build** the motif's only role is a **period filter** — drop
  homopolymers (period 1) and motifs above the SSR period ceiling (period > 6).
  After that filter, a locus is just a located, TRF-flagged region.
- The **repeat count** is a derived annotation, never identity — the model works
  from observed **sequences**, not counts.
- The motif **returns later** as a Stage-2 stutter covariate (it defines the
  whole-unit step, §6) — a property the *model* uses, still never locus identity.

The Mark-2 consequence: the model never rebuilds a tract from `motif × count` — it
works from observed **sequences** over a located reference frame, exactly what
Mark-1's `Allele::to_sequence` ([types.rs](../../src/ssr/types.rs), rebuilding
rungs from `motif × count`) got wrong.

## 2. The Mark-2 shift, in one statement (agreed 2026-06-17)

> **Candidate alleles are *observed*, not *invented*.** The set of real alleles at
> a locus is assembled from the reads' observed repeat-region sequences (identity
> = byte-equality), not generated as whole-unit steps off the reference. The
> reference's *only* roles are (a) locating the loci and the motif (Stage 0) and
> (b) providing a coordinate frame to delimit each read's repeat region from its
> flanks. It makes **no claim about which alleles are real.**

Consequences that fall straight out:

- **There is no on-ladder / off-ladder label in Stage 1.** Every read contributes
  its observed repeat-region sequence, uniformly. Pure and impure tracts are the
  same kind of thing — a sequence.
- **The on/off-ladder *concept* does not vanish — it relocates to Stage 2's
  stutter model**, as a *relationship between candidate alleles*: "is candidate B
  reachable from candidate A by whole-unit slippage?" (stutter moves in whole
  motif units). That relationship is what the stutter kernel needs, and Stage 2 is
  its right home (§6).
- **The reference bias is gone.** A population shifted off the reference is simply
  a set of observed sequences whose peak isn't the reference — assembled correctly,
  with no allele mislabelled.

---

## 3. Stage 0 — catalog (reference-only, unchanged)

Run TRF on the reference to locate SSR loci; per locus store `(chrom, start, end,
motif, ref_seq + flanks)` exactly as today (spec §3, [ssr_catalog.md](ssr_catalog.md)).
The catalog stays **reference-only** — no cohort dependency, works for a new
species (a stated design goal). It supplies the **motif** and the **coordinate
frame** (flanks for border-anchoring), nothing about real alleles.

---

## 4. Stage 1 — per-sample pileup (delimit + observe)

Per locus, per read — no likelihoods, no candidate scoring here:

1. **Border alignment (delimit — Option A, decided §7-Q2).** Align the *whole read*
   against the single locus reference frame (`left_flank + ref_tract + right_flank`)
   with our pair-HMM in **Viterbi + traceback** mode (unbanded — cheap at locus
   sizes), and read the `left_border | repeat_region | right_border` split off the
   two flank-junction columns. One alignment, no likelihood, no candidate sweep
   (delimit only); the flanks anchor both boundaries jointly and the read's length
   difference sits as an indel in the repeat region. **Determinism** (the
   cross-sample identity guarantee): fixed traceback tie-break
   `Match > Deletion > Insertion`, and a junction indel is assigned to the
   **preceding (5′) block** (HipSTR's rule — so a left-junction indel joins the
   flank, a right-junction indel joins the repeat). A read missing a border (the
   repeat runs off the read end) is **not spanning** → QC only ("allele ≥ read
   length", spec §1.4). *(The repeat region is known only after this step — hence
   the quality gate comes next, not first.)*
2. **Quality gate on the repeat region.** With the repeat region now delimited,
   compute the **first quartile of its base qualities** and drop the read if that
   quartile is below a threshold (assembly-style). Survivors are confident enough to
   be treated as **uniform, good-enough quality** — which is what lets Stage 1 store
   sequences without per-base quals (§7-Q1). *Caveat to measure, not assume:*
   repeat-region quality tends to be lower at the longest alleles, so this is an
   allele-length-dependent dropout — tolerable under precision-first (spec §1.4) but
   to be quantified.
3. **Per-sample empirical distribution.** Tally the surviving reads' repeat-region
   **sequences** in a `HashMap<sequence, count>` (sequence-keyed, so the full
   observed sequence is preserved; the count is the observation weight).
4. **Store to `.ssr.psp`.** Per locus: the distinct observed repeat-region
   **sequences + their counts** (the per-sample observed ladder) + the QC scalars.
   **No quality, no per-read rows**: the gate (step 2) already made the survivors
   uniform-quality, so the signal is the counts (the Stage-2 flat-error HMM supplies
   the error rate, §7-S3). This replaces the Mark-1 `(length, log-lik)` profiles and
   the off-ladder columns with one uniform "observed sequences + counts" column set.

Stage 1 stays **per-sample and embarrassingly parallel** (no cohort context).

---

## 5. Stage 2 — cohort (assemble candidates, then score)

1. **Pool** the per-sample distributions across the cohort (coordinate-merge scan,
   spec §2) → the cohort-aggregate observed distribution over sequences.
2. **Assemble the candidate allele set** from that distribution: the abundant
   sequences are the real-allele candidates, with stutter-aware pruning of the
   satellites (the generalisation of spec §5.1's peak+adjacent rule from the
   integer ladder to sequences; §7-S1) and the hard cap `MAX_CANDIDATE_ALLELES`.
3. **Likelihoods (HMM, here not in Stage 1).** Score each observed sequence against
   each candidate with the pair-HMM under a **flat (uniform-quality) error model**
   → `Qᵣ`. Identical observed sequences collapse across the whole cohort, so this
   is `|distinct sequences| × |candidates|` HMM runs per locus, once (θ-independent),
   reused across EM iterations.
4. **Stutter + EM + genotype + VCF** as spec §5 — over the empirical candidate set
   rather than the reference ladder. The stutter convolution `P(read|a) = Σ
   Qᵣ·S_θ` runs over the stutter-reachable candidates (§6).

Stage 2 is heavier than Mark-1 (it holds observed sequences and runs the HMM) —
**accepted** (2026-06-17). Worth a working-set sanity-check against the
memory/cohort-scaling thesis, but not a blocker.

---

## 6. Where on/off-ladder went

It becomes a **Stage-2 relationship between candidate alleles**, used only by the
stutter kernel: candidate B is a *stutter-reachable* (was: "on-ladder") neighbour
of candidate A iff B is A with whole motif units added/removed and the rest of the
sequence (any interruption) preserved; otherwise B is a *distinct* (was:
"off-ladder") allele, not a stutter product of A. This is exactly the
discriminator the stutter model already needs (spec §5.2: slippage is whole-unit,
in-frame). So the concept is not lost — it is **defined between observed
candidates, not against the reference**, which is what removes the bias. (How to
decide "whole units added/removed with interruption preserved" for two arbitrary
sequences is §7-S2.)

---

## 7. Open details to iron out (live agenda)

Quality / Stage 1:
- **Q1 — the region quality gate.** Method **decided (2026-06-17): the first
  quartile of the repeat-region base qualities, dropped below a threshold.** Open:
  the threshold value (repurpose `MIN_BASE_QUAL`?) and the long-allele-dropout
  measurement.
- **Q2 — border-alignment method + determinism. DECIDED (2026-06-17): Option A.**
  One unbanded global pair-HMM (Viterbi + traceback) of the whole read against the
  single reference frame `left_flank + ref_tract + right_flank`, delimit-only (no
  likelihood, no candidate sweep); boundary = the flank-junction columns of the
  traceback. Determinism (the cross-sample identity guarantee — same molecule →
  same bytes in every sample): fixed tie-break `Match > Deletion > Insertion`,
  junction indels assigned to the **preceding (5′) block** (HipSTR's rule).
  Grounding: HipSTR aligns + tracebacks (per-haplotype); GangSTR falls back to SSW;
  **neither treats determinism as designed — we do.**
- **Q3 — what to store. DECIDED (2026-06-17): distinct `(sequence, count)` per
  locus, no quality.** The gate (§4 step 2) already removed the low-Q reads, so the
  survivors are uniform-quality and the signal lives entirely in the counts.
  Container shape: a CSR of repeat-region byte blobs + a parallel count column, per
  locus, + the QC scalars.

Stage 2 / candidates:
- **S1 — candidate assembly rule.** Generalise peak+adjacent from the integer
  ladder to sequences: which abundant sequences are kept, how stutter satellites
  are pruned without pruning real adjacent alleles, the `MAX_CANDIDATE_ALLELES`
  cap + no-call.
- **S2 — stutter reachability between sequences.** The rule deciding "B = A ± whole
  units, interruption preserved" for two observed sequences (the relocated
  on/off-ladder test, §6) — and how the stutter kernel's `δ` (whole-unit step) is
  defined when an allele is impure.
- **S3 — the flat-error HMM model.** Single per-base error rate: where it comes
  from (the gate threshold? a global estimate?), and whether a pair-HMM is even
  needed vs a simpler affine-gap score now that quality is uniform.

Cross-cutting:
- **X1 — VCF REF/ALT.** REF = reference tract sequence; ALT = observed candidate
  sequences. Repeat-count fields (`REPCN`) derived from sequence length / motif
  (fractional for impure), as spec §5.9.
- **X2 — depth cap / reservoir.** Still per-locus in Stage 1 (determinism gate).
- **X3 — what to retire.** The Mark-1 `Allele::OnLadder/OffLadder` enum,
  `build_rungs`, the reference-ladder window, and the off-ladder columns all go;
  `types.rs`/`candidate_generation.rs`/`registry_ssr.rs`/the spec need amending.
