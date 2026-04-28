# Phase chains — what they are, what they answer, why we use them

**Status:** Explanatory document, 2026-04-28. Companion to
[per_sample_caller.md](per_sample_caller.md) and
[calling_pipeline_architecture.md](calling_pipeline_architecture.md).
This is not a design spec — those already commit to the phase-chain
mechanism. This document exists to make the *concept* concrete:
what a phase chain represents, why the project needs one, how it
relates to the two existing mainstream approaches (GATK and
freebayes), and what shape Stage 1 ultimately writes to disk. A
reader new to the project should be able to read this end-to-end
and walk away with a clean mental model.

The order:

1. The problem phasing is supposed to solve.
2. How GATK solves it (assembly-mediated, with phase output).
3. How freebayes solves it (per-read linkage internally, no phase output).
4. Why we adopt neither approach as-is.
5. What our "phase chain" is, in precise terms.
6. The slot encoding, with a worked example.
7. A side-by-side comparison.
8. What downstream stages actually do with the chains.

A single toy example threads through sections 1, 6, and 8 so the
abstract definitions never lose their grounding.

## 1. The problem: per-position calling cannot tell compound haplotypes from independent variants

A variant caller emits one record per reference position (or one
record per anchor for indels). The biological reality, however, is
per-haplotype: each chromosome is a single connected sequence, and
a sample's two haplotypes (for a diploid) carry their respective
sets of variants. Two variants near each other can be **on the
same haplotype** (a compound allele) or **on different haplotypes**
(two independent heterozygous variants), and the right model for
the sample's evidence is different in each case.

### A small worked example

Say the reference is `AGTC` at positions 100-103 and the sample's
reads, after alignment, look like this:

```
position:    100   101   102   103
reference:    A     G     T     C

read R1:      A     G     T     C    (REF everywhere)
read R2:      A     G     T     C    (REF everywhere)
read R3:      A     T     -     -    (SNP G→T at 101, DEL covering 102-103)
read R4:      A     T     -     -    (SNP G→T at 101, DEL covering 102-103)
```

A naive per-position view sees:

- At 101: alleles `G` (from R1, R2) and `T` (from R3, R4) — looks
  heterozygous.
- At 102 (deletion anchored at 101 by VCF convention; the
  deletion's evidence sits there, but the per-position record at
  102 still shows): allele `T` (from R1, R2) and "no
  observation" from R3, R4 — looks heterozygous for a deletion.

A naive joint caller would emit two separate heterozygous
variants. But that is not what is actually in this sample.
R3 and R4 show that the SNP at 101 and the deletion at 102 came
from the same physical molecule — i.e. they sit on the same
haplotype. The other haplotype is plain reference. So this sample
is heterozygous for a single **compound allele** that combines
"SNP + DEL" on one chromosome, not heterozygous for two
independent variants.

The compound model and the independent-variants model differ in
their predicted read counts, their predicted likelihoods, and
ultimately their genotype probabilities. The joint stage needs
some way to tell them apart.

### What information is actually needed

The fact that distinguishes the two models is one specific
question, asked one position-pair at a time:

> *In this sample, do the alleles at p₁ and p₂ co-occur on the
> same physical DNA molecule, or do they come from different
> molecules?*

The phase information the caller needs is exactly this. Not "which
chromosome (maternal or paternal)?", not "what is the global
haplotype across the genome?" — just per-pair co-occurrence on a
shared molecule.

Anything that answers this question per pair is sufficient. The
question of *how* to record that information is what the rest of
this document is about.

## 2. How GATK answers the question: assembly-mediated phasing

GATK's HaplotypeCaller solves it by reconstructing the haplotypes
themselves and then reading the linkage off the reconstruction. It
does not record per-read linkage as a primary data structure.

The flow ([AssemblyBasedCallerUtils.java:728-900](../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java#L728-L900)
in the local GATK clone):

1. **Active region detection.** GATK windows the genome into
   regions where variants are likely. Most of the genome is
   inactive and is skipped for phasing purposes.
2. **Local re-assembly.** Inside each active region, GATK assembles
   candidate haplotypes from the reads using a de Bruijn graph
   approach. This produces a small set of plausible long sequences
   spanning the entire region (typically 1-2 kb).
3. **EventMap per haplotype.** Each assembled haplotype is
   compared to the reference and decorated with an `EventMap` —
   the list of variant events (SNPs, indels) that haplotype
   carries.
4. **Read scoring.** Every read in the region is scored against
   every assembled haplotype via PairHMM, producing a per-read
   per-haplotype likelihood. Critically, GATK stores these
   likelihoods in `AlleleLikelihoods<sample, allele, read>` — a
   3-D matrix with no cross-position linkage stored per read.
5. **Variant calling per position.** Each position's genotype is
   called by marginalising over haplotypes.
6. **Phasing post-processing.** This is the step that produces the
   PID/PGT/PS tags in the output VCF
   ([AssemblyBasedCallerUtils.java:790-900](../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java#L790-L900)):

   - For each called variant, collect the set of assembled
     haplotypes whose `EventMap` contains it.
   - For each pair of called variants in the region, compare their
     haplotype sets:
     - **Identical sets** (variants always co-occur on the same
       haplotypes) → phase them together with PGT `0|1`.
     - **Mutually exclusive sets** (variants partition the
       haplotypes) → opposite phase, PGT `1|0`.
     - Otherwise unphased.
   - Connected components in this pairwise relation become phase
     blocks. Each phase block gets a single PID equal to
     `pos_ref_alt` of its first variant
     ([AssemblyBasedCallerUtils.java:972](../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java#L972)).
   - PS is the numeric position of the same first variant
     ([AssemblyBasedCallerUtils.java:997](../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/AssemblyBasedCallerUtils.java#L997)).

### Applied to the toy example

Inside the active region containing positions 100-103:

- Assembly produces two haplotypes:
  - `H_REF`: `AGTC`. EventMap = {} (no variants).
  - `H_ALT`: `AT--` (SNP at 101 + DEL spanning 102-103). EventMap
    = {SNP@101 G→T, DEL@101 starting one base after the anchor}.
- Reads R1, R2 score best against `H_REF`; R3, R4 best against
  `H_ALT`.
- Phasing step:
  - SNP at 101 appears on `{H_ALT}`. DEL at 101+ appears on
    `{H_ALT}`. Same set → phased together.
  - PID = `101_G_T`, PS = 101, PGT = `0|1`.

This works **even when no single read spans both events** — the
linkage came from assembly, not from individual reads. That is
GATK's strength.

### The cost

Assembly + PairHMM are expensive: every read against every
candidate haplotype, with a numerical inner loop. Worse, assembly
quality depends on coverage. With 30× reads and a clean active
region, GATK reconstructs cleanly. With 2-10× reads — the
project's target — the de Bruijn graph is sparse, candidate
haplotypes are unreliable, and the resulting "phasing" is at best
the same as direct read evidence and at worst noise from spurious
assemblies.

## 3. How freebayes answers the question: per-read linkage in adaptive windows, used internally only

Freebayes takes a third route, distinct from GATK's full assembly
and from the proposal in this document. It records per-read allele
linkage as a primary internal data structure, uses it during
likelihood evaluation, and then discards it before writing the
output. Its windows are dynamic and locus-focused rather than
fixed-size active regions.

The flow ([AlleleParser.cpp:3214-3324](../../freebayes/src/AlleleParser.cpp#L3214-L3324)
and the surrounding code, plus the supporting headers):

1. **Per-read allele extraction.** When a BAM alignment is parsed,
   every allele observation it contributes is bundled into a
   `RegisteredAlignment` ([AlleleParser.h:47-80](../../freebayes/src/AlleleParser.h#L47-L80)).
   The struct carries the read's name, alignment span, and a
   vector of all `Allele` objects extracted from that read. Each
   `Allele` also carries a `readID` backpointer
   ([Allele.h:123](../../freebayes/src/Allele.h#L123)). This is
   the per-read linkage: reads that contribute to multiple
   positions have all their alleles in one place, indexed by read.
2. **Adaptive haplotype window.** At each candidate locus,
   freebayes computes a window that grows from 1 bp until it
   accommodates the longest observed allele plus any
   indel-in-tandem-repeat structure ([AlleleParser.cpp:3223-3310](../../freebayes/src/AlleleParser.cpp#L3223-L3310)).
   The window is *not* a fixed active region; it adapts to local
   complexity. The `--haplotype-length` /
   `--max-complex-gap` parameter (default 3 bp,
   [Parameters.cpp:205](../../freebayes/src/Parameters.cpp#L205))
   controls how close mismatches must be to be merged into a
   single complex allele.
3. **Allele enumeration.** All alleles within the window — across
   reads — are collected ([AlleleParser.cpp:3282](../../freebayes/src/AlleleParser.cpp#L3282)),
   grouped by sequence equivalence, and enumerated as candidate
   haplotype-alleles ([AlleleParser.cpp:3285](../../freebayes/src/AlleleParser.cpp#L3285)).
   No de Bruijn assembly: candidates are observed allele
   combinations only, never invented ones.
4. **Likelihood and genotyping.** The Bayesian model evaluates
   each candidate genotype using the per-read evidence assembled
   in step 1.
5. **Output.** The best genotype per sample per position is
   emitted to VCF with GT, GQ, depth, and allele counts. **No
   PID, no PGT, no PS, no phase information about which variants
   share a haplotype** appears in the FORMAT or INFO fields
   ([ResultData.cpp:8-637](../../freebayes/src/ResultData.cpp#L8-L637),
   exhaustive scan).

### Two distinct senses of "phase" inside freebayes

Freebayes uses the word "phase" in two different ways and only one
of them is what this document is about:

- **Within a single complex allele.** When two adjacent SNPs sit on
  the same read and within the `--max-complex-gap` window,
  freebayes merges them into a single complex MNP/indel allele
  (e.g. ref `AC` → alt `GT` represented as one allele, not two
  phased SNPs). Inside that single allele, the two changes are
  implicitly phased because they share the allele string. This is
  what
  [Parameters.cpp:44-48](../../freebayes/src/Parameters.cpp#L44-L48)
  calls "phased in this way." It is allele construction, not
  cross-record phase emission, and our pipeline does the same
  thing via the architecture doc's "extend the anchor REF" rule.
- **Across records.** None. Each VCF position is emitted
  independently; the per-read linkage in `RegisteredAlignment` is
  internal scaffolding for the likelihood calculation and is
  dissolved once the per-position genotype is called.

So a SNP at position 100 and an indel at position 200 supported
by the same reads are correctly used during freebayes' likelihood
calculation, but the output VCF gives no signal that those
variants co-occurred on shared molecules. A downstream consumer
sees only the per-position genotype calls.

### Applied to the toy example

Same example as §1: SNP at 101 and DEL anchored at 101 spanning
102-103.

- Reads R3 and R4 each produce one `RegisteredAlignment` whose
  `alleles` vector contains both the SNP and the DEL.
- The haplotype window at position 101 expands to cover the
  deletion (the DEL allele's reference-side length forces it).
- The candidate alleles enumerated in the window include the
  compound `T-` allele (SNP+DEL together).
- The likelihood is computed using R3 and R4 as supporting
  evidence for the compound allele, R1 and R2 for the reference
  allele.
- The output VCF emits one record at position 101 with the
  compound allele as a single ALT entry. No PID is set; no
  separate phase relationship is exposed for downstream consumers
  to consult.

The within-record handling here is essentially what our pipeline
also does at this single anchor record. The difference shows up
only when a question crosses record boundaries — which is exactly
the question Stage 5 has to ask, and which freebayes' output does
not carry the answer to.

### Why per-read linkage is internal-only in freebayes

Freebayes is a single-shot caller: BAMs in, VCF out, in one
process invocation. There is no separation between per-sample
evidence collection and joint posterior computation, so the
per-read linkage is generated, consumed, and discarded inside the
same run. The output VCF format does not need a per-read linkage
field because freebayes already used the linkage to compute the
posterior.

This is the structural reason freebayes has no need for something
like a phase chain id in its output. Our pipeline, in contrast,
splits the work into a per-sample artefact (the `.psf`) and a
later joint stage that can be re-run cheaply across many cohort
calls. The per-read linkage that freebayes computes-and-throws-away
has to survive in the `.psf` for Stage 5 to consult — which is the
job phase chain ids exist to do.

## 4. Why we adopt none of the obvious alternatives as-is

Three reasonable alternatives almost work, and none is taken as-is.
The first two are the existing reference designs (GATK and
freebayes); the third is a tempting shortcut at the cohort level
that does not actually substitute for read-level evidence.

### Why not GATK's assembly-mediated approach

Already in
[calling_pipeline_architecture.md §"Per-read likelihood quality"](calling_pipeline_architecture.md)
as the rejection of Option C, summarised here:

- **Coverage.** At 2-10× per sample, a single sample's reads are
  not enough to reassemble haplotypes reliably. With 3-5 reads
  covering a region, most "reassembled" candidate haplotypes are
  artefacts.
- **Reassembly is really a joint operation.** It pays off when
  reads from *many samples* can be pooled — that does not fit a
  per-sample Stage 1 topology, where each sample is processed
  independently.
- **Implementation cost.** De Bruijn graph + PairHMM is a
  multi-month engineering effort with significant correctness
  surface.
- **Speed.** Per-sample work would be 10-30× slower than the
  freebayes-shaped likelihood the project picked.

The architectural payoff of GATK's approach — phasing across gaps
that no single read spans — is real but expensive, and the
project's coverage target is in the regime where the payoff is
unreliable anyway.

### Why not freebayes' approach as-is

Freebayes' per-read linkage is structurally very close to what we
want: it records, per read, every allele the read contributes to.
Our chains are essentially the same idea. But freebayes' approach
is missing the one thing this pipeline absolutely needs:
**persistence**.

Two specific gaps:

- **The linkage is not exported.** Freebayes' `RegisteredAlignment`
  lives only in memory during one invocation. Re-running the
  cohort with a new sample added forces freebayes to re-process
  every BAM from scratch. The architecture doc lists "read each
  BAM exactly once" and "cohort re-callable without re-doing
  per-sample work" as binding constraints; reusing freebayes
  unmodified breaks both.
- **The linkage is locus-scoped, not read-scoped.** Freebayes uses
  the linkage only inside the immediate haplotype window when
  evaluating a candidate genotype. It does not keep linkage
  information alive across record boundaries, so even if we were
  willing to re-run freebayes on every cohort change, there is no
  cheap way to carry the linkage out.

So the project takes a hybrid: keep freebayes' direct per-read
linkage idea, but make it persistent and read-scoped (not
window-scoped) so Stage 5 can ask cross-position questions later
without re-touching any BAM. That is what §5 defines.

### Why not lean on cohort haplotype frequencies instead?

A tempting shortcut at the cohort level: with thousands of
samples, some individuals will be homozygous at every constituent
position of a candidate compound. Their genotypes pin down the
haplotypes that exist in the population. For heterozygous
individuals, why not assume their two haplotypes are drawn from
that observed homozygote pool — using cohort haplotype frequencies
to *infer* compound calls instead of carrying per-read chains
through the pipeline?

This is structurally what classical population-phasing tools
(fastPHASE, BEAGLE, IMPUTE, SHAPEIT) do for whole-genome haplotype
reconstruction. It is real signal, and at thousands of samples it
can be very strong. But it cannot replace the chain mechanism for
the specific job the chain does, for one fundamental and three
practical reasons.

**The fundamental reason: it confuses likelihood with prior.**
Per-sample inference combines, in Bayes' rule, a *per-sample
likelihood* `P(reads | genotype)` with a *cohort-level prior*
`P(genotype)`. These enter at different points and cannot
substitute for each other.

The chain's job is to make the *likelihood* correct for compound
genotypes. Computing `P(reads | compound)` requires knowing, **per
read**, whether that read shows alt at every constituent position
or only at some — because a read that shows the compound directly
contributes very different evidence from a read that shows only
one constituent. Cohort haplotype frequencies are a population
fact: they tell you how often a given haplotype exists across
samples. They do not tell you, for sample s, which of s's reads
support the compound. Substituting them into the per-read
likelihood is a category error that biases every per-sample call
in proportion to how often the per-read independence assumption is
violated.

**Three practical consequences:**

- **Per-sample likelihood becomes systematically biased.** Without
  the chain, the only available per-sample likelihood treats
  constituent variants as independent. For a sample heterozygous
  at two compound constituents, the independent model overstates
  support for the "ref-alt + alt-ref" haplotype pair when the
  truth is "compound + ref" (or vice versa). That bias propagates
  into the genotype calls fed back to cohort-frequency estimation
  — a corrupted base for the iteration.

- **Bootstrapping is circular and HWE-fragile.** Cohort haplotype
  frequencies depend on per-sample genotype calls, which depend
  on per-sample likelihoods, which (without the chain) are wrong
  in the way above. EM can iterate, but only under Hardy-Weinberg.
  HWE breaks under inbreeding, admixture, selection at the locus,
  population structure, and LD with selected loci — exactly the
  conditions under which interesting compound mutations often
  appear.

- **Rare compounds vanish entirely.** Hom-alt frequency under HWE
  is `q²`. With compound frequency `q = 0.05` and N = 1000,
  expected hom-alts ≈ 2.5; with `q = 0.01`, expected hom-alts
  ≈ 0.1. A compound with no homozygous anchor in the cohort
  cannot be inferred from cohort frequencies, regardless of N.
  Per-read chains, by contrast, work for a single sample of any
  rarity.

So cohort haplotype frequencies are not statistically substitutable
for chains in the per-sample likelihood. They *are* a useful
**complementary** signal — adding them as a prior in Stage 6
strengthens the phase-broken fallback path, where the chain
returns no information for some samples (long-range compounds,
mate-missing pairs, low-coverage constituents). That extension
belongs in the posterior engine, not in Stage 1, and is described
in
[calling_pipeline_architecture.md §"Stage 6 — posterior engine"](calling_pipeline_architecture.md).
The chain mechanism stays load-bearing for everything below the
prior layer.

## 5. The phase chain: read-level direct linkage

A **phase chain** is the set of allele observations contributed by
**one read** (or one read pair, when both mates pass the filter)
across all the reference positions it touches.

Three pieces of that definition are load-bearing:

- **Per read or per pair.** One physical DNA molecule, one chain.
  When two paired-end mates of the same fragment overlap or sit
  near each other on the genome, they are evidence about the same
  molecule and share a chain.
- **Across positions.** A chain spans every reference position the
  read or pair touches with `M`/`=`/`X` CIGAR ops, plus any
  anchored indels. Insertions and deletions belong to the same
  chain as the surrounding bases.
- **Allele observations, not haplotypes.** A chain is a *bag of
  observations*, not an assembled sequence. Each observation is a
  (position, allele) pair tagged with the chain. The chain says
  "these observations all came from the same molecule"; it does
  not say "this is what the molecule's full sequence is."

### What the chain answers

Exactly the question identified in §1:

> *In this sample, do the alleles at p₁ and p₂ co-occur on the
> same physical DNA molecule?*

If both observations carry the same chain identifier, yes. If
they carry different identifiers, no. If one of them does not
appear at all (the chain has no observation at that position),
the chain doesn't link them — there is no evidence either way.

### What the chain does *not* provide

- **No maternal/paternal assignment.** A chain says two
  observations came from one molecule, not which of the two
  parental chromosomes that molecule is. PGT-style `0|1` vs
  `1|0` is *not* something Stage 1 can produce; that is a
  population-genetic or trio-based decision that lives elsewhere.
- **No statistical / population phasing.** There is no LD-based
  reasoning, no reference-panel imputation, no read-aware phaser
  like WhatsHap or SHAPEIT. Those add long-range linkage from
  outside the reads; chains add nothing beyond what the reads
  themselves directly observed.
- **No linkage across non-overlapping reads.** When a read ends
  at p₂ and the next read starts at p₃ > p₂ with no overlap, no
  chain links them. (For paired-end fragments the gap between the
  two mates is bridged by the shared fragment identifier, but only
  within that one pair.)
- **No haplotype reconstruction.** Stage 1 never tries to figure
  out the consensus sequence of the molecule a chain came from.
  The five per-allele scalars are aggregated on top of the chain
  membership; the chain itself is just an identifier.
- **No "the read carries the variant" claim.** A chain at slot k
  means a read contributed observations at every position the
  chain appears in; one of those observations might be the alt
  allele, another might be the ref. Stage 5 looks at the
  per-position allele to decide.

### Chain length in practice

Bounded by the read or pair span:

- Single short read: ~150 bp → chain spans up to ~150 reference
  positions (minus any deletions, which the chain still covers
  via their anchor).
- Paired-end pair, typical 200-500 bp fragment: chain spans the
  union of the two mates' alignments, often with a small gap in
  the unaligned middle.
- Paired-end pair, 1 kb fragment: chain reaches up to ~1 kb,
  again with a middle gap.

So chains are *short-lived*. At any given reference position, the
number of chains "currently active" is bounded by per-position
depth times a small factor (chain spans ÷ position spacing), which
is to say: roughly per-position depth itself. At 30×, that is ~30;
at 2-10× (the project target), 5-15. This bound is what makes the
slot encoding cheap, and is the subject of the next section.

## 6. The slot encoding, with a worked example

Stage 1 assigns each chain a **slot id** — a small integer (`u8`
in practice; `u16` if we ever want to be defensive). Slots are
recycled: when a chain ends and a new one starts, the new chain
can occupy a freed slot. To prevent consumers from confusing the
old chain with the new one in a recycled slot, every per-position
record carries two **lifecycle markers**:

- `new_chains`: slot ids that started since the previous emitted
  position (a chain begins when its read enters the active set).
- `expired_chains`: slot ids that ended since the previous
  emitted position (a chain ends when its read or pair has fully
  left the active set).

The conceptual model is the same as VCF's `|` / `/` phase markers:
within an open phase block, the slot id is the linkage; the
markers declare where blocks open and close.

### Per-record shape

```
record at position p:
    delta_pos:        distance to previous record
    new_chains:       [slot id, ...]
    expired_chains:   [slot id, ...]
    alleles:
      [
        { sequence: "A",
          num_obs: 5,
          q_sum: ...,
          fwd: 3, placed_left: 2, placed_start: 1,
          chain_slots: [0, 2] },          // contributed by chains 0 and 2
        { sequence: "T",
          ...
          chain_slots: [1] },
      ]
```

`chain_slots` per allele is the small list of slots that
contributed observations of that allele at this position. It is
compact in practice (typically 1-3 entries; bounded by depth).

### A trace through the toy example

Take the same sample as in §1, but flesh out the read positions:

```
read R1:  starts at 100, ends at 101  (covers 100, 101)
read R2:  starts at 100, ends at 103  (covers 100, 101, 102, 103)
read R3:  starts at 101, ends at 103  (covers 101, 102, 103) — SNP+DEL on the same molecule
read R4:  starts at 101, ends at 103  (covers 101, 102, 103) — SNP+DEL on the same molecule
read R5:  starts at 102, ends at 103  (covers 102, 103) — REF
```

Stage 1 walks the merge in coordinate order. Slot allocation:

```
position 100:
  R1 enters active set → assigned slot 0
  R2 enters active set → assigned slot 1
  emitted record:
    new_chains:     [0, 1]
    expired_chains: []
    alleles:
      { "A", num_obs=2, chain_slots=[0, 1] }      # both reads see REF "A"

position 101:
  R3 enters active set → assigned slot 2
  R4 enters active set → assigned slot 3
  R1, R2 still active.
  emitted record:
    new_chains:     [2, 3]
    expired_chains: []
    alleles:
      { "G", num_obs=2, chain_slots=[0, 1] },     # R1 and R2 saw REF "G"
      { "T", num_obs=2, chain_slots=[2, 3] }      # R3 and R4 saw the SNP

  Note: the deletion that R3 and R4 carry is anchored at 101 by
  VCF convention, so the same record gets a third allele:
      { "TGT" → "T" (DEL), num_obs=2, chain_slots=[2, 3] }
  (i.e. the SNP and DEL allele at 101 are *separate* alleles in
  the same record, both supported by chains 2 and 3.)

  R1's last covered position is 101; R1 expires *after* this
  record is emitted (the active set walker only knows R1 has
  ended once it advances past 101).

position 102:
  R5 enters active set → assigned slot 4 (slot 0 is freed by R1's
                                          expiration but we don't
                                          reuse it within the same
                                          tick for clarity; the
                                          algorithm could also
                                          recycle slot 0 here)
  R1 leaves active set (already past its end).
  R2, R3, R4 still active.
  emitted record:
    new_chains:     [4]
    expired_chains: [0]
    alleles:
      { "T", num_obs=2, chain_slots=[1, 4] }      # R2 and R5 saw REF "T"
                                                  # R3, R4 do not contribute here:
                                                  # the deletion already accounted for them
                                                  # at the anchor position 101.

position 103:
  R2 leaves active set after this record.
  R3, R4, R5 still active here.
  emitted record:
    new_chains:     []
    expired_chains: []                            # nothing has *just* expired yet
    alleles:
      { "C", num_obs=3, chain_slots=[1, 4, ?] }
      # R5 contributes; R3 and R4 do not (still inside their deletion span).
```

The crucial observation is in the record at position 101:

- The SNP allele `T` carries `chain_slots = [2, 3]`.
- The DEL allele `TGT → T` carries `chain_slots = [2, 3]`.
- Both alleles share the same chains.

That sharing is the answer to the per-position-pair question: for
this sample, the SNP and the DEL co-occur on chains 2 and 3,
i.e., on two distinct physical molecules. Stage 5 reads this
directly — it does not need to assemble haplotypes.

### Why the markers are needed

Slot 0 is freed when R1 expires and could be reused later — say a
new read R6 starts at position 200 and gets assigned slot 0.
Without lifecycle markers, a downstream consumer scanning records
might see "slot 0 at position 100" and "slot 0 at position 200"
and incorrectly conclude they came from the same molecule. The
`expired_chains: [0]` marker between them tells the consumer:
"slot 0's previous chain has ended; any future appearance of
slot 0 is a new chain." That cleanly separates the two
incarnations of slot 0.

This is exactly the role `|` plays in VCF: within a phase block,
positions linked by `|` share the genotype-phase-set ID; the
opening/closing of phase blocks demarcates which `|`-linked
positions are in the same phase set.

## 7. GATK, freebayes, and our slot id, side by side

| Axis | GATK PID/PGT | freebayes (internal only) | Our slot id |
|---|---|---|---|
| Mechanism | Local re-assembly → haplotype EventMap → set comparison | Adaptive haplotype window + observed-allele enumeration; per-read linkage in `RegisteredAlignment` | Direct per-read tagging at scan time |
| Range | Active region (typically 1-2 kb) | Per-locus haplotype window (typically a handful of bp; longer at indel-in-repeat sites) | Read or read-pair span (~150 bp – ~1 kb) |
| Linkage source | "Variants on the same assembled haplotype" | "Alleles co-observed on the same read within the window" (per-read backpointer `readID` on each `Allele`) | "Observations from the same physical molecule" |
| Requires direct read overlap | No (assembly can bridge gaps) | Yes (only what reads see in the window) | Yes (chain ends when the read ends) |
| Tells you which chromosome (mat/pat) | Yes (`0\|1` vs `1\|0` via PGT) | No | No |
| Persistence in output | Yes (PID/PGT/PS in the VCF) | No (used only during calling, then discarded) | Yes (slot ids + lifecycle markers in the `.psf`) |
| Cost per locus | Heavy: PairHMM × #reads × #haplotypes | Light-medium: enumeration + per-window likelihood | Light: small slot table, integer comparisons |
| Behaviour at 2-10× coverage | Unreliable; assembly under-powered | Works | Works |
| Implementation size | ~thousands of lines (de Bruijn + PairHMM) | ~hundreds of lines (the relevant linkage parts) | ~hundreds of lines (active set + slot allocator) |
| Stored on disk per record | One PID + one PGT per phased genotype | — (linkage never reaches disk) | Two short slot lists (`new_chains`, `expired_chains`) plus a per-allele slot list |
| Downstream consumer's work | Read PID/PGT directly per genotype | N/A — consumer never sees the linkage | Walk records; intersect slot sets between two records |

The three approaches solve overlapping problems with very different
trade profiles:

- **GATK's** strength is reach beyond direct read overlap. Its cost
  is everything that goes with assembly: implementation surface,
  CPU per locus, and unreliability at low coverage.
- **Freebayes'** strength is its directness — per-read linkage as a
  primary internal data structure — and its low cost. Its
  limitation is that the linkage is for internal use during a
  single calling pass and is never made available to a separate
  joint stage.
- **Our slot id** keeps freebayes' direct per-read approach but
  makes it *persistent* and read-scoped (not window-scoped) so the
  linkage survives in the `.psf` for Stage 5 to consult, without
  paying anything close to assembly cost.

For the project's coverage target, scale, and split-stage
architecture, the slot mechanism takes the right pieces from both
references.

## 8. What Stage 5 does with the chains

The compound-haplotype check
([calling_pipeline_architecture.md §"Compound haplotype alleles: the phase-chain check"](calling_pipeline_architecture.md))
is exactly the operation §1 set up: given two specific (position,
allele) observations, ask whether they co-occur on the same
molecule in this sample, and use the answer to choose between the
compound-allele likelihood and the independent-variants
likelihood.

### When the chain is consulted (and when it isn't)

It is worth being precise about the scope, because the chain
mechanism is narrower than it looks. The chain is read by exactly
one operation in the whole pipeline: Stage 5's compound-haplotype
check inside an `OverlappingVariantGroup`. Outside that one
operation, the chain is never consulted.

**When the chain is *informative*.** When the sample is
**heterozygous at two or more constituent positions** of a candidate
compound allele. This is the case where "compound on one haplotype"
versus "two independent variants on opposite haplotypes" is
genuinely ambiguous from per-position observations alone, and where
the chain's answer materially changes the likelihood. The toy
example in §1 hits exactly this case: het at 101 (`G`/`T`) and het
at the deletion anchored at 101.

**When the chain is *degenerate but cheap*.** The check still runs
in less ambiguous configurations, and just confirms whatever the
per-position observations already imply:

| Sample's per-position state | What the chain check does |
|---|---|
| Het at both constituent positions | Most informative: "compound on same molecule?" answered directly |
| Hom-alt at both | Chain confirms (in a normal diploid, both haplotypes carry both events); no shift in posterior |
| Het at one, hom-ref at the other | The compound is structurally impossible for this sample; the check just confirms an empty intersection |
| Hom-ref at both | The compound's likelihood for this sample does not depend on the chain (no alt reads to link); skipped |

**When the chain does *nothing*.**

- Variants that are **not in the same `OverlappingVariantGroup`**.
  No compound, no question.
- **Single-position variants standing alone** in their group. No
  compound, no question.
- Variants in the same group but **far apart, beyond read/pair
  span** (typically because the group was extended by a long
  deletion). The chain check returns an empty intersection in every
  sample; Stage 5 falls back to the phase-broken approximation
  everywhere and flags `CA=1` on the record. The chain is
  *consulted* in this case, but the answer is "no" universally and
  adds no new information.

**Why this scope matters.** Phase chains are deliberately limited
to read/pair span (~150 bp – ~1 kb) because that is the size of
the compound haplotypes the joint stage will ever construct. We
do not need genome-wide phasing — that would be a different
problem (population genetics) requiring different machinery
(statistical phasing, reference panels, trios). For the project's
actual question — *"in this sample, do the alleles at p₁ and p₂
co-occur on the same physical molecule?"* — the narrow rule above
is the entire scope.

### How Stage 5 walks the records

Sequentially, maintaining a tiny per-sample active-slot view:

1. On each per-position record:
   - Apply `expired_chains`: remove those slot ids from the
     active view.
   - Apply `new_chains`: add those slot ids to the active view.
   - For each allele, read its `chain_slots` list to know which
     active chains support that allele at this position.

2. To answer "do alleles X at p₁ and Y at p₂ co-occur in sample s
   on the same molecule":

   - Walk to record at p₁; read `chain_slots` for allele X →
     get the set S₁.
   - Walk to record at p₂; read `chain_slots` for allele Y →
     get the set S₂.
   - Compute |S₁ ∩ S₂|. That is the number of physical
     molecules in this sample that support both X and Y.

3. Use the count to choose the right likelihood model. A
   positive intersection is direct evidence of the compound
   haplotype; an empty intersection means the sample's reads
   only support the constituents separately, in which case
   Stage 5 falls back to a phase-broken approximation and
   flags the record (the `CA` per-sample FORMAT field in the
   final VCF — see the architecture doc).

### Applied to the toy example

Recall: at position 101, chains 2 and 3 supported both the SNP
allele `T` and the DEL allele. Suppose Stage 5 wants to know: does
sample s have the compound allele "SNP at 101 + DEL spanning 102"?

- Read at p=101 (SNP allele "T"): chain_slots = [2, 3].
- Read at p=101 (DEL allele "TGT→T"): chain_slots = [2, 3].
  (Both events live at the same anchor record; the question is
  cross-allele within one record, but the same intersection
  applies.)
- |{2, 3} ∩ {2, 3}| = 2.

So this sample has 2 reads that support the compound allele on
the same molecule. Stage 5 evaluates the compound-allele
likelihood directly, with that observation count feeding the
five-scalar reconstruction described in
[freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md).

### What happens when the chains *don't* line up

If sample t had R3 and R4 replaced by:

```
read R3': starts at 101, ends at 101 (SNP only, mate not seen)
read R4': starts at 102, ends at 103 (DEL only, mate not seen)
```

then in t:

- At 101 SNP allele `T`: chain_slots = [slot for R3'].
- At 101 DEL allele "TGT→T": chain_slots = [slot for R4'].
- Intersection is empty.

Stage 5 sees that neither read supports both events together. It
falls back to the phase-broken approximation: under the
hypothesis that sample t carries the compound, each constituent's
likelihood is computed separately and combined under an
independence assumption that overstates support. Stage 5 then
sets `CA = 1` on this sample's final FORMAT field so a downstream
consumer that cares about the approximation can filter on it.

## Recap

Phase chains are an explicit, per-read record of which allele
observations came from the same physical DNA molecule. They
answer one question only — "do these two observations co-occur on
the same molecule?" — and they answer it cheaply and reliably at
the project's coverage target. They are not haplotype assignments,
not population-phasing, and not assemblies; they are just shared
identifiers attached to the observations a single read or pair
contributes.

In spirit, the chains are closest to **freebayes' internal
`RegisteredAlignment` per-read linkage — but persisted**. Freebayes
computes that linkage, uses it for one round of variant calling,
and discards it. Our pipeline separates per-sample evidence
collection from joint posterior computation, so Stage 1 has to
keep the linkage in a form Stage 5 can consume later, possibly
across many cohort recalls. The slot encoding does that without
inventing a global chain-id namespace, without paying for assembly
the way GATK does, and without any extra burden on Stage 2 beyond
"write these small integers."

The on-disk encoding is small: per-record `new_chains` and
`expired_chains` lifecycle markers plus a per-allele list of slot
ids. Stage 2's job is to write the integers verbatim. Stage 5's
job is to walk the records, maintain its own active-slot view,
and intersect slot sets when it needs to evaluate a compound
haplotype.

The two design decisions that make the chain mechanism cheap and
correct are spelled out in
[per_sample_caller.md §"Phase chain identifiers"](per_sample_caller.md):
slot recycling bounded by lifecycle markers, and per-pair (rather
than per-mate) chain assignment via a small in-flight QNAME map.
This document explains the *why*; that one specifies the *what*
in implementation terms.
