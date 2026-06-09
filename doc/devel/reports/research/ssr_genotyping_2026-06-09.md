# Research summary — genotyping SSRs/STRs from Illumina short reads

**Date:** 2026-06-09
**Scope:** algorithms and software for **population/germline length genotyping**
of short tandem repeats (STRs / microsatellites / SSRs) from Illumina
short-read data — i.e. calling the **repeat-allele length (copy number)** per
individual at known or discovered loci, for population genetics, diversity,
GWAS and breeding. **Out of scope:** pathogenic repeat-expansion detection,
forensic panel typing. Method-agnostic across organisms (human tools are the
methodological state of the art; plant/polyploid relevance flagged).
**Goal of this line of work:** implement a tool that genotypes SSRs **directly
from BAM/CRAM files**, ideally reusing this project's cohort infrastructure.

## 0. Provenance & confidence

This summary merges (a) an automated deep-research run (5 search angles, 20
primary sources, 99 candidate claims, adversarial 3-vote verification) whose
**verification + synthesis phase was cut short by a Claude session rate-limit**,
and (b) expert synthesis of the primary algorithm papers. Confidence tags:

- **[V]** — confirmed by the run's adversarial vote (2–3 of 3 verifiers agreed).
- **[S]** — taken from a fetched primary source, but the verification vote did
  not complete (rate-limited; abstained ≠ refuted). Plausible, not independently
  re-checked this run.
- **[B]** — established background/domain knowledge filling structural gaps the
  run did not reach. Confirm against the cited paper before implementing.

One claim was genuinely **refuted**: that HipSTR has *no* length ceiling beyond
read length (1–2 vote) — it *is* read-length-bound, consistent with [V] below.

## 1. Why SSRs need their own algorithm (and why our SNP pipeline fails)

An SSR genotype is a **repeat count (allele length)**, not a base substitution.
Two physical facts dominate:

- **PCR/replication stutter** produces a ladder of products differing from the
  true allele by ±1 (occasionally ±2–3) motif units. Stutter is
  **contraction-biased** (deletions > insertions), **grows with allele length**,
  and **shrinks with motif unit length** [S, NAR 47/5/2436]. **Dinucleotide
  motifs have the highest genotyping error — > 2× any other motif length**
  [V, PMC10984476].
- **Read-length ceiling**: a read can only directly measure an allele it fully
  spans.

Why standard SNP/indel callers (GATK HaplotypeCaller, freebayes, **our own
caller**) underperform on SSRs:

1. **No stutter model** — they treat each ladder product as an independent indel
   allele instead of collapsing the ladder into one length call; evidence
   fragments, genotypes go wrong. [B]
2. **Multi-allelic reality** — population STR loci are often >2 alleles and
   individuals heterozygous for two *non-reference* lengths; biallelic indel
   logic mishandles this. [B]
3. **Normalization / left-alignment ambiguity** — an indel in a repeat has many
   equivalent placements; representation diverges between caller and truth, so
   allele match (POS+REF+ALT) fails even when the locus *is* called. **Measured
   in our human benchmark: 57 of our residual indel false-negatives were
   representation mismatches at STR loci.** [B + in-house]
4. **Low-complexity masking** — DUST-style masks delete the repeat outright.
   **Measured in-house: DUST dropped 222 of our 311 missed truth indels.** The
   literature corroborates that low-complexity regions (~2% of the genome)
   concentrate a large fraction of false-positive indel calls [S, PMC4271055].

Conclusion: SSR genotyping is a **locus-oriented, per-read, realignment**
problem — structurally different from our per-position columnar SNP path.

## 2. The design space (decisions every SSR genotyper makes)

Every tool is a choice on five axes:

1. **Locus model** — fixed catalog (BED of known/TRF loci) vs de-novo per region.
2. **Read evidence** — spanning only, or + flanking / in-repeat / insert-size.
3. **Allele-length extraction** — *how a read becomes a length*: trust the
   aligner CIGAR, remap flanks, k-mer count, or **realign to candidate alleles**.
   **The single biggest accuracy lever.**
4. **Error/likelihood model** — `P(observed | true allele)`: none → fixed ±1 →
   **per-locus learned stutter**.
5. **Genotype inference** — tally/vote → Bayesian → ML+EM → population/joint;
   and **ploidy** (diploid vs N).

### Read evidence classes (set the read-length ceiling)

| class | tells you | precision |
|---|---|---|
| **Spanning/enclosing** (flank–STR–flank in one read) | exact length ≤ read len | high |
| **Flanking/anchored** (one end in unique flank) | lower-bounds length | medium |
| **In-repeat reads (IRR/FRR)** (read/mate fully inside repeat) | allele ≫ read len | low (count) |
| **Insert-size** (paired-end fragment shift) | indel size | low |

HipSTR uses **no** mate-pair info → **read-length-bound**; GangSTR &
ExpansionHunter add insert-size + spanning reads to genotype **beyond** read
length [V, PMC10984476]. For common within-read SSRs, **spanning evidence is the
precise one**; the extra classes mainly buy the long tail.

## 3. The core algorithmic paradigms (advantages / limitations)

**P1 — CIGAR / alignment tally.** Read allele length straight from the existing
BAM alignment at the locus (RepeatSeq; STRsensor primary path; naive callers).
- ✅ Trivial, fast, uses base qualities.
- ❌ **Inherits the aligner's reference bias** (reads far from the reference
  length soft-clip / mis-map → large alleles undercounted); CIGAR placement in a
  repeat is ambiguous. *This is what our current pipeline effectively does.*

**P2 — Flank-anchored / split remapping.** Map only the unique flanks; the
spanned gap = allele (STR-FM; lobSTR's flank alignment).
- ✅ Removes whole-read reference bias; mapper-agnostic.
- ❌ Needs unique flanks both sides; spanning-only; sensitive to flank SNPs/errors.

**P3 — Realignment to candidate alleles under a stutter HMM.** Enumerate
candidate alleles (observed read lengths + reference ± SNPs), realign each read
to each candidate with an HMM encoding **stutter + base error**, score
`P(read | haplotype)`, pick the ML diploid pair; **learn per-locus stutter by
EM** across samples (HipSTR).
- ✅ **Best accuracy** [V: HipSTR 95.2% vs lobSTR 88.2%, RepeatSeq 57.8% vs
  capillary; 98.9% after filtering the least-confident 10%, nmeth.4267].
  Reference-bias-free, principled stutter, integrates **phased SNPs + population
  data** [V].
- ❌ Heaviest compute; **spanning-only** (read-length ceiling); needs coverage to
  learn per-locus stutter (else a shared prior); diploid by default.

**P4 — Repeat-graph realignment + multi-class likelihood.** Build a sequence
graph (flanks + a repeat node traversable *k* times), realign reads — including
mismapped IRRs found by k-mer content — and combine spanning + flanking + IRR +
insert-size into one likelihood (ExpansionHunter = graph; GangSTR = explicit
multi-class ML).
- ✅ Genotypes **beyond read length**; robust to mismapping; one framework for
  short + long.
- ❌ More machinery; IRR/insert-size evidence is low-resolution → **less precise
  on short common STRs**; insert-size calibration.

**P5 — k-mer / alignment-free.** Count flank-anchored / motif k-mers (STRsensor
fallback; assembly methods). STRsensor uses locally-unique flank k-mers
(typically **k=8**) in hash tables, preferring CIGAR and falling back to k-mers
[V, bbae637].
- ✅ Works when alignment fails; fast.
- ❌ Coarse length; needs unique flank k-mers.

**P6 — Copy-number / ploidy-aware enumeration (orthogonal).** For a given
per-locus copy number, enumerate all genotypes (multisets of lengths) and pick
the best-fitting; user can override CN for aneuploid/polyploid regions
(ConSTRain).
- ✅ **The only paradigm built for polyploids** (demonstrated on triploid banana;
  trisomy-21 sim) [S, Commun Biol 2025] — directly relevant to plant work; fast.
- ❌ Combinatorial in CN; still bounded by the feeding extraction method's ceiling.

## 4. Tools → paradigm map

| Tool | Year | Extraction | Stutter model | Genotype inference | > read len | Ploidy | Status |
|---|---|---|---|---|---|---|---|
| lobSTR | 2012 | P2 + FFT/entropy read-sensing | simple learned noise | allelotype from length dist | no | 2 | superseded [B] |
| RepeatSeq | 2013 | P1 | learned empirical | Bayesian posterior | no | 2 | legacy [V baseline] |
| STR-FM | 2015 | P2 | empirical profiles | tally | no | 2 | legacy [S] |
| **HipSTR** | 2017 | **P3 (HMM realign)** | **per-locus EM, in/out-of-frame** | **ML+EM, SNP-phased, population** | no | 2 | **standard for common STRs**, maintained [V] |
| popSTR/popSTR2 | 2017/2020 | P1/P2 | per-marker logistic, attribute-based | likelihood + **population pooling** | partial | 2 | cohort-scale [S] |
| GangSTR | 2019 | P4 (multi-class ML) | within spanning class | **joint ML over classes** | **yes** | 2 | maintained [V] |
| ExpansionHunter (v5) | 2019→ | **P4 (graph)** | within model | graph-evidence inference | **yes** | 2 | actively maintained (Illumina/DRAGEN); highest per-allele accuracy in one TR benchmark [S, PMC10197592] |
| STRsensor | 2024 | P1 + **P5 fallback** | **2-param (u,d) MLE** | **MAP** | no | 2 | new [V] |
| ConSTRain | 2025 | P1 + **P6** | yes | **CN-enumeration best-fit** | no | **N** | new; HG002 98.28%, faster than GangSTR/HipSTR [S] |
| *contrast:* STRetch | 2018 | decoy refs | — | — | expansion | — | expansion-only (out of scope) |
| *contrast:* Tandem-genotypes | 2019 | **long reads** | — | — | yes | — | not Illumina |
| *meta:* EnsembleTR | 2023 | merges call sets | — | — | 2 | — | builds consensus panels [S] |

**Historical arc:** extraction **tally (P1) → flank-remap (P2) → haplotype-HMM
(P3) → graph/multi-class (P4)**; stutter **none → fixed → per-locus-EM**;
inference **vote → Bayesian → ML/EM → population/joint**.

### Locus discovery (step 0, needed for non-model genomes)

- **Tandem Repeats Finder (TRF, 1999)** — the standard reference scanner. [B]
- **MISA / Krait / pal_finder / GMATA** — SSR mining & marker design, heavily
  used in the plant world; Krait appears actively available [S, krait.biosv.com].

## 5. Key quantitative findings (verified)

- HipSTR 95.2% / lobSTR 88.2% / RepeatSeq 57.8% vs capillary; HipSTR → 98.9%
  after filtering the least-confident 10% [V, nmeth.4267].
- >85% of common-STR calls identical across HipSTR/GangSTR/ExpansionHunter; ≥89%
  of those match the reference [V, PMC10984476].
- Dinucleotide STRs: highest error, >2× other motifs [V, PMC10984476].
- HipSTR is read-length-bound (no mate-pair); GangSTR/EH exceed read length
  [V, PMC10984476].
- HipSTR: per-locus stutter EM + HMM realignment + phased-SNP/population
  integration [V, nmeth.4267].
- STRsensor: CIGAR + k=8 k-mer fallback; 2-param (u,d) MLE stutter; MAP
  genotype [V, bbae637].
- [S, not re-verified] ExpansionHunter highest per-allele accuracy across
  motifs/sizes on one genome-wide TR truth set, where HipSTR no-called 21.4% of
  loci (PMC10197592); ConSTRain HG002 ~1.39M/1.7M loci at 98.28% in ~20 min/32
  threads, faster than GangSTR (97.69%, 14.9 h) and HipSTR (97.74%, 12.0 h)
  (Commun Biol 2025); PCR-free libraries cut STR error severalfold (Genome Res
  25/5).

## 6. Recommended architecture for a BAM-input genotyper (this project)

Our goal (BAM → per-individual SSR length genotype, **cohort**, plant-relevant)
points to a clear synthesis: **the HipSTR realignment paradigm (P3) as the
core**, with two grafts that play to our strengths — **population/joint
inference** (we are already a cohort caller; pooling demonstrably helps) and
**optional ploidy-awareness (P6, ConSTRain-style)** for polyploid crops.

1. **Locus catalog.** Input a BED of STR loci. For non-model/plant, build it with
   **TRF** (or MISA/Krait) on the assembly — no curated catalog exists, so
   catalog-building is step 0.
2. **Per-locus read collection** from the BAM: all reads overlapping
   `[locus − flank, locus + flank]`. A **window/locus** data model — *not* our
   per-position columnar one.
3. **Allele extraction = realignment (P3), not CIGAR (P1).** Build candidate
   alleles from observed spanning-read lengths + reference; realign each read to
   each candidate with a stutter-aware HMM; emit `P(read | allele)`. **This
   directly fixes both failure modes we measured** (DUST masking → genotype the
   locus instead of masking it; representation ambiguity → emit a *repeat-count*
   genotype, not an indel allele).
4. **Stutter model:** per-locus, learned by **EM across the cohort**; pool to a
   **motif-class prior** (condition on motif length & allele length; di-nuc worst
   [V]) when a locus is low-coverage.
5. **Genotype inference:** ML/MAP diploid with a length prior; emit the **full
   per-length likelihood vector** + a quality. Add **ConSTRain-style CN
   enumeration (P6)** as the ploidy generalization for polyploid samples.
6. **Output:** VCF using the **repeat-unit / copy-number representation** (GA4GH
   STR conventions) — sidesteps normalization ambiguity by construction.
7. **Later extensions:** GangSTR-style FRR + insert-size (P4) only if alleles
   beyond read length matter. (**SNP-phasing deliberately excluded** — see
   "Independence" in §6b.)

### Critical caveat for our codebase

SSR genotyping is **locus/window-oriented and needs per-read evidence across the
whole repeat** (realignment over a haplotype). Our pipeline is **per-position
columnar** (`pileup` → `.psp` → per-position merge), which structurally cannot do
P3. So this is a **parallel module that reads the BAM directly** (or a
locus-aware `.psp` extension retaining per-read spanning sequences over a
window), **not** a tweak to the SNP path. The realignment HMM is the part to get
right — read the HipSTR source/paper for the exact transition structure first.

### Implementation pitfalls

- Candidate-allele generation **must include lengths not seen in the reference**,
  or reference bias is reintroduced.
- The stutter EM needs enough loci/coverage — design the **motif-class fallback
  prior** up front.
- **Flank uniqueness fails** in tandem-of-tandems and segmental duplications —
  filter the catalog.
- Di-/homopolymer loci are the permanent **error floor** [V] — expose a per-class
  confidence filter rather than chasing them.

## 6b. Design v0 — constraints-driven architecture (2026-06-09)

Five constraints set by the PM, and how each lands in the design:

1. **No catalog for some species** → must **discover SSR loci de novo from the
   reference**. (DUST/sdust suggested — see below: useful as a *pre-filter*, not
   a catalog.)
2. **Posteriors preferred, likelihoods acceptable** → a **generative likelihood**
   `P(reads | genotype)` that we turn into **genotype posteriors** with a
   population prior. Emit `GL`/`PL` and `GP`/`GQ`.
3. **Learn the best lesson from every prior tool** → see the table below.
4. **Independent tool, or at minimum a fully independent module** that runs
   standalone and emits **its own VCF coded with allele lengths**.
5. **Population-first design** → joint genotyping: pool the cohort to learn
   stutter + allele frequencies, set HWE priors, emit population fields.

### Best lesson to steal from each prior tool/algorithm (constraint 3)

| Source | Lesson to adopt |
|---|---|
| lobSTR | **Detect the motif/period from the data** (autocorrelation/FFT), don't trust the reference to define it; enables a fast SSR-read sensor |
| STR-FM | **Anchor on unique flanks** to escape whole-read reference-length bias |
| RepeatSeq | **Bayesian posterior** + base-quality-aware error model |
| **HipSTR** | **Per-locus stutter learned by EM + HMM realignment to candidate alleles** + optional SNP-phase/population integration; in-frame vs out-of-frame stutter; confidence filtering. *The core.* |
| popSTR2 | **Population pooling** for the error model + scalable **incremental sample addition** |
| GangSTR | **One ML likelihood unifying heterogeneous read classes** (spanning / flanking / FRR / insert-size) |
| ExpansionHunter | **Repeat-graph** representation of variable copy number; recover mismapped **in-repeat reads by k-mer content** |
| STRsensor | **Pragmatic dual extraction** (CIGAR primary + k-mer fallback) and a **minimal viable 2-param stutter model** |
| ConSTRain | **Explicit per-locus copy number / ploidy** → enumerate genotypes for the CN (polyploids); speed |
| EnsembleTR / **TRtools** | **Conform to a TRtools-compatible VCF** so the popgen ecosystem (mergeSTR / statSTR / qcSTR / dumpSTR / associaTR-GWAS) works out of the box |
| MISA / TRF / Krait | **Motif-aware tandem scanning** to build the catalog from the reference |

### Resulting architecture: three stages, two of them mirroring our SNP pipeline

The PM's "population-first" + "independent module" constraints fit our existing
**two-stage (per-sample extract → cohort joint-genotype)** pattern exactly — the
`pileup → .psp → var-calling` philosophy, reused for a locus-oriented problem.
That gives population priors + stutter pooling "for free".

**Stage 0 — Catalog builder (`ssr-catalog`, constraint 1).** Scan the reference
with a **motif-aware tandem detector** (MISA/Krait-style: motif lengths 1–6,
per-motif min-copy thresholds, bounded imperfection; handle compound/interrupted
SSRs), emitting a catalog: `chrom, start, end, repeat_unit, ref_copies, purity`.
- **DUST/sdust** (already vendored in this repo) is **only a fast pre-filter** to
  prioritise windows — it is *not* motif-aware and over-masks (homopolymers,
  AT-rich), so it cannot define motif, period, or precise boundaries. Use it to
  cut the search space, then run the motif-aware scan; optionally cross-check
  against TRF.

**Stage 1 — Per-sample evidence extraction from BAM (the "STR pileup",
constraints 4 + 5).** For each catalog locus, pull reads over
`[locus − flank, locus + flank]` (noodles BAM/CRAM), classify into spanning /
flanking / in-repeat (GangSTR's classes), **realign spanning reads to candidate
alleles under a stutter+base-error model** (HipSTR), anchoring on flanks
(STR-FM) and confirming the motif from the data (lobSTR). Emit a compact
per-sample, per-locus **evidence artifact** (per-read length-likelihoods +
read-class counts) — the STR analogue of our `.psp`/gVCF (memory-efficient,
cohort-scaling).

**Stage 2 — Joint / population genotyping (constraints 2 + 5).** Pool all samples
per locus:
- **EM #1**: learn per-locus **stutter** parameters (pool = more data; fall back
  to a **motif-class prior** at low coverage; di-nuc worst [V]).
- **EM #2**: estimate **population allele frequencies** → **HWE genotype prior**.
- Per sample: `P(reads | G)` for each candidate genotype `G` (multiset of lengths
  sized by **ploidy** — ConSTRain generalisation for polyploids), × prior →
  **posterior** `P(G | reads)`. Emit `GL`/`PL` **and** `GP`/`GQ`.
- Output a **VCF coded with allele lengths / copy numbers**, **TRtools-compatible**
  (INFO: `RU`, ref copies, `AF`/`AC`/`AN` per length; FORMAT: `GT`, `REPCN` e.g.
  `12/15`, `GL`/`PL`, `GP`, `GQ`). Following the GangSTR/EnsembleTR `REPCN`
  convention unlocks the whole popgen toolchain.

**Independence (constraint 4) — two independent callers from one BAM.** The SNP/
indel caller and the SSR caller are **fully independent pipelines that share only
the raw input data** (BAM/CRAM), each emitting **its own VCF** (an allele-base
VCF and an allele-length VCF). Rationale: their statistical models are
fundamentally different (per-position substitution/indel likelihoods vs
per-locus repeat-length + stutter), and keeping the callsets independent lets
population hypotheses be analysed **independently** on each marker type. Concrete
consequences:
- A **standalone binary/crate** with its own pipeline and VCF; shares only
  low-level BAM/CRAM I/O (noodles) and sdust if useful. The realignment +
  per-read-spanning-sequence-over-a-window data model does **not** fit our
  per-position `.psp`, so extraction is independent (not a `pileup` reuse).
  Same-repo separate binary vs a new repo is a later call.
- **No SNP-phasing / no coupling to the SNP callset.** HipSTR improves accuracy by
  physically phasing the STR against nearby heterozygous SNPs, but that needs a
  read co-spanning the STR *and* an informative SNP — rare here, because SSRs are
  sparse and far apart — and it would couple the two callers, defeating the
  independence goal. **Excluded by design**, not deferred. (The SSR caller never
  reads SNP calls; the only shared substrate is the alignments.)

**Read-length scope.** Start **spanning-only** (precise for common, short plant
SSRs); add GangSTR-style **FRR + insert-size** later only if alleles beyond read
length matter.

### Open design decisions (to resolve before/with the spike)

- **Imperfect/compound SSR model** — how much fuzziness in the catalog scan and
  in realignment (perfect-only is simplest; bounded mismatches/interruptions add
  recall but complexity).
- **Stutter parameterisation** — start STRsensor-simple (2-param up/down +
  in/out-of-frame), allow per-locus EM à la HipSTR; geometric step size optional.
- **Ploidy** — per-sample fixed vs per-locus CN (ConSTRain) for polyploids/CNV.
- **VCF spec** — lock to GangSTR/EnsembleTR `REPCN` for TRtools compatibility
  (recommended) vs HipSTR `GB`.
- **Evidence-artifact format** — design the STR `.psp` equivalent (its on-disk
  schema is what makes the population stage scale; learn from our `.psp`).

## 6c. Stage 0 — Catalog builder (detailed)

*Goal: from a reference FASTA (no prior catalog), emit a locus catalog
`chrom, start, end, repeat_unit, ref_copies, purity` for genotyping. FP-averse:
a spurious locus wastes effort and can produce spurious "variants".*

### Algorithms used by other tools

| Family | Tools | How it works | Advantages | Limitations |
|---|---|---|---|---|
| **Exact / threshold scan** | MISA, PERF, Krait/pytrf, GMATA | Find maximal runs of ≥ k *perfect* copies of a 1–6 bp motif (per-motif-length thresholds, e.g. mono≥10, di≥6, tri–hexa≥5); compound SSRs = two within a max gap | Deterministic, **very fast** (genome-wide in seconds–minutes; PERF/Krait are C-fast), simple, reproducible, no parameters beyond thresholds; plant-world standard (MISA) | **Perfect repeats only** — misses interrupted/degenerate SSRs; thresholds arbitrary; boundaries are exact-run edges (may clip a real-but-impure locus) [B] |
| **Approximate / statistical alignment** | **TRF** (Benson 1999), phobos, **ULTRA**, mreps | TRF: k-tuple seed matches at offset d → statistical candidate detection → wraparound DP alignment to infer consensus motif + copy number; models % identity and indels between copies | **Handles imperfect/fuzzy repeats** (the real biology); any period; well-validated; **the de-facto catalog source** (HipSTR/GangSTR/EH reference catalogs are all TRF-derived + filtered) | Parameter-sensitive (match/mismatch/indel weights, min score); produces **overlapping/redundant** calls needing post-merge; slower; TRF licensing/integration friction (C, registration) [B] |
| **Combinatorial maximal-repetition** | mreps | Finds approximate tandem repeats with formal guarantees (theory of maximal repetitions) | Fast, principled, catches approximate repeats | Less used for SSR catalogs; tuning of the "resolution" param [B, confirm] |
| **Satellite/long-period (contrast)** | TRASH, TideHunter | Plant satellite arrays / long periods | Good for satellites | Not aimed at 1–6 bp SSRs [B, confirm] |
| **Curated community catalogs** | EnsembleTR / GA4GH TR / Gymrek-lab refs | Population-polymorphic TR sets | Best for **human** (skip building) | None exist for most non-model species → must build [S] |

**The field-standard recipe for a genotyping catalog** is *run TRF → filter
(period ≤ 6, sufficient purity, drop overlaps/segmental-dups, require unique
flanks) → format*. MISA/PERF are the fast exact alternative when perfect-SSR
recall suffices.

### How DUST/sdust could be used as a pre-filter

DUST/sdust (already vendored, [[project_reference_codebases]]) scores
**low-complexity** via short-word (triplet) over-representation and emits masked
intervals. SSRs are a *subset* of low-complexity, so:

- **Use:** run sdust genome-wide (cheap), then run the expensive motif-aware
  scan **only inside masked windows + a margin**. Cuts the search space.
- **Limitations / risks:** sdust is **not motif-aware** (gives no period, copy
  number, or precise boundaries), and it **over-masks** (AT-rich / compositionally
  skewed but non-tandem regions) and can **mis-bound or miss short perfect SSRs**
  that don't trip the complexity threshold. So it is a **recall-oriented coarse
  prefilter, never the catalog**.
- **Decision rule:** only worth it if the motif-aware scanner is TRF-class
  (slow). If we use a fast exact scanner (PERF/MISA-class), it runs genome-wide
  directly and **DUST adds nothing** — measure scanner runtime first. If DUST is
  used, tune its threshold to favour recall and **measure DUST window recall vs a
  full genome-wide scan** so we know how many true loci the prefilter drops.

### How to test / measure catalog accuracy

Catalog errors: missed loci (FN), spurious loci (FP — non-tandem low-complexity),
wrong motif/period, wrong boundaries (which set the genotyping flank window).

1. **Cross-tool concordance** — run TRF, PERF/MISA, and ours on the same
   reference; report locus precision/recall, motif agreement, boundary-offset
   distribution. (For human, compare to the EnsembleTR/GA4GH curated catalog.)
2. **Randomized-sequence FP test** — run on shuffled / Markov-matched random
   sequence; any "SSR" found above threshold is a false positive → directly
   measures the FP rate and calibrates thresholds (FP-averse).
3. **Simulation recall** — implant SSRs of known motif/copy/purity, measure
   recall by motif length × copy number × purity, and boundary accuracy.
4. **Polymorphism enrichment (functional)** — fraction of catalog loci that yield
   ≥ 2 alleles in a pilot cohort (a useful catalog is enriched for polymorphic
   loci).
5. **★ Tomato capillary markers as ground truth** — published *S. lycopersicum*
   SSR marker panels have **capillary-validated motif + allele sizes**; check the
   catalog recovers those loci with the right motif and that genotyping later
   reproduces the known sizes. This is the strongest species-specific test.

### Decision (2026-06-09) — TRF, optionally accelerated by a DUST pre-search

**Chosen approach: Tandem Repeats Finder (TRF) as the catalog detector**, because
it handles imperfect/degenerate SSRs (the real biology) and is the validated,
de-facto source of every established STR catalog. The fast exact scanners
(MISA/PERF) are rejected as the primary because perfect-repeat-only recall would
silently drop interrupted loci.

**Acceleration: optional DUST/sdust pre-search.** Run sdust genome-wide first and
restrict TRF to masked windows (+ margin) to cut TRF's runtime. This is kept
**only if** it is shown to be near-lossless *and* genome-wide TRF is too slow on
the target genomes — sdust is not motif-aware and can miss/mis-bound SSRs, so it
must not silently drop true loci.

**We must measure accuracy.** Two distinct measurements, both required:

1. **Catalog accuracy of TRF itself** (per §6c "test / measure"): randomized-
   sequence FP rate, simulation recall by motif × copy × purity, boundary
   accuracy, and — the anchor — recovery of **tomato published capillary SSR
   markers** (known motif + sizes).
2. **DUST-prefilter recall** — run TRF **genome-wide** vs **TRF-within-sdust-
   windows** on the same reference and report the fraction of genome-wide TRF loci
   the prefilter drops (by motif length / purity / locus length), plus the
   wall-time saved. The prefilter is adopted only if recall loss is negligible
   (FP-averse, but we also don't want to lose real polymorphic loci) for a
   worthwhile speedup. Tune the sdust threshold toward recall; add a margin around
   masked windows.

**Practical / integration notes for the Stage-0 plan:**
- TRF is C with a registration/licence and `.dat`/`.html` output — decide
  **vendor-and-wrap vs subprocess**, and write a parser to the catalog schema.
- Post-process TRF output: keep **period ≤ 6** (SSR scope), filter by **purity /
  score**, **merge overlapping/redundant** calls, drop loci without **unique
  flanks** (mappability) — these filters are themselves accuracy knobs to sweep.
- Parameters to pin/sweep: TRF match/mismatch/indel weights, min score, max
  period; sdust threshold + window margin.
- Evaluate on **both** a human reference (cross-check vs the EnsembleTR/GA4GH
  curated catalog) **and** the tomato reference (capillary-marker truth).

→ This decision seeds the **Stage 0 implementation plan** (catalog builder +
accuracy harness).

### Catalog format & schema (final, 2026-06-09)

- **One self-describing file** — a bgzip+tabix BED-like TSV with a VCF/GFF-style
  `##` metadata header (reference path + md5, TRF params, filters, tool/version,
  date), then a `#`-prefixed column header, then rows. **No sidecar** (metadata
  travels with the data; tabix skips comment lines). *Unifying principle: every
  artifact is self-describing — catalog `##` header, evidence Parquet footer
  metadata, VCF `##` header.*
- **Minimal schema — only non-derivable columns** (PM: don't store redundant
  data, QA at the test level):

  ```
  chrom   start   end   motif   purity
  ```

  Dropped as derivable: `period` = len(motif); `ref_copies` = (end−start)/period;
  `class` (perfect/imperfect) = threshold(purity), perfect ⟺ purity = 1.0;
  `locus_id` = f(chrom,start,motif). **Cross-file linkage is positional** (the
  evidence file is one row per catalog locus in catalog order); the VCF `ID` is
  constructed at output. `purity` is the one judgment call — keep iff Stage 2 uses
  it for confidence weighting (noisier low-purity loci → FP-aversion); drop if it
  is only a build-time filter.
- **Imperfect single-motif loci: kept and genotyped.** **Compound loci: split**
  into their component single-motif sub-loci (each a normal perfect/imperfect
  locus) → the catalog builder needs a TRF-output **compound-splitting** post-step.
  - ⚠ **Stage-1 implication (captured for that stage):** a split sub-locus's
    *inner* flank is the adjacent motif's repeat — **not unique sequence**. Flank-
    anchoring / realignment must treat that boundary as *known repeat structure*,
    not a random unique flank, or extraction misbehaves at compound junctions.

## 6d. Stage 1 — Per-sample evidence extraction (detailed)

*Goal: for each catalog locus × sample, read the BAM and emit a compact evidence
artifact (the "STR `.psp`"). **What data is needed depends on the Stage-2
algorithm** — this is the key coupling.*

### What each genotyping algorithm consumes (so we know what to extract)

| Stage-2 approach | Minimal per-locus×sample evidence needed | Artifact size |
|---|---|---|
| **Length-histogram** (lobSTR/RepeatSeq/STRsensor) | multiset `{allele_length → count}` from spanning reads (+ per-read mapq, base-qual summary) | tiny |
| **Per-read realignment likelihoods** (HipSTR) | for each spanning read, `log P(read | allele a)` vs a candidate-allele set (from realignment under base-error+stutter) | medium |
| **Multi-class likelihood** (GangSTR) | above **+** flanking-read length bounds, FRR counts, insert-size deviations | medium |
| **CN best-fit** (ConSTRain) | observed length distribution (histogram) + per-locus copy number | small |

Since we want HipSTR-grade accuracy + posteriors, the artifact should be the
**per-read length-likelihood** form (richer than a bare histogram, so Stage 2 can
do realignment-quality genotyping and learn stutter), but stored compactly.

### Extraction algorithm (spanning-only scope)

1. **Pull reads** overlapping `[locus − flank, locus + flank]` (noodles BAM/CRAM).
2. **Use the FULL read sequence, not just the aligned portion** — the aligner
   **soft-clips the bases that don't match the reference allele**, but those
   clipped bases carry the true allele. Naive CIGAR parsing of the aligned span
   is the classic failure mode. Anchor on the original alignment, then re-examine
   the whole read. [B — critical]
3. **Locate the repeat by flank match** (STR-FM lesson): require ≥ `min_flank`
   matching bases on **both** sides ⇒ the read is *spanning*; else classify as
   flanking (one-sided) or in-repeat (no flank).
4. **Realign / count** the repeat span against candidate allele lengths under a
   base-error model (motif confirmed from the read data, lobSTR lesson); emit
   `log P(read | length)` per candidate. k-mer fallback if flanks don't map
   (STRsensor).
5. **Record per locus × sample:** per spanning read → {observed length (units &
   bp), per-candidate log-lik or a compact mismatch/indel profile, mapq, strand,
   purity flag, base-qual summary}; plus **read-class counts** (#spanning,
   #flanking+bounds, #IRR, depth, #filtered). Plus locus context (motif, ref
   copies, flank seqs).

### Per-sample evidence schema (what to store — the engineering crux)

Determined entirely by the Stage-2 math
`P(reads|G,θ) = Πᵣ (1/ploidy) Σ_{a∈G} Σ_L Qᵣ(L)·S_θ(L|a)` plus the cohort EM for
`θ` (stutter) and `π` (allele freqs). The on-disk schema is what makes the cohort
stage scale (directly analogous to our `.psp`).

**Non-negotiable rule: store the *stutter-free* per-read length likelihood
`Qᵣ(L)` (sequencing/alignment only), NOT `P(read|allele)`.** Baking stutter in at
Stage 1 would make the cohort stutter model unlearnable in Stage 2. Stage 1 does
the heavy per-read realignment once; Stage 2 convolves with `S_θ` and runs the
light EM.

**Per locus × sample, store:**

1. **Per-spanning-read length evidence (θ-free) — the core.** Per read: ML
   observed length `L*` (motif units) + a compact `Qᵣ(L)` profile (log-liks at
   `L*` and a few neighbours), over a window wide enough to cover the cohort
   candidate alleles ± stutter; plus `mapq`, `strand`, base-quality summary, and a
   **purity flag** (read repeat matches motif vs interrupted/SNP).
   - **Compression (cohort-scaling form):** aggregate confident reads into a
     **`(length × quality-class)` count histogram**; keep explicit `Qᵣ(L)`
     profiles only for *ambiguous* reads (e.g. bimodal at a homopolymer edge).
     Columnar, `.psp`-style.
2. **Read-class & QC counts:** depth; `n_spanning`, `n_flanking` (+ length lower
   bounds), `n_FRR`; `n_filtered` (low-mapq/dup/clip); `n_stutter_artifact`,
   `n_flank_indel` (HipSTR's ≤10% gate); per-length support counts (candidate
   `≥2 reads & ≥20%` rule + the `≥20%` support filter); mapped-read count
   (ConSTRain normalized-depth = reads / CN).
3. **Locus reference by id** into the catalog (motif, ref length, flanks) — *not*
   duplicated per sample.
4. **Deferred / optional** (only if spanning-only is lifted): flanking length
   bounds, FRR counts, insert-size deviations (GangSTR classes). *No SNP-phasing
   evidence* — the SSR caller is independent of the SNP callset by design (§6b).

**Why each piece (→ Stage 2):** `Qᵣ(L)` → genotype likelihood **and** stutter EM
(`u,d,ρ` from posterior-weighted observed-vs-assigned length diffs); length
histogram → allele-frequency EM + candidate set; counts → the precision-first
filters and call quality. Candidate alleles are shared across samples by union'ing
each sample's observed-length support in Stage 2, padded by the stutter range —
so the per-sample `Qᵣ(L)` window must extend a few stutter units beyond the
locus's observed lengths.

### Evidence file format (final, A2 — 2026-06-09)

**Parquet, one file per sample, a single table with exactly one row per catalog
locus in catalog (genomic) order — including no-coverage loci as zero/empty
rows** (they compress to ~nothing via RLE/dictionary). Two invariants at once:
- **Synchronized scan:** row `i` is the same locus in every sample's file →
  Stage 2 reads row `i` across all N files, no join/merge logic. (The catalog
  defines the row universe + order; each per-sample file is a parallel column over
  it.)
- **Position query:** rows sorted by `(chrom, start)` → Parquet per-row-group
  **min/max statistics** + the **page index** prune to relevant row groups; no
  tabix needed (DuckDB/pyarrow/polars get range filtering for free). Optional
  `contig → row-group range` map in the footer for O(1) seek.

Single table (not normalized); variable-length evidence goes in **nested list /
CSR columns** (`.psp`-style), so it stays one file per sample.

| column | type | notes |
|---|---|---|
| `chrom` | dict&lt;string&gt; | contig; row-group min/max → position pruning |
| `start` | int32 | 0-based tract start; **rows sorted by (chrom,start)** |
| `end` | int32 | tract end (exclusive) |
| `depth` | int32 | total reads at locus |
| `n_spanning` | int32 | usable spanning reads |
| `n_flanking` | int32 | (length lower-bounds carried later if scope lifts) |
| `n_frr` | int32 | fully-repetitive |
| `n_filtered` | int32 | low-mapq / dup / clipped |
| `n_stutter_artifact` | int32 | HipSTR ≤10% gate |
| `n_flank_indel` | int32 | |
| `mapped_reads` | int32 | ConSTRain normalized depth |
| `hist_lengths` | list&lt;int16&gt; | distinct observed allele lengths (**repeat units**), ascending |
| `hist_counts` | list&lt;int32&gt; | confident-read count per length (parallel to `hist_lengths`) |
| `hist_weight` | list&lt;float32&gt; | *optional* base-qual aggregate weight per length |
| `amb_read_offsets` | list&lt;int32&gt; | CSR prefix offsets for ambiguous reads (len = n_amb + 1) |
| `amb_lengths` | list&lt;int16&gt; | flattened per-read candidate lengths |
| `amb_logliks` | list&lt;float32&gt; | flattened **stutter-free** log-liks (parallel to `amb_lengths`) |

**The ambiguous `Qᵣ(L)` profiles (the `amb_*` CSR columns) are in v1** (PM
decision — not deferred): confident reads collapse to the `(length → count)`
histogram; genuinely bimodal reads carry an explicit sparse profile. (A later
refinement could make the confident block 2-D `(length × quality-class)` via
parallel `hist_lengths/hist_qualbins/hist_counts`.)

**Footer key-value metadata (self-describing, no sidecar):** `schema_version,
sample_name, reference_md5, catalog_md5, n_loci, ploidy(if fixed),
extraction_params{min_mapq, min_flank_bp, min_base_qual, ambiguity_threshold},
tool_version` + a `contig` name↔id table.

**Invariants:** `int16` repeat-unit lengths; **all stored likelihoods are
stutter-free**; every per-sample file has exactly `len(catalog)` rows in catalog
order. **Write knobs:** sort by `(chrom,start)`; row-group ≈ a few × 10⁴ loci;
enable the page/column index; ZSTD.

### Testing

- Per-read length-extraction accuracy vs simulated reads of known allele length
  (recall of correct length by allele size, stutter level, flank-SNP density).
- Spanning-classification correctness (does the soft-clip handling recover large
  alleles the aligner clipped?).
- Compare extracted length histograms to HipSTR's internal evidence on a shared
  BAM+catalog.

## 6e. Stage 2 — Cohort genotyping (detailed)

*Goal: from the per-sample artifacts, jointly genotype the cohort with
**population-informed posteriors**, biased toward **precision over recall**.*

### Statistical approaches used by other tools

| Tool | Statistical model | Output | Notes |
|---|---|---|---|
| RepeatSeq | Bayesian posterior over diploid genotypes, learned error model | MAP GT + posterior | per-sample |
| lobSTR | ML allele-pair under per-locus noise model | GT | per-sample |
| STRsensor | 2-param stutter MLE → **MAP** | GT | per-sample [V] |
| **HipSTR** | **EM: alternate (stutter params) ↔ (per-sample genotypes)** over HMM-realignment read likelihoods; optional population AF + phased-SNP prior | GT + posterior Q | **inherently cohort** (pools to learn stutter) [V] |
| popSTR2 | per-marker logistic (attribute) error model + population AF; per-sample genotyping | GT + qual | scales to 100k+ genomes [S] |
| GangSTR | joint **ML** over read-class likelihoods | GT + **REPCI** (conf. interval) | per-sample [V] |
| ExpansionHunter | graph-evidence likelihood | GT + CI | per-sample [S] |
| ConSTRain | enumerate genotypes for the locus **copy number**, best length-distribution fit | GT | ploidy-aware [S] |

### Accuracy (verified where marked)

- HipSTR **95.2%** vs capillary, **98.9% after dropping the least-confident 10%**
  [V, nmeth.4267] — *confidence filtering is the biggest precision lever*.
- > 85% concordance across HipSTR/GangSTR/EH; **di-nucleotides worst (>2×)** [V].
- ConSTRain ~**98.28%** filtered on HG002 [S]. General pattern: tetra/tri reliable,
  di/homopolymer are the error floor.

### Proposed statistical model (fleshed out, 2026-06-09)

**HipSTR-like in spirit** (cohort-pooled stutter learning + realignment
likelihoods + a population allele-frequency prior = proven SOTA), with three
deliberate adaptations for our goals: a clean two-stage split, a *hierarchical*
stutter prior, and posteriors + an explicit precision-first operating point.

**Conceptual clarification — two objects, two distinct roles (not both priors):**
- **Allele-length distribution `π_ℓ`** (population frequencies at locus ℓ) → the
  **prior** over genotypes, via an **inbreeding-adjusted** (Wright / HWE-with-`F`)
  rule — *not* strict HWE (see "Fixation index" below).
- **Stutter model `S_θ(L | a)`** → part of the **likelihood**: the *generative*
  process by which a *true* allele `a` produces *observed* lengths `L` (the
  PCR/replication slippage ladder). It is **not** a prior.

  `P(G | reads) ∝ P(reads | G, θ) · P(G | π, F)` — stutter in the first factor,
  the allele distribution + inbreeding in the second. The cohort EM learns `π`
  (prior), `θ` (a likelihood parameter), and optionally `F` by pooling all samples
  at the locus.

**Fixation index `F` in the prior (mandatory for plants — they are far from HWE).**
Most plant species (tomato is highly selfing) carry a large excess of homozygotes,
so a strict-HWE prior is badly miscalibrated. We **reuse the SNP caller's exact
prior** ([`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs)):
the **multiallelic Wright / IBD-mixture** genotype prior with a **per-sample
fixation index `F ∈ [0,1]`** (`0` = outcrossing/HWE, `1` = full inbreeding):

- homozygote allele `i`: `P(ii | π, F) = F·π_i + (1 − F)·π_i²`
- heterozygote `i ≠ j`:  `P(ij | π, F) = (1 − F)·2·π_i·π_j`

i.e. with prob `F` the two alleles are identical-by-descent (one draw from `π`),
else two independent HWE draws. Carries over verbatim from SNPs (per-sample
`fixation_index_default` + `fixation_index_overrides`, pseudocounted `π̂`, cohort
estimate `f̂_C`, homogeneous-`F` fast path); STRs differ only in that the allele
set is repeat lengths and is highly multiallelic.

**Synergy with FP-aversion:** in a selfing line a high `F` strongly down-weights
*heterozygous* genotypes a priori — which is exactly the spurious call stutter
tends to produce (a stutter rung read as a second allele). So the inbreeding
prior is a *second* structural defence against stutter-driven false heterozygotes,
on top of the stutter model itself. For polyploids, `F` generalises to the
polysomic IBD coefficient (defer; note double-reduction).

**Generative model (per locus ℓ — what EM inverts):**
1. Population allele frequencies `π = (π_a)` over candidate allele lengths.
2. Sample `s` draws genotype `G_s` = multiset of `ploidy` alleles from `π` under
   the **inbreeding-adjusted prior `P(G_s | π, F_s)`** (Wright/HWE-with-`F`, above).
3. Read `r` in `s`: pick one `a ∈ G_s` uniformly (1/ploidy); the repeat region is
   generated from `a` by **stutter** `S_θ` (true→observed length) then
   **sequencing/alignment error** (observed length→read bases).
   ⇒ `P(reads_s | G_s) = Π_r [ (1/ploidy) Σ_{a∈G_s} P(read_r | a, θ) ]`.

**Where the two error sources live (the key split):**
- **Stage 1 (per-sample, stutter-free):** emit per read
  `Q_r(L) = P(read_r | observed length L)` — the *sequencing/alignment* term
  only; peaked at the length the read's repeat best matches. Compact (short vector
  near best `L`, or `(best L, quality)`). **θ-independent** ⇒ Stage 1 never needs
  the stutter model.
- **Stage 2 (cohort):** convolve with stutter:
  `P(read_r | a, θ) = Σ_L Q_r(L) · S_θ(L | a)`. So the **entire stutter model +
  its EM live in Stage 2**, exactly where the cohort pooling is.

**Cohort EM at locus ℓ (params `π, θ`; latent `G_s`):**
- **E-step:** `γ_s(G) = P(G | reads_s, π, θ, F_s) ∝ P(reads_s | G, θ) · P(G | π, F_s)`.
- **M-step (allele freqs):** `π_a ← (1/(ploidy·N)) Σ_s Σ_G γ_s(G)·count_a(G)`
  (standard AF-EM → AC/AN/AF). **Dirichlet/pseudocount smoothing** so plausible
  **rare alleles are not hard-zeroed**.
- **M-step (fixation index, optional):** if `F` is not supplied (e.g. a known
  selfing rate `F ≈ s/(2−s)`), estimate it from the posterior heterozygote
  deficit — cohort `f̂_C = 1 − Σ_s E[het_s] / Σ_s HWE-expected-het`, or per-sample.
  Mirrors the SNP caller's `f̂_C`.
- **M-step (stutter):** re-fit `θ` from posterior-weighted observed-vs-true length
  differences over all reads × allele-assignments. **Hierarchical shrinkage:**
  `θ_ℓ` drawn from a **motif-class prior** (keyed by motif length × reference
  allele length — the two factors stutter scales with); pool genome-wide for the
  class hyperparameters, shrink each locus toward its class. *(Improvement over
  HipSTR's per-locus-only stutter — robust at low coverage and on the long tail.)*
- Iterate to convergence.

**Output per sample:** MAP genotype + `GP`/`GQ` from `γ_s`; `GL`/`PL` from
`P(reads_s | G, θ̂)`; population `AF`/`AC`/`AN` from `π̂` and the estimated/used
fixation index `f̂` (per locus or cohort). **TRtools-compatible** VCF coded with
allele lengths/copy numbers (FORMAT `GT:REPCN:GL:GP:GQ`, INFO `RU`,`AF`/`AC`/`AN`,`F`).

### Output VCF format (final, A3 — 2026-06-09)

Decided against the **live TRtools source** (cloned at `TRTools/`, gitignored;
verified in `TRTools/trtools/utils/tr_harmonizer.py`): **emit a GangSTR-format-
compatible VCF** — the cleanest TRtools target for a copy-number/length
genotyper (native `REPCN` integer copy numbers; no HipSTR START/END flank-offset
gymnastics).

- **TRtools detection:** the harmonizer types a VCF as GangSTR iff the header
  contains both `command=` *and* the substring `gangstr` (case-insensitive). We
  trigger this **honestly** via a `##source`/`##command` line stating the VCF is
  *GangSTR-compatible* (the substring `GangSTR` legitimately appears); users can
  also force it with `--vcftype gangstr`. **No TRtools source change needed.**
- **REF/ALT = actual repeat-tract sequences** (not symbolic `<STR>`): TRtools
  derives lengths as `len(allele)/period`, and `REPCN` gives the copy number.
- **Required (GangSTR-compatible):** INFO `END`, `RU` (motif), `PERIOD`, `REF`
  (reference copy number, Float); FORMAT `GT`, `DP`, `Q` (genotype posterior 0–1 —
  TRtools' quality field), `REPCN` (per-allele copy numbers, e.g. `12,15`).
- **Our model maps onto GangSTR's own fields:** learned stutter `(u,d,ρ)` →
  INFO `STUTTERUP`/`STUTTERDOWN`/`STUTTERP`; per-allele CI (if bootstrapped) →
  FORMAT `REPCI`, bootstrap SE → `STDERR`.
- **Our extras (TRtools ignores them, harmless):** FORMAT `GP` (posteriors),
  `GL`/`PL`, `GQ` (phred) — carry our full posterior + Beagle/imputation compat;
  INFO `AF`/`AC`/`AN` (per-length population), `NS`, `F` (used/estimated fixation
  index). FILTER `PASS` + `lowGQ`/`lowDepth`/`segdup`/`lowSupport` (§6e/§6f gates).
  *Note TRtools reads only the scalar `Q`, not `GL`/`PL`/`GP` — so `Q` must carry
  the called genotype's posterior.*
- **Ploidy caveat:** GangSTR `REPCN` is `Number=2` (diploid). Diploid cohorts are
  fully TRtools-harmonizable; **polyploid** VCFs stay valid (REPCN with ploidy
  entries + ploidy-general `GT`) but the TRtools harmonizer is human-diploid-
  centric, so polyploid harmonization may need `--vcftype` care or fall outside
  TRtools.

Supersedes the earlier "GA4GH conventions" placeholder. Reference repos cloned
for this (gitignored, per the freebayes/gatk/... convention): `TRTools/`,
`HipSTR/`, `GangSTR/`.

**Differences from HipSTR (deliberate):** (1) clean Stage-1/Stage-2 split with a
θ-free per-sample artifact (mirrors `pileup→.psp→var-calling`, lets the cohort
stage scale); (2) **hierarchical motif-class stutter** vs per-locus-only;
(3) **posteriors + explicit precision-first operating point** as first-class
outputs; (4) **polyploid-general** genotype space (ConSTRain). **We deliberately do NOT
adopt HipSTR's SNP-phasing** — SSRs are too sparse for reads to co-span an STR and
an informative SNP, and coupling to the SNP callset would break the
two-independent-callers design (§6b).

**FP-aversion is baked in (precision ≫ recall):**
- The **stutter model** is the first defence — it explains the ladder as noise
  instead of calling each rung an allele.
- The **HWE/population prior** down-weights genotypes implausible given the cohort
  → suppresses spurious stutter-heterozygotes.
- A **QUAL-like likelihood ratio** (best genotype vs most-likely homozygous-common
  genotype); **no-call** when posterior mass is diffuse.
- **`GQ`/posterior threshold + the precision/call-rate sweep** (reuse the SNP
  QUAL-sweep tooling) sets the operating point; mappability + min-spanning-depth
  hard guards.

**Incremental build (per [[feedback_incremental_steps]]):**
- **v1:** fixed stutter by motif-class (empirical/literature default, no per-locus
  refinement) + per-locus AF-EM (Dirichlet-smoothed) + per-sample diploid
  **inbreeding-adjusted (HWE-with-`F`) posteriors** with a supplied per-sample `F`
  (default + overrides, reusing the SNP path) + GQ filter & sweep. Tractable;
  testable on tomato capillary markers.
- **v2:** hierarchical per-locus stutter shrinkage; polyploid genotypes;
  likelihood-ratio QUAL.
- **v3:** incremental cohort addition (popSTR2-style); beyond-read-length classes
  (GangSTR) if needed. *(No SNP-phasing — excluded by design, §6b.)*

### FP-aversion (the PM's stated priority — precision ≫ recall)

- **Posterior threshold filtering** (HipSTR's 10%-drop → 95.2→98.9 lesson) — emit
  a tunable GQ/posterior cutoff; **prefer no-call (`./.`) over a low-confidence
  call**.
- **Posterior-threshold precision/recall sweep** — *reuse the exact ROC-style
  sweep methodology we just built for SNP QUAL*: plot precision vs call-rate over
  the genotype posterior, pick the FP-averse operating point.
- **Hard guards:** min spanning-read depth; flank-uniqueness / mappability filter
  (drop segmental dups); per-class confidence (di/homopolymer auto-downweighted).
- **Population sanity:** flag extreme **HWE deviation** / implausible allele
  spectra (free once we have cohort AF) — strong FP signal.

### Read-length scope — agree, with one caveat

Spanning-only is the right call **for precision**: spanning reads give exact
lengths; FRR/insert-size evidence is coarse and adds FP — and the PM prioritises
precision, so excluding it is *consistent with the stated goal*. The only cost is
**recall on alleles longer than the read** (systematic FN for long loci), which
is acceptable under a precision-first mandate, and most plant SSR markers are
short (di/tri, < 60 bp). **One cheap extension that stays in the spanning
paradigm:** merge overlapping read pairs (fastp/bbmerge) into ~250–300 bp
fragments before extraction → raises the spanning ceiling without any FRR/insert
machinery. Recommend: **spanning-only + optional read-pair merging; defer
FRR/insert-size** unless long-locus recall becomes a requirement.

## 6f. Lessons from deep-reading HipSTR, GangSTR, ConSTRain (2026-06-09)

Full texts read: HipSTR (PMC5482724), GangSTR (PMC6735967), ConSTRain
(PMC12504596). Extraction was via a fast model over the PMC full text —
**high-confidence on structure/formulas, but confirm exact constants in the
supplements + source before coding.** Concrete refinements to the §6e model:

### Stutter model — adopt HipSTR's exact 3-parameter geometric form

This *is* our `S_θ`. Stutter `δ` measured in **whole repeat units**; per-locus
params `u` (insertion prob), `d` (deletion prob), `ρ` (geometric step size):

```
S_θ(δ) =  1 − u − d              if δ = 0
          u · ρ · (1−ρ)^(δ−1)    if δ > 0   (gain of δ units)
          d · ρ · (1−ρ)^(−δ−1)   if δ < 0   (loss of |δ| units)
```
with `L = a + δ·motif`. Only **in-frame** (whole-unit) stutter is modelled here;
**out-of-frame** changes are absorbed by the base-error term (our Stage-1
`Q_r(L)`), not the stutter model. HipSTR **EM M-step** (gives us the stutter-EM
update directly): `u` = posterior-weighted fraction of reads whose called allele
is *longer* than the assigned true allele; `d` = fraction *shorter*; `ρ` = 1 /
mean-weighted-step-size (the geometric MLE); `f_j` = allele frequency.
- Contrast: **GangSTR uses a cruder fixed `Geometric(P=0.9)`** on enclosing-read
  counts (no per-locus EM); **ConSTRain models no stutter at all** (leans on dup
  marking + depth filters). Both confirm our choice of HipSTR-style per-locus
  (hierarchical, motif-class-shrunk) stutter as the precision lever they lack.

### Stage-1 per-read likelihood `Q_r(L)` — HipSTR's realignment, minus stutter

HipSTR uses a **flanking-sequence HMM** (Match/Ins/Del states, emission
`Q(b,h) = q_b` if match else `(1−q_b)/3`, base-error per **Albers et al. 2011 /
Dindel**) plus an **STR model** that assumes the read differs from a haplotype by
**≤ 1 indel of magnitude D** (a multiple of the motif), inserted bases assumed to
be periodic copies of the motif. Our Stage-1 `Q_r(L)` = exactly this realignment
likelihood **without** the stutter term (we moved stutter to Stage-2). Adopt the
Dindel emission model + the ≤1-indel STR alignment.

### Candidate-allele set — adopt the support thresholds

- HipSTR: a spanning read needs **≥10 bp exact match on both ends** and **no
  longer exact match 15 bp up/down**; a length becomes a candidate if seen in
  **≥2 reads AND ≥20% of a sample's reads**; then **iterate** (re-derive
  candidates from stutter-corrected alignments after a genotyping round).
- GangSTR: candidates from enclosing reads with **support ≥2**. ConSTRain:
  spanning = **≥5 bp flank each side**.
- → **Our rule:** candidate length set = {lengths in ≥2 reads & ≥20% of a sample}
  ∪ reference, **union'd across the cohort**, padded by the stutter range; one
  refinement iteration. (Resolves the §6d "candidate set shared across samples"
  open question.)

### Polyploid genotype space — adopt ConSTRain's integer-partition enumeration

ConSTRain's key tractability trick: enumerate genotypes as **integer partitions
of the copy number `c`** (alleles in descending abundance), **not** weak
compositions — `c = 20` ⇒ **627 partitions vs 6.9×10¹⁰ compositions**. This makes
our polyploid posteriors tractable. Reuse ConSTRain's "**each allele contributes
≈ total_reads / c reads**" as the *mean* of our read mixture, but **keep the
proper likelihood + stutter + posterior** instead of ConSTRain's crude
`argmin ‖D_i − O‖₁` (Manhattan-distance fit, which has no stutter and no
posteriors).

### Quality + FP-aversion — adopt the concrete filter thresholds

- **HipSTR:** genotype posterior **Q ≥ 0.9**; **≥10 spanning reads**; **≤10% reads
  with stutter/flank-indel artifacts**; **≥20% allele support**; "drop the 10%
  least-confident → 95.2→98.9%". These are ready-made FP-averse defaults.
- **GangSTR:** `Q = −10·log₁₀(1 − L_ML/L)` (a uniform-prior posterior; we compute
  the analogue from our *real* population posterior → `GQ`); bootstrap-resample
  reads for a CI (optional, → `REPCI`).
- **ConSTRain:** **normalized depth = reads / copy_number**, keep within
  `[--min-norm-depth=1.0, drop top/bottom 2.5%]`; **exclude segmental
  duplications**. Adopt normalized-depth bounds + segdup/mappability exclusion.
- **Call-rate caveat (precision↔recall):** at matched ~98% accuracy on HG002,
  loci-called were **HipSTR 69% vs GangSTR 81% vs ConSTRain 82%** — HipSTR
  no-calls far more (stricter gate). So our precision-first stance will cost call
  rate; make the trade explicit via the GQ sweep.

### Cohort mechanics

HipSTR genotypes **200 samples per batch**, learning allele frequencies `f_j` in
the EM — direct confirmation of our cohort-pooled design. ConSTRain is **Rust +
multithreaded** (architecturally identical to us; HG002 in ~19.5 min/32 threads).

### On file — GangSTR's beyond-read-length read-class likelihoods

If we ever lift the spanning-only scope (§6e), GangSTR's exact per-class terms:
**Enclosing** `Geometric(P=0.9)` on repeat count; **Spanning** `Normal(μ = L_frag
− A·m, σ)` on fragment length; **FRR** count `Poisson(λ)` with
`λ = C_v·r·[u(A·m−r)+u(B·m−r)]/m`; **Flanking** `Uniform[1, A]`. Joint
`LL(A,B) = Σ_pairs log[0.5·P(·|A)+0.5·P(·|B)] + log Poisson(|FRR|; λ)` — the same
0.5/0.5 diploid mixture we use, plus the Poisson FRR term.

**Net effect on the design:** §6e's `S_θ` is now concretely the `(u,d,ρ)`
geometric form with HipSTR's EM updates; Stage-1 `Q_r(L)` is the Dindel-emission
≤1-indel realignment; the polyploid genotype space uses integer-partition
enumeration; and FP-aversion has ready-made thresholds. No change to the
two-stage split or the inbreeding-adjusted prior.

## 7. Open questions / next steps

- [ ] **Deep-read the 3 core papers** to pin exact formulations before coding:
      HipSTR (realignment HMM + stutter EM), GangSTR (multi-class likelihood),
      ConSTRain (ploidy enumeration).
- [ ] **Re-run the verification + synthesis** of the deep-research workflow after
      the session-limit reset to upgrade the **[S]** claims to **[V]** and add the
      tools the run didn't reach (lobSTR internals, popSTR2, GATK PCR-indel model
      specifics, STR-FM, error-rate-vs-length curves).
- [ ] Decide the **locus catalog** strategy for the target genomes (tomato:
      TRF/MISA on SL reference).
- [ ] Decide BAM-direct module vs **locus-aware `.psp` extension** (retain
      per-read spanning sequence over a window) — affects whether `pileup` is
      reused.
- [ ] Prototype the **realignment + per-locus stutter EM** on a single locus as a
      feasibility spike; validate vs HipSTR on a shared catalog.

## 8. Sources (from the deep-research run)

Primary (verified-source-backed claims):
- HipSTR — Willems et al., *Nat Methods* 2017. https://www.nature.com/articles/nmeth.4267
- Genome-wide HipSTR/GangSTR/ExpansionHunter comparison — https://pmc.ncbi.nlm.nih.gov/articles/PMC10984476/
- STRsensor — *Brief Bioinform* 2024. https://academic.oup.com/bib/article/26/1/bbae637/7922198
- ConSTRain — *Commun Biol* 2025. https://www.nature.com/articles/s42003-025-08837-8
- EnsembleTR / TR benchmark — https://pmc.ncbi.nlm.nih.gov/articles/PMC10197592/
- STR error rates / STR-FM — *Genome Res* 25/5 2015. https://genome.cshlp.org/content/25/5/736.full
- In-vitro stutter / mutability model — *NAR* 47/5 2019. https://academic.oup.com/nar/article/47/5/2436/5304323
- Low-complexity regions & variant-call error — Li 2014. https://pmc.ncbi.nlm.nih.gov/articles/PMC4271055/
- HipSTR docs/repo — https://hipstr-tool.github.io/HipSTR/ , https://github.com/HipSTR-Tool/HipSTR
- Krait (SSR mining) — http://krait.biosv.com/
- Additional fetched: bioinformatics 36/7/2269; PMC7327730; pubmed 34260828;
  frontiers fgene 1474611; bioinformatics 33/24/4041.
  *(Note: the run had mislabelled `PMC6262739` as GangSTR — it is actually
  Šarhanová et al. SSR-seq, unrelated; corrected below.)*

**Deep-read in full (PMC full text, 2026-06-09 — basis for §6f):**
- HipSTR — Willems et al., *Nat Methods* 2017 — https://pmc.ncbi.nlm.nih.gov/articles/PMC5482724/
- GangSTR — Mousavi et al., *NAR* 2019;47(15):e90 — https://pmc.ncbi.nlm.nih.gov/articles/PMC6735967/
- ConSTRain — *Commun Biol* 2025 — https://pmc.ncbi.nlm.nih.gov/articles/PMC12504596/
  (bioRxiv: 10.1101/2024.12.13.628141)

Run stats: 5 angles, 20 sources, 99 claims, 25 verified (8 confirmed, 1 refuted,
16 abstained/rate-limited). Full machine output retained in the session task
log (`tasks/wom9kpb1d.output`).
