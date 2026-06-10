# SSR/STR genotyping — specification

**Status:** design specification, 2026-06-09. Distilled from the research report
[`doc/devel/reports/research/ssr_genotyping_2026-06-09.md`](../reports/research/ssr_genotyping_2026-06-09.md),
which carries the literature evidence, alternatives considered, and the
decision trail. This file is the authoritative, implementation-facing summary;
where the two disagree, this file wins and the research report should be updated.

This specification follows the project's
[`design_principles.md`](design_principles.md) (clarity paramount; reasons and
trade-offs spelled out; no magic numbers — every non-trivial literal is a named
`const` with units and source).

---

## Glossary — acronyms & terms

Defined once here; used unqualified thereafter.

**Domain / biology**

- **SSR** — Simple Sequence Repeat (microsatellite); a tandem repeat with a short
  motif. Used interchangeably with **STR**.
- **STR** — Short Tandem Repeat. Synonym of SSR in this document.
- **VNTR** — Variable Number Tandem Repeat; the broader family (minisatellites and
  up), explicitly out of scope here (period ≤ 6 only).
- **period** — length of the repeat unit (motif) in bp; period 1 = mono-, …, 6 =
  hexanucleotide.
- **motif / RU** — the repeat unit itself (e.g. `CAG`). RU = Repeat Unit (the VCF
  field name).
- **stutter** — PCR/replication slippage that shifts the observed repeat count by
  whole units; the dominant STR artifact, modelled in Stage 2.
- **FRR** — Fully Repetitive Read; a read lying entirely inside the repeat tract
  (no flank), used as a length lower bound for long alleles.
- **purity** — fraction of the tract matching a perfect motif tiling; 1.0 =
  perfect repeat, < 1.0 = imperfect / interrupted.

**Tools / file formats**

- **TRF** — Tandem Repeats Finder; the catalog detector (Stage 0).
- **DUST / sdust** — low-complexity sequence masker; optional TRF prefilter.
- **BAM / CRAM** — aligned-read containers (input).
- **VCF** — Variant Call Format (output). **GFF / BED / TSV** — annotation/tabular
  formats used for the catalog.
- **Parquet** — columnar storage format for the per-sample evidence table.
- **psp** — the project's per-sample pileup/evidence artifact ("STR `.psp`"); here
  the Stage-1 evidence file.
- **md5** — checksum used to bind artifacts to their reference/catalog.

**Methods / statistics**

- **CSR** — Compressed Sparse Row: a ragged array stored as one flat values array
  plus an offsets array marking each row's slice. Here: `amb_read_offsets`
  (offsets) + `amb_lengths`/`amb_logliks` (values) hold a variable number of
  candidate lengths per ambiguous read.
- **HMM** — Hidden Markov Model; the banded pair-HMM is the slow-path realigner.
- **pair-HMM** — a 3-state (Match/Insertion/Deletion) HMM scoring two sequences.
- **BAQ** — Base Alignment Quality; the existing engine whose banded-forward
  pattern the pair-HMM borrows.
- **EM** — Expectation-Maximization; the Stage-2 cohort inference loop.
- **HWE** — Hardy-Weinberg Equilibrium; the random-mating genotype prior, here
  relaxed by the inbreeding index `F`.
- **IBD** — Identity By Descent; two alleles inherited from a common ancestor
  (the `F` component of the prior).
- **SMM** — Stepwise Mutation Model; the reference-centred allele-length prior.
- **`θ` (stutter)** — the stutter-kernel parameters `(u, d, ρ)` of §5.2. Distinct
  from `θ_pop`.
- **`θ_pop` = `4·Nₑ·μ`** — population mutation parameter (`Nₑ` effective size, `μ`
  per-generation mutation rate); sets the SMM allele-spread variance in §5.5. Not
  the stutter `θ`.
- **PCR** — Polymerase Chain Reaction; its absence/presence (PCR-free vs PCR+)
  drives the stutter kernel.
- **FP / FN** — False Positive / False Negative. **QC** — Quality Control.

**VCF / genotype fields**

- **GT** — genotype; **DP** — depth; **GQ** — genotype quality; **GP** — genotype
  posterior probabilities; **GL / PL** — genotype likelihoods (log / Phred).
- **AF / AC / AN** — allele frequency / allele count / allele number; **NS** —
  number of samples with data.
- **REPCN** — per-allele repeat copy numbers (GangSTR field).
- **F** — fixation index (inbreeding coefficient), per sample.

---

## 1. Overview, scope, goals

A tool that **genotypes short tandem repeats (SSRs / microsatellites / STRs)
from aligned short reads (BAM/CRAM)**, producing **per-individual repeat-allele
lengths (copy numbers)** across a cohort, for **population genetics, diversity,
GWAS and breeding**.

It is a **length genotyper**: the genotype at a locus is a multiset of integer
repeat counts, not a base substitution.

### 1.1 Two independent callers from one BAM (foundational)

The SSR caller and the project's SNP/indel caller are **fully independent
pipelines that share only the raw alignments**, each emitting its own VCF (an
allele-base VCF and an allele-length VCF). Their statistical models are
fundamentally different — per-position substitution/indel likelihoods vs
per-locus repeat-length + stutter — and keeping the call sets independent lets
population hypotheses be analysed independently on each marker type. Concretely:

- A **standalone pipeline** (own binary/crate, own VCF); it shares only low-level
  BAM/CRAM I/O (noodles) and the vendored `sdust`.
- It does **not** read the SNP call set, and is **not** a tweak to the
  per-position `.psp` path — the SSR data model is locus/window-oriented.
- **No SNP-phasing.** HipSTR gains accuracy by physically phasing an STR against
  nearby heterozygous SNPs, but SSRs are too sparse for reads to co-span an STR
  and an informative SNP, and it would couple the two callers. Excluded by
  design.

### 1.2 Priorities

- **Precision ≫ recall (FP-averse).** Prefer a no-call over a low-confidence call.
  Quality is a first-class output and the operating point is tunable and measured.
- **Population-first.** The cohort is pooled to learn stutter and allele
  frequencies and to set genotype priors.
- **Memory-efficient, cohort-scaling**, in the project's house style (columnar
  intermediates; the two-stage extract→genotype pattern mirrors
  `pileup → .psp → var-calling`).

### 1.3 Non-goals (explicit)

- Pathogenic repeat-**expansion** detection (alleles ≫ read length) and forensic
  panel typing are out of scope.
- **No aneuploidy / no per-locus or per-contig copy number / no mixed ploidy.** A
  single uniform ploidy per run; chromosomes of differing ploidy are run
  separately.
- **Spanning-reads only** (see §1.4).
- No coupling to CNV calling; no SNP-phasing.

### 1.4 Read-length scope

**Spanning reads only** — reads that fully cross the repeat with anchored flanks
on both sides give exact lengths; this is the precise evidence and the
precision-first choice. Alleles longer than the read are systematically missed
(acceptable under a precision-first mandate; most plant SSR markers are short).
A cheap in-paradigm extension — **merging overlapping read pairs** into longer
fragments — is allowed to raise the spanning ceiling. The beyond-read-length read
classes (flanking / in-repeat / insert-size, à la GangSTR) are **deferred**
(§12).

---

## 2. Architecture

Three pipeline stages plus a test simulator. Stages 1–2 reproduce the project's
two-stage `pileup → .psp → var-calling` pattern: heavy per-sample work once,
summarised to a columnar artifact, then a light cohort math stage.

```
reference FASTA ──► [Stage 0: catalog builder] ──► catalog (.ssr_catalog.bed.gz)
                                                        │
BAM/CRAM (per sample) ─┐                                ▼
                       └─► [Stage 1: per-sample extract] ──► evidence (.parquet, one per sample)
                                                        │
                                                        ▼
                              [Stage 2: cohort genotyping] ──► VCF (allele lengths, GangSTR-compatible)
```

**Self-describing artifacts (unifying rule):** every artifact carries its own
metadata — the catalog in a `##` header, the evidence in the Parquet footer, the
VCF in `##` headers. No sidecar files.

**Locus identity** is positional: the catalog defines the row universe and order;
each per-sample evidence file has exactly one row per catalog locus in catalog
order, so the cohort stage is a synchronized N-file scan.

---

## 3. Stage 0 — Catalog builder

Build the locus catalog de novo from the reference (no prior catalog is assumed
for new species).

### 3.1 Algorithm

- **Tandem Repeats Finder (TRF)** is the detector — it handles imperfect /
  degenerate repeats (the real biology) and is the validated source of every
  established STR catalog. (Fast exact scanners like MISA/PERF are rejected as
  primary because perfect-only recall silently drops interrupted loci.)
- **Optional DUST/sdust pre-search** to accelerate: run sdust genome-wide and
  restrict TRF to masked windows + a margin. Kept **only if** shown near-lossless
  *and* genome-wide TRF is too slow on the target genome (sdust is not
  motif-aware and can miss/mis-bound loci). Must be measured (§3.4).
- **Post-process** TRF output: keep **period ≤ 6** (SSR scope), filter by
  purity/score, **merge overlapping/redundant** calls, drop loci without **unique
  flanks** (mappability). These filters are accuracy knobs (§3.4).
- **Imperfect single-motif loci are kept and genotyped.** **Compound loci are
  split** into their component single-motif sub-loci (each a normal
  perfect/imperfect locus).
  - **Inner-flank consequence (carries into Stage 1):** a split sub-locus's
    *inner* flank is the adjacent motif's repeat — *not* unique sequence. Stage 1
    must anchor on the **outer** flank and treat the inner boundary as
    catalog-known repeat structure, not a random unique flank.

### 3.2 Catalog format

One **self-describing bgzip+tabix BED-like TSV**: a VCF/GFF-style `##` metadata
header (reference path + md5, TRF params, filters, tool/version, date), then a
`#`-prefixed column header, then rows. Tabix skips comment lines. No sidecar.

**Minimal schema (only non-derivable columns):**

```
chrom   start   end   motif   purity
```

- `start`/`end` are 0-based half-open; `end − start` is the reference tract length
  (the reference allele); `motif` gives the period.
- Dropped as derivable: `period` = len(motif); `ref_copies` = (end−start)/period;
  `class` (perfect/imperfect) = threshold(purity), perfect ⇔ purity = 1.0;
  `locus_id` = f(chrom,start,motif). Cross-file linkage is **positional**; the VCF
  `ID` is constructed at output.
- `purity` is retained only because Stage 2 may use it for confidence weighting;
  it may be dropped if it ends up a build-time filter only.

### 3.3 Decisions

- TRF primary; DUST prefilter conditional and measured.
- Compounds split; imperfect kept.
- No redundant columns; QA at the test level.

### 3.4 Accuracy harness (catalog)

Two distinct measurements:
1. **Catalog accuracy** — randomized-sequence FP rate (any "SSR" in shuffled
   sequence is a false positive → calibrates thresholds); simulation recall by
   motif × copy × purity; boundary accuracy; and **recovery of tomato published
   capillary SSR markers** (known motif + sizes) as the species anchor.
2. **DUST-prefilter recall** (if used) — TRF genome-wide vs TRF-within-sdust-
   windows: fraction of true loci dropped, by motif/purity, vs wall-time saved.
   Adopt the prefilter only if the recall loss is negligible.

---

## 4. Stage 1 — Per-sample evidence extraction

For each catalog locus, read the BAM and emit a compact, **stutter-free**
per-read length-likelihood summary — the "STR `.psp`".

### 4.1 Read handling

- Pull reads overlapping `[locus − flank, locus + flank]` (noodles).
- **Use the full read sequence including soft-clipped bases.** The aligner
  soft-clips exactly the bases carrying a non-reference allele length; naive CIGAR
  parsing of the aligned span silently undercounts large alleles. Anchor on the
  original alignment, then re-examine the whole read.
- A read is **spanning** if it crosses the tract with ≥ `MIN_FLANK_BP` clean
  matched bases on both sides; otherwise it is flanking / in-repeat and is counted
  but (spanning-only scope) not used for the length likelihood in v1.

### 4.2 Qᵣ(L) — the per-read length likelihood (two-tier)

**Reads fall into two natural classes, which is why the likelihood is computed
two ways.** For many reads the repeat count is unambiguous: the read clearly
spans the tract with clean unique flanks on both sides, and the tract itself is a
clean integer tiling of the motif — you can simply *count* the repeats and read
off a single confident length. Other reads are ambiguous, for any of several
reasons: an interior sequencing indel or substitution blurs the motif boundaries,
the repeat is impure/interrupted so the count is not integer-obvious, low base
quality near a boundary makes the last unit uncertain, or the read only partially
covers the tract. For these, no single length is defensible — the honest answer
is a *distribution* over plausible lengths. The fast path handles the first class
by direct counting; the slow path handles the second by full probabilistic
realignment.

`Qᵣ(L)` is the likelihood the read came from a molecule of repeat length `L`
units, under **sequencing/alignment error only** — stutter is *not* applied here
(it is learned and applied in Stage 2). Both tiers realign against candidate-length
haplotypes `H_L = outer_flank + (motif × L) + inner_context` (flanks/motif from
catalog + reference) to **escape reference bias**.

- **Fast path — flank-anchored exact motif count.** A read qualifies iff it spans
  with ≥ `MIN_FLANK_BP` clean flanks on both sides, its tract is a pure integer
  motif tiling with no interior sequencing indel, and boundary base-quality ≥
  `MIN_BASE_QUAL`. → a confident length `L*` (+ weight). Handles the bulk in
  O(read length).
- **Slow path — banded pair-HMM forward.** For impure / ambiguous reads (anything
  failing the fast-path test, bundled as `AMBIGUITY_THRESHOLD`), compute a forward
  (sum-over-alignments) score against `H_L` for `L ∈ [count − W, count + W]`
  (`W = STUTTER_WINDOW_UNITS`). → a real `Qᵣ(L)` distribution. **Forward, not
  max** (true likelihood, honest ambiguity).
  - **Base-error model:** Dindel/Albers — a 3-state pair-HMM (Match/Insertion/
    Deletion), per-base-quality emissions
    (`match = 1 − 10^(−Q/10)`, `mismatch = (10^(−Q/10))/3`), affine gaps for
    *sequencing* indels (distinct from stutter). Implemented as a **bespoke banded
    pair-HMM in Rust** (no `rust-bio` dependency), learning the banded-forward
    pattern from `baq_engine` but not coupling to it.
  - **Compound inner-flank:** anchor on the outer unique flank; build the inner
    context of `H_L` from the catalog neighbour (`motif × ref_copies`) via
    coordinate adjacency. A sub-locus sandwiched between two repeats → flag
    low-confidence.

**How `Qᵣ(L)` is stored — the two tiers map onto two storage regimes (§4.3),
not a uniform per-read profile.** A read produces *either* a single length or a
sparse distribution, never a dense vector:

- **Fast-path reads store one length, and not per-read.** A confident `L*`
  carries no distribution — `Qᵣ` is a point mass at `L*`. These reads are *not*
  written individually; they collapse into the locus-level histogram
  (`hist_lengths` / `hist_counts` / `hist_weight`). The per-read identity is
  discarded; only the tally survives. This is the bulk of reads.
- **Slow-path reads store many likelihoods, per-read.** The forward is
  *evaluated* over the `2W + 1` candidate lengths `L ∈ [count − W, count + W]`
  (`W = STUTTER_WINDOW_UNITS`), but the *stored* profile is **sparse**: drop
  candidates whose log-lik is more than `AMB_LL_DROP` below the per-read max,
  then renormalize over what remains. So the per-read entry count is **variable
  and usually ≪ 2W + 1** (typically the 1–3 lengths with real mass). These are
  the only reads written individually, in the `amb_*` CSR columns
  (`amb_read_offsets` gives each read's variable-length slice).

So: one stored length for the confident majority (aggregated, not per-read), and
a short sparse `(length, log-lik)` list for each genuinely ambiguous read.

> **Naming note:** `STUTTER_WINDOW_UNITS` sizes the candidate-`L` range even
> though Stage 1 is stutter-free — the window must be wide enough to cover
> plausible stutter excursions so Stage 2 has the support points to convolve the
> stutter kernel against. The window bounds *which lengths get a likelihood*, not
> a stutter computation here.

Performance is negligible (sparse loci, most reads fast-path) next to the
whole-genome alignment already paid for the BAM.

### 4.3 Evidence format (per-sample Parquet)

**Parquet, one file per sample, a single table, exactly one row per catalog
locus in catalog (genomic) order — including no-coverage loci as zero/empty
rows** (they compress to ~nothing). Two invariants:

- **Synchronized scan:** row `i` is the same locus in every sample → Stage 2 reads
  row `i` across all files, no join/merge.
- **Position query:** rows sorted by `(chrom, start)` → Parquet row-group min/max
  + page index prune to row groups (no tabix). Optional `contig → row-group range`
  footer map for O(1) seek.

| column | type | notes |
|---|---|---|
| `chrom` | dict&lt;string&gt; | contig |
| `start` | int32 | 0-based tract start; rows sorted by (chrom,start) |
| `end` | int32 | tract end (exclusive) |
| `depth` | int32 | total reads at locus |
| `n_spanning` | int32 | usable spanning reads |
| `n_flanking` | int32 | count only in v1; per-read length *lower bounds* (allele ≥ X) would be added here if scope lifts beyond spanning-only to call long alleles |
| `n_frr` | int32 | fully-repetitive |
| `n_filtered` | int32 | low-mapq / dup / clipped |
| `n_flank_indel` | int32 | |
| `mapped_reads` | int32 | for normalized-depth QC |
| `hist_lengths` | list&lt;int16&gt; | distinct observed allele lengths (repeat units), ascending |
| `hist_counts` | list&lt;int32&gt; | confident-read count per length (parallel) |
| `hist_weight` | list&lt;float32&gt; | optional base-qual aggregate per length |
| `amb_read_offsets` | list&lt;int32&gt; | CSR prefix offsets for ambiguous reads (len = n_amb + 1) |
| `amb_lengths` | list&lt;int16&gt; | flattened per-read candidate lengths |
| `amb_logliks` | list&lt;float32&gt; | flattened **stutter-free** log-liks (parallel) |

Confident reads collapse to the `(length → count)` histogram; genuinely bimodal
reads carry an explicit sparse `Qᵣ(L)` profile in the `amb_*` CSR columns (in
v1, not deferred). Footer key-value metadata: `schema_version, sample_name,
reference_md5, catalog_md5, n_loci, ploidy, extraction_params{MIN_FLANK_BP,
MIN_BASE_QUAL, AMBIGUITY_THRESHOLD, STUTTER_WINDOW_UNITS, AMB_LL_DROP},
tool_version` + a
contig name↔id table. Write knobs: sort by `(chrom,start)`; row-group ≈ a few ×
10⁴ loci; page index; ZSTD. All stored likelihoods are **stutter-free**; lengths
are `int16` repeat units.

---

## 5. Stage 2 — Cohort genotyping

From the per-sample evidence, jointly genotype the cohort with population-informed
posteriors, biased toward precision.

### 5.1 Generative model (per locus ℓ)

1. **Population allele frequencies `π = (π_a)` over candidate allele lengths.** At
   this locus an allele *is* a repeat length (10 units, 11, 12, …); the candidate
   lengths are the plausible ones. `π_a` is the fraction of chromosomes in the
   cohort carrying length `a`, so `π` is a probability distribution over those
   lengths (`π_a ≥ 0`, `Σ_a π_a = 1`). It is one of the two unknowns Stage-2 EM
   estimates (§5.4) and the root of the hierarchy: `π` → genotype (2) → reads (3).
2. Sample `s` draws genotype `G_s` = multiset of `ploidy` alleles under the
   **inbreeding-adjusted prior** `P(G_s | π, F_s)` (§5.3).
3. **Read `r` in `s` is born in three steps.** *(a)* It comes off one physical
   homolog, so it picks one allele `a ∈ G_s` uniformly — probability `1/ploidy`
   (½ for diploid); this is what makes a heterozygote's reads a 50/50 mixture of
   its two alleles. *(b)* **Stutter `S_θ`** then maps the true allele length `a` to
   a **molecule length `L`** by adding/dropping whole repeat units (PCR slippage);
   `L` is a *hidden intermediate*, never observed. *(c)* **Sequencing error** reads
   that length-`L` molecule imperfectly — exactly the Stage-1 `Qᵣ(L)` (§4.2). The two
   noise sources are ordered and kept apart on purpose: stutter (b) is a shared
   chemistry property *learned in Stage 2*; sequencing error (c) is per-read and
   was computed *stutter-free in Stage 1* (§4.2). Summing over the hidden `L` gives the
   per-allele read likelihood below.

`P(reads_s | G_s) = Πᵣ [ (1/ploidy)·Σ_{a∈G_s} P(read_r | a, θ) ]`, with
`P(read_r | a, θ) = Σ_L Qᵣ(L)·S_θ(L | a)` — the convolution of the stutter-free
Stage-1 likelihood with the stutter kernel. **Stutter lives entirely in Stage 2.**

### 5.2 Stutter model `S_θ`

**What stutter is.** When a tandem repeat is replicated — by the polymerase in
PCR, or in vivo — the nascent and template strands can transiently dissociate and
**re-anneal out of register by one or more whole repeat units** (replication
slippage). The result is a molecule whose repeat count differs from the true
allele by an integer number of units, *before* a single base is sequenced. This
is the dominant STR artifact: it manufactures a "ladder" of minor length variants
flanking the true allele and is what turns a clean homozygote into something that
*looks* heterozygous. `S_θ(L | a)` is the probability model for it — given the
true allele length `a` (in units), the distribution over the **molecule length
`L`** that actually enters sequencing.

**The kernel, term by term.** Write the slippage as a whole-unit step
`δ = (L − a)/period` (`δ = 0` no slip, `δ > 0` gain/expansion, `δ < 0`
loss/contraction). HipSTR's 3-parameter geometric form:

```
S_θ(δ) =  1 − u − d              δ = 0   (faithful copy — most reads)
          u · ρ · (1−ρ)^(δ−1)    δ > 0   (gain of δ units)
          d · ρ · (1−ρ)^(−δ−1)   δ < 0   (loss of |δ| units)
```

- **`u`** = total probability the molecule gained units, **`d`** = total
  probability it lost units (each summed over all step sizes); `1 − u − d` =
  probability of a faithful copy. They are **separate** because slippage is
  biased — STRs contract more often than they expand, so typically `d > u`
  ("contraction bias").
- **`ρ`** = the geometric decay: *given* a slip occurred, the chance it is exactly
  one more unit before stopping. So a ±1-unit slip is the most likely, ±2 is
  `(1−ρ)×` as likely, and so on — large jumps fall off geometrically (mean step
  size `1/ρ`). Real STR stutter is overwhelmingly ±1 unit, which this captures.
- The product `u·ρ(1−ρ)^(δ−1)` is just "a gain happened (`u`) **and** it was
  exactly `δ` units (geometric)"; the `δ < 0` branch is the mirror image with `d`.
  The three branches sum to 1 over all `δ` — it is a proper distribution.

**Only in-frame.** The kernel moves probability *only* in whole-unit steps,
because slippage re-registers on the repeat period. A length change that is **not**
a multiple of the motif (e.g. a 1 bp indel inside a trinucleotide tract) is **not**
stutter — it is a sequencing/replication base error, and was already accounted for
by the Stage-1 `Qᵣ(L)` term (§4.2). Splitting them this way keeps each noise source
in exactly one place.

**`θ` is covariate-parameterised**, predicted per allele:
`θ = θ(period, allele_length, motif, purity)`. This is biology, not bookkeeping —
stutter rate is not constant: it rises sharply with **allele length** (longer
tracts slip more), is highest for **short periods** (mononucleotide runs are the
worst, hexamers mild), varies by **motif** (`AT` vs `GC` content), and falls with
**purity** (interruptions arrest slippage). So `(u, d, ρ)` are *functions* of these
covariates rather than one global triple.

**Inference — joint mixture-deconvolution (D3).** We never get to see `a`
directly; we see a pile of read lengths per (sample, locus) that is a **blend of
`ploidy` stutter ladders, one centred on each true allele**:
`O_{s,ℓ}(x) = (1/ploidy)·Σ_{a∈G_{s,ℓ}} K_{θ(cov)}(x − a)`. Fitting means
simultaneously (i) placing the allele centres `a` (the genotype) and (ii) learning
the kernel shape `K_θ`. Two design choices make this identifiable:

- **The kernel is shared, the alleles are local.** `K_θ` is a property of the
  chemistry + motif, **identical across all samples of a library** — so pooling
  reads across many (sample, locus) units of the same covariate cell pins it down,
  and it **auto-adapts to PCR vs PCR-free** (PCR-free libraries simply fit a much
  smaller `u`, `d`). The allele centres, by contrast, are **per-(sample,locus)
  latent** — estimated locally.
- **Stutter is strictly within-(sample,locus).** Raw read lengths are **never**
  pooled across samples to estimate allele positions — doing so would confuse
  one sample's stutter ladder with genuine **population** length variation (that
  is `π`'s job, §5.1). Only the *kernel parameters* pool; the *allele calls* stay
  local.
- **Heterozygotes need no special-casing.** A het is just a ≥2-component mixture;
  modelling the mixture directly means we never need a separate "is this sample
  homozygous?" classifier — the data decide how many centres are supported.

The covariate kernel is fit by **regression + hierarchical shrinkage**: covariate
cells with abundant data (common motifs) get their own estimate; sparse cells
**shrink toward the broader period-level estimate** rather than overfitting a
handful of reads. Reads pool across all (sample, locus) units that fall in a
covariate cell.

### 5.3 Genotype prior — inbreeding-adjusted (Wright / HWE-with-`F`)

Plants are far from HWE (tomato is highly selfing), so the prior is **not** strict
HWE. **Reuse the SNP caller's exact multiallelic IBD-mixture prior**
([`src/var_calling/posterior_engine.rs`](../../../src/var_calling/posterior_engine.rs))
with a **per-sample fixation index `F ∈ [0,1]`** (`0` = outcrossing/HWE, `1` =
full inbreeding):

- homozygote `i`: `P(ii | π, F) = F·π_i + (1 − F)·π_i²`
- heterozygote `i ≠ j`: `P(ij | π, F) = (1 − F)·2·π_i·π_j`

(With prob `F` the two alleles are identical-by-descent; else two independent HWE
draws.) STRs differ from SNPs only in that the allele set is repeat lengths and is
highly multiallelic. High `F` in a selfing line a-priori down-weights
heterozygotes — a *second* structural defence against stutter-driven false hets.

### 5.4 Cohort EM (per locus; parameters `π`, `θ`; latent `G_s`)

- **E-step:** `γ_s(G) = P(G | reads_s, π, θ, F_s) ∝ P(reads_s | G, θ)·P(G | π, F_s)`.
- **M-step (allele freqs):** `π_a ← (1/(ploidy·N))·Σ_s Σ_G γ_s(G)·count_a(G)`
  (standard AF-EM → AC/AN/AF). **Dirichlet-smoothed** with pseudocount `α` and a
  **non-uniform base measure** (§5.5) so plausible rare alleles are never
  hard-zeroed and stutter-only candidates self-prune.
- **M-step (stutter):** update the shared covariate kernel from
  responsibility-weighted observed-vs-assigned length residuals (the deconvolution
  update); HipSTR's estimators (`u` = posterior-weighted fraction longer, `d` =
  fraction shorter, `ρ` = geometric MLE).
- **M-step (`F`, optional):** estimate a cohort/per-sample `f̂` from the
  heterozygote deficit only above an effective-sample threshold; else use the
  supplied default. Per-locus `F` is not estimated.

**Structure — initialisation by a confident-homozygote pre-pass.** EM is
non-convex and the stutter kernel and the genotypes are mutually dependent (learn
the kernel ⇄ know which reads are stutter), so a good starting kernel matters. The
pre-pass exploits the one case where the genotype is knowable *without* a stutter
model: a **confident homozygote**, selected by two criteria — **deep** spanning
coverage (the length histogram is well-resolved) **and** a **single dominant
length** (one length holds the overwhelming majority of reads, the rest a small
ladder around it). For such a (sample, locus) there is exactly one true allele
`a` = the modal length, so **every non-modal read is, by construction, a stutter
product of that single allele** — a labelled stutter observation with known step
`δ = (L − a)/period`, no genotype inference required. Pooling these across many
confident-homozygote loci within a covariate cell yields a clean initial
`(u, d, ρ)` via the same estimators as the M-step (`u` = fraction longer,
`d` = fraction shorter, `ρ` = geometric MLE of the ladder decay). The joint
deconvolution EM then **refines** that seed using *all* loci, now including
heterozygotes — where the two alleles' ladders overlap and the genotype is
genuinely latent. This homozygote init is the spec's chosen defence against local
optima (see Identifiability below); multi-restart is deferred.
**Convergence:** penalised-log-lik / parameter tolerance + max-iter; assert
non-decreasing log-lik.

**Identifiability:** within one locus, adjacent-het vs hom+stutter is ambiguous;
broken by (a) the **shared kernel** across many loci (it cannot be simultaneously
small and large) and (b) **population recurrence** (a real allele recurs at stable
frequency; a stutter satellite's mass tracks its parent's read count). No
label-switching (alleles are ordered lengths). Local optima → the homozygote init;
multi-restart deferred.

### 5.5 Allele-length prior (the Dirichlet base measure)

**What it is and where it plugs in.** The §5.4 M-step estimates allele
frequencies `π` as **Dirichlet-smoothed** counts. A Dirichlet prior needs two
things: a **strength** `α` (pseudocount) and a **base measure** `G₀` — a prior
distribution over candidate lengths answering "before any data, where do we expect
allele-frequency mass to sit?". The smoothed update is the blend

```
π_a  ←  ( α·G₀(a)  +  Σ_s Σ_G γ_s(G)·count_a(G) )  /  ( α + ploidy·N )
```

so `G₀` is the fallback `π` leans on when data are thin. This section specifies
`G₀`.

**Not uniform.** A uniform `G₀` would call every candidate length equally
plausible a priori — contradicting microsatellite biology. Under the **stepwise
mutation model (SMM)**, within-population allele-size distributions are
**unimodal**: a modal length with frequency decaying as you move away in whole
repeat units (Valdes, Slatkin & Freimer 1993). A length 8 units off the mode is a
priori far less likely than ±1; a uniform prior would, at small N, let far-out or
stutter-spawned lengths hold real frequency.

**The shape — reference-centred, unimodal** in the signed offset
`Δ = (L − ref)/period` (units from the reference allele):

- **v1 — geometric**, `G₀ ∝ p^|Δ|` (`p < 1`): mass falls by a constant factor per
  unit away from `ref`. One parameter, robust.
- **upgrade — discretized Gaussian** centred at `ref`, variance a hyperparameter.
  The variance may optionally be set **per-locus** from SMM theory: the
  within-population variance of allele size ≈ `θ_pop / 2`, where
  `θ_pop = 4·Nₑ·μ` is the **population mutation parameter** (the SMM equilibrium
  link). ⚠ This `θ_pop` is **not** the stutter kernel `θ` of §5.2 — distinct
  quantity, distinct symbol.

**Deliberately weak — its job is small N.** `G₀` enters with weight `α` against
`ploidy·N` observations, so at useful cohort sizes the **empirical data swamp it**
and `π` ≈ the observed allele frequencies (as it should). Its whole purpose is the
**small-N / single-sample** regime, where the data cannot pin `π` down: there `G₀`
regularises `π` toward a biologically plausible unimodal shape around `ref` instead
of something degenerate. A guardrail that fades as evidence accumulates.

**Population spread ≠ stutter decay.** Both produce ladders of nearby lengths, but
they differ in kind: population spread (`G₀`/`π`) is **real, heritable** allele
variation (evolutionary, SMM); stutter (`S_θ`, §5.2) is a **non-heritable PCR /
measurement artifact**. They are kept as **separate parameters** precisely so the
model can attribute a minor length either to a true rare allele or to slippage —
collapsing them would reintroduce the confusion §5.4 works to avoid.

### 5.6 Small-N / single-sample

Single-sample is a **first-class use case, by construction, no N=1 special-case**:
stutter stays estimable via the genome-wide motif-class kernel; the
Dirichlet-smoothed `π` self-interpolates from the base measure (small N) to the
empirical population AF (large N); `F` is supplied at low N (default 0). Selfing
species supply abundant homozygotes for the kernel init.

### 5.7 Ploidy & `F` input

- **Ploidy:** a single uniform `--ploidy` per run (default `2`); the genotype space
  is the **integer partitions of that ploidy** over the candidate alleles (the
  ConSTRain enumeration trick). No aneuploidy / per-locus CN / mixed ploidy.
  Polyploid supported (uniform); polyploid `F` uses a single-parameter polysomic
  approximation (ignores double reduction — a v-later refinement).
- **`F`:** per-sample `fixation_index_default` (default `0`) + optional per-sample
  overrides; optionally estimated (§5.4).

### 5.8 FP-aversion (precision ≫ recall)

- The **stutter model** explains the ladder as noise (first defence); the
  **inbreeding prior** suppresses spurious hets (second defence).
- **Posterior threshold + precision/call-rate sweep** (reuse the SNP QUAL-sweep
  tooling) sets the operating point; prefer no-call when posterior mass is diffuse.
  HipSTR-derived defaults: posterior `Q ≥ MIN_POSTERIOR`, ≥ `MIN_SPANNING_READS`,
  ≤ `MAX_STUTTER_FRAC` stutter/flank-indel reads, ≥ `MIN_ALLELE_SUPPORT_FRAC`
  allele support; normalized-depth bounds; segdup/mappability exclusion.
- **Stutter-artifact fraction is computed here, not in Stage 1.** Whether a read
  is a stutter artifact is defined relative to the called allele and the fitted
  kernel `S_θ` — both unknown until this stage. After EM, each read's
  stutter-vs-true responsibility (the `S_θ(L | a)` mass in §5.1's convolution)
  gives the per-call artifact fraction directly; the `MAX_STUTTER_FRAC` gate
  applies to that. Stage 1 stays stutter-free and stores only the raw
  `hist_lengths`/`hist_counts` this is derived from — it carries **no**
  `n_stutter_artifact` column.

### 5.9 Output VCF — GangSTR-format-compatible

Maximally compatible with the **TRtools** ecosystem (dumpSTR/mergeSTR/statSTR/
qcSTR/associaTR), verified against the vendored `TRTools/trtools/utils/
tr_harmonizer.py`.

- **REF/ALT = actual repeat-tract sequences** (not symbolic); per-allele copy
  numbers in **`REPCN`**.
- **Required:** INFO `END`, `RU` (motif), `PERIOD`, `REF` (reference copy number);
  FORMAT `GT`, `DP`, `Q` (genotype posterior, 0–1 — TRtools' quality field),
  `REPCN`.
- **Our model → GangSTR fields:** stutter `(u,d,ρ)` → INFO
  `STUTTERUP`/`STUTTERDOWN`/`STUTTERP`.
- **Our extras (TRtools ignores; carry our full posterior + Beagle compat):**
  FORMAT `GP`, `GL`/`PL`, `GQ`; INFO `AF`/`AC`/`AN`, `NS`, `F`. FILTER `PASS` +
  `lowGQ`/`lowDepth`/`segdup`/`lowSupport`. *(TRtools reads only the scalar `Q`, so
  `Q` must carry the called genotype's posterior.)*
- **Detection:** the harmonizer types GangSTR on the `gangstr` header substring —
  triggered honestly via a `##source` line stating "GangSTR-compatible" (or
  `--vcftype gangstr`). No TRtools change needed.
- **Ploidy caveat:** diploid is fully harmonizable; polyploid VCFs are valid but
  TRtools is diploid-centric.

---

## 6. Parameters (named constants — units & source)

Defaults are starting points; the precision-critical ones are swept (§7).

| constant | role | default / source |
|---|---|---|
| `MIN_FLANK_BP` | clean flank bases for a spanning read | ~10 bp (HipSTR) |
| `MIN_BASE_QUAL` | boundary base-quality for the fast path | ~Q20 |
| `AMBIGUITY_THRESHOLD` | bundles fast-path gate (flanks, 0 interior indels, base-qual) | — |
| `STUTTER_WINDOW_UNITS` (`W`) | Qᵣ(L) candidate window ±units | 3 |
| `ploidy` | uniform genotype size | 2 |
| `fixation_index_default` (`F`) | inbreeding prior | 0.0 (SNP caller) |
| `alpha` (`α`) | Dirichlet pseudocount on `π` | — |
| base-measure decay `p` | reference-centred allele-length prior | — (SMM) |
| `MIN_POSTERIOR` (`Q`) | call-quality gate | 0.9 (HipSTR) |
| `MIN_SPANNING_READS` | depth gate | ~10 (HipSTR) |
| `MAX_STUTTER_FRAC` | artifact-read gate | 0.10 (HipSTR) |
| `MIN_ALLELE_SUPPORT_FRAC` | per-allele support | 0.20 (HipSTR) |

---

## 7. Validation & testing (two buckets)

**Bucket 1 — synthetic data for code tests (build critical path).** A bespoke
STR-aware simulator emitting at two levels: **evidence-level** (synthetic Parquet
`Qᵣ(L)` → tests Stage 2 statistics in isolation, before extraction exists) and
**read/BAM-level** (synthetic BAM + truth genotype table → tests the whole
pipeline incl. the pair-HMM). **Anti-tautology rule:** generate data whose
stutter/error *deviates* from the caller's `(u,d,ρ)` assumptions, and keep the
simulator's model definition separate from the caller's code. Tests: deterministic
unit fixtures + property/statistical tests (recovers injected `π/θ/F`; calibrated
posteriors; no false hets on monomorphic loci; di-nuc / long-allele / low-coverage
/ polyploid stress).

**Bucket 2 — comparison vs reference tools (post-build benchmark; data off the
critical path).** Human gold standard: HG002 → GIAB-TR v1.0 truth, compared with
**Truvari bench+refine**, against **HipSTR/GangSTR** (vendored) — reusing the SNP
cross-caller benchmark + dashboard methodology. Tomato best-effort: matched
CE+WGS if findable, else Mendelian/replicate/cross-tool concordance + HWE/`F`
sanity. Truth-free QC always on (Mendelian error on the GIAB trio, replicate
concordance, HWE/spectrum sanity, call-rate). Acceptance target (FP-averse form):
genotype precision at a chosen posterior operating point with call-rate reported,
per-motif-length breakdown, precision ≥ HipSTR/GangSTR at matched call-rate.

---

## 8. CLI surface (provisional — see §11/E)

Mirrors `pileup → var-calling`:

- `ssr-catalog` — reference FASTA → catalog.
- `ssr-extract` — BAM/CRAM + reference + catalog → per-sample evidence Parquet.
- `ssr-genotype` — N evidence files + catalog → cohort VCF.
- `ssr-simulate` — test/dev: inject genotypes+stutter → synthetic BAM and/or
  evidence + truth table.

Repo/crate placement (same-repo separate binary vs new repo) is the open
structural decision (E), to be settled before the implementation plan.

---

## 9. Lessons adopted from prior tools

HipSTR — per-locus stutter EM + realignment to candidate alleles (the core);
in-frame stutter; confidence filtering. GangSTR — multi-class likelihood (its
beyond-read-length classes are on file for §12). ConSTRain — integer-partition
genotype enumeration for ploidy. STR-FM/lobSTR — flank anchoring; motif from data.
popSTR2 — population pooling; covariate (regression) error model. RepeatSeq —
Bayesian posterior + reference-informed prior. EnsembleTR/TRtools — emit a
compatible VCF so the ecosystem works.

---

## 10. References

Papers (full texts read; PMC): HipSTR (PMC5482724); GangSTR (PMC6735967);
ConSTRain (PMC12504596); RepeatSeq (PMC3592458); Valdes/Slatkin/Freimer SMM
(PMC1205356); GIAB-TR v1.0 (PMC11952744). Vendored reference repos (gitignored):
`TRTools/`, `HipSTR/`, `GangSTR/`. Tooling: Truvari (`github.com/ACEnglish/truvari`).

---

## 11. Open structural decision (E)

Repo/crate placement and final CLI naming. Everything else in this spec is
settled; once E is decided, each stage (0, 1, 2) and the simulator can be turned
into its own implementation plan, in data-flow order, with Bucket-1 synthetic
tests on the critical path.

---

## 12. Deferred / future

- Beyond-read-length classes (GangSTR Geometric/Normal-insert/Poisson-FRR/
  Uniform-flank) + read-pair merging.
- Polyploid `F` with partial IBD / double reduction.
- Multi-restart for hard-locus local optima; a second stutter-refinement round
  (full joint hierarchical EM) if measurements demand it.
- Discretized-Gaussian (θ-linked) allele-length base measure if geometric proves
  insufficient.
