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
  motif (period ≤ 6). Used interchangeably with **STR**. **We say SSR, not the
  broader "TR", deliberately** — see the *SSR vs STR vs TR* note below.
- **STR** — Short Tandem Repeat. Synonym of SSR in this document. The
  human-genetics/forensics name for the same short-motif class.
- **TR** — Tandem Repeat; the *umbrella* term (any period, homopolymer →
  satellite) used by the modern long-read ecosystem (GIAB-TR, TRtools, TRGT).
  **Broader than our scope** — we appropriate it only at the interop seam (§5.9
  VCF, GIAB-TR benchmark), never as our own feature word.
- **VNTR** — Variable Number Tandem Repeat; the broader family (minisatellites and
  up), explicitly out of scope here (period ≤ 6 only).

  > **SSR vs STR vs TR — a deliberate naming choice, not an oversight.** The
  > short-motif tandem repeat has several names that differ by field and by
  > scope. *SSR / microsatellite* (period ~1–6 bp) is the dominant term in
  > plant, agricultural, ecological and population genetics — our domain and
  > audience. *STR* is the same class under its human-genetics/forensics name
  > (CODIS markers are STRs); we treat it as a synonym. *TR (tandem repeat)* is
  > the **superset** covering every period from homopolymers through
  > kilobase VNTRs and pathogenic expansions. Modern sequencing tools
  > (GIAB-TR, TRtools, TRGT, LongTR) standardised on "TR" for *principled*
  > reasons — long reads dissolved the assay-defined size boundaries, the old
  > motif-length cutoffs were always arbitrary, and expansion disorders broke
  > the "short" in STR — but those reasons are driven by a scope that spans the
  > whole continuum. **Ours does not:** period ≤ 6, spanning short reads,
  > pop-gen markers; VNTRs and expansions are explicit non-goals (§1.3–§1.4).
  > For *that* scope, SSR is the more precise term and the one our readers use,
  > so we keep it. (And note: **TRF ≠ TR** — TRF is a *tool*, see below.)
- **period** — length of the repeat unit (motif) in bp; period 1 = mono-, …, 6 =
  hexanucleotide.
- **motif / RU** — the repeat unit itself (e.g. `CAG`). RU = Repeat Unit (the VCF
  field name).
- **stutter** — PCR/replication slippage that shifts the observed repeat count by
  whole units; the dominant STR artifact, modelled in Stage 2.
- **FRR** — Fully Repetitive Read; a read lying entirely inside the repeat tract
  (no flank), used as a length lower bound for long alleles.
- **purity** — fraction of the tract matching a perfect motif tiling; 1.0 =
  perfect repeat, < 1.0 = imperfect / interrupted. "purity" is the *concept*
  used in prose; the concrete catalog column / `Locus` field is named
  **`purity_fraction`** (the `_fraction` suffix marks it a degree in [0, 1], not
  a boolean).
- **allele (here)** — a repeat allele's *identity* is its actual tract
  **sequence**, not a bare number. Two molecules that differ by a
  non-motif-multiple amount (e.g. 12 clean repeats vs 12 repeats + 1 bp) are
  **distinct alleles**. The integer **repeat count** is a *derived coordinate*
  read off the sequence (used by stutter and the prior), never the identity. See
  on-ladder / off-ladder.
- **on-ladder / off-ladder** — an allele is **on-ladder** if its length is a
  whole-motif multiple of the reference tiling — a clean rung `ref ± k` units; its
  sequence is fully reconstructible from `(reference tract, k)`, so it is the same
  allele in every sample by construction. An allele is **off-ladder** if it differs
  from a clean tiling by a non-motif amount (a 1 bp indel, a partial unit, a
  *variable* interruption); its sequence is **not** reconstructible from a rung
  number, so it is carried explicitly. Stutter (§5.2) moves only between on-ladder
  rungs; an off-ladder allele is a genuine distinct allele, never a stutter
  product. (A *fixed* interruption shared by the whole population is not off-ladder
  — it is pinned into the reference tract, §3.1.)

**Tools / file formats**

- **TRF** — Tandem Repeats Finder (Benson 1999); the catalog **detector tool**
  run in Stage 0. A *tool name*, never a name for the feature itself — a locus
  is an SSR, not "a TRF". Do not confuse with **TR** (the feature umbrella,
  above).
- **DUST / sdust** — low-complexity sequence masker; optional TRF prefilter.
- **BAM / CRAM** — aligned-read containers (input).
- **VCF** — Variant Call Format (output). **GFF / BED / TSV** — annotation/tabular
  formats used for the catalog.
- **Parquet** — a generic columnar storage format; considered for the per-sample
  evidence but **rejected** in favour of the project's `.psp` columnar-block
  container (§2 / §4.3), which gives positional cross-sample block alignment and
  finer memory control.
- **`.psp` / `.snp.psp` / `.ssr.psp`** — the project's bespoke columnar-block
  evidence container. One shared format (genomic-window blocks, per-column zstd, CSR
  ragged columns, tail block index, TOML header), two schemas: `.snp.psp` for the
  SNP/indel caller, `.ssr.psp` for this one. The `kind` header field is the
  authoritative schema tag; the extension is convenience.
- **md5** — checksum binding artifacts to their reference/catalog. `reference_md5`
  is the md5 of the reference's **concatenated, upper-cased sequence** (per-contig in
  contig order — the same content convention as CRAM's per-contig `M5`), *not* the
  file bytes, so it is invariant to FASTA line-wrapping and (de)compression.
  `catalog_md5` is over the catalog's data rows (post-`#` header).

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

It is a **length genotyper**: the genotype at a locus is a multiset of repeat
**alleles**, not a base substitution. An allele's *identity* is its actual tract
**sequence** — so two molecules differing by a non-motif-multiple amount (12 clean
repeats vs 12 repeats + 1 bp) are distinct alleles, the way HipSTR represents them.
Its integer **repeat count** is a *derived coordinate* read off that sequence,
used by the stutter model (§5.2) and the population prior (§5.5) but never as the
allele's identity. The reason is cross-sample comparison: Stage 2 must decide when
two samples carry "the same" allele, and the integer count is not a faithful key
(it collapses 12 and 12+1 bp onto the same number) — only the sequence is. Most
alleles are **on-ladder** (clean rungs of the integer ladder) and are stored
compactly as a rung number whose sequence is reconstructible from the catalog's
reference tract; the minority that are **off-ladder** carry their sequence
explicitly (the *hybrid* representation, §4.2/§4.3). See the glossary entries for
*allele (here)* and *on-ladder / off-ladder*.

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
fragments — is allowed to raise the spanning ceiling, but it is **not free**: the
overlap must be merged into a consistent base-quality model (combined quals where
mates agree, conflict handling where they disagree), and a chimeric/mis-paired merge
*inside the repeat* fabricates a false length — so merged fragments need the same
slow-path scrutiny (§4.2) as soft-clipped reads, not a blind concatenation. The
beyond-read-length read classes (flanking / in-repeat / insert-size, à la GangSTR)
are **deferred** (§12).

---

## 2. Architecture

Three pipeline stages plus a test simulator. Stages 1–2 reproduce the project's
two-stage `pileup → .snp.psp → var-calling` pattern: heavy per-sample work once,
summarised to a columnar artifact, then a light cohort math stage.

```
reference FASTA ──► [Stage 0: catalog builder] ──► catalog (.ssr_catalog.bed.gz)
                                                        │
BAM/CRAM (per sample) ─┐                                ▼
                       └─► [Stage 1: per-sample extract] ──► evidence (.ssr.psp, one per sample)
                                                        │
                                                        ▼
                              [Stage 2: cohort genotyping] ──► VCF (allele lengths, GangSTR-compatible)
```

**Reuse the `.psp` columnar-block container (not Parquet).** The per-sample evidence
uses the **same bespoke columnar-block format the SNP/indel caller's `.psp` already
uses**, abstracted into a shared container that both callers ride on. The SNP/indel
schema becomes `.snp.psp`, the SSR schema `.ssr.psp` — same container, two column
registries (the `kind` header field carries the authoritative schema identity; the
extension is human/glob convenience). We chose this over Parquet deliberately: it
gives **genomic-window-aligned blocks across samples** (Parquet's byte-sized row
groups can't align positionally), full control of the decode unit and memory
(`--block-window-bp`), and reuse of battle-tested read/write/index/compression code.
The cost — no third-party-tool interop on an internal intermediate, and one extra
column registry to maintain — is acceptable. (See §4.3; structural placement in §11.)

**Self-describing artifacts (unifying rule):** every artifact carries its own
metadata — the catalog in a `##` header, the evidence in the `.psp` TOML header,
the VCF in `##` headers. No sidecar files.

**Locus identity** is positional **by coordinate**, not by row index: the catalog
defines the locus universe; each per-sample evidence file is **sparse**
(no-coverage loci are simply absent) and is keyed by `(chrom, start)` through the
`.psp` block index. The cohort stage is a synchronized scan that **merges by
coordinate** across the N files — the genomic-window block grid makes every
sample's blocks cover the same intervals, so the merge is block-aligned with no
dense empty rows.

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
header (reference path + md5, TRF params, filters, **`flank_bp`** margin,
tool/version, date), then a `#`-prefixed column header, then rows. Tabix skips
comment lines. No sidecar.

**Schema — non-derivable columns, plus the embedded local reference (§"Self-contained"):**

```
chrom   start   end   motif   purity_fraction   ref_seq_start   ref_seq
```

- `start`/`end` are 0-based half-open; `end − start` is the reference tract length
  (the reference allele); `motif` gives the period.
- `motif` is the **verbatim** repeat unit — reference-strand, phase-faithful (TRF's
  consensus unit as it tiles the tract), *not* canonicalized: rotating it would
  break tiling, and reconstruction reads phase-correct bytes from `ref_seq` anyway.
  Its canonical *class* (lexicographically-min rotation over the motif + its
  reverse complement) is **derived** as the §5.2 stutter "motif-or-GC class"
  pooling key, never stored.
- `purity_fraction` is the fraction of the tract that is a clean motif tiling, in
  [0, 1] — a *degree*, not a flag (1.0 = perfect, < 1.0 = interrupted). The
  `_fraction` suffix is deliberate: the bare name reads like a boolean.
- `ref_seq` is the **upper-cased reference bases of the tract plus a `flank_bp`
  margin on each side** (clamped at contig ends); `ref_seq_start` is the 0-based
  genomic coordinate of `ref_seq[0]`, so the tract is
  `ref_seq[(start − ref_seq_start) .. (end − ref_seq_start)]` and the flanks are
  what remain either side. (See §"Self-contained" for why this lives in the
  catalog.)
- Dropped as derivable: `period` = len(motif); `ref_copies` = (end−start)/period;
  `class` (perfect/imperfect) = threshold(purity_fraction), perfect ⇔
  `purity_fraction` = 1.0; `locus_id` = f(chrom,start,motif). Cross-file linkage
  is **positional**; the VCF `ID` is constructed at output.
- **Reference allele identity is the reference tract sequence** (`ref_seq`'s tract
  slice), not `ref_copies`. For an **imperfect** locus `(end − start)` is generally
  *not* a multiple of `period`, so `ref_copies` is **fractional** — it is a derived
  annotation only (rounded for the VCF integer field, §5.9), never the reference
  allele's identity. The integer ladder Stage 1/2 work on is anchored on this
  reference tract; the fractional remainder is carried as fixed reference-tract
  context, not as a ladder rung (§4.2).
- `purity_fraction` is retained only because Stage 2 may use it for confidence
  weighting; it may be dropped if it ends up a build-time filter only.

**Self-contained — the catalog embeds the local reference so the SSR *algorithm*
in Stages 1–2 never reads the FASTA.** `ref_seq` is the *one deliberate
exception* to "only non-derivable columns": the tract+flank bases **are**
derivable from the reference, but storing them is the point. The reference FASTA
is a fixed, multi-hundred-MB resident cost; embedding the few bytes each locus
actually needs means the SSR math — candidate-ladder reconstruction (§4.2),
off-ladder normalization, VCF REF/ALT (§5.9) — sources every reference base from
the catalog. A direct expression of the project's memory-for-scaling thesis, and
it makes the catalog a single self-contained input; the only thing that opens the
reference *for the SSR algorithm* is this Stage-0 builder, which already has it
for TRF.

- ⚠ **CRAM caveat (I/O layer, not the algorithm).** CRAM stores reads as deltas
  against the reference, so **decoding a CRAM** still needs the reference at the
  noodles/htslib layer regardless of this embedding. So: **BAM input is fully
  FASTA-free**; **CRAM input still requires the reference for decoding** (passed
  to `ssr-pileup` for that purpose only — the SSR math never consults it). The
  embedding removes the *algorithm's* reference dependency; it cannot remove
  CRAM's *decoder* dependency.
- **Costs.** The catalog grows by `ref_seq` (~tract + 2·`flank_bp` per locus,
  ≈ tens of MB bgzip-compressed for a plant genome — small beside a genome), and
  it is a **denormalization** of reference-derivable data, kept honest by the
  header's `reference_md5` (any mismatch between `ref_seq` and the declared
  reference is detectable).
- `flank_bp` is sized to Stage 1's read-anchoring + pair-HMM band need (set when
  building `ssr-pileup`); Stage 2 needs only the tract slice, so Stage 1 is the
  binding constraint.

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
per-read length-likelihood summary — the per-sample `.ssr.psp` evidence file.

### 4.1 Read handling

- Pull reads overlapping `[locus − flank, locus + flank]` (noodles).
- **Use the full read sequence including soft-clipped bases.** When a read carries
  an allele longer than the *reference*, the aligner cannot place the extra repeat
  units and often **soft-clips** them — pushing one flank *into the clip*, so naive
  CIGAR parsing of the aligned span silently undercounts large alleles. Anchor on the
  matched side, then **realign the clipped side** to recover where the tract ends and
  the flank resumes. This recovery is a **slow-path** operation (§4.2), not a count
  off the CIGAR — and it doubles as the test that the clip is a real long allele (it
  realigns to extra units + a clean flank) rather than junk (adapter / low-quality
  end, which yields no clean flank).
- **Spanning is a property of the read's recovered *content*, not its CIGAR.** A read
  is **spanning** if its sequence carries ≥ `MIN_FLANK_BP` clean flank on **both**
  sides of the full tract — whether those flank bases were *aligned* or *recovered
  from a soft-clip*. So a soft-clipped read whose clip yields a clean flank **is**
  spanning and carries a real long allele; only when the clip yields **no** clean
  flank does the allele run off the read end ("allele ≥ read length") and the read is
  counted but not used for the length likelihood — the §1.4 scope boundary.
  Otherwise-non-spanning reads are flanking / in-repeat: counted, not used in v1.

### 4.2 Qᵣ — the per-read allele likelihood (two-tier; on-ladder + off-ladder)

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
realignment. A subset of the ambiguous class is not fuzzy at all but cleanly
**off-ladder** — the read supports a definite sequence that simply is not a clean
rung (an in-frame count plus a 1 bp indel, a partial unit, a variable
interruption); the slow path emits that sequence as its own allele rather than
forcing it onto the nearest rung.

`Qᵣ(·)` is the likelihood the read came from a given **candidate allele**, under
**sequencing/alignment error only** — stutter is *not* applied here (it is learned
and applied in Stage 2). Candidates are *sequences*, in two families:

- **on-ladder** — the integer ladder of clean tilings `H_L = outer_flank + (motif ×
  L) + inner_context` for `L` near the observed count (flanks/motif from the
  catalog's embedded `ref_seq` + `motif`, §3.2 — no FASTA access). The ladder is
  catalog-defined, hence **identical across all samples**,
  and is the spine every read is scored against. An on-ladder candidate is
  identified by its rung `L` alone (its sequence is reconstructible from the
  reference tract).
- **off-ladder** — an *actual* tract sequence a read supports that is **not** a
  clean rung. It is identified by its sequence, carried as a **normalized** delta
  from the reference tract (canonical left-aligned form, reusing the SNP caller's
  indel-normalization discipline so the *same* off-ladder allele in two samples
  gets the *same* key — this is what makes the §5.1 cohort union work).

Both families are realigned against to **escape reference bias**. We do **not**
assemble haplotypes across reads (HipSTR does): an off-ladder allele is captured
only when at least one read supports its sequence cleanly enough to emit it. An
off-ladder allele present in the whole cohort *only* in error-buried form is
missed — an accepted precision-first **sensitivity** limit, not an
identity/comparison limit (every *emitted* off-ladder sequence is reconciled across
samples in §5.1).

- **Fast path — flank-anchored exact motif count.** A read qualifies iff both flanks
  are **cleanly aligned** (matched, ≥ `MIN_FLANK_BP`, no soft-clip on either side),
  its tract is a pure integer motif tiling with no interior sequencing indel, and
  boundary base-quality ≥ `MIN_BASE_QUAL`. → a confident on-ladder length `L*`
  (+ weight). A fast-path read is on-ladder by construction (a pure tiling is a clean
  rung). Handles the bulk in O(read length). A **soft-clipped read never qualifies**
  (a clipped flank is not matched), so it falls to the slow path — which is where its
  length is recovered (§4.1).
- **Slow path — banded pair-HMM forward.** For impure / ambiguous reads **and
  soft-clip-recovered spanning reads** (anything failing the fast-path test, bundled
  as `FAST_PATH_GATE`), compute a forward (sum-over-alignments) score against
  `H_L` for `L ∈ [count − W, count + W]` (`W = STUTTER_WINDOW_UNITS`), where `count`
  is the *recovered* repeat count — for a soft-clipped read, the realignment of the
  clip establishes it, so the window centres on the true (possibly large) length, not
  the reference. → a real `Qᵣ(L)` distribution. **Forward, not max** (true likelihood,
  honest ambiguity); a cleanly-long soft-clipped read simply yields a sharply-peaked
  `Qᵣ`. If the read's best alignment is **off-ladder** (its consensus tract is a
  definite non-rung sequence), it *also* emits that normalized off-ladder sequence as
  a candidate with its forward score, so Stage 2 treats it as a distinct allele.
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
- **Off-ladder reads store a sequence — symmetrically to the on-ladder regimes.**
  When a slow-path read supports an off-ladder allele, its normalized sequence
  delta enters the locus-level off-ladder list (`offl_seqs`); confident off-ladder
  reads collapse to an `(offl_seq → count)` tally (`offl_counts`, the off-ladder
  analogue of the histogram), and a read genuinely ambiguous *between* an off-ladder
  sequence and ladder rungs carries its off-ladder leg in the `offl_amb_*` CSR.
  These columns are **empty at the vast majority of loci**, so the cost is paid only
  where real off-ladder signal exists.

So: one stored length for the confident on-ladder majority (aggregated, not
per-read), a short sparse `(length, log-lik)` list for each genuinely ambiguous
read, and a normalized sequence for the rare off-ladder allele.

> **Naming note:** `STUTTER_WINDOW_UNITS` sizes the candidate-`L` range even
> though Stage 1 is stutter-free — the window must be wide enough to cover
> plausible stutter excursions so Stage 2 has the support points to convolve the
> stutter kernel against. The window bounds *which lengths get a likelihood*, not
> a stutter computation here.

**Cost — expected cheap, but to be measured, not asserted.** The relevant baseline
is Stage 1's *own* wall time (and the per-sample SNP pileup it parallels) — **not**
the whole-genome alignment, which was paid once upstream and never re-paid, so it is
no yardstick. Cost = `BAM/CRAM decode over covered regions` + per-read work
(`O(read length)` counting on the fast path; a banded pair-HMM forward over `2W+1`
haplotypes × band × read length on the slow path). We *expect* it cheap — catalog
loci are sparse (a small fraction of the genome, unlike the every-position pileup)
and most reads fast-path — but the **slow-path fraction is the driver**, and it rises
on repeat-rich / impure genomes (the very target of an STR caller) and with soft-clip
recovery (every soft-clipped long-allele read is slow-path, §4.1). So Stage-1
throughput is a claim to **benchmark on a repeat-rich genome** — reporting Stage-1
wall, the fast/slow read split, and the slow-path share — not to assert.

### 4.3 Evidence format (per-sample `.ssr.psp`)

**The shared `.psp` columnar-block container (§2), one file per sample, with the
SSR column registry below.** One **record per covered locus** — the file is
**sparse**: no-coverage loci are simply absent (no dense empty rows). Records are
grouped into **blocks on the genomic-window grid** (`--block-window-bp`), blocks
never cross a chromosome, and a tail-mounted **block index**
`(chrom, first_start, last_end, n_records, offset)` makes the file seekable. Two
invariants:

- **Coordinate-merge scan:** Stage 2 opens the N files and merges records by
  `(chrom, start)`. Because every sample uses the same window grid, blocks cover the
  same intervals, so the merge is **block-aligned** — the sparse analogue of a
  synchronized scan, with no join index and no empty rows.
- **Region query:** binary-search the block index for the first block overlapping
  `[start, end)`, decode forward until past it (`ssr-call --regions`). The index
  keys on **interval** extent — `last_end`, not `last_start` — so a locus that
  starts before a query window but extends into it is not missed (the point-vs-
  interval difference from the SNP `.psp`, whose records are points).

The columns below are the SSR registry; they map onto the container's encodings —
scalars as fixed/varint columns, the per-locus lists (`hist_*`, `amb_*`, `offl_*`)
as the container's **CSR ragged columns**, and `chrom`/`offl_seqs` strings via its
dictionary/bytes columns. Types shown are logical. Each column is independently
zstd-framed, so a reader decodes only the columns a pass needs (column-selective
decode).

| column | type | notes |
|---|---|---|
| `chrom` | dict&lt;string&gt; | contig |
| `start` | int32 | 0-based tract start; records in `(chrom, start)` order (the block-grid key) |
| `end` | int32 | tract end (exclusive) |
| `depth` | int32 | total reads at locus |
| `n_spanning` | int32 | usable spanning reads |
| `n_flanking` | int32 | count only in v1; per-read length *lower bounds* (allele ≥ X) would be added here if scope lifts beyond spanning-only to call long alleles |
| `n_frr` | int32 | fully-repetitive |
| `n_filtered` | int32 | low-mapq / dup / clipped |
| `n_flank_indel` | int32 | |
| `mapped_reads` | int32 | for normalized-depth QC |
| `hist_lengths` | list&lt;uint16&gt; | distinct observed **on-ladder** allele lengths (repeat units, non-negative), ascending |
| `hist_counts` | list&lt;int32&gt; | confident-read count per on-ladder length (parallel) |
| `hist_weight` | list&lt;float32&gt; | optional base-qual aggregate per on-ladder length |
| `amb_read_offsets` | list&lt;int32&gt; | CSR prefix offsets for ambiguous reads (len = n_amb + 1) |
| `amb_lengths` | list&lt;uint16&gt; | flattened per-read candidate **on-ladder** lengths (non-negative) |
| `amb_logliks` | list&lt;float32&gt; | flattened **stutter-free** log-liks (parallel) |
| `offl_seqs` | list&lt;string&gt; (dict) | distinct **off-ladder** allele sequences at this locus, as normalized deltas vs the ref tract, canonical order; empty when none |
| `offl_counts` | list&lt;int32&gt; | confident off-ladder read count per sequence (parallel to `offl_seqs`) |
| `offl_weight` | list&lt;float32&gt; | optional base-qual aggregate per off-ladder sequence (parallel) |
| `offl_amb_offsets` | list&lt;int32&gt; | CSR offsets for ambiguous reads carrying off-ladder mass (len = n_offl_amb + 1) |
| `offl_amb_idx` | list&lt;int32&gt; | flattened index into `offl_seqs` per ambiguous-read off-ladder candidate |
| `offl_amb_logliks` | list&lt;float32&gt; | flattened **stutter-free** log-liks (parallel) |

The `offl_*` block mirrors the on-ladder block exactly — `offl_seqs`/`offl_counts`
↔ `hist_lengths`/`hist_counts` (confident tally) and `offl_amb_*` ↔ `amb_*`
(per-read sparse profile) — the only difference being that an off-ladder allele is
keyed by a **normalized sequence** rather than an integer rung. The exact column
split is an implementation detail; the invariant is that every read's likelihood is
expressible over {its on-ladder window rungs} ∪ {the off-ladder sequences it
supports}.

Confident on-ladder reads collapse to the `(length → count)` histogram; genuinely
bimodal reads carry an explicit sparse `Qᵣ(L)` profile in the `amb_*` CSR columns
(in v1, not deferred); off-ladder alleles are stored symmetrically in the `offl_*`
columns and are empty at the vast majority of loci. **TOML header** (human-readable;
the container's, not a binary footer): `kind = "ssr"`, `container_version`,
`schema_version`, `sample_name`, `reference_md5`, `catalog_md5`, `ploidy`,
`extraction_params{MIN_FLANK_BP, MIN_BASE_QUAL, FAST_PATH_GATE,
STUTTER_WINDOW_UNITS, AMB_LL_DROP}`, `tool_version`, and a contig name↔id table.
`kind`/`schema_version` are authoritative (the filename is convenience);
`container_version` versions the shared format **independently** of the SSR schema,
so it and the `.snp.psp` schema evolve without lockstep. Write knobs:
`--block-window-bp` (the decode unit and the primary memory lever, as on the SNP
path); per-column ZSTD. All stored likelihoods are **stutter-free**; on-ladder
lengths are `uint16` repeat units (non-negative; the only signed length
quantity, the ref offset Δ, is derived, never stored — §5.5).

---

## 5. Stage 2 — Cohort genotyping

From the per-sample evidence, jointly genotype the cohort with population-informed
posteriors, biased toward precision.

### 5.1 Generative model (per locus ℓ)

1. **Population allele frequencies `π = (π_a)` over the candidate allele set.** At
   this locus an allele is a **sequence** — in practice a repeat length on the
   integer ladder (10 units, 11, 12, …) for the on-ladder majority, plus any
   **off-ladder** sequences the cohort emitted (§4.2). `π_a` is the fraction of
   chromosomes in the cohort carrying allele `a`, so `π` is a probability
   distribution over the candidate set (`π_a ≥ 0`, `Σ_a π_a = 1`). It is one of the
   two unknowns Stage-2 EM estimates (§5.4) and the root of the hierarchy: `π` →
   genotype (2) → reads (3). The candidate set is assembled before EM — see below.
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
The sum runs over the molecule lengths `L` reachable from allele `a` by whole-unit
slippage, *in `a`'s own frame*: if `a` is on-ladder, those `L` are on-ladder rungs
(`amb_*`); if `a` is off-ladder, they are the off-ladder sequences sharing `a`'s
remainder (`offl_*`). An off-ladder candidate a read has no support for contributes
negligibly to that read's likelihood automatically — no special-casing.

**Candidate allele set `A_ℓ` — the cohort union (assembled before EM).** The
per-locus candidate set is built by the same synchronized N-file scan the cohort
stage already does, from the **cohort-aggregate** on-ladder length profile (per-rung
sum of all samples' `hist_*`/`amb_*` support) plus the off-ladder sequences:

```
A_ℓ  =  { on-ladder rungs that are a LOCAL MAXIMUM of the cohort-aggregate profile
          (support ≥ both neighbours),  OR  adjacent (±1 unit) to such a maximum }
     ∪  { supported off-ladder sequences,
          unioned across samples by their normalized key (§4.2) }
```

Only evidenced rungs are eligible; a rung no read anywhere supports is never a
candidate (there is no assembly — the precision-first sensitivity limit of §4.2).
The off-ladder union by normalized key makes one biological off-ladder allele one
candidate, not one-per-sample; off-ladder alleles are kept whenever supported (they
cannot be in-frame stutter of an on-ladder allele, so they are their own peaks).
This is the Stage-2 analogue of HipSTR's cross-sample candidate assembly, by
sequence-key union rather than pooling raw reads — forbidden (§5.2) because it would
confuse one sample's stutter ladder with population length variation. `π`, the
genotype space (integer partitions of ploidy over `A_ℓ`), and the convolution above
are all defined over this `A_ℓ`.

**The peak + adjacent rule — pruning only the *unambiguous* stutter.** A length
profile is peaks (true alleles) on decaying stutter shoulders. The rule keeps every
peak and every rung within one unit of a peak, pruning the rest — *without ever
having to tell a real allele from a stutter satellite*, which is impossible from the
profile shape alone (a skewed adjacent het `a`/`a+1` puts a real allele at `a+1`
that is **not** a local maximum — topologically identical to a +1 stutter satellite).
The two clauses divide the labour:

- **Local-maximum clause** keeps a real allele ≥2 rungs from a bigger one — it forms
  its *own* peak (a real `a+2` is above `a+1` and `a+3`).
- **±1-adjacent clause** keeps the one case the peak clause structurally cannot: a
  real allele exactly one repeat from a bigger one. It is kept **unconditionally**,
  and whether it is a rare real allele or a stutter satellite is left to the EM,
  which *can* resolve it (shared kernel + population recurrence, §5.4) where the
  profile shape cannot.
- **Everything pruned** is then a non-peak rung ≥2 from every peak — an
  *unambiguous* stutter tail or inter-peak valley, where a real allele could not live
  without either peaking or being adjacent to a peak.

Crucially, pruning a rung drops it as a genotype *candidate* but **not its reads**: a
read at a pruned `a+2` rung still enters `P(read | a) = Σ_L Qᵣ(L)·S_θ(L | a)` as a +2
slip of the kept peak, so the stutter kernel still sees the full ladder — only the
disbelieved genotype hypothesis is removed. The single real-allele class this can
lose is one that is *both* not a peak *and* ≥2 rungs from every peak anywhere in the
cohort — permanently buried in a stutter tail at the sensitivity floor: a small,
documented precision-first recall trade (an extension of the §4.2 no-assembly limit).
Peak detection runs on the cohort-*aggregate* profile so **population recurrence**
protects real alleles; an allele carried by very few samples is the sensitivity edge
(admitting per-sample peaks is the escape hatch if data show it matters — deferred).

**Bounding cost — hard cap, precision-first.** The genotype space is the integer
partitions of ploidy over `|A_ℓ|` (§5.7), which blows up at a hypervariable locus ×
high ploidy — and there the many candidates are genuine *peaks*, so the peak rule
does not shrink them. The backstop is a hard cap `MAX_CANDIDATE_ALLELES`: if `|A_ℓ|`
exceeds it, **no-call the locus and log it** (never a silent top-K truncation — a
genotype over a clipped allele set is not trustworthy, and hypervariable loci are
where genotyping is least reliable anyway).

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

**Only in-frame — and what a non-multiple change means.** The kernel moves
probability *only* in whole-unit steps, because slippage re-registers on the repeat
period (HipSTR additionally carries a small *out-of-frame*, bp-step term; we fold
that into the cases below rather than model it as stutter — a deliberate
simplification, see §12). A length change that is **not** a multiple of the motif —
e.g. a 1 bp indel inside a trinucleotide tract — is **not** stutter, but it splits
into two distinct things that must not be confused:

- **A per-read, non-recurrent non-multiple change is a sequencing/replication base
  error.** It is already accounted for by the Stage-1 `Qᵣ` term (§4.2) and never
  becomes an allele.
- **A recurrent non-multiple change is a real, heritable off-ladder allele.** It is
  its own candidate in `A_ℓ` (§5.1), with its own in-frame stutter ladder in its own
  frame. It is *not* a stutter satellite of the nearest on-ladder allele.

The discriminator is **recurrence**: an error appears sporadically in single reads;
an off-ladder allele recurs across reads and samples at a stable frequency (the same
population-vs-stutter logic the identifiability argument uses, §5.4). Splitting them
this way keeps each phenomenon in exactly one place: transient error in Stage-1 `Qᵣ`,
in-frame slippage in `S_θ`, and heritable off-ladder variation in `π` over `A_ℓ`.

**`θ` is covariate-parameterised**, predicted per allele:
`θ = θ(period, allele_length, motif, purity)`. This is biology, not bookkeeping —
stutter rate is not constant: it rises sharply with **allele length** (longer
tracts slip more), is highest for **short periods** (mononucleotide runs are the
worst, hexamers mild), varies by **motif** (`AT` vs `GC` content), and falls with
**purity** (interruptions arrest slippage). So `(u, d, ρ)` are *functions* of these
covariates rather than one global triple.

**Pooling is a deliberate tradeoff, strongest for *pure* loci.** Three prior tools
bracket the design space (verified against the vendored sources): **GangSTR** uses
one fixed global triple (`0.05, 0.05, 0.9`) with an optional per-locus override
table — robust but blind to real variation; **HipSTR** estimates a *separate* triple
per locus by EM — captures everything locus-specific but is data-hungry and noisy at
low depth; **popSTR2** (our approach) regresses the triple on covariates and pools
across loci — borrows strength, robust at low depth, but **only as good as the
covariates**. We bet the four covariates capture stutter well for **pure** loci
(period + length + motif/GC essentially determine slippage there). The bet is
weakest exactly where the covariate is weakest — **impure / interrupted loci**: a
single scalar `purity` cannot encode *interruption structure* (one mid-tract
interruption vs scattered mismatches vs a near-compound tract all read as the same
`purity` yet arrest slippage differently). So impure loci are the regime to watch —
hence the per-locus relief valve and the open data question below.

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

**Fitting the covariate kernel — proposed model (provisional; the exact form is the
one genuinely open core of Stage 2, settled on real data — not pre-judged here).**

*The invariant — valid-probability links.* However the covariate dependence is
shaped, the three outcomes (no-slip / up / down) must always form a valid split and
`ρ` must stay in `(0,1)`. So the parameters are modelled through link functions: a
**multinomial-logit (softmax)** over `(1−u−d, u, d)` and a **logit** for `ρ`. Every
prediction is then a proper distribution regardless of the covariates.

*The structure — two options, choice deferred to data:*
- **Option 1 (v1) — discretized cells + empirical-Bayes shrinkage.** Bin the
  covariate space (period × `log(allele_length)` bin × motif-or-GC class × purity
  class), estimate `(u,d,ρ)` per cell by the M-step estimators, and **shrink sparse
  cells toward their period-level parent** (James–Stein / hierarchical) so a thin
  cell inherits the period rate instead of overfitting. Simple, transparent,
  trivially fits inside EM (responsibility-weighted counts per cell + shrinkage),
  validates per-cell. Con: bin-edge artifacts, no smooth length extrapolation.
- **Option 2 (upgrade) — continuous hierarchical GLM.** The softmax/logit above as
  *smooth* functions of `log(allele_length)` (period, motif/GC as factors, purity
  continuous), with **partial pooling** via a hierarchical prior on the
  coefficients. Smooth, extrapolates, no arbitrary bins. Con: a constrained weighted
  GLM re-fit each EM iteration — more machinery and robustness care.

Recommendation: **Option 1 as v1** (it *is* the "cells shrink toward period-level"
idea, and the cheapest honest thing that works), **Option 2 as the documented
upgrade** if cell artifacts prove material — mirroring the v1/upgrade pattern of the
base measure (§5.5). Reads pool across all (sample, locus) units in a cell; only the
*kernel parameters* pool, never the allele calls (above).

**Per-locus relief valve + impure-locus policy.** Pooling is the default, but a
**high-depth locus that clearly disagrees with its covariate cell** may carry a
per-locus deviation on top of the cell mean (partial pooling at the locus level —
HipSTR-like where the data support it). This matters most for **impure loci**, where
the policy is *per-locus when depth allows, else **no-call*** — we do not force an
impure locus onto a covariate-cell rate we don't trust (precision-first). A cell or
locus with neither enough data nor a trustworthy pooled rate falls back to the
**fixed default kernel** at the bottom of the §5.4 seed hierarchy (GangSTR's
`0.05/0.05/0.9`, ideally a per-period default table).

⚠ **Open — settle on real data before committing.** How well covariate-pooling holds
across **repeat size, GC content, and purity/impurity** is empirical and not
pre-judged: measure per-cell stutter fit and residual locus-to-locus variation on
real cohorts, then decide (a) Option 1 vs 2, (b) the impure-locus per-locus/no-call
depth cutoff, and (c) whether a per-locus relief term is needed at all. This is the
statistical core; §11 lists it as open, not settled.

### 5.3 Genotype prior — inbreeding-adjusted (Wright / HWE-with-`F`)

Plants are far from HWE (tomato is highly selfing), so the prior is **not** strict
HWE. **Reuse the SNP caller's exact multiallelic IBD-mixture prior**
([`src/var_calling/posterior_engine.rs`](../../../src/var_calling/posterior_engine.rs))
with a **per-sample fixation index `F ∈ [0,1]`** (`0` = outcrossing/HWE, `1` =
full inbreeding). With probability `F` the individual is **fully autozygous** (all
`ploidy` copies IBD from one allele); else its `ploidy` copies are independent HWE
draws from `π`. For a genotype `G` with allele counts `(k_a)`, `Σ_a k_a = ploidy`:

- **fully homozygous** (one allele `i`, `k_i = ploidy`):
  `P(G | π, F) = F·π_i + (1 − F)·π_i^ploidy`
- **≥2 distinct alleles** (incl. partial polyploid hets like `AABB`, `AAAB`):
  `P(G | π, F) = (1 − F)·(ploidy! / Π_a k_a!)·Π_a π_a^{k_a}` — the multinomial HWE
  term, with **no** `F` contribution (here `F` means *all* copies IBD, which forces
  full homozygosity).

Diploid is the familiar special case: `P(ii) = F·π_i + (1−F)·π_i²`,
`P(ij) = (1−F)·2·π_i·π_j`. This is exactly what the SNP engine computes (verified: it
enumerates all non-decreasing `ploidy`-tuples and uses the `π^ploidy` homozygous
term), reused verbatim. It is the **single-parameter polysomic** approximation: `F`
captures only full autozygosity, so **partial IBD** (e.g. 2 of 4 tetraploid copies
IBD) and **double reduction** are ignored — the refinement deferred in §5.7/§12.

STRs differ from SNPs only in that the allele set is repeat alleles (on-ladder
lengths plus off-ladder sequences, §5.1) and is highly multiallelic. High `F` in a
selfing line a-priori down-weights heterozygotes — a *second* structural defence
against stutter-driven false hets.

### 5.4 Cohort EM (per locus; parameters `π`, `θ`; latent `G_s`)

- **E-step:** `γ_s(G) = P(G | reads_s, π, θ, F_s) ∝ P(reads_s | G, θ)·P(G | π, F_s)`.
- **M-step (allele freqs):** `π_a ← (1/(ploidy·N))·Σ_s Σ_G γ_s(G)·count_a(G)`
  (standard AF-EM → AC/AN/AF). **Dirichlet-smoothed** with pseudocount `α` and a
  **non-uniform base measure** (§5.5) so plausible rare alleles are never
  hard-zeroed. This is also where the **kept-but-ambiguous candidates self-prune**:
  a kept ±1 shoulder (§5.1) whose apparent reads are all assigned to its peak via
  `S_θ` converges to `π_a → 0` and falls below `MIN_ALLELE_SUPPORT_FRAC` — the
  *model* makes that call, not a pre-EM heuristic (§5.1 already excludes the
  *unambiguous* stutter — non-peak rungs ≥2 from any peak; the ±1 shoulders it keeps
  are the ambiguous ones the EM resolves here).
- **M-step (stutter):** update the shared covariate kernel from
  responsibility-weighted observed-vs-assigned length residuals (the deconvolution
  update), following HipSTR's estimators — but note the normalisation: `u` and `d`
  are each the posterior-weighted count of gain / loss events divided by **all**
  events (no-slip included), not longer-vs-shorter among slips only; `ρ` is the
  geometric MLE of the step-size decay. (`em_stutter_genotyper.cpp` normalises over
  the full event total.)
- **M-step (`F`, optional):** estimate a cohort/per-sample `f̂` from the
  heterozygote deficit only above an effective-sample threshold; else use the
  supplied default. Per-locus `F` is not estimated.

**Structure — initialisation by a confident-homozygote pre-pass.** EM is
non-convex and the stutter kernel and the genotypes are mutually dependent (learn
the kernel ⇄ know which reads are stutter), so a good *starting* kernel matters —
but only as a start: the M-step re-estimates `(u,d,ρ)` from all loci every
iteration, so the seed shifts **which basin EM starts in**, not where it converges.
The pre-pass exploits the one case where the genotype is *nearly* knowable without a
stutter model: a **confident homozygote**, selected by two named gates — spanning
depth ≥ `SEED_MIN_DEPTH` (the length histogram is well-resolved) **and** a single
dominant length holding fraction ≥ `SEED_MIN_DOMINANCE_FRAC` (~0.9) of spanning
reads, the rest a small ladder around it. At such a (sample, locus) the modal
length is taken as the single true allele `a`, so each non-modal read is treated as
a stutter product with known step `δ = (L − a)/period` — a labelled stutter
observation, no genotype inference. Pooling these across many confident-homozygote
loci within a covariate cell yields an initial `(u, d, ρ)` via the M-step
estimators.

*This is a seed, not ground truth.* The selection is not perfectly clean: a
**skewed adjacent-unit heterozygote** (true `a`/`a+1` with, say, 85/15 balance) can
pass the dominance gate and masquerade as a homozygote, feeding its real minor
allele into the stutter estimate. The resulting bias **inflates `u`/`d`** — it errs
toward *more* stutter, which suppresses hets: a **precision-safe** direction,
consistent with the mandate. `SEED_MIN_DOMINANCE_FRAC` is the primary control (at
0.9 it excludes the 85/15 masquerader outright); stricter → cleaner seed but fewer
seed loci, and genome-wide pooling supplies plenty, so lean strict. And because the
seed only sets the starting basin (above), the residual contamination is
low-stakes — its one real risk is starting EM in a bad basin, which the diagnostic
below flags and the deferred multi-restart ultimately backstops.

**Measured cut-off (per covariate cell).** Rather than fix the gate blindly, build a
**per-cell** distribution over the seed-candidate loci and detect the het/stutter
confusion directly (stutter rate varies by cell, so the cut-off must too). The
discriminating statistic is **not** raw dominance — true homozygotes and skewed
adjacent hets both sit at high dominance and overlap there — but the **one-sided
adjacent-neighbour excess**: how much the ±1-unit mass exceeds, and is more
lop-sided than, a plausible stutter ladder (a true homozygote's ladder is two-sided
and decays per `u/d`; a masquerader's is a one-sided spike at the het partner).
Bimodality in that statistic is the confusion signature: its valley **recommends**
`SEED_MIN_DOMINANCE_FRAC` for the cell, and its strength is a **warning** flagging
the regime where the otherwise-low-stakes seed could mislead the basin. Exposed as
`--seed-dominance auto` (apply the per-cell valley) vs a fixed conservative default.
⚠ **Open — to be validated on real data:** whether a clean per-cell cut-off can be
*reliably* detected is an empirical question (the bimodality may be weak or smeared
in practice). So `auto` ships **available but off by default**, and is promoted to
the default **only if** that validation is convincing.

**Seed fallback hierarchy** (so the kernel is always seedable, even where confident
homozygotes are scarce — e.g. hypervariable cells in a diverse outcrosser): a cell
with enough confident-homozygote observations uses its own seed; sparse cells
**shrink toward the period-level pooled seed**; cells with neither fall back to a
**fixed default kernel** (literature-derived per period). This removes any hard
dependence on confident homozygotes existing in every cell.

The joint deconvolution EM then **refines** the seed using *all* loci, now including
heterozygotes — where the two alleles' ladders overlap and the genotype is genuinely
latent. This homozygote init is the chosen defence against local optima (see
Identifiability below); multi-restart is deferred.
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

**Off-ladder alleles.** `G₀` above lives on the integer `Δ` ladder. An off-ladder
allele (§4.2) sits at a fractional offset, so it takes the `G₀` mass of its
**nearest on-ladder rung** (its integer repeat count) scaled by a small constant
`OFFLADDER_PRIOR_FACTOR` (< 1) — encoding that an off-ladder variant is a-priori
rarer than the clean rung at the same length. This keeps off-ladder alleles
representable in `π` without letting a thin off-ladder ladder masquerade as a major
allele at small N; like the rest of `G₀`, it is swamped by data at useful cohort
sizes.

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
empirical population AF (large N); `F` is supplied at low N (default 0). The kernel
seeds from a *single* sample **regardless of mating system**, because
**homozygosity is per-locus**: even an outcrossing individual is homozygous at a
large fraction of its loci, so its genome-wide confident-homozygote loci (§5.4)
seed the kernel without needing other samples (selfing species simply supply more).
Single-sample is, however, the **weakest identifiability regime**: with one sample
there is no **population recurrence** (§5.4) to separate a real rare allele from a
stutter satellite, so identifiability rests on the shared kernel + the base-measure
prior alone — and the precision-first no-call (§5.8) is the safety valve when that
is not enough.

### 5.7 Ploidy & `F` input

- **Ploidy:** a single uniform `--ploidy` per run (default `2`); the genotype space
  is the **integer partitions of that ploidy** over the candidate alleles `A_ℓ` (the
  ConSTRain enumeration trick). Its size grows steeply with `|A_ℓ|` × ploidy, so
  `MAX_CANDIDATE_ALLELES` (§5.1) is what bounds Stage-2 cost at hypervariable loci.
  No aneuploidy / per-locus CN / mixed ploidy. Polyploid supported (uniform);
  polyploid `F` uses a single-parameter polysomic approximation (ignores double
  reduction — a v-later refinement).
- **`F`:** per-sample `fixation_index_default` (default `0`) + optional per-sample
  overrides; optionally estimated (§5.4).

### 5.8 FP-aversion (precision ≫ recall)

- The **stutter model** explains the ladder as noise (first defence); the
  **inbreeding prior** suppresses spurious hets (second defence).
- **Known bias — heterozygosity runs low (document, don't hide).** Those two het
  defences systematically under-call **adjacent-unit heterozygotes** (true
  `a`/`a+1` → called `a`/`a`), so **observed heterozygosity Ho is biased downward** —
  a real cost on a headline popgen statistic, not a free precision win. This is a
  property of the **converged** kernel + these thresholds (distinct from, and not
  caused by, the §5.4 *seed*). It is **worst at small N** and recovers as cohort
  size grows, because **population recurrence** (§5.4) re-establishes a recurring
  `a+1` as a genuine allele and lets its hets be called. Reported, swept alongside
  the operating point, and checked in validation (§7).
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

- **REF/ALT = actual repeat-tract sequences** (not symbolic) — the allele's
  identity (§1); on- and off-ladder alleles round-trip losslessly (12 vs 12+1 bp
  are distinct ALT strings). Copy number is a **derived** annotation, never the
  identity: TRtools recomputes it as `len(allele)/len(motif)` and **tolerates
  fractional values** (verified in `tr_harmonizer.py` — it does not even read our
  `REPCN`), so an imperfect / off-ladder allele whose length is not a clean motif
  multiple is represented honestly via its sequence rather than rounded into a
  neighbour. We still emit per-allele **`REPCN`** (integer copies, rounded) for
  consumers that expect it, plus a per-ALT **`BPDIFFS`** (bp difference from REF,
  HipSTR-style) so exact bp-level detail survives.
- **Fields, by who actually needs them** (verified against the gangstr harmonizer
  path, not assumed):
  - *Harmonizer-mandatory* (the path raises without them): INFO `RU` (motif);
    FORMAT `GT`; and a `##FORMAT=<ID=Q,…>` header declaration for
    `HasQualityScore`. Must **not** emit INFO `VID` (AdVNTR marker) or `VARID` (EH
    marker) — the gangstr path raises on either.
  - *GangSTR-conventional, emitted for fidelity / other consumers* (TRtools neither
    reads nor requires them): INFO `END`, `PERIOD`, `REF` (reference copy number);
    FORMAT `DP`, `REPCN`. `Q` (genotype posterior, 0–1) is TRtools' quality field —
    read when its header is declared and the record is not Beagle-imputed.
- **Our model → GangSTR fields:** stutter `(u,d,ρ)` → INFO
  `STUTTERUP`/`STUTTERDOWN`/`STUTTERP`.
- **Our extras (TRtools ignores; carry our full posterior + Beagle compat):**
  FORMAT `GP`, `GL`/`PL`, `GQ`; INFO `AF`/`AC`/`AN`, `NS`, `F`. FILTER `PASS` +
  `lowGQ`/`lowDepth`/`segdup`/`lowSupport`. *(TRtools reads only the scalar `Q`, so
  `Q` must carry the called genotype's posterior.)*
- **Detection (honest, collision-free).** TRtools infers `gangstr` when the
  lowercased raw header contains **both** the substrings `command=` **and**
  `gangstr` — two tokens, anywhere, possibly on different lines
  (`tr_harmonizer.py:210`). We satisfy both truthfully: a real
  `##command=ssr-call …` line (our actual invocation supplies `command=`) plus a
  `##source=ssr-call <ver> (GangSTR-compatible output)` line (truthfully supplies
  `gangstr`). **No spoofed GangSTR command** — we never claim GangSTR ran. To stay a
  *unique* match (TRtools hard-errors on multi-match) we avoid the other detectors'
  triggers: no `hipstr`/`longtr` tokens anywhere, no `source=advntr`/`source=popstr`,
  no symbolic `##ALT=<ID=STR\d+>` (we emit real sequence ALTs regardless).
  - ⚠ `--vcftype gangstr` is **not** a bypass: inference still runs and the flag
    only succeeds if the header is *already* gangstr-detectable
    (`tr_harmonizer.py:237`). So the two header tokens are required either way — the
    flag only disambiguates, it cannot rescue a non-detectable header. No TRtools
    change needed.
- **Verified by a live test, not by assertion.** A test asserts
  `TRRecordHarmonizer(our.vcf)` infers `gangstr` *uniquely* and round-trips a record
  (REF/ALT-derived copy numbers, `Q`, `GT`) — guarding both our header construction
  and against TRtools changing its detection rules.
- **Ploidy caveat:** diploid is fully harmonizable; polyploid VCFs are valid but
  TRtools is diploid-centric.

---

## 6. Parameters (named constants — units & source)

Defaults are starting points; the precision-critical ones are swept (§7).

| constant | role | default / source |
|---|---|---|
| `MIN_FLANK_BP` | clean flank bases for a spanning read | ~10 bp (HipSTR) |
| `MIN_BASE_QUAL` | boundary base-quality for the fast path | ~Q20 |
| `FAST_PATH_GATE` | fast-path qualification **predicate** (cleanly-aligned flanks, 0 interior indels, base-qual) — a gate, not a tunable scalar | n/a (predicate) |
| `STUTTER_WINDOW_UNITS` (`W`) | Qᵣ(L) candidate window ±units | 3 |
| `ploidy` | uniform genotype size | 2 |
| `fixation_index_default` (`F`) | inbreeding prior | 0.0 (SNP caller) |
| `SEED_MIN_DEPTH` | min spanning depth for a kernel-seed homozygote (§5.4) | — |
| `SEED_MIN_DOMINANCE_FRAC` | min modal-length fraction for a seed homozygote; per-cell `auto` or fixed (§5.4) | ~0.9 |
| `alpha` (`α`) | Dirichlet pseudocount on `π` | — |
| base-measure decay `p` | reference-centred allele-length prior | — (SMM) |
| `OFFLADDER_PRIOR_FACTOR` | base-measure down-weight for off-ladder vs nearest rung | — (< 1) |
| `MAX_CANDIDATE_ALLELES` | cap on `|A_ℓ|`; exceed → no-call + log the locus (§5.1) | — |
| `DEFAULT_STUTTER_(U,D,RHO)` | fixed fallback kernel when a cell/locus lacks data (§5.2/§5.4) | 0.05 / 0.05 / 0.9 (GangSTR) |
| `MIN_PERLOCUS_STUTTER_DEPTH` | depth to trust a per-locus kernel over the cell (impure loci; §5.2) | — |
| `MIN_POSTERIOR` (`Q`) | call-quality gate | 0.9 (HipSTR) |
| `MIN_SPANNING_READS` | depth gate | ~10 (HipSTR) |
| `MAX_STUTTER_FRAC` | artifact-read gate | 0.10 (HipSTR) |
| `MIN_ALLELE_SUPPORT_FRAC` | per-allele support | 0.20 (HipSTR) |

---

## 7. Validation & testing (two buckets)

**Bucket 1 — synthetic data for code tests (build critical path).** A bespoke
STR-aware simulator emitting at two levels: **evidence-level** (synthetic
`.ssr.psp` `Qᵣ(L)` → tests Stage 2 statistics in isolation, before extraction
exists) and
**read/BAM-level** (synthetic BAM + truth genotype table → tests the whole
pipeline incl. the pair-HMM). **Anti-tautology rule:** generate data whose
stutter/error *deviates* from the caller's `(u,d,ρ)` assumptions, and keep the
simulator's model definition separate from the caller's code. Tests: deterministic
unit fixtures + property/statistical tests (recovers injected `π/θ/F`; calibrated
posteriors; no false hets on monomorphic loci; di-nuc / long-allele / low-coverage
/ polyploid / **off-ladder-allele** (e.g. a 12 vs 12+1 bp cohort: must call two
distinct alleles, not collapse them, and reconcile the off-ladder allele across
samples) / **imperfect-locus** stress). Plus an **Ho/He-vs-cohort-size**
calibration: quantify the adjacent-unit-het under-call (§5.8) and confirm it shrinks
with N as population recurrence kicks in; and a **seed-cut-off detection** check —
on data with injected adjacent-het contamination, does the per-cell bimodality
diagnostic (§5.4) recover a clean cut-off? (gates whether `--seed-dominance auto`
can become the default). And a **covariate-kernel-fit** study (the §5.2 open core):
on real cohorts, measure per-cell stutter fit and residual locus-to-locus variation
across repeat size / GC / purity — does covariate-pooling hold for pure loci, and
how badly does it fail on impure ones? — to settle Option 1 vs 2 and the
impure-locus per-locus/no-call cutoff.

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

## 8. CLI surface

Three subcommands on the existing binary, named to mirror the SNP caller's
`pileup → var-calling` roles one-for-one (architecture doc §2/§3.3):

- `ssr-catalog` — reference FASTA → catalog (Stage 0). The only stage that reads
  the FASTA; it embeds the local reference into the catalog (`ref_seq`, §3.2).
- `ssr-pileup` — BAM/CRAM + catalog → per-sample evidence `.ssr.psp` (Stage 1;
  the SSR analog of `pileup`). Takes **no reference for the SSR algorithm** (it
  reads `ref_seq` from the catalog); a reference is needed **only to decode CRAM
  input**, not for BAM (§3.2 CRAM caveat).
- `ssr-call` — N evidence files + catalog → cohort VCF (Stage 2; the SSR
  analog of `var-calling`). Knobs include `--seed-dominance {fixed|auto}`
  (§5.4; `fixed` default until the cut-off detection is validated on real
  data).

**The simulator is not a subcommand.** Injecting genotypes+stutter →
synthetic BAM and/or evidence + truth table is **test/dev scaffolding**, so it
lives as a **crate module the test suite calls directly** (`src/ssr/simulate/`),
not on the production CLI (architecture doc §2.1).

Decision E (repo/crate placement, the shared `.psp` container split, CLI
names) is **settled** in the architecture doc: same crate / new `src/ssr/`
module tree; a generic schema-agnostic `psp` core + `snp`/`ssr` schemas (flat
in `src/psp/`, submodules deferred); the three subcommand names above.

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

## 11. Structural decision (E) — SETTLED

Repo/crate placement, CLI naming, **and where the shared `.psp`
columnar-block container lives** are settled in the architecture doc
([`../architecture/ssr_genotyping_architecture.md`](../architecture/ssr_genotyping_architecture.md),
§2–§5, 2026-06-11):

- **Placement:** same crate, new `src/ssr/` module tree, `ssr-*` subcommands on
  the existing binary (not a new repo, not a workspace split).
- **CLI:** three subcommands — `ssr-catalog`, `ssr-pileup`, `ssr-call` — named
  to mirror the SNP `pileup`/`var-calling` roles. The simulator is **not** a
  subcommand; it is a crate/test module (§8).
- **Container:** a generic schema-agnostic core plus two thin schemas (`snp`,
  `ssr`) in `src/psp/` — both `.snp.psp` and `.ssr.psp` ride the same core. The
  three-layer split is **conceptual, enforced by file naming**; the directory
  **stays flat** (~13–14 files) and grows `core/`/`snp/`/`ssr/` submodules only
  on a crowding trigger (architecture doc §3.2). Extracted *from the two
  concrete consumers* as SSR is built, not designed up front; **pre-alpha, so no
  backwards-compat** — free to rev the format and restructure the SNP schema,
  the regression gate being the SNP caller's end-to-end tests, not byte-identity
  to the old format.

With E settled, each stage (0, 1, 2) and the simulator becomes its own
implementation plan, in data-flow order, with Bucket-1 synthetic tests on the
critical path.

**Two genuinely-open *modelling* cores** (distinct from E, which is structural) —
both to be settled on real data during the Stage-2 plan, not pre-judged here:
- **The covariate stutter-kernel form (§5.2):** Option 1 (cells + shrinkage) vs
  Option 2 (continuous GLM); the impure-locus per-locus/no-call policy; whether a
  per-locus relief term is needed. This is the statistical core of the caller.
- **`--seed-dominance auto` (§5.4):** whether the per-cell cut-off can be reliably
  detected.

"Settled" elsewhere means the architecture and the generative model; several model
defaults remain to be swept — see the §6 parameter table and §7.

---

## 12. Deferred / future

- Beyond-read-length classes (GangSTR Geometric/Normal-insert/Poisson-FRR/
  Uniform-flank) + read-pair merging.
- Polyploid `F` with partial IBD / double reduction.
- Multi-restart for hard-locus local optima; a second stutter-refinement round
  (full joint hierarchical EM) if measurements demand it.
- **Out-of-frame stutter** as its own model term. We currently fold non-unit-multiple
  slippage into sequencing error (Stage-1 `Qᵣ`) or a heritable off-ladder allele
  (§5.2) — a deliberate simplification. HipSTR carries a separate bp-step out-of-frame
  term; add one here if validation shows genuine out-of-frame slippage being
  misattributed.
- Discretized-Gaussian (θ-linked) allele-length base measure if geometric proves
  insufficient.
