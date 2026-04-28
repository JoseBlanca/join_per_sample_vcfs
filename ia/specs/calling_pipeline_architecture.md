# Multi-sample calling pipeline — proposed architecture

**Status:** Design proposal, 2026-04-22. Intended to replace the gVCF-in /
VCF-out design described in `general_project_specification.md` and to guide
the rewrite. Supersedes and updates
[snp_calling_from_bam.md](../feature_implementation_plans/snp_calling_from_bam.md)
and [per_sample_pileup_format.md](per_sample_pileup_format.md) — both of which
anticipated this design but left key data-flow questions open.

## Design principles

A few commitments shape both this specification and the code that will
implement it. They take precedence when a concrete design decision is
ambiguous.

1. **Clarity and readability are paramount.** They apply equally to
   this specification and to the code that implements it — not one more
   than the other. A reader new to the project should be able to
   understand both without needing to consult additional sources.
   Cleverness that hurts readability is never a good trade.
   - *In the spec:* where a choice is not obvious, the reason is spelled
     out; where a trade-off exists, both sides are named.
   - *In the code:* identifiers name what they are; control flow reads
     top to bottom; abstractions are introduced only when they earn
     their complexity; decisions that can't be understood from the
     names alone get a short comment explaining *why* (never *what* —
     the code says that already).

2. **No guessing at user intent — every default and decision is
   explicit.** Every parameter default, every fallback behaviour, every
   assumption about the input is named in this document. Hidden
   heuristics and "sensible" values that only exist in the source
   accumulate into systems nobody can debug later. If a behaviour isn't
   written here, the implementation should not silently invent it.

3. **Errors must not pass silently.** Corrupt inputs, unexpected
   formats, out-of-range values, missing files — each produces a clear,
   actionable error message and halts. The default posture of every
   component is "reject and report," never "guess and proceed." It is
   far better to fail loudly at Stage 1 than to emit subtly wrong
   genotypes at the final stage.

## Constraints this design has to satisfy

1. **Scale to thousands of samples.** freebayes' joint-BAM approach does not
   — memory grows with cohort size. GATK's `GenotypeGVCFs` does, but it only
   joins, it does not truly merge.
2. **Merge genotypes across samples**, not merely join per-site records. In
   particular, correctly reconcile overlapping / complex events (a deletion in
   one sample that overlaps a SNP in another) into unified haplotype alleles.
3. **Read each BAM exactly once.** Per-sample analysis must produce a reusable
   artefact.
4. **Support proper Bayesian posteriors** (GT + GQ + QUAL) on the merged
   allele set, not just per-sample-local posteriors carried forward naively.
5. **Cohort re-callable without re-doing per-sample work.** Adding a new
   sample should not require re-running the per-sample stage for existing
   samples.

See [freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md) and
[posterior_gt_probs.md](posterior_gt_probs.md) for the two approaches this
design draws on.

## Pipeline overview

```
BAM_i  ──► Stage 1: per-sample caller ──► sample_i.psf (Stage 2 format) ┐
                                                                          │
BAM_j  ──► Stage 1: per-sample caller ──► sample_j.psf (Stage 2 format) ├─► multi-way per-position iterator
                                                                          │             │
BAM_k  ──► Stage 1: per-sample caller ──► sample_k.psf (Stage 2 format) ┘             ▼
                                                                           Stage 3: DUST filter ◄── reference FASTA
                                                                                        │
                                                                                        ▼
                                                                           Stage 4: grouping ──► OverlappingVarGroup stream
                                                                                                   │
                                                                                                   ▼
                                                                           Stage 5: per-group processing (rayon parallel)
                                                                                                   │
                                                                                                   ▼
                                                                           Stage 6: posterior engine ──► multi-sample VCF
```

Six stages total:

- **Stage 1 — per-sample caller.** Runs once per sample, reads the BAM,
  writes one `.psf` (per-sample file) artefact.
- **Stage 2 — per-sample file contract.** Specifies the `.psf` format
  that Stage 1 produces. Not a runtime step — it is the interface
  between the per-sample and cohort sides of the pipeline.
- **Stage 3 — low-complexity filter.** A streaming per-position filter.
  Computes DUST scores from the reference on the fly and silently drops
  per-position records whose reference context is low-complexity. No
  intermediate mask file is written.
- **Stage 4 — grouping.** Walks the filtered per-position stream in
  genomic order and bundles overlapping / adjacent candidate positions
  (e.g. a deletion spanning several positions, a SNP falling inside a
  deletion) into an `OverlappingVarGroup`. Sequential by nature, but
  emits a stream of independent groups that can be processed in
  parallel downstream.
- **Stage 5 — per-group processing.** Takes one `OverlappingVarGroup`
  at a time and produces a merged multi-sample variant record: allele
  unification, per-sample likelihood reconstruction from the stored
  scalars, compound-haplotype consistency via phase chain ids. Groups
  are independent of each other, so Stage 5 runs in parallel across
  groups via rayon.
- **Stage 6 — posterior engine.** Runs EM over the merged records and
  emits the final VCF.

Only Stage 1 touches BAMs. Stages 3–6 operate on the reference FASTA
and on the `.psf` files alone, which is what makes cohort recall cheap
and BAM-free. Only Stage 4 is inherently sequential; Stages 1 and 5
parallelise cleanly (across samples and across groups respectively).

## Stage 1 — per-sample caller

### Purpose

One BAM → one `.psf`. Encodes everything downstream stages might need to
compute a likelihood for this sample against any allele that may later be
introduced into the merged variant set by other samples.

Stage 1 is built in-tool rather than wrapped around an existing
per-sample caller (freebayes, bcftools mpileup). The output format
is the custom `.psf` (see Stage 2), which no existing caller
emits; the quality adjustment we want is BAQ applied inline and
parallelised across reads (see §"Why BAQ in-process" below), which
`samtools calmd` cannot do in-process; and the per-read likelihood
model we need is the simple freebayes-shaped one (Option B in
§"Per-read likelihood quality" below) rather than HaplotypeCaller's
PairHMM. Wrapping an existing caller would force us to write the
same BAM-walking + scalar-accumulation code anyway, and would add
an external dependency and an IO round-trip for no saved work.

### What it does, conceptually

1. Streams the BAM, position by position, accumulating per-sample per-read
   evidence at each reference locus.
2. **Applies BAQ (Base Alignment Quality)** per read, computing the
   posterior probability that each base is correctly aligned and capping
   effective base quality at `min(BQ, BAQ)`. This is done in-process,
   parallelised across reads — **not** by shelling out to `samtools
   calmd`. See "Why BAQ in-process" below.
3. Emits **one record per covered reference position**, pileup-style.
   Each record carries the five freebayes-sufficient scalars for every
   distinct allele observed at that position, computed from BAQ-adjusted
   qualities. When all reads at a position agree on REF, the record has
   one allele (REF); as soon as any read disagrees, the record has two
   or more. Stage 1 does not distinguish "variant" from "non-variant"
   positions — that decision only happens at the joint stage. Phase
   chain identifiers link alleles that co-occur on the same haplotype
   across positions. No per-genotype PLs are stored — they are fully
   derivable from the scalars when the joint stage needs them; see
   Stage 2 §"Why PLs are *not* stored".
4. **Indel alleles are anchored at the position one before the variable
   region** (VCF convention). A deletion is stored as one allele at the
   anchor position; interior positions of the deletion appear in the
   `.psf` only if they are covered by non-deleting reads. See Stage 2
   §"Indel anchoring" for the exact REF/ALT shape.
5. **Indel BQ proxy.** A deleted base has no sequencing quality of its
   own, so for a supporting read we take the window of `l + 2`
   BAQ-adjusted base qualities centred on the indel (where `l` is the
   indel length, with edge-clamping at read ends as in freebayes
   [AlleleParser.cpp:1626](../../freebayes/src/AlleleParser.cpp#L1626))
   and use the **minimum** BQ in that window as the per-read
   contribution to the indel allele's `Σ max(ln_BQ, ln_MQ)` scalar.
   This matches freebayes' `--useMinIndelQuality` mode. Reads reporting
   an indel as their first or last CIGAR operation are rejected — the
   lack of flanking evidence on both sides makes the placement
   untrustworthy (same rule as freebayes).
6. **Sub-threshold observations are kept.** A single read supporting an
   alt at an otherwise-REF position is recorded as a second allele
   entry for that position, not filtered out. At low coverage the
   cohort-level signal often lives in these single-read observations.
7. **Uncovered positions produce no record.** Positions where the
   sample has zero reads are simply absent from the `.psf`. Records
   are separated by a `delta_pos` field (distance from the previous
   record), so gaps are represented implicitly and the absence of a
   record unambiguously means "no data."

### Per-read likelihood quality: BAQ vs. PairHMM vs. trust-the-aligner

Three options exist for how trustworthy per-read observations should
be made before they enter the likelihood. This design picks the middle
option (BAQ) and rejects the others for the reasons below.

**Option A — trust the aligner.** Take reads as placed by the aligner,
use the BAM's reported base qualities (possibly post-BQSR) without any
local re-alignment or BAQ. This is essentially freebayes' approach.

- *Rejected because* it relies entirely on upstream BAM processing being
  clean. The single biggest known failure mode — false SNP calls
  adjacent to indels due to aligner misplacement — is unaddressed. At
  low-to-moderate coverage the resulting false positives are expensive
  to filter out later.

**Option B — BAQ (chosen).** Run a local HMM per read that computes
`P(base is correctly aligned)` per base; cap effective BQ at `min(BQ,
BAQ)`. Then compute per-read likelihoods in freebayes' style.

- Fixes the indel-adjacent misplacement signal at modest cost.
- Per-read and embarrassingly parallel — unlike `samtools calmd`,
  which is single-threaded and slow enough to be a pipeline bottleneck
  in practice. Implementing BAQ directly inside our per-sample caller
  lets us parallelise it across reads natively (rayon), amortise it
  with BAM decoding, and avoid the calmd roundtrip.
- Keeps the per-read likelihood simple: BAQ-adjusted BQ flows into
  `max(ln_BQ, ln_MQ)` as it would in freebayes.
- Matches what bcftools mpileup / bcftools call do by default, so the
  calibration is well-studied.

**Option C — local reassembly + PairHMM (HaplotypeCaller-style).**
Assemble candidate haplotypes from reads in an active window, then
score each read against each haplotype with a pair HMM.

- *Rejected because*:
  - At the target coverage (2–10× per sample) and one sample at a
    time, there are not enough reads per haplotype for reassembly to
    propose reliable novel haplotypes. With 3–5 reads covering a
    region, most "reassembled" haplotypes are supported by one or
    zero reads and only add noise.
  - Reassembly really pays off when reads from **many samples** can
    be pooled to build consensus haplotypes — that's a joint
    operation that does not fit the per-sample Stage 1 topology.
    Doing it well is essentially a separate project (cross-sample
    haplotype reassembly + PairHMM scoring).
  - Implementation cost is large: de Bruijn graph assembly + PairHMM
    is a multi-month engineering effort with significant correctness
    surface.
  - The per-sample work would be ~10–30× slower than Option B.

A future extension could add cross-sample reassembly as an optional
post-processing step, leaving Stage 1 as BAQ-based.

### Why BAQ in-process rather than requiring `samtools calmd -r`

`samtools calmd` is single-threaded and does not parallelise across
reads. In practice it becomes a bottleneck at scale — a per-sample
preprocessing cost that has to be paid serially before our tool can
start. Since the BAQ algorithm itself is embarrassingly parallel (a
local HMM per read, independent of every other read), implementing it
inside our per-sample caller lets us:

- parallelise across reads with rayon;
- pipeline it with BAM decompression and CIGAR walking;
- skip writing a BAQ-decorated BAM back to disk (`calmd` writes a new
  BAM, our pipeline just carries the adjusted BQ in memory into the
  likelihood computation);
- have one tool to run per sample instead of two.

This is a modest implementation win (BAQ is ~a few hundred lines
following Heng Li's 2011 paper) that removes a real-world performance
bottleneck users have hit with `calmd` in the past.

### Related existing draft

The binary layout for the per-sample artefact is sketched in
[per_sample_pileup_format.md](per_sample_pileup_format.md) (PSP v1).
That draft stores per-read records (3 bytes each). This architecture
rejects per-read storage in favour of per-allele scalar summaries —
sufficient to reproduce freebayes' posterior calculation exactly (see
[freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md),
section "The data freebayes actually consumes") at constant size in
depth. Per-read storage would buy flexibility for alternative
per-read likelihood models (PairHMM-style, novel error models) that
this pipeline does not plan to support, at a cost linear in coverage
— a trade we explicitly decline. PSP v1 needs to be revised
accordingly; the BAQ-adjusted BQ values should also be reflected in
whatever quality summary the revised format records.

## Stage 2 — per-sample file contract (`.psf`)

The `.psf` is a stream of **per-position records**, pileup-style. There
are no separate record types for variant vs non-variant positions and no
gVCF-style reference-confidence blocks — every covered position gets
the same record structure. Uncovered positions produce no record; gaps
are represented implicitly via `delta_pos`.

The contents of a per-position record:

| Field | Purpose |
|---|---|
| `delta_pos` | distance from the previous record (variable-length encoded). Absence of a record means the sample has no coverage at that position. |
| chromosome id (when the chromosome changes) | locus identification |
| `n_alleles` | number of distinct alleles observed at this position. `1` when all reads at this position agree on REF; `≥2` as soon as any read disagrees — even a single read supporting an alt. Stage 1 does not make variant / non-variant decisions; it just records what the reads show. |
| per observed allele: allele sequence + the **five scalars** below | the single source of truth for per-sample evidence at this position |
| per observed allele: phase chain identifier | links alleles across positions that co-occur on the same haplotype in this sample; lets the merger reason about compound haplotypes without needing any extra sequence context |

**No explicit haplotype-context field.** A natural-looking
alternative is to store a "local consensus sequence" on each allele
— the surrounding bases the read carries around the variant — so
the merger can reconstruct compound haplotypes by sequence
alignment. The record above deliberately does not carry one,
because every piece of information such a field would provide is
already in the record or reachable from it: the reference sequence
is in the FASTA, the sample's alt alleles sit in the per-position
records themselves, and the linkage between alleles at different
positions (which is the only cross-position fact the merger
actually needs) comes from matching phase chain ids. Adding a
consensus-sequence field would duplicate data already present
elsewhere without enabling any reasoning the phase chain ids do
not already support.

### The five per-allele scalars

The per-local-allele observation summary is the fixed set of scalars
sufficient to reproduce freebayes' likelihood and observation-bias
priors exactly. Five scalars per allele per sample:

| Scalar | Drives |
|---|---|
| observation count (reads supporting this allele) | allele-balance multinomial in the likelihood ([DataLikelihood.cpp:155](../../freebayes/src/DataLikelihood.cpp#L155)) + `observations` field in the combo-level `AlleleCounter` |
| `Σ max(ln_BQ, ln_MQ)` over supporting reads | `prodQout` contribution in the likelihood ([DataLikelihood.cpp:34](../../freebayes/src/DataLikelihood.cpp#L34)) — the nonlinear per-read combination of base quality and mapping quality, pre-summed |
| forward-strand count | strand-bias prior ([Genotype.cpp:1560](../../freebayes/src/Genotype.cpp#L1560)) |
| placed-left count | read-placement bias prior ([Genotype.cpp:1561](../../freebayes/src/Genotype.cpp#L1561)) |
| placed-start count | fragment-position bias prior ([Genotype.cpp:1562](../../freebayes/src/Genotype.cpp#L1562)) |

Total: ~24 B per allele per sample (five 4-byte or 8-byte fields with a
small header), constant in depth. At 2 alleles per site this is ~50 B
per sample per site.

**Note on the quality scalar.** Storing `sum(BQ)` and `sum(MQ)`
*separately* does **not** work — freebayes combines BQ and MQ per read
via `max` in log space before aggregating, which is nonlinear in the
pair. The per-allele scalar therefore has to be the pre-combined sum
`Σ max(ln_BQ, ln_MQ)`, not two marginal sums. See
[freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md#compactness-implications--important-subtlety-about-bq-and-mq).

The BQ values that feed into this scalar are the **BAQ-adjusted** base
qualities produced by Stage 1 — i.e. `min(BQ, BAQ)` per base rather
than the BAM's raw BQ. See Stage 1 §"Per-read likelihood quality".

### Allele and record conventions

A few simple rules make every allele in a `.psf` record
self-describing. Anyone reading the file can identify the type of each
allele — SNP, deletion, insertion, MNP — from its REF and ALT strings
alone, without needing any extra type tags or metadata.

**Allele strings are literal sequences.** Uppercase nucleotide strings
over `{A, C, G, T, N}` — no symbolic alleles (`<DEL>`, `<INS>`,
`<NON_REF>`, `<*>`, `*`) are ever written. Indels are expressed
through the REF/ALT length difference, and the need for `<NON_REF>`
is eliminated by the scalar-based likelihood reconstruction (see §"Why
PLs are *not* stored").

**REF is always `alleles[0]`.** Downstream consumers can rely on this
unconditionally.

**`ref_span = len(alleles[0])`.** The number of reference positions
the record covers is defined by the REF sequence length; no separate
field is stored. Every alt in the record replaces the same reference
stretch.

**Variant shape is derivable from REF/ALT lengths alone:**

- **SNP:** `len(REF) == len(ALT) == 1`.
- **Deletion:** `len(REF) > len(ALT)`; REF = anchor base + deleted
  bases, ALT = anchor base alone.
- **Insertion:** `len(REF) < len(ALT)`; REF = anchor base alone, ALT
  = anchor base + inserted bases.
- **MNP / complex:** `len(REF) == len(ALT) > 1`, or length mismatch
  with bases differing at multiple positions — represented as a
  single allele pair, no decomposition into component events at
  Stage 1.

**Indel anchoring follows VCF convention** — indels are anchored at
the position one before the variable region.

- Deletion of bases at reference positions X..X+l−1:
  - record at position X−1
  - REF = anchor base + deleted bases (length `l+1`)
  - DEL allele = anchor base (length 1, same `ref_span` as REF)
  - interior positions X..X+l−1 only produce records if covered by
    non-deleting reads; the deletion's evidence lives only at the
    anchor record
- Insertion of bases after reference position X:
  - record at position X
  - REF = base at X (length 1)
  - INS allele = base at X + inserted bases (longer than REF, same
    `ref_span = 1`)

**Overlapping events extend the anchor REF.** When a SNP and an indel
co-occur in the window of a single anchor (e.g. a deletion spanning
positions 101–102 co-occurring with a SNP at position 101 in some
reads), Stage 1 extends the record's REF sequence to cover the widest
event and expresses all observed alts against that extended span —
freebayes' haplotype-allele convention. The invariant that Stage 2
records must maintain is: **every allele in a record spans the same
reference stretch.**

**Alt ordering within a record is Stage-1-internal** and carries no
cohort-level meaning. The merger re-indexes at join time.

Consequence: given a record, a consumer can determine every allele's
type from the REF/ALT sequences alone, and can compute `ref_span`
from `len(alleles[0])`.

### Per-position record layout and compression

The `.psf` is pileup-style: one record per covered reference position.
This design accepts a larger file per sample than GATK-style gVCF
blocks in exchange for strictly lossless preservation of per-position
information — including alt observations supported by even a single
read, exact depth, exact quality sums, and strand / placement detail.
Stage 1 applies no thresholding; every observed allele is kept
regardless of how many reads support it.

Three design decisions shape the format:

1. **All observed alleles are kept**, down to single-read support.
   Stage 1 has no "minimum reads per allele" or similar filter —
   anything the BAM shows, the `.psf` records. At the project's
   target coverage (2–10×) these single-read observations are often
   the only per-sample trace of real cohort-level signal, and the
   decision to accept or reject an allele is deferred entirely to
   the joint stage.
2. **Uncovered positions produce no record.** Each record carries a
   `delta_pos` field (distance from the previous record, variable-length
   encoded). Absence of a record at a position unambiguously means
   "this sample had zero coverage there." This makes exome and targeted
   data naturally compact: off-target regions generate no records.
3. **Per-position overhead is accepted**, even at very low depth
   where the record header is comparable in size to the data itself.
   A coarser format for ultra-low-depth regions is a possible future
   optimisation but not needed now.

Compression, stacked:

- **Implicit position** via `delta_pos` (typically 1 for consecutive
  covered positions).
- **Columnar storage** within a block of records — each scalar lives
  in its own column so that adjacent-position similarity compresses
  well.
- **zstd** applied per block.
- Optionally, **lossless run-length encoding** for runs of positions
  whose entire record (alleles + scalars) is bit-identical. These
  runs are common in uniformly-covered non-variant regions and
  compress extremely hard. The decompressor expands them back to
  one record per position — no information is lost.

Note: lossless RLE is *not* a gVCF block. A gVCF block merges adjacent
positions with *similar* data into a single summary; lossless RLE
merges adjacent positions with *identical* records and preserves the
full data on decompression.

### Why custom binary rather than BCF

The `.psf` is an internal pipeline artefact consumed only by
Stages 3–6 of this tool; no external consumer reads it. That
removes the main reason to adopt BCF (tool interop), and the data
shape — per-position pileup records with per-allele scalars and
phase chain ids — does not map naturally onto BCF's per-variant
FORMAT model. Scalars would have to be shoehorned into
non-standard FORMAT fields that other BCF consumers could not
interpret correctly anyway, so the nominal interop advantage
would be lost in practice.

Custom binary gives, in exchange for giving up bcftools interop:

- **native fit to the per-position record shape** — scalars and
  phase chain ids live as first-class fields rather than
  overloaded FORMAT columns;
- the full compression stack described above (`delta_pos`
  implicit position, columnar per-block layout, zstd, optional
  lossless RLE for identical runs) — realistic ~2–3× smaller
  than an equivalent BCF;
- **streaming write and streaming read** with no back-patched
  headers and no CSI index needed (Stages 3–6 scan entire files
  in order);
- **one file per sample as the natural granularity** — which is
  what cheap cohort recall requires, and which BCF tooling is
  not built around (BCF's natural granularity is multi-sample
  cohort files).

The costs are real but bounded: no bcftools interop, a small
`psf dump` / `psf head` utility for debugging (~a few hundred
lines), and explicit version management via a magic number +
version word at the file head. Readers reject unknown versions
rather than silently guessing.

A `psf → BCF`/VCF exporter is not planned now. It is a
straightforward data-shape transformation that can be added as
its own module if and when an external consumer requires it; no
design decision here commits us either way.

### Why PLs are *not* stored

A natural instinct — and how GATK's gVCF works — is to pre-compute and
store per-genotype PL values in the per-sample file. We deliberately do
not, because in this design **the PLs are fully derivable from the five
scalars**. freebayes' per-sample likelihood is

```
L(G) = prodQout(G) + multinomialSamplingProbLn(allele_probs(G), obs_counts)

     = [ Σ over alleles a ∉ G of S_a ]
     + multinomial( allele_probs(G), obs_counts )
```

where `S_a = Σ max(ln_BQ, ln_MQ)` and `obs_counts` are two of the stored
scalars, and `allele_probs(G)` is determined by `G` itself. Storing PLs
as well would be redundant: it pre-computes, for a specific subset of
the possible genotypes (those over the sample-local allele set), a
quantity that is trivially recomputable from data we already have.

Worse, a stored PL is pinned to the sample-local allele set. The
hardest case for the merger — a genotype involving an allele the
sample did not locally observe — is *not* covered by any stored PL; it
has to be computed from the scalars anyway. So pre-stored PLs only help
in the easy case, where recomputation from scalars is trivially cheap.

The reasons to store PLs that apply to GATK's gVCF do not apply here:

| Possible reason to store PLs | Applies to `.psf`? |
|---|---|
| Locks in calibration if the likelihood formula changes | no — we've committed to freebayes-style; a formula change would be a new file version |
| BCF / gVCF interoperability | no — `.psf` is a custom format; BCF export is a separate concern |
| Saves joint-stage computation | marginally — reconstruction is a sum plus one multinomial call |
| Supports per-sample analysis standalone | no — a consumer can compute PLs from scalars just as easily |
| Encodes a formula the scalars cannot | no — scalars are sufficient for freebayes' likelihood |

The scalars are the single durable source of truth for per-sample
evidence.

### Why a stored `<NON_REF>` / `<OTHER>` PL is also not needed

For the same reason: the likelihood of a genotype involving an allele
not in the sample's local set is computable from the stored scalars.
Every local read disagrees with such an allele by construction, so the
`prodQout` contribution is the sum of `S_a` across all the sample's
local alleles, scaled by the number of copies of the unknown allele in
the genotype. No separate stored `<NON_REF>` PL is required.

GATK's gVCF needs `<NON_REF>` because its stored object *is* the PL; it
has no alternative source of information. Our stored object is the
observation summary, which is strictly more information — it functions
as a generalised `<NON_REF>` that works for any new allele, including
compound haplotype alleles that `<NON_REF>` cannot represent.

## Stage 3 — low-complexity filter

### Purpose

Skip reference positions where short-read alignment and genotyping are
known to be unreliable — principally tandem repeats, homopolymers, and
other low-complexity sequence. Enabled by default. This is the single
most effective simple filter against pipeline false positives.

### Shape: a streaming filter, not a batch stage

Stage 3 is **not** a separate pass that produces a mask artefact on
disk. It is a streaming per-position filter that sits in the iterator
chain, between the multi-way merge of `.psf` files and the grouping
in Stage 4. When the DUST scorer flags a reference position as
low-complexity, the per-position records at that position are
silently dropped from the stream — they never reach Stage 4, so no
`OverlappingVarGroup` is ever constructed for them, and no merged
record is emitted.

Conceptually:

```
.psf_i ─┐
.psf_j ─┼─► multi-way per-position iterator ─► [Stage 3 DUST filter] ─► Stage 4 grouping ─► ...
.psf_k ─┘                                              ▲
                                                       │
                                        reference FASTA ─► DUST scorer
```

No mask file is written between stages. The DUST scorer is a small
component that sits alongside the iterator and is consulted in-line
as each position flows past. Implementation-wise this is just
another iterator adaptor in the Rust iterator chain, analogous to
`Iterator::filter`.

### Why this filter exists

Short-read aligners place reads ambiguously inside tandem repeats and
long homopolymers: the same read can be positioned at multiple
offsets with identical alignment scores, and the aligner picks one
essentially arbitrarily. Neither BAQ nor BQSR fixes this — they
adjust per-base confidence, not repeat-wide placement ambiguity. Both
freebayes and GATK routinely emit false-positive variant calls in
these regions, and the cheapest reliable mitigation across the whole
field is to drop the affected reference positions up front.

Filtering at the iterator level (rather than at Stage 1) keeps the
per-sample `.psf` files faithful to what the BAM actually showed: the
filter is a cohort-level policy, not a property of the raw
observations. The same `.psf` files can be re-run later with a
different complexity setting — or none at all — without regeneration.

The other reason for keeping the filter out of Stage 1 is that
pushing it in would not actually save meaningful disk space.
Low-complexity regions are typically short stretches (usually tens to
a few hundred bp each) and account for a small fraction of the total
reference. Dropping them at `.psf`-write time would shrink per-sample
files only marginally, and in exchange we would be committing every
`.psf` to a single DUST threshold chosen at the time of writing.
That trade — small space saving in exchange for locking in one
complexity parameter forever — is clearly not worth it, so the
filter lives at Stage 3 where the threshold can be revisited freely.

### Algorithm — BLAST DUST

The same algorithm used by NCBI BLAST's DUST filter (Tatusov and
Lipman). For each sliding window of `w` reference bases:

1. Count each triplet (`AAA`, `AAC`, ..., `TTT`) within the window.
2. Score the window as `s = Σ over triplets t of f_t · (f_t − 1) / 2`,
   where `f_t` is the count of triplet `t` in the window.
3. If `s` exceeds the threshold `T`, every position in the window is
   low-complexity.

Defaults (matching the standard DUST parameters): `w = 64`, `T = 20`.

Because DUST uses a sliding window, the filter only needs to maintain a
small moving buffer of the reference (~`w` bp each side of the current
position). It does not need a whole-genome mask in memory, and
computation is incremental as the iterator advances.

### How the filter behaves at the iterator

At each position in the multi-way per-position stream:

1. Advance the DUST scorer's window to cover the position in question.
2. If the DUST score at that window flags the position as
   low-complexity, the filter silently drops this position's records
   and advances to the next.
3. Otherwise, the records for that position are yielded downstream to
   Stage 4's grouping iterator unchanged.

The filter does not aggregate, reorder, or buffer records; it just
decides "pass or skip" per incoming position. Order is preserved.
Memory is constant (bounded by the DUST window size).

### Parameters and opt-out

Following Design principle 2 (no silent defaults, every decision
explicit), all filter behaviour is user-visible:

| Flag | Default | Effect |
|---|---|---|
| *(default — filter on)* | on | Apply DUST filter with `w = 64`, `T = 20`. |
| `--no-complexity-filter` | — | Disable the filter entirely. Every position from the upstream iterator is passed through. |
| `--complexity-window N` | 64 | Override DUST window size. |
| `--complexity-threshold T` | 20 | Override DUST score threshold. |

Opt-out is explicit (`--no-complexity-filter`), never silent — per
Design principle 3, if the pipeline is running without the filter the
user must have asked for it.

### What the filter does *not* do

- It does not modify `.psf` files. Those remain faithful records of
  what the BAM showed.
- It does not adjust per-read qualities. BAQ and the `max(ln_BQ,
  ln_MQ)` scalar remain responsible for per-base confidence.
- It does not use sample data. The filter depends only on the
  reference sequence — by design.
- It does not mask *soft*. There is no "low-confidence" tier between
  dropped and kept. A position is either skipped or processed.
- It does not produce any intermediate file between Stage 3 and
  Stage 4. The filter is purely a transformation on the iterator
  stream.

## Stage 4 — grouping

### Purpose

Walk the filtered per-position stream in genomic order and bundle
positions whose variants overlap into a single `OverlappingVarGroup`.
The output is a stream of independent groups; everything downstream
(Stage 5 onward) can treat each group in isolation.

### What "overlapping" means here

Positions are grouped together when their alleles share reference
coverage in at least one sample. The canonical cases:

- A deletion anchored at position P with `ref_span = k` forces
  positions P..P+k−1 into the same group with P.
- A SNP at position Q that falls inside another sample's deletion
  spanning positions P..P+k−1 is drawn into the group at P.
- Two indels anchored a few bp apart whose spans touch get merged
  into one group.
- Chains of the above propagate transitively: if positions P–Q and
  Q–R each overlap, P, Q and R are all in one group.

Positions with no overlap across samples form trivial groups of size
one. Most of the genome falls in this category.

### What it produces

An `OverlappingVarGroup` contains:

- the group's chromosome and `[start, end]` range
- every `.psf` per-position record from every sample that falls in
  that range
- the per-allele scalars and phase chain ids from those records,
  carried through verbatim

Everything Stage 5 needs for a group is inside the group itself. No
cross-group context is required, which is what makes Stage 5
trivially parallelisable.

### Why grouping is its own stage

Grouping is **sequential by nature**: you cannot decide a group's
right boundary without walking positions in order, because the next
position might extend an in-progress deletion into the current group.
Any attempt to parallelise grouping directly would require solving
that boundary-detection problem with speculation or chunking.

Splitting grouping out as its own stage isolates the sequential work
to the minimum that actually has to be sequential. The downstream
stream of independent groups is then free for parallel consumption.

### Order and streaming

Stage 4 emits groups in genomic order, one at a time, as soon as each
is closed. Memory is bounded by the widest open group (typically a
few hundred bp, rarely more). It does not buffer the whole cohort.

The implementation is conceptually similar to the existing
[`VariantGroupIterator`](../../src/variant_grouping.rs) —
multi-way-merging the per-sample streams and accumulating records into
a group until the group closes — now operating on `.psf` input and on
the filtered per-position stream produced by Stage 3.

## Stage 5 — per-group processing

### Purpose

For each `OverlappingVarGroup`, produce one merged multi-sample
variant record: unify alleles across samples, reconstruct per-sample
likelihoods against the merged allele set, and reconcile compound
haplotype alleles where needed.

### What it does, per group

1. **Allele unification.** Across all samples in the group, build the
   union of observed alleles. Handle compound haplotype alleles
   spanning multiple original per-sample variant records (e.g.
   deletion + nearby SNP co-occurring on the same haplotype in some
   samples).
2. **Per-sample likelihood reconstruction.** For each sample, for
   each genotype `G` in the merged allele set, compute the likelihood
   directly from the stored scalars using freebayes' standard
   formula:

   ```
   L(G) = [ Σ over alleles a ∉ G of S_a ]
        + multinomial( allele_probs(G), obs_counts )
   ```

   where `S_a` and `obs_counts` come from the stored scalars. Alleles
   in `G` that the sample never locally observed contribute `S_a = 0`
   and `obs_count = 0`, so they fall out of the multinomial; the
   sample's reads supporting its own local alleles contribute to
   `prodQout` for any genotype that doesn't include those alleles.
   The generalised `<NON_REF>`-equivalent behaviour is automatic.
3. **Compound haplotype consistency.** For a merged genotype `G`
   containing a compound haplotype allele, use the phase chain ids
   across the group's `.psf` records to decide whether the sample's
   own reads carry any support for the full compound haplotype. If
   the constituent alleles of the compound share a phase chain id in
   this sample, the compound has support in this sample; otherwise
   the reads supporting any constituent are treated as disagreeing
   with `G` (standard `prodQout` contribution).
4. **Emit a merged variant record** ready for Stage 6.

Sites where compound-allele reasoning cannot be evaluated cleanly
(e.g. a compound haplotype spanning more than one read's worth of
genome, with no phasing evidence in this sample) fall back to the
standard formula and are flagged in the merged record so Stage 6
knows the likelihood is approximate.

### The statistics, explained

This subsection expands the likelihood formula in step 2 above into
the full statistical story, written for geneticists not specialised
in statistics. Every formula is given an intuitive gloss. For
source-level detail on the same quantities, see
[freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md),
which walks through freebayes' original C++ implementation, and
[gatk_vs_freebayes_comparison.md](gatk_vs_freebayes_comparison.md)
for the same terms in GATK.

#### The question being asked

For each sample `s` at each site, we want the posterior:

```
P(G_s | reads from every sample in the cohort)
```

"given everything every sample's reads say, what is the probability
that sample `s`'s genotype is `G`?". Bayes' theorem decomposes this
into two pieces:

```
P(G_s | reads)  ∝   P(reads_s | G_s)    ×   P(G_s | cohort-level information)
                    └── likelihood ──┘      └─────────── prior ───────────┘
```

Stage 5 computes the **likelihood** — `P(reads_s | G)` — for every
merged candidate genotype `G` and every sample `s`. Stage 6 combines
that likelihood with the prior and iterates to a final posterior.

#### What the likelihood actually measures

`P(reads_s | G)` asks: "if sample `s` really had genotype `G`, how
likely are the reads we observed?". Two conditions have to hold for
the reads to be likely under `G`:

1. Any read supporting an allele **not in** `G` has to be explained
   as something other than sample `s`'s germline genotype —
   principally a sequencing error or a mismapping, though at low
   levels contaminating DNA from another source can also account for
   such reads (see the contamination caveat below). The BQ/MQ-covered
   part of that cost is what appears in the error-cost term.
2. Among reads supporting alleles **in** `G`, the split across alleles
   must look like what `G` predicts. A heterozygous `A/T` is expected
   to produce roughly 50/50 A:T reads; a 20/0 split is unlikely under
   `A/T`. A homozygous `A/A` is expected to produce ≈100% A reads; a
   10/10 split is unlikely under `A/A`.

These two checks are exactly the two terms in the formula:

```
L(G)  =   Σ over alleles a ∉ G of S_a              ←  "error-cost" term
        + multinomial(allele_probs(G), obs_counts_in_G)   ←  "allele-balance" term
```

Reads supporting alleles not in `G` flow *only* through the
error-cost term. They are taken out of the multinomial so that the
multinomial is not multiplied by a zero expected probability.

#### The error-cost term: `Σ_{a ∉ G} S_a`

`S_a` is the per-allele quality scalar stored in the `.psf` for
sample `s`:

```
S_a  =  Σ over reads supporting a of  max( ln(BQ_error), ln(MQ_error) )
```

where `BQ_error = 10^(-Q/10)` is the per-base sequencing-error
probability implied by the Phred base quality and `MQ_error =
10^(-MQ/10)` the per-read mismapping probability from MAPQ. Per read
the **larger** error dominates — the read is only as trustworthy as
its weakest link — so the maximum of the two log-error-probs is
taken before summing.

Under the hypothesis that sample `s` has genotype `G`, every read
supporting an allele `a ∉ G` has to be explained as a sequencing
error (BQ) or a mismapping (MQ). The total log-probability "all
these reads are BQ/MQ-explained" is exactly `Σ_{a ∉ G} S_a`. A
third mechanism — contamination — is not captured by
`max(ln_BQ, ln_MQ)`; see the subsection below.

**Geneticist intuition.** Twenty A-reads at Q=30 accumulate into
`S_A ≈ 20 × -6.9 = -138`. If we evaluate `G = T/T`, the full `-138`
is applied — effectively zero likelihood. This is how the error-cost
term crushes genotypes that are grossly contradicted by the reads.

#### When a read "not in G" isn't actually an error: contamination

The error-cost term above models out-of-`G` reads as BQ/MQ errors
only. In reality such reads can also arise from:

- **Cross-sample contamination.** A fraction of reads may derive
  from another individual's DNA — library-prep carry-over, index
  hopping on the sequencer, or an outright sample swap. These reads
  are real DNA, not sequencing artefacts, so their BQ and MQ are
  fine; they simply do not belong to sample `s`.
- **Somatic mosaicism.** Reads carrying a truly present but
  non-germline variant in a small subset of cells.

The `max(ln_BQ, ln_MQ)` scalar models neither case explicitly, so
the base error-cost formula overpenalises genotypes in samples
with real contamination. In well-prepared libraries contamination
is typically <1 % and the mis-modelling is noise-level, but at
higher levels — FFPE tumour material, environmental samples,
deliberately pooled libraries — it biases the pipeline toward
calling contaminating alt reads as real low-frequency variants.
Because low-level contamination is common enough in practice to
matter, **the pipeline models it explicitly**, as follows.

**The mixture model.** Each read in sample `s` is treated as
coming from a two-component mixture:

- with probability `1 − c_s`, the read is from `s`'s own DNA (the
  standard per-read model — consistent with `G_s` at expected
  allele fractions, modulated by BQ/MQ errors);
- with probability `c_s`, the read is from a random cohort member,
  whose allele at this site is distributed as the population
  allele frequency `p` (one read = one chromosome drawn under HWE).

Per-read probability of supporting allele `a`:

```
P(read → a | G_s, c_s, p)
   =  (1 − c_s) · P_own(a | G_s, BQ, MQ)     ← existing own-DNA term
   +  c_s       · p_a                         ← contamination term
```

The genetics insight behind this form is that contaminating reads
do not support random alleles — they support alleles that exist
in *other samples*. A hom-A sample whose "odd reads" at a site
with population AFs `(p_A = 0.3, p_B = 0.7)` preferentially
support `B` rather than `C` or `G` is showing the fingerprint of
contamination, not random sequencing error (which would be
roughly unbiased across the three non-`A` bases). The mixture
term `c_s · p_a` captures this formally.

**Estimating `c_s`.** `c_s` is one more parameter estimated by
the Stage 6 EM loop, jointly with the population allele
frequencies `p`:

- **E-step** — per-sample genotype posteriors given the current
  `(p, c_s)`.
- **M-step on `p`** — as before (posterior-weighted allele
  counts + Dirichlet pseudocounts).
- **M-step on `c_s`** — a 1D maximisation per sample: choose the
  `c_s` that best explains the sample's odd-read pattern given
  `p`.

Most of the information about `c_s` comes from sites where the
sample is confidently non-heterozygous and yet carries odd reads
matching population-common alt alleles. No external estimate or
reference panel is required. A user-supplied seed value (as in
freebayes' `--contamination-estimates` or the output of
VerifyBamID2 / GATK `CalculateContamination`) is still accepted
and overrides the EM estimate when present.

**Scalar-storage approximation.** The mixture log does not
factor cleanly over the stored scalars:

```
log L_a = Σ over reads of log[ (1 − c_s) · P_err(read) + c_s · p_a ]
```

The `log` wraps a sum *inside* the per-read product, so each
read needs its individual quality — a per-read value that `S_a`
threw away when Stage 1 pre-summed. The pipeline therefore uses
the **homogeneous-quality approximation**: all reads supporting
allele `a` are treated as sharing the mean quality
`S_a / count_a`, so the mixture is evaluated once per
(sample, allele) and multiplied by `count_a`. No extra scalar is
added to the `.psf`.

The approximation is exact when all reads supporting `a` happen
to have identical quality, and degrades as per-read qualities
spread. At the 1–10 % contamination regime this matters for, the
contamination term in the mixture tends to dominate each
per-read log and the approximation stays tight. A follow-up
option — adding one more per-allele scalar that makes the
mixture reconstruction exact — is deferred until we see
approximation error in practice.

Relationship to freebayes: freebayes' own optional correction
([DataLikelihood.cpp:116-130](../../freebayes/src/DataLikelihood.cpp#L116-L130),
[freebayes_posterior_gt_probs.md §2d](freebayes_posterior_gt_probs.md#2d-contamination-and-reference-bias))
takes a different shape — rescaling the reference-allele
sampling probability rather than the per-read mixture — because
freebayes operates on reads in memory and does not face the
scalar-reconstruction constraint. It also requires the user to
supply a contamination fraction externally; this pipeline infers
`c_s` from the cohort as a side effect of the EM.

#### The allele-balance term

`obs_counts_in_G` are the observation counts for the alleles in `G`
and `allele_probs(G)` is the vector of expected allele fractions
under `G` — ploidy-aware:

- `A/A` diploid: `(A, T) → (1, 0)`
- `A/T` diploid: `(A, T) → (0.5, 0.5)`
- `A/A/T` triploid: `(A, T) → (2/3, 1/3)`

The multinomial PMF of the observed counts given those expected
fractions is one standard probability calculation. Its role is to
**penalise genotypes whose predicted balance does not match the
observed counts**.

**Geneticist intuition.** The discriminating power of this term
grows with depth. At 20× coverage:

- 20 A + 0 T under `A/A`: multinomial = 1 (log 0). Under `A/T`:
  multinomial = `0.5²⁰ ≈ 10⁻⁶` — roughly a Phred-60 penalty. So the
  multinomial alone can strongly rule out `A/T` from a clean
  homozygous pileup.
- 10 A + 10 T under `A/T`: multinomial peaks at the 50/50
  expectation. Under `A/A` the T-reads leave the multinomial and
  contribute to the error-cost term instead (~`10 × -6.9 = -69`),
  which is how `A/A` is penalised in this case.

At 5× the same comparison is much gentler (`0.5⁵ ≈ 0.03`, ~Phred
15). Low-depth samples therefore rely more heavily on the
error-cost term and on the cohort-level prior to tell genotypes
apart — one of the main reasons joint calling across many samples
adds value over per-sample calling.

#### Why the five `.psf` scalars are enough

| Scalar | Role |
|---|---|
| observation count | `obs_counts` in the multinomial (likelihood) |
| `S_a = Σ max(ln_BQ, ln_MQ)` | `S_a` in the error-cost term (likelihood) |
| forward-strand count | strand-bias prior (Stage 6) |
| placed-left count | read-placement bias prior (Stage 6) |
| placed-start count | fragment-position bias prior (Stage 6) |

Both likelihood terms are written purely in terms of `obs_counts`
and `S_a`, so the Stage 5 reconstruction is **exact** — bit-for-bit
identical to what freebayes would compute with the full BAM in
memory. The three remaining scalars feed observation-bias priors in
Stage 6.

A natural question is whether `sum(BQ)` and `sum(MQ)` separately
could replace `S_a`. They cannot: the `max` in log space is
nonlinear and has to be applied **per read** before summing. See
[freebayes_posterior_gt_probs.md §Compactness implications](freebayes_posterior_gt_probs.md#compactness-implications--important-subtlety-about-bq-and-mq).

#### Alleles the sample never locally observed

A sample's `.psf` lists only alleles its own reads support. At the
cohort level, the merged genotype `G` for sample `s` may include an
allele `x` that `s` never observed — e.g. another sample
contributed a novel ALT. The formula handles this automatically:

- `obs_count(x) = 0` and `S_x = 0` for sample `s`, so `x`
  contributes nothing in either term.
- The sample's **own** locally observed alleles `a` that are not in
  `G` still contribute their full `S_a` to the error-cost term —
  those reads have to be errors under `G`.

This reproduces freebayes' `<NON_REF>`-like behaviour for free. No
separate `<NON_REF>` PL is stored, and the same mechanism
generalises to compound haplotype alleles that `<NON_REF>` alone
cannot represent (see §"Why a stored `<NON_REF>` / `<OTHER>` PL is
also not needed" above).

#### Compound haplotype alleles: the phase-chain check

For most genotypes `G` the formula above is exact. The one case
that needs extra machinery is a **compound haplotype allele** — an
allele that is itself a combination of per-position variants
co-occurring on the same haplotype (e.g. "SNP at 101 + deletion
spanning 102–105" on one chromosome).

A per-position `.psf` record stores evidence for each variant
separately, but a compound is only "present" in sample `s` if `s`'s
reads carry the pieces *together*, on the same read or read pair.
Stage 1 therefore records a **phase chain identifier** on each
allele whose support is linked to another position's allele by a
shared read. At Stage 5:

- If the compound's constituents share a phase chain id in sample
  `s`, the sample's reads do carry the compound → evaluate the
  likelihood normally, using the compound's joint observation count
  and joint quality scalar.
- Otherwise, the sample's reads carry only the constituents
  separately. Under the hypothesis that `s` has the compound, those
  reads have to be explained away. Stage 5 falls back to the
  standard formula and flags the record as "compound-approximate"
  so Stage 6 can interpret the resulting GQ conservatively.

For every non-compound site — and for compound sites where the
evidence is clean in every sample — the reconstruction is exact.
For the fraction of compound sites where some samples have only
partial phasing, the likelihood is approximate and explicitly
labelled.

**VCF encoding of the flag.** The approximation is surfaced in
the final multi-sample VCF as a single boolean per-sample FORMAT
field, `CA` (compound-approximate): `1` on records where this
sample's likelihood was reconstructed under the phase-broken
fallback, `0` otherwise. The field is absent entirely on
non-compound records, to keep simple-SNP records uncluttered.
No tiered scheme and no per-genotype granularity: the flag
addresses a small minority of sites (compound haplotype +
partial phase in some sample) and does not need a richer
taxonomy to be useful. Downstream consumers that care about
likelihood provenance filter on `CA`; everyone else ignores it.

#### From likelihood to posterior (Stage 6, in one picture)

Stage 5 outputs a per-sample, per-genotype likelihood table. Stage
6 combines it with two pieces of population-genetic prior
information:

1. **Allele frequency in the cohort — Dirichlet prior with tiny
   alt pseudocounts.** Rare alleles are *a priori* less common
   than common ones; a single-sample homozygous-alt call out of
   a thousand samples is implausible unless the likelihood is
   very strong. Stage 6 estimates `p` at each site by EM,
   iteratively tightening a point estimate from soft per-sample
   posteriors. Full derivation in
   [posterior_gt_probs.md](posterior_gt_probs.md).

   The prior on `p` is a **Dirichlet distribution** with
   GATK-style pseudocounts: `α_ref = 10`, `α_alt = 0.01` for
   SNPs, `α_alt = 0.00125` for indels. User-overridable via
   `--ref-pseudocount`, `--snp-alt-pseudocount`, and
   `--indel-alt-pseudocount`. The M-step posterior mean is
   closed-form:

   ```
   p̂_k  ∝  α_k  +  E[n_k]
   ```

   where `E[n_k]` is the posterior-weighted expected count of
   allele `k` across the cohort. The pseudocounts act as a
   pre-observation belief (≈99.9 % ref, ≈0.1 % any given alt):
   at small cohorts they shrink weak alt evidence toward zero;
   at large cohorts real alt signal overwhelms them and `p̂`
   tracks the data. This is the same prior shape used by GATK's
   `GenotypeGVCFs`.

   **Why Dirichlet rather than Ewens' sampling formula.**
   freebayes' equivalent prior is Ewens' sampling formula
   ([freebayes_posterior_gt_probs.md §Component 2](freebayes_posterior_gt_probs.md#component-2--pallele-frequency-ewens-sampling-formula)),
   parameterised by a population mutation rate `θ`. Both priors
   encode the same genetics — rare alleles are *a priori* more
   plausible than common ones, 50/50 splits are unusual — but
   they differ in mathematical shape and, crucially, in whether
   they fit EM. We pick Dirichlet for three reasons, recorded
   here so a future redesign knows what has to be re-argued
   before switching:

   - **Conjugacy with the M-step.** Dirichlet-multinomial is
     the conjugate pair for allele-frequency updates. The
     M-step reduces to a closed-form sum — one addition per
     allele per iteration. Ewens has no conjugacy with anything
     in this pipeline; mixing it into the EM loop would require
     a numerical maximisation every round, giving up the main
     reason EM is the right algorithm at this scale. This is
     the same reason GATK uses Dirichlet, and we follow.
   - **Fits the pipeline's data structures.** Ewens scores the
     *joint partition* of allele counts across samples (e.g. a
     "one major + one rare" partition is scored differently
     from a "two equal" partition). This partition is the
     natural object in freebayes' combo search, where each
     combo has one concrete partition. Our pipeline never
     builds partitions of this form — it tracks per-sample
     genotype posteriors and a single `p` vector. The object
     Ewens scores is simply not represented in our data flow.
   - **At our target scale, the choice is invisible.** At
     thousands of samples the data dominates either prior and
     `p̂` converges to the same value regardless. The
     theoretical gap between Ewens and Dirichlet matters only
     at small cohorts with marginal variants — which is exactly
     the regime where the EM point estimate of `p` is itself
     suspect (see §"Why EM rather than freebayes' combo search"
     below). So there is no scenario in which Dirichlet loses
     meaningfully to Ewens at target scale; at small scale the
     entire EM-vs-combo-search decision reopens, and the
     prior-shape choice rides along with it rather than
     standing alone.

   **Consequence for a future redesign.** Swapping Dirichlet
   for Ewens in isolation is not a well-posed change: Ewens
   only pays off in combination with combo-search
   marginalisation (freebayes' full approach), and combo search
   does not scale to thousands of samples. A future pipeline
   that wants Ewens-level calibration at small cohorts should
   revisit the EM-vs-combo-search decision first; the prior
   choice follows from that, not the other way around.
2. **Hardy–Weinberg proportions given `p`, adjusted by an
   inbreeding coefficient `F`.** Under pure HWE (`F = 0`) the
   expected genotype frequencies in the cohort are `(1 − p)²`,
   `2 p (1 − p)`, `p²` (diploid biallelic). With `F > 0`:

   ```
   P(AA) = (1 − p)²  +  F · p (1 − p)
   P(AB) = 2 p (1 − p) · (1 − F)
   P(BB) = p²        +  F · p (1 − p)
   ```

   — homozygotes are enriched and heterozygotes depleted by a
   factor of `F`. `F = 0` recovers pure HWE; `F = 1` means fully
   homozygous lines (e.g. doubled haploids or long-inbred
   material); negative `F` is allowed and encodes excess
   heterozygosity (e.g. heterotic F1 populations).

   `F` is exposed as the `--inbreeding` CLI parameter — a single
   cohort-wide scalar, **default `0`**. The project targets
   plant populations, where self-pollination, clonal
   propagation, and inbred-line breeding make `F > 0` the norm
   rather than the exception, and where assuming strict HWE
   would systematically overcall heterozygotes. Typical values:
   ≈0.98 for pure inbred lines or doubled haploids, ≈0.5 for a
   few generations of selfing from an outcrossing ancestor, 0
   for open-pollinated species or wild panels.

   Internally, `F` is stored as a **per-sample vector `F_s`**,
   initialised from the single user-supplied scalar (all samples
   get the same value). The HWE prior then reads `F_s` rather
   than a scalar `F` throughout. Keeping the internal
   representation per-sample from day one costs essentially
   nothing (one extra vector allocation) but leaves the upgrade
   path open: adding per-sample `F` later (via a user-supplied
   per-sample file) becomes a user-interface change, not an
   algorithmic change. The user-facing interface stays simple —
   one scalar — until there is a real case that needs more.
   Beyond a user-supplied per-sample file, data-driven inference
   of `F_s` from excess homozygosity in the sample's own calls
   is a further possible extension, deferred for the same
   reason: no concrete need today.

Putting it together:

```
P(G_s | reads_cohort)  ∝   L_s(G_s)  ×  HWE(G_s | p̂)  ×  bias_priors
```

normalised so each sample's per-genotype posterior sums to 1. From
there:

- **GT** = argmax genotype.
- **GQ** = Phred of `1 − P(best_genotype_s)`.
- **QUAL** (site-level) = Phred of `Π_s P(hom-ref)_s`, i.e. the
  posterior probability that the site is actually variant in at
  least one sample.

The observation-bias priors use the three bias scalars from the
`.psf`. Conceptually: a candidate whose alt-supporting reads are
all on one strand, or bunched at one end of the fragment, is more
likely to be a systematic artefact than a real variant, so its
posterior is down-weighted. Structural form as in
[freebayes_posterior_gt_probs.md §Component 4](freebayes_posterior_gt_probs.md#component-4--observation-level-bias-priors-optional);
unlike GATK, which handles these as post-hoc annotations (SOR, FS),
this pipeline keeps them inside the posterior as freebayes does.

#### Why EM rather than freebayes' combo search

freebayes does not estimate a single `p̂`; it enumerates (searches)
joint genotype *combos* across all samples and marginalises back to
per-sample posteriors. The full three-way comparison (freebayes
combo search, GATK EM, full MCMC) is in
[freebayes_posterior_gt_probs.md §Three ways to solve the same problem](freebayes_posterior_gt_probs.md#three-ways-to-solve-the-same-problem-freebayes-gatk-mcmc)
and [gatk_vs_freebayes_comparison.md §Breaking the chicken-and-egg](gatk_vs_freebayes_comparison.md#breaking-the-chicken-and-egg).
Short version:

- **Small cohorts (tens to a few hundred samples):** combo search
  is better calibrated because it implicitly integrates over
  uncertainty in `p`. EM collapses `p` to a point estimate and can
  be slightly overconfident at rare variants.
- **Target scale of this pipeline (thousands of samples):** the
  posterior on `p` is so tight that `p̂` effectively equals the
  integral. EM and combo search converge to the same numbers, and
  EM wins on scaling (linear in sample count; combo search grows
  much faster and runs out of budget past a few hundred samples).

So EM here is a cheap approximation at small N and *the* correct
answer at large N. See §"Stage 6 — posterior engine §Why EM, not
combo-search marginalisation" below for the corresponding
discussion at the Stage-6 level.

### Parallelism

Groups are independent of each other: each one contains everything it
needs, and nothing in the processing of group `g` depends on the
processing of group `g − 1` or `g + 1`. Stage 5 processes the
`OverlappingVarGroup` stream with rayon's parallel iterators, across
worker threads, with an output reorder buffer to preserve genomic
order for Stage 6. Memory is bounded by the in-flight group count
(worker count × a few groups each), not by cohort size or genome size.

This is the headline parallelism of the cohort side of the pipeline.
Stage 1 parallelises across samples; Stage 5 parallelises across
groups. Stage 4 (grouping) and Stage 6 (posterior EM) remain
sequential but are cheap relative to Stage 5's per-group work.

### No in-tool chunking

The pipeline does not tile the genome into chunks for internal
parallelism. The iterator chain already bounds memory at every
stage: per-sample `.psf` readers hold one decompressed block each;
the multi-way merger holds one next-record slot per sample
(O(N_samples) × ~50 B); Stage 4 holds the currently-open group
(a few hundred bp × cohort-wide records); Stage 5 holds only the
groups currently being processed by rayon workers. None of those
grow with `N_samples × genome_length`, so a single machine handles
thousands of samples end-to-end without chunking.

Two scenarios that *would* motivate chunking are explicitly not
current requirements:

- **Distribution across a cluster.** Achievable externally by
  invoking the tool per chromosomal region — every node needs
  only the same `.psf` files and reference, which are already
  distributable.
- **Region-selective re-processing** after adding samples. Same
  external-invocation story.

Both scenarios additionally require random-access seek into each
`.psf` (a sidecar index with file offsets per chromosomal chunk),
which contradicts the streaming-only posture of the `.psf` format
(see Stage 2 §"Why custom binary rather than BCF"). If a real
in-tool need later emerges — e.g. a cluster-native deployment
that wants a single invocation across many nodes — the index
plus tile-boundary handling (for compound variants that straddle
tile edges) can be added against that concrete requirement
rather than speculatively now.

### Relationship to the existing code

This is conceptually the same algorithm as
[`create_variant_for_region`](../../src/genotype_merging.rs#L20),
generalised to:

- operate on `.psf` per-position records rather than parsed gVCF rows
- use the five per-allele scalars to reconstruct likelihoods, rather
  than relying on carried-over PLs
- reason about compound haplotypes via phase chain ids
- run under rayon across groups rather than the current
  double-buffered batch-of-1000 arrangement

## Stage 6 — posterior engine

### Algorithm

EM over allele frequency, extended with a per-sample contamination
parameter. The baseline algorithm is described in
[posterior_gt_probs.md](posterior_gt_probs.md); the contamination
extension is in Stage 5 §"When a read 'not in G' isn't actually an
error: contamination".

1. Initialise flat allele frequencies `p` and `c_s = 0` for every
   sample.
2. **E-step**: for each sample, compute per-genotype posteriors
   using the mixture likelihood (current `p` and `c_s`) + the
   HWE-with-`F` prior (see Stage 5 §"From likelihood to
   posterior"). `F` is cohort-wide, user-supplied via
   `--inbreeding`, default `0`.
3. **M-step on `p`**: update allele frequencies from
   posterior-weighted allele counts + Dirichlet pseudocounts
   (`α_ref = 10`, `α_alt = 0.01` SNP, `α_alt = 0.00125` indel;
   overridable via `--ref-pseudocount` /
   `--snp-alt-pseudocount` / `--indel-alt-pseudocount`). Rationale
   for Dirichlet over freebayes' Ewens prior is in Stage 5
   §"From likelihood to posterior".
4. **M-step on `c_s`**: 1D per-sample update that maximises the
   sample's likelihood at its soft-assigned genotypes, using the
   pattern of odd reads at confidently-non-heterozygous sites as
   the main signal. Skipped when the user has supplied a fixed
   `c_s` externally.
5. Iterate until convergence. Typically 3–5 rounds.
6. QUAL = Phred of `Π_s P(hom-ref)_s`.
7. Assign GT = argmax genotype. GQ = Phred of `1 − P(best_genotype)`.

### Why EM, not combo-search marginalisation

At thousands of samples the posterior on allele frequency `p` is tight
enough that EM's point-estimate approximation and freebayes'
full-marginalisation approach agree to many decimal places. EM is linear
in sample count; combo-search is much worse. Standing on EM is the
scalability win. See [freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md)
for the three approaches compared.

A future upgrade path, if better calibration at rare variants is ever
needed, is to swap the point estimate of `p` for a Dirichlet posterior on
`p` (variational Bayes). This stays linear in sample count and requires
only local changes to the posterior engine.

### Future extension: cohort haplotype-frequency priors for compound alleles

The current Stage 6 estimates per-position allele frequencies `p` and
uses them as priors on per-sample genotypes. A natural and load-bearing
extension is to **also estimate per-compound haplotype frequencies** and
feed them as priors on compound genotypes for samples where the
phase-chain check returns no information — the phase-broken fallback
path that Stage 5 currently flags as `CA = 1`.

The motivation is that thousands of samples carry strong cohort-level
signal about which compound haplotypes actually exist in the
population. Homozygous-for-the-compound individuals pin the haplotype
down unambiguously; heterozygotes' calls then have a meaningful prior
to lean on when read-level evidence (the phase chain) is insufficient.
This is structurally what classical population-phasing tools
(fastPHASE, BEAGLE, IMPUTE, SHAPEIT) do for whole-genome haplotype
reconstruction; here we'd apply the same idea narrowly inside an
`OverlappingVariantGroup`, leveraging the cohort to strengthen
specifically the cases where the chain is silent.

**Where it goes in the inference.** Strictly at the *prior* layer, on
top of the chain-based likelihood — it does **not** substitute for the
chain in `P(reads | genotype)` (see
[phase_chain.md §"Why not lean on cohort haplotype frequencies instead?"](phase_chain.md)
for why a substitution would be statistically unsound). Concretely:

- **E-step extension.** For each compound allele `C` in a group, the
  current E-step computes a per-sample posterior over per-position
  genotypes; extend it to compute a per-sample posterior over the
  *compound* genotypes too, using whatever per-sample likelihood is
  available (chain-based when the chain has evidence, phase-broken
  fallback otherwise).
- **M-step on `f_C`.** Add a new M-step that updates the cohort
  frequency `f_C` of each compound haplotype as the posterior-weighted
  count, with Dirichlet pseudocounts (small alt pseudocount, e.g.
  `α_compound = 0.001`, tighter than per-position alt because compound
  haplotypes are conditionally rarer).
- **Feedback into samples with empty chain intersection.** For
  samples that fell into the phase-broken fallback at Stage 5,
  replace the independent-constituents prior with the cohort-derived
  `f_C`. Samples with chain evidence keep using the chain-based
  likelihood unchanged — `f_C` is just an updated prior for them too.
- **HWE with `F`.** The same Wright's-`F` extension already used at
  the per-position level applies directly to compound haplotype
  pairs, so inbreeding / selfing populations remain handled.

**What the extension buys.**

- **Long-range compounds.** Compounds whose constituents lie beyond
  read/pair span (chains can't reach) get cohort-anchored calls
  whenever the cohort contains enough homozygotes.
- **Mate-missing and low-coverage cases.** Per-sample fallbacks
  currently flagged `CA = 1` get a much stronger prior, so their
  GQ becomes informative rather than conservatively-down-weighted.
- **Common compounds at scale.** With N in the thousands and
  compound frequency above ~5%, the cohort prior is tight enough
  to dominate noise in any single sample's fallback likelihood,
  recovering most of the calibration the chain would have given.

**What it does not buy.**

- **Rare compounds (`q < ~0.01`).** Hom-alt frequency under HWE is
  `q²`; below this threshold even N = 1000 produces ~0 homozygotes,
  and the extension has no anchor. Chain evidence remains the only
  way to call these — exactly the regime where the chain mechanism
  is load-bearing.
- **Replacement for the chain.** This is a prior, not a likelihood.
  Samples with chain evidence still need the chain-based per-read
  likelihood; cohort frequencies cannot be substituted there
  without breaking Bayesian soundness.

**Implementation cost.** Local to the posterior engine: one extra
inner loop in the E-step (per-compound posterior), one extra
M-step (per-compound frequency), and a small change to the
fallback path so it consults `f_C` instead of constituent-product
priors. Memory grows by O(#compounds_per_group), bounded by the
group's allele structure. The chain machinery is unaffected.

This extension is **not implemented in the v1 pipeline**. It is
recorded here so a future implementer has a complete picture of
how the cohort signal can strengthen the existing chain-based
inference without compromising the per-sample likelihood path.

## Properties of this design vs. the alternatives

| Property | freebayes joint | GATK GenotypeGVCFs | This design |
|---|---|---|---|
| Scales to thousands of samples | ✗ | ✓ | ✓ |
| BAMs needed at joint-calling time | ✓ | ✗ | ✗ |
| Per-sample artefacts reusable across cohorts | ✗ | ✓ | ✓ |
| Correct joining of simple alleles | ✓ | ✓ | ✓ |
| Correct **merging** of compound haplotype alleles | ✓ | ✗ | approximate, principled |
| Read-level observation-bias priors | ✓ | partial | partial (via observation summary) |
| Full Bayesian posterior on allele frequency | ✓ | no (point estimate) | no (EM point estimate) |

The approximation row deserves a note: for haplotype-compound merges,
neither the observation summary nor the reconstructed likelihood can
perfectly match what freebayes computes with the reads in memory. In the
small fraction of sites where this matters, Stage 5 flags the affected
records so downstream consumers know the likelihood was approximated.
For the overwhelming majority of sites — simple SNPs, simple indels,
and merges where the samples carrying the compound alleles are the ones
with direct per-sample evidence — the tiered rule gives exact or
near-exact likelihoods.

## Relationship to existing specs in this directory

- [general_project_specification.md](general_project_specification.md)
  describes the current gVCF-in / VCF-out tool. This architecture is the
  intended replacement; the current spec will need to be rewritten or
  retired once this design is firm.
- [per_sample_pileup_format.md](per_sample_pileup_format.md) is the PSP v1
  binary format draft. Its per-read record layout is superseded by the
  per-allele scalar summary decided here (see Stage 1 §"Related
  existing draft" and Stage 2 §"The five per-allele scalars"). PSP
  v1 needs to be revised: the per-read 3-byte record is replaced by
  five per-allele scalars (observation count, `Σ max(ln_BQ, ln_MQ)`,
  forward-strand count, placed-left count, placed-start count).
- [genotype_joining_specification.md](genotype_joining_specification.md)
  documents the current allele-joining algorithm. The merge step in this
  architecture generalises that algorithm to operate on `.psf` instead of
  gVCF; most of the algorithmic content carries over.
- [posterior_gt_probs.md](posterior_gt_probs.md) and
  [freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md)
  describe the two statistical approaches to posterior calculation. This
  architecture adopts the GATK-style EM approach (described in the first)
  for scalability reasons argued in the second.
- [snp_calling_from_bam.md](../feature_implementation_plans/snp_calling_from_bam.md)
  is an earlier sketch of this same pipeline with different emphasis.
  Open questions raised there (quality assignment for indels, use of
  rust-htslib, per-read vs summary storage) are unchanged and still open.

