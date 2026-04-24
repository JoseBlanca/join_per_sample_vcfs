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
depth. PSP v1 needs to be revised accordingly; the BAQ-adjusted BQ
values should also be reflected in whatever quality summary the
revised format records.

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

EM over allele frequency, exactly as described in
[posterior_gt_probs.md](posterior_gt_probs.md):

1. Initialise flat allele frequencies.
2. **E-step**: for each sample, compute per-genotype posteriors using the
   likelihood assigned above + HWE prior (with optional inbreeding
   coefficient `F`).
3. **M-step**: update allele frequencies from posterior-weighted allele
   counts + Dirichlet pseudocounts.
4. Iterate until convergence. Typically 3–5 rounds.
5. QUAL = Phred of `Π_s P(hom-ref)_s`.
6. Assign GT = argmax genotype. GQ = Phred of `1 − P(best_genotype)`.

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

## Open decisions

1. **Per-read records vs. per-allele summaries in `.psf`.** — **Resolved: per-allele summaries.**
   The five scalars listed in Stage 2 are sufficient to reproduce
   freebayes' likelihood and priors exactly. Per-read storage (the PSP v1
   approach) buys flexibility for future alternative likelihood models we
   have no plans for, at a linear cost in coverage. PSP v1 needs to be
   revised to match.

2. **Haplotype context representation.** — **Resolved: not stored;
   phase chain identifiers are sufficient.** The per-position allele
   records + phase chain ids linking alleles across positions provide
   everything the merger needs to reason about compound haplotypes.
   A separate "local consensus sequence" field on each allele was
   considered and rejected as redundant: the reference comes from the
   FASTA, the sample's alt alleles are already in the per-position
   records, and the "which reads carry which combination" linkage
   comes from matching phase chain ids across positions.

3. **Serialisation: BCF vs. custom binary.**
   - BCF gives tool interop, existing indexing, and a well-understood
     FORMAT extension mechanism.
   - Custom binary (PSP v1 route) is ~2–3× smaller and avoids BCF overhead
     for fields we use in a non-standard way.
   - The two can coexist: start custom, add a BCF exporter later.

4. **Build Stage 1 or leverage existing callers.** — **Resolved: build
   Stage 1 ourselves**, with BAQ computed inline in the per-sample
   caller. Rationale: the per-read likelihood model we need is
   freebayes-shaped and simple (Option B in the Stage 1 discussion),
   but `samtools calmd` is single-threaded and has been a real-world
   performance bottleneck. Implementing BAQ directly lets us
   parallelise it with rayon and pipeline it with BAM decoding. See
   Stage 1 §"Per-read likelihood quality" and §"Why BAQ in-process"
   for the full argument, including why Option C (PairHMM + local
   reassembly) is rejected at target coverage and left as a possible
   future project.

5. **Reference-confidence blocks: format and density.** — **Resolved:
   no blocks; per-position records instead.** The `.psf` stores one
   record per covered reference position with the full five-scalar
   summary for every observed allele, down to single-read support.
   Uncovered positions produce no record; gaps are implicit via
   `delta_pos`. Compression is achieved by columnar storage + zstd,
   plus optional lossless RLE for runs of identical adjacent records
   — but not by the lossy summarisation that gVCF blocks rely on,
   because at the target-coverage regime (2–10×) single-read
   observations can be the only per-sample trace of real cohort-level
   signal and must not be discarded at Stage 1. See Stage 2
   §"Per-position record layout and compression" for the full
   rationale.

6. **Chunking strategy for Stages 4 + 5.**
   - Cohort-wide merging at thousands of samples will be memory-heavy if
     done in one pass. Chunking by genomic region (e.g. 1 Mb tiles) and
     merging tiles independently, writing in parallel, is the natural
     path. Needs a careful story for tile boundaries that fall inside a
     compound variant.

7. **Approximation flagging.**
   - Which merged-record FORMAT field records "this sample's likelihood
     for this genotype was approximated via tier 3 / tier 4 fallback"?
   - Downstream consumers (GQ interpretation, filtering) should be able to
     distinguish exact from approximate per-sample likelihoods.

## Relationship to existing specs in this directory

- [general_project_specification.md](general_project_specification.md)
  describes the current gVCF-in / VCF-out tool. This architecture is the
  intended replacement; the current spec will need to be rewritten or
  retired once this design is firm.
- [per_sample_pileup_format.md](per_sample_pileup_format.md) is the PSP v1
  binary format draft. Its per-read record layout is superseded by the
  per-allele scalar summary decided here (Open decision 1). PSP v1 needs
  to be revised: the per-read 3-byte record is replaced by five per-allele
  scalars (observation count, `Σ max(ln_BQ, ln_MQ)`, forward-strand count,
  placed-left count, placed-start count).
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

