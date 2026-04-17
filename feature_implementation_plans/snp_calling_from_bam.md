# SNP Calling from BAM — First Draft

**Status:** Early draft. Open questions are flagged with **[QUESTION]** and **[DECISION]** tags; they are for us to discuss together in the next iteration.

## Goal

Replace the current gVCF-in / VCF-out workflow with an end-to-end variant caller that takes per-sample BAM files and produces a multi-sample VCF. The new pipeline computes genotype likelihoods from read evidence directly (instead of trusting gVCF PL/GT), merges alleles across samples (full merge, not GATK-style per-position joining), runs EM across the cohort, and writes final calls with posterior-informed GT/GQ/PL.

Inspiration (to be studied, not copied): GATK HaplotypeCaller / GenotypeGVCFs, bcftools/samtools mpileup+call, freebayes.

## High-level pipeline

```
[BAM sample 1] ─┐
[BAM sample 2] ─┼─ (stage 1: pileup per sample) ──► [binary.psp per sample]
[BAM sample N] ─┘                                              │
                                                               │
                      (cached, reusable across runs)           │
                                                               ▼
                                    (stage 2: multi-sample variant grouping)
                                                               │
                                                               ▼
                                      OverlappingVarGroup iterator
                                                               │
                                                               ▼
                                    (stage 3: allele merging across samples)
                                                               │
                                                               ▼
                                   (stage 4: per-sample genotype likelihoods
                                     computed from read evidence)
                                                               │
                                                               ▼
                                         (stage 5: EM + posteriors)
                                                               │
                                                               ▼
                                      (stage 6: genotype assignment,
                                         missing / partial handling)
                                                               │
                                                               ▼
                                                         [multi-sample VCF]
```

## Stage 1 — Per-sample BAM analysis (pileup + candidate sites)

**Purpose:** Read one BAM, emit a sorted stream of candidate variant sites with all per-allele evidence needed for downstream likelihood computation.

Is this step similar to samtools pileup?[QUESTION]
Could we use or should we use samtools rust-htslib?[QUESTION]


### What information do we need per candidate site?

Minimum required:
- `chrom`, `pos`, `ref_allele` (from reference FASTA)
- For each observed allele at this site (including ref):
  - read count (total depth)
  - count on forward strand / reverse strand (for strand-bias detection later)
  - list of base qualities (or a summary sufficient to compute likelihoods)
  - list of mapping qualities (or summary)

Probably also:
- total depth at the site (including filtered reads, for DP field)
- count of reads filtered out (low MAPQ, low BQ, duplicate, secondary, etc.)

**[QUESTION]** Do we store per-read values (bq, mq) as flat lists, or do we summarize upfront (e.g., precomputed log-likelihood contributions per allele)? Flat lists give us flexibility later at the cost of file size; summaries are compact but lock in modeling choices.

### Handling deletions

For a read supporting a deletion of length `k` starting at position `p`:
- At position `p`, the read contributes the deletion allele (ref starts with base at `p`, alt is just that base).
- At positions `p+1` .. `p+k`, the read is "inside" the deletion; no base observed.

Per the user's instruction: the quality assigned to those interior deleted positions should come from **flanking base qualities** in the alignment (e.g., average or min of the base before and the base after the deletion).

**[QUESTION]** Flanking quality: min, mean, or something else? GATK uses min for BQSR-like contexts. Let's pick one and document why.

### Handling insertions

For a read supporting an insertion `+XY` after position `p`:
- At position `p`, the read contributes `ref_base + XY` (extended allele).
- No positions are "skipped."
- Base quality for the inserted bases comes from the read itself (directly in the BAM).

### Variant discovery strategy

**[DECISION]** Which discovery model?

1. **Naive pileup (bcftools/samtools style):** At every position, if any read shows a non-ref base, consider it a candidate. Simple, fast, works well for SNPs and short indels.
2. **Active-region + local reassembly (GATK HaplotypeCaller style):** Detect noisy regions, locally reassemble haplotypes, realign reads. Much more accurate for indels and complex regions, much more complex to implement.
3. **Haplotype-based with partial local realignment (freebayes style):** Enumerate candidate haplotypes over short windows, evaluate jointly.

Naive pileup is the natural starting point; (2) or (3) are realistic later upgrades. Recommend starting with (1) and designing the binary format to support richer information later.
At this first point, having information of only one sample, Jose Blanca thinks that naive pileup is enough, local realigments or reasemblies will be more effective once we are in the merging step with all samples.

### Filters (per read, before counting)

Standard filters to consider:
- Skip unmapped, secondary, supplementary alignments
- Skip duplicates (flag set)
- Minimum MAPQ (configurable, common default: 20)
- Minimum BQ (configurable, common default: 13 or 20)
- Optional: maximum read depth (cap to bound memory / runtime)
- Optional: BAQ (base alignment quality) recalibration — bcftools/samtools does this; GATK doesn't because it does local realignment

**[DECISION]** Do we implement BAQ? It helps a lot near indels for callers without local realignment, but it's extra work. For a first cut: skip BAQ, document the limitation.

### BAM library choice

**[DECISION]** Rust has two mature BAM/SAM libraries:
- **rust-htslib** — bindings to htslib. Proven, fast, standard. Requires C toolchain.
- **noodles** — pure Rust, supports a wide range of bioinformatics formats. More idiomatic Rust, might be slower than htslib for some operations.

Recommend rust-htslib for battle-testing; noodles if we want zero C dependencies.

### Reference FASTA access

Need indexed FASTA (`.fai`). Both rust-htslib and noodles support faidx.

### Output of stage 1

A sorted stream of candidate sites, written to the per-sample binary file (stage 2).

## Stage 2 — Per-sample binary file format

**Purpose:** Cache the pileup result so that it's computed once per sample and reused across all cohort calls the user does in the future.

### Requirements

- Sorted by (chrom, pos) — ideally following the BAM header's chromosome order, which is also the variant-grouping order.
- Streaming-friendly: we read sequentially, no need to load whole file into memory.
- Optional random access by region (chrom:start-end) — useful for parallelism and `--chroms` flag compatibility. Implies some kind of index. (Having indexes inside the binary file to be able to parallelize by chrom could be useful.)
- Compact: binary, compressed.
- Schema-evolvable: we want to be able to add fields later without invalidating old files (or at least detect version mismatch and re-run stage 1). The binary files should include a version number to be future proof.

### Format options

**[DECISION]** Pick one:

1. **Custom binary with zstd compression** — full control, minimum dependencies, but we implement framing/versioning ourselves.
2. **Apache Arrow IPC / Parquet** — columnar, great compression, standard tooling, but a heavy dependency and overkill if no one else reads it.
3. **Cap'n Proto / FlatBuffers** — schema-driven, fast, zero-copy reads.
4. **bgzf-framed records with a tabix-compatible index** — reuses ecosystem tooling; we can essentially write a custom binary record format but use the same indexing as VCF/BAM.

Given the user's framing ("no one else reads this file, we want it fast"), option **1 (custom binary + zstd)** is probably the right default. Option **4** is nice if we want to leverage existing indexing.

### Suggested record schema (draft)

```
file header:
  magic bytes (e.g. "PSP\0")
  format version
  reference name + md5 (to detect FASTA mismatches)
  BAM sample name
  chromosome list in order

per record (one per candidate site):
  chrom_id: u32
  pos: u32 (1-based)
  ref_base: u8 (or ref_allele_len + bytes, for multi-base refs from indels)
  n_alleles: u8
  for each allele:
    allele_bytes: variable-length
    fwd_count: u16
    rev_count: u16
    bq_summary: (either list of u8, or min/mean/max/sum_of_squares for moment-based likelihood)
    mq_summary: similar
  total_depth_raw: u32
  total_depth_filtered: u32
```

Not final — just to make the conversation concrete.

### Naming

**[QUESTION]** What do we call the file extension? Candidates: `.psp` (per-sample pileup), `.alnev` (alignment evidence), `.sce` (sample candidate evidence). The user named the project "join_per_sample_vcfs" so there's a naming theme to extend.

## Stage 3 — Multi-sample variant grouping

**Purpose:** Take N per-sample binary files, merge their candidate sites into `OverlappingVarGroup`s.

The existing `variant_grouping.rs` does exactly this for `VarIterator<Variant>`. We need a new iterator that yields the same `OverlappingVarGroup` shape but reads from the binary pileup files.

**[DECISION]** Two options:
**Introduce a new grouped type** (e.g. `OverlappingPileupGroup`), inspired by the current OverlappingVarGroup, that carries richer evidence per sample/allele. Cleaner, but duplicates the grouping logic.

## Stage 4 — Full allele merging across samples

**Purpose:** Given an `OverlappingVarGroup`, produce one merged variant covering the whole region, with a unified allele list that correctly represents all samples' observations — including complex overlapping deletions.

This is closer to what the current `genotype_merging.rs` does (one output per group, haplotype-level string concatenation) than to the GATK joining approach (one output per position with `*` spanning alleles).

In this step we could consider reassembling.

### Differences from current `genotype_merging.rs`

- Current code starts from per-sample GT calls. The new code starts from per-sample *allele evidence* (read counts per allele). No GTs have been called yet — that happens in stages 5-6.
- Current code assembles haplotypes from called alleles. The new code needs to consider all observed alleles (each with its evidence), not just the top-called allele.
- Complex overlapping deletions: when sample A has a 3bp deletion and sample B has a SNP at position 2, the merged allele list needs to contain both the full-length deletion allele (from A) and the SNP-extended-to-match-length allele (from B), and evidence from each sample maps correctly.

### Pseudocode sketch

```
fn merge_alleles(group) -> MergedSite:
    collect all distinct alleles across samples, normalized to the longest ref length in the group
    unify equivalent alleles (e.g. an A->G SNP at the first base of a 3bp ref becomes "GNN" where NN is ref suffix)
    for each sample:
        remap the sample's observed-allele evidence to the unified allele list
    return MergedSite { alleles, per_sample_evidence }
```

### Open questions

When sample A has a deletion covering positions 1-3 and sample B has a SNP at position 2, does the merged variant emit:
- A single multi-base record with alleles `[AGC, A, AAC]` (deletion + SNP-extended)? (merge style)

**[QUESTION]** How do we handle reads that only partially span the merged region (e.g., reads that cover position 1 of a 3bp deletion but end before position 3)? These can support the deletion allele (if they show the deletion at pos 1) but don't give evidence for ref across the whole region. Options: (a) count them as supporting their observed allele fragment, (b) discard them for this site.

## Stage 5 — Genotype likelihoods from read evidence

**Purpose:** Per sample, per possible genotype at the merged site, compute P(reads | genotype). These replace the PL values we currently read from gVCFs.

### Model (standard, GATK / bcftools-style)

For a diploid sample with observed reads R = {r_1, r_2, ...} and a candidate genotype G = (a1, a2) with alleles a1, a2:

```
P(R | G) = ∏_r P(r | G)
         = ∏_r [0.5 * P(r | a1) + 0.5 * P(r | a2)]

P(r | allele a) =
    (1 - err_r)    if r supports a
    err_r / 3      otherwise   (per-base error distributed uniformly across 3 other bases)
```

where `err_r = 10^(-BQ_r / 10)` is the per-base error probability.

For ploidy > 2, the mixing weights are multinomial coefficients.

For indel alleles, per-base error no longer captures the full story (gap-open / gap-extend matter). A simple first cut: treat indel alleles with a fixed per-read indel error rate (e.g. 10^-4) when computing `P(r | indel_allele)`; refine later.

### What gets stored vs computed on the fly

Because we kept per-allele read counts and BQ summaries in the binary format, we can compute P(reads|G) for any genotype at any site on the fly. This gives us PLs for all enumerated genotypes, directly from read evidence.

**[QUESTION]** Summaries vs raw lists. If we stored `(fwd_count, rev_count, list_of_BQs)` for each allele, we get exact GATK-equivalent likelihoods. If we stored only moment summaries (mean BQ, count), we can approximate. Raw lists are more flexible but balloon file size on high-depth samples. Reasonable compromise: store the counts at each BQ bin (histogram by BQ, e.g., 0-40 binned into 5-unit buckets).

## Stage 6 — EM and posteriors across the cohort

The existing `genotype_posteriors::estimate_posteriors` already does this. It takes flat PLs per sample, iterates allele frequencies, emits posteriors + QUAL. We reuse it essentially unchanged.

## Stage 7 — Final genotype assignment and missing-data handling

For each sample at each merged site:
- Best posterior → GT.
- GQ = phred-scaled gap to second-best posterior (cap 99).
- If depth is below a threshold OR best posterior is below a threshold, emit missing: `./.` (or partial missing `0/.` when only one haplotype has sufficient info — do we want partial missing? phase-aware?).
- QUAL = site-level phred-scaled probability of being non-variant.

**[QUESTION]** Per-haplotype missing (`0/.`) is tricky without phased reads. For unphased pileups, we probably only emit fully-missing (`./.`) or fully-called (`0/1`, `1/1`, etc.) — not partial. For phased read data (linked reads, long reads) we could infer per-haplotype, but that's a much later feature.

## Stage 8 — VCF output

Extend the existing `VcfWriter` to emit:
- `GT` (already supported)
- `GQ` (phred-scaled confidence)
- `DP` (total depth)
- `AD` (allele depths, per allele)
- `PL` (genotype likelihoods, rescaled to min-zero)
- QUAL column populated from `SitePosteriors.qual`
- Per-site INFO fields (optional, later): `AC`, `AN`, `AF`, `MQ`, strand bias metrics.

## Reuse inventory from current codebase

| Current code                              | Reuse? | Notes                                                                 |
|-------------------------------------------|--------|-----------------------------------------------------------------------|
| `variant_grouping.rs`                     | Yes    | Works unchanged if we yield `Variant` from the new pileup iterator.   |
| `genotype_posteriors.rs` (EM core)        | Yes    | Same EM, just better PL input.                                        |
| `vcf_writer.rs`                           | Extend | Needs GQ/PL/DP/AD/QUAL.                                               |
| `genotype_merging.rs`                     | Inspire| Allele-merging logic is conceptually similar but starts from evidence, not GTs. Likely a rewrite, informed by this.      |
| `genotype_joining.rs` (planned, paused)   | Drop   | Superseded by this plan.                                              |
| `gvcf_parser.rs`                          | Keep   | Still useful for optional gVCF-input mode; not on the critical path.  |
| `decompression_pool.rs`, `threaded_reader.rs` | Keep   | Generic I/O infrastructure, reusable.                                 |

## Suggested module layout (first pass)

```
src/
  bam_pileup.rs           — stage 1: BAM → per-site evidence stream
  sample_evidence.rs      — stage 2: binary format reader/writer, schema, version
  variant_grouping.rs     — stage 3: unchanged
  allele_merging.rs       — stage 4: full-merge across samples, evidence remap
  genotype_likelihoods.rs — stage 5: P(reads | G) computation
  genotype_posteriors.rs  — stage 6: existing EM, unchanged
  variant_calling.rs      — stage 7: orchestrator; GT/GQ/missing logic
  vcf_writer.rs           — stage 8: extended with GQ/PL/DP/AD
  pipeline.rs             — top-level pipeline wiring
  main.rs                 — CLI: two subcommands, one for stage 1 (per-sample),
                            one for stages 3-8 (cohort call)
```

CLI suggestion:
- `tool pileup --bam sample1.bam --fasta ref.fa --out sample1.psp`
- `tool call sample1.psp sample2.psp ... --fasta ref.fa > out.vcf`

## Things to study from GATK / freebayes

Clone and read (not copy). Useful references:
- **GATK HaplotypeCaller:** `src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/` — active regions, PairHMM, local assembly.
- **GATK GenotypeGVCFs:** `GenotypeGVCFs.java` and `GenotypeGVCFsEngine.java` — the joining logic we've already studied.
- **GATK AlleleLikelihoods / GenotypeLikelihoodCalculator:** how PLs are computed from reads.
- **freebayes:** `src/` — Bayesian model, haplotype windows. Concise C++.
- **bcftools mpileup + call:** `bam2bcf.c`, `ccall.c` — classic per-site pileup + per-site calling.
- **noodles examples:** for idiomatic Rust BAM iteration.

**[ACTION]** Clone into `third_party/` (gitignored) for reference:
```
git clone --depth 1 https://github.com/broadinstitute/gatk.git
git clone --depth 1 https://github.com/freebayes/freebayes.git
git clone --depth 1 https://github.com/samtools/bcftools.git
```

## Risks and unknowns

1. **Indel calling is hard.** Naive pileup does SNPs well, indels less well. If indel accuracy matters, we'll need local realignment or haplotype assembly, which is a big subproject.
2. **Performance.** Per-site P(reads|G) over millions of sites × thousands of samples can be slow. We'll need a clear perf target early (e.g., "1000-sample cohort, whole-genome, in < N hours").
3. **BAM I/O is expensive.** We should benchmark the BAM pileup stage early — it may dominate runtime, in which case parallelism / chunking matters a lot.
4. **Reference FASTA footprint.** Need streaming access, not load-all; both rust-htslib and noodles support this.
5. **Binary format evolution.** If we change the schema after collecting many `.psp` files, users will need to regenerate. Version field + clear error message is essential.
6. **Scope creep.** This is effectively a new variant caller. Want to scope v1 tightly — e.g., "diploid, SNPs + short indels, naive pileup, no BAQ, no active regions" — and iterate.

## Open decisions for the next discussion

Tagged in the doc above as **[DECISION]** and **[QUESTION]**. Top ones to resolve before drafting a detailed plan:

- [DECISION] BAM library: rust-htslib vs noodles.
- [DECISION] Discovery model: naive pileup vs active regions vs haplotype windows.
- [DECISION] Binary format: custom / Arrow / Cap'n Proto / bgzf+tabix.
- [DECISION] BAQ: yes/no for v1.
- [QUESTION] Evidence summaries vs raw lists in the binary format.
- [QUESTION] Merge-style output for overlapping deletions — confirm this is desired over join-style.
- [QUESTION] Naming: file extension, CLI subcommand names.
- [QUESTION] Scope for v1: SNPs only, or SNPs + indels, or full?
- [QUESTION] Partial missing genotypes (`0/.`) — support or not?

## Next step after this draft

Review together, resolve the decisions above, then I'll produce a detailed per-stage plan (on the same shape as the previous `genotype_joining.md`: function signatures, pseudocode, test strategy) for the stages we agree to tackle first.
