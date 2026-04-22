# Per-Sample Pileup (PSP) — Binary Format Specification

**Status:** Draft derived from the 2026-04-21 design discussion.
This draft should be completely reconsidered once we decide which data we need to calculate the genotype posterior probabilities.

## Purpose

Defines the on-disk format produced by Stage 1 of the SNP calling pipeline (see `snp_calling_from_bam.md`). One PSP file per sample caches the candidate sites and per-read evidence extracted from a single BAM so that cohort-level variant calling can be re-run without revisiting the BAM.

## Design goals

- Compact on low-coverage data. Typical target: 2–10× per sample.
- Streamable sequentially — the common cohort-call pattern walks the whole genome in chrom/pos order.
- Seekable by genomic region for parallelism and `--chroms` filtering.
- Schema-versioned so new fields can be added without silently breaking old files.
- Written once per sample, read many times.

## Scope

In scope: site-level and per-read evidence needed for downstream genotype likelihoods, standard GATK-style rank-sum filters (MQRankSum, ReadPosRankSum, BaseQRankSum), strand bias filters, and cross-sample allele merging.

Out of scope for v1: haplotype blocks, long-read phase, per-read linkage tags.

## Data model

Granularity is three nested levels: site → allele → read.

The format stores only independent data — anything derivable from the per-read records (depth, strand counts, sum of MQ², RMS MQ, etc.) is computed on read, not stored. At the project's target coverage, iterating reads is cheap and downstream stages already iterate them for likelihood computation.

### Site-level fields

One record per candidate site.

| Field      | Type      | Notes                                                                         |
|------------|-----------|-------------------------------------------------------------------------------|
| chrom_id   | u32       | Index into the file header's chromosome table.                                |
| pos        | u32       | 1-based reference start of the site's ref allele.                             |
| ref_allele | var-len   | 1-byte length prefix + sequence bytes. One base for SNPs, longer for indel-extended refs. |
| n_alleles  | u8        | Number of per-allele records that follow.                                     |

Pre-filter depth and filter-reason counters are deliberately not stored. Dropped reads (duplicates, low-MAPQ, secondary, etc.) are not recoverable from the file, but none of the calling stages use that information. If QC reporting later needs pre-filter counts, it can be added behind a feature flag without a breaking version bump.

### Per-allele fields

One entry per observed allele at the site. The first entry is always the reference allele (including when no read supports a variant).

| Field      | Type     | Notes                                                                      |
|------------|----------|----------------------------------------------------------------------------|
| allele_seq | var-len  | 1-byte length prefix + sequence bytes.                                     |
| ref_span   | u8       | Reference bases consumed by this allele (≥ 1).                             |
| n_reads    | u16      | Number of per-read records that follow.                                    |

Strand counts (fwd/rev), allele depth, and other per-allele summaries are computed from the per-read records on read. Partial-span reads, if tracked separately, would be a format extension — see open decisions.

### Per-read fields

One 3-byte packed record per read supporting this allele:

| Byte | Bits | Field         | Notes                                                 |
|------|------|---------------|-------------------------------------------------------|
| 0    | 8    | BQ            | Base quality at the site (Phred 0–93).                |
| 1    | 8    | MQ            | Read mapping quality (Phred; usually 0–60).           |
| 2    | 7    | pos_in_read   | Position-in-read, clamped or binned to 0–127.         |
| 2    | 1    | strand        | 0 = forward, 1 = reverse.                             |

Strand is duplicated with the per-allele fwd/rev counts. Keeping it in the per-read record makes records self-contained and simplifies downstream per-read iteration.

### Size estimate at 5× coverage

- Site header: ~11 B.
- Per-allele header: ~6 B (SNP case).
- Per-read record: 3 B.

Typical site with 2 alleles and 5 total reads: `11 + 2·6 + 5·3 ≈ 38 B` pre-compression.

### Indel handling

**Deletions.** For a read supporting a deletion of length k starting at position p:
- The read contributes a per-read record to the deletion allele at position p.
- The read contributes no evidence at positions p+1 … p+k.
- BQ for the deletion record: **[OPEN]** min or mean of the base qualities immediately flanking the deletion, vs a fixed indel-error proxy.

**Insertions.** For a read supporting an insertion of bases XY after position p:
- The read contributes a per-read record to the allele `ref_base + XY` at position p.
- BQ: **[OPEN]** min or mean of the BQs across the inserted bases.

## File structure

```
+-------------------+
| file header       |
+-------------------+
| block 0           |
+-------------------+
| block 1           |
+-------------------+
| ...               |
+-------------------+
| block N-1         |
+-------------------+
| block index       |
+-------------------+
| footer (20 B)     |
+-------------------+
```

### File header

| Field              | Type            | Notes                                                               |
|--------------------|-----------------|---------------------------------------------------------------------|
| magic              | 4 B             | `PSP\0`.                                                            |
| format_version     | u16             | Starts at 1. Bumped on breaking schema change.                      |
| flags              | u16             | Reserved; 0 in v1. Used for optional feature bits.                  |
| sample_name        | var-len utf-8   | Length-prefixed.                                                    |
| reference_name     | var-len utf-8   | Length-prefixed.                                                    |
| reference_md5      | 16 B            | MD5 of the FASTA used, for mismatch detection at cohort-call time.  |
| chrom_table        | var-len         | u32 count, then per chrom: u32 name length, name bytes, u32 length_bp. |
| creation_timestamp | u64             | Unix seconds.                                                       |
| writer_version     | var-len utf-8   | Tool version that wrote the file.                                   |

### Block layout

Each block is a standalone zstd frame. Records never cross block boundaries, and blocks never cross chromosome boundaries.

| Field             | Type   | Notes                                              |
|-------------------|--------|----------------------------------------------------|
| compressed_size   | u32    | Bytes of zstd payload that follow.                 |
| uncompressed_size | u32    | Expected size after decompression.                 |
| zstd_payload      | bytes  | Zstd frame containing a concatenation of site records. |

**Target uncompressed block size:** 4 MB as the starting default, configurable. Larger blocks are preferred when RAM is not scarce (the common case); they improve zstd ratio and have negligible seek cost for any realistic region query. Test 16 MB before freezing the default.

**Compression level:** **[OPEN]** Default zstd 9 as a starting point. Files are written once per sample and read many times, so a slower-write / better-ratio tradeoff is justified. Benchmark 3 / 9 / 19 on representative samples.

**Dictionary:** not used in v1. Revisit only if small-block compression ratio turns out to be poor.

### Block index

Written immediately after the last data block. Prefixed by `u64 n_blocks`.

| Field             | Type | Notes                                           |
|-------------------|------|-------------------------------------------------|
| chrom_id          | u32  |                                                 |
| first_pos         | u32  |                                                 |
| last_pos          | u32  |                                                 |
| file_offset       | u64  | Byte offset of this block's `compressed_size` field. |
| compressed_size   | u32  |                                                 |
| uncompressed_size | u32  |                                                 |
| n_sites           | u32  |                                                 |

Followed by `u32 index_crc` (CRC32C over the index bytes) for corruption detection.

### Footer

Fixed 20 bytes at end of file:

| Field        | Type | Notes                                            |
|--------------|------|--------------------------------------------------|
| index_offset | u64  | Byte offset of the block index start.            |
| index_size   | u64  | Index length in bytes.                           |
| magic        | 4 B  | `PSP\0` again, for tail-based format detection.  |

A reader pread's the last 20 bytes, validates magic, reads the index, then seeks to the blocks it needs.

## Access patterns

### Sequential read (full cohort call)

1. Open file, read header.
2. Read blocks in order; decompress each and emit its site records.
3. Stop at the start of the block index.

The block index is not needed for sequential reads. This matches the common case.

### Random access (region query)

1. Open file, read header and footer.
2. Seek to `index_offset`, decode the block index.
3. Binary-search for blocks overlapping the requested `(chrom, start, end)`.
4. Decompress only those blocks; filter records inside.

### Parallel write

Stage 1 naturally pipelines as:

- BAM decode
- Pileup producing site records into an in-memory block buffer
- Compression pool running zstd on completed blocks
- Writer appending compressed blocks sequentially

Block independence means the compression pool has no ordering constraints beyond final write order.

## Design rationale

- **Raw per-read records, not histograms.** At the project's target coverage (2–10×), per-site read counts are small. Raw records (~3 B/read) come out smaller than fixed-size histograms (~50+ B per allele) and preserve exact per-read correlations between BQ/MQ for rank-sum filters. Quantization error from binning is not a concern the project cares about; size at low coverage is.
- **Normalized format — no derived fields.** Depth, strand counts, RMS MQ, and other summaries are recomputed from the per-read records on demand. This keeps the file smaller and removes a class of consistency bugs (summary vs raw disagreeing). Recomputation is cheap because downstream stages iterate the reads anyway to compute likelihoods.
- **Blocked zstd, not BGZF or zstd seekable format.** Our queries are always genomic-coord queries. A custom coord-native index is simpler and smaller than layering (chrom,pos) → byte-offset on top of zstd seekable. zstd replaces gzip for roughly 2× better ratio on our data.
- **Large blocks (≥ 4 MB).** Typical usage walks large ranges; RAM is abundant. Larger blocks give zstd more context and better ratio, with negligible decompression latency on modern CPUs.
- **Chromosome boundaries enforced.** Simplifies region queries and avoids spurious cross-chrom decompression at chrom transitions.
- **In-file tail index.** Atomic with content; no sidecar file to manage or lose. Reader locates it from the fixed footer.

## Open decisions

To be resolved before format v1 is frozen:

- **[OPEN]** Indel BQ: flanking-base stat (min or mean) for deletions, inserted-base stat (min or mean) for insertions, or a fixed indel-error proxy per read.
- **[OPEN]** Partial-span reads: record them as a distinct allele bucket, mark them with a flag in the per-read record, or drop them entirely.
- **[OPEN]** Reference context caching at each site (flanking bases, for microsatellite detection), or re-read from FASTA at cohort-call time.
- **[OPEN]** Default compression level (3 / 9 / 19), to be chosen by benchmarking.
- **[OPEN]** Block-level checksum: rely on zstd's built-in frame checksum, or add a dedicated u32 per block.
- **[OPEN]** Default block size (4 MB vs 16 MB), to be chosen by benchmarking.
- **[OPEN]** mmap vs pread on the read path (implementation detail, not format).

## Versioning policy

- Bump `format_version` on any breaking schema change.
- Readers must refuse to open files with a major version they do not understand.
- New optional fields can be gated behind bits in the header's `flags` u16 without bumping the major version, provided older readers ignoring the flag still decode correctly.

## Relation to the broader pipeline

This format is consumed by Stage 3 (`variant_grouping`) and later stages. Stage 3 reads PSP files in genomic order, aligning sites across samples into `OverlappingVarGroup`s. Stage 4 consumes the per-allele and per-read evidence to perform cross-sample allele merging. Stage 5 computes per-sample genotype likelihoods from the per-read BQ/MQ records. The format is independent of the BAM library choice in Stage 1 (rust-htslib vs noodles) — both libraries can produce the same record stream.
