---
name: Project specification and architecture
description: Complete specification of the merge_per_sample_vcfs Rust project — structure, modules, data flow, key types, algorithms, and test coverage
type: reference
---

# merge_per_sample_vcfs — Project Specification

## Purpose

A high-performance Rust CLI tool that merges per-sample gVCF files into a multi-sample VCF, replicating GATK's GenotypeGVCFs functionality. Designed for thousands of input files with streaming I/O and parallel processing.

## CLI Interface

```
merge_per_sample_vcfs --chroms <chrom1,chrom2,...> [--threads <N>] <vcf1.gz> <vcf2.gz> ...
```

- `--chroms` (required): comma-separated chromosome names defining processing order
- `--threads` (optional, default: all cores): thread count for rayon and decompression pool
- Input files must be gzipped gVCF files (one per sample)
- Output: VCFv4.2 written to stdout

## Pipeline Architecture

```
gVCF files -> DecompressionPool -> Parsing -> Grouping -> Joining -> VCF output
```

### 5-Stage Flow

1. **Decompression** (`decompression_pool.rs`): N worker threads shared across all files. 64KB chunks, 4-chunk bounded buffers per reader. Avoids 1-thread-per-file explosion.
2. **Parsing** (`gvcf_parser.rs`): Streaming `VariantIterator<B>` per file. Zero-alloc fast paths for 0/0, 0|0. LRU cache for FORMAT fields. 256KB I/O buffers. Parses GT, PL, phase, alleles. Strips `<NON_REF>` and `.` from alleles.
3. **Grouping** (`variant_grouping.rs`): `VariantGroupIterator` merges multiple streams by position. Phase A: parallel peek via rayon. Phase B: sequential consumption. Produces `OverlappingVariantGroup` bins. Non-variable groups (all samples hom-ref) are skipped.
4. **Joining** (`genotype_joining.rs`): Double-buffered — background thread collects groups into batches of 1000 while main thread processes previous batch with rayon. Handles deletion spanning, phase tracking, allele unification, EM-based genotype posteriors.
5. **Output** (`vcf_writer.rs`): `VcfWriter` writes VCFv4.2 to files (plain/gzip) or stdout. Silently handles broken pipes.

Orchestrated by `pipeline.rs` (`merge_alleles_and_genotypes()`).

## Source Modules

### src/main.rs
CLI entry point. Parses args, creates `DecompressionPool`, builds `VariantIterator` per file, calls `merge_alleles_and_genotypes()`.

### src/lib.rs
Declares public modules: `decompression_pool`, `errors`, `genotype_merging`, `genotype_posteriors`, `gvcf_parser`, `pipeline`, `threaded_reader`, `utils_magic`, `variant_grouping`, `vcf_writer`.

### src/pipeline.rs
`merge_alleles_and_genotypes<B>(var_iters, sorted_chromosomes, writer, prior)` — orchestrates grouping, merging, and writing.

### src/gvcf_parser.rs
- **`Variant`** struct: core data type. Fields: `chrom`, `pos`, `alleles` (ref first), `ref_allele_len`, `qual`, `genotypes` (flat `Vec<i8>`, -1=missing), `phase` (per sample bool), `pls` (flat `Vec<f64>`), `pls_per_sample`, `n_samples`.
- **`VariantIterator<B: BufRead>`**: streaming iterator with peek support. Constructors: `from_reader`, `from_gzip_reader`, `from_gzip_path`, `from_gzip_path_threaded`, `from_gzip_path_pooled`.
- Uses `CommonGenotypePatterns` for fast-path matching at any ploidy.
- `parse_format_field()`: extracts GT and PL indices from FORMAT field, cached via LRU.
- `parse_genotypes()`: handles fast paths (0/0, 0|0, ./.), general parsing with ploidy validation.
- `parse_pls()`: parses PL values from sample fields.

### src/variant_grouping.rs
- **`OverlappingVariantGroup`**: `chrom`, `start`, `end`, `variants: Vec<Variant>`, `source_var_iter_idxs: Vec<usize>`.
- **`VariantGroupIterator<B>`**: merges multiple streams. `compute_next_span_seed()` uses rayon parallel peek. `group_overlapping_vars()` sequential consumption. Validates chromosome order and position sorting. Detects duplicate sample names.
- **`VariantIteratorInfo`**: stores sample names per iterator.

### src/genotype_joining.rs
To be written.

### src/genotype_posteriors.rs
Implements GATK-style Bayesian genotype calling:
- **`PriorConfig`**: `snp_heterozygosity` (0.001), `indel_heterozygosity` (1.25e-4), `heterozygosity_std_dev` (0.01), `fixation_index` (0.0, Wright's F).
- **`Genotype`**: `allele_counts`, `pl_index`, `log10_multinomial_coefficient`.
- **`SitePosteriors`**: `allele_frequencies`, `genotype_posteriors` (flat), `num_genotypes`, `qual`.
- `enumerate_genotypes(num_alleles, ploidy)`: VCF PL-order enumeration. Works for any ploidy.
- `synthetic_pls_from_gt()`: generates synthetic PLs when VCF has no PL field (99% confidence on called GT).
- `estimate_posteriors()`: full EM algorithm. Dirichlet prior pseudocounts (ref=10, alt=0.01). Up to 50 iterations, convergence at max_change < 1e-6. Incorporates fixation index F for selfing species.
- `compute_genotype_posteriors()` (E-step): combines log10 likelihood (from PL), Hardy-Weinberg prior with F, multinomial coefficient.
- `genotype_log10_prior()`: mixture model: (1-F)*P_hw + F*P_inbred.
- `update_allele_frequencies()` (M-step): posterior-weighted allele counts + Dirichlet pseudocounts.
- `compute_qual()`: QUAL = -10 * log10(product of P(hom-ref) across all samples).

### src/decompression_pool.rs
- `DecompressionPool`: fixed N worker threads service all readers cooperatively via shared job queue. Each worker reads one 64KB chunk at a time, enqueues result, moves to next reader.
- `PooledReader`: implements `Read`, pulls from per-reader 4-chunk bounded buffer.
- Cancellation via `Drop`. Proper shutdown on pool `Drop`.

### src/threaded_reader.rs
To be removed.
Simpler 1-thread-per-file alternative (predates pool). `ThreadedReader` wraps a `Read` with background thread, 64KB chunks, bounded channel of 4.

### src/errors.rs
`VcfParseError` enum with 23 variants via thiserror. Includes `NotVariable` (used for filtering).

### src/utils_magic.rs
Gzip magic byte detection (`0x1f 0x8b`). `file_is_gzipped()` reads first 4 bytes.

### src/vcf_writer.rs
`VcfWriter`: writes VCFv4.2 header + data lines. Supports stdout, plain file, gzip file, or any `Write`. Handles broken pipes silently. FORMAT field is always just `GT`.

## Key Data Types

### Variant (central data structure)
```
chrom: String           // chromosome name
pos: u32                // 1-based position
alleles: Vec<String>    // [ref, alt1, alt2, ...] — no <NON_REF>
ref_allele_len: u8      // length of ref allele
qual: f32               // QUAL score (NaN if missing)
genotypes: Vec<i8>      // flat: sample_i * ploidy + haplotype_j. -1 = missing
phase: Vec<bool>        // per sample: true = phased (|), false = unphased (/)
pls: Vec<f64>           // flat PL values per sample per genotype
pls_per_sample: usize   // number of PL values per sample
n_samples: usize
```

### OverlappingVariantGroup
```
chrom: String
start: u32              // 1-based inclusive
end: u32                // 1-based inclusive
variants: Vec<Variant>  // all variants in the group
source_var_iter_idxs: Vec<usize>  // which iterator each variant came from
```

## Nomenclature (project-specific terms)
- **vallele**: an allele found in a variant (e.g., AC, A, AG are three valleles)
- **gallele**: an allele found in a genotype for a sample (one per chromosome)


## Test Coverage

### Unit Tests (in-module)
- `genotype_posteriors.rs`: 15 tests — genotype enumeration, multinomial coefficients, prior pseudocounts, EM convergence, QUAL scores, fixation index effects, synthetic PLs, tetraploid support.
- `vcf_writer.rs`: 14 tests — header format, genotype output, phasing, missing GTs, multi-sample, multi-variant, file output (plain/gz), broken pipe handling.
- `decompression_pool.rs`: 8 tests — single/multi/concurrent readers, cancellation, errors, gzip decompression.

### Integration Tests (tests/)
- `gvcf_parser_test.rs`: 15 tests — parsing, peek, gzip, record fields, allele filtering, sample parsing, genotype parsing.
- `variant_group_test.rs`: 7 tests — binning, deletion spanning, ordering violations, duplicate samples.
- `genotype_merging_test.rs`: 18 tests — simple merge, insertion, deletion (various), overlapping deletions, phase tracking (7 phase-specific tests), missing alleles, non-variable filtering.
- `integration_test.rs`: 12 end-to-end pipeline tests — two samples, all-ref, deletion+SNP, missing GT, multi-chromosome, overlapping deletions, phase problems.
- `utils_magic_test.rs`: 6 tests — gzip detection.

### Benchmark
- `benches/gvcf_perf.rs`: criterion benchmark parsing a large gVCF file.

## Dependencies
- `flate2`: gzip compression/decompression
- `rayon`: data parallelism
- `ahash`: fast hashing for LRU cache
- `lru`: FORMAT field index cache
- `thiserror`: error type derivation
- `tempfile`: test temp files
- `criterion` (dev): benchmarks

## Rust Edition
2024 (Cargo.toml `edition = "2024"`)

## Documentation Files
- `README.md`: project overview and architecture summary
- `genotype_joining_specification.md`: joining algorithm spec with examples, nomenclature, and the EM-based re-genotyping steps
- `information_about_the_problem_to_solve.md`: background on gVCF format, PL meaning, GQ, GATK joint VCF, * allele
- `posterior_gt_probs.md`: detailed description of GATK's EM algorithm for genotype posteriors, including call flow, math, intuitive explanation

## GATK repository

As a reference you have the GATK repository cloned in the directory ~/devel/gatk/
