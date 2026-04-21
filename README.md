# merge_per_sample_vcfs

A high-performance tool for merging per-sample gVCF files into a multi-sample VCF. Designed for thousands of input files with streaming I/O, without creating an intermediate variant database, and parallel processing.

## Usage

```
$ merge_per_sample_vcfs --chroms <chrom1,chrom2,...> [--threads <N>] <vcf1.gz> <vcf2.gz> ...
```

### Options

- `--chroms <list>` — Comma-separated chromosome names in processing order (required)
- `--threads <N>` — Number of threads for parallel merging (default: all cores)

### Arguments

- `<vcf.gz>` — One or more gzipped gVCF files (one per sample)

## Development Container

A [Containerfile](Containerfile) is provided for reproducible development with [podman](https://podman.io/). The image pins Rust (currently `1.95-bookworm`) and includes `bcftools`/`tabix` for inspecting VCF output, plus the [Claude Code](https://claude.com/claude-code) CLI for in-container agent sessions.

### Getting started

```
./scripts/dev.sh                 # interactive shell inside the container
./scripts/dev.sh cargo test      # run a one-off command
./scripts/dev.sh cargo build --release
```

The first invocation builds the image (~4 min). Subsequent runs reuse the cached image.

### Running Claude Code in the container

```
./scripts/claude.sh              # new Claude Code session
./scripts/claude.sh --continue   # resume last session
```

The container acts as the sandbox: Claude can freely read/write/execute inside the mounted project directory but cannot touch anything else on the host, so `--dangerously-skip-permissions` is passed by default. Host `~/.claude.json` and `~/.claude/` are mounted so login and per-project auto-memory persist.

### Design notes

- The project directory is bind-mounted at the **same absolute path** inside and outside the container. This keeps paths consistent for tools that embed `cwd` into state (notably Claude auto-memory).
- Container builds write to `target-container/` instead of `target/` so host and container don't invalidate each other's incremental compilation state.
- Rootless podman maps container UID 0 to the host user, so files created inside the container appear with the correct ownership on the host.

## Architecture

The pipeline flows through 5 stages:

```
gVCF files -> Decompression Pool -> Parsing -> Grouping -> Merging -> VCF output
```

### Modules

| Module | Role |
|---|---|
| `main.rs` | CLI entry point. Parses `--chroms`, `--threads`, and input paths |
| `pipeline.rs` | `merge_alleles_and_genotypes()` — orchestrates the full pipeline |
| `gvcf_parser.rs` | `VariantIterator<B>` — streaming parser with zero-alloc fast paths for common genotypes (`0/0`, `0|0`), LRU cache for FORMAT fields, 256KB I/O buffers |
| `variant_grouping.rs` | `VariantGroupIterator` — merges multiple streams by position. Phase A: parallel peek via rayon to find next position. Phase B: sequential consumption of overlapping variants into `OverlappingVariantGroup` |
| `genotype_merging.rs` | `merge_vars_in_groups()` — double-buffered: background thread collects groups into batches while main thread processes previous batch with rayon. Handles deletion spanning, phase tracking, allele unification |
| `vcf_writer.rs` | `VcfWriter` — writes VCFv4.2 to files (plain/gzip) or stdout. Silently handles broken pipes |
| `decompression_pool.rs` | `DecompressionPool` — N worker threads shared across all input files (avoids 1-thread-per-file explosion). 64KB chunks, 4-chunk bounded buffers per reader |
| `threaded_reader.rs` | `ThreadedReader` — simpler 1-thread-per-file alternative (predates pool) |
| `errors.rs` | `VcfParseError` — 23 variants via thiserror |
| `utils_magic.rs` | Gzip magic byte detection for input validation |

### Key Data Structures

- **`Variant`** — flat layout: `genotypes: Vec<i8>` (sample `i`, ploidy `p` at `[i*p..(i+1)*p]`), `-1` for missing
- **`OverlappingVariantGroup`** — all variants across samples that overlap a genomic span, with `source_var_iter_idxs` tracking which iterator each came from

### Parallelism Strategy

1. **Decompression**: Shared thread pool (N workers), bounded 4-chunk buffers per file
2. **Grouping**: Rayon parallel peek across all iterators, sequential consumption
3. **Merging**: Double-buffered — collector thread fills batch N+1 while rayon processes batch N (batch size: 1000 groups, channel capacity: 1)

### Genomic Algorithm Highlights

- **Allele unification**: Per sample alleles are merged and deduplicated into a shared allele index
- **Deletion handling**: Builds merged alleles covering overlapping deletions and SNPs to create complex variants
- **Phase tracking**: Detects broken phase chains (unphased het after prior het marks sample as missing)
- **Non-variable filtering**: Groups where all samples are hom-ref are discarded

### Dependencies

- `flate2` — gzip compression/decompression
- `rayon` — data parallelism for batch processing
- `ahash` — fast hashing for LRU cache
- `lru` — FORMAT field index cache
- `thiserror` — error type derivation
- `tempfile` — temporary file handling for tests
