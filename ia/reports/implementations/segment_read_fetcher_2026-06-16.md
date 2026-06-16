# Segment read fetcher — implementation report (2026-06-16)

Plan: [segment_read_fetcher.md](../../../doc/devel/implementation_plans/segment_read_fetcher.md).
Branch/worktree: `segment-read-fetcher`.

## 1. Plan

Built the standalone, thread-safe indexed-segment read source in
[segment_reader.rs](../../../src/bam/segment_reader.rs) (plan increments #1–#3):
`AlignmentFile` enum (`BamFile` / `CramFile`) with a pooled-reader design,
`get_reads_from_segment`, and `from_input`. Reuses the existing per-format
decode helpers; adds only reader ownership + pooling. The existing
`OwnedIndexed{Bam,Cram}Records` scanners are left untouched (they still back
`AlignmentMergedReader::query` until increment #5).

## 2. Assumptions (silent choices surfaced)

1. **Increments #4 (SSR wiring) and #5 (SNP `--regions` retrofit) deferred.**
   The SSR module does not exist on this branch (`src/ssr/` absent), so the
   primitive ships standalone with its own test suite, exactly as the plan's
   "standalone, independently testable primitive" intends. #4 is deferred for
   lack of a consumer here; #5 is its own measured change per the plan.
2. **No `FilterCounts` in the primitive.** It applies the pre-decode flag/MAPQ
   filter (`classify_pre_decode`) but does not tally dropped reads — no
   consumer needs counts here, and the SNP path keeps its own tallies.
3. **`from_input` is single-file, not `PileupInputs`-shaped.** It takes one
   input's already-validated handles (`path`, `Arc<sam::Header>`,
   `AlignmentIndex`, optional `fasta::Repository`, `cfg`, `source_file_index`);
   the caller extracts these from `PileupInputs` per input. This matches the
   plan's "one file → its reads" non-merge contract.
4. **Per-record pipeline** (mirrors what `AlignmentMergedReader` emits per
   read, minus the cross-file merge): target-contig check → early-stop once
   `alignment_start > segment.end` → `ContigInterval::overlaps_record` →
   `classify_pre_decode` → `min_read_length` → `record_buf_to_mapped_read`.
   The downstream G2/F1/F3 read processing is intentionally *not* run here (it
   is the consumer's stage), same as the merged reader.
5. **`#![allow(dead_code)]` on the module**, documented, because the only live
   consumers are the tests until #4/#5 land (same convention as
   `FilterBucket::BaqRejected`). Remove when the consumer is wired.

## 3. Changes made

- **New** [src/bam/segment_reader.rs](../../../src/bam/segment_reader.rs):
  - `AlignmentFile { Bam(BamFile), Cram(CramFile) }` — `Sync`, shared by `&`.
  - `BamFile` / `CramFile` — `path`, `Arc`-shared index, `Arc<HeaderRefMap>`
    (bundles the `Arc<sam::Header>` the decoder needs with a `name → ref_id`
    map), `cfg`, `source_file_index`, and a `Mutex<Vec<Handle>>` reader pool
    (CRAM also holds the `fasta::Repository`).
  - `from_input(...)` — dispatches on file extension × index variant; typed
    errors for format/index mismatch and a CRAM with no repository.
  - `get_reads_from_segment(chrom, start, end)` — resolves the ref-id,
    (BAM) runs `BinningIndex::query`, eagerly borrows a pooled reader, returns
    `MappedReadsInSegment<'_>`. On `Drop` the iterator returns its reader.
  - `BamSegmentReads` / `CramSegmentReads` — the per-call streaming iterators
    (chunk-walk / `.crai` container-walk) with the per-record pipeline above.
- **Visibility lifts** so the primitive reuses existing logic (no fork):
  - [alignment_input.rs](../../../src/bam/alignment_input.rs): `classify_pre_decode`,
    `PreDecodeOutcome`, `record_buf_to_mapped_read` → `pub(super)`.
  - [bam_input.rs](../../../src/bam/bam_input.rs): `query_interval`,
    `open_bam_reader_with_header` → `pub(super)`.
  - [cram_input.rs](../../../src/bam/cram_input.rs): `open_cram_reader_with_header`
    → `pub(super)`.
- **New errors** in [errors.rs](../../../src/bam/errors.rs): `MissingCramReference`,
  `InvalidSegment`.
- **Module registration** in [mod.rs](../../../src/bam/mod.rs).

## 4. Tests added (18, all in the module)

- `bam_returns_exactly_overlapping_reads_with_inclusive_edges` + source-file
  stamping; `bam_segment_boundaries_are_one_based_inclusive` (4↔4 in, 5↔5 out,
  1↔1 in).
- `bam_read_spanning_two_segments_is_yielded_by_both_whole`.
- `bam_empty_segment_yields_nothing_and_returns_handle`.
- `bam_pool_opens_once_and_returns_to_resting_size` (pool 0 fresh → 0 during
  borrow → 1 after drop, across 5 calls; never grows).
- `bam_parallel_segments_match_sequential` (`par_iter` over a shared
  `&AlignmentFile` equals the sequential run — the `Sync` design proof).
- `bam_target_contig_filter_excludes_other_contigs`,
  `bam_low_mapq_read_is_filtered_under_default_config` (MAPQ gate isolated),
  `bam_invalid_segment_is_rejected`, `bam_unknown_contig_is_rejected`.
- CRAM analogues: overlap, spanning, empty+handle-return, pool reuse, parallel
  determinism (`.crai` built via `cram::fs::index`).
- `cram_without_repository_is_rejected`, `from_input_rejects_format_index_mismatch`,
  `alignment_file_is_sync` (compile-time `Sync` assertion).

## 5. Validation results

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib --tests` — 1056 lib + all integration tests pass (the 18
  new tests included); 1 pre-existing ignored.
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib` — clean.

## 6. Tradeoffs and follow-ups

- **Chunk-walk duplication.** `BamSegmentReads` / `CramSegmentReads` repeat the
  chunk/container-walk shape of `OwnedIndexed{Bam,Cram}Records` rather than
  refactoring them, because those still back `AlignmentMergedReader::query`.
  Increment #5 (the SNP `--regions` retrofit) is where the two converge.
- **No CRAM container cache** (plan §8) — each segment re-decodes the
  containers it touches. The reader is still pooled (no re-open). Revisit only
  if profiling on real SSR loci says so.
- **CRAM early-stop is container-granular** — the `.crai` walk skips containers
  disjoint from the segment but does not early-terminate the contig scan; BAM
  has the per-record `alignment_start > end` early-stop. Acceptable for tiny
  SSR loci; tighten if a hot path appears.
- **Deferred:** #4 SSR `fetch_locus_reads` wiring; #5 SNP `--regions` retrofit
  (+ its measurement of the `--regions` tax improvement).
