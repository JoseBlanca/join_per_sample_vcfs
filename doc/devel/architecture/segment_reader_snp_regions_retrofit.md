# SNP `--regions` retrofit onto the segment reader — architecture sketch

**Status:** sketch for alignment, 2026-06-16 (branch `segment-read-fetcher`).
Pre-implementation outline of the types, enums, and functions. The "why"
and the increment sequence live in the plan
([segment_reader_snp_regions_retrofit.md](../implementation_plans/segment_reader_snp_regions_retrofit.md));
this doc is just the shapes. **Step 1 (filter-drop counts in the fetcher) is
already done** on `main` @ `3e5237e` — this sketch reflects the shipped API,
not the originally-proposed `RecordVerdict`.

## The idea, in one breath

Today the pileup re-opens each input file (and re-walks its index) **per
region** via `AlignmentMergedReader::query`. Replace that with **one pooled,
re-seekable [`AlignmentFile`](../../../src/bam/segment_reader.rs) per input,
opened once** and held for the whole run; per region, **k-way-merge the
per-input `MappedReadsInSegment` streams at the `MappedRead` level**. The
merge (coordinate order + cross-file dedup + within-file order check +
counts) survives — it just moves from `RecordBuf` to `MappedRead`, because
every key it needs is already a `MappedRead` field and the cheap filter
already ran inside the fetcher. Then `OwnedIndexed{Bam,Cram}Records` and the
old `query` body are deleted.

## Topology

```
PileupInputs (built once at startup)
  ├── files: Vec<AlignmentFile>     one pooled reader per input (Sync, open-once)
  ├── contigs, sample_name
  └── (headers/indexes consumed into the AlignmentFiles at build time)

per region  (serial loop in run_pileup, unchanged):
  PileupInputs::query_region(chrom, [start,end])
        │
        │  for each input file:  AlignmentFile::get_reads_from_segment(...)
        ▼
  SegmentMergedReads<'a>            ← NEW: the MappedRead k-way merge
   • peek each per-input stream
   • argmin by (ref_id, pos)        coordinate order
   • within-file order check        OutOfOrderRead + fuse
   • at-locus cross-file dedup      DuplicateReadAcrossFiles
   • (cheap filter already applied inside each fetcher)
   • counts() = Σ per-fetcher FilterCounts
        │  Iterator<Item = Result<MappedRead, AlignmentInputError>>
        ▼
  existing BAQ / read_processor / walker chain   (unchanged)
        │
        ▼  drop SegmentMergedReads → each MappedReadsInSegment returns its
           reader to its AlignmentFile's pool (re-seeked next region)
```

The region loop is serial, so each `AlignmentFile`'s pool holds ~one reader;
the merge borrows one segment iterator per file per region.

## Types

### Owning bundle (replaces today's handle vectors)

```rust
// src/bam/alignment_input.rs (or a new src/bam/pileup_inputs.rs)
pub struct PileupInputs {
    pub sample_name: String,
    pub contigs: ContigList,
    /// One pooled, re-seekable reader per input file, in input order.
    /// Built once by `load_pileup_inputs` (header + index + repository +
    /// SegmentReadFilter folded into each via `AlignmentFile::from_input`).
    files: Vec<AlignmentFile>,
    paths: Vec<PathBuf>,            // for merge error messages
}
```

### The merge (the only genuinely new struct)

```rust
// src/bam/segment_merge.rs  (new module; or beside the fetcher)
pub(crate) struct SegmentMergedReads<'a> {
    /// One peekable per-input segment stream — already cheap-filtered,
    /// coordinate-sorted, yielding `MappedRead`.
    streams: Vec<Peekable<MappedReadsInSegment<'a>>>,
    paths: &'a [PathBuf],

    /// Within-file coordinate-order guard, parallel to `streams`.
    per_file_prev_locus: Vec<Option<Locus>>,
    /// Cross-file dedup buffer for the current locus (cleared on advance).
    current_locus: Option<Locus>,
    current_locus_fingerprints: Vec<ReadFingerprintWithSourceFile>,

    fused: bool,                    // fuse-on-error
}

impl Iterator for SegmentMergedReads<'_> {
    type Item = Result<MappedRead, AlignmentInputError>;
    // argmin head → order check → dup check → emit; fuse on first Err.
}
```

Note: **no `FilterCounts` field on the merge** — the cheap filter already
ran in the fetchers, and dedup/order are *errors*, not drops. Counts are
summed from the streams (see below).

### Counts in the fetcher — already on `main` (`3e5237e`)

Plan step 1 is **done**. The fetcher already tallies cheap-filter drops; the
merge only needs to **sum** them. The shipped shape (not the `RecordVerdict`
I'd sketched — they threaded a `&mut FilterCounts` out-param instead, less
churn):

```rust
// src/bam/segment_reader.rs  (as merged on main)
fn classify_segment_record(
    record: &RecordBuf, target: usize, segment: &ContigInterval,
    filter: &SegmentReadFilter,
    counts: &mut FilterCounts,            // ← bumped on a filter drop
    source_file_index: usize, path: &Path,
) -> Result<Option<MappedRead>, AlignmentInputError>;
//   classify_pre_decode drop → counts.record_drop(bucket)
//   min_read_length  drop → counts.too_short += 1
//   out-of-segment (wrong contig / overlap miss) → Ok(None), UNCOUNTED

// each iterator owns a `filter_counts: FilterCounts` field; exposed as:
impl MappedReadsInSegment<'_> {
    pub(crate) fn filter_counts(&self) -> &FilterCounts;   // read after draining
}

// src/bam/alignment_input.rs  (as merged on main)
impl FilterCounts {
    pub(super) fn record_drop(&mut self, bucket: FilterBucket);  // central mapping
}
```

So the merge just aggregates:

```rust
impl SegmentMergedReads<'_> {
    /// Σ streams[i].filter_counts() — no merge-level drops exist
    /// (dedup/order are errors, not drops).
    pub(crate) fn filter_counts(&self) -> FilterCounts;
}
```

**Count parity is free.** Out-of-segment drops are uncounted, which is
exactly what the old reader did (its `OwnedIndexed*Records` pre-filter
ref_id+overlap *before* the merge counts), so the per-segment totals match
the old path by construction.

## Functions (the seam the driver calls)

Replace the 9-argument `AlignmentMergedReader::query` with a method on the
bundle that owns everything:

```rust
impl PileupInputs {
    pub fn load(
        alignment_files: &[PathBuf],
        fasta: &Path,
        build_index_if_missing: bool,
        filter: SegmentReadFilter,
    ) -> Result<Self, AlignmentInputError>;     // builds the pooled AlignmentFiles

    /// Reads overlapping `[start,end]` on `contig_name`, k-way-merged across
    /// this sample's inputs, coordinate-sorted. Borrows the pooled readers;
    /// the returned iterator releases them on drop.
    pub fn query_region<'a>(
        &'a self,
        contig_name: &str,
        region: ContigInterval,
    ) -> Result<SegmentMergedReads<'a>, AlignmentInputError>;
}
```

Single-file fast path (inside `query_region`): when `files.len() == 1`,
return a thin wrapper over one `MappedReadsInSegment` instead of the k-way
machinery — but still run the within-file order check (the fetcher trusts
the coordinate sort and does not).

## Reused / lifted types

- `Locus = (usize, u64)`, `ReadFingerprint`, `ReadFingerprintWithSourceFile`
  — today private in `alignment_input`; lift to `pub(crate)` (or co-locate
  with the merge). The argmin/dedup/order logic is a near-verbatim port of
  `AlignmentMergedReader`'s, minus the `classify_pre_decode`/`RecordBuf`
  conversion (already done in the fetcher).
- `FilterCounts`, `MappedRead`, `ContigInterval`, `AlignmentFile`,
  `MappedReadsInSegment`, `SegmentReadFilter` — unchanged.

## Key ideas & subtleties (the byte-identity gate)

1. **Count source moved (done), totals must not.** The cheap-filter buckets
   are tallied in the fetcher (`3e5237e`) and just summed by the merge; the
   reconciliation with the downstream `read_processor` G2/F1 tallies in
   `stage1_pipeline` is unchanged. Parity holds by construction (out-of-segment
   drops uncounted on both old and new paths) — still pin it with a
   counts-parity test.

2. **Order-check vs fetcher-filter ordering.** The old merge runs the
   within-file order check *before* dropping a too-short read; the new fetcher
   drops it *first* (and now also *counts* it as `too_short` there). So a read
   that is **both** too-short **and** a coordinate regression would error under
   the old path but be silently dropped (and counted) under the new one.
   Edge-only (normal data unaffected); narrow to exactly the *too-short* case —
   `classify_pre_decode` drops never reached the old order check either (the
   old `refill_heads` drops them before argmin), so those already agree.
   **Decision (2026-06-16): accept the divergence.** A read that is *both*
   too-short *and* a coordinate regression is dropped (and counted as
   `too_short`) rather than erroring — a corrupt-input edge case with no effect
   on well-formed sorted data, and not worth complicating the fetcher to
   reproduce the old error. The within-file order check stays (it still catches
   a regression among *kept* reads). Noted at the order-check call site in
   `segment_merge`.

3. **Dedup buffer membership.** Too-short reads never entered the old dedup
   buffer (dropped before the fingerprint push); they don't reach the new
   merge at all → same end state. Confirm with a multi-file dup test.

4. **Argmin tie-breaking** must match the old merge (lowest input index wins
   at equal locus) so emit order is identical.

5. **Pooled re-seek across contig transitions.** Regions arrive sorted; on a
   contig change the same pooled reader is re-seeked (BGZF virtual position /
   CRAM container offset). No re-open. The FASTA repo cache-clear on contig
   transition is unchanged.

## Open questions

1. Land `query_region` as a `PileupInputs` method (favoured — kills the
   9-arg signature) vs a reworked `AlignmentMergedReader::query`.
2. New module `src/bam/segment_merge.rs` vs extending `segment_reader.rs`.
3. Keep the old `query` path behind a flag for one release as the
   byte-identity oracle, or delete once the in-tree oracle test lands.
4. ~~`RecordVerdict` vs out-param~~ — **resolved** on `main` (`3e5237e`): a
   `&mut FilterCounts` out-param threaded into `classify_segment_record`.
