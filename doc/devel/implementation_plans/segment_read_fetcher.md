# Segment read fetcher — a shared indexed-segment read source (implementation plan)

**Status:** shape sketch for alignment, 2026-06-16. A standalone, independently
testable primitive, to be wired into the **SSR** Stage-1 fetcher first and
**retrofitted into the SNP `--regions`** path later (a separate, measurable
change). Built/owned in `src/bam/` because both callers share it.

**Goal.** Given a sorted+indexed BAM/CRAM and a genomic segment, return the
reads overlapping that segment — efficiently and thread-safely — by reusing the
fact that the file is coordinate-sorted: index-seek to the segment, stream reads
until past its end. Calling it for ~10⁶ tiny SSR loci must not re-open the file
or re-parse the header per call.

Cross-refs: the analysis that motivated this (per-`query` re-open + CRAM
container re-decode), `AlignmentMergedReader`/`OwnedIndexed{Bam,Cram}Records`
(the existing per-format scanners we refactor), the SNP `--regions` tax TODO.

---

## 1. The algorithm (general terms)

For one open reader + a segment `[start, end]` on a contig:

1. **Index → candidate chunks.** The binning index (`BinningIndex::query`) +
   linear index give a set of BGZF/CRAI virtual-offset ranges that may contain
   overlapping reads — and a lower-bound offset to skip everything before the
   segment. (No per-record bisect exists or is possible — variable-length,
   block-compressed records; linear scan within index-narrowed chunks *is* the
   standard.)
2. **Seek + linear-scan + overlap-filter.** Seek to the first chunk; read
   records in coordinate order; for each: read **ends before** `start` → skip;
   read **starts after** `end` → **stop** (sort order guarantees nothing later
   overlaps); else → **yield**.

A read overlapping two segments is found independently by each segment's query
and **yielded whole to both** — the *consumer* decides what to do with the part
outside its segment (no clipping in the primitive).

---

## 2. Interface (the names settled with the PM)

```rust
// src/bam/segment_reader.rs

/// A BAM/CRAM file, immutable + Sync, shared across threads. Holds the path,
/// the Arc-shared parsed index, the header ref→ContigList map, the read-filter
/// config, and a POOL of idle open readers (so per-segment calls don't re-open).
pub(crate) enum AlignmentFile {
    Bam(BamFile),
    Cram(CramFile),
}

impl AlignmentFile {
    /// Reads overlapping `[start, end]` (1-based inclusive) on `chrom`.
    /// Borrows a pooled reader, index-seeks it to the segment, and returns an
    /// iterator that streams until past the segment end. Thread-safe: callable
    /// concurrently from many threads for different segments.
    pub(crate) fn get_reads_from_segment(
        &self, chrom: &str, start: u32, end: u32,
    ) -> Result<MappedReadsInSegment<'_>, AlignmentInputError>;

    pub(crate) fn from_input(/* path + Arc index + header map + cfg */) -> Result<Self, _>;
}

/// The per-call iterator. `(chrom,start,end)` are fixed at creation; the only
/// mutable state is the reader cursor. Borrows `&'a AlignmentFile` (no 'static).
/// On Drop it returns its borrowed reader to the pool.
pub(crate) enum MappedReadsInSegment<'a> {
    Bam(BamSegmentReads<'a>),
    Cram(CramSegmentReads<'a>),
}
impl Iterator for MappedReadsInSegment<'_> {
    type Item = Result<MappedRead, AlignmentInputError>;
}
```

`BamFile` / `CramFile` each:
```rust
pub(crate) struct BamFile {
    path: PathBuf,
    index: Arc<BamIndex>,                 // BAI/CSI, parsed once (shared)
    ref_map: Arc<HeaderRefMap>,           // header ref-id ↔ name ↔ ContigList id
    cfg: ReadFilterConfig,                // mapq/dup/flag pre-decode filter
    readers_pool: Mutex<Vec<BamReaderHandle>>, // idle open readers (the pool)
}
struct BamReaderHandle { reader: bam::io::Reader<bgzf::io::Reader<File>> } // seekable
```

`get_reads_from_segment` returns a `Result` because borrowing (pop-or-open) +
the initial index-seek can fail; iteration errors surface as `Some(Err(_))`.

---

## 3. The reader pool (in v1, not deferred)

- `BamFile`/`CramFile` hold `Mutex<Vec<Handle>>`.
- **Borrow** (`get_reads_from_segment`): lock, `pop()` an idle handle, or — if
  empty — **open a fresh reader** (the only file-open; happens ~once per thread
  on warmup); unlock; seek the handle to the segment's first chunk.
- **Return** (`MappedReadsInSegment::drop`): lock, `push()` the handle back.
- The lock is held only for the pop/push, **never during iteration**, so
  contention is near-zero (≈ one borrow per thread at a time, per the PM:
  "only one segment per thread usually"). The pool naturally grows to the
  concurrency level and is reused thereafter.
- Thread-safety: `AlignmentFile: Sync` (Mutex pool + Arc index/map + immutable
  cfg); shared by `&` across rayon workers; each `get_reads_from_segment` is
  independent. No shared cursor, no lock held across work.

This is the project's established "per-worker owned reader, brief uncontended
lock" pattern (cf. the ref-fetcher review preferring it over a held `Mutex`).

---

## 4. Internals per format (reuse existing scanners)

Most of this is **refactoring** the existing per-query owned scanners into
pool-backed, re-seekable handles:

- **BAM** (`OwnedIndexedBamRecords` today): `BinningIndex::query(ref_id,
  interval) → Vec<Chunk>`; seek the handle to `chunk.start()`; read records
  until `virtual_position() >= chunk.end`; advance chunks; per record:
  `ContigInterval::overlaps_record` filter + the flag/mapq pre-decode filter
  (`classify_pre_decode`); **early-stop** once `alignment_start > end`
  (optimization; the chunk set already bounds it). Records →
  `record_buf_to_mapped_read` (existing) → `MappedRead`.
- **CRAM** (`OwnedIndexedCramRecords` today): `.crai` container walk; seek to a
  container's `offset()`; decode its slices; filter records by overlap.
  **No container cache in v1** (PM: thousands of records is small vs the genome;
  revisit only if profiling says so) — but the reader IS pooled (no re-open).
- Both reuse: the `MappedRead` type, the canonical-`ContigList` ref mapping, the
  pre-decode flag/mapq filter, and `ContigInterval` (1-based inclusive) +
  `overlaps_record`.

**noodles dependency (confirmed public API — no reimplementation).** This
primitive is entirely noodles-backed; we are *not* building a noodles-independent
reader. Every method the seek + linear-scan + (CRAM) container-decode flow needs
is part of noodles' **public** surface, so the work is thin ownership/pooling
wrappers over noodles, not a fork of its decoders:
- BAM: `csi::BinningIndex::query` → `Chunk::{start,end}`; `bgzf::io::Reader::{seek,
  virtual_position}`; `bam::io::Reader::read_record_buf`.
- CRAM (the complex path — decode logic stays *inside* noodles): `crai::Record::
  {reference_sequence_id, alignment_start, alignment_span, offset}`;
  `cram::io::Reader::{seek, read_container}`; `Container::{compression_header,
  slices}` → `Slice::{decode_blocks, records}` → `RecordBuf::
  try_from_alignment_record`. All verified `pub` in noodles-cram 0.93.
The existing `OwnedIndexed{Bam,Cram}Records` already mirror noodles' own `Query`
iterators; we only add reader ownership + the pool.

**No multi-input merge here** (PM): the primitive is *one file → its reads for a
segment*. If a sample has several files, the driver queries each and combines
downstream — and SSR needs no ordered merge at all (a locus's reads go into the
order-independent reservoir), so it just concatenates per file.

---

## 5. Coordinate conventions (pin with fixtures)

- Public API: `(chrom, start, end)` **1-based inclusive** (matches
  `ContigInterval`). The SSR driver converts the catalog `Locus` (0-based
  half-open) + `flank_bp` to this; the SNP driver passes its `ContigInterval`
  directly.
- Internally → noodles `Region`/`Interval` for the index query (the existing
  `query_interval` conversion).

---

## 6. Increment sequence (each its own commit + tests)

1. **BAM path** — `BamFile` + pool + `BamSegmentReads` + `get_reads_from_segment`
   (refactor `OwnedIndexedBamRecords` into the pooled, re-seekable shape).
2. **CRAM path** — `CramFile` + pool + `CramSegmentReads` (refactor
   `OwnedIndexedCramRecords`).
3. **`AlignmentFile` enum** + `from_input` construction (build from the existing
   `PileupInputs` per-input handles: path, `Arc` index, header map, cfg).
4. **Wire the SSR fetcher** onto it (`fetch_locus_reads` → one
   `get_reads_from_segment` per locus into the reservoir).
5. **(Separate, later) SNP `--regions` retrofit** — replace the per-BED-interval
   re-`query` with this primitive; measure the `--regions` tax improvement.

Build #1–#4 now; #5 is its own measured change.

---

## 7. Tests (the safety net)

Synthetic indexed BAM (via `tests/common` BAM builders) with reads at known
coordinates + a `.bai`/`.csi`:
- **Overlap correctness** — a segment returns exactly the reads overlapping it;
  reads strictly before/after are excluded; reads touching each boundary are
  handled (pin the 1-based-inclusive edges).
- **Read spanning two segments** — a long read straddling two queried segments
  is yielded by **both** calls (whole, unclipped).
- **Pool reuse** — N sequential `get_reads_from_segment` calls open the file
  **once** (assert via a counter / handle identity), and the pool returns to its
  resting size after each iterator drops.
- **Thread-safety** — `par_iter` over many segments on a shared `&AlignmentFile`
  returns the same reads as the sequential run (determinism), with no data race
  (the test is the proof the `Sync` design holds).
- **Empty segment** — no overlapping reads → empty iterator, handle returned.
- **CRAM** — the same overlap + pool tests against a synthetic indexed CRAM
  (reuse the cram_input test fixtures).

---

## 8. Deferred (explicitly not v1)

- **CRAM container cache** (decode-once across adjacent same-container segments).
- **Multi-input merge** (downstream / not needed for SSR).
- **SNP `--regions` retrofit** (increment #5, separate measured change).

---

## 9. Open items to pin while coding

1. Exact `BamReaderHandle` type + how noodles' bam reader is seeked to a
   `bgzf::VirtualPosition` and read record-by-record while owning the reader
   (the existing `OwnedIndexedBamRecords` already does this — lift it).
2. The `HeaderRefMap` shape (chrom name → index ref-id → canonical ContigList
   id) and where it's built (once, in `from_input`, from `PileupInputs`).
3. Whether `get_reads_from_segment` borrows the handle eagerly (at call) or
   lazily (first `next`). Lean **eager** (simpler; the SSR driver drains each
   iterator before the next).
4. Error surface: open/seek failure as the `Result` of `get_reads_from_segment`
   vs the first `next()`. Lean **`Result` on the call** for open/seek; per-record
   decode errors as `Some(Err(_))`.
