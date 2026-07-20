# ng — the sample's reads: types & interfaces

*Status: architecture draft (2026-07-20), companion to the spec
[`../spec/sample_reads.md`](../spec/sample_reads.md) (the design and its rationale). Stands on
[`alignment_file.md`](alignment_file.md) (this doc's `AlignmentFile` and `RegionReads` come from
there) and on the shared [`ng_step_interfaces.md`](ng_step_interfaces.md) /
[`module_layout.md`](module_layout.md). Naming per
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md). Signatures are
illustrative; the **contract** is the deliverable. See the spec for every "why".*

## Module home

`src/ng/read/input/`, shared with `alignment_file.md`. This doc owns:

```
src/ng/read/input/
├── mod.rs    – SampleReads, SampleRegionReads, IngestError
└── merge.rs  – MergedRegionReads: the argmin k-way merge + the same-file-twice check
```

The `Single | Merged` enum lives in `mod.rs` with the entry point, not in `merge.rs`: the single-file
arm is the *absence* of a merge, not a variant of one (spec §3.4).

## 1. Types

### 1.1 `GenomePosition` — seeded in the shared vocabulary, not here

The merge's sort key is *a position in the genome*. That is a general value, not a merge-specific
one, so it belongs in `ng::types` beside `GenomeRegion` — **not** as a local `ReadKey`, which would
name the value after one of its uses and guarantee the next consumer mints a duplicate. The same
pair is already threaded as two parameters through a dozen `ref_seq.rs` signatures
([`ref_seq.rs:43`](../../../../src/ng/ref_seq.rs) onward), and every later consumer that says "where
along the reference" wants it.

```rust
/// One base, genome-wide: which contig, and where in it. `Position` alone is only
/// meaningful once a contig is known — this is that pair given a name, and the
/// point-shaped sibling of `GenomeRegion`.
///
/// **The derived `Ord` is genome order** — contig index first, then position — which is
/// why it serves directly as a sort key. That ordering is only meaningful because every
/// file's `ref_id` was proved equal to `ContigId` at open (`../spec/alignment_file.md`
/// §3.1); without that, contig indices would not be comparable across files at all.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct GenomePosition {
    pub contig: ContigId,
    pub position: Position,
}
```

Added to `ng_step_interfaces.md` §1; this module is its first consumer, not its owner.

### 1.2 Local types

```rust
/// One sample's files, opened and cross-validated. Built once, queried per region.
///
/// Not `Clone` (it owns k `AlignmentFile`s, each owning a reader pool).
pub struct SampleReads {
    files: Vec<AlignmentFile>,
    /// The one name all k files agree on (spec §3.1).
    sample_name: String,
}
```

No wrapper around `MappedRead`, and no per-file struct beyond `AlignmentFile`: this module adds
ordering, not data.

## 2. Errors

```rust
/// Sample-level failures. The two own variants are the cross-file checks; everything else
/// arrives from one file and is wrapped with the index of the file it came from, so a
/// message can always name the culprit.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum IngestError {
    /// The k files do not name the same sample. Fires at `open`, before any read
    /// (spec §3.1).
    SampleNameMismatch { files: Vec<PathBuf>, names: Vec<String> },
    /// The identical read surfaced from two files — the caller passed the same file twice,
    /// or two files with overlapping content (spec §3.2).
    DuplicateReadAcrossFiles { qname: Vec<u8>, key: GenomePosition, files: (usize, usize) },
    /// Anything one file raised, at open or mid-stream.
    File { source_file_index: usize, source: AlignmentFileError },
}
```

**Wrap, don't flatten** — the spec left this open (§4). Flattening would duplicate every
`AlignmentFileError` variant and lose which file raised it; wrapping keeps one source of truth and
carries the index for free.

## 3. Opening and serving

```rust
impl SampleReads {
    /// Open every file through `AlignmentFile::open` (which validates each one), then check
    /// they all name the same sample. Both happen before any read flows.
    pub fn open(
        paths: &[PathBuf],
        reference: &ReferenceInfo,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
    ) -> Result<Self, IngestError>;

    /// Every read of the sample overlapping `region`, coordinate-ordered across all files.
    ///
    /// `&self` for the same reason `AlignmentFile::reads_in_region` is (`alignment_file.md`
    /// §5): the per-file pools make concurrent queries possible without a signature change.
    /// `reference` is the caller's accessor, passed through per file.
    pub fn reads_in_region<R: RawRefSeq + Clone>(
        &self,
        region: GenomeRegion,
        reference: R,
    ) -> Result<SampleRegionReads<'_>, IngestError>;

    pub fn sample_name(&self) -> &str;
    /// Per input file, indexed to match `MappedRead::source_file_index`. Never summed —
    /// a per-experiment drop rate is the signal (spec §3.3).
    pub fn counts(&self) -> Vec<ReadFilterCounts>;
    /// Per input file; forwarded for the deferred assembly check (spec §3.4).
    pub fn sq_md5s(&self) -> Vec<&[Option<[u8; 16]>]>;
}

/// The served stream. An **enum, never `Box<dyn Iterator>`** — dynamic dispatch would block
/// inlining on the hottest loop in the module (spec §3.4).
pub enum SampleRegionReads<'a> {
    /// One file: the per-file chain, verbatim. No merge exists.
    Single(RegionReads<'a>),
    Merged(MergedRegionReads<'a>),
}

impl Iterator for SampleRegionReads<'_> {
    type Item = Result<MappedRead, IngestError>;
}
impl std::iter::FusedIterator for SampleRegionReads<'_> {}
```

**Contract.** Coordinate-ordered, lazy, fused. Every read carries the `source_file_index` of the file
it came from. The module drops nothing — its rejections are errors. It re-checks neither within-file
order nor contig identity: both are `AlignmentFile`'s guarantees (spec §2, preconditions 1 and 2).

## 4. The merge

```rust
/// Argmin k-way merge over the sample's per-file streams.
///
/// **Keys are held beside the heads, not read through them.** `keys[i]` mirrors `heads[i]`
/// and is refreshed only when that head is refilled, so the argmin scans a small contiguous
/// array of `GenomePosition`s in cache rather than dereferencing k `MappedRead`s per step. This is
/// the per-read budget of spec §3.2, expressed in the layout.
pub struct MergedRegionReads<'a> {
    streams: Vec<RegionReads<'a>>,
    /// `None` = that stream is exhausted.
    heads: Vec<Option<MappedRead>>,
    /// `keys[i]` is `heads[i]`'s key; `None` in lockstep.
    keys: Vec<Option<GenomePosition>>,
    done: bool,
}
```

**Contract, and the budget it is held to.** Each `next()`:

1. Scans `keys` for the minimum, ties to the **lowest index** — deterministic output order, which
   matters because per-experiment files tie constantly (spec §3.2).
2. **On a tie only**, checks the tied heads for the same read: `flag` (a `u16`) first, `qname` (a
   short string) only if flags also match. Two reads at different positions cannot be the same read,
   so a non-tie costs nothing — the check rides on a comparison the argmin already made. A match is
   `DuplicateReadAcrossFiles`.
3. Emits the winner with `Option::take` and refills that slot from its stream, recomputing only that
   one key.

**Per read: O(k) `GenomePosition` comparisons, one `MappedRead` move, zero allocations, and no `clone()`.**
A `clone()` in this loop is a defect, not a slow path — `MappedRead` owns its sequence, qualities and
CIGAR ([`alignment_input.rs:78`](../../../../src/bam/alignment_input.rs)). Deliberately **not** built
on `Peekable`: `peek()` hands out a `&Result<MappedRead, _>` and the emit path then moves out of the
peek slot, which is easy to write and quietly doubles the work.

Fused: the first `Err` is yielded once, then `None`.

## 5. Design decisions — decided

- **The sort key is the shared `GenomePosition`, not a local `ReadKey` — decided.** It needs to be a
  type rather than a `(ContigId, Position)` tuple because it is compared, copied and stored per read,
  and a transposed pair would silently mis-order the stream instead of failing to compile. But it must
  not be *this module's* type: "a position in the genome" is general, and naming it after one use
  (`ReadKey`) invites the next consumer to mint a duplicate. It is therefore added to `ng::types`
  beside `GenomeRegion` and to the shared vocabulary (`ng_step_interfaces.md` §1); this module is its
  first consumer, not its owner (§1.1).
- **Keys stored beside heads — decided.** The argmin must not chase pointers into `MappedRead`s; the
  key is copied once per refill (spec §3.2's budget).
- **Enum over the two arms, not `Box<dyn Iterator>` — decided.** Dynamic dispatch is opaque to the
  optimiser, so the whole per-read chain stops being inlined at the boundary — on the module's
  hottest loop (spec §3.4).
- **`IngestError` wraps `AlignmentFileError` with the file index — decided** (§2; spec §4 left the
  shape open).
- **Linear argmin, not a heap — decided.** For a few files a linear scan of a contiguous array beats
  `O(log k)` on constants and locality; production reached the same conclusion (spec §3.2).
- **The merge is concrete, not generic — decided.** No `Sortable` trait, no generic k-way merger. The
  one plausible second caller (the cohort layer) probably wants *group-by-position*, not *interleave*,
  so the shared abstraction cannot yet be identified; and the same-file-twice check does not
  generalise. If it is generalised later: a key-extractor `Fn(&T) -> K`, never a trait on the item.
  Full reasoning and the reopen trigger: **spec §5**.
- **Rebuild rather than reuse `SegmentMergedReads` — decided.** Production's merge re-checks per-file
  order (redundant here), has no merge-free single-file arm, reports counts in the wrong shape, and is
  unaudited against the budget above (spec §6, O3-merge).

## 6. Reconciliation with existing code

Every row read at the cited line.

| ng name | existing code | action |
|---|---|---|
| `AlignmentFile`, `RegionReads`, `AlignmentFileError` | [`alignment_file.md`](alignment_file.md) §1–§4 | consume; this module adds no reader logic |
| `GenomePosition` | no such type exists — `ContigId` [`ng/types.rs:11`](../../../../src/ng/types.rs) and `Position` (`:32`) are separate, and `ref_seq.rs` threads them as two parameters ([`:43`](../../../../src/ng/ref_seq.rs) onward) | **new, in `ng::types`** beside `GenomeRegion` (`:55`); built per head-refill from `MappedRead.ref_id`/`.pos` [`bam/alignment_input.rs:82-84`](../../../../src/bam/alignment_input.rs) |
| `ContigId`, `Position` | [`ng/types.rs:11`, `:32`](../../../../src/ng/types.rs) | reuse as-is |
| `MappedRead` | [`bam/alignment_input.rs:78`](../../../../src/bam/alignment_input.rs) — `qname:79`, `flag:80`, `source_file_index:102` | reuse as-is; the same-file-twice check reads `qname`/`flag` |
| `ReadFilterCounts` | [`ng/read/filtering.rs:117`](../../../../src/ng/read/filtering.rs) | reuse as-is, one per file |
| `ReferenceInfo` | [`ng/reference_info.rs:71`](../../../../src/ng/reference_info.rs) | pass through to `AlignmentFile::open` |
| `MergedRegionReads` | `SegmentMergedReads` [`bam/segment_merge.rs`](../../../../src/bam/segment_merge.rs) | **model, not reuse** (§5) |

## 7. Open items

- **`counts()` return shape** — `Vec<ReadFilterCounts>` copies k small structs per call; a slice would
  need the per-file tallies to live contiguously here rather than inside each `AlignmentFile`. Pin
  alongside `AlignmentFile::counts` (`alignment_file.md` §7), which has the same question.
- **`R: RawRefSeq + Clone`** — the k per-file chains each need an accessor. `Clone` is the simplest
  route; if the impls make cloning costly, take `&R` and have the chains share by reference, or accept
  a small factory. Confirm against `WindowedRefSeq`'s cost when coding.

## Test & bench shape (spec §7)

- Unit tests beside `merge.rs`: interleaving with correct `source_file_index` (T6), tie-break to the
  lower index plus run-to-run identity, the same-file-twice error at the **first** collision (T7), and
  same-position-different-`qname` reads both surviving — the case the cheap-first comparison order
  must not get wrong.
- `SampleReads`-level: cross-file sample-name mismatch at open (T12b), and the single-file arm
  yielding exactly what a two-file merge with one empty input yields (T11), asserted through the same
  `SampleRegionReads` type so the arms are provably indistinguishable.
- **One bench (T14)** — the merge's per-read budget over a synthetic k-file stream: zero per-read
  allocations, no `MappedRead` clone. Worth the setup because this failure is silent: a stray `clone()`
  costs throughput without changing an output byte, so no correctness test would catch it.
