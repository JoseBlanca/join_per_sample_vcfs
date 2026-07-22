# ng — the alignment file: types & interfaces

*Status: architecture draft (2026-07-20), companion to the spec
[`../spec/alignment_file.md`](../spec/alignment_file.md) (the design and its rationale) and to the
shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary + step traits) and
[`module_layout.md`](module_layout.md) (the `src/ng/` tree). Sits directly under
[`read_filtering.md`](read_filtering.md) — this module builds the `RecordSource` that step 1's
`ReadFilter` consumes. Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): domain nouns for types,
verbs for functions, newtypes for domain scalars. Signatures are illustrative; the **contract** is
the deliverable. See the spec for the "why" behind every decision here.*

## Module home

`src/ng/read/input/`, a folder beside `read/filtering.rs` (spec §6). A folder rather than a file
because it holds two independent concerns with their own tests, and it shares `read/` with steps 1+2
per `module_layout.md` principle 1b — this is step 1's input edge, not a new step.

```
src/ng/read/input/
├── mod.rs          – re-exports; AlignmentFileError; check_assembly (§4)
├── open_bam.rs     – AlignmentFile::open — the validate-on-open gate (§3)
└── region_query.rs – BamRegionSource / CramRegionSource + OrderVerified (§4)
```

`sample_reads.md` owns `mod.rs`'s sample-level surface and `merge.rs`. No trait ceremony here: there
is one region reader per container format, not competing implementations of one idea, so
`BamRegionSource`/`CramRegionSource` are two `RecordSource` impls and there is no new trait.

## 1. Types

### 1.1 Reused, not re-minted

The module introduces **no scalar newtypes**, and one shared type it consumes rather than owns:
`GenomePosition` (`{ contig, position }`, §2), seeded in `ng::types` — see
[`sample_reads.md`](sample_reads.md) §1.1. Region queries speak `GenomeRegion`
([`types.rs:55`](../../../../src/ng/types.rs)) — contig + 1-based inclusive range, already the ng
region currency — and reads come back as `MappedRead`
([`alignment_input.rs:78`](../../../../src/bam/alignment_input.rs)) unchanged. Per-contig digests are
`Option<[u8; 16]>`, matching `ContigInfo::md5`
([`reference_info.rs:66`](../../../../src/ng/reference_info.rs)); ng has no `Md5` newtype today and
this module does not invent one.

### 1.2 Local types

```rust
/// One opened, **validated** alignment file: the gate of spec §3.1 has passed, so
/// `ref_id == ContigId` holds for every record it will ever yield, and both the reader
/// pool and the parsed index are owned here for the life of the run.
///
/// Not `Clone`: it owns a reader pool and an index: sharing is `&` (or `Arc`), never a copy.
pub struct AlignmentFile {
    path: PathBuf,
    header: sam::Header,
    /// Parsed once, at open — never re-read per query (spec §3.3).
    index: AlignmentIndex,
    /// From `@RG SM`; the k-file agreement check is `sample_reads`' (spec §3.1).
    sample_name: String,
    /// The `@SQ M5` tags, indexed by `ContigId`. `None` where the file carries no `M5`.
    /// Captured at open; compared only later, by `check_assembly` (§4).
    sq_md5s: Vec<Option<[u8; 16]>>,
    /// Idle readers + their scratch, borrowed per query and returned on `Drop`.
    /// A `Mutex` so `reads_in_region` can take `&self` (spec §3.3).
    readers: Mutex<Vec<ReaderHandle>>,
    filter_config: ReadFilterConfig,
}

/// One idle reader positioned past the header, **plus the scratch its next query will
/// reuse**. Pooling the buffers with the reader is what keeps a region query
/// allocation-free: a fresh `ReadFilter` per query would otherwise allocate a record
/// buffer and a reference-fetch buffer ~10⁶ times (spec §3.3).
struct ReaderHandle {
    reader: ReaderKind,           // Bam(..) | Cram(..)
    record_buf: NoodlesRawRecord, // the reused record buffer (filtering.rs:344)
    ref_buf: Vec<u8>,             // reused scratch for filter #8's reference fetch
}

/// What `check_assembly` found. Returned even when nothing was comparable, so a caller
/// can report "verified 18 of 24 contigs" rather than implying a guarantee it lacks
/// (spec §3.1). `compared == 0` is normal, not an error.
pub struct AssemblyCheck {
    pub compared: usize,
    pub total: usize,
}
```

## 2. Errors

```rust
/// Everything that can go wrong for one file. Split by *when*: the first five are
/// returned by `open` before any read flows; `OutOfOrderRead` and `Filter` arrive in the
/// item stream, after which the iterator fuses.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum AlignmentFileError {
    /// The file could not be opened, or its header could not be parsed.
    Open { path: PathBuf, source: io::Error },
    /// `@HD SO` is absent or not `coordinate` (spec §3.1 check 1).
    NotCoordinateSorted { path: PathBuf },
    /// The `@SQ` list is not the reference's contig table. `detail` is
    /// `first_disagreement`'s message, naming the first differing field and index.
    ContigReconcile { path: PathBuf, detail: String },
    /// The index is missing or unparseable (spec §3.1 check 3).
    Index { path: PathBuf, source: AlignmentIndexError },
    /// The `@RG` records name more than one sample (spec §3.1 check 4).
    MultipleSampleNames { path: PathBuf, names: Vec<String> },
    /// A record regressed in `(ref_id, pos)` — the file claims `SO:coordinate` and is not
    /// sorted. Carries both keys so the message can say where it breaks (spec §3.2).
    OutOfOrderRead { path: PathBuf, previous: GenomePosition, current: GenomePosition },
    /// Step 1 hit a fatal condition (source read, decode, or filter #8's reference fetch).
    Filter(#[source] ReadFilterError),
    /// The requested region is invalid, or names a contig absent from the reference.
    Region { region: GenomeRegion },
}

/// Raised by `check_assembly` only. A **separate** type, not a variant above: it fires
/// long after the stream is done, so it can never appear as a stream item (spec §4).
#[derive(Debug, thiserror::Error)]
pub struct AssemblyMismatch {
    pub path: PathBuf,
    pub contig: ContigId,
    pub expected: [u8; 16],
    pub observed: [u8; 16],
}
```

`GenomePosition` is `{ contig: ContigId, position: Position }` — a **shared** ng type seeded in
`ng::types` beside `GenomeRegion`, not this module's. Its derived `Ord` is genome order, which is what
makes it both the order guard's comparison here and the merge's sort key in
[`sample_reads.md`](sample_reads.md) §1.1, where it is introduced.

## 3. Opening — the gate

```rust
impl AlignmentFile {
    /// Open one indexed BAM/CRAM and **validate it or fail**: `@HD SO == coordinate`, `@SQ`
    /// equals `reference.contig_list()` exactly (order included), the index parses, and the
    /// `@RG` records name exactly one sample. No handle exists in an unvalidated state.
    pub fn open(
        path: &Path,
        reference: &ReferenceInfo,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
    ) -> Result<Self, AlignmentFileError>;

    pub fn sample_name(&self) -> &str;
    pub fn sq_md5s(&self) -> &[Option<[u8; 16]>];
}
```

**Contract.** `open` performs all I/O it will ever need for validation; after it returns `Ok`, the
invariant `ref_id == ContigId` holds for every record the file can yield, and no later code re-checks
it. The `@SQ` comparison delegates to `ContigList::first_disagreement`
([`fasta/mod.rs:69`](../../../../src/fasta/mod.rs)) against `reference.contig_list()`
([`reference_info.rs:91`](../../../../src/ng/reference_info.rs)), so the MD5 arm is automatically a
wildcard when either side lacks a digest (`ContigEntry`'s `PartialEq`,
[`fasta/mod.rs:43-55`](../../../../src/fasta/mod.rs)) — which is exactly the `.fai`-only behaviour
the spec wants, with no extra branch here.

```rust
/// Compare the `@SQ M5` tags captured at open against a **verified** `ReferenceInfo` — the
/// one carrying real per-contig digests, available once the caller joins
/// `read_reference_verifying_or_creating_fai`'s background handle. Pure, and a free
/// function rather than a method: it must run when the file handles are long gone.
///
/// A contig is compared only when **both** sides carry a digest; contigs without one are
/// skipped and counted, never an error and never a warning (spec §3.1).
pub fn check_assembly(
    path: &Path,
    observed: &[Option<[u8; 16]>],
    verified: &ReferenceInfo,
) -> Result<AssemblyCheck, AssemblyMismatch>;
```

## 4. Serving a region

Two `RecordSource` impls, then step 1's filter, then the order guard. The first seam is *not* an
iterator — see `read_filtering.md` §4 for why the source fills a reused buffer.

```rust
/// BAM: `BinningIndex::query` on the already-parsed index → chunks; seek the pooled reader
/// to each; drop records missing the interval (chunk-edge slop — a reader concern, uncounted);
/// stop once a record on the target contig starts past `end`.
struct BamRegionSource<'a> { /* borrowed handle, &index, region, chunk cursor */ }

/// CRAM: walk the `.crai` for containers on the target contig, decode each, drop
/// non-overlapping, stop once a container starts past `end`. **Carries a cursor into the
/// `.crai`** rather than rescanning from entry 0 per query (spec §3.3).
struct CramRegionSource<'a> { /* borrowed handle, &index, region, crai cursor */ }

impl RecordSource for BamRegionSource<'_> { type Record = NoodlesRawRecord; /* … */ }
impl RecordSource for CramRegionSource<'_> { type Record = NoodlesRawRecord; /* … */ }

/// Hard-errors `OutOfOrderRead` when a read's `(contig, pos)` is **strictly less** than the
/// previous one's. Equal keys are legal. State lives here, in the per-region iterator — never
/// on `AlignmentFile` — so querying region B then region A is a new forward scan, not a
/// regression (spec §3.2).
struct OrderVerified<I> { inner: I, last: Option<GenomePosition>, path: PathBuf }

/// The served stream: source → ReadFilter → order-verify, returning the pooled handle on Drop.
pub struct RegionReads<'a> { /* … */ }

impl Iterator for RegionReads<'_> {
    type Item = Result<MappedRead, AlignmentFileError>;
}
impl std::iter::FusedIterator for RegionReads<'_> {}

impl AlignmentFile {
    /// Every read overlapping `region`, coordinate-ordered and step-1-filtered.
    ///
    /// Takes **`&self`**: a handle is borrowed from the pool and returned when the returned
    /// iterator drops, so N threads may query one file concurrently (spec §3.3).
    ///
    /// `reference` is the caller's own reference accessor, passed per query — see §5.
    pub fn reads_in_region<R: RawRefSeq>(
        &self,
        region: GenomeRegion,
        reference: R,
    ) -> Result<RegionReads<'_>, AlignmentFileError>;

    /// This file's step-1 tally for the queries served so far.
    pub fn counts(&self) -> ReadFilterCounts;
}
```

**Contract.** Lazy and forward-only; one region query is one forward scan, and re-entry means a new
query. Records the index over-returns are dropped **uncounted** (reader concern, not a filter drop);
every read the filter drops is charged to exactly one `DropReason`, per `read_filtering.md` §1. A
fatal error is yielded once as `Some(Err(_))` and then `None` — the `ReadFilter` convention
([`filtering.rs:550`](../../../../src/ng/read/filtering.rs)) carried outward. The pooled handle is
returned on `Drop`, including on the error path.

## 5. Design decisions — decided

The spec argued each of these; the record carries the code shape plus the cross-ref.

- **The reference accessor is a per-query parameter, not a field — decided.** The spec left this
  open (§3.4). `RawRefSeq` impls are *stateful readers* — `WindowedRefSeq`
  ([`ref_seq.rs:500`](../../../../src/ng/ref_seq.rs)) holds an open per-contig reader — so a single
  accessor stored on `AlignmentFile` would need a mutex on the hot path the moment queries run
  concurrently, and `&self` queries would fight over one file cursor. Passing it per query gives each
  caller (and later, each worker) its own, which is the same shape production uses for its per-worker
  CRAM cache. Cost: the caller constructs one; that is a line at the call site, not a design burden.
- **No probe in the per-query `ReadFilter` — decided.** `ReadFilter::new`
  ([`filtering.rs:596`](../../../../src/ng/read/filtering.rs)) validates that every `@SQ` contig
  resolves, by fetching each one. That is correct for a whole-file pass but is **O(contigs) reference
  fetches per region query**, i.e. ~10⁶ × 24 on the STR path. The open gate (§3) already proves a
  strictly stronger property, so the probe is redundant here, exactly as spec §3.1 anticipated. Add a
  probe-free constructor to `filtering.rs` (ng-owned, so this is an extension, not a production edit)
  and use it. `ReadFilter::new` stays as-is for whole-file callers.
- **Reader pool built now, and it carries the scratch — decided.** The pool is what makes
  `reads_in_region(&self)` possible, and the signature is the part that is expensive to retrofit
  (spec §3.3). Pooling the record buffer and `ref_buf` *with* the reader is what keeps a query
  allocation-free once a fresh `ReadFilter` is built per region.
- **Two `RecordSource` impls, no new trait — decided.** BAM and CRAM are two containers for one
  idea, not competing implementations to bake off, so there is no trait ceremony
  (`module_layout.md` principle 1a applied within a folder).
- **`check_assembly` is a free function returning its own error — decided.** It runs after the
  streams are gone and cannot appear as a stream item, so a method on a dropped handle would be the
  wrong shape and an `AlignmentFileError` variant would be unreachable in the stream (spec §3.1, §4).
- **ng extracts the header fields itself — decided (corrects the spec's reuse map).** The spec §5
  listed `alignment_input::extract_single_sample_name` as reusable. It is **module-private**
  ([`alignment_input.rs:378`](../../../../src/bam/alignment_input.rs); `extract_header` likewise at
  `:292`), so ng cannot call it, and making it `pub(crate)` would be a production edit the freeze
  forbids. ng reads `@HD SO`, `@SQ` and `@RG SM` off the noodles `sam::Header` directly — a few lines,
  and it lets the errors be ng's own rather than `AlignmentInputError`'s.

## 6. Reconciliation with existing code

Every row read at the cited line.

| ng name | existing code | action |
|---|---|---|
| `GenomeRegion` (query type) | [`ng/types.rs:55`](../../../../src/ng/types.rs) | reuse as-is — contig + 1-based inclusive range |
| `ContigId`, `Position` | [`ng/types.rs:11`, `:32`](../../../../src/ng/types.rs) | reuse as-is |
| `GenomePosition` | does not exist yet | **new, in `ng::types`** — introduced by [`sample_reads.md`](sample_reads.md) §1.1; used here by the order guard |
| `MappedRead` (item) | [`bam/alignment_input.rs:78`](../../../../src/bam/alignment_input.rs) — `flag:80`, `ref_id:82` (`usize`), `pos:84` (`u64`, 1-based), `source_file_index:102` | reuse as-is |
| index load | `load_alignment_index` [`bam/index_preflight.rs:120`](../../../../src/bam/index_preflight.rs) → `AlignmentIndex` (`:95`, `Crai`/`BamCsi`/`BamBai`, each `Arc`-wrapped) | call as-is; store the result on `AlignmentFile` |
| index policy | `preflight_alignment_indexes` [`bam/index_preflight.rs:199`](../../../../src/bam/index_preflight.rs) (takes `build_if_missing`) | call as-is; also rejects mixed BAM+CRAM |
| `@SQ` comparison | `ContigList::first_disagreement` [`fasta/mod.rs:69`](../../../../src/fasta/mod.rs) — `pub(crate)`, returns `Result<(), String>` | call as-is; the `String` becomes `ContigReconcile.detail` |
| MD5 wildcard rule | `ContigEntry`'s `PartialEq` [`fasta/mod.rs:43-55`](../../../../src/fasta/mod.rs) — absent MD5 matches anything | rely on it; no extra branch |
| canonical contig table | `ReferenceInfo::contig_list()` [`ng/reference_info.rs:91`](../../../../src/ng/reference_info.rs) | call per open |
| verified digests | `ContigInfo::md5` [`ng/reference_info.rs:66`](../../../../src/ng/reference_info.rs) — `Option<[u8; 16]>` | read in `check_assembly` |
| the join point | `VerificationHandle::join` [`ng/reference_info.rs:931`](../../../../src/ng/reference_info.rs), from `read_reference_verifying_or_creating_fai` (`:1011`) | the caller's, not ours — we only supply `check_assembly` |
| `RecordSource` / `RawRecord` | [`ng/read/filtering.rs:310`, `:285`](../../../../src/ng/read/filtering.rs) | implement; do **not** redefine |
| reused record buffer | `NoodlesRawRecord` [`ng/read/filtering.rs:344`](../../../../src/ng/read/filtering.rs) | reuse as the pooled `record_buf` |
| whole-file siblings | `BamRecordSource` [`:379`](../../../../src/ng/read/filtering.rs), `CramRecordSource` [`:429`](../../../../src/ng/read/filtering.rs) | **model, not reuse** — the region impls are index-driven siblings |
| the filter | `ReadFilter` [`ng/read/filtering.rs:581`](../../../../src/ng/read/filtering.rs) + `ReadFilterConfig` (`:47`), `ReadFilterCounts` (`:117`), `ReadFilterError` (`:550`) | compose per region; **extend** with a probe-free constructor (§5) |
| reference access | `RawRefSeq` [`ng/ref_seq.rs:180`](../../../../src/ng/ref_seq.rs), `RefSeqError` (`:39`) | reuse; passed per query (§5) |
| region-read mechanics | `BamFile`/`CramFile` [`bam/segment_reader.rs:457`, `:679`](../../../../src/bam/segment_reader.rs); `BinningIndex::query` (`:485`); `borrow_handle` (`:508`) | **model, not reuse** (spec §5, O3-file → rebuild) |
| `@HD SO` / `@SQ` / `@RG SM` extraction | `extract_header` [`bam/alignment_input.rs:292`](../../../../src/bam/alignment_input.rs), `extract_single_sample_name` (`:378`) — both **module-private** | **cannot reuse**; ng reads the `sam::Header` itself (§5) |

## 7. Open items

Impl-time confirmations, not design questions.

- **`ReaderKind` shape** — whether the pooled handle holds an enum of BAM/CRAM readers or the pool is
  split per format. Production splits (`BamFile`/`CramFile`); one `AlignmentFile` with an enum is
  simpler for a caller that does not care. Confirm when the CRAM container cache lands, since spec
  §3.3 notes that cache is per-worker and non-pooled.
- **Probe-free constructor name** — `ReadFilter::with_validated_contigs` vs a `Contigs` typestate
  parameter. Prefer the plain constructor; the guarantee is documented, not encoded.
- **`counts()` under `&self`** — with a pooled, `&self` API the tally needs interior mutability
  (an atomic set, or per-handle counts summed on read). Pin when the pool is written.
- `OPEN:` **naming** — `open_bam.rs` also opens CRAM (spec §6). Kept on the owner's call; the
  alternative is `open_alignment_file.rs`.

## Test & bench shape (spec §7)

- Unit tests beside the code: the gate's four rejections (T1–T3, T12a), the order guard's four cases
  (T4a–T4d), region-query overlap/early-stop (T5), and `check_assembly`'s three outcomes (T2b) — the
  last is pure, so it needs no fixture at all.
- Fixture-driven: tiny in-memory BAM/CRAM via noodles writers plus an in-crate CSI/CRAI build, over a
  small multi-contig reference indexed by `reference_info`. BAM/CRAM parity (T8) and the full-filter
  proof (T9) are the regression anchors.
- One instrumented test for T13 (the index is parsed once across many queries) — a counting reader
  wrapper, since this is the guarantee the whole per-query cost model rests on.
- No `bench/`: no competing implementations. The merge's per-read budget is benched in
  `sample_reads.md`.
