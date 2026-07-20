//! `AlignmentFile::open` — the validate-on-open gate, and the handle it
//! produces.
//!
//! Every check lives inside the function that opens the file, so a file is
//! either opened *and* validated or it is an `Err`: there is no window in which
//! an unvalidated handle exists. The checks run fail-fast in this order —
//! `@HD SO`, `@SQ`↔reference, the index, `@RG SM`
//! (`doc/devel/ng/spec/alignment_file.md` §3.1).
//!
//! **The invariant this establishes** is what every later layer leans on: once
//! the gate passes, `ref_id == ContigId` holds by construction for every record
//! the file can yield. That is what makes it sound for the merge one layer up
//! to compare positions across files without re-checking anything.
//!
//! ng reads `@HD SO`, `@SQ` and `@RG SM` off the noodles `sam::Header` itself
//! rather than reusing production's extractors: those are module-private, and
//! making them visible would be a production edit the ng freeze forbids
//! (`doc/devel/ng/arch/alignment_file.md` §5).
//!
//! The name says BAM but the module opens CRAM too — "BAM" in the everyday
//! sense of "the alignment file" (spec §6).

use std::fs::File;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use crate::bam::alignment_input::{MappedRead, build_fasta_repository};
use crate::bam::index_preflight::{
    AlignmentFileKind, AlignmentIndex, load_alignment_index, preflight_alignment_indexes,
};
use crate::fasta::{ContigEntry, ContigList};
use crate::ng::read::filtering::{
    NoodlesRawRecord, ReadFilter, ReadFilterBuffers, ReadFilterConfig, ReadFilterCounts,
};
use crate::ng::ref_seq::RawRefSeq;
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::types::GenomeRegion;

use super::AlignmentFileError;
use super::region_query::{
    BamRegionSource, CramRegionPlan, CramRegionSource, OrderVerified, RegionPlan, RegionSource,
};

/// Which planner ran, and its result. Paired with the pooled reader's own
/// format immediately below; both come from the same path, so they always
/// agree.
enum QueryPlan {
    Bam(RegionPlan),
    Cram(CramRegionPlan),
}

/// One opened, **validated** alignment file.
///
/// The gate of [`AlignmentFile::open`] has passed, so `ref_id == ContigId`
/// holds for every record this file will ever yield — the invariant the region
/// query, the order guard and the cross-file merge all lean on without
/// re-checking. There is no way to hold one of these unvalidated: `open` either
/// returns a validated handle or an error.
///
/// The parsed index is owned here and lives for the whole run: a region query
/// is an in-memory lookup plus a seek, never a file open or an index parse
/// (spec §3.3). Not `Clone` — it owns an index and a reader pool; sharing is
/// by reference.
pub struct AlignmentFile {
    /// `Arc` so the per-query order guard can hold it for its error message
    /// without an allocation per query.
    path: Arc<Path>,
    /// Kept for the region query, which resolves a contig to a `ref_id` and
    /// hands the header to the record source. Read from C2 onwards.
    header: sam::Header,
    /// Parsed once, at open — never re-read per query (spec §3.3). Queried from
    /// C2 onwards, which is the guarantee the whole per-query cost model rests
    /// on: a query is an in-memory lookup plus a seek.
    index: AlignmentIndex,
    /// From `@RG SM`. The k-file agreement check belongs to `SampleReads`.
    sample_name: String,
    /// The `@SQ M5` tags, indexed by `ContigId` — which is sound precisely
    /// because the gate just proved this file's `@SQ` order *is* the
    /// reference's. `None` where the file carries no usable `M5`.
    ///
    /// Captured at open and compared much later, by `check_assembly` (D1),
    /// once the caller has joined `reference_info`'s background verification.
    sq_md5s: Vec<Option<[u8; 16]>>,
    /// Handed to the per-query `ReadFilter` from C4 onwards. Held on the file
    /// rather than passed per query because the filtering policy is the file's
    /// for the whole run, not the caller's per region.
    filter_config: ReadFilterConfig,
    /// Idle readers and their scratch, borrowed per query and returned on
    /// `Drop`. A `Mutex` so `reads_in_region` can take `&self` (spec §3.3).
    ///
    /// The lock only ever guards a `Vec` pop or push — never a read, and never
    /// a file open (see `borrow_reader`) — so it is held for a few instructions
    /// and never across iteration.
    readers: Mutex<Vec<ReaderHandle>>,
    /// How many readers have actually been opened, ever.
    ///
    /// The pool's whole purpose is that this stays tiny however many queries
    /// run, so it is worth being able to state rather than assume — T13 asserts
    /// it directly.
    readers_opened: AtomicUsize,
    /// This CRAM's `.crai` entries, grouped by contig — `crai_by_contig[i]`
    /// holds contig `i`'s entries in file order. Empty for a BAM.
    crai_by_contig: Vec<Arc<[cram::crai::Record]>>,
    /// The reference bases CRAM slice decoding consults, built once at open.
    ///
    /// `None` for a BAM, which stores its own sequences and needs no
    /// reference to decode.
    reference_repository: Option<fasta::Repository>,
    /// This file's step-1 tally, summed as each region stream ends.
    ///
    /// The `&self` API's answer to `counts()`: a per-read running total would
    /// need a lock on the hot path or an atomic per counter, so each stream
    /// keeps its own tally and folds it in here on `Drop` — one lock per query,
    /// beside the one the pool already takes.
    counts: Mutex<ReadFilterCounts>,
    /// Stamped onto every `MappedRead` this file yields.
    source_file_index: usize,
}

/// One idle reader positioned past the header, **plus the scratch its next
/// query will reuse**.
///
/// Pooling the buffers with the reader is what keeps a region query
/// allocation-free: a fresh `ReadFilter` per query would otherwise allocate a
/// record buffer and a reference-fetch buffer ~10⁶ times (spec §3.3). The index
/// and the header are deliberately *not* here — they stay on `AlignmentFile`,
/// shared, so no pooled reader ever re-parses them.
struct ReaderHandle {
    reader: ReaderKind,
    buffers: ReadFilterBuffers<NoodlesRawRecord>,
}

/// A reader for whichever container this file is.
///
/// One pool holding an enum, rather than production's split `BamFile` /
/// `CramFile` (`segment_reader.rs`): ng has a single `AlignmentFile` for both
/// formats because a caller asking for a region does not care which it is, so
/// splitting the pool would mean splitting the type above it too. Arch §7 left
/// this open and flagged it for revisiting when the CRAM container cache lands
/// — that cache is per-worker and *not* pooled, so it will sit beside this
/// rather than inside it.
pub(crate) enum ReaderKind {
    Bam(bam::io::Reader<bgzf::io::Reader<File>>),
    Cram(cram::io::Reader<File>),
}

/// A reader borrowed from the pool, returned when this drops.
///
/// The handle is an `Option` purely so `Drop` can move it out — `drop` gets
/// `&mut self` and cannot take ownership otherwise. Returning on `Drop` rather
/// than at the end of a successful read is what makes the error and early-exit
/// paths safe: a caller that abandons a region stream half-way still gives the
/// reader back.
struct BorrowedReader<'a> {
    file: &'a AlignmentFile,
    handle: Option<ReaderHandle>,
}

impl BorrowedReader<'_> {
    /// Take the handle's parts, leaving the borrow to return an empty slot.
    ///
    /// Used when the parts are moved into a region stream that will hand them
    /// back itself; the `Drop` below then has nothing to return.
    ///
    /// **Must be the last, infallible step.** Taking the handle severs it from
    /// this borrow, so the obligation to return it transfers to whatever
    /// receives it. A `?` between here and the point that obligation is
    /// re-established loses a reader silently — the pool would just open
    /// another, and only `readers_opened` climbing would show it. Do every
    /// fallible part of building a region stream *before* calling this.
    fn take(mut self) -> ReaderHandle {
        self.handle
            .take()
            .expect("a borrowed reader always holds its handle until taken")
    }
}

impl Drop for BorrowedReader<'_> {
    fn drop(&mut self) {
        if let Some(handle) = self.handle.take() {
            self.file.return_handle(handle);
        }
    }
}

/// One region's reads: source → step-1 filter → order guard, coordinate-ordered
/// and lazy.
///
/// Forward-only. A region query is one forward scan, and re-entry means a new
/// query. A fatal error is yielded once and then the iterator is done — the
/// `ReadFilter` convention carried outward, so `let read = read?;` surfaces it
/// and it cannot be mistaken for a clean end of input.
///
/// **The pooled reader goes back on `Drop`, including on the error path and
/// when a caller abandons the stream half-way**, along with this query's drop
/// tally. An `Option` because `drop` gets `&mut self` and has to move the
/// stream out to unwrap it.
pub struct RegionReads<'a, R: RawRefSeq> {
    file: &'a AlignmentFile,
    stream: Option<OrderVerified<ReadFilter<RegionSource<'a>, R>>>,
}

impl<R: RawRefSeq> Iterator for RegionReads<'_, R> {
    type Item = Result<MappedRead, AlignmentFileError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.stream.as_mut()?.next()
    }
}

impl<R: RawRefSeq> std::iter::FusedIterator for RegionReads<'_, R> {}

impl<R: RawRefSeq> Drop for RegionReads<'_, R> {
    fn drop(&mut self) {
        let Some(stream) = self.stream.take() else {
            return;
        };

        // Unwrap the chain in the order it was built: guard → filter → source
        // → reader. `into_parts` also yields the buffers and the tally, which
        // is why the filter has to be taken apart rather than just dropped —
        // the drops this query recorded would vanish with it.
        let (source, buffers, counts) = stream.into_inner().into_parts();
        self.file.return_handle(ReaderHandle {
            reader: source.into_reader(),
            buffers,
        });
        self.file.add_counts(&counts);
    }
}

impl AlignmentFile {
    /// Open one indexed BAM/CRAM and **validate it or fail**.
    ///
    /// Four checks, fail-fast **in this order** (spec §3.1): `@HD SO` is
    /// `coordinate`; the `@SQ` list equals `reference.contig_list()` exactly,
    /// order included; the index loads; the `@RG` records name exactly one
    /// sample. The order is cheapest-and-most-fundamental first, so a file that
    /// is wrong in several ways reports the most basic fault rather than
    /// whichever check happens to run.
    ///
    /// A CRAM then needs one more thing before it can be read at all — the
    /// reference *bases*, to decode against — so the repository is built after
    /// those four, and a reference that cannot supply one is rejected here
    /// rather than at the first query.
    ///
    /// **The `@SQ` check is the permutation fix.** Comparing *in order* against
    /// the reference is what distinguishes this from the resolves-only check it
    /// replaces: a file whose contig list is a re-ordering of the reference's
    /// resolves every index and then fetches the wrong contig for every read.
    ///
    /// With `build_index_if_missing`, an absent index is built next to the file
    /// rather than rejected; that is a caller policy, not this module's.
    pub fn open(
        path: &Path,
        reference: &ReferenceInfo,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
    ) -> Result<Self, AlignmentFileError> {
        Self::open_as(path, reference, filter_config, build_index_if_missing, 0)
    }

    /// `open`, with the file's index within its sample — the value stamped onto
    /// every `MappedRead` as `source_file_index`. `SampleReads` (E1) passes the
    /// position in its list; a lone file is 0.
    pub fn open_as(
        path: &Path,
        reference: &ReferenceInfo,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
        source_file_index: usize,
    ) -> Result<Self, AlignmentFileError> {
        let header = read_header(path)?;

        // 1. @HD SO.
        let observed_sort_order = sort_order(&header);
        if observed_sort_order.as_deref() != Some("coordinate") {
            return Err(AlignmentFileError::NotCoordinateSorted {
                path: path.to_path_buf(),
                sort_order: observed_sort_order,
            });
        }

        // 2. @SQ vs the reference's contig table — names, lengths, order, and
        //    digests where both sides carry one. `first_disagreement` applies
        //    the absent-digest-is-a-wildcard rule itself (`fasta/mod.rs`, its
        //    `if let (Some, Some)` arm) rather than going through
        //    `ContigEntry`'s `PartialEq`, which encodes the same rule
        //    separately — so a `.fai`-only reference compares on name and
        //    length alone with no branch needed here.
        //
        //    The file is the left operand so the message reads "file value vs
        //    reference value", the direction a user needs: the reference is the
        //    authority and the file is the thing that is wrong.
        let file_contigs =
            contig_list(&header).map_err(|bad| AlignmentFileError::MalformedMd5 {
                path: path.to_path_buf(),
                contig: bad.contig,
                detail: bad.detail,
            })?;
        file_contigs
            .first_disagreement(&reference.contig_list())
            .map_err(|detail| AlignmentFileError::ContigReconcile {
                path: path.to_path_buf(),
                detail,
            })?;

        // 3. The index — parsed here, once, and held for the life of the run.
        //
        // Building a missing index reaches `noodles_cram::fs::index`, which
        // decodes *multi-reference* slices against an empty
        // `fasta::Repository` (marked `// TODO` in noodles 0.93) and therefore
        // **panics** on a CRAM whose reads are stored as differences from the
        // reference and whose slices span contigs — plausible on a fragmented
        // assembly. Building the index outside ng avoids it. Noted here
        // because this is the only ng call that can reach that path.
        if build_index_if_missing {
            preflight_alignment_indexes(&[path.to_path_buf()], true).map_err(|source| {
                AlignmentFileError::Index {
                    path: path.to_path_buf(),
                    source,
                }
            })?;
        }
        let index = load_alignment_index(path).map_err(|source| AlignmentFileError::Index {
            path: path.to_path_buf(),
            source,
        })?;

        // 4. Exactly one @RG SM.
        let sample_name = match sample_names(&header) {
            SampleNames::One(name) => name,
            SampleNames::Several(names) => {
                return Err(AlignmentFileError::MultipleSampleNames {
                    path: path.to_path_buf(),
                    names,
                });
            }
            SampleNames::MissingTag { read_group_id } => {
                return Err(AlignmentFileError::MissingSampleName {
                    path: path.to_path_buf(),
                    read_group: Some(read_group_id),
                });
            }
            SampleNames::NoReadGroups => {
                return Err(AlignmentFileError::MissingSampleName {
                    path: path.to_path_buf(),
                    read_group: None,
                });
            }
        };

        // 5. A CRAM needs the reference *bases* to decode at all, so the
        //    repository is built here — once — and a reference that cannot
        //    supply one is a hard error now rather than a mystery at the first
        //    query.
        let crai_by_contig = match &index {
            AlignmentIndex::Crai(crai) => group_crai_by_contig(crai, file_contigs.entries.len()),
            _ => Vec::new(),
        };

        let reference_repository = match AlignmentFileKind::from_path(path) {
            Some(AlignmentFileKind::Cram) => {
                let fasta = reference.fasta_path.as_deref().ok_or_else(|| {
                    AlignmentFileError::CramNeedsReferenceFasta {
                        path: path.to_path_buf(),
                    }
                })?;
                Some(
                    build_fasta_repository(fasta).map_err(|source| AlignmentFileError::Open {
                        path: fasta.to_path_buf(),
                        source: std::io::Error::other(source),
                    })?,
                )
            }
            _ => None,
        };

        // Free, and only sound now: check 2 proved this order is the
        // reference's, so position i really is `ContigId(i)`.
        let sq_md5s = file_contigs.entries.iter().map(|entry| entry.md5).collect();

        Ok(Self {
            path: Arc::from(path),
            header,
            index,
            sample_name,
            sq_md5s,
            filter_config,
            // Empty: the first query opens the first reader. `open` itself
            // never leaves one behind, so a file that is opened and never
            // queried costs no descriptor.
            readers: Mutex::new(Vec::new()),
            readers_opened: AtomicUsize::new(0),
            crai_by_contig,
            reference_repository,
            counts: Mutex::new(ReadFilterCounts::default()),
            source_file_index,
        })
    }

    /// The single sample this file's `@RG` records name.
    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }

    /// The `@SQ M5` tags, indexed by `ContigId`, for the deferred assembly
    /// check. `None` where the file carried no usable digest.
    pub fn sq_md5s(&self) -> &[Option<[u8; 16]>] {
        &self.sq_md5s
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Every read overlapping `region`, coordinate-ordered and step-1-filtered.
    ///
    /// The chain is source → [`ReadFilter`] → order guard, and it is the
    /// *complete* product of this module: the sample layer either uses it as-is
    /// or merges k of them, and this module is indifferent to which.
    ///
    /// Takes **`&self`**, not `&mut self` — the load-bearing signature choice
    /// (spec §3.3). A reader comes from the internal pool and the returned
    /// iterator gives it back on `Drop`, so N threads may query one file
    /// concurrently without touching a single call site.
    ///
    /// `reference` is the caller's own accessor, passed **per query** rather
    /// than stored: `RawRefSeq` impls are stateful readers, so one shared
    /// accessor would need a mutex on the hot path the moment queries run
    /// concurrently, and `&self` queries would fight over one file cursor
    /// (arch §5).
    ///
    /// The filter is built **probe-free**: the open gate already proved
    /// something strictly stronger than the per-contig resolve probe, and
    /// paying that probe per query would mean ~10⁶ × the contig count in
    /// reference fetches.
    pub fn reads_in_region<R: RawRefSeq>(
        &self,
        region: GenomeRegion,
        reference: R,
    ) -> Result<RegionReads<'_, R>, AlignmentFileError> {
        // **Everything fallible happens first**, before a reader leaves the
        // pool: the format check, resolving the region, and querying the index
        // all touch no reader, so any of them failing costs nothing to recover
        // from — whereas a `?` after the handle is taken would lose that reader
        // silently (`BorrowedReader::take`).
        //
        // The format decides which planner runs, and it has to be decided
        // first: a CRAM's index is a `.crai`, which the BAM planner cannot read.
        let plan = match AlignmentFileKind::from_path(&self.path) {
            Some(AlignmentFileKind::Cram) => QueryPlan::Cram(CramRegionSource::plan(
                &self.header,
                &self.crai_by_contig,
                region,
            )?),
            _ => QueryPlan::Bam(BamRegionSource::plan(
                &self.header,
                &self.index,
                region,
                &self.path,
            )?),
        };

        let borrowed = self.borrow_reader()?;
        let handle = borrowed.take();
        let ReaderHandle { reader, buffers } = handle;

        let source = match (reader, plan) {
            (ReaderKind::Bam(reader), QueryPlan::Bam(plan)) => RegionSource::Bam(
                BamRegionSource::new(reader, &self.header, plan, self.source_file_index),
            ),
            (ReaderKind::Cram(reader), QueryPlan::Cram(plan)) => {
                RegionSource::Cram(CramRegionSource::new(
                    reader,
                    &self.header,
                    self.reference_repository
                        .clone()
                        .expect("open builds a repository for every CRAM"),
                    plan,
                    self.source_file_index,
                ))
            }
            // Both the reader and the plan are chosen from the same path's
            // extension, so they cannot disagree.
            _ => unreachable!("the reader's format and the plan's are both the path's"),
        };

        let filter =
            ReadFilter::with_validated_contigs(source, reference, self.filter_config, buffers);

        Ok(RegionReads {
            file: self,
            stream: Some(OrderVerified::new(filter, Arc::clone(&self.path))),
        })
    }

    /// This file's step-1 tally, summed over the region queries served so far.
    ///
    /// **Per query, not per read**: each stream adds its own tally here when it
    /// ends, so a stream still running is not yet reflected. That is the price
    /// of `&self` — a running total updated per read would need either a lock
    /// on the hot path or an atomic per counter, and neither buys anything a
    /// caller reading counts *after* draining actually wants.
    ///
    /// A stream that is never dropped — leaked, or held for the process
    /// lifetime — never contributes at all, and takes its pooled reader with
    /// it. That is inherent to returning things on `Drop`.
    pub fn counts(&self) -> ReadFilterCounts {
        self.counts
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .clone()
    }

    /// Take an idle reader from the pool, opening one only if the pool is
    /// empty.
    ///
    /// This is the only place the file is ever opened after `AlignmentFile::open`,
    /// so in a single-threaded run it happens exactly once however many regions
    /// are queried — the guarantee the whole per-query cost model rests on
    /// (spec §3.3). With N threads it settles at N.
    ///
    /// A freshly opened reader is positioned past the header. A *returned* CRAM
    /// reader sits wherever its last query stopped, which is harmless only
    /// because `CramRegionSource` seeks to a container before every read and
    /// never reads without seeking — the BAM source does the same per chunk.
    /// Neither source may assume an incoming reader's position.
    fn borrow_reader(&self) -> Result<BorrowedReader<'_>, AlignmentFileError> {
        // The `let` ends the guard's scope, so the lock is released before
        // `open_reader` runs. Written as a `match` scrutinee it would still be
        // held inside the `None` arm — temporaries live to the end of the whole
        // `match` — and every other thread would block for the length of a file
        // open plus a header read, serializing at exactly the point the pool
        // exists to keep parallel. Production takes the same care
        // (`segment_reader.rs` `borrow_handle`, via an early return).
        let pooled = self.lock_pool().pop();
        let handle = match pooled {
            Some(handle) => handle,
            None => self.open_reader()?,
        };
        Ok(BorrowedReader {
            file: self,
            handle: Some(handle),
        })
    }

    /// Fold a finished stream's tally into the file's.
    fn add_counts(&self, counts: &ReadFilterCounts) {
        let mut total = self
            .counts
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        total.add(counts);
    }

    fn return_handle(&self, handle: ReaderHandle) {
        self.lock_pool().push(handle);
    }

    /// Open a reader and skip its header, leaving it where a returned one would
    /// be. The header is read and discarded: this file already holds the parsed
    /// one, and the reader has to be advanced past it either way.
    fn open_reader(&self) -> Result<ReaderHandle, AlignmentFileError> {
        let open_error = |source: std::io::Error| AlignmentFileError::Open {
            path: self.path.to_path_buf(),
            source,
        };

        let reader = match AlignmentFileKind::from_path(&self.path) {
            Some(AlignmentFileKind::Bam) => {
                let mut reader = bam::io::reader::Builder
                    .build_from_path(&self.path)
                    .map_err(open_error)?;
                reader.read_header().map_err(open_error)?;
                ReaderKind::Bam(reader)
            }
            Some(AlignmentFileKind::Cram) => {
                let file = File::open(&self.path).map_err(open_error)?;
                let mut reader = cram::io::Reader::new(file);
                reader.read_header().map_err(open_error)?;
                ReaderKind::Cram(reader)
            }
            // `open` already rejected any other extension, so this is
            // reachable only if the path changed kind underneath us. Reported
            // the same way `read_header` reports it, rather than as a CRAM
            // magic-number failure.
            None => {
                return Err(open_error(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    "expected a '.bam' or '.cram' extension",
                )));
            }
        };

        self.readers_opened.fetch_add(1, Ordering::Relaxed);
        Ok(ReaderHandle {
            reader,
            buffers: ReadFilterBuffers::default(),
        })
    }

    /// Lock the pool, recovering the guard if a panic elsewhere poisoned the
    /// mutex.
    ///
    /// The lock only ever guards a `Vec` pop or push, neither of which can
    /// unwind, so this type's own code cannot poison it. Recovering rather than
    /// panicking keeps a poison introduced elsewhere from cascading into every
    /// later query — and, on the return path, from silently dropping the reader.
    fn lock_pool(&self) -> std::sync::MutexGuard<'_, Vec<ReaderHandle>> {
        self.readers
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    /// How many readers this file has opened, ever — one per concurrent caller,
    /// not one per query.
    fn readers_opened(&self) -> usize {
        self.readers_opened.load(Ordering::Relaxed)
    }

    #[cfg(test)]
    fn pooled_readers(&self) -> usize {
        self.lock_pool().len()
    }
}

/// Hand-written so the reader pool does not have to be `Debug` (noodles'
/// readers are not), and so the output says what identifies the file rather
/// than dumping a parsed index.
impl std::fmt::Debug for AlignmentFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AlignmentFile")
            .field("path", &self.path)
            .field("sample_name", &self.sample_name)
            .field("contigs", &self.sq_md5s.len())
            .field("readers_opened", &self.readers_opened())
            .finish_non_exhaustive()
    }
}

/// Bucket a `.crai`'s entries by contig, keeping each contig's entries in file
/// order. One O(n) pass, at open.
///
/// **This replaces a binary search, and the reason is correctness rather than
/// speed.** Searching assumes the `.crai` is ordered by contig with unplaced
/// entries last — and that assumption is false for an index noodles itself
/// writes: within a multi-reference slice it emits entries sorted by
/// `Option<usize>`, and `None < Some(0)` in Rust, so the unplaced entry comes
/// *first*. A `partition_point` over an unpartitioned slice returns an
/// unspecified index, and the walk would then find a foreign entry and report
/// end-of-input — **losing every read of the region, with no error**.
///
/// Grouping assumes nothing about the order between contigs. It also makes the
/// per-query lookup O(1) rather than O(log n), which matters at one query per
/// STR locus. Within a contig the file order is kept, because the
/// container-level early stop does rely on containers ascending by start —
/// which follows from the file being coordinate-sorted, the same assumption the
/// BAM early stop makes.
///
/// Unplaced entries are dropped: they can never overlap a region.
pub(crate) fn group_crai_by_contig(
    index: &cram::crai::Index,
    contig_count: usize,
) -> Vec<Arc<[cram::crai::Record]>> {
    let mut by_contig: Vec<Vec<cram::crai::Record>> = vec![Vec::new(); contig_count];

    for record in index.iter() {
        if let Some(contig) = record.reference_sequence_id()
            && contig < contig_count
        {
            by_contig[contig].push(record.clone());
        }
    }

    by_contig.into_iter().map(Arc::from).collect()
}

/// Read just the SAM header, dispatching on the file's extension.
///
/// The reader is opened and dropped here: the gate needs the header, and the
/// readers a query will use come from the pool (step C1).
fn read_header(path: &Path) -> Result<sam::Header, AlignmentFileError> {
    let open_error = |source: std::io::Error| AlignmentFileError::Open {
        path: path.to_path_buf(),
        source,
    };

    let file = File::open(path).map_err(open_error)?;

    match AlignmentFileKind::from_path(path) {
        Some(AlignmentFileKind::Bam) => bam::io::Reader::new(file).read_header(),
        Some(AlignmentFileKind::Cram) => cram::io::Reader::new(file).read_header(),
        // `load_alignment_index` rejects this too, but the header read comes
        // first, so the extension has to be understood here.
        None => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "expected a '.bam' or '.cram' extension",
        )),
    }
    .map_err(open_error)
}

/// The `@HD SO` value, or `None` when the header carries no `@HD` line or no
/// sort-order tag at all.
///
/// The value is returned rather than compared here so the caller owns the
/// policy: "missing" and "queryname" are the same rejection, but they are not
/// the same *diagnosis*, and only the caller knows whether its error can say
/// which it was.
pub(crate) fn sort_order(header: &sam::Header) -> Option<String> {
    use sam::header::record::value::map::header::tag::SORT_ORDER;

    let raw = header.header()?.other_fields().get(&SORT_ORDER)?;
    Some(String::from_utf8_lossy(raw.as_ref()).into_owned())
}

/// The header's `@SQ` list as a [`ContigList`] — name, length, and the `M5`
/// digest where the file carries one.
///
/// The order is the header's own, which is the whole point: the gate compares
/// this list against the reference's **including order**, and that is what
/// catches a permutation. Building it in header order is therefore load-bearing
/// rather than incidental.
///
/// **A missing `M5` is fine; a malformed one is fatal.** Spec §3.1 settles the
/// missing case: never an error, never a warning, because the digest check is
/// opportunistic and refusing a file for not offering a bonus would punish the
/// common case. A tag that is *present* but not 32 hex characters is a
/// different thing — a header its writer got wrong — and reading it as absent
/// would pass an error silently. The caller turns this into
/// `MalformedMd5`, naming the contig.
pub(crate) fn contig_list(header: &sam::Header) -> Result<ContigList, MalformedMd5Tag> {
    let entries = header
        .reference_sequences()
        .iter()
        .map(|(name, reference_sequence)| {
            let name = String::from_utf8_lossy(name.as_ref()).into_owned();
            Ok(ContigEntry {
                length: usize::from(reference_sequence.length()) as u64,
                md5: md5_tag(reference_sequence).map_err(|detail| MalformedMd5Tag {
                    contig: name.clone(),
                    detail,
                })?,
                name,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(ContigList { entries })
}

/// An `@SQ` line whose `M5` tag could not be read as a digest. Carries the
/// facts; the path and the phrasing belong to the caller's error.
#[derive(Debug)]
pub(crate) struct MalformedMd5Tag {
    pub(crate) contig: String,
    pub(crate) detail: String,
}

/// The `M5` tag of one `@SQ` entry, decoded from hex.
///
/// `Ok(None)` = no tag, which is ordinary. `Err(detail)` = a tag that is there
/// but unreadable, which is not.
fn md5_tag(
    reference_sequence: &sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
) -> Result<Option<[u8; 16]>, String> {
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let Some(hex) = reference_sequence.other_fields().get(&MD5_CHECKSUM) else {
        return Ok(None);
    };
    let hex = hex.as_ref();

    decode_md5_hex(hex).map(Some).ok_or_else(|| {
        if hex.len() != 32 {
            format!("expected 32 hex characters, got {}", hex.len())
        } else {
            "contains a non-hex character".to_string()
        }
    })
}

/// Decode 32 hex characters into the 16 digest bytes they spell. `None` for any
/// other length or a non-hex character.
fn decode_md5_hex(hex: &[u8]) -> Option<[u8; 16]> {
    let hex: &[u8; 32] = hex.try_into().ok()?;

    let mut digest = [0u8; 16];
    for (byte, pair) in digest.iter_mut().zip(hex.chunks_exact(2)) {
        *byte = (hex_nibble(pair[0])? << 4) | hex_nibble(pair[1])?;
    }
    Some(digest)
}

fn hex_nibble(character: u8) -> Option<u8> {
    match character {
        b'0'..=b'9' => Some(character - b'0'),
        b'a'..=b'f' => Some(character - b'a' + 10),
        b'A'..=b'F' => Some(character - b'A' + 10),
        _ => None,
    }
}

/// What the `@RG` records say the file's sample is.
///
/// The gate requires **exactly one** (spec §3.1, check 4), so each way of
/// failing that is a separate variant carrying the *facts* — never a rendered
/// message. Phrasing belongs to the caller's error type, which is also the only
/// layer that knows the file's path; and cross-file agreement is a different
/// check again, one layer further up, which needs the name rather than an
/// error.
pub(crate) enum SampleNames {
    /// Every `@RG` carries an `SM`, and they all agree.
    One(String),
    /// Two or more distinct `SM` values, in first-seen order.
    Several(Vec<String>),
    /// An `@RG` record carries no `SM` tag.
    MissingTag { read_group_id: String },
    /// The file has no `@RG` records at all.
    NoReadGroups,
}

/// Collect the distinct `@RG SM` values, in first-seen order.
///
/// First-seen rather than sorted, so a message lists the names the way the file
/// does and a reader can find them.
///
/// **The first fault in header order is the one reported**, matching
/// production: a file that both names two samples and has an `@RG` with no `SM`
/// stops at whichever comes first. Both are fatal at open and both are fixed by
/// correcting the header, so ranking them would add a rule to explain without
/// changing what the user does about it.
pub(crate) fn sample_names(header: &sam::Header) -> SampleNames {
    use sam::header::record::value::map::read_group::tag::SAMPLE;

    let mut distinct: Vec<String> = Vec::new();

    for (read_group_id, read_group) in header.read_groups() {
        let Some(raw) = read_group.other_fields().get(&SAMPLE) else {
            return SampleNames::MissingTag {
                read_group_id: String::from_utf8_lossy(read_group_id.as_ref()).into_owned(),
            };
        };
        let sample = String::from_utf8_lossy(raw.as_ref()).into_owned();
        if !distinct.contains(&sample) {
            distinct.push(sample);
        }
        if distinct.len() > 1 {
            return SampleNames::Several(distinct);
        }
    }

    match distinct.into_iter().next() {
        Some(name) => SampleNames::One(name),
        None => SampleNames::NoReadGroups,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::input::test_fixtures::{
        FIXTURE_CONTIGS, bam_header, fixture_reference, header, indexed_bam, matching_contigs,
        one_read, read_named, unindexed_bam,
    };

    // --- @HD SO ---

    #[test]
    fn sort_order_reads_the_tag_when_present() {
        assert_eq!(
            sort_order(&header(Some("coordinate"), &[], &[])).as_deref(),
            Some("coordinate")
        );
        assert_eq!(
            sort_order(&header(Some("queryname"), &[], &[])).as_deref(),
            Some("queryname"),
            "a wrong value is reported, not silently normalised — the gate \
             wants to name what it found"
        );
    }

    #[test]
    fn sort_order_is_none_when_the_tag_or_the_hd_line_is_absent() {
        assert_eq!(sort_order(&header(None, &[], &[])), None, "@HD with no SO");
        assert_eq!(
            sort_order(&sam::Header::default()),
            None,
            "no @HD line at all"
        );
    }

    // --- @SQ ---

    #[test]
    fn contig_list_preserves_header_order_with_names_and_lengths() {
        let contigs = contig_list(&header(
            Some("coordinate"),
            &[("chr2", 200, None), ("chr1", 100, None)],
            &[],
        ))
        .expect("no M5 tags to be malformed");

        let names: Vec<&str> = contigs.entries.iter().map(|e| e.name.as_str()).collect();
        assert_eq!(
            names,
            vec!["chr2", "chr1"],
            "header order is kept as-is — comparing it against the reference \
             *in order* is what catches a permutation"
        );
        assert_eq!(contigs.entries[0].length, 200);
        assert_eq!(contigs.entries[1].length, 100);
    }

    #[test]
    fn contig_list_decodes_an_m5_tag_and_tolerates_its_absence() {
        let contigs = contig_list(&header(
            Some("coordinate"),
            &[
                ("chr1", 100, Some("0123456789abcdef0123456789ABCDEF")),
                ("chr2", 200, None),
            ],
            &[],
        ))
        .expect("both tags are well formed");

        assert_eq!(
            contigs.entries[0].md5,
            Some([
                0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab,
                0xcd, 0xef
            ]),
            "decoded from hex, and upper case is accepted alongside lower"
        );
        assert_eq!(
            contigs.entries[1].md5, None,
            "a contig without M5 is ordinary input, not a fault"
        );
    }

    /// A malformed digest is a **hard error**, not a silently-absent one.
    /// Reading it as absent would pass a real defect silently — and a file
    /// whose tags were all malformed would look untagged, switching
    /// wrong-assembly detection off with nobody told.
    #[test]
    fn a_malformed_m5_is_rejected_naming_the_contig_and_the_reason() {
        let wrong_length = contig_list(&header(
            Some("coordinate"),
            &[("chr1", 100, Some("abc"))],
            &[],
        ))
        .expect_err("a three-character digest is malformed");
        assert_eq!(wrong_length.contig, "chr1");
        assert_eq!(wrong_length.detail, "expected 32 hex characters, got 3");

        let not_hex = contig_list(&header(
            Some("coordinate"),
            &[("chr2", 200, Some("zzzz56789abcdef0123456789abcdef0"))],
            &[],
        ))
        .expect_err("a non-hex character is malformed");
        assert_eq!(not_hex.contig, "chr2");
        assert_eq!(not_hex.detail, "contains a non-hex character");

        for bad in ["", &"a".repeat(33)] {
            assert!(
                contig_list(&header(
                    Some("coordinate"),
                    &[("chr1", 100, Some(bad))],
                    &[]
                ))
                .is_err(),
                "malformed M5 {bad:?} must be rejected"
            );
        }
    }

    /// The same fault through the gate, so the path and the phrasing are
    /// covered too — and so it is unmistakable that this rejects the *file*,
    /// not merely the tag.
    #[test]
    fn the_gate_rejects_a_file_with_a_malformed_m5() {
        let contigs = vec![("chr1", 100, Some("not-a-digest")), ("chr2", 200, None)];

        let error = open_fixture(&contigs, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::MalformedMd5 { contig, detail, .. } => {
                assert_eq!(contig, "chr1");
                assert_eq!(detail, "expected 32 hex characters, got 12");
            }
            other => panic!("expected MalformedMd5, got {other:?}"),
        }
        assert!(
            error
                .to_string()
                .contains("@SQ 'chr1' has a malformed M5 tag"),
            "{error}"
        );
    }

    // --- @RG SM ---

    #[test]
    fn one_sample_name_across_several_read_groups_is_one_name() {
        let names = sample_names(&header(
            Some("coordinate"),
            &[],
            &[("rg1", Some("NA12878")), ("rg2", Some("NA12878"))],
        ));
        assert!(matches!(names, SampleNames::One(name) if name == "NA12878"));
    }

    /// The single-file half of the sample check: two `@RG` records naming
    /// different samples is a file that cannot say whose reads it holds.
    #[test]
    fn two_distinct_sample_names_are_reported_in_first_seen_order() {
        let names = sample_names(&header(
            Some("coordinate"),
            &[],
            &[
                ("rg1", Some("NA12892")),
                ("rg2", Some("NA12878")),
                ("rg3", Some("NA12892")),
            ],
        ));
        match names {
            SampleNames::Several(names) => assert_eq!(
                names,
                vec!["NA12892", "NA12878"],
                "distinct, de-duplicated, and in the order the file lists them"
            ),
            _ => panic!("expected two distinct names"),
        }
    }

    #[test]
    fn a_file_with_no_read_groups_names_no_sample() {
        let names = sample_names(&header(Some("coordinate"), &[("chr1", 100, None)], &[]));
        assert!(matches!(names, SampleNames::NoReadGroups));
    }

    #[test]
    fn a_read_group_without_an_sm_tag_is_reported_by_id() {
        let names = sample_names(&header(
            Some("coordinate"),
            &[],
            &[("rg1", Some("NA12878")), ("rg2", None)],
        ));
        assert!(matches!(
            names,
            SampleNames::MissingTag { read_group_id } if read_group_id == "rg2"
        ));
    }

    /// Two faults at once: the **first in header order** is reported, matching
    /// production. Both orderings are pinned so the rule is the loop's stated
    /// behaviour rather than an accident of where the checks sit.
    #[test]
    fn the_first_sample_fault_in_header_order_is_the_one_reported() {
        let second_name_comes_first = sample_names(&header(
            Some("coordinate"),
            &[],
            &[("rg1", Some("A")), ("rg2", Some("B")), ("rg3", None)],
        ));
        match second_name_comes_first {
            SampleNames::Several(names) => assert_eq!(names, vec!["A", "B"]),
            _ => panic!("the second distinct name appears before the untagged group"),
        }

        let missing_tag_comes_first = sample_names(&header(
            Some("coordinate"),
            &[],
            &[("rg1", Some("A")), ("rg2", None), ("rg3", Some("B"))],
        ));
        assert!(
            matches!(
                missing_tag_comes_first,
                SampleNames::MissingTag { read_group_id } if read_group_id == "rg2"
            ),
            "the untagged group appears before the second distinct name"
        );
    }

    #[test]
    fn contig_list_is_empty_for_a_header_with_no_reference_sequences() {
        // B2 compares this against the reference, so an empty list is rejected
        // on length rather than on anything subtler.
        assert!(
            contig_list(&header(Some("coordinate"), &[], &[]))
                .expect("no tags at all")
                .entries
                .is_empty()
        );
    }

    // -----------------------------------------------------------------
    // The reader pool (C1)
    // -----------------------------------------------------------------

    /// **The guarantee the per-query cost model rests on.** Ten sequential
    /// borrows must open the file *once*: each returns its reader on drop, and
    /// the next takes that one back rather than opening another. If this
    /// regressed, ~10⁶ STR region queries would become ~10⁶ file opens, and
    /// nothing about the *answers* would change — which is why it is asserted
    /// rather than assumed.
    #[test]
    fn sequential_borrows_open_the_file_once_and_reuse_the_reader() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("opens");

        assert_eq!(
            file.readers_opened(),
            0,
            "opening the file opens no reader — a file that is never queried \
             costs no descriptor"
        );
        assert_eq!(file.pooled_readers(), 0);

        for _ in 0..10 {
            let borrowed = file.borrow_reader().expect("borrow");
            assert_eq!(
                file.pooled_readers(),
                0,
                "while borrowed, the pool is empty — the handle is out on loan"
            );
            drop(borrowed);
            assert_eq!(file.pooled_readers(), 1, "and comes back on drop");
        }

        assert_eq!(
            file.readers_opened(),
            1,
            "ten borrows, one open — the pool did its job"
        );
    }

    /// Two borrows held *at once* need two readers, because a handle out on
    /// loan cannot be lent again. Both come back, so the pool ends up holding
    /// both — which is what lets N threads query one file concurrently.
    #[test]
    fn concurrent_borrows_each_get_their_own_reader_and_all_return() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("opens");

        let first = file.borrow_reader().expect("borrow");
        let second = file.borrow_reader().expect("borrow");
        assert_eq!(
            file.readers_opened(),
            2,
            "the pool was empty for the second"
        );

        drop(first);
        drop(second);
        assert_eq!(file.pooled_readers(), 2);

        // A third borrow now reuses rather than opening a third.
        let third = file.borrow_reader().expect("borrow");
        assert_eq!(file.readers_opened(), 2);
        drop(third);
    }

    /// The pool exists so `reads_in_region(&self)` can be called from several
    /// threads at once, and none of the tests above actually does that. Each
    /// thread must get a reader of its own, and every one must come back: a
    /// handle lost under contention would leak a descriptor per query, and a
    /// lock held across the file open would show up here as a hang rather than
    /// a failure.
    #[test]
    fn borrows_from_many_threads_each_get_a_reader_and_all_return() {
        const THREADS: usize = 8;

        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("opens");

        std::thread::scope(|scope| {
            for _ in 0..THREADS {
                scope.spawn(|| {
                    let borrowed = file.borrow_reader().expect("borrow");
                    // Hold it long enough that the borrows genuinely overlap.
                    std::thread::yield_now();
                    drop(borrowed);
                });
            }
        });

        let opened = file.readers_opened();
        assert!(
            (1..=THREADS).contains(&opened),
            "at most one reader per concurrent borrow, at least one overall; \
             got {opened}"
        );
        assert_eq!(
            file.pooled_readers(),
            opened,
            "every reader opened came back — none lost on any path"
        );
    }

    /// The CRAM arm of `open_reader`, which no other test reaches. It also
    /// pins the interchangeability claim the borrow rests on: a freshly opened
    /// reader must be positioned past the header, exactly where a returned one
    /// would be, or the two could not be pooled together.
    #[test]
    fn a_cram_file_pools_readers_the_same_way() {
        use crate::pileup::per_sample::cram_files::{HeaderOverrides, build_cram};

        let specs: Vec<ContigSpec> = FIXTURE_CONTIGS
            .iter()
            .map(|(name, length)| ContigSpec {
                name: (*name).to_string(),
                length: *length as u64,
            })
            .collect();
        let (_fasta_dir, fasta) = build_fasta(&specs).expect("build fasta");
        let (_cram_dir, cram_path) = build_cram(
            &fasta,
            &specs,
            &HeaderOverrides {
                read_groups: vec![("rg1".to_string(), Some("NA12878".to_string()))],
                ..HeaderOverrides::default()
            },
            &[read_named("read-1", 0, 1)],
        )
        .expect("build cram");

        // The gate requires an index. This test is about pooling, not querying,
        // so an empty `.crai` is enough — C5 builds real ones.
        cram::crai::fs::write(
            cram_path.with_extension("cram.crai"),
            &cram::crai::Index::default(),
        )
        .expect("write an empty crai");

        // The `Fasta` arm, because a CRAM can only be decoded against the
        // bases — see `a_cram_against_a_fai_only_reference_is_refused_at_open`.
        let reference = read_reference_info(ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        })
        .expect("read reference");
        let file = AlignmentFile::open(&cram_path, &reference, ReadFilterConfig::default(), false)
            .expect("a matching CRAM opens");

        let borrowed = file.borrow_reader().expect("borrow a CRAM reader");
        assert!(
            matches!(
                borrowed.handle.as_ref().map(|h| &h.reader),
                Some(ReaderKind::Cram(_))
            ),
            "a .cram path yields a CRAM reader, not the BAM fallback"
        );
        drop(borrowed);

        assert_eq!(file.readers_opened(), 1, "the pool still works for CRAM");
        assert_eq!(file.pooled_readers(), 1);
    }

    /// A failed open must not count, or `readers_opened` stops meaning
    /// "readers that exist" and T13's assertion becomes unfalsifiable.
    #[test]
    fn a_failed_open_does_not_count_as_an_opened_reader() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("opens");

        // Delete the BAM out from under the handle, then ask for a reader.
        std::fs::remove_file(file.path()).expect("remove the bam");

        assert!(file.borrow_reader().is_err(), "the file is gone");
        assert_eq!(
            file.readers_opened(),
            0,
            "the counter tracks readers that exist, not attempts"
        );
    }

    /// `take()` is how a region stream adopts the parts; the borrow must then
    /// return nothing, or the same reader would be pushed back twice and handed
    /// to two callers at once.
    #[test]
    fn a_taken_handle_is_not_also_returned_by_the_borrow() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("opens");

        let handle = file.borrow_reader().expect("borrow").take();
        assert_eq!(
            file.pooled_readers(),
            0,
            "the borrow returned nothing — the stream owns the handle now"
        );

        file.return_handle(handle);
        assert_eq!(file.pooled_readers(), 1, "and the stream hands it back");
    }

    // -----------------------------------------------------------------
    // The composed chain — reads_in_region (C4): T9, T10, T13
    // -----------------------------------------------------------------

    use noodles_sam::alignment::RecordBuf;

    use crate::ng::read::input::test_fixtures::read_named_with_length;
    use crate::ng::ref_seq::InMemoryRefSeq;

    /// An all-`A` reference matching the fixture contigs, so a read of `A`s
    /// matches perfectly and a read of `C`s mismatches at every base.
    fn reference_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        )
    }

    /// A read whose bases all mismatch the all-`A` reference — filter #8's
    /// business, and only reachable if the *whole* filter is composed in.
    fn all_mismatching_read(qname: &str, start: usize) -> RecordBuf {
        let mut record = read_named_with_length(qname, 0, start, 30);
        record.sequence_mut().as_mut().fill(b'C');
        record
    }

    fn opened_over(records: &[RecordBuf]) -> (TempDir, TempDir, AlignmentFile) {
        let (reference_dir, reference) = fixture_reference(false);
        let (bam_dir, path) = indexed_bam(&bam_header(&matching_contigs()), records);
        let file = AlignmentFile::open(&path, &reference, ReadFilterConfig::default(), false)
            .expect("the fixture matches");
        (reference_dir, bam_dir, file)
    }

    fn whole_first_contig() -> GenomeRegion {
        GenomeRegion {
            contig: crate::ng::types::ContigId(0),
            start: crate::ng::types::Position(1),
            end: crate::ng::types::Position(FIXTURE_CONTIGS[0].1 as u64),
        }
    }

    /// **T9 — the full step-1 filter runs, not just the cheap subset.**
    ///
    /// A read that mismatches the reference at every base is filter #8's
    /// business, and #8 is the *reference-dependent* filter — the one a region
    /// reader that applied only flag/MAPQ checks would miss. Its drop being
    /// charged to `high_mismatch_fraction` is what proves `ReadFilter` is
    /// composed into the chain rather than a subset of it. This is the property
    /// that decided the rebuild over reusing production's reader.
    #[test]
    fn t9_the_served_stream_runs_the_reference_dependent_filter() {
        let (_reference_dir, _bam_dir, file) = opened_over(&[
            read_named_with_length("clean", 0, 1, 30),
            all_mismatching_read("mismatching", 40),
        ]);

        let reads: Vec<MappedRead> = file
            .reads_in_region(whole_first_contig(), reference_bases())
            .expect("query")
            .collect::<Result<_, _>>()
            .expect("no fatal error");

        assert_eq!(reads.len(), 1, "only the clean read survives");
        assert_eq!(reads[0].qname, b"clean");

        let counts = file.counts();
        assert_eq!(
            counts.high_mismatch_fraction, 1,
            "the drop is charged to filter #8, so the whole cascade ran"
        );
        assert_eq!(counts.kept, 1);
    }

    /// Chunk-edge slop is dropped **uncounted**: a read the index over-returned
    /// is not a read the filter rejected, and charging it to a `DropReason`
    /// would make the tally mean something different for an indexed read than
    /// for a whole-file one.
    #[test]
    fn records_outside_the_region_are_dropped_without_being_counted() {
        let (_reference_dir, _bam_dir, file) = opened_over(&[
            read_named_with_length("inside", 0, 1, 30),
            read_named_with_length("outside", 0, 60, 30),
        ]);

        let region = GenomeRegion {
            contig: crate::ng::types::ContigId(0),
            start: crate::ng::types::Position(1),
            end: crate::ng::types::Position(35),
        };
        let reads: Vec<MappedRead> = file
            .reads_in_region(region, reference_bases())
            .expect("query")
            .collect::<Result<_, _>>()
            .expect("no fatal error");

        assert_eq!(reads.len(), 1);
        let counts = file.counts();
        assert_eq!(counts.kept, 1);
        assert_eq!(
            counts.duplicate
                + counts.low_mapq
                + counts.supplementary
                + counts.secondary
                + counts.unmapped
                + counts.qc_fail
                + counts.too_short
                + counts.high_mismatch_fraction
                + counts.bad_cigar,
            0,
            "the out-of-region read is charged to no drop reason at all"
        );
    }

    /// **T13 — the file is opened once and the index parsed once**, however
    /// many regions are queried. The index cannot be re-parsed by construction
    /// (it is owned by the handle), so what is asserted is the observable half:
    /// many queries, one reader.
    ///
    /// This is the guarantee the whole per-query cost model rests on, and the
    /// one that decides whether ~10⁶ STR queries are affordable. A regression
    /// changes no answers — only cost — so nothing but an explicit assertion
    /// would catch it.
    #[test]
    fn t13_many_region_queries_open_the_file_once() {
        let (_reference_dir, _bam_dir, file) =
            opened_over(&[read_named_with_length("r", 0, 1, 30)]);

        for _ in 0..25 {
            let reads: Vec<MappedRead> = file
                .reads_in_region(whole_first_contig(), reference_bases())
                .expect("query")
                .collect::<Result<_, _>>()
                .expect("no fatal error");
            assert_eq!(reads.len(), 1);
        }

        assert_eq!(
            file.readers_opened(),
            1,
            "25 queries, one reader — the pool served every one after the first"
        );
        assert_eq!(file.pooled_readers(), 1, "and it is back in the pool");
    }

    /// **T10 — a fatal mid-stream error is yielded once, then the stream is
    /// done.** A truncated file must not look like a short region.
    #[test]
    fn t10_a_truncated_file_yields_one_error_and_then_nothing() {
        let (_reference_dir, reference) = fixture_reference(false);
        // Enough records to span several BGZF blocks: truncating a
        // single-block file would destroy the *header*, which is a different
        // fault (the gate's) from the mid-stream one this test is about.
        let mut records: Vec<RecordBuf> = Vec::new();
        for start in 1..=60 {
            for copy in 0..120 {
                records.push(read_named_with_length(
                    &format!("r{start}_{copy}"),
                    0,
                    start,
                    30,
                ));
            }
        }
        let (_bam_dir, path) = indexed_bam(&bam_header(&matching_contigs()), &records);

        let file = AlignmentFile::open(&path, &reference, ReadFilterConfig::default(), false)
            .expect("opens");

        // Truncate *after* opening, so the gate and the index are intact and
        // the fault can only appear mid-stream.
        let full = std::fs::metadata(&path).expect("stat").len();
        std::fs::OpenOptions::new()
            .write(true)
            .open(&path)
            .expect("reopen")
            .set_len(full * 3 / 4)
            .expect("truncate");

        let mut reads = file
            .reads_in_region(whole_first_contig(), reference_bases())
            .expect("the query itself still plans");

        let mut reads_before_error = 0;
        let mut errors = 0;
        let mut after_error = 0;
        while let Some(item) = reads.next() {
            match item {
                Ok(_) => reads_before_error += 1,
                Err(_) => {
                    errors += 1;
                    // Everything the iterator yields after the first error.
                    after_error = reads.by_ref().count();
                    break;
                }
            }
        }

        // Without this the test would pass against a chain that yielded
        // *nothing* — it would still reach the truncation and still fuse. That
        // is exactly the shape of the fixture bug this module already hit once.
        assert!(
            reads_before_error > 0,
            "reads must flow before the truncation is reached"
        );
        assert_eq!(errors, 1, "the truncation surfaced as a fatal error");
        assert_eq!(after_error, 0, "and the stream fused rather than resuming");

        // The reader still comes back on the error path — the case where its
        // state is least obvious.
        drop(reads);
        assert_eq!(
            file.pooled_readers(),
            1,
            "a stream that ended in an error still returns its reader"
        );
    }

    /// The signature exists so N threads can query one file at once, and
    /// nothing so far has actually done that through `reads_in_region` — the
    /// pool tests stop at `borrow_reader`. A `Rc` or `RefCell` creeping into
    /// the chain would break this and surface only at the first parallel call
    /// site, a step or two away.
    #[test]
    fn a_region_stream_is_send_so_threads_can_share_one_file() {
        fn assert_send<T: Send>() {}
        assert_send::<RegionReads<'_, InMemoryRefSeq>>();
        assert_send::<&AlignmentFile>();
    }

    /// And the same thing dynamically: eight concurrent queries, every reader
    /// returned, **every tally folded in**. The counts assertion is the
    /// valuable half — it is what would catch a fold lost under contention.
    #[test]
    fn concurrent_region_queries_each_get_a_reader_and_bank_every_tally() {
        let (_reference_dir, _bam_dir, file) =
            opened_over(&[read_named_with_length("r", 0, 1, 30)]);

        std::thread::scope(|scope| {
            for _ in 0..8 {
                scope.spawn(|| {
                    let reads: Vec<MappedRead> = file
                        .reads_in_region(whole_first_contig(), reference_bases())
                        .expect("query")
                        .collect::<Result<_, _>>()
                        .expect("no fatal error");
                    assert_eq!(reads.len(), 1);
                });
            }
        });

        assert_eq!(
            file.counts().kept,
            8,
            "every stream folded its tally in — none lost to the race"
        );
        assert_eq!(
            file.pooled_readers(),
            file.readers_opened(),
            "and every reader opened came back"
        );
    }

    /// **T4d, re-sited from C3.** The order guard's state lives in the
    /// per-region iterator, so querying a later region and then an earlier one
    /// **on the same handle** is a new forward scan, not a regression. C3 could
    /// only assert this with two separate guards, which could not fail if the
    /// state ever migrated to the handle. This can.
    #[test]
    fn t4d_a_later_region_then_an_earlier_one_on_one_handle_is_not_a_regression() {
        let (_reference_dir, _bam_dir, file) = opened_over(&[
            read_named_with_length("early", 0, 1, 30),
            read_named_with_length("late", 0, 60, 30),
        ]);

        let later = GenomeRegion {
            contig: crate::ng::types::ContigId(0),
            start: crate::ng::types::Position(55),
            end: crate::ng::types::Position(100),
        };
        let earlier = GenomeRegion {
            contig: crate::ng::types::ContigId(0),
            start: crate::ng::types::Position(1),
            end: crate::ng::types::Position(40),
        };

        let first: Vec<MappedRead> = file
            .reads_in_region(later, reference_bases())
            .expect("query")
            .collect::<Result<_, _>>()
            .expect("the later region streams");
        assert_eq!(first.len(), 1);

        let second: Vec<MappedRead> = file
            .reads_in_region(earlier, reference_bases())
            .expect("query")
            .collect::<Result<_, _>>()
            .expect("an earlier region on the same handle is a new forward scan");
        assert_eq!(second.len(), 1);
    }

    /// A query that fails to plan must not cost a reader — the pool is only
    /// touched once everything fallible has succeeded.
    #[test]
    fn a_failed_query_returns_no_reader_because_it_never_took_one() {
        let (_reference_dir, _bam_dir, file) =
            opened_over(&[read_named_with_length("r", 0, 1, 30)]);

        let nonexistent_contig = GenomeRegion {
            contig: crate::ng::types::ContigId(9),
            start: crate::ng::types::Position(1),
            end: crate::ng::types::Position(10),
        };
        assert!(matches!(
            file.reads_in_region(nonexistent_contig, reference_bases()),
            Err(AlignmentFileError::Region { .. })
        ));
        assert_eq!(
            file.readers_opened(),
            0,
            "planning failed before the pool was touched"
        );
    }

    /// Abandoning a stream half-way still returns the reader, and still banks
    /// the drops recorded so far — both happen in `Drop`, so neither depends on
    /// the caller draining.
    #[test]
    fn abandoning_a_stream_returns_the_reader_and_banks_its_counts() {
        let (_reference_dir, _bam_dir, file) = opened_over(&[
            all_mismatching_read("bad", 1),
            read_named_with_length("good", 0, 40, 30),
        ]);

        {
            let mut reads = file
                .reads_in_region(whole_first_contig(), reference_bases())
                .expect("query");
            // Pull exactly one read, then walk away.
            let _ = reads.next();
            assert_eq!(file.pooled_readers(), 0, "the reader is out on loan");
        }

        assert_eq!(file.pooled_readers(), 1, "returned on drop");
        assert_eq!(
            file.counts().high_mismatch_fraction,
            1,
            "the drop seen before abandoning was banked, not lost"
        );
    }

    // -----------------------------------------------------------------
    // BAM/CRAM parity (C5) — T8
    // -----------------------------------------------------------------

    use crate::ng::read::input::test_fixtures::indexed_cram;

    /// A spread of reads with pile-ups, so parity means more than "both
    /// returned nothing".
    ///
    /// One contig, because the CRAM fixture cannot hold more — see
    /// `test_fixtures::indexed_cram` for why (a noodles indexing limitation,
    /// not a choice).
    fn parity_records() -> Vec<RecordBuf> {
        let mut records = Vec::new();
        let mut start = 1;
        while start + 30 < FIXTURE_CONTIGS[0].1 {
            records.push(read_named_with_length(&format!("r{start}"), 0, start, 30));
            if start % 20 == 1 {
                records.push(read_named_with_length(&format!("r{start}_b"), 0, start, 30));
            }
            start += 7;
        }
        records
    }

    fn names_from(file: &AlignmentFile, region: GenomeRegion) -> Vec<String> {
        file.reads_in_region(region, reference_bases())
            .expect("query")
            .map(|read| String::from_utf8_lossy(&read.expect("no fatal error").qname).into_owned())
            .collect()
    }

    /// **T8 — the same reads written as BAM and as CRAM produce the same
    /// ordered stream.**
    ///
    /// The two containers share nothing below the `RecordSource` seam: BAM
    /// reads one record at a time from bgzf chunks, CRAM decodes whole
    /// containers against the reference and walks a `.crai`. Everything above —
    /// the filter, the order guard — is the same code. So an disagreement here
    /// is a container-reader bug, which is exactly what this is looking for,
    /// and BAM is the oracle because it is the simpler reader and was verified
    /// first against a linear scan (T5).
    #[test]
    fn t8_a_cram_yields_the_same_ordered_reads_as_the_same_bam() {
        let records = parity_records();

        let (_reference_dir, bam_reference) = fixture_reference(false);
        let (_bam_dir, bam_path) = indexed_bam(&bam_header(&matching_contigs()), &records);
        let bam = AlignmentFile::open(
            &bam_path,
            &bam_reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("the BAM opens");

        let (_cram_dir, cram_path, _fasta_dir, fasta) = indexed_cram(&records);
        let cram_reference = read_reference_info(ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        })
        .expect("read reference");
        let cram = AlignmentFile::open(
            &cram_path,
            &cram_reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("the CRAM opens");

        let regions = [
            (0u32, 1u64, 100u64),
            (0, 1, 30),
            (0, 45, 55),
            (0, 60, 62),
            (0, 90, 100),
        ];

        let mut total = 0;
        for (contig, start, end) in regions {
            let region = GenomeRegion {
                contig: crate::ng::types::ContigId(contig),
                start: crate::ng::types::Position(start),
                end: crate::ng::types::Position(end),
            };

            let from_bam = names_from(&bam, region);
            let from_cram = names_from(&cram, region);

            assert_eq!(
                from_cram, from_bam,
                "BAM and CRAM disagreed for contig {contig} [{start}, {end}]"
            );
            total += from_bam.len();
        }

        assert!(
            total > 20,
            "the fixture must actually cover these regions — {total} reads is \
             too few for the comparison to mean anything"
        );
    }

    /// **The owner's rule, enforced at open.** A CRAM stores sequences as
    /// differences from the reference, so it cannot be decoded from a `.fai`,
    /// which holds only geometry. That is a hard error the moment the file is
    /// opened — not a mystery at the first query, and not a silently empty
    /// stream.
    #[test]
    fn a_cram_against_a_fai_only_reference_is_refused_at_open() {
        let (_cram_dir, cram_path, _fasta_dir, fasta) = indexed_cram(&one_read());

        let fai_only = read_reference_info(ReferenceSource::Fai(
            crate::ng::reference_info::sibling_fai_path(&fasta),
        ))
        .expect("read reference");

        let error = AlignmentFile::open(&cram_path, &fai_only, ReadFilterConfig::default(), false)
            .expect_err("a CRAM needs the bases");
        assert!(matches!(
            error,
            AlignmentFileError::CramNeedsReferenceFasta { .. }
        ));
        assert!(
            error.to_string().contains("supply the reference FASTA"),
            "the message must say what to do about it: {error}"
        );
    }

    /// A BAM is unaffected — it stores its own sequences, so a `.fai`-only
    /// reference is perfectly ordinary input for it.
    #[test]
    fn a_bam_against_a_fai_only_reference_opens_normally() {
        let (_reference_dir, _bam_dir, file) =
            opened_over(&[read_named_with_length("r", 0, 1, 30)]);
        assert_eq!(file.sample_name(), "NA12878");
    }

    /// **The `.crai` walk, with more than one entry to walk.**
    ///
    /// The parity fixture fits in a single container, so it exercises the
    /// container decode but none of the walk's control flow. This one spans
    /// several containers over a long contig, so a query for a slice of it must
    /// visit some entries, skip others, and **stop early** rather than decoding
    /// every container to the end — which is the difference between a seek and
    /// a whole-file scan, ~10⁶ times over.
    #[test]
    fn a_multi_container_cram_walks_its_crai_and_stops_early() {
        use crate::ng::read::input::test_fixtures::multi_container_cram;

        const CONTIG_LENGTH: usize = 400_000;
        // Comfortably over noodles' 10240 records per container.
        const READS: usize = 30_000;

        let (_cram_dir, cram_path, _fasta_dir, fasta) = multi_container_cram(CONTIG_LENGTH, READS);
        let reference = read_reference_info(ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        })
        .expect("read reference");
        let file = AlignmentFile::open(&cram_path, &reference, ReadFilterConfig::default(), false)
            .expect("opens");

        let bases = || InMemoryRefSeq::from_contigs(vec![vec![b'A'; CONTIG_LENGTH]]);
        let near_the_start = GenomeRegion {
            contig: crate::ng::types::ContigId(0),
            start: crate::ng::types::Position(1),
            end: crate::ng::types::Position(200),
        };

        let reads: Vec<MappedRead> = file
            .reads_in_region(near_the_start, bases())
            .expect("query")
            .collect::<Result<_, _>>()
            .expect("no fatal error");

        assert!(
            !reads.is_empty(),
            "the region is covered, so the walk must find its container"
        );
        for read in &reads {
            assert!(
                read.pos <= 200 && read.pos + 29 >= 1,
                "a read outside the region survived: pos {}",
                read.pos
            );
        }

        // Deep in the covered range, to prove the walk reaches later containers
        // rather than only ever decoding the first. (Well inside where reads
        // actually are — the spread stops short of the contig end.)
        let late_start = (CONTIG_LENGTH as u64 * 3) / 4;
        let near_the_end = GenomeRegion {
            contig: crate::ng::types::ContigId(0),
            start: crate::ng::types::Position(late_start),
            end: crate::ng::types::Position(late_start + 200),
        };
        let late: Vec<MappedRead> = file
            .reads_in_region(near_the_end, bases())
            .expect("query")
            .collect::<Result<_, _>>()
            .expect("no fatal error");
        assert!(
            !late.is_empty(),
            "the walk must reach the contig's last containers too"
        );
        assert!(
            late[0].pos >= late_start - 30,
            "and those are genuinely different reads from the first query"
        );
    }

    // --- the hex decoder's edges ---

    /// The four malformed cases above put the bad character early or vary the
    /// length; none proves the **last** pair is inspected. A 32-character string
    /// that is valid hex except for its final character catches a decoder that
    /// stops one pair short — the off-by-one the `zip`/`chunks_exact` pairing
    /// could plausibly hide.
    #[test]
    fn decode_md5_hex_rejects_a_non_hex_character_in_the_final_position() {
        let mut hex = [b'0'; 32];
        assert!(decode_md5_hex(&hex).is_some(), "all-zeros is valid hex");

        hex[31] = b'z';
        assert!(decode_md5_hex(&hex).is_none(), "last nibble is inspected");

        let mut hex = [b'0'; 32];
        hex[30] = b'z';
        assert!(
            decode_md5_hex(&hex).is_none(),
            "high nibble of the last byte is inspected too"
        );
    }

    /// Range patterns fail at their edges, and a `z`-only test cannot see it.
    /// These six characters each sit one step outside a valid range.
    // -----------------------------------------------------------------
    // The gate — `AlignmentFile::open` (T1, T2, T3, T12a)
    // -----------------------------------------------------------------
    use tempfile::TempDir;

    use crate::ng::reference_info::{ReferenceSource, read_reference_info};
    use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};

    /// An opened fixture **plus the temp dirs its files live in**.
    ///
    /// The dirs are returned rather than dropped at helper exit because they
    /// own the BAM and the reference on disk: drop them and the handle points
    /// at deleted files. Today's assertions would survive that (the index is
    /// already parsed in memory), but the first test that actually *queries*
    /// the file would fail somewhere far from the cause.
    struct OpenedFixture {
        file: Result<AlignmentFile, AlignmentFileError>,
        _reference_dir: TempDir,
        _bam_dir: TempDir,
    }

    fn open_fixture(
        contigs: &[(&str, usize, Option<&str>)],
        reference_has_digests: bool,
    ) -> OpenedFixture {
        open_fixture_with_header(&bam_header(contigs), reference_has_digests)
    }

    fn open_fixture_with_header(
        header: &sam::Header,
        reference_has_digests: bool,
    ) -> OpenedFixture {
        let (reference_dir, reference) = fixture_reference(reference_has_digests);
        let (bam_dir, path) = indexed_bam(header, &one_read());
        OpenedFixture {
            file: AlignmentFile::open(&path, &reference, ReadFilterConfig::default(), false),
            _reference_dir: reference_dir,
            _bam_dir: bam_dir,
        }
    }

    /// The `ContigReconcile` detail for a header that should fail check 2.
    fn reconcile_detail(
        contigs: &[(&str, usize, Option<&str>)],
        reference_has_digests: bool,
    ) -> String {
        match open_fixture(contigs, reference_has_digests)
            .file
            .expect_err("must not open")
        {
            AlignmentFileError::ContigReconcile { detail, .. } => detail,
            other => panic!("expected ContigReconcile, got {other:?}"),
        }
    }

    #[test]
    fn a_matching_file_opens_and_exposes_its_sample_and_digests() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("the file matches");

        assert_eq!(file.sample_name(), "NA12878");
        assert_eq!(
            file.sq_md5s().len(),
            2,
            "one slot per contig, indexed by ContigId"
        );
        assert!(
            file.sq_md5s().iter().all(Option::is_none),
            "this fixture's @SQ lines carry no M5, which is ordinary input"
        );
    }

    /// **T1 — the permutation hole, closed.**
    ///
    /// `@SQ` lists the reference's contigs in the wrong order. Every `ref_id`
    /// still *resolves*, which is why the check this replaces let it through
    /// and then fetched the wrong contig for every read. Order-aware equality
    /// catches it on the first transposed name.
    ///
    /// Mutation-verified below by
    /// `a_resolves_only_check_accepts_the_permutation_this_gate_rejects`.
    #[test]
    fn t1_a_permuted_sq_list_is_rejected_naming_the_first_transposed_contig() {
        let permuted = vec![("chr2", 200, None), ("chr1", 100, None)];

        let detail = reconcile_detail(&permuted, false);

        assert_eq!(
            detail, "name disagreement at index 0 ('chr2' vs 'chr1')",
            "the first transposed position, file value before reference value"
        );
    }

    /// **Why T1 needs an ordered comparison — the superseded probe, recorded.**
    ///
    /// `ReadFilter::new`'s check asks only "does every `@SQ` index resolve in
    /// the reference?", fetching each contig and accepting if the fetch
    /// succeeds. Run here against the permuted header, it **passes** — order is
    /// never consulted, so the file goes on to fetch the wrong contig for every
    /// read.
    ///
    /// **This is documentation, not the mutation guard.** It does not call
    /// `AlignmentFile::open`, so gutting the gate's comparison would leave it
    /// green; T1 is what catches that. The real mutation was performed against
    /// the gate itself during B2 — the `first_disagreement` call was replaced
    /// with a resolves-only length check and T1 duly failed — and this test is
    /// the standing record of *what* the old check did, so the reason the gate
    /// is stricter does not have to be taken on trust.
    #[test]
    fn the_superseded_resolves_only_probe_cannot_see_order() {
        use crate::ng::ref_seq::{InMemoryRefSeq, RefSeq};
        use crate::ng::types::ContigId;

        let (_reference_dir, reference) = fixture_reference(false);
        let permuted = bam_header(&[("chr2", 200, None), ("chr1", 100, None)]);

        // The superseded check: fetch every contig the file's `@SQ` list can
        // name, and accept if each one resolves. Order is never consulted.
        let reference_bases = InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        );
        let mut probe = Vec::new();
        let every_contig_resolves = (0..permuted.reference_sequences().len()).all(|index| {
            reference_bases
                .fetch_into(ContigId(index as u32), 1, 0, &mut probe)
                .is_ok()
        });
        assert!(
            every_contig_resolves,
            "every @SQ index of the permuted file resolves — so the old probe \
             accepts it, and then every read fetches the wrong contig"
        );

        // The gate's own comparison rejects that very same header.
        assert!(
            contig_list(&permuted)
                .expect("no M5 tags")
                .first_disagreement(&reference.contig_list())
                .is_err(),
            "the ordered comparison must reject what the probe accepted"
        );
    }

    /// **T2 — name, length and count mismatches**, each naming the right field
    /// and index.
    #[test]
    fn t2_name_length_and_count_mismatches_are_each_named() {
        let cases = [
            (
                vec![("chrX", 100, None), ("chr2", 200, None)],
                "name disagreement at index 0",
            ),
            (
                vec![("chr1", 100, None), ("chr2", 999, None)],
                "length disagreement at index 1",
            ),
            (
                vec![("chr1", 100, None)],
                "@SQ list length differs (1 vs 2)",
            ),
        ];

        for (contigs, expected) in cases {
            let detail = reconcile_detail(&contigs, false);
            assert!(
                detail.contains(expected),
                "expected {expected:?}, got {detail:?}"
            );
        }
    }

    /// **T2, the digest half.** A wrong `@SQ M5` is caught at open **only**
    /// when the `ReferenceInfo` in hand carries digests — the `Fasta` arm. Under
    /// a `.fai`-only table the comparison is a wildcard and the same file opens
    /// cleanly, which is the behaviour production ships and what makes the
    /// deferred `check_assembly` (D1) necessary rather than redundant.
    #[test]
    fn t2_a_wrong_m5_is_caught_only_when_the_reference_carries_digests() {
        let wrong_digest = "ffffffffffffffffffffffffffffffff";
        let contigs = vec![
            ("chr1", 100, Some(wrong_digest)),
            ("chr2", 200, Some(wrong_digest)),
        ];

        let detail = reconcile_detail(&contigs, true);
        assert!(
            detail.contains("md5 disagreement at index 0"),
            "got: {detail}"
        );

        assert!(
            open_fixture(&contigs, false).file.is_ok(),
            "against a .fai-only reference the digest is a wildcard, so the \
             same file opens — the gap check_assembly closes later"
        );
    }

    /// **T3 — `SO` wrong or missing**, rejected before anything else, and the
    /// message says which it was.
    #[test]
    fn t3_a_file_that_is_not_coordinate_sorted_is_rejected_at_open() {
        let contigs = matching_contigs();

        let queryname = header(Some("queryname"), &contigs, &[("rg1", Some("NA12878"))]);
        let error = open_fixture_with_header(&queryname, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::NotCoordinateSorted { sort_order, .. } => assert_eq!(
                sort_order.as_deref(),
                Some("queryname"),
                "carries the observed value as a fact, not pre-quoted prose"
            ),
            other => panic!("expected NotCoordinateSorted, got {other:?}"),
        }
        assert!(
            error.to_string().contains("@HD SO is 'queryname'"),
            "and the rendered message quotes it: {error}"
        );

        let no_sort_order = header(None, &contigs, &[("rg1", Some("NA12878"))]);
        let error = open_fixture_with_header(&no_sort_order, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::NotCoordinateSorted { sort_order, .. } => {
                assert_eq!(sort_order.as_deref(), None)
            }
            other => panic!("expected NotCoordinateSorted, got {other:?}"),
        }
        assert!(error.to_string().contains("@HD SO is missing"), "{error}");
    }

    /// The gate's fail-fast order is spec §3.1's — `SO` → `@SQ` → index → `SM`
    /// — and it is only observable through a file that fails two checks at
    /// once. Each boundary gets its own case, so re-ordering any adjacent pair
    /// fails a test rather than quietly changing which fault a user is told
    /// about.
    #[test]
    fn the_four_checks_run_in_the_specified_order() {
        // SO before @SQ.
        let bad_sort_and_bad_contigs = header(
            Some("queryname"),
            &[("chrX", 1, None)],
            &[("rg1", Some("NA12878"))],
        );
        assert!(matches!(
            open_fixture_with_header(&bad_sort_and_bad_contigs, false).file,
            Err(AlignmentFileError::NotCoordinateSorted { .. })
        ));

        // @SQ before the index: bad contigs on a BAM with no index at all.
        let bad_contigs = bam_header(&[("chrX", 1, None)]);
        assert!(matches!(
            open_unindexed_with_header(&bad_contigs).1,
            Err(AlignmentFileError::ContigReconcile { .. })
        ));

        // The index before @RG SM: no index and no sample.
        let no_sample = header(Some("coordinate"), &matching_contigs(), &[]);
        assert!(matches!(
            open_unindexed_with_header(&no_sample).1,
            Err(AlignmentFileError::Index { .. })
        ));
    }

    /// Write a BAM with **no** index beside it and try to open it. Returns the
    /// temp dirs so the file outlives the call.
    fn open_unindexed_with_header(
        header: &sam::Header,
    ) -> (
        (TempDir, TempDir),
        Result<AlignmentFile, AlignmentFileError>,
    ) {
        let (reference_dir, reference) = fixture_reference(false);
        let records = if header.reference_sequences().is_empty() {
            Vec::new()
        } else {
            one_read()
        };
        let (dir, path) = unindexed_bam(header, &records);

        let opened = AlignmentFile::open(&path, &reference, ReadFilterConfig::default(), false);
        ((reference_dir, dir), opened)
    }

    /// **T12a — one file naming two samples.** The cross-file half (two files
    /// naming different samples) is `SampleReads`' T12b.
    #[test]
    fn t12a_a_file_whose_read_groups_name_two_samples_is_rejected_at_open() {
        let contigs = matching_contigs();
        let two_samples = header(
            Some("coordinate"),
            &contigs,
            &[("rg1", Some("NA12878")), ("rg2", Some("NA12892"))],
        );

        match open_fixture_with_header(&two_samples, false)
            .file
            .expect_err("must not open")
        {
            AlignmentFileError::MultipleSampleNames { names, .. } => {
                assert_eq!(names, vec!["NA12878", "NA12892"])
            }
            other => panic!("expected MultipleSampleNames, got {other:?}"),
        }
    }

    /// The gate's mapping of `SampleNames::MissingTag`, and the `Some` branch
    /// of `MissingSampleName`'s message — the only non-trivial formatting in
    /// the new variant, and until now exercised by nothing.
    #[test]
    fn a_read_group_with_no_sm_tag_is_rejected_naming_that_read_group() {
        let contigs = matching_contigs();
        let untagged = header(Some("coordinate"), &contigs, &[("rg1", None)]);

        let error = open_fixture_with_header(&untagged, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::MissingSampleName { read_group, .. } => {
                assert_eq!(read_group.as_deref(), Some("rg1"))
            }
            other => panic!("expected MissingSampleName, got {other:?}"),
        }
        assert!(
            error.to_string().contains("@RG 'rg1' has no SM tag"),
            "{error}"
        );
    }

    /// The other half of "exactly one sample": a file with no `@RG` at all
    /// cannot say whose reads it holds, so the sample layer would have nothing
    /// to agree on.
    #[test]
    fn a_file_naming_no_sample_is_rejected_at_open() {
        let contigs = matching_contigs();
        let no_read_groups = header(Some("coordinate"), &contigs, &[]);

        match open_fixture_with_header(&no_read_groups, false)
            .file
            .expect_err("must not open")
        {
            AlignmentFileError::MissingSampleName { read_group, .. } => {
                assert_eq!(read_group, None)
            }
            other => panic!("expected MissingSampleName, got {other:?}"),
        }
    }

    /// A missing index is an error, not a silent whole-file scan — and with the
    /// build flag it is repaired instead.
    #[test]
    fn a_missing_index_is_rejected_unless_the_caller_asks_for_one() {
        let header = bam_header(&matching_contigs());
        let (dirs, opened) = open_unindexed_with_header(&header);

        assert!(
            matches!(opened, Err(AlignmentFileError::Index { .. })),
            "an unindexed file is an error, never a silent whole-file scan"
        );

        // The same file in the same place — only the caller's policy differs.
        let (_reference_dir, bam_dir) = &dirs;
        let (_fresh_reference_dir, reference) = fixture_reference(false);
        let path = bam_dir.path().join("sample.bam");
        assert!(
            AlignmentFile::open(&path, &reference, ReadFilterConfig::default(), true).is_ok(),
            "with build_index_if_missing the index is created next to the file"
        );
    }

    #[test]
    fn hex_nibble_rejects_the_characters_bordering_each_valid_range() {
        for character in [b'/', b':', b'`', b'g', b'@', b'G'] {
            assert_eq!(
                hex_nibble(character),
                None,
                "{:?} borders a valid range but is not hex",
                character as char
            );
        }
        assert_eq!(hex_nibble(b'0'), Some(0));
        assert_eq!(hex_nibble(b'9'), Some(9));
        assert_eq!(hex_nibble(b'a'), Some(10));
        assert_eq!(hex_nibble(b'f'), Some(15));
        assert_eq!(hex_nibble(b'A'), Some(10));
        assert_eq!(hex_nibble(b'F'), Some(15));
    }
}
