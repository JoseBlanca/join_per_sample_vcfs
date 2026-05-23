//! Production `RefSeqFetcher` impl that bounds memory by evicting
//! the previous chromosome's bases when the walker advances to a
//! new one.
//!
//! Background: `noodles_fasta::Repository` keeps a `HashMap<Vec<u8>,
//! Arc<Sequence>>` cache that has no eviction policy. A naive impl
//! that builds one Repository at startup and shares it across
//! chromosomes would accumulate the bases of every visited contig
//! for the lifetime of the run — roughly 3 GB across the 24 human
//! chromosomes. The walker is sequential single-pass per chromosome
//! (per `ia/specs/pileup_walker.md` §"Chromosome boundaries"), so
//! steady-state we only need the *current* contig in memory; the
//! previous one's `Arc<Sequence>` can be dropped the moment the
//! walker advances. See finding S6 in
//! `ia/reviews/pileup_samtools_comparison_2026-05-07.md`.

use std::cell::Cell;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use noodles_core::Region;
use noodles_fasta as fasta;
use noodles_fasta::fai;

use super::cram_input::ContigList;
use super::pileup::RefSeqFetcher;

/// Sliding-buffer size used by [`StreamingChromRefFetcher`]. 100× the
/// `--var-group-max-span` default (10 KB), so under the
/// monotonic-forward access pattern of the per-group merger a buffer
/// refill fires only once per ~100 group fetches in dense regions.
/// Per-worker peak overhead is therefore ~1 MB regardless of contig
/// size — vs Phase B's full-contig footprint (~91 MB on tomato ch01,
/// ~250 MB on human chr1).
pub const STREAMING_REF_BUFFER_BYTES: usize = 1024 * 1024;

/// FASTA-read buffer used by [`StreamingChromRefFetcher::refill`] for
/// raw on-disk bytes (sequence + newlines) before they're stripped
/// and uppercased into the streamer's `buf`. 64 KiB matches the
/// streaming MD5 verify in
/// [`crate::pop_var_caller::common::compute_contig_md5_streaming`]
/// and the project-wide buffered I/O capacity.
const STREAMING_REF_FILE_READ_CHUNK: usize = 64 * 1024;

/// `RefSeqFetcher` that holds a single `noodles_fasta::Repository`
/// for the whole run and clears its cache whenever the walker moves
/// to a new chromosome. Steady-state cache size is exactly one
/// contig regardless of how many chromosomes have been visited.
pub struct ChromBoundaryRefFetcher {
    repository: fasta::Repository,
    contigs: ContigList,
    /// Last `chrom_id` seen by `fetch`. `Cell` (not `RefCell`) is
    /// enough because it only holds a `Copy` `u32`. The walker is
    /// single-threaded per call (see S6 §"Risk"), so non-`Sync`
    /// interior mutability is fine.
    current_chrom: Cell<Option<u32>>,
}

impl ChromBoundaryRefFetcher {
    /// Open the `.fa` (and its sibling `.fai`) at `fasta_path` and
    /// build a fetcher that resolves `chrom_id` → contig name via
    /// `contigs`.
    pub fn new(fasta_path: &Path, contigs: ContigList) -> io::Result<Self> {
        let indexed_reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
        let repository = fasta::Repository::new(adapter);
        Ok(Self {
            repository,
            contigs,
            current_chrom: Cell::new(None),
        })
    }

    /// Number of contigs currently held in the underlying
    /// `Repository`'s cache. With the chromosome-boundary eviction
    /// in place, this is 0 before the first fetch and 1 thereafter.
    /// Exposed so tests can pin the eviction invariant.
    pub fn cached_contig_count(&self) -> usize {
        self.repository.len()
    }
}

impl RefSeqFetcher for ChromBoundaryRefFetcher {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        if self.current_chrom.get() != Some(chrom_id) {
            // Walker moved to a new chromosome (or this is the
            // first fetch of the run). Drop the previous contig's
            // bases before loading the next.
            self.repository.clear();
            self.current_chrom.set(Some(chrom_id));
        }

        fetch_from_repository(
            &self.repository,
            &self.contigs,
            chrom_id,
            start_1based,
            length,
        )
    }
}

// ---------------------------------------------------------------------
// SyncRefFetcher — Sync-safe, no eviction. For the BAQ stage.
// ---------------------------------------------------------------------

/// `RefSeqFetcher` variant for the rayon-parallel BAQ stage, which
/// requires a `Sync` fetcher to share across worker threads.
///
/// Trades eviction for thread-safety: under the hood it is the same
/// `noodles_fasta::Repository` (already `Sync` thanks to its
/// internal `Arc<RwLock<...>>`) but without the per-fetch
/// chromosome-boundary check that makes [`ChromBoundaryRefFetcher`]
/// non-`Sync`. The pipeline uses both:
///
/// - BAQ stage → [`SyncRefFetcher`] (parallel; non-evicting cache).
/// - Pileup walker → [`ChromBoundaryRefFetcher`] (sequential;
///   one-chrom-resident cache).
///
/// Memory: this fetcher's cache grows to hold every chromosome the
/// run visits (~3 GB on a 24-chrom human reference). For Stage 1's
/// per-sample whole-CRAM scan that is the dominant memory user, but
/// it is bounded and predictable. A future slice can introduce a
/// thread-safe chrom-boundary fetcher if the doubled cache footprint
/// becomes a problem; today the two-fetcher split is the simpler
/// design.
pub struct SyncRefFetcher {
    repository: fasta::Repository,
    contigs: ContigList,
}

impl SyncRefFetcher {
    pub fn new(fasta_path: &Path, contigs: ContigList) -> io::Result<Self> {
        let indexed_reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
        let adapter = fasta::repository::adapters::IndexedReader::new(indexed_reader);
        let repository = fasta::Repository::new(adapter);
        Ok(Self {
            repository,
            contigs,
        })
    }
}

impl RefSeqFetcher for SyncRefFetcher {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        fetch_from_repository(
            &self.repository,
            &self.contigs,
            chrom_id,
            start_1based,
            length,
        )
    }
}

/// Shared body of both fetchers — read a window from a
/// noodles `Repository` after validating `chrom_id` and 1-based
/// coordinates. Promoted to a free function so the two fetchers
/// agree on the error shapes (test invariants checked against
/// `ChromBoundaryRefFetcher` cover `SyncRefFetcher` by reuse).
fn fetch_from_repository(
    repository: &fasta::Repository,
    contigs: &ContigList,
    chrom_id: u32,
    start_1based: u32,
    length: u32,
) -> Result<Vec<u8>, io::Error> {
    let entry = contigs.entries.get(chrom_id as usize).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "chrom_id {chrom_id} out of range (have {} contigs)",
                contigs.entries.len()
            ),
        )
    })?;

    let seq_arc = match repository.get(entry.name.as_bytes()) {
        None => {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("contig {} not in FASTA", entry.name),
            ));
        }
        Some(Err(e)) => return Err(e),
        Some(Ok(seq)) => seq,
    };
    let bytes = seq_arc.as_ref().as_ref();

    if start_1based == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "start_1based must be >= 1",
        ));
    }
    let start_idx = (start_1based - 1) as usize;
    let end_idx = start_idx
        .checked_add(length as usize)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"))?;
    if end_idx > bytes.len() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!(
                "fetch [{}, {}) past contig {} length {}",
                start_1based,
                start_1based + length,
                entry.name,
                bytes.len()
            ),
        ));
    }
    // Soft-masked FASTAs (Ensembl/Gencode default) encode repeat
    // regions as lowercase `acgtn`. The mask carries no information
    // used downstream, and the PSP writer rejects anything outside
    // A/C/G/T/N; uppercase here so every consumer sees canonical bases.
    Ok(bytes[start_idx..end_idx]
        .iter()
        .map(|b| b.to_ascii_uppercase())
        .collect())
}

// ---------------------------------------------------------------------
// SingleChromRefFetcher — per-worker, owns one contig's bytes.
// ---------------------------------------------------------------------

/// `RefSeqFetcher` bound to **one** `chrom_id`. Each cohort
/// `var-calling` per-chromosome worker constructs its own instance,
/// uses it through the run, and drops it when the worker finishes —
/// at which point the contig's `Vec<u8>` is freed.
///
/// No shared mutable state: each worker owns its `SingleChromRefFetcher`
/// outright. No `Arc`-of-cache shared across threads, no `Mutex`,
/// no lease ref-counting. This is the simplest shape that delivers
/// per-chrom eviction under per-chrom outer parallelism — and the
/// shape with zero possibility of cross-thread interference.
///
/// Construction reads the contig once via a one-shot
/// `noodles_fasta::io::IndexedReader::query` (each instance opens its
/// own `File` — noodles supports concurrent read-only indexed readers
/// over the same FASTA) and uppercases the bytes in place. Subsequent
/// `fetch` calls slice into the owned `Vec<u8>` and copy out — no
/// per-fetch uppercase pass, no `Repository` cache.
///
/// Memory shape under `var-calling` parallelism: peak resident is
/// `min(threads, n_chroms) × max_chrom_size`, not `Σ contig.length`.
/// Phase B of the
/// [`reference_fasta_streaming`](../../../doc/devel/implementation_plans/reference_fasta_streaming.md)
/// plan.
pub struct SingleChromRefFetcher {
    /// The single chrom_id this fetcher is bound to. `fetch` calls
    /// with any other chrom_id return `InvalidInput` — catches
    /// "wrong worker fetcher passed to wrong chrom" bugs early.
    chrom_id: u32,
    /// Contig name for error messages only; not used in the fetch
    /// hot path.
    contig_name: String,
    /// Uppercased contig bases. Owned by-value; freed when this
    /// fetcher is dropped.
    bytes: Vec<u8>,
}

impl SingleChromRefFetcher {
    /// Load the contig identified by `(chrom_id, contig_name)` from
    /// the FASTA at `fasta_path`. Opens its own indexed reader,
    /// queries the whole contig, uppercases in place, drops the
    /// reader. Returns ready-to-fetch.
    pub fn new(fasta_path: &Path, chrom_id: u32, contig_name: String) -> io::Result<Self> {
        let mut reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(fasta_path)?;
        let region = Region::new(contig_name.as_bytes().to_vec(), ..);
        let record = reader.query(&region)?;
        let mut bytes: Vec<u8> = record.sequence().as_ref().to_vec();
        bytes.make_ascii_uppercase();
        Ok(Self {
            chrom_id,
            contig_name,
            bytes,
        })
    }

    /// The chrom_id this fetcher is bound to. Test/diagnostic only.
    #[cfg(test)]
    pub fn chrom_id(&self) -> u32 {
        self.chrom_id
    }

    /// Owned contig length in bytes. Test/diagnostic only.
    #[cfg(test)]
    pub fn contig_len(&self) -> usize {
        self.bytes.len()
    }
}

impl RefSeqFetcher for SingleChromRefFetcher {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        if chrom_id != self.chrom_id {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "fetch on chrom_id {chrom_id} but this fetcher is bound to chrom_id {} ({})",
                    self.chrom_id, self.contig_name
                ),
            ));
        }
        if start_1based == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "start_1based must be >= 1",
            ));
        }
        let start_idx = (start_1based - 1) as usize;
        let end_idx = start_idx
            .checked_add(length as usize)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"))?;
        if end_idx > self.bytes.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "fetch [{}, {}) past contig {} length {}",
                    start_1based,
                    start_1based + length,
                    self.contig_name,
                    self.bytes.len()
                ),
            ));
        }
        // Bytes are pre-uppercased at construction time, so this is
        // a straight slice copy (one of the L4 wins in the
        // 2026-05-20 perf review).
        Ok(self.bytes[start_idx..end_idx].to_vec())
    }
}

// ---------------------------------------------------------------------
// StreamingChromRefFetcher — per-worker, sliding-buffer reader.
// ---------------------------------------------------------------------

/// `RefSeqFetcher` bound to one `chrom_id` that reads its bytes from
/// disk through a small sliding buffer. Designed for the cohort
/// `var-calling` per-chromosome workers, where access to the contig
/// is sequential and monotonically non-decreasing in position:
///
/// - DUST mask construction does **one** forward sequential pass over
///   the contig (see [`Self::bases`]).
/// - `PerGroupMerger` fetches small windows (≤ `--var-group-max-span`,
///   default 10 KB) at non-decreasing `start` positions.
///
/// Trade-off compared to [`SingleChromRefFetcher`] (Phase B):
///
/// - **Memory:** `STREAMING_REF_BUFFER_BYTES` per worker (~1 MB)
///   instead of one full contig (~91 MB on tomato ch01, ~250 MB on
///   human chr1). Under per-chrom outer parallelism at T workers,
///   peak fetcher contribution is `T × 1 MB`, not `T × max_chrom_size`.
/// - **I/O:** a handful of `read(2)` syscalls per contig instead of
///   one big read. With a warm page cache this is microseconds; with
///   a cold cache the kernel readahead does the same work either way.
/// - **Contract:** the same as `SingleChromRefFetcher` — `fetch` only
///   answers for the bound `chrom_id`. The fetcher tolerates
///   non-monotonic access (refills on any miss), but is sized for the
///   monotonic-forward case; backwards seeks degrade to "extra refill
///   per backward jump".
///
/// Phase C of the
/// [`reference_fasta_streaming`](../../../doc/devel/implementation_plans/reference_fasta_streaming.md)
/// plan.
pub struct StreamingChromRefFetcher {
    /// Bound chromosome id. `fetch` calls with any other id are
    /// rejected with `InvalidInput`.
    chrom_id: u32,
    /// Contig name for error messages only.
    contig_name: String,
    /// `.fai` record for this contig — sequence start offset + line
    /// layout. Captured at construction so we don't re-parse the
    /// `.fai` on every refill.
    fai: ContigFai,
    /// All mutable streamer state. `Mutex` because the trait alias
    /// `SharedRefFetcher = Arc<dyn RefSeqFetcher + Send + Sync>`
    /// requires `Sync`, but each per-chrom worker owns its fetcher
    /// outright, so contention is zero by construction.
    inner: Mutex<StreamState>,
}

/// Pared-down `.fai` record — just the fields the streamer needs.
/// Avoids carrying around a `noodles::fai::Record` (which is harder
/// to construct from raw values in the in-memory test fixture).
#[derive(Debug, Clone, Copy)]
struct ContigFai {
    /// Byte offset of the first base in the FASTA file.
    seq_offset: u64,
    /// Total bases in the contig.
    length: u32,
    /// Bases per text line (excluding the trailing newline).
    line_bases: u32,
    /// Line width including the trailing newline.
    line_width: u32,
}

impl ContigFai {
    /// File-byte offset of the (1-based) base coordinate `b`.
    fn base_to_file_offset(&self, b: u32) -> u64 {
        let zero_based = (b - 1) as u64;
        let line_idx = zero_based / self.line_bases as u64;
        let in_line = zero_based % self.line_bases as u64;
        self.seq_offset + line_idx * self.line_width as u64 + in_line
    }
}

/// Source of FASTA bytes for the streamer. Production uses a real
/// `File`; tests use an in-memory `Vec<u8>` (`Cursor`) so the suite
/// stays fast and hermetic.
enum Source {
    File(File),
    #[cfg_attr(not(test), allow(dead_code))]
    Memory(io::Cursor<Vec<u8>>),
}

impl Source {
    fn seek_to(&mut self, offset: u64) -> io::Result<()> {
        match self {
            Source::File(f) => {
                f.seek(SeekFrom::Start(offset))?;
            }
            Source::Memory(c) => {
                c.seek(SeekFrom::Start(offset))?;
            }
        }
        Ok(())
    }
    fn read_chunk(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Source::File(f) => f.read(buf),
            Source::Memory(c) => c.read(buf),
        }
    }
}

struct StreamState {
    source: Source,
    /// Uppercased, newline-stripped bases currently resident.
    /// Length ≤ `STREAMING_REF_BUFFER_BYTES`.
    buf: Vec<u8>,
    /// 1-based contig coordinate of `buf[0]`. Bases held in `buf`
    /// cover contig positions `[buf_start_base, buf_start_base + buf.len())`.
    /// Zero when no buffer has been filled yet.
    buf_start_base: u32,
}

impl StreamingChromRefFetcher {
    /// Open the FASTA at `fasta_path`, parse the sibling `.fai`, look
    /// up `contig_name`, and return a fetcher ready to serve `fetch`
    /// calls for `chrom_id`. The contig's bytes are **not** loaded
    /// here — that happens lazily on the first fetch.
    pub fn new(fasta_path: &Path, chrom_id: u32, contig_name: String) -> io::Result<Self> {
        let mut fai_pathbuf: OsString = fasta_path.as_os_str().to_os_string();
        fai_pathbuf.push(".fai");
        let fai_path = PathBuf::from(fai_pathbuf);
        let index = fai::fs::read(&fai_path)?;
        let record = index
            .as_ref()
            .iter()
            .find(|r| AsRef::<[u8]>::as_ref(r.name()) == contig_name.as_bytes())
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("contig {contig_name} not in FASTA index"),
                )
            })?;
        let length = u32::try_from(record.length()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("contig {} length {} exceeds u32::MAX", contig_name, record.length()),
            )
        })?;
        let line_bases = u32::try_from(record.line_bases()).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, ".fai line_bases exceeds u32::MAX")
        })?;
        let line_width = u32::try_from(record.line_width()).map_err(|_| {
            io::Error::new(io::ErrorKind::InvalidData, ".fai line_width exceeds u32::MAX")
        })?;
        let file = File::open(fasta_path)?;
        Ok(Self {
            chrom_id,
            contig_name,
            fai: ContigFai {
                seq_offset: record.offset(),
                length,
                line_bases,
                line_width,
            },
            inner: Mutex::new(StreamState {
                source: Source::File(file),
                buf: Vec::with_capacity(STREAMING_REF_BUFFER_BYTES),
                buf_start_base: 0,
            }),
        })
    }

    /// Test-only constructor: stream bytes out of an in-memory FASTA
    /// blob instead of a real file. The blob must be a single
    /// well-formed contig (one `>name\n` header line followed by
    /// `line_bases`-wrapped sequence lines, ending with the contig's
    /// final line). Used by the unit tests to keep the suite
    /// hermetic + fast.
    #[cfg(test)]
    fn from_memory(
        fasta_blob: Vec<u8>,
        chrom_id: u32,
        contig_name: String,
        seq_offset: u64,
        length: u32,
        line_bases: u32,
        line_width: u32,
    ) -> Self {
        Self {
            chrom_id,
            contig_name,
            fai: ContigFai {
                seq_offset,
                length,
                line_bases,
                line_width,
            },
            inner: Mutex::new(StreamState {
                source: Source::Memory(io::Cursor::new(fasta_blob)),
                buf: Vec::with_capacity(STREAMING_REF_BUFFER_BYTES),
                buf_start_base: 0,
            }),
        }
    }

    /// Bound chromosome id. Test-only.
    #[cfg(test)]
    pub fn chrom_id(&self) -> u32 {
        self.chrom_id
    }

    /// Forward sequential iterator over every uppercased base of the
    /// contig. Pulls bytes through the streamer's sliding buffer one
    /// at a time, refilling on demand. Used by the DUST mask
    /// construction in
    /// [`crate::var_calling::dust_filter::DustFilter::ensure_mask_for`].
    pub fn bases(&self) -> StreamingBaseIter<'_> {
        StreamingBaseIter {
            fetcher: self,
            next_base: 1,
            done: false,
        }
    }
}

impl RefSeqFetcher for StreamingChromRefFetcher {
    fn iter_bases<'a>(
        &'a self,
        chrom_id: u32,
        length: u32,
    ) -> Result<Box<dyn Iterator<Item = io::Result<u8>> + 'a>, io::Error> {
        if chrom_id != self.chrom_id {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "iter_bases on chrom_id {chrom_id} but this fetcher is bound to chrom_id {} ({})",
                    self.chrom_id, self.contig_name
                ),
            ));
        }
        if length != self.fai.length {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "iter_bases asked for {length} bases but contig {} has length {}",
                    self.contig_name, self.fai.length
                ),
            ));
        }
        // Override: stream through the sliding buffer instead of
        // materialising the contig.
        Ok(Box::new(self.bases()))
    }

    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        if chrom_id != self.chrom_id {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "fetch on chrom_id {chrom_id} but this fetcher is bound to chrom_id {} ({})",
                    self.chrom_id, self.contig_name
                ),
            ));
        }
        if start_1based == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "start_1based must be >= 1",
            ));
        }
        let end_1based = start_1based.checked_add(length).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow")
        })?;
        if end_1based.saturating_sub(1) > self.fai.length {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "fetch [{}, {}) past contig {} length {}",
                    start_1based, end_1based, self.contig_name, self.fai.length
                ),
            ));
        }

        let mut state = self
            .inner
            .lock()
            .expect("StreamingChromRefFetcher mutex poisoned");
        // Check if the requested window already lies inside the
        // current buffer; refill otherwise.
        let buf_covers = !state.buf.is_empty()
            && start_1based >= state.buf_start_base
            && (end_1based as u64)
                <= state.buf_start_base as u64 + state.buf.len() as u64;
        if !buf_covers {
            // Refill at least the requested span, capped at the
            // configured buffer size. If the request exceeds the
            // buffer size (rare — would mean a single fetch larger
            // than `STREAMING_REF_BUFFER_BYTES`), grow the buffer for
            // this one call. The capacity stays at the larger value
            // afterward, but typical access keeps it at 1 MB.
            let want = (length as usize).max(STREAMING_REF_BUFFER_BYTES);
            refill(&mut state, &self.fai, start_1based, want)?;
        }
        let start_idx = (start_1based - state.buf_start_base) as usize;
        let end_idx = start_idx + length as usize;
        Ok(state.buf[start_idx..end_idx].to_vec())
    }
}

/// Iterator over uppercased contig bases, in 1..=length order. Yields
/// `io::Result<u8>` so the caller can surface refill errors. Used by
/// the DUST mask construction; one logical pass per chrom.
pub struct StreamingBaseIter<'a> {
    fetcher: &'a StreamingChromRefFetcher,
    /// 1-based coordinate of the next base to yield.
    next_base: u32,
    /// Latches on the first error so subsequent `next()` calls return
    /// `None` (matches the convention for fused error-returning
    /// iterators).
    done: bool,
}

impl Iterator for StreamingBaseIter<'_> {
    type Item = io::Result<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done || self.next_base > self.fetcher.fai.length {
            return None;
        }
        let mut state = match self.fetcher.inner.lock() {
            Ok(s) => s,
            Err(_) => {
                self.done = true;
                return Some(Err(io::Error::other(
                    "StreamingChromRefFetcher mutex poisoned in bases()",
                )));
            }
        };
        let need_refill = state.buf.is_empty()
            || self.next_base < state.buf_start_base
            || (self.next_base as u64)
                >= state.buf_start_base as u64 + state.buf.len() as u64;
        if need_refill
            && let Err(e) = refill(
                &mut state,
                &self.fetcher.fai,
                self.next_base,
                STREAMING_REF_BUFFER_BYTES,
            )
        {
            self.done = true;
            return Some(Err(e));
        }
        let idx = (self.next_base - state.buf_start_base) as usize;
        let byte = state.buf[idx];
        self.next_base += 1;
        Some(Ok(byte))
    }
}

/// Seek to `target_base_1based` in the FASTA and fill `state.buf`
/// with up to `want` uppercased non-newline bases (capped at the
/// contig's remaining length). Resets `state.buf_start_base` to the
/// new origin.
fn refill(
    state: &mut StreamState,
    fai: &ContigFai,
    target_base_1based: u32,
    want: usize,
) -> io::Result<()> {
    let offset = fai.base_to_file_offset(target_base_1based);
    state.source.seek_to(offset)?;
    state.buf.clear();
    state.buf_start_base = target_base_1based;

    // Bound by what's left on the contig — never read past the last
    // base, even if the caller asked for more. Refill error reporting
    // is left to fetch / next() so the caller sees the right kind.
    let remaining = (fai.length - target_base_1based + 1) as usize;
    let target_len = want.min(remaining);

    // Read 64 KiB raw windows of FASTA bytes; strip `\n`/`\r` and
    // uppercase the rest into `state.buf` until we've harvested
    // `target_len` bases (or the source runs out, which is an error).
    let mut read_buf = [0u8; STREAMING_REF_FILE_READ_CHUNK];
    while state.buf.len() < target_len {
        let n = state.source.read_chunk(&mut read_buf)?;
        if n == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "FASTA ended with {} base(s) still expected at refill from base {} (target_len {})",
                    target_len - state.buf.len(),
                    target_base_1based,
                    target_len,
                ),
            ));
        }
        for &b in &read_buf[..n] {
            if state.buf.len() >= target_len {
                break;
            }
            if b == b'\n' || b == b'\r' {
                continue;
            }
            state.buf.push(b.to_ascii_uppercase());
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    use crate::per_sample_pileup::cram_files::{ContigSpec, build_fasta};
    use crate::per_sample_pileup::cram_input::{ContigEntry, ContigList};

    fn contig_list(entries: &[(&str, u64)]) -> ContigList {
        ContigList {
            entries: entries
                .iter()
                .map(|(n, l)| ContigEntry {
                    name: (*n).into(),
                    length: *l,
                    md5: None,
                })
                .collect(),
        }
    }

    #[test]
    fn fetch_returns_bases_at_correct_offset() {
        // build_fasta writes 'A' for every base, so the test is
        // about offsets / range handling, not byte values.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 100,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 100)])).expect("fetcher");

        let bases = fetcher.fetch(0, 5, 4).expect("fetch");
        assert_eq!(bases, b"AAAA");
        assert_eq!(bases.len(), 4);
    }

    #[test]
    fn fetch_uppercases_soft_masked_bases() {
        // Soft-masked FASTAs (Ensembl/Gencode default) encode repeat
        // regions as lowercase `acgtn`. The fetcher must canonicalise
        // to uppercase so the PSP writer (which rejects lowercase
        // bytes) never sees a soft-masked REF byte.
        use std::fs::File;
        use std::io::Write;

        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("soft.fa");
        let fai_path = dir.path().join("soft.fa.fai");

        let header = b">chr0\n";
        let seq = b"ACGTacgtNn";
        {
            let mut fa = File::create(&fasta_path).expect("fa");
            fa.write_all(header).expect("hdr");
            fa.write_all(seq).expect("seq");
            fa.write_all(b"\n").expect("nl");
        }
        {
            let mut fai = File::create(&fai_path).expect("fai");
            writeln!(
                fai,
                "chr0\t{}\t{}\t{}\t{}",
                seq.len(),
                header.len(),
                seq.len(),
                seq.len() + 1
            )
            .expect("fai");
        }

        let fetcher =
            ChromBoundaryRefFetcher::new(&fasta_path, contig_list(&[("chr0", seq.len() as u64)]))
                .expect("fetcher");

        let bases = fetcher.fetch(0, 1, seq.len() as u32).expect("fetch");
        assert_eq!(bases, b"ACGTACGTNN");
    }

    #[test]
    fn cache_evicts_previous_chromosome_on_chrom_change() {
        // The eviction invariant: cache size stays at 1 across
        // chrom changes, regardless of how many chromosomes the
        // walker has visited. Without `Repository::clear()` it
        // would grow unboundedly (a regression S6 was added to
        // prevent).
        let specs = vec![
            ContigSpec {
                name: "chr0".into(),
                length: 100,
            },
            ContigSpec {
                name: "chr1".into(),
                length: 100,
            },
        ];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 100), ("chr1", 100)]))
                .expect("fetcher");

        assert_eq!(
            fetcher.cached_contig_count(),
            0,
            "fresh fetcher should have an empty cache"
        );

        fetcher.fetch(0, 1, 4).expect("fetch chr0");
        assert_eq!(
            fetcher.cached_contig_count(),
            1,
            "chr0 must be cached after first fetch"
        );

        fetcher.fetch(1, 1, 4).expect("fetch chr1");
        assert_eq!(
            fetcher.cached_contig_count(),
            1,
            "chr0 must have been evicted when fetcher saw chr1; \
             cache should hold the single current contig only"
        );

        // A second fetch on the same chromosome must reuse the
        // cached entry — no new eviction, no growth.
        fetcher.fetch(1, 5, 4).expect("re-fetch chr1");
        assert_eq!(fetcher.cached_contig_count(), 1);

        // Going back to chr0 also re-loads (chr0 was evicted) and
        // does not double-cache.
        fetcher.fetch(0, 1, 4).expect("re-fetch chr0");
        assert_eq!(fetcher.cached_contig_count(), 1);
    }

    #[test]
    fn fetch_past_contig_end_returns_unexpected_eof() {
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 10)])).expect("fetcher");

        let err = fetcher.fetch(0, 8, 5).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn unknown_chrom_id_returns_invalid_input() {
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 10)])).expect("fetcher");

        let err = fetcher.fetch(99, 1, 4).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn start_1based_zero_is_rejected() {
        // 1-based coordinates: 0 is a contract violation, not a
        // valid offset.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 10)])).expect("fetcher");

        let err = fetcher.fetch(0, 0, 1).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn contig_name_in_fasta_must_match_contig_list() {
        // contig_list says the chrom is "chr0", but the FASTA on
        // disk only knows "different_name". The Repository surfaces
        // the missing contig as an `io::Error`; the exact kind is
        // noodles' choice (currently `InvalidInput` from the
        // indexed-reader adapter), so we assert only that fetch
        // errors out — not which kind.
        let specs = vec![ContigSpec {
            name: "different_name".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher =
            ChromBoundaryRefFetcher::new(&path, contig_list(&[("chr0", 10)])).expect("fetcher");

        fetcher.fetch(0, 1, 4).expect_err("must fail");
    }

    // ---------------------------------------------------------------
    // SingleChromRefFetcher tests
    // ---------------------------------------------------------------

    #[test]
    fn single_chrom_fetcher_loads_contig_and_serves_bases() {
        // Positive path: build a 50-base 'AAA…' FASTA, construct a
        // fetcher bound to chrom_id 0, fetch a window. The bytes are
        // uppercased at construction so the slice is verbatim.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 50,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher = SingleChromRefFetcher::new(&path, 0, "chr0".into()).expect("fetcher");
        assert_eq!(fetcher.chrom_id(), 0);
        assert_eq!(fetcher.contig_len(), 50);

        let bytes = fetcher.fetch(0, 5, 4).expect("fetch");
        assert_eq!(bytes, b"AAAA");
    }

    #[test]
    fn single_chrom_fetcher_uppercases_soft_masked_bases() {
        // Soft-masked FASTAs encode repeats as lowercase `acgtn`. The
        // constructor must uppercase so `fetch` returns canonical
        // SAM-spec bytes — matches the long-standing `SyncRefFetcher`
        // contract.
        use std::fs::File;
        use std::io::Write;

        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("soft.fa");
        let fai_path = dir.path().join("soft.fa.fai");

        let header = b">chr0\n";
        let seq = b"ACGTacgtNn";
        {
            let mut fa = File::create(&fasta_path).expect("fa");
            fa.write_all(header).expect("hdr");
            fa.write_all(seq).expect("seq");
            fa.write_all(b"\n").expect("nl");
        }
        {
            let mut fai = File::create(&fai_path).expect("fai");
            writeln!(
                fai,
                "chr0\t{}\t{}\t{}\t{}",
                seq.len(),
                header.len(),
                seq.len(),
                seq.len() + 1
            )
            .expect("fai");
        }

        let fetcher = SingleChromRefFetcher::new(&fasta_path, 0, "chr0".into()).expect("fetcher");
        let bytes = fetcher.fetch(0, 1, seq.len() as u32).expect("fetch");
        assert_eq!(bytes, b"ACGTACGTNN");
    }

    #[test]
    fn single_chrom_fetch_on_wrong_chrom_id_errors() {
        // Contract: this fetcher answers only for the bound
        // `chrom_id`. Passing any other id is an `InvalidInput`
        // — catches "wrong worker passed wrong fetcher" bugs early.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher = SingleChromRefFetcher::new(&path, 0, "chr0".into()).expect("fetcher");

        let err = fetcher.fetch(7, 1, 4).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn single_chrom_fetch_past_contig_end_returns_unexpected_eof() {
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher = SingleChromRefFetcher::new(&path, 0, "chr0".into()).expect("fetcher");

        let err = fetcher.fetch(0, 8, 5).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn single_chrom_fetch_start_zero_is_rejected() {
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher = SingleChromRefFetcher::new(&path, 0, "chr0".into()).expect("fetcher");

        let err = fetcher.fetch(0, 0, 1).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn single_chrom_construct_missing_contig_errors() {
        // The supplied contig name is not in the FASTA. Construction
        // must fail (rather than build a zero-length fetcher).
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");

        // `expect_err` would require Debug on the Ok variant; just
        // assert via `is_err()`.
        assert!(SingleChromRefFetcher::new(&path, 7, "chr_missing".into()).is_err());
    }

    #[test]
    fn single_chrom_fetchers_are_independent_across_threads() {
        // Two threads, each constructs its own SingleChromRefFetcher
        // on a distinct contig from the same FASTA. No shared state,
        // no Arc, no Mutex — sanity that the "each worker owns its
        // own" design composes under real concurrency.
        let specs = vec![
            ContigSpec {
                name: "chr0".into(),
                length: 32,
            },
            ContigSpec {
                name: "chr1".into(),
                length: 32,
            },
        ];
        let (dir, path) = build_fasta(&specs).expect("fasta");
        let path = std::sync::Arc::new(path);
        let _dir_keep = dir; // tempdir lifetime

        let mut handles = Vec::new();
        for (cid, name) in [(0u32, "chr0"), (1u32, "chr1")] {
            let path = std::sync::Arc::clone(&path);
            handles.push(std::thread::spawn(move || {
                let fetcher = SingleChromRefFetcher::new(&path, cid, name.into()).expect("fetcher");
                let bytes = fetcher.fetch(cid, 1, 4).expect("fetch");
                assert_eq!(bytes, b"AAAA");
            }));
        }
        for h in handles {
            h.join().expect("thread");
        }
    }

    // ---------------------------------------------------------------
    // StreamingChromRefFetcher tests
    // ---------------------------------------------------------------

    /// Build an in-memory FASTA blob with a single contig wrapped at
    /// `line_bases` chars per line, plus the `.fai`-equivalent
    /// numbers (offset, length, line_bases, line_width). Used to
    /// drive [`StreamingChromRefFetcher::from_memory`] without
    /// touching disk.
    fn build_memory_fasta(
        name: &str,
        seq: &[u8],
        line_bases: u32,
    ) -> (Vec<u8>, u64 /* seq_offset */, u32 /* line_width */) {
        let header = format!(">{name}\n");
        let mut blob = header.as_bytes().to_vec();
        let seq_offset = blob.len() as u64;
        let lb = line_bases as usize;
        let mut i = 0;
        while i < seq.len() {
            let end = (i + lb).min(seq.len());
            blob.extend_from_slice(&seq[i..end]);
            blob.push(b'\n');
            i = end;
        }
        let line_width = line_bases + 1;
        (blob, seq_offset, line_width)
    }

    #[test]
    fn streaming_fetcher_returns_bytes_for_bound_chrom() {
        // 50-base contig wrapped at 50/line (single line). Fetch a
        // window in the middle and assert the right bytes come back.
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAN";
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", seq, 50);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            50,
            line_width,
        );
        assert_eq!(fetcher.chrom_id(), 0);

        let bytes = fetcher.fetch(0, 5, 4).expect("fetch");
        // 1-based start=5, length=4 → bytes at indices 4..8.
        assert_eq!(bytes, &seq[4..8]);
    }

    #[test]
    fn streaming_fetcher_serves_line_wrapped_contig() {
        // 250-base contig wrapped at 60 chars/line (Ensembl default).
        // Multiple fetches across line boundaries; every one must
        // return the correct bases without picking up newlines.
        let seq: Vec<u8> = (0..250).map(|i| b"ACGT"[i % 4]).collect();
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", &seq, 60);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            60,
            line_width,
        );

        // Spans across line boundary.
        let bytes = fetcher.fetch(0, 58, 5).expect("fetch");
        assert_eq!(bytes, seq[57..62]);
        // Wholly inside one line.
        let bytes = fetcher.fetch(0, 70, 10).expect("fetch");
        assert_eq!(bytes, seq[69..79]);
    }

    #[test]
    fn streaming_fetcher_uppercases_soft_masked_bases() {
        // Lowercase input → uppercased output. Same contract as the
        // other fetchers.
        let seq = b"ACGTacgtNnacgtACGT".to_vec();
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", &seq, 18);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            18,
            line_width,
        );
        let bytes = fetcher.fetch(0, 1, seq.len() as u32).expect("fetch");
        assert_eq!(bytes, b"ACGTACGTNNACGTACGT");
    }

    #[test]
    fn streaming_fetcher_buffer_refill_on_forward_jump() {
        // Two fetches separated by a gap larger than the buffer.
        // Validates the refill path: second fetch must re-seek and
        // return the right bytes. We use a small synthetic buffer
        // size by jumping past `STREAMING_REF_BUFFER_BYTES` worth of
        // bases; cheaper to assert correctness with two short
        // fetches at distant positions on a 2 MB contig.
        let seq: Vec<u8> = (0..2_000_000_u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", &seq, 80);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            80,
            line_width,
        );
        // Fetch near the start.
        let a = fetcher.fetch(0, 100, 16).expect("fetch a");
        assert_eq!(a, seq[99..115]);
        // Then jump well past the 1 MB buffer.
        let b = fetcher.fetch(0, 1_500_000, 16).expect("fetch b");
        assert_eq!(b, seq[1_499_999..1_500_015]);
    }

    #[test]
    fn streaming_fetcher_buffer_refill_on_backward_jump() {
        // Backward jump after a forward fetch. The streamer should
        // detect the buffer miss and refill from the new (earlier)
        // position. Documents the slower path; production access
        // never hits this (PerGroupMerger is monotonic), but tests
        // pin the behaviour.
        let seq: Vec<u8> = (0..2_000_000_u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", &seq, 80);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            80,
            line_width,
        );
        let _ = fetcher.fetch(0, 1_500_000, 16).expect("fetch forward");
        let back = fetcher.fetch(0, 100, 16).expect("fetch back");
        assert_eq!(back, seq[99..115]);
    }

    #[test]
    fn streaming_fetcher_fetch_past_contig_end_returns_unexpected_eof() {
        let seq = b"ACGTACGTAC".to_vec();
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", &seq, 10);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            10,
            line_width,
        );
        let err = fetcher.fetch(0, 8, 5).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn streaming_fetcher_fetch_on_wrong_chrom_id_errors() {
        let seq = b"ACGT".to_vec();
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", &seq, 4);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            4,
            line_width,
        );
        let err = fetcher.fetch(7, 1, 4).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn streaming_fetcher_fetch_start_zero_is_rejected() {
        let seq = b"ACGT".to_vec();
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", &seq, 4);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            4,
            line_width,
        );
        let err = fetcher.fetch(0, 0, 1).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn streaming_fetcher_construct_missing_contig_errors() {
        // Real-file path: build a FASTA with one contig and try to
        // open a fetcher for a contig that's not in the .fai.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        assert!(StreamingChromRefFetcher::new(&path, 7, "chr_missing".into()).is_err());
    }

    #[test]
    fn streaming_fetcher_bases_iterator_yields_uppercased_contig_in_order() {
        // DUST mask path: bases() iterator must yield every contig
        // base in 1..=length order, uppercased, exactly once. No
        // newlines, no headers.
        let seq: Vec<u8> = (0..200).map(|i| b"acgtACGTN"[i % 9]).collect();
        let mut expected: Vec<u8> = seq.clone();
        expected.make_ascii_uppercase();
        let (blob, seq_offset, line_width) = build_memory_fasta("chr0", &seq, 50);
        let fetcher = StreamingChromRefFetcher::from_memory(
            blob,
            0,
            "chr0".into(),
            seq_offset,
            seq.len() as u32,
            50,
            line_width,
        );
        let got: Vec<u8> = fetcher.bases().map(|b| b.expect("base")).collect();
        assert_eq!(got, expected);
    }

    #[test]
    fn streaming_fetcher_against_real_file() {
        // Smoke test: the production constructor (real File + real
        // .fai) returns the same bytes as the existing
        // SingleChromRefFetcher on the same fixture.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 64,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let streaming = StreamingChromRefFetcher::new(&path, 0, "chr0".into()).expect("fetcher");
        let single = SingleChromRefFetcher::new(&path, 0, "chr0".into()).expect("fetcher");

        for (start, length) in [(1u32, 16u32), (8, 32), (50, 14), (64, 1)] {
            let a = streaming.fetch(0, start, length).expect("streaming");
            let b = single.fetch(0, start, length).expect("single");
            assert_eq!(a, b, "mismatch at start={start} length={length}");
        }
    }
}
