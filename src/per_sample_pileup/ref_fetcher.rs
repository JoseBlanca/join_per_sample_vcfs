//! Reference-FASTA fetchers. One concrete impl —
//! [`StreamingChromRefFetcher`] — backed by a sliding buffer that
//! reads from disk on demand. Two trait-bridging adapters let it
//! drop into both the new single-contig [`ChromRefFetcher`] API
//! (used by cohort var-calling per-chrom workers) and the legacy
//! multi-chrom `RefSeqFetcher` API (used by the Stage 1 walker, BAQ,
//! and `var_calling_from_bam`).
//!
//! ## Why one fetcher + two adapters, not one fetcher per access shape
//!
//! Earlier iterations had a separate fetcher type per access pattern
//! (chrom-boundary-evicting for the walker, accumulating Sync cache
//! for BAQ, single-chrom owned-bytes for cohort workers). The
//! `unified_chrom_ref_fetcher` plan
//! (`doc/devel/implementation_plans/unified_chrom_ref_fetcher.md`)
//! collapsed that menu: the access pattern is encoded by *which
//! adapter* a caller picks, not by which fetcher *type* exists. The
//! sliding-buffer streamer underneath is the same in every case.
//!
//! ## Layout
//!
//! - [`StreamingChromRefFetcher`] — the workhorse. Bound to one
//!   contig at construction; serves `fetch(start, length)` and
//!   `iter_bases()` through a 1 MB sliding buffer. Implements both
//!   the new [`ChromRefFetcher`] trait and (via the adapters below)
//!   the legacy `RefSeqFetcher` trait.
//! - [`WalkerLegacyAdapter`] — wraps a swappable
//!   `StreamingChromRefFetcher` and exposes the legacy multi-chrom
//!   `RefSeqFetcher` API for consumers that traverse multiple
//!   contigs (Stage 1 walker, Stage 1 BAQ shared instance,
//!   `var_calling_from_bam`). Rebuilds the inner streamer on chrom
//!   transition. `Sync` via internal `Mutex`.
//!
//! Both adapters are transitional in spirit (the plan's step 3
//! would migrate every legacy consumer to the new trait directly),
//! but they're load-bearing today and would only go away if every
//! multi-chrom consumer were restructured to per-chrom processing.
//! That restructure is non-trivial; we ship the unified shape and
//! revisit if the adapter layer becomes a real friction point.

use std::cell::RefCell;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

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
/// Per worker the resident footprint is `STREAMING_REF_BUFFER_BYTES`
/// (~1 MB) regardless of contig size. Under per-chrom outer
/// parallelism at T workers the peak fetcher contribution is
/// `T × 1 MB`, not `T × max_chrom_size`.
///
/// I/O cost: a handful of `read(2)` syscalls per contig instead of
/// one big read. With a warm page cache this is microseconds; with a
/// cold cache the kernel readahead does the same work either way.
///
/// **Contract.** `fetch` only answers for the bound `chrom_id`. The
/// fetcher tolerates non-monotonic access today (refills on any
/// miss), but is sized for the monotonic-forward case; a backwards
/// seek beyond the sliding buffer triggers an extra refill (slow but
/// correct). The follow-up
/// [`unified_chrom_ref_fetcher`](../../../doc/devel/implementation_plans/unified_chrom_ref_fetcher.md)
/// plan replaces this with an explicit `OutOfPattern` error.
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
    /// All mutable streamer state. `RefCell` (not `Mutex`) because
    /// each per-chrom rayon worker owns its fetcher outright — no
    /// cross-thread sharing — so we don't pay atomics on every
    /// `ChromRefBaseIter::next`. The cohort `SharedRefFetcher` alias
    /// is `Arc<dyn ChromRefFetcher + Send>` (not `Send + Sync`),
    /// matching the per-worker ownership invariant documented at the
    /// cohort driver's worker construction.
    inner: RefCell<StreamState>,
}

// M14 (2026-05-23 code review): explicit compile-time assertion of
// `Send + !Sync`. Today this holds because `RefCell<StreamState>`
// auto-derives both bounds correctly — but the cohort concurrency
// model documented at the `inner` field is load-bearing, and a
// future field addition (e.g. an Arc<something> or a Mutex) could
// silently flip the bounds and let the cohort driver share the
// fetcher across threads via accident. The Send check below is a
// hard fence; the !Sync requirement is enforced by `RefCell` itself
// (it's `!Sync` by definition), so any field change that breaks the
// !Sync property has to first remove the `RefCell`, which the
// reviewer would catch.
const _: () = {
    const fn assert_send<T: Send>() {}
    assert_send::<StreamingChromRefFetcher>();
};

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
    /// Validate a parsed `.fai` record. The noodles parser accepts
    /// values that would crash the byte-offset math:
    ///
    /// - `line_bases = 0` → `base_to_file_offset` divides by zero.
    /// - `line_width < line_bases` → wrong file offsets (the trailing
    ///   newline can't have negative width).
    /// - `line_width = 0` (with non-zero `line_bases`) is rejected by
    ///   the first check above transitively.
    ///
    /// B1 of the 2026-05-23 code review: the `--reference` path is
    /// attacker-influenced; this validation guards every constructor
    /// from a panic-on-construction surface.
    fn validate(&self, contig_name: &str) -> Result<(), io::Error> {
        if self.line_bases == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "malformed .fai for contig {contig_name}: line_bases = 0 \
                     (would divide-by-zero in offset arithmetic)"
                ),
            ));
        }
        if self.line_width < self.line_bases {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "malformed .fai for contig {contig_name}: line_width ({}) \
                     < line_bases ({}) — line_width must include the trailing newline",
                    self.line_width, self.line_bases,
                ),
            ));
        }
        Ok(())
    }

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
///
/// M8 (2026-05-23 code review): the `Memory` variant is
/// `#[cfg(test)]`-gated rather than always-present, so production
/// builds carry no discriminant overhead for the test-only path
/// and the `Source` enum is a single-variant production enum
/// (which LLVM optimises out entirely).
enum Source {
    File(File),
    #[cfg(test)]
    Memory(io::Cursor<Vec<u8>>),
}

impl Source {
    fn seek_to(&mut self, offset: u64) -> io::Result<()> {
        match self {
            Source::File(f) => {
                f.seek(SeekFrom::Start(offset))?;
            }
            #[cfg(test)]
            Source::Memory(c) => {
                c.seek(SeekFrom::Start(offset))?;
            }
        }
        Ok(())
    }
    fn read_chunk(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            Source::File(f) => f.read(buf),
            #[cfg(test)]
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
                format!(
                    "contig {} length {} exceeds u32::MAX",
                    contig_name,
                    record.length()
                ),
            )
        })?;
        let line_bases = u32::try_from(record.line_bases()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                ".fai line_bases exceeds u32::MAX",
            )
        })?;
        let line_width = u32::try_from(record.line_width()).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                ".fai line_width exceeds u32::MAX",
            )
        })?;
        let fai = ContigFai {
            seq_offset: record.offset(),
            length,
            line_bases,
            line_width,
        };
        // B1: reject malformed .fai (line_bases = 0, etc.) before
        // any fetch can panic on the offset arithmetic.
        fai.validate(&contig_name)?;
        let file = File::open(fasta_path)?;
        Ok(Self {
            chrom_id,
            contig_name,
            fai,
            inner: RefCell::new(StreamState {
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
            inner: RefCell::new(StreamState {
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
        let end_1based = start_1based
            .checked_add(length)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"))?;
        if end_1based.saturating_sub(1) > self.fai.length {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "fetch [{}, {}) past contig {} length {}",
                    start_1based, end_1based, self.contig_name, self.fai.length
                ),
            ));
        }

        let mut state = self.inner.borrow_mut();
        // Check if the requested window already lies inside the
        // current buffer; refill otherwise.
        let buf_covers = !state.buf.is_empty()
            && start_1based >= state.buf_start_base
            && (end_1based as u64) <= state.buf_start_base as u64 + state.buf.len() as u64;
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
        let mut state = self.fetcher.inner.borrow_mut();
        let need_refill = state.buf.is_empty()
            || self.next_base < state.buf_start_base
            || (self.next_base as u64) >= state.buf_start_base as u64 + state.buf.len() as u64;
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
            // M26 (2026-05-23 code review): canonicalise to
            // `A`/`C`/`G`/`T`/`N` so the trait contract holds for
            // IUPAC-ambiguity and soft-masked input bytes.
            state.buf.push(canonicalise(b));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------
// New unified ChromRefFetcher API (in-flight migration).
// ---------------------------------------------------------------------
//
// The trait below is the target shape for the
// `unified_chrom_ref_fetcher` plan. It coexists with the legacy
// `RefSeqFetcher` trait during migration; consumers move one at a
// time. Once every consumer is on the new trait, the old one (and
// the old fetcher types that only implement it) get deleted.
//
// Differences from the legacy trait:
//
// - **No `chrom_id` parameter.** The fetcher is bound to one contig
//   at construction; passing `chrom_id` to every call was redundant
//   and validated-but-never-useful (single-chrom impls returned
//   `InvalidInput` on mismatch).
// - **Typed error.** `ChromRefFetchError` distinguishes "I/O
//   failure", "out of contig bounds", and "out of access pattern"
//   so consumers can route them differently. The legacy trait uses
//   `io::Error` for everything.
// - **`OutOfPattern` is the load-bearing contract.** A streaming
//   impl errors when a `fetch` start lies before the current sliding
//   buffer's origin. The contract is "monotonic non-decreasing
//   `start` across `fetch` calls within a phase"; `iter_bases()`
//   defines its own phase (resets the buffer state on construction
//   and on Drop). If the error fires in real code, that consumer
//   needs a random-access impl rather than a streaming one.
// - **`length()`** is exposed on the trait so consumers can use
//   the contig's declared length without threading it through
//   constructor calls.

/// Why a [`ChromRefFetcher`] call failed.
// M2 (2026-05-23 code review): #[non_exhaustive] so adding a future
// variant (e.g. an InvalidIndex variant for B1's validation if it
// graduates from io::Error to a typed variant) doesn't break callers.
// The enum lives behind a `pub` surface but the crate is unpublished;
// in-tree matchers must accept the marker today so they pay the
// add-a-variant cost upfront.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ChromRefFetchError {
    /// The requested `[start, start+length)` window exceeds the
    /// bound contig's total length.
    #[error("fetch [{start}, {end}) past contig {contig_name} length {contig_length}")]
    OutOfBounds {
        contig_name: String,
        contig_length: u32,
        start: u32,
        end: u32,
    },
    /// `start_1based` was 0 — 1-based coordinate contract violated.
    #[error("start_1based must be >= 1")]
    InvalidStart,
    /// Streaming impl only. The requested `start_1based` lies
    /// before the current sliding buffer's origin; the consumer's
    /// access pattern is not monotonic non-decreasing. Production
    /// callers (DUST + PerGroupMerger) satisfy the contract by
    /// construction; this error firing in real code means a
    /// non-monotonic consumer that needs a random-access fetcher
    /// impl with the same trait.
    #[error(
        "fetch at base {requested_start} lies before the streamer's current buffer (origin = {buffer_origin}). \
         StreamingChromRefFetcher requires monotonic non-decreasing access"
    )]
    OutOfPattern {
        requested_start: u32,
        buffer_origin: u32,
    },
    /// Underlying I/O failure: file open, seek, read, or `.fai`
    /// parse.
    #[error("I/O failure on contig {contig_name}: {source}")]
    Io {
        contig_name: String,
        #[source]
        source: io::Error,
    },
}

/// Sealed-marker module for [`ChromRefFetcher`]. Visibility is
/// `pub(crate)` so in-tree modules and tests can add new
/// `ChromRefFetcher` impls (each one must also `impl
/// sealed::Sealed`); external crates cannot. See M9 in the
/// 2026-05-23 ref_fetcher code review.
pub(crate) mod sealed {
    pub trait Sealed {}
}

/// Reference-FASTA fetcher bound to a single contig at construction
/// time. Implementations differ in the access patterns they accept:
///
/// - [`StreamingChromRefFetcher`] (the only impl today) requires
///   monotonic non-decreasing `start` across `fetch` calls within a
///   phase. Backward jumps beyond the current sliding buffer return
///   [`ChromRefFetchError::OutOfPattern`]. `iter_bases()` defines
///   its own phase: the buffer is reset on iter construction and on
///   Drop, so a `fetch` after an `iter_bases` walk starts fresh.
/// - A future random-access impl (built only if a real consumer
///   needs it) would accept any access pattern with a different
///   memory cost. Same trait, drop-in.
///
/// Consumers should not interleave `iter_bases()` and `fetch()` on
/// the same fetcher within a phase — use one or the other. Our
/// pipeline runs `iter_bases()` (DUST mask construction) once at
/// chrom load, then switches to `fetch()` (PerGroupMerger window
/// reads).
///
/// M9 (2026-05-23 code review): the trait is **sealed** —
/// `sealed::Sealed` is `pub(crate)` so in-tree modules and tests can
/// add new impls but external crates cannot. The trait carries a
/// non-trivial cross-call contract (`iter_bases` Drop-reset; the
/// monotonic-forward `fetch` invariant) that a downstream impl
/// could violate silently, and the API is still evolving (the
/// `fetch_into` default impl was added in a perf-review fix).
/// Sealing keeps both moving parts in-tree.
pub trait ChromRefFetcher: sealed::Sealed {
    /// Total bases in the bound contig.
    fn length(&self) -> u32;

    /// Return `length` uppercased bases starting at 1-based
    /// `start`. Bytes are SAM-spec canonical (`A`/`C`/`G`/`T`/`N`)
    /// regardless of how the FASTA is masked on disk.
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::InvalidStart`] if `start_1based == 0`.
    /// - [`ChromRefFetchError::OutOfBounds`] if
    ///   `start_1based + length - 1` exceeds the contig length.
    /// - [`ChromRefFetchError::OutOfPattern`] if `start_1based` is
    ///   before the current sliding-buffer window (impls that enforce
    ///   monotonic-forward access).
    /// - [`ChromRefFetchError::Io`] if the underlying FASTA read fails.
    fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError>;

    /// Forward sequential iterator over every uppercased base of
    /// the bound contig in 1..=length order. Used by DUST's mask
    /// construction. Yields per-byte `Result<u8, ChromRefFetchError>`
    /// so refill / I/O failures propagate without panicking the
    /// iterator early.
    ///
    /// Implementation contract: constructing this iterator resets
    /// the fetcher's internal sliding buffer; dropping it resets
    /// the buffer again. This makes `iter_bases()` followed by
    /// `fetch()` from an arbitrary `start` legal (the next `fetch`
    /// starts a fresh buffer phase).
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::Io`] if the initial seek to the contig
    ///   start (or the first refill) fails. Per-byte refill errors are
    ///   surfaced through the iterator's `Item = Result<u8, _>` rather
    ///   than from this constructor.
    fn iter_bases<'a>(
        &'a self,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>;

    /// Append `length` bases starting at 1-based `start` into the
    /// caller-supplied `dst`. The buffer is cleared first; on success
    /// `dst.len() == length`. Allocation-free in the hot path when
    /// `dst` is reused across calls — that's the reason this method
    /// exists alongside `fetch`, which returns an owned `Vec<u8>`.
    ///
    /// Default impl forwards to `fetch` and `mem::swap`s the
    /// returned `Vec<u8>` into `dst` (transferring ownership of the
    /// fresh buffer; the caller's previous `dst` is dropped). The
    /// production [`StreamingChromRefFetcher`] override copies
    /// straight from its sliding buffer, avoiding the per-call
    /// `Vec::with_capacity(length)` / drop pair entirely.
    ///
    /// M10 (2026-05-23 code review): the default impl was previously
    /// `dst.clear() + extend_from_slice(&v)`, which is a fresh
    /// allocation followed by a memcpy — exactly the per-call alloc
    /// the method exists to avoid. The `mem::swap` shape inherits
    /// the freshly-allocated buffer instead of memcpying through.
    ///
    /// # Errors
    ///
    /// Same failure modes as [`Self::fetch`]: `InvalidStart`,
    /// `OutOfBounds`, `OutOfPattern`, and `Io`.
    fn fetch_into(
        &self,
        start_1based: u32,
        length: u32,
        dst: &mut Vec<u8>,
    ) -> Result<(), ChromRefFetchError> {
        let mut v = self.fetch(start_1based, length)?;
        std::mem::swap(dst, &mut v);
        Ok(())
    }
}

/// Forwarding impl so callers may pass either an owned fetcher or a
/// shared reference to a generic `F: ChromRefFetcher`. The
/// `drive_cohort_pipeline` call site is the motivating consumer:
/// it holds a `SharedRefFetcher = Arc<dyn ChromRefFetcher + Send>`
/// and passes `&*fetcher` (a `&(dyn ChromRefFetcher + Send)`) into
/// `DustFilter::new`, which is generic over `F: ChromRefFetcher`.
impl<T: ChromRefFetcher + ?Sized> sealed::Sealed for &T {}
impl<T: ChromRefFetcher + ?Sized> ChromRefFetcher for &T {
    fn length(&self) -> u32 {
        (**self).length()
    }

    fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError> {
        (**self).fetch(start_1based, length)
    }

    fn iter_bases<'a>(
        &'a self,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
    {
        (**self).iter_bases()
    }

    fn fetch_into(
        &self,
        start_1based: u32,
        length: u32,
        dst: &mut Vec<u8>,
    ) -> Result<(), ChromRefFetchError> {
        (**self).fetch_into(start_1based, length, dst)
    }
}

impl StreamingChromRefFetcher {
    /// New-API constructor bound to one contig by name. Internally
    /// stores the same state as the legacy [`Self::new`] but with
    /// `chrom_id = 0` (the new trait doesn't expose `chrom_id`, so
    /// the stored value is meaningless for new-API consumers).
    ///
    /// Allocates a sliding buffer of [`STREAMING_REF_BUFFER_BYTES`]
    /// (= 1 MiB by default) — the buffer is filled lazily on the
    /// first fetch, not at construction. Tests can override the
    /// buffer size via the `#[cfg(test)]`-only
    /// `for_contig_with_buffer_bytes` constructor.
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::Io`] wrapping an
    ///   [`io::ErrorKind::NotFound`] error if `contig_name` is not
    ///   present in the sibling `.fai` index.
    /// - [`ChromRefFetchError::Io`] wrapping an
    ///   [`io::ErrorKind::InvalidData`] error if the `.fai` is
    ///   malformed (line_bases = 0 or line_width < line_bases — see
    ///   [`ContigFai::validate`]), or if any of the contig's
    ///   `length` / `line_bases` / `line_width` fields exceed
    ///   `u32::MAX`.
    /// - [`ChromRefFetchError::Io`] wrapping the underlying
    ///   [`io::Error`] from `noodles_fasta::fai::fs::read` or
    ///   `File::open` on a missing / unreadable FASTA or `.fai`
    ///   file.
    pub fn for_contig(fasta_path: &Path, contig_name: &str) -> Result<Self, ChromRefFetchError> {
        Self::for_contig_with_fai_path(fasta_path, None, contig_name)
    }

    /// Variant of [`Self::for_contig`] for non-standard `.fai`
    /// locations. `None` falls back to `<fasta_path>.fai`.
    ///
    /// # Errors
    ///
    /// Same as [`Self::for_contig`].
    pub fn for_contig_with_fai_path(
        fasta_path: &Path,
        fai_path: Option<&Path>,
        contig_name: &str,
    ) -> Result<Self, ChromRefFetchError> {
        Self::for_contig_internal(
            fasta_path,
            fai_path,
            contig_name,
            STREAMING_REF_BUFFER_BYTES,
        )
    }

    /// Test-only constructor allowing a custom buffer size so the
    /// `OutOfPattern` contract can be exercised with small fixtures
    /// (production code only ever sees `STREAMING_REF_BUFFER_BYTES`).
    #[cfg(test)]
    pub fn for_contig_with_buffer_bytes(
        fasta_path: &Path,
        contig_name: &str,
        buffer_bytes: usize,
    ) -> Result<Self, ChromRefFetchError> {
        Self::for_contig_internal(fasta_path, None, contig_name, buffer_bytes)
    }

    fn for_contig_internal(
        fasta_path: &Path,
        fai_path: Option<&Path>,
        contig_name: &str,
        buffer_bytes: usize,
    ) -> Result<Self, ChromRefFetchError> {
        let mut fai_pathbuf: OsString;
        let fai_path_owned: PathBuf;
        let fai_path_resolved: &Path = match fai_path {
            Some(p) => p,
            None => {
                fai_pathbuf = fasta_path.as_os_str().to_os_string();
                fai_pathbuf.push(".fai");
                fai_path_owned = PathBuf::from(fai_pathbuf);
                &fai_path_owned
            }
        };
        let index = fai::fs::read(fai_path_resolved).map_err(|source| ChromRefFetchError::Io {
            contig_name: contig_name.to_string(),
            source,
        })?;
        let record = index
            .as_ref()
            .iter()
            .find(|r| AsRef::<[u8]>::as_ref(r.name()) == contig_name.as_bytes())
            .ok_or_else(|| ChromRefFetchError::Io {
                contig_name: contig_name.to_string(),
                source: io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("contig {contig_name} not in FASTA index"),
                ),
            })?;
        let length = u32::try_from(record.length()).map_err(|_| ChromRefFetchError::Io {
            contig_name: contig_name.to_string(),
            source: io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig {} length {} exceeds u32::MAX",
                    contig_name,
                    record.length()
                ),
            ),
        })?;
        let line_bases =
            u32::try_from(record.line_bases()).map_err(|_| ChromRefFetchError::Io {
                contig_name: contig_name.to_string(),
                source: io::Error::new(
                    io::ErrorKind::InvalidData,
                    ".fai line_bases exceeds u32::MAX",
                ),
            })?;
        let line_width =
            u32::try_from(record.line_width()).map_err(|_| ChromRefFetchError::Io {
                contig_name: contig_name.to_string(),
                source: io::Error::new(
                    io::ErrorKind::InvalidData,
                    ".fai line_width exceeds u32::MAX",
                ),
            })?;
        let fai = ContigFai {
            seq_offset: record.offset(),
            length,
            line_bases,
            line_width,
        };
        // B1: reject malformed .fai (line_bases = 0, etc.) before
        // any fetch can panic on the offset arithmetic.
        fai.validate(contig_name)
            .map_err(|source| ChromRefFetchError::Io {
                contig_name: contig_name.to_string(),
                source,
            })?;
        let file = File::open(fasta_path).map_err(|source| ChromRefFetchError::Io {
            contig_name: contig_name.to_string(),
            source,
        })?;
        Ok(Self {
            chrom_id: 0, // unused under the new API
            contig_name: contig_name.to_string(),
            fai,
            inner: RefCell::new(StreamState {
                source: Source::File(file),
                buf: Vec::with_capacity(buffer_bytes),
                buf_start_base: 0,
            }),
        })
    }
}

/// Forward sequential iterator returned by
/// [`ChromRefFetcher::iter_bases`]. Holds a borrow on the fetcher;
/// the borrow guarantees no concurrent `fetch()` can race with the
/// buffer-state resets the iter performs on construction and on
/// Drop.
pub struct ChromRefBaseIter<'a> {
    fetcher: &'a StreamingChromRefFetcher,
    /// 1-based coordinate of the next base to yield. Starts at 1
    /// after the construction-time reset.
    next_base: u32,
    /// Latches on the first error so subsequent `next()` calls
    /// return `None` (fused-error convention).
    done: bool,
}

impl Iterator for ChromRefBaseIter<'_> {
    type Item = Result<u8, ChromRefFetchError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done || self.next_base > self.fetcher.fai.length {
            return None;
        }
        let mut state = self.fetcher.inner.borrow_mut();
        let need_refill = state.buf.is_empty()
            || self.next_base < state.buf_start_base
            || (self.next_base as u64) >= state.buf_start_base as u64 + state.buf.len() as u64;
        if need_refill {
            // We don't know the buffer size from this scope — use the
            // production constant; tests that need a custom size go
            // through `fetch` rather than the iter.
            if let Err(source) = refill(
                &mut state,
                &self.fetcher.fai,
                self.next_base,
                STREAMING_REF_BUFFER_BYTES,
            ) {
                self.done = true;
                return Some(Err(ChromRefFetchError::Io {
                    contig_name: self.fetcher.contig_name.clone(),
                    source,
                }));
            }
        }
        let idx = (self.next_base - state.buf_start_base) as usize;
        let byte = state.buf[idx];
        self.next_base += 1;
        Some(Ok(byte))
    }
}

impl Drop for ChromRefBaseIter<'_> {
    fn drop(&mut self) {
        // Reset the sliding buffer so the next `fetch()` starts a
        // fresh access-pattern phase regardless of where the iter
        // left off.
        //
        // M15 (2026-05-23 code review): `try_borrow_mut` (not
        // `borrow_mut`) so a panic that escaped from `next()` while
        // holding its own `borrow_mut` doesn't double-panic and
        // abort the process. The buffer reset is best-effort:
        // skipping it is safe because the next `iter_bases()` call
        // resets the buffer at the top, and any post-panic
        // `fetch()` would see the stale buffer and either return
        // OutOfPattern (correct) or refill (correct on monotonic
        // access).
        if let Ok(mut state) = self.fetcher.inner.try_borrow_mut() {
            state.buf.clear();
            state.buf_start_base = 0;
        }
    }
}

impl sealed::Sealed for StreamingChromRefFetcher {}
impl ChromRefFetcher for StreamingChromRefFetcher {
    fn length(&self) -> u32 {
        self.fai.length
    }

    fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError> {
        // fetch is just fetch_into into a fresh Vec; the heavy
        // lifting lives in fetch_into.
        let mut dst = Vec::with_capacity(length as usize);
        self.fetch_into(start_1based, length, &mut dst)?;
        Ok(dst)
    }

    fn fetch_into(
        &self,
        start_1based: u32,
        length: u32,
        dst: &mut Vec<u8>,
    ) -> Result<(), ChromRefFetchError> {
        if start_1based == 0 {
            return Err(ChromRefFetchError::InvalidStart);
        }
        let end_1based =
            start_1based
                .checked_add(length)
                .ok_or_else(|| ChromRefFetchError::Io {
                    contig_name: self.contig_name.clone(),
                    source: io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"),
                })?;
        if end_1based.saturating_sub(1) > self.fai.length {
            return Err(ChromRefFetchError::OutOfBounds {
                contig_name: self.contig_name.clone(),
                contig_length: self.fai.length,
                start: start_1based,
                end: end_1based,
            });
        }

        let mut state = self.inner.borrow_mut();

        // OutOfPattern check: if the buffer is non-empty and the
        // request is before its origin, the consumer is violating
        // monotonic non-decreasing access.
        if !state.buf.is_empty() && start_1based < state.buf_start_base {
            return Err(ChromRefFetchError::OutOfPattern {
                requested_start: start_1based,
                buffer_origin: state.buf_start_base,
            });
        }

        let buf_covers = !state.buf.is_empty()
            && start_1based >= state.buf_start_base
            && (end_1based as u64) <= state.buf_start_base as u64 + state.buf.len() as u64;
        if !buf_covers {
            // Refill forward (or initial fill on an empty buffer).
            // The forward direction is guaranteed by the OutOfPattern
            // check above plus the bounds check; we never call refill
            // with a target before the current origin.
            let want = (length as usize).max(state.buf.capacity());
            refill(&mut state, &self.fai, start_1based, want).map_err(|source| {
                ChromRefFetchError::Io {
                    contig_name: self.contig_name.clone(),
                    source,
                }
            })?;
        }
        let start_idx = (start_1based - state.buf_start_base) as usize;
        let end_idx = start_idx + length as usize;
        dst.clear();
        dst.extend_from_slice(&state.buf[start_idx..end_idx]);
        Ok(())
    }

    fn iter_bases<'a>(
        &'a self,
    ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
    {
        // Reset the buffer so the iter starts a fresh phase.
        let mut state = self.inner.borrow_mut();
        state.buf.clear();
        state.buf_start_base = 0;
        drop(state);
        Ok(Box::new(ChromRefBaseIter {
            fetcher: self,
            next_base: 1,
            done: false,
        }))
    }
}

// ---------------------------------------------------------------------
// WalkerLegacyAdapter — multi-chrom, single-thread bridge.
// ---------------------------------------------------------------------
//
// Bridge for consumers that traverse multiple chromosomes in
// sequence with a single thread — the Stage 1 pileup walker today.
// On each chrom transition the inner [`StreamingChromRefFetcher`]
// is rebuilt for the new contig; fetches within a chrom hit the
// same sliding buffer.
//
// Transitional. A future migration replaces this adapter with a
// direct `ChromRefFetcher` consumer in the Stage 1 walker.

/// Adapter that maintains *one* [`StreamingChromRefFetcher`] at a
/// time, swapping it when [`RefSeqFetcher::fetch`] is called with a
/// different `chrom_id` than the currently-bound contig. Implements
/// the legacy `RefSeqFetcher` trait so it drops into existing
/// walker plumbing.
///
/// `Sync` via internal `Mutex<Option<...>>`. The walker itself is
/// single-threaded so the Mutex is uncontended; `var_calling_from_bam`
/// also uses this adapter and passes it through `Arc<dyn
/// RefSeqFetcher + Send + Sync>`, so the `Sync` bound is real.
/// Stage 1's BAQ uses a separate per-rayon-worker
/// [`ManualEvictChromRefFetcher`] (built inside `BaqStream`'s
/// `map_init`); BAQ does **not** share a fetcher with the walker
/// — the access patterns are too different for one abstraction
/// to serve both.
///
/// Behaviour on chrom transition: the existing inner is dropped
/// (releasing its sliding buffer), a fresh
/// [`StreamingChromRefFetcher::for_contig`] is built for the new
/// chrom (open file + parse .fai), and the call delegates via the
/// legacy *tolerant* `RefSeqFetcher::fetch` path (silent backward
/// refill). The walker doesn't need the strict `OutOfPattern`
/// contract — its single-thread sequential access pattern doesn't
/// produce backward jumps in the first place, and tolerating the
/// rare one is simpler than adding monotonicity to its contract.
pub struct WalkerLegacyAdapter {
    fasta_path: PathBuf,
    contigs: ContigList,
    inner: Mutex<Option<(u32, StreamingChromRefFetcher)>>,
}

impl WalkerLegacyAdapter {
    /// Build an adapter bound to a FASTA. No inner fetcher is
    /// constructed yet — the first `fetch` triggers the lazy load.
    pub fn new(fasta_path: PathBuf, contigs: ContigList) -> Self {
        Self {
            fasta_path,
            contigs,
            inner: Mutex::new(None),
        }
    }
}

impl RefSeqFetcher for WalkerLegacyAdapter {
    fn fetch(&self, chrom_id: u32, start_1based: u32, length: u32) -> Result<Vec<u8>, io::Error> {
        let mut slot = self.inner.lock().unwrap_or_else(|e| e.into_inner());
        let needs_rebuild = match slot.as_ref() {
            Some((current, _)) => *current != chrom_id,
            None => true,
        };
        if needs_rebuild {
            let entry = self.contigs.entries.get(chrom_id as usize).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "WalkerLegacyAdapter: chrom_id {chrom_id} out of range (have {} contigs)",
                        self.contigs.entries.len()
                    ),
                )
            })?;
            // Drop the old fetcher (releases its sliding buffer)
            // before building the new one to keep peak memory at
            // one contig's working set, not two.
            *slot = None;
            let fresh = StreamingChromRefFetcher::for_contig(&self.fasta_path, &entry.name)
                .map_err(legacy_io_error)?;
            *slot = Some((chrom_id, fresh));
        }
        let (_, inner) = slot.as_ref().expect("post-rebuild slot is populated");
        // Use the legacy *tolerant* `RefSeqFetcher::fetch` on the
        // inner (silent backward refill), not the new strict
        // `ChromRefFetcher::fetch` with the `OutOfPattern` check.
        // The strict contract suits the cohort var-calling path
        // (monotonic by construction); the walker just needs a
        // working multi-chrom fetcher and doesn't care about the
        // contract.
        //
        // The inner's `chrom_id` field is 0 (set by `for_contig`);
        // the legacy `fetch`'s `chrom_id == self.chrom_id` check
        // accepts 0.
        RefSeqFetcher::fetch(inner, 0, start_1based, length)
    }

    // iter_bases falls back to the trait default. The walker doesn't
    // use it (only DUST / PerGroupMerger in the var-calling pipeline
    // do), so defining it explicitly would be dead code.
}

/// Map [`ChromRefFetchError`] back into the legacy `io::Error`
/// shape. The legacy trait returns `io::Result`; until the
/// consumer-side migration in step 3 we have to flatten the typed
/// error here. Each variant picks the closest `io::ErrorKind` so
/// the legacy error-handling code still routes correctly.
fn legacy_io_error(e: ChromRefFetchError) -> io::Error {
    match e {
        ChromRefFetchError::Io { source, .. } => source,
        ChromRefFetchError::InvalidStart => {
            io::Error::new(io::ErrorKind::InvalidInput, "start_1based must be >= 1")
        }
        ChromRefFetchError::OutOfBounds { .. } => {
            io::Error::new(io::ErrorKind::UnexpectedEof, e.to_string())
        }
        ChromRefFetchError::OutOfPattern { .. } => {
            io::Error::new(io::ErrorKind::InvalidInput, e.to_string())
        }
    }
}

// ---------------------------------------------------------------------
// ManualEvictChromRefFetcher — caller-managed buffer lifecycle.
// ---------------------------------------------------------------------
//
// Designed for consumers whose access pattern is genuinely random
// access *within a known region* — Stage 1 BAQ is the motivating
// example. Each rayon worker owns its own instance via `map_init`
// and processes reads that may visit positions in any order within
// the current chunk's coord range. The worker calls `evict_before`
// at deterministic points to keep memory bounded.
//
// Contract:
//
// - `fetch(start, length)` returns a borrowed `&[u8]` slice. If the
//   range isn't covered, the buffer extends in whichever direction
//   is needed — forward by appending fresh bytes, backward by
//   prepending. Never shrinks on its own.
// - `evict_before(pos)` is the *only* shrinking operation. The
//   caller decides when bytes are no longer needed. CRAM is
//   coordinate-sorted, so calling `evict_before(read.pos)` after
//   each read keeps the buffer at "current-position window plus
//   whatever the next read extends forward."
//
// Trade-off vs [`StreamingChromRefFetcher`]:
//
// - Streaming forces monotonic forward access (small look-back
//   within the buffer is fine; large backward jumps are
//   `OutOfPattern`). Good for cohort var-calling.
// - ManualEvict tolerates any access pattern within the resident
//   range and trusts the caller to keep memory bounded. Good for
//   rayon-parallel random-access-within-chunk consumers. Worst case
//   (no `evict_before` calls + adversarial input) is one whole
//   contig resident — acceptable because the consumer controls it.
//
// Not `Sync`. Each thread holds its own instance; concurrent
// fetches from multiple threads on a single instance are forbidden.
// `Send` because every field is `Send`.

/// Reference fetcher with caller-managed buffer lifecycle. See the
/// module-level section comment above for the contract.
pub struct ManualEvictChromRefFetcher {
    contig_name: String,
    fai: ContigFai,
    file: File,
    /// Uppercased, newline-stripped bases currently resident.
    buf: Vec<u8>,
    /// 1-based contig coordinate of `buf[0]`. Meaningful only when
    /// `!buf.is_empty()`; treated as `1` when the buffer is empty.
    buf_start_base: u32,
}

impl ManualEvictChromRefFetcher {
    /// Open the FASTA at `fasta_path`, parse the sibling `.fai`,
    /// look up `contig_name`, return a fetcher with an empty
    /// buffer ready to serve `fetch` calls.
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::Io`] with `ErrorKind::NotFound` if
    ///   `contig_name` is absent from the FASTA index.
    /// - [`ChromRefFetchError::Io`] with `ErrorKind::InvalidData` if
    ///   the `.fai` entry is malformed (zero `line_bases`,
    ///   `line_width < line_bases`, lengths overflowing `u32`).
    /// - [`ChromRefFetchError::Io`] propagating the underlying
    ///   `io::Error` if the `.fai` read or `File::open(fasta_path)`
    ///   fails.
    pub fn for_contig(fasta_path: &Path, contig_name: &str) -> Result<Self, ChromRefFetchError> {
        Self::for_contig_with_fai_path(fasta_path, None, contig_name)
    }

    /// Variant of [`Self::for_contig`] for non-standard `.fai`
    /// locations. `None` falls back to `<fasta_path>.fai`.
    ///
    /// # Errors
    ///
    /// Same as [`Self::for_contig`].
    pub fn for_contig_with_fai_path(
        fasta_path: &Path,
        fai_path: Option<&Path>,
        contig_name: &str,
    ) -> Result<Self, ChromRefFetchError> {
        let mut fai_pathbuf: OsString;
        let fai_path_owned: PathBuf;
        let fai_path_resolved: &Path = match fai_path {
            Some(p) => p,
            None => {
                fai_pathbuf = fasta_path.as_os_str().to_os_string();
                fai_pathbuf.push(".fai");
                fai_path_owned = PathBuf::from(fai_pathbuf);
                &fai_path_owned
            }
        };
        let index = fai::fs::read(fai_path_resolved).map_err(|source| ChromRefFetchError::Io {
            contig_name: contig_name.to_string(),
            source,
        })?;
        let record = index
            .as_ref()
            .iter()
            .find(|r| AsRef::<[u8]>::as_ref(r.name()) == contig_name.as_bytes())
            .ok_or_else(|| ChromRefFetchError::Io {
                contig_name: contig_name.to_string(),
                source: io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("contig {contig_name} not in FASTA index"),
                ),
            })?;
        let length = u32::try_from(record.length()).map_err(|_| ChromRefFetchError::Io {
            contig_name: contig_name.to_string(),
            source: io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "contig {} length {} exceeds u32::MAX",
                    contig_name,
                    record.length()
                ),
            ),
        })?;
        let line_bases =
            u32::try_from(record.line_bases()).map_err(|_| ChromRefFetchError::Io {
                contig_name: contig_name.to_string(),
                source: io::Error::new(
                    io::ErrorKind::InvalidData,
                    ".fai line_bases exceeds u32::MAX",
                ),
            })?;
        let line_width =
            u32::try_from(record.line_width()).map_err(|_| ChromRefFetchError::Io {
                contig_name: contig_name.to_string(),
                source: io::Error::new(
                    io::ErrorKind::InvalidData,
                    ".fai line_width exceeds u32::MAX",
                ),
            })?;
        let fai = ContigFai {
            seq_offset: record.offset(),
            length,
            line_bases,
            line_width,
        };
        // B1: reject malformed .fai (line_bases = 0, etc.) before
        // any fetch can panic on the offset arithmetic.
        fai.validate(contig_name)
            .map_err(|source| ChromRefFetchError::Io {
                contig_name: contig_name.to_string(),
                source,
            })?;
        let file = File::open(fasta_path).map_err(|source| ChromRefFetchError::Io {
            contig_name: contig_name.to_string(),
            source,
        })?;
        Ok(Self {
            contig_name: contig_name.to_string(),
            fai,
            file,
            buf: Vec::new(),
            buf_start_base: 1,
        })
    }

    /// Total bases in the bound contig.
    pub fn length(&self) -> u32 {
        self.fai.length
    }

    /// Return uppercased bases `[start_1based, start_1based+length)`
    /// as a borrowed slice. If the buffer doesn't cover the
    /// requested range, extends in whichever direction is needed
    /// (forward by appending, backward by prepending).
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::InvalidStart`] if `start_1based == 0`.
    /// - [`ChromRefFetchError::OutOfBounds`] if
    ///   `start_1based + length - 1 > self.length()`.
    /// - [`ChromRefFetchError::Io`] if a refill from the underlying
    ///   FASTA fails, or if the requested range overflows `u32`.
    pub fn fetch(&mut self, start_1based: u32, length: u32) -> Result<&[u8], ChromRefFetchError> {
        if start_1based == 0 {
            return Err(ChromRefFetchError::InvalidStart);
        }
        let end_exclusive =
            start_1based
                .checked_add(length)
                .ok_or_else(|| ChromRefFetchError::Io {
                    contig_name: self.contig_name.clone(),
                    source: io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"),
                })?;
        // [start_1based, end_exclusive); valid only if end <= length+1.
        if end_exclusive.saturating_sub(1) > self.fai.length {
            return Err(ChromRefFetchError::OutOfBounds {
                contig_name: self.contig_name.clone(),
                contig_length: self.fai.length,
                start: start_1based,
                end: end_exclusive,
            });
        }

        if self.buf.is_empty() {
            // Empty buffer: read [start, end) from disk fresh.
            self.read_into_buffer_at(start_1based, length as usize)?;
            return Ok(&self.buf[..length as usize]);
        }

        let buf_end_exclusive = self.buf_start_base as u64 + self.buf.len() as u64;
        let want_start = start_1based as u64;
        let want_end_exclusive = end_exclusive as u64;
        let buf_start = self.buf_start_base as u64;

        // Extend forward if requested end exceeds buffer end.
        if want_end_exclusive > buf_end_exclusive {
            let extra_bases = (want_end_exclusive - buf_end_exclusive) as usize;
            self.append_forward(extra_bases)?;
        }
        // Extend backward if requested start precedes buffer start.
        if want_start < buf_start {
            let extra_bases = (buf_start - want_start) as usize;
            self.prepend_backward(start_1based, extra_bases)?;
        }

        let local_start = (start_1based - self.buf_start_base) as usize;
        let local_end = local_start + length as usize;
        Ok(&self.buf[local_start..local_end])
    }

    /// Drop buffer bytes whose 1-based base coordinate is `< pos`.
    /// Compacts in place; the `Vec`'s capacity is retained so the
    /// next forward fetch can reuse the freed slots without
    /// reallocating.
    ///
    /// No-op if `pos <= buf_start_base` or the buffer is empty.
    /// If `pos` would drop the entire buffer, clears it
    /// (and updates `buf_start_base` to `pos`).
    pub fn evict_before(&mut self, pos: u32) {
        if self.buf.is_empty() || pos <= self.buf_start_base {
            return;
        }
        let buf_end_exclusive = self.buf_start_base as u64 + self.buf.len() as u64;
        if (pos as u64) >= buf_end_exclusive {
            self.buf.clear();
            self.buf_start_base = pos;
            return;
        }
        let drop_count = (pos - self.buf_start_base) as usize;
        // `Vec::drain(..n)` shifts the remainder to the front; the
        // backing allocation is preserved.
        self.buf.drain(..drop_count);
        self.buf_start_base = pos;
    }

    /// Current resident buffer length in bytes. Test/diagnostic.
    #[cfg(test)]
    pub fn buf_len(&self) -> usize {
        self.buf.len()
    }

    /// Current 1-based buffer origin. Test/diagnostic.
    #[cfg(test)]
    pub fn buf_start_base(&self) -> u32 {
        self.buf_start_base
    }

    /// Read `n_bases` uppercased, newline-stripped bases starting
    /// at 1-based `start` directly into `self.buf` (which is empty
    /// on entry). Sets `buf_start_base = start`.
    fn read_into_buffer_at(
        &mut self,
        start_1based: u32,
        n_bases: usize,
    ) -> Result<(), ChromRefFetchError> {
        debug_assert!(self.buf.is_empty());
        let offset = self.fai.base_to_file_offset(start_1based);
        self.file
            .seek(SeekFrom::Start(offset))
            .map_err(|source| ChromRefFetchError::Io {
                contig_name: self.contig_name.clone(),
                source,
            })?;
        let remaining = (self.fai.length - start_1based + 1) as usize;
        let target_len = n_bases.min(remaining);
        self.buf.reserve(target_len);
        self.buf_start_base = start_1based;
        read_uppercased_bases(&mut self.file, &mut self.buf, target_len).map_err(|source| {
            ChromRefFetchError::Io {
                contig_name: self.contig_name.clone(),
                source,
            }
        })
    }

    /// Append `extra_bases` more bases to the end of the buffer by
    /// reading from disk just past the current end.
    fn append_forward(&mut self, extra_bases: usize) -> Result<(), ChromRefFetchError> {
        let next_base = self.buf_start_base + self.buf.len() as u32;
        let offset = self.fai.base_to_file_offset(next_base);
        self.file
            .seek(SeekFrom::Start(offset))
            .map_err(|source| ChromRefFetchError::Io {
                contig_name: self.contig_name.clone(),
                source,
            })?;
        let remaining = (self.fai.length - next_base + 1) as usize;
        let take = extra_bases.min(remaining);
        self.buf.reserve(take);
        read_uppercased_bases(&mut self.file, &mut self.buf, take).map_err(|source| {
            ChromRefFetchError::Io {
                contig_name: self.contig_name.clone(),
                source,
            }
        })
    }

    /// Prepend `extra_bases` bases to the front of the buffer by
    /// reading from disk at `new_start` and memmove-shifting the
    /// existing buffer content to make room.
    fn prepend_backward(
        &mut self,
        new_start: u32,
        extra_bases: usize,
    ) -> Result<(), ChromRefFetchError> {
        debug_assert!(new_start < self.buf_start_base);
        debug_assert_eq!((self.buf_start_base - new_start) as usize, extra_bases,);

        // Build the prefix bytes in a scratch Vec, then splice into
        // the front of `self.buf`. `Vec::splice` does the memmove
        // for us and reuses capacity when possible.
        let offset = self.fai.base_to_file_offset(new_start);
        self.file
            .seek(SeekFrom::Start(offset))
            .map_err(|source| ChromRefFetchError::Io {
                contig_name: self.contig_name.clone(),
                source,
            })?;
        let mut prefix = Vec::with_capacity(extra_bases);
        read_uppercased_bases(&mut self.file, &mut prefix, extra_bases).map_err(|source| {
            ChromRefFetchError::Io {
                contig_name: self.contig_name.clone(),
                source,
            }
        })?;
        self.buf.splice(..0, prefix);
        self.buf_start_base = new_start;
        Ok(())
    }
}

/// Read `n_bases` uppercased, newline-stripped bases from `reader`
/// (currently positioned at the file offset for those bases) and
/// append them to `dst`. Returns an `UnexpectedEof` if the source
/// runs out before `n_bases` bases are gathered. Shared by both
/// streaming and manual-evict fetchers.
/// Canonicalise a single FASTA byte to SAM-spec uppercase `ACGT`/`N`.
///
/// M26 (2026-05-23 code review): the `ChromRefFetcher::fetch` trait
/// docs promise output is "uppercased bases, SAM-spec canonical
/// (`A`/`C`/`G`/`T`/`N`) regardless of how the FASTA is masked on
/// disk". The previous implementation used only
/// `b.to_ascii_uppercase()`, which leaves IUPAC ambiguity codes
/// (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`), the RNA `U`,
/// gap characters (`-`, `.`), and any other ASCII byte unchanged.
/// Downstream DUST + BAQ + PerGroupMerger then saw raw IUPAC bytes
/// that they were not built to handle.
///
/// This function folds any non-`ACGT` byte to `N` AFTER
/// uppercasing — so `a`/`c`/`g`/`t` (soft-masked) map to `A`/`C`/
/// `G`/`T`, and everything else (`R`, `n`, `N`, `U`, `-`, `.`,
/// stray whitespace bytes that survived the newline filter) maps
/// to `N`.
#[inline]
fn canonicalise(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => b.to_ascii_uppercase(),
        _ => b'N',
    }
}

fn read_uppercased_bases(reader: &mut File, dst: &mut Vec<u8>, n_bases: usize) -> io::Result<()> {
    let want_total = dst.len() + n_bases;
    let mut read_buf = [0u8; STREAMING_REF_FILE_READ_CHUNK];
    while dst.len() < want_total {
        let n = reader.read(&mut read_buf)?;
        if n == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "FASTA ended with {} base(s) still expected (want_total {})",
                    want_total - dst.len(),
                    want_total,
                ),
            ));
        }
        for &b in &read_buf[..n] {
            if dst.len() >= want_total {
                break;
            }
            if b == b'\n' || b == b'\r' {
                continue;
            }
            dst.push(canonicalise(b));
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
    ) -> (
        Vec<u8>,
        u64, /* seq_offset */
        u32, /* line_width */
    ) {
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

        let bytes = RefSeqFetcher::fetch(&fetcher, 0, 5, 4).expect("fetch");
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
        let bytes = RefSeqFetcher::fetch(&fetcher, 0, 58, 5).expect("fetch");
        assert_eq!(bytes, seq[57..62]);
        // Wholly inside one line.
        let bytes = RefSeqFetcher::fetch(&fetcher, 0, 70, 10).expect("fetch");
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
        let bytes = RefSeqFetcher::fetch(&fetcher, 0, 1, seq.len() as u32).expect("fetch");
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
        let seq: Vec<u8> = (0..2_000_000_u32)
            .map(|i| b"ACGT"[(i % 4) as usize])
            .collect();
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
        let a = RefSeqFetcher::fetch(&fetcher, 0, 100, 16).expect("fetch a");
        assert_eq!(a, seq[99..115]);
        // Then jump well past the 1 MB buffer.
        let b = RefSeqFetcher::fetch(&fetcher, 0, 1_500_000, 16).expect("fetch b");
        assert_eq!(b, seq[1_499_999..1_500_015]);
    }

    #[test]
    fn streaming_fetcher_buffer_refill_on_backward_jump() {
        // Backward jump after a forward fetch. The streamer should
        // detect the buffer miss and refill from the new (earlier)
        // position. Documents the slower path; production access
        // never hits this (PerGroupMerger is monotonic), but tests
        // pin the behaviour.
        let seq: Vec<u8> = (0..2_000_000_u32)
            .map(|i| b"ACGT"[(i % 4) as usize])
            .collect();
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
        let _ = RefSeqFetcher::fetch(&fetcher, 0, 1_500_000, 16).expect("fetch forward");
        let back = RefSeqFetcher::fetch(&fetcher, 0, 100, 16).expect("fetch back");
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
        let err = RefSeqFetcher::fetch(&fetcher, 0, 8, 5).expect_err("must fail");
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
        let err = RefSeqFetcher::fetch(&fetcher, 7, 1, 4).expect_err("must fail");
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
        let err = RefSeqFetcher::fetch(&fetcher, 0, 0, 1).expect_err("must fail");
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
        // .fai) returns the expected uppercase 'A' bytes that
        // `build_fasta` writes for a 64-base contig.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 64,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let streaming = StreamingChromRefFetcher::new(&path, 0, "chr0".into()).expect("fetcher");

        for (start, length) in [(1u32, 16u32), (8, 32), (50, 14), (64, 1)] {
            let bytes = RefSeqFetcher::fetch(&streaming, 0, start, length).expect("streaming");
            let expected = vec![b'A'; length as usize];
            assert_eq!(bytes, expected, "mismatch at start={start} length={length}");
        }
    }

    // ---------------------------------------------------------------
    // ChromRefFetcher (new API) tests
    // ---------------------------------------------------------------
    //
    // The new trait coexists with the legacy `RefSeqFetcher` impl on
    // the same struct. These tests exercise the new-API methods:
    // `for_contig` construction, `length()`, `fetch(start, length)`
    // (without chrom_id), `iter_bases()` (with buffer reset
    // semantics), and the `OutOfPattern` contract on backward
    // jumps beyond the buffer.

    use std::fs::File as StdFile;
    use std::io::Write as _;

    /// Build a real FASTA on disk with a single contig of
    /// `seq.len()` bases wrapped at `line_bases` per line. Returns
    /// the tempdir (kept alive for the test) and the FASTA path.
    fn build_line_wrapped_fasta(
        name: &str,
        seq: &[u8],
        line_bases: usize,
    ) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        let mut fa = StdFile::create(&fasta_path).expect("fa");
        let header = format!(">{name}\n");
        fa.write_all(header.as_bytes()).expect("hdr");
        let seq_offset = header.len() as u64;
        let mut i = 0;
        while i < seq.len() {
            let end = (i + line_bases).min(seq.len());
            fa.write_all(&seq[i..end]).expect("seq");
            fa.write_all(b"\n").expect("nl");
            i = end;
        }
        let mut fai = StdFile::create(&fai_path).expect("fai");
        writeln!(
            fai,
            "{}\t{}\t{}\t{}\t{}",
            name,
            seq.len(),
            seq_offset,
            line_bases,
            line_bases + 1
        )
        .expect("fai entry");
        (dir, fasta_path)
    }

    #[test]
    fn chrom_ref_for_contig_returns_bytes_for_bound_contig() {
        // Positive path: build a 50-base 'A' FASTA, construct via
        // the new `for_contig`, fetch a window. length() reports the
        // bound contig length.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 50,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");
        assert_eq!(ChromRefFetcher::length(&fetcher), 50);

        let bytes = ChromRefFetcher::fetch(&fetcher, 5, 4).expect("fetch");
        assert_eq!(bytes, b"AAAA");
    }

    #[test]
    fn chrom_ref_serves_line_wrapped_contig() {
        // 60-col wrap (Ensembl default). Multiple fetches across
        // line boundaries; every one returns the correct bases.
        let seq: Vec<u8> = (0..250).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 60);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");

        let bytes = ChromRefFetcher::fetch(&fetcher, 58, 5).expect("fetch");
        assert_eq!(bytes, seq[57..62]);
        let bytes = ChromRefFetcher::fetch(&fetcher, 70, 10).expect("fetch");
        assert_eq!(bytes, seq[69..79]);
    }

    #[test]
    fn chrom_ref_fetch_uppercases_soft_masked() {
        // Lowercase input → uppercased output.
        let seq = b"ACGTacgtNn".to_vec();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 18);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");

        let bytes = ChromRefFetcher::fetch(&fetcher, 1, seq.len() as u32).expect("fetch");
        assert_eq!(bytes, b"ACGTACGTNN");
    }

    #[test]
    fn chrom_ref_fetch_past_contig_end_returns_out_of_bounds() {
        let seq = b"ACGTACGTAC".to_vec();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 10);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");

        let err = ChromRefFetcher::fetch(&fetcher, 8, 5).expect_err("must fail");
        match err {
            ChromRefFetchError::OutOfBounds {
                contig_length,
                start,
                end,
                ..
            } => {
                assert_eq!(contig_length, 10);
                assert_eq!(start, 8);
                assert_eq!(end, 13);
            }
            other => panic!("expected OutOfBounds, got {other:?}"),
        }
    }

    #[test]
    fn chrom_ref_fetch_start_zero_is_rejected() {
        let seq = b"ACGT".to_vec();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 4);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");

        let err = ChromRefFetcher::fetch(&fetcher, 0, 1).expect_err("must fail");
        match err {
            ChromRefFetchError::InvalidStart => {}
            other => panic!("expected InvalidStart, got {other:?}"),
        }
    }

    #[test]
    fn chrom_ref_construct_missing_contig_errors() {
        // The supplied contig name is not in the FASTA. Construction
        // must fail (rather than build a zero-length fetcher).
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 10,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");

        let result = StreamingChromRefFetcher::for_contig(&path, "chr_missing");
        match result {
            Err(ChromRefFetchError::Io {
                contig_name,
                source,
            }) => {
                assert_eq!(contig_name, "chr_missing");
                assert_eq!(source.kind(), io::ErrorKind::NotFound);
            }
            Ok(_) => panic!("expected Err, got Ok"),
            Err(other) => panic!("expected Io, got {other:?}"),
        }
    }

    #[test]
    fn chrom_ref_buffer_refill_on_forward_jump() {
        // Large forward jump triggers a refill; the second fetch
        // returns the correct bytes from the new buffer window.
        let seq: Vec<u8> = (0..2_000_000_u32)
            .map(|i| b"ACGT"[(i % 4) as usize])
            .collect();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 80);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");

        let a = ChromRefFetcher::fetch(&fetcher, 100, 16).expect("fetch a");
        assert_eq!(a, seq[99..115]);
        let b = ChromRefFetcher::fetch(&fetcher, 1_500_000, 16).expect("fetch b");
        assert_eq!(b, seq[1_499_999..1_500_015]);
    }

    #[test]
    fn chrom_ref_small_backward_within_buffer() {
        // Backward jump that stays inside the current sliding
        // buffer: should serve from the buffer without error.
        let seq: Vec<u8> = (0..2_000_000_u32)
            .map(|i| b"ACGT"[(i % 4) as usize])
            .collect();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 80);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");

        let _first = ChromRefFetcher::fetch(&fetcher, 100, 16).expect("forward 1");
        let _second = ChromRefFetcher::fetch(&fetcher, 100_000, 16).expect("forward 2");
        // Small backward inside the 1 MB buffer that opened at base
        // 100_000: any position in [100_000, 100_000 + 1MB) is
        // valid. Pick 100_100.
        let backward = ChromRefFetcher::fetch(&fetcher, 100_100, 16).expect("backward in buffer");
        assert_eq!(backward, seq[100_099..100_115]);
    }

    #[test]
    fn chrom_ref_backward_beyond_buffer_returns_out_of_pattern() {
        // Build a fetcher with a deliberately small buffer (4 KB),
        // advance past it via a forward jump, then issue a backward
        // fetch that lies before the buffer's current origin.
        // OutOfPattern must fire with the right field values —
        // pins the contract for tests.
        let seq: Vec<u8> = (0..20_000_u32).map(|i| b"ACGT"[(i % 4) as usize]).collect();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 80);
        let fetcher =
            StreamingChromRefFetcher::for_contig_with_buffer_bytes(&path, "chr0", 4 * 1024)
                .expect("fetcher");

        // Forward jump: loads buffer at start=10_000, capacity 4 KB.
        let _ = ChromRefFetcher::fetch(&fetcher, 10_000, 16).expect("forward");
        // Backward to start=100: before the buffer's origin (10_000),
        // so OutOfPattern.
        let err = ChromRefFetcher::fetch(&fetcher, 100, 16).expect_err("must fail");
        match err {
            ChromRefFetchError::OutOfPattern {
                requested_start,
                buffer_origin,
            } => {
                assert_eq!(requested_start, 100);
                assert_eq!(buffer_origin, 10_000);
            }
            other => panic!("expected OutOfPattern, got {other:?}"),
        }
    }

    #[test]
    fn chrom_ref_iter_bases_yields_uppercased_contig_in_order() {
        // The iter_bases iter visits every contig base 1..=length in
        // order, uppercased.
        let seq: Vec<u8> = (0..200).map(|i| b"acgtACGTN"[i % 9]).collect();
        let mut expected: Vec<u8> = seq.clone();
        expected.make_ascii_uppercase();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 50);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");

        let bytes: Vec<u8> = ChromRefFetcher::iter_bases(&fetcher)
            .expect("iter")
            .map(|b| b.expect("base"))
            .collect();
        assert_eq!(bytes, expected);
    }

    #[test]
    fn chrom_ref_fetch_after_iter_bases_starts_fresh_phase() {
        // The DUST + PerGroupMerger composition: iter_bases walks
        // the full contig forward, then fetch() is called near the
        // start of the contig. The iter_bases Drop must reset the
        // sliding buffer so this composition doesn't trip
        // OutOfPattern.
        let seq: Vec<u8> = (0..2_000_000_u32)
            .map(|i| b"ACGT"[(i % 4) as usize])
            .collect();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 80);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("fetcher");

        // Drive iter_bases to completion (consume + drop).
        {
            let it = ChromRefFetcher::iter_bases(&fetcher).expect("iter");
            let count = it.filter_map(Result::ok).count();
            assert_eq!(count, 2_000_000);
        }
        // Now a fetch near the contig's start must succeed without
        // OutOfPattern — iter_bases reset the buffer on Drop.
        let bytes = ChromRefFetcher::fetch(&fetcher, 100, 16).expect("post-iter fetch");
        assert_eq!(bytes, seq[99..115]);
    }

    #[test]
    fn walker_adapter_serves_each_chrom_via_swap() {
        // Build a two-chrom FASTA, construct the walker adapter,
        // alternate fetches across chroms. Every fetch returns the
        // right bytes; the adapter rebuilds the inner streamer on
        // each transition.
        let specs = vec![
            ContigSpec {
                name: "chr0".into(),
                length: 16,
            },
            ContigSpec {
                name: "chr1".into(),
                length: 24,
            },
        ];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let contigs = contig_list(&[("chr0", 16), ("chr1", 24)]);
        let adapter = WalkerLegacyAdapter::new(path, contigs);

        assert_eq!(
            RefSeqFetcher::fetch(&adapter, 0, 1, 4).expect("chr0"),
            b"AAAA"
        );
        assert_eq!(
            RefSeqFetcher::fetch(&adapter, 1, 1, 4).expect("chr1"),
            b"AAAA"
        );
        // Switch back: forces another swap.
        assert_eq!(
            RefSeqFetcher::fetch(&adapter, 0, 5, 4).expect("chr0 again"),
            b"AAAA"
        );
    }

    #[test]
    fn walker_adapter_rejects_unknown_chrom_id() {
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 8,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let contigs = contig_list(&[("chr0", 8)]);
        let adapter = WalkerLegacyAdapter::new(path, contigs);

        let err = RefSeqFetcher::fetch(&adapter, 99, 1, 4).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn walker_adapter_holds_one_chrom_at_a_time() {
        // Sequential fetches across two chroms; only one inner
        // StreamingChromRefFetcher is resident at any time. We
        // verify the swap happened by checking that fetching from
        // an out-of-bounds offset on the previously-loaded chrom
        // (after the swap) errors — confirming the inner is now
        // bound to the new chrom.
        let specs = vec![
            ContigSpec {
                name: "chr0".into(),
                length: 8,
            },
            ContigSpec {
                name: "chr1".into(),
                length: 32,
            },
        ];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let contigs = contig_list(&[("chr0", 8), ("chr1", 32)]);
        let adapter = WalkerLegacyAdapter::new(path, contigs);

        // Load chr0.
        RefSeqFetcher::fetch(&adapter, 0, 1, 4).expect("chr0 load");
        // Switch to chr1 — long fetch only works on chr1.
        let bytes = RefSeqFetcher::fetch(&adapter, 1, 1, 32).expect("chr1 full");
        assert_eq!(bytes.len(), 32);
        // Fetching chr1's full length on chr0 (8 bases) should error
        // (the inner is now bound to chr1, but legacy fetch chrom_id
        // = 0 forces a swap back to chr0 which then errors on the
        // out-of-bounds 32-base request).
        let err = RefSeqFetcher::fetch(&adapter, 0, 1, 32).expect_err("must fail");
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    // ---------------------------------------------------------------
    // ManualEvictChromRefFetcher tests
    // ---------------------------------------------------------------

    /// Helper for `ManualEvictChromRefFetcher` tests: build a small
    /// line-wrapped FASTA on disk with one contig of arbitrary bytes
    /// (already SAM-spec uppercase ACGT/N) and return a fetcher
    /// pointing at it plus the tempdir guard.
    fn manual_evict_fetcher_from_bytes(
        seq: &[u8],
        line_bases: usize,
    ) -> (tempfile::TempDir, ManualEvictChromRefFetcher) {
        let (dir, path) = build_line_wrapped_fasta("chr0", seq, line_bases);
        let fetcher =
            ManualEvictChromRefFetcher::for_contig(&path, "chr0").expect("manual-evict fetcher");
        (dir, fetcher)
    }

    #[test]
    fn manual_evict_fetch_grows_buffer_forward() {
        // Two forward fetches at increasing positions. The buffer
        // grows to cover the union of the two ranges; no eviction.
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 60);

        let first = fetcher.fetch(10, 8).expect("forward 1").to_vec();
        assert_eq!(first, &seq[9..17]);
        assert_eq!(fetcher.buf_start_base(), 10);
        let len_after_first = fetcher.buf_len();
        assert!(len_after_first >= 8);

        let second = fetcher.fetch(80, 8).expect("forward 2").to_vec();
        assert_eq!(second, &seq[79..87]);
        // Buffer grew forward to cover up to base 87.
        assert!(fetcher.buf_len() >= 87 - 10 + 1 - 1);
        // Buffer origin unchanged (no eviction).
        assert_eq!(fetcher.buf_start_base(), 10);
    }

    #[test]
    fn manual_evict_fetch_grows_buffer_backward() {
        // Fetch at a high position first, then at a lower position.
        // The buffer must prepend bytes to cover the lower range.
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 60);

        let high = fetcher.fetch(120, 8).expect("high").to_vec();
        assert_eq!(high, &seq[119..127]);
        assert_eq!(fetcher.buf_start_base(), 120);

        let low = fetcher.fetch(20, 8).expect("low").to_vec();
        assert_eq!(low, &seq[19..27]);
        // Buffer origin moved backward to 20.
        assert_eq!(fetcher.buf_start_base(), 20);
        // And still covers the high range (no eviction).
        assert!(fetcher.buf_len() >= 127 - 20 + 1 - 1);
    }

    #[test]
    fn manual_evict_evict_before_compacts_buffer() {
        // After a few forward fetches the buffer covers a wide
        // range. `evict_before` drops bytes before the given
        // coordinate and updates `buf_start_base`.
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 60);

        let _ = fetcher.fetch(10, 8).expect("fetch 1");
        let _ = fetcher.fetch(120, 8).expect("fetch 2");
        let before_evict = fetcher.buf_len();

        fetcher.evict_before(100);
        assert_eq!(fetcher.buf_start_base(), 100);
        // Buffer length shrank by 90 (bases 10..99 dropped).
        assert_eq!(before_evict - fetcher.buf_len(), 90);

        // Subsequent fetch in the retained range still works.
        let bytes = fetcher.fetch(120, 8).expect("post-evict fetch");
        assert_eq!(bytes, &seq[119..127]);
    }

    #[test]
    fn manual_evict_evict_before_below_origin_is_noop() {
        // `evict_before(pos)` when pos <= buf_start_base must not
        // touch the buffer (no negative drain).
        let seq: Vec<u8> = (0..50).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 50);

        let _ = fetcher.fetch(20, 8).expect("fetch");
        let origin = fetcher.buf_start_base();
        let len_before = fetcher.buf_len();

        fetcher.evict_before(origin);
        assert_eq!(fetcher.buf_start_base(), origin);
        assert_eq!(fetcher.buf_len(), len_before);

        fetcher.evict_before(origin - 1);
        assert_eq!(fetcher.buf_start_base(), origin);
        assert_eq!(fetcher.buf_len(), len_before);
    }

    #[test]
    fn manual_evict_evict_before_past_end_clears_buffer() {
        // Edge case: evict_before pointing past the buffer's end
        // clears the buffer entirely and pins `buf_start_base` to
        // the evict point.
        let seq: Vec<u8> = (0..50).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 50);

        let _ = fetcher.fetch(10, 8).expect("fetch");

        fetcher.evict_before(40);
        assert_eq!(fetcher.buf_len(), 0);
        assert_eq!(fetcher.buf_start_base(), 40);

        // Fresh fetch after the buffer was cleared.
        let bytes = fetcher.fetch(40, 8).expect("post-clear fetch");
        assert_eq!(bytes, &seq[39..47]);
    }

    #[test]
    fn manual_evict_fetch_uppercases_soft_masked() {
        // Soft-masked input → uppercased output, same as the other
        // fetchers.
        let seq = b"ACGTacgtNn".to_vec();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 10);
        let bytes = fetcher.fetch(1, seq.len() as u32).expect("fetch");
        assert_eq!(bytes, b"ACGTACGTNN");
    }

    #[test]
    fn manual_evict_fetch_past_contig_end_returns_out_of_bounds() {
        let seq = b"ACGTACGT".to_vec();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 8);

        let err = fetcher.fetch(6, 5).expect_err("must fail");
        match err {
            ChromRefFetchError::OutOfBounds { .. } => {}
            other => panic!("expected OutOfBounds, got {other:?}"),
        }
    }

    #[test]
    fn manual_evict_fetch_start_zero_is_rejected() {
        let seq = b"ACGT".to_vec();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 4);
        let err = fetcher.fetch(0, 1).expect_err("must fail");
        match err {
            ChromRefFetchError::InvalidStart => {}
            other => panic!("expected InvalidStart, got {other:?}"),
        }
    }

    #[test]
    fn manual_evict_full_lifecycle_simulates_baq_chunk_flow() {
        // Simulate the BAQ access pattern: a "chunk" of reads at
        // positions 10..50, processed in pseudo-random order (rayon
        // work-stealing would scramble them similarly). The
        // fetcher's buffer grows to cover the chunk's range; after
        // each "read" we call evict_before. End state: buffer ≤
        // chunk_max - chunk_min in size.
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 60);

        // Out-of-order reads in [10..=50].
        let read_positions = [30u32, 50, 12, 40, 22, 18, 48];

        for pos in read_positions {
            let bytes = fetcher.fetch(pos, 6).expect("fetch").to_vec();
            assert_eq!(bytes, &seq[(pos - 1) as usize..(pos - 1 + 6) as usize]);
            fetcher.evict_before(pos);
        }

        // After the full chunk: buffer origin advanced to the last
        // evicted position (48), buffer holds at most a small
        // window.
        assert_eq!(fetcher.buf_start_base(), 48);
        assert!(fetcher.buf_len() <= 8);
    }

    // ---------------------------------------------------------------
    // B1 (2026-05-23 code review): .fai validation tests
    // ---------------------------------------------------------------
    //
    // The .fai file format permits attacker-influenced values
    // (noodles_fasta::fai::Record parses without validating
    // `line_bases > 0` or `line_width >= line_bases`). The
    // `--reference` CLI argument is a real attacker surface.
    // Pre-B1 a malformed .fai with `line_bases = 0` would parse,
    // and the first `fetch` call would panic inside
    // ContigFai::base_to_file_offset on integer division by zero.

    /// Write `<fa>` and `<fai>` with caller-supplied .fai fields.
    /// Used to inject malformed .fai values that the noodles
    /// parser accepts but that the fetcher should reject.
    fn write_fixture_with_fai_line(
        contig_name: &str,
        contig_seq_len: u64,
        seq_offset: u64,
        line_bases: u64,
        line_width: u64,
    ) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("tempdir");
        let fa_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        // The .fa is mostly placeholder; the validation under
        // test fires at construction time (before any read).
        let mut fa = StdFile::create(&fa_path).expect("fa");
        write!(fa, ">{contig_name}\n").expect("fa header");
        fa.write_all(&vec![b'A'; contig_seq_len as usize])
            .expect("fa seq");
        fa.write_all(b"\n").expect("fa nl");
        let mut fai = StdFile::create(&fai_path).expect("fai");
        writeln!(
            fai,
            "{contig_name}\t{contig_seq_len}\t{seq_offset}\t{line_bases}\t{line_width}"
        )
        .expect("fai entry");
        (dir, fa_path)
    }

    #[test]
    fn streaming_for_contig_returns_io_error_on_zero_line_bases() {
        let (_dir, fa_path) =
            write_fixture_with_fai_line("chr0", 10, b">chr0\n".len() as u64, 0, 0);
        let result = StreamingChromRefFetcher::for_contig(&fa_path, "chr0");
        match result {
            Ok(_) => panic!("expected Io error, got Ok"),
            Err(ChromRefFetchError::Io { .. }) => {}
            Err(other) => panic!("expected Io error, got {other:?}"),
        }
    }

    #[test]
    fn streaming_for_contig_returns_io_error_on_line_width_less_than_line_bases() {
        // line_bases=10 line_width=5 — line_width must be >= line_bases
        // (line_width includes the trailing newline).
        let (_dir, fa_path) =
            write_fixture_with_fai_line("chr0", 100, b">chr0\n".len() as u64, 10, 5);
        let result = StreamingChromRefFetcher::for_contig(&fa_path, "chr0");
        match result {
            Ok(_) => panic!("expected Io error, got Ok"),
            Err(ChromRefFetchError::Io { .. }) => {}
            Err(other) => panic!("expected Io error, got {other:?}"),
        }
    }

    #[test]
    fn manual_evict_for_contig_returns_io_error_on_zero_line_bases() {
        let (_dir, fa_path) =
            write_fixture_with_fai_line("chr0", 10, b">chr0\n".len() as u64, 0, 0);
        let result = ManualEvictChromRefFetcher::for_contig(&fa_path, "chr0");
        match result {
            Ok(_) => panic!("expected Io error, got Ok"),
            Err(ChromRefFetchError::Io { .. }) => {}
            Err(other) => panic!("expected Io error, got {other:?}"),
        }
    }

    // ---------------------------------------------------------------
    // M26 (2026-05-23 code review): IUPAC + soft-mask canonicalisation
    // ---------------------------------------------------------------

    #[test]
    fn streaming_fetch_canonicalises_iupac_and_softmask_to_n_and_acgt() {
        // Mix of soft-masked, IUPAC-ambiguity, gap, RNA, and N bytes.
        // Expected: lowercase ACGT → uppercase; everything else → N.
        // Input bytes:   a c g t  R Y S W K M  B D H V N  -  .  U  X
        // Expected:      A C G T  N N N N N N  N N N N N  N  N  N  N
        let input = b"acgtRYSWKMBDHVN-.UX";
        let expected: Vec<u8> = vec![
            b'A', b'C', b'G', b'T', // lowercase ACGT folded
            b'N', b'N', b'N', b'N', b'N', b'N', // IUPAC 2-codes
            b'N', b'N', b'N', b'N', // IUPAC 3-codes
            b'N', // N stays N
            b'N', b'N', // gap chars
            b'N', // RNA U
            b'N', // any other ASCII
        ];
        let (_dir, fa_path) = build_line_wrapped_fasta("chr0", input, 8);
        let fetcher =
            StreamingChromRefFetcher::for_contig(&fa_path, "chr0").expect("for_contig");
        let bytes =
            ChromRefFetcher::fetch(&fetcher, 1, input.len() as u32).expect("fetch");
        assert_eq!(bytes, expected);
    }

    #[test]
    fn manual_evict_fetch_canonicalises_iupac_to_n() {
        let input = b"acgtRYSWKMBDHVN";
        let expected: Vec<u8> = b"ACGTNNNNNNNNNNN".to_vec();
        let (_dir, fa_path) = build_line_wrapped_fasta("chr0", input, 8);
        let mut fetcher =
            ManualEvictChromRefFetcher::for_contig(&fa_path, "chr0").expect("for_contig");
        let bytes = fetcher.fetch(1, input.len() as u32).expect("fetch").to_vec();
        assert_eq!(bytes, expected);
    }
}
