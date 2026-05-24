//! Reference-FASTA fetchers. Three types live here:
//!
//! - [`StreamingChromRefFetcher`] — per-contig, sliding-buffer reader.
//!   The workhorse for monotonic-forward access (cohort var-calling,
//!   PerGroupMerger). Implements the single-contig typed-error
//!   [`ChromRefFetcher`] trait. Bound to one contig at construction;
//!   serves `fetch(start, length)` and `iter_bases()` through a 1 MB
//!   sliding buffer.
//! - [`ManualEvictChromRefFetcher`] — per-contig, caller-managed
//!   buffer lifecycle. Designed for genuinely random-access patterns
//!   within a known region (Stage 1 BAQ per-worker instance).
//! - [`MultiChromStreamingRefFetcher`] — multi-contig fetcher used by
//!   the Stage 1 pileup walker. Wraps a swappable
//!   `StreamingChromRefFetcher` and rebuilds it on chrom transition.
//!   Implements the multi-contig typed-error
//!   [`super::pileup::MultiChromRefFetcher`] trait. `Sync` via an
//!   internal `Mutex` (the only `Send + Sync` fetcher in the module).
//!
//! M12 (2026-05-23 code review): the legacy `RefSeqFetcher` trait
//! (`io::Error`-returning, multi-chrom) was retired. All consumers
//! migrated to the typed-error
//! [`super::pileup::MultiChromRefFetcher`] trait.

use std::cell::RefCell;
use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use noodles_fasta::fai;

use super::{ContigList, MultiChromRefFetcher};

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

/// Single-contig sliding-buffer FASTA reader. Implements the typed-
/// error [`ChromRefFetcher`] trait. Designed for the cohort
/// `var-calling` per-chromosome workers, where access to the contig
/// is sequential and monotonically non-decreasing in position:
///
/// - DUST mask construction does **one** forward sequential pass over
///   the contig (see [`ChromRefFetcher::iter_bases`]).
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
/// **Contract.** Access must be monotonically non-decreasing in
/// position. A backward jump beyond the sliding buffer returns
/// [`ChromRefFetchError::OutOfPattern`]. Production callers
/// (PerGroupMerger, DUST) satisfy this by construction; the
/// [`MultiChromStreamingRefFetcher`] (the walker-side multi-contig
/// wrapper) rebuilds the inner fetcher on chrom transition, so
/// chrom-boundary backward jumps don't violate the contract either.
pub struct StreamingChromRefFetcher {
    /// Contig name for error messages only.
    chrom_name: String,
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
    fn validate(&self, chrom_name: &str) -> Result<(), io::Error> {
        if self.line_bases == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "malformed .fai for contig {chrom_name}: line_bases = 0 \
                     (would divide-by-zero in offset arithmetic)"
                ),
            ));
        }
        if self.line_width < self.line_bases {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "malformed .fai for contig {chrom_name}: line_width ({}) \
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

struct StreamState {
    /// Production reads from a `File`; tests reach the same code path
    /// through real on-disk fixtures (`build_fasta` /
    /// `build_line_wrapped_fasta` in the test module), so we don't
    /// carry an in-memory variant. The previous `enum Source { File,
    /// Memory }` was retired when M12 deleted the legacy
    /// `from_memory` test helper.
    source: File,
    /// Uppercased, newline-stripped bases currently resident.
    /// Length ≤ `STREAMING_REF_BUFFER_BYTES`.
    buf: Vec<u8>,
    /// 1-based contig coordinate of `buf[0]`. Bases held in `buf`
    /// cover contig positions `[buf_start_base, buf_start_base + buf.len())`.
    /// Zero when no buffer has been filled yet.
    buf_start_base: u32,
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
    state.source.seek(SeekFrom::Start(offset))?;
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
        let n = state.source.read(&mut read_buf)?;
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
// The trait below is the single-contig typed-error fetcher API
// used by per-chrom consumers (DUST, PerGroupMerger). The
// multi-contig sibling is [`super::pileup::MultiChromRefFetcher`]
// (used by the Stage 1 walker). Both use [`ChromRefFetchError`].
//
// Design choices:
//
// - **No `chrom_id` parameter.** The fetcher is bound to one contig
//   at construction; the multi-contig wrapper
//   [`MultiChromStreamingRefFetcher`] rebuilds the inner fetcher on
//   chrom transition rather than passing chrom ids per call.
// - **Typed error.** `ChromRefFetchError` distinguishes "I/O
//   failure", "out of contig bounds", and "out of access pattern"
//   so consumers can route them differently.
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
    #[error("fetch [{start}, {end}) past contig {chrom_name} length {chrom_length}")]
    OutOfBounds {
        chrom_name: String,
        chrom_length: u32,
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
    #[error("I/O failure on contig {chrom_name}: {source}")]
    Io {
        chrom_name: String,
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
    /// Constructor bound to one contig by name. Opens the FASTA,
    /// parses the sibling `.fai`, looks up the contig, and returns
    /// a fetcher ready to serve `fetch` calls. The contig's bytes
    /// are **not** loaded here — the sliding buffer fills lazily on
    /// the first fetch.
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
    ///   [`io::ErrorKind::NotFound`] error if `chrom_name` is not
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
    pub fn for_contig(fasta_path: &Path, chrom_name: &str) -> Result<Self, ChromRefFetchError> {
        Self::for_contig_with_fai_path(fasta_path, None, chrom_name)
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
        chrom_name: &str,
    ) -> Result<Self, ChromRefFetchError> {
        Self::for_contig_internal(fasta_path, fai_path, chrom_name, STREAMING_REF_BUFFER_BYTES)
    }

    /// Test-only constructor allowing a custom buffer size so the
    /// `OutOfPattern` contract can be exercised with small fixtures
    /// (production code only ever sees `STREAMING_REF_BUFFER_BYTES`).
    #[cfg(test)]
    pub fn for_contig_with_buffer_bytes(
        fasta_path: &Path,
        chrom_name: &str,
        buffer_bytes: usize,
    ) -> Result<Self, ChromRefFetchError> {
        Self::for_contig_internal(fasta_path, None, chrom_name, buffer_bytes)
    }

    fn for_contig_internal(
        fasta_path: &Path,
        fai_path: Option<&Path>,
        chrom_name: &str,
        buffer_bytes: usize,
    ) -> Result<Self, ChromRefFetchError> {
        let (fai, file) = open_contig(fasta_path, fai_path, chrom_name)?;
        Ok(Self {
            chrom_name: chrom_name.to_string(),
            fai,
            inner: RefCell::new(StreamState {
                source: file,
                buf: Vec::with_capacity(buffer_bytes),
                buf_start_base: 0,
            }),
        })
    }
}

/// M13: shared `.fai`-open + index-load + record-lookup +
/// `u32::try_from` + validation + `File::open` path used by both
/// fetcher constructors. Returning `(ContigFai, File)` lets each
/// caller wrap into its own state struct.
///
/// `fai_path = None` falls back to `<fasta_path>.fai`. The B1
/// validation (`line_bases > 0`, `line_width >= line_bases`,
/// no `u32::try_from` overflow) lives here, so every constructor
/// inherits it.
fn open_contig(
    fasta_path: &Path,
    fai_path: Option<&Path>,
    chrom_name: &str,
) -> Result<(ContigFai, File), ChromRefFetchError> {
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
    let wrap_io = |source: io::Error| ChromRefFetchError::Io {
        chrom_name: chrom_name.to_string(),
        source,
    };
    let index = fai::fs::read(fai_path_resolved).map_err(wrap_io)?;
    let record = index
        .as_ref()
        .iter()
        .find(|r| AsRef::<[u8]>::as_ref(r.name()) == chrom_name.as_bytes())
        .ok_or_else(|| {
            wrap_io(io::Error::new(
                io::ErrorKind::NotFound,
                format!("contig {chrom_name} not in FASTA index"),
            ))
        })?;
    let length = u32::try_from(record.length()).map_err(|_| {
        wrap_io(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "contig {} length {} exceeds u32::MAX",
                chrom_name,
                record.length()
            ),
        ))
    })?;
    let line_bases = u32::try_from(record.line_bases()).map_err(|_| {
        wrap_io(io::Error::new(
            io::ErrorKind::InvalidData,
            ".fai line_bases exceeds u32::MAX",
        ))
    })?;
    let line_width = u32::try_from(record.line_width()).map_err(|_| {
        wrap_io(io::Error::new(
            io::ErrorKind::InvalidData,
            ".fai line_width exceeds u32::MAX",
        ))
    })?;
    let fai = ContigFai {
        seq_offset: record.offset(),
        length,
        line_bases,
        line_width,
    };
    // B1: reject malformed .fai (line_bases = 0, etc.) before any
    // fetch can panic on the offset arithmetic.
    fai.validate(chrom_name).map_err(wrap_io)?;
    let file = File::open(fasta_path).map_err(wrap_io)?;
    Ok((fai, file))
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
                    chrom_name: self.fetcher.chrom_name.clone(),
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
                    chrom_name: self.chrom_name.clone(),
                    source: io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"),
                })?;
        if end_1based.saturating_sub(1) > self.fai.length {
            return Err(ChromRefFetchError::OutOfBounds {
                chrom_name: self.chrom_name.clone(),
                chrom_length: self.fai.length,
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
                    chrom_name: self.chrom_name.clone(),
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
        // Mi5: scope the borrow_mut so the RefCell guard is released
        // before we return the iterator (which the caller will then
        // poll, each .next() taking a fresh borrow_mut). Avoids an
        // accidental long-lived borrow if this function ever grows
        // post-reset work.
        {
            let mut state = self.inner.borrow_mut();
            state.buf.clear();
            state.buf_start_base = 0;
        }
        Ok(Box::new(ChromRefBaseIter {
            fetcher: self,
            next_base: 1,
            done: false,
        }))
    }
}

// ---------------------------------------------------------------------
// MultiChromStreamingRefFetcher — multi-chrom, single-thread bridge.
// ---------------------------------------------------------------------
//
// Bridge for consumers that traverse multiple chromosomes in
// sequence — the Stage 1 pileup walker today. On each chrom
// transition the inner [`StreamingChromRefFetcher`] is rebuilt for
// the new contig; fetches within a chrom hit the same sliding
// buffer.

/// Multi-contig FASTA fetcher used by the Stage 1 pileup walker.
/// Maintains *one* [`StreamingChromRefFetcher`] at a time, rebuilding
/// it on the first `fetch` call for a different `chrom_id`.
/// Implements [`MultiChromRefFetcher`], so it drops into the walker's
/// generic plumbing.
///
/// `Sync` via internal `Mutex<Option<...>>`. The walker itself is
/// single-threaded so the Mutex is uncontended; `var_calling_from_bam`
/// also uses this fetcher and passes it through `Arc<dyn
/// MultiChromRefFetcher + Send + Sync>`, so the `Sync` bound is real.
/// Stage 1's BAQ uses a separate per-rayon-worker
/// [`ManualEvictChromRefFetcher`] (built inside `BaqStream`'s
/// `map_init`); BAQ does **not** share a fetcher with the walker
/// — the access patterns are too different for one abstraction
/// to serve both.
///
/// Behaviour on chrom transition: a fresh
/// [`StreamingChromRefFetcher::for_contig`] is built for the new
/// chrom *outside* the lock (M16 of the 2026-05-23 code review —
/// the inner-fetcher build does the .fai parse + `File::open` work,
/// which is millisecond-scale and shouldn't serialise other walker
/// threads). The old fetcher is dropped when the new one is
/// installed; both are briefly resident across the swap.
/// Within-chrom fetches delegate to the strict
/// [`ChromRefFetcher::fetch`] on the inner; the walker's
/// single-threaded sequential access pattern satisfies the
/// monotonic-forward contract.
pub struct MultiChromStreamingRefFetcher {
    fasta_path: PathBuf,
    contigs: ContigList,
    inner: Mutex<Option<(u32, StreamingChromRefFetcher)>>,
}

impl MultiChromStreamingRefFetcher {
    /// Build an adapter bound to a FASTA. No inner fetcher is
    /// constructed yet — the first `fetch` triggers the lazy load.
    pub fn new(fasta_path: PathBuf, contigs: ContigList) -> Self {
        Self {
            fasta_path,
            contigs,
            inner: Mutex::new(None),
        }
    }

    /// M4: surface poison as a typed error rather than swallowing it.
    /// The walker is single-threaded inside a chrom worker; the only
    /// way to poison this mutex is for the inner
    /// `StreamingChromRefFetcher::for_contig` rebuild to panic (which
    /// panic-frees the critical section, leaving `inner = None`).
    /// Returning here lets the caller terminate the chrom worker
    /// rather than silently continuing with a possibly-inconsistent
    /// inner state.
    fn lock_inner(
        &self,
    ) -> Result<
        std::sync::MutexGuard<'_, Option<(u32, StreamingChromRefFetcher)>>,
        ChromRefFetchError,
    > {
        self.inner.lock().map_err(|_poison| ChromRefFetchError::Io {
            chrom_name: String::from("<walker-adapter>"),
            source: io::Error::other(
                "MultiChromStreamingRefFetcher mutex poisoned (prior panic during rebuild)",
            ),
        })
    }
}

// M23 (2026-05-23 code review): MultiChromStreamingRefFetcher is the only
// `Send + Sync` fetcher in the module — Stage 1 reader threads share
// it through `Arc<dyn MultiChromRefFetcher + Send + Sync>`. Compile-
// time assertion guards both bounds; a future field change that
// flipped either would silently break the walker.
const _: () = {
    const fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<MultiChromStreamingRefFetcher>();
};

impl MultiChromRefFetcher for MultiChromStreamingRefFetcher {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError> {
        // M16: rebuild outside the lock. The lock is taken twice (once
        // to check, once to install), but each critical section is
        // tiny (a pointer comparison or a single move). Previously the
        // lock was held across `for_contig` (file open + .fai parse +
        // contig lookup), serialising every walker thread through a
        // millisecond-scale operation.
        let needs_rebuild = {
            let slot = self.lock_inner()?;
            !matches!(slot.as_ref(), Some((current, _)) if *current == chrom_id)
        };

        if needs_rebuild {
            let entry = self.contigs.entries.get(chrom_id as usize).ok_or_else(|| {
                ChromRefFetchError::Io {
                    chrom_name: String::from("<walker-adapter>"),
                    source: io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "MultiChromStreamingRefFetcher: chrom_id {chrom_id} out of range \
                             (have {} contigs)",
                            self.contigs.entries.len()
                        ),
                    ),
                }
            })?;
            // Build outside the lock (the expensive step: .fai parse
            // + File::open + ContigFai::validate). The old fetcher
            // for the previous chrom is still resident in the slot
            // here, so peak memory briefly covers two contig
            // working sets — accepted to keep the lock-holding
            // window minimal.
            let fresh = StreamingChromRefFetcher::for_contig(&self.fasta_path, &entry.name)?;
            // Install the new fetcher; the old one drops here.
            let mut slot = self.lock_inner()?;
            *slot = Some((chrom_id, fresh));
        }

        let slot = self.lock_inner()?;
        // PANIC-FREE: when `needs_rebuild` was true, the block above
        // installed `Some((chrom_id, _))` before releasing the lock;
        // when it was false, the slot was already
        // `Some((chrom_id, _))` from the matches! check. A concurrent
        // racer could have rebuilt for a *different* chrom_id between
        // our install and this re-lock — in that case we'd fall
        // through, but the walker is single-threaded per chrom so
        // this can't happen in production. Returning a typed error
        // rather than panicking on the race is defence in depth.
        let (_, inner) = match slot.as_ref() {
            Some(pair) if pair.0 == chrom_id => pair,
            _ => {
                return Err(ChromRefFetchError::Io {
                    chrom_name: String::from("<walker-adapter>"),
                    source: io::Error::other(format!(
                        "MultiChromStreamingRefFetcher: slot raced — expected chrom_id {chrom_id} \
                         but found a different binding"
                    )),
                });
            }
        };
        ChromRefFetcher::fetch(inner, start_1based, length)
    }

    // iter_bases falls back to the trait default. The walker doesn't
    // use it (only DUST / PerGroupMerger in the var-calling pipeline
    // do), so defining it explicitly would be dead code.
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
    chrom_name: String,
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
    /// look up `chrom_name`, return a fetcher with an empty
    /// buffer ready to serve `fetch` calls.
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::Io`] with `ErrorKind::NotFound` if
    ///   `chrom_name` is absent from the FASTA index.
    /// - [`ChromRefFetchError::Io`] with `ErrorKind::InvalidData` if
    ///   the `.fai` entry is malformed (zero `line_bases`,
    ///   `line_width < line_bases`, lengths overflowing `u32`).
    /// - [`ChromRefFetchError::Io`] propagating the underlying
    ///   `io::Error` if the `.fai` read or `File::open(fasta_path)`
    ///   fails.
    pub fn for_contig(fasta_path: &Path, chrom_name: &str) -> Result<Self, ChromRefFetchError> {
        Self::for_contig_with_fai_path(fasta_path, None, chrom_name)
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
        chrom_name: &str,
    ) -> Result<Self, ChromRefFetchError> {
        let (fai, file) = open_contig(fasta_path, fai_path, chrom_name)?;
        Ok(Self {
            chrom_name: chrom_name.to_string(),
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
                    chrom_name: self.chrom_name.clone(),
                    source: io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"),
                })?;
        // [start_1based, end_exclusive); valid only if end <= length+1.
        if end_exclusive.saturating_sub(1) > self.fai.length {
            return Err(ChromRefFetchError::OutOfBounds {
                chrom_name: self.chrom_name.clone(),
                chrom_length: self.fai.length,
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
                chrom_name: self.chrom_name.clone(),
                source,
            })?;
        let remaining = (self.fai.length - start_1based + 1) as usize;
        let target_len = n_bases.min(remaining);
        self.buf.reserve(target_len);
        self.buf_start_base = start_1based;
        read_uppercased_bases(&mut self.file, &mut self.buf, target_len).map_err(|source| {
            ChromRefFetchError::Io {
                chrom_name: self.chrom_name.clone(),
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
                chrom_name: self.chrom_name.clone(),
                source,
            })?;
        let remaining = (self.fai.length - next_base + 1) as usize;
        let take = extra_bases.min(remaining);
        self.buf.reserve(take);
        read_uppercased_bases(&mut self.file, &mut self.buf, take).map_err(|source| {
            ChromRefFetchError::Io {
                chrom_name: self.chrom_name.clone(),
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
                chrom_name: self.chrom_name.clone(),
                source,
            })?;
        let mut prefix = Vec::with_capacity(extra_bases);
        read_uppercased_bases(&mut self.file, &mut prefix, extra_bases).map_err(|source| {
            ChromRefFetchError::Io {
                chrom_name: self.chrom_name.clone(),
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

    use crate::fasta::{ContigEntry, ContigList};
    use crate::per_sample_pileup::cram_files::{ContigSpec, build_fasta};

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

    // ---------------------------------------------------------------
    // ChromRefFetcher (new API) tests
    // ---------------------------------------------------------------
    //
    // These tests exercise the new-API methods: `for_contig`
    // construction, `length()`, `fetch(start, length)` (without
    // chrom_id), `iter_bases()` (with buffer reset semantics), and
    // the `OutOfPattern` contract on backward jumps beyond the
    // buffer.

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
                chrom_length,
                start,
                end,
                ..
            } => {
                assert_eq!(chrom_length, 10);
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
            Err(ChromRefFetchError::Io { chrom_name, source }) => {
                assert_eq!(chrom_name, "chr_missing");
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
        let adapter = MultiChromStreamingRefFetcher::new(path, contigs);

        assert_eq!(
            MultiChromRefFetcher::fetch(&adapter, 0, 1, 4).expect("chr0"),
            b"AAAA"
        );
        assert_eq!(
            MultiChromRefFetcher::fetch(&adapter, 1, 1, 4).expect("chr1"),
            b"AAAA"
        );
        // Switch back: forces another swap.
        assert_eq!(
            MultiChromRefFetcher::fetch(&adapter, 0, 5, 4).expect("chr0 again"),
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
        let adapter = MultiChromStreamingRefFetcher::new(path, contigs);

        let err = MultiChromRefFetcher::fetch(&adapter, 99, 1, 4).expect_err("must fail");
        match err {
            ChromRefFetchError::Io { source, .. } => {
                assert_eq!(source.kind(), io::ErrorKind::InvalidInput);
            }
            other => panic!("expected Io(InvalidInput), got {other:?}"),
        }
    }

    #[test]
    fn walker_adapter_holds_one_chrom_at_a_time() {
        // Sequential fetches across two chroms; only one inner
        // StreamingChromRefFetcher is resident at any time. After a
        // swap to chr1 the adapter's `inner` slot is bound to chr1.
        // Fetching chr0 again forces another rebuild; an out-of-bounds
        // request on chr0 (8 bases) for 32 must surface OutOfBounds.
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
        let adapter = MultiChromStreamingRefFetcher::new(path, contigs);

        // Load chr0.
        MultiChromRefFetcher::fetch(&adapter, 0, 1, 4).expect("chr0 load");
        // Switch to chr1 — long fetch only works on chr1.
        let bytes = MultiChromRefFetcher::fetch(&adapter, 1, 1, 32).expect("chr1 full");
        assert_eq!(bytes.len(), 32);
        // After switching back to chr0, a 32-base fetch must error
        // (chr0 is only 8 bases).
        let err = MultiChromRefFetcher::fetch(&adapter, 0, 1, 32).expect_err("must fail");
        assert!(
            matches!(err, ChromRefFetchError::OutOfBounds { .. }),
            "expected OutOfBounds, got {err:?}",
        );
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
        chrom_name: &str,
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
        writeln!(fa, ">{chrom_name}").expect("fa header");
        fa.write_all(&vec![b'A'; contig_seq_len as usize])
            .expect("fa seq");
        fa.write_all(b"\n").expect("fa nl");
        let mut fai = StdFile::create(&fai_path).expect("fai");
        writeln!(
            fai,
            "{chrom_name}\t{contig_seq_len}\t{seq_offset}\t{line_bases}\t{line_width}"
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
        let fetcher = StreamingChromRefFetcher::for_contig(&fa_path, "chr0").expect("for_contig");
        let bytes = ChromRefFetcher::fetch(&fetcher, 1, input.len() as u32).expect("fetch");
        assert_eq!(bytes, expected);
    }

    #[test]
    fn manual_evict_fetch_canonicalises_iupac_to_n() {
        let input = b"acgtRYSWKMBDHVN";
        let expected: Vec<u8> = b"ACGTNNNNNNNNNNN".to_vec();
        let (_dir, fa_path) = build_line_wrapped_fasta("chr0", input, 8);
        let mut fetcher =
            ManualEvictChromRefFetcher::for_contig(&fa_path, "chr0").expect("for_contig");
        let bytes = fetcher
            .fetch(1, input.len() as u32)
            .expect("fetch")
            .to_vec();
        assert_eq!(bytes, expected);
    }

    // ---------------------------------------------------------------
    // B3 + M22 + M24 + M25 (2026-05-23 code review): contract /
    // edge-case tests that were missing.
    // ---------------------------------------------------------------

    /// Stub [`ChromRefFetcher`] for M22: only overrides the required
    /// methods (`length`, `fetch`, `iter_bases`). Inherits the
    /// **default** `fetch_into` impl from the trait — the whole
    /// point of the M22 tests below is to exercise that default.
    struct DefaultFetchIntoStub {
        payload: Vec<u8>,
    }
    impl sealed::Sealed for DefaultFetchIntoStub {}
    impl ChromRefFetcher for DefaultFetchIntoStub {
        fn length(&self) -> u32 {
            self.payload.len() as u32
        }
        fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError> {
            if start_1based == 0 {
                return Err(ChromRefFetchError::InvalidStart);
            }
            let start = (start_1based - 1) as usize;
            let end = start + length as usize;
            Ok(self.payload[start..end].to_vec())
        }
        fn iter_bases<'a>(
            &'a self,
        ) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>
        {
            Ok(Box::new(self.payload.iter().copied().map(Ok)))
        }
        // No fetch_into override: tests below pin the trait's default
        // impl (mem::swap, dst-clearing).
    }

    #[test]
    fn fetch_into_default_impl_clears_dst() {
        // M22: the trait's default fetch_into impl must clear `dst`
        // before writing. A regression that prepended stale bytes
        // would silently corrupt every consumer relying on the
        // default impl.
        let stub = DefaultFetchIntoStub {
            payload: b"ACGTACGTAC".to_vec(),
        };
        let mut dst = b"GARBAGE_PRE_EXISTING_BYTES".to_vec();
        stub.fetch_into(3, 4, &mut dst).expect("fetch_into");
        assert_eq!(dst, b"GTAC");
        assert_eq!(dst.len(), 4);
    }

    #[test]
    fn fetch_into_through_ref_forwarding_clears_dst() {
        // M22: same contract via the `&T` blanket impl. The cohort
        // driver passes `&*fetcher` into DustFilter::new, which is
        // generic over `F: ChromRefFetcher`, so this seam is on the
        // hot path.
        let stub = DefaultFetchIntoStub {
            payload: b"ACGTACGTAC".to_vec(),
        };
        let by_ref: &dyn ChromRefFetcher = &stub;
        let mut dst = b"STALE".to_vec();
        by_ref.fetch_into(1, 4, &mut dst).expect("fetch_into");
        assert_eq!(dst, b"ACGT");
        assert_eq!(dst.len(), 4);
    }

    #[test]
    fn streaming_fetch_into_returns_exact_length_spanning_refill() {
        // M22: the StreamingChromRefFetcher override copies straight
        // from the slab. Test that the slab path also keeps the
        // `dst.len() == length` invariant when the request spans the
        // initial fill (no warm cache).
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, path) = build_line_wrapped_fasta("chr0", &seq, 60);
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("for_contig");
        let mut dst = b"GARBAGE".to_vec();
        ChromRefFetcher::fetch_into(&fetcher, 50, 32, &mut dst).expect("fetch_into");
        assert_eq!(dst.len(), 32);
        assert_eq!(dst, &seq[49..81]);
    }

    #[test]
    fn streaming_fetch_at_production_buffer_size_rejects_backward_jump() {
        // B3: build a contig wider than the production sliding
        // buffer (1 MiB) and exercise the OutOfPattern contract at
        // the real buffer size, not the test-only 4 KiB. A refactor
        // that silently turned OutOfPattern into a refill would
        // shrink throughput on the cohort var-caller hot path but
        // pass every other test in this file.
        let specs = vec![ContigSpec {
            name: "chr0".into(),
            length: 4 * 1024 * 1024,
        }];
        let (_dir, path) = build_fasta(&specs).expect("fasta");
        let fetcher = StreamingChromRefFetcher::for_contig(&path, "chr0").expect("for_contig");
        // Forward fetch at 2 MiB — primes the buffer past its origin.
        ChromRefFetcher::fetch(&fetcher, 2 * 1024 * 1024, 16).expect("forward fetch");
        // Backward jump beyond the buffer's origin.
        let err = ChromRefFetcher::fetch(&fetcher, 1, 16).expect_err("backward must fail");
        match err {
            ChromRefFetchError::OutOfPattern { .. } => {}
            other => panic!("expected OutOfPattern, got {other:?}"),
        }
    }

    #[test]
    fn manual_evict_fetch_in_buffer_hit_does_not_refill() {
        // M25: the warm-cache happy path for BAQ. After the first
        // fetch primes the buffer, a request whose range is already
        // covered must not touch `buf_start_base` and must not
        // re-read from disk (verified by observing buf_len + origin
        // are unchanged).
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 60);

        // Prime: fetch a 32-byte window starting at base 50.
        let primed = fetcher.fetch(50, 32).expect("prime").to_vec();
        assert_eq!(primed, &seq[49..81]);
        let origin_after_prime = fetcher.buf_start_base();
        let len_after_prime = fetcher.buf_len();

        // In-buffer hit: fetch a sub-range covered by the primed
        // window. Must not change buffer state.
        let hit = fetcher.fetch(60, 8).expect("in-buffer").to_vec();
        assert_eq!(hit, &seq[59..67]);
        assert_eq!(fetcher.buf_start_base(), origin_after_prime);
        assert_eq!(fetcher.buf_len(), len_after_prime);
    }

    #[test]
    fn manual_evict_evict_before_past_contig_end_clears_and_pins_buffer() {
        // M24: `evict_before(chrom_length + 1)` is a legitimate
        // caller intent (the BAQ chunk ran off the end of the
        // contig and the caller wants to free the buffer). The
        // buffer must clear and `buf_start_base` must pin to the
        // requested coordinate; a subsequent fetch from anywhere
        // valid in the contig still works.
        let seq: Vec<u8> = (0..50).map(|i| b"ACGT"[i % 4]).collect();
        let (_dir, mut fetcher) = manual_evict_fetcher_from_bytes(&seq, 50);

        fetcher.fetch(20, 8).expect("prime");
        assert!(fetcher.buf_len() > 0);

        let past_end = (seq.len() as u32) + 1; // == chrom_length + 1
        fetcher.evict_before(past_end);
        assert_eq!(fetcher.buf_len(), 0);
        assert_eq!(fetcher.buf_start_base(), past_end);

        // Buffer is empty: next fetch starts fresh from disk,
        // ignoring the stale `buf_start_base = past_end`.
        let bytes = fetcher.fetch(10, 8).expect("post-evict fetch").to_vec();
        assert_eq!(bytes, &seq[9..17]);
    }
}
