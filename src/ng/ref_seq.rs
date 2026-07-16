//! Reference-sequence access — the [`RefSeq`] trait and its implementations. The shared
//! source of reference-genome bases for the pileup, BAQ, DUST, left-alignment, and read
//! filtering. Design spec: `doc/devel/ng/spec/ref_seq.md`.
//!
//! [`RefSeq`] is the **universal** surface: canonical `{A,C,G,T,N}` fetch, which every
//! implementation provides. Raw bytes are the [`RawRefSeq`] sub-trait; buffer eviction is
//! an inherent method on the windowed impl ([`WindowedRefSeq`]). Both fetch surfaces write
//! into a caller-owned buffer (`&self`, alloc-free when the buffer is reused) rather than
//! returning a borrowed slice — that keeps every impl `&self`-shareable and avoids the
//! lazy-load/borrow-escape problem the file-backed impls would otherwise hit.
//!
//! This file holds the traits, the error type, and all three implementations: the
//! synthetic [`InMemoryRefSeq`], the whole-contig FASTA-backed [`ResidentRefSeq`], and the
//! streaming, caller-evictable [`WindowedRefSeq`]. Canonicalisation reuses the production
//! [`crate::fasta::fetcher::canonicalise`], so canonical bytes are byte-identical to the
//! production fetchers by construction.
//!
//! `noodles`' FASTA crate is referred to by its full name `noodles_fasta` here (only the
//! `Repository` type is named); `crate::fasta` is our own module. Keeping the two spelled
//! differently avoids the "same token, two crates" ambiguity.

use crate::fasta::fetcher::canonicalise;
use crate::fasta::{ChromRefFetchError, ContigEntry, ContigList, ManualEvictChromRefFetcher};
use crate::ng::types::ContigId;
use noodles_fasta::Repository;
use std::cell::RefCell;
use std::io;
use std::ops::Range;
use std::path::PathBuf;

/// Errors from a reference fetch. `#[non_exhaustive]` so matchers must accept future
/// variants; the shape mirrors the production `fasta::fetcher::ChromRefFetchError`.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RefSeqError {
    /// The requested `[start, start + length)` window exceeds the contig's length.
    #[error("fetch [{start}, {end}) past {contig:?} (length {contig_length})")]
    OutOfBounds {
        contig: ContigId,
        contig_length: u64,
        start: u64,
        end: u64,
    },
    /// `start_1based` was 0 — the 1-based coordinate contract was violated.
    #[error("start_1based must be >= 1")]
    InvalidStart,
    /// No contig with this id exists in the reference contig table.
    #[error("unknown {0:?}")]
    UnknownContig(ContigId),
    /// A reference read failed: either a genuine FASTA I/O error, or a contig that is in
    /// the table but absent from the FASTA repository. Both fold into one variant to
    /// mirror the production `ChromRefFetchError::Io` (kept for port-back parity); the
    /// `source` distinguishes them (a synthesised `NotFound` for the missing-contig case).
    #[error("reference read failed on {contig:?}")]
    Io {
        contig: ContigId,
        #[source]
        source: io::Error,
    },
    // NOTE: a strict monotonic-forward (auto-slide) fetcher would add an `OutOfPattern`
    // variant here (cf. production `ChromRefFetchError`). ng's windowed impl is
    // caller-evictable and bounds memory by eviction, so it does not need one;
    // `#[non_exhaustive]` keeps adding it later non-breaking. (ref_seq.md, Decision 6.)
}

/// Validate a 1-based `[start_1based, start_1based + length)` request against a contig of
/// `contig_len` bytes, returning the 0-based [`Range`] to slice. The single home for the
/// coordinate/bounds contract, shared by every impl so they cannot drift. The caller must
/// have already established that `contig` exists and that `start_1based >= 1`.
fn validate_window(
    contig: ContigId,
    contig_len: usize,
    start_1based: u64,
    length: u64,
) -> Result<Range<usize>, RefSeqError> {
    debug_assert!(
        start_1based >= 1,
        "caller must reject start_1based == 0 first"
    );
    let start0 = (start_1based - 1) as usize;
    // `checked_add` guards a 32-bit `usize`, where `start0 + length` could wrap and slip
    // past the bounds check into a slice-index panic; on 64-bit it never fails.
    match start0.checked_add(length as usize) {
        Some(end0) if end0 <= contig_len => Ok(start0..end0),
        _ => Err(RefSeqError::OutOfBounds {
            contig,
            // `u64`, so this reports the contig's real length. It used to be
            // `u32::try_from(contig_len).unwrap_or(u32::MAX)` — a >4 Gb contig
            // **silently clamped**, in built ng code, and the error then lied about
            // the very number it existed to report (spec §4). Widening deletes the
            // clamp rather than guarding it.
            contig_length: contig_len as u64,
            start: start_1based,
            end: start_1based.saturating_add(length),
        }),
    }
}

/// Narrow ng's `u64` request to the `u32` production's `ChromRefFetcher` takes.
///
/// **The one place ng's width meets production's, and it fails rather than folds.**
/// `src/fasta/` is frozen (spec Revision), so ng cannot widen it; but a silent
/// narrowing here would reintroduce exactly the bug B2 removed a few lines up —
/// the `unwrap_or(u32::MAX)` that turned a >4 Gb contig into a wrong number with
/// no error. So a request this fetcher cannot express is `OutOfBounds`, reported
/// against the contig's real (`u64`) length.
///
/// Unreachable on any real assembly — the largest known contig is ~250 Mb, and
/// `u32` holds 4.29 Gb — which is *why* it must be an error and not an assert: it
/// is a property of the reference, not a caller bug.
fn narrow_for_fetcher(
    contig: ContigId,
    entry: &ContigEntry,
    start_1based: u64,
    length: u64,
) -> Result<(u32, u32), RefSeqError> {
    match (u32::try_from(start_1based), u32::try_from(length)) {
        (Ok(s), Ok(l)) => Ok((s, l)),
        _ => Err(RefSeqError::OutOfBounds {
            contig,
            contig_length: entry.length,
            start: start_1based,
            end: start_1based.saturating_add(length),
        }),
    }
}

/// Access to reference-genome bases by (contig, 1-based range). The universal surface
/// every implementation provides: canonical `{A,C,G,T,N}` fetch. Raw bytes are the
/// [`RawRefSeq`] capability; buffer eviction is the windowed impl's inherent
/// `evict_before` — neither is a method here.
///
/// Coordinates are bare `u64` (1-based `start_1based`, `length`) rather than newtypes, to
/// match the production `MultiChromRefFetcher::fetch` signature so a winning impl ports
/// back with no signature change (ref_seq.md, Decision 3).
pub trait RefSeq {
    /// Write the canonical `{A,C,G,T,N}` bases for the `length` bases starting at 1-based
    /// `start_1based` on `contig` into `dst`. On success, `dst`'s previous contents are
    /// **replaced**; on error `dst` is left unchanged. Alloc-free when `dst` is reused.
    fn fetch_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError>;

    /// Owned convenience over [`Self::fetch_into`].
    fn fetch(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
    ) -> Result<Vec<u8>, RefSeqError> {
        // Start empty rather than `with_capacity(length)`: `length` is caller-controlled
        // and only validated inside `fetch_into`, so reserving up front would let an
        // invalid (huge) request allocate before it is rejected.
        let mut out = Vec::new();
        self.fetch_into(contig, start_1based, length, &mut out)?;
        Ok(out)
    }
}

/// Raw, un-canonicalised reference bytes — the left-alignment / mismatch-fraction path,
/// which needs the verbatim reference bytes (no `{A,C,G,T,N}` folding). A capability of
/// impls that can serve raw bytes; an impl whose only representation is already
/// canonicalised simply does not implement it, so "no raw here" is a compile-time fact
/// rather than a runtime error.
///
/// Like [`RefSeq::fetch_into`], raw bytes are **written into the caller's `dst`** (`&self`,
/// alloc-free when reused) rather than returned as a borrowed slice: a file-backed impl
/// loads its contig lazily, and a borrowed return would force `&mut self` (a resident
/// cache) — the copy into `dst` keeps the whole trait `&self`. (ref_seq.md, Decision 4.)
pub trait RawRefSeq: RefSeq {
    /// Write the raw bytes for the window into `dst` (replacing its contents on success;
    /// `dst` unchanged on error).
    fn fetch_raw_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError>;
}

/// A synthetic, fully in-memory reference: each contig's bytes held directly, indexed by
/// [`ContigId`] (`contigs[i]` is `ContigId(i)`). For tests — it needs no FASTA on disk,
/// so the mismatch filter and the pileup can be exercised against a hand-built reference.
/// Stored bytes are kept verbatim (raw); [`RefSeq::fetch_into`] folds them to canonical
/// on the fly, exactly as the file-backed impls do.
pub struct InMemoryRefSeq {
    contigs: Vec<Vec<u8>>,
}

impl InMemoryRefSeq {
    /// Build from one byte vector per contig, in [`ContigId`] order.
    pub fn from_contigs(contigs: Vec<Vec<u8>>) -> Self {
        Self { contigs }
    }

    /// Resolve a `(contig, start_1based, length)` request to the raw stored slice, or the
    /// matching [`RefSeqError`]. Both fetch paths go through it.
    fn resolve_range(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
    ) -> Result<&[u8], RefSeqError> {
        if start_1based == 0 {
            return Err(RefSeqError::InvalidStart);
        }
        let bytes = self
            .contigs
            .get(contig.get() as usize)
            .ok_or(RefSeqError::UnknownContig(contig))?;
        let range = validate_window(contig, bytes.len(), start_1based, length)?;
        Ok(&bytes[range])
    }
}

impl RefSeq for InMemoryRefSeq {
    fn fetch_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError> {
        let raw = self.resolve_range(contig, start_1based, length)?;
        dst.clear();
        dst.extend(raw.iter().copied().map(canonicalise));
        Ok(())
    }
}

impl RawRefSeq for InMemoryRefSeq {
    fn fetch_raw_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError> {
        let raw = self.resolve_range(contig, start_1based, length)?;
        dst.clear();
        dst.extend_from_slice(raw);
        Ok(())
    }
}

/// A FASTA-backed reference over a shared noodles [`Repository`] (the whole-contig cache
/// the production CRAM/pileup path already uses). Stateless per call: each fetch pulls the
/// contig's resident `Arc<Sequence>` from the repository, slices the window, and copies it
/// into the caller's `dst` — so the impl is `&self` and holds no per-fetch state (a clean
/// [`RefSeq`] + [`RawRefSeq`] pair). Mirrors the production `RepositoryRefFetcher`
/// (canonical) + `RawContigRefCache` (raw). One contig is resident at a time;
/// [`Self::clear`] drops it at a contig transition (the `--regions` memory discipline).
pub struct ResidentRefSeq {
    repository: Repository,
    contigs: ContigList,
}

impl ResidentRefSeq {
    /// Build over a shared repository and its contig table. The `repository` is an
    /// `Arc`-backed handle (cheap clone); pass the same instance the reader uses so the
    /// resident contigs are shared. `contigs` maps [`ContigId`] → contig name.
    pub fn new(repository: Repository, contigs: ContigList) -> Self {
        Self {
            repository,
            contigs,
        }
    }

    /// Drop the resident contig(s) from the underlying repository cache — call at a contig
    /// transition to keep at most one contig resident (the `--regions` memory rule). A
    /// subsequent fetch transparently reloads.
    pub fn clear(&self) {
        self.repository.clear();
    }

    /// Pull the contig's resident bytes, validate the window, and hand the raw slice to
    /// `consume` while the `Arc<Sequence>` is alive (the slice cannot escape the `Arc`, so
    /// both fetch paths consume it in place — one canonicalising, one copying verbatim).
    fn with_window<R>(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        consume: impl FnOnce(&[u8]) -> R,
    ) -> Result<R, RefSeqError> {
        if start_1based == 0 {
            return Err(RefSeqError::InvalidStart);
        }
        let entry = self
            .contigs
            .entries
            .get(contig.get() as usize)
            .ok_or(RefSeqError::UnknownContig(contig))?;
        let seq = match self.repository.get(entry.name.as_bytes()) {
            Some(Ok(seq)) => seq,
            Some(Err(source)) => return Err(RefSeqError::Io { contig, source }),
            None => {
                return Err(RefSeqError::Io {
                    contig,
                    source: io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("contig {} not in FASTA repository", entry.name),
                    ),
                });
            }
        };
        // Arc<Sequence> -> &Sequence -> &[u8] (the resident bytes, verbatim).
        let raw: &[u8] = AsRef::<[u8]>::as_ref(seq.as_ref());
        let range = validate_window(contig, raw.len(), start_1based, length)?;
        Ok(consume(&raw[range]))
    }
}

impl RefSeq for ResidentRefSeq {
    fn fetch_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError> {
        self.with_window(contig, start_1based, length, |raw| {
            dst.clear();
            dst.extend(raw.iter().copied().map(canonicalise));
        })
    }
}

impl RawRefSeq for ResidentRefSeq {
    fn fetch_raw_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError> {
        self.with_window(contig, start_1based, length, |raw| {
            dst.clear();
            dst.extend_from_slice(raw);
        })
    }
}

/// Map the single-contig fetcher's error into the ng error, tagging the contig.
fn map_chrom_error(contig: ContigId, err: ChromRefFetchError) -> RefSeqError {
    match err {
        ChromRefFetchError::InvalidStart => RefSeqError::InvalidStart,
        ChromRefFetchError::OutOfBounds {
            chrom_length,
            start,
            end,
            ..
        } => RefSeqError::OutOfBounds {
            contig,
            // Widening from production's `u32` error: lossless, and the only
            // direction that is.
            contig_length: u64::from(chrom_length),
            start: u64::from(start),
            end: u64::from(end),
        },
        ChromRefFetchError::Io { source, .. } => RefSeqError::Io { contig, source },
        // `OutOfPattern` (streaming-only) and any future variant: the manual-evict fetcher
        // never produces them, so fold defensively into Io.
        other => RefSeqError::Io {
            contig,
            source: io::Error::other(other.to_string()),
        },
    }
}

/// A FASTA-backed reference that keeps only a **sub-range window** of the current contig
/// resident, extending it on demand in either direction and shrinking it on the caller's
/// command via [`Self::evict_before`] — the memory-bounded, any-access-within-the-window
/// alternative to [`ResidentRefSeq`]'s whole-contig residency. Canonical-only (its buffer
/// holds `{A,C,G,T,N}` bytes), so it does **not** implement [`RawRefSeq`]. It wraps the
/// production `ManualEvictChromRefFetcher` (one per resident contig, rebuilt on a contig
/// change), so its buffer / canonicalisation / eviction logic — and its bytes — are
/// identical to the production streaming path. A `RefCell` keeps `fetch_into(&self)` while
/// the inner fetcher mutates its buffer; `evict_before` is the inherent `&mut self`
/// capability. `Send` but not `Sync` (per-worker ownership, like the production fetchers).
pub struct WindowedRefSeq {
    fasta_path: PathBuf,
    contigs: ContigList,
    current: RefCell<Option<(ContigId, ManualEvictChromRefFetcher)>>,
}

impl WindowedRefSeq {
    /// Build over a FASTA path (its sibling `<path>.fai` is used) and the contig table. No
    /// contig is resident until the first fetch.
    pub fn new(fasta_path: PathBuf, contigs: ContigList) -> Self {
        Self {
            fasta_path,
            contigs,
            current: RefCell::new(None),
        }
    }

    /// The reference's contig table — names and lengths, in `@SQ` / `.fai` order.
    ///
    /// **Provenance, not convenience** (typed_regions.md §2.6, §3). The walk must
    /// tell `admit` how long a contig *really* is, and that number cannot be
    /// derived from the window in hand — a slice mistaken for a chromosome is the
    /// silent bug §2.6 exists to kill. `admit` guards what arithmetic can catch,
    /// but its check is inherently incomplete: a caller passing the window's own
    /// end as `contig_len` is arithmetically legal, and only *reading the length
    /// from here* rules it out. That is why this accessor is part of the walk's
    /// substrate rather than a getter.
    pub fn contigs(&self) -> &ContigList {
        &self.contigs
    }

    /// Release buffered bytes before `pos` on the currently-resident contig — the
    /// caller-driven memory bound (production's `ManualEvictChromRefFetcher::evict_before`:
    /// drain `[.., pos)`, keep capacity). No-op when no contig is resident. Correctness is
    /// preserved: a later fetch of an evicted position simply re-reads it.
    pub fn evict_before(&mut self, pos: u32) {
        if let Some((_, fetcher)) = self.current.get_mut() {
            fetcher.evict_before(pos);
        }
    }
}

impl RefSeq for WindowedRefSeq {
    fn fetch_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError> {
        let entry = self
            .contigs
            .entries
            .get(contig.get() as usize)
            .ok_or(RefSeqError::UnknownContig(contig))?;
        let mut current = self.current.borrow_mut();
        let needs_rebuild = !matches!(&*current, Some((resident, _)) if *resident == contig);
        if needs_rebuild {
            let fetcher = ManualEvictChromRefFetcher::for_contig(&self.fasta_path, &entry.name)
                .map_err(|e| map_chrom_error(contig, e))?;
            *current = Some((contig, fetcher));
        }
        let (_, fetcher) = current.as_mut().expect("current set above");
        // **ng speaks `u64`; production's fetcher speaks `u32`** (spec §4 — `src/fasta/`
        // is frozen, so the narrowing lives here, at the seam, and nowhere else).
        //
        // It is an **error, not a clamp**. That distinction is the whole of B2: the
        // code this replaced wrote `unwrap_or(u32::MAX)` and a >4 Gb contig silently
        // became a wrong number. ng can now *represent* such a coordinate; this
        // fetcher simply cannot *serve* it, and saying so is the honest answer.
        let (start_u32, len_u32) = narrow_for_fetcher(contig, entry, start_1based, length)?;
        // The fetcher returns canonical {A,C,G,T,N} bytes (it reuses `canonicalise` on
        // refill), so they can be copied out verbatim — byte-identical to ResidentRefSeq.
        let bytes = fetcher
            .fetch(start_u32, len_u32)
            .map_err(|e| map_chrom_error(contig, e))?;
        dst.clear();
        dst.extend_from_slice(bytes);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::fetcher::RepositoryRefFetcher;
    use crate::fasta::{
        ChromRefFetcher, ContigEntry, MultiChromRefFetcher, StreamingChromRefFetcher,
    };
    use crate::pileup::per_sample::read_processor::RawContigRefCache;
    use std::io::Write;

    // ----- InMemoryRefSeq -------------------------------------------------

    /// contig 0 mixes case, an ambiguity code, a gap, and a stray byte to exercise
    /// canonicalisation; contig 1 is a clean `ACGT` for coordinate tests.
    fn in_memory() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![b"acgtNR-x".to_vec(), b"ACGT".to_vec()])
    }

    /// Convenience: raw bytes as an owned Vec (the trait writes into a caller buffer).
    fn raw(r: &impl RawRefSeq, contig: ContigId, start: u64, len: u64) -> Vec<u8> {
        let mut dst = Vec::new();
        r.fetch_raw_into(contig, start, len, &mut dst).unwrap();
        dst
    }

    #[test]
    fn fetch_canonical_uppercases_acgt_and_folds_the_rest_to_n() {
        let r = in_memory();
        assert_eq!(r.fetch(ContigId(0), 1, 8).unwrap(), b"ACGTNNNN");
    }

    #[test]
    fn canonicalise_folds_rna_uracil_to_n() {
        let r = InMemoryRefSeq::from_contigs(vec![b"UuAa".to_vec()]);
        assert_eq!(r.fetch(ContigId(0), 1, 4).unwrap(), b"NNAA");
    }

    #[test]
    fn fetch_raw_returns_stored_bytes_verbatim() {
        let r = in_memory();
        assert_eq!(raw(&r, ContigId(0), 1, 8), b"acgtNR-x");
    }

    #[test]
    fn fetch_into_replaces_dst_contents() {
        let r = in_memory();
        let mut dst = b"stale-contents".to_vec();
        r.fetch_into(ContigId(1), 1, 4, &mut dst).unwrap();
        assert_eq!(dst, b"ACGT");
    }

    #[test]
    fn fetch_into_leaves_dst_unchanged_on_error() {
        let r = in_memory();
        let mut dst = b"KEEP".to_vec();
        assert!(r.fetch_into(ContigId(1), 1, 99, &mut dst).is_err());
        assert_eq!(dst, b"KEEP");
    }

    #[test]
    fn fetch_matches_fetch_into() {
        let r = in_memory();
        let owned = r.fetch(ContigId(0), 2, 3).unwrap();
        let mut dst = Vec::new();
        r.fetch_into(ContigId(0), 2, 3, &mut dst).unwrap();
        assert_eq!(owned, dst);
    }

    #[test]
    fn sub_range_fetch_is_1_based() {
        let r = in_memory();
        assert_eq!(r.fetch(ContigId(1), 2, 3).unwrap(), b"CGT");
        assert_eq!(raw(&r, ContigId(1), 2, 3), b"CGT");
    }

    #[test]
    fn fetch_last_base_of_contig_returns_tail() {
        let r = in_memory();
        assert_eq!(r.fetch(ContigId(1), 4, 1).unwrap(), b"T");
        assert_eq!(raw(&r, ContigId(1), 4, 1), b"T");
    }

    #[test]
    fn zero_length_fetch_is_empty_not_error() {
        let r = in_memory();
        assert_eq!(r.fetch(ContigId(1), 1, 0).unwrap(), b"");
        assert_eq!(r.fetch(ContigId(1), 5, 0).unwrap(), b"");
    }

    #[test]
    fn fetch_on_empty_contig_returns_empty_for_zero_length() {
        let r = InMemoryRefSeq::from_contigs(vec![Vec::new()]);
        assert_eq!(r.fetch(ContigId(0), 1, 0).unwrap(), b"");
        assert!(matches!(
            r.fetch(ContigId(0), 1, 1),
            Err(RefSeqError::OutOfBounds {
                contig_length: 0,
                ..
            })
        ));
    }

    #[test]
    fn start_zero_is_invalid_start() {
        let r = in_memory();
        assert!(matches!(
            r.fetch(ContigId(1), 0, 1),
            Err(RefSeqError::InvalidStart)
        ));
    }

    #[test]
    fn unknown_contig_id_is_reported() {
        let r = in_memory();
        assert!(matches!(
            r.fetch(ContigId(9), 1, 1),
            Err(RefSeqError::UnknownContig(ContigId(9)))
        ));
    }

    #[test]
    fn window_past_contig_end_is_out_of_bounds() {
        let r = in_memory();
        match r.fetch(ContigId(1), 1, 6) {
            Err(RefSeqError::OutOfBounds {
                contig,
                contig_length,
                start,
                end,
            }) => {
                assert_eq!(contig, ContigId(1));
                assert_eq!(contig_length, 4);
                assert_eq!(start, 1);
                assert_eq!(end, 7);
            }
            other => panic!("expected OutOfBounds, got {other:?}"),
        }
    }

    /// `contigs()` hands back the reference's own table — the **provenance** the
    /// walk needs to tell `admit` a contig's true length (typed_regions.md §2.6).
    /// `admit`'s own guard cannot catch a caller who passes the window's end as
    /// `contig_len`; reading the length from here is what rules that out.
    /// **B2's payoff: a coordinate past `u32` is reported, not clamped.**
    ///
    /// The line this pins used to read
    /// `contig_length: u32::try_from(contig_len).unwrap_or(u32::MAX)`, so on a
    /// contig above 4 Gb the error *silently lied about the very number it
    /// existed to report*, in built ng code. Nothing could catch it, because
    /// nothing could represent the right answer.
    ///
    /// Now `RefSeqError` is `u64` end to end: a request beyond `u32::MAX` comes
    /// back with its own coordinates intact, so a >4 Gb assembly (they exist:
    /// lungfish, some conifers) gets a true error rather than a wrong number.
    /// Unreachable through `InMemoryRefSeq`'s tiny fixture in the *length*, but
    /// exactly reachable in the *request*, which is the half that used to clamp.
    #[test]
    fn a_request_beyond_u32_is_reported_verbatim_not_clamped() {
        let r = in_memory();
        let past_u32 = u64::from(u32::MAX) + 1_000;
        match r.fetch(ContigId(1), past_u32, 10) {
            Err(RefSeqError::OutOfBounds {
                contig_length,
                start,
                end,
                ..
            }) => {
                // The request survives the round trip at full width — no `u32::MAX`
                // anywhere.
                assert_eq!(start, past_u32, "start reported verbatim");
                assert_eq!(end, past_u32 + 10, "end reported verbatim");
                assert!(start > u64::from(u32::MAX), "the point of the test");
                // And the contig's real length, not a clamp of it.
                assert_eq!(contig_length, 4);
            }
            other => panic!("expected OutOfBounds, got {other:?}"),
        }
    }

    #[test]
    fn one_base_past_the_end_is_out_of_bounds_but_exact_fit_succeeds() {
        let r = in_memory();
        assert!(matches!(
            r.fetch(ContigId(1), 1, 5),
            Err(RefSeqError::OutOfBounds { .. })
        ));
        assert_eq!(r.fetch(ContigId(1), 1, 4).unwrap(), b"ACGT");
    }

    #[test]
    fn raw_and_canonical_differ_only_by_folding() {
        let r = in_memory();
        let raw_bytes = raw(&r, ContigId(0), 1, 8);
        let canon = r.fetch(ContigId(0), 1, 8).unwrap();
        let expect: Vec<u8> = raw_bytes.iter().copied().map(canonicalise).collect();
        assert_eq!(canon, expect);
    }

    // ----- ResidentRefSeq -------------------------------------------------

    /// Write a multi-contig, single-line-per-contig FASTA + `.fai` to a tempdir and return
    /// the guard, the path, and the matching `ContigList` (mirrors the production
    /// `fetcher.rs` test fixtures). Bytes must be FASTA-legal (letters/IUPAC/`N`); the
    /// contig strings here exercise soft-masking + ambiguity-code canonicalisation.
    fn build_fasta(contigs: &[(&str, &[u8])]) -> (tempfile::TempDir, PathBuf, ContigList) {
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");
        let mut fa = std::fs::File::create(&fasta_path).expect("create fasta");
        let mut fai_txt = String::new();
        let mut offset = 0usize;
        let mut entries = Vec::new();
        for (name, bases) in contigs {
            let header = format!(">{name}\n");
            fa.write_all(header.as_bytes()).expect("hdr");
            fa.write_all(bases).expect("seq");
            fa.write_all(b"\n").expect("nl");
            let seq_offset = offset + header.len();
            // name\tlength\toffset\tline_bases\tline_width
            fai_txt.push_str(&format!(
                "{name}\t{len}\t{seq_offset}\t{len}\t{width}\n",
                len = bases.len(),
                width = bases.len() + 1,
            ));
            offset = seq_offset + bases.len() + 1;
            entries.push(ContigEntry {
                name: (*name).to_string(),
                length: bases.len() as u64,
                md5: None,
            });
        }
        std::fs::write(&fai_path, fai_txt).expect("write fai");
        (dir, fasta_path, ContigList { entries })
    }

    fn repository_for(fasta_path: &std::path::Path) -> Repository {
        let indexed = noodles_fasta::io::indexed_reader::Builder::default()
            .build_from_path(fasta_path)
            .expect("indexed fasta");
        let adapter = noodles_fasta::repository::adapters::IndexedReader::new(indexed);
        Repository::new(adapter)
    }

    /// Fixture: chr0 exercises soft-masking + IUPAC folding; chr1 is clean.
    const FASTA_CONTIGS: &[(&str, &[u8])] = &[("chr0", b"acgtNRYK"), ("chr1", b"ACGTACGT")];

    fn resident() -> (tempfile::TempDir, ResidentRefSeq, ContigList) {
        let (dir, path, contigs) = build_fasta(FASTA_CONTIGS);
        let resident = ResidentRefSeq::new(repository_for(&path), contigs.clone());
        (dir, resident, contigs)
    }

    #[test]
    fn resident_canonical_matches_production_repository_fetcher() {
        let (_dir, path, contigs) = build_fasta(FASTA_CONTIGS);
        let resident = ResidentRefSeq::new(repository_for(&path), contigs.clone());
        let production = RepositoryRefFetcher::new(repository_for(&path), contigs);

        for (chrom_id, (_name, bases)) in FASTA_CONTIGS.iter().enumerate() {
            let len = bases.len() as u64;
            for start in 1..=len {
                for length in 0..=(len - start + 1) {
                    let ours = resident
                        .fetch(ContigId(chrom_id as u32), start, length)
                        .unwrap();
                    // ng is `u64`, production's fetcher is `u32` (spec §4; `src/fasta/`
                    // is frozen). Narrowing here is what makes the two comparable at
                    // all, and the fixture is a few dozen bases.
                    let prod = MultiChromRefFetcher::fetch(
                        &production,
                        chrom_id as u32,
                        start as u32,
                        length as u32,
                    )
                    .expect("production fetch");
                    assert_eq!(
                        ours, prod,
                        "canonical mismatch chrom {chrom_id} start {start} len {length}"
                    );
                }
            }
        }
    }

    #[test]
    fn resident_raw_matches_production_raw_cache() {
        let (_dir, path, contigs) = build_fasta(FASTA_CONTIGS);
        let resident = ResidentRefSeq::new(repository_for(&path), contigs.clone());
        let mut production = RawContigRefCache::new(repository_for(&path), contigs);

        for (chrom_id, (_name, bases)) in FASTA_CONTIGS.iter().enumerate() {
            let len = bases.len() as u64;
            for start in 1..=len {
                for length in 1..=(len - start + 1) {
                    let ours = raw(&resident, ContigId(chrom_id as u32), start, length);
                    let prod = production
                        .fetch_raw_slice(chrom_id, start, length as u32)
                        .expect("production raw slice");
                    assert_eq!(
                        ours.as_slice(),
                        prod,
                        "raw mismatch chrom {chrom_id} start {start} len {length}"
                    );
                }
            }
        }
    }

    #[test]
    fn resident_agrees_with_in_memory_on_canonical_and_raw() {
        let (_dir, resident, _contigs) = resident();
        let in_mem =
            InMemoryRefSeq::from_contigs(FASTA_CONTIGS.iter().map(|(_, b)| b.to_vec()).collect());

        for (chrom_id, (_name, bases)) in FASTA_CONTIGS.iter().enumerate() {
            let id = ContigId(chrom_id as u32);
            let len = bases.len() as u64;
            for start in 1..=len {
                for length in 1..=(len - start + 1) {
                    assert_eq!(
                        resident.fetch(id, start, length).unwrap(),
                        in_mem.fetch(id, start, length).unwrap(),
                        "canonical chrom {chrom_id} start {start} len {length}"
                    );
                    assert_eq!(
                        raw(&resident, id, start, length),
                        raw(&in_mem, id, start, length),
                        "raw chrom {chrom_id} start {start} len {length}"
                    );
                }
            }
        }
    }

    #[test]
    fn resident_raw_is_verbatim_fasta_bytes() {
        let (_dir, resident, _contigs) = resident();
        assert_eq!(raw(&resident, ContigId(0), 1, 8), b"acgtNRYK");
        assert_eq!(resident.fetch(ContigId(0), 1, 8).unwrap(), b"ACGTNNNN");
    }

    #[test]
    fn resident_fetch_into_replaces_stale_dst_contents() {
        let (_dir, resident, _contigs) = resident();
        let mut dst = b"stale".to_vec();
        resident.fetch_into(ContigId(1), 1, 4, &mut dst).unwrap();
        assert_eq!(dst, b"ACGT");
    }

    #[test]
    fn resident_leaves_dst_unchanged_on_error() {
        let (_dir, resident, _contigs) = resident();
        let mut dst = b"KEEP".to_vec();
        assert!(resident.fetch_into(ContigId(1), 1, 99, &mut dst).is_err());
        assert_eq!(dst, b"KEEP");
        let mut raw_dst = b"KEEP-RAW".to_vec();
        assert!(
            resident
                .fetch_raw_into(ContigId(1), 1, 99, &mut raw_dst)
                .is_err()
        );
        assert_eq!(raw_dst, b"KEEP-RAW");
    }

    #[test]
    fn resident_error_cases() {
        let (_dir, path, mut contigs) = build_fasta(FASTA_CONTIGS);
        // Add a contig present in the table but absent from the FASTA -> Io on fetch.
        contigs.entries.push(ContigEntry {
            name: "chrGhost".to_string(),
            length: 4,
            md5: None,
        });
        let resident = ResidentRefSeq::new(repository_for(&path), contigs);

        assert!(matches!(
            resident.fetch(ContigId(0), 0, 1),
            Err(RefSeqError::InvalidStart)
        ));
        assert!(matches!(
            resident.fetch(ContigId(9), 1, 1),
            Err(RefSeqError::UnknownContig(ContigId(9)))
        ));
        assert!(matches!(
            resident.fetch(ContigId(2), 1, 1),
            Err(RefSeqError::Io {
                contig: ContigId(2),
                ..
            })
        ));
        assert!(matches!(
            resident.fetch(ContigId(1), 1, 99),
            Err(RefSeqError::OutOfBounds { .. })
        ));
    }

    #[test]
    fn resident_clear_reloads_transparently() {
        let (_dir, resident, _contigs) = resident();
        let before = resident.fetch(ContigId(0), 1, 8).unwrap();
        resident.clear();
        let after = resident.fetch(ContigId(0), 1, 8).unwrap();
        assert_eq!(before, after);
    }

    // ----- WindowedRefSeq -------------------------------------------------

    fn windowed() -> (tempfile::TempDir, WindowedRefSeq) {
        let (dir, path, contigs) = build_fasta(FASTA_CONTIGS);
        (dir, WindowedRefSeq::new(path, contigs))
    }

    /// `contigs()` hands back the reference's own table — the **provenance** the
    /// walk needs to tell `admit` a contig's true length (typed_regions.md §2.6).
    /// `admit`'s guard cannot catch a caller who passes the *window's* end as
    /// `contig_len`, because that is arithmetically legal; reading the length from
    /// here is what rules it out. Hence an accessor on the walk's substrate, not a
    /// getter for its own sake.
    #[test]
    fn windowed_exposes_the_contig_table_as_provenance() {
        let (_dir, path, contigs) = build_fasta(FASTA_CONTIGS);
        let windowed = WindowedRefSeq::new(path, contigs.clone());

        assert_eq!(windowed.contigs(), &contigs, "the table, verbatim");
        // The number the walk actually reaches for — the CONTIG's length, and
        // never any window's. Already `u64` on production's table: only ng's
        // *fetch* surface ever narrowed it (the clamp B2 deleted).
        let entry = &windowed.contigs().entries[0];
        assert_eq!(entry.length, contigs.entries[0].length);
        assert!(!windowed.contigs().entries.is_empty());
    }

    #[test]
    fn windowed_canonical_matches_resident_across_all_windows() {
        let (_dir, path, contigs) = build_fasta(FASTA_CONTIGS);
        let resident = ResidentRefSeq::new(repository_for(&path), contigs.clone());
        let windowed = WindowedRefSeq::new(path, contigs);
        for (chrom_id, (_name, bases)) in FASTA_CONTIGS.iter().enumerate() {
            let id = ContigId(chrom_id as u32);
            let len = bases.len() as u64;
            for start in 1..=len {
                for length in 1..=(len - start + 1) {
                    assert_eq!(
                        windowed.fetch(id, start, length).unwrap(),
                        resident.fetch(id, start, length).unwrap(),
                        "chrom {chrom_id} start {start} len {length}"
                    );
                }
            }
        }
    }

    #[test]
    fn windowed_matches_production_streaming_on_a_forward_walk() {
        let (_dir, path, contigs) = build_fasta(FASTA_CONTIGS);
        let windowed = WindowedRefSeq::new(path.clone(), contigs);
        let streaming = StreamingChromRefFetcher::for_contig(&path, "chr1").unwrap();
        // chr1 = ACGTACGT (len 8); walk 2-base windows forward (monotonic → streaming-legal).
        for start in 1..=7u64 {
            let ours = windowed.fetch(ContigId(1), start, 2).unwrap();
            let prod = ChromRefFetcher::fetch(&streaming, start as u32, 2).unwrap();
            assert_eq!(ours, prod, "forward walk start {start}");
        }
    }

    #[test]
    fn windowed_supports_within_window_random_and_backward_access() {
        let (_dir, windowed) = windowed();
        // chr1 = ACGTACGT (len 8). Fetch forward, then jump backward — the manual-evict
        // buffer allows any access within the resident range (unlike strict streaming).
        assert_eq!(windowed.fetch(ContigId(1), 5, 4).unwrap(), b"ACGT");
        assert_eq!(windowed.fetch(ContigId(1), 1, 3).unwrap(), b"ACG"); // backward jump
        assert_eq!(windowed.fetch(ContigId(1), 3, 4).unwrap(), b"GTAC"); // overlaps both
    }

    #[test]
    fn windowed_eviction_preserves_correctness() {
        let (_dir, mut windowed) = windowed();
        let before = windowed.fetch(ContigId(1), 1, 8).unwrap();
        windowed.evict_before(5); // drop resident bytes before position 5
        // Re-fetching an evicted position simply re-reads it; the bytes are unchanged.
        assert_eq!(windowed.fetch(ContigId(1), 1, 8).unwrap(), before);
        // evict on a fresh contig / with nothing resident is a harmless no-op.
        assert_eq!(windowed.fetch(ContigId(0), 1, 8).unwrap(), b"ACGTNNNN");
    }

    #[test]
    fn windowed_rebuilds_on_contig_transition() {
        let (_dir, windowed) = windowed();
        assert_eq!(windowed.fetch(ContigId(0), 1, 8).unwrap(), b"ACGTNNNN");
        assert_eq!(windowed.fetch(ContigId(1), 1, 4).unwrap(), b"ACGT");
        assert_eq!(windowed.fetch(ContigId(0), 5, 4).unwrap(), b"NNNN"); // back to chr0
    }

    #[test]
    fn windowed_error_cases() {
        let (_dir, path, mut contigs) = build_fasta(FASTA_CONTIGS);
        contigs.entries.push(ContigEntry {
            name: "chrGhost".to_string(),
            length: 4,
            md5: None,
        });
        let windowed = WindowedRefSeq::new(path, contigs);
        assert!(matches!(
            windowed.fetch(ContigId(0), 0, 1),
            Err(RefSeqError::InvalidStart)
        ));
        assert!(matches!(
            windowed.fetch(ContigId(9), 1, 1),
            Err(RefSeqError::UnknownContig(ContigId(9)))
        ));
        assert!(matches!(
            windowed.fetch(ContigId(2), 1, 1),
            Err(RefSeqError::Io {
                contig: ContigId(2),
                ..
            })
        ));
        assert!(matches!(
            windowed.fetch(ContigId(1), 1, 99),
            Err(RefSeqError::OutOfBounds { .. })
        ));
    }
}
