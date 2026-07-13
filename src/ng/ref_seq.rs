//! Reference-sequence access — the [`RefSeq`] trait and its implementations. The shared
//! source of reference-genome bases for the pileup, BAQ, DUST, left-alignment, and read
//! filtering. Design spec: `doc/devel/ng/spec/ref_seq.md`.
//!
//! [`RefSeq`] is the **universal** surface: canonical `{A,C,G,T,N}` fetch, which every
//! implementation provides. Capabilities only some impls have live on the impls that
//! have them — never as trait methods that silently do nothing on the rest: raw bytes
//! are the [`RawRefSeq`] sub-trait; buffer eviction is an inherent method on the windowed
//! impl (a later milestone). This file holds the traits, the error type, and the
//! synthetic [`InMemoryRefSeq`]; the file-backed impls (`ResidentRefSeq`,
//! `WindowedRefSeq`) land in later milestones.
//!
//! Canonicalisation reuses the production [`crate::fasta::fetcher::canonicalise`], so
//! canonical bytes are byte-identical to the production fetchers by construction.

use crate::fasta::fetcher::canonicalise;
use crate::ng::types::ContigId;

/// Errors from a reference fetch. `#[non_exhaustive]` so matchers must accept future
/// variants; the shape mirrors the production `fasta::fetcher::ChromRefFetchError`.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RefSeqError {
    /// The requested `[start, start + length)` window exceeds the contig's length.
    #[error("fetch [{start}, {end}) past {contig:?} (length {contig_length})")]
    OutOfBounds {
        contig: ContigId,
        contig_length: u32,
        start: u32,
        end: u32,
    },
    /// `start_1based` was 0 — the 1-based coordinate contract was violated.
    #[error("start_1based must be >= 1")]
    InvalidStart,
    /// No contig with this id exists in the reference.
    #[error("unknown {0:?}")]
    UnknownContig(ContigId),
    /// Underlying reference I/O failure (file-backed impls only).
    #[error("reference I/O failure on {contig:?}: {source}")]
    Io {
        contig: ContigId,
        #[source]
        source: std::io::Error,
    },
    // NOTE: a strict monotonic-forward (auto-slide) fetcher would add an `OutOfPattern`
    // variant here (cf. production `ChromRefFetchError`). ng's windowed impl is
    // caller-evictable and bounds memory by eviction, so it does not need one;
    // `#[non_exhaustive]` keeps adding it later non-breaking. (ref_seq.md, Decision 6.)
}

/// Access to reference-genome bases by (contig, 1-based range). The universal surface
/// every implementation provides: canonical `{A,C,G,T,N}` fetch. Raw bytes and buffer
/// eviction are capabilities on the impls that support them ([`RawRefSeq`]; the windowed
/// impl's inherent `evict_before`), not methods here — see the module docs.
///
/// Coordinates are bare `u32` (1-based `start_1based`, `length`) rather than newtypes, to
/// match the production `MultiChromRefFetcher::fetch` signature so a winning impl ports
/// back with no signature change (ref_seq.md, Decision 3).
pub trait RefSeq {
    /// Write the canonical `{A,C,G,T,N}` bases for the `length` bases starting at 1-based
    /// `start_1based` on `contig` into `dst`. On success, `dst`'s previous contents are
    /// **replaced**; on error `dst` is left unchanged. The allocation-free hot path.
    fn fetch_into(
        &self,
        contig: ContigId,
        start_1based: u32,
        length: u32,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError>;

    /// Owned convenience over [`Self::fetch_into`]. Canonical fetch returns owned bytes
    /// because canonicalisation produces fresh bytes; the raw read ([`RawRefSeq`]) can
    /// borrow instead.
    fn fetch(
        &self,
        contig: ContigId,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, RefSeqError> {
        // Start empty rather than `with_capacity(length)`: `length` is caller-controlled
        // and only validated inside `fetch_into`, so reserving up front would let an
        // invalid (huge) request allocate before it is rejected. `fetch_into` grows the
        // buffer only after the range is validated.
        let mut out = Vec::new();
        self.fetch_into(contig, start_1based, length, &mut out)?;
        Ok(out)
    }
}

/// Raw, un-canonicalised reference bytes (borrowed) — the left-alignment /
/// mismatch-fraction path. A capability of impls that keep bytes resident; an impl whose
/// buffer is already canonicalised simply does not implement it, so "no raw here" is a
/// compile-time fact rather than a runtime error.
pub trait RawRefSeq: RefSeq {
    fn fetch_raw(
        &self,
        contig: ContigId,
        start_1based: u32,
        length: u32,
    ) -> Result<&[u8], RefSeqError>;
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
    /// matching [`RefSeqError`]. The single place the coordinate/bounds contract is
    /// enforced; both `fetch_into` (which canonicalises the result) and `fetch_raw`
    /// (which returns it as-is) go through it.
    fn resolve_range(
        &self,
        contig: ContigId,
        start_1based: u32,
        length: u32,
    ) -> Result<&[u8], RefSeqError> {
        if start_1based == 0 {
            return Err(RefSeqError::InvalidStart);
        }
        let bytes = self
            .contigs
            .get(contig.get() as usize)
            .ok_or(RefSeqError::UnknownContig(contig))?;
        let start0 = (start_1based - 1) as usize;
        // `checked_add` guards a 32-bit `usize`, where `start0 + length` could wrap and
        // slip past the bounds check into a slice-index panic; on 64-bit it never fails.
        let end0 = match start0.checked_add(length as usize) {
            Some(end0) if end0 <= bytes.len() => end0,
            _ => {
                return Err(RefSeqError::OutOfBounds {
                    contig,
                    contig_length: bytes.len() as u32,
                    start: start_1based,
                    end: start_1based.saturating_add(length),
                });
            }
        };
        Ok(&bytes[start0..end0])
    }
}

impl RefSeq for InMemoryRefSeq {
    fn fetch_into(
        &self,
        contig: ContigId,
        start_1based: u32,
        length: u32,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError> {
        let raw = self.resolve_range(contig, start_1based, length)?;
        dst.clear();
        dst.extend(raw.iter().copied().map(canonicalise));
        Ok(())
    }
}

impl RawRefSeq for InMemoryRefSeq {
    fn fetch_raw(
        &self,
        contig: ContigId,
        start_1based: u32,
        length: u32,
    ) -> Result<&[u8], RefSeqError> {
        self.resolve_range(contig, start_1based, length)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// contig 0 mixes case, an ambiguity code, a gap, and a stray byte to exercise
    /// canonicalisation; contig 1 is a clean `ACGT` for coordinate tests.
    fn reference() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![b"acgtNR-x".to_vec(), b"ACGT".to_vec()])
    }

    #[test]
    fn fetch_canonical_uppercases_acgt_and_folds_the_rest_to_n() {
        let r = reference();
        // acgt -> ACGT; N stays N; R, -, x all fold to N.
        assert_eq!(r.fetch(ContigId(0), 1, 8).unwrap(), b"ACGTNNNN");
    }

    #[test]
    fn canonicalise_folds_rna_uracil_to_n() {
        // `U`/`u` are named in canonicalise's contract; pin that they fold to N.
        let r = InMemoryRefSeq::from_contigs(vec![b"UuAa".to_vec()]);
        assert_eq!(r.fetch(ContigId(0), 1, 4).unwrap(), b"NNAA");
    }

    #[test]
    fn fetch_raw_returns_stored_bytes_verbatim() {
        let r = reference();
        assert_eq!(r.fetch_raw(ContigId(0), 1, 8).unwrap(), b"acgtNR-x");
    }

    #[test]
    fn fetch_into_replaces_dst_contents() {
        let r = reference();
        let mut dst = b"stale-contents".to_vec();
        r.fetch_into(ContigId(1), 1, 4, &mut dst).unwrap();
        assert_eq!(dst, b"ACGT");
    }

    #[test]
    fn fetch_into_leaves_dst_unchanged_on_error() {
        let r = reference();
        let mut dst = b"KEEP".to_vec();
        // An overrunning request errors from resolve_range before dst.clear() runs.
        assert!(r.fetch_into(ContigId(1), 1, 99, &mut dst).is_err());
        assert_eq!(dst, b"KEEP");
    }

    #[test]
    fn fetch_matches_fetch_into() {
        let r = reference();
        let owned = r.fetch(ContigId(0), 2, 3).unwrap();
        let mut dst = Vec::new();
        r.fetch_into(ContigId(0), 2, 3, &mut dst).unwrap();
        assert_eq!(owned, dst);
    }

    #[test]
    fn sub_range_fetch_is_1_based() {
        let r = reference();
        // positions 2..=4 (1-based) of "ACGT" are "CGT".
        assert_eq!(r.fetch(ContigId(1), 2, 3).unwrap(), b"CGT");
        assert_eq!(r.fetch_raw(ContigId(1), 2, 3).unwrap(), b"CGT");
    }

    #[test]
    fn fetch_last_base_of_contig_returns_tail() {
        let r = reference();
        // start at the last base (1-based 4 of "ACGT"), length 1.
        assert_eq!(r.fetch(ContigId(1), 4, 1).unwrap(), b"T");
        assert_eq!(r.fetch_raw(ContigId(1), 4, 1).unwrap(), b"T");
    }

    #[test]
    fn zero_length_fetch_is_empty_not_error() {
        let r = reference();
        assert_eq!(r.fetch(ContigId(1), 1, 0).unwrap(), b"");
        // even at the one-past-the-end boundary (start = length + 1).
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
        let r = reference();
        assert!(matches!(
            r.fetch(ContigId(1), 0, 1),
            Err(RefSeqError::InvalidStart)
        ));
    }

    #[test]
    fn unknown_contig_id_is_reported() {
        let r = reference();
        assert!(matches!(
            r.fetch(ContigId(9), 1, 1),
            Err(RefSeqError::UnknownContig(ContigId(9)))
        ));
    }

    #[test]
    fn window_past_contig_end_is_out_of_bounds() {
        let r = reference();
        // "ACGT" has length 4; 6 bases from position 1 overruns.
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

    #[test]
    fn one_base_past_the_end_is_out_of_bounds_but_exact_fit_succeeds() {
        let r = reference();
        // "ACGT" is length 4: length 5 from pos 1 overruns by exactly one...
        assert!(matches!(
            r.fetch(ContigId(1), 1, 5),
            Err(RefSeqError::OutOfBounds { .. })
        ));
        // ...while the exact-fit neighbour (length 4) is fine.
        assert_eq!(r.fetch(ContigId(1), 1, 4).unwrap(), b"ACGT");
    }

    #[test]
    fn raw_and_canonical_differ_only_by_folding() {
        let r = reference();
        let raw = r.fetch_raw(ContigId(0), 1, 8).unwrap().to_vec();
        let canon = r.fetch(ContigId(0), 1, 8).unwrap();
        let expect: Vec<u8> = raw.iter().copied().map(canonicalise).collect();
        assert_eq!(canon, expect);
    }
}
