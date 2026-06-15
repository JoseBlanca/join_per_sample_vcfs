//! Read triage — classify reads near a locus and route usable ones to the
//! fast/slow path (arch §3).
//!
//! Everything downstream trusts this stage's sorting: anchor each read to the
//! catalog flanks, recover soft-clipped tract ends, and classify it as spanning
//! / flanking / in-repeat (arch §3.3); for a spanning read, decide its observed
//! repeat count and whether it takes the fast path (direct count) or the slow
//! path (pair-HMM), recording *why* for the slow path (arch §2).
//!
//! **Built so far: the content pre-probe ([`find_longest_stretch`]).** It is the
//! dependency-free estimator at the heart of soft-clip recovery (arch §3.2),
//! operating on raw read bytes. The read-consuming parts — anchoring, spanning
//! classification, the fast/slow gate, and the `TriageResult`/`SpanningRead`
//! types — are added next, over the shared `MappedRead`
//! ([`crate::bam::alignment_input`]) the SNP pileup already yields (arch §3.1).

use crate::ssr::types::Motif;

/// Result of the content pre-probe: a cheap `O(read len)` scan of a read's own
/// repeat content (arch §3.2). Used **only** when a flank lives in a soft-clip,
/// to estimate the (possibly large) repeat count without trusting the CIGAR.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ProbeHit {
    /// Longest run of **contiguous** (adjacent, non-overlapping) motif copies.
    /// This centres the pair-HMM window (`count ± W`) — the robust count at short
    /// periods, where incidental copies elsewhere inflate the total but not the
    /// contiguous run (arch §3.2).
    pub(crate) longest_run: u32,
    /// Total motif copies anywhere in the read (non-overlapping). A junk bound /
    /// slow-path trigger (against `MIN_MOTIF_RUN_BY_PERIOD`, §10), **not** the
    /// window centre.
    pub(crate) total_copies: u32,
}

/// Scan `seq` for its repeat content against `motif`, GangSTR-style — the
/// longest contiguous run of motif copies and the total copy count
/// (port of GangSTR `realignment.cpp` `find_longest_stretch`, arch §3.2).
///
/// Copies are counted **non-overlapping** (a match at `i` advances the scan by
/// `period`); a run is a maximal stretch of copies that are immediately
/// adjacent. The longest run centres the window; the total is the junk bound.
///
/// **One deliberate deviation from the C++ source:** GangSTR's loop bound
/// (`i < len − period`) misses a copy ending exactly at the read's end; we scan
/// every valid start (`i + period ≤ len`) so the final copy is counted. The
/// thresholds this feeds (`MIN_MOTIF_RUN_BY_PERIOD`) are uncalibrated placeholders
/// (arch §14), so the off-by-one is not a parity concern.
pub(crate) fn find_longest_stretch(seq: &[u8], motif: &Motif) -> ProbeHit {
    let unit = motif.as_bytes();
    let period = unit.len();
    if period == 0 || seq.len() < period {
        return ProbeHit {
            longest_run: 0,
            total_copies: 0,
        };
    }

    let mut longest_run = 0u32;
    let mut current_run = 0u32;
    let mut total_copies = 0u32;
    let mut i = 0usize;
    while i + period <= seq.len() {
        if &seq[i..i + period] == unit {
            total_copies += 1;
            current_run += 1;
            longest_run = longest_run.max(current_run);
            i += period; // non-overlapping; the next copy must abut to extend the run
        } else {
            current_run = 0;
            i += 1;
        }
    }

    ProbeHit {
        longest_run,
        total_copies,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn motif(bytes: &[u8]) -> Motif {
        Motif::new(bytes).unwrap()
    }

    #[test]
    fn pure_repeat_run_counts_all_copies() {
        let hit = find_longest_stretch(b"CACACACA", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 4,
                total_copies: 4
            }
        );
    }

    #[test]
    fn flanked_repeat_counts_only_the_tract() {
        // GGG | CACACA | TTT — three contiguous CA copies, nothing in the flanks.
        let hit = find_longest_stretch(b"GGGCACACATTT", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 3,
                total_copies: 3
            }
        );
    }

    #[test]
    fn interruption_splits_the_run_but_total_counts_both_halves() {
        // CACA | GG | CACA — two runs of 2, four copies total.
        let hit = find_longest_stretch(b"CACAGGCACA", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 2,
                total_copies: 4
            }
        );
    }

    #[test]
    fn incidental_copies_inflate_total_not_the_contiguous_run() {
        // The arch §3.2 point at short periods: isolated CA copies bump the total
        // but the longest contiguous run stays 1, so centring on the run is robust.
        let hit = find_longest_stretch(b"CAGGGGCAGGGGCA", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 1,
                total_copies: 3
            }
        );
    }

    #[test]
    fn homopolymer_counts_each_base() {
        let hit = find_longest_stretch(b"AAAAA", &motif(b"A"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 5,
                total_copies: 5
            }
        );
    }

    #[test]
    fn trinucleotide_run() {
        let hit = find_longest_stretch(b"CAGCAGCAGCAG", &motif(b"CAG"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 4,
                total_copies: 4
            }
        );
    }

    #[test]
    fn read_shorter_than_period_finds_nothing() {
        assert_eq!(
            find_longest_stretch(b"C", &motif(b"CA")),
            ProbeHit {
                longest_run: 0,
                total_copies: 0
            }
        );
        assert_eq!(
            find_longest_stretch(b"", &motif(b"CA")),
            ProbeHit {
                longest_run: 0,
                total_copies: 0
            }
        );
    }

    #[test]
    fn single_copy_at_the_read_end_is_counted() {
        // The deliberate fix vs GangSTR's `i < len - period` bound: a lone copy
        // exactly filling the read is found (GangSTR would miss it).
        let hit = find_longest_stretch(b"CA", &motif(b"CA"));
        assert_eq!(
            hit,
            ProbeHit {
                longest_run: 1,
                total_copies: 1
            }
        );
    }
}
