//! Per-contig post-processing: parsed TRF-mod records → catalog [`Locus`]es
//! (architecture [`ssr_catalog.md`](../../../doc/devel/architecture/ssr_catalog.md) §4).
//!
//! We **drop** compounds and bundles — we do not split (the GangSTR reference
//! pipeline does exactly this). Order (§4):
//!
//! 1. **Period 2..=6** (+ an early TRF `score` floor) — cheap scope/volume cut.
//!    Period-1 homopolymers are dropped (GangSTR/HipSTR convention; also stops
//!    poly-A/T runs from bundle-dropping adjacent SSRs in step 3).
//! 2. **Drop compound-motif loci** — a motif that is itself internally periodic
//!    (`ATAT` = `(AT)²`). Port of GangSTR `minimal_trim.py::is_compound`.
//! 3. **Drop bundles** — any locus within `bundle_threshold` bp of another is
//!    discarded with its whole cluster. Port of GangSTR `remove_bundles.py`.
//!    With `bundle_threshold ≥ flank_bp`, every survivor has clean unique
//!    flanks, so there is no inner-flank case downstream.
//! 4. **End-trim** partial motifs to clean whole-motif boundaries (GangSTR
//!    `minimal_trim`), then apply the per-period **copy-number floor**.
//! 5. **Recompute purity** from the trimmed tract vs a perfect motif tiling
//!    (spec §3.2 — *not* TRF's `fracMatch`), then drop below the floor.
//!    **Imperfect single-motif loci are kept** (our one divergence from
//!    GangSTR's perfect-only `remove_messy`).
//! 6. **Embed `ref_seq`** (trimmed tract + `flank_bp` each side, clamped at the
//!    contig ends), and **drop any locus whose flank clamped to zero** on either
//!    side — a tract abutting position 0 of the contig, or ending on its last
//!    base. The Stage-1 delimiter anchors on both flank junctions, so a survivor
//!    is guaranteed a non-empty flank each side (contig-boundary microsatellites
//!    are not genotypeable and are dropped here, not crashed downstream).
//!
//! Two adaptations from the GangSTR scripts: (a) GangSTR reads the repeat
//! sequence from a TRF column; we slice it from the resident contig
//! (`contig_seq[start..end]`). (b) GangSTR uses TRF's reported motif; our motif
//! is the first `period` bases of the reference tract (verbatim, phase-faithful,
//! types §5) — which, after `minimal_trim` aligns the tract to a motif boundary,
//! is exactly `trimmed_tract[0..period]`.

use super::CatalogParams;
use super::trf::TrfRecord;
use crate::ssr::types::{Locus, Motif};

/// Per-period minimum copy number a tract must reach to survive (GangSTR
/// `minimal_trim.py` `thresholds = {1:10, 2:5, 3:4, 4:3, 5:3, 6:3}`; default 3
/// for any other period). Period 1 is filtered out upstream ([`MIN_PERIOD`]),
/// so its floor is retained only for parity with the source table.
const fn copy_number_floor(period: usize) -> u32 {
    match period {
        1 => 10,
        2 => 5,
        3 => 4,
        4 => 3,
        5 => 3,
        6 => 3,
        _ => 3,
    }
}

/// The narrowest SSR period the catalog keeps. Period-1 **homopolymers are
/// excluded** (architecture §4): the standard GangSTR/HipSTR drop — error-prone
/// for STR genotyping and not the di/tri/tetra-nucleotide target — and dropping
/// them *before* bundling stops a long poly-A/T run from bundle-dropping an
/// adjacent real SSR.
const MIN_PERIOD: u16 = 2;

/// The widest SSR period the catalog keeps (architecture §4 step 1).
const MAX_PERIOD: u16 = 6;

/// Post-process one contig's TRF-mod records into start-sorted catalog loci.
/// `chrom` is the contig name; `contig_seq` is the full contig (any case — the
/// tract, motif, and `ref_seq` are upper-cased here for case-stable identity).
pub(crate) fn build_loci(
    recs: Vec<TrfRecord>,
    chrom: &str,
    contig_seq: &[u8],
    p: &CatalogParams,
) -> Vec<Locus> {
    // 1. scope + score gate, then 2. compound-motif drop. Both need the tract,
    //    so we slice the (upper-cased) prefix motif here and filter on it.
    let mut kept: Vec<TrfRecord> = recs
        .into_iter()
        .filter(|r| {
            r.period >= MIN_PERIOD
                && r.period <= MAX_PERIOD
                && r.score >= p.min_score
                && r.end > r.start
                && (r.end as usize) <= contig_seq.len()
        })
        .filter(|r| {
            let period = r.period as usize;
            let tract = &contig_seq[r.start as usize..r.end as usize];
            // Need at least one full motif to form the prefix.
            if tract.len() < period {
                return false;
            }
            let motif = upper(&tract[..period]);
            !is_compound(&motif)
        })
        .collect();

    // 3. drop bundles (on the raw, pre-trim coordinates). Records must be
    //    start-sorted for the streaming clustering.
    kept.sort_by_key(|r| (r.start, r.end));
    let kept = drop_bundles(kept, p.bundle_threshold);

    // 4-5. per record: end-trim + copy floor, recompute purity + floor, embed.
    let mut out = Vec::with_capacity(kept.len());
    for r in &kept {
        if let Some(locus) = finish_locus(r, chrom, contig_seq, p) {
            out.push(locus);
        }
    }
    out
}

/// Steps 4-5 for one record: end-trim, copy-number floor, motif, purity floor,
/// `ref_seq` embed. `None` if the record fails any gate.
fn finish_locus(r: &TrfRecord, chrom: &str, contig_seq: &[u8], p: &CatalogParams) -> Option<Locus> {
    let period = r.period as usize;
    let raw_tract = upper(&contig_seq[r.start as usize..r.end as usize]);
    let motif_bytes = raw_tract.get(..period)?.to_vec();

    // End-trim to clean whole-motif boundaries (GangSTR minimal_trim).
    let (st, en) = minimal_trim(&raw_tract, &motif_bytes)?;
    let new_start = r.start + st as u32;
    let new_end = r.start + en as u32;
    let trimmed = &raw_tract[st..en];

    // Copy-number floor — GangSTR computes copies from the ORIGINAL TRF span
    // (integer division), as an accept-gate after trimming.
    let ref_copy = (r.end - r.start) / r.period as u32;
    if ref_copy < copy_number_floor(period) {
        return None;
    }

    // After minimal_trim the tract starts on a motif boundary, so
    // `trimmed[..period] == motif_bytes`; use it as the phase-faithful motif.
    let motif = Motif::new(&motif_bytes).ok()?;

    // Recompute purity from the trimmed tract vs a perfect motif tiling.
    let purity = recompute_purity(trimmed, &motif_bytes);
    if purity < p.min_purity {
        return None;
    }

    // Embed ref_seq: trimmed tract + flank each side, clamped at contig ends.
    let ref_start = new_start.saturating_sub(p.flank_bp);
    let ref_end = (new_end + p.flank_bp).min(contig_seq.len() as u32);

    // Drop a locus whose flank clamped to nothing on either side — a tract
    // abutting position 0 of the contig (empty left flank) or ending on the
    // contig's last base (empty right flank). The empirical-candidate delimiter
    // (`ssr::pileup::alignment`) anchors the repeat region on *both* flank
    // junctions; a zero-length flank leaves nothing to anchor against, so the
    // tract is not genotypeable and must not reach Stage 1.
    if ref_start == new_start || ref_end == new_end {
        return None;
    }

    let ref_bytes = upper(&contig_seq[ref_start as usize..ref_end as usize]);

    Locus::new(
        chrom.to_string().into_boxed_str(),
        new_start,
        new_end,
        motif,
        purity,
        ref_bytes.into_boxed_slice(),
        ref_start,
    )
    .ok()
}

/// ASCII upper-case copy of `bytes`.
fn upper(bytes: &[u8]) -> Vec<u8> {
    bytes.iter().map(|b| b.to_ascii_uppercase()).collect()
}

/// Count non-overlapping greedy occurrences of `motif` in `repeat` (GangSTR
/// `count_motif`): scan left-to-right, advancing by `motif.len()` on a match
/// and by 1 otherwise.
fn count_motif(repeat: &[u8], motif: &[u8]) -> usize {
    let m = motif.len();
    let mut c = 0;
    let mut s = 0;
    while s < repeat.len() {
        if s + m <= repeat.len() && &repeat[s..s + m] == motif {
            c += 1;
            s += m;
        } else {
            s += 1;
        }
    }
    c
}

/// `true` if `motif` is itself internally periodic — a non-fundamental period
/// such as `ATAT = (AT)²` (GangSTR `is_compound`, threshold 0.8). A motif whose
/// shorter prefix `sub` tiles more than 80% of it is compound.
fn is_compound(motif: &[u8]) -> bool {
    let l = motif.len();
    const THRESHOLD: f64 = 0.8;
    for i in 1..=(l / 2) {
        let sub = &motif[..i];
        if (count_motif(motif, sub) * i) as f64 > l as f64 * THRESHOLD {
            return true;
        }
    }
    false
}

/// End-trim a tract to clean whole-motif boundaries (GangSTR `minimal_trim`):
/// find the smallest `start_offset` where two motif copies (`motif*2`) appear,
/// and the largest `end_offset` (exclusive) where they appear before it.
/// Returns `(start_offset, end_offset)` into `rep`, or `None` if the tract has
/// no clean motif boundary within `min(3·|motif|, |rep|/2)` of each end.
fn minimal_trim(rep: &[u8], motif: &[u8]) -> Option<(usize, usize)> {
    let m = motif.len();
    let ll = m * 2;
    let max_trim_len = (m * 3).min(rep.len() / 2);

    let two_motifs = |slice: &[u8]| -> bool {
        slice.len() == ll && &slice[..m] == motif && &slice[m..] == motif
    };

    // Smallest start offset whose next 2·m bytes are motif*2. GangSTR checks
    // the bound first and bails (not continues) — replicate.
    let mut start = None;
    for so in 0..=max_trim_len {
        if so + ll >= rep.len() {
            return None;
        }
        if two_motifs(&rep[so..so + ll]) {
            start = Some(so);
            break;
        }
    }
    let st = start?;

    // Largest end offset (exclusive) whose preceding 2·m bytes are motif*2.
    let mut end = None;
    for eo in (rep.len().saturating_sub(max_trim_len)..=rep.len()).rev() {
        if eo < ll {
            return None;
        }
        if two_motifs(&rep[eo - ll..eo]) {
            end = Some(eo);
            break;
        }
    }
    let en = end?;
    if st >= en {
        return None;
    }
    Some((st, en))
}

/// Fraction of `tract` matching a perfect tiling of `motif` from phase 0 (spec
/// §3.2). With `motif == tract[..period]`, the first repeat always matches, so
/// a perfect tract scores 1.0 and interruptions lower it proportionally.
fn recompute_purity(tract: &[u8], motif: &[u8]) -> f32 {
    if tract.is_empty() || motif.is_empty() {
        return 0.0;
    }
    let m = motif.len();
    let matches = tract
        .iter()
        .enumerate()
        .filter(|&(i, &b)| b == motif[i % m])
        .count();
    matches as f32 / tract.len() as f32
}

/// Two records are "close" (bundle candidates) if any of their start/end
/// coordinates are within `thresh` bp (GangSTR `is_close`, `check_motif=False`;
/// chrom equality is implicit — these are one contig's records).
fn is_close(a: &TrfRecord, b: &TrfRecord, thresh: u32) -> bool {
    a.start.abs_diff(b.start) < thresh
        || a.start.abs_diff(b.end) < thresh
        || b.start.abs_diff(a.end) < thresh
        || b.end.abs_diff(a.end) < thresh
}

/// Drop bundles: discard every maximal run of consecutive close records in
/// their entirety (GangSTR `remove_bundles.py`). `recs` must be start-sorted.
fn drop_bundles(recs: Vec<TrfRecord>, thresh: u32) -> Vec<TrfRecord> {
    let mut out = Vec::new();
    let n = recs.len();
    let mut i = 0;
    while i < n {
        if i + 1 < n && is_close(&recs[i], &recs[i + 1], thresh) {
            // Advance through the whole close cluster; drop all of it.
            let mut j = i;
            while j + 1 < n && is_close(&recs[j], &recs[j + 1], thresh) {
                j += 1;
            }
            i = j + 1;
        } else {
            out.push(recs[i].clone());
            i += 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> CatalogParams {
        CatalogParams {
            min_purity: 0.8,
            min_score: 0,
            flank_bp: 5,
            bundle_threshold: 50,
        }
    }

    // ---- is_compound -------------------------------------------------------

    #[test]
    fn is_compound_flags_atat_but_not_at_or_atc() {
        assert!(is_compound(b"ATAT"), "ATAT = (AT)^2 is compound");
        assert!(is_compound(b"ATATAT"), "ATATAT = (AT)^3 is compound");
        assert!(!is_compound(b"AT"), "fundamental AT is not compound");
        assert!(!is_compound(b"ATC"), "ATC is not internally periodic");
        assert!(!is_compound(b"ATCG"), "ATCG is not internally periodic");
    }

    // ---- recompute_purity --------------------------------------------------

    #[test]
    fn purity_is_one_for_perfect_tiling() {
        assert_eq!(recompute_purity(b"ATATATAT", b"AT"), 1.0);
        assert_eq!(recompute_purity(b"CAGCAGCAG", b"CAG"), 1.0);
    }

    #[test]
    fn purity_drops_with_one_interruption() {
        // One mismatch in 8 bases → 7/8.
        assert_eq!(recompute_purity(b"ATATCTAT", b"AT"), 7.0 / 8.0);
    }

    // ---- minimal_trim ------------------------------------------------------

    #[test]
    fn minimal_trim_strips_partial_motifs_at_both_ends() {
        // "T" + (AT)*5 + "A": a trailing half-motif at each end. motif = "TA"
        // is the prefix here; instead drive it with the clean interior.
        let rep = b"GAGATCGATCGATCGATCAG"; // not a clean repeat; expect a trim window
        // Use a clean case: 6 copies of "AT" with a 1bp junk prefix/suffix.
        let clean = b"GATATATATATATC";
        let (st, en) = minimal_trim(clean, b"AT").expect("clean AT run trims");
        assert_eq!(&clean[st..en], b"ATATATATATAT");
        // The junk-laden one should still trim to whole AT copies or fail
        // cleanly (never panic).
        let _ = minimal_trim(rep, b"AT");
    }

    #[test]
    fn minimal_trim_returns_none_without_two_clean_copies() {
        assert!(minimal_trim(b"ACGT", b"AT").is_none());
    }

    // ---- drop_bundles ------------------------------------------------------

    #[test]
    fn drop_bundles_discards_whole_clusters_keeps_isolated() {
        // A,B,C are a close cluster (within 50 bp); D is isolated (far away).
        let recs = vec![
            TrfRecord::for_test(100, 130, 2, 100, b"AT"),
            TrfRecord::for_test(150, 180, 3, 100, b"CAG"),
            TrfRecord::for_test(200, 230, 2, 100, b"AT"),
            TrfRecord::for_test(5000, 5030, 2, 100, b"AT"),
        ];
        let kept = drop_bundles(recs, 50);
        assert_eq!(kept.len(), 1, "only the isolated D survives");
        assert_eq!(kept[0].start, 5000);
    }

    #[test]
    fn drop_bundles_keeps_all_when_none_close() {
        let recs = vec![
            TrfRecord::for_test(100, 130, 2, 100, b"AT"),
            TrfRecord::for_test(1000, 1030, 2, 100, b"AT"),
            TrfRecord::for_test(2000, 2030, 2, 100, b"AT"),
        ];
        assert_eq!(drop_bundles(recs, 50).len(), 3);
    }

    // ---- build_loci end-to-end --------------------------------------------

    /// A synthetic contig with one clean (AT)*8 tract, isolated, period 2.
    /// Built by concatenation so the coordinates can't drift from the bytes.
    #[test]
    fn build_loci_emits_a_clean_locus_with_flanks() {
        let left = b"CGCGC"; // 5 bp left flank
        let tract = b"ATATATATATATATAT"; // 16 bp = 8 copies of AT
        let right = b"GCGCG"; // 5 bp right flank
        let contig = [left.as_ref(), tract.as_ref(), right.as_ref()].concat();
        let tract_start = left.len() as u32; // 5
        let tract_end = tract_start + tract.len() as u32; // 21
        let recs = vec![TrfRecord::for_test(tract_start, tract_end, 2, 100, b"AT")];
        let loci = build_loci(recs, "chr1", &contig, &params());
        assert_eq!(loci.len(), 1);
        let l = &loci[0];
        assert_eq!(l.chrom(), "chr1");
        assert_eq!(l.start(), tract_start);
        assert_eq!(l.end(), tract_end);
        assert_eq!(l.motif().as_bytes(), b"AT");
        assert_eq!(l.purity_fraction(), 1.0);
        assert_eq!(l.ref_tract(), tract);
        // 5 bp flank each side, clamped — here both flanks fit.
        assert_eq!(l.left_flank(), left);
        assert_eq!(l.right_flank(), right);
    }

    #[test]
    fn build_loci_drops_period_over_six_and_compound() {
        let contig = b"AAAACGATATATATATATATCCCC";
        // period 7 (> MAX_PERIOD) → dropped by the scope filter.
        let too_long = TrfRecord::for_test(0, 14, 7, 100, b"ACGATAT");
        // compound AT-as-period-4 (ATAT) → dropped by is_compound.
        let compound = TrfRecord::for_test(6, 20, 4, 100, b"ATAT");
        let loci = build_loci(vec![too_long, compound], "chr1", contig, &params());
        assert!(loci.is_empty(), "both records are dropped, got {loci:?}");
    }

    #[test]
    fn build_loci_drops_below_copy_number_floor() {
        // period 2 needs >= 5 copies; (AT)*3 = 6 bp is only 3 copies → dropped.
        let contig = b"CCCCCATATATCCCCC";
        let recs = vec![TrfRecord::for_test(5, 11, 2, 100, b"AT")];
        assert!(build_loci(recs, "chr1", contig, &params()).is_empty());
    }

    #[test]
    fn build_loci_clamps_flanks_at_contig_ends() {
        // Tract (AT)*8 starts at position 1; left flank (5 bp) is clamped to 1.
        // A 1 bp flank is still a (weak) anchor, so the locus is kept.
        let contig = b"GATATATATATATATATCGCGCGCGCG";
        let recs = vec![TrfRecord::for_test(1, 17, 2, 100, b"AT")];
        let loci = build_loci(recs, "chr1", contig, &params());
        assert_eq!(loci.len(), 1);
        let l = &loci[0];
        assert_eq!(l.ref_bytes_start(), 0, "left flank clamped to contig start");
        assert_eq!(l.left_flank(), b"G", "only 1 bp of left flank available");
    }

    #[test]
    fn build_loci_drops_locus_with_empty_left_flank() {
        // Tract (AT)*8 abuts position 0 — the left flank clamps to nothing, so
        // the delimiter would have no left junction to anchor: dropped.
        let contig = b"ATATATATATATATATCGCGC";
        let recs = vec![TrfRecord::for_test(0, 16, 2, 100, b"AT")];
        assert!(
            build_loci(recs, "chr1", contig, &params()).is_empty(),
            "a tract at contig position 0 has no left flank and is not genotypeable"
        );
    }

    #[test]
    fn build_loci_drops_locus_with_empty_right_flank() {
        // Tract (AT)*8 ends on the contig's last base — the right flank clamps
        // to nothing: dropped for the same reason as the empty-left case.
        let contig = b"CGCGCATATATATATATATAT";
        let recs = vec![TrfRecord::for_test(5, 21, 2, 100, b"AT")];
        assert!(
            build_loci(recs, "chr1", contig, &params()).is_empty(),
            "a tract ending on the contig's last base has no right flank"
        );
    }

    /// Period-1 homopolymers are dropped (MIN_PERIOD = 2), and — because the
    /// drop happens before bundling — a poly-A run adjacent to a real SSR no
    /// longer bundle-drops it.
    #[test]
    fn build_loci_drops_period_one_homopolymer_and_spares_the_neighbour_ssr() {
        // 20 bp poly-A, then (CAG)*10, then a right flank so the CAG locus is
        // not dropped for abutting the contig end (orthogonal to this test).
        let mut contig = vec![b'A'; 20];
        for _ in 0..10 {
            contig.extend_from_slice(b"CAG");
        }
        contig.extend_from_slice(b"TTTTT");
        let homopolymer = TrfRecord::for_test(0, 20, 1, 100, b"A");
        let cag = TrfRecord::for_test(20, 50, 3, 100, b"CAG");
        let loci = build_loci(vec![homopolymer, cag], "chr1", &contig, &params());
        assert_eq!(loci.len(), 1, "the homopolymer is gone, the CAG survives");
        assert_eq!(loci[0].motif().as_bytes(), b"CAG");
        assert_eq!(loci[0].period(), 3);
    }
}
