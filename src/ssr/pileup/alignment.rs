//! Per-read delimitation — the Mark-2 core (arch `ssr_pileup_mark2.md` §3, model
//! Q2/Q1).
//!
//! For one read region (the locus-window slice from
//! [`super::footprint::extract_region`], clips included), align it against the
//! locus reference frame `left_flank + ref_tract + right_flank` with a **per-Q
//! Viterbi (max-path) pair-HMM + traceback**, and read the repeat region off the
//! two flank-junction columns. This is *delimitation only* — no likelihood, no
//! candidate scoring (that is Stage 2). The emission/transition model is the
//! Mark-1 pair-HMM's (per-base-quality Dindel emissions, affine gaps); the new
//! part is the max-path DP with a full backpointer matrix and the traceback.
//!
//! Determinism (the cross-sample identity guarantee — two reads of one molecule
//! must delimit to identical bytes): the traceback tie-break is the fixed priority
//! **Match > Deletion > Insertion**, and a junction indel is assigned to the
//! preceding (5′) block (so a left-junction insertion joins the flank, a
//! right-junction insertion joins the repeat).

use std::ops::Range;
use std::sync::LazyLock;

use crate::ssr::types::Locus;

/// DP-cell state indices into a `[_; 3]` cell: Match / Insertion / Deletion.
const M: usize = 0;
const I: usize = 1;
const D: usize = 2;
const NEG: f64 = f64::NEG_INFINITY;

/// First-quartile repeat-region quality cutoff (model Q1/P1): drop a read whose
/// repeat-region base-quality first quartile is below this. First value (Phred
/// 15) — calibrate on real data.
pub(crate) const MIN_REGION_Q1: u8 = 15;

/// Gap-open probability (Match → Insertion or → Deletion). Dindel's base value
/// for short homopolymer runs (Mark-1 pair-HMM, arch §5.4).
const GAP_OPEN_PROB: f64 = 2.9e-5;

/// Gap-extension probability: Dindel's fixed `e^-1 ≈ 0.368`.
static GAP_EXTEND_PROB: LazyLock<f64> = LazyLock::new(|| (-1.0f64).exp());

/// Emission log-probabilities per Phred quality, `(match_ln, mismatch_ln)`:
/// `match = ln(1 − ε)`, `mismatch = ln(ε / 3)` with `ε = 10^(−Q/10)` — the
/// Dindel base-quality model. Q0 (ε = 1) is floored so a Q0 match is merely
/// very unlikely, not an annihilating `-∞`.
static EMISSION_LN: LazyLock<[(f64, f64); 256]> = LazyLock::new(|| {
    let mut table = [(0.0f64, 0.0f64); 256];
    for (q, slot) in table.iter_mut().enumerate() {
        let eps = 10f64.powf(-(q as f64) / 10.0);
        let match_ln = (1.0 - eps).max(f64::MIN_POSITIVE).ln();
        let mismatch_ln = (eps / 3.0).ln();
        *slot = (match_ln, mismatch_ln);
    }
    table
});

/// Insertion emission: an inserted read base scored against a uniform base
/// distribution, `ln(1/4)`.
static INS_EMIT_LN: LazyLock<f64> = LazyLock::new(|| 0.25f64.ln());

/// The pair-HMM's log-space transition probabilities (Mark-1's, constant gap-open).
#[derive(Debug, Clone, Copy)]
pub(crate) struct HmmModel {
    ln_mm: f64,
    ln_gap_open: f64,
    ln_gap_close: f64,
    ln_gap_extend: f64,
}

impl Default for HmmModel {
    fn default() -> Self {
        let extend = *GAP_EXTEND_PROB;
        Self {
            ln_mm: (1.0 - 2.0 * GAP_OPEN_PROB).ln(),
            ln_gap_open: GAP_OPEN_PROB.ln(),
            ln_gap_close: (1.0 - extend).ln(),
            ln_gap_extend: extend.ln(),
        }
    }
}

impl HmmModel {
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

/// Per-worker Viterbi scratch: rolling score rows (max, not log-sum-exp) + a full
/// backpointer matrix for the traceback. Grow-and-keep so the hot path never
/// reallocates.
#[derive(Debug, Default)]
pub(crate) struct ViterbiScratch {
    prev: Vec<[f64; 3]>,
    cur: Vec<[f64; 3]>,
    /// Per cell `(i, j)` and per state, the predecessor state that won the max —
    /// flat `(m+1) × (n+1)` matrix, row-major with stride `n + 1`.
    back: Vec<[u8; 3]>,
    /// Reused buffer for the quality-gate quantile select.
    qual_buf: Vec<u8>,
}

impl ViterbiScratch {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    fn resize(&mut self, m: usize, n: usize) {
        let row = n + 1;
        if self.prev.len() < row {
            self.prev.resize(row, [NEG; 3]);
            self.cur.resize(row, [NEG; 3]);
        }
        let cells = (m + 1) * row;
        if self.back.len() < cells {
            self.back.resize(cells, [0u8; 3]);
        }
    }
}

/// Outcome of delimiting one read region against a locus.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Delimited {
    /// The repeat region is `region[range]` (read-coordinate, relative to the
    /// slice handed in). Both flank junctions landed inside the read.
    Region(Range<usize>),
    /// A flank ran off the read end (allele ≥ read length) — counted, not used.
    BorderOffEnd,
}

/// Max of the candidates, keeping the **first on ties** — so the caller encodes
/// the tie-break by passing candidates in priority order.
#[inline]
fn pick(cands: &[(f64, u8)]) -> (f64, u8) {
    let mut best = cands[0];
    for &c in &cands[1..] {
        if c.0 > best.0 {
            best = c;
        }
    }
    best
}

/// Delimit one read region (arch §3). Aligns `region_seq` (with `region_qual`)
/// against the locus reference frame `locus.ref_bytes()` and returns the
/// read-coordinate span of the repeat region, or [`Delimited::BorderOffEnd`] when
/// a flank is absent from the read.
pub(crate) fn delimit_read(
    region_seq: &[u8],
    region_qual: &[u8],
    locus: &Locus,
    model: &HmmModel,
    scratch: &mut ViterbiScratch,
) -> Delimited {
    debug_assert_eq!(
        region_seq.len(),
        region_qual.len(),
        "region seq/qual length mismatch"
    );
    let hap = locus.ref_bytes();
    let left_len = locus.left_flank().len();
    let right_len = locus.right_flank().len();
    let m = region_seq.len();
    let n = hap.len();
    if n == 0 {
        // No reference frame — not a real locus; nothing to delimit.
        return Delimited::BorderOffEnd;
    }

    scratch.resize(m, n);
    let nn = n + 1;
    let emission = &*EMISSION_LN;
    let ins_emit = *INS_EMIT_LN;

    let ViterbiScratch {
        prev, cur, back, ..
    } = scratch;

    // Row 0 (no read base consumed): begin in Match at the corner; the only
    // reachable cells are deletions of leading haplotype bases.
    prev[0] = [0.0, NEG, NEG];
    back[0] = [M as u8, M as u8, M as u8];
    for j in 1..=n {
        let (d, dp) = pick(&[
            (model.ln_gap_open + prev[j - 1][M], M as u8),
            (model.ln_gap_extend + prev[j - 1][D], D as u8),
        ]);
        prev[j] = [NEG, NEG, d];
        back[j] = [M as u8, M as u8, dp];
    }

    // Rows 1..=m: consume one read base per row.
    for i in 1..=m {
        let read_base = region_seq[i - 1];
        let (match_ln, mismatch_ln) = emission[region_qual[i - 1] as usize];
        let row = i * nn;

        // Column 0: only an insertion (read base before any haplotype base).
        let (ic, ip) = pick(&[
            (model.ln_gap_open + prev[0][M], M as u8),
            (model.ln_gap_extend + prev[0][I], I as u8),
        ]);
        cur[0] = [NEG, ins_emit + ic, NEG];
        back[row] = [M as u8, ip, M as u8];

        for j in 1..=n {
            let emit = if read_base == hap[j - 1] {
                match_ln
            } else {
                mismatch_ln
            };
            // Match: from the diagonal cell, any state. Priority M > D > I.
            let (mm, mp) = pick(&[
                (model.ln_mm + prev[j - 1][M], M as u8),
                (model.ln_gap_close + prev[j - 1][D], D as u8),
                (model.ln_gap_close + prev[j - 1][I], I as u8),
            ]);
            // Insertion: from the cell above (read consumed, hap not).
            let (ii, ipp) = pick(&[
                (model.ln_gap_open + prev[j][M], M as u8),
                (model.ln_gap_extend + prev[j][I], I as u8),
            ]);
            // Deletion: from the cell left (hap consumed, read not).
            let (dd, dpp) = pick(&[
                (model.ln_gap_open + cur[j - 1][M], M as u8),
                (model.ln_gap_extend + cur[j - 1][D], D as u8),
            ]);
            cur[j] = [emit + mm, ins_emit + ii, dd];
            back[row + j] = [mp, ipp, dpp];
        }
        std::mem::swap(prev, cur);
    }

    // Final cell (m, n): best terminal state, tie-break M > D > I.
    let last = prev[n];
    let (_, final_state) = pick(&[(last[M], M as u8), (last[D], D as u8), (last[I], I as u8)]);

    // Traceback. `r(k)` = read bases consumed before haplotype base `k` is
    // consumed; insertions are counted as soon as they occur, so they fall to
    // the block on their 5′ side (the preceding-block rule). The repeat region is
    // `read[r(left_len) .. r(n - right_len)]`.
    let left_junction = left_len; // first tract hap base
    let right_junction = n - right_len; // first right-flank hap base
    let mut tract_start = 0usize; // default: no left flank consumed
    let mut tract_end = m; // default: no right flank (runs to read end)

    let mut i = m;
    let mut j = n;
    let mut state = final_state as usize;
    while i != 0 || j != 0 {
        let pred = back[i * nn + j][state];
        match state {
            M => {
                let k = j - 1; // hap base this M consumes; r(k) = i - 1
                if k == left_junction {
                    tract_start = i - 1;
                }
                if right_len > 0 && k == right_junction {
                    tract_end = i - 1;
                }
                i -= 1;
                j -= 1;
            }
            D => {
                let k = j - 1; // hap base this D consumes; r(k) = i
                if k == left_junction {
                    tract_start = i;
                }
                if right_len > 0 && k == right_junction {
                    tract_end = i;
                }
                j -= 1;
            }
            _ => {
                // Insertion: consumes a read base, no haplotype base.
                i -= 1;
            }
        }
        state = pred as usize;
    }

    // A flank is "off the read end" when no read base sits on its side of the
    // tract (relies on Stage 0's clean unique flanks + the upstream reach gate).
    let left_off = left_len > 0 && tract_start == 0;
    let right_off = right_len > 0 && tract_end == m;
    if left_off || right_off || tract_start > tract_end {
        Delimited::BorderOffEnd
    } else {
        Delimited::Region(tract_start..tract_end)
    }
}

/// The Q1 quality gate: the **first quartile** of `quals` ≥ `threshold`. An empty
/// region (a zero-unit tract) passes vacuously.
pub(crate) fn passes_quality_gate(
    quals: &[u8],
    threshold: u8,
    scratch: &mut ViterbiScratch,
) -> bool {
    if quals.is_empty() {
        return true;
    }
    let buf = &mut scratch.qual_buf;
    buf.clear();
    buf.extend_from_slice(quals);
    let k = (buf.len() - 1) / 4;
    let (_, q1, _) = buf.select_nth_unstable(k);
    *q1 >= threshold
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::types::Motif;

    /// A perfect dinucleotide locus: GGG | CACACA | TTT, tract ref [13, 19),
    /// embedded window ref [10, 22). left_flank GGG (3), right_flank TTT (3).
    fn ca_locus() -> Locus {
        Locus::new(
            "chr1".into(),
            13,
            19,
            Motif::new(b"CA").unwrap(),
            1.0,
            (*b"GGGCACACATTT").into(),
            10,
        )
        .unwrap()
    }

    /// Delimit a region of all-Q40 bases.
    fn delimit(region: &[u8]) -> Delimited {
        let locus = ca_locus();
        let model = HmmModel::new();
        let mut s = ViterbiScratch::new();
        let quals = vec![40u8; region.len()];
        delimit_read(region, &quals, &locus, &model, &mut s)
    }

    #[test]
    fn clean_read_delimits_the_reference_tract() {
        // left GGG + tract CACACA (3 units) + right TTT.
        match delimit(b"GGGCACACATTT") {
            Delimited::Region(r) => {
                assert_eq!(&b"GGGCACACATTT"[r.clone()], b"CACACA");
                assert_eq!(r, 3..9);
            }
            other => panic!("expected Region, got {other:?}"),
        }
    }

    #[test]
    fn longer_allele_is_extracted_between_the_flanks() {
        // tract CACACACA (4 units): the extra CA is an insertion vs the ref tract,
        // but the flanks anchor, so the whole 8-base tract is extracted.
        match delimit(b"GGGCACACACATTT") {
            Delimited::Region(r) => assert_eq!(&b"GGGCACACACATTT"[r], b"CACACACA"),
            other => panic!("expected Region, got {other:?}"),
        }
    }

    #[test]
    fn impure_tract_is_extracted_verbatim() {
        // An interruption (the extra A) makes the tract non-tiling: CACAACA. The
        // flanks still anchor, so the off-ladder content is extracted exactly.
        match delimit(b"GGGCACAACATTT") {
            Delimited::Region(r) => assert_eq!(&b"GGGCACAACATTT"[r], b"CACAACA"),
            other => panic!("expected Region, got {other:?}"),
        }
    }

    #[test]
    fn shorter_allele_is_extracted() {
        // tract CACA (2 units): a deletion vs the ref tract; flanks anchor.
        match delimit(b"GGGCACATTT") {
            Delimited::Region(r) => assert_eq!(&b"GGGCACATTT"[r], b"CACA"),
            other => panic!("expected Region, got {other:?}"),
        }
    }

    #[test]
    fn missing_right_flank_is_border_off_end() {
        // The right flank ran off the read end (allele ≥ read length).
        assert_eq!(delimit(b"GGGCACACA"), Delimited::BorderOffEnd);
    }

    #[test]
    fn missing_left_flank_is_border_off_end() {
        // The read starts inside the tract — left flank absent.
        assert_eq!(delimit(b"CACACATTT"), Delimited::BorderOffEnd);
    }

    #[test]
    fn delimitation_is_deterministic() {
        // Same input twice → identical span (the traceback tie-break is fixed).
        let r1 = delimit(b"GGGCACAACATTT");
        let r2 = delimit(b"GGGCACAACATTT");
        assert_eq!(r1, r2);
    }

    #[test]
    fn quality_gate_passes_high_q_and_fails_low_q() {
        let mut s = ViterbiScratch::new();
        assert!(passes_quality_gate(
            &[40, 40, 40, 40],
            MIN_REGION_Q1,
            &mut s
        ));
        assert!(!passes_quality_gate(
            &[10, 10, 10, 10],
            MIN_REGION_Q1,
            &mut s
        ));
        // Empty region passes vacuously.
        assert!(passes_quality_gate(&[], MIN_REGION_Q1, &mut s));
    }

    #[test]
    fn quality_gate_keys_on_the_first_quartile_not_the_min() {
        let mut s = ViterbiScratch::new();
        // One low base among many high ones: first quartile stays high → pass.
        let mostly_high = [5u8, 40, 40, 40, 40, 40, 40, 40];
        assert!(passes_quality_gate(&mostly_high, MIN_REGION_Q1, &mut s));
        // A quarter low: first quartile drops below the threshold → fail.
        let quarter_low = [5u8, 5, 40, 40, 40, 40, 40, 40];
        assert!(!passes_quality_gate(&quarter_low, MIN_REGION_Q1, &mut s));
    }
}
