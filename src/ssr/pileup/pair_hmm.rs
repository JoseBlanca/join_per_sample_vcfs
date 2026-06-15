//! Slow path — the pair-HMM forward scorer (arch §5).
//!
//! For a read that fails the fast-path gate, this scores the read against a
//! candidate haplotype `H_L = left_flank + (motif × L) + right_flank` and
//! returns `Qᵣ(L)` — the log-likelihood the read's bases arose from a molecule
//! whose tract is `L` units long, under **sequencing/alignment error only** (no
//! stutter; that is Stage 2's `S_θ`). It is a pure scorer (`read × haplotype →
//! f64`): a forward (sum over alignments), never a traceback, so we keep only a
//! rolling two-row scratch and never reconstruct an alignment (arch §5.7).
//!
//! The model is HipSTR's flank pair-HMM minus its stutter marginalization
//! (arch §5.1): a 3-state (Match/Insertion/Deletion) forward in **log-space**
//! with log-sum-exp, emission from the Dindel/BAQ per-Q model (arch §5.3), and
//! affine gaps (arch §5.4). We reuse the BAQ engine's *patterns* — the per-Q
//! lookup, the grow-and-keep scratch — not its code (arch §5.6); the math here
//! is bespoke and simpler (no backward pass, no posterior decoding, no htslib
//! parity literals).
//!
//! **Deliberately deferred (arch §5.4/§5.5/§14, calibration not structure):**
//! - the **homopolymer-run-indexed gap-open** (Dindel's per-run-length table) —
//!   this uses a single constant gap-open for now; the table is a parameter
//!   refinement that swaps the constant for a position-indexed lookup;
//! - **banding** — this computes the full `O(read × hap)` DP. The band
//!   (`PAIR_HMM_BAND_BP`, arch §5.5) is an optimization to add once the slow
//!   path is shown to bind; an unbanded forward is the correct `W → ∞` limit.

use std::sync::LazyLock;

use crate::ssr::pileup::candidate_generation::CandidateAllele;
use crate::ssr::types::Allele;

/// DP-cell state indices into a `[f64; 3]` cell: Match / Insertion / Deletion.
const M: usize = 0;
const I: usize = 1;
const D: usize = 2;

/// Gap-open probability (Match → Insertion or → Deletion). Dindel's base value
/// for short homopolymer runs (`AlignmentModel.cpp`, run length ≤ 4, arch §5.4).
/// A single constant for now — the homopolymer-run-indexed table is deferred
/// (see the module doc).
const GAP_OPEN_PROB: f64 = 2.9e-5;

/// Gap-extension probability (Insertion → Insertion, Deletion → Deletion):
/// Dindel's fixed `e^-1 ≈ 0.368` (arch §5.4). `gap → Match` is `1 - e^-1`.
static GAP_EXTEND_PROB: LazyLock<f64> = LazyLock::new(|| (-1.0f64).exp());

/// Emission log-probabilities per Phred quality, `(match_ln, mismatch_ln)`:
/// `match = ln(1 − ε)`, `mismatch = ln(ε / 3)` with `ε = 10^(−Q/10)` — the
/// Dindel base-quality model (arch §5.3). Built once, process-wide (the BAQ
/// `Q2P` lookup *pattern*). Q0 (ε = 1) is pathological — its match probability
/// would be `ln 0`; we floor it to the smallest positive `f64` so a Q0 match is
/// merely extremely unlikely, not an impossible `-∞` that annihilates the row.
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

/// Insertion emission: an inserted read base is scored against a uniform base
/// distribution, `ln(1/4)` (Durbin et al.'s background `qₐ`). Constant per
/// inserted base, so it only shifts a candidate's score when candidates differ
/// in their insertion count — exactly when it should.
static INS_EMIT_LN: LazyLock<f64> = LazyLock::new(|| 0.25f64.ln());

/// The pair-HMM's log-space transition probabilities, precomputed once. With a
/// constant gap-open these are run-independent; the homopolymer-indexed
/// refinement (deferred) would make the open terms position-dependent.
#[derive(Debug, Clone, Copy)]
pub(crate) struct HmmModel {
    /// Match → Match: `ln(1 − 2·gap_open)`.
    ln_mm: f64,
    /// Match → Insertion and Match → Deletion (gap open): `ln(gap_open)`.
    ln_gap_open: f64,
    /// Insertion → Match and Deletion → Match (gap close): `ln(1 − gap_extend)`.
    ln_gap_close: f64,
    /// Insertion → Insertion and Deletion → Deletion (gap extend):
    /// `ln(gap_extend)`.
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
    /// The default sequencing-error model (constant gap-open, Dindel emission +
    /// gap-extend). The only model for now; a constructor taking calibrated
    /// parameters can be added when the residual-validation pass tunes them.
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

/// Per-worker forward scratch: the previous and current DP rows, each
/// `hap_len + 1` cells of `[M, I, D]` log-probabilities. Grow-and-keep
/// (`resize_for`) so the hot path never allocates (arch §5.7, the BAQ
/// scratch.rs ethos). No full matrix — the forward sum needs only the prior
/// row, so there is nothing to trace back.
#[derive(Debug, Default)]
pub(crate) struct PairHmmScratch {
    prev: Vec<[f64; 3]>,
    cur: Vec<[f64; 3]>,
}

impl PairHmmScratch {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Ensure both rows hold at least `hap_len + 1` cells. Keeps capacity
    /// across calls; only grows.
    fn resize_for(&mut self, hap_len: usize) {
        let needed = hap_len + 1;
        if self.prev.len() < needed {
            self.prev.resize(needed, [f64::NEG_INFINITY; 3]);
            self.cur.resize(needed, [f64::NEG_INFINITY; 3]);
        }
    }
}

/// `ln(e^a + e^b)`, numerically stable, treating `-∞` as a true zero
/// probability.
#[inline]
fn ln_sum_exp2(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// `ln(e^a + e^b + e^c)`.
#[inline]
fn ln_sum_exp3(a: f64, b: f64, c: f64) -> f64 {
    ln_sum_exp2(ln_sum_exp2(a, b), c)
}

/// Forward log-likelihood `ln P(read | candidate_seq)` summed over all
/// alignments, under the sequencing-error model `model`. `read` and `quals` are
/// the read bases and their matching Phred scores (`read.len() == quals.len()`);
/// `candidate_seq` is the full candidate haplotype the read is scored against
/// (upper-cased A/C/G/T bytes). Pure — no traceback, no allocation beyond the
/// reused `scratch`.
///
/// The recursion is the standard global pair-HMM forward with `Match[0][0]`
/// as the begin state; leading/trailing insertions and deletions are reached
/// through the I/D boundary cells. Returns `f64::NEG_INFINITY` only for the
/// degenerate empty-haplotype-and-empty-read corner via the normal recursion.
pub(crate) fn forward(
    read: &[u8],
    quals: &[u8],
    candidate_seq: &[u8],
    scratch: &mut PairHmmScratch,
    model: &HmmModel,
) -> f64 {
    debug_assert_eq!(read.len(), quals.len(), "read and quals length mismatch");
    let m = read.len();
    let n = candidate_seq.len();
    scratch.resize_for(n);
    let ins_emit = *INS_EMIT_LN;
    let emission = &*EMISSION_LN;

    // Row 0 (no read base consumed): begin in Match at the corner; the only
    // reachable cells along the row are deletions of leading haplotype bases.
    {
        let prev = &mut scratch.prev;
        prev[0] = [0.0, f64::NEG_INFINITY, f64::NEG_INFINITY];
        for j in 1..=n {
            let d = ln_sum_exp2(
                model.ln_gap_open + prev[j - 1][M],
                model.ln_gap_extend + prev[j - 1][D],
            );
            prev[j] = [f64::NEG_INFINITY, f64::NEG_INFINITY, d];
        }
    }

    // Rows 1..=m: consume one read base per row.
    for i in 1..=m {
        let read_base = read[i - 1];
        let (match_ln, mismatch_ln) = emission[quals[i - 1] as usize];

        // Column 0: only an insertion (read base before any haplotype base).
        let cur0_i = ins_emit
            + ln_sum_exp2(
                model.ln_gap_open + scratch.prev[0][M],
                model.ln_gap_extend + scratch.prev[0][I],
            );
        scratch.cur[0] = [f64::NEG_INFINITY, cur0_i, f64::NEG_INFINITY];

        for j in 1..=n {
            let emit = if read_base == candidate_seq[j - 1] {
                match_ln
            } else {
                mismatch_ln
            };
            let m_cell = emit
                + ln_sum_exp3(
                    model.ln_mm + scratch.prev[j - 1][M],
                    model.ln_gap_close + scratch.prev[j - 1][I],
                    model.ln_gap_close + scratch.prev[j - 1][D],
                );
            let i_cell = ins_emit
                + ln_sum_exp2(
                    model.ln_gap_open + scratch.prev[j][M],
                    model.ln_gap_extend + scratch.prev[j][I],
                );
            let d_cell = ln_sum_exp2(
                model.ln_gap_open + scratch.cur[j - 1][M],
                model.ln_gap_extend + scratch.cur[j - 1][D],
            );
            scratch.cur[j] = [m_cell, i_cell, d_cell];
        }

        std::mem::swap(&mut scratch.prev, &mut scratch.cur);
    }

    let last = &scratch.prev[n];
    ln_sum_exp3(last[M], last[I], last[D])
}

/// Score a read against a pre-built candidate set, returning the **dense** `Qᵣ`:
/// one `(allele, forward log-likelihood)` per candidate, in candidate order.
///
/// This is the per-read join of candidate generation (arch §6) and the forward
/// (arch §5): each candidate's full sequence is handed to [`forward`] and its raw
/// log-likelihood read back, reusing the one `scratch` across the set. Pruning to
/// a sparse profile (`AMB_LL_DROP`) and renormalization are the **aggregator's**
/// job (`locus_record`, arch §11), so this returns the raw scores untouched — the
/// dense distribution over the window the aggregator then sparsifies.
pub(crate) fn score_candidates(
    read: &[u8],
    quals: &[u8],
    candidates: &[CandidateAllele],
    scratch: &mut PairHmmScratch,
    model: &HmmModel,
) -> Vec<(Allele, f64)> {
    candidates
        .iter()
        .map(|c| {
            let loglik = forward(read, quals, &c.candidate_seq, scratch, model);
            (c.allele.clone(), loglik)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A model and a fresh scratch for a test.
    fn setup() -> (HmmModel, PairHmmScratch) {
        (HmmModel::new(), PairHmmScratch::new())
    }

    fn assert_close(got: f64, want: f64) {
        assert!(
            (got - want).abs() < 1e-9,
            "got {got}, want {want} (diff {})",
            (got - want).abs()
        );
    }

    #[test]
    fn single_base_exact_match_is_hand_computable() {
        let (model, mut s) = setup();
        // read "A" @ Q40 vs hap "A": the only surviving path is begin → Match,
        // so the score is ln(1 − ε₄₀) + ln(1 − 2·gap_open).
        let got = forward(b"A", &[40], b"A", &mut s, &model);
        let eps = 10f64.powf(-4.0);
        let want = (1.0 - eps).ln() + (1.0 - 2.0 * GAP_OPEN_PROB).ln();
        assert_close(got, want);
    }

    #[test]
    fn single_base_mismatch_uses_the_mismatch_emission() {
        let (model, mut s) = setup();
        let got = forward(b"A", &[40], b"C", &mut s, &model);
        let eps = 10f64.powf(-4.0);
        let want = (eps / 3.0).ln() + (1.0 - 2.0 * GAP_OPEN_PROB).ln();
        assert_close(got, want);
    }

    #[test]
    fn a_match_scores_far_higher_than_a_mismatch() {
        let (model, mut s) = setup();
        let m = forward(b"A", &[40], b"A", &mut s, &model);
        let mm = forward(b"A", &[40], b"C", &mut s, &model);
        assert!(m > mm, "match {m} should beat mismatch {mm}");
    }

    #[test]
    fn identical_haplotype_beats_a_one_unit_longer_one() {
        let (model, mut s) = setup();
        let read = b"CACACA"; // 3 units of CA
        let quals = [40u8; 6];
        let exact = forward(read, &quals, b"CACACA", &mut s, &model);
        let longer = forward(read, &quals, b"CACACACA", &mut s, &model); // 4 units
        assert!(
            exact > longer,
            "exact-length hap {exact} should beat the 2-bp-longer hap {longer}"
        );
    }

    #[test]
    fn longer_exact_match_carries_more_certainty() {
        // A longer all-match alignment multiplies in more (1 − ε) match factors,
        // so its log-likelihood is more negative than a short one — but both are
        // dominated by the match path; we check the per-length ordering holds.
        let (model, mut s) = setup();
        let short = forward(b"ACGT", &[40; 4], b"ACGT", &mut s, &model);
        let long = forward(b"ACGTACGT", &[40; 8], b"ACGTACGT", &mut s, &model);
        // More bases ⇒ more (1 − ε) factors ⇒ smaller (more negative) log-lik.
        assert!(short > long, "short {short} vs long {long}");
        // Both remain valid finite log-probabilities (≤ 0).
        assert!(short <= 0.0 && long <= 0.0);
    }

    #[test]
    fn low_quality_mismatch_is_less_penalized_than_high_quality() {
        // A mismatch at low Q is more forgivable (ε large ⇒ ε/3 larger) than at
        // high Q, so its score is higher.
        let (model, mut s) = setup();
        let lowq = forward(b"A", &[5], b"C", &mut s, &model);
        let highq = forward(b"A", &[40], b"C", &mut s, &model);
        assert!(lowq > highq, "low-Q mismatch {lowq} vs high-Q {highq}");
    }

    #[test]
    fn empty_read_scores_an_all_deletion_path() {
        // P(empty read | hap) is the all-deletion path: a finite, ≤ 0 log-prob.
        let (model, mut s) = setup();
        let got = forward(b"", &[], b"ACGT", &mut s, &model);
        assert!(got.is_finite() && got < 0.0, "empty-read score {got}");
    }

    #[test]
    fn scratch_is_reused_across_reads_of_growing_haplotypes() {
        // resize_for only grows; scoring a short hap after a long one must still
        // be correct (stale cells beyond `n` are never read).
        let (model, mut s) = setup();
        let _long = forward(b"ACGTACGT", &[40; 8], b"ACGTACGTACGT", &mut s, &model);
        let short = forward(b"A", &[40], b"A", &mut s, &model);
        let eps = 10f64.powf(-4.0);
        let want = (1.0 - eps).ln() + (1.0 - 2.0 * GAP_OPEN_PROB).ln();
        assert_close(short, want);
    }

    /// A perfect (CA) locus: GGG | CACACA | TTT, tract at [13, 19).
    fn ca_locus() -> crate::ssr::types::Locus {
        use crate::ssr::types::{Locus, Motif};
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

    #[test]
    fn score_candidates_returns_one_dense_score_per_candidate_in_order() {
        use crate::ssr::pileup::candidate_generation::build_rungs;
        let locus = ca_locus();
        let mut cands = Vec::new();
        build_rungs(&locus, 3, 2, &mut cands); // L ∈ 1..=5 → 5 candidates
        let (model, mut s) = setup();
        // Any read; we only check the shape + ordering here.
        let scored = score_candidates(b"GGGCACACATTT", &[40; 12], &cands, &mut s, &model);
        assert_eq!(scored.len(), cands.len());
        let alleles: Vec<&Allele> = scored.iter().map(|(a, _)| a).collect();
        let want: Vec<&Allele> = cands.iter().map(|c| &c.allele).collect();
        assert_eq!(alleles, want, "dense scores follow candidate order");
    }

    #[test]
    fn score_candidates_ranks_the_matching_rung_highest() {
        use crate::ssr::pileup::candidate_generation::build_rungs;
        let locus = ca_locus();
        let mut cands = Vec::new();
        build_rungs(&locus, 3, 2, &mut cands); // L ∈ 1..=5
        let (model, mut s) = setup();
        // A high-quality read equal to the L=3 rung (== ref_bytes): it must win.
        let scored = score_candidates(b"GGGCACACATTT", &[40; 12], &cands, &mut s, &model);
        let best = scored
            .iter()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();
        assert_eq!(best.0, Allele::OnLadder { units: 3 });
    }
}
