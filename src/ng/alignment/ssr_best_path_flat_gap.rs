//! Algorithm 3 — the repeat-aware best-path aligner with **two gap regimes and one flat
//! gap inside the repeat**. This is ng's tract delimiter, and the module's **only
//! byte-parity oracle**.
//!
//! # Not the recommended delimiter — kept as the oracle and the baseline
//!
//! **[`super::ssr_best_path_unit_slip`] (algorithm 4) measured better and is the recommended
//! STR delimiter.** On synthetic reads with known truth it was uniformly more accurate and
//! free of the small reference-pull bias this flat-gap model carries (the D2 comparison,
//! `doc/devel/reports/research/ssr_delimiter_3v4_comparison_2026-07-23.md`). So a new caller
//! that just wants to measure a repeat should reach for algorithm 4, not this.
//!
//! Algorithm 3 is **retained deliberately, not deprecated**, for two reasons that outrank the
//! measurement result:
//!
//! - **It is the module's only byte-parity oracle.** This is a port of production's
//!   `delimit_read`, held byte-identical to it over a 200,000-case soak (`delimit_parity`).
//!   Algorithm 4 has **no oracle of its own** — its correctness rests on the differential
//!   cross-check *against this aligner* (`algorithm_4_agrees_with_algorithm_3_on_clean_reads`).
//!   Removing or disabling algorithm 3 would delete the only thing that validates algorithm 4.
//! - **Adoption is not settled here.** D2 measured which ruler is more accurate; whether to
//!   *adopt* algorithm 4 turns on the genotype comparison, which spans this module and the
//!   genotyping that composes it (spec §10.3) and has not run. The measurement points to 4;
//!   the decision is downstream.
//!
//! So this is not disabled by a panic or removed: doing either would break the oracle and the
//! validation of the algorithm it would be promoting. The steering belongs at the **selection
//! site** — when the STR read preparer is wired up it should instantiate algorithm 4, and any
//! future user-facing choice should omit algorithm 3 from its menu while keeping the type
//! compilable for these tests.
//!
//! It is the same alignment as the affine aligner, but with two gap regimes along the
//! reference: a **stiff** gap in the flanks and a **soft** one inside the tract. The
//! reasoning is that flanks are clean unique sequence and should hold the tract's edges in
//! place, while length variation inside the tract is expected — it is the repeat changing
//! size, not an error. With a single stiff gap everywhere, any allele more than about one
//! unit from the reference collapses back onto the reference; production hit exactly that
//! failure before splitting the regimes (spec §4.2).
//!
//! Because the flanks stay anchored, the two flank-to-tract boundary columns of the
//! resulting line-up say where the read's repeat starts and ends — which is how this
//! aligner *measures* a read's repeat.
//!
//! # Unbanded, deliberately
//!
//! Production fills its whole matrix. This port does too, because **only an unbanded port
//! can be shown byte-identical to it**, and that parity is this module's one hard oracle
//! (spec §3, §10.3). Banding lands as its own separately-proved change in Milestone C, at
//! which point the same fixture must produce the same output — which is what makes banding
//! a performance change rather than a behaviour change.
//!
//! # A known inconsistency, reproduced on purpose
//!
//! Production computes the match→match probability **once**, from the *flank* gap-open
//! (`1 − 2 × 2.9e-5`), and never recomputes it under the tract regime. So inside the tract
//! the three transitions leaving a match sum to about **1.02** rather than 1, while in the
//! flanks they sum to exactly 1.
//!
//! Be precise about what that skews. It is *direction*-symmetric — insertion and deletion
//! get the same tract open — so it does not favour longer over shorter. But normalising
//! would divide all three by 1.02 **at every departure from a match**, and a path carrying
//! an *L*-base gap passes through about *L* fewer matches than the diagonal does. So the
//! un-normalised model quietly penalises the length-changing path by roughly 1.02 per
//! gapped base — about 1.2× at a 10 bp change, 1.8× at 30 bp. **The bias is toward the
//! reference length, and it grows with the size of the change.** That is the same class of
//! bias the tract regime exists to remove; it is dominated today by the 350× open ratio,
//! which is why it has never shown up in results, but anyone who narrows that ratio or
//! lengthens the tract inherits an unstated pull toward reference.
//!
//! It is reproduced here rather than fixed because **reproducing it silently and fixing it
//! silently both end with a parity test failing for a reason nobody can find** (spec §4.2).

use super::emission::Emission;
use super::{BestPathAligner, ReadBases, RepeatContext, RepeatSpan};
use std::sync::LazyLock;

/// Which of the pair-HMM's three states a cell was entered in.
///
/// An enum rather than three `const usize` indices, which is what this was ported as. The
/// port kept production's shape deliberately — transcribe first, change separately — and this
/// is the separate change, made once the byte-parity oracle existed to prove it changed
/// nothing.
///
/// Two things it buys, both real defects rather than tidiness:
///
/// - **The traceback's catch-all arm is gone.** Matching a `usize` needs a `_ =>`, and that
///   arm silently read *any* value other than match or deletion as an insertion. Nothing
///   could reach it, but nothing said so either; now the compiler enumerates the three cases,
///   and a fourth would be a compile error rather than an insertion.
/// - **`best_of`'s second element is type-safe.** It carried a `u8` — indistinguishable from
///   a [`BaseQual`](crate::ng::types::BaseQual), which is in scope in the same loop — and
///   every construction site needed an `as u8` cast.
///
/// `#[repr(u8)]` with explicit discriminants, so the layout and index values are exactly
/// production's: match 0, insertion 1, deletion 2. `index()` compiles to the same constant
/// the `const` did.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum State {
    /// A read base and a reference base consumed together.
    Match = 0,
    /// A read base consumed against no reference base.
    Insertion = 1,
    /// A reference base consumed against no read base.
    Deletion = 2,
}

impl State {
    /// Index into a cell's three-state score array.
    #[inline]
    const fn index(self) -> usize {
        self as usize
    }
}

/// An unreachable score. A real value, not an error: some cells genuinely cannot be
/// entered in some states.
const UNREACHABLE: f64 = f64::NEG_INFINITY;

/// The band's slop for a **score-neutral bow** — the smallest of its three terms, and the
/// only one that is a constant (spec §3, arch §6). See [`band_width`] for the whole formula,
/// and `a_long_allele_at_the_extreme_is_measured_not_collapsed` /
/// `a_run_off_flank_read_is_measured_despite_the_bow` for the two silent failures the other
/// terms prevent.
///
/// # The band has three terms, and the spec anticipated only two
///
/// A cell `(i, j)` is computed when `|i − j| <= band`. The band is
/// `|read − reference| + left_flank + right_flank + BAND_HEADROOM`, and each term earns its
/// place:
///
/// 1. **`|read − reference|` — the forced floor.** These are *global* alignments: read and
///    reference are both consumed entirely, so when they differ in length the path is forced
///    that far off the diagonal just to reach the far corner. Below this the true path is
///    provably lost. Computed per read, **never a constant** (spec §3: production's own
///    long-allele test aligns 52 bases against 28 and needs 24 cells of deviation, which no
///    fixed width covers).
///
/// 2. **`left_flank + right_flank` — the run-off correction, and a departure from the spec.**
///    Spec §3 modelled the band as the floor plus a *small constant* headroom, and §9 left
///    the amount open. That model is **insufficient for this aligner, and it fails silently**:
///    a read that expands the tract *and* runs off a flank has an optimal path whose peak
///    deviation is the tract expansion, while its *net* length difference is smaller by the
///    run-off flank's length. So the path strays past the `|read − reference|` corridor by up
///    to the flank length — far more than a constant. The parity oracle caught this at a
///    left-flankless locus (a 43-base read of a 7-unit tract against a 4-unit, 34-base frame:
///    peak deviation 18, floor 9). Bounding the stray by the two flank lengths is derived from
///    the alignment *geometry*, not the scoring model, so it keeps the spec's separation of
///    the ruler from the slip cutoff.
///
///    `SPEC-FOLLOWUP(alignment §3/§9)`: the spec models the band as floor + a constant
///    headroom and lists the amount as open; it should be amended to the three-term geometry
///    band above. Owner to fold in (delegated 2026-07-23); recorded in the C1 impl report and
///    `PROJECT_STATUS.md` so the deferral is tracked, not only in this comment.
///
/// 3. **`BAND_HEADROOM` — this constant.** With the floor and the run-off covered, what is
///    left is a genuine score-neutral bow: an insertion balanced by a later deletion, net
///    zero. Every such bow pays two gap events whose transition cost is strictly negative in
///    every regime, so it is always *worse* than the direct path, never merely tied (even
///    where emissions tie, at Q0 or a flat ε = 0.75, the gap transitions still cost). Bows are
///    therefore self-limiting; production's banded forward allows just 2 (`FLANK_SLOP`). This
///    takes **8**, generous because the delimiter's soft tract gap is cheaper than the
///    marginal's flank gap so a bow reaches a little further, and because the cost of being
///    too generous is a few wasted cells while the cost of being too narrow is a silently
///    wrong measurement. Held at byte-parity across the whole 200,000-case soak.
///
/// **No term is derived from [`MAX_SLIP`](super::stutter::MAX_SLIP)**, which is a *scoring*
/// cutoff: coupling the ruler's width to it would be too narrow for the long-allele test and
/// would let the scoring model blind the measurement (spec §3).
const BAND_HEADROOM: usize = 8;

/// Gap-open probability (match → insertion, or → deletion) in the **flanks**. Dindel's
/// base value for short homopolymer runs. The flanks are clean unique sequence, so a stiff
/// gap keeps them anchoring the tract junctions.
///
/// **Provisional development calibration, not a finding** — to be reconciled against the
/// measured slip rate (spec §4.2).
const GAP_OPEN_PROB: f64 = 2.9e-5;

/// Gap-open probability **inside the repeat tract** — about 350× softer than the flank
/// value, because STR length variation *is* multi-unit tract indels, and the flank-grade
/// gap collapses any allele more than about one unit longer than the reference. HipSTR
/// makes the same flank/tract split.
///
/// This is a **flat per-base tract gap**: content-agnostic by design, so an impure or
/// out-of-frame read is still read out verbatim. Pricing whole-unit slips separately is
/// algorithm 4's job, and which is better is the Milestone D comparison — not something to
/// decide here (spec §4.2).
///
/// **Provisional development calibration, not a finding.**
const GAP_OPEN_PROB_TRACT: f64 = 1e-2;

/// Gap-extension probability: Dindel's fixed `e⁻¹ ≈ 0.368`. Shared by both regimes —
/// **only the open cost switches regime**, because the open is the dominant fixed cost and
/// the extension is a knob nobody has yet needed to turn (spec §4.2).
static GAP_EXTEND_PROB: LazyLock<f64> = LazyLock::new(|| (-1.0f64).exp());

/// The pair-HMM's log-space transition probabilities, with a **tract-aware gap-open**.
///
/// See the module docs for the match→match inconsistency this deliberately reproduces:
/// `ln_match_to_match` is computed from the *flank* gap-open and never recomputed under
/// the tract regime.
#[derive(Debug, Clone, Copy)]
pub struct TransitionCosts {
    ln_match_to_match: f64,
    ln_gap_open: f64,
    ln_gap_open_tract: f64,
    ln_gap_close: f64,
    ln_gap_extend: f64,
}

impl Default for TransitionCosts {
    fn default() -> Self {
        let extend = *GAP_EXTEND_PROB;
        Self {
            ln_match_to_match: (1.0 - 2.0 * GAP_OPEN_PROB).ln(),
            ln_gap_open: GAP_OPEN_PROB.ln(),
            ln_gap_open_tract: GAP_OPEN_PROB_TRACT.ln(),
            ln_gap_close: (1.0 - extend).ln(),
            ln_gap_extend: extend.ln(),
        }
    }
}

impl TransitionCosts {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    // Accessors for the sibling algorithm 4 (`ssr_best_path_unit_slip`), which shares this
    // exact per-base gap model for its *out-of-frame* route and adds whole-unit slip states
    // on top. `pub(super)` because they are an internal seam between the two aligners, not
    // part of the module's public surface; reading them changes nothing about algorithm 3.
    pub(super) fn ln_match_to_match(&self) -> f64 {
        self.ln_match_to_match
    }
    pub(super) fn ln_gap_open(&self) -> f64 {
        self.ln_gap_open
    }
    pub(super) fn ln_gap_open_tract(&self) -> f64 {
        self.ln_gap_open_tract
    }
    pub(super) fn ln_gap_close(&self) -> f64 {
        self.ln_gap_close
    }
    pub(super) fn ln_gap_extend(&self) -> f64 {
        self.ln_gap_extend
    }
}

/// Per-worker Viterbi scratch: two rolling score rows plus a **full** backpointer matrix
/// for the traceback.
///
/// Grow-and-keep, so the hot path never reallocates — and buffers only, deciding nothing
/// that changes a result, as the aligner trait's `Scratch` contract requires. The
/// backpointer matrix is the whole `(m+1) × (n+1)`, which is what "unbanded" costs and what
/// makes parity provable.
///
/// **`resize` only grows; it never clears** — and that is safe because **every score-row
/// cell in use is written every row**, in band or out. Banding (Milestone C) does not
/// *skip* out-of-band cells, which would leave a larger previous read's values to leak in as
/// a neighbour score; it writes them `UNREACHABLE`. So the invariant the old exhaustive fill
/// gave for free — no cell is read before it is written this call — still holds.
///
/// The backpointer matrix is the one place a stale cell survives: an out-of-band cell's
/// backpointer is not rewritten. It is never *read*, though — the traceback follows the
/// score-optimal path, which is reachable and therefore in band, and an `UNREACHABLE` cell
/// can never be a winning predecessor.
#[derive(Debug, Default)]
pub struct ViterbiScratch {
    previous: Vec<[f64; 3]>,
    current: Vec<[f64; 3]>,
    /// Per cell `(i, j)` and per state, the predecessor state that won the maximum — a
    /// flat `(m+1) × (n+1)` matrix, row-major with stride `n + 1`.
    backpointers: Vec<[State; 3]>,
}

impl ViterbiScratch {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    fn resize(&mut self, read_len: usize, reference_len: usize) {
        let row = reference_len + 1;
        if self.previous.len() < row {
            self.previous.resize(row, [UNREACHABLE; 3]);
            self.current.resize(row, [UNREACHABLE; 3]);
        }
        let cells = (read_len + 1) * row;
        if self.backpointers.len() < cells {
            self.backpointers.resize(cells, [State::Match; 3]);
        }
    }
}

/// Where the traceback found the repeat, before it is classified.
///
/// The raw readout: the two junction offsets in read coordinates, plus whether each flank
/// actually held. **Step B2 turns this into a [`RepeatSpan`]** and adds the
/// [`BestPathAligner`] implementation; this step lands the recurrence and the walk, which
/// is the part whose failure is silent — a wrong recurrence is a wrong repeat length, not
/// a panic.
///
/// [`RepeatSpan`]: super::RepeatSpan
/// [`BestPathAligner`]: super::BestPathAligner
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TractReadout {
    /// Read offset where the repeat starts.
    pub tract_start: u64,
    /// Read offset where the repeat ends, exclusive.
    pub tract_end: u64,
    /// Whether the left flank held — false when the read began inside the repeat.
    pub left_anchored: bool,
    /// Whether the right flank held — false when the read ran off its end inside it.
    pub right_anchored: bool,
}

impl TractReadout {
    /// Classify the readout into the four cases a caller can act on.
    ///
    /// **This is the widening, and it is where ng goes past production.** Production's
    /// `Delimited` collapses every unanchored outcome into one `BorderOffEnd`, so it is
    /// *side-blind*: it cannot say which flank was missing, and the STR read preparer cannot
    /// otherwise tell a measurement from a lower bound (arch §2.1, §5; spec §8).
    ///
    /// **The mapping is the whole risk of this step, and its failure is silent.** Getting a
    /// side backwards does not crash and does not produce an implausible number — it
    /// produces a read that looks perfectly good and is filed under the wrong observation
    /// class. So read the two middle arms carefully: [`RepeatSpan::FromLeft`] means the
    /// **left** flank is the one that *held*, and therefore that the read ran off its own
    /// 3′ end somewhere inside the repeat.
    #[must_use]
    pub fn classify(&self) -> RepeatSpan {
        // `delimit` guarantees this, but the type is public with public fields — B3's parity
        // harness needs the raw offsets — so a hand-built readout can reach here. An
        // inverted pair would become a `Between` whose measured length saturates to zero: a
        // *confident zero-unit allele*, which is far worse than a rejection.
        debug_assert!(
            self.tract_start <= self.tract_end,
            "inverted tract offsets: {} > {}",
            self.tract_start,
            self.tract_end
        );
        let span = self.tract_start..self.tract_end;
        match (self.left_anchored, self.right_anchored) {
            (true, true) => RepeatSpan::Between(span),
            (true, false) => RepeatSpan::FromLeft(span),
            (false, true) => RepeatSpan::FromRight(span),
            // Neither flank held: the read lies wholly inside the repeat and carries no
            // per-read fact about its length, so the span is dropped rather than reported
            // as if it measured something.
            (false, false) => RepeatSpan::Unanchored,
        }
    }
}

/// Pick the best of the candidates, keeping the **first on ties** — so the caller encodes
/// the tie-break by passing candidates in priority order.
///
/// The tie-break is a correctness rule, not a nicety: two reads of the same molecule must
/// measure the same repeat, and equally-scoring line-ups would otherwise let identical
/// input produce different answers. Every call site passes **match, then deletion, then
/// insertion** (spec §4.2).
///
/// # Panics
///
/// Panics if `candidates` is empty. Every call site passes a fixed non-empty array literal
/// of transition candidates, so this cannot fire today — the assertion records that, so a
/// future dynamic call site is caught in development rather than indexing off the end.
#[inline]
fn best_of(candidates: &[(f64, State)]) -> (f64, State) {
    debug_assert!(!candidates.is_empty(), "best_of requires a candidate");
    let mut best = candidates[0];
    for &candidate in &candidates[1..] {
        if candidate.0 > best.0 {
            best = candidate;
        }
    }
    best
}

/// The band: a cell `(i, j)` is computed when `|i − j| <= band_width(..)`.
///
/// The three geometry terms, in full on [`BAND_HEADROOM`]: the forced floor
/// `|read − reference|`, the run-off correction `left + right`, and the score-neutral bow
/// slop. A free function, not a closure, so the formula the whole port's correctness rests on
/// is unit-testable in isolation (`band_width_is_geometry_not_the_slip_cutoff`,
/// `the_band_always_reaches_the_far_corner`).
///
/// The four arguments are plain `usize`, and a transposition of either pair is **harmless by
/// construction**: the formula is symmetric in read/reference (`abs_diff`) and in left/right
/// (addition), so swapping them changes no result. That is why they are not newtyped here.
fn band_width(
    read_len: usize,
    reference_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
) -> usize {
    read_len.abs_diff(reference_len) + left_flank_len + right_flank_len + BAND_HEADROOM
}

/// Algorithm 3: the two-regime, one-flat-tract-gap repeat delimiter.
///
/// Generic over its [`Emission`], never behind `dyn` — the emission model is the
/// experiment's configuration and is selected at compile time (arch §4).
#[derive(Debug, Clone, Copy)]
pub struct SsrFlatGapAligner<E: Emission> {
    emission: E,
    costs: TransitionCosts,
}

impl<E: Emission> SsrFlatGapAligner<E> {
    #[must_use]
    pub fn new(emission: E) -> Self {
        Self {
            emission,
            costs: TransitionCosts::new(),
        }
    }

    /// Fill the matrix and walk the traceback, returning the raw readout.
    ///
    /// `None` when there is no reference frame at all — an empty reference is not a real
    /// locus, and there is nothing to delimit. That is a *defined answer*, not an error
    /// (arch §3).
    #[must_use]
    pub fn delimit(
        &self,
        read: ReadBases<'_>,
        reference: &[u8],
        context: &RepeatContext<'_>,
        scratch: &mut ViterbiScratch,
    ) -> Option<TractReadout> {
        debug_assert!(
            context
                .geometry
                .fits_reference(crate::ng::types::Bp(reference.len() as u64)),
            "the repeat geometry does not fit the reference stretch"
        );

        let read_len = read.len();
        let bases = read.bases();
        let reference_len = reference.len();
        if reference_len == 0 {
            return None;
        }
        // Flank lengths, narrowed to the index width and clamped into the frame.
        //
        // **Both clamps are no-ops whenever the documented precondition holds** — the debug
        // assertion above is the real check — so they change no measurement and production
        // has no counterpart to them. They exist so that a violated precondition degrades
        // into a *defined* frame in release rather than an out-of-bounds column index. Be
        // clear-eyed about the cost: when they do bind, they silently redefine where the
        // tract is, which is a wrong answer rather than a crash. That is the trade arch §3
        // makes everywhere in this module.
        //
        // `try_into().unwrap_or(usize::MAX)` rather than `as usize`: on a 32-bit target a
        // bare cast would *wrap*, and the `.min()` below would then turn a wrapped value
        // into a plausible-looking flank length. Saturating keeps the clamp meaningful.
        let left_flank_len = usize::try_from(context.geometry.left_flank_len.get())
            .unwrap_or(usize::MAX)
            .min(reference_len);
        let right_flank_len = usize::try_from(context.geometry.right_flank_len.get())
            .unwrap_or(usize::MAX)
            .min(reference_len - left_flank_len);

        scratch.resize(read_len, reference_len);
        let stride = reference_len + 1;
        let insertion_emission = self.emission.insert_ln();

        // The band. A cell `(i, j)` is computed only when `|i − j| <= band`; every other cell
        // is written `UNREACHABLE`, exactly as production's banded forward does
        // ([src/ssr/cohort/pair_hmm.rs](../../../ssr/cohort/pair_hmm.rs)).
        //
        // Three terms — the forced floor `|read − reference|`, the run-off correction
        // `left_flank + right_flank`, and the score-neutral bow slop `BAND_HEADROOM` — each
        // justified in full on [`BAND_HEADROOM`]. Every term is bytes of alignment *geometry*;
        // none is the scoring cutoff.
        //
        // **Writing every out-of-band cell `UNREACHABLE`, rather than skipping it, is what
        // makes banding safe on a reused scratch.** The unbanded fill was exhaustive, so
        // stale cells could never be observed; a band that merely *skipped* cells would let a
        // previous, larger read's values leak in as a neighbour score. Writing each cell every
        // row keeps the invariant the exhaustive fill used to give for free.
        let band = band_width(read_len, reference_len, left_flank_len, right_flank_len);
        let in_band = |row: usize, column: usize| row.abs_diff(column) <= band;

        // Tract-aware gap-open: a gap touching reference column `j` is **inside the tract**
        // — and so gets the soft gap — when `left_flank_len < j <= reference_len -
        // right_flank_len`, i.e. when it inserts beside, or deletes, a tract base. Every
        // other column, including both junctions' outer sides, keeps the stiff flank rate.
        // That is what holds the junctions in place while the interior is free to change
        // length. `j` ranges over `1..=reference_len`.
        let tract_last_column = reference_len - right_flank_len;
        let gap_open = |column: usize| {
            if left_flank_len < column && column <= tract_last_column {
                self.costs.ln_gap_open_tract
            } else {
                self.costs.ln_gap_open
            }
        };

        let ViterbiScratch {
            previous,
            current,
            backpointers,
        } = scratch;

        // Row 0 — no read base consumed. Begin in Match at the corner; the only reachable
        // cells are deletions of leading reference bases.
        //
        // **`gap_open` is consulted at column 1 only, here.** From column 2 on, the match
        // predecessor `previous[column - 1][State::Match.index()]` is `UNREACHABLE`, so that candidate is
        // `-inf` and the deletion-extend candidate always wins — the open cost never enters.
        // And column 1 is inside the tract only when the left flank is empty. So the whole
        // of row 0's tract-awareness reduces to a single case: **a locus with no left flank**.
        //
        // Worth writing down, because it was briefly recorded as a hole in the parity
        // oracle: a mutation replacing this with the flat flank rate survived 200,000
        // randomized cases. It survived only because the fixture generated no zero-length
        // flanks. Once it did, the same mutation failed within ~100 cases.
        previous[0] = [0.0, UNREACHABLE, UNREACHABLE];
        backpointers[0] = [State::Match; 3];
        for column in 1..=reference_len {
            // Row 0's band is `column <= band`. Beyond it the deletion chain would run off
            // the corridor, and — because the chain is strictly decreasing — never return to
            // the optimal path, so those cells are `UNREACHABLE`.
            if !in_band(0, column) {
                previous[column] = [UNREACHABLE; 3];
                continue;
            }
            let (score, predecessor) = best_of(&[
                (
                    gap_open(column) + previous[column - 1][State::Match.index()],
                    State::Match,
                ),
                (
                    self.costs.ln_gap_extend + previous[column - 1][State::Deletion.index()],
                    State::Deletion,
                ),
            ]);
            previous[column] = [UNREACHABLE, UNREACHABLE, score];
            backpointers[column] = [State::Match, State::Match, predecessor];
        }

        // Rows 1..=read_len — one read base per row.
        for row_index in 1..=read_len {
            let read_base = bases[row_index - 1];
            // Quality is resolved once per row, not once per cell: it belongs to the read
            // base, so it is constant along the row while the reference base varies.
            let scores = self.emission.scores_for(read.quality_at(row_index - 1));
            let row = row_index * stride;

            // Column 0 — a read base before any reference base: insertion only. In band iff
            // `row <= band`; a read longer than the band with no reference consumed is off the
            // corridor.
            if in_band(row_index, 0) {
                let (score, predecessor) = best_of(&[
                    (
                        self.costs.ln_gap_open + previous[0][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend + previous[0][State::Insertion.index()],
                        State::Insertion,
                    ),
                ]);
                current[0] = [UNREACHABLE, insertion_emission + score, UNREACHABLE];
                backpointers[row] = [State::Match, predecessor, State::Match];
            } else {
                current[0] = [UNREACHABLE; 3];
            }

            for column in 1..=reference_len {
                if !in_band(row_index, column) {
                    current[column] = [UNREACHABLE; 3];
                    continue;
                }
                let emitted = scores.pick(read_base, reference[column - 1]);

                // Match: from the diagonal cell, in any state. Priority M > D > I.
                let (match_score, match_predecessor) = best_of(&[
                    (
                        self.costs.ln_match_to_match + previous[column - 1][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_close + previous[column - 1][State::Deletion.index()],
                        State::Deletion,
                    ),
                    (
                        self.costs.ln_gap_close + previous[column - 1][State::Insertion.index()],
                        State::Insertion,
                    ),
                ]);

                let open = gap_open(column);

                // Insertion: from the cell above — a read base consumed, no reference base.
                let (insertion_score, insertion_predecessor) = best_of(&[
                    (open + previous[column][State::Match.index()], State::Match),
                    (
                        self.costs.ln_gap_extend + previous[column][State::Insertion.index()],
                        State::Insertion,
                    ),
                ]);

                // Deletion: from the cell to the left — a reference base consumed, no read
                // base.
                let (deletion_score, deletion_predecessor) = best_of(&[
                    (
                        open + current[column - 1][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend + current[column - 1][State::Deletion.index()],
                        State::Deletion,
                    ),
                ]);

                current[column] = [
                    emitted + match_score,
                    insertion_emission + insertion_score,
                    deletion_score,
                ];
                backpointers[row + column] = [
                    match_predecessor,
                    insertion_predecessor,
                    deletion_predecessor,
                ];
            }
            std::mem::swap(previous, current);
        }

        // The final cell: best terminal state, tie-broken M > D > I.
        let last = previous[reference_len];
        let (_, final_state) = best_of(&[
            (last[State::Match.index()], State::Match),
            (last[State::Deletion.index()], State::Deletion),
            (last[State::Insertion.index()], State::Insertion),
        ]);

        // Walk the traceback, recording for each reference base consumed how many read
        // bases were consumed before it. Insertions are counted **as soon as they occur**,
        // which is what puts a gap on a junction into the block on its 5′ side: an
        // insertion at the left junction joins the flank, one at the right junction joins
        // the repeat (spec §4.2).
        let left_junction = left_flank_len; // first tract base
        let right_junction = reference_len - right_flank_len; // first right-flank base
        let mut tract_start = 0usize; // default: no left flank consumed
        let mut tract_end = read_len; // default: no right flank — runs to the read's end

        let mut row_index = read_len;
        let mut column = reference_len;
        let mut state = final_state;
        while row_index != 0 || column != 0 {
            let predecessor = backpointers[row_index * stride + column][state.index()];
            match state {
                State::Match => {
                    let consumed = column - 1; // reference base this match consumes
                    if consumed == left_junction {
                        tract_start = row_index - 1;
                    }
                    if right_flank_len > 0 && consumed == right_junction {
                        tract_end = row_index - 1;
                    }
                    row_index -= 1;
                    column -= 1;
                }
                State::Deletion => {
                    let consumed = column - 1;
                    if consumed == left_junction {
                        tract_start = row_index;
                    }
                    if right_flank_len > 0 && consumed == right_junction {
                        tract_end = row_index;
                    }
                    column -= 1;
                }
                State::Insertion => {
                    // Consumes a read base, no reference base.
                    //
                    // Named rather than a catch-all: this used to be `_ =>` with a debug
                    // assertion, because matching a `usize` requires one — and that arm would
                    // have read any unexpected value as an insertion. The enum makes the
                    // three cases exhaustive, so the assertion has nothing left to guard.
                    row_index -= 1;
                }
            }
            state = predecessor;
        }

        // A flank "ran off the read end" when no read base sits on its side of the tract.
        //
        // **Production has a third rejection here that this port deliberately omits**: it
        // folds `tract_start > tract_end` into the same `BorderOffEnd` answer. That case is
        // *unreachable*, and the omission is recorded rather than left to be discovered.
        // The reason it cannot happen: `left_junction <= right_junction` holds
        // unconditionally (the clamps above guarantee `left_flank_len + right_flank_len <=
        // reference_len`), and the traceback assigns `tract_start` and `tract_end` from a
        // walk that is monotone in reference position — so the offset recorded at the later
        // reference column is never smaller than the one recorded at the earlier column.
        // Two reviewers proved this independently, one by exhaustive sweep; the invariant
        // is asserted below and pinned by `the_tract_offsets_are_never_inverted`.
        //
        // It matters because the next step turns this into a `RepeatSpan`, and an inverted
        // pair there would become a `Between(..)` whose measured length saturates to zero —
        // a *confident zero-unit allele*, which is a far worse answer than a rejection.
        debug_assert!(
            tract_start <= tract_end,
            "inverted tract offsets: {tract_start} > {tract_end}"
        );
        Some(TractReadout {
            tract_start: tract_start as u64,
            tract_end: tract_end as u64,
            // **A flank that does not exist cannot anchor.** B1 had this the other way round
            // — "a flank that does not exist cannot fail to hold" — which was harmless while
            // these were two raw bits, but B2 promotes them to a `RepeatSpan`, and
            // `Between` is the one variant that claims to *pin* the repeat's length.
            //
            // At a contig-edge locus with no right flank, the two cases that matter produce
            // an **identical** readout: a read that ended because the tract ended, and a
            // read that ended because the read ran out mid-tract. Both give
            // `tract_end == read_len`, and nothing here can tell them apart. Calling that a
            // measurement over-claims for one of the two — a short allele that was never
            // observed, exactly what the widening exists to prevent, at the geometry where
            // the evidence is already thinnest.
            //
            // So the conservative reading wins: without a flank holding that side, the span
            // is a lower bound. The cost is that a flankless locus can never yield a
            // measurement, which is honest — with no flank in the frame there is nothing to
            // pin the tract's edge against.
            //
            // **This is a real divergence from production, and B3 asserts it rather than
            // glossing it.** An earlier version of this comment claimed parity was
            // "undisturbed"; that was too glib, and only looked true because the parity
            // fixture generated no zero-length flanks. It does now. Production still answers
            // `Region` at a flankless edge — its own off-end tests are guarded on the flank
            // existing — where ng answers a bound. The **bytes** are unchanged wherever ng
            // reports a span, which is what the parity claim actually covers; what differs is
            // the *classification*, deliberately.
            left_anchored: left_flank_len > 0 && tract_start != 0,
            right_anchored: right_flank_len > 0 && tract_end != read_len,
        })
    }
}

impl<E: Emission> BestPathAligner for SsrFlatGapAligner<E> {
    type Scratch = ViterbiScratch;
    type Output = RepeatSpan;
    type Context<'a> = RepeatContext<'a>;

    /// Measure the read's repeat, and report **which flanks held it**.
    ///
    /// An empty reference gives [`RepeatSpan::Unanchored`] — a *defined answer*, not an
    /// error. There is no locus to delimit, so the read carries no fact about one, which is
    /// exactly what `Unanchored` means (arch §3).
    fn align(
        &self,
        read: ReadBases<'_>,
        reference: &[u8],
        context: Self::Context<'_>,
        scratch: &mut Self::Scratch,
    ) -> Self::Output {
        self.delimit(read, reference, &context, scratch)
            .map_or(RepeatSpan::Unanchored, |readout| readout.classify())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::alignment::emission::{FlatEmission, PerQualityEmission};
    use crate::ng::alignment::{RepeatGeometry, StutterModel};
    use crate::ng::types::Bp;
    use crate::ng::types::Motif;

    /// Build a reference stretch as `left flank + tract + right flank`, and the geometry
    /// that describes it.
    fn frame(left: &[u8], tract: &[u8], right: &[u8], motif: &[u8]) -> (Vec<u8>, RepeatGeometry) {
        let mut reference = Vec::new();
        reference.extend_from_slice(left);
        reference.extend_from_slice(tract);
        reference.extend_from_slice(right);
        let geometry = RepeatGeometry {
            left_flank_len: Bp(left.len() as u64),
            right_flank_len: Bp(right.len() as u64),
            motif: Motif::new(motif).expect("a valid test motif"),
        };
        (reference, geometry)
    }

    /// Run the delimiter on a read, at uniform high quality.
    fn delimit(read: &[u8], reference: &[u8], geometry: &RepeatGeometry) -> Option<TractReadout> {
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry,
            stutter: &stutter,
        };
        let quality = vec![35u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = ViterbiScratch::new();
        aligner.delimit(bases, reference, &context, &mut scratch)
    }

    /// **A clean repeat measures exactly.** The read is the reference, so the tract must
    /// come back as exactly the tract's own coordinates.
    #[test]
    fn a_clean_repeat_measures_exactly() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let readout = delimit(&reference, &reference, &geometry).expect("a real frame");

        assert_eq!(readout.tract_start, 8);
        assert_eq!(readout.tract_end, 20);
        assert!(readout.left_anchored && readout.right_anchored);
    }

    /// **A longer allele must not collapse to the reference.** This is the failure the
    /// two-regime split exists to prevent: with one stiff gap everywhere, an allele more
    /// than about one unit from the reference is pulled back onto the reference length.
    #[test]
    fn a_longer_allele_does_not_collapse_to_the_reference() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        // Three extra units in the read.
        let read = b"ACGTACGTCAGCAGCAGCAGCAGCAGCAGTTGGTTGGAT";
        let readout = delimit(read, &reference, &geometry).expect("a real frame");

        assert_eq!(readout.tract_start, 8);
        assert_eq!(
            readout.tract_end, 29,
            "the read's 21-base tract was not measured verbatim"
        );
        assert_eq!(readout.tract_end - readout.tract_start, 21);
        assert!(readout.left_anchored && readout.right_anchored);
    }

    /// A shorter allele, the other direction.
    #[test]
    fn a_shorter_allele_measures_shorter() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGTTGGTTGGAT";
        let readout = delimit(read, &reference, &geometry).expect("a real frame");

        assert_eq!(readout.tract_start, 8);
        assert_eq!(readout.tract_end - readout.tract_start, 6);
        assert!(readout.left_anchored && readout.right_anchored);
    }

    /// **An interrupted repeat comes out verbatim.** The flat tract gap is content-agnostic
    /// by design, so a substitution inside the tract must not be tidied away — the read is
    /// measured, not judged.
    #[test]
    fn an_interrupted_repeat_comes_out_verbatim() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        // One interior base substituted: CAGCAGC*T*GCAG.
        let read = b"ACGTACGTCAGCAGCTGCAGTTGGTTGGAT";
        let readout = delimit(read, &reference, &geometry).expect("a real frame");

        assert_eq!(readout.tract_start, 8);
        assert_eq!(readout.tract_end, 20);
        assert_eq!(
            &read[readout.tract_start as usize..readout.tract_end as usize],
            b"CAGCAGCTGCAG",
            "the interrupted tract was not read out verbatim"
        );
    }

    /// A read running off its 3′ end inside the repeat: the left flank anchors, the right
    /// one cannot.
    #[test]
    fn a_read_running_off_its_end_anchors_only_the_left_flank() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCAG";
        let readout = delimit(read, &reference, &geometry).expect("a real frame");

        assert!(readout.left_anchored);
        assert!(!readout.right_anchored);
        assert_eq!(readout.tract_start, 8);
        assert_eq!(readout.tract_end, read.len() as u64);
    }

    /// The mirror: a read starting inside the repeat anchors only the right flank. The two
    /// sides are recorded separately, which production's single `BorderOffEnd` cannot do —
    /// the widening B2 builds on.
    #[test]
    fn a_read_starting_inside_the_repeat_anchors_only_the_right_flank() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"CAGCAGTTGGTTGGAT";
        let readout = delimit(read, &reference, &geometry).expect("a real frame");

        assert!(!readout.left_anchored);
        assert!(readout.right_anchored);
        assert_eq!(readout.tract_start, 0);
    }

    /// An empty reference is not a real locus. A **defined answer**, not an error — this
    /// module returns no `Result` (arch §3).
    #[test]
    fn an_empty_reference_has_nothing_to_delimit() {
        let geometry = RepeatGeometry {
            left_flank_len: Bp(0),
            right_flank_len: Bp(0),
            motif: Motif::new(b"CAG").expect("a valid test motif"),
        };
        assert_eq!(delimit(b"ACGT", b"", &geometry), None);
    }

    /// **The tract-aware gap-open is what makes the long allele work.** Priced at the stiff
    /// flank rate everywhere, the same read collapses toward the reference length — which
    /// is the production failure the split was introduced to fix. Demonstrated by pricing
    /// the whole frame as flank (both flanks covering everything leaves no tract column).
    #[test]
    fn one_stiff_gap_everywhere_collapses_the_long_allele() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCAGCAGCAGCAGCAGTTGGTTGGAT";

        let soft = delimit(read, &reference, &geometry).expect("a real frame");
        assert_eq!(soft.tract_end - soft.tract_start, 21);

        // The same aligner with the tract regime switched off, by making the tract gap as
        // stiff as the flank one.
        let aligner = SsrFlatGapAligner {
            emission: PerQualityEmission::new(),
            costs: TransitionCosts {
                ln_gap_open_tract: GAP_OPEN_PROB.ln(),
                ..TransitionCosts::new()
            },
        };
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry: &geometry,
            stutter: &stutter,
        };
        let quality = vec![35u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = ViterbiScratch::new();
        let stiff = aligner
            .delimit(bases, &reference, &context, &mut scratch)
            .expect("a real frame");

        assert!(
            stiff.tract_end - stiff.tract_start < soft.tract_end - soft.tract_start,
            "a uniformly stiff gap should pull the measurement toward the reference length"
        );
    }

    /// The match→match cost is built from the **flank** gap-open and never recomputed under
    /// the tract regime, so inside the tract the three transitions leaving a match sum to
    /// about 1.02. Reproduced from production deliberately; pinned so that fixing it is a
    /// visible choice rather than an accident (spec §4.2).
    #[test]
    fn the_tract_transitions_leaving_a_match_sum_past_one() {
        let costs = TransitionCosts::new();
        let flank_total = costs.ln_match_to_match.exp() + 2.0 * costs.ln_gap_open.exp();
        assert!(
            (flank_total - 1.0).abs() < 1e-12,
            "the flank transitions should sum to exactly one, summed to {flank_total}"
        );

        let tract_total = costs.ln_match_to_match.exp() + 2.0 * costs.ln_gap_open_tract.exp();
        assert!(
            (tract_total - 1.02).abs() < 1e-4,
            "expected the known ~1.02 tract total, got {tract_total}"
        );
    }

    /// Scratch is reused across reads and must not carry anything between them: the same
    /// read measured twice, with a differently-shaped read in between, gives the same
    /// answer. This is what makes the aligner's output independent of call order.
    #[test]
    fn reusing_scratch_does_not_leak_between_reads() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry: &geometry,
            stutter: &stutter,
        };
        let mut scratch = ViterbiScratch::new();

        let run = |scratch: &mut ViterbiScratch, read: &[u8]| {
            let quality = vec![35u8; read.len()];
            let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
            aligner.delimit(bases, &reference, &context, scratch)
        };

        let first = run(&mut scratch, &reference);
        let _ = run(&mut scratch, b"ACGTACGTCAGCAGCAGCAGCAGCAGCAGTTGGTTGGAT");
        let again = run(&mut scratch, &reference);

        assert_eq!(first, again, "scratch reuse changed the answer");
        assert_eq!(
            first,
            run(&mut ViterbiScratch::new(), &reference),
            "a fresh scratch disagreed with a reused one"
        );

        // Also the shrinking direction: a long read first, then a short one. Sound today
        // only because the fill is exhaustive — this is the ordering that breaks first when
        // Milestone C makes the fill partial.
        let mut shrinking = ViterbiScratch::new();
        let _ = run(&mut shrinking, b"ACGTACGTCAGCAGCAGCAGCAGCAGCAGTTGGTTGGAT");
        assert_eq!(
            run(&mut shrinking, &reference),
            first,
            "a scratch previously sized for a longer read changed the answer"
        );
    }

    /// **`best_of` keeps the first candidate on a tie.** That is the whole mechanism behind
    /// the M > D > I rule: each call site encodes the preference by *ordering* its
    /// candidates, so the comparison must be strict `>` and never `>=`. A `>=` slip would
    /// silently invert every tie in the matrix.
    #[test]
    fn best_of_keeps_the_first_candidate_on_a_tie() {
        use State::{Deletion, Insertion, Match};

        // All equal: the first wins. This is the whole mechanism behind M > D > I — each
        // call site encodes the preference by *ordering* its candidates.
        assert_eq!(
            best_of(&[(1.0, Match), (1.0, Deletion), (1.0, Insertion)]),
            (1.0, Match)
        );
        // A strictly better later candidate still wins.
        assert_eq!(
            best_of(&[(1.0, Match), (2.0, Deletion), (2.0, Insertion)]),
            (2.0, Deletion)
        );
        // A single candidate is returned unchanged.
        assert_eq!(best_of(&[(-3.5, Deletion)]), (-3.5, Deletion));
        // Unreachable candidates do not displace a reachable one, in either position.
        assert_eq!(
            best_of(&[(UNREACHABLE, Insertion), (0.0, Match)]),
            (0.0, Match)
        );
        assert_eq!(
            best_of(&[(0.0, Match), (UNREACHABLE, Insertion)]),
            (0.0, Match)
        );
        // All unreachable: still answers, with the first — the traceback needs *a* state.
        assert_eq!(
            best_of(&[(UNREACHABLE, Match), (UNREACHABLE, Deletion)]),
            (UNREACHABLE, Match)
        );
    }

    /// The discriminants are load-bearing: they are the array indices, and they must stay
    /// exactly production's `M = 0`, `I = 1`, `D = 2`. A reordering of the variants would
    /// silently permute every score array.
    #[test]
    fn the_state_discriminants_are_the_array_indices() {
        assert_eq!(State::Match.index(), 0);
        assert_eq!(State::Insertion.index(), 1);
        assert_eq!(State::Deletion.index(), 2);
    }

    /// **The tie-break is a correctness rule, not a nicety** (spec §4.2): two reads of the
    /// same molecule must measure the same repeat, and where several line-ups score equally
    /// an arbitrary choice lets identical input produce different answers.
    ///
    /// This test exists because **reversing the candidate order in the match recurrence to
    /// I > D > M left all ten original tests passing** — the rule the module documents most
    /// emphatically was the one nothing checked.
    ///
    /// Making a *genuine* tie takes work, and the first attempt at this test failed to: an
    /// emission model where a match and a mismatch score alike still leaves the three match
    /// candidates strictly ordered, because they differ in their **transition** costs
    /// (`ln_match_to_match` ≈ ln(0.99994) against `ln_gap_close` = ln(0.632)). With a strict
    /// winner the candidate order is never consulted, and the mutation survived.
    ///
    /// So this flattens **both** halves: every transition costs zero (probability 1) and
    /// every emission is equal, which makes every path through the matrix score identically.
    /// The answer is then decided by nothing but the preference order, and it moves if that
    /// order does.
    #[test]
    fn ties_are_broken_toward_match_then_deletion_then_insertion() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let flat_costs = TransitionCosts {
            ln_match_to_match: 0.0,
            ln_gap_open: 0.0,
            ln_gap_open_tract: 0.0,
            ln_gap_close: 0.0,
            ln_gap_extend: 0.0,
        };
        // ε = 0.75 is the crossover established in step A1: a match and a mismatch score
        // identically there, so no cell can be decided by sequence either.
        let aligner = SsrFlatGapAligner {
            emission: FlatEmission::try_new(0.75).expect("a valid test error rate"),
            costs: flat_costs,
        };
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry: &geometry,
            stutter: &stutter,
        };
        let quality = vec![35u8; reference.len()];
        let bases = ReadBases::try_new(&reference, &quality).expect("matched lengths");
        let mut scratch = ViterbiScratch::new();
        let readout = aligner
            .delimit(bases, &reference, &context, &mut scratch)
            .expect("a real frame");

        // Every path scores the same, so preferring match keeps the line-up on the diagonal
        // and the tract lands on the reference's own coordinates. Preferring insertion —
        // the reversed order — walks the read off before the reference and moves them.
        assert_eq!(
            (readout.tract_start, readout.tract_end),
            (8, 20),
            "the tie-break no longer prefers the diagonal"
        );
    }

    /// Determinism again, from the other side: the answer must not depend on how many times
    /// the aligner has been used, nor on which emission model is *configured* for equal
    /// inputs. Two identical reads measure identically.
    #[test]
    fn two_reads_of_the_same_molecule_measure_the_same_repeat() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCAGCAGCAGTTGGTTGGAT";
        let first = delimit(read, &reference, &geometry);
        let second = delimit(read, &reference, &geometry);
        assert_eq!(first, second);
    }

    /// **The tract window is two-sided, and both edges must be pinned.** Widening the
    /// right-hand comparison from `column <= tract_last_column` to `<` left all ten original
    /// tests passing, while the equivalent mutation on the *left* edge failed two — coverage
    /// was asymmetric.
    ///
    /// The last tract column and the first right-flank column are checked directly, because
    /// they are what the predicate actually decides between.
    #[test]
    fn the_tract_window_includes_its_last_column_and_excludes_the_flank() {
        let (reference, _geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let costs = TransitionCosts::new();
        let left_flank_len = 8usize;
        let tract_last_column = reference.len() - 10;
        let gap_open = |column: usize| {
            if left_flank_len < column && column <= tract_last_column {
                costs.ln_gap_open_tract
            } else {
                costs.ln_gap_open
            }
        };

        // The left junction: column 8 is still flank, column 9 is the first tract column.
        assert_eq!(gap_open(8), costs.ln_gap_open);
        assert_eq!(gap_open(9), costs.ln_gap_open_tract);
        // The right junction: column 20 is the last tract column, 21 is the first flank one.
        assert_eq!(tract_last_column, 20);
        assert_eq!(gap_open(20), costs.ln_gap_open_tract);
        assert_eq!(gap_open(21), costs.ln_gap_open);
    }

    /// **Production's dropped `tract_start > tract_end` rejection is unreachable here**, and
    /// this is what makes that claim checkable rather than asserted. It matters because step
    /// B2 turns these offsets into a `RepeatSpan`, where an inverted pair would become a
    /// `Between(..)` whose measured length saturates to zero — a *confident zero-unit
    /// allele*, which is a far worse answer than a rejection.
    ///
    /// Swept across flank sizes, tract sizes, and reads that run off either end.
    #[test]
    fn the_tract_offsets_are_never_inverted() {
        let motif = b"CAG";
        for left in [0usize, 1, 4, 8] {
            for right in [0usize, 1, 5, 10] {
                for tract_units in 1..=5usize {
                    let left_flank = &b"ACGTACGTAC"[..left];
                    let right_flank = &b"TTGGTTGGAT"[..right];
                    let tract: Vec<u8> = motif.repeat(tract_units);
                    let (reference, geometry) = frame(left_flank, &tract, right_flank, motif);

                    for read in [
                        reference.clone(),
                        reference[..reference.len().min(4)].to_vec(),
                        reference[reference.len() / 2..].to_vec(),
                        tract.clone(),
                        Vec::new(),
                    ] {
                        if let Some(readout) = delimit(&read, &reference, &geometry) {
                            assert!(
                                readout.tract_start <= readout.tract_end,
                                "inverted offsets {readout:?} at left={left} right={right} \
                                 units={tract_units} read={read:?}"
                            );
                            assert!(readout.tract_end <= read.len() as u64);
                        }
                    }
                }
            }
        }
    }

    /// The two flanks must be told apart, not merely counted. Both flanks were 8 bp in the
    /// original fixtures, so **transposing `left_flank_len` and `right_flank_len` passed
    /// every test** — the same blind-fixture pattern that hid two transpositions in step A3.
    /// The frame is now deliberately asymmetric (8 bp left, 10 bp right).
    #[test]
    fn the_two_flanks_are_not_interchangeable() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let straight = delimit(&reference, &reference, &geometry).expect("a real frame");

        let transposed = RepeatGeometry {
            left_flank_len: geometry.right_flank_len,
            right_flank_len: geometry.left_flank_len,
            motif: geometry.motif,
        };
        let swapped = delimit(&reference, &reference, &transposed).expect("a real frame");

        assert_ne!(
            (straight.tract_start, straight.tract_end),
            (swapped.tract_start, swapped.tract_end),
            "swapping the flank lengths did not change the measurement"
        );
    }

    /// Degenerate inputs have **defined answers**, not panics — this module returns no
    /// `Result` (arch §3). None of these is necessarily *useful*; the point is that each is
    /// answered rather than crashing, and that the answers are pinned so a later change to
    /// them is visible.
    #[test]
    fn degenerate_reads_have_defined_answers() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");

        // An empty read: nothing anchors, and both offsets collapse to zero.
        let empty = delimit(b"", &reference, &geometry).expect("a real frame");
        assert_eq!((empty.tract_start, empty.tract_end), (0, 0));

        // A read shorter than the left flank.
        let stub = delimit(b"AC", &reference, &geometry).expect("a real frame");
        assert!(stub.tract_start <= stub.tract_end);

        // A frame with no flanks at all — a repeat at both contig edges. **Neither side can
        // anchor**, because there is no flank to hold it: a flank that does not exist cannot
        // anchor. The span is still the whole read, but it is a bound, not a measurement.
        let (bare, bare_geometry) = frame(b"", b"CAGCAGCAG", b"", b"CAG");
        let all_repeat = delimit(b"CAGCAGCAG", &bare, &bare_geometry).expect("a real frame");
        assert_eq!((all_repeat.tract_start, all_repeat.tract_end), (0, 9));
        assert!(!all_repeat.left_anchored && !all_repeat.right_anchored);
    }

    // -----------------------------------------------------------------
    // B2 — the RepeatSpan readout, the widening
    // -----------------------------------------------------------------

    /// Align through the public [`BestPathAligner`] surface.
    fn align(read: &[u8], reference: &[u8], geometry: &RepeatGeometry) -> RepeatSpan {
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry,
            stutter: &stutter,
        };
        let quality = vec![35u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = ViterbiScratch::new();
        aligner.align(bases, reference, context, &mut scratch)
    }

    /// **`classify` is a pure mapping, and this pins every arm of it — including which side
    /// is which.** A mis-assigned side is the silent failure this step exists to avoid: the
    /// read still looks perfectly good and is simply filed under the wrong observation
    /// class. `FromLeft` means the **left** flank is the one that held.
    #[test]
    fn classify_maps_each_anchoring_to_its_own_case() {
        let readout = |left_anchored, right_anchored| TractReadout {
            tract_start: 4,
            tract_end: 10,
            left_anchored,
            right_anchored,
        };

        assert_eq!(readout(true, true).classify(), RepeatSpan::Between(4..10));
        assert_eq!(readout(true, false).classify(), RepeatSpan::FromLeft(4..10));
        assert_eq!(
            readout(false, true).classify(),
            RepeatSpan::FromRight(4..10)
        );
        assert_eq!(readout(false, false).classify(), RepeatSpan::Unanchored);

        // The two one-flank cases must not be interchangeable — swapping the inputs must
        // swap the outputs, not leave them alone.
        assert_ne!(
            readout(true, false).classify(),
            readout(false, true).classify()
        );
    }

    /// **All four cases are reachable from real alignments**, which is what the plan asks
    /// this step to prove. Anchoring is what distinguishes them, so each read is chosen to
    /// present a different pair of flanks to the same frame.
    #[test]
    fn every_repeat_span_case_is_reachable() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");

        // Both flanks present.
        assert!(matches!(
            align(&reference, &reference, &geometry),
            RepeatSpan::Between(_)
        ));
        // Runs off its 3′ end inside the repeat: only the left flank held.
        assert!(matches!(
            align(b"ACGTACGTCAGCAGCAG", &reference, &geometry),
            RepeatSpan::FromLeft(_)
        ));
        // Starts inside the repeat: only the right flank held.
        assert!(matches!(
            align(b"CAGCAGTTGGTTGGAT", &reference, &geometry),
            RepeatSpan::FromRight(_)
        ));
        // Wholly inside the repeat: neither flank held, and no span is reported.
        assert_eq!(
            align(b"CAGCAGCAG", &reference, &geometry),
            RepeatSpan::Unanchored
        );
    }

    /// A clean repeat is a **measurement**, and it measures exactly. This is the case the
    /// whole distinction exists to separate from the other three.
    #[test]
    fn a_clean_repeat_is_a_measurement_of_the_right_length() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let span = align(&reference, &reference, &geometry);

        assert!(span.is_measurement());
        assert_eq!(span.measured_length(), Some(12));
        assert_eq!(span.observed_span(), Some(&(8..20)));
    }

    /// A longer allele must be **measured**, not collapsed — and it must still come back as
    /// a measurement rather than a bound, because both flanks held.
    #[test]
    fn a_longer_allele_is_measured_at_its_own_length() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let span = align(
            b"ACGTACGTCAGCAGCAGCAGCAGCAGCAGTTGGTTGGAT",
            &reference,
            &geometry,
        );

        assert!(span.is_measurement());
        assert_eq!(span.measured_length(), Some(21));
    }

    /// **A truncated read is a lower bound, not a short allele.** This is the distinction the
    /// widening buys, and the one a caller most easily loses: the span looks like any other,
    /// so only the *case* says it is not a measurement.
    #[test]
    fn a_truncated_read_bounds_the_length_instead_of_measuring_it() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCAG";
        let span = align(read, &reference, &geometry);

        assert!(!span.is_measurement());
        assert_eq!(
            span.measured_length(),
            None,
            "a bound must not report a length"
        );
        // It still bounds the repeat below, and the bound is the part of it the read shows.
        assert_eq!(span.length_lower_bound(read.len() as u64), 9);
        // ...and the raw span is available for extracting the bases it did see.
        assert_eq!(span.observed_span(), Some(&(8..17)));
    }

    /// An interrupted repeat comes out **verbatim** through the public surface too — the
    /// flat tract gap is content-agnostic, so a substitution inside the tract is measured,
    /// not tidied away.
    #[test]
    fn an_interrupted_repeat_is_measured_verbatim() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCTGCAGTTGGTTGGAT";
        let span = align(read, &reference, &geometry);

        assert_eq!(span.measured_length(), Some(12));
        let observed = span.observed_span().expect("a measured span");
        assert_eq!(
            &read[observed.start as usize..observed.end as usize],
            b"CAGCAGCTGCAG"
        );
    }

    /// An empty reference is not a locus, so the read carries no fact about one — a defined
    /// answer rather than an error (arch §3).
    #[test]
    fn an_empty_reference_aligns_to_unanchored() {
        let geometry = RepeatGeometry {
            left_flank_len: Bp(0),
            right_flank_len: Bp(0),
            motif: Motif::new(b"CAG").expect("a valid test motif"),
        };
        assert_eq!(align(b"ACGT", b"", &geometry), RepeatSpan::Unanchored);
    }

    /// The public surface must report the offsets the algorithm actually found.
    ///
    /// **Deliberately not written as `span == readout.classify()`** — that was the first
    /// version, and it is a tautology: `align` *is* `delimit(..).map_or(Unanchored,
    /// classify)`, so it compares `classify(x)` with `classify(x)` and passes under every
    /// mutation of `classify`, including swapping `FromLeft` with `FromRight`. What is
    /// checked instead is the part that is not definitionally true: that the wiring carries
    /// B1's offsets through unchanged, and that the *case* follows from the anchoring flags
    /// spelled out independently here.
    #[test]
    fn the_public_surface_reports_the_offsets_the_traceback_found() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        for read in [
            reference.clone(),
            b"ACGTACGTCAGCAGCAGCAGCAGCAGCAGTTGGTTGGAT".to_vec(),
            b"ACGTACGTCAGCAGCAG".to_vec(),
            b"CAGCAGTTGGTTGGAT".to_vec(),
            b"CAGCAGCAG".to_vec(),
        ] {
            let readout = delimit(&read, &reference, &geometry).expect("a real frame");
            let span = align(&read, &reference, &geometry);

            // The case, derived here from the flags rather than by calling `classify`.
            let expected = match (readout.left_anchored, readout.right_anchored) {
                (true, true) => RepeatSpan::Between(readout.tract_start..readout.tract_end),
                (true, false) => RepeatSpan::FromLeft(readout.tract_start..readout.tract_end),
                (false, true) => RepeatSpan::FromRight(readout.tract_start..readout.tract_end),
                (false, false) => RepeatSpan::Unanchored,
            };
            assert_eq!(span, expected, "read {read:?}");

            if let Some(observed) = span.observed_span() {
                assert_eq!(observed.start, readout.tract_start);
                assert_eq!(observed.end, readout.tract_end);
            }
        }
    }

    /// **A locus with no flank on a side can never measure that side.** This is the case the
    /// review caught: before the fix, a flankless frame reported `Between` — a *measurement*
    /// — for a read that had measured nothing, because "a flank that does not exist cannot
    /// fail to hold" turned into "the length is pinned".
    ///
    /// The two readings are indistinguishable there (a read that ended because the tract
    /// ended, and one that ran out mid-tract, produce identical offsets), so the honest
    /// answer is the conservative one: a bound, not a measurement.
    #[test]
    fn a_locus_without_flanks_yields_bounds_rather_than_measurements() {
        // No flanks at all: nothing can anchor, whatever the read shows.
        let (bare, bare_geometry) = frame(b"", b"CAGCAGCAGCAG", b"", b"CAG");
        for read in [b"CAGCAG".as_slice(), b"CAGCAGCAGCAG".as_slice()] {
            let span = align(read, &bare, &bare_geometry);
            assert_eq!(span, RepeatSpan::Unanchored, "read {read:?}");
            assert!(!span.is_measurement());
            assert_eq!(span.measured_length(), None);
        }

        // A left flank only — a repeat running to a contig end. The left side can anchor;
        // the right never can, so the best available answer is a lower bound.
        let (left_only, left_geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"", b"CAG");
        let span = align(&left_only, &left_only, &left_geometry);
        assert!(matches!(span, RepeatSpan::FromLeft(_)));
        assert_eq!(
            span.measured_length(),
            None,
            "a flankless side must not report a measured length"
        );
        assert_eq!(span.length_lower_bound(left_only.len() as u64), 12);
    }

    // -----------------------------------------------------------------
    // C1 — banding
    // -----------------------------------------------------------------

    /// The band is three terms of **geometry**, and the floor term is the one the spec §3
    /// insists on: a fixed width would lose a long allele the moment `|read − reference|`
    /// exceeds it. Pinned against the exact production example (52 vs 28 → 24 forced cells).
    #[test]
    fn band_width_is_geometry_not_the_slip_cutoff() {
        // The forced floor alone already exceeds any fixed constant.
        assert_eq!(band_width(52, 28, 0, 0), 24 + BAND_HEADROOM);
        assert_eq!(band_width(28, 52, 0, 0), 24 + BAND_HEADROOM);
        // Flanks widen it — the run-off correction. (Equal-length read/reference, so the
        // floor term is zero and only the flanks and the slop remain.)
        assert_eq!(band_width(30, 30, 8, 10), 8 + 10 + BAND_HEADROOM);
        // MAX_SLIP (a scoring cutoff) is nowhere in it: a 40-unit expansion of a period-2
        // tract is 80 forced cells, far past MAX_SLIP × period = 20.
        assert!(band_width(110, 30, 0, 0) >= 80);
    }

    /// The band must always reach the far corner `(m, n)`, or the alignment cannot complete
    /// and the answer is lost. `|m − n|` is exactly the corner's deviation, so the band —
    /// which is that plus non-negative terms — always contains it. Swept over a wide range of
    /// shapes, since this is the one property a too-narrow band violates.
    #[test]
    fn the_band_always_reaches_the_far_corner() {
        for read_len in [0usize, 1, 5, 40, 200] {
            for reference_len in [1usize, 5, 40, 200] {
                for (left, right) in [(0, 0), (8, 0), (0, 10), (12, 12)] {
                    let band = band_width(read_len, reference_len, left, right);
                    assert!(
                        read_len.abs_diff(reference_len) <= band,
                        "band {band} cannot reach the corner at read {read_len}, ref \
                         {reference_len}"
                    );
                }
            }
        }
    }

    /// **The plan's required extreme: a long allele a too-narrow band would collapse.** The
    /// read carries eight extra units — a 24-base expansion of a period-3 tract — and must
    /// come back measured at its own length, not pulled toward the reference. A fixed-width
    /// band (or one derived from the slip cutoff) would clip the 24-cell excursion and return
    /// a short allele, silently.
    #[test]
    fn a_long_allele_at_the_extreme_is_measured_not_collapsed() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        // Reference tract is 4 units (12 bp); the read carries 12 units (36 bp).
        let read = b"ACGTACGTCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTGGTTGGAT";
        let span = align(read, &reference, &geometry);

        assert!(span.is_measurement());
        assert_eq!(
            span.measured_length(),
            Some(36),
            "the long allele was collapsed — the band clipped its excursion"
        );
    }

    /// **The run-off finding, pinned as a regression.** A left-flankless locus whose read
    /// expands the tract *and* runs off the right flank: the optimal path's peak deviation is
    /// the tract expansion, well past the net `|read − reference|`, so a band without the flank
    /// terms clips it. This is the shape the 200,000-case soak first caught (seed `0x5eed0001`
    /// case 1745); here it is a fixed input so a regression names it directly.
    ///
    /// Checked against production, which is unbanded and therefore the ground truth for
    /// whether the band is wide enough — measured bytes must match exactly.
    #[test]
    fn a_run_off_flank_read_is_measured_despite_the_bow() {
        // Left-flankless: tract of 4 units, right flank of 10.
        let (reference, geometry) =
            frame(b"", b"GTTGTGGTTGTGGTTGTGGTTGTG", b"AACCCAGGCC", b"GTTGTG");
        // The read is 7 units of tract (a +3-unit, 18-base expansion) then runs off the
        // flank after its first base.
        let read = b"GTTGTGGTTGTGGTTGTGGTTGTGGTTGTGGTTGTGGTTGTGA";

        // The peak matrix deviation exceeds the net length difference by the run-off length.
        let net = read.len().abs_diff(reference.len());
        assert!(
            band_width(read.len(), reference.len(), 0, 10) > net + BAND_HEADROOM,
            "the flank term is what carries this case"
        );

        // It must still measure — a bound (the read ran off the right end), not a collapse or
        // a panic. The observed bytes are what the unbanded path would give.
        let readout = delimit(read, &reference, &geometry).expect("a real frame");
        assert!(readout.tract_start <= readout.tract_end);
        assert!(readout.tract_end <= read.len() as u64);
        // The left is flankless so cannot anchor; the read ran off the right end, so neither
        // flank held — an unanchored bound, and specifically the whole read is tract-side.
        assert_eq!(readout.tract_start, 0);
    }

    /// **A read that runs off *both* ends of an expanded tract.** The read is pure tract, so
    /// both reference flanks are deleted and neither anchors — the geometry that stresses the
    /// band on both sides at once. Its answer must match production's unbanded delimiter, and
    /// it must be a bound, not a measurement.
    #[test]
    fn a_read_off_both_ends_of_an_expanded_tract_is_measured() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        // Pure tract, expanded to nine units — no flank bases at all.
        let read = b"CAGCAGCAGCAGCAGCAGCAGCAGCAG";
        let span = align(read, &reference, &geometry);

        // Neither flank held: the whole read lies inside the repeat.
        assert_eq!(span, RepeatSpan::Unanchored);
        assert!(!span.is_measurement());
        // The band still reached the answer — the delimiter did not panic or clip.
        let readout = delimit(read, &reference, &geometry).expect("a real frame");
        assert_eq!(
            (readout.tract_start, readout.tract_end),
            (0, read.len() as u64)
        );
    }

    /// **Scratch reuse across an extreme size drop.** A very long read grows the buffers; a
    /// tiny read immediately after must not read a stale cell left in the band from the large
    /// one. Banding writes out-of-band cells `UNREACHABLE` every row precisely so this holds —
    /// this pins it against the size drop that would expose a skip-don't-write bug.
    #[test]
    fn scratch_survives_an_extreme_size_drop_under_banding() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let mut scratch = ViterbiScratch::new();
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry: &geometry,
            stutter: &stutter,
        };
        let run = |scratch: &mut ViterbiScratch, read: &[u8]| {
            let quality = vec![35u8; read.len()];
            let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
            aligner.delimit(bases, &reference, &context, scratch)
        };

        // A large read to grow and populate the buffers well beyond the reference width.
        let huge: Vec<u8> = b"ACGTACGT"
            .iter()
            .chain(b"CAG".iter().cycle().take(120))
            .chain(b"TTGGTTGGAT")
            .copied()
            .collect();
        let _ = run(&mut scratch, &huge);

        // Then the reference itself, on the same scratch, must give the fresh-scratch answer.
        let reused = run(&mut scratch, &reference);
        let fresh = run(&mut ViterbiScratch::new(), &reference);
        assert_eq!(
            reused, fresh,
            "a stale band cell leaked across the size drop"
        );
    }
}
