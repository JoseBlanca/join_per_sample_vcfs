//! Algorithm 4 — the repeat-aware best-path aligner with **two penalties**: a whole-unit
//! slip priced from the [`StutterModel`], everything else priced as an ordinary gap.
//!
//! Algorithm 3 ([`super::ssr_best_path_flat_gap`]) charges one flat per-base gap anywhere
//! inside the tract, so a three-unit expansion costs three units' worth of per-base gaps.
//! This aligner instead prices a **whole-unit** slip by the stutter model — direction-
//! asymmetric, geometric in the number of units — and keeps the flat gap only for the
//! *out-of-frame* changes that are not a whole number of units (spec §4.2). It is the
//! HipSTR-shaped model, and it is the more faithful of the two; which one measures better is
//! the Milestone-D comparison, and **neither is settled** (spec §10.3).
//!
//! **This algorithm has no parity oracle.** Algorithm 3 is a port of production's delimiter
//! and can be shown byte-identical to it; algorithm 4 is new (arch §5: "no production
//! counterpart"), so it is *measured*, not verified. Its tests establish that it computes
//! what it claims — the three properties the spec makes mandatory, plus the degenerate cases
//! — not that it matches an external truth. What confidence it has comes from the differential
//! against algorithm 3 (the byte-parity-validated sibling), which is exactly why algorithm 3
//! is kept rather than deleted (see [`super::ssr_best_path_flat_gap`]).
//!
//! **This is the recommended STR delimiter.** It measured better than algorithm 3 in the D2
//! comparison (`doc/devel/reports/research/ssr_delimiter_3v4_comparison_2026-07-23.md`) — more
//! accurate, and without the small reference-pull bias the flat-gap model carries. Whether to
//! formally *adopt* it turns on the genotype comparison, which lives with the genotyping
//! (spec §10.3) and has not run; the measurement evidence points here.
//!
//! # The recurrence, and the three decisions behind it
//!
//! Algorithm 3's M / I / D matrix, unchanged, **plus two whole-unit slip states** inside the
//! tract: `SlipIns` (the read carries an extra unit the reference does not) and `SlipDel`
//! (the reference carries a unit the read does not). A slip transition jumps `period` cells
//! at once, which is why the score scratch keeps a small ring of the last `period + 1` rows
//! rather than the two rolling rows algorithm 3 needs.
//!
//! 1. **The geometric decomposes onto affine slip transitions.** The stutter probability of
//!    an `n`-unit expansion is `in_up · in_geom · (1 − in_geom)^(n − 1)` (spec §5.2), which in
//!    log space is `[ln(in_up · in_geom)] + (n − 1)·[ln(1 − in_geom)]` — an affine *open* plus
//!    `n − 1` *extends*. So a run of unit slips is priced exactly like an affine gap, one
//!    open and then extends, and the direction split (`in_up` vs `in_down`) is the asymmetry
//!    the spec demands.
//!
//! 2. **The slip cost is relative to no-slip.** A best-path aligner *maximises*, and only
//!    relative scores matter in an argmax — so the `equal` mass (the stutter probability of
//!    Δ = 0) is folded into the baseline: no slip scores `0`, and an `n`-unit slip scores
//!    `ln(P(n)) − ln(equal)`. This is the only place the `equal` mass can live in a Viterbi,
//!    since Δ = 0 corresponds to *not taking* a slip transition and so has nothing to attach a
//!    cost to. The open cost therefore carries the `− ln(equal)` term; the extend does not.
//!
//! 3. **A slipped unit's bases are scored against the motif** (decided with the owner,
//!    2026-07-23): an inserted unit that really is a copy of the repeat is rewarded, and one
//!    whose bases are wrong pays a mismatch — which is production's behaviour (its `q_r`
//!    resizes the candidate by tiling the motif and scores composition, `hipstr.rs`), and it
//!    is what makes the arithmetic in-frame test catch a wrong-base insertion downstream as a
//!    base error rather than as stutter (spec §4.2). A deletion consumes reference bases and
//!    so scores nothing.
//!
//! **Placement multiplicity: best placement, no divide** (decided with the owner). Production
//! *sums* over the reachable placements of a slip and divides by their count (`reach_variants`,
//! `1/|placements|`) — a marginal operation. A best-path aligner cannot: it maximises, so it
//! picks the single best-scoring placement. Here that falls out for free — the Viterbi may
//! start a slip at *any* tract column, and its max over those columns is the best placement.
//! The gap between this max and production's sum is a real difference between a ruler and a
//! score, and it is one of the things the Milestone-D comparison exists to measure.
//!
//! **This aligner is, by design, a ruler that prices toward tidy in-frame lengths.** Making a
//! whole-unit slip cheap pulls the measurement toward a whole number of units — which is the
//! §4.2 risk (does it round an interruption away?), and precisely the hypothesis the
//! comparison tests. It is intended, not a defect.

use super::emission::Emission;
use super::stutter::StutterModel;
use super::{BestPathAligner, ReadBases, RepeatContext, RepeatSpan};

// Reuse algorithm 3's tract-aware per-base gap transitions and the enum's shape; the
// difference is entirely in the two extra slip states, which algorithm 3 does not have.
use super::ssr_best_path_flat_gap::TransitionCosts;

/// The five states a cell can be entered in. The first three are algorithm 3's; the last two
/// are the whole-unit slip states this algorithm adds.
///
/// `#[repr(u8)]` with explicit discriminants so the values are stable array indices, matching
/// the pattern algorithm 3 uses for its three states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum State {
    /// A read base and a reference base consumed together.
    Match = 0,
    /// One read base against no reference base — an *out-of-frame* insertion (the flat gap).
    Insertion = 1,
    /// One reference base against no read base — an *out-of-frame* deletion (the flat gap).
    Deletion = 2,
    /// A whole unit of read bases the reference does not have — an *in-frame* expansion,
    /// priced by the stutter model. Consumes `period` read bases at once.
    SlipInsertion = 3,
    /// A whole unit of reference bases the read does not have — an *in-frame* contraction.
    /// Consumes `period` reference bases at once.
    SlipDeletion = 4,
}

/// The number of states — the width of each cell's score array and backpointer array.
const STATES: usize = 5;

impl State {
    #[inline]
    const fn index(self) -> usize {
        self as usize
    }
}

/// An unreachable score. A real value: some cells cannot be entered in some states.
const UNREACHABLE: f64 = f64::NEG_INFINITY;

/// The whole-unit slip transition costs, in log space, derived **from the shared
/// [`StutterModel`]** — not a second copy of its parameters (spec §4.2 requires the model be
/// written once and shared).
///
/// See the module docs for the derivation: `open` carries the affine open plus the
/// `− ln(equal)` baseline shift; `extend` is the geometric's per-extra-unit factor.
#[derive(Debug, Clone, Copy)]
struct SlipCosts {
    /// Log cost to open a run of **expansion** units (the read gains units), relative to
    /// no slip: `ln(in_up · in_geom) − ln(equal)`.
    open_expansion: f64,
    /// Log cost to open a run of **contraction** units (the read loses units), relative to
    /// no slip: `ln(in_down · in_geom) − ln(equal)`.
    open_contraction: f64,
    /// Log cost of each unit **after** the first, either direction: `ln(1 − in_geom)`.
    extend: f64,
}

impl SlipCosts {
    /// Derive the slip costs from the stutter model. Reads the model's parameters; keeps no
    /// copy of them.
    fn from_model(model: &StutterModel) -> Self {
        let ln_equal = model.equal().ln();
        Self {
            open_expansion: (model.in_up() * model.in_geom()).ln() - ln_equal,
            open_contraction: (model.in_down() * model.in_geom()).ln() - ln_equal,
            extend: (1.0 - model.in_geom()).ln(),
        }
    }
}

/// Per-worker scratch for algorithm 4.
///
/// Unlike algorithm 3's two rolling rows, a whole-unit insertion reaches cell `(i, j)` from
/// `(i − period, j)`, so the last `period + 1` score rows must be live at once — a small ring
/// (`period` is at most 6). The backpointer matrix is the full `(m + 1) × (n + 1)`, as in
/// algorithm 3; each entry is `[State; STATES]`, the winning predecessor per state.
///
/// Grow-and-keep; buffers only, deciding nothing that changes a result.
#[derive(Debug, Default)]
pub struct UnitSlipScratch {
    /// A ring of score rows, each `reference_len + 1` cells of `[f64; STATES]`. Indexed by
    /// `row_index % ring_len`.
    rows: Vec<Vec<[f64; STATES]>>,
    /// The winning predecessor state per cell and per state — flat, row-major, stride
    /// `reference_len + 1`.
    backpointers: Vec<[State; STATES]>,
}

impl UnitSlipScratch {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Size the ring to hold `ring_len` rows of `reference_len + 1` cells, and the full
    /// backpointer matrix. Grows only.
    fn resize(&mut self, read_len: usize, reference_len: usize, ring_len: usize) {
        let width = reference_len + 1;
        if self.rows.len() < ring_len {
            self.rows.resize_with(ring_len, Vec::new);
        }
        for row in &mut self.rows {
            if row.len() < width {
                row.resize(width, [UNREACHABLE; STATES]);
            }
        }
        let cells = (read_len + 1) * width;
        if self.backpointers.len() < cells {
            self.backpointers.resize(cells, [State::Match; STATES]);
        }
    }
}

/// Pick the best of the candidates, keeping the **first on ties** — the same tie-break rule
/// and encoding algorithm 3 uses (spec §4.2): candidates are passed in priority order, so the
/// caller encodes "match, then deletion, then insertion" (and, here, the slip states last) by
/// their order.
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

/// Algorithm 4: the two-penalty repeat delimiter. Generic over its [`Emission`], never behind
/// `dyn` (arch §4).
#[derive(Debug, Clone, Copy)]
pub struct SsrUnitSlipAligner<E: Emission> {
    emission: E,
    costs: TransitionCosts,
}

impl<E: Emission> SsrUnitSlipAligner<E> {
    #[must_use]
    pub fn new(emission: E) -> Self {
        Self {
            emission,
            costs: TransitionCosts::new(),
        }
    }

    /// Measure the read's repeat, pricing whole-unit slips from the stutter model.
    ///
    /// `None` when there is no reference frame — a defined answer, not an error (arch §3).
    #[must_use]
    pub fn delimit(
        &self,
        read: ReadBases<'_>,
        reference: &[u8],
        context: &RepeatContext<'_>,
        scratch: &mut UnitSlipScratch,
    ) -> Option<super::ssr_best_path_flat_gap::TractReadout> {
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

        let left_flank_len = usize::try_from(context.geometry.left_flank_len.get())
            .unwrap_or(usize::MAX)
            .min(reference_len);
        let right_flank_len = usize::try_from(context.geometry.right_flank_len.get())
            .unwrap_or(usize::MAX)
            .min(reference_len - left_flank_len);
        let tract_last = reference_len - right_flank_len;

        let motif = context.geometry.motif.as_bytes();
        let period = motif.len();
        debug_assert_eq!(period, context.geometry.motif.period());

        let slip = SlipCosts::from_model(context.stutter);
        let ring_len = period + 1;
        scratch.resize(read_len, reference_len, ring_len);
        let stride = reference_len + 1;
        let insertion_emission = self.emission.insert_ln();

        // Tract-aware gap-open for the *out-of-frame* per-base route (identical to algorithm 3).
        let gap_open = |column: usize| {
            if left_flank_len < column && column <= tract_last {
                self.costs.ln_gap_open_tract()
            } else {
                self.costs.ln_gap_open()
            }
        };
        // A reference column is inside the tract when a gap touching it is a tract gap.
        let column_in_tract = |column: usize| left_flank_len < column && column <= tract_last;
        // The `k`-th base of a whole unit inserted at tract column `column`: the motif,
        // phase-aligned so a unit inserted at a unit boundary begins at `motif[0]`.
        let motif_base = |column: usize, k: usize| motif[((column - left_flank_len) + k) % period];

        let UnitSlipScratch { rows, backpointers } = scratch;

        // Row 0 — no read base consumed. Reachable states: Match at the corner, and the two
        // deletion routes (out-of-frame `D` and whole-unit `SlipDeletion`) across the reference.
        {
            // Row 0 only reads and writes itself: a whole-unit deletion here (`i − period`
            // would be negative) chains from earlier columns of the same row, read as
            // `row0[column - period]`. So a plain `&mut rows[0]` suffices — no split needed.
            let row0 = &mut rows[0];
            row0[0] = [0.0, UNREACHABLE, UNREACHABLE, UNREACHABLE, UNREACHABLE];
            backpointers[0] = [State::Match; STATES];
            for column in 1..=reference_len {
                // Out-of-frame single-base deletion.
                let (d, d_pred) = best_of(&[
                    (
                        gap_open(column) + row0[column - 1][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend() + row0[column - 1][State::Deletion.index()],
                        State::Deletion,
                    ),
                ]);
                // Whole-unit deletion — only when the deleted unit lies wholly in the tract.
                let (sd, sd_pred) =
                    if column >= period && (column - period + 1..=column).all(column_in_tract) {
                        best_of(&[
                            (
                                slip.open_contraction + row0[column - period][State::Match.index()],
                                State::Match,
                            ),
                            (
                                slip.extend + row0[column - period][State::SlipDeletion.index()],
                                State::SlipDeletion,
                            ),
                        ])
                    } else {
                        (UNREACHABLE, State::Match)
                    };
                row0[column] = [UNREACHABLE, UNREACHABLE, d, UNREACHABLE, sd];
                backpointers[column] = [State::Match, State::Match, d_pred, State::Match, sd_pred];
            }
        }

        // Rows 1..=read_len.
        for row_index in 1..=read_len {
            let read_base = bases[row_index - 1];
            let scores = self.emission.scores_for(read.quality_at(row_index - 1));
            let back_row = row_index * stride;

            // `previous` is row i−1; `slip_src` is row i−period (for whole-unit insertions).
            // They live in distinct ring slots from the row being written.
            let prev_slot = (row_index - 1) % ring_len;
            let slip_slot = row_index.checked_sub(period).map(|r| r % ring_len);
            let cur_slot = row_index % ring_len;

            // Column 0 — read base before any reference base: out-of-frame insertion only.
            let (ins0, ins0_pred) = best_of(&[
                (
                    self.costs.ln_gap_open() + rows[prev_slot][0][State::Match.index()],
                    State::Match,
                ),
                (
                    self.costs.ln_gap_extend() + rows[prev_slot][0][State::Insertion.index()],
                    State::Insertion,
                ),
            ]);
            rows[cur_slot][0] = [
                UNREACHABLE,
                insertion_emission + ins0,
                UNREACHABLE,
                UNREACHABLE,
                UNREACHABLE,
            ];
            backpointers[back_row] = [
                State::Match,
                ins0_pred,
                State::Match,
                State::Match,
                State::Match,
            ];

            for column in 1..=reference_len {
                let emit = scores.pick(read_base, reference[column - 1]);

                // Match: from the diagonal, any state (a match may follow a completed slip).
                let (m, m_pred) = best_of(&[
                    (
                        self.costs.ln_match_to_match()
                            + rows[prev_slot][column - 1][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_close()
                            + rows[prev_slot][column - 1][State::Deletion.index()],
                        State::Deletion,
                    ),
                    (
                        self.costs.ln_gap_close()
                            + rows[prev_slot][column - 1][State::Insertion.index()],
                        State::Insertion,
                    ),
                    // A completed slip closes to a match at no extra cost — the geometric has
                    // already priced the whole run (module docs, decision 2).
                    (
                        rows[prev_slot][column - 1][State::SlipInsertion.index()],
                        State::SlipInsertion,
                    ),
                    (
                        rows[prev_slot][column - 1][State::SlipDeletion.index()],
                        State::SlipDeletion,
                    ),
                ]);

                let open = gap_open(column);

                // Out-of-frame single-base insertion (read base, no reference base).
                let (ins, ins_pred) = best_of(&[
                    (
                        open + rows[prev_slot][column][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend()
                            + rows[prev_slot][column][State::Insertion.index()],
                        State::Insertion,
                    ),
                ]);

                // Out-of-frame single-base deletion (reference base, no read base).
                let (del, del_pred) = best_of(&[
                    (
                        open + rows[cur_slot][column - 1][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend()
                            + rows[cur_slot][column - 1][State::Deletion.index()],
                        State::Deletion,
                    ),
                ]);

                // Whole-unit insertion: `period` read bases forming a unit, scored against the
                // motif, at a tract column. Reachable from row i−period, same column.
                //
                // **Each inserted base is scored at its own quality**, not the row's shared
                // `scores`: the unit spans read rows `row_index − period + 1 ..= row_index`,
                // and base `k` sits at read index `row_index − period + k` with quality
                // `quality_at` of that index. Using the current row's quality for all `period`
                // bases skews the slip score whenever intra-unit qualities vary.
                let (sins, sins_pred) = if let Some(slip_slot) = slip_slot {
                    if column_in_tract(column) {
                        let unit_emit: f64 = (0..period)
                            .map(|k| {
                                let idx = row_index - period + k;
                                self.emission
                                    .scores_for(read.quality_at(idx))
                                    .pick(bases[idx], motif_base(column, k))
                            })
                            .sum();
                        let (score, pred) = best_of(&[
                            (
                                slip.open_expansion + rows[slip_slot][column][State::Match.index()],
                                State::Match,
                            ),
                            (
                                slip.extend + rows[slip_slot][column][State::SlipInsertion.index()],
                                State::SlipInsertion,
                            ),
                        ]);
                        (unit_emit + score, pred)
                    } else {
                        (UNREACHABLE, State::Match)
                    }
                } else {
                    (UNREACHABLE, State::Match)
                };

                // Whole-unit deletion: `period` reference bases the read lacks, priced as a
                // contraction, when the deleted unit lies wholly in the tract. No emissions.
                let (sdel, sdel_pred) =
                    if column >= period && (column - period + 1..=column).all(column_in_tract) {
                        best_of(&[
                            (
                                slip.open_contraction
                                    + rows[cur_slot][column - period][State::Match.index()],
                                State::Match,
                            ),
                            (
                                slip.extend
                                    + rows[cur_slot][column - period][State::SlipDeletion.index()],
                                State::SlipDeletion,
                            ),
                        ])
                    } else {
                        (UNREACHABLE, State::Match)
                    };

                rows[cur_slot][column] = [emit + m, insertion_emission + ins, del, sins, sdel];
                backpointers[back_row + column] =
                    [m_pred, ins_pred, del_pred, sins_pred, sdel_pred];
            }
        }

        // Final cell (m, n): best terminal state. Tie-break order M > D > I > SlipIns > SlipDel.
        let last = &rows[read_len % ring_len][reference_len];
        let (_, final_state) = best_of(&[
            (last[State::Match.index()], State::Match),
            (last[State::Deletion.index()], State::Deletion),
            (last[State::Insertion.index()], State::Insertion),
            (last[State::SlipInsertion.index()], State::SlipInsertion),
            (last[State::SlipDeletion.index()], State::SlipDeletion),
        ]);

        let geometry = MatrixGeometry {
            stride,
            read_len,
            reference_len,
            left_flank_len,
            right_flank_len,
            period,
        };
        Some(trace_back(final_state, backpointers, geometry))
    }
}

/// The matrix shape the traceback needs — bundled so the walk takes one `Copy` argument
/// instead of six positional `usize`s, several of them interchangeable (a transposition of
/// `read_len`/`reference_len` or the two flank lengths would be a silent wrong measurement).
#[derive(Debug, Clone, Copy)]
struct MatrixGeometry {
    stride: usize,
    read_len: usize,
    reference_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
    period: usize,
}

/// Walk the traceback and read the tract off the two flank junctions, exactly as algorithm 3
/// does — the slip states simply consume `period` cells at a step instead of one.
fn trace_back(
    final_state: State,
    backpointers: &[[State; STATES]],
    geometry: MatrixGeometry,
) -> super::ssr_best_path_flat_gap::TractReadout {
    let MatrixGeometry {
        stride,
        read_len,
        reference_len,
        left_flank_len,
        right_flank_len,
        period,
    } = geometry;
    let left_junction = left_flank_len; // first tract reference base
    let right_junction = reference_len - right_flank_len; // first right-flank base
    let mut tract_start = 0usize;
    let mut tract_end = read_len;

    let mut i = read_len;
    let mut j = reference_len;
    let mut state = final_state;
    while i != 0 || j != 0 {
        let pred = backpointers[i * stride + j][state.index()];
        match state {
            State::Match => {
                let consumed = j - 1;
                if consumed == left_junction {
                    tract_start = i - 1;
                }
                if right_flank_len > 0 && consumed == right_junction {
                    tract_end = i - 1;
                }
                i -= 1;
                j -= 1;
            }
            State::Deletion => {
                let consumed = j - 1;
                if consumed == left_junction {
                    tract_start = i;
                }
                if right_flank_len > 0 && consumed == right_junction {
                    tract_end = i;
                }
                j -= 1;
            }
            State::Insertion => {
                i -= 1;
            }
            State::SlipInsertion => {
                // A whole unit of read bases inserted — no reference consumed.
                i -= period;
            }
            State::SlipDeletion => {
                // A whole unit of reference bases deleted — no read consumed. The deleted
                // 0-based reference indices are `j - period ..= j - 1`. A contraction at the
                // very start of the tract deletes the first tract base, which **is** the
                // left junction — so the crossing must be recorded, exactly as an ordinary
                // deletion records it (no read base consumed, so the offset is the current
                // read row). The right junction is the first right-flank base and lies
                // beyond any tract-interior deletion, so it cannot be consumed here; the
                // check is kept defensively.
                let deleted = (j - period)..j;
                if deleted.contains(&left_junction) {
                    tract_start = i;
                }
                if right_flank_len > 0 && deleted.contains(&right_junction) {
                    tract_end = i;
                }
                j -= period;
            }
        }
        state = pred;
    }

    let left_anchored = left_flank_len > 0 && tract_start != 0;
    let right_anchored = right_flank_len > 0 && tract_end != read_len;
    super::ssr_best_path_flat_gap::TractReadout {
        tract_start: tract_start as u64,
        tract_end: tract_end as u64,
        left_anchored,
        right_anchored,
    }
}

impl<E: Emission> BestPathAligner for SsrUnitSlipAligner<E> {
    type Scratch = UnitSlipScratch;
    type Output = RepeatSpan;
    type Context<'a> = RepeatContext<'a>;

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
    use crate::ng::alignment::emission::PerQualityEmission;
    use crate::ng::alignment::stutter::StutterRates;
    use crate::ng::alignment::{RepeatGeometry, StutterModel};
    use crate::ng::types::{Bp, Motif};

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

    /// Contraction-biased parameters — HipSTR's fitted values are, and the direction-asymmetry
    /// property is the whole point of algorithm 4.
    fn contraction_biased() -> StutterModel {
        StutterModel::new(StutterRates {
            in_up: 0.03,
            in_down: 0.07,
            in_geom: 0.9,
            out_up: 0.004,
            out_down: 0.012,
            out_geom: 0.8,
        })
    }

    fn measure(
        read: &[u8],
        reference: &[u8],
        geometry: &RepeatGeometry,
        model: &StutterModel,
    ) -> RepeatSpan {
        let aligner = SsrUnitSlipAligner::new(PerQualityEmission::new());
        let context = RepeatContext {
            geometry,
            stutter: model,
        };
        let quality = vec![35u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = UnitSlipScratch::new();
        aligner.align(bases, reference, context, &mut scratch)
    }

    /// **Mandatory property 1: a clean read of the reference measures the reference length.**
    /// The baseline sanity check — with no slip, the aligner must not invent one.
    #[test]
    fn a_clean_read_measures_the_reference_tract() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let span = measure(&reference, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(12));
    }

    /// A genuine whole-unit expansion (the read carries the motif verbatim) is measured at its
    /// own length, not collapsed to the reference — the slip route must be usable.
    #[test]
    fn a_whole_unit_expansion_is_measured() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        // Two extra units.
        let read = b"ACGTACGTCAGCAGCAGCAGCAGCAGTTGGTTGGAT";
        let span = measure(read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(18));
    }

    /// **Mandatory property 2: out-of-frame changes keep a route.** An interrupted repeat —
    /// a single-base substitution inside the tract, not a whole-unit event — must still be
    /// measured verbatim, not forced onto a tidy unit count.
    #[test]
    fn an_out_of_frame_change_still_has_a_route() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        // One interior base substituted — same length, out of frame in composition.
        let read = b"ACGTACGTCAGCAGCTGCAGTTGGTTGGAT";
        let span = measure(read, &reference, &geometry, &contraction_biased());
        // Same length as the reference tract, and the interrupted bases come out verbatim.
        assert_eq!(span.measured_length(), Some(12));
        let observed = span.observed_span().expect("a measured span");
        assert_eq!(
            &read[observed.start as usize..observed.end as usize],
            b"CAGCAGCTGCAG"
        );
    }

    /// **Mandatory property 3 — the shared model, and direction asymmetry.** Algorithm 4 reads
    /// the passed `StutterModel`; a contraction-biased model makes a one-unit contraction
    /// *cheaper to enter* than a one-unit expansion. Verified at the cost level, since that is
    /// where the asymmetry lives.
    #[test]
    fn a_contraction_is_cheaper_to_open_than_an_expansion() {
        let slip = SlipCosts::from_model(&contraction_biased());
        assert!(
            slip.open_contraction > slip.open_expansion,
            "a contraction-biased model must make contraction the cheaper slip to open"
        );
    }

    /// The slip costs reconstruct the stutter model's own probabilities exactly — the affine
    /// open + extends must sum back to `ln(P(n)) − ln(equal)`, or the aligner is pricing a
    /// different distribution than the one it was handed.
    #[test]
    fn the_slip_costs_reconstruct_the_stutter_probability() {
        let model = contraction_biased();
        let slip = SlipCosts::from_model(&model);
        let period = std::num::NonZeroU8::new(3).unwrap();
        let ln_equal = model.equal().ln();

        for n in 1..=5i64 {
            // Expansion of n units: open + (n−1) extends, relative to equal.
            let reconstructed = slip.open_expansion + (n - 1) as f64 * slip.extend;
            let expected = model.probability(n * 3, period).ln() - ln_equal;
            assert!(
                (reconstructed - expected).abs() < 1e-12,
                "expansion of {n} units diverged: {reconstructed} vs {expected}"
            );
            // Contraction likewise.
            let reconstructed = slip.open_contraction + (n - 1) as f64 * slip.extend;
            let expected = model.probability(-n * 3, period).ln() - ln_equal;
            assert!(
                (reconstructed - expected).abs() < 1e-12,
                "contraction of {n} units diverged"
            );
        }
    }

    /// Measure the same read with algorithm 3 (the flat-gap delimiter).
    fn measure_flat(
        read: &[u8],
        reference: &[u8],
        geometry: &RepeatGeometry,
        model: &StutterModel,
    ) -> RepeatSpan {
        use crate::ng::alignment::ssr_best_path_flat_gap::{SsrFlatGapAligner, ViterbiScratch};
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let context = RepeatContext {
            geometry,
            stutter: model,
        };
        let quality = vec![35u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = ViterbiScratch::new();
        aligner.align(bases, reference, context, &mut scratch)
    }

    /// **The strongest correctness check available for a no-oracle algorithm: cross-check
    /// against algorithm 3, which is byte-parity-validated against production.**
    ///
    /// On a *clean* read — the reference verbatim, or a clean whole-unit expansion or
    /// contraction whose flanks unambiguously anchor the tract — the two aligners must measure
    /// the **same span**. The slip pricing changes the *score* of the interior alignment, but
    /// with the flanks pinning both junctions it cannot move where they land. Where the two
    /// *do* diverge is on ambiguous or interrupted reads, and measuring that divergence is
    /// exactly the Milestone-D comparison's job — not something to assert here.
    #[test]
    fn algorithm_4_agrees_with_algorithm_3_on_clean_reads() {
        let model = contraction_biased();
        // Motifs of several periods, including a homopolymer (period 1) and a period-6, since
        // the slip stride is the period and the traceback junction handling is period-sensitive
        // (a clean-read differential first caught a contraction-at-the-tract-start bug here).
        let motifs: &[&[u8]] = &[b"A", b"CA", b"CAG", b"CAGT", b"CAGGA", b"GTTGTG"];
        let left_flanks: &[&[u8]] = &[b"", b"ACGT", b"ACGTACGT"];
        let right_flanks: &[&[u8]] = &[b"", b"TTGG", b"TTGGTTGGAT"];

        for motif in motifs {
            let period = motif.len();
            for left in left_flanks {
                for right in right_flanks {
                    let ref_units = 4usize;
                    let ref_tract: Vec<u8> = motif
                        .iter()
                        .cycle()
                        .take(period * ref_units)
                        .copied()
                        .collect();
                    let (reference, geometry) = frame(left, &ref_tract, right, motif);

                    // The reference itself, plus clean whole-unit expansions and contractions.
                    for delta in [-3i32, -2, -1, 0, 1, 2, 4] {
                        let read_units = (ref_units as i32 + delta).max(0) as usize;
                        let mut read = left.to_vec();
                        read.extend(motif.iter().cycle().take(period * read_units).copied());
                        read.extend_from_slice(right);

                        let four = measure(&read, &reference, &geometry, &model);
                        let three = measure_flat(&read, &reference, &geometry, &model);
                        assert_eq!(
                            four,
                            three,
                            "algorithm 4 and 3 disagreed: motif {:?} left {} right {} read_units {read_units}",
                            std::str::from_utf8(motif).unwrap(),
                            left.len(),
                            right.len()
                        );
                    }
                }
            }
        }
    }

    /// Algorithm 4 must be usable through the priced slip route without panicking or
    /// diverging on a real unit expansion — the route the whole algorithm exists for.
    #[test]
    fn a_priced_slip_expansion_is_measured_at_its_length() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCAGCAGCAGTTGGTTGGAT"; // +1 unit (15-base tract)
        let span = measure(read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(15));
    }

    /// **The load-bearing discriminants.** `State` writes score arrays by position and reads
    /// them by `index()`, so reordering the variants would silently permute every cell — the
    /// same guard algorithm 3 carries for its three states.
    #[test]
    fn the_state_discriminants_are_the_array_indices() {
        assert_eq!(State::Match.index(), 0);
        assert_eq!(State::Insertion.index(), 1);
        assert_eq!(State::Deletion.index(), 2);
        assert_eq!(State::SlipInsertion.index(), 3);
        assert_eq!(State::SlipDeletion.index(), 4);
    }

    /// **A slipped unit's bases are each scored at their own quality**, not the current row's.
    /// A whole-unit insertion spans `period` read rows, and each base sits at read index
    /// `row_index − period + k`; scoring all `period` of them at the row's single quality
    /// skews the slip score whenever intra-unit qualities vary. The fix indexes each base's
    /// own quality by its own position — correct by construction.
    ///
    /// A pure black-box assertion is weak here (a score-level shift rarely flips the measured
    /// span), so this exercises the slip route under a **non-uniform** quality profile —
    /// which no other test does — and pins that the call stays sound and still measures the
    /// expansion. It is the varied-quality coverage the fix needs; the mechanism itself is
    /// verified by inspection of the `quality_at(idx)` / `bases[idx]` pairing.
    #[test]
    fn a_slipped_unit_is_scored_under_non_uniform_quality() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCAGCAGCAGTTGGTTGGAT"; // +1 clean unit
        let aligner = SsrUnitSlipAligner::new(PerQualityEmission::new());
        let model = contraction_biased();
        let context = RepeatContext {
            geometry: &geometry,
            stutter: &model,
        };

        // A jagged quality profile, low across the inserted unit's rows.
        let mut quality = vec![35u8; read.len()];
        for (offset, q) in quality.iter_mut().enumerate() {
            if (20..23).contains(&offset) {
                *q = 3;
            }
        }
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = UnitSlipScratch::new();
        let span = aligner.align(bases, &reference, context, &mut scratch);
        assert_eq!(span.measured_length(), Some(15));
    }

    /// A read running off its 3′ end anchors only the left flank — the four `RepeatSpan` cases
    /// must be reachable from algorithm 4 too, since it shares the classification.
    #[test]
    fn a_read_off_its_end_anchors_only_the_left_flank() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let span = measure(
            b"ACGTACGTCAGCAGCAG",
            &reference,
            &geometry,
            &contraction_biased(),
        );
        assert!(matches!(span, RepeatSpan::FromLeft(_)));
    }

    /// An empty reference is a defined answer, not an error.
    #[test]
    fn an_empty_reference_is_unanchored() {
        let geometry = RepeatGeometry {
            left_flank_len: Bp(0),
            right_flank_len: Bp(0),
            motif: Motif::new(b"CAG").expect("valid"),
        };
        assert_eq!(
            measure(b"ACGT", b"", &geometry, &contraction_biased()),
            RepeatSpan::Unanchored
        );
    }

    /// Scratch reuse across a size drop must not leak — the ring is grown, never cleared, so a
    /// smaller read after a larger one must give the fresh-scratch answer.
    #[test]
    fn scratch_reuse_does_not_leak_across_periods() {
        let aligner = SsrUnitSlipAligner::new(PerQualityEmission::new());
        let model = contraction_biased();

        // Loci of assorted periods and sizes — the ring length is `period + 1`, so it *changes*
        // between calls, which is the hard case: a slot written by a longer, higher-period
        // alignment must never be read stale by a shorter, lower-period one.
        let loci: &[(&[u8], &[u8])] = &[
            (b"GTTGTG", b"GTTGTGGTTGTGGTTGTGGTTGTGGTTGTG"),
            (b"A", b"AAAA"),
            (b"CAG", b"CAGCAGCAGCAG"),
            (b"CA", b"CACACACACACA"),
        ];

        let mut reused_scratch = UnitSlipScratch::new();
        for &(motif, tract) in loci {
            let (reference, geometry) = frame(b"ACGTACGT", tract, b"TTGGTTGGAT", motif);
            let quality = vec![35u8; reference.len()];
            let context = RepeatContext {
                geometry: &geometry,
                stutter: &model,
            };
            let bases = ReadBases::try_new(&reference, &quality).expect("matched lengths");
            let reused = aligner.align(bases, &reference, context, &mut reused_scratch);
            let bases = ReadBases::try_new(&reference, &quality).expect("matched lengths");
            let fresh = aligner.align(bases, &reference, context, &mut UnitSlipScratch::new());
            assert_eq!(
                reused,
                fresh,
                "scratch reuse leaked across periods at motif {:?}",
                std::str::from_utf8(motif).unwrap()
            );
        }
    }

    /// Multiple periods work — a period-1 (homopolymer) and a period-2 locus, since the slip
    /// stride is the period and period 1 is where the indel deficit lives (spec §4.2).
    #[test]
    fn homopolymer_and_dinucleotide_loci_are_measured() {
        // Period 1.
        let (reference, geometry) = frame(b"ACGTACGT", b"AAAAAA", b"TTGGTTGGAT", b"A");
        assert_eq!(
            measure(&reference, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(6)
        );
        // A one-base expansion of the homopolymer, a genuine unit slip at period 1.
        let read = b"ACGTACGTAAAAAAATTGGTTGGAT";
        assert_eq!(
            measure(read, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(7)
        );

        // Period 2.
        let (reference, geometry) = frame(b"ACGTACGT", b"CACACACA", b"TTGGTTGGAT", b"CA");
        assert_eq!(
            measure(&reference, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(8)
        );
    }
}
