# ng alignment — Milestone D: algorithm 4, the comparison entry

**Date:** 2026-07-23
**Plan:** [alignment_best_path.md](../../ng/impl_plan/alignment_best_path.md), Milestone D
**Design authority:** spec [alignment.md](../../ng/spec/alignment.md) §4.2, §5.2, §10.3; arch
[alignment.md](../../ng/arch/alignment.md) §5
**Process:** plan-driven-implementation — implement → review → apply-fixes → commit.

Milestone D is **research, off the STR generator's critical path** (which was unblocked at Checkpoint
B), and its two steps are different in kind from A–C: **algorithm 4 is genuinely new code with no
parity oracle** (arch §5: "no production counterpart"; algorithms 2, 4, 5, 6 are "measured, not
verified"), and D2 is a comparison *report*. Because the design left the recurrence's concrete matrix
embedding open, the model was **settled with the owner before building** (see below).

---

## Step D1 — the two-penalty aligner

**Status:** shipped (reviewed, fixes applied).
**Review:** [ng_alignment_d1_2026-07-23.md](../reviews/ng_alignment_d1_2026-07-23.md) — no Blockers,
1 Major, several Minors; all applied.

### 1. The design, settled with the owner

The spec specifies algorithm 4's *properties* (direction-asymmetric, out-of-frame route, shares the
`StutterModel`) but not the Viterbi embedding, and a best-path aligner cannot do what production's
candidate-scoring model does (sum over placements, divide by count). Four decisions were taken with the
owner before writing code:

1. **Slipped-unit bases scored against the motif tiling** (production-aligned: rewards a real unit,
   catches a wrong base as a mismatch).
2. **Placement multiplicity: best placement, no divide** — a max aligner picks the best-scoring
   placement; the max-vs-sum gap is a documented difference for D2 to measure, not a bug.
3. **Slip cost relative to no-slip** (a derivation, not a free choice): only relative scores matter in
   an argmax, and Δ = 0 has no transition to hang the `equal` mass on, so `slip_open = ln(direction ·
   geom) − ln(equal)`, baseline 0. This decomposes the stutter geometric exactly onto affine
   open/extend transitions.
4. The aligner is, by design, a **ruler that prices toward tidy in-frame lengths** — the §4.2 risk D2
   exists to measure, intended not defective.

### 2. What shipped

[ssr_best_path_unit_slip.rs](../../../../src/ng/alignment/ssr_best_path_unit_slip.rs) — algorithm 3's
M/I/D matrix plus two whole-unit slip states (`SlipInsertion`, `SlipDeletion`), reachable only inside
the tract, each jumping `period` cells. The score scratch keeps a ring of the last `period + 1` rows
(insertions look back `period` rows). `SlipCosts` derives the affine slip costs from the shared
`StutterModel` (accessed by reference, not copied). It shares algorithm 3's out-of-frame gap model via
`pub(super)` accessors.

### 3. Validation — a differential against the one trusted aligner

With no oracle, the strongest check is the **differential against algorithm 3** (byte-parity-validated
against production): on clean reads the two must measure the same span. Shipped as
`algorithm_4_agrees_with_algorithm_3_on_clean_reads` (378 cases); the review broadened it to 3666 and
found **zero disagreements** on clean junctions, the only divergences being the designed ruler-vs-score
gap. Plus 12 property tests: the three mandatory behaviours (direction asymmetry, out-of-frame route,
shared model), the slip costs reconstructing the stutter probability exactly, the four `RepeatSpan`
cases, degenerate inputs, cross-period scratch reuse, and per-base slip quality.

### 4. The bugs the process caught

- **The differential caught a real traceback bug** during implementation: a contraction at the very
  start of the tract deletes the first tract base, which *is* the left junction, and the arm was not
  recording `tract_start`. No property test would have found it; the cross-check against algorithm 3
  did, immediately.
- **The review caught a per-base-quality bug** (D1-1): a slipped unit's `period` bases were all scored
  at the current row's quality instead of each base's own — invisible to uniform-quality tests, a wrong
  slip score when qualities vary. Fixed.

Both are the same lesson as the rest of the plan: on a no-oracle algorithm, correctness comes from
cross-checking against the trusted sibling and from a review that reasons about the recurrence — not
from the code looking right.

### 5. Validation results

Container: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -D warnings` exit 0;
`cargo test --lib` **2220 passed, 0 failed, 4 ignored**; algorithm 3's 200,000-case parity soak still
byte-identical (the accessors added to its file did not disturb it).

### 6. D2 — deferred, to reassess with the owner

D2 is the 3-vs-4 comparison on synthetic reads — a research *report* measuring **calibration and
accuracy**, split by period, and at period 1 split by whether the indel is the repeat's own base (spec
§10.3). It is a distinct deliverable (a simulator + a measurement harness + a written finding), and it
is explicitly off the critical path. **Reassess with the owner before starting it.**
