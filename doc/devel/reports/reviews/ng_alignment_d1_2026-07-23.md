# Code Review: ng_alignment_d1
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step D1 of `doc/devel/ng/impl_plan/alignment_best_path.md` — algorithm 4, the two-penalty aligner
**Status:** Approve-with-changes (all applied — see §Author response)

---

### 1. Scope

The uncommitted diff for D1 (HEAD `f96b213`): new
[ssr_best_path_unit_slip.rs](../../../../src/ng/alignment/ssr_best_path_unit_slip.rs) (algorithm 4),
`pub(super)` accessors added to `TransitionCosts` in
[ssr_best_path_flat_gap.rs](../../../../src/ng/alignment/ssr_best_path_flat_gap.rs) (additive; algorithm
3's parity soak still passes), and a one-line `mod.rs` change.

Two reviewers: one attacking the **recurrence's correctness by experiment** (this algorithm has **no
parity oracle** — arch §5: algorithms 2, 4, 5, 6 are "measured, not verified"), one on code quality.

### 2. Verdict

**Approve-with-changes.** No Blockers. **One Major** (found by the correctness reviewer, applied), the
rest Minor/Nit. The recurrence survived every correctness attack.

### 3. The recurrence is sound — established by experiment against the trusted sibling

Algorithm 4 has no external oracle, so the strongest evidence is the **differential cross-check against
algorithm 3**, which *is* byte-parity-validated against production. On a clean read the two must
measure the same span (the flanks pin both junctions), and the reviewer confirmed this at scale:

- **Broadened differential** (8 motifs including period 1 and 6, 5×5 flanks, tracts of 3/5/8 units,
  deltas −8…+9): **3666 clean-junction cases agreed, 0 disagreements.**
- The one disagreement class that surfaced is the **designed** ruler-vs-score divergence, not a bug:
  it occurs only when a flank base extends the repeat run so the boundary can legally slide — exactly
  what the module doc and the shipped differential's own comment call out. A non-colliding flank agrees.
- **Hard scratch reuse** (jagged size sequences on one scratch, period-6→period-1 and back): every
  reused-scratch result equalled a fresh-scratch result — **no stale-ring leak.**
- **Runs-off-both-ends and Q0 reads** across periods 1/2/3/6: no panic, correct `Unanchored` — **no
  period-jump underflow.**
- A **direction-swap mutation** (`open_expansion`/`open_contraction`) is caught by the cost-level
  tests; correctly *not* caught by the clean-read differential, since clean lengths are direction-blind.

### 4. Findings

#### Major (applied)

**D1-1: a slipped unit's bases were all scored at one shared quality, not each base's own.**
A whole-unit insertion spans `period` read rows, but the emission summed `scores.pick(bases[row_index −
period + k], …)` using `scores` — resolved once from the *current* row's quality — for all `period`
bases. Invisible to every test (all use uniform quality), but it silently skews the slip-route score
whenever intra-unit qualities vary: a wrong measurement, with no oracle to catch it.
**Fix applied:** each inserted base is now scored at its own quality, `scores_for(quality_at(idx))`
with `idx = row_index − period + k` — correct by construction (the base and its quality share the same
index). A varied-quality slip test was added for coverage; a purely black-box assertion is weak here
(a score-level shift rarely flips the measured span), and the report says so plainly rather than
over-claiming.

#### Minor / Nits (applied)

- **The `State` discriminants were unpinned** — the 5-variant enum writes score arrays by position and
  reads by `index()`, so a reorder would silently permute them. Algorithm 3 guards this;
  `the_state_discriminants_are_the_array_indices` now does too.
- **`trace_back` took 8 positional arguments** under `#[allow(clippy::too_many_arguments)]`, six of them
  interchangeable `usize`s (a `read_len`/`reference_len` or flank transposition would be a silent wrong
  measurement). Replaced with a `Copy` `MatrixGeometry` struct destructured at the top — the `allow` is
  gone and the transposition hazard with it. It also became a free function (it used no `self`).
- **A clumsy `split_at_mut(1)` + `let _ = &from` row-0 workaround** replaced by a plain `&mut rows[0]`.
- **The `fits_reference` precondition assert** (present in algorithm 3) was missing; restored.
- **The scratch-reuse test was single-period/single-drop**; broadened to `scratch_reuse_does_not_leak_across_periods`,
  which is the harder case since the ring length is `period + 1` and changes between calls.

#### Assessed and accepted with reason

- **The two same-named `State` enums** (3-variant in algorithm 3, 5-variant here) are genuinely
  different state spaces, both private — correct, not confusing.
- **`best_of` and `UNREACHABLE` are duplicated** between the two aligners. `best_of` cannot be trivially
  shared because the two `State` types differ (it would need to be generic); left as a small, isolated
  duplicate rather than reaching into the committed algorithm-3 file for a generic lift this step does
  not need. Recorded as a candidate cleanup.
- **`TransitionCosts`/`TractReadout` living in algorithm 3's file** while both aligners consume them —
  the `pub(super)` accessor seam's *visibility* is right; a future lift to the module root is the
  natural home, deferred with the port-first convention.
- **`delimit`'s length** is justified as a flat DP fill, consistent with algorithm 3's sibling.
- **The `Vec<Vec<[f64; STATES]>>` ring** is a nested allocation, less flat than the sibling's rolling
  rows; correct as written, flagged as a possible flattening but not reworked (reworking correct DP
  risks a regression for a small perf gain).

### 5. Missing tests added

`the_state_discriminants_are_the_array_indices`, `a_slipped_unit_is_scored_under_non_uniform_quality`,
`scratch_reuse_does_not_leak_across_periods`, and the broadened
`algorithm_4_agrees_with_algorithm_3_on_clean_reads` (378 cases: 6 motifs × 3 left × 3 right flanks × 7
deltas).

### 6. What's good

- The module doc derives the recurrence and states all four design decisions with their owner sign-off,
  and is candid that the algorithm is *measured, not verified* — the reader knows exactly what the tests
  do and do not establish.
- The differential against algorithm 3 is the right validation shape for a no-oracle algorithm, and it
  **caught a real traceback bug during implementation** (a contraction at the tract start deletes the
  first tract base, which is the left junction — the arm was not recording it).
- `SlipCosts::from_model` reads the shared model rather than copying its parameters, satisfying the
  spec's write-once-and-share requirement, and `the_slip_costs_reconstruct_the_stutter_probability`
  pins that the affine decomposition reproduces the model's own distribution.

### 7. Commands to re-verify

```
./scripts/dev.sh cargo test --lib
./scripts/dev.sh env PVC_PARITY_CASES=50000 cargo test --release --lib ng::alignment::delimit_parity
```

### Author response

The Major and every actioned Minor **fixed in this step's commit**; the accepted items are recorded
above with their reasons. The recurrence's correctness rests on the differential against algorithm 3,
which the review exercised at 3666 cases rather than trusting.

Per-category audit trail at `tmp/review_2026-07-23_ng-alignment-d1/`.
