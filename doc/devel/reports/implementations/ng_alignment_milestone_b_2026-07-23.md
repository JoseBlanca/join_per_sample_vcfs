# ng alignment — Milestone B: algorithm 3, the parity anchor

**Date:** 2026-07-23
**Plan:** [alignment_best_path.md](../../ng/impl_plan/alignment_best_path.md), Milestone B
**Design authority:** spec [alignment.md](../../ng/spec/alignment.md) §4.2, §10.3; arch
[alignment.md](../../ng/arch/alignment.md) §5
**Process:** plan-driven-implementation — one implement → review → apply-fixes → commit loop per step.

**This is the milestone that unblocks the STR locus generator.** Its `align_read` is what
[locus_generation_ssr.md](../../ng/impl_plan/locus_generation_ssr.md) Milestone D has been waiting
on, and it is the alignment module's **only byte-parity oracle**, so every later comparison gets a
fixed reference point from it.

---

## Step B1 — the two-regime matrix and traceback, unbanded

**Status:** shipped (reviewed, fixes applied).
**Review:** [ng_alignment_b1_2026-07-23.md](../reviews/ng_alignment_b1_2026-07-23.md) — **3 Blockers,
2 Majors**, all applied. None of them in the transcription.

### 1. Plan

Port `delimit_read`: match/insertion/deletion states, affine gaps, and the **tract-aware gap-open**
(stiff in the flanks, soft inside the repeat) keyed by which reference column a gap touches. Fill the
whole matrix, as production does — **unbanded is deliberate** and is what makes parity provable. Carry
the two determinism rules. *Source:* spec §4.2, arch §5.

### 2. Decisions

- **The port is faithful, and that is now demonstrated rather than asserted.** A reviewer read
  production and the port side by side across nine checkpoints, then built a temporary in-crate
  differential harness — production's `delimit_read` is `pub(crate)`, so this is possible — and ran
  **~800,000 randomized cases across four seeds, debug and release: zero divergences**, ~570,000 of
  them production `Region` results measured to the same two offsets.
- **`delimit` returns a raw `TractReadout`, and there is no `BestPathAligner` impl yet.** The plan
  splits B1 (recurrence) from B2 (the `RepeatSpan` widening), and a stub `align` would have been
  dead code carrying an `unimplemented!()` — a release panic waiting for its first caller. B2 adds
  the impl.
- **The match→match inconsistency is reproduced, not fixed** (spec §4.2 requires the choice to be
  deliberate). Inside the tract the three transitions leaving a match sum to ~1.02, biasing toward
  the reference length by roughly 1.02 per gapped base. Reproducing it silently and fixing it
  silently both end with a parity test failing for a reason nobody can find.
- **Production's `tract_start > tract_end` rejection is omitted, and the omission is now argued.**
  Two reviewers independently established it is unreachable — one structurally, one by exhaustive
  sweep. Documented as a lemma with a `debug_assert!` and a sweeping test rather than reinstated as a
  redundant guard.
- **The state-index-to-enum refactor was declined for now.** A genuine improvement (it would remove
  15 casts and close the traceback's silent `_ =>` arm), but it is a refactor of a port whose
  byte-identity was just established over 800,000 cases, and this plan's rule is transcribe first,
  change separately. A `debug_assert_eq!` closes the immediate hazard.

### 3. Changes made

New [ssr_best_path_flat_gap.rs](../../../../src/ng/alignment/ssr_best_path_flat_gap.rs) —
`SsrFlatGapAligner`, `TransitionCosts`, `ViterbiScratch`, `TractReadout`, `best_of`, `delimit`; plus
the module declaration in `mod.rs`.

### 4. Tests — 17

**The lesson of this step is that the tests were the weak part, not the code.** All three Blockers
were found by *mutating the source and re-running*, and each mutation survived the original ten
tests: reversing the match tie-break to I > D > M; transposing the two flank lengths; widening the
tract window's right edge.

**And the first fix for the tie-break did not work either.** A test built on an emission model where
a match and a mismatch score alike still failed to catch the mutation, because the three match
candidates stay strictly ordered by their *transition* costs — with a strict winner, the candidate
order is never consulted. The shipped test flattens both the transitions and the emissions, so every
path scores identically and only the preference order decides. **Verified by re-running the mutation
against it: it fails under the mutation and passes without it.**

### 5. Validation results

Container: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -D warnings` exit
0; `cargo test --lib` **2185 passed, 0 failed, 4 ignored**.

### 6. Tradeoffs and follow-ups

- **One unreproduced anomaly, recorded rather than dismissed.** Two early runs of the fidelity
  reviewer's small sweep reported disagreements; five of the reported cases replay as *agreeing* when
  pinned, and 18 later runs including six at the identical seed were clean. Stale incremental
  artifacts are the plausible explanation. **B3's permanent parity harness is where this would
  resurface if it is real.**
- **Scratch's grow-without-clear is safe only because the fill is exhaustive.** Documented, because
  **Milestone C's banding makes it a live hazard** — once cells go unwritten, the leftovers become
  reachable.
- The state enum refactor, deferred above.
