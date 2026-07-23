# ng alignment — Milestone C: banding

**Date:** 2026-07-23
**Plan:** [alignment_best_path.md](../../ng/impl_plan/alignment_best_path.md), Milestone C
**Design authority:** spec [alignment.md](../../ng/spec/alignment.md) §3, §9, §10.3; arch
[alignment.md](../../ng/arch/alignment.md) §4, §6
**Process:** plan-driven-implementation — implement → review → apply-fixes → commit.

Milestone C is a single step, C1, then Checkpoint C.

---

## Step C1 — band the delimiter

**Status:** shipped (reviewed, fixes applied).
**Review:** [ng_alignment_c1_2026-07-23.md](../reviews/ng_alignment_c1_2026-07-23.md) — no Blockers, no
Majors; the band survived every mutation attack. All findings applied.

### 1. Plan

Restrict the delimiter's matrix to a per-read band — a floor of `|read − reference|` plus headroom —
as a **performance change that provably alters no output**. Isolated because its failure is silent: a
too-narrow band loses long alleles with no error. *Source:* spec §3, §10.3, arch §4, §6.

### 2. The step's finding — the band is wider than the spec's model, and the spec is incomplete

**This step raised a design question mid-implementation and it was taken to the owner (Checkpoint-C
decision, resolved "adopt the geometry band").**

Spec §3 models the band as `|read − reference|` (a forced floor) plus a *small constant* headroom, and
§9 lists the amount as open. **That model is insufficient for this aligner, and it fails silently.**
The parity oracle caught it: a read that expands the tract *and* runs off a flank has an optimal path
whose **peak** matrix deviation is the tract expansion, while its **net** `|read − reference|` is
smaller by the run-off flank length — so the path strays past the floor by up to the flank length, far
more than a constant. (Case seed `0x5eed0001`/1745: a left-flankless 43-base read of a 7-unit tract
against a 4-unit, 34-base frame; peak deviation 18, floor 9.)

The correct band is therefore three terms:

```
band = |read − reference|            (the forced floor — spec §3, mandatory)
     + left_flank + right_flank      (the run-off correction — this step's finding)
     + BAND_HEADROOM (= 8)           (the score-neutral bow slop — a named constant)
```

Every term is bytes of alignment **geometry**; none is the scoring cutoff, so the spec's separation of
the ruler from `MAX_SLIP` holds.

**`SPEC-FOLLOWUP(alignment §3/§9)` — for the owner.** The spec should be amended to the three-term
band. The owner chose (Checkpoint C) to *note* the correction rather than have this step edit the
design doc; a grep-able marker sits on `BAND_HEADROOM`, and this is the record of it.

### 3. Decisions

- **Out-of-band cells are written `UNREACHABLE`, not skipped.** This mirrors production's banded
  forward, and it is what makes banding safe on a *reused* scratch: the unbanded fill was exhaustive,
  so a skip-don't-write band would let a larger previous read's values leak in as a neighbour score.
  This closes the grow-without-clear note B1 raised for exactly this milestone.
- **`BAND_HEADROOM = 8`, a 4× margin.** The review established the empirical minimum is 2 (a headroom
  of 1 passes 12k but fails the 200k soak at case 28,307). 8 is deliberate slack against a silent
  failure, narrowable only against a measured perf need.
- **`band_width` is a free function**, extracted from an inline closure so the formula the port's
  correctness rests on is unit-testable in isolation.

### 4. Tests

The plan's required pair — the parity fixture unchanged (the committed 12k/200k differential oracle,
which is what caught the run-off gap) and a long-allele extreme
(`a_long_allele_at_the_extreme_is_measured_not_collapsed`) — plus the case-1745 finding pinned as a
fixed-input regression, the band-formula unit tests, and the two the reviewer proposed (both-ends
run-off; scratch reuse across an extreme size drop).

### 5. Validation

Container: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -D warnings` exit 0;
`cargo test --lib` **2207 passed, 0 failed, 4 ignored**; the 200,000-case parity soak
(`PVC_PARITY_CASES=50000`, release) byte-identical to the unbanded delimiter.

### 6. The lesson, consistent with the rest of the plan

The band was validated by the same oracle that validated the port: **byte-parity with production's
unbanded delimiter.** The one gross error I made (a constant headroom) was caught by that oracle, not
by inspection — the reason the plan built the oracle first. And the review's confidence came from
*mutating each band term and watching the oracle fail*, not from reading the formula.
