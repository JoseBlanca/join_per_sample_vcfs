# Code Review: ng_alignment_b3
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step B3 of `doc/devel/ng/impl_plan/alignment_best_path.md` — byte-parity against production, the port anchor
**Status:** Approve-with-changes (all applied — see §Author response)

---

### 1. Scope

The uncommitted diff for B3 (B2 committed as `e1e676d`): new test-only
[delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs), plus a `#[cfg(test)] mod`
declaration and a module-doc refresh in `mod.rs`.

One reviewer, scoped to the question that actually matters for a parity harness: **can it fail?** A
parity test that passes vacuously is worse than none, because it manufactures confidence exactly
where the plan says confidence matters most.

### 2. Verdict

**Approve-with-changes.** **3 Majors**, no Blockers. The oracle is real — but it had two reachable
blind spots and one generator hole, all three found by *running experiments*, not by reading.

### 3. The oracle is real: 5 of 7 mutations caught, 2 slipped through

The reviewer mutated the aligner in ways a genuine porting error would look like, ran each in the
container, and reverted:

| mutation | caught? |
|---|---|
| swap tract/flank gap-open | **yes**, case 0 |
| tract window right edge `<=`→`<` | **yes**, case 31 |
| tract window left edge `<`→`<=` | **yes**, case 7 |
| traceback `tract_start` off by one | **yes**, case 0 |
| reverse the match `best_of` order | **yes** |
| **transpose `left_anchored`/`right_anchored`** | **NO — all three tests passed** |
| **row-0 init loses the tract-aware gap-open** | **NO — passed at 200,000 cases too** |

### 4. Findings

#### Major

**B3-1: the harness could not see a transposed side** — the very failure the aligner's own docs call
"the whole risk of this step, and its failure is silent". It is invisible to *parity* by
construction: production's `BorderOffEnd` discarded which flank was missing, so swapping the two
anchors satisfies every parity assertion.
**Fix applied:** a **self-consistency** check, which is not a tautology because the side is
independently derivable from ng's own offsets — a read that ran off its 3′ end has its span reaching
the read's end; one that began inside the repeat has its span starting at zero.
**And the first version of that check was wrong**, caught by its own test: both implications hold
only *when the flank exists*. A side reports "not anchored" for two different reasons — the read ran
past that junction, **or** there is no flank there to anchor against — and only the first constrains
the offsets. Now conditioned accordingly. **Verified: the transposition mutation now fails both
parity tests at the first generated case.**

**B3-2: zero-length flanks were never generated**, and that is precisely where ng's B2 widening
*deliberately* disagrees with production. Confirmed by hand: `left=0, right=0` gives production
`Region(0..12)` against ng `Unanchored`, which would have made the harness **panic as a parity
failure** on an intentional divergence.
**This also falsified a claim in B2's own source comment** — that the flankless behaviour "does not
disturb the parity oracle". That only looked true because the fixture never produced the case.
**Fix applied:** zero-length flanks are now generated (~1 in 8 per side); the divergence is
*asserted* rather than tripped over (ng must not claim a measurement; where it does report a span the
bytes must still match production exactly); a `flankless_divergences > 50` floor keeps the assertion
from going dead; and the B2 comment is corrected to say what is actually true.

**B3-3: a fresh `ViterbiScratch` per case** made the only randomized, size-varying driver in the
module structurally blind to the grow-but-never-clear stale-cell hazard that
`ssr_best_path_flat_gap` names as *the* Milestone-C risk.
**Fix applied:** one scratch hoisted across the whole run — free, and now thousands of
differently-sized reads and frames flow through the same buffers.

#### Minor (applied)

- **The `≥2 of 3` off-end shape assertion was looser than reality** — the measured distribution
  reaches all three (`FromLeft` 1,759 / `FromRight` 1,255 / `Unanchored` 252 over 12,000 cases).
  Tightened to require all three, so generator drift is caught.
- **No way to run the big sweep without editing source.** The port was justified by ~800,000 cases;
  the permanent harness runs 12,000 as a *sentinel* — defensible, since every recurrence mutation
  died within 40 cases, but the soak should be one command away. Added `PVC_PARITY_CASES`.
  **Re-run at 50,000/seed (200,000 cases, release): green.**

#### Verified sound, no change needed

- **Production is driven correctly.** `Locus::new(.., 1.0, left+tract+right, 0)` yields exactly the
  `left_flank()`/`right_flank()`/`ref_bytes()` the delimiter reads; `purity_fraction` only feeds
  `is_perfect()`, which `delimit_read` never calls, so `1.0` makes nothing easier. Random flanks
  occasionally extending the tract is not a shared blind spot — both implementations see the
  identical frame and read, so it is a legitimate shared input.
- **Reproducibility works** — the reviewer replayed seed `0x5eed0001` index 7 under a mutation and
  reproduced the reported divergence exactly.
- **The earlier unreproducible divergence** (B1's review) now has ~1M cases of clean evidence behind
  the stale-incremental-artifact explanation.

#### Known remaining gap, recorded

**The row-0 initialisation mutation still slips through** even at 200,000 cases. Row 0 is the
all-deletion prefix, and losing its tract-aware gap-open evidently never changes a measured offset on
generated input. Recorded rather than papered over: it is a genuine hole in what this oracle can
detect, and the honest statement is that parity covers the recurrence body, not every cell of the
initialisation.

### 5. Commands to re-verify

```
./scripts/dev.sh cargo test --lib ng::alignment::delimit_parity
./scripts/dev.sh env PVC_PARITY_CASES=50000 cargo test --release --lib ng::alignment::delimit_parity
```

### Author response

All three Majors and both Minors **fixed in this step's commit**; the row-0 gap is recorded above and
in `PROJECT_STATUS.md`. Audit trail at `tmp/review_2026-07-23_ng-alignment-b3/`.
