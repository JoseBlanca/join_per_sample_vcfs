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

---

## Step B2 — the `RepeatSpan` readout, the widening

**Status:** shipped (reviewed, fixes applied).
**Review:** [ng_alignment_b2_2026-07-23.md](../reviews/ng_alignment_b2_2026-07-23.md) — 3 Majors, no
Blockers; all applied.

### 1. Plan

Walk the traceback and read the repeat off the two flank-junction columns, **reporting which flanks
anchored** — the part production's side-blind `Delimited::BorderOffEnd` cannot express, and the part
the STR preparer needs to tell a measurement from a lower bound. *Source:* spec §4.2, arch §2.1, §5.

### 2. Decisions

- **`TractReadout::classify` maps the 2×2 of anchoring onto the four `RepeatSpan` cases**, and the
  trait impl wires it. Verified by mutation, not argument: swapping `FromLeft`/`FromRight` is caught
  by 2 tests, swapping the anchor initialisers by 3, transposing the endpoints by 6.
- **The trait's `Context` became a generic associated type passed by value** (`type Context<'a>: Copy`).
  The arch sketch's plain `type Context` taken by reference cannot name `RepeatContext<'a>` without
  forcing `RepeatContext<'static>`, which would make every caller own its locus data forever — the
  opposite of arch §2.2's reason for borrowing it. Both reviewers assessed the alternatives and upheld
  the GAT; one verified by compilation that `type Context<'a> = ()` still works for the gated affine
  aligner.
- **A flank that does not exist cannot anchor** — the review's main finding, and a correction to B1.

### 3. Tests — 8 added

Four `RepeatSpan` cases all reachable from real alignments; a clean repeat measures exactly; a long
allele is measured at its own length; a truncated read **bounds** rather than measures; an interrupted
repeat verbatim; the flankless cases.

### 4. Validation

Container: fmt clean, clippy clean, `cargo test --lib` **2194 passed, 0 failed**.

### 5. The finding worth carrying forward

**A flankless geometry turned an unmeasured read into a measurement**, and the cause was a convention
I wrote in B1: "a flank that does not exist cannot fail to hold". Harmless while it was two raw bits;
**B2 is the step that promotes those bits to a claim about length**, and `Between` is the one variant
that says the length is pinned. Demonstrated live — frame `("", "CAGCAGCAGCAG", "")`, read `CAGCAG` →
`measured_length() == Some(6)`, a short allele nothing observed, at contig-edge loci where the
evidence is already thinnest.

Inverted: **a flank that does not exist cannot anchor.** At a flankless edge, a read that ended
because the tract ended and one that ran out mid-tract produce an *identical* readout, so the
measurement claim over-claims for one of them. A flankless locus can therefore never yield a
measurement — honest, since there is no flank to pin the edge against. **Parity is undisturbed:**
production's `Region` makes no measurement claim, and B3 compares the measured bytes.

The general lesson, now three steps running: **the code was right and the test was not.** One of this
step's own tests (`the_classified_span_agrees_with_the_raw_readout`) was a tautology — `align` *is*
`classify ∘ delimit`, so comparing against `readout.classify()` compares `classify(x)` with
`classify(x)`, and it passed under the very mutation it was written to catch.

---

## Step B3 — byte-parity against production, the port anchor

**Status:** shipped (reviewed, fixes applied). **Milestone B complete.**
**Review:** [ng_alignment_b3_2026-07-23.md](../reviews/ng_alignment_b3_2026-07-23.md) — 3 Majors, all
applied.

### 1. Plan

On the shared fixture, every read whose production result is a `Delimited::Region` must measure the
**same bytes**; reads production called `BorderOffEnd` must land in a one-flank or no-flank case,
checked by count since that reclassification has no production oracle. *Source:* spec §10.3, arch
§Test & bench shape.

### 2. What shipped

[delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs) — a **permanent, `#[cfg(test)]`**
differential harness. A hand-written SplitMix64 over four fixed seeds generates loci (random motifs
of period 1–6, independently-drawn flank lengths **including zero**) and reads derived from them by
unit gains and losses, out-of-frame indels, substitutions and end truncations. Both implementations
run on identical input; ng's shipping code depends on nothing in production.

**Result: 12,000 cases per run, zero divergences** — and a soak at `PVC_PARITY_CASES=50000`
(200,000 cases, release) is likewise green.

### 3. The review's contribution: the harness could not fail in three ways

The reviewer's job was to ask whether the oracle is real, and answered it by mutating the aligner:
**5 of 7 mutations were caught, 2 were not.** The fixes:

- **A transposed anchor was invisible** — the failure the aligner's own docs call the whole risk of
  the widening. Invisible to *parity* by construction, since production's `BorderOffEnd` discarded
  which flank was missing. Fixed with a **self-consistency** check against ng's own offsets, which
  is not a tautology. **My first version of that check was itself wrong** and its own test caught it:
  the implication holds only when the flank exists, because "not anchored" has two causes.
- **Zero-length flanks were never generated**, which is exactly where B2's widening deliberately
  disagrees with production — and where the harness would have **panicked on an intentional
  divergence**. This falsified a claim in B2's own source comment, now corrected. The divergence is
  asserted rather than tripped over.
- **A fresh scratch per case** made the only size-varying driver blind to the grow-but-never-clear
  hazard that becomes live at Milestone C. Hoisted.

### 4. Validation

Container: fmt clean, clippy clean, `cargo test --lib` **2197 passed, 0 failed**. Transposition
mutation re-run: **fails both parity tests at the first generated case.**

### 5. The recorded gap, and its closure

At the time of the review **one mutation still slipped through**: losing the tract-aware gap-open in
the row-0 initialisation survived 200,000 cases. It was stated plainly rather than papered over.

**It closed itself.** Row 0's `gap_open` is consulted at **column 1 only** — from column 2 on the
match predecessor is `UNREACHABLE`, so the deletion-extend candidate always wins — and column 1 is
inside the tract only when the left flank is empty. So row 0's tract-awareness reduces to a single
case: a locus with **no left flank**. That is precisely what the generator did not produce, and
precisely what finding B3-2 made it produce. Re-running the mutation now fails within ~100 cases.

Two lessons, both already familiar from this plan: the gap was in the *fixture*, not the code; and
a coverage hole is worth understanding structurally rather than accepting, because the structure
said exactly which input would close it.

---

## Post-checkpoint follow-ups (owner-delegated)

After Checkpoint B the owner delegated the accumulated open decisions. Four landed as their own
commits, each validated in the container:

- **`Motif` lifted to `types.rs`** — closing the peer→stage back-reference that put
  `region_typing` (ng step 3) on `alignment`'s public surface, contradicting its own module doc.
  `segment_criteria` re-exports, so step 3's callers are untouched.
- **`FlatEmission::try_new`** — the check is now real in every build profile, not compiled out of
  release. Arch §3's `Result` ban is justified by per-*cell* cost, which does not reach a
  once-per-run constructor; arch §3's own escape clause names this remedy.
- **`ReadBases`** — the read/quality-length precondition is dissolved rather than documented. The
  three interchangeable `&[u8]` arguments are now two, and the invariant is true by construction.
  The parity oracle passing unchanged is the evidence the signature change altered no behaviour.
- **The architecture doc reconciled** with what shipped, and its `Motif`-allocation error removed.

Still open, and correctly so: the `ViterbiScratch` grow-without-clear note, which is **Milestone C's**
to handle — banding is what creates the hazard, so the fix belongs with it.
