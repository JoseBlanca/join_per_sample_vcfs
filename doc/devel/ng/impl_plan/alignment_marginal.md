# ng alignment — marginal aligners: implementation plan (2 of 3)

**Status:** draft, 2026-07-23. The build order for the **aligners that return a summed probability
rather than a line-up**: `LogProb`, the `MarginalAligner` interface, the sequence-versus-sequence
marginal that production already uses, and the whole-read forward it is compared against. Design is
settled in [`../spec/alignment.md`](../spec/alignment.md) (§5, §7, §9, §10) and
[`../arch/alignment.md`](../arch/alignment.md). This turns that design into build order; it is **not**
a place for new design.

Follows [`alignment_best_path.md`](alignment_best_path.md) (plan 1); followed by
[`alignment_normalization.md`](alignment_normalization.md) (plan 3).

> **Read this before starting: the module owns less of the marginal than it looks.** Production's
> read-likelihood model computes *P(observed repeat | candidate allele)* in two steps — how likely
> the length change is (stutter), times how well the letters match once the candidate has been
> stretched to the observed length. **Only the second step is alignment**, and only it is built here.
> The stutter half belongs to the genotyping, which composes this module (spec §5.1). Porting the
> whole of `HipstrModel` into this plan would put a genetics model inside a module that knows only
> about alignment.

---

## Scope

**In:** `LogProb` landed in `types.rs`; the `MarginalAligner` trait; **algorithm 5** (one sequence
against another, gaps confined to the ends — production's `align_subst`, returned in logarithms);
**algorithm 6** (a forward over the whole read against flank-repeat-flank, built with fixed synthetic
qualities).

**Out (later plan / other owners):**

- **`AlignmentNormalizer` and the three normalizers** → plan 3,
  [`alignment_normalization.md`](alignment_normalization.md).
- **The stutter distribution's *use* in a read likelihood** — the genotyping's. This plan builds no
  stutter term; `StutterModel` already exists from plan 1 and is not needed here (spec §5.1).
- **The 5-versus-6 bake-off.** The two are not head-to-head rivals — 5 scores a measured repeat
  against a candidate, 6 scores a whole read — so the comparison that matters spans this module *and*
  the genotyping that composes it. Its harness lives with the genotyping; this plan proves only that
  each algorithm computes what it claims (spec §10.3).
- **Carrying more read sequence, or per-read qualities, to genotyping time** — the storage change
  algorithm 6 would imply if it won. Deliberately not paid for here (spec §9).
- **A forward over a repeat-loop graph** — out of the initial set, with its trigger recorded
  (spec §10.1).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **Reuse over rewrite.** Algorithm 5 is a port of one function, `align_subst`. Its restrictions are
  ported with it, not tidied away (arch §5).
- **Isolate the silently-wrong steps.** Two things here fail without a crash: the linear-to-logarithm
  conversion at the port boundary (a transposition returns a plausible wrong number), and dropping
  the interior-gap restriction (which destroys identifiability and changes nothing visible). Each is
  its own commit with its test.
- **Test the path the source never exercises.** Production's caller always makes the two sequences
  equal-length, so the unequal-length forward is **unreachable from it**. A test built by imitating
  that caller would never touch the forward pass at all; it must be driven deliberately (spec §5.1).
- **Cheapest honest comparison.** Algorithm 6 is built with **fixed synthetic qualities**, which
  settles whether it wins before anyone pays for the per-read storage it would otherwise require
  (spec §9).
- **Incremental, with pauses.** One milestone, stop for review.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; native host build at the end.

## Preconditions (already in place)

- **Plan 1 is done** — `src/ng/alignment/` exists with its skeleton, and the `Emission` component is
  available (algorithm 5's flat-rate scoring reuses it).
- **The design is settled**: spec §5.1 (the two length cases, the interior-gap restriction, the
  normalisation defect, the boundary against the genotyping), §7 (the logarithm contract, the
  per-call context), §9 (algorithm 6 and its storage consequence); arch §1, §3, §5.
- **`LogProb` is already specified** — declared in
  [`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) §1 and assigned to `types.rs` by
  [`../arch/module_layout.md`](../arch/module_layout.md), but **not yet written**: `grep` over `src/`
  finds none. This plan is its first user, not its author (arch §1).
- **The reuse target:** `align_subst` ([pair_hmm.rs:52](../../../../src/ssr/cohort/pair_hmm.rs#L52)),
  its `banded_forward` ([:83](../../../../src/ssr/cohort/pair_hmm.rs#L83)), its `HmmScratch`
  ([:32](../../../../src/ssr/cohort/pair_hmm.rs#L32)) and its `FLANK_SLOP`
  ([:28](../../../../src/ssr/cohort/pair_hmm.rs#L28)).

---

## The steps

### Milestone A — `LogProb` and the interface (no logic)

**✅ A1. Land `LogProb` in `types.rs`.**
`LogProb(pub f64)` — a probability held as its natural logarithm, where negative infinity means
impossible and is a **value, not an error**. Land it where the shared vocabulary already put it,
beside `Bp`. Its point is that the compiler refuses to mix it with an ordinary probability, which is
the transposition this module is most exposed to, since the code being ported returns linear
probabilities. *Depends:* plan 1. *Source:* arch §1, §5.

**✅ A2. The `MarginalAligner` trait.**
`Scratch` and `Context` associated types; `marginal_probability(read, reference, context, scratch)`
returning `LogProb`. `Context` is `()` for a sequence-versus-sequence marginal and the repeat context
for one that needs the geometry. No implementations. *Depends:* A1. *Source:* arch §3.

> **Checkpoint A:** `LogProb` exists in the shared vocabulary and the trait compiles with no
> implementations. Pause for review.

### Milestone B — algorithm 5, the sequence marginal

**✅ B1. The equal-length case.**
When the two sequences are the same length there is exactly one way to line them up, so the sum has a
single term and the result is a base-by-base comparison under a **single flat error rate** — not
per-base qualities. Tests: an exact match scores `(1−ε)^len`; one substitution costs exactly a factor
`(ε/3)/(1−ε)`. *Depends:* A2. *Source:* spec §5.1, arch §5.

**✅ B2. The unequal-length forward, with interior gaps forbidden. Own commit; do not bundle.**
The banded forward summing the few ways a length difference can be absorbed by gaps **within a couple
of bases of either end**. **The interior-gap restriction is load-bearing and must be ported, not
tidied away:** an indel in the middle of the repeat is what the *stutter* model explains, and letting
this algorithm explain it too makes base error and slippage indistinguishable. **Isolated because its
failure is silent:** a general forward runs fine and destroys identifiability with no symptom. Test
that an interior indel scores far below an end gap of the same size — the assertion that pins the
restriction. **Drive this path deliberately:** production's own caller never reaches it (see
Principles). *Depends:* B1. *Source:* spec §5.1, arch §5.

**✅ B3. The logarithm boundary. Own commit; do not bundle.**
Convert at the interface: the ported function returns a plain probability, this trait returns its
logarithm. **Isolated because its failure is silent:** a missed or doubled conversion returns a
plausible number, not an error. Tests: parity against the ported function, compared in whichever
space is convenient; and an unreachable line-up returning negative infinity rather than zero, which
is the distinction the contract exists for. Record the known **normalisation defect** — over
differing-length outputs the total mass slightly exceeds one — in the implementation's doc comment,
so it is reproduced knowingly. *Depends:* B2. *Source:* spec §5.1, §7, arch §5.

> **Checkpoint B:** algorithm 5 matches the ported function on equal-length inputs, sums correctly
> over end gaps on unequal-length ones, refuses interior gaps, and returns logarithms. Pause for
> review.

### Milestone C — algorithm 6, the whole-read forward

**✅ C1. The whole-read forward.**
A forward over the whole read against a full flank-repeat-flank sequence, scoring with **one flat
error rate** and **fixed synthetic qualities** — following production, which pairs a flat rate at
scoring with a base-quality gate at read ingestion. Building it this way is deliberate: it settles
whether the design wins **before** anyone pays for the per-read storage it would otherwise require,
and it keeps the per-locus observation table a tally rather than per-read records (spec §9).
*Depends:* B3. *Source:* spec §9, §10.1.

**☐ C2. Each algorithm computes what it claims.**
Not a bake-off — that spans this module and the genotyping (see Scope). Prove instead, on synthetic
reads with known truth, that each returns the probability it says it returns: a read generated from
an allele scores higher against that allele than against its neighbours, and the ordering is stable
under the flat error rate. *Depends:* C1. *Source:* spec §10.3.

> **Checkpoint C:** both marginals run on synthetic truth and each scores its own generating allele
> highest. **Plan 2 is complete.** Pause for review, and hand the 5-versus-6 comparison to the
> genotyping, which owns the harness that can run it.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | `LogProb` compiles in `types.rs`; the trait compiles with no implementations |
| B | exact-match and single-substitution closed forms; **an interior indel scores far below an end gap** (the identifiability guard); parity against `align_subst`; an impossible line-up returns −∞, not 0 |
| C | synthetic reads with known truth — each algorithm scores its own generating allele highest, stably |

## Out of scope (next plan)

- **`AlignmentNormalizer` and the three normalizers** → plan 3,
  [`alignment_normalization.md`](alignment_normalization.md).
- **The 5-versus-6 comparison** → the genotyping's harness; it spans two modules (spec §10.3).
- **The storage change algorithm 6 implies**, if it wins — read-preparation and locus-generation
  concerns, deliberately unpaid here (spec §9).
