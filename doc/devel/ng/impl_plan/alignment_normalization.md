# ng alignment — normalization: implementation plan (3 of 3)

**Status:** draft, 2026-07-23. The build order for the **last third of the alignment module**: the
`AlignmentNormalizer` interface, its three implementations, the property test that grades them
against a definition, and the screen that says whether they differ at all. Design is settled in
[`../spec/alignment.md`](../spec/alignment.md) (§6, §10) and
[`../arch/alignment.md`](../arch/alignment.md). This turns that design into build order; it is **not**
a place for new design.

Follows [`alignment_best_path.md`](alignment_best_path.md) (plan 1) and
[`alignment_marginal.md`](alignment_marginal.md) (plan 2).

> **Why this runs last, and when to pull it forward.** Nothing depends on it: no other algorithm, no
> other plan, and no blocked work. The two plans before it unblock the STR locus generator and
> reproduce production's read likelihood, which is why they come first. But this plan is **short,
> self-contained, and answers a live question cheaply** — Milestone D's screen can close or confirm
> the hypothesis that indel *placement* explains part of production's indel deficit, for the cost of
> one run. If that question becomes urgent, this plan can be pulled ahead of plan 2 without
> disturbing anything; it needs only the `Alignment` type from plan 1.

---

## Scope

**In:** the `AlignmentNormalizer` trait; three implementations —
`left_align_structured.rs` (1a, production's), `left_align_repeated.rs` (1b, freebayes'), and the
fixpoint wrapper (1c); the leftmost **property checker**; the differ-at-all screen.

**Out (other owners):**

- **Which normalizer the generic read preparer calls** → read preparation's own plan. This module
  supplies algorithms and knows no caller (spec §1).
- **Whether normalization explains the indel gap** — Milestone D produces the evidence; acting on it
  is a calling-path decision. And note what the evidence can and cannot say: freebayes having **no
  stutter model at all** shows a stutter model is not *necessary* to beat production, but it does not
  identify production's defect. Normalization is one candidate; the reference-favouring genotype
  prior and the filters are others, and they live outside this module (spec §6).
- **Everything from the first two plans** — the aligners, the stutter model, `LogProb`.

## Principles (how the order was chosen)

- **The oracle before the implementations.** The leftmost property checker (Milestone A) is built and
  tested against hand-built alignments *before* any normalizer exists, so every implementation is
  graded from its first commit. This is the module's **only** comparison against a definition rather
  than against a rival implementation — it earns being built first, and it outranks agreement between
  the three (spec §6).
- **Reuse over rewrite.** 1a is a port: it calls the existing `left_align_indels` and supplies only
  the trait wrapper. Nothing about left-alignment is re-derived (arch §5).
- **Isolate the silently-wrong step.** A misplaced indel is a wrong variant, not a crash. 1a lands as
  **its own commit** with the property test green before and after, so `git bisect` can find it if
  placement moves.
- **Cheapest discriminating measurement first.** The screen (D) costs one run and can close the whole
  normalization hypothesis; it precedes any larger comparison (spec §6).
- **Incremental, with pauses.** One milestone, stop for review.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; native host build at the end.

## Preconditions (already in place)

- **Plan 1 is done** — `src/ng/alignment/` exists with its skeleton and the `Alignment` type. Plan 2
  is *not* a dependency; this plan uses nothing from it.
- **The design is settled**: spec §6 (the three algorithms, the property test, the screen, and what
  the freebayes comparison does and does not establish); arch §3 (the trait), §5 (the reconciliation
  rows), §Test & bench shape.
- **The reuse targets and the port oracle:** `left_align_indels`
  ([indel_norm.rs:408](../../../../src/pileup/walker/indel_norm.rs#L408)), `left_align_cigar`
  ([:254](../../../../src/pileup/walker/indel_norm.rs#L254)), `normalize_alleles`
  ([norm_seqs.rs:108](../../../../src/norm_seqs.rs#L108)).
- **freebayes is vendored read-only** for 1b's shape: `freebayes/src/LeftAlign.cpp:385`
  (`stablyLeftAlign`, cap 20 at `LeftAlign.h:118`). Implement from the description, not by
  transliteration (spec §8).

---

## The steps

### Milestone A — the interface and the oracle, before any implementation

**☐ A1. The `AlignmentNormalizer` trait.**
`normalize(&mut Alignment, read, reference)`. It takes the whole `Alignment` and not just its
operations because left-alignment can drop a leading deletion and move the alignment's start, which
an operations-only signature could not express. The read is required too: how far a gap may shift is
bounded by matching on **both** sequences. No implementations. *Depends:* plan 1. *Source:* arch §3,
§5.

**☐ A2. The leftmost property checker.**
A test-support function: given an `Alignment`, its read and its reference, assert that **no indel can
shift one base left and still represent the same pair of sequences**. Not an implementation of
left-alignment — a checker for one. Unit-tested against hand-built alignments that are
known-leftmost and known-not-leftmost, including an indel inside a homopolymer and one inside a
period-2 repeat, where shifting is possible and must be detected. *Depends:* A1. *Source:* spec §6,
arch §Test & bench shape.

> **Checkpoint A:** the checker accepts leftmost alignments and rejects shiftable ones, with no
> normalizer written yet. Every later step is graded by it. Pause for review.

### Milestone B — 1a, the port

**☐ B1. `left_align_structured.rs`. Own commit; do not bundle.**
An `AlignmentNormalizer` calling the existing `left_align_indels`. **Port its two load-bearing
behaviours, they are the point:** it merges consecutive indels, and it propagates an indel across
alignment blocks rather than stopping at the first — this is *not* a naive single pass, and the
failure it is commonly assumed to have is the case those two behaviours exist to handle (arch §5).
**Isolated because its failure is silent:** a misplaced indel is a wrong variant, not a panic.
Verified by A2 on a fixture, plus byte-parity of the resulting operations against calling
`left_align_indels` directly. *Depends:* A2. *Source:* spec §6, arch §5.

> **Checkpoint B:** 1a normalizes a fixture, the property test passes on its output, and its
> operations match the production function byte for byte. Pause for review.

### Milestone C — 1b and 1c

**☐ C1. `left_align_repeated.rs` — repeated simple passes.**
Re-run a simple shift pass until nothing moves, bounded (freebayes uses 20). **Do not copy
freebayes' handling of exhaustion:** it returns a failure flag its own caller ignores, so
non-convergent results ship silently (arch §5). Here, exhaustion is reported. Verified by A2.
*Depends:* B1. *Source:* spec §6, arch §5.

**☐ C2. 1c — the fixpoint wrapper.**
1a applied repeatedly until nothing moves, with the iteration cap as a safety net that **fails
loudly** rather than returning a half-normalised alignment. A thin wrapper over 1a, not a third
implementation of the shifting. Verified by A2, plus a test that a deliberately-capped run surfaces
the failure instead of swallowing it. *Depends:* B1. *Source:* spec §6.

> **Checkpoint C:** all three normalizers exist and all three pass the property test. Where any of
> them *fails* it, that is a finding to record — **the property test outranks agreement between
> them.** Pause for review.

### Milestone D — the screen

**☐ D1. The differ-at-all screen.**
Run all three over the same real reads and count outputs that disagree, reported per algorithm pair.
**This is the cheap discriminating measurement:** if the count is near zero, normalization cannot
explain any difference in calling and that avenue closes for one run (spec §6). Land it as a small
example binary alongside the module, following `examples/ng_ssr_loci_dump.rs`'s shape. *Depends:* C2.
*Source:* spec §6, §10.3.

> **Checkpoint D:** the screen runs on a real fixture and reports disagreement counts. **Plan 3 is
> complete, and with it the alignment module.** Pause for review, and record the screen's result — it
> is an input to whether the calling path adopts a different normalizer.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the property checker accepts hand-built leftmost alignments and rejects shiftable ones, in a homopolymer and a period-2 repeat |
| B | **the property test on 1a's output, plus byte-parity of its operations against `left_align_indels`** — the port anchor |
| C | the property test on 1b and 1c; 1c surfaces cap-exhaustion instead of swallowing it |
| D | the screen runs on real reads and reports per-pair disagreement counts |

## Out of scope (other plans)

- **Adopting a normalizer on the generic read-preparation path** — read preparation's plan; gated on
  D1's result.
- **The remaining open questions** the module leaves behind: the band's headroom, the default emission
  for the repeat-aware best-path aligner, whether an affine marginal is worth building, and what
  triggers re-alignment on the generic path (spec §9, and `../spec/read_preparation.md` §9). All are
  settled by measurement, and none blocks this plan.
