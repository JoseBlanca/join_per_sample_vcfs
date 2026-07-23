# ng alignment — Milestone A: module skeleton and the aligner types (no logic)

**Date:** 2026-07-23
**Plan:** [alignment_best_path.md](../../ng/impl_plan/alignment_best_path.md) (plan 1 of 3), Milestone A
**Design authority:** spec [alignment.md](../../ng/spec/alignment.md), arch [alignment.md](../../ng/arch/alignment.md)
**Process:** plan-driven-implementation — one implement → review → apply-fixes → commit loop per step.

This report accumulates across Milestone A's four steps (A0–A3). Each step's own code review is
saved separately under [`../reviews/`](../reviews/).

---

## Step A0 — scaffold `src/ng/alignment/` and define `Alignment`

**Status:** shipped (reviewed, fixes applied).
**Review:** [ng_alignment_a0_2026-07-23.md](../reviews/ng_alignment_a0_2026-07-23.md) — no Blockers, no Majors, 7 Minors, all applied.

### 1. Plan

Create the module folder and the one type the plan assigns to this step; no algorithms.
*Source:* arch §Module home, §2.1.

- `src/ng/alignment/mod.rs` — a **folder, not a file**, because the module holds competing
  implementations that get compared (five aligners and three normalizers land in sibling files).
- `Alignment { reference_offset, cigar }`, reusing production's `CigarOp` rather than minting a
  parallel operation type.
- Wire `pub mod alignment;` into `src/ng/mod.rs`.

**A0 was run as its own loop, not merged with A1.** The plan permits sharing an iteration when a
step would leave a stage near-empty; that latitude was not needed — A0 lands a public type with a
real contract, and the review found seven Minor findings on it, so the stage was not empty.

### 2. Assumptions and decisions

- **`CigarOp` reuse costs production nothing — checked, not assumed.** The arch doc mandates reuse,
  but ng's own precedent cuts both ways: `region_typing` had to *port* `Motif` because production's
  is `pub(crate)` and would have leaked from a `pub` ng signature. Verified that `CigarOp`
  ([walker/mod.rs:43](../../../../src/pileup/walker/mod.rs#L43)) is fully `pub`, so that objection
  does not apply here.
- **`reference_offset` is `u64`, not the arch sketch's `usize`.** ng speaks one width for
  coordinates and lengths; the closest analogue, `RepeatInterval.start` — also a 0-based offset into
  the slice its producer was handed — is `u64` for that reason
  ([tandem_repeat.rs:198-212](../../../../src/ng/tandem_repeat.rs#L198-L212)). Raised by the review
  at High confidence and verified directly.
- **The field is `cigar`, not the arch sketch's `ops`.** Every other site in the crate qualifies the
  name (`cigar` / `cigar_ops`, seven sites). Licensed by arch's own preamble: "Signatures are
  illustrative; the **contract** is the deliverable."
- **Which `CigarOp` variants an affine aligner emits is NOT decided here.** The enum carries nine,
  including three (`Skip`, `HardClip`, `Padding`) that describe an input record rather than a
  computed placement. The producer is Milestone E, which is *gated*, and arch §6 still carries
  `OPEN:` on the affine aligner's output shape. The field records the obligation on the producing
  step rather than guessing an answer — declining the review's proposed concrete subset.
- **No `try_new` and no `Result`.** Consistent with arch §3 (nothing in this module returns a
  `Result`) and with the project's newtype convention: `Alignment`'s well-formedness is uncheckable
  from the struct alone, and `AlignmentNormalizer::normalize` mutates in place, so a
  constructor-enforced invariant could be broken without passing through the constructor.

### 3. Changes made

| file | change |
|---|---|
| [src/ng/alignment/mod.rs](../../../../src/ng/alignment/mod.rs) | **new** — module doc (folder rationale, "not a pipeline step") and `Alignment` |
| [src/ng/mod.rs](../../../../src/ng/mod.rs) | `pub mod alignment;` + the landed-modules doc sentence (repunctuated to a semicolon list; it had grown to three coordinating "and"s) |

### 4. Tests added

- `alignment_distinguishes_offset_and_operation_order` — a deliberate **compile-time anchor**, since
  A0 lands no algorithm. It asserts whole-struct inequality against a changed `reference_offset` and
  against permuted operations, so it fails if a future hand-written `PartialEq` ignored either field
  — which field-by-field assertions on a struct literal could not catch.
- `TODO(Milestone E)` records the two tests the producing step owes:
  `align_returns_offset_relative_to_the_reference_stretch` (the stretch-versus-genome confusion is a
  wrong answer, not a panic) and `align_returns_operations_in_read_order`.

### 5. Validation results

All in the container via `./scripts/dev.sh`, after the review fixes:

- `cargo fmt --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `2126 passed; 0 failed; 4 ignored`, including the new test.
- `cargo doc --no-deps` — red, **pre-existing**; the three intra-doc links this step adds
  (`Bp`, `RepeatInterval`, `Motif`) all resolve, and the failing set is unchanged by this diff.

**Two pre-existing failures, verified by `git stash -u` and re-running on the clean tree:**
`benches/psp_writer_perf.rs:386` panics under `cargo test --all-targets --all-features`, and
`cargo doc --no-deps` fails on 11 unresolved intra-doc links across `src/ng/` and `src/ssr/`.
Neither is caused by this work; both make a project-standard validation command red for every step
of this plan. Flagged at Checkpoint A.

### 6. Tradeoffs and follow-ups

- **Two divergences from the arch doc's §2.1 sketch** (`ops` → `cigar`, `usize` → `u64`). Both are
  inside the latitude arch grants its own signatures, and both are recorded on the fields. **The
  owner may want to fold them into arch §2.1** — raised at Checkpoint A; this skill does not edit
  design docs.
- **`CigarOp`'s home is a pre-existing misplacement.** It is crate-wide CIGAR vocabulary living
  inside a pipeline stage, and A0 puts that on ng's *public* surface — `ng::alignment::Alignment`
  cannot be named without naming `pileup::walker`. Lifting it to `src/cigar.rs` is a production
  edit and production is frozen, so it is recorded as known debt at the reuse site; the port-back of
  this module is the natural moment.
- **`Hash` is not derivable** on `Alignment` because `CigarOp` does not derive it. Recorded in code
  so the next person finds the blocker rather than the symptom.
