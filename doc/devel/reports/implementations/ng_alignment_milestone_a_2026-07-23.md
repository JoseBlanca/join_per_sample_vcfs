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

---

## Step A1 — `Emission` and its two implementations

**Status:** shipped (reviewed, fixes applied).
**Review:** [ng_alignment_a1_2026-07-23.md](../reviews/ng_alignment_a1_2026-07-23.md) — no Blockers, **5 Majors**, ~12 Minors, all applied.

### 1. Plan

The `Emission` trait plus the per-base-quality implementation (porting `EMISSION_LN` **including
its quality-zero floor**) and the flat-rate one. *Source:* spec §4.2, arch §2.3, §5.

### 2. Assumptions and decisions

- **The flat model's `insert_ln` is `ln(1/4)` — a decision, not a port**, as arch §2.3 warned it
  would have to be. Production's per-quality path uses `ln(1/4)`; its flat path has no emission for
  an inserted base at all, only a *transition* cost (`gap = eps`,
  [pair_hmm.rs](../../../../src/ssr/cohort/pair_hmm.rs)). Reasoning for choosing `ln(1/4)`: `ε` is
  the chance of misreading a base whose true identity the reference supplies, and an inserted base
  has no such base — so there is nothing for `ε` to describe, and both models make the same uniform
  assumption. Taking production's `gap = eps` instead would put a transition cost in an emission
  slot and **price the same event twice** once the aligner's gap model lands on top.
- **Arch §6's deferred signature question is now settled: quality resolves per *row*.** Arch listed
  "whether quality arrives per call or as a pre-resolved row" as an impl-time confirmation that
  "resolves when the first two implementations exist" — they now do. `scores_for(quality) ->
  BaseScores` is the trait's primary method, with `emit_ln` a provided convenience over it. The
  reason is structural, not measured (**nothing is built yet to measure**): a quality belongs to a
  *read base*, so it is constant along a matrix row while the reference base varies along it, and
  production's own loop hoists it per row. Matching that loop shape matters here beyond speed —
  Milestone B has to reproduce production byte for byte.
- **The bases are compared by raw byte equality, and canonicalizing them is the caller's job.**
  This is the arch §3 precondition shape, and it did not need a design decision: ng already offers
  both forms deliberately — `RefSeq::fetch` canonicalizes, `RawRefSeq` is verbatim with soft-mask
  intact. Upper-casing inside `emit_ln` was rejected as a scoring-model change smuggled into a
  component, and it would cost per cell.
- **`ε = 1` is the treatment for a non-finite rate** — the no-information end of the scale, so
  garbage input cannot yield confident output.
- **A `try_new` was recommended and deliberately not applied** — see Tradeoffs.

### 3. Changes made

| file | change |
|---|---|
| [src/ng/alignment/emission.rs](../../../../src/ng/alignment/emission.rs) | **new** — `Emission` (with `Sized`), `BaseScores`, `PerQualityEmission`, `FlatEmission` |
| [src/ng/alignment/mod.rs](../../../../src/ng/alignment/mod.rs) | `pub mod emission;` + re-exports + doc correction |

### 4. Tests added — 15

Highlights: `per_quality_table_is_bit_exact` (exact `f64` equality across all 256 qualities, guarding
the parity contract against a *reformulation* that a tolerance test would wave through);
`emission_is_finite_at_every_quality` and `flat_emission_stays_total_for_rates_outside_the_contract`
(the totality contract, the latter on the release path specifically);
`a_match_outscores_a_mismatch_above_the_quality_one_crossover` and its flat counterpart;
`mismatch_floor_never_binds_over_the_quality_domain`; `emission_compares_bases_by_raw_byte_equality`.

**One test wrote itself into a correction.** The first version asserted "a match always outscores a
mismatch" and **failed at Q1**. The test was wrong, not the code: a match beats a mismatch exactly
when `1 − ε > ε/3`, i.e. `ε < 0.75`, i.e. `Q > 1.249`, so at Q0 and Q1 the ported model genuinely
prefers the mismatch — at `ε ≈ 0.79` the base is nearly noise, and there are three ways to disagree
against one to agree. Inherited from production, not introduced here. The test now pins both sides
of the crossover instead of the tidy claim.

### 5. Validation results

Container, after fixes: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -D
warnings` exit 0; `cargo test --lib` **2141 passed, 0 failed, 4 ignored**. The two pre-existing red
commands are unchanged (see A0 §5 and Standing items).

### 6. Tradeoffs and follow-ups

- **The review's strongest argument was not applied, deliberately.** `FlatEmission::new` should
  arguably be a checked `try_new`: arch §3 bans `Result` on *hot-path* cost grounds, but this rate
  is per-run configuration consumed once into two `f64` fields, so checking it costs nothing per
  cell — and arch §3's own escape clause names the checked constructor on the *context type* as the
  remedy. Two things held it back: it needs a new `DomainError` variant in `src/ng/types.rs`, which
  this plan's preconditions say it "adds nothing to"; and arch's clause is conditioned on the value
  being reachable from **untrusted input**, which it is not today (the only callers are tests).
  Instead the totality hole is closed by clamping, and **the question goes to the owner at
  Checkpoint A.** It becomes mandatory the moment a CLI flag or config field feeds this rate.
- **Two further recorded-not-applied items**, both for the same reason (they touch `types.rs` or
  arch §2.3): a `Base` newtype for the two `u8` base parameters, which are positionally adjacent and
  transposable; and renaming `insert_ln`.
- **The hot-path claim is structural and unmeasured.** No bench exists for this module; the first
  one belongs with the first algorithm, not with a component nothing calls yet.
- **`from_unchecked_rate` exists for testability** — `new` minus the debug assertion, because the
  assertion fires in the test profile and the release path would otherwise be untestable without
  compiling the test out of the build anyone runs.

---

## Step A2 — `RepeatSpan`, `RepeatGeometry`, and the `BestPathAligner` trait

**Status:** shipped (reviewed, fixes applied).
**Review:** [ng_alignment_a2_2026-07-23.md](../reviews/ng_alignment_a2_2026-07-23.md) — no Blockers, **8 Majors**, several Minors; six Majors applied, two carried to Checkpoint A.

### 1. Plan and the one deviation

The plan's A2 is "`RepeatSpan`, `RepeatGeometry`, `RepeatContext`, and the `BestPathAligner` trait".
**`RepeatContext` moved to A3, recorded.** It bundles a `&StutterModel`, and that type does not exist
until A3 — so A2 could only have landed it by inventing a placeholder type or by shipping a struct
that gains its second field one commit later. Moving one type across the boundary between two
adjacent types-only steps preserves both steps' intent and keeps each commit coherent. *Source:*
arch §2.1, §2.2, §3.

### 2. Assumptions and decisions

- **`RepeatSpan`'s four variants are the widening**, and the reason is worth restating: production's
  `Delimited` collapses every unanchored case into one **side-blind** `BorderOffEnd`, so it cannot
  say which flank was missing. B2 is the step that produces these, and the plan flags it as
  silently-failing precisely because a mis-assigned side is a wrong observation class on a read that
  still looks fine.
- **The accessors are split by what they can honestly claim** — see the review's M1. `measured_length`
  answers only where a measurement exists; `length_lower_bound` answers everywhere; `observed_span`
  hands back bases, not a length.
- **`Range<u64>`, not ng's existing `RegionSpan`** — `RegionSpan` is the region seam's
  *reference*-coordinate vocabulary and this module knows no coordinates. Costs `Copy`.
- **`Bp` for lengths, bare `u64` for offsets.** The distinction is real: `Bp` is documented as ng's
  length currency, while an offset into a caller-supplied slice is the `RepeatInterval` shape.
- **A factual error in the architecture doc, not repeated in the code.** Arch §2.2 justifies the
  per-call geometry partly on "`Motif` is an allocation". ng's `Motif` is `Copy` over an inline
  `[u8; 6]` and allocates nothing. The decision stands on its other ground — the geometry changes at
  every locus — and the code now says only that. **Arch §2.2 still carries the claim.**

### 3. Changes made

[src/ng/alignment/mod.rs](../../../../src/ng/alignment/mod.rs) — `RepeatSpan` (+4 accessors),
`RepeatGeometry` (+2 methods), `BestPathAligner`, module-doc update.

### 4. Tests added — 12

The ones that carry weight: `only_a_two_flank_span_yields_a_measured_length` (which also asserts the
trap directly — the two variants' raw spans are *equal*, so nothing reading the span alone can tell
them apart); `fitting_is_exactly_having_a_repeat_length` (sweeps every flank pair against every
reference length in a small window, and asserts the three parts reconstruct the whole);
`repeat_len_reports_no_answer_rather_than_underflowing` (both subtractions, including the mirror case
that is the only route to the second). `TODO(Milestone B)` records the four tests the trait contract
owes once an implementation exists.

### 5. Validation results

Container: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -D warnings` exit
0; `cargo test --lib` **2153 passed, 0 failed, 4 ignored**.

### 6. Tradeoffs and follow-ups — three for the owner at Checkpoint A

- **Lift `Motif` to `types.rs`?** The import is a genuine peer→stage back-reference: `region_typing`
  is ng step 3, and `alignment` states two lines above the import that it knows no steps. The `pub
  motif` field puts it on ng's public surface. Recorded as debt in code rather than fixed, because
  the move reaches beyond this step and touches `types.rs`.
- **A `ReadBases` newtype in `align`'s signature?** `read`, `quality` and `reference` are three
  positionally interchangeable `&[u8]`, and the read/quality-length precondition exists exactly
  because two of them must correspond. Bundling read and quality dissolves that precondition with no
  `Result` — the check moves to a once-per-read constructor, off the hot path. Not applied: it
  changes the trait the arch doc specifies and ripples into every later aligner. **Cheapest now,
  while there are no implementors.**
- **Arch §2.2's `Motif`-allocation claim** should be corrected in the doc.
