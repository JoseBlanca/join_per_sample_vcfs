# Code Review: ng_alignment_a0
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A0 of `doc/devel/ng/impl_plan/alignment_best_path.md` — the `src/ng/alignment/` scaffold and the `Alignment` type
**Status:** Approve-with-changes (all applied — see §Author response)

---

### 1. Scope

- **What was reviewed:** the uncommitted working-tree diff for plan step A0.
- **Reviewed against:** branch `ng-locus-evidence`, working tree at the time of review (parent commit `f81c56b`).
- **In-scope files:**
  - [mod.rs](../../../../src/ng/alignment/mod.rs) — new
  - [mod.rs](../../../../src/ng/mod.rs) — module declaration + one doc sentence
- **Deliberately out of scope:** production `src/ssr/` and `src/pileup/` (read-only oracles — ng never edits production, owner 2026-07-16); later plan steps A1–A3 and B1–B3, which add the `Emission` trait, the repeat-aware types, the stutter model and the algorithms. A0 is types-first with no algorithms **by plan design**, so the absence of logic, error paths and callers is intent, not a gap — every sub-agent was told so explicitly.
- **Categories dispatched (8):** `reliability` (always), `errors` (always), `naming` (always), `idiomatic` (always), `refactor_safety` (always), `smells` (always), `module_structure` (scope spans two files, and the folder-vs-file decision is the step's main structural claim), `defaults` (scope adds public API).
- **Categories skipped:** `unsafe_concurrency` (no `unsafe`, no `Arc`/`Mutex`/atomics/channels/`async`/threads in scope), `tooling` (no `Cargo.toml` change; scope is a diff, not a crate), `extras` (no parser, validator, security boundary, untrusted input or hot-path code — A0 contains no executable logic at all).

### 2. Verdict

**Approve-with-changes.** No Blocker and no Major findings. Seven Minor findings, of which two are real code changes and five are contract-documentation gaps that are cheapest to close before the first producer and consumer exist.

### 3. Execution status

All commands run in the project container via `./scripts/dev.sh`.

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib` | `test result: ok. 2126 passed; 0 failed; 4 ignored; 0 measured; 0 filtered out; finished in 39.70s` |
| `cargo test --all-targets --all-features` | **one failure, pre-existing** — `thread 'main' panicked at benches/psp_writer_perf.rs:386:60`. Verified by `git stash -u` and re-running on the clean tree: identical failure without this diff. |
| `cargo doc --no-deps` | **fails, pre-existing** — 11 unresolved intra-doc links, in `src/ng/locus_generation/ssr.rs`, `src/ng/region_typing/{segment_criteria.rs,mod.rs}`, `src/ng/tandem_repeat.rs` and several `src/ssr/` files. **None in the in-scope files**; re-checked after the fixes, and the three intra-doc links this step adds (`Bp`, `RepeatInterval`, `Motif`) all resolve. |
| `cargo audit` | not run — no dependency change in scope. |

Findings labelled "Needs verification": 0. The one Medium-confidence finding (Mi3) is explicitly about a step that does not exist yet, and was resolved by recording the obligation rather than guessing the answer.

### 4. Open questions and assumptions

1. **Does renaming a field away from the architecture doc's sketch need owner sign-off?** Affects Mi1. Resolved in the code: [arch/alignment.md:9](../../ng/arch/alignment.md) states "Signatures are illustrative; the **contract** is the deliverable", which places a field rename inside implementer latitude. The arch doc's §2.1 sketch still reads `ops`; **the owner may want to fold the rename into the doc** — raised at Checkpoint A.
2. **Which `CigarOp` variants can the affine aligner emit?** Affects Mi3. Genuinely undecided — the producer is Milestone E, which is *gated*, and arch §6 still carries `OPEN:` on the affine aligner's output shape. Not answered; the obligation is recorded on the field instead.
3. **Is `usize` or `u64` right for a slice offset in ng?** Affects Mi2. Resolved against project precedent, verified directly: [tandem_repeat.rs:198-212](../../../../src/ng/tandem_repeat.rs#L198-L212) types `RepeatInterval.start` — the same shape, a 0-based offset into the slice its producer was handed — as `u64`, with the comment "`u64` is the width regardless".

### 5. Top 3 priorities

1. **Mi2** — `reference_offset: usize` violates ng's documented one-width rule; zero callers today makes this the cheapest it will ever be to fix.
2. **Mi1** — `ops` is the crate's only unqualified name for a `Vec<CigarOp>`; the divergence hardens as the folder fills with aligners.
3. **Mi5/Mi6** — the type's contract (may `cigar` be empty? do the operation lengths sum to the read?) exists only in the author's head; getting it wrong downstream produces wrong likelihoods, not a panic.

### 6. Findings

#### Minor

**Mi1: src/ng/alignment/mod.rs — `ops` is the crate's odd name for a `Vec<CigarOp>`**
**Categories:** naming. **Confidence:** High (inconsistency), Medium (that changing it is right).
Every other binding of this concept in the crate qualifies it: `cigar_ops` at `src/pileup/per_sample/record_specs.rs:30`; `cigar` at `src/pileup/walker/mod.rs:246`, `src/pileup/walker/indel_norm.rs:37`, `src/bam/alignment_input.rs:86`; and as a parameter in `baq_engine.rs`, `indel_norm.rs:408`, `src/ng/read/filtering.rs:1023`. `Alignment.ops` was the single site dropping the qualifier. The aligner's output will meet production's CIGAR-speaking code — `left_align_indels(cigar: &mut Vec<CigarOp>, …)` is named in arch §5 as a reuse site — so a differing name costs a mental translation at each boundary.
**Fix:** rename to `cigar`.

**Mi2: src/ng/alignment/mod.rs — `reference_offset: usize` diverges from ng's one-width rule**
**Categories:** idiomatic (Minor, High), naming (Nit — newtype trigger), errors (Nit — `usize`/`u32` mismatch against `CigarOp`'s payloads).
ng fixed one width for coordinates and lengths: `Bp(pub u64)` is documented as `u64` "since B2 (spec §4): ng speaks one width, so nothing narrows, nothing is checked, and no off-by-width bug is possible" ([types.rs:142-151](../../../../src/ng/types.rs#L142-L151)), and the closest analogue `RepeatInterval.start` is `u64` for exactly that reason. A `usize` here forces a cast at every boundary where it meets a `Bp`, a `RepeatInterval` or a genome position.
**Fix:** `pub reference_offset: u64`.
*Convergent note (errors, Nit):* if Milestone E ever converts an `Alignment` back to a production CIGAR, flag the `u64 → u32` direction — an `as` truncation there is a wrong placement, not a crash.

**Mi3: src/ng/alignment/mod.rs — `Vec<CigarOp>` admits operations an affine aligner can never emit**
**Categories:** idiomatic (filed), smells (cross-category note). **Confidence:** Medium.
`CigarOp` ([walker/mod.rs:43-53](../../../../src/pileup/walker/mod.rs#L43-L53)) has nine variants, including `Skip`, `HardClip` and `Padding` — production's *parser* vocabulary, describing an input record rather than a placement this module computed. Reuse is right per arch §5, but nothing states which subset is inhabited, so every future consumer must match all nine arms or guess.
**Fix:** documentation. **Note the reviewer's proposed doc text asserted the subset** (`Match`/`Insertion`/`Deletion`, plus `SoftClip` when local); that answer is not available yet — see open question 2 — so the applied fix records the *obligation* on the producing step instead of inventing the answer.

**Mi4: src/ng/alignment/mod.rs — the sole test is tautological**
**Categories:** reliability (Minor), smells (Nit), refactor_safety (Nit).
The test built a struct literal and asserted the public fields equal the literals just written; with `pub` plain-data fields and no constructor those assertions are true by language construction. Its real power is compile-time only. Its doc comment also claimed it pinned that the offset "is stored verbatim rather than folded into the first operation" — no code exists that could fold it. `refactor_safety` adds that field-by-field assertions would let a future third field be forced into the literal but never into an assertion.
**Fix:** whole-struct `assert_ne!` comparisons that fail if a hand-written `PartialEq` ignored either field, plus an honest doc comment naming the test a compile-time anchor.

**Mi5: src/ng/alignment/mod.rs — two documented invariants have no test and no owner**
**Categories:** reliability. The stretch-relative (not genomic) offset and the read-order operations are untestable at A0 by design, but nothing recorded the debt — and the stretch-versus-genome confusion is precisely a wrong-answer-without-panicking failure mode.
**Fix:** a `TODO(Milestone E)` naming the two tests the producing step owes.

**Mi6: src/ng/alignment/mod.rs — the struct's relational invariants live only in the author's head**
**Categories:** smells. **Confidence:** High.
The docs said what each field *is* and never what a consumer may rely on: whether an empty `cigar` is legal, and whether the read-consuming lengths must sum to the read's length. With `pub` fields and no `try_new`, nothing enforces either.
**Fix:** state the contract on the type, and say who owns it. Cross-checked against `errors`, which independently concluded `Alignment` is **not** owed a `try_new`: its well-formedness is uncheckable from the struct alone, and `AlignmentNormalizer::normalize` mutates in place, so a constructor-enforced invariant could be broken without passing through the constructor — the same case the project already settled for `GenomeRegion` with public fields.

**Mi7: src/ng/alignment/mod.rs:17 — `CigarOp` is a cross-stage interchange type imported from a pipeline-stage module**
**Categories:** module_structure. **Confidence:** High.
`use crate::pileup::walker::CigarOp` is a back-reference (algorithm module → pipeline-stage module), and `pub struct Alignment` puts it on ng's public surface: `crate::ng::alignment::Alignment` cannot be named without naming `crate::pileup::walker`, which tells every future reader the alignment module has something to do with the pileup walker — false, per spec §1. The reuse decision itself is sound; the misplacement is one level up and pre-existing (`CigarOp` is consumed by `bam/`, `ssr/pileup/`, `pop_var_caller/cli/`, two benches, two examples and now ng).
**Fix:** documentation only — the lift to a top-level `src/cigar.rs` is a **production edit and production is frozen**. Recorded as known debt at the reuse site.

#### Nits

- `Hash` is not derived. **Not applicable, and now recorded in the code:** `CigarOp` does not derive `Hash`, so `#[derive(Hash)]` on `Alignment` would not compile; adding it to `CigarOp` is a production edit.
- `pub struct Alignment` has no caller outside its module — accepted, matching every other ng module surface (`types`, `ref_seq`, `tandem_repeat`), which are `pub` because ng is a research crate consumed by examples and benches.
- Absence of `#[non_exhaustive]` is the **right** call and was noted as such by `refactor_safety`: arch §6 anticipates the type widening, and omitting it means the compiler flags every literal when it does.
- Prose: the module doc said "the traits live here in `mod.rs`" in the present tense while A0 lands no trait (raised by `idiomatic` and `module_structure`); and the `src/ng/mod.rs` landed-modules sentence had grown to a four-item list joined by three coordinating "and"s (raised by `naming`, `idiomatic`, `module_structure` — the most convergent finding in the review). Both fixed.

#### Clean categories

`defaults` — `No findings.` (no `Default` impl, no constructor, no `Option` parameter, no literal, no configuration value). `errors` — no Blocker/Major/Minor; confirmed the type does not foreclose the module's no-`Result` policy, since an empty `cigar` at offset 0 is constructible, so "nothing lined up" is expressible as a value. `refactor_safety` — no Blocker/Major/Minor; the single struct literal spells both fields with no `..`, all four traits are derived, no slice indexing, no `..Default::default()`, no partial destructure.

### 7. Out of scope observations

- **`benches/psp_writer_perf.rs:386` panics** under `cargo test --all-targets --all-features`. Pre-existing and verified independent of this diff. Follow-up: separate fix; it makes the project's own stated validation command red for every step of this plan.
- **`cargo doc --no-deps` is red** on 11 pre-existing unresolved intra-doc links across `src/ng/` and `src/ssr/`. Same follow-up class: it is one of the review skill's five verification commands and currently cannot pass.
- **`CigarOp`'s home.** See Mi7 — the lift out of `pileup::walker` is blocked on production being frozen; the port-back of the alignment module is the natural moment.

### 8. Missing tests to add now

None at A0 — the step lands no executable behaviour, and `reliability` was explicit that property/fuzz, concurrency, error-variant and doc-example tests are all not-applicable here rather than omitted. The tests the step *owes* are recorded as `TODO(Milestone E)` in the test module: `align_returns_offset_relative_to_the_reference_stretch` and `align_returns_operations_in_read_order`.

### 9. What's good

- The folder-over-flat-file decision states its reason **in the code** and cites the arch section that enumerates the ten sibling files — which is what earns the "imminent submodules" exemption on that rule's own terms, rather than by assertion (`module_structure` checked this and filed no finding).
- The `CigarOp` reuse justification names the *counter*-precedent it does not fall under — ng's `Motif`, ported because production's is `pub(crate)` — so a reader learns why this case differs instead of wondering whether it was considered.
- "It is not a pipeline step, and it knows none" is stated in the module doc with the spec cross-reference, so the module's most unusual structural property is the first thing a reader meets.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib
```

### Author response

All seven Minor findings and both prose nits **fixed in this step's commit**, before it was made.

| finding | resolution |
|---|---|
| Mi1 | fixed — field renamed `ops` → `cigar`, with the reason and the arch-doc divergence recorded on the field. Raised at Checkpoint A for the owner to fold into arch §2.1. |
| Mi2 | fixed — `reference_offset` is now `u64`, citing `Bp` and `RepeatInterval`. |
| Mi3 | fixed **as an obligation, not an answer** — the field documents that the affine aligner inhabits only part of the enum and that Milestone E's producing step must state the subset. The reviewer's concrete subset was declined as premature (open question 2). |
| Mi4 | fixed — the test now asserts whole-struct inequality on both a changed offset and permuted operations, and its doc comment names it a compile-time anchor. Renamed `alignment_distinguishes_offset_and_operation_order`. |
| Mi5 | fixed — `TODO(Milestone E)` naming the two owed tests. |
| Mi6 | fixed — contract documented; no `try_new`, per the convergent `errors` reasoning. |
| Mi7 | fixed (documentation) — the back-reference and the deferred lift to `src/cigar.rs` are recorded at the reuse site. The arch §5 half of the reviewer's suggested fix is **not** applied: this skill does not edit design docs. Raised at Checkpoint A. |
| Nits | `Hash` non-derivability recorded in code; prose fixed in both files. |

Per-category audit trail left in place at `tmp/review_2026-07-23_ng-alignment-a0/`.
