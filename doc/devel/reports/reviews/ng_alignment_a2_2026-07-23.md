# Code Review: ng_alignment_a2
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A2 of `doc/devel/ng/impl_plan/alignment_best_path.md` — `RepeatSpan`, `RepeatGeometry`, `BestPathAligner`
**Status:** Approve-with-changes (all applied except two recorded debts — see §Author response)

---

### 1. Scope

- **What was reviewed:** the uncommitted working-tree diff for A2 (A0/A1 committed as `bb8cf20`, `66327c1`). All three reviewers independently confirmed scope via `git diff`: one modified file, `src/ng/alignment/mod.rs`, +273/−3.
- **In-scope file:** [mod.rs](../../../../src/ng/alignment/mod.rs) — the A2 additions only.
- **Out of scope:** the committed `Alignment` and `emission.rs`; production `src/ssr/`, `src/pileup/`; steps A3, B1–B3.
- **Categories dispatched (8):** reliability, errors, naming, idiomatic, smells, refactor_safety, module_structure, defaults. **Skipped:** `unsafe_concurrency` (nothing applicable), `tooling` (no manifest change), `extras` (types only — no algorithm, no untrusted input, nothing to measure).

### 2. Verdict

**Approve-with-changes.** No Blockers. **8 Majors**, spread across five categories — the highest count of the run so far, and appropriate for a step whose entire job is to shape distinctions that later code cannot recover if the shape is wrong.

### 3. Execution status

Container, after fixes: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -- -D warnings` exit 0; `cargo test --lib` **2153 passed, 0 failed, 4 ignored** (12 new tests). The two pre-existing red commands (`--all-targets` bench panic; `cargo doc`) are unchanged and unrelated — see Standing items in `PROJECT_STATUS.md`.

Findings labelled "Needs verification": 0. One reviewer explicitly noted it ran no build commands and did not restate the supplied output as its own.

### 4. Open questions and assumptions

1. **Should `Motif` move to `src/ng/types.rs`?** Affects M6. This is the fix for a genuine peer→stage back-reference, but it reaches beyond A2 and touches `types.rs`, which the plan says this plan does not. **Recorded as debt in code; raised at Checkpoint A.**
2. **The architecture doc contains a factual error.** Arch §2.2 justifies the per-call geometry partly on "`Motif` is an allocation". Verified false for ng's `Motif` — `Copy` over an inline `[u8; 6]`, never heap-allocated ([segment_criteria.rs:147](../../../../src/ng/region_typing/segment_criteria.rs#L147)). The *decision* stands on its other ground (the geometry changes at every locus). Two reviewers found this independently. **The code no longer repeats the claim; arch §2.2 still carries it** — raised at Checkpoint A, since this skill does not edit design docs.

### 5. Top 3 priorities

1. **M1** — `span()` made the wrong answer a one-liner: it returned a range for all three anchored variants, so a length taken off it is a fabricated short allele for a truncated read.
2. **M2** — `is_measurement()`'s `matches!` hid an implicit `_ => false`, the one place a wildcard would absorb exactly the change that matters.
3. **M8** — `Scratch: Default` let an implementation choose a result-changing value that no call site can see or report.

### 6. Findings

#### Major

**M1: `span()` makes the type's central distinction easy to lose**
**Categories:** errors (filed Major), smells (convergent). **Confidence:** High.
All three anchored variants carry an identical span, and `span()` returned it for all three. So `span().map(|s| s.end - s.start)` yields a plausible allele length for a read that **ran off its own end** — a short allele that was never observed. The correct version differs from the wrong one by an `if`, the wrong one is shorter, and it leaves nothing a grep can find. `is_measurement()` existing is not enough: it has to be *remembered*.
**Fix applied:** the accessor is split by what it can honestly claim. `measured_length()` returns `Some` only for `Between`; `length_lower_bound(read_len)` is defined for every case (including `Unanchored`, where a read lying wholly inside the repeat still proves a bound); and the raw span survives as `observed_span()`, documented for **extracting bases** — which is legitimate for a lower bound, and is what an interrupted repeat needs, since its sequence must come out verbatim. A test asserts the trap directly: `Between(4..10).observed_span() == FromLeft(4..10).observed_span()`.

**M2: `is_measurement()` used `matches!`, absorbing the one change that matters**
**Categories:** refactor_safety. **Confidence:** High.
`matches!(self, Self::Between(_))` compiles an implicit `_ => false`. It was the only consumer in the file a fifth variant would *not* break — `observed_span()` has no wildcard and would fail loudly, and the tests enumerate variants by name. A future measurement-bearing variant (an interrupted-repeat case) would silently read as "lower bound" → **systematically short alleles**, with the compiler saying nothing.
**Fix applied:** exhaustive `match`, with the reasoning recorded at the site.

**M3: nothing pinned `fits_reference` and `repeat_len` to the same predicate**
**Categories:** reliability, errors, smells (three-way convergent). **Confidence:** High.
They computed one predicate by two different arithmetic routes (`checked_add` + `<=` versus two chained `checked_sub`s). They agreed — one reviewer worked the proof — but no test would have caught `<=`→`<` or `checked_sub`→`saturating_sub`.
**Fix applied:** `fits_reference` is now *defined* as `repeat_len(..).is_some()`, so they cannot drift; and `fitting_is_exactly_having_a_repeat_length` sweeps every flank pair and reference length in a small window, additionally asserting `left + right + repeat == reference_len` whenever it fits.

**M4: the second `checked_sub` was never exercised**
**Categories:** reliability. **Confidence:** High.
The arithmetic *was* complete, but coverage was not symmetric: `geometry(u64::MAX, 1)` short-circuits at the **first** `?`. The mirror case — left flank fits, right flank overruns at full width — was never tested, so a regression in the second subtraction would have passed.
**Fix applied:** the mirror case, plus a small exact-fit pair.

**M5: the trait contract's testable properties were owed to nobody**
**Categories:** reliability. **Confidence:** High.
The contract states four properties no test can reach until an implementation exists (the tie-break order, the 5′ junction rule, empty reference → `Unanchored`, and the debug-asserted preconditions). A0 set the precedent with an explicit `TODO(Milestone E)`; A2 carried no equivalent.
**Fix applied:** `TODO(Milestone B)` naming all four owed tests.

**M6: `Motif` is imported from a pipeline-stage module into a module that says it knows none**
**Categories:** module_structure. **Confidence:** High.
`region_typing` calls itself ng step 3; `alignment`'s own header says two lines above the import that it is not a pipeline step and knows none. The `pub motif` field puts that back-reference on ng's public surface. `Motif` is STR domain vocabulary with consumers in three modules across stage boundaries, and `module_layout.md` principle 3 already assigns such types to the shared vocabulary.
**Not fixed — recorded as debt in code.** The fix (lift to `types.rs`, re-export from `segment_criteria` so step 3 is untouched) reaches beyond this step and touches `types.rs`. See open question 1.

**M7: `align` takes three positionally interchangeable `&[u8]`**
**Categories:** idiomatic. **Confidence:** Medium.
`read`, `quality` and `reference` are mutually transposable at every call site, and the read/quality-length precondition exists precisely because two of them must correspond. A `ReadBases` newtype bundling read and quality would dissolve that precondition **without introducing a `Result`** — the check moves to a once-per-read constructor, off the hot path.
**Not applied.** It changes the trait signature the arch doc specifies (§3) and would ripple into every Milestone B/D/E implementation; the plan puts the trait here and the aligners there. **Raised at Checkpoint A** — it is cheapest now, while there are no implementors.

**M8: `Scratch: Default` lets an implementation choose the answer invisibly**
**Categories:** defaults. **Confidence:** High.
The bound is the only way a caller can build a scratch, and the otherwise-thorough contract said nothing about what `Default` may decide. For a **banded** aligner — Milestone C, the next thing after the parity port — the natural scratch contents include a band half-width or a matrix cap, and a too-narrow band **silently loses long alleles**. That would be a result-changing value no call site can see, choose, or report in a bake-off. The reviewer also noted `Default` takes no `&self`, so the bound has *already* decided grow-on-demand for every implementation without writing it down.
**Fix applied:** the bound stays (it is right for what scratch is); the contract now states that `Default` must decide nothing that changes a result — scratch is buffers only, empty, grown on demand, inert with respect to output — and names the escape hatch (`fn scratch(&self) -> Self::Scratch`) for an aligner that ever needs a configured one.

#### Minor (applied)

- **`fits_reference`/`repeat_len` took bare `u64` while the struct's own fields are `Bp`** (naming). Using `Bp` adds nothing to `types.rs`, so the plan constraint held. Fixed. The reviewer correctly distinguished these *lengths* from the *offsets* (`reference_offset`, `RepeatSpan`'s payload) that stay bare.
- **5′/3′ vocabulary inverts for reverse-strand reads** (naming). The `FromLeft`/`FromRight` docs explained slice-end facts in strand terms, but reads are held reference-forward, so the words would mean the opposite of the fact for half the data. Rewritten in terms of read start and end.
- **No producer invariant on the span payload** (errors). An inverted range reads as a zero-length tract — the same short-allele failure the enum was widened to prevent. Fixed: documented as the producer's invariant, and both length accessors saturate so a violation cannot manufacture an enormous length. Tested (the literal `10..4` is rejected by clippy's `reversed_empty_ranges`, which is itself evidence the invariant is one the language's tooling treats as a mistake).
- **Uncovered degenerate cases** (reliability): every `RepeatSpan` test used the same `4..10`; the empty reference had no zero-flank geometry test. Both added.

#### Argued and accepted as-is

- **`#[non_exhaustive]` correctly omitted** (refactor_safety): it constrains only *other* crates, so it would be inert here — and it would force the very `_` arm that M2 is about.
- **The enum beats a product type** (smells, pressed deliberately): a "span plus two anchored booleans" encoding makes the `(false, false)` row change whether the payload *exists*, so it needs `Option<Span>` and yields 8 representable states for 4 meaningful ones.
- **`Range<u64>` over ng's existing `RegionSpan`** (idiomatic, Nit): `RegionSpan` is the region seam's *reference*-coordinate vocabulary, and this module knows no coordinates. The `Copy` loss is real; the decision is now recorded rather than merely made.
- `RepeatGeometry` field addition breaks every literal (no `Default`, no constructor, all-`pub` fields, no `..`); trait method addition breaks implementors loudly (no default bodies).

### 7. Out of scope observations

The two pre-existing red validation commands, and arch §2.2's `Motif`-allocation error (open question 2).

### 8. Missing tests added now

`only_a_two_flank_span_yields_a_measured_length`, `every_case_bounds_the_length_below`, `span_lengths_handle_empty_inverted_and_extreme_payloads`, `an_empty_reference_fits_only_a_geometry_with_no_flanks`, `fitting_is_exactly_having_a_repeat_length`, plus the mirror underflow case. `TODO(Milestone B)` records the four the trait contract owes.

### 9. What's good

- The four-variant enum is the right shape and was defended on its merits when challenged, not by appeal to the arch doc.
- Every doc comment on the new types states the *failure mode*, not just the rule — "counts a truncated read as a short allele" is what makes the distinction stick.
- `repeat_len` returning `None` rather than a saturated zero, because zero is indistinguishable from a genuinely empty tract.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib
```

### Author response

Six of eight Majors and every Minor **fixed in this step's commit**. Two carried to Checkpoint A as owner decisions, both because they reach beyond the step: **M6** (lift `Motif` to `types.rs`) and **M7** (a `ReadBases` newtype in the trait signature). Also carried: arch §2.2's factual error about `Motif` allocating.

Per-category audit trail at `tmp/review_2026-07-23_ng-alignment-a2/`.
