# Code Review: ng_alignment_a3
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A3 of `doc/devel/ng/impl_plan/alignment_best_path.md` — the stutter model (and `RepeatContext`)
**Status:** Approve-with-changes (all Blocker/Major/Minor findings applied — see §Author response)

---

### 1. Scope

- **What was reviewed:** the uncommitted working-tree diff for A3 (A0–A2 committed as `bb8cf20`, `7692f43`, `3309e7d`). All three reviewers confirmed scope independently: `src/ng/alignment/mod.rs` modified (a small docs + re-export hunk), `src/ng/alignment/stutter.rs` new.
- **Categories dispatched (9):** reliability, errors, naming, idiomatic, smells, refactor_safety, module_structure, defaults, **extras** (numerical stability; the model's output feeds a shared likelihood). **Skipped:** `unsafe_concurrency`, `tooling`.

### 2. Verdict

**Approve-with-changes.** **1 Blocker, 5 Majors**, plus Minors. The step is almost entirely arithmetic, so reviewers were asked to work the formulas by hand rather than trust the tests — a test that shares a formula's mistake passes happily.

**The arithmetic itself came back correct.** One reviewer re-derived every branch against spec §5.2 and production's `stutter_pmf`: the in-frame and out-of-frame terms match term by term; `e = Δ − Δ/period` is right and symmetric for negatives because Rust truncates toward zero (Δ = −4, p = 3 → −3); the two-scale cutoff numbers check out; and `unsigned_abs()` throughout means `i64::MIN` is safe, where production's `-units - 1` form would have overflowed at `i32::MIN`. `MAX_SLIP`, `GEOM_MIN` and `GEOM_MAX` all equal production's, including production's reuse of `GEOM_MIN` as the `equal` floor.

**Every finding below is about what the tests could not see, or about a contract the code did not keep — not about a wrong formula.**

### 3. Execution status

Container, after fixes: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -- -D warnings` exit 0; `cargo test --lib` **2168 passed, 0 failed, 4 ignored**. The two pre-existing red commands are unchanged. No benchmark exists for this code and none was run; no performance number is claimed anywhere.

### 4. Open questions and assumptions

1. **Is the clamping question here the same as step A1's?** Explicitly asked, and the answer was **no**, with a distinction worth keeping: `FlatEmission::new` already had a loud `debug_assert!` *plus* a release clamp *plus* a documented trigger for becoming `try_new`, so A1's open question is only about the release half. `StutterModel::new` had **no check of any kind** — not a `Result`, not a clamp on four of its six arguments, not an assertion. And the clamp is defensible here in a way it is not for `MismatchFraction`, because it is **ported** behaviour: production does exactly this. So the finding was the missing assertion and the missing mass sanitizing, both of which need no `Result` and add nothing to `types.rs`. **Applied; A1's separate question still stands for the owner.**

### 5. Top 3 priorities

1. **B1** — the out-of-frame direction split was untested; swapping `out_up`/`out_down` passed the whole suite.
2. **M1** — `new` clamped two of six arguments while the type doc claimed "the clamps **are** the contract".
3. **M2** — six positional `f64`s with fixtures that made two separate transpositions invisible.

### 6. Findings

#### Blocker

**B1: the out-of-frame direction split was untested**
**Categories:** reliability. **Confidence:** High.
Both fixtures set `out_up == out_down == 0.01`, and **every directed out-of-frame assertion used a positive Δ**. So swapping the two out-of-frame direction masses passed the entire suite. In-frame direction *was* covered; out-of-frame was the hole. This matters because direction asymmetry is the model's defining structure — stutter is contraction-biased — and a swap here is silent.
**Fix applied:** an `all_distinct()` fixture whose six rates are all different (`out_up: 0.004`, `out_down: 0.012`, `in_geom: 0.95`, `out_geom: 0.8`), the out-of-frame formula test extended to negative Δ, and `a_contraction_outscores_an_expansion_of_the_same_size_in_both_regimes`.

#### Major

**M2: six positional `f64` arguments, with two transpositions invisible to every test**
**Categories:** naming. **Confidence:** High.
`new(in_up, in_down, in_geom, out_up, out_down, out_geom)`. The sharpest evidence was in the tests themselves: `shipped_defaults()` was `0.05, 0.05, 0.95, 0.01, 0.01, 0.95`, so an `in_up`↔`in_down` swap was invisible — and **every** fixture in the file used `0.95` for *both* geometrics, so an `in_geom`↔`out_geom` swap was invisible in all twelve tests. HipSTR keeps the two geometrics genuinely independent (spec §5.2), so the pairing is real.
**Fix applied:** a `StutterRates` named-field struct, plus the all-distinct fixture (which also closes B1).

**M1: `new` clamped two of six arguments while the doc claimed the clamps were the contract**
**Categories:** errors, defaults, reliability (convergent). **Confidence:** High — *verified by execution*.
The geometrics were clamped; the four direction masses were not. A reviewer ran the expressions: a negative mass gives `-0.0475`, a mass of `2.0` gives `1.9`, and `NaN` gives `NaN`. Two further holes the docs actively denied: **`f64::clamp` passes `NaN` through**, so "geometrics held strictly inside (0, 1)" was false; and **`f64::max` absorbs `NaN`**, so a `NaN` mass produced `equal() == 0.01` — the constructor reporting a healthy floor while every score was `NaN`.
**Fix applied:** a `debug_assert!` on all six rates (the real check), plus release-mode sanitizing — non-finite → `0` for a mass, `GEOM_MIN` for a geometric — and `equal` now `clamp`ed rather than `max`ed, so it cannot exceed 1 either. The doc now separates the two jobs explicitly: the geometric clamps and the `equal` floor are *ported model behaviour*; the mass sanitizing is a *backstop for a violated precondition*.

**M3: the "must stay equal to production's `MAX_SLIP`" invariant had no test**
**Categories:** errors, module_structure (convergent). Copying rather than importing is right (ng does not depend on production), but the doc asserted an equality that nothing enforced.
**Fix applied:** `the_copied_cutoff_still_equals_productions`, a **test-only** reference to production's constant — so shipping ng code still depends on nothing there.

**M4: `RepeatContext` was in the wrong file**
**Categories:** module_structure. **Confidence:** High.
Arch §2.2 puts it beside `RepeatGeometry`, and arch §Module home says `mod.rs` holds "the three traits + the shared types". The stated reason for putting it in `stutter.rs` — "it names a `StutterModel`" — **is not a Rust constraint**: `mod.rs` can simply `use` the type.
**Fix applied:** moved to `mod.rs`. (It still lands in step A3 rather than A2, which remains the recorded deviation — that part *is* forced, since the type it names arrives here.)

**M5: `ReferenceAlleleLen` was a dead alias that enforced nothing**
**Categories:** smells, module_structure (convergent). Publicly re-exported with zero consumers, and being a *transparent* alias it gave none of the reference-versus-candidate enforcement its doc claimed.
**Fix applied:** deleted. The distinction it gestured at is real and is documented on `StutterModel` itself; a type that would actually enforce it belongs with the adapting constructor, when that lands.

#### Minor (applied)

- **The `period − 1` error** (extras). The doc said out-of-frame changes are cut off "roughly `period` times sooner"; the arithmetic gives `period − 1`, so at period 2 the claim was off by 2× — the effect actually *vanishes* there. Corrected, with period 2 now pinned by a test. The neighbouring numbers (13 bp, 40 bp, 14 → 11) were all verified correct.
- **The cutoff was pinned only on the zero side**, so a `>` / `>=` slip survived (reliability). Now pinned on both sides, in both regimes and both directions.
- **Two near-duplicate regime branches** (smells) — collapsed into one `regime()` helper, which also makes the "size is at least one" precondition provable in one place instead of two.
- **`Copy` on a 56-byte struct** (idiomatic), buying nothing since the context holds it by reference. Dropped.
- **Tests recomputed expectations from the constants they were meant to pin** (refactor_safety) — the same trap A0 hit. The clamp assertions now use the contractual literals `0.01` / `0.99`.
- **The linear-versus-log boundary was never stated** (naming), though HipSTR's own fields are logs. Now stated in the module doc, with the `_ln` convention referenced.
- **Nothing mechanical stopped the two HipSTR parameter rows being mixed** (defaults) — a hazard spec §5.2 records the spec itself falling into. Now `hipstr_shipped()` and `hipstr_em_start()` named constructors, pinned by a test.

#### Accepted with reason

- **`NonZeroU8` for the period** — recorded by two reviewers as arch §3's `debug_assert!` rule done *better*: the illegal state is unrepresentable, so there is no release path to get wrong.
- **Seven accessors are the right price** for private fields: the `_private: ()` alternative leaves the `pub` fields writable, so a caller with a `mut` model could assign past the clamps.
- **No `Default`** — correct; a default stutter set would be a hidden behaviourally-significant value, and the two-parameter-set hazard makes it worse than usual.
- **`stutter.rs` placement and the absent `ssr_` prefix** — both correct per arch §Module home (a prefix should *disambiguate*, and there is no non-repeat stutter model to disambiguate from). No back-reference introduced.

### 7. Out of scope observations

The two pre-existing red validation commands. Also: three deferred decisions live only in rustdoc with no greppable tracker — they are carried in this plan's follow-up list and in `PROJECT_STATUS.md` instead.

### 8. Missing tests added now

`the_out_of_frame_branch_reproduces_the_published_formula_in_both_directions`, `a_contraction_outscores_an_expansion_of_the_same_size_in_both_regimes`, `ill_formed_rates_still_yield_probabilities` (all six slots × six ill-formed values), `the_two_hipstr_parameter_sets_are_kept_as_matched_rows`, `the_copied_cutoff_still_equals_productions`, plus negative-rank and both-sides-of-the-cutoff coverage.

### 9. What's good

- The re-indexing rationale states what it is **not** (double-counting) before what it is (rank compression), which is the misconception a reader most likely arrives with.
- `a_single_unit_slip_outweighs_a_larger_one` was confirmed to actually catch the inverted-geometric trap — and the reviewer established *why*: monotonicity alone would not catch it, since at geom = 0.05 the sequence still decreases. The quantitative ratio assertion is what does the work, and it is correct at 1/(1−0.95) = 20.
- The `equal`-floor contract is documented as a reason **not** to write an obvious test, which is a hard thing to record and easy to lose.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib
```

### Author response

The Blocker and all five Majors **fixed in this step's commit**, along with every Minor. Nothing carried forward from A3 itself; step A1's `try_new` question remains open for the owner and is distinct from this step's (open question 1).

Per-category audit trail at `tmp/review_2026-07-23_ng-alignment-a3/`.
