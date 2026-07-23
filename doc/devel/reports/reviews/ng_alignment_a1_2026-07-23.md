# Code Review: ng_alignment_a1
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step A1 of `doc/devel/ng/impl_plan/alignment_best_path.md` — the `Emission` component and its two implementations
**Status:** Approve-with-changes (all Majors and Minors applied — see §Author response)

---

### 1. Scope

- **What was reviewed:** the uncommitted working-tree diff for plan step A1 (A0 was already committed as `bb8cf20`).
- **In-scope files:** [emission.rs](../../../../src/ng/alignment/emission.rs) (new — the whole step); [mod.rs](../../../../src/ng/alignment/mod.rs) (module declaration, re-exports, one doc sentence).
- **Deliberately out of scope:** production `src/ssr/` and `src/pileup/` (read-only oracles); later steps A2, A3, B1–B3.
- **Categories dispatched (9):** reliability, errors, naming, idiomatic, refactor_safety, smells, module_structure, defaults, **extras** (newly applicable — this code is on the module's hottest path and its output feeds a byte-identity guarantee).
- **Skipped:** `unsafe_concurrency` (no `unsafe`/`Arc`/`Mutex`/atomics/channels/`async`/threads), `tooling` (no `Cargo.toml` change).

### 2. Verdict

**Approve-with-changes.** No Blockers. **Five Majors**, three of which are genuine defects rather than documentation gaps — and two of those were found by *executing* the code, not by reading it.

### 3. Execution status

Container, via `./scripts/dev.sh`, after fixes:

| command | result |
|---|---|
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib` | `test result: ok. 2141 passed; 0 failed; 4 ignored` (15 new tests) |
| `cargo test --all-targets --all-features` | one **pre-existing** failure (`benches/psp_writer_perf.rs:386`), verified on a clean stashed tree |
| `cargo doc --no-deps` | **pre-existing** red (11 unresolved intra-doc links elsewhere); none in the in-scope files |

Findings labelled "Needs verification": 0. Both hot-path claims are labelled structural, **not measured** — no benchmark exists for this code and none was run.

### 4. Open questions and assumptions

1. **Should `FlatEmission::new` become a checked `try_new`?** Affects M1. The `errors` reviewer argued yes and made the strongest case in the review: arch §3's ban on `Result` is justified by *hot-path* cost, but this rate is per-run configuration consumed once into two `f64` fields, so a checked constructor costs nothing per cell — and arch §3's own escape clause ("it becomes a checked constructor on the *context* type") describes this position exactly. **Not applied**, because it needs a new `DomainError` variant in `src/ng/types.rs` and this plan's preconditions state it "adds nothing to it". Arch's clause is also conditioned on the value being *reachable from untrusted input*, which it is not today — the only callers are tests. **Raised at Checkpoint A as the owner's call**; it becomes mandatory the moment a CLI flag or config field feeds this.
2. **Should `read_base`/`reference_base` be a `Base` newtype?** Affects a naming Minor. Same blocker (it would live in `types.rs`) and the same disposition — recorded, not applied.
3. **Should `insert_ln` be renamed** (`inserted_base_ln`)? Declined to limit divergence from arch §2.3's names, given the trait is already being reshaped for M2. Recorded.

### 5. Top 3 priorities

1. **M1** — `FlatEmission` breaks its own totality contract in release for out-of-contract rates; `ε = ∞` yields a **`+∞` mismatch score**.
2. **M2** — `emit_ln` bills per *cell* what production bills per *row*, in the module's hottest loop.
3. **M3** — tolerance assertions (`< 1e-12`) guarding a table that must stay **bit**-identical to production.

### 6. Findings

#### Major

**M1: emission.rs — `FlatEmission` violates the totality contract once the debug assertion is compiled out**
**Categories:** reliability, errors (convergent). **Confidence:** High — *verified by execution*, not inference.
The floor is a **lower** clamp only, so an out-of-range rate escapes upward. The `errors` reviewer ran it: `ε = -0.5` gives a match score of **`+0.405`** — a log-probability above 1; `ε = +inf` gives a **`+inf` mismatch score**; `NaN` is swallowed silently because `f64::max` returns the non-NaN operand, yielding two finite, ordinary-looking scores. A `+inf` is worse than the `-inf` the floor exists to prevent: it makes the event it scores *infinitely preferred*, so the model is not skewed but inverted, with nothing in the output to show for it. `debug_assert!` does not help — it compiles out of the release build this project runs.
**Fix applied:** clamp the rate into `[0, 1]` (non-finite → `1`, the no-information end) before computing, so totality holds unconditionally. The assertion stays as the real check in development.

**M2: emission.rs — the signature bills per matrix cell what production bills per row**
**Categories:** extras. **Confidence:** High (structural). **Not measured** — nothing is built yet to measure.
`emit_ln` took a quality per call, so `PerQualityEmission` deref'd a `LazyLock` (an initialisation check plus an atomic acquire the compiler need not hoist) and loaded a 16-byte pair **on every cell**. Production hoists both once per row and leaves a compare-and-select in the inner loop ([alignment.rs](../../../../src/ssr/pileup/alignment.rs)). A quality belongs to a *read base*, so it is constant along a matrix row while the reference base varies along it — the per-cell lookup is work the structure of the problem does not require. This is squarely arch §6's flagged **impl-time confirmation** ("whether quality arrives per call or as a pre-resolved row… resolves when the first two implementations exist"), and they now exist.
**Fix applied:** `scores_for(quality) -> BaseScores` is now the trait's primary method; `emit_ln` remains as a provided convenience over it. `row_resolved_scores_agree_with_the_per_call_form` pins that the two agree at every quality.

**M3: emission.rs — tolerance assertions guarding a bit-identity contract**
**Categories:** extras. **Confidence:** High.
Every check on the ported table used `< 1e-12` — about a million ulps on a log-space value near −6.9. The repeat-aware aligner built on this table must reproduce production's measured repeats **byte for byte** (spec §10.3), so reformulating the arithmetic (`ln_1p`, `powi`, splitting the division) would pass the suite while moving every downstream score.
**Fix applied:** `per_quality_table_is_bit_exact` asserts exact `f64` equality across all 256 qualities, against an expectation re-derived from the published model (so it stays portable, with no hardcoded platform bits). The hand-checkable decimal cases keep a tolerance and now say why.

**M4: emission.rs — nothing tested `FlatEmission`'s match-versus-mismatch ordering**
**Categories:** reliability. **Confidence:** High.
Swapping the two structurally identical lines that build `match_ln` and `mismatch_ln` passed the entire suite. This matters more for the flat model than the quality table, because its `ε` is a **constructor parameter**, so the `ε = 0.75` crossover is reachable by configuration rather than fixed at Q0/Q1.
**Fix applied:** `flat_emission_orders_match_above_mismatch_below_the_crossover`, asserting both sides of the crossover across six rates.

**M5: emission.rs — base comparison is raw byte equality with no stated precondition**
**Categories:** reliability. **Confidence:** High.
`b'a' != b'A'`, so a **soft-masked reference scores as a mismatch at every base** — and soft-masking marks repeats, which is exactly where the repeat-aware aligner works. Separately, `N` against `N` scores a full-confidence *match*. Neither was documented or tested.
**Resolution:** this is a caller precondition, not a bug, and did not need a design decision: ng deliberately offers both shapes — `RefSeq::fetch` canonicalizes (upper-cases ACGT, folds the rest to `N`, [ref_seq.rs:676](../../../../src/ng/ref_seq.rs#L676)) while `RawRefSeq` returns bases verbatim with soft-mask intact ([raw_chrom_reader.rs:179](../../../../src/ng/raw_chrom_reader.rs#L179)). Which one feeds the aligner is step 2's choice. **Fix applied:** documented as a precondition on the trait, in the shape arch §3 prescribes, plus `emission_compares_bases_by_raw_byte_equality` pinning both behaviours. Silently upper-casing inside the hot path was rejected — it would be a scoring-model change smuggled into a component, and it would cost per cell.

#### Minor (all applied)

- **Port fidelity, raised independently by three categories** (naming, smells, extras): the table floors the *mismatch* term where production floors only the match ([alignment.rs:64](../../../../src/ssr/pileup/alignment.rs#L64)), under a doc line claiming a faithful port. Verified inert — the smallest `ε/3` a `u8` quality yields is ≈`1.05e-26`, some 282 orders of magnitude above `f64::MIN_POSITIVE` — so the tables are bit-identical. **Fixed:** the divergence and its inertness are now documented, and `mismatch_floor_never_binds_over_the_quality_domain` asserts the claim over the whole domain rather than assuming it.
- **`insert_ln` paid a `LazyLock` atomic deref per call** (smells) — inconsistent with the module's own per-cell-cost argument. **Fixed:** `UNIFORM_BASE_LN` is a `const`, with `uniform_base_ln_is_ln_of_a_quarter` pinning the literal against `0.25f64.ln()`.
- **`#[derive(Default)]` on `PerQualityEmission`** (refactor_safety) — the one construction path that would survive the struct gaining a field, silently zero-filling it, where `new()` would fail to compile. **Fixed:** written out by hand, calling `new()`.
- **The floor tests recomputed their expectation from `PROBABILITY_FLOOR` itself** (refactor_safety) — editing the constant moved both sides of the assertion, and the only independent guard (`< -700.0`) passes for any floor below ~`1e-304`. **Fixed:** the constant is now checked against its literal value.
- **`Emission` was dyn-compatible** (refactor_safety), so arch §4's "never `dyn`" was documentation only. **Fixed:** a `Sized` supertrait makes it a compile error, at zero cost.
- **Undocumented `pub fn new`** (idiomatic), **missing `#[must_use]`** (idiomatic), **no `# Panics` section for the `debug_assert!`** (extras), **provenance cited by bare identifier** across two different production files (extras), **`&f64` array over `Copy` scalars** (idiomatic), **`emission` naming both a model and an `f64`** (naming), **`mod.rs` still said "the traits will live here in `mod.rs`"** (module_structure) — all fixed.
- **`FlatEmission` correctly has no `Default`** (defaults) — a default error rate would be exactly the hidden behaviourally-significant value the rules forbid. Its absence is a decision; now documented as one.

#### Clean

`module_structure` — one Nit only; layout matches arch §Module home on every checked point (named sibling file, no `ssr_` prefix — correctly, since `emission.rs` is shared by both algorithm families — no back-references, knows no callers, and A1 resolves A0's `mod.rs`-alone-in-its-directory smell). `defaults` — no Blocker/Major. The `_ln` suffix audit came back **clean**: every log/linear crossing is an explicit `.ln()`, and nothing linear escapes an `_ln` name, which was the transposition hazard worth checking given production's counterparts return linear probabilities.

### 7. Out of scope observations

- The two pre-existing red validation commands (see A0's review and Standing items in `PROJECT_STATUS.md`).
- No `benches/` entry guards the emission hot path. Arch's test-and-bench shape calls for a `bench/` in this module; it belongs with the first algorithm (Milestone B/D), not with a component that nothing calls yet.

### 8. Missing tests added now

`per_quality_table_is_bit_exact`, `mismatch_floor_never_binds_over_the_quality_domain`, `flat_emission_stays_total_for_rates_outside_the_contract`, `flat_emission_orders_match_above_mismatch_below_the_crossover`, `emission_compares_bases_by_raw_byte_equality`, `uniform_base_ln_is_ln_of_a_quarter`, `row_resolved_scores_agree_with_the_per_call_form`.

One implementation change was needed to make a test possible: `FlatEmission::from_unchecked_rate` is `new` **minus** the debug assertion — that is, exactly what a release build runs. Without it the totality test could not reach the path it exists to check, because the assertion fires in the test profile; the alternative, `#[cfg(not(debug_assertions))]`, would compile the test out of the build anyone actually runs.

### 9. What's good

- `insert_ln`'s `ln(1/4)` is documented at all three levels a reader might arrive from, and states the **rejected** alternative (production's `gap = eps` transition, which would double-price the event once a gap model lands) rather than only the choice — called "exemplary" by the `errors` reviewer.
- The Q0 floor's doc explains the *failure mode* (`-inf` annihilates every path through the cell), not just the mechanism, so the next reader knows why removing it would be catastrophic rather than merely different.
- Totality is discharged **exhaustively** over all 256 quality bytes rather than by sampling — the domain is small enough that there is no reason to sample.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib
```

### Author response

All five Majors and every Minor **fixed in this step's commit**, before it was made. Three items are deliberately **not** applied and are carried to Checkpoint A as owner decisions: the `try_new`/`DomainError` question (open question 1), the `Base` newtype (2), and the `insert_ln` rename (3) — all three would touch `src/ng/types.rs` or arch §2.3, which this loop does not edit.

Per-category audit trail at `tmp/review_2026-07-23_ng-alignment-a1/`.
