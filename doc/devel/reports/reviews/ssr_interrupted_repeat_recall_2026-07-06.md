# Code Review: ssr-interrupted-repeat-recall
**Date:** 2026-07-06
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** the `ssr-interruptions` branch feature diff (`0771462..HEAD`, 11 commits) — Phase 1 sequence-keyed SSR alleles + Phase 2 purity→stutter-level covariate
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** a feature-branch diff (11 commits, `0771462..HEAD`), 12 changed files under `src/ssr/cohort/` plus a new `examples/ssr_slip_dump.rs` (+1764 / −129 lines).
- **Reviewed against:** branch `ssr-interruptions` at HEAD, base `0771462`.
- **In-scope files:** [rung_ladder.rs](../../../../src/ssr/cohort/rung_ladder.rs), [candidate_set.rs](../../../../src/ssr/cohort/candidate_set.rs), [em_init.rs](../../../../src/ssr/cohort/em_init.rs), [attribution.rs](../../../../src/ssr/cohort/attribution.rs), [em.rs](../../../../src/ssr/cohort/em.rs), [vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs), [driver.rs](../../../../src/ssr/cohort/driver.rs), [prepass.rs](../../../../src/ssr/cohort/prepass.rs), [param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs), [inbreeding.rs](../../../../src/ssr/cohort/inbreeding.rs) (import + one field), [bakeoff.rs](../../../../src/ssr/cohort/bakeoff.rs) (import + one field), [examples/ssr_slip_dump.rs](../../../../examples/ssr_slip_dump.rs).
- **Deliberately out of scope:** `src/var_calling/`, `src/vcf/`, `src/ssr/cohort/read_model/`, and unchanged `src/ssr/cohort/` modules — these carry PRE-EXISTING clippy/test/doc failures unrelated to this feature (see §7).
- **Categories dispatched (11):** reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, unsafe_concurrency (determinism focus), smells, tooling, extras (stable-output + hot-path + diff-matches-intent). Each ran as a parallel sub-agent; audit trail in `tmp/review_2026-07-06_ssr-interrupted-repeat/`.

## 2. Verdict

**Approve-with-changes.** No Blocker. One Major (a bounded, silent statistical shift + hot-path regression sharing one root cause and one fix). The remainder are Minors — mostly documentation, test-coverage, and tuning-knob formalism. The two hard contracts — cross-thread byte-identity and panic-freedom — hold up under adversarial tracing, and the diff matches the spec with the seven "already conformant" functions verified untouched.

## 3. Execution status

Run in the project container (`./scripts/dev.sh`):
- `cargo fmt --check` — in-scope files **clean** (the only diffs are pre-existing, in `src/paralog/` and the untracked `examples/ssr_psp_seqdump.rs`).
- `cargo clippy --all-targets --all-features -- -D warnings` — **FAILS**, exit 101, on exactly one PRE-EXISTING lint `src/vcf/writer.rs:384` (needless_lifetimes, SNP path). In-scope code is clippy-clean.
- `cargo test --lib` — `1592 passed; 1 failed; 3 ignored`. The single failure `var_calling::posterior_engine::tests::larger_ref_pseudocount_cannot_increase_p_alt` is deterministic and PRE-EXISTING (SNP path, untouched). All 206 `ssr::cohort` tests pass. (3 ignored includes the intentional P2.0 diagnostic.)
- `cargo doc --no-deps` — **FAILS** only on pre-existing links in `read_model/mod.rs` + `var_calling/diversity.rs`. No in-scope file has a doc error.
- `cargo audit` — **not run** (cargo-audit not installed).
- **"Needs verification" findings:** 0 — every finding cites a location the reviewing agent read.

## 4. Open questions and assumptions

1. **Is the composition tie-break for equidistant *different-length* alleles intended?** (affects **M1**.) The spec claimed the slip attribution is "effect-neutral" only for same-length ties; the code also lets `align_subst` decide different-length equidistant ties. Confirm whether that is a wanted improvement (then doc + test it) or an accidental divergence (then restrict it).
2. **Are the Phase-2 purity-fit gates (`MIN_PURITY_FIT_READS`, `PURITY_FACTOR_FLOOR`) meant to be swept in F2?** (affects **Mi7**, **Mi9**.) They are hardcoded consts while the sibling §5.2 knobs are on `CandidateCfg` for sweeping. Decide "tunable → lift to a `PurityFitCfg`" vs "frozen → say so".
3. **Does `interruption_count`'s first-unit-phase precondition actually hold for every candidate/peak tract?** (affects **Mi1**.) It assumes the tract's first repeat unit carries the motif phase; nothing in-scope enforces it.

## 5. Top 3 priorities

1. **M1** — `nearest_called_by_sequence` runs `align_subst` on every length-nearest allele: a per-slip-read `banded_forward` DP added to the whole-cohort θ_locus refit (perf), and for composition-asymmetric *equidistant different-length* reads it can flip the slip `Δ` sign vs the old `nearest_parent` (silent θ_locus shift, untested). One fix addresses both.
2. **Mi2** — `allele_balance` silently returns `1.0` (FP defence off) when `allele_support` is empty, untested; add a `debug_assert!` + regression test so a future data-wiring bug fails loud instead of dropping the filter.
3. **Mi1** — `interruption_count` over-scores a first-unit interruption (→ collapses that allele's modelled stutter toward 0); guard the phase precondition or pin it with a test. Latent while Phase 2 is inert, live once a cohort activates the factor.

## 6. Findings

### Major

**M1: src/ssr/cohort/attribution.rs:56 — `nearest_called_by_sequence` scores composition on the length-nearest allele unconditionally: a hot-path DP regression, and a slip-sign flip on equidistant different-length ties**
**Categories:** refactor_safety, extras (convergent) · **Confidence:** High (perf), Medium (correctness)

Two facets of one root cause — the `filter(|(_, &(_, units))| (units−read_units).abs() == min_dist).map(align_subst).reduce(...)` shape:

- **Hot path (High):** when a single allele is length-nearest (every pure locus, every length-separated het), `align_subst` still runs once; for a *slip* read (`obs.len() != variant.len()`) it falls through to `banded_forward` ([pair_hmm.rs:59](../../../../src/ssr/cohort/pair_hmm.rs#L59)), an `O(m·n)` DP with a scratch resize. The old `nearest_parent` used integer length distance. `attribute_locus` runs this inside the refit loop for *every admitted locus*, so every stutter read now pays a DP where it paid a compare — including the pure-locus majority the feature is meant to leave unchanged.
- **Correctness (Medium):** the `== min_dist` filter also admits two *different-length* alleles equidistant from the read (het `(6,8)` with a read at 7 units). The old `nearest_parent` broke that tie to the lowest index (`nearest_parent(8,&[6,10]) == Some((0,2))`, [attribution.rs:106](../../../../src/ssr/cohort/attribution.rs#L106)); the new code lets `align_subst` decide. For composition-asymmetric equidistant alleles (an interruption sibling — this feature's data), `banded_forward` scores are asymmetric and can select the other allele, flipping `Δ = read_units − called[best].1` (e.g. `+1 → −1`), moving counts between `fit.profile.up`/`down` in `attribute_locus` ([em.rs:635](../../../../src/ssr/cohort/em.rs#L635)) and shifting `θ_locus`. It is deterministic (no byte-identity break) — a *silent statistical shift*. The doc calls the single-nearest call "ignored" and the change "effect-neutral"; both are false for slip / equidistant reads. Untested (the new tests cover only same-length and non-equidistant-length cases). No demonstrated miscall (the benchmark concordance held), so impact is bounded to a nuisance parameter on a rare read class — but it contradicts a documented claim and is unverified.

**Fix (one, covers both):** short-circuit the single-nearest case, and scope the composition tie-break to genuine same-length ties (matching `nearest_parent` for different-length equidistant), then either restrict as below or, if the composition tie-break on equidistant-different-length is *intended*, correct the doc and add a test pinning it:
```rust
let mut nearest = called.iter().enumerate()
    .filter(|&(_, &(_, u))| (i32::from(u) - read_units).abs() == min_dist);
let (first_idx, _) = nearest.next()?;
let best_idx = if nearest.next().is_none() {
    first_idx // single nearest → align_subst is irrelevant, skip the DP
} else {
    // ≥2 tie: max align_subst over the tied set, lowest index on tie
    // (optionally restrict the tied set to same-length alleles to match nearest_parent)
    /* existing reduce over the tied set */
};
```

### Minor

**Mi1: src/ssr/cohort/param_estimation.rs:110 — `interruption_count` over-scores an interruption in the first repeat unit** (reliability, Medium)
`seq[i] != seq[i % period]` compares each base to the *first unit*; if the interruption is in unit 1 (e.g. period-2 `CA` tract `TACACACACACA`) every later position mismatches the corrupted `seq[0]`, yielding ~`units−1` interruptions instead of 1. That inflated `k` feeds `PurityLevel::level_factor` as `factor^k`, collapsing the allele's modelled stutter toward 0 — a wrong result, no panic. The doc names the precondition; nothing enforces it. **Fix:** add a `debug_assert!`/documented invariant at the call sites (`cand_interruptions` in em.rs, `peak_interruptions` in prepass.rs) that tracts are phase-aligned, or score against the cohort motif rather than `seq[0..period]`; at minimum add `interruption_count_overcounts_a_first_unit_interruption` (§8).

**Mi2: src/ssr/cohort/vcf_out.rs:142 — `allele_balance` silently disables the FP defence when `allele_support` is empty, untested** (reliability + refactor_safety convergent, Medium)
`apply_fp_control` no longer sees the locus reads; `allele_balance` returns `1.0` (skip the test) when `support.len() < 2`. Safe today (`fill_allele_support` always runs before return, [em.rs:385](../../../../src/ssr/cohort/em.rs#L385)), but the coupling is implicit: any future path building a `LocusCall` with populated `allele_indices` but empty `allele_support` loses the depth-inflated-false-het defence with no error. **Fix:** `debug_assert!(!call.allele_support.is_empty())` for a non-empty het in `apply_fp_control`, plus `allele_balance_returns_one_when_support_absent` (§8).

**Mi3: src/ssr/cohort/em.rs:427 — `candidate_level` recomputed per `(obs, candidate)` though independent of `obs`** (reliability + extras convergent, efficiency, Low)
In `compute_data_ll` the level is recomputed for every read row against every candidate, but it depends only on `(level, cand_units[c], multiplier, purity, cand_interruptions[c])` — constant across reads. Hoisting to a per-candidate `Vec<f64>` once per sample removes `k × n_obs` redundant `powi`/`clamp`. (Pre-existing inline pattern; the purity factor added one more `powi` per row.) On the hot path, worth doing alongside M1.

**Mi4: src/ssr/cohort/em.rs:633 — the swapped `.expect` lacks the module's `// PANIC-FREE:` marker** (errors, High)
Every other new panic site in the module (`sample_chemistry`, `final_calls`, `format_vcf_record`, `genotype_locus`) documents its invariant; the `nearest_called_by_sequence(...).expect(...)` does not. It is provably guarded (the loop skips empty `allele_indices`; `called` is parallel), but a future edit letting `genotype_units`/`allele_indices` diverge would turn it into a live panic with no local record. **Fix:** add the marker comment.

**Mi5: src/ssr/cohort/em_init.rs:86 — `sample_eps`'s decision-E panic path is untested** (reliability, High)
`sample_eps` duplicates `sample_chemistry`'s two `.expect`-guarded lookups but, unlike `sample_chemistry` (which has `sample_chemistry_panics_on_a_sample_missing_its_group`), has no regression test. The two can drift; without a test the em_init copy could be silently weakened to a default-ε fallback — the exact silent-default that decision-E forbids. **Fix:** add `sample_eps_panics_on_a_sample_missing_its_group` (§8).

**Mi6: src/ssr/cohort/prepass.rs:123 — `MIN_PURITY_FIT_READS = 50` doc cites a figure (100/side) that neither matches the value nor the gated quantity** (defaults, Medium)
The doc says "≥100/side" but the const is 50, applied to the *cell total* `f + s` (not per-side). Halving is unexplained and "each side" is wrong. **Fix:** state the real rationale for 50 and align the wording with the gated quantity (cell total).

**Mi7: src/ssr/cohort/prepass.rs:126 — `PURITY_FACTOR_FLOOR = 0.05` documents the direction but not the source of the magnitude** (defaults, Medium)
No spec/measurement source for `0.05` vs `0.1`/`0.01`. **Fix:** record the origin ("conservative guard, no measurement yet; revisit in F2").

**Mi8: src/ssr/cohort/candidate_set.rs:78 — the §5.2 defaults are literals in `dev_default()` *and* restated in each field's prose → drift risk** (defaults, Medium)
`min_same_length_{reads,samples,fraction}` values live in two places per field, and these three are exactly the precision/recall knob slated for a P1.5/F2 sweep — the value most likely to be edited in one place only. `dev_default()`'s own doc does not enumerate them. **Fix:** hoist to named `const`s referenced by both the field doc and `dev_default()`, or drop the numeric value from the field prose. (Module-wide `*Cfg::dev_default` convention — a convention choice.)

**Mi9: src/ssr/cohort/prepass.rs:123 — the purity-fit gates are non-sweepable module consts, asymmetric with the sibling `CandidateCfg` knobs** (defaults, Medium)
The §5.2 bar is on `CandidateCfg` (tunable, `Debug`-inspectable, "swept in P1.5"); the two equally precision/recall-affecting purity gates are baked into `fit_purity_level`, needing a recompile to move and invisible at runtime. **Fix:** lift into a small `PurityFitCfg { min_fit_reads, factor_floor }` (mirroring `G0FitCfg`), or document them as deliberately frozen.

**Mi10: src/ssr/cohort/driver.rs:1138 — the ~185-line P2.0 diagnostic buried in `mod tests` re-clones `run`'s two-pass pipeline; its `#[ignore]` gives no removal condition** (module_structure + smells convergent, Low)
`AlleleSlip` / `accumulate_allele_slips` / `measure_allele_slips` / the `#[ignore]`d `p20_dump_per_allele_slip_stats` are a file-writing, env-var-driven measurement tool, and `measure_allele_slips` hand-duplicates the burn-in + per-locus loop from `run`. The repo's own pattern (`examples/`) is one directory over. **Fix:** promote to an example or a dedicated `#[cfg(test)]` module that calls the existing driver passes instead of re-deriving them; add a removal condition to the `#[ignore]` reason (the gate is settled — P2.0 is done).

**Mi11: cross-module — duplicated test fixtures across six test modules** (smells, Low)
`ca_seq`/`ca`, `pure6`/`interrupted6`/`noise6`, the `s[5]=b'T'` interruption construction, and the byte-sorting sample builder are copy-pasted across rung_ladder / candidate_set / em_init / em / vcf_out / driver test modules. **Fix:** one `#[cfg(test)]` `test_fixtures` helper module.

### Nits

- **Naming:** `purity_slip` ([prepass.rs:59](../../../../src/ssr/cohort/prepass.rs#L59)) conveys neither its `(length,purity)` keys nor `(faithful,slipped)` values → `slip_counts_by_length_and_purity`; `PurityLevel` reads as a measure but holds a *factor* → consider `PurityLevelFactor`; `RungSeq.reads`/`.samples` are `u32` counts named like collections; `called` → `called_alleles`; `interruption_count` is a noun-phrase for a fn (verb `count_interruptions` — minor).
- **Casts:** prefer `f64::from`/`u64::from` over `as` where the module already does, at `em.rs:479` (`*count as f64`), `param_estimation.rs` (`interruptions as i32` in `powi`, `.count() as u32`). (The `u64→f64` casts in `fit_purity_level` are *not* nits — no `f64: From<u64>`.)
- **Duplication:** the `align_subst`+`reduce` lowest-index tie-break is written twice (em_init.rs + attribution.rs); `purity_slip`'s `(u64,u64)` tuple-of-tuples invites a small named struct.
- **examples/ssr_slip_dump.rs:** the committed "reads-only proxy" the in-caller `measure_allele_slips` supersedes; overlaps the untracked `examples/ssr_psp_seqdump.rs`. Keep one, point to the other. Its `chroms[chrom_id]` index (line 41) and lack of `--help`/arg validation are acceptable for a dev binary.
- **idiomatic (test-only):** `accumulate_allele_slips` clones a `Box<[u8]>` as a BTreeMap key on every read even on the entry-exists path (`acc.entry(...clone())`); `get_mut`-then-insert clones only on miss.
- **defaults:** the `F2`/`P1.5` milestone codes on `dev_default()` docs are opaque at the call site with no pointer to where they're defined.

## 7. Out of scope observations

- **`src/vcf/writer.rs:384`** — `needless_lifetimes` clippy error fails the crate-wide `clippy -D warnings` CI gate (and, since the lib won't compile under `-D`, blocks clippy from fully checking in-scope targets). Pre-existing, SNP path. Follow-up: elide the lifetime in a separate housekeeping change.
- **`src/var_calling/posterior_engine.rs`** — `larger_ref_pseudocount_cannot_increase_p_alt` fails deterministically (a monotonicity property broken by the recent SFS-prior work per PROJECT_STATUS). Pre-existing, out of scope. Follow-up: investigate on the SNP side.
- **`src/ssr/cohort/read_model/mod.rs:20`** (`[ClassicStutterModel]`, a `#[cfg(test)]` item) **and `src/var_calling/diversity.rs:6-7`** — unresolved intra-doc links fail `cargo doc`. Pre-existing. Follow-up: escape or fix the links.
- **Untracked working-tree files:** `examples/ssr_psp_seqdump.rs` and `proptest-regressions/` are present but not committed — confirm intentional.

## 8. Missing tests to add now

From the `reliability` challenge-test pass (grouped by function):

- **`interruption_count`** — `interruption_count_returns_zero_for_a_tract_shorter_than_one_period` (empty / `"C"` / exactly-one-period → 0; guards the `period..len` range on a 0-unit deletion allele); `interruption_count_overcounts_a_first_unit_interruption` (pins the Mi1 phase precondition — `interruption_count(b"TACACACACACA",2) > 1`; update if a phase guard lands).
- **`nearest_called_by_sequence`** — `nearest_called_by_sequence_breaks_an_equidistant_length_tie_by_composition` (read equidistant from two *different-length* alleles; pins the M1 semantics — the current suite tests only same-length and clean-winner cases). Add both the 2-unit-slip case (`align_subst`≈0 → lowest index) and the **1-unit composition-asymmetric** case (the one M1 can flip).
- **`candidate_for_sequence`** — `candidate_for_sequence_returns_none_when_no_candidate_sits_at_the_length` (the untested `None` branch).
- **`fill_allele_support`** — `fill_allele_support_skips_a_read_matching_no_called_allele` (the `Σ Qᵣ ≤ 0` junk-read skip; asserts support sums to the depth of *matching* reads).
- **`fit_purity_level`** — `fit_purity_level_floors_an_extreme_contrast` (an extreme ratio clamps to `PURITY_FACTOR_FLOOR`); `fit_purity_level_fits_multiple_interruption_levels_through_origin` (k=1 and k=2 at one length → compounding `factor≈0.5`, pins the through-origin WLS end-to-end).
- **`allele_balance`** — `allele_balance_returns_one_when_support_absent` (Mi2 defensive branch).
- **`em_init::sample_eps`** — `sample_eps_panics_on_a_sample_missing_its_group` (Mi5 decision-E invariant, `#[should_panic]`).
- **e2e determinism** — `run_is_byte_identical_across_threads_with_a_purity_contrast`: the existing T1-vs-T4 test uses a clean fixture that leaves `purity_level` at `none()`, so the Phase-2 fit + `candidate_level` scaling are never covered by a cross-thread test.

## 9. What's good

- **`candidate_level` (em.rs)** genuinely dedupes the per-candidate stutter-level formula so the data likelihood and the allele-balance deconvolution can never disagree on the level — a single source of truth extracted at exactly the right seam.
- **`fit_purity_level` (prepass.rs)** sorts its cells into canonical `(length, k)` order *before* the float `sxy/sxx` reduce, keeping the byte-identity contract that a HashMap-iteration-order sum would have broken — the determinism instinct applied to genuinely new statistics.
- **Panic discipline:** nearly every new `.expect`/`panic!` carries a `// PANIC-FREE:` invariant or a deliberate loud-fail doc (the decision-E contract), and every `Option`/`None`/zero-evidence path degrades to a documented safe fallback rather than a silent wrong answer.
- **The Q-I1 decoupling** (compute deconvolved support in `em::fill_allele_support`, consume it in `vcf_out::allele_balance`) *paid off*: `vcf_out.rs` dropped its `nearest_parent`/`SampleEvidence` imports and now reads a field — a real boundary simplification, not a reshuffle.
- **Phase-2 collapse-to-Phase-1** is airtight and verified: `PurityLevel::none()` → factor 1, a pure allele → factor 1, and `fit_purity_level` → `none()` with no contrast, so the guarded feature is a byte-for-byte no-op until the data activates it.

## 10. Commands to re-verify

Run in the container (`./scripts/dev.sh`):
- `cargo test --lib ssr::cohort` (expect 206 passed; the branch's own suite).
- `cargo fmt --check` (in-scope clean; ignore pre-existing `src/paralog/` drift).
- After M1/Mi-fixes: `cargo clippy --lib` on the touched files, and the new tests in §8.
- Determinism spot-check: `cargo test --lib ssr::cohort::driver::tests::same_length_interruption_locus_emits_a_variable_row_identically_across_threads`.

### Author response convention
Address each finding by identifier (M1, Mi1…) with `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Resolve the §4 open questions first (they gate M1, Mi1, Mi7, Mi9).
