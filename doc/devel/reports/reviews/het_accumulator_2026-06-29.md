# Code Review: het accumulator + HetCounts amendment (B3)

**Date:** 2026-06-29
**Reviewer:** rust-code-review skill (orchestrator + 1 consolidated stats-focused sub-agent)
**Scope:** feature B3 — `src/sample_summary/het.rs` (new) + `HetCounts` two→four-count amendment in `mod.rs` + spec/architecture amendments
**Status:** Approve-with-changes

---

## 1. Scope
- Reviewed: B3 working diff on branch `tomato2-paralog-filter`.
- In-scope: `src/sample_summary/het.rs` (new), `mod.rs` (`HetCounts`),
  `coverage.rs` (one test literal). Companion doc amendments:
  spec `hidden_paralog_filter.md` §3, architecture Premise 1b/2/4.
- Categories (consolidated, stats-focused): reliability, errors,
  idiomatic, naming, smells, refactor_safety.

## 2. Verdict
**Approve-with-changes.** The statistics are correct — the reviewer
numerically verified the binomial-coefficient cancellation is valid for
both comparisons (variant gate and het/hom split are differences, so
`ln C(n,k)` cancels), the `>`/`<`/`≤` margins are internally consistent,
and the documented depth-aware cases reproduce exactly. Findings about
the input-precondition surface and test coverage.

## 3. Execution status
- `cargo test --lib sample_summary` → 24 passed (pre-fix incl. 8 het).
- `cargo fmt --check` → exit 0. `cargo clippy --lib` → clean on scope.

## 4. Findings

### Major
- **M1 (errors / reliability / idiomatic — one defect, three angles):
  `observe_site(k_alt, n_total)` only debug-asserted `k ≤ n`.** In release
  a `k > n` argument makes `ref_obs` negative and silently miscounts the
  site as confident hom-alt (verified numerically), polluting `Hobs` with
  no signal; the two same-typed args are also transposable. Fix: take a
  **`SiteCounts { ref_obs, alt_obs }`** input and derive `n = ref + alt`,
  so `alt ≤ total` holds **by construction** — the bad state is
  unrepresentable and the transposition risk is gone. *(Applied.)*
- **M2 (reliability): no exact-margin boundary test** (`logLR == +M`,
  gate `== M`). The strict `>`/`<` vs `≤` band is the classifier's whole
  contract; a comparator flip would pass every existing test. *(Applied —
  `site_at_exact_gate_margin_is_rejected`,
  `site_at_exact_split_margin_is_ambiguous`.)*
- **M3 (reliability): no `k==n`/`k==0` test at small `n` near the gate**
  (the `ref_obs == 0` arithmetic edge, where collapsed-paralog artefacts
  live). *(Applied — `pure_hom_alt_at_min_depth_is_hom_alt`.)*

### Minor
- **Mi1 (reliability): `min_depth` inclusion boundary untested.**
  *(Applied — `site_at_exactly_min_depth_is_scored`.)*
- **Mi2 (reliability): no property test for the accumulator's four-count
  identity.** *(Applied — `accumulator_maintains_count_identity` proptest
  over arbitrary `(ref, alt)` sequences.)*
- **Mi3 (errors): `+= 1` / `finish()` sum unchecked while the consumer
  `validate()` uses `checked_add`.** *(Applied — comment naming the
  site-count bound `<< u64::MAX`, making the discipline explicit.)*
- **Mi4 (idiomatic): `max`/comparison NaN-freedom depended on an
  undocumented invariant.** *(Applied — comment ties finiteness to the
  constructor asserts + the by-construction `alt ≤ total`.)*

### Nits
- ln(10) margin literal repeated across tests → use `std::f64::consts::LN_10`.
  *(Applied.)*

### Clean categories
naming, refactor_safety — no findings (every `HetCounts` literal spells
all fields, which is why the new fields correctly forced the `coverage.rs`
test update).

## 9. What's good
- Binomial-coefficient cancellation makes the classifier a few multiplies
  with no log-gamma — "rough" in cost, principled in shape.
- Depth-awareness is a tested property, not a claim (5/6@6× ambiguous,
  59/60@60× hom-alt).
- The `SiteCounts` newtype makes `alt ≤ total` a type invariant.

## 10. Commands to re-verify
- `./scripts/dev.sh cargo test --lib sample_summary::het`
- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings -A clippy::doc_lazy_continuation`
