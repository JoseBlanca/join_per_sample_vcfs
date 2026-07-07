# Fix Application Report: rough_snp_het_caller_2026-07-07.md

**Date:** 2026-07-07
**Source review:** `doc/devel/reports/reviews/rough_snp_het_caller_2026-07-07.md`
**Source state reviewed against:** commit `87cc19b`, branch `rough-snp-heuristics`
**Execution mode:** interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 6 (M1–M6)
- Minors: 7 (Mi1–Mi7)
- Nits: 3 (grouped)

### User decisions (up front)
- **Q1 (M1/M2):** remove the `PVC_HET_STRAND_Z` env override as scaffolding. → M1 resolved by removal; M2 still applied (persist `strand_bias_z`).
- **Q2 (M5/Mi1/Mi4):** minimal + shared constructor — keep fields public, add `SiteCounts::from_record`, `debug_assert`+`saturating_sub` guard.

### Outcome totals
- Applied: 11 (M1, M3, M4, M5, M6, Mi1, Mi2, Mi3, Mi4, Mi5, Mi6)
- Applied with adaptation: 1 (M2 — required-field + schema version bump, not `serde(default)`)
- Already fixed: 0
- Deferred: 2 (Mi7; the `Option<f64>`-sentinel nit)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 1 (the fmt nit → M3)
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → **in-scope clean.** My files are formatted; the only remaining diff is pre-existing drift in the untouched `src/var_calling/vcf_writer.rs` (restored deliberately, out of scope).
- `cargo clippy --lib --all-features -- -D warnings` → **exit 0, clean** (the changed code is fully linted after M4). `--all-targets` still aborts on a pre-existing `sort_by_key` lint in the untouched `examples/ssr_psp_seqdump.rs` (out of scope; previously masked by the M4 error).
- `cargo test --all-targets --all-features` → **1633 passed, 1 failed.** +12 new tests, all passing. The 1 failure is the pre-existing, out-of-scope `var_calling::posterior_engine::tests::larger_ref_pseudocount_cannot_increase_p_alt`.
- `cargo doc --no-deps` → **no new errors.** The one link I introduced (`het::is_strand_biased`, a private fn) was fixed; 4 pre-existing unresolved links remain in untouched ssr/var_calling modules.
- `cargo audit` → not run (`cargo-audit` not installed; no dependency changes).
- Performance check → **skipped.** See §9.

### Unresolved high-priority findings
- None. All Majors Applied and validated.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | `PVC_HET_STRAND_Z` env override unvalidated | Apply (remove) | Applied | Q1 | cli.rs | Pass | No |
| M2 | Major | `strand_bias_z` not persisted in `.psp` | Apply | Applied with adaptation | Q1 | mod.rs, het.rs, diversity.rs, coverage.rs, inbreeding.rs, prepass.rs, test_support.rs | Pass | No |
| M3 | Major | Diff not `cargo fmt`-clean | Apply | Applied | No | het_baseline.rs, pileup_to_psp.rs, (+all touched) | Pass | No |
| M4 | Major | New code never clippy-linted | Apply | Applied | No | vcf/writer.rs | Pass | Pre-existing example lint remains (OOS) |
| M5 | Major | `strand_biased` trusts unenforced `fwd ≤ obs` | Apply | Applied | Q2 | het.rs | Pass | No |
| M6 | Major | Constructor guards untested | Apply | Applied | No | het.rs (tests) | Pass | No |
| Mi1 | Minor | Duplicated feed + non-exhaustive ALT build | Apply | Applied | Q2 | het.rs, pileup_to_psp.rs, paralog_score_parity.rs | Pass | No |
| Mi2 | Minor | Stale `HetAccumulator::new` panic doc | Apply | Applied | No | het.rs | Pass | No |
| Mi3 | Minor | `strand_biased` → `is_strand_biased` | Apply | Applied | No | het.rs | Pass | No |
| Mi4 | Minor | `AlleleGroupStats.obs` → `num_obs` | Apply | Applied | Q2 | het.rs, pileup_to_psp.rs | Pass | No |
| Mi5 | Minor | `DEFAULT_HET_*` family split | Apply | Applied | No | mod.rs, het.rs | Pass | No |
| Mi6 | Minor | `min_depth ≥ 1` invariant unenforced | Apply | Applied | No | het.rs | Pass | No |
| Mi7 | Minor | `SiteCounts` `Copy` at 48 bytes | Defer | Deferred | No | None | N/A | Low-priority |
| Nit-var | Nit | Unreachable `variance<=0` comment | Apply | Applied | No | het.rs | Pass | No |
| Nit-opt | Nit | `Option<f64>` disable-sentinel | Defer | Deferred | No | None | N/A | Ripple cost |
| Nit-fmt | Nit | fmt (subsumed by M3) | Superseded | Superseded | No | — | — | By M3 |

## 3. Questions asked and answers

1. **M1/M2 (Q1)** — Keep `PVC_HET_STRAND_Z` as a shipped knob, or remove it as scaffolding?
   - **Answer:** Remove it. (M1 resolved by removal; M2 still applied.)
2. **M5/Mi1/Mi4 (Q2)** — Keep fields public with a shared constructor + defensive guard, or make fields private behind a validating constructor?
   - **Answer:** Minimal + shared constructor (keep public).

## 4. Per-finding log

### M1 — `PVC_HET_STRAND_Z` env override unvalidated
- **Severity:** Major · **Initial:** Apply (remove) · **Final:** Applied
- **Reasoning:** User chose to remove the override (session scaffolding). Removal eliminates both the silent-fallback and the panic-on-bad-value paths — the finding dissolves.
- **Implementation:** `cli.rs` now sets `strand_bias_z: DEFAULT_HET_STRAND_BIAS_Z` directly (the `std::env::var(...)` chain deleted).
- **Files changed:** `src/pop_var_caller/cli.rs`
- **Validation:** builds + clippy `--lib` clean.
- **Residual risk:** None.

### M2 — `strand_bias_z` not persisted in the `.psp`
- **Severity:** Major · **Initial:** Apply · **Final:** Applied with adaptation
- **Reasoning:** The review suggested adding the field "beside the three existing knobs." The three existing knobs are *required* v-fields; the crate's documented policy (`SAMPLE_SUMMARY_VERSION` doc) is that summary fields are required and old documents are regenerated, not `serde(default)`-backfilled. I adapted accordingly.
- **Implementation:** added required `strand_bias_z` to `HetCounts`; bumped `SAMPLE_SUMMARY_VERSION` 3 → 4 with a v4 doc line; `finish()` populates it; `validate()` rejects `NaN`/`≤ 0`; the 2 exhaustive `HetCounts` destructures in `diversity.rs` and all literal construction sites updated (compiler-guided). Tests: `rejects_strand_bias_z_out_of_range`, `document_missing_strand_bias_z_is_rejected`, and `strand-bias-z` added to `wire_keys_are_kebab_case`.
- **Review suggestion used verbatim?:** No.
- **Adaptation:** required-field + version bump instead of `serde(default)`, to match the crate's on-disk-schema convention. (Het counts are unchanged; only a recorded provenance field is added.)
- **Files changed:** `src/sample_summary/mod.rs`, `src/sample_summary/het.rs`, `src/sample_summary/coverage.rs`, `src/var_calling/diversity.rs`, `src/paralog/inbreeding.rs`, `src/var_calling/paralog_filter/{prepass,test_support}.rs`
- **Tests added:** `rejects_strand_bias_z_out_of_range`, `document_missing_strand_bias_z_is_rejected`
- **Validation:** lib tests pass (round-trip, wire-keys, rejection).
- **Residual risk:** existing v3 `.psp` are rejected and must be regenerated by re-running `pileup` — the crate's stated pre-alpha policy.

### M3 — diff not `cargo fmt`-clean
- **Severity:** Major · **Final:** Applied
- **Implementation:** `cargo fmt` over the touched files; the 4 pre-existing-drift files (`ssr_psp_seqdump.rs`, `paralog/locus_score.rs`, `paralog_filter/write_pass.rs`, `var_calling/vcf_writer.rs`) `git checkout`-restored so they stay out of this diff.
- **Validation:** `cargo fmt --check` shows no diff in any in-scope file.
- **Residual risk:** the crate-wide fmt gate remains red on the untouched `vcf_writer.rs` (pre-existing, §7).

### M4 — new code never clippy-linted
- **Severity:** Major · **Final:** Applied
- **Reasoning:** the clippy gate aborted on a pre-existing `needless_lifetimes` in `src/vcf/writer.rs:384` before reaching the changed lib. Fixed that one line (compiler-suggested elision) so clippy compiles through the lib; `cargo clippy --lib` is now clean, confirming the new ε̂/strand code lints cleanly.
- **Implementation:** `fn table<'a>(cache: &'a mut …) -> &'a …` → elided lifetimes.
- **Files changed:** `src/vcf/writer.rs` (a single out-of-scope-but-necessary line, to unblock validation of the in-scope code).
- **Validation:** `cargo clippy --lib --all-features -- -D warnings` → exit 0.
- **Residual risk:** `--all-targets` clippy now surfaces a *different* pre-existing lint in the untouched `examples/ssr_psp_seqdump.rs` (§7). Not introduced here.

### M5 — `strand_biased` trusts unenforced `fwd ≤ obs`
- **Severity:** Major · **Final:** Applied
- **Implementation:** `is_strand_biased` now `debug_assert!(rf <= nr && af <= na, …)` and uses `saturating_sub` for the reverse cells, so a violating direct construction panics in debug and degrades to "no veto" (reverse cell 0 → margin guard rejects) in release rather than wrapping. Test `strand_veto_rejects_fwd_exceeding_num_obs`.
- **Files changed:** `src/sample_summary/het.rs`
- **Validation:** the `#[should_panic]` test passes; full het suite green.

### M6 — constructor guards untested
- **Severity:** Major · **Final:** Applied
- **Implementation:** added `#[should_panic]` tests `error_rate_at_ceiling_panics`, `nan_strand_bias_z_panics`, `zero_strand_bias_z_panics` (and `zero_min_depth_panics` for Mi6).
- **Files changed:** `src/sample_summary/het.rs`
- **Validation:** all pass.

### Mi1 — duplicated feed + non-exhaustive ALT build
- **Severity:** Minor · **Final:** Applied
- **Implementation:** added `SiteCounts::from_record(&PileupRecord) -> Option<SiteCounts>` (REF via `AlleleGroupStats::from_support`; ALT pooled as an **exhaustive** struct literal). Both feed sites (`pileup_to_psp.rs`, `paralog_score_parity.rs`) call it, removing the duplicate loop; a new `AlleleGroupStats` field is now a compile error in the ALT build. Test `from_record_pools_alts_and_skips_pure_ref`.
- **Files changed:** `src/sample_summary/het.rs`, `src/pileup/per_sample/pileup_to_psp.rs`, `examples/paralog_score_parity.rs`
- **Adaptation:** `from_record` lives on `SiteCounts` (het.rs now imports `pileup_record`, the universal data model — an acceptable edge).
- **Validation:** feed byte-behavior unchanged (het tests + from_record test pass).

### Mi2 — stale `HetAccumulator::new` panic doc
- **Final:** Applied. Rewrote the `# Panics` section to the current `(0, MAX_SITE_ERROR_RATE)` bound + `strand_bias_z` + `min_depth`.

### Mi3 — `strand_biased` → `is_strand_biased`
- **Final:** Applied (predicate naming; call site updated).

### Mi4 — `AlleleGroupStats.obs` → `num_obs`
- **Final:** Applied (matches the record's canonical `num_obs`; all uses updated).

### Mi5 — `DEFAULT_HET_*` family split
- **Final:** Applied. Moved the `3.0` literal (with its doc) to `mod.rs` as `DEFAULT_HET_STRAND_BIAS_Z`, matching the three sibling literals; deleted `het::DEFAULT_STRAND_BIAS_Z`; the test helper uses the `3.0` literal like its neighbours.

### Mi6 — `min_depth ≥ 1` invariant unenforced
- **Final:** Applied. `assert!(params.min_depth >= 1, …)` in the constructor + `zero_min_depth_panics` test.

### Mi7 — `SiteCounts` `Copy` at 48 bytes
- **Final:** Deferred. No correctness impact; `Copy` is not load-bearing but dropping it is a pure-style change with no forcing reason. Carried forward.

### Nit — unreachable `variance<=0` comment
- **Final:** Applied. Comment now states the branch is unreachable after the margin guard (defense-in-depth against a future-refactor NaN).

### Nit — `Option<f64>` disable-sentinel
- **Final:** Deferred (ripple across params + all construction sites; the `INFINITY` sentinel is documented and tested).

## 5. Deferred findings to carry forward
- **Mi7** — `SiteCounts` `Copy` at 48 bytes (style; drop `Copy` when convenient).
- **Nit** — `Option<f64>` disable-sentinel for `strand_bias_z` (deferred: ripple cost).

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** No.
- **Baseline saved:** No — edits began before a baseline could be captured (the run started directly from the review).
- **Outcome:** Skipped. The `observe_site` hot-path changes are perf-neutral by inspection: the `debug_assert` compiles out of release; `saturating_sub` is the same cost as `-`; `SiteCounts::from_record` relocates the identical per-record summation (same arithmetic, same allele-index order); the veto rename is a no-op. No arithmetic on the scoring path changed.
- **Notes:** per the skill, a baseline is not back-filled after editing; recorded as skipped rather than fabricated.

## 10. Commands run
- `cargo fmt`
- `cargo fmt --check`
- `cargo clippy --lib --all-features -- -D warnings`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo doc --no-deps`

## 11. Command results
- `cargo fmt --check` → in-scope clean (pre-existing `vcf_writer.rs` drift only).
- `cargo clippy --lib …` → exit 0, clean.
- `cargo clippy --all-targets …` → aborts on pre-existing `examples/ssr_psp_seqdump.rs` `sort_by_key` (OOS).
- `cargo test --all-targets …` → 1633 passed, 1 pre-existing failure.
- `cargo doc --no-deps` → no new errors (4 pre-existing unresolved links remain).

## 12. Notes
- **M2 schema bump:** `SAMPLE_SUMMARY_VERSION` is now 4; existing `.psp` must be regenerated (`pileup`). The tomato measurement `.psp` under `benchmarks/tomato1/results/` are session scratch and can be regenerated if needed; the het *counts* are unaffected by M2 (only a recorded field is added).
- **Out-of-scope touch:** M4 edited one line of `src/vcf/writer.rs` (a pre-existing clippy error) solely to unblock linting of the in-scope code. The crate-wide clippy/fmt/doc/test gates retain other pre-existing reds (§7 of the review), none introduced here.
