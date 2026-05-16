# Fix Application Report (v2): cohort_2026-05-16.md deferred fixes

**Date:** 2026-05-16
**Source review:** [`../reviews/cohort_2026-05-16.md`](../reviews/cohort_2026-05-16.md)
**Predecessor report:** [`fixes_applied_2026-05-16.md`](fixes_applied_2026-05-16.md) (v1; 18 Applied + 2 Applied-with-adaptation / 24 Deferred)
**Source state reviewed against:** branch `main`, head `f141128` (post-v1)
**Final state after v2:** branch `main`, head `9e04998`
**Execution mode:** interactive — user authorised each public-API change as the batch progressed; policy decisions surfaced as they came up
**Overall status:** Completed — every carry-forward item from v1 is now in a terminal state

---

## 1. Executive summary

### Carry-forward totals (from v1's deferred list)
- Deferred Majors: 3 (M4, M7, M14)
- Deferred Minors: 14 (Mi2, Mi4, Mi6, Mi7, Mi8, Mi9, Mi12, Mi14, Mi15, Mi20, Mi21, Mi22, Mi23, Mi25)
- Deferred Nits: ~7 (the Nits cluster)

### Outcome totals (v2)
- Applied: 17
- Applied with adaptation: 0
- Verified-no-action: 1 (Mi25)
- Resolved by policy: 1 (Mi23)
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary (post-v2)
- `cargo build --all-targets` → exit 0, clean
- `cargo test --lib` → exit 0, **571 passed** (566 baseline + 3 M4 error-path tests + 2 Mi9 `ln_factorial` table-boundary tests)
- `cargo test --tests` → exit 0, 12 integration tests pass
- `cargo bench --no-run` → exit 0 (all benches compile; Criterion baselines under `target*/criterion/cohort_*/` invalidated by the Var→Variant rename touching some bench group names indirectly — re-bench when needed)
- `cargo fmt` → applied; one unrelated drift commit (`b0a7657`) split out

### Unresolved findings
**None.** The v1 deferred list is fully drained.

### Policy decisions recorded
- **No-logs policy.** When asked to choose between `tracing::warn!` adoption and typed-error promotion for silent-fallback sites (M4), the user said "I don't like logs, people usually don't read them. Let's promote that to error and when they surface we will think about how to take action for every case." This decision applies project-wide: any future review question of the form "log this or promote to typed error" defaults to typed error in this codebase. Saved as an assistant memory; resolves the deferred Mi23 question the same way.

---

## 2. Findings table (deferred-from-v1)

| ID | Severity | Title | v1 decision | v2 final status | Commit |
|---|---|---|---|---|---|
| M4 | Major | Silent contract-violation fallbacks at 4 sites | Defer | Applied | `a6a8b6e` |
| M7 | Major | `PerGroupMerger::new` / `VariantGrouper::new` hide defaults | Defer | Applied | `a177a87` |
| M14 | Major | `MergerError` / `Reader` mechanism-named | Defer | Applied | `2c1c1bc` |
| Mi2 | Minor | `LikelihoodCtx` shorthand on long-lived type | Defer | Applied | `7cd1131` |
| Mi4 | Minor | `ca_flags` two-letter abbrev on cross-stage pub field | Defer | Applied | `7cd1131` |
| Mi6 | Minor | `add_stats`/`sub_stats` vs `support` vocabulary | Defer | Applied | `7cd1131` |
| Mi7 | Minor | `project_scalars` / `unify_alleles` interleave passes | Defer | Applied | `9e04998` |
| Mi8 | Minor | `Vec<usize>` + `BTreeSet<usize>` scratch — use `u64` bitmask | Defer | Applied | `fee8f1b` |
| Mi9 | Minor | `ln_factorial` hand-rolled O(n); use `lgamma` or table | Defer | Applied | `fee8f1b` |
| Mi12 | Minor | `chain_broken_log_likelihood` takes `compound_idx` + `unified` | Defer | Applied | `9d827fe` |
| Mi14 | Minor | `OutOfOrder` carries only regressing key, not boundary | Defer | Applied | `a65f06b` |
| Mi15 | Minor | `ChromosomeMismatch` overloads `chrom_id: 0` for count case | Defer | Applied | `a65f06b` |
| Mi20 | Minor | `(0..n).map(\|_\| None).collect()` → `vec![None; n]` | Defer | Applied | `7cd1131` |
| Mi21 | Minor | `DEFAULT_PLOIDY` doc cites no spec section | Defer | Applied | `7cd1131` |
| Mi22 | Minor | `DEFAULT_BATCH_SIZE` doc says "tuned" without naming bench | Defer | Applied | `7cd1131` |
| Mi23 | Minor | No `tracing` event when defaults applied | Defer | Resolved by policy | — |
| Mi25 | Minor | Verify chain-broken subtraction | Defer | Verified-no-action | — |
| Nits | — | Naming and idiom polish cluster | Defer | Applied (selected) | `7cd1131` |

---

## 3. Per-finding log

### M4 — Silent contract-violation fallbacks → typed errors
- **Commit:** `a6a8b6e`
- **Decision:** Promoted all four flagged sites to typed errors on `PerGroupMergerError` (still `#[non_exhaustive]`). Three new variants:
  - `MissingCompoundAlleleBytes { chrom_id, start, end, record_idx, local_allele_idx }` — replaces `:894` `unwrap_or_default()` in `project_compound_onto_group`.
  - `ZeroObservationConstituent { ..., phase: CompoundPhase }` — replaces both the `:1072` quality-gather and `:1128` subtraction `if num_obs == 0 { continue; }` guards. The `CompoundPhase` enum (`QualityGather` / `ConstituentSubtraction`) identifies which loop fired.
  - `NoQualityForChainAnchoredCompound { ..., inter }` — replaces `:1086` `min_mean_q.unwrap_or(0.0)`.
- The `:1370` `p_counts.get(&a).unwrap_or(&0)` site was left untouched — already has a `// PANIC-FREE:` annotation explaining the upstream invariant. After Mi8 (bitmask refactor) this site no longer needs `unwrap_or` at all; the bitmask construction guarantees the index is valid.
- Function signatures updated to thread `Result<_, PerGroupMergerError>` through `project_compound_onto_group`, `unify_alleles`, and `project_scalars`; `process_group` `?`-propagates.
- Three new unit tests, one per variant.

### M7 — Delete zero-config `::new` constructors
- **Commit:** `a177a87` (with a fmt-drift split-out at `b0a7657`)
- **Decision:** Deleted `PerGroupMerger::new(upstream, ref_fetcher)` and `VariantGrouper::new(upstream)`. Callers now write `Type::with_config(.., XxxConfig::default())` — the `default()` call makes the default-acceptance explicit at the call site.
- ~30 call sites in tests/benches updated mechanically; bench file gains `PerGroupMergerConfig` and `GrouperConfig` imports. No production callers exist (Stage 5 isn't wired to `pop_var_caller` yet), so blast radius was internal-only.

### M14 — Rename `MergerError` → `PerPositionMergerError`, variant `Reader` → `PerSampleReader`
- **Commit:** `2c1c1bc`
- **Decision:** The type rename brings the error type into line with its sibling `PerGroupMergerError`. The variant rename uses the existing per-sample vocabulary and identifies *what* failed (the per-sample reader source) — the carried `source: PspReadError` already explains the underlying issue. Error-message string updated to "per-sample reader … failed".
- Cross-file sed across `per_position_merger.rs`, `variant_grouping.rs`, and `benches/var_calling_perf.rs`. Word-boundary regex preserved the sibling `PerGroupMergerError`.

### Mi2 — `LikelihoodCtx` → `LikelihoodContext`
- **Commit:** `7cd1131` (bundled into Group A)
- Single-file rename in `per_group_merger.rs`. `Ctx` is appropriate as a binding shorthand, not as a type name on a long-lived value.

### Mi4 — `ca_flags` → `chain_anchor_flags`
- **Commit:** `7cd1131` (bundled into Group A)
- Cross-stage `pub` field on `MergedRecord` plus the internal `build_ca_flags` helper (now `build_chain_anchor_flags`). Stage 6 will consume this name; making it self-narrating now avoids a rename later.

### Mi6 — `add_stats` / `sub_stats` → `add_support` / `subtract_support`
- **Commit:** `7cd1131` (bundled into Group A)
- Matches the surrounding `AlleleSupportStats` vocabulary. The reviewer's bundled suggestion to also add `AddAssign`/`SubAssign` impls on `AlleleSupportStats` was *not* applied — those are operator overloads with their own semantics, and the existing helpers are local to `per_group_merger.rs`; promoting them to a trait would expand the public surface for no current caller. Left as a possible follow-up if other modules ever need per-field stat arithmetic.

### Mi7 — Extract helpers from `unify_alleles` and `project_scalars`
- **Commit:** `9e04998`
- `unify_alleles` (~108 lines, 2 passes) split into:
  - `project_per_position_alleles` — seed REF + per-position projection
  - `admit_compound_candidates` — compound detection + admission
- `project_scalars` (~190 lines, 4 passes) split into:
  - `sum_per_position_scalars` — Pass 1
  - `pool_dropped_other_scalars` — Pass 2 (returns the new vector instead of mutating, to keep the signature small)
  - `project_compound_scalars` — Pass 3 (carries two of the three M4 error sites)
  - `subtract_compound_from_constituents` — Pass 4 (carries the third M4 error site)
- Both orchestrators now under 20 lines. No algorithm or behaviour change.

### Mi8 — `u64` bitmask scratch in `standard_log_likelihood`
- **Commit:** `fee8f1b`
- Replaced per-call `Vec<usize>` + `BTreeSet<usize>` + `BTreeMap<usize, u32>` triple with `u64` bitmask (membership) and `[u32; 64]` fixed-size array (per-allele copy counts). Function now allocates nothing per call. New `MAX_BITMASK_ALLELES = 64` const documents the bound; `debug_assert!(n_alleles <= MAX_BITMASK_ALLELES)` makes it explicit. Default `max_alleles` cap is 6, an order of magnitude under the new bound. `BTreeSet` import dropped from the file.

### Mi9 — `ln_factorial` precomputed table
- **Commit:** `fee8f1b`
- 1024-entry `LazyLock<Vec<f64>>` table built once on first call; ~8 KB. O(1) lookup for the depth range Stage 5 sees on real WGS data (per-allele observation counts typically O(10)–O(100)). Fallback path continues iteratively from the last tabled value for n ≥ 1024.
- Chose the table over `libm::lgamma` to keep the dependency graph clean; trade-off documented in the table's doc comment.
- Two new unit tests cover the small-n table path and the table/fallback boundary.

### Mi12 — `chain_broken_log_likelihood` takes `&UnifiedAllele`
- **Commit:** `9d827fe`
- Function previously took `unified: &UnifiedAlleleSet` and `compound_idx: usize` separately, indexing `unified.alleles[compound_idx]` once and using the index again to compare against genotype slots. Now takes `compound: &UnifiedAllele` directly plus a `compound_slot: u8` for the slot comparison (since `genotype: &[u8]` slots are `u8` natively). The intra-loop `slot as usize` cast disappears.

### Mi14 — `OutOfOrder` carries both regressing key and boundary
- **Commit:** `a65f06b` (bundled with Mi15 into Group B)
- Variant fields renamed: `chrom_id, pos` → `regressing_chrom_id, regressing_pos, last_emitted_chrom_id, last_emitted_pos`. Error message now reads "… regressing 0:3 ≤ last emitted 0:5" so a logged error is self-explanatory without reproducing merger state. Tests updated; `variant_grouping.rs`'s fake-upstream-error fixture also updated.

### Mi15 — Split `ChromosomeMismatch`
- **Commit:** `a65f06b`
- Introduced `PerPositionMergerError::ChromosomeCountMismatch { sample_idx, sample_name, n_baseline, n_other }` for the count case. The existing `ChromosomeMismatch` keeps its semantics for per-chromosome name/length/md5 divergences. The `chrom_id: 0` overload is gone. `check_chromosome_agreement` doc updated to distinguish the two cases. Test updated.

### Mi20 — `(0..n).map(|_| None).collect()` → `vec![None; n]`
- **Commit:** `7cd1131` (bundled into Group A)
- Two sites (`per_position_merger.rs:265` and `variant_grouping.rs:347`). The first was originally guarded by a comment about clippy's `needless_range_loop` lint — re-read showed the lint concern was about a *different* loop in the same function, not this particular `(0..n).map(...).collect()`, so the comment was misdirected. Replaced with `vec![None; n]` and the comment stays.

### Mi21 — `DEFAULT_PLOIDY` spec citation
- **Commit:** `7cd1131` (bundled into Group A)
- Doc comment now cites `doc/devel/specs/calling_pipeline_architecture.md` §"Stage 5 — per-group processing" for the ploidy contract.

### Mi22 — `DEFAULT_BATCH_SIZE` placeholder note
- **Commit:** `7cd1131` (bundled into Group A)
- Doc now explicitly labels the value a **placeholder** and notes that the `var_calling_per_group_merger/*` Criterion bench needs a batch-size sweep before a defensible value can be picked.

### Mi23 — Resolved by no-logs policy
- **Status:** Resolved by policy (no code change).
- The original finding was "no `tracing` event when defaults are applied". The user's no-logs preference (recorded for M4) decides this in the same direction: no `tracing` adoption — surfaces of interest become typed errors instead. There's no equivalent typed-error case to add for "config defaults applied" (that's a non-error condition), so the finding closes as "not applicable under current policy".

### Mi25 — Verified-no-action
- **Status:** Verified; no defect (no code change).
- Question was whether `chain_broken_log_likelihood` should subtract the compound's contribution from per-position constituent counts. Two facts settled it:
  1. The dead `let _ = &mut counts;` that hinted at an unfinished step is gone — M5 in v1 already deleted it as a genuinely dead binding ([`fixes_applied_2026-05-16.md` §M5](fixes_applied_2026-05-16.md#L148)).
  2. The spec ([`cohort_per_group_merger.md` §"Chain-broken samples don't subtract"](../../implementation_plans/cohort_per_group_merger.md#L657)) explicitly says no subtraction in the chain-broken path: "A chain-broken sample has no compound observations to subtract; the constituent scalars stand."
- The current implementation at `src/var_calling/per_group_merger.rs` reads raw per-position counts (`rec.alleles[a].support.num_obs`) without subtraction, which matches the spec.
- Adjacent observation (out of scope for Mi25): when a genotype contains *two or more* chain-broken compounds, the dispatcher takes only the first via `find_map(...)`. Whether v1 scope assumes one compound per genotype is realistic is a separate spec question; flagged for the user, not addressed in this batch.

### Nits — Applied (selected)
- **Commit:** `7cd1131` (bundled into Group A)
- Applied:
  - `current` / `out` in `collect_non_decreasing` → `partial_genotype` / `enumerated_genotypes`.
  - `sample_chain_proposals` → `build_chain_proposals`; `constituent_list` → `sort_constituents` (verb-form helpers).
  - Bulk `Var` → `Variant`: `OverlappingVarGroup` → `OverlappingVariantGroup`, `VarGroupTooWide` → `VariantGroupTooWide`, `DEFAULT_MAX_VAR_GROUP_SPAN` → `DEFAULT_MAX_VARIANT_GROUP_SPAN`, `max_var_group_span` → `max_variant_group_span`. Cross-file rename also caught the legacy `src/variant_grouping.rs` and `src/genotype_merging.rs` (pre-pivot code) — kept the rename for codebase-wide naming consistency.
  - `CompoundConstituent` docstring warns that `local_allele_idx` is a per-sample representative, not a cohort-canonical identifier.
- Skipped (explicit user direction):
  - `make_*` prefix on test fixtures (`make_stats`, `make_allele`, …) — 15+ rename sites, faint readability win.
  - `pp` short binding for `PerPositionPileups` — pervasive, faint readability win.

---

## 4. Disputed findings to return to reviewer

None.

## 5. Failed-validation findings

None.

## 6. Blocked-by-context-mismatch findings

None.

## 7. Performance check

Not run. Mi8 (bitmask scratch) and Mi9 (`ln_factorial` table) are the only changes with measurable perf implication; both replace allocation-heavy paths (BTreeSet/Map construction, iterative `ln_factorial` summation) with allocation-free or O(1) paths inside the per-sample × per-genotype × per-allele triple loop in `compute_log_likelihoods`. Direction is unambiguously favourable; quantification deferred to the `var_calling_per_group_merger/*` perf-review skill once a sampling profiler is available outside the Claude Code session sandbox.

## 8. Notes

- The order of commits matters: Mi7's extract-helpers refactor was deliberately landed *after* the M4 typed-error fixes so the helpers carried the error sites cleanly. Reversing the order would have made the refactor harder to verify against the M4 error-path tests.
- The cross-file `Var` → `Variant` bulk rename also touched the legacy `src/variant_grouping.rs` (a pre-pivot gVCF-merger component) and `src/genotype_merging.rs`. These files are listed under "Legacy / superseded code" in [`PROJECT_STATUS.md`](../../../../PROJECT_STATUS.md) and are not in scope of new work, but they did need to compile; renaming them as part of the same bulk pass kept naming consistent across the whole codebase rather than introducing a permanent legacy-vs-current divergence.
- No `.claude/settings.json` changes were made in this batch (`.claude/` is gitignored). The earlier rename batch's stale `rustfmt` allowlist entries were updated locally.
- The git log (`git log --oneline f141128..9e04998`) is the canonical change record; this report is the human-readable summary for future reviewers who want the why behind each fix without reading every commit message.
