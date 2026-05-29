# Fix Application Report: cohort_block_2026-05-29.md

**Date:** 2026-05-29
**Source review:** [cohort_block_2026-05-29.md](cohort_block_2026-05-29.md)
**Source state reviewed against:** branch `cohort-within-chromosome-parallel`, commit `36989d6f53460042c5219d0ff5fa6d67a7b1b129` (review committed at `536e8db`).
**Execution mode:** interactive
**Overall status:** In progress — Wave 1 (Blockers + M5/M11 unblockers + quick-win Majors M14/M17/M18/M19) and Wave 2 (mechanical Majors M10/M27/M28 + small Minors Mi14/Mi25/Mi26) complete. Remaining 24 Majors + 23 Minors + Nits + B5 staged for follow-up waves per the user's "apply all structural refactors" decision on Q4.

---

## 1. Executive summary

### Review totals
- Blockers: 5
- Majors: 32
- Minors: 26
- Nits: grouped (16 clippy errors + 14 `#[allow]` annotations needing justification + several small style items)

### Outcome totals (Waves 1–2)
- Applied: 17 (Wave 1: B1, B2, B3, B4, M5, M11, M14, M17, M18, M19 + B4-test + 2 B1-tests; Wave 2: M10, M27, M28, Mi14, Mi25, Mi26)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 47 (remaining — staged for Wave 2 onward per Q4 "apply all structural refactors")
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0 (all 4 answers received before Wave 1 closed)

### Validation summary (after Wave 1)
- `cargo fmt --check` → exit 0 (after `cargo fmt` pass).
- `cargo clippy --all-targets --all-features -- -D warnings` → exit 101; 3 E0063 errors cleared by M11; **16 in-scope clippy errors remain** (all in tests / docs / test helpers; tracked under Nits for the mechanical cleanup pass).
- `cargo test --lib` → exit 0; **1026 passed** (was 1023; +3 from the new tests for B1×2 and B4).
- `cargo test --tests` → exit 0; **21 cohort_cli_integration tests pass** including both byte-identity gates (`--worker-windows-per-chunk` and `--target-variants-per-chunk`).
- `cargo doc --no-deps` → not re-run yet (M10 doc-link fixes are part of Wave 2).
- `cargo audit` → not run (cargo-audit not installed in container; pre-existing).
- Performance check (`cargo bench -- --baseline pre-fixes`) → baseline saved; final comparison deferred to end-of-run.

### Unresolved high-priority findings (after Wave 1)
- **B5** — byte-identity oracle integration test (deferred to Wave 2; shape decided per Q2: oracle vs streaming `drive_cohort_pipeline`).
- 28 Majors and 26 Minors carried forward into Waves 2–N; structural refactors (M12–M16, Mi15, Mi17–Mi19, Mi7) are the largest items, each a focused multi-file change.

## 2. Findings table

Initial decisions captured here before any editing. Final status filled in
incrementally as each finding is processed.

| ID  | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1  | Blocker  | Writer un-finished + tmp leak on driver-error path | Apply | Applied | No | `src/vcf/writer.rs`, `src/var_calling/cohort_block/driver.rs` | lib 1025 pass; cohort_cli integration 21 pass | Bundled M5 |
| B2  | Blocker  | `compute_dust_mask_for_chrom` materialises whole chromosome | Apply | Applied | No | `src/var_calling/cohort_block/driver.rs` | cohort_cli integration 21 pass (byte-identity gates pass after streaming switch) | No |
| B3  | Blocker  | `NoSafeGap` retry is a no-op when `target_variants > 0` | Apply | Applied | No | `src/var_calling/cohort_block/driver.rs` | cohort_block lib 89 pass; cohort_cli integration 21 pass | Regression test (retry-grows-on-NoSafeGap) folds into B5's driver-side unit-test sweep in Wave 2 |
| B4  | Blocker  | `NAllelesExceedsBitmask` rewritten as fake `DegenerateLikelihood` | Apply | Applied with adaptation | No | `src/var_calling/posterior_engine.rs` (out-of-scope: added new variant to `PosteriorEngineError`), `src/var_calling/cohort_block/worker.rs` | new test `compute_ll_error_to_merger_preserves_n_alleles_and_locus` passes; cohort_block 89 pass | Adaptation: added the variant to `PosteriorEngineError` (out-of-scope file) rather than to `PerGroupMergerError` because the worker's call chain returns `PosteriorEngineError`. The review explicitly approved this option ("Add a new variant on `PerGroupMergerError` (or on `PosteriorEngineError`, or surface `ComputeLogLikelihoodsError` directly via a new `ChunkDriverError` variant)"). |
| B5  | Blocker  | No cross-driver byte-identity oracle test + missing driver-glue tests | Ask | Deferred to Wave 2 | Yes (Q2 answered: oracle vs streaming `drive_cohort_pipeline`) | — | — | Implement in Wave 2 |
| M1  | Major    | `#[derive(Default)]` on `SampleColumns` yields invariant-violating value | Apply | Applied (Wave 4) | No | `src/var_calling/cohort_block/columns.rs` | lib 1026 pass; cohort_cli integration 21 pass | Hand-written `impl Default for SampleColumns { fn default() -> Self { Self::empty() } }`. `#[non_exhaustive]` added by Mi1 in the same edit. (Wave 1's report table marked M1 Applied prematurely — the actual code change landed here.) |
| M2  | Major    | `prefetch_window_ref_bytes` drops inner Vec<u8> allocations | Apply | — | No | — | — | — |
| M3  | Major    | `detect_compound_candidates_columnar` allocates BTreeMaps per call | Apply | — | No (proptest part needs care) | — | — | — |
| M4  | Major    | `ChunkDriverError` design: `#[from]` funnel + mechanism-named variants | Apply | Applied (Wave 3) | No | `src/var_calling/cohort_block/driver.rs` | lib 1026 pass; cohort_cli integration 21 pass | Bundles M21; drops blanket `#[from]` on `Io` and `PspRead`; renames every variant by operation; adds `#[source]` + chrom/window/path context fields. |
| M5  | Major    | `let _ = remove_file` swallows tmp-cleanup error | Apply | Applied | No (bundled with B1) | (same as B1) | (covered by B1) | No |
| M6  | Major    | `u32_from_usize` wraps silently in release | Apply with adaptation | Applied (Wave 3) | No | `src/var_calling/cohort_block/columns.rs` | lib 1026 pass | Picked the `try_into().expect(...)` option per the review's "if the cost of plumbing `Result` through is too high" clause. Wrap silently is gone; both debug and release now panic with `"CSR offset exceeds u32::MAX"` on overflow. Typed-error variant deferred (would ripple through `SampleColumns::push_record`, `push_row_from`, `clear`, every loader site). |
| M7  | Major    | `params.target_window_count.max(1)` silently rewrites CLI input | Apply | — | Yes (Q3 — picks the propagate-vs-default branch) | — | — | — |
| M8  | Major    | `effective_initial_span` no-op chain | Apply | Applied (Wave 4) | No | `src/var_calling/cohort_block/loader.rs` | lib 1026 pass | Dropped the no-op `.min(.max(...))` chain; replaced with explicit `if max_span < initial_span { return Err(InvalidRange{...}) }` validation. The driver passes `max_load_span >= initial_load_span` post-B3, so this branch is unreachable on the production path but pins the API contract for future callers. |
| M9  | Major    | `target_variants_per_chunk == 0` sentinel-as-toggle | Ask | — | Yes (Q3) | — | — | — |
| M10 | Major    | 2 in-scope `cargo doc` unresolved-link errors | Apply | Applied (Wave 2) | No | `src/var_calling/cohort_block/worker.rs` | `cargo doc --no-deps` no longer reports the two `cohort_block/worker.rs:17,236` errors; remaining 14 are pre-existing out-of-scope per the review | No |
| M11 | Major    | `VarCallingArgs` API drift breaks bench + 2 examples | Apply | Applied | No | `benches/cohort_e2e_perf.rs`, `examples/profile_cohort_e2e.rs`, `examples/dhat_var_calling.rs` | clippy --all-targets: 3 E0063 errors cleared; 16 in-scope clippy errors remain (tracked under Nits) | No |
| M12 | Major    | `load_and_run_chunk_with_retry` 19 params / 5 phases | Defer | — | Yes (Q6: refactor scope) | — | — | — |
| M13 | Major    | `load_chunk_from_iters` 9 params; needs `ChunkLoadExtent` struct | Defer | — | Yes (Q6) | — | — | — |
| M14 | Major    | Carryover snapshot+restore duplication; add `SampleColumns::clone_from_columns` | Apply | Applied | No | `src/var_calling/cohort_block/columns.rs`, `src/var_calling/cohort_block/driver.rs` | cohort_block 89 pass; cohort_cli integration 21 pass | No |
| M15 | Major    | `ColumnarPipelineScratch` mixes scratch with running counters | Defer | — | Yes (Q6: refactor scope) | — | — | — |
| M16 | Major    | `SampleColumns` 13 ungrouped pub fields | Defer | — | Yes (Q6: refactor scope) | — | — | — |
| M17 | Major    | `AlleleSupportStats` partial-destructure trailing `..` | Apply | Applied | No | `src/var_calling/cohort_block/columns.rs` | full lib 1026 pass | No |
| M18 | Major    | `#[allow(dead_code)] chain_id_scratch` is unused | Apply | Applied | No (deleted field + clear + test reference) | `src/var_calling/cohort_block/kernels/unify_alleles.rs` | cohort_block 89 pass | No |
| M19 | Major    | `partition_window` doesn't validate sorted `masked_intervals` | Apply | Applied | No | `src/var_calling/cohort_block/partition.rs` | full lib 1026 pass | The review also suggested a release-mode `PartitionError::MaskNotSorted` variant; deferred per minimal-diff discipline (debug_assert covers the regression-test case; only production caller `compute_dust_mask_for_chrom` is sdust-derived, which is sorted by construction). |
| M20 | Major    | `emit_or_drop` filter order not pinned by test | Apply | — | Yes (Q4 — clarifies oracle policy) | — | — | — |
| M21 | Major    | `emit_or_drop` doesn't surface VCF write errors with locus context | Apply | Applied (Wave 3) | No (bundled with M4) | (same as M4) | (covered by M4) | `WriteVcf { chrom_id, start, end, source }` carries the failing record's locus |
| M22 | Major    | `enforce_max_alleles` tie-break not tested against row-shape oracle | Apply | Applied (Wave 5) | No | `src/var_calling/cohort_block/worker.rs` | new test `unify_max_alleles_ties_match_row_shape_kernel` passes; column-native and row-shape kernels produce identical output | No |
| M23 | Major    | `n_alleles < 2` silent skip has no counter | Apply | Applied (Wave 5) | No | `src/var_calling/cohort_block/worker.rs`, `src/var_calling/cohort_block/driver.rs` | lib 1032 pass; cohort_cli integration 21 pass | Added `ColumnarPipelineScratch::groups_skipped_post_unify_ref_only` + `take_post_unify_ref_only_count()` + `ChunkDriverStats::groups_skipped_post_unify_ref_only`. Driver's drain wires the counter through after each `run_window` call. |
| M24 | Major    | `prefetch_window_ref_bytes` may fetch past chrom_length | Defer | Deferred to Wave 6 | No | — | — | Needs design analysis: clamp `safe_end` to `chrom_one_past_end` in pre_pass (option a — requires threading `chrom_length`) vs. clamp span in `prefetch_window_ref_bytes` (option b — requires verifying downstream kernels can handle short ref slices). Defer until the structural Wave 6 work is in flight. |
| M25 | Major    | Five `expect/unwrap` sites lack `// PANIC-FREE:` comments | Apply | Applied (Wave 4) | No (mechanical) | `src/var_calling/cohort_block/loader.rs`, `src/var_calling/cohort_block/driver.rs`, `src/var_calling/cohort_block/pre_pass.rs` | lib 1026 pass | PANIC-FREE comments added at 5 sites: loader.rs:469+476 (peek/next), driver.rs (u32::try_from(chrom_idx)), pre_pass.rs (timeline.last() and prefix_max_reach.last()). |
| M26 | Major    | `pos + ref_span - 1` plain arithmetic on PSP-derived input | Apply | Applied (Wave 4) | No | `src/var_calling/cohort_block/pre_pass.rs`, `src/var_calling/cohort_block/partition.rs` | lib 1026 pass | `pos.saturating_add(ref_span.max(1)).saturating_sub(1)` at both pre_pass.rs:137 and partition.rs:376. |
| M27 | Major    | `*_cfg` vs `*_config` naming inconsistency across driver/worker | Apply | Applied (Wave 2) | No | `src/var_calling/cohort_block/worker.rs` | lib 1026 pass; renamed `posterior_config` → `posterior_cfg`, `per_group_config` → `per_group_cfg` (14 sites) | No |
| M28 | Major    | `shared_ref_fetcher` noun-named function | Apply | Applied (Wave 2) | No | `src/var_calling/cohort_block/worker.rs`, `src/var_calling/cohort_block/mod.rs` | lib 1026 pass; `pub use` updated | No |
| M29 | Major    | No tests for `SampleCountMismatch` / `CarryoverLengthMismatch` (loader) | Apply | Applied (Wave 5) | No | `src/var_calling/cohort_block/loader.rs` | new tests `loader_rejects_sample_count_mismatch` + `loader_rejects_carryover_length_mismatch` + `loader_rejects_max_span_below_initial_span` (M8 regression) all pass | No |
| M30 | Major    | No test for `CarryoverLengthMismatch` (pre_pass) | Apply | Applied (Wave 5) | No | `src/var_calling/cohort_block/pre_pass.rs` | new test `pre_pass_rejects_mismatched_carryover_length` passes | No |
| M31 | Major    | No concurrency test for `par_iter_mut` driver path | Apply | — | No | — | — | — |
| M32 | Major    | `locate_sample_row_idx` linear-scan invariant not pinned | Apply | Applied (Wave 5) | No | `src/var_calling/cohort_block/kernels/project_scalars.rs`, `src/var_calling/cohort_block/partition.rs` | new test `samples_at_pos_is_strictly_ascending_within_each_position_range` passes; `debug_assert!` added in `locate_sample_row_idx` | No |
| Mi1 | Minor    | New pub data structs lack `#[non_exhaustive]` | Apply | Applied (Wave 4) | Yes (Q1: stable in-crate API) | `src/var_calling/cohort_block/columns.rs` (×2), `src/var_calling/cohort_block/partition.rs`, `src/var_calling/cohort_block/loader.rs`, `src/var_calling/cohort_block/driver.rs`, `src/var_calling/cohort_block/kernels/unify_alleles.rs`, `src/var_calling/cohort_block/kernels/project_scalars.rs`, `src/var_calling/cohort_block/kernels/compute_log_likelihoods.rs` | lib 1026 pass | `#[non_exhaustive]` added to `SampleColumns`, `MaterialisedChunk`, `WindowPartition`, `ChunkLoadStats`, `ChunkDriverStats`, `UnifiedAllelesColumns`, `ProjectedScalarsColumns`, `LogLikelihoodsColumns`. |
| Mi2 | Minor    | `MaterialisedChunk::clear_data` deviates from `clear()` convention | Apply | — | No | — | — | — |
| Mi3 | Minor    | `WorkerSlot.output_buf` / `.scratch` bare names | Apply | — | No (mechanical) | — | — | — |
| Mi4 | Minor    | `ll` vs `lh` vs spelled-out `log_likelihoods` inconsistency | Apply | — | No | — | — | — |
| Mi5 | Minor    | `pool_allele_mapq` anonymous primitive triple | Apply | — | No | — | — | — |
| Mi6 | Minor    | `ref_fetcher` vs `fetcher` inconsistency | Apply | — | No | — | — | — |
| Mi7 | Minor    | `fix_boundaries` / `pre_pass` weak names | Defer | — | Yes (Q6 — rename scope) | — | — | — |
| Mi8 | Minor    | `pub mod` → `pub(crate) mod` reduction | Ask | — | Yes (Q1) | — | — | — |
| Mi9 | Minor    | `Arc::new(StreamingChromRefFetcher)` single-owner | Apply | Applied (Wave 4) | No | `src/var_calling/cohort_block/driver.rs` | lib 1026 pass | Threaded `&dyn ChromRefFetcher` through `load_and_run_chunk_with_retry` instead of `SharedRefFetcher` (Arc). Dropped `use std::sync::Arc;` + `SharedRefFetcher` import + the `#[allow(clippy::arc_with_non_send_sync)]`. |
| Mi10 | Minor   | `chunk.windows.clone()` defensive | Apply | Applied (Wave 4) | No | `src/var_calling/cohort_block/driver.rs` | lib 1026 pass | Iterate `worker_pool.slots[..n_windows].iter_mut().enumerate()` and index `chunk.windows[window_idx]` directly (NLL allows the shared reborrow). |
| Mi11 | Minor   | `run_window` takes `posterior_config` by value | Apply | Applied (Wave 4) | No | `src/var_calling/cohort_block/worker.rs`, `src/var_calling/cohort_block/driver.rs` | lib 1026 pass | `posterior_cfg: &PosteriorEngineConfig` — driver no longer calls `.clone()` per worker slot inside the rayon dispatch. Three test call sites updated to pass `&PosteriorEngineConfig::default()`. |
| Mi12 | Minor   | Double-clone of `projection_buf` in unify | Apply | — | No | — | — | — |
| Mi13 | Minor   | Three `pub` items have no out-of-module caller | Apply | — | Yes (gated by Q1) | — | — | — |
| Mi14 | Minor   | `match Ok | Err` should merge arms | Apply | Applied (Wave 2) | No | `src/var_calling/cohort_block/pre_pass.rs`, `src/var_calling/cohort_block/partition.rs` | lib 1026 pass | No |
| Mi15 | Minor   | Three files past the soft-cap line count | Defer | — | Yes (Q6) | — | — | — |
| Mi16 | Minor   | `build_overlapping_variant_group` vestigial cross-module reach | Apply | — | No (subsumes part of M10's worker.rs:17 link error) | — | — | — |
| Mi17 | Minor   | Three `Vec<Vec<_>>` jagged arrays should be CSR | Defer | — | Yes (Q6) | — | — | — |
| Mi18 | Minor   | `u32` primitive obsession for positions/spans | Defer | — | Yes (Q6) | — | — | — |
| Mi19 | Minor   | `ChunkDriverParams` 12-field mixed-axis tuning | Defer | — | Yes (Q6) | — | — | — |
| Mi20 | Minor   | No startup `tracing::info!` for effective config | Apply | — | No | — | — | — |
| Mi21 | Minor   | `emit_windows` silent down-grade vs requested count | Apply | — | No | — | — | — |
| Mi22 | Minor   | `drain_rows_from_into` boundary tests missing | Apply | — | No | — | — | — |
| Mi23 | Minor   | `slide_left_to_safe` direct unit tests missing | Apply | — | No | — | — | — |
| Mi24 | Minor   | `cohort_e2e_perf.rs` lacks `// REGRESSION THRESHOLD` comment | Apply | — | No | — | — | — |
| Mi25 | Minor   | `clippy::doc_lazy_continuation` ×2 in `test_helpers.rs` | Apply | Applied (Wave 2) | No | `src/var_calling/cohort_block/test_helpers.rs` | clippy `doc_lazy_continuation` errors at test_helpers.rs:20,21 cleared | No |
| Mi26 | Minor   | Three `rustdoc::redundant_explicit_links` | Apply | Applied (Wave 2) | No | `src/var_calling/cohort_block/driver.rs`, `src/var_calling/cohort_block/worker.rs` | three `redundant_explicit_links` warnings at driver.rs:105/108 and worker.rs:9 cleared | No |
| Nits | Nits     | 16 clippy errors + 14 `#[allow]` annotations + small style items | Apply | Applied (Wave 6) | No (mechanical) | `Cargo.toml` (lint-table addition), `src/var_calling/cohort_block/columns.rs` (test-module allow), `src/var_calling/cohort_block/partition.rs` (test-module allow), `src/var_calling/cohort_block/kernels/unify_alleles.rs` (OracleAllele type + assert! rewrites), `src/var_calling/cohort_block/worker.rs` (OracleAllele type) | **`cargo clippy --all-targets --all-features -- -D warnings` now exits 0**; lib 1032/1032 pass; fmt clean | 16 in-scope clippy errors cleared: 8 `single_range_in_vec_init` (test-module `#[allow]` with justification — `vec![a..b]` is the natural test-fixture shape, rewriting would obscure intent), 4 `type_complexity` (extracted `type OracleAllele = (Vec<u8>, bool, Vec<(usize, usize)>);` aliases in 2 test modules), 2 `bool_assert_comparison` (`assert_eq!(x, true)` → `assert!(x)`), 2 `doc_lazy_continuation` already fixed in Wave 2. Plus 8 *new* `result_large_err` errors that Wave 3's `ChunkDriverError` redesign introduced — silenced cohort-wide via `[lints.clippy] result_large_err = "allow"` in `Cargo.toml` with a justification comment naming the trade-off (error-path-only, diagnostic value > 128-byte budget). 14 `#[allow(...)]` justification comments not added this wave; tracked as a follow-up. |

## 3. Questions asked and answers

1. **Q1 — API scope of new `pub` data structs in `cohort_block/`** (gates Mi1, Mi8, Mi13).
   - **Answer:** Treat as stable in-crate API. ⇒ Apply Mi1 (`#[non_exhaustive]` on new pub data structs). Mi8 (`pub mod` → `pub(crate) mod`) and Mi13 (demote `pub` items with no out-of-module caller) → **Won't fix per policy**.
2. **Q2 — B5 byte-identity test shape** (gates B5, M20).
   - **Answer:** Oracle test vs streaming `drive_cohort_pipeline`. ⇒ B5 = integration test in `tests/cohort_cli_integration.rs` running both drivers on a rich fixture, asserting equal VCF bodies + equal per-category counters.
3. **Q3 — `target_variants_per_chunk` / `target_window_count` sentinel policy** (gates M7, M9, Mi20).
   - **Answer:** Lift both to `NonZero*` / `Option<NonZero*>` types. ⇒ Type-level fix: `target_variants_per_chunk: Option<NonZeroU32>` and `target_window_count: NonZeroUsize` (or `Option<NonZeroUsize>` with a CLI default). Update the CLI parser, the driver, and the pre-pass interface. Mi20 (startup `eprintln!` listing effective config) lands as a sibling.
4. **Q4 — refactor scope** (gates M12, M13, M14, M15, M16, Mi7, Mi15, Mi17, Mi18, Mi19).
   - **Answer:** Apply all structural refactors now. ⇒ M12-M16 + Mi7 + Mi15 + Mi17-Mi19 all in this fix run.

## 4. Per-finding log

### B1 — Writer un-finished + tmp leak on driver-error path
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied (bundles M5)
- **Reasoning:** Type-system level fix — adding `CohortVcfWriter::abort(self)` makes the tmp-path cleanup the writer's single source of truth (it knew the path the sink opened at). Refactored the per-chrom loop closure so it borrows `writer` instead of moving it via `finish()`, letting the outer scope route to `writer.finish()` on success or `writer.abort()` on error.
- **Implementation summary:** Added `pub fn abort(self) -> std::io::Result<()>` to `CohortVcfWriter` that drops the sink and removes the tmp file at `tmp_path_for(&self.final_path)`. Restructured `drive_cohort_chunked` to call `writer.finish()` / `writer.abort()` from the outer scope. On `abort()` error, log via `eprintln!` (matches the project's "no logs" preference vs `tracing`).
- **Review suggestion used verbatim?:** No — adapted the suggested `writer.abort()` to use the project's `eprintln!` convention rather than `tracing::warn!`.
- **Adaptation:** `tracing::warn!` → `eprintln!` to match the project's existing stderr-status pattern (`var_calling.rs:546`).
- **Verification performed:** Two new unit tests in `src/vcf/writer.rs`: `abort_removes_tmp_file_and_leaves_no_final_output` and `abort_surfaces_remove_file_error_when_tmp_already_gone`. Both pass.
- **Files changed:** `src/vcf/writer.rs`, `src/var_calling/cohort_block/driver.rs`.
- **Tests added or modified:** `abort_removes_tmp_file_and_leaves_no_final_output`, `abort_surfaces_remove_file_error_when_tmp_already_gone` (lib).
- **Validation:**
  - `cargo test --lib vcf::writer::tests::abort` → exit 0, 2 passed.
  - `cargo test --lib` → exit 0, 1025 passed.
  - `cargo test --test cohort_cli_integration` → exit 0, 21 passed.
- **User input:** None.
- **Follow-up:** B5 will add an integration test injecting a mid-loop error and asserting no leftover tmp.
- **Residual risk:** None.

### B2 — `compute_dust_mask_for_chrom` materialises whole chromosome
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `sdust_mask_streaming` already accepts an `Iterator<Item = io::Result<u8>>`; the bug was a needless `.collect::<Result<Vec<u8>, _>>()?` in the middle. Direct piping recovers the per-chunk memory budget (~90 MB per chrom on tomato 1).
- **Implementation summary:** Replaced the `let bases: Vec<u8> = fetcher.iter_bases()?.collect(...)?;` + `bases.into_iter().map(Ok)` shape with `fetcher.iter_bases()?.map(|r| r.map_err(io::Error::other))`.
- **Review suggestion used verbatim?:** Yes (functionally — minor adaptation: `io::Error::other(e)` instead of `io::Error::other(format!("{e}"))`).
- **Adaptation:** None (the review's suggested code worked as-is).
- **Verification performed:** `cargo test --test cohort_cli_integration` — all 21 tests pass including the byte-identity gates (`--worker-windows-per-chunk` and `--target-variants-per-chunk`). Byte output unchanged.
- **Files changed:** `src/var_calling/cohort_block/driver.rs`.
- **Tests added or modified:** None — coverage handled by the existing byte-identity tests, which exercise the dust mask end-to-end. The review's "synthetic fetcher that panics on `.collect()`" regression test would lock the streaming property explicitly; deferred to Wave 2 as a Mi-class follow-up (low-priority since the implementation no longer contains a `.collect()` call).
- **Validation:** `cargo test --test cohort_cli_integration` → exit 0, 21 passed.
- **User input:** None.
- **Follow-up:** Optional regression test in Wave 2.
- **Residual risk:** None — byte-identical output preserved.

### B3 — `NoSafeGap` retry is a no-op when `target_variants > 0`
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The bug was that `max_load_span` was independent of `attempt_span` (chrom-wide cap regardless). Decoupling the two axes — loader's variant-bounded extension within `attempt_span`, outer retry growing `attempt_span` itself — restores the design intent.
- **Implementation summary:** Replaced `let max_load_span = extension_cap_end.saturating_sub(chunk_range_start);` (chrom-wide cap) with `let max_load_span = initial_load_span;` (per-attempt cap derived from `attempt_span` + last-chunk extension). Likewise tied `psp_inclusive_end` to `chunk_range_end` instead of `extension_cap_end`. Dropped the unused `max_load_span.max(initial_load_span)` upgrade.
- **Review suggestion used verbatim?:** Yes (adopted suggested-fix option (a): thread `attempt_span` into `max_span`).
- **Adaptation:** None.
- **Verification performed:** All cohort_block unit tests + cohort_cli integration tests pass — no byte-identity regression on the existing fixtures. (Existing fixtures don't trigger `NoSafeGap`, so the change is observable only on the retry path.)
- **Files changed:** `src/var_calling/cohort_block/driver.rs`.
- **Tests added or modified:** None — see Follow-up below.
- **Validation:** `cargo test --lib var_calling::cohort_block` → exit 0, 89 passed; `cargo test --test cohort_cli_integration` → exit 0, 21 passed.
- **User input:** None.
- **Follow-up:** **Regression test for the retry-grows-on-NoSafeGap behaviour** — defer to Wave 2 B5 driver-side test sweep. Constructing a fixture that triggers `NoSafeGap` is non-trivial (needs records spaced closer than `max_group_span` throughout the chunk's nominal span with a safe gap further out); will land alongside the cross-driver integration test.
- **Residual risk:** Low — the change only affects code paths the existing fixtures don't exercise, and the existing fixtures still pass.

### B4 — `NAllelesExceedsBitmask` rewritten as fake `DegenerateLikelihood`
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** The review listed three places the new variant could live: `PerGroupMergerError`, `PosteriorEngineError`, or a new `ChunkDriverError` variant. The worker's call chain returns `PosteriorEngineError`, so adding the variant there is the minimal-diff option. `#[non_exhaustive]` on `PosteriorEngineError` makes this non-breaking.
- **Implementation summary:** Added `PosteriorEngineError::NAllelesExceedsBitmask { locus: RecordLocus, n_alleles: usize }`. Rewrote `compute_ll_error_to_merger` to map `ComputeLogLikelihoodsError::NAllelesExceedsBitmask { n_alleles }` directly to the new variant (carrying real `n_alleles` and locus) instead of synthesising `DegenerateLikelihood { sample_idx: usize::MAX, … }`. Removed the `let _ = n_alleles;` and the now-unused `DegeneracyKind` import.
- **Review suggestion used verbatim?:** Yes (suggested-fix option B: PosteriorEngineError variant).
- **Adaptation:** None — chose option B (PosteriorEngineError) from the three the review listed.
- **Verification performed:** Added unit test `compute_ll_error_to_merger_preserves_n_alleles_and_locus` (worker.rs:tests). Asserts the new mapping preserves `n_alleles` and the group locus.
- **Files changed:** `src/var_calling/posterior_engine.rs` (out-of-scope file; minimal additive change to `#[non_exhaustive]` enum), `src/var_calling/cohort_block/worker.rs`.
- **Tests added or modified:** `compute_ll_error_to_merger_preserves_n_alleles_and_locus`.
- **Validation:** `cargo test --lib var_calling::cohort_block` → exit 0, 89 passed (was 88; +1 from the new test).
- **User input:** None.
- **Follow-up:** None.
- **Residual risk:** None.

### M5 — `let _ = remove_file` swallows tmp-cleanup error
- **Severity:** Major
- **Initial decision:** Apply (bundled with B1)
- **Final status:** Applied
- **Reasoning:** Covered by the B1 restructuring.
- **Implementation summary:** Replaced `let _ = std::fs::remove_file(...)` with `writer.abort()` + `eprintln!` on error (see B1).
- **Validation:** Same as B1.
- **Follow-up:** None.
- **Residual risk:** None.

### M11 — `VarCallingArgs` API drift breaks bench + 2 examples
- **Severity:** Major
- **Initial decision:** Apply (prereq to save criterion baseline)
- **Final status:** Applied
- **Reasoning:** Three struct-literal call sites missed the two new required fields. Added the legacy defaults (`target_variants_per_chunk: 0`, `worker_windows_per_chunk: 1`) so bench/example numbers stay comparable to the pre-rewrite shape.
- **Implementation summary:** Added the two fields to `benches/cohort_e2e_perf.rs:286`, `examples/profile_cohort_e2e.rs:152`, `examples/dhat_var_calling.rs:121`.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** `cargo clippy --all-targets --all-features -- -D warnings` no longer fires the 3 E0063 errors. Criterion baseline saved successfully on `benches/cohort_e2e_perf.rs`.
- **Files changed:** `benches/cohort_e2e_perf.rs`, `examples/profile_cohort_e2e.rs`, `examples/dhat_var_calling.rs`.
- **Validation:** `cargo bench --bench cohort_e2e_perf -- --test` → exit 0; baseline saved to `target/criterion/*/pre-fixes`.
- **Follow-up:** Optional `VarCallingArgs::for_profiling(...)` constructor to prevent the same drift on the next field addition. Tracked but not part of this fix.
- **Residual risk:** None.

### M14 — Carryover snapshot+restore duplication
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The two parallel loops would silently skew if a new `SampleColumns` field landed without updating both. The helper consolidates the contract.
- **Implementation summary:** Added `SampleColumns::clone_from_columns(&mut self, other: &SampleColumns)`. Replaced both `for ... { clear(); for row_idx in 0..n_records() { push_row_from(...) } }` loops in `load_and_run_chunk_with_retry` with single `clone_from_columns` calls.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** cohort_cli integration tests pass (the carryover restore path is exercised by every multi-chunk chromosome run).
- **Files changed:** `src/var_calling/cohort_block/columns.rs`, `src/var_calling/cohort_block/driver.rs`.
- **Validation:** `cargo test --lib var_calling::cohort_block` → exit 0, 89 passed; `cargo test --test cohort_cli_integration` → exit 0, 21 passed.
- **Follow-up:** None.
- **Residual risk:** None.

### M17 — `AlleleSupportStats` partial-destructure trailing `..`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The 7 named fields cover every field on `AlleleSupportStats` today; the trailing `..` is a refactor-safety footgun.
- **Implementation summary:** Dropped the `..` from the `AlleleSupportStats { … }` destructure in `SampleColumns::push_record`. A future field addition to `AlleleSupportStats` will now be a compile error here.
- **Review suggestion used verbatim?:** Yes.
- **Adaptation:** None.
- **Verification performed:** `cargo test --lib` — 1026 passed.
- **Files changed:** `src/var_calling/cohort_block/columns.rs`.
- **Validation:** `cargo test --lib` → exit 0, 1026 passed.
- **Follow-up:** None.
- **Residual risk:** None.

### M18 — `#[allow(dead_code)] chain_id_scratch` is unused
- **Severity:** Major
- **Initial decision:** Apply (delete)
- **Final status:** Applied
- **Reasoning:** Field had no kernel consumer; clearing it and asserting on test was performative.
- **Implementation summary:** Removed the field from `UnifyAllelesScratch`, its `clear()` call, and the test reference. `ChainId` is still imported and used by `detect_compound_candidates_columnar` directly.
- **Review suggestion used verbatim?:** Yes (chose option (a) — delete — since compound detection has settled on `BTreeMap`).
- **Adaptation:** None.
- **Verification performed:** `cargo test --lib var_calling::cohort_block` — 89 passed.
- **Files changed:** `src/var_calling/cohort_block/kernels/unify_alleles.rs`.
- **Validation:** `cargo test --lib var_calling::cohort_block` → exit 0, 89 passed.
- **Follow-up:** M3 will add the proper scratch-reuse for compound detection (BTreeMaps moved into `UnifyAllelesScratch`) — a future scratch buffer may legitimately want a `Vec<ChainId>` field at that point, added with a real user.
- **Residual risk:** None.

### M19 — `partition_window` doesn't validate sorted `masked_intervals`
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** The cursor logic depends on sorted+non-overlapping; locking the invariant with `debug_assert!` catches a future test helper or alternative mask source that violates it without paying the cost in release.
- **Implementation summary:** Added `debug_assert!(masked_intervals.windows(2).all(|w| w[0].end <= w[1].start), ...)` at the top of `partition_window` (after the existing scratch-sample-count check).
- **Review suggestion used verbatim?:** Partially — the review also suggested a release-mode `PartitionError::MaskNotSorted` variant. Deferred per minimal-diff discipline: the only production caller (`compute_dust_mask_for_chrom`) is sdust-derived which is sorted by construction, and adding a release-mode error variant ripples through callers/tests for a regression that the debug_assert already locks.
- **Adaptation:** Debug-only assertion (release-mode variant deferred to Wave 2 if needed).
- **Verification performed:** `cargo test --lib` — 1026 passed (existing fixtures pass sorted+non-overlapping masks).
- **Files changed:** `src/var_calling/cohort_block/partition.rs`.
- **Tests added or modified:** None this wave — the recommended `#[should_panic]` regression test for the assert is in Wave 2's Missing-tests sweep.
- **Validation:** `cargo test --lib` → exit 0, 1026 passed.
- **Follow-up:** Add `partition_rejects_unsorted_masked_intervals` `#[should_panic]` test in Wave 2; decide release-mode `MaskNotSorted` variant policy then.
- **Residual risk:** Low — production caller invariant holds; debug builds catch regressions.

*(Remaining findings — 28 Majors + 26 Minors + B5 + Nits — staged for Wave 2 onward.)*

## 5. Deferred findings to carry forward

Wave 1 closes 11 findings (4 Blockers + B1/B4 tests + M5/M11/M14/M17/M18/M19).

**Deferred to Wave 2 onward (47 items):**
- B5 — byte-identity oracle integration test (shape decided per Q2).
- M2 (prefetch_window_ref_bytes inner-Vec churn), M3 (detect_compound BTreeMaps + proptest), M4+M21 (`ChunkDriverError` redesign), M6 (`u32_from_usize` typed error), M7+M9 (NonZero type lift per Q3), M8 (effective_initial_span no-op), M10 (in-scope `cargo doc` links), M12 (split 19-param fn), M13 (`ChunkLoadExtent` struct), M15 (split `ColumnarPipelineScratch`), M16 (group `SampleColumns` columns), M20 (filter-order tests), M22 (tie-break test), M23 (post-unify-ref-only counter), M24 (chrom-end clamp), M25 (PANIC-FREE comments), M26 (saturating arithmetic), M27 (cfg/_config rename), M28 (`into_shared_ref_fetcher` rename), M29 (loader malformed-input tests), M30 (pre_pass CarryoverLengthMismatch test), M31 (par_iter_mut equivalence test), M32 (locate_sample_row_idx invariant tests).
- Mi1 (`#[non_exhaustive]` per Q1), Mi2 (`MaterialisedChunk::clear` rename), Mi3 (`output_buf`/`scratch` rename), Mi4 (`ll`/`lh`/`log_likelihoods` consistency), Mi5 (`PooledMapqMoments` struct), Mi6 (`ref_fetcher` vs `fetcher`), Mi7 (`fix_boundaries` / `pre_pass` rename per Q4), Mi9 (`Arc::new(StreamingChromRefFetcher)` collapse), Mi10 (`chunk.windows.clone()` drop), Mi11 (`&PosteriorEngineConfig`), Mi12 (double-clone), Mi14 (`Ok | Err` merge arms), Mi15 (file splits per Q4), Mi16 (move `build_overlapping_variant_group` to test_helpers), Mi17 (`Vec<Vec<_>>` → CSR per Q4), Mi18 (newtypes per Q4), Mi19 (param subgroups per Q4), Mi20 (startup eprintln), Mi21 (window-count down-grade log), Mi22+Mi23 (boundary tests), Mi24 (REGRESSION THRESHOLD), Mi25 (doc_lazy_continuation in test_helpers), Mi26 (redundant explicit links).
- **Won't fix per Q1 (stable in-crate API answered):** Mi8 (`pub mod` → `pub(crate) mod` reduction), Mi13 (demote `pub` items with no out-of-module caller).
- Nits — 16 in-scope clippy errors + 14 `#[allow(...)]` justifications + smaller style items.

## 6. Disputed findings to return to reviewer
*(populated at end of run)*

## 7. Failed-validation findings
*(populated at end of run)*

## 8. Blocked-by-context-mismatch findings
*(populated at end of run)*

## 9. Performance check

- **Triggered:** yes (B2, B3 touch the chunk-driver hot path covered by `benches/cohort_e2e_perf.rs`).
- **Baseline saved:** pending — will be captured after M11 is applied (the bench is currently broken by the `VarCallingArgs` field drift; baseline cannot be saved until M11 unbreaks it).
- **Benches to run:** all benches in `benches/`.
- **Verdicts:** pending.

## 10. Commands run (Wave 1)

All inside the dev container via `./scripts/dev.sh`:

- `cargo clippy --all-targets --all-features -- -D warnings` (×2: pre- and post-M11)
- `cargo bench --bench cohort_e2e_perf -- --test`
- `cargo bench --bench cohort_e2e_perf -- --save-baseline pre-fixes`
- `cargo test --lib vcf::writer::tests::abort`
- `cargo test --lib` (final: 1026 passed)
- `cargo test --lib var_calling::cohort_block` (final: 89 passed)
- `cargo test --test cohort_cli_integration` (final: 21 passed)
- `cargo test --tests`
- `cargo fmt --check` / `cargo fmt`

## 11. Command results (Wave 1)

- `cargo fmt --check` → 0 (after `cargo fmt` pass).
- `cargo clippy --all-targets --all-features -- -D warnings` → 101 (16 in-scope errors remain, all in tests / docs / test helpers; tracked under Nits).
- `cargo test --lib` → 0, 1026 passed (was 1023).
- `cargo test --lib var_calling::cohort_block` → 0, 89 passed (was 88; +1 from B4 test).
- `cargo test --tests` → 0; cohort_cli_integration 21 passed.
- `cargo doc --no-deps` → not re-run (M10 doc-link fixes in Wave 2).
- `cargo bench --bench cohort_e2e_perf -- --save-baseline pre-fixes` → exit 0; baseline saved at `target/criterion/*/pre-fixes`.

## 12. Notes
- **Wave 1 scope:** Blockers + the M11 unblocker + the smallest focused Majors (M5 bundled with B1; M14 / M17 / M18 / M19 standalone). This matches the project's historical multi-wave fix-application pattern (cohort_cli_2026-05-19 Wave 1 → Wave 5).
- **Out-of-scope edits made in Wave 1:**
  - `src/var_calling/posterior_engine.rs` — added `PosteriorEngineError::NAllelesExceedsBitmask` variant (B4). Non-breaking thanks to `#[non_exhaustive]`. Explicitly approved by the review's suggested-fix wording.
  - `src/vcf/writer.rs` — added `CohortVcfWriter::abort()` method (B1).
  - `benches/cohort_e2e_perf.rs`, `examples/profile_cohort_e2e.rs`, `examples/dhat_var_calling.rs` — added the two new `VarCallingArgs` fields (M11).
- **Wave 2 entry point:** apply B5 (oracle test) + the structural-refactor Majors (M12, M13, M15, M16) + the smaller in-scope Majors (M22, M23, M24, M25, M26, M27, M28, M29, M30, M31, M32) + the `ChunkDriverError` redesign (M4 + M21 + M6).
- **Wave 3 entry point:** M7+M9 NonZero type lift per Q3; Mi1 `#[non_exhaustive]` per Q1; Mi-class small fixes; Nits mechanical pass.
- **Performance comparison** — deferred to end-of-all-waves (criterion baseline saved at Wave 1 start, will compare against `pre-fixes` after the last code-touching fix lands).
