# Code Review: cohort_block
**Date:** 2026-05-29
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** PR-equivalent review of the new `src/var_calling/cohort_block/` module (12 files, ~8 400 LoC) on branch `cohort-within-chromosome-parallel` at commit `36989d6`.
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the chunk-based rewrite of the cohort var-calling pipeline, as a PR-equivalent diff against `main`.
- **Reviewed against:** branch `cohort-within-chromosome-parallel`, commit `36989d6f53460042c5219d0ff5fa6d67a7b1b129`.
- **In-scope files:**
  - [mod.rs](../../../src/var_calling/cohort_block/mod.rs) (41 LoC)
  - [columns.rs](../../../src/var_calling/cohort_block/columns.rs) (601 LoC)
  - [loader.rs](../../../src/var_calling/cohort_block/loader.rs) (954 LoC)
  - [pre_pass.rs](../../../src/var_calling/cohort_block/pre_pass.rs) (628 LoC)
  - [partition.rs](../../../src/var_calling/cohort_block/partition.rs) (824 LoC)
  - [driver.rs](../../../src/var_calling/cohort_block/driver.rs) (650 LoC)
  - [worker.rs](../../../src/var_calling/cohort_block/worker.rs) (1166 LoC)
  - [test_helpers.rs](../../../src/var_calling/cohort_block/test_helpers.rs) (135 LoC, `#[cfg(test)]`)
  - [kernels/mod.rs](../../../src/var_calling/cohort_block/kernels/mod.rs) (24 LoC)
  - [kernels/unify_alleles.rs](../../../src/var_calling/cohort_block/kernels/unify_alleles.rs) (1766 LoC)
  - [kernels/project_scalars.rs](../../../src/var_calling/cohort_block/kernels/project_scalars.rs) (921 LoC)
  - [kernels/compute_log_likelihoods.rs](../../../src/var_calling/cohort_block/kernels/compute_log_likelihoods.rs) (720 LoC)
- **Deliberately out of scope:**
  - `var_calling/per_group_merger.rs`, `per_position_merger.rs`, `posterior_engine.rs` — row-shape modules the columnar pipeline calls into for type definitions only.
  - `var_calling/dust_filter.rs`, `var_calling/variant_grouping.rs` — narrow-API consumers.
  - The 17 pre-existing rustup 1.94→1.95 clippy errors outside `cohort_block/` (Phase A impl report deferred-work list).
- **Categories dispatched:**
  - `reliability` — always.
  - `errors` — always.
  - `naming` — always.
  - `defaults` — scope contains a public `ChunkDriverParams` + many `Default`-deriving structs.
  - `idiomatic` — always.
  - `refactor_safety` — always.
  - `module_structure` — multi-file module + nested `kernels/`.
  - `unsafe_concurrency` — rayon dispatch on shared `&chunk` + `Send`-bound choices for the fetcher.
  - `smells` — always.
  - `tooling` — newly-introduced `#[allow(...)]` annotations + broken bench/example call sites.
  - `extras` — hot path, stable byte-identical output, public-API surface, diff-vs-plan check.

## 2. Verdict

**Request-changes.** Five Blocker findings (one parallel-driver byte-identity test gap, one writer-resource leak on the error path, one `O(chrom_length)` materialisation that defeats the rewrite's memory contract, one retry that is silently a no-op when the operator-recommended knob is on, and one error that gets silently rewritten as a different error kind with `usize::MAX` placeholder fields). The clippy/doc gates also fail (Major, not Blocker — does not affect correctness, but blocks CI). Architecture is sound: parallel-section soundness is statically enforced by the `Send + !Sync` typedef and the `par_iter_mut`/`try_for_each` shape is clean (`unsafe_concurrency` returned no findings).

## 3. Execution status

All verification commands ran inside the dev container via `./scripts/dev.sh`. Verbatim output saved at [tmp/review_2026-05-29_cohort_block/verification.md](../../../tmp/review_2026-05-29_cohort_block/verification.md).

- `cargo fmt --check` — exit 0; clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — **exit 101**; 16 new clippy errors in `cohort_block/` (8× `single_range_in_vec_init`, 4× `type_complexity`, 2× `bool_assert_comparison`, 2× `doc_lazy_continuation`) + 3 `E0063` compile errors in `benches/cohort_e2e_perf.rs:286`, `examples/profile_cohort_e2e.rs:152`, `examples/dhat_var_calling.rs:121` from added required fields on `VarCallingArgs`.
- `cargo test --lib` — exit 0; **1023 passed** (88 of which are in `var_calling::cohort_block::`).
- `cargo test --tests` — exit 0; all integration binaries pass.
- `cargo doc --no-deps` — **exit 101**; 17 `rustdoc::broken_intra_doc_links` errors crate-wide, **2 in scope** ([worker.rs:17](../../../src/var_calling/cohort_block/worker.rs#L17) `build_overlapping_variant_group`; [worker.rs:236](../../../src/var_calling/cohort_block/worker.rs#L236) `PosteriorEngine`); 4 `redundant_explicit_links` warnings crate-wide, **3 in scope** ([driver.rs:105](../../../src/var_calling/cohort_block/driver.rs#L105), [driver.rs:108](../../../src/var_calling/cohort_block/driver.rs#L108), [worker.rs:9](../../../src/var_calling/cohort_block/worker.rs#L9)).
- `cargo audit` — **not run**: `error: no such command: audit` (cargo-audit not installed in the container).

Findings labeled "Needs verification": 0. Findings at Medium confidence: 9.

## 4. Open questions and assumptions

The author should resolve these before responding to individual findings.

1. **Is the in-module pub surface (`SampleColumns`, `MaterialisedChunk`, `WindowPartition`, `ChunkLoadStats`, `ChunkDriverStats`, `UnifiedAllelesColumns`, `ProjectedScalarsColumns`, `LogLikelihoodsColumns`, the `kernels::*` submodules) intended as a stable external API, or strictly internal scaffolding that downstream `pop_var_caller/` consumes?** — Gates Mi3 (`#[non_exhaustive]` on data structs) and Mi8 (`pub mod` → `pub(crate) mod` reduction). The fact that all error enums in the same scope *are* `#[non_exhaustive]` and `ChunkDriverParams` is `#[non_exhaustive]` argues "yes-stable"; the absence of out-of-crate consumers argues "internal-only".

2. **Is the streaming driver (`drive_cohort_pipeline` in `src/pop_var_caller/cohort_driver.rs`) staying around long-term as the byte-identity oracle, or scheduled for removal once `drive_cohort_chunked` ships?** — Gates the *shape* of B5's missing test: "vs streaming driver oracle in CI" vs "snapshot test against a checked-in golden VCF". If the streaming driver stays, the oracle approach is cheaper; if it goes, the golden-snapshot approach is the only durable option.

3. **Is `target_variants_per_chunk == 0` ("variant-bounded extension disabled") a long-lived operator choice, or only a transition mode while Phase B stabilises?** — Gates M9 (sentinel-as-toggle): if long-lived, the toggle should be lifted into the type system (`Option<NonZeroU32>`); if transitional, a `pub const TARGET_VARIANTS_DISABLED: u32 = 0` named constant + a debug-log on startup is enough.

4. **Filter-order divergence between drivers — is it byte-equivalent by construction?** — The streaming `drive_cohort_pipeline` applies `min_alt_obs_per_sample` *before* `PosteriorEngine` (between merger and EM); `emit_or_drop` applies it *after* the EM (driver.rs:271–275 acknowledges this). Per-category counters must still match byte-for-byte under every parameter combination, including `min_alt_obs = 0`, the LH-cap skip (`Ok(None)` at worker.rs:351), and the `n_alleles < 2` skip (worker.rs:354). Gates B5, M19, M20.

## 5. Top 3 priorities

1. **B5 — Add a byte-identity oracle integration test for the new driver.** Three of the other four Blockers (B1 partial-VCF leak, B2 whole-chrom alloc, B3 no-op retry) would be caught by a single rich-enough fixture run through both drivers with per-category counter assertions. Without it the only safety net is the 3-tomato manual run, which the impl report quotes but CI does not enforce.
2. **B3 — Fix the `NoSafeGap` retry no-op.** Today `load_and_run_chunk_with_retry`'s "double the span" path does no useful work whenever `target_variants_per_chunk > 0` (the operator-recommended mode) because `max_load_span` is the chromosome-wide cap regardless of `attempt_span`. Pathological real-data inputs that the retry was designed to recover from will surface as hard `NoSafeGap` failures. Find at [driver.rs:521-583](../../../src/var_calling/cohort_block/driver.rs#L521-L583).
3. **B2 — Stream the FASTA bases into `sdust_mask_streaming` instead of materialising the whole chromosome.** `compute_dust_mask_for_chrom` currently does `fetcher.iter_bases()?.collect::<Result<Vec<u8>, _>>()?` then re-streams the result. On tomato chrom 1 (~90 Mb) that's a 90 MB allocation per chromosome — defeating the rewrite's "per-chunk memory bound by `target_variants_per_chunk × n_samples`" contract. Find at [driver.rs:348-373](../../../src/var_calling/cohort_block/driver.rs#L348-L373).

## 6. Findings

### Blocker

- `src/var_calling/cohort_block/driver.rs:206-226` — **[Blocker]** B1: Driver-error path leaves the `CohortVcfWriter` un-`finish`ed and the partial tmp file behind
- **Categories:** reliability
- **Confidence:** High
- **Problem:** On any `Err` from the per-chrom loop the driver does `let _ = std::fs::remove_file(tmp_path_for(output));` and returns. The writer is dropped without `finish()` (or any abort), and the `tmp_path_for(output)` deletion *assumes* it names the same path the writer used internally — there is no debug-mode assertion linking the two. There is no test covering the failure path. A mid-run failure can silently leave a partial uncompressed/compressed VCF that looks indistinguishable from a successful output.
- **Why it matters:** Confusable partial output is the canonical "wrong result without a panic" failure mode. The 3-tomato byte-identity contract does not exercise this path.
- **Suggested fix:** Make `tmp_path_for` the single source of truth via an explicit `writer.abort()` method that takes the tmp path it actually used and removes it; call it from the error branch instead of the bare `remove_file`. Add an integration test that injects an error mid-loop (e.g. a PSP reader that errors on its second `region_records` call) and asserts no leftover tmp on disk.
  ```rust
  if let Err(e) = driver_result {
      let _ = writer.abort(); // best-effort, deletes the exact tmp path
      return Err(e);
  }
  ```

- `src/var_calling/cohort_block/driver.rs:348-373` — **[Blocker]** B2: `compute_dust_mask_for_chrom` materialises the whole chromosome's REF bases into a `Vec<u8>` before "streaming" to sdust
- **Categories:** reliability, smells (cross-cat), extras (cross-cat), errors (cross-cat)
- **Confidence:** High
- **Problem:** The mask helper does `let bases: Vec<u8> = fetcher.iter_bases()?.collect::<Result<Vec<u8>, _>>()?;`, then passes `bases.into_iter().map(Ok::<_, io::Error>)` to `sdust_mask_streaming`. On tomato chrom 1 (~90 Mb) that's a 90 MB allocation per chromosome — alive for the duration of the per-chrom call. The driver's own header comment promises *"one chunk × N samples, not T_chrom × N samples"*. The actual operating profile is `O(chunk × N_samples) + O(chrom_length)` per chromosome.
- **Why it matters:** Phase B's memory budget is the rewrite's headline win; this path nullifies it for small-N cohorts and constrains scaling for larger genomes (tomato max ~90 Mb is small; mammalian or plant chromosomes 200–300 Mb each are the next-up case). Byte-identity is preserved (the mask intervals are unchanged) but the design promise is broken silently.
- **Suggested fix:** Stream directly — `sdust_mask_streaming` already takes an `Iterator<Item = io::Result<u8>>`.
  ```rust
  let intervals = sdust_mask_streaming(
      fetcher.iter_bases()?.map(|r| r.map_err(io::Error::other)),
      chrom_length,
      dust_cfg.window(),
      dust_cfg.threshold(),
  )?;
  ```
  Add a regression test using a synthetic `ChromRefFetcher` whose `iter_bases` panics if `.collect::<Vec<_>>()`-ed, to lock the streaming property.

- `src/var_calling/cohort_block/driver.rs:483-583` — **[Blocker]** B3: `load_and_run_chunk_with_retry`'s "double the span" retry is a no-op when `target_variants_per_chunk > 0`
- **Categories:** reliability
- **Confidence:** High
- **Problem:** On `NoSafeGap` the helper doubles `attempt_span` and retries. But every retry calls `load_chunk_from_iters(..., initial_span=attempt_span, target_variants, max_span=max_load_span, ...)` where `max_load_span = extension_cap_end - chunk_range_start` is the chromosome-end cap, **independent of `attempt_span`**. When `target_variants_per_chunk > 0` the loader's internal extension runs up to `max_load_span` on the first attempt, so the chunk's actual span is fixed by `target_variants_per_chunk` and the data, not by `attempt_span`. The retry produces an identical chunk → identical `NoSafeGap`. The loop terminates because `attempt_span < max_span` eventually fails — but it does zero useful work.
- **Why it matters:** The operator-recommended mode (`target_variants_per_chunk > 0`) is exactly the one where the retry doesn't help. Dense variant clusters that need a wider chunk to find a safe gap will surface as hard failures even though a wider load could have processed them. The 3-tomato fixture does not exercise this path.
- **Suggested fix:** Thread `attempt_span` into the loader's `max_span` so each retry can actually load more data, or push the retry inside the loader as "keep extending until a safe gap is found or chrom cap hit". Add a test that constructs a chunk where the first attempt satisfies `target_variants_per_chunk` but the resulting chunk has no safe gap, and assert the retry actually grows the chunk.
  ```rust
  let max_load_span = attempt_span.min(extension_cap_end.saturating_sub(chunk_range_start));
  // …or pass attempt_span as both initial_span and max_span on retry.
  ```

- `src/var_calling/cohort_block/worker.rs:571,586` — **[Blocker]** B4: `compute_ll_error_to_merger` silently rewrites `NAllelesExceedsBitmask { n_alleles }` as `DegenerateLikelihood { sample_idx: usize::MAX, genotype_idx: usize::MAX, kind: NaN }` and discards `n_alleles` via `let _ = n_alleles;`
- **Categories:** errors, reliability (M-rel-5, convergent), smells (M-smell-4, convergent)
- **Confidence:** High
- **Problem:** `ComputeLogLikelihoodsError::NAllelesExceedsBitmask` is an **upstream-invariant-break** condition (the unifier's cap was misconfigured: `cfg.max_alleles > MAX_BITMASK_ALLELES`). The conversion at [worker.rs:571-595](../../../src/var_calling/cohort_block/worker.rs#L571-L595) drops the only context field (`n_alleles`) via `let _ = n_alleles;` and synthesises a `DegenerateLikelihood` with `usize::MAX` placeholder sample/genotype indices. The doc comment acknowledges the placeholder — but the result is a *misleading* operator-facing diagnostic that points at an impossible sample/genotype and hides the actual fault category. Three categories (errors, reliability, smells) flagged it independently.
- **Why it matters:** When this fires on real input, the operator sees `sample usize::MAX, genotype usize::MAX, NaN` — internally contradictory, no `n_alleles` value, no signal of the actual misconfiguration. Reliability rule: errors must carry the context the caller needs to debug; rewriting an internal-bug condition as a math-degenerate condition violates "errors are typed".
- **Suggested fix:** Add a dedicated variant carrying the lost context.
  ```rust
  // PerGroupMergerError (or a new ChunkDriverError variant):
  #[error("n_alleles={n_alleles} exceeds bitmask cap at {chrom_id}:{group_start}..{group_end}; \
           max_alleles_lh_calc misconfigured")]
  NAllelesExceedsBitmask {
      chrom_id: u32, group_start: u32, group_end: u32, n_alleles: usize,
  },
  ```
  Map directly into that variant. Drop the `let _ = n_alleles;`. Add a test that triggers the path (`PerGroupMergerConfig::new(2, MAX_BITMASK_ALLELES + 1, 64, 32)`) and asserts the new variant.

- Multiple files — **[Blocker]** B5: No in-tree byte-identity test for the new driver vs the streaming driver oracle, and no unit tests for `drive_cohort_chunked` / `drive_one_chrom_generic` / `load_and_run_chunk_with_retry` / `emit_or_drop` / `compute_dust_mask_for_chrom` / `passes_min_alt_obs` / `record_fails_mapq_diff_t`
- **Categories:** reliability, extras (M-ext-1, convergent)
- **Confidence:** High
- **Problem:** The chunk-driver glue carries the integration risk for the whole module. Today its only end-to-end coverage is the **external** 3-tomato benchmark quoted in the impl report (records_emitted=8329, hom_ref=8308, low_qual=1596, low_alt_obs=233990, low_mapq_diff_t=1134). The 88 in-tree tests under `var_calling::cohort_block::*` exercise sub-modules in isolation but never the assembled driver. Three existing integration tests (`var_calling_emits_deterministic_vcf_across_runs`, `var_calling_byte_identical_across_worker_windows_per_chunk`, `var_calling_byte_identical_across_target_variants_per_chunk`) check **run-vs-run** determinism and **knob-vs-knob** invariance, but **not** byte-identity vs the still-extant streaming `drive_cohort_pipeline` oracle on a fixture rich enough to exercise every drop category. Filter order even diverges between drivers (streaming: `min_alt_obs` pre-EM; chunked: `min_alt_obs` post-EM; see Open Question 4) — the per-category counters must remain byte-equal but no checked-in test pins it.
- **Why it matters:** The task brief calls byte-identical VCFs vs `main` the **hard correctness contract**. As soon as someone touches `emit_or_drop`, `build_posterior_record_columnar`, or the kernels, only the manual 3-tomato run would catch a regression. Adding an integration test now (while both drivers exist) costs the team one fixture + one assertion file; deferring it to "after streaming driver is removed" makes the test impossible because the oracle is gone.
- **Suggested fix:** Add one integration test in [tests/cohort_cli_integration.rs](../../../tests/cohort_cli_integration.rs) (or a new file) that (1) builds a multi-position fixture covering ≥3 samples, ≥1 MNP, ≥1 multi-allelic site forcing the LH cap, ≥1 hom-REF-only group, ≥1 site below `min_alt_obs`, ≥1 site failing `qual_phred`, ≥1 site failing `mapq_diff_t`; (2) runs both `drive_cohort_chunked` and `drive_cohort_pipeline` on identical input; (3) asserts equal VCF bodies (stripping `##source=` / `##commandline=`) **and** field-by-field equal counters (`records_written`, `records_dropped_hom_ref`, `records_dropped_low_qual`, `records_dropped_low_alt_obs`, `records_dropped_low_mapq_diff_t`, `records_unconverged`, `lh_cap_groups_skipped`, `lh_cap_alleles_in_skipped`). Both stat structs already have compatible shapes (per the comment at [driver.rs:117-120](../../../src/var_calling/cohort_block/driver.rs#L117-L120)). The §8 "Missing tests to add now" section enumerates the unit-level gaps that fall out (per-category-counter assertions, the chrom-end clamp, the retry-wider-on-NoSafeGap regression, the par_iter_mut vs sequential equivalence test).

### Major

- `src/var_calling/cohort_block/columns.rs:47` — **[Major]** M1: `#[derive(Default)]` on `SampleColumns` produces a structurally invalid value
- **Categories:** defaults, idiomatic — convergent
- **Confidence:** High
- **Problem:** `SampleColumns` derives `Default` (columns.rs:47) and also exposes a hand-written `empty()` (columns.rs:68) that seeds the three CSR offset columns with their `vec![0]` sentinel so `len(offsets) == n + 1` holds for `n = 0`. `Default::default()` produces empty `Vec`s without the sentinel; calling `n_alleles_total()` on a defaulted value panics ([columns.rs:96](../../../src/var_calling/cohort_block/columns.rs#L96): `expect("CSR sentinel always present")`), and `push_record` paths underflow when reading `allele_offsets[record_idx + 1]`. The loader's own code at [loader.rs:229](../../../src/var_calling/cohort_block/loader.rs#L229) uses `SampleColumns::empty` precisely because the author knows `Default` is wrong; sibling types `WindowPartition::default` and `UnifiedAllelesColumns::default` delegate to `Self::empty()` for the same reason.
- **Why it matters:** Two ways to construct one type, one of which silently produces a broken state. Future `..Default::default()` or `#[derive(Default)]` on a wrapping struct picks up the broken shape. Canonical "make illegal states unrepresentable" failure.
- **Suggested fix:** Remove `Default` from the derive and add a hand-written impl that delegates to `empty()`.
  ```rust
  #[derive(Debug, Clone, PartialEq)]
  pub struct SampleColumns { /* … */ }
  impl Default for SampleColumns { fn default() -> Self { Self::empty() } }
  ```

- `src/var_calling/cohort_block/worker.rs:284` — **[Major]** M2: `prefetch_window_ref_bytes` drops every per-group `Vec<u8>` allocation on every call
- **Categories:** idiomatic, reliability (Minor convergence)
- **Confidence:** High
- **Problem:** The function does `out.clear()` on `Vec<Vec<u8>>` — which drops every inner `Vec<u8>` and frees its capacity — then `Vec::with_capacity(span as usize)` for each group. The `WorkerSlot.pre_fetched_ref_bytes` field is documented at [worker.rs:164-166](../../../src/var_calling/cohort_block/worker.rs#L164-L166) as a reusable per-slot buffer, but the outer-clear pattern defeats that. With one slot per worker × `n_groups` per window × chunks-per-chrom, this is thousands to millions of `Vec::with_capacity` + `drop` pairs per run. This is the canonical "load / use / clear / reload" anti-pattern the project memory ([feedback_scratch_memory](../../../.claude/projects/-Users-jose-devel-pop-var-caller/memory/feedback_scratch_memory.md)) names.
- **Why it matters:** Allocator churn on the hottest path of the rewrite, in a function that was specifically designed to amortize across chunks. Directly contradicts the module-wide scratch-reuse design.
- **Suggested fix:** Resize the outer Vec without dropping inner Vecs, then clear each inner Vec in place.
  ```rust
  let n_groups = partition.n_groups();
  if out.len() < n_groups {
      out.resize_with(n_groups, Vec::new);
  } else {
      out.truncate(n_groups);
  }
  for (g, bytes) in out.iter_mut().enumerate().take(n_groups) {
      bytes.clear();
      let group_start = partition.group_starts[g];
      let span = partition.group_ends[g] - group_start + 1;
      bytes.reserve(span as usize);
      ref_fetcher.fetch_into(group_start, span, bytes)?;
  }
  ```

- `src/var_calling/cohort_block/kernels/unify_alleles.rs:770` — **[Major]** M3: `detect_compound_candidates_columnar` allocates two `BTreeMap`s plus per-key `Vec`s on every group
- **Categories:** idiomatic, reliability — convergent
- **Confidence:** High
- **Problem:** `detect_compound_candidates_columnar` (lines 764–798) and its helper `build_chain_proposals_columnar` allocate a fresh `BTreeMap<Vec<(usize,usize)>, CompoundCandidate>` and a `BTreeMap<ChainId, Vec<CompoundConstituent>>` on every call. `by_chain.clear()` runs per-sample but the maps, their per-key `Vec<(usize,usize)>` keys, and the per-key `Vec<CompoundConstituent>` values all allocate per call. `UnifyAllelesScratch` exists to amortize per-group buffers; these two are exactly the kind of state it was built for. The row-shape merger had the same churn — this rewrite was meant to fix it.
- **Why it matters:** Compound detection is on the hot path of every variant group with chain ids. `BTreeMap` insertion allocates per key; with N samples × C candidates per group × G groups per chunk, this is the next-biggest allocator caller after PSP's `decode_bytes_split` per recent perf profiling. Also: no proptest covers the order-invariance law of the compound-aggregation algorithm (different sample / chain-id orderings must produce the same candidate set).
- **Suggested fix:** Move both maps and their per-key buffers into `UnifyAllelesScratch`; pass `&mut UnifyAllelesScratch` into the helper; clear in place. Add a proptest that shuffles sample/chain-id orderings and asserts the candidate set is invariant.

- `src/var_calling/cohort_block/driver.rs:143-162` — **[Major]** M4: `ChunkDriverError` design — variant names describe mechanism not operation, `#[from]` funnels multiple origins through single variants (esp. `Io` and `PspRead`), `#[error("…: {0}")]` flattens inner causes into the parent message, public `Io(#[from] io::Error)` exposes a foreign type
- **Categories:** errors
- **Confidence:** High
- **Problem:** The enum has 9 variants — `ChunkLoad`, `FixBoundaries`, `Partition`, `Posterior`, `PerGroupMerger`, `RefFetch`, `VcfWrite`, `PspRead`, `Io`. `PspReadError` is reachable directly (`PspRead`) **and** inside `ChunkLoadError<PspReadError>` (`ChunkLoad`); `io::Error` is reachable directly (`Io`) **and** inside `ChromRefFetchError` (`RefFetch`) **and** inside `VcfWriteError`. Every `?` over an `io::Result` (file open, fs::remove_file, reader read) collapses into `Io` regardless of what the driver was attempting. Variant names describe subsystem (`Posterior`, `RefFetch`) rather than operation (`OpenPsp`, `FetchRefBases`, `WriteVcf`). The `: {0}` interpolation flattens the inner cause into the parent message, duplicating `source()` chain rendering.
- **Why it matters:** Operator-facing log lines say `io: ...` without naming the operation; triage requires reading the source to recover which `?` site fired. With many variants flattening the same way, the error stream loses operator-facing structure. The 17 pre-existing reviews flagged the same shape; reuse the same fix here.
- **Suggested fix:** Drop the `#[from]` impls on `Io` and `PspRead`; rename variants by operation; remove the `: {0}` interpolation so `source()` carries the cause.
  ```rust
  #[error("failed to load chunk on chromosome {chrom_id}")]
  LoadChunk { chrom_id: u32, #[source] source: ChunkLoadError<PspReadError> },
  #[error("failed to open PSP file {path}")]
  OpenPsp { path: PathBuf, #[source] source: PspReadError },
  #[error("failed to write cohort VCF at {chrom_id}:{start}..{end}")]
  WriteVcf { chrom_id: u32, start: u32, end: u32, #[source] source: VcfWriteError },
  // …
  ```

- `src/var_calling/cohort_block/driver.rs:228-230` — **[Major]** M5: `let _ = std::fs::remove_file(tmp_path_for(output));` swallows the cleanup error with no log
- **Categories:** errors, smells, defaults (cross-cat) — convergent
- **Confidence:** High
- **Problem:** On the error path the driver does `let _ = std::fs::remove_file(...)` with the comment "Best-effort cleanup of the writer's tmp file" and no `tracing` event. The rule (errors.md §24 + smells.md §14) requires either propagation or a structured log event naming the failure mode being discarded. When the driver fails *and* leaves a tmp file behind (permissions, stale path, disk full), the operator has no breadcrumb.
- **Why it matters:** Best-effort cleanup is fine; *silent* best-effort cleanup is not. Pair with B1's `writer.abort()` fix.
- **Suggested fix:**
  ```rust
  if let Err(remove_err) = std::fs::remove_file(tmp_path_for(output)) {
      tracing::warn!(
          error = %remove_err,
          path = %tmp_path_for(output).display(),
          "failed to remove tmp VCF during error cleanup",
      );
  }
  ```

- `src/var_calling/cohort_block/columns.rs:438-448` — **[Major]** M6: `u32_from_usize` wraps silently in release builds; the cast is in `push_record` / `push_row_from`'s hot path
- **Categories:** reliability, errors (Minor convergence)
- **Confidence:** High
- **Problem:** `value as u32` in release truncates; only a `debug_assert!` catches. Doc comment names the wrap as deliberate ("release builds wrap silently — same contract as PSP's CSR writers"). CSR offset overflow corrupts every downstream index → unreadable columns → silently wrong likelihoods on the affected chunk.
- **Why it matters:** PSP files come from disk and are not provably bounded by the type system. A pathological / corrupted PSP could overflow this and the worker would compute wrong outputs without any operator-visible signal. Filed Major (not Blocker) because High confidence the PSP cap makes this unreachable on real tomato data today; promote to Blocker if PSP becomes a less-trusted boundary.
- **Suggested fix:** Surface as a typed error.
  ```rust
  pub(crate) fn u32_from_usize(value: usize) -> Result<u32, ChunkLoadError<I>> {
      u32::try_from(value).map_err(|_| ChunkLoadError::CsrOffsetOverflow { value })
  }
  ```
  Or, if the cost of plumbing `Result` through is too high, use `try_into().expect(...)` to panic loudly with a named invariant rather than wrap silently.

- `src/var_calling/cohort_block/driver.rs:564` — **[Major]** M7: `params.target_window_count.max(1)` silently rewrites an out-of-range CLI input
- **Categories:** defaults
- **Confidence:** High
- **Problem:** The pre-pass rejects `target_window_count = 0` with `FixBoundariesError::ZeroTargetWindowCount`; the driver bypasses that validation by clamping with `.max(1)` at the call site. An operator passing `--worker-windows-per-chunk 0` gets silent single-window behaviour with no log explaining the rewrite. `ChunkDriverParams::target_window_count` has no documented `Default`.
- **Why it matters:** Either `0` is valid (then it must be (a) named `DEFAULT_WORKER_WINDOWS_PER_CHUNK`, (b) doc'd at the call site, (c) logged when the fallback fires), or `0` is invalid (then propagate the pre-pass error).
- **Suggested fix:** Pick one of:
  ```rust
  // (a) document the fallback + log
  let target = if params.target_window_count == 0 {
      tracing::debug!("target_window_count = 0 → falling back to {}", DEFAULT_WORKER_WINDOWS_PER_CHUNK);
      DEFAULT_WORKER_WINDOWS_PER_CHUNK
  } else { params.target_window_count };
  fix_boundaries(chunk, carryover, fix_scratch, max_group_span, target)?;
  // (b) propagate
  fix_boundaries(chunk, carryover, fix_scratch, max_group_span, params.target_window_count)?;
  ```

- `src/var_calling/cohort_block/loader.rs:220-274` — **[Major]** M8: `effective_initial_span` no-op + `max_span < initial_span` silently treated as `max_span = initial_span`
- **Categories:** reliability
- **Confidence:** High
- **Problem:** `let effective_initial_span = initial_span.min(max_span.max(initial_span));` is structurally `initial_span` (because `max(max_span, initial_span) >= initial_span` ⇒ `min(initial_span, …) == initial_span`). The next line similarly upgrades `max_span < initial_span` to `initial_span`. A caller passing `max_span < initial_span` does not get the cap they asked for; the function silently widens it. The driver currently passes `max_load_span >= initial_load_span` (so the bug doesn't fire today), but the API contract is wrong and would bite the next caller.
- **Why it matters:** A future caller reading the API doc ("hard cap on the chunk's BP span across all extension iterations. The loader never pulls records past range_start + max_span.") sees behavior that violates the contract.
- **Suggested fix:** Either reject at entry with a new `ChunkLoadError::MaxSpanBelowInitial` variant, or delete the no-op chain and document the invariant explicitly. Cross-category note: `as u32` saturations at `loader.rs:261` and `pre_pass.rs:163` should use `u32::try_from` for the same reason.

- `src/var_calling/cohort_block/driver.rs:103` — **[Major]** M9: `target_variants_per_chunk == 0` is a sentinel-as-toggle, not a named default
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Both `ChunkDriverParams::target_variants_per_chunk` (driver.rs:96–102) and `load_chunk_from_iters` (loader.rs:150–153) document `0` as "disables the variant-bounded extension". The off-switch and the threshold share one `u32` field. No named constant, no startup `tracing` event recording "extension disabled because target = 0". Reading `drive_cohort_chunked(..., params)` at a call site does not reveal that `0` is special — only the field doc does.
- **Why it matters:** Per the defaults rule ("every behavioral default lives in a named `pub const`"), the disabling sentinel deserves a named constant and the call site deserves an explicit "off" representation. Gates Open Question 3.
- **Suggested fix:** Either name the sentinel (`pub const TARGET_VARIANTS_DISABLED: u32 = 0;` doc-commented next to `DEFAULT_CHUNK_GENOMIC_SPAN`) and `tracing::debug!("extension disabled")` when applied, or lift the toggle into the type:
  ```rust
  pub target_variants_per_chunk: Option<std::num::NonZeroU32>,
  ```

- `src/var_calling/cohort_block/worker.rs:17,236` — **[Major]** M10: Two unresolved-link errors in scope (`build_overlapping_variant_group`, `PosteriorEngine`) fail the crate's `broken_intra_doc_links = "deny"` policy
- **Categories:** tooling, extras (Minor convergence), reliability (Nit convergence)
- **Confidence:** High
- **Problem:** [worker.rs:17](../../../src/var_calling/cohort_block/worker.rs#L17) links to `build_overlapping_variant_group`, which is `#[cfg(test)]` only — rustdoc's default cfg can't reach it. [worker.rs:236](../../../src/var_calling/cohort_block/worker.rs#L236) links to `PosteriorEngine`, which is not imported into `worker.rs`'s `use` list (only `PosteriorEngineConfig` is). `Cargo.toml` sets `broken_intra_doc_links = "deny"`, so these are hard `cargo doc --no-deps` errors.
- **Why it matters:** `cargo doc` is one of the CI gates; the crate's own deny-policy blocks merge on these.
- **Suggested fix:** (a) For [worker.rs:17](../../../src/var_calling/cohort_block/worker.rs#L17): drop the link target, leave the prose `` `build_overlapping_variant_group` `` as a code span. (b) For [worker.rs:236](../../../src/var_calling/cohort_block/worker.rs#L236): spell the path so rustdoc can resolve without an import.
  ```rust
  ///   tuning passed straight to
  ///   [`PosteriorEngine`](crate::var_calling::posterior_engine::PosteriorEngine).
  ```

- `benches/cohort_e2e_perf.rs:286`, `examples/profile_cohort_e2e.rs:152`, `examples/dhat_var_calling.rs:121` — **[Major]** M11: `VarCallingArgs` gained two new required fields without updating call sites; `cargo clippy --all-targets -- -D warnings` exits 101
- **Categories:** tooling, refactor_safety (positive observation — the safety mechanism worked)
- **Confidence:** High
- **Problem:** The driver added required `target_variants_per_chunk: u32` and `worker_windows_per_chunk: usize` to `VarCallingArgs` ([src/pop_var_caller/var_calling.rs:105,116](../../../src/pop_var_caller/var_calling.rs#L105)). Three workspace call sites still use the pre-rewrite struct literal and fail with `E0063: missing field`. The whole `--all-targets` CI step exits 101 before clippy runs. Refactor-safety side: this is exactly the safety mechanism working — the compiler flagged every stale call site — but the breakage was left unfixed.
- **Why it matters:** Examples and the bench are how perf for this module is reproduced. Leaving them broken on the feature branch means the bench regression detector goes dark exactly when the rewrite lands.
- **Suggested fix:** Add the two new fields to each of the three call sites with the same defaults the CLI uses (`target_variants_per_chunk: 0`, `worker_windows_per_chunk: 1`) so the bench/example numbers stay comparable to the pre-rewrite shape. A small follow-up `VarCallingArgs::for_profiling(...)` constructor would prevent the same drift on the next field addition.

- `src/var_calling/cohort_block/driver.rs:483` — **[Major]** M12: `load_and_run_chunk_with_retry` is 19 positional parameters across ~165 lines, mixing five distinct phases
- **Categories:** smells
- **Confidence:** High
- **Problem:** Carryover snapshot, load+retry loop, worker-pool sizing, sequential partition+prefetch, parallel rayon math, sequential drain — all share one parameter list under `#[allow(clippy::too_many_arguments)]` with no justification comment. Every new per-stage knob lands as another positional `u32`; a swap at the call site is silent.
- **Why it matters:** Highest-friction surface in the module. Refactoring opportunity that the bench/example E0063 errors (M11) already validate: the compiler is the wrong place to catch positional `u32` swaps.
- **Suggested fix:** Promote per-iteration mutable state into a `ChunkLoopState<'_>` bundle constructed in `drive_one_chrom_generic`. Split into three private helpers all taking `&mut ChunkLoopState`: `load_chunk_with_safe_boundary_retry`, `run_chunk_windows_parallel`, `drain_window_outputs`. Drop the `#[allow]`.

- `src/var_calling/cohort_block/loader.rs:186` — **[Major]** M13: `load_chunk_from_iters` 9 positional parameters with no extent-grouping struct
- **Categories:** smells
- **Confidence:** High
- **Problem:** `scratch`, `out`, `chrom_id`, `range_start`, `initial_span`, `target_variants`, `max_span`, `per_sample_iters`, `carryover` — four of these (`range_start`, `initial_span`, `max_span`, `target_variants`) always travel together as the load-extent policy. Tests build them positionally; the `#[allow(clippy::too_many_arguments)]` has no justification.
- **Why it matters:** Same shape problem as M12. Swapping `initial_span` and `target_variants` at a call site is the `u32`/`u32` mix-up the compiler will never catch.
- **Suggested fix:**
  ```rust
  pub struct ChunkLoadExtent {
      pub chrom_id: u32,
      pub range_start: u32,
      pub initial_span: u32,
      pub target_variants: u32,
      pub max_span: u32,
  }
  pub fn load_chunk_from_iters<I, E>(
      scratch: &mut ChunkLoadScratch, out: &mut MaterialisedChunk,
      extent: ChunkLoadExtent,
      per_sample_iters: Vec<I>, carryover: &mut [SampleColumns],
  ) -> Result<ChunkLoadStats, ChunkLoadError<E>>
  ```

- `src/var_calling/cohort_block/driver.rs:511-519, 570-575` — **[Major]** M14: Carryover snapshot and restore are two near-identical loops, missing a `SampleColumns::clone_from_columns` helper
- **Categories:** smells
- **Confidence:** High
- **Problem:** Two read-write loops walking the same parallel-array structure invite skew. A new `SampleColumns` field would have to update both by hand; either side regressing would silently flatten the carryover on retry and surface as wrong VCFs.
- **Suggested fix:** Add `SampleColumns::clone_from_columns(&mut self, other: &SampleColumns)` that does `self.clear(); for row_idx in 0..other.n_records() { self.push_row_from(other, row_idx); }`. Both driver sites collapse to one line each. Add a unit test.

- `src/var_calling/cohort_block/worker.rs:51-96` — **[Major]** M15: `ColumnarPipelineScratch` mixes per-group scratch buffers with running per-window counters
- **Categories:** smells
- **Confidence:** High
- **Problem:** Eight fields are per-group scratch (`unify`, `unified`, `project`, `projection`, `log_likelihoods`, `chain_anchor_flags`, `record_scratch`, `math`); two (`lh_cap_groups_skipped`, `lh_cap_alleles_in_skipped`) are running counters needing `take_lh_cap_stats` to extract. Different methods touch disjoint subsets of fields. Counter lifetime (per-`run_window`) and scratch lifetime (per-group) are bundled into one struct.
- **Suggested fix:** Move counters into a sibling field of `WorkerSlot` (`pub stats: WindowRunStats`) or return them from `run_window`. `ColumnarPipelineScratch` becomes pure scratch, and the slot's interaction with the driver becomes explicit.

- `src/var_calling/cohort_block/columns.rs:48-62` — **[Major]** M16: `SampleColumns` exposes 13 ungrouped `pub` column fields with all `push_record` / `push_row_from` / `clear` / `truncate` paths touching every one by name
- **Categories:** smells
- **Confidence:** Medium (this is a deliberate columnar layout choice — flagging the *surface*, not the layout)
- **Problem:** Every field is `pub`, every field is a `Vec`. New per-allele scalar fields (the next column to land in `AlleleObservation`) propagate to 7+ hand-edited sites: the struct, four maintenance loops (`push_record`, `push_row_from`, `clear`, `truncate`), the test fixtures, the carryover snapshot/restore (M14).
- **Why it matters:** The 13 fields are the columnar promise — that's the design. The *all-`pub` open-access surface with no grouping* is what invites the shotgun-edit smell.
- **Suggested fix:** Group along CSR boundaries: `PerAlleleFixed` (7 fixed-width per-allele scalars), `PerAlleleSeq` (offsets + bytes), `PerAlleleChainIds` (offsets + ids). `push_record` writes one line per group instead of per column. The visibility stays `pub`; downstream readers' column-iteration access patterns are preserved.

- `src/var_calling/cohort_block/columns.rs:116` — **[Major]** M17: Partial destructure of `AlleleSupportStats` with trailing `..` silently absorbs new fields
- **Categories:** refactor_safety
- **Confidence:** High
- **Problem:** `SampleColumns::push_record` destructures `AlleleSupportStats { num_obs, q_sum, fwd, placed_left, placed_start, mapq_sum, mapq_sum_sq, .. }`. `AlleleSupportStats` is `#[non_exhaustive]` at `src/pileup_record.rs:43` — that attribute is a *cross-crate* signal; inside the same crate the catch-all `..` is purely a refactor-safety footgun. The 7 named fields already cover the whole struct, so the `..` adds zero information today. A future `AlleleSupportStats` field would be silently dropped on the floor for every record flowing through `push_record` / `push_row_from`.
- **Why it matters:** Drops would not necessarily change current byte output (the existing columns don't carry the field) but would silently break byte-identity for any downstream code that adds the field through `materialise_record`. Exactly the refactor-safety failure mode the rule targets.
- **Suggested fix:** Drop the trailing `..`. The 7 named fields already cover the struct.
  ```rust
  let AlleleSupportStats {
      num_obs, q_sum, fwd, placed_left, placed_start, mapq_sum, mapq_sum_sq,
  } = allele.support;
  ```

- `src/var_calling/cohort_block/kernels/unify_alleles.rs:267` — **[Major]** M18: `#[allow(dead_code)] pub(crate) chain_id_scratch` is read only by one scratch-clearing unit test
- **Categories:** refactor_safety, tooling (Minor convergence), smells (Minor convergence)
- **Confidence:** High
- **Problem:** `UnifyAllelesScratch.chain_id_scratch: Vec<ChainId>` is `clear()`-ed in `clear` and asserted in `scratch_clear_drops_contents_preserves_capacity` (line 1758) but never read or written by the kernel — `detect_compound_candidates_columnar` (M3) uses local `BTreeMap`s instead. The annotation hides a real signal: the field is a no-op carrying allocator capacity.
- **Why it matters:** A reader cannot tell whether this is "waiting to be wired up" or "obsolete". The cleared field gives a false sense of scratch-reuse hygiene.
- **Suggested fix:** Either delete (field, clear call, test references) — the column-native compound-detection code walks chain ids directly — or wire it into the M3 fix as the scratch buffer the doc-comment promises.

- `src/var_calling/cohort_block/partition.rs:261-421` — **[Major]** M19: `partition_window` doesn't validate that `masked_intervals` is sorted and non-overlapping; the cursor logic silently misbehaves if the contract is violated
- **Categories:** reliability
- **Confidence:** High
- **Problem:** `partition_window` checks only the *current* mask interval (`masked_intervals[mask_cursor]`). If the slice is unsorted (e.g. `[10..20, 5..7, 30..40]`) or has overlaps, the cursor advances past the wrong interval and emit/skip decisions become silently wrong. The function doc says "sorted, non-overlapping"; no debug-assert, no error variant.
- **Why it matters:** Today the only caller (`compute_dust_mask_for_chrom`) returns sorted+non-overlapping intervals because sdust does. A future test helper or alternative mask source could pass malformed input and get silent wrong filtering. Other entry points in the module validate their contracts explicitly.
- **Suggested fix:**
  ```rust
  debug_assert!(
      masked_intervals.windows(2).all(|w| w[0].end <= w[1].start),
      "masked_intervals must be sorted and non-overlapping",
  );
  ```
  Plus a release-mode error variant `PartitionError::MaskNotSorted` for the public-API entry point.

- `src/var_calling/cohort_block/driver.rs:240-268` — **[Major]** M20: `emit_or_drop` filter order is asserted in code but not pinned to the streaming pipeline's order by any test
- **Categories:** reliability, extras (cross-cat)
- **Confidence:** High
- **Problem:** Order is `min_alt_obs` → `is_variant_call` (hom_ref) → `qual_phred` → `mapq_diff_t` → `unconverged` counter → write. The comment claims "same order" as streaming, but the streaming driver applies `min_alt_obs` *pre*-EM (between merger and EM at `cohort_driver.rs:286-310`). If a record would trip multiple filters, the per-category counter assignment depends on order — a swap would silently move records between categories and the 3-tomato fixture's reported numbers would diverge. No unit-level test pins each filter against a synthetic single-record fixture covering each rejection case. Ties to Open Question 4.
- **Why it matters:** Per-category counts are reported to operators; an order swap is silently wrong. B5 fixes the cross-driver equality; M20 fixes the per-filter unit equivalence.
- **Suggested fix:** Add five unit tests in a new `mod tests` block inside `driver.rs`, each fixture tripping exactly one of the four filters; assert the matching counter increments and others stay at zero. Add a fifth fixture that *would* trip multiple filters and document which one wins. See §8.

- `src/var_calling/cohort_block/driver.rs:240-268` — **[Major]** M21: `emit_or_drop` does not surface VCF write errors with chrom/range context
- **Categories:** reliability, errors (cross-cat)
- **Confidence:** High
- **Problem:** `writer.write_record(&record)?;` propagates `VcfWriteError` directly via `From`; the parent `ChunkDriverError::VcfWrite` carries no chrom/range fields. Real-world I/O failures (full disk mid-write) need the failing record's coordinates for triage.
- **Suggested fix:** Reshape `ChunkDriverError::VcfWrite` to `ChunkDriverError::WriteVcf { chrom_id: u32, start: u32, end: u32, #[source] source: VcfWriteError }` (pairs with M4).

- `src/var_calling/cohort_block/kernels/unify_alleles.rs:425-491` — **[Major]** M22: `enforce_max_alleles_columnar` tie-break (stable sort by `cohort_count`, original index breaks ties) is not pinned by any test against the row-shape kernel's tie-break
- **Categories:** reliability
- **Confidence:** Medium
- **Problem:** Existing `unify_max_alleles_cap_*` tests use distinct cohort_counts so ties are never exercised. The 3-tomato byte-identity fixture may not hit ties either. A real-data input with `cohort_count` ties at the cap boundary could pick a different ALT vs `main`, silently violating byte-identity.
- **Suggested fix:** Add a 4-ALT fixture where two have equal counts and `max_alleles` cuts between them; assert column-native and row-shape kernels pick the same ALT. Document the tie-break rule in the function's doc string. If the two differ, add a secondary tie-break (lex-smallest seq) and lock it.

- `src/var_calling/cohort_block/worker.rs:301-473` — **[Major]** M23: `build_posterior_record_columnar` returns `Ok(None)` for `n_alleles < 2` silently — no counter, while the sibling `n_alleles > max_alleles_lh_calc` skip *does* increment `lh_cap_groups_skipped`
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Two adjacent silent-skip paths; one has an operator-visible counter (`lh_cap_*`), the other doesn't. Operators auditing why no record was emitted at a known variant locus have no breadcrumb for the REF-only-post-unify case.
- **Suggested fix:** Add `groups_skipped_post_unify_ref_only: u64` to `ChunkDriverStats` and a corresponding `take_*` API mirroring `take_lh_cap_stats`.

- `src/var_calling/cohort_block/worker.rs:284-299` — **[Major]** M24: `prefetch_window_ref_bytes` may fetch past `chrom_length` on the last chunk
- **Categories:** reliability
- **Confidence:** Medium
- **Problem:** Group ends are bounded by `chunk.safe_end ≤ chunk.range.end`. The last chunk's `range.end == chrom_one_past_end + last_chunk_logical_extension`. The pre-pass doesn't know `chrom_length` and doesn't clamp `safe_end` to `chrom_one_past_end`; a group whose `group_end > chrom_length` then asks the fetcher for bytes past the contig end. Behaviour depends on the fetcher (out-of-scope to verify here).
- **Suggested fix:** Clamp `safe_end` to `chrom_one_past_end` in the pre-pass for the last chunk, or have `prefetch_window_ref_bytes` clamp `span` to `chrom_length - group_start + 1`. Add a regression test on a chunk whose last group sits at chrom-end.

- `src/var_calling/cohort_block/loader.rs:469,476` and `src/var_calling/cohort_block/driver.rs:208` and `src/var_calling/cohort_block/pre_pass.rs:204,205` — **[Major]** M25: `.expect(...)` / `.unwrap()` on five sites without `// PANIC-FREE:` comments
- **Categories:** errors
- **Confidence:** Medium (the invariants hold today; rule 11 is about future-proofing)
- **Problem:** Five `expect`/`unwrap` sites — two in `pull_records_with_pos_under` (peek/next discipline), one in `drive_cohort_chunked` (`u32::try_from(chrom_idx)`), two in `pick_safe_end` (`scratch.timeline.last().unwrap()` / `prefix_max_reach.last().unwrap()`). Each panic string names the invariant but rule 11 wants a `// PANIC-FREE:` comment plus the call-site condition that maintains it.
- **Suggested fix:** Add `// PANIC-FREE:` comments naming the invariant + the maintaining site, or restructure with `let-else` to make the discipline visible.
  ```rust
  let Some(&(last_pos, _)) = scratch.timeline.last() else {
      return Ok(chunk.range.end);
  };
  ```

- `src/var_calling/cohort_block/pre_pass.rs:137` and `src/var_calling/cohort_block/partition.rs:365` — **[Major]** M26: Plain arithmetic `pos + ref_span - 1` on PSP-derived input would underflow on a zero-length REF
- **Categories:** errors
- **Confidence:** Medium
- **Problem:** `PileupRecord`'s doc says `alleles[0]` is REF with `seq.len() >= 1`, so `ref_span` is ≥ 1 by invariant. But that's a documented invariant of `PileupRecord`, not a type-system guarantee. A malformed PSP (or future PSP-format extension allowing zero-length REF) underflows the expression (debug panic, release wrap to `u32::MAX`). PSP is not fully trusted.
- **Suggested fix:** Saturate explicitly.
  ```rust
  let reach = pos.saturating_add(sample.ref_span_at(row_idx).max(1)).saturating_sub(1);
  ```

- `src/var_calling/cohort_block/driver.rs:82-85` ⇿ `worker.rs:247-248` — **[Major]** M27: `*_cfg` vs `*_config` collides for the same payload across the driver/worker boundary
- **Categories:** naming
- **Confidence:** High
- **Problem:** `ChunkDriverParams` names its fields `dust_cfg`, `grouper_cfg`, `per_group_cfg`, `posterior_cfg`. The worker's `run_window` / `build_posterior_record_columnar` name the same payloads `per_group_config` / `posterior_config`. The driver even rebinds them (`let posterior_cfg = &params.posterior_cfg;` / `let per_group_cfg = params.per_group_cfg;`) before passing into the worker's `*_config` parameters.
- **Why it matters:** Every reader has to mentally rename; greps for one form miss the other.
- **Suggested fix:** Pick one form crate-wide. `cfg` is the form `ChunkDriverParams` exposes; rename the worker-side parameters to match. Same change in `build_posterior_record_columnar` and call sites.

- `src/var_calling/cohort_block/worker.rs:648` — **[Major]** M28: `shared_ref_fetcher` is a noun-named function (it's a constructor)
- **Categories:** naming
- **Confidence:** High
- **Problem:** `pub fn shared_ref_fetcher<F>(fetcher: F) -> SharedRefFetcher` wraps an owned `ChromRefFetcher` into an `Arc<dyn …>`. Bare-noun function names are reserved for field-like accessors in Rust convention; converting constructors are `from_X` / `into_X` / `to_X`.
- **Suggested fix:** Rename:
  ```rust
  pub fn into_shared_ref_fetcher<F>(fetcher: F) -> SharedRefFetcher
  ```
  Update `mod.rs:41`.

- `src/var_calling/cohort_block/loader.rs:186-305` — **[Major]** M29: No test exercises `ChunkLoadError::SampleCountMismatch` or `ChunkLoadError::CarryoverLengthMismatch`
- **Categories:** reliability
- **Confidence:** High
- **Problem:** Both variants are declared and returned at the top of `load_chunk_from_iters`; the existing tests only assert `InvalidRange`, `UnexpectedChromosome`, and `UpstreamRead`. Future refactor dropping these checks (e.g. moving to debug-assert) would silently pass; the next mis-wiring would trip an `index out of range` panic deep in the loop instead of a clean error.
- **Suggested fix:** Add `loader_rejects_sample_count_mismatch` and `loader_rejects_carryover_length_mismatch` with the matching variant assertions (see §8).

- `src/var_calling/cohort_block/pre_pass.rs:113-185` — **[Major]** M30: No test covers `fix_boundaries`'s `CarryoverLengthMismatch` error path
- **Categories:** reliability
- **Confidence:** High
- **Suggested fix:** Add `pre_pass_rejects_mismatched_carryover_length` (see §8).

- `src/var_calling/cohort_block/driver.rs:170-235, 623-635` — **[Major]** M31: No concurrency test for the rayon `par_iter_mut().try_for_each(|slot| run_window(...))` path
- **Categories:** reliability
- **Confidence:** High
- **Problem:** The end-to-end byte-identity test exercising the parallel dispatch is the external 3-tomato fixture; all unit tests in `worker.rs` are single-window. A regression breaking per-slot disjointness (e.g. someone adds a `static mut` cache, or `Rc` inside a slot) would not be caught by the existing tests. The `unsafe_concurrency` agent confirmed soundness today is statically enforced by the `Send + !Sync` typedef; what's missing is a regression test that pins the parallel-vs-sequential equivalence.
- **Suggested fix:** Add `drive_one_chrom_par_iter_matches_sequential_run_window` — same chunk processed with `target_window_count = 4` (parallel) vs `target_window_count = 1` (sequential); assert emitted records byte-equal. Use `ThreadPoolBuilder` or `RAYON_NUM_THREADS` to force the parallel dispatch deterministically. See §8.

- `src/var_calling/cohort_block/kernels/project_scalars.rs:531-548` — **[Major]** M32: `locate_sample_row_idx` relies on `partition_window`'s `samples_at_pos` ordering invariant; neither side has a test linking them
- **Categories:** reliability
- **Confidence:** Medium
- **Problem:** `locate_sample_row_idx` does a linear scan and returns the first match. The function depends on `samples_at_pos[lo..hi]` being sorted-deduped by `sample_idx` (per [partition.rs:295](../../../src/var_calling/cohort_block/partition.rs#L295) comment). If a future partition change emitted a sample twice at the same position (dedup regression), the first hit would be returned silently wrong.
- **Suggested fix:** Add a `debug_assert!` in `locate_sample_row_idx` that the returned sample_idx is unique in `sample_range`. Add a proptest in `partition.rs` over random fixtures asserting `samples_at_pos[lo..hi]` is strictly ascending.

### Minor

- `src/var_calling/cohort_block/columns.rs:48`, `:378`; `partition.rs:89`; `loader.rs:81`; `driver.rs:120`; `kernels/unify_alleles.rs:86`; `kernels/project_scalars.rs:44`; `kernels/compute_log_likelihoods.rs:52` — **[Minor]** Mi1: Newly-public data structs with all-`pub` fields lack `#[non_exhaustive]`
- **Categories:** refactor_safety
- **Confidence:** Medium — gated on Open Question 1.
- **Problem:** Error enums in the same scope *are* `#[non_exhaustive]` and `ChunkDriverParams` is `#[non_exhaustive]`; the matching armour on data structs is missing. Adding fields silently breaks struct-literal construction at every external call site (the M11 bench/example breakage is exactly the failure mode `#[non_exhaustive]` on the *data* side would prevent for in-crate consumers).
- **Suggested fix:** Add `#[non_exhaustive]` to each. Replace the inline test struct literal at columns.rs:516 with the factory + field-assignment pattern.

- `src/var_calling/cohort_block/columns.rs:426` — **[Minor]** Mi2: `MaterialisedChunk::clear_data` deviates from the crate-wide `clear()` convention used by every other persistent buffer
- **Categories:** naming
- **Confidence:** High
- **Suggested fix:** Rename to `clear`; update the single in-tree caller at `loader.rs:224`.

- `src/var_calling/cohort_block/worker.rs:164,167,170` — **[Minor]** Mi3: `WorkerSlot.output_buf` is a bare adjective + generic `buf`; `WorkerSlot.scratch` is a bare generic noun where the type carries the domain
- **Categories:** naming
- **Confidence:** High
- **Suggested fix:** Rename `output_buf` → `posterior_records`, `scratch` → `pipeline_scratch` (mirrors the sibling `partition_scratch`).

- `src/var_calling/cohort_block/worker.rs:499-547` — **[Minor]** Mi4: Asymmetric abbreviation across sibling helpers: `unify_error_to_merger` / `project_scalars_error_to_merger` / `compute_ll_error_to_merger` (`ll` is the abbreviation; the kernel module spells it `log_likelihoods`); the same file uses `lh_cap_*` for the same concept
- **Categories:** naming
- **Confidence:** High
- **Suggested fix:** Rename to `unify_alleles_error_to_merger` / `project_scalars_error_to_merger` / `compute_log_likelihoods_error_to_merger`. Separately consider renaming `lh_cap_*` to `lh_calc_cap_*` to match the `max_alleles_lh_calc` config name.

- `src/var_calling/cohort_block/driver.rs:330` — **[Minor]** Mi5: `pool_allele_mapq` returns an anonymous `(u64, u64, u128)` primitive triple
- **Categories:** naming
- **Confidence:** Medium
- **Suggested fix:** Introduce a named return:
  ```rust
  struct PooledMapqMoments { n: u64, sum: u64, sum_of_squares: u128 }
  ```

- `src/var_calling/cohort_block/worker.rs:284, 682` ⇿ `driver.rs:349` — **[Minor]** Mi6: `ref_fetcher` vs `fetcher` for the same concept
- **Categories:** naming
- **Confidence:** High
- **Suggested fix:** Pick `ref_fetcher` everywhere (more domain-explicit; this codebase has PSP fetchers and FASTA fetchers).

- `src/var_calling/cohort_block/pre_pass.rs:113` (and the module name `pre_pass`) — **[Minor]** Mi7: `fix_boundaries` is a vague verb; `pre_pass` is a layer-pattern module name, not a domain-concept name
- **Categories:** naming
- **Confidence:** Medium
- **Suggested fix:** Rename `pre_pass.rs` → `chunk_boundaries.rs`; rename `fix_boundaries` → `finalise_chunk_boundaries` (or `commit_chunk_boundaries`). Propagate to `FixBoundariesError` / `FixBoundariesScratch`.

- `src/var_calling/cohort_block/mod.rs:23-31` — **[Minor]** Mi8: Every submodule declared `pub mod` while `mod.rs` only re-exports a curated surface
- **Categories:** module_structure
- **Confidence:** High — gated on Open Question 1.
- **Problem:** All seven submodules under `cohort_block/` are `pub mod`, so the entire submodule namespace becomes part of the crate's public API and surfaces in rustdoc. The `pub use` block then re-exports a narrower set, but items the re-export omits are still reachable via the deep path. The kernels' ~40 `pub(crate)` / `pub` items are all in the public docs. Call sites in `pop_var_caller/contamination_chunked_stream.rs:43-44` and `pop_var_caller/estimate_contamination.rs:42` actually use deep paths (`cohort_block::columns::{...}`, `cohort_block::loader::{...}`) even though `mod.rs` re-exports those names — only `pop_var_caller/var_calling.rs:42-45` uses the curated surface.
- **Suggested fix:** `pub mod` → `pub(crate) mod` for every submodule (one-line change per module). Update the three external call sites to import via the re-exports for a single canonical path per item.

- `src/var_calling/cohort_block/driver.rs:408` — **[Minor]** Mi9: `Arc::new(StreamingChromRefFetcher)` for single-owner sequential use, with an `#[allow(clippy::arc_with_non_send_sync)]` that documents the smell
- **Categories:** idiomatic, smells (Nit convergence)
- **Confidence:** High
- **Problem:** The fetcher is `!Sync` and is never handed to rayon (the parallel block only reads `slot.pre_fetched_ref_bytes`). There is exactly one owner (`drive_one_chrom_generic`), exactly one consumer (`prefetch_window_ref_bytes` in the sequential phase), and no shared-ownership requirement. The `Arc` adds a heap allocation per chromosome, a refcount bump per `clone()`, and a lint suppression.
- **Suggested fix:** Thread `&StreamingChromRefFetcher` (or `&dyn ChromRefFetcher`) through `load_and_run_chunk_with_retry`; delete the `#[allow]`. The `unsafe_concurrency` audit confirmed nothing in the parallel section needs the Arc.

- `src/var_calling/cohort_block/driver.rs:586` — **[Minor]** Mi10: `let windows = chunk.windows.clone();` is an unnecessary defensive clone
- **Categories:** idiomatic
- **Confidence:** Medium
- **Problem:** Inside the loop body, `chunk` is borrowed only immutably (`partition_window(chunk, …)`, `chunk.n_samples()`). The mutable borrow lives on `worker_pool.slots`, disjoint from `chunk`. Iterating `chunk.windows.iter()` directly is a shared reborrow that compiles under modern Rust's NLL.
- **Suggested fix:** Drop the clone; iterate `chunk.windows.iter()` directly. Verify with `cargo check`.

- `src/var_calling/cohort_block/worker.rs:243` — **[Minor]** Mi11: `run_window` takes `posterior_config: PosteriorEngineConfig` by value, forcing the driver to clone per worker slot
- **Categories:** idiomatic
- **Confidence:** High
- **Problem:** `PosteriorEngineConfig` is not `Copy` (it owns `Option<Vec<f64>>` for per-sample fixation overrides). Nothing inside `run_window` needs ownership — it passes through as `&posterior_config` to `build_posterior_record_columnar`, which already takes by-ref.
- **Suggested fix:** Change to `posterior_config: &PosteriorEngineConfig`. Drop the `posterior_cfg.clone()` at the rayon site.

- `src/var_calling/cohort_block/kernels/unify_alleles.rs:614, 739` — **[Minor]** Mi12: Double-clone of `projection_buf` for the dedup-map insertion
- **Categories:** idiomatic
- **Confidence:** High
- **Problem:** The pattern `let seq_copy = scratch.projection_buf.clone(); scratch.byte_index.insert(seq_copy.clone(), idx); push_allele_into_scratch(scratch, &seq_copy, …);` clones twice where one move suffices.
- **Suggested fix:**
  ```rust
  let key = scratch.projection_buf.clone();
  push_allele_into_scratch(scratch, &key, false, false, n_samples);
  scratch.byte_index.insert(key, idx); // move, no second clone
  ```

- `src/var_calling/cohort_block/worker.rs:648, 678; kernels/unify_alleles.rs:534` — **[Minor]** Mi13: Three `pub` items have no caller outside their module
- **Categories:** idiomatic
- **Confidence:** High — gated on Open Question 1.
- **Problem:** `worker::shared_ref_fetcher` is `pub use`'d in `mod.rs:41` but has no consumer in `src/`, `tests/`, `benches/`, or `examples/`. `worker::unified_alleles_for_group_columnar` + `worker::WorkerUnifyError` are only called from `worker.rs:914` inside `#[cfg(test)]`. `kernels::unify_alleles::project_per_position_alleles_columnar` is called only from its own test sites.
- **Suggested fix:** Demote to `pub(crate)` (or `pub(super)`). Drop the `pub use shared_ref_fetcher` from `mod.rs` if no caller materialises.

- `src/var_calling/cohort_block/pre_pass.rs:172` and `src/var_calling/cohort_block/partition.rs:296` — **[Minor]** Mi14: `match Ok(idx) => idx, Err(idx) => idx` on `Result<usize, usize>` should be `Ok(idx) | Err(idx) => idx`
- **Categories:** idiomatic
- **Confidence:** High

- `src/var_calling/cohort_block/kernels/unify_alleles.rs:1`, `worker.rs:1`, `loader.rs:1` — **[Minor]** Mi15: Three files past the soft-cap line count (`unify_alleles.rs` 1766 LoC, `worker.rs` 1166 LoC, `loader.rs` 954 LoC); sub-step boundaries are already commented but stay in one file
- **Categories:** smells, module_structure (Nit convergence)
- **Confidence:** High
- **Suggested fix:** Promote `kernels/unify_alleles.rs` to `kernels/unify_alleles/{mod.rs, per_position.rs, compound.rs, cap.rs, serialize.rs, tests.rs}`. Split `worker.rs` into `worker/{mod.rs, pool.rs, run.rs, error_shims.rs, test_adapter.rs}`. Move `loader.rs`'s test module into `loader/tests.rs`.

- `src/var_calling/cohort_block/worker.rs:609` — **[Minor]** Mi16: `#[cfg(test)] pub(crate) build_overlapping_variant_group` is reached from 5 test sites in 4 modules — vestigial cross-module reach back into `worker`
- **Categories:** smells
- **Confidence:** High
- **Problem:** Test-only oracle's home is wrong. Moving to `test_helpers.rs` (already `#[cfg(test)] pub(crate)`) drops the cross-module reach and resolves M10 (worker.rs:17 broken doc link).
- **Suggested fix:** Move to `test_helpers.rs`; update the 5 import sites; drop the `#[cfg(test)] use` block at `worker.rs:45-49`.

- `src/var_calling/cohort_block/kernels/unify_alleles.rs:228, 289`; `worker.rs:164, 287` — **[Minor]** Mi17: Three `Vec<Vec<_>>` jagged arrays inside a module that preaches CSR
- **Categories:** smells
- **Confidence:** Medium
- **Problem:** `WorkingAllele.sources_per_sample: Vec<Vec<(u32,u32)>>`, `UnifyAllelesScratch.other_per_sample: Vec<Vec<(u32,u32)>>`, `WorkerSlot.pre_fetched_ref_bytes: Vec<Vec<u8>>` — each inner `Vec` is a separate allocation that never compacts. The CSR shape these are paired with is the deliberate alternative.
- **Suggested fix:** Convert to flat-`values` + `offsets` CSR. Pair with M2 (which fixes the prefetch case specifically).

- `src/var_calling/cohort_block/columns.rs:49` and throughout — **[Minor]** Mi18: `u32` for 1-based positions / chrom_ids / spans is primitive obsession
- **Categories:** smells
- **Confidence:** Medium
- **Suggested fix:** Introduce `OneBasedPos(u32)`, `OneBasedRange(Range<u32>)`, `ChromId(u32)` newtypes. Adopt incrementally at the highest-payoff sites: `MaterialisedChunk.range`, `chunk.windows`, `partition.group_starts/group_ends`, loader's `range_start`/`max_span`/`attempt_end`, and the DUST mask's 0-based → 1-based conversion.

- `src/var_calling/cohort_block/driver.rs:80-113` — **[Minor]** Mi19: `ChunkDriverParams` is 12 fields with mixed-axis tuning
- **Categories:** smells
- **Confidence:** Medium
- **Suggested fix:** Group along stage boundaries: `pub sizing: ChunkSizingParams` (chunk_genomic_span, target_variants_per_chunk, target_window_count), `pub downstream: DownstreamFilterParams` (min_qual_phred, min_alt_obs, no_mapq_diff_filter, min_mapq_diff_t).

- `src/var_calling/cohort_block/driver.rs:119` and startup — **[Minor]** Mi20: No startup `tracing::info!` listing effective `ChunkDriverParams` values; operators can't recover "what config did this run use?"
- **Categories:** defaults
- **Confidence:** Medium
- **Suggested fix:**
  ```rust
  tracing::info!(
      chunk_genomic_span = params.chunk_genomic_span,
      target_variants_per_chunk = params.target_variants_per_chunk,
      target_window_count = params.target_window_count,
      min_qual_phred = params.min_qual_phred,
      no_complexity_filter = params.no_complexity_filter,
      "chunk driver starting"
  );
  ```

- `src/var_calling/cohort_block/pre_pass.rs:245` — **[Minor]** Mi21: `emit_windows` may silently produce fewer windows than `target_window_count` requested
- **Categories:** defaults
- **Confidence:** High
- **Problem:** Documented on the function ("the partition may shrink when safe positions are sparse, but never errors") but invisible at the call site. An operator setting `--worker-windows-per-chunk 8` and observing no speed-up has no signal whether chunks were partitioned into 8 or silently down-graded to 1.
- **Suggested fix:** Add a `tracing::debug!` in `emit_windows` (or `fix_boundaries`) when `chunk.windows.len() < target_window_count`. Add a `chunks_with_fewer_windows_than_requested: u64` counter on `ChunkDriverStats`.

- `src/var_calling/cohort_block/columns.rs:316-321` — **[Minor]** Mi22: `drain_rows_from_into` has no test for `split_row_idx == 0` and `split_row_idx == n_records()` boundaries
- **Categories:** reliability

- `src/var_calling/cohort_block/pre_pass.rs:296-315` — **[Minor]** Mi23: `slide_left_to_safe` has no direct unit test for the `desired > safe_end` clamp path or the `min_open == desired` no-iteration path
- **Categories:** reliability

- `benches/cohort_e2e_perf.rs` — **[Minor]** Mi24: Hot-path bench exists but carries no `// REGRESSION THRESHOLD: N%` comment
- **Categories:** extras
- **Confidence:** High
- **Suggested fix:** Add `// REGRESSION THRESHOLD: 10%` above each `criterion_group!` / `bench_function!`, and wire a `cargo bench --bench cohort_e2e_perf --baseline main` step into the verify skill.

- `src/var_calling/cohort_block/test_helpers.rs:20-21` — **[Minor]** Mi25: `clippy::doc_lazy_continuation` ×2 in `allele()`'s rustdoc — `-D warnings` clippy job ships a deny error in a `pub(crate)` doc
- **Categories:** tooling
- **Suggested fix:** Re-flow the doc so the `+` isn't at the start of a line.

- `src/var_calling/cohort_block/driver.rs:105, 108`; `worker.rs:9` — **[Minor]** Mi26: Three `rustdoc::redundant_explicit_links` warnings; `cargo doc -D warnings` would gate on these
- **Categories:** tooling
- **Suggested fix:** Drop the explicit URL, let rustdoc resolve from context.

### Nits

A single mechanical pass clears the remaining lint surface; do not enumerate individually.

- **16 `cargo clippy --all-targets -- -D warnings` errors in scope**, all in tests / docs / test helpers:
  - `clippy::single_range_in_vec_init` ×8 at `columns.rs:520`, `partition.rs:711,727,742,771,789,793,818`.
  - `clippy::type_complexity` ×4 at `kernels/unify_alleles.rs:1361,1396` and `worker.rs:926,956` — return type is the `(Vec<u8>, bool, Vec<(usize, usize)>)` oracle tuple; either add `#[allow(clippy::type_complexity)]` with a one-line justification or introduce `type OracleAllele = (Vec<u8>, bool, Vec<(usize, usize)>);` in the test module.
  - `clippy::bool_assert_comparison` ×2 at `kernels/unify_alleles.rs:1608,1610` — rewrite `assert_eq!(out.cap_protected[0], true)` as `assert!(out.cap_protected[0])`.
  - `clippy::doc_lazy_continuation` ×2 at `test_helpers.rs:20,21` — see Mi25.
- **Twelve `#[allow(clippy::too_many_arguments)]` annotations** across the module (`driver.rs:171,379,483`, `loader.rs:186`, `kernels/compute_log_likelihoods.rs:123,287`, `kernels/project_scalars.rs:168`, `kernels/unify_alleles.rs:858,950`, `worker.rs:242,310,677`) lack the per-call-site justification comment the tooling rule mandates. Add a one-line "scratch-threaded columnar kernel; arg list mirrors the columnar inputs the kernel reads" each.
- `#[allow(clippy::arc_with_non_send_sync)]` at `driver.rs:408` and `#[allow(clippy::needless_range_loop)]` at `kernels/project_scalars.rs:263` likewise lack justification comments.
- `driver.rs:472` — `MAX_CHUNK_SPAN_GROWTH: u32 = 8` is declared between two functions; group with `DEFAULT_CHUNK_GENOMIC_SPAN` at file top.
- `kernels/unify_alleles.rs:262` — `byte_index: AHashMap<Vec<u8>, usize>` should name the mapping (`seq_to_allele_idx`), not the storage shape.
- `pre_pass.rs:147-155` — hand-rolled dedup-with-merge loop; the standard `Vec::dedup_by` covers it.
- `driver.rs:380` — `drive_one_chrom_generic` encodes the type parameter into the name; the function is private with no non-generic sibling, drop `_generic`.
- `kernels/unify_alleles.rs:454` — `to_remove_sorted` reads as the action, not the value; `dropped_indices` matches surrounding `dropped` / `prunable` / `protected` vocabulary.

## 7. Out of scope observations

- **Streaming driver (`drive_cohort_pipeline` in `src/pop_var_caller/cohort_driver.rs`)**: ties to Open Question 2. If staying long-term, B5's oracle test is straightforward. If being removed, B5 needs to be a checked-in golden VCF snapshot test instead — and the removal should not land before that snapshot is captured.
- **Pre-existing 17 clippy errors elsewhere in the crate** (rustup 1.94 → 1.95 bump): explicitly out of scope per the brief; tracked under the Phase A impl report's deferred-work list.
- **PSP reader trust boundary** (`u32_from_usize` silent wrap, M6 / Mi-class): if PSP is reclassified as a less-trusted boundary, M6 may need to promote to Blocker; today the High-confidence Reliability rule downgrade-to-Major holds.

## 8. Missing tests to add now

Grouped by function. Names follow `function_returns_expected_on_condition`.

### `drive_cohort_chunked` / `drive_one_chrom_generic` / `load_and_run_chunk_with_retry` (driver.rs)

These would also satisfy B5's cross-driver byte-identity check.

- `drive_cohort_chunked_matches_streaming_oracle_on_multi_category_fixture` — runs both `drive_cohort_chunked` and `drive_cohort_pipeline` on the same fixture (≥3 samples, ≥1 MNP, ≥1 LH-cap-tripping site, ≥1 hom-REF group, ≥1 below-`min_alt_obs` site, ≥1 below-`qual_phred` site, ≥1 above-`mapq_diff_t` site). Asserts VCF bodies byte-equal (stripping headers) **and** field-by-field equal counter sets. Closes B5.
- `drive_cohort_chunked_leaves_no_tmp_file_when_psp_reader_errors_mid_loop` — injects an error from the second `region_records` call; asserts no leftover tmp file at the expected path. Closes B1.
- `drive_cohort_chunked_zero_length_chromosome_yields_no_records` — input class: `ParsedChromosome.length == 0`. Catches early-return regression at `driver.rs:398-400`.
- `drive_one_chrom_par_iter_matches_sequential_run_window` — same chunk processed with `target_window_count = 4` (parallel) vs `target_window_count = 1` (sequential), forced via `ThreadPoolBuilder`. Asserts byte-equal emitted records. Closes M31.
- `load_and_run_chunk_with_retry_actually_widens_chunk_on_no_safe_gap` — a chunk that satisfies `target_variants_per_chunk` on first attempt but has no safe gap. Documents the B3 bug as a known-fail until fixed.

### `load_chunk_from_iters` (loader.rs)

- `loader_rejects_sample_count_mismatch` — `per_sample_iters.len() != scratch.n_samples()`. Closes M29.
  ```rust
  let mut scratch = ChunkLoadScratch::with_n_samples(2);
  let mut out = MaterialisedChunk::with_n_samples(2);
  let mut carry = vec![SampleColumns::empty(); 2];
  let iters: Vec<_> = vec![Vec::<PileupRecord>::new()] // 1 iter, scratch sized for 2
      .into_iter().map(|rs| rs.into_iter().map(Ok::<_, std::convert::Infallible>)).collect();
  let err = load_chunk_from_iters(&mut scratch, &mut out, 0, 10, 90, 0, 90, iters, &mut carry).unwrap_err();
  assert!(matches!(err, ChunkLoadError::SampleCountMismatch { expected: 2, got: 1 }));
  ```
- `loader_rejects_carryover_length_mismatch` — symmetric on `carryover.len()`. Closes M29.
- `loader_validates_max_span_below_initial_span` — documents (or asserts after-fix) the M8 contract.
- `loader_panics_in_debug_on_duplicate_position_in_one_sample_iter` — `#[should_panic]` pinning the `SampleColumns::push_record` monotonic-position invariant.

### `fix_boundaries` / `pick_safe_end` (pre_pass.rs)

- `pre_pass_rejects_mismatched_carryover_length` — `carryover.len() != chunk.n_samples()`. Closes M30.
- `pre_pass_safe_end_clamps_to_chrom_one_past_end_on_last_chunk` — chunk where `range.end > chrom_one_past_end`. Closes M24.
- `slide_left_to_safe_returns_none_when_desired_equals_min_open` — boundary case.
- `slide_left_to_safe_clamps_above_safe_end` — boundary case.

### `partition_window` (partition.rs)

- `partition_rejects_unsorted_masked_intervals` — `#[should_panic]` (or asserts the new `PartitionError::MaskNotSorted`) on an unsorted mask. Closes M19.
- `samples_at_pos_is_strictly_ascending_within_each_position_range` — proptest over random fixtures. Locks the M32 invariant.

### `prefetch_window_ref_bytes` (worker.rs)

- `prefetch_window_ref_bytes_reuses_inner_vec_capacities_across_calls` — call twice on same partition; capture inner `Vec<u8>.capacity()` after first call; assert ≥ after second. Closes M2.

### `emit_or_drop` and per-category counters (driver.rs)

- `emit_or_drop_increments_low_alt_obs_when_record_fails_min_alt_obs`
- `emit_or_drop_increments_hom_ref_when_is_variant_call_returns_false`
- `emit_or_drop_increments_low_qual_when_qual_phred_below_threshold`
- `emit_or_drop_increments_low_mapq_diff_t_when_filter_triggers`
- `emit_or_drop_writes_record_and_increments_records_written_on_pass`
- `emit_or_drop_filter_order_pins_hom_ref_winning_over_low_qual` — record that *would* trip both; documents which wins. Closes M20.

### `compute_dust_mask_for_chrom` (driver.rs)

- `compute_dust_mask_for_chrom_streams_without_materialising_full_chrom_bases` — fetcher whose `iter_bases` panics if `collect::<Vec<_>>()`-ed; passes after B2 fix.
- `compute_dust_mask_for_chrom_translates_zero_based_to_one_based` — assert 1-based half-open output.

### `passes_min_alt_obs`, `record_fails_mapq_diff_t`, `pool_allele_mapq` (driver.rs)

- `passes_min_alt_obs_returns_true_when_min_obs_is_zero` — boundary.
- `passes_min_alt_obs_returns_true_when_record_has_only_ref_allele` — `n_alleles == 1`.
- `passes_min_alt_obs_returns_true_when_one_alt_meets_threshold_in_one_sample` — max-across-samples wins.
- `record_fails_mapq_diff_t_returns_false_when_threshold_is_nan` — boundary.
- `record_fails_mapq_diff_t_returns_false_when_n_ref_below_minimum` — boundary.
- `record_fails_mapq_diff_t_skips_low_count_alt_and_filters_via_high_count_alt` — multi-ALT.
- `record_fails_mapq_diff_t_matches_streaming_oracle_across_thresholds` — calls both this and `cohort_driver::record_fails_mapq_diff_t_for_test` on the same `PosteriorRecord` for several thresholds (incl. NAN, INFINITY, 0.0); asserts equality. Closes the M20 unit-level oracle gap.

### `enforce_max_alleles_columnar` (kernels/unify_alleles.rs)

- `enforce_max_alleles_ties_match_row_shape_kernel` — 4 ALTs with equal `cohort_count` and `max_alleles = 3`; asserts cn and row-shape pick the same ALT. Closes M22.
- `enforce_max_alleles_drops_everything_when_all_prunable_and_max_is_one` — boundary.
- `enforce_max_alleles_keeps_only_protected_when_max_below_protected_count` — exercise `budget_for_prunable saturating_sub`.

### `detect_compound_candidates_columnar` (kernels/unify_alleles.rs)

- `detect_compound_candidates_invariant_under_sample_permutation` — proptest: random per-sample chain proposals with shuffled orderings; assert candidate set is invariant. Closes the M3 missing-proptest gap.

### `compute_log_likelihoods_columnar` (kernels/compute_log_likelihoods.rs)

- `compute_log_likelihoods_returns_n_alleles_exceeds_bitmask_when_passed_too_many_alleles` — `n_alleles = MAX_BITMASK_ALLELES + 1`. Pairs with B4's new variant.
- `compute_log_likelihoods_returns_degenerate_nan_error_on_synthetic_input` — pin the existing path.
- `chain_broken_log_likelihood_zero_when_no_constituents` — empty `constituent_offsets[lo..hi]`.

### `SampleColumns` (columns.rs)

- `sample_columns_default_returns_csr_sentinels_matching_empty` — pinning M1's `impl Default { Self::empty() }` fix.
- `sample_columns_truncate_to_zero_resets_csr_sentinels` — boundary.
- `sample_columns_push_record_handles_record_with_zero_alleles` — boundary on the CSR invariant.
- `u32_from_usize_panics_on_overflow_in_debug` — `#[test] #[should_panic]`. Pins the M6 debug-assert.

### `WindowPartition` (partition.rs)

- `window_partition_default_matches_empty_helper` — assert `Default::default()` equals `WindowPartition::empty()`.

### `WorkerPool::ensure_capacity` (worker.rs)

- `worker_pool_ensure_capacity_grows_monotonically_and_preserves_existing_slot_buffers` — `ensure_capacity(2)` → `ensure_capacity(4)` → `ensure_capacity(3)`; assert slots 0/1 buffer capacities survive.

### `ColumnarPipelineScratch::take_lh_cap_stats` (worker.rs)

- `take_lh_cap_stats_returns_and_clears_accumulator` — straight unit test.

## 9. What's good

Five specific, transferable patterns worth keeping:

1. **`Send + !Sync` typedef makes the parallel-section soundness statically enforced.** `SharedRefFetcher = Arc<dyn ChromRefFetcher + Send>` (in `per_group_merger.rs`, out of scope but consumed here) means a future edit attempting to capture the fetcher inside the rayon closure would fail to compile, not just fail review. The `unsafe_concurrency` audit returned `No findings.` on the basis of this design choice — that is a deliberate type-system win, worth repeating across the codebase wherever a per-worker single-owner resource exists.
2. **Driver-side sequential pre-fetch (`prefetch_window_ref_bytes`) before rayon dispatch** routes around the `!Sync` fetcher cleanly: the parallel section never touches it. The two-phase shape (sequential partition+prefetch / parallel math / sequential drain) at [driver.rs:596-647](../../../src/var_calling/cohort_block/driver.rs#L596-L647) is the right pattern for "compute that needs shared state but math that doesn't".
3. **Sequential drain in window order** at [driver.rs:640-647](../../../src/var_calling/cohort_block/driver.rs#L640-L647) preserves per-record emit order across parallel execution. This is exactly the byte-identity-friendly shape — rayon does the math in any order, the driver re-serialises emission. Worth pointing at when teaching parallel-pipeline determinism in this codebase.
4. **Every error enum in scope is `#[non_exhaustive]`** (`ChunkLoadError`, `FixBoundariesError`, `PartitionError`, `ChunkDriverError`, `UnifyAllelesError`, `ProjectScalarsError`, `ComputeLogLikelihoodsError`, `WorkerUnifyError`). Consistent error-side armour across a fresh module makes future variant additions safe without follow-up findings.
5. **88 unit tests across the cohort_block sub-modules** exercise the loader/pre-pass/partition/kernels at the function level with a clear `function_returns_expected_on_condition` naming convention. The gaps the review surfaced (B5, M20, M22, etc.) are concentrated at the driver glue + cross-module byte-identity layer, not at the kernels themselves — the unit-test discipline at the leaves is sound.

## 10. Commands to re-verify

Commands the reviewer ran (re-run to confirm they still pass):
- `./scripts/dev.sh cargo fmt --check` (current: exit 0)
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` (current: exit 101 — fixing M11 + the 16 in-scope clippy errors gets this to exit 0; expect that to be a single PR)
- `./scripts/dev.sh cargo test --lib` (current: 1023 passed; B5 + Missing-tests section adds ~25 unit tests, expect ~1048 passed)
- `./scripts/dev.sh cargo test --tests` (current: all pass; B5 adds ≥1 cross-driver integration test)
- `./scripts/dev.sh cargo doc --no-deps` (current: exit 101 — fixing M10 + Mi26 gets this to exit 0)
- `./scripts/dev.sh cargo audit` (not installable in the dev container today; out of band)

New commands introduced by this review:
- `./scripts/dev.sh cargo test --lib var_calling::cohort_block::driver::tests` once the driver-side test module is added.
- `./scripts/dev.sh cargo bench --bench cohort_e2e_perf --baseline main` once Mi24's `// REGRESSION THRESHOLD` comments + a CI baseline exist.

### Author response convention

Address each finding by its identifier (e.g., `B2`, `M14`) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the four open questions in section 4 before responding to the findings they gate.
