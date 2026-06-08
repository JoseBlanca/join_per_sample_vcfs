# Code Review: var_calling

**Date:** 2026-06-08
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** `src/var_calling/` subtree (cohort `.psp` → multi-sample VCF pipeline), excluding the `posterior_engine` module (deferred by PM)
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the `src/var_calling/` module subtree — the re-architected record-streaming cohort `.psp` → VCF pipeline merged to `main` on 2026-06-07.
- **Reviewed against:** branch `main` @ `35a6b67` (working tree carried one uncommitted, out-of-scope change in `pipeline.rs` — the `PVC_QUEUE_DEPTH` env override — explicitly excluded per the requester).
- **In-scope files (14):**
  - [mod.rs](../../../../src/var_calling/mod.rs)
  - [types.rs](../../../../src/var_calling/types.rs)
  - [pipeline.rs](../../../../src/var_calling/pipeline.rs)
  - [sample_reader.rs](../../../../src/var_calling/sample_reader.rs)
  - [cohort_integration.rs](../../../../src/var_calling/cohort_integration.rs)
  - [pileup_overlaps.rs](../../../../src/var_calling/pileup_overlaps.rs)
  - [em_posterior_calc.rs](../../../../src/var_calling/em_posterior_calc.rs)
  - [vcf_writer.rs](../../../../src/var_calling/vcf_writer.rs)
  - [per_position_merger.rs](../../../../src/var_calling/per_position_merger.rs) *(verbatim kernel — spot-checked)*
  - [variant_grouping.rs](../../../../src/var_calling/variant_grouping.rs) *(verbatim kernel — spot-checked)*
  - [dust_filter.rs](../../../../src/var_calling/dust_filter.rs) *(verbatim kernel — spot-checked)*
  - [contamination_estimation.rs](../../../../src/var_calling/contamination_estimation.rs) *(verbatim kernel — spot-checked)*
  - [per_group_merger.rs](../../../../src/var_calling/per_group_merger.rs) *(verbatim kernel — spot-checked)*
  - [test_helpers.rs](../../../../src/var_calling/test_helpers.rs)
- **Deliberately out of scope:**
  - `src/var_calling/posterior_engine.rs` + `posterior_engine/{backends,interp,shape}.rs` — PM deferred to a later review.
  - The uncommitted `PVC_QUEUE_DEPTH` block in `pipeline.rs` (~lines 257-265) — ongoing work from another session.
- **Categories dispatched (10):** reliability, errors, naming, defaults, idiomatic, refactor_safety, module_structure, unsafe_concurrency, smells, extras. `tooling` was skipped (subtree review; the build config — `lto=fat`/`codegen-units=1`/`panic=abort`/pinned `target-cpu` — was assessed in the prior perf reviews and is unchanged).

The rebuilt-structure files (the new integration surface) were reviewed thoroughly; the six numeric kernels are carried verbatim from the prior architecture and were each individually reviewed mid-May, so they were spot-checked at the seams for Major+ only.

## 2. Verdict

**Approve-with-changes.** 0 Blockers, 8 Major, 13 Minor, grouped Nits.

No correctness Blockers in the emitted calls: all gates are green, 1052 lib + integration tests pass, and `unsafe_concurrency` re-confirmed the parallel topology is sound (no `unsafe`; deadlock-free bounded channels; order-independent parallel fold; `VariantCaller: Sync`). The byte-identical-VCF contract is verified out-of-tree (vs the prior architecture + the GIAB benchmark) per the standing PM decision. The verdict is driven by quality/maintainability debt that should be paid down now while the module is fresh — dead public API at the central producer→caller seam (M1), an over-broad public surface (M2), a typed-error-flattening regression (M3), a stage-coupling back-reference (M4), and stale migration-era documentation whose "copied verbatim from <module>" provenance now points at deleted code (M7) — plus a real test-coverage gap on the byte-identity contract and its release-level guards (M5, M8). None blocks correctness; hence Approve-with-changes rather than Request-changes.

## 3. Execution status

All commands run in the project dev container (`./scripts/dev.sh`) per `CLAUDE.md`:

- `cargo fmt --check` → **exit 0**, clean.
- `cargo clippy --all-targets --all-features -- -D warnings` → **exit 0**, no warnings.
- `cargo doc --no-deps` → **exit 0**, clean (no broken intra-doc links).
- `cargo test --lib --tests` (the CI gate) → **exit 0**, **1052 passed, 0 failed** across lib + integration binaries.
- `cargo test --all-targets --all-features` → **fails**, but only at the `psp_writer_perf` bench (`benches/psp_writer_perf.rs:386` panic). This is the **pre-existing, out-of-scope** failure already documented in `PROJECT_STATUS.md` ("breaks `cargo test --all-targets`; not the CI gate, which uses `--lib --tests`"). Not attributable to this subtree.
- `cargo audit` → not run.

Findings labeled "Needs verification": 1 (M7 — the verbatim-copy claims cannot be diffed in-tree because their source modules are deleted; verification must go against git history / the out-of-tree GIAB diff).

## 4. Open questions and assumptions

1. **Is the deferred-until-P7 documentation debt ready to be resolved now?** The `mod.rs` header, `types.rs`, and several kernel doc comments still describe the subtree as "built alongside the shipping `var_calling` until the one-commit swap (§P7)". The swap has happened — this *is* the only `var_calling`. Resolving this unblocks M7, Mi4, Mi5, Mi6 and the `Variant` type-alias decision. *(Affects M7, Mi4, Mi5.)*
2. **Should `PileupCohortChunk` be deleted or is it a deliberate not-yet-wired target?** It is dead today (no constructor/consumer). The design notes treat the columnar `RawCohortChunk` path as settled, which implies deletion. *(Affects M1, Mi5.)*
3. **Is the out-of-tree byte-identity oracle considered sufficient, or should an in-tree worker-count A/B test be added?** Prior PM decision (deferred M7 on the 2026-06-05 review) said out-of-tree verification suffices. M5/M8 ask whether at least the writer-reorder mechanism and the release-level guards should get cheap in-tree unit coverage. *(Affects M5, M8.)*
4. **Should the five rebuilt-structure modules be `pub(crate)`?** They have zero external consumers today; the only blocker is `vcf_writer::WriterStats` named in a CLI struct field. *(Affects M2.)*

## 5. Top 3 priorities

1. **M3 — restore typed config errors in `pipeline.rs`.** `run_var_calling` funnels every typed config-builder error through `cfg_err = |e| PipelineError::Config(e.to_string())`, flattening the cause and breaking the `source()` chain — a diagnosability regression. Cheap, mechanical fix (operation-named variants with `#[source]`).
2. **M1 — delete the dead `PileupCohortChunk` type.** A `pub` type documented as "the producer→caller work-unit" that is never constructed or consumed sits at the exact seam the re-architecture set out to clarify, and creates the `CohortPileupRecord`/`PileupCohortChunk` near-anagram. Convergently flagged by four categories.
3. **M5 — add a byte-identical-across-worker-counts test.** The subtree's hardest contract ("byte-identical for any worker count") has no test that varies `threads`; every integration test runs `threads: None`. The writer's whole reorder machine is unguarded at the dimension it exists for.

## 6. Findings

### Major

**M1: [types.rs:93-107](../../../../src/var_calling/types.rs#L93-L107) — `PileupCohortChunk` is dead public API at the central producer→caller seam**
**Categories:** idiomatic, smells, module_structure, naming (convergent)
**Confidence:** High.
`PileupCohortChunk` is a `pub struct` documented as "the producer→caller work-unit", but a crate-wide grep (`src/`, `tests/`, `benches/`, `examples/`) finds **only** its definition plus three doc-comment cross-references ([types.rs:10](../../../../src/var_calling/types.rs#L10), [:126](../../../../src/var_calling/types.rs#L126), [:173](../../../../src/var_calling/types.rs#L173)) — no constructor, no consumer, no test. The work-unit actually shipped over the channel is the columnar `RawCohortChunk` ([pipeline.rs:272](../../../../src/var_calling/pipeline.rs#L272), [cohort_integration.rs:923](../../../../src/var_calling/cohort_integration.rs#L923), [em_posterior_calc.rs:112](../../../../src/var_calling/em_posterior_calc.rs#L112)). A reader sees two structs both described as "the producer→caller work-unit" and must grep to learn which is live; the `CallStats`/`CalledChunk` docs even cross-reference the dead one as if it were the source type.
**Fix:** Delete `PileupCohortChunk` and retarget the three doc references at `RawCohortChunk` (its `chunk_order` is the field they describe). This also dissolves the Mi5 naming clash. If it is genuinely a near-term target shape, gate it behind an explicit `// UNUSED: pre-columnar variant` marker so it cannot rot silently.

**M2: [mod.rs:39-53](../../../../src/var_calling/mod.rs#L39-L53) — five rebuilt-structure modules are `pub` with zero external consumers (extends prior M8)**
**Categories:** module_structure, idiomatic (convergent)
**Confidence:** High.
All 13 modules are declared `pub mod`. An external-consumer audit (`grep crate::var_calling::<mod>` outside the subtree) shows `types`, `sample_reader`, `cohort_integration`, `pileup_overlaps`, and `em_posterior_calc` have **zero** references outside `src/var_calling/` — they are used only by `pipeline` and by each other. Making them `pub` puts `CohortChunkIntegrator`, `VariantCaller`, `SamplePspChunk`/`TwoPhaseSegment`, and every interchange struct on the crate's public API surface (and into rustdoc) for no caller. Within `sample_reader` specifically, several getters/constructors (`records_for`, `take_scalar`/`take_seq`/`take_chain_ids`, `n_alleles_at`, `RefSpan::empty`) are reachable only from `#[cfg(test)]`.
**Fix:** Demote the five internal modules to `pub(crate) mod`. `pipeline` stays `pub` (external entry); `vcf_writer` stays `pub` only because `pop_var_caller/var_calling.rs` names `crate::var_calling::vcf_writer::WriterStats` in a field — consider re-exporting `WriterStats` from `pipeline` and demoting `vcf_writer` too. The six kernels keep `pub` (all have external consumers). Narrow the test-only `sample_reader` getters to `pub(crate)` (see also M6).

**M3: [pipeline.rs:149-176](../../../../src/var_calling/pipeline.rs#L149-L176),[:197](../../../../src/var_calling/pipeline.rs#L197) — config-construction errors flattened to `PipelineError::Config(String)`, dropping the typed cause and `source()` chain**
**Categories:** errors
**Confidence:** High.
Every per-stage config builder returns a typed error (`GrouperConfigError`, `PerGroupMergerConfigError`, `PosteriorEngineConfigError`, the `--regions` BED error, the rayon pool-build error), but all are funnelled through `let cfg_err = |e: &dyn std::fmt::Display| PipelineError::Config(e.to_string());`. The typed cause is stringified into the parent's `Display` and never attached via `#[source]`, so the error chain is broken and a caller cannot match on the underlying validation variant. `Config(String)` is a mechanism-named `Other(String)`-shaped catch-all that several unrelated failure modes collapse into. A bad `--em-convergence-threshold`, a bad `--var-group-max-span`, a malformed `--regions` BED, and a pool-build failure all surface as one opaque stringified `Config`.
**Fix:** Give `PipelineError` operation-named variants carrying each cause via `#[source]` (`GrouperConfig(#[source] GrouperConfigError)`, `MergerConfig(...)`, `PosteriorConfig(...)`, `Regions(...)`, `ProducerPool(#[source] rayon::ThreadPoolBuildError)`), and map each `?` site to its specific variant instead of the shared stringifier.

**M4: [em_posterior_calc.rs:31](../../../../src/var_calling/em_posterior_calc.rs#L31),[:120](../../../../src/var_calling/em_posterior_calc.rs#L120) — the caller stage imports `merge_compacted_samples` from the producer module (caller→producer back-reference)**
**Categories:** module_structure
**Confidence:** High.
`em_posterior_calc` (the **caller**) imports and calls `merge_compacted_samples`, a `pub fn` defined in [cohort_integration.rs:438](../../../../src/var_calling/cohort_integration.rs#L438) (the **producer**). The runtime topology is producer → caller → writer, so this is a backward edge — and the caller's own doc says this is "the record-building work moved off the single producer thread onto this parallel caller", i.e. the work conceptually belongs to the caller. It lives in the producer module only because the helper types were authored there. It is the sole edge keeping the two stage modules cyclically coupled, and blocks demoting `cohort_integration` to a narrower surface (M2).
**Fix:** Move `merge_compacted_samples` (plus `ok_record` / `KeptRecordIter`, [cohort_integration.rs:387-426](../../../../src/var_calling/cohort_integration.rs#L387-L426)) into the caller. It depends only on `SamplePspChunk::records_all`, `PerPositionMerger`, and `CohortPileupRecord` — all importable from the caller. After the move, the only cross-edge between the stage modules is the `RawCohortChunk` type (correctly in `types`).

**M5: [tests/cohort_cli_integration.rs](../../../../tests/cohort_cli_integration.rs) / [vcf_writer.rs:102](../../../../src/var_calling/vcf_writer.rs#L102) — the byte-identical-for-any-worker-count contract is untested at the dimension it is defined over**
**Categories:** extras, reliability (convergent)
**Confidence:** High.
The contract is stated plainly ("The output is byte-identical for any worker count", [pipeline.rs:16](../../../../src/var_calling/pipeline.rs#L16)) and the writer's entire `BTreeMap` reorder + gapless-`chunk_order` drain ([vcf_writer.rs:102-115](../../../../src/var_calling/vcf_writer.rs#L102-L115)) exists only to honour it. But every integration test that drives the pipeline sets `threads: None` — none varies `W` and diffs the resulting VCFs. There is a two-run determinism test and a target-variants-per-chunk byte-identity test, but the caller→writer reorder under genuine out-of-order completion is never driven concurrently. A reorder/cursor regression (off-by-one in `next_expected`, emit-in-arrival-order) would pass the whole suite and surface only as a non-deterministic field diff in production. `vcf_writer.rs` has zero in-module tests.
**Fix:** Add `var_calling_byte_identical_across_worker_counts` driving the same cohort at `threads ∈ {Some(1), Some(2), Some(8)}` and asserting each VCF body equals the `Some(1)` run (the natural extension of the existing per-chunk-target byte-identity test); and a direct `VcfWriter` unit test submitting `CalledChunk`s `2,0,1` and asserting emit order `0,1,2`. Byte-identity is verified out-of-tree, so this is defense-in-depth, but it is cheap and targets the exact rewrite-sensitive mechanism.

**M6: [sample_reader.rs:447](../../../../src/var_calling/sample_reader.rs#L447),[:381](../../../../src/var_calling/sample_reader.rs#L381),[:514-584](../../../../src/var_calling/sample_reader.rs#L514-L584),[:22-32](../../../../src/var_calling/sample_reader.rs#L22-L32) — dead / mislabeled `SamplePspChunk` decode API; the module doc misdirects on the live decode path**
**Categories:** smells, idiomatic, reliability (convergent)
**Confidence:** High.
The module doc presents the move-out take-getters (`take_seq`/`take_chain_ids`/`take_scalar`) as the production "Phase 5 column-selective seam" and `from_block`/`records_for` as the "Phase 1 simple decode" production path. In the shipped wiring, production reads **exclusively** via `next_two_phase` → `TwoPhaseSegment::set_variable_rows` → `records_all`. A usage scan shows `append_kept` ([:447](../../../../src/var_calling/sample_reader.rs#L447) — confirmed: only its own definition + one doc link) and `SamplePspChunk::chrom_id()` have **zero** callers; `take_*`/`records_for`/`n_alleles_at`/`from_block`/`next_chunk` are **test-only** (the eager path now serves as the byte-identity oracle). So the doc describes a production decode path the code no longer uses, over a large `pub` surface that must stay compiling/tested for no production benefit.
**Fix:** (1) Delete the genuinely unused `append_kept` and `SamplePspChunk::chrom_id()`. (2) For the test-only oracle methods, relabel honestly — gate behind `#[cfg(test)]` where callers are tests only, and rewrite the module doc to say the eager `from_block` path is the **byte-identity oracle** and `next_two_phase` is the production decode. Verify with a non-test `cargo build` after gating.

**M7: [cohort_integration.rs](../../../../src/var_calling/cohort_integration.rs) / [mod.rs:1-36](../../../../src/var_calling/mod.rs#L1-L36) — stale migration-era docs; the "copied verbatim from <module>" provenance now points at deleted code and the in-tree byte-identity safety net is unreachable**
**Categories:** refactor_safety, module_structure, naming, smells (convergent)
**Confidence:** High (verbatim sources confirmed deleted); the byte-identity-equivalence claim itself is **Needs verification**.
Many byte-identity-critical functions carry "copied verbatim from X" doc comments citing the *old* package — `two_pass::{reach,WindowSummary::*}`, `loader::find_block_cut`, `driver::{merge_block_ranges,emit_or_drop,record_fails_mapq_diff_t,dust_mask_for_interval,restrict_intervals_to_regions}`, `worker::columnar_passes_min_alt_obs`. None of those modules exist under `src/` anymore (confirmed: `find` for `two_pass.rs`/`worker.rs`/`loader.rs` → none; only the unrelated `src/pileup/walker/driver.rs`). The P7 swap the `mod.rs` header still frames as future has already happened — `lib.rs` exposes exactly one `var_calling`, and `[var_calling](crate::var_calling)` is a self-link. So the only in-repo defense against silent drift in these math seams — "this is a verbatim copy of <source>" — can no longer be diffed against anything, and there is no in-tree A/B oracle (byte-identity is confirmed out-of-tree only).
**Fix:** (1) Replace each dangling "copied verbatim from `two_pass`/`driver`/`worker`/`loader`" provenance with the **commit hash** the body was copied at, so the claim stays auditable against git history. (2) Rewrite the `mod.rs`/`types.rs`/kernel headers from migration future-tense ("built alongside … until the P7 swap", "Phase N builds this", "becomes the real type at the P7 swap") to present-tense descriptions of the production pipeline, and resolve the deferred-until-swap decisions (the `Variant` alias, the stale `doc/devel/...` doc links vs `ia/specs/`). (3) Consider freezing the seam functions behind a fixed-fixture regression test that pins their output, giving an in-tree drift signal.

**M8: [vcf_writer.rs:146](../../../../src/var_calling/vcf_writer.rs#L146),[:120](../../../../src/var_calling/vcf_writer.rs#L120) / [cohort_integration.rs:878](../../../../src/var_calling/cohort_integration.rs#L878),[:784](../../../../src/var_calling/cohort_integration.rs#L784) / [pileup_overlaps.rs:39](../../../../src/var_calling/pileup_overlaps.rs#L39) / [pipeline.rs:436](../../../../src/var_calling/pipeline.rs#L436),[:480](../../../../src/var_calling/pipeline.rs#L480) — release-level guards and byte-identity helpers have no targeted tests**
**Categories:** reliability
**Confidence:** High.
A cluster of byte-identity-critical / data-loss-guarding paths in the rebuilt files have no isolating test (the modules `pipeline.rs`, `vcf_writer.rs`, `pileup_overlaps.rs` have zero in-module tests):
- `WriterError::MissingChunks` — the gapless-invariant guard against silent VCF truncation; a regression returning `Ok` on a non-empty reorder buffer would truncate output silently.
- `ProducerError::StalledCut` — the anti-spin guard; a regression that lets `cut` not advance would *hang* (not panic, not fail a test).
- `emit_or_drop` — the load-bearing filter order (hom-ref → QUAL → MAPQ-t) and four `WriterStats` counters, plus the unconverged-but-written branch; only `record_fails_mapq_diff_t` is tested in isolation.
- `overlapping_groups` — the only function in `pileup_overlaps.rs`, the `CohortPileupRecord`→`PerPositionPileups` re-wrap feeding the grouper.
- `restrict_intervals_to_regions` + `dust_mask_for_interval` — the `--regions` clip arithmetic and the sub-span DUST coalesce, both claiming byte-identity to a single whole-interval scan.
- `rebuild_fold` — the parallel reduce whose byte-identity rests on `CohortSpanFold::merge` being associative+commutative; the existing `merge` test uses distinct values, so a tie-break regression would slip through.
- `compact_samples` straddler (`Pending` partial-inflate) branch — whether a straddling segment's partial-inflate columns are byte-identical to the eventual `Ready` finalisation is the crux and is not isolated.
**Fix:** Add the targeted unit tests listed in §8. These are cheap and pin the exact mechanisms the rewrite is most likely to regress; the reliability sub-agent supplied concrete specs for each.

### Minor

**Mi1: [pipeline.rs:252-271](../../../../src/var_calling/pipeline.rs#L252-L271) — startup log records resolved knob values but not their source.** The `eprintln!` emits `producer_threads`/`workers`/`queue_cap`/`target_variants_per_chunk` but an operator cannot tell whether a value came from a `PVC_*` env var, `--threads`, or the built-in default — the precise question the env-dependent-defaults rule says logs must answer. A mistyped `PVC_PRODUCERTHREADS` would look identical to an honoured one, making perf-sweep attribution unreliable. *(Tag each value with its provenance; have `resolve_thread_split` return the source alongside the value.)* The `defaults` sub-agent rated this Major; filed at Minor as it is observability, not behavior.

**Mi2: [pipeline.rs:43-44](../../../../src/var_calling/pipeline.rs#L43-L44) — `BUFFERED_IO_CAPACITY` const-via-comment drift.** Documented to "match `pop_var_caller::common::DEFAULT_BUFFERED_IO_CAPACITY`" but hard-codes `64 * 1024`; that const is `pub(crate)` and importable, so nothing forces them equal. **Categories:** smells, defaults, idiomatic. *(Reference the canonical const directly.)*

**Mi3: [vcf_writer.rs:26-31](../../../../src/var_calling/vcf_writer.rs#L26-L31),[:155-160](../../../../src/var_calling/vcf_writer.rs#L155-L160) — co-dependent filter fields.** `DownstreamFilters.no_mapq_diff_filter: bool` gates whether `min_mapq_diff_t: f32` means anything; the state `(no_mapq_diff_filter=true, min_mapq_diff_t=<any>)` carries a dead field. **Categories:** smells. *(Collapse to `enum MapqDiffFilter { Off, On { min_t: f32 } }`, byte-identity-neutral.)*

**Mi4: [em_posterior_calc.rs:1](../../../../src/var_calling/em_posterior_calc.rs#L1) — module named after a delegated sub-step.** The file hosts the `VariantCaller` worker (grouping + per-group merge + per-position merge + EM), not just the EM. **Categories:** naming, smells (re-confirms prior Mi11). *(Rename to `variant_caller.rs`; update `mod.rs:48` + importers.)*

**Mi5: [types.rs:38](../../../../src/var_calling/types.rs#L38),[:102](../../../../src/var_calling/types.rs#L102) — `CohortPileupRecord` / `PileupCohortChunk` near-anagram.** Same three tokens permuted on the two central data types; the single-position-vs-list distinction rides on the easily-missed final token. **Categories:** naming (re-confirms prior Mi12). Largely resolved by deleting the dead `PileupCohortChunk` (M1); otherwise rename the work-unit so its kind leads.

**Mi6: [types.rs:64-70](../../../../src/var_calling/types.rs#L64-L70) — `RefSpan::empty()` dead constructor with a fictional recycling doc.** Called only by its own unit test; the doc describes a "buffer recycling (the producer refills … before shipping)" protocol the producer does not implement (it builds a fresh `RefSpan` per chunk). **Categories:** smells. *(Delete it — callers can use `RefSpan::default()` — or rewrite the doc to drop the fictional narrative.)*

**Mi7: [vcf_writer.rs:130-136](../../../../src/var_calling/vcf_writer.rs#L130-L136) — `CallStats`→`WriterStats` roll-up is field-access, not exhaustive destructure.** Correct today (all four counters covered), but adding a fifth `CallStats` counter would compile and silently never reach the run summary. **Categories:** refactor_safety. *(Destructure `let CallStats { .. } = stats;` so a new field is compiler-flagged.)*

**Mi8: [pipeline.rs:363](../../../../src/var_calling/pipeline.rs#L363),[:393](../../../../src/var_calling/pipeline.rs#L393),[:419](../../../../src/var_calling/pipeline.rs#L419) / [cohort_integration.rs:1019](../../../../src/var_calling/cohort_integration.rs#L1019),[:974](../../../../src/var_calling/cohort_integration.rs#L974) — release-path `.expect()`s lack `// PANIC-FREE:` / `// UNREACHABLE:` comments.** The `scope` join `.expect("… thread panicked")`, the `ref_fetch`/`dust_mask_for` `expect("just set")` cache re-lookups, the `fetch_ref_span` `binary_search(...).expect("variable position must be in the fold")`, and the `new_column_decompressor().expect("… infallible …")` are all genuinely unreachable/justified but undocumented as such. **Categories:** errors, refactor_safety, unsafe_concurrency, idiomatic, reliability (convergent). *(Add one-line invariant comments; consider returning the `&mut` from the `get_or_insert` branch to drop the re-lookup.)*

**Mi9: [pipeline.rs:362-367](../../../../src/var_calling/pipeline.rs#L362-L367) — a caller/writer-thread panic aborts the run via `.expect` instead of surfacing a typed `PipelineError`.** Not a deadlock (the panicking thread's channel handles drop on unwind, so survivors drain and join — verified), but a single data-dependent panic in one chunk becomes a non-recoverable crash with a generic message and no chunk context, bypassing the "first error wins" graceful path. **Categories:** unsafe_concurrency. *(Either document it as fatal-by-design + deadlock-free, or capture the panic into `first_err` via `match handle.join()`.)*

**Mi10: [pipeline.rs:205](../../../../src/var_calling/pipeline.rs#L205) — `SamplePspReader::new(r, 0, 1, 1)` passes undocumented magic placeholder args.** The `(chrom_id=0, region_start=1, region_end=1)` literals are inert placeholders reset by `begin_interval` before any read, but read as behavioural defaults at the call site. **Categories:** defaults. *(Add a one-line comment, or a `new_unpositioned(reader)` constructor.)*

**Mi11: [sample_reader.rs:729](../../../../src/var_calling/sample_reader.rs#L729) — single-letter look-alike scratch buffers `sd`/`so`/`cd`/`co` in `set_variable_rows`.** Four same-shaped `Vec`s named by two-letter initialisms spanning the function's longest block; a transposed pair would be a silent bug. **Categories:** naming. *(Qualify to `seq_data`/`seq_offsets`/`chain_data`/`chain_offsets`.)*

**Mi12: [benches/cohort_var_calling_perf.rs](../../../../benches/cohort_var_calling_perf.rs) — hot-path bench has no encoded regression threshold or CI baseline gate.** Well-documented and asserts `records_written`, but no `// REGRESSION THRESHOLD: N%` marker and no automated comparison against `main`, so a producer-fold/compaction regression (the exact code these reviews keep tuning) lands silently. **Categories:** extras. *(Add a threshold marker per group and a CI baseline step, or document why the manual+memory workflow is deliberate.)*

**Mi13: [em_posterior_calc.rs:147-152](../../../../src/var_calling/em_posterior_calc.rs#L147-L152) — `passes_min_alt_obs` called with positional `MergedRecord` fields.** The `n_alleles == record.alleles.len()` coupling that the flat `scalars[sample*n_alleles+allele]` index depends on is convention-only; a layout change would still typecheck. **Categories:** refactor_safety. *(Move the predicate onto `MergedRecord` so the coupling lives next to the layout.)*

### Nits

- [cohort_integration.rs:582](../../../../src/var_calling/cohort_integration.rs#L582) — bare-letter field `n` for sample count; `n_samples` matches the domain term used elsewhere.
- [cohort_integration.rs:806](../../../../src/var_calling/cohort_integration.rs#L806) — `.reduce(CohortSpanFold::new, …)` where `new` is just `default()`; `CohortSpanFold::default` (or dropping the `new` wrapper) is less ceremony.
- [pipeline.rs:149](../../../../src/var_calling/pipeline.rs#L149) — the per-site `|e| cfg_err(&e)` closures could be a generic `fn cfg_err<E: Display>(e: E) -> PipelineError` used directly as `.map_err(cfg_err)` (folds into M3's rework).
- [sample_reader.rs:66-73](../../../../src/var_calling/sample_reader.rs#L66-L73) — `len()`/`is_empty()` column accessors could carry `#[must_use]` if the project enables `must_use_candidate`.
- [sample_reader.rs:539-583](../../../../src/var_calling/sample_reader.rs#L539-L583) — `take_scalar`/`take_seq`/`take_chain_ids` share a keep-walk duplicated a third time in `set_variable_rows`; factor the keep→range computation into one helper (folds into M6's cleanup).
- [pipeline.rs:126-373](../../../../src/var_calling/pipeline.rs#L126-L373) — `run_var_calling` is ~247 lines; the three config-construction blocks are natural extract-to-helper boundaries that would bring it under ~100 lines.
- [pipeline.rs:76-86](../../../../src/var_calling/pipeline.rs#L76-L86) — `resolve_thread_split`'s `PVC_PRODUCER_THREADS`/`PVC_CALLER_THREADS` experimental knobs have a stated removal condition but no tracking issue/owner reference.
- Recurring `dz`/`cs`/`ds` scratch-pair locals could match the parameter names (`compressed_scratch`/`decompressed_scratch`) they bind to.

## 7. Out of scope observations

- **[dust_filter.rs:891-903](../../../../src/var_calling/dust_filter.rs#L891-L903)** *(verbatim kernel)* — `is_masked` does `debug_assert!(false, …)` then returns `false` in release if the mask is unloaded, which would pass a record through unmasked (silently wrong output) on a release-build contract violation. Worth a correctness follow-up on the verbatim kernel; out of this review's seam scope.
- **[sample_reader.rs:930](../../../../src/var_calling/sample_reader.rs#L930)** (`advance_pos`) — `u32::try_from(delta)` errors on `u64` overflow but the subsequent `saturating_add` caps a position at `u32::MAX` rather than erroring. This mirrors the canonical row reader verbatim (`src/psp/reader.rs` Mi5), so it is intended behaviour, not a divergence — flagged only to confirm the saturate-vs-error choice is desired for the new path.
- **`benches/psp_writer_perf.rs:386`** — the bench panic that breaks `cargo test --all-targets` is pre-existing and already tracked in `PROJECT_STATUS.md`; not attributable to this subtree.
- **`posterior_engine` module** — deferred by the PM to a later review; not assessed.

## 8. Missing tests to add now

Grouped by function under test. (The `reliability` sub-agent supplied full specs; condensed here.)

**`VcfWriter` (vcf_writer.rs — currently zero in-module tests)**
- `handle_emits_in_chunk_order_when_chunks_arrive_out_of_order` — submit `CalledChunk`s with `chunk_order` 2, 0, 1; assert records reach the inner writer in 0,1,2 order. Catches: emit-in-arrival-order / wrong-cursor reorder regression (silent VCF mis-ordering). *(M5/M8)*
- `finish_errors_with_missing_chunks_on_gap` — `handle(0)`, `handle(2)`, then `finish()` → `Err(WriterError::MissingChunks { count: 1 })`. Catches: silent VCF truncation guard returning `Ok` on a gap. *(M8)*
- `emit_or_drop_applies_filters_in_order_and_counts_each_bucket` — one record per drop reason (hom-ref, sub-QUAL, MAPQ-t-failing) + one clean + one clean-but-unconverged; assert `records_written == 2` and each `records_dropped_*`/`records_unconverged` counter. Catches: reordered gates, dropped-instead-of-written unconverged record, mis-counted buckets. *(M8)*

**`run_var_calling` (pipeline.rs — integration)**
- `var_calling_byte_identical_across_worker_counts` — same cohort at `threads ∈ {Some(1), Some(2), Some(8)}`; assert each VCF body equals the `Some(1)` run. Catches: any worker-count-dependent divergence in the reorder/gapless machinery. *(M5)*

**`pileup_overlaps::overlapping_groups`**
- `overlapping_groups_bridges_dropped_positions_and_splits_on_gap` — overlapping reaches + a gap past `var_group_max_span` + a pure-REF position inside an open group; assert the `OverlappingVariantGroup` boundaries and that the pure-REF position folds in. Catches: bad re-wrap / mis-passed `GrouperConfig`. *(M8)*

**`pipeline::restrict_intervals_to_regions` + `dust_mask_for_interval`**
- `restrict_intervals_to_regions_clips_at_boundaries` — region inside an interval / straddling two / touching the exclusive end / empty regions (table-driven). Catches: off-by-one in the 1-based-inclusive→half-open conversion. *(M8)*
- `dust_mask_for_interval_subspan_invariant` — a masked run straddling a `DUST_SUBSPAN` boundary; assert small-subspan result equals the `subspan = u32::MAX` single scan. Catches: a coalescing bug splitting one run at the seam. *(M8)*

**`cohort_integration` (producer)**
- `rebuild_fold_independent_of_reduce_grouping` — ≥3 folds with *tied* `max_ref_span`/`max_nonref_obs` at shared positions reduced as `((a∘b)∘c)∘d` vs `(a∘b)∘(c∘d)`; assert identical. Catches: a non-commutative `merge` tie-break the distinct-value `merge` test misses. *(M8)*
- `compact_samples_straddler_matches_finalized` — one segment straddling a cut (force via `target_variants = 1`); assert the straddling chunk's `per_sample` records equal the keep-all finalisation restricted to the same range. Catches: partial-inflate diverging from `Ready` finalisation. *(M8)*
- `produce_chunk_with_zero_samples_yields_none` — empty reader vec; assert `Ok(None)`. Catches: a future watermark/`all_exhausted` default change turning n=0 into a spin/panic.
- A `StalledCut` boundary test (hand-built fold where `find_cut` returns `next_chunk_start`) asserting `Err(ProducerError::StalledCut)` rather than a hang. *(M8)*

**`types::RefSpan`**
- `ref_span_slice_below_start_panics` — `#[should_panic]` slicing with `group_start < genomic_start` (documents the release-mode precondition-violation panic).

## 9. What's good

- **Channel topology is provably deadlock-free and panic-safe.** Both bounded hand-offs have `cap ≥ 1`, handle ownership/close-ordering is clean (main drops its handles post-spawn, producer drops `chunk_tx` on return, callers break on `send().is_err()`), and even a producer/worker panic out of `producer_pool.install` unwinds `chunk_tx` before the join — verified by the `unsafe_concurrency` pass with no findings. ([pipeline.rs:276-372](../../../../src/var_calling/pipeline.rs#L276-L372))
- **Byte-identity-under-parallelism is structurally enforced, not hoped for.** The writer's `BTreeMap`+`next_expected` reorder makes emit order independent of worker finish order, and the producer's parallel fold reduces via an associative+commutative integer-`max`/position-union `merge`, so any rayon reduce tree yields the same summary. ([vcf_writer.rs:102-115](../../../../src/var_calling/vcf_writer.rs#L102-L115), [cohort_integration.rs:205-243](../../../../src/var_calling/cohort_integration.rs#L205-L243))
- **`VariantCaller` is stateless-by-construction and the type system carries it.** `call_chunk`/`call_records` take `&self`, build all per-call state locally, and route the merge through the free `merge_group_with_ref` (shared-ref tables) rather than the `Arc`/atomic-bearing `PerGroupMerger`, so `&caller` is soundly shared across threads. ([em_posterior_calc.rs:78-185](../../../../src/var_calling/em_posterior_calc.rs#L78-L185))
- **Release-load-bearing invariants were correctly promoted from `debug_assert!` to typed errors.** `WriterError::MissingChunks` (anti-truncation) and `ProducerError::StalledCut` (anti-spin) are real release guards, not debug-only — exactly addressing the prior review's M6. ([vcf_writer.rs:120-128](../../../../src/var_calling/vcf_writer.rs#L120-L128), [cohort_integration.rs:878-883](../../../../src/var_calling/cohort_integration.rs#L878-L883))
- **Error-type design across the new surface is strong.** `ProducerError`, `CallerError`, `WriterError`, `PerPositionMergerError` use operation-named variants with single-origin `#[from]`/`#[source]`, `#[non_exhaustive]` on the kernel publics, and no foreign-type leaks (the producer's `RefFetchError = Box<dyn Error>` is an internal closure-return type). The lone exception is the config-flattening in `pipeline.rs` (M3).

## 10. Commands to re-verify

Reviewer ran (in the dev container; re-run to confirm they still pass):
```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo doc --no-deps
./scripts/dev.sh cargo test --lib --tests          # CI gate: 1052 passed, 0 failed
```
New commands the review introduces (after applying fixes):
```
# M5 / M8 — the new byte-identity + guard tests
./scripts/dev.sh cargo test --lib --tests var_calling
./scripts/dev.sh cargo test --test cohort_cli_integration var_calling_byte_identical_across_worker_counts
# M7 verification (out-of-tree, against git history)
git log --oneline -- src/var_calling/   # confirm the commit the verbatim kernels were copied at
# Mi12 (optional) — encode + run the hot-path bench baseline
./scripts/dev.sh cargo bench --bench cohort_var_calling_perf
```

### Author response convention
Address each finding by its identifier (e.g., "M3", "Mi8") with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the §4 open questions first — they gate M1, M2, M5, M7, Mi4, Mi5.
