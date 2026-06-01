# Code Review: var_calling

**Date:** 2026-06-01
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Full module review of `src/var_calling/` (Stages 3–6 of the cohort SNP caller), 23 files / 27,489 LoC
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the whole `src/var_calling/` subtree — the DUST filter, per-position merger, variant grouper, per-group merger, posterior engine (+ `backends`/`interp`/`shape`), contamination estimator, and the `from_psp/` chunk-at-a-time cohort driver (driver, loader, columns, column_span_reader, partition, worker, kernels).
- **Reviewed against:** commit `3c9ebf2` on `main`.
- **Review style:** *full but prioritized* — prioritized the never-deeply-reviewed `from_psp/` integration seams (driver↔loader↔worker↔partition handoffs, the post-reorg public surface, the `from_bam`-removal fallout) and spot-checked the stage internals that were each individually reviewed in mid-May 2026 (`posterior_engine`, `per_group_merger`, `per_position_merger`, `dust_filter`, `variant_grouping`, `contamination_estimation`). The `from_psp/` subtree itself was reviewed 2026-05-29 under its old name `cohort_block`.
- **In-scope files:** [mod.rs](../../../../src/var_calling/mod.rs), [dust_filter.rs](../../../../src/var_calling/dust_filter.rs), [per_position_merger.rs](../../../../src/var_calling/per_position_merger.rs), [variant_grouping.rs](../../../../src/var_calling/variant_grouping.rs), [per_group_merger.rs](../../../../src/var_calling/per_group_merger.rs), [posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs), [posterior_engine/backends.rs](../../../../src/var_calling/posterior_engine/backends.rs), [posterior_engine/interp.rs](../../../../src/var_calling/posterior_engine/interp.rs), [posterior_engine/shape.rs](../../../../src/var_calling/posterior_engine/shape.rs), [contamination_estimation.rs](../../../../src/var_calling/contamination_estimation.rs), [from_psp/mod.rs](../../../../src/var_calling/from_psp/mod.rs), [from_psp/columns.rs](../../../../src/var_calling/from_psp/columns.rs), [from_psp/column_span_reader.rs](../../../../src/var_calling/from_psp/column_span_reader.rs), [from_psp/loader.rs](../../../../src/var_calling/from_psp/loader.rs), [from_psp/partition.rs](../../../../src/var_calling/from_psp/partition.rs), [from_psp/driver.rs](../../../../src/var_calling/from_psp/driver.rs), [from_psp/worker.rs](../../../../src/var_calling/from_psp/worker.rs), [from_psp/test_helpers.rs](../../../../src/var_calling/from_psp/test_helpers.rs), [from_psp/kernels/mod.rs](../../../../src/var_calling/from_psp/kernels/mod.rs), [from_psp/kernels/unify_alleles.rs](../../../../src/var_calling/from_psp/kernels/unify_alleles.rs), [from_psp/kernels/compute_log_likelihoods.rs](../../../../src/var_calling/from_psp/kernels/compute_log_likelihoods.rs), [from_psp/kernels/project_scalars.rs](../../../../src/var_calling/from_psp/kernels/project_scalars.rs).
- **Deliberately out of scope:** all code outside `src/var_calling/` (`src/psp/`, `src/bam/`, `src/fasta/`, `src/vcf/`, `src/cli/`, `src/pop_var_caller/`) — read for caller/callee context only; the `benches/psp_writer_perf.rs` bench panic (psp module).
- **Categories dispatched:** reliability (always; test-gap focus on the integration seams), errors (always; `from_psp` error enums + propagation), naming (always; post-reorg vocabulary), defaults (public config + `DEFAULT_*` consts), idiomatic (always), refactor_safety (always; reorg fallout + doc gate), module_structure (multi-file scope), unsafe_concurrency (`Arc`/`Mutex`/`Condvar`/`Atomic`/`rayon` in the driver), smells (always), extras (hot path + stable VCF output + parser-like columnar decode). `tooling` skipped — module-level scope, no `Cargo.toml` ownership.

## 2. Verdict

**Request-changes.** No correctness Blockers: the concurrency review verified the parallel-worker soundness (the `Send + !Sync` REF fetcher never crosses a thread boundary; `seq_idx`-ordered emit and additive stats back the byte-identity claim; no `unsafe`, no lock-ordering hazard). But the `cargo doc` CI gate is **red** (a real, in-scope build-target failure → at least Major), the rewrite's *hard* correctness contract (byte-identical VCFs) is **unguarded by any test**, and the `from_bam` removal left a trail of **dead code, dead doc links, and an inert operator knob**. These are must-fix before the subtree is shippable, hence Request-changes rather than Approve-with-changes — even though none is a behavioral bug today.

## 3. Execution status

All commands run inside the dev container per `CLAUDE.md`. Output quoted verbatim.

- `cargo fmt --check` → **exit 0** (clean).
- `cargo clippy --all-targets --all-features -- -D warnings` → **exit 0** (clean).
- `cargo test --lib var_calling` → **exit 0**: `test result: ok. 404 passed; 0 failed; 1 ignored; 0 measured; 655 filtered out`.
- `cargo test --all-targets --all-features` → **exit 101** — the only failure is an **out-of-scope** benchmark panic: `thread 'main' panicked at benches/psp_writer_perf.rs:385:60: index out of bounds: the len is 3300000 but the index is 3300000`. No in-scope test failed. (See Out of scope observations.)
- `cargo doc --no-deps` → **exit 101**: `pop_var_caller (lib doc) generated 2 warnings` → `error: could not document pop_var_caller`. Both warnings are **in scope** (see M1).

Findings labeled "Needs verification": **0** — every cited location was read; the two doc-gate warnings and the dead-code/inert-knob claims were re-verified against source at synthesis.

## 4. Open questions and assumptions

Resolve these first; several gate the findings below.

1. **Is `chunk_genomic_span` a reserved future knob, or dead?** It is set, logged, and documented (citing a nonexistent `MAX_CHUNK_SPAN_GROWTH`), but the production loader is driven solely by `target_variants_per_chunk`. Gates **M5**: delete vs. document-as-reserved.
2. **Is `from_psp` intended as crate-public API, or implementation detail of the cohort CLI?** Nothing in `tests/`/`benches/` consumes its `pub` surface except one `_for_test` shim, and the two non-CLI consumers bypass the curated facade. Gates **M8**: `pub(crate)` sweep vs. keep-public.
3. **What is the byte-identity oracle now that the streaming `drive_cohort_pipeline` was deleted?** The contract names "the prior streaming pipeline," but that code is gone, so the only candidates are (a) a committed golden VCF, (b) serial-vs-parallel self-equivalence. Gates **M2** and **M9**.
4. **Are both `drive_blocks_serial` and `drive_blocks_parallel` permanent, or is always-parallel the intent?** Which runs is an environment property (rayon pool size), so both ship today. Gates **M2** and the stale "the serial path no longer exists" comment.

## 5. Top 3 priorities

1. **M1 — fix the red `cargo doc` gate** (two redundant explicit link targets). One-line edits; CI cannot pass until done.
2. **M2 — add a serial-vs-parallel byte-identity test.** The rewrite's hard contract is held today only by a static-reasoning argument; the test that would have caught a divergence was deliberately morphed into a weaker run-vs-run determinism check.
3. **M3 / M4 / M7 — clean up the `from_bam`-removal debris**: the dead `PerGroupMerger` error variant + `From` impl, the dead `into_shared_ref_fetcher`, and the dead `from_bam` intra-doc links that silently render as broken URLs. Convergent across many categories; cheap to fix; they actively mislead.

## 6. Findings

### Blocker

None.

### Major

**M1: [from_psp/column_span_reader.rs:3](../../../../src/var_calling/from_psp/column_span_reader.rs#L3) + [posterior_engine.rs:925](../../../../src/var_calling/posterior_engine.rs#L925) — `cargo doc` gate is red on two redundant explicit intra-doc link targets**
**Categories:** refactor_safety (+ cross-noted by naming, idiomatic, smells, module_structure, extras — convergent). **Confidence:** High.
`cargo doc --no-deps` denies warnings and exits 101 on exactly these two: `[`BlockColumnReader`](crate::psp::BlockColumnReader)` and `[`genotype_order`](crate::var_calling::per_group_merger::genotype_order)`. In each, the explicit `(path)` duplicates the target rustdoc already resolves from the link text. **Fix:** drop the explicit target — `//! Wraps a [`BlockColumnReader`] and …` / `/// The canonical [`genotype_order`] …`. These are the only two in-scope doc warnings; both still resolve to the same items.

**M2: [from_psp/driver.rs:471](../../../../src/var_calling/from_psp/driver.rs#L471), [driver.rs:536](../../../../src/var_calling/from_psp/driver.rs#L536), [driver.rs:412-432](../../../../src/var_calling/from_psp/driver.rs#L412-L432) — the byte-identical serial-vs-parallel contract is unguarded by any test**
**Categories:** reliability (filed Blocker), refactor_safety (Major), extras (Major) — convergent. **Confidence:** High. *Downgraded from the reliability sub-agent's Blocker to Major: the concurrency review verified the `seq_idx` reorder + additive `roll_window_stats` are sound by static reading, so this is a missing-regression-test gap, not a confirmed defect.*
`drive_cohort_chunked` dispatches on `rayon::current_num_threads()`: `drive_blocks_serial` (one block at a time) when ≤1, else `drive_blocks_parallel` (producer thread + bounded `BlockQueue` + N workers + a collector that re-imposes `seq_idx` order). The doc at [driver.rs:521-535](../../../../src/var_calling/from_psp/driver.rs#L521-L535) asserts the parallel path is "byte-identical to `drive_blocks_serial`," and byte-identity is the rewrite's hard contract — but `var_calling_emits_deterministic_vcf_across_runs` ([tests/cohort_cli_integration.rs:240](../../../../tests/cohort_cli_integration.rs#L240)) runs both invocations through the *same* rayon pool, so it never compares the two collectors' emit order/stats. A `seq_idx`-ordering, `pending`-drain, or stats-rollup regression would produce a wrong VCF only under multi-threading and pass every test. The integration-test comment even claims "the serial path no longer exists" — but it is live at `driver.rs:471` whenever the pool is single-threaded. **Fix:** drive the same fixture cohort through both paths (force a 1-thread and an N-thread `rayon::ThreadPoolBuilder`) and `assert_eq!` on emitted VCF bytes *and* every `ChunkDriverStats` field. Test sketch in §8.

**M3: [from_psp/driver.rs:296-323](../../../../src/var_calling/from_psp/driver.rs#L296-L323) — dead `ChunkDriverError::PerGroupMerger` variant + `From<PerGroupMergerError>` impl (from_bam debris, provenance hazard)**
**Categories:** errors, smells — convergent. **Confidence:** High (verified: `driver.rs:52` imports `PerGroupMergerError` solely to feed this variant; the only construction site is the `From` impl itself).
The variant's own doc justifies it by "the `var-calling-from-bam` shim," deleted in `d7a97d0`. No `?` site in the driver produces a bare `PerGroupMergerError` — the worker maps all per-group errors into `PosteriorEngineError`, surfaced via `ChunkDriverError::RunWindow`. Worse, `PosteriorEngineError` already owns `From<PerGroupMergerError>`, so a `PerGroupMergerError` now has two conversion routes into `ChunkDriverError`; a future inference-driven `?` could collapse the same failure into either variant — the exact provenance loss the module's own comment (lines 308-311) warns against. **Fix:** delete the variant, the `From` impl, and `PerGroupMergerError` from the `use` at line 52.

**M4: [from_psp/worker.rs:702-715](../../../../src/var_calling/from_psp/worker.rs#L702-L715) + [from_psp/mod.rs:41](../../../../src/var_calling/from_psp/mod.rs#L41) — `into_shared_ref_fetcher` is dead `pub` code with a misleading doc**
**Categories:** smells. **Confidence:** High (verified: zero call sites in `src/`/`tests/`/`benches/` — only the definition and the re-export).
The production REF path builds its `Arc` fetchers directly (`StreamingChromRefFetcher::for_contig`, `ManualEvictChromRefFetcher::for_contig`); workers consume pre-fetched byte buffers via `prefetch_window_ref_bytes` and take no `SharedRefFetcher`. The function's doc describes a `run_window` signature that no longer exists. **Fix:** delete the function and its `pub use` entry; drop the now-unused `SharedRefFetcher` import if nothing else in the non-test module uses it.

**M5: [from_psp/driver.rs:83-101](../../../../src/var_calling/from_psp/driver.rs#L83-L101) — `ChunkSizingParams::chunk_genomic_span` is an announced, logged, but inert knob; its doc cites a nonexistent const**
**Categories:** defaults (+ smells angle). **Confidence:** High (verified: read only at its decl, the startup `eprintln!`, and the two sites that set it to `DEFAULT_CHUNK_GENOMIC_SPAN`; the production `StreamingBlockLoader::fill_block` takes only `target_variants`; `MAX_CHUNK_SPAN_GROWTH` exists nowhere but this doc comment).
The startup log prints `chunk_genomic_span=100000`, telling an operator the run used a 100 kb span — but changing it changes nothing; byte-identity is held by `target_variants_per_chunk` + clean group-boundary cuts. **Fix (gated by Q1):** remove the field + const + log line + both construction sites (preferred), or rewrite the doc to mark it explicitly reserved (drop the `MAX_CHUNK_SPAN_GROWTH` and "first pull attempt" claims) and stop logging it as effective.

**M6: [from_psp/driver.rs:352-366](../../../../src/var_calling/from_psp/driver.rs#L352-L366) — startup config log omits every per-stage default**
**Categories:** defaults. **Confidence:** High.
The startup `eprintln!` logs 7 chunk-sizing/downstream-filter values but none of the behaviorally significant per-stage defaults in the same `ChunkDriverParams`: DUST window/threshold, `max_variant_group_span`, ploidy/max-alleles/batch-size, and the entire `PosteriorEngineConfig` (pseudocounts, convergence, max-iterations, max-GQ, contamination on/off). `print_run_summary` reports only counters. An operator cannot recover "what ploidy / pseudocounts / DUST threshold did this run use?" from logs — and a partial log gives false confidence it is complete. **Fix:** append the configs' `Debug` renderings (all derive `Debug`) or emit one structured event per resolved value.

**M7: [from_psp/driver.rs:23](../../../../src/var_calling/from_psp/driver.rs#L23), [driver.rs:150](../../../../src/var_calling/from_psp/driver.rs#L150), [per_group_merger.rs:25](../../../../src/var_calling/per_group_merger.rs#L25) — dangling intra-doc links to the removed `from_bam` module render as dead URLs (doc tooling does *not* catch them)**
**Categories:** refactor_safety, errors, naming, smells, module_structure, idiomatic, defaults, reliability — heavily convergent. **Confidence:** High (verified: `from_bam` is gone everywhere; these are the only surviving references).
All three use `[text](crate::var_calling::from_bam::pipeline::…)` explicit-target syntax. Because the URL part contains `::`, rustdoc treats it as a verbatim URL, emits a dead hyperlink, and runs *no* intra-doc resolution — which is why they survived the removal and don't appear in the `cargo doc` warning list. `per_group_merger.rs:25` is the worst: it points a reader at the *current* parallelism seam, which has moved to `from_psp::driver`, so it actively misleads. ([mod.rs:33](../../../../src/var_calling/mod.rs#L33) is past-tense prose, not a link — harmless.) **Fix:** repoint at the surviving driver (`crate::var_calling::from_psp::driver::{drive_cohort_chunked, ChunkDriverStats}`) or convert to prose; prefer plain `[`Foo`]` intra-doc form going forward so the doc gate catches the *next* such removal.

**M8: [from_psp/mod.rs:20-41](../../../../src/var_calling/from_psp/mod.rs#L20-L41) — the entire `from_psp` public surface could be `pub(crate)`, and its curated facade is bypassed**
**Categories:** module_structure. **Confidence:** High (verified: zero `tests/`/`benches/` consumers of the `pub` surface except the `record_fails_mapq_diff_t_for_test` shim).
Two related problems: (a) all submodules are `pub mod` *and* there is a `pub use` re-export layer, but consumers split — `var_calling.rs` uses the flat re-exports while `contamination_chunked_stream.rs`/`estimate_contamination.rs` reach in via submodule paths, and `DEFAULT_CHUNK_GENOMIC_SPAN` is imported two different ways; (b) the chunk-driver/columnar internals are exposed as crate API that nobody outside the crate uses, blocking dead-code analysis and presenting internals as stable rustdoc API. **Fix (gated by Q2):** downgrade the submodules + re-exports to `pub(crate)` and route the two non-CLI consumers through the facade; at minimum settle `DEFAULT_CHUNK_GENOMIC_SPAN` on one path.

**M9: [tests/cohort_cli_integration.rs:240-341](../../../../tests/cohort_cli_integration.rs#L240-L341) — the external byte-identity contract has no committed golden VCF (only self-consistency tests)**
**Categories:** extras. **Confidence:** High.
The two byte-identity tests compare the pipeline against *itself* (run-vs-run determinism; `target_variants_per_chunk` invariance). There is no `tests/golden/` file, no `expected*.vcf`, no snapshot. A *uniform* drift (QUAL/GQ formatting, INFO order, record ordering) would pass both and break the external contract. The named oracle ("the prior streaming pipeline") no longer exists in the tree (Q3). **Fix:** commit a small golden VCF from a fixed multi-sample fixture (SNP + MNP + compound + homref-at-variant), assert `run_var_calling`'s body (volatile `##source=`/`##commandline=` stripped, as the existing tests do) equals it byte-for-byte, regenerated only behind an explicit `BLESS=1`-style gate.

**M10: [benches/var_calling_perf.rs:812](../../../../benches/var_calling_perf.rs#L812) — no performance guard on the `from_psp` hot path**
**Categories:** extras. **Confidence:** High.
The only `var_calling` criterion group benches the *streaming* stage internals (`per_position_merger`, `variant_grouper`, `per_group_merger`, `posterior_engine`). None of `drive_cohort_chunked` / `StreamingBlockLoader::fill_block` / `ColumnSpanReader::read_span` / the three columnar kernels — the actual production hot path and the RAM-for-scaling thesis — has a benchmark, and no `// REGRESSION THRESHOLD: N%` markers exist. A regression in the columnar load/compact or kernel inner loops lands silently. **Fix:** add at least one end-to-end `bench_drive_cohort_chunked` (N ∈ {1, 8, 24}) and one `bench_compute_log_likelihoods_columnar`, each with a documented input and a regression-threshold comment.

**M11: [from_psp/driver.rs:684-748](../../../../src/var_calling/from_psp/driver.rs#L684-L748) — `BlockQueue` (hand-rolled Mutex + 2 Condvars) has zero concurrency tests**
**Categories:** reliability. **Confidence:** High. (Code verified *sound* by the concurrency review — this is a missing-test gap, but a high-value one.)
`push` blocks-while-full then `Err(block)` after `close`; `pop` drains then `None` after `close`; `close` is idempotent and wakes all. None tested — while the sibling `DustAheadPool` *is* thoroughly tested ([driver.rs:2067-2181](../../../../src/var_calling/from_psp/driver.rs#L2067-L2181)). A dropped `not_full.notify_one()` in `pop` would deadlock the producer under back-pressure with no test catching it (a lost wakeup is a hang, not a panic). **Fix:** mirror the `DustAheadPool` suite — see §8.

**M12: [from_psp/driver.rs:1253](../../../../src/var_calling/from_psp/driver.rs#L1253) — `DustAheadPool::shutdown` swallows worker panics (`let _ = handle.join()`)**
**Categories:** unsafe_concurrency. **Confidence:** Medium.
A DUST worker that panics after stashing its last result (e.g. on drop, with no further lock) is observed by no one; in the mid-run case the panic surfaces only as a secondary `"DustAheadPool mutex poisoned"` `expect` with no link to origin. **Fix:** capture/log the join result, or have `recv_next` translate a poisoned lock into the already-existing typed `ChunkDriverError::DustAheadGone` ([driver.rs:285-286](../../../../src/var_calling/from_psp/driver.rs#L285-L286)) instead of `expect`-panicking.

**M13: [from_psp/driver.rs:550-551](../../../../src/var_calling/from_psp/driver.rs#L550-L551) — unbounded `mpsc` result/recycle channels, memory-bounded only by an undocumented transitive invariant**
**Categories:** unsafe_concurrency. **Confidence:** Medium.
Both `result_tx` (owns `Vec<PosteriorRecord>`) and `recycle_tx` (owns a whole `MaterialisedChunk` × N) are unbounded. The real bound (`≤ cap + n_workers` in-flight, because the producer blocks on `BlockQueue` once `cap = n_workers*2` blocks are outstanding, and each block yields one result) is correct but written nowhere; the collector's `pending: BTreeMap` ([driver.rs:627](../../../../src/var_calling/from_psp/driver.rs#L627)) accumulates up to that bound. A future change to the cap or worker fan-out could silently defeat the back-pressure the RAM-for-scaling thesis depends on. **Fix:** document the bound at the channel-creation site, or use `sync_channel(cap + n_workers)` so the type enforces it.

### Minor

**Mi1: [per_group_merger.rs:410-413](../../../../src/var_calling/per_group_merger.rs#L410-L413), [from_psp/worker.rs:730](../../../../src/var_calling/from_psp/worker.rs#L730) — `Display` concatenates the source's message into the parent.** *errors, High.* `Upstream`/`RefFetch`/`WorkerUnifyError::RefFetch` interpolate `{0}`/`{source}` while also carrying the cause via `#[from]`/`#[source]`, so the cause double-renders under chain printing. Drop the interpolation, keep operation context (the shape `posterior_engine.rs:574` already uses).

**Mi2: [from_psp/driver.rs:104-124](../../../../src/var_calling/from_psp/driver.rs#L104-L124) — `DownstreamFilterParams`/`ChunkSizingParams` have no `Default` and no doc binding to the `DEFAULT_*` consts.** *defaults, High.* The const↔field link is enforced only on the CLI path; a library caller gets no announced default, unlike the exemplary per-stage configs. Add `/// Defaults to [`DEFAULT_MIN_QUAL_PHRED`]`-style field docs and consider a `with_project_defaults()` constructor.

**Mi3: [from_psp/driver.rs:109-124](../../../../src/var_calling/from_psp/driver.rs#L109-L124) — co-dependent `no_mapq_diff_filter: bool` + `min_mapq_diff_t: f32` with two disable channels.** *idiomatic, Medium.* When the bool is set the threshold is silently inert; `record_fails_mapq_diff_t` *also* treats a non-finite threshold as disabled ([driver.rs:807](../../../../src/var_calling/from_psp/driver.rs#L807)) — two ways to express "off." Collapse to `enum MapqDiffFilter { Off, On { min_t: f32 } }` with predicate accessors, or document the inert combination on the field.

**Mi4: [from_psp/driver.rs:68](../../../../src/var_calling/from_psp/driver.rs#L68), [driver.rs:74](../../../../src/var_calling/from_psp/driver.rs#L74) — two constants duplicated to "avoid cross-module reach."** *module_structure, Medium.* `DRIVER_PSP_BUFFER_BYTES` copies an already-`pub(crate)` const, so the stated justification no longer holds; `MAPQ_FILTER_MIN_READS_PER_SIDE = 3` duplicates a magic number that has *no* named source in `per_position_merger`. Both diverge silently. Import the buffer const directly; give `per_position_merger` a named `pub(crate) const` and import it (or add an equality assertion test).

**Mi5: `from_psp/` whole subtree — "chunk" and "block" alias the same materialised unit.** *naming, High.* `MaterialisedChunk`/`ChunkLoadScratch`/`load_chunk_from_iters`/`chunk_genomic_span`/`MIN_CHUNK_SPAN` vs. `StreamingBlockLoader`/`ReadyBlock`/`BlockIterator`/`BlockQueue`/`merge_block_ranges`/`drive_blocks_serial` — often one line apart ([loader.rs:546](../../../../src/var_calling/from_psp/loader.rs#L546); [driver.rs:1460](../../../../src/var_calling/from_psp/driver.rs#L1460) literally "The chunk/block generator"). Pick one term (MEMORY's vocab note argues "block = materialised output"), do a rename pass, and add a glossary entry to `from_psp/mod.rs`.

**Mi6: [from_psp/partition.rs:1-10](../../../../src/var_calling/from_psp/partition.rs#L1-L10), [columns.rs:717-720](../../../../src/var_calling/from_psp/columns.rs#L717-L720), [driver.rs:8-9](../../../../src/var_calling/from_psp/driver.rs#L8-L9) — stale "window" vocabulary + a dangling `chunk.windows` doc reference.** *naming, smells, High/Medium.* `driver.rs:1471` states "Genomic windows are gone: one block = one unit of work," yet `partition.rs` is built around `WindowPartition`/`partition_window`/"worker window" and its doc names a `chunk.windows` field `MaterialisedChunk` does not have; `columns.rs:717-720`'s `#[allow(clippy::single_range_in_vec_init)]` justification cites the same nonexistent field. Update the module docs to the post-reorg "one block = one partition" model and either rename the window-flavoured surface to block-flavoured or pin "window" in a glossary; fix the dangling `chunk.windows` references.

**Mi7: [from_psp/loader.rs:596-598](../../../../src/var_calling/from_psp/loader.rs#L596-L598) — `fill_block` doc says `Ok(false)`/`Ok(true)` but returns `Ok(u32)`.** *smells, High.* Stale doc from a former `bool` signature; the caller correctly treats the value as a kept-position count. Reword the "Returns" sentence to describe the count (`0` = exhausted/empty).

**Mi8: [per_group_merger.rs:217-222](../../../../src/var_calling/per_group_merger.rs#L217-L222) — `Relaxed` atomic counters lack the required ordering-justification comment.** *unsafe_concurrency, High.* `Relaxed` is in fact correct (independent monotonic counters read only after the producing chain joins), but the rule requires a comment naming the synchronization relationship, and the "read only after join" precondition is easy to break with a mid-run progress read. Add the comment (snippet in the category file).

**Mi9: [from_psp/loader.rs:716](../../../../src/var_calling/from_psp/loader.rs#L716) — `find_block_cut` (byte-identity-critical block boundary) has no direct tests.** *reliability, High.* Covered only transitively by the `streaming_matches_batch_*` equivalence tests, which can't localise an off-by-one at the `group_end_pos > watermark` comparison or exercise the empty-union / reach-equals-watermark boundaries. Add direct tests — see §8.

**Mi10: [from_psp/column_span_reader.rs:333-461](../../../../src/var_calling/from_psp/column_span_reader.rs#L333-L461), [columns.rs:707-715](../../../../src/var_calling/from_psp/columns.rs#L707-L715) — parser-like decode paths lack malformed-PSP tests.** *extras, Medium.* The decode is well-defended (typed `VarintOverflow` in `advance_pos`; loud `u32_from_usize` panic; upstream CSR validation in `psp/block.rs`), but every test feeds well-formed PSP from the project's own writer, so nothing asserts the typed-error / panic contracts actually fire. A refactor downgrading `advance_pos` to a wrapping add would pass CI. Add the two targeted tests in §8.

### Nits

- [posterior_engine/shape.rs:48](../../../../src/var_calling/posterior_engine/shape.rs#L48) — `#[allow(dead_code)]` on `GenotypeShape::n_genotypes` names *why* it's kept but no removal condition; add one or expose via a test-only accessor. *(smells)*
- [from_psp/worker.rs:574-583](../../../../src/var_calling/from_psp/worker.rs#L574-L583) — hand-rolled running-max in `columnar_passes_min_alt_obs` where `.map().max()` fits; hot path, so confirm identical lowering before changing. *(idiomatic)*
- [from_psp/columns.rs:546](../../../../src/var_calling/from_psp/columns.rs#L546) — `clone_from_columns` shadows the `Clone::clone_from` convention with different semantics; rename to `replace_with`/`clone_rows_from`. *(idiomatic)*
- [from_psp/column_span_reader.rs:168](../../../../src/var_calling/from_psp/column_span_reader.rs#L168) — `read_span_inner` returns an unnamed positional 4-tuple; borderline, only worth a struct if the block grows. *(idiomatic)*
- [from_psp/kernels/unify_alleles.rs:533](../../../../src/var_calling/from_psp/kernels/unify_alleles.rs#L533) — doc-link uses `super::super::partition::partition_window`; prefer the crate-absolute path (fragile under restructuring). *(module_structure/refactor_safety)*
- [from_psp/driver.rs:444-454](../../../../src/var_calling/from_psp/driver.rs#L444-L454) — abort-cleanup failure logs via `eprintln!` rather than structured `tracing`; observability nit (the failure is printed, not silently swallowed). *(errors)*
- 19 `..Default::default()` uses, all confined to `#[cfg(test)]` fixtures — low risk; no production `..Default::default()` in scope. *(refactor_safety)*

## 7. Out of scope observations

- **`benches/psp_writer_perf.rs:385` index-out-of-bounds panic** breaks `cargo test --all-targets --all-features` (exit 101). Unrelated to `var_calling` (psp writer bench). Follow-up: fix or `#[ignore]` it so the all-targets gate is green — separate from this review.
- **`genotype_order` / `MergedAllele` placement.** The `vcf` peer module ([src/vcf/writer.rs:30](../../../../src/vcf/writer.rs#L30), [record_encode.rs:604](../../../../src/vcf/record_encode.rs#L604)) imports `genotype_order` (defined in `per_group_merger.rs:513`), as do `posterior_engine`, `from_psp/kernels`, and four `src/vcf/` files. This cohort-wide genotype-combinatorics primitive arguably belongs in a top-level peer (`src/genotype_order.rs`) rather than inside a pipeline stage. Most consumers are out of scope (`src/vcf/`); raise as a cross-module placement decision in a separate PR.

## 8. Missing tests to add now

Grouped by function under test. The reliability challenge-tests pass feeds this directly; full bodies/sketches are in `tmp/review_2026-06-01_var_calling/{reliability,extras}.md`.

**`drive_blocks_serial` / `drive_blocks_parallel` (driver.rs:471/536) — M2**
- `drive_blocks_serial_and_parallel_emit_identical_records_and_stats` — input: a multi-block in-memory PSP cohort dense enough to exercise `seq_idx` ordering, run through a 1-thread and an N-thread (`ThreadPoolBuilder::new().num_threads(4)`) pool. Catches: a `seq_idx`-reorder, `pending`-drain, or stats-rollup regression that diverges only under parallelism. Assert byte-equal VCF bodies *and* every `ChunkDriverStats` field.

**`BlockQueue` (driver.rs:684-748) — M11** (mirror the `DustAheadPool` suite)
- `block_queue_push_errs_after_close` — `close()` then `push` returns `Err(block)`.
- `block_queue_pop_drains_then_returns_none_after_close` — push 2, close, pop/pop = `Some`, pop = `None`.
- `block_queue_push_blocks_until_pop_frees_a_slot` — cap=1, fill, spawn a blocked pusher, probe it's still blocked, pop once, join. Catches a lost `not_full` wakeup (hang).
- `block_queue_close_wakes_all_waiters` — spawn a blocked `push` + a blocked `pop`, `close()`, both return.

**`find_block_cut` (loader.rs:716) — Mi9**
- `find_block_cut_returns_watermark_plus_one_for_empty_union` — `assert_eq!(find_block_cut(&[], &[], 100), 101)`.
- `find_block_cut_closes_group_with_reach_equal_to_watermark` — `union=[100], ref_span=[1]`, watermark 100 → 101.
- `find_block_cut_returns_open_group_start_when_reach_exceeds_watermark` — `union=[100,105], ref_span=[10,1]`, watermark 102 → 100. Catches an off-by-one at `group_end_pos > watermark` that splits a group across blocks.

**`prefetch_window_ref_bytes` (worker.rs:296) — reliability**
- `prefetch_window_ref_bytes_reuses_buffer_across_changing_n_groups` — call into the *same* `out` with 3-, then 1-, then 2-group partitions; assert `out.len() == n_groups` and per-group bytes are correct each time. Catches a `truncate`→`clear` regression feeding a stale REF buffer into layer 1.

**DUST-plan construction (driver.rs:920/1004/1110) — reliability**
- `build_dust_plans_skips_empty_and_uncovered_contigs` and `flatten_to_jobs_is_chrom_major_interval_order` — catch a producer↔pool lockstep desync that applies the wrong DUST mask (release: wrong mask; debug: `debug_assert_eq!` trip).

**`advance_pos` (column_span_reader.rs:234) — Mi10**
- `advance_pos_returns_column_decode_error_on_u32_overflow` — `advance_pos(0, u64::from(u32::MAX) + 1, 3)` returns `Err(PspReadError::ColumnElementDecode { .. })`, not a silent wrap.

**Columnar decode malformed-input (Mi10)**
- `read_span` over a hand-crafted block whose delta-pos varint exceeds `u32::MAX` → asserts `PspReadError::ColumnElementDecode { column: "delta-pos", source: VarintOverflow, .. }`.
- `#[should_panic(expected = "…")]` driving `u32_from_usize` past `u32::MAX`.

**External byte-identity (M9)**
- Committed `tests/golden/` VCF + a test asserting `run_var_calling`'s body (volatile lines stripped) equals it byte-for-byte, regenerated only behind a `BLESS=1` gate.

**Hot-path benches (M10)**
- `bench_drive_cohort_chunked` (N ∈ {1, 8, 24}) and `bench_compute_log_likelihoods_columnar`, each with a `// REGRESSION THRESHOLD: N%` comment.

## 9. What's good

- **Parallel-worker soundness is correct by construction.** The `Send + !Sync` `StreamingChromRefFetcher` ([fasta/fetcher.rs:110-113](../../../../src/fasta/fetcher.rs#L110-L113) `assert_send` fence) stays on the producer thread; workers receive only pre-fetched `Vec<u8>` REF bytes, so the `!Sync` fetcher never crosses a boundary — the central post-reorg claim holds.
- **Loader data-parallelism is race-free by disjointness:** every `par_iter_mut().zip(...)` in [loader.rs](../../../../src/var_calling/from_psp/loader.rs) operates over disjoint per-sample buffers with shared *immutable* `position_union`/`is_kept` — deterministic regardless of thread count.
- **Lock discipline:** `BlockQueue` and `DustAheadPool` each hold a single lock at a time, never nested, with loop-guarded condvar waits and `notify_all` on shutdown — no lock-ordering or missed-wakeup hazard.
- **The per-stage config structs are an exemplar** ([posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs), [per_group_merger.rs](../../../../src/var_calling/per_group_merger.rs), [dust_filter.rs](../../../../src/var_calling/dust_filter.rs)): named `pub const DEFAULT_*` with unit/source docs, `# Defaults` tables, and validated constructors — the bar the `from_psp` driver knobs (Mi2) should meet.
- **The four manual `Debug` impls** ([per_position_merger.rs:165](../../../../src/var_calling/per_position_merger.rs#L165), [per_group_merger.rs:621](../../../../src/var_calling/per_group_merger.rs#L621), [posterior_engine.rs:1196](../../../../src/var_calling/posterior_engine.rs#L1196), [variant_grouping.rs:173](../../../../src/var_calling/variant_grouping.rs#L173)) use exhaustive named-field destructures with comments — so a new field forces a compile-time decision.

## 10. Commands to re-verify

Re-run (inside the dev container per `CLAUDE.md`):

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib var_calling`
- `./scripts/dev.sh cargo doc --no-deps` — **currently red (M1)**; should pass once the two redundant-link targets are dropped.

New invocations this review would add:
- the serial-vs-parallel and `BlockQueue`/`find_block_cut`/decode tests (§8): `cargo test --lib var_calling::from_psp`.
- the golden-VCF integration test: `cargo test --test cohort_cli_integration`.
- the new `from_psp` benches: `cargo bench --bench var_calling_perf`.

### Author response convention

Address each finding by identifier (`M1`, `Mi3`, …) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the §4 open questions first — Q1 gates M5, Q2 gates M8, Q3 gates M2/M9, Q4 gates M2.
