# Code Review: re-architecture-pipeline
**Date:** 2026-06-05
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** The re-architected cohort `.psp` → VCF pipeline (record-streaming producer → W callers → writer), swapped into production on branch `re-architect`
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** Diff/PR — the re-architecture of the cohort `var-calling` pipeline. The old columnar pipeline (`driver.rs` / `worker.rs` / `loader.rs` / `columns.rs` / `partition.rs` / `two_pass.rs` / `kernels/`, ~14 k LoC) was deleted and replaced by a flat record-streaming topology, then swapped into production (Phase 7).
- **Reviewed against:** branch `re-architect` @ `ef93b67`, vs `main` (`git diff main...re-architect`).
- **In-scope files:**
  - New: [types.rs](../../../../src/var_calling/types.rs), [sample_reader.rs](../../../../src/var_calling/sample_reader.rs), [cohort_integration.rs](../../../../src/var_calling/cohort_integration.rs), [pileup_overlaps.rs](../../../../src/var_calling/pileup_overlaps.rs), [em_posterior_calc.rs](../../../../src/var_calling/em_posterior_calc.rs), [vcf_writer.rs](../../../../src/var_calling/vcf_writer.rs), [pipeline.rs](../../../../src/var_calling/pipeline.rs), [mod.rs](../../../../src/var_calling/mod.rs).
  - Modified production: [per_group_merger.rs](../../../../src/var_calling/per_group_merger.rs) (the `merge_group_with_ref` / `MergeGroupOutcome` extraction), [pop_var_caller/var_calling.rs](../../../../src/pop_var_caller/var_calling.rs) (CLI rewiring), [pop_var_caller/estimate_contamination.rs](../../../../src/pop_var_caller/estimate_contamination.rs) (contamination port), [pop_var_caller/mod.rs](../../../../src/pop_var_caller/mod.rs), [Cargo.toml](../../../../Cargo.toml).
- **Deliberately out of scope:** the copied-verbatim numeric kernels (`posterior_engine`, `variant_grouping`, `per_position_merger`, `dust_filter`, `contamination_estimation`) — byte-for-byte unchanged from `main`, already reviewed; deleted files (referenced only as the equivalence oracle); markdown design docs; `test_helpers.rs`; `benches/psp_writer_perf.rs` (pre-existing, see §7).
- **Categories dispatched:** all 11 — `reliability` (always), `errors` (always), `naming` (always; clarity is the project's stated goal), `defaults` (CLI/config surface), `idiomatic` (always), `refactor_safety` (large behavior-preserving rewrite with a byte-identity contract), `module_structure` (multi-file restructure), `unsafe_concurrency` (crossbeam channels + `thread::scope` + rayon `par_iter_mut`), `smells` (always), `tooling` (crate with `Cargo.toml`; new dep), `extras` (`.psp` parser / hot path / stable VCF output / PR-intent).

## 2. Verdict

**Request-changes.**

The structural rewrite is well-executed and the byte-identity-sensitive math checks out line-by-line against `main` (see §9). There are **no correctness Blockers in the emitted calls**: `fmt`/`clippy` are clean, lib + integration tests pass, and `unsafe_concurrency` returned a clean topology (no `unsafe`, sound channel/close ordering, order-independent parallel fold). The verdict is driven by **a red `cargo doc` CI gate** (M1, 7 unresolved-link errors), **a silent VCF-provenance regression** (M2, `##commandline` hard-coded), **a malformed-input panic regression** (M3), and the fact that the rewrite **deleted its own byte-identity oracle** while leaving the parallel-ordering guarantor (the writer reorder buffer) and several release-only invariants untested (M6, M7).

## 3. Execution status

All commands run in the project dev container (Apple `container`, `scripts/dev.sh`).

| Command | Exit | Result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib` | 0 | **996 passed**, 0 failed, 1 ignored |
| `cargo test --test '*'` | 0 | all integration green (cohort_cli 12, cohort_vcf_writer 6, contamination 5, pileup_cli 9, posterior 7, psp_to_pileup 4, …) |
| `cargo test --all-targets --all-features` | 101 | **red, but only** from the pre-existing out-of-scope `benches/psp_writer_perf.rs:386` panic (`index out of bounds: the len is 3300000 but the index is 3300000`). CI's test step is `cargo test --lib --tests` (excludes benches), so the CI test gate is **green**. See §7. |
| `cargo doc --no-deps` | 101 | **RED — in-scope.** 7 `unresolved link` errors + 5 `redundant explicit link target` warnings → `could not document pop_var_caller`. See M1. |
| `cargo audit` | — | not run (no advisory-DB change relevant to this diff; the only new dep is `crossbeam-channel = "0.5"`). |

Findings labeled "Needs verification": 0 (one Major, M3, is filed at Medium confidence with an explicit verification step).

## 4. Open questions and assumptions

1. **Is the `##commandline` header change intentional?** The pipeline hard-codes `command_line: "var-calling"` (was `current_command_line()`). Gates **M2**. Assumed *unintentional* (the import was removed as part of the rewiring, not by design).
2. **Should `--target-variants-per-chunk` remain a silent no-op for `estimate-contamination`?** The stream port stopped reading it. Gates **Mi2**.
3. **Now that `main`'s columnar path is deleted, where does the byte-identity oracle live?** The in-repo A/B test it provided is gone; equivalence currently rests on self-consistency tests. Gates **M7** and the §8 missing tests. (Was equivalence captured as an out-of-tree golden/diff, per the var_calling 2026-06-01 review's M9 decision?)
4. **Is the crate-wide `#![allow(dead_code)]` meant to persist post-P7?** Its own comment says "removed at the P7 swap". Gates **M5**.

## 5. Top 3 priorities

1. **M1 — Red `cargo doc` gate.** 7 unresolved intra-doc links (dead references to the deleted `driver` / `drive_cohort_chunked` / `contamination_chunked_stream`, plus broken `GrouperError` / `CallStats` / `PosteriorEngineConfig::contamination`) block merge under CI's `RUSTDOCFLAGS=-D warnings`. Mechanical, but must be fixed and all 5 warnings cleared too.
2. **M2 — `##commandline` provenance regression.** Every emitted VCF now records `var-calling` instead of the real argv. One-line fix; restore `current_command_line()`.
3. **M3 — Zero-allele `.psp` record panics / mis-decodes `from_block`.** A crafted/corrupt `.psp` that the row reader handles cleanly makes the new producer panic (trailing record) or silently emit a wrong `ref_span` (interior record). Robustness regression vs the path it mirrors.

## 6. Findings

### Blocker

*(none)*

### Major

**M1: [src/var_calling/types.rs:111](../../../../src/var_calling/types.rs#L111) (+ 6 more sites) — Red `cargo doc` CI gate: 7 unresolved intra-doc links + 5 redundant-explicit-link warnings**
**Categories:** tooling, module_structure, smells, extras (convergent)
**Confidence:** High
`cargo doc --no-deps` exits 101. CI sets `RUSTDOCFLAGS: -D warnings` on the `cargo doc --no-deps --lib --all-features` step (`.github/workflows/ci.yml:45-46`), so this blocks merge. The full set (only partially visible in any single sub-agent's run — the author must run `cargo doc` to see all of them):

- **7 `error: unresolved link`** — dead references left by the swap:
  - `src/pop_var_caller/var_calling.rs:11` → `GrouperError`; `:12` → `crate::var_calling::driver::drive_cohort_chunked`; `:19` → `PosteriorEngineConfig::contamination`; `:223` → `drive_cohort_chunked`
  - `src/var_calling/contamination_estimation.rs:633` → `crate::pop_var_caller::contamination_chunked_stream::ContaminationStreamError`
  - `src/var_calling/per_group_merger.rs:25` → `crate::var_calling::driver`
  - `src/var_calling/vcf_writer.rs:36` → `CallStats` (needs the full path / in-scope import)
- **5 `warning: redundant explicit link target`** — at `src/var_calling/types.rs:111`, `src/var_calling/cohort_integration.rs:30/38/66`, `src/pop_var_caller/estimate_contamination.rs:353`.

**Fix:** repoint or delete each dead link (the `driver` / `drive_cohort_chunked` / `contamination_chunked_stream` targets no longer exist; `GrouperError` / `CallStats` / `PosteriorEngineConfig::contamination` need a correct in-scope path), and drop the redundant explicit targets (e.g. ``[`PosteriorRecord`](crate::…::PosteriorRecord)`` → ``[`PosteriorRecord`]``). Then re-run `cargo doc --no-deps` to confirm green. This is the same class of gate the var_calling 2026-06-01 review flagged (M1) — worth a crate-wide `cargo doc` pass before the next merge.

**M2: [src/var_calling/pipeline.rs:179](../../../../src/var_calling/pipeline.rs#L179) — Hard-coded `command_line: "var-calling"` silently drops VCF provenance**
**Categories:** refactor_safety, defaults, extras, reliability (convergent)
**Confidence:** High
The writer's `CohortMetadata` is built with `command_line: "var-calling".to_string()`. On `main` this was `current_command_line()` (the space-joined `args_os()`); the CLI diff removed exactly that import and assignment. The field is written into the VCF header as `##commandline=…` ([src/vcf/header.rs:75](../../../../src/vcf/header.rs#L75)), so every emitted VCF now records `var-calling` regardless of `--reference`, `--regions`, `--threads`, output path, ploidy, etc.

**This is *not* a byte-identity-contract violation** — the project's oracle strips `##source=` / `##commandline=` before diffing (confirmed at `tests/cohort_cli_integration.rs:380,436`), and the contract is "drop `^##`, md5". But it is a real, user-visible **provenance regression**: "how was this VCF produced?" is no longer answerable from the file, and the fallback is invisible to and un-overridable by the caller. It fails the PR's stated intent ("byte-identical output vs `main`").
**Fix:** restore the real command line — `command_line: crate::pop_var_caller::common::current_command_line()` (it is `pub(crate)`, already reachable from the pipeline module), or thread it in from the CLI wrapper so all provenance fields are filled at one layer.

**M3: [src/var_calling/sample_reader.rs:339-344](../../../../src/var_calling/sample_reader.rs#L339-L344) — Zero-allele `.psp` record panics / mis-decodes `from_block`; the row reader is robust**
**Category:** extras
**Confidence:** Medium (with a stated verification step)
The block decoder validates `sum(n_alleles) == n_total_alleles` and `n_total_alleles >= n_records`, but **not** per-record `n_alleles[r] >= 1` — only the *writer* rejects zero-allele records (`InvalidRecordKind::ZeroAlleles`). So a crafted/corrupt `.psp` with `n_alleles = [0, 2]` survives decode. For each kept record `from_block` unconditionally reads `cols.allele_seq_offsets[lo + 1] - cols.allele_seq_offsets[lo]` as the REF span (line 344):
- a **trailing** zero-allele record (`lo == n_total_alleles`) indexes `allele_seq_offsets[lo + 1]`, one past the `n_total_alleles + 1`-length array → **index-out-of-bounds panic**;
- an **interior** zero-allele record silently reads the *next* record's offset → a wrong `ref_span` that perturbs grouping and REF-fetch.

The row decoder (`reader.rs:706-760`) is robust here (its per-allele loop is simply empty), so this diff *regresses* untrusted-input robustness vs the path it mirrors.
**Fix:** preferably add a per-record `n_alleles[r] >= 1` invariant in `decode_block_payload` returning a typed error (fixes both readers, matches the writer's rule); or guard `from_block`:
```rust
let ref_span = if hi > lo { cols.allele_seq_offsets[lo + 1] - cols.allele_seq_offsets[lo] } else { 0 };
```
Add a malformed-`.psp` regression test (craft the block bytes — the writer can't emit a zero-allele record) asserting a typed `PspReadError`, not a panic. **Verification step (Medium → High):** confirm no layer between `decode_block_payload` and `from_block` rejects per-record zero alleles (the decoder and `block.rs` invariants do not).

**M4: [src/var_calling/cohort_integration.rs:407-409](../../../../src/var_calling/cohort_integration.rs#L407-L409) — REF-fetch error context flattened into `ProducerError::Ref(String)`**
**Category:** errors
**Confidence:** High
The producer's REF-fetch failure is a `String`. The closure that feeds it (`pipeline.rs::ref_fetch`, lines 312-328) stringifies a typed `ChromRefFetchError` *twice* (`.map_err(|e| e.to_string())` on both `for_contig` and `fetch`), collapsing the root cause (`OutOfPattern { … }`, `OutOfBounds { … }`, an `io::Error`) into a flat message with no `source()` chain — even though `PipelineError` already has a typed `RefFetch(#[from] ChromRefFetchError)` variant (pipeline.rs:65) that this path bypasses. The producer's most diagnostically important failure mode (a non-monotonic streaming-fetch bug, a contig-length mismatch) loses its typed cause exactly where it is needed.
**Fix:** carry the typed error — make `ProducerError::Ref(#[source] Box<dyn std::error::Error + Send + Sync>)` and have `ref_fetch` return the boxed typed error instead of `e.to_string()`; or thread `ChromRefFetchError` through the closure's error type.

**M5: [src/var_calling/mod.rs:42](../../../../src/var_calling/mod.rs#L42) — Crate-subtree `#![allow(dead_code)]` persists post-P7, masking dead public surface**
**Categories:** module_structure, smells, tooling, idiomatic, reliability (convergent)
**Confidence:** Medium-High
The blanket `#![allow(dead_code)]` is justified by a comment as "Build-phase only … Removed at the P7 swap" — but P7 has shipped (`1d34f85`) and the package is now production. The allow now weakens the dead-code lint across the whole `var_calling` subtree on live code. `module_structure` identified verified-dead `pub` items it masks (e.g. `SamplePspReader::into_reader`, `CohortSpanFold::n_positions`, `CohortSpanFold::chunk_cuts`). (Note: a sub-agent that rebuilt `--all-targets` after removing the allow saw 0 warnings, because the test-only `take_*`/`RefSpan::empty`/`merge` items are used under `#[cfg(test)]`; run a plain `cargo build` to surface the lib-dead items.)
**Fix:** remove the blanket allow; prune the genuinely-dead items or tighten their visibility to `pub(crate)`; pinpoint-`#[allow(dead_code)]` only the few intentionally-retained verbatim-kernel items.

**M6: [src/var_calling/cohort_integration.rs:752](../../../../src/var_calling/cohort_integration.rs#L752) & [src/var_calling/vcf_writer.rs:112](../../../../src/var_calling/vcf_writer.rs#L112) — Release-load-bearing invariants guarded only by `debug_assert!`**
**Categories:** errors, reliability, unsafe_concurrency (convergent)
**Confidence:** Medium
Two correctness backstops are compiled out in release:
- `produce_chunk` asserts `cut > self.next_chunk_start` only via `debug_assert!`. The `variant.is_empty()` branch does `self.next_chunk_start = cut; continue;`, so if `find_cut` ever returned `cut == next_chunk_start` in release (a degenerate fold at the chunk boundary), the loop **spins forever** — a release-only hang, worse than a panic, on the critical path. I did not find a reachable counterexample, but the guard that would catch one is gone.
- `VcfWriter::finish` enforces the gapless `chunk_order` invariant (one `CalledChunk` per chunk, the data-completeness backstop) via `debug_assert!` only; in release a lost chunk silently truncates the VCF.

**Fix:** promote both to release-level guards that return a typed error — a new `ProducerError::StalledCut { next_chunk_start, cut }` and a `WriterError::MissingChunks { count }` — so the failure is diagnosable rather than a hang / silent truncation.

**M7: [src/var_calling/vcf_writer.rs:94](../../../../src/var_calling/vcf_writer.rs#L94) — The parallel-ordering guarantor and the verbatim filter chain are untested at the unit level**
**Category:** reliability
**Confidence:** High
The rewrite deleted `main`'s columnar path — the live A/B oracle — so byte-identity now rests on self-consistency tests that cannot catch a bug the deleted and new paths don't share. The highest-value gaps:
- `VcfWriter` has **no `#[cfg(test)] mod tests`**. Its `handle` → `BTreeMap` reorder → `next_expected` drain is the *sole* guarantor of genomic VCF order under W>1 workers, yet the full-pipeline tests feed chunks already in `chunk_order` (the producer emits them ordered), so the out-of-order insert-ahead/drain-later branch is **never exercised**. A drain bug (wrong key, off-by-one, failure to flush a now-contiguous backlog) passes every test yet reorders/drops records whenever two workers finish out of order.
- `emit_or_drop`'s filter *order* (hom-ref → QUAL → MAPQ-diff t) and the `WriterStats` counter wiring are hand-reassembled here and asserted by no unit test.
- `read_samples`' rayon `par_iter_mut` decode is never run under an actual multi-thread pool against the one-shot reference, nor with a corrupt sample to check error propagation.

**Fix:** add the unit tests specified in §8 (permuted-`chunk_order` reorder, buffered-future-chunk, per-gate `emit_or_drop`, multi-thread `read_samples` vs reference, decode-error surfacing). For the deleted oracle, record (in the impl report) the out-of-tree commit/run that confirmed byte-identity vs `main`, per the var_calling 2026-06-01 M9 decision.

**M8: [src/var_calling/pipeline.rs:48,145-149,199-207](../../../../src/var_calling/pipeline.rs#L145-L149) — Behavioral defaults hidden behind sentinels/magic numbers, undocumented and unlogged**
**Category:** defaults
**Confidence:** High
Three compounding issues:
- `--target-variants-per-chunk` is `default_value_t = 0`; the pipeline maps `0 → DEFAULT_TARGET_VARIANTS (1024)`. The const is private to `pipeline.rs`, the CLI literal `0` doesn't reference it, `main`'s named `DEFAULT_DESIRED_VARIANTS_PER_BLOCK = 1024` was deleted (so the two `1024`s can drift behind a comment), and the **CLI arg doc still describes deleted behavior** ("keep the legacy BP-only loop … `--chunk-genomic-span`"). The effective value appears in no startup log or run summary.
- `n_workers` (`available_parallelism().unwrap_or(1)`) and `cap = (2 * n_workers).max(1)` — concurrency width and queue back-pressure depth (a peak-RSS knob) — default implicitly; the `2×` is a bare magic number with no named const, and neither resolved value is logged or summarised. The `--threads` doc doesn't mention it now also sizes the channel topology.

**Fix:** make the CLI default reference a shared `DEFAULT_TARGET_VARIANTS` const (drop the `0` sentinel, or document the mapping); extract `const QUEUE_DEPTH_PER_WORKER: usize = 2` with a rationale doc; rewrite the stale arg doc; emit the resolved `target_variants` / `n_workers` / `cap` in a startup `tracing` event or in `print_run_summary`.

### Minor

- **Mi1** — [pop_var_caller/var_calling.rs:409](../../../../src/pop_var_caller/var_calling.rs#L409): the `chunks_loaded` / `avg_variants_per_chunk` stderr run-summary line was dropped with the old driver (`WriterStats` has no chunk count). Stderr-only, no VCF impact; re-thread a counter or document the intentional removal. *(Categories: refactor_safety, defaults.)*
- **Mi2** — [pop_var_caller/estimate_contamination.rs:92](../../../../src/pop_var_caller/estimate_contamination.rs#L92): `--target-variants-per-chunk` is now an accepted no-op for `estimate-contamination` (the `PerPositionMerger` port ignores it; it never affected output on `main` either). Drop the flag or note in `--help` that it no longer affects the side-pass.
- **Mi3** — contamination port: the "byte-identical to the chunk-loader stream" claim has no in-repo equivalence test (the deleted stream is gone; the surviving test only proves invariance to the now-inert `target_variants_per_chunk`). `refactor_safety` verified equivalence by `git show` diff and found the memory consequence *favorable* (a streaming k-way merge holds ≈N peeked heads, not a chunk × N buffer) — so this is a coverage gap, not a defect. Add a golden-artefact regression or document the external verification.
- **Mi4** — [per_group_merger.rs:803](../../../../src/var_calling/per_group_merger.rs#L803): error-precedence change — `main`'s `process_group` checked `ploidy == 0` before `end < start`; the extraction moved `ploidy == 0` into `merge_group_with_ref` (checked *after* the inverted-range guard). A pathological `ploidy=0` + inverted-range fixture now surfaces a different typed error. Both are error-path-only and unreachable on validated input; pin the precedence with a one-line test or restore the order.
- **Mi5** — [pipeline.rs:296,300](../../../../src/var_calling/pipeline.rs#L296-L300): `handle.join().expect("caller thread panicked")` / `expect("writer thread panicked")` re-raise a worker panic (correct, and `unsafe_concurrency` confirmed it is not swallowed), but lack the mechanical `// PANIC-FREE:` annotation the errors checklist requires. Add a one-line comment, or `std::panic::resume_unwind(handle.join().unwrap_err())` to forward the original payload.
- **Mi6** — [cohort_integration.rs:859-862](../../../../src/var_calling/cohort_integration.rs#L859-L862): `binary_search(...).expect("variable position must be in the fold")` in `fetch_ref_span` relies on a cross-method data-flow invariant (`variable ⊆ positions()`), not a type guarantee. Add a `// PANIC-FREE:` comment naming it, or carry the fold index alongside `variable` and drop the search.
- **Mi7** — [cohort_integration.rs:405](../../../../src/var_calling/cohort_integration.rs#L405): `ProducerError::Merge(Box<PerPositionMergerError>)` interpolates the inner error via `{0}` rather than marking it `#[source]`, so the cause isn't in the `source()` chain. Use `Merge(#[source] Box<…>)` and drop `: {0}`.
- **Mi8** — stale module docs contradict the shipped code: [pipeline.rs:10-18](../../../../src/var_calling/pipeline.rs#L10-L18) ("**Phase 4: single-threaded** … the writer's reorder is a no-op here") describes a superseded phase — the body *is* the crossbeam parallel topology, and the byte-identity argument now depends on the reorder actually running; [sample_reader.rs:13-14](../../../../src/var_calling/sample_reader.rs#L13-L14) ("Used single-threaded … `!Send`-by-ownership") contradicts the `par_iter_mut` parallel decode that requires `R: Send`. *(Categories: extras, reliability, module_structure.)*
- **Mi9** — [pipeline.rs:45-51](../../../../src/var_calling/pipeline.rs#L45-L51): `BUFFERED_IO_CAPACITY` / `DUST_SUBSPAN` / `DEFAULT_TARGET_VARIANTS` re-type `main`'s values behind "matches X" comments instead of `use`-ing the source const; `DEFAULT_BUFFERED_IO_CAPACITY` is `pub(crate)` and importable, and the `DEFAULT_DESIRED_VARIANTS_PER_BLOCK` comment now points at a deleted const.
- **Mi10** — [cohort_integration.rs:516](../../../../src/var_calling/cohort_integration.rs#L516): `CohortChunkIntegrator::new` clamps `target_variants.max(1)` (`0 ⇒ 1`) while the pipeline maps `0 ⇒ 1024`, so the same field has two divergent zero conventions. Take an already-resolved `NonZeroU32` so the `0 ⇒ ?` decision lives at one boundary (as `main`'s `ChunkSizingParams` did).
- **Mi11** — module/type name: [em_posterior_calc.rs](../../../../src/var_calling/em_posterior_calc.rs) is named after one sub-step (EM/posterior) but owns the whole caller section (group → merge → `min_alt_obs` filter → EM, the `VariantCaller`). Every sibling names its concept (`sample_reader`, `cohort_integration`, `pileup_overlaps`, `vcf_writer`); consider `variant_caller`. *(Categories: naming, module_structure.)*
- **Mi12** — [types.rs](../../../../src/var_calling/types.rs): `CohortPileupRecord` vs `PileupCohortChunk` — a word-order inversion on the two hottest data-path types, distinguished only by the trailing word. Consider `CohortPileupChunk` so the container shares the element's stem.
- **Mi13** — [vcf_writer.rs:26-31](../../../../src/var_calling/vcf_writer.rs#L26-L31): `DownstreamFilters` pairs `no_mapq_diff_filter: bool` + `min_mapq_diff_t: f32` — the threshold is dead when the filter is off. A `MapqDiffFilter` enum (off / on(threshold)) with a predicate method makes the dependency type-enforced (byte-identity-neutral).
- **Mi14** — [pipeline.rs:29](../../../../src/var_calling/pipeline.rs#L29): the pipeline imports `VarCallingArgs` from the `pop_var_caller` CLI stage — a cross-stage back-reference the code itself flagged "temporary until P7", now permanent. Consider a stage-local config struct the CLI builds.
- **Mi15** — [Cargo.toml:64](../../../../Cargo.toml#L64): the `crossbeam-channel` comment points at `src/var_calling_new/pipeline.rs`, a path renamed away in the swap.
- **Mi16** — [cohort_integration.rs:864](../../../../src/var_calling/cohort_integration.rs#L864): `let len = max_reach - first + 1;` is unguarded `u32` arithmetic on position-derived values (cannot underflow today — `max_reach >= first` by construction), while the neighbouring `reach`/`find_cut`/`merge_block_ranges` deliberately use `saturating_*`. Use `saturating_sub`/`saturating_add` for consistency, or assert the invariant.

### Nits

Grouped: `RefSpan::slice` debug-asserts the lower bound but lets the upper bound panic via raw index (add a symmetric `debug_assert!(hi <= self.bytes.len())`); `CohortChunkIntegrator.n` field duplicates `readers.len()` and uses a bare abbreviation against the crate's `n_samples` convention; watermark threaded as bare `w` and sample index as bare `s` where the crate uses `watermark`/`sample_idx`; the `take_*` getters' "call-once" contract is enforced only by emptying the column (a second call silently returns nothing); the six `range.clone()`s in `PerAlleleFixed::extend_from_range`; `VariantCaller::new` is an undocumented `pub fn`; a test helper `fn fold` shadows the domain concept.

## 7. Out of scope observations

- **`benches/psp_writer_perf.rs:386`** panics (`index out of bounds: the len is 3300000 but the index is 3300000`), making local `cargo test --all-targets` red. This is **pre-existing** (PROJECT_STATUS records it as the "out-of-scope `psp_writer_perf` bench panic") and does not affect the CI test gate (`cargo test --lib --tests` excludes benches). Suggested follow-up: separate issue to fix the off-by-one in the bench harness.
- **Hot-path allocations on the producer's serial critical path** (perf is an explicit non-goal here, so informational): `build_records` ([cohort_integration.rs:818-824](../../../../src/var_calling/cohort_integration.rs#L818-L824)) allocates a fresh `Vec<bool>` keep-mask per buffered chunk per sample and does a `variable.binary_search(&p)` per record; `fetch_ref_span` re-`binary_search`es each `variable` position back into `fold.positions()` that `derive_is_kept` already walked. Both `variable` and the chunk positions are sorted, so a single merge-walk would drop the searches and the mask `Vec`. This matches the known "producer record-building bottleneck" memory note — worth folding into that deferred perf lever, not this review.

## 8. Missing tests to add now

Grouped by function under test (specs from the `reliability` challenge pass):

- **`VcfWriter::handle` / reorder** — `vcf_writer_emits_records_in_genomic_order_when_chunks_arrive_permuted`: feed `CalledChunk`s with `chunk_order` permuted (2,0,1,4,3), each with records at known ascending positions; `finish`; assert emitted positions strictly ascending and count == sum of inputs. `vcf_writer_handle_buffers_future_chunk_until_predecessors_arrive`: after `handle(2)` assert `records_written == 0`; after `handle(0)`,`handle(1)` assert all flushed. Catches a drain bug invisible to the in-order pipeline tests.
- **`emit_or_drop`** — `emit_or_drop_drops_hom_ref_before_qual_and_counts_correctly`: hom-ref `Variant` with `qual_phred` below threshold → `records_dropped_hom_ref == 1`, `records_dropped_low_qual == 0`, nothing written. Catches a filter-order swap. Plus `emit_or_drop_counts_unconverged_but_still_writes`.
- **`read_samples`** — `read_samples_parallel_pool_matches_one_shot_reference`: run `streaming(...)` inside `rayon::ThreadPoolBuilder::new().num_threads(8).build().unwrap().install(...)`; assert `== reference(...)`. `produce_chunk_surfaces_decode_error_from_one_sample`: one corrupt `.psp` buffer → `Err(ProducerError::Decode(_))`, no hang.
- **`from_block`** — `from_block_returns_none_when_all_records_below_floor`; `from_block_clamps_inclusive_region_end_exactly` (region_end == a record pos keeps it; one less drops it); `from_block_single_allele_record_has_zero_nonref_obs_and_ref_span_eq_seq_len`; and the **M3** malformed case `from_block_returns_typed_error_on_zero_allele_record`.
- **`produce_chunk` / `fill_to_target`** — `produce_chunk_makes_progress_on_degenerate_zero_span_group`: a fold where reaches collapse at the boundary; assert termination within a bounded iteration count (guards M6).
- **`merge_group_with_ref`** — `merge_group_with_ref_rejects_inverted_range_and_zero_ploidy_consistently`: `OverlappingVariantGroup{start:50,end:10}` → `Err`; `ploidy:0` config → `Err`; pins Mi4 and confirms no `end - start + 1` underflow panic.
- **`run_var_calling` (pipeline)** — `run_var_calling_surfaces_caller_error_without_deadlock`: a fixture forcing a `CallerError` (group exceeding `var_group_max_span`); assert `Err(PipelineError::Caller(_))` within a timeout (first-error-wins + no channel deadlock).

## 9. What's good

- **Disciplined "rebuild the skeleton, transplant the organs."** `refactor_safety` verified line-by-line against the `main` originals that `reach` / `derive_is_kept` / `merge` / `find_cut` / `count_kept_below` / `chunk_cuts` / `merge_block_ranges` ([cohort_integration.rs](../../../../src/var_calling/cohort_integration.rs)), `emit_or_drop` + `record_fails_mapq_diff_t` ([vcf_writer.rs](../../../../src/var_calling/vcf_writer.rs)), `passes_min_alt_obs` ([em_posterior_calc.rs:163](../../../../src/var_calling/em_posterior_calc.rs#L163)), and the `merge_group_with_ref` extraction ([per_group_merger.rs](../../../../src/var_calling/per_group_merger.rs)) are byte-identical to the deleted columnar path — the byte-identity-sensitive math was genuinely copied, not paraphrased.
- **Sound concurrency with zero `unsafe`.** `unsafe_concurrency` confirmed the `thread::scope` channel-close ordering terminates all three thread kinds on every path (including a dead caller), the parallel-decode-then-serial-fold is order-independent (`max`/union aggregation), `first-error-wins` loses no error/panic, and the only interior mutability reachable across callers is a `thread_local!` `SHAPE_CACHE`.
- **The reader's typed move-out getters** (`take_seq` / `take_chain_ids` / `take_fixed`, [sample_reader.rs:430-475](../../../../src/var_calling/sample_reader.rs#L430-L475)) are already shaped (`&mut self`, move-out, once-only) for the deferred Phase-5 column-selective decode — the structure earns the future win without paying for it now.
- **Strong in-module equivalence testing.** `streaming_matches_reference_across_chunk_sizes` / `run_walks_chromosomes_and_intervals` ([cohort_integration.rs](../../../../src/var_calling/cohort_integration.rs)) cross-check the streaming producer against a one-shot reference across `target_variants ∈ {1,3,17,100000}` and `min_alt_obs ∈ {1,2}`, asserting `chunk_order` gaplessness and global sort — exactly the chunk-size-independence the design rests on.
- **Clear domain vocabulary** (the project's stated goal): `chunk_order` is consistent end-to-end (producer stamp → `PileupCohortChunk` → `CalledChunk` → writer `next_expected`); `CallStats` (per-chunk, caller-side) vs `WriterStats` (run-level) is a deliberate, documented split, not accidental overlap.

## 10. Commands to re-verify

In the dev container (`./scripts/dev.sh bash -c '…'`):

- `cargo fmt --check` — expect clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — expect clean.
- `cargo test --lib && cargo test --test '*'` — expect 996 lib + all integration green.
- `cargo doc --no-deps` — **currently red (M1)**; re-run after fixing the 7 unresolved links + 5 redundant-link warnings; expect clean.

New invocations the review introduces (once the §8 tests land): the `VcfWriter` reorder unit tests, the multi-thread `read_samples` test, the `from_block` boundary/malformed tests, and the `run_var_calling` error-path test.

### Author response convention
Address each finding by identifier (B/M/Mi) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the §4 open questions first (they gate M2, M5, M7, Mi2).
