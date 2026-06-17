# Code Review: ssr_pileup (perf branch)
**Date:** 2026-06-17
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** branch `ssr-pileup-review` vs `main` (`git diff main..HEAD`, merge-base `b8ad8b4`, HEAD `9c47399`) — SSR Stage-1 pileup performance work (pair-HMM math rewrites, default window 10→6, per-worker CRAM container cache)
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** PR/branch diff (8 commits) — the changed `*.rs` + `Cargo.toml` surface only.
- **Reviewed against:** `main` (merge-base `b8ad8b43`), branch HEAD `9c47399`.
- **In-scope files:**
  - [src/ssr/pileup/pair_hmm.rs](../../../../src/ssr/pileup/pair_hmm.rs) — `ln_sum_exp2`/`ln_sum_exp3` rewrites + shared-prefix DP (`score_candidates`, `longest_common_prefix_len`, `fill_row`, `PairHmmScratch.seam`) + equivalence test
  - [src/bam/segment_reader.rs](../../../../src/bam/segment_reader.rs) — **new** `CachingCramReader` / `WorkerReader` / `ContainerCache` + gate test (new code only; ~L916-1040, test ~L1716+)
  - [src/ssr/pileup/driver.rs](../../../../src/ssr/pileup/driver.rs) — `DEFAULT_WINDOW` 10→6; `WorkerReader` `map_init` wiring; ignored `ssr_psp_concordance` test; `LocusScratch::parts_mut`
  - [src/ssr/pileup/fetch_reads.rs](../../../../src/ssr/pileup/fetch_reads.rs) — `fetch_locus_reads` now takes `&mut [WorkerReader]`
  - [src/ssr/pileup/bench_harness.rs](../../../../src/ssr/pileup/bench_harness.rs) — new `#[doc(hidden)]` bench seam
  - [src/ssr/pileup/mod.rs](../../../../src/ssr/pileup/mod.rs) — `pub mod bench_harness`
  - [src/ssr/catalog/io.rs](../../../../src/ssr/catalog/io.rs) — new ignored `filter_catalog_to_regions` fixture tool only
  - [benches/ssr_pileup_perf.rs](../../../../benches/ssr_pileup_perf.rs), [examples/profile_ssr_pileup.rs](../../../../examples/profile_ssr_pileup.rs), [Cargo.toml](../../../../Cargo.toml) — new bench target
- **Deliberately out of scope:** the design/perf docs (context); the SNP pileup path; the pre-existing `segment_reader.rs` machinery that predates this branch (`BamFile`/`CramFile`/`CramSegmentReads`/`classify_segment_record`/the pooled readers — pre-existing Blocker-class only); noodles and vendored deps; SSR catalog Stage-0 logic; `read_analysis.rs` (unchanged — consistency-checked only).
- **Categories dispatched:** reliability (always; test-gap heavy), errors (always), naming (always), defaults (window + cache-cap defaults changed), idiomatic (always), refactor_safety (the byte-identity contract is the crux), module_structure (multi-file + new `pub mod`), unsafe_concurrency (`map_init` + per-worker readers + `Arc` shares), smells (always), tooling (`Cargo.toml` + bench target + doc gate), extras (hot path + stable on-disk output + PR intent-match).

## 2. Verdict

**Request-changes.** Two things gate merge: a verified `cargo doc` gate break the diff introduced (**B1**), and the branch's central claim — that the pair-HMM math rewrites are byte-identical to `main` — is **asserted but never verified against `main`** (**M1**); the one in-tree bit-identity test compares the new code against itself. Everything else is sound: the shared-prefix DP is provably bit-identical (and well-gated), the per-worker reader concurrency is correct, and the window-default change is well-documented. Fix B1 (trivial), resolve M1 with the concordance run the branch already scaffolds, and address the `decoded==0` walk divergence (M2) + its missing test (M3); the rest are quality.

## 3. Execution status

Run in the dev container (`./scripts/dev.sh cargo …`), verbatim:

- `cargo fmt --check` — **exit 0** (clean).
- `cargo clippy --all-targets --all-features -- -D warnings` — **exit 0** (clean).
- `cargo doc --no-deps` — **FAILED**:
  ```
  error: unresolved link to `AlignmentFile`
   --> src/ssr/pileup/fetch_reads.rs:9:9
  error: could not document `pop_var_caller`
  ```
- `cargo test --lib` — **exit 0**: `test result: ok. 1167 passed; 0 failed; 3 ignored; 0 measured; 0 filtered out; finished in 33.33s`
- `cargo test --all-targets --all-features` — **FAILED only at** the pre-existing out-of-scope bench: `thread 'main' panicked at benches/psp_writer_perf.rs:386:60` / `error: test failed, to rerun pass --bench psp_writer_perf`. This is the known pre-existing `--all-targets` failure documented repeatedly in PROJECT_STATUS; not introduced by this branch. (Output was tailed, so the per-suite lib counts above come from the dedicated `--lib` run.)
- `cargo audit` — **not run** (cargo-audit historically unavailable in-container per PROJECT_STATUS; no new deps this branch, so nothing in scope to audit).

**Findings labeled "Needs verification": 2** (M1 — byte-identity of the math vs `main`; M2 — `decoded==0` mid-walk reachability for real CRAM).

## 4. Open questions and assumptions

1. **(gates M1)** Has anyone diffed a `.ssr.psp` produced on `main` against one produced on this branch *with the window pinned to 10* (so only the math differs), over a real catalog+CRAM? The ignored `ssr_psp_concordance` test is built exactly for this but I could not run it (needs a real CRAM fixture). The "byte-identical" claim rests entirely on this.
2. **(gates M2)** Can a real, well-formed CRAM ever have a `.crai` record pointing at a container that `read_container` decodes to zero records *before* a later in-segment container? If never (zero-decode only at the EOF marker, which the `.crai` never indexes), M2 is latent-only and can drop to Minor.
3. **(affects M4)** Is `DEFAULT_MAX_CACHED_CONTAINERS` intended to stay a hardcoded constant, or become tunable? If it stays fixed it is still not recoverable from a produced artefact (not in the `.ssr.psp` header) — see M4.
4. **(affects Mi5)** Should the two `#[ignore]` env-var manual tools (`ssr_psp_concordance`, `filter_catalog_to_regions`) live as committed tests, move to `examples/`, or be deleted after the investigation?

## 5. Top 3 priorities

1. **B1 — fix the broken `cargo doc` gate** ([fetch_reads.rs:9](../../../../src/ssr/pileup/fetch_reads.rs#L9)). Trivial: re-point the module-doc link from `AlignmentFile` to `WorkerReader` and refresh the now-stale prose. The project sets `broken_intra_doc_links = "deny"`, so this is a hard CI failure.
2. **M1 — verify (or downgrade) the math byte-identity claim** ([pair_hmm.rs:150-187](../../../../src/ssr/pileup/pair_hmm.rs#L150-L187)). `ln_1p` and the single-pass reduction are ULP-different from `main`; the only test compares new-vs-new. Run the shipped concordance diff with the window pinned to 10 and require 0% profile diff, or restate the claim as "algebraically equivalent" and re-baseline.
3. **M2 + M3 — close the CRAM cache gaps**: mirror `refill`'s `decoded==0` walk-stop semantics ([segment_reader.rs:1010-1037](../../../../src/bam/segment_reader.rs#L1010-L1037)) and add a multi-container / cache-eviction case to the single-container gate test.

## 6. Findings

### Blocker

**B1: [src/ssr/pileup/fetch_reads.rs:9](../../../../src/ssr/pileup/fetch_reads.rs#L9) — `cargo doc` gate broken by an unresolved intra-doc link**
**Categories:** reliability, errors, naming, idiomatic, refactor_safety, module_structure, smells, tooling, extras (convergent — surfaced by 9 of 11 sub-agents)
**Confidence:** High (reproduced — see §3 verbatim).
The module doc still reads `[`AlignmentFile`]s (the pooled indexed-segment reader)`, but the import was changed from `AlignmentFile` to `WorkerReader` and `fetch_locus_reads` now takes `&mut [WorkerReader<'_>]`. With the crate's `-D rustdoc::broken-intra-doc-links` this is a hard error, so `cargo doc --no-deps` fails to build — a committed gate the project has previously treated as release-blocking (cf. the `ssr_types` B1). The surrounding prose (lines 9, 16, 146) also now mis-describes the input type.
**Fix:** re-link to `[`WorkerReader`]` (or `[`super::super::bam::segment_reader::WorkerReader`]`) and reword lines 9/16 to describe the per-worker reader; verify with `./scripts/dev.sh cargo doc --no-deps`.

### Major

**M1: [src/ssr/pileup/pair_hmm.rs:150-187](../../../../src/ssr/pileup/pair_hmm.rs#L150-L187) — the `ln_sum_exp2`/`ln_sum_exp3` rewrite is ULP-different from `main`, and its byte-identity-to-`main` is untested**
**Categories:** refactor_safety, extras, reliability (cross), smells (cross) — convergent.
**Confidence:** Medium (Blocker-class risk filed at Major per the severity/confidence rule).
`main`'s `ln_sum_exp2` is `m + ((a-m).exp() + (b-m).exp()).ln()`; the branch's is `m + (other-m).exp().ln_1p()` ([pair_hmm.rs:161](../../../../src/ssr/pileup/pair_hmm.rs#L161)). `(1.0 + x).ln()` and `x.ln_1p()` are *not* the same f64 bits — `ln_1p` is deliberately more accurate near 0. `ln_sum_exp3` also changed from nested `ln_sum_exp2(ln_sum_exp2(a,b),c)` to a single-pass `for x in [a,b,c]` reduction, reordering the floating-point sum. A direct probe of the two `ln_sum_exp2` formulas over 2,000,000 random `(a,b)` ∈ `[-60,0]` (run by the refactor_safety sub-agent) gave `f64-bit-diffs=94828 (4.741%) f32-diffs=0`: ~4.7% of single evaluations differ in f64 bits. The forward DP chains thousands of these per read, so the *final* loglik differs from `main` at the ULP scale. The stored value is `(ll - z) as f32` after `prune_and_renormalize`, and the single-op probe shows no per-op f32 flip — but two end-to-end decision points are *not* bounded by it: (1) the keep/drop test `*ll >= max - AMB_LL_DROP` flips if a candidate sits within ULP-drift of the boundary; (2) accumulated DP drift before the `as f32` cast can exceed a single op's drift. The in-tree `score_candidates_is_bit_identical_to_per_candidate_forward` test compares `score_candidates` against `forward` **using the new math on both sides** — it cannot detect new-vs-`main` divergence. The "byte-identical" claim (module doc / PROJECT_STATUS) is therefore asserted, not verified, for the only comparison that matters.
**Why it matters:** the branch contract is byte-identity of `.ssr.psp` to `main` for everything except the window. If one locus flips a pruned length or an f32 bit, the contract is broken and the commit is not the byte-identical refactor it claims.
**Fix:** build `.ssr.psp` on `main` and on the branch with the window pinned to 10 over a real catalog+CRAM, then `PVC_PSP_A=main.ssr.psp PVC_PSP_B=branch.ssr.psp cargo test --release ssr_psp_concordance -- --ignored --nocapture` and require `profile differs: 0 (0.00%)`. If non-zero, either revert the math to `main`'s exact ops or drop the "byte-identical" wording, document the change as ULP-affecting, and re-baseline downstream goldens.

**M2: [src/bam/segment_reader.rs:1010-1037](../../../../src/bam/segment_reader.rs#L1010-L1037) — `CachingCramReader` `.crai` walk diverges from `refill` on `decoded==0` (continues vs. stops)**
**Categories:** refactor_safety, reliability, errors — convergent.
**Confidence:** Medium (Needs verification — see Q2).
`CramSegmentReads::refill` treats a container that decodes to zero records as end-of-walk ([segment_reader.rs:838-840](../../../../src/bam/segment_reader.rs#L838-L840): `if decoded == 0 { return Ok(true); }`), stopping `next()` entirely. The new `decode_container` instead returns `Ok(Vec::new())` ([segment_reader.rs:1075-1077](../../../../src/bam/segment_reader.rs#L1075-L1077)) and the `for idx in 0..len` loop **continues to the next index record**. For well-formed CRAM the two agree (zero-decode only at the genuine EOF marker, which the `.crai` doesn't index), but if any `.crai` record points at a zero-decode container *before* a later in-segment container, the oracle yields nothing further while the new reader yields the later container — not byte-identical. The gate test (`caching_cram_reader_matches_per_call_path`) uses a single-container fixture and never exercises a mid-walk `decoded==0`. (Note: the new path is arguably the *more* correct one — but the contract is "identical to the per-call path", so they must agree.)
**Fix:** mirror `refill` — `let decoded = self.get_or_decode(record.offset())?; if decoded.is_empty() { break; }` — and add a two-container fixture (a zero-record container followed by an overlapping one) to the gate test.

**M3: [src/bam/segment_reader.rs:928-1040](../../../../src/bam/segment_reader.rs#L928-L1040) — `ContainerCache` eviction + multi-container path is untested**
**Categories:** reliability.
**Confidence:** High.
The gate test fixture builds a single CRAM container, so with `cap=3` the FIFO eviction (`pop_front`), the dedup-on-insert (`if … any(|(o,_)| *o == offset) { return; }`), the `cap.max(1)` floor, and the `Arc`-outlives-eviction rationale **never execute**. The cache's correctness-relevant policy has zero direct coverage; a bug in eviction order or dedup would pass CI.
**Fix:** add a unit test on `ContainerCache` directly (insert > cap distinct offsets, assert the oldest is evicted and a re-inserted offset is a no-op) and a multi-container `fetch_mapped_reads` case with `cap` smaller than the number of overlapping containers, diffed against `get_reads_from_segment` (closes M2's fixture need too).

**M4: [src/bam/segment_reader.rs:926](../../../../src/bam/segment_reader.rs#L926) — `DEFAULT_MAX_CACHED_CONTAINERS` is not inspectable at runtime**
**Categories:** defaults.
**Confidence:** High.
Unlike `window`/`reservoir_cap` (clap `[default: …]` + stamped into the `.ssr.psp` header parameters), the cache cap is hardcoded in `worker_reader`, exposed by no CLI/env/accessor, and is **not** written to the `.ssr.psp` header. Its effective value is unrecoverable from a running instance or a produced artefact — "read the source" is the only way to learn it, failing the defaults "inspectable/logged when applied" rule.
**Fix:** record it in the `.ssr.psp` header parameters alongside `window`/`reservoir_cap` (no CLI flag needed — it's a cache-sizing knob, not result-affecting). Cheap symmetry with the SNP header convention.

**M5: [src/bam/segment_reader.rs:1009-1110](../../../../src/bam/segment_reader.rs#L1009-L1110) — the `.crai` walk + container decode is duplicated between `fetch_mapped_reads`/`decode_container` and `CramSegmentReads::refill`**
**Categories:** smells, refactor_safety (cross), module_structure (cross), extras (cross) — convergent.
**Confidence:** High.
Two copies of correctness-critical container-overlap arithmetic (the `reference_sequence_id`/`alignment_start > segment.end` break/`container_end < segment.start` continue cascade) and the per-slice decode loop (`compression_header` → `slices()` → `decode_blocks` → `records` → `try_from_alignment_record`) now exist side by side. M2 is the concrete cost of this duplication: the two copies have *already* drifted (`decoded==0` handling). The author's own comment flags the parallel.
**Fix:** extract a shared `overlapping_container_offsets(&index, target, &segment) -> impl Iterator<Item=&crai::Record>` (with the early-break) and a `decode_container_records(reader, repository, header, offset) -> Result<Vec<RecordBuf>>`; have both `refill` and `decode_container` call them so they cannot diverge.

### Minor

- **Mi1: stale `window = 10` across the new perf scaffolding.** [bench_harness.rs:44-45](../../../../src/ssr/pileup/bench_harness.rs#L44-L45) ("the `ssr-pileup` default is 10"), [benches/ssr_pileup_perf.rs](../../../../benches/ssr_pileup_perf.rs) (comments + the hardcoded/swept `10`), [examples/profile_ssr_pileup.rs:14](../../../../examples/profile_ssr_pileup.rs#L14) (`window=10` default + doc). All now misdescribe the shipped default of 6, and the bench/profile no longer measure at the production operating point. **Categories:** defaults, smells, reliability, extras (convergent). Fix: update the prose to 6 and set the bench/example baseline window to `DEFAULT_WINDOW` (import it).
- **Mi2: the `bench_harness` seam widens the *public* API and compiles into the production library.** [mod.rs:16](../../../../src/ssr/pileup/mod.rs#L16) `pub mod bench_harness` exposes the only genuinely-public items under `ssr/pileup/` (`SyntheticLocusWorkload`, `build_synthetic_workload`, `analyze_workload`); `#[doc(hidden)]` hides them from rustdoc but not from `use`, contradicting the module's own "not part of the supported surface" doc. **Categories:** module_structure, tooling (convergent). For an internal lab CLI this is acceptable as a *deliberately flagged* widening, but the clean fix is `#[cfg(any(test, feature = "bench-internals"))]` on the module so it does not ship in normal builds.
- **Mi3: [src/bam/segment_reader.rs:1086](../../../../src/bam/segment_reader.rs#L1086) — `.expect("opened by ensure_open")` lacks the project-required `// PANIC-FREE:` invariant comment.** The invariant holds (`ensure_open()?` populates `reader` first), but the rule mandates the named comment or a restructure. **Category:** errors. Cleaner fix: have `ensure_open` *return* `&mut Reader`, eliminating the `expect`.
- **Mi4: [src/ssr/pileup/pair_hmm.rs:597-660](../../../../src/ssr/pileup/pair_hmm.rs#L597-L660) — the bit-identity test omits a single-candidate / `n==lcp` case.** That is the corner where Pass 2's tail loop is empty and the final read comes straight from the seam (logic verified correct, but least-exercised). **Categories:** refactor_safety, reliability. Fix: add a `build_rungs(&locus, obs, 0, …)` single-rung case asserting `to_bits()` equality vs `forward`.
- **Mi5: two `#[ignore]` env-var manual tools committed to the tree.** `ssr_psp_concordance` ([driver.rs:789](../../../../src/ssr/pileup/driver.rs#L789)) and `filter_catalog_to_regions` ([catalog/io.rs](../../../../src/ssr/catalog/io.rs)) are near-identical PVC_*-driven manual fixtures masquerading as tests. **Categories:** smells, reliability. Fix: move to `examples/` (or delete after the window investigation closes); if kept, give them the project's required removal condition. (Note: `ssr_psp_concordance` is also the verification vehicle for M1, so keep it at least until M1 is resolved.)
- **Mi6: [src/bam/segment_reader.rs:1089-1100](../../../../src/bam/segment_reader.rs#L1089-L1100) — `decode_container` spells out `AlignmentInputError::Io { path: self.path.clone(), source }` inline twice instead of using the `io_error` helper, and the justifying comment misstates the borrow reason.** **Category:** smells. Fix: the seek/read errors don't need disjoint-field borrows (no `&self` method call mid-expression) — route them through `self.io_error(e)` like the rest, and drop the misleading comment.
- **Mi7: [src/ssr/pileup/pair_hmm.rs:296](../../../../src/ssr/pileup/pair_hmm.rs#L296) — `#[allow(clippy::too_many_arguments)]` on `fill_row` has no justification comment.** **Category:** smells. Fix: add a one-line rationale (hot-path row kernel; a params struct would add indirection) or bundle `(read_base, match_ln, mismatch_ln, ins_emit)` into a small `RowEmission` struct.
- **Mi8: [src/bam/segment_reader.rs:934](../../../../src/bam/segment_reader.rs#L934) — `ContainerCache` holds `Arc<Vec<RecordBuf>>` for data frozen after `decode_container`; prefer `Arc<[RecordBuf]>`.** Drops one indirection on the hot fetch path; the `Arc` itself is correct (served windows must outlive eviction). **Category:** idiomatic. Fix: `Arc::from(out)` at insert, type `Arc<[RecordBuf]>`.
- **Mi9: no compile-time `Send` guard on `WorkerReader`.** The bound rayon actually requires is currently checked only implicitly by the driver compiling; the sibling `AlignmentFile` has an explicit `alignment_file_is_sync` test. **Category:** unsafe_concurrency. Fix: add a `fn worker_reader_is_send<T: Send>(){} ` style static assert so a future non-`Send` field fails at a named test, not a deep generic rayon error.
- **Mi10: `ContainerCache.cap` vs `Reservoir.capacity` — same concept ("max items held"), two spellings, shipping together.** **Category:** naming. Fix: `capacity` for both.
- **Mi11: the new hot-path bench measures the ~54% win but does not guard it** ([benches/ssr_pileup_perf.rs](../../../../benches/ssr_pileup_perf.rs)). No `main`-baseline / regression threshold, so a future regression won't fail CI. **Category:** extras. Fix: note the expected throughput in a comment or wire a criterion baseline check (acknowledging criterion benches aren't CI-gated here).
- **Mi12: [src/ssr/pileup/driver.rs:399-411](../../../../src/ssr/pileup/driver.rs#L399-L411) — the batch-fill `?` propagates a per-locus error with no positional context** (which locus/index failed). **Category:** errors. Downgradeable if `SsrPileupError`/`CatalogError` already pin the failing record. Fix: wrap with the locus coordinate if not.

### Nits

- [segment_reader.rs:1064-1067](../../../../src/bam/segment_reader.rs#L1064-L1067) — the doc comment above `ensure_open` actually describes `decode_container` (the decode-equivalence claim landed on the wrong function). (idiomatic, reliability)
- [segment_reader.rs:963-973](../../../../src/bam/segment_reader.rs#L963-L973) — `CachingCramReader` doc says index/header/**repository/filter** are "`Arc`-shared", but `repository` is a clone-shared `Repository` newtype and `filter` is `Copy`. Wording drift; sound either way. (unsafe_concurrency)
- [segment_reader.rs:928-961](../../../../src/bam/segment_reader.rs#L928-L961) — `ContainerCache`'s linear-scan `get`/`insert` is *correct and faster than a map* at `cap≤3`, but undocumented; a one-line "linear scan is intentional at this cap" note prevents a future "optimize to HashMap" regression. (idiomatic)
- [pair_hmm.rs:357](../../../../src/ssr/pileup/pair_hmm.rs#L357), `LocusScratch.cands` — abbreviation; the crate spells it `candidates` elsewhere. `out` in `decode_container` and `v`/`value` in `ContainerCache` are generic but short-scoped. (naming)
- [benches/ssr_pileup_perf.rs](../../../../benches/ssr_pileup_perf.rs) — the `alloc-mimalloc` doc-comment example still names `psp_writer_perf`. (tooling)

## 7. Out of scope observations

- **`benches/psp_writer_perf.rs:386` panics under `cargo test --all-targets`** — pre-existing, documented repeatedly in PROJECT_STATUS as the only `--all-targets` failure; not introduced here. Follow-up: separate fix (gate the bench's debug-build panic or move the fixture).
- **`AlignmentInputError::Io` is mechanism-named and funnels seek/read/decode/convert origins** — pre-existing; the new code adds six more `Io` sites but does not create the issue. Follow-up: a typed CRAM-decode error variant, repo-wide.
- **`fasta::Repository` is `Arc<RwLock<AdapterCache>>` internally** — concurrent per-worker CRAM fetches synchronize on one `RwLock`, a potential `--threads` scaling ceiling (correctness is fine). Pre-existing shared-repository design; flag for the next perf pass.
- **`get_reads_from_segment` `get_`-prefix naming (pre-existing Mi3 from the segment-reader review)** — the new `fetch_mapped_reads` is the better target name; converge when the SNP `--regions` retrofit (increment #5) lands.

## 8. Missing tests to add now

Grouped by unit under test. (The reliability + refactor_safety sub-agents' challenge passes feed this section.)

**`ContainerCache` (segment_reader.rs)**
- `container_cache_evicts_oldest_when_over_capacity` — insert `cap+2` distinct offsets into a `cap`-sized cache; assert the two oldest miss and the newest `cap` hit. Catches an eviction-order/`pop_front` bug (M3).
- `container_cache_insert_is_idempotent_on_repeat_offset` — insert the same offset twice with different `Arc`s; assert the first wins and len is unchanged. Catches a dedup regression (M3).

**`CachingCramReader::fetch_mapped_reads` (segment_reader.rs)**
- `caching_reader_matches_per_call_across_multiple_containers` — a multi-container CRAM fixture with `cap` < overlapping-container count; diff `(reads, FilterCounts)` against `get_reads_from_segment` over a window spanning ≥2 containers. Catches eviction-under-load + cross-container boundary divergence (M2/M3).
- `caching_reader_stops_walk_on_zero_decode_like_refill` — a fixture (or seam) with a zero-record container before an in-segment one; assert the caching path matches `refill` (after the M2 fix this should yield the `refill` result, not the extra container).

**`score_candidates` (pair_hmm.rs)**
- `score_candidates_single_candidate_is_bit_identical_to_forward` — one-rung set (`n==lcp`), assert `to_bits()` equality vs `forward` (Mi4).

**End-to-end byte-identity (the M1 gate)**
- `ssr_psp_is_byte_identical_to_main_at_window_10` — automate the concordance: produce `.ssr.psp` on `main` and on the branch with window=10 over a checked-in small real CRAM fixture; assert 0 profile diffs. This is the test that actually backs the "byte-identical" claim; the ignored `ssr_psp_concordance` is its manual ancestor.

## 9. What's good

- **The shared-prefix DP is genuinely bit-identical and proven, not hoped.** `score_candidates` reconstructs `forward`'s exact operand order from the saved seam column; the edge cases (lcp==0, n==lcp, m==0, stale-cells-below-lcp) all reduce correctly, and the `to_bits()` test ([pair_hmm.rs:597](../../../../src/ssr/pileup/pair_hmm.rs#L597)) gates it. This is the right way to ship a memoization optimization.
- **Per-worker reader concurrency is correct by construction.** Each `CachingCramReader` owns its own `File` handle + cache (no cross-thread cursor aliasing), shares only immutable `Arc`s, and `map_init` preserves `--threads`-invariance — backed by `run_output_is_thread_count_invariant` ([driver.rs:848](../../../../src/ssr/pileup/driver.rs#L848)).
- **The `DEFAULT_WINDOW` 10→6 change is documented exactly right** ([driver.rs:54-69](../../../../src/ssr/pileup/driver.rs#L54-L69)): the measurement source, the 0.48%/0.06% output impact, and the rung-count/DP saving are all at the const, and the value flows to clap `[default: 6]` so a user sees it.
- **The caching path reuses `classify_segment_record` as the single filter chokepoint**, so the `(reads, FilterCounts)` membership and ordering are identical to the per-call path by construction rather than by parallel re-implementation.
- **The cache copies the `Arc` out before serving** so an evicted container's records outlive eviction — defensive correctness even though the current consumer drains within the call.

## 10. Commands to re-verify

Reviewer ran (re-run to confirm):
- `./scripts/dev.sh cargo fmt --check` — exit 0.
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `./scripts/dev.sh cargo test --lib` — 1167 passed; 0 failed; 3 ignored.
- `./scripts/dev.sh cargo doc --no-deps` — **FAILS** (B1); must pass after the fix.

Introduced by this review (run after fixes):
- The M1 gate: `PVC_PSP_A=main.ssr.psp PVC_PSP_B=branch-window10.ssr.psp ./scripts/dev.sh cargo test --release ssr_psp_concordance -- --ignored --nocapture` → require 0% profile diff.
- New tests in §8 once added.

### Author response convention
Address each finding by identifier (B1, M1…, Mi1…) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the §4 open questions first (they gate M1/M2/M4/Mi5).
