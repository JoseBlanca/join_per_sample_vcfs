# Performance Review: var_calling (re-architected cohort pipeline)
**Date:** 2026-06-06
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** cohort `var-calling` pipeline (`.psp` per-sample files → multi-sample VCF), branch `re-architect` @ `37e02c2` — the record-streaming redesign (producer → caller workers → writer)
**Verdict:** Run experiments
**Hot-path evidence:** CPU sampling profile (macOS `sample`, T=1 30 s + T=8 6 s, N=50 real tomato cohort) — collected this run; quoted in [_profile_evidence.md](../../../../tmp/perf_review_2026-06-06_var-calling-cohort/_profile_evidence.md). DHAT attribution reused from prior work (N=50/T=8/tvpc=256).

---

## 1. Scope and constraints

- **What was reviewed:** the cohort `var-calling` module subtree on the re-architected record-streaming pipeline. Not the whole crate, not a diff — the production path from N `.psp` files to one multi-sample VCF.
- **Reviewed against:** branch `re-architect`, commit `37e02c2`.
- **Targets / input sizes / hardware:** memory-efficient multi-sample SNP caller (trades RAM for sample-count scaling — a competitor's OOM is a headline result, so peak RSS is a first-class metric). Real workload: **N=50 tomato samples, whole genome (~12 chromosomes, ~800 Mbp), run as one `--threads T` process** (cohort joint calling is not per-sample-parallel). Production defaults: `--target-variants-per-chunk 256`, hand-off queue depth `2×threads`. Profiling host: macOS Apple Silicon; production builds: Linux container. **Goal: close the wall-time gap vs `main` at T≥2 without giving up the memory advantage, and surface further memory wins.**
- **Hard constraint:** byte-identity of calls — GT/GQ/AD/AF/AC/FILTER must stay identical to `main`. **QUAL is exempt** (float-noise accepted). Every finding that risks call output is flagged and held to a correctness bar.
- **Hot-path evidence available:** a real CPU sampling profile (the load-bearing input the prior reviews lacked at this commit). Benchmark wall/RSS matrix and DHAT attribution from the handoff are real (measured 2026-06-06 / prior). Per the skill's rubric, code-level findings are promoted to **Hot-path** only where the `sample` profile names the site.
- **In-scope files:** `src/var_calling/pipeline.rs`, `cohort_integration.rs`, `sample_reader.rs`, `per_position_merger.rs`, `per_group_merger.rs`, `variant_grouping.rs`, `pileup_overlaps.rs`, `posterior_engine.rs`, `em_posterior_calc.rs`, `vcf_writer.rs`; `src/psp/reader.rs`, `src/psp/block.rs`; `src/pop_var_caller/var_calling.rs`, `common.rs`; `Cargo.toml`.
- **Deliberately out of scope:** `contamination_estimation.rs`, `dust_filter.rs` (the sdust pre-pass — separately profiled, cold on this fixture), the `pileup` / `psp-to-pileup` subcommands, vendored reference codebases under `/Users/jose/devel/pop_var_caller/{freebayes,gatk,bcftools,htslib,...}`, and the zstd C internals (cost noted, not rewritten).
- **Categories dispatched:** `methodology` (always — build config + the missing bench), `concurrency` (the T≥2 scaling gap), `allocations` (RSS is the headline metric), `hot_loops` (EM + decode self-time), `data_layout` (columnar layout + the dead column), `io_and_syscalls` (zstd decode + VCF write).

## 2. Verdict

**Run experiments.** A real sampling profile now exists, and it points cleanly at two dominant, independent levers plus several smaller ones:

- The **T≥2 wall gap** is a *scheduling* problem, not a per-record cost: the producer's rayon compaction pool and the crossbeam caller pool are **two independent pools, each sized to `--threads`**, live at the same time → oversubscription (**H1**). The profile pins it — the producer thread is **64 % blocked** waiting on rayon at T=8, with context-switch/futex frames elevated vs T=1, exactly tracking the on-par-at-T=1 / +6/+12/+15 %-at-T=2/4/8 benchmark shape.
- The **next memory win** is the **chain-id dead-weight** (**H2**): REF-allele chain_ids are 65 % of the in-flight payload and ~96.6 % of all chain_ids, yet the merger provably never reads them (`per_group_merger.rs:1304` skips allele 0). Three categories converged here.
- One **cheap, byte-identical, apply-now** codegen win exists (**H3** `ln_factorial` inlining), confirmed by `cargo asm`.

H1 and H2 need their experiments run before merge (H1 to choose among pool-unification variants and confirm the gap closes; H2's deep variant needs the byte-identity oracle), so the verdict is *Run experiments* rather than *Apply* — but H3 and the cheapest variants of H1/H2 are low-risk and can land immediately behind the byte-identity diff. **The single biggest process gap is the absence of an end-to-end cohort criterion bench** (deleted with the from-bam removal); building it is the first item in the measurement plan because it makes every code finding here rankable.

## 3. Measurement plan

In dependency order:

1. **Build the missing cohort end-to-end bench (unblocks everything else).** There is no committed cohort `var-calling` bench (`benches/` has only psp_reader/psp_writer/baq/pileup_walker/freebayes; the `cohort_e2e_perf` bench was deleted). Add `benches/cohort_var_calling_perf.rs` (`harness = false`, criterion). Drive `run_var_calling` (or the `pipeline.rs` producer→caller→writer scope over pre-opened readers, to keep file IO out of the timed region) over a synthetic N-replicate cohort built the way `examples/profile_cohort_e2e.rs` already replicates one `.psp` into N samples. Sweep `N ∈ {1,8,50}` × `threads ∈ {1,2,8}` at the **production `target_variants_per_chunk = 256`**; `black_box` every input; build the fixture **once outside `b.iter`**; assert on returned `WriterStats` so a chunking/fixture regression trips the bench. Size the fixture to one chromosome / fixed window so a sample completes in seconds. Pair with the `dhat_var_calling` heap run at matching N so the memory-vs-speed trade is rankable. **Act threshold for any code finding: a cross-commit median shift ≥ 5 % that survives a revert experiment** (per the methodology rule that one within-run CI is not trustworthy).
2. **Fix the maintained profiling driver's fidelity.** `examples/profile_cohort_e2e.rs:161` hard-codes `target_variants_per_chunk: 0` (legacy single-pull) — the opposite of the production chunked/streaming shape this review targets. Default it to `256` (or add the `--target-variants-per-chunk` flag `dhat_var_calling.rs` already has) and confirm the producer-blocked-on-rayon frame (`pipeline.rs:296`) reproduces.
3. **H1 oversubscription sweep** (the wall-gap experiment). On the N=50 workload at T=2/4/8, wall to last VCF byte, best-of-3 median, **diff GT/GQ/AD/AF/AC/FILTER vs `main` (QUAL exempt)**:
   - Variant 1 (cheapest): cap the rayon pool below `--threads` so rayon + callers ≈ `T` total.
   - Variant 3: serialize `compact_samples` (drop its `par_iter_mut`) — test whether per-sample compaction is even worth parallelizing once the callers carry the EM in parallel.
   - Variant 2 (structural, highest byte-identity risk): move `compact_samples` onto the caller threads (ship a still-compressed chunk), collapsing the two pools into one. Defer behind 1/3.
   - Confirm via a fresh T=8 `sample` capture that `swtch_pri` / `__ulock_wait2` / `__psynch_mutexwait` fall toward T=1 levels.
4. **H2 chain-id drop** — DHAT before/after with the command in the brief (`--n-samples 50 --threads 4 --target-variants-per-chunk 256`). Step 1 (REF slot → empty `Vec` in `records_all`/`records_for`): expect `records_all`'s 262 MB and `append_range`'s 662 MB / 20.8 GB churn to fall by the REF chain-id share. Step 2 (CSR stores non-REF alleles only): expect the inflate/copy traffic itself to drop. **Byte-identity gate: zero VCF diff** + the existing `sample_reader.rs` round-trip tests pass.
5. **H3 `ln_factorial` inline** — apply, re-run `cargo asm` to confirm the four `bl …ln_factorial` become inline table loads, then re-`sample` and confirm the `ln_factorial` leaf disappears and `compute_log_likelihoods` self-time drops.
6. **Secondary sweeps** (after H1 lands, since it moves the floor): queue depth `QUEUE_DEPTH_PER_WORKER ∈ {1,2,4}` tracking wall **and** peak RSS; allocator A/B (`--features alloc-mimalloc`) at N=50/T=8.

## 4. Build / toolchain configuration

The release profile is already well-tuned — **do not change it as part of code-level work**:

- [Cargo.toml:45-49](../../../../Cargo.toml#L45-L49): `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`, `debug = "line-tables-only"`, `opt-level` inherits 3. LTO/codegen-units/panic are optimal; there is no whole-program-inlining lever left to pull. `[profile.bench]` inherits release with `debug = true`, so benches are profilable.
- **Allocator A/B is the one unstruck build lever.** The profile shows ~3.5k self-time samples in the macOS system allocator (`_xzm_free 1438`, `_xzm_xzone_malloc_tiny 1268`) under a multithreaded pipeline — the shape `mimalloc`/`jemalloc` per-thread caches target. The `alloc-mimalloc` feature already exists ([Cargo.toml:107](../../../../Cargo.toml#L107)) but is never A/B-measured on the cohort path (no cohort bench to measure against). Measure once the bench (item 1) lands; one hypothesis per measurement; track RSS too (memory-efficiency thesis). **Speculative**, gated.
- **`[profile.profiling]` with full `debug = true` is absent** despite project memory describing one; release ships `debug = "line-tables-only"`. It did not block this `sample` capture (function names resolved), but degrades DHAT inlined-frame attribution. Cheap to add (`inherits = "release"`, `debug = true`); build profiling drivers with `--profile profiling`. **Speculative.**

## 5. Code-level findings

### Hot-path

**H1: [src/pop_var_caller/var_calling.rs:238](../../../src/pop_var_caller/var_calling.rs#L238) + [src/var_calling/cohort_integration.rs:745](../../../src/var_calling/cohort_integration.rs#L745),[926](../../../src/var_calling/cohort_integration.rs#L926) + [src/var_calling/pipeline.rs:212-243](../../../src/var_calling/pipeline.rs#L212-L243) — Two independent thread pools, each sized to `--threads`, run concurrently → oversubscription (the T≥2 wall-gap driver)**
- **Confidence:** High
- **Hot-path evidence:** _profile_evidence.md, T=8 producer (main) thread, 4406 samples: **2806 (64 %) BLOCKED** in `pipeline.rs:296` → `produce_chunk` → rayon `bridge_producer_consumer` → `in_worker_cold` → `LockLatch::wait_and_reset` → `__psynch_cvwait`. Topology note: "up to 8 rayon + 8 caller + 1 writer + producer-main contend for the cores … two independent pools, each sized to `--threads`, running concurrently → oversubscription." Oversubscription signatures elevated vs T=1: `swtch_pri 750`, `__ulock_wait2 309`, `__psynch_mutexwait 156`, `__ulock_wake 80`. Benchmark matrix: on-par/faster at T=1, wall +6/+12/+15 % at T=2/4/8 — the loss grows with T exactly as oversubscription predicts.
- **Pattern matched:** rayon fork-join overhead + general contention — `configure_rayon_pool(args.threads)` sizes the global rayon pool to `--threads`, which drives the producer's `read_samples`/`compact_samples` `par_iter_mut`; `pipeline.rs` *separately* spawns `--threads` crossbeam caller threads + a writer in the same `thread::scope`. That is up to `2·T + 2` runnable CPU-bound threads on `T` cores.
- **Mechanism:** the scheduler must time-slice ~`2·T` threads onto `T` cores. Every producer rayon burst wakes/sleeps `T` rayon workers while `T` caller threads are mid-EM; the extra context switches / futex waits are pure scheduling overhead absent at T=1 (where one rayon worker ≈ the main thread, no contention). This is the parallel-scaling ceiling.
- **Container caveat (the highest-risk default).** When `--threads` is *omitted*, both `configure_rayon_pool(None)` and `n_workers` fall back to `available_parallelism()` — which on Linux reads the **cgroup CPU limit, not the request/quota the process actually owns** (this was PostHog's literal root cause in [*Untangling Tokio and Rayon in production*](https://posthog.com/blog/untangling-rayon-and-tokio)). Under a Kubernetes/container CPU quota the 2× oversubscription is then layered on top of a core count the scheduler will *throttle*, so the no-`--threads` path is more exposed than the explicit-`--threads` sweeps below. Both the asymmetric "starve the parked side" recipe and the rayon-cap variant assume both pools are CPU-bound here (producer = zstd decode + compaction, callers = EM) — neither side is "mostly parked," so the goal is strictly **rayon + callers ≈ owned-core-count**, and the pool size should derive from the owned quota, not `available_parallelism()`'s view of the limit. When validating on the Linux container, watch cgroup `cpu.stat` `nr_throttled`/`throttled_time` (treat >30 % throttling as "still oversubscribed") alongside wall.
- **Measurement plan:** §3 item 3 — the three variants, cheapest first, each gated on byte-identical calls and a fresh T=8 `sample` confirming the oversubscription frames drop.
- **Complexity cost:** Variant 1 (cap rayon): trivial (`min`/subtraction at the call site + a justifying comment); risk is under-feeding decode if the reserve is too large. Variant 3 (serialize compaction): small (delete `par_iter_mut`, serial `map`). Variant 2 (compaction onto callers): largest — `RawCohortChunk` would carry per-sample `TwoPhaseSegment`s instead of finalized `SamplePspChunk`s, the caller must reproduce the producer's exact `kept_all` row selection, **high byte-identity risk**; defer.
- **Suggested experiment / fix:** start with variant 1.
  ```rust
  // src/pop_var_caller/var_calling.rs — cap the producer's rayon pool so it does
  // not contend 1:1 with the crossbeam caller pool. EXPERIMENT: sweep the reserve.
  let rayon_threads = args.threads.map(|t| t.saturating_sub(t / 2).max(1));
  configure_rayon_pool(rayon_threads)
      .map_err(|_| VarCallingCliError::RayonAlreadyConfigured)?;
  ```
  (The exact split is what the T=2/4/8 sweep determines; the point is rayon + callers ≈ `T`, not `2·T`.)

**H2: [src/var_calling/sample_reader.rs:472-513](../../../src/var_calling/sample_reader.rs#L472-L513) (`records_all`/`records_for`) + [set_variable_rows](../../../src/var_calling/sample_reader.rs#L751-L759), consumed at [per_group_merger.rs:1303-1307](../../../src/var_calling/per_group_merger.rs#L1303-L1307) — REF-allele chain_ids are decoded, inflated, copied, and materialized end-to-end but never read (chain-id dead-weight)**
- **Confidence:** High (converged independently by `allocations`, `data_layout`, and `hot_loops`)
- **Hot-path evidence:** DHAT (N=50/T=8/tvpc=256, peak 1.352 GB): `append_range` 662 MB (49 %, 20.8 GB churned), `records_all` 262 MB (19 %); **65 % of the in-flight payload is chain_ids; ~96.6 % of chain_ids are REF-allele ids.** CPU profile: `decode_list_column_csr 3303` (#2 leaf), `PerAlleleChainIds::extend_from_range 218`, the `records_all`-side `to_vec()`, and the ~7.3k memmove/memset bucket. **Verified directly:** the only merger consumer is `per_group_merger.rs:1307` (`for &chain_id in &allele.chain_ids`), unreachable for REF because of the explicit `if local_allele_idx == 0 { continue; }` at `:1304`; `PileupRecord::ref_span()` reads `alleles[0].seq`, not chain_ids; `variant_grouping.rs` and `vcf_writer.rs` never touch chain_ids.
- **Pattern matched:** dead data copied through the hot path — REF chain_id cells are copied 3–4× (inflate → `extend_from_range` → `append_range` → `to_vec` per allele in `records_all`) for a column no consumer reads.
- **Mechanism:** every chain-id cell is a `u64`; REF cells are 96.6 % of them. Dropping them removes ~96.6 % of the chain-id memmove traffic and the matching per-REF-allele `Vec<ChainId>` allocations, shrinking the chain-id component (65 % of the in-flight payload).
- **Measurement plan:** §3 item 4. Merge threshold: peak-live drop ≥ ~15 % **or** churn drop ≥ ~30 %, with **zero VCF diff**.
- **Complexity cost:** Step 1 is cheap and self-contained (pass `Vec::new()` for the REF allele — allocation-free until pushed-to — instead of `slice_at(lo).to_vec()` in `records_all`/`records_for`); new invariant "rebuilt REF alleles carry an empty `chain_ids`," to be documented at the carrier and re-checked if any future consumer reads `alleles[0].chain_ids`. Step 2 (CSR stores non-REF alleles only, so the 96.6 % are never inflated/copied) is a real layout change touching `from_block`/`set_variable_rows`/`extend_from_range`/`records_all` with byte-exact offset arithmetic — gate on the same byte-identity tests, do only if step-1 DHAT shows the inflated buffers still dominate. No `unsafe`, no dependency.
- **Suggested experiment / fix:**
  ```rust
  // Step 1 (cheap, no layout change): in records_all / records_for, do not clone
  // REF chain_ids — the merger skips allele 0 (per_group_merger.rs:1304).
  for a in lo..hi {
      let chain_ids = if a == lo { Vec::new() } else { self.chain_ids.slice_at(a).to_vec() };
      alleles.push(AlleleObservation::new(self.seq.slice_at(a).to_vec(), self.fixed.support_at(a), chain_ids));
  }
  // Step 2 (layout, gated): store PerAlleleChainIds for non-REF alleles only.
  ```

**H3: [src/var_calling/per_group_merger.rs:2151](../../../src/var_calling/per_group_merger.rs#L2151) — `ln_factorial` is not `#[inline]`; it stays an out-of-line `bl` on the per-(sample × genotype × allele) likelihood loop**
- **Confidence:** High
- **Hot-path evidence:** profile leaves `ln_factorial 448` + `compute_log_likelihoods 1512`. `cargo asm --lib --simplify` of `compute_log_likelihoods` shows four `bl …per_group_merger::ln_factorial` per genotype evaluation (called, not inlined). Smaller share than H1/H2 (~1 % of non-idle self-time) but High-confidence, zero-risk, trivial.
- **Pattern matched:** `#[inline]` on a small, very-frequently-called function — `ln_factorial`'s hot path (`n < 1024`) is one bounds-checked table load, but carries no attribute.
- **Mechanism:** each call is a branch-and-link + register save/restore around what should be one `LDR d, [table, n*8]`; paid `n_samples · n_genotypes · (kept_alleles+1)` times per record. `#[inline]` lets LLVM fold the table load into the caller and keep the accumulator in registers; the cold `n ≥ 1024` fallback splits into a `#[cold] #[inline(never)]` tail.
- **Measurement plan:** §3 item 5 — `ln_factorial` leaf disappears, `compute_log_likelihoods` self-time drops measurably; cheap `cargo asm` A/B confirms inline loads.
- **Complexity cost:** one attribute + an optional cold-tail fn. No `unsafe`, no new invariant. Pure codegen — byte-identical (does not touch GT/AD/QUAL math).
- **Suggested experiment / fix:** `#[inline]` on the table-lookup fast path, `#[cold] #[inline(never)]` on the iterative `n ≥ 1024` tail (snippet in the category file).

### Likely

**L1: [src/var_calling/cohort_integration.rs:827-897](../../../src/var_calling/cohort_integration.rs#L827-L897) + [pipeline.rs:264-308](../../../src/var_calling/pipeline.rs#L264-L308) — Producer is fully serial on the main thread between sends; chunk N's shaping cannot overlap chunk N−1's EM**
- **Confidence:** High (structural) / Medium (that it is the residual floor after H1)
- **Hot-path evidence:** profile: producer ~700 inline fold samples at `pipeline.rs:296` (the non-blocked remainder after the 64 % rayon wait); design note `pipeline.rs:46-52` "the producer is the wall floor"; project memory `producer_record_building_bottleneck` ("serial records_for + per-position merge on the producer thread caps wall-scaling past ~T=4").
- **Pattern matched:** the serial fraction of a producer/consumer pipeline — one producer thread runs the light-column fold + safe-gap cut + dust mask + `kept_all`/`variable` collects + REF fetch serially *between* the rayon bursts, then `send`s.
- **Mechanism:** the producer cannot build chunk N+1 while callers consume chunk N; by Amdahl this serial section caps `W`-way speedup. `target-variants-per-chunk=256` mitigates (smaller work-units feed callers sooner) but does not remove the floor.
- **Measurement plan:** **measure-only, gated behind H1** — after the pool fix, re-`sample` the producer (the BLOCKED-in-rayon share should fall; whatever remains as non-blocked serial self-time is the floor) and take the wall-vs-T curve at T=1,2,4,8. Act if caller-thread idle (high `__psynch_cvwait` on `chunk_rx`) exceeds ~25 % at the knee.
- **Complexity cost:** removing it means pipelining the producer's per-chunk shaping against the channel drain (a dedicated producer thread, or splitting "shape" from "ship" so chunk N+1's shaping starts before chunk N's `send` returns) — a topology change. High cost; byte-identity holds as long as `chunk_order` stamping stays monotonic (the writer already reorders). **Do not attempt before H1** — H1 accounts for the dominant 64 % and may move the floor enough that this churn isn't worth it.

**L2: [src/var_calling/posterior_engine.rs:2693-2719](../../../src/var_calling/posterior_engine.rs#L2693-L2719) — `e_step_simd` lays SIMD lanes across the sample axis via 4-way gather; for biallelic-diploid (`n_genotypes == 3`) the SIMD inner trip count is 3 and every lane-load is a strided gather**
- **Confidence:** Medium
- **Hot-path evidence:** profile leaf `e_step_simd 2578` (T=1) / `2492` (T=8) — a top project leaf, the named hot region after the prior QUAL convolution went cold. `cargo asm` confirms the `wide::f64x4` ln/exp **does** lower to NEON (260 `*.2d` ops) — so this is a loop-*shape* issue, not a vectorization failure: the batch is built across the sample axis (stride `n_genotypes`) over row-major sample-major data, so each `f64x4` is four non-contiguous loads + a `to_array()` scatter back.
- **Pattern matched:** autovectorization over non-contiguous data — SIMD laid across the wrong axis for this shape; the contiguous genotype axis (length 3) is scalar-iterated.
- **Mechanism:** gather/scatter cost more than contiguous `ldr q` and the `to_array()` round-trip materializes lanes to the stack, amortized over only a 3-trip inner loop. A `(2,3)` scalar fast path (mirroring the existing `accumulate_expected_counts` `(2,3)` precedent) may beat the SIMD batch outright for the dominant shape; keep SIMD for wide-`n_genotypes` (multiallelic/polyploid).
- **Measurement plan:** `examples/profile_posterior_engine.rs` under `sample` at the biallelic-diploid shape; benchmark a `(2,3)` scalar specialization vs the gather batch. Threshold ≥ 5 % drop in `e_step_simd` self-time; confirm with `cargo asm` that gathers become contiguous loads.
- **Complexity cost:** a shape-specialized inner loop (one extra path) or a lane-axis transpose. **Risk: feeds `posteriors` → GT/GQ/AD — must stay byte-identical**; preserve the exact `exp(log_post_unnorm - log_z)` per-cell evaluation order; gate on the byte-identity harness. **This is a restructure experiment, not a guaranteed win** (the gain is bounded by gather/scatter + short-trip overhead, not by adding SIMD).

**L3: [src/var_calling/sample_reader.rs:472-513](../../../src/var_calling/sample_reader.rs#L472-L513) (`records_all`/`records_for`) — per-allele `Vec<u8>` + `Vec<ChainId>` re-allocation when the caller rebuilds records**
- **Confidence:** Medium
- **Hot-path evidence:** DHAT `records_all` 262 MB (19 %); record-build memmove inside the ~7.3k bucket; allocator `_xzm_xzone_malloc_tiny 1268` / `_xzm_free 1438`.
- **Pattern matched:** per-allele `slice_at(a).to_vec()` allocates one fresh `Vec` per allele per record per chunk on the caller thread.
- **Mechanism:** `records_all` hands the byte-identity-critical `PerPositionMerger` the owned `PileupRecord` shape it consumes, so the owned-`Vec` allocations are inherent to that boundary. **The largest removable slice is the same REF chain-id dead-weight (H2) resurfacing at materialize time** — H2 step 1 collapses it for free (REF slot → allocation-free `Vec::new()`). The residual (ALT chain_ids + all seqs) is genuine merger input, not removable without changing the merger's ownership contract.
- **Measurement plan:** same DHAT command; isolate the `records_all` line before/after H2 — expect it to drop ~proportionally to the REF chain-id share, no separate change needed. Only if a residual stays top-5, consider a `bumpalo` arena over a per-chunk record batch.
- **Complexity cost:** zero for the bundled effect (rides H2). A standalone arena adds a `bumpalo` dependency + a merger-input lifetime — not worth it unless DHAT still flags `records_all` after H2.

**L4: [src/var_calling/cohort_integration.rs:992](../../../src/var_calling/cohort_integration.rs#L992) + [src/fasta/fetcher.rs:697](../../../src/fasta/fetcher.rs#L697) — REF span fetched into a fresh `Vec::with_capacity` per chunk instead of a reused producer-side buffer**
- **Confidence:** Medium
- **Hot-path evidence:** profile `read_uppercased_bases 256` (T=1) / 203 (T=8 producer breakdown) — on the single producer thread the T=8 profile shows is the wall floor.
- **Pattern matched:** per-chunk owned-`Vec` on a streaming read path where a scratch buffer would do (project's load/use/clear/reload scratch pattern).
- **Mechanism:** `fetch_ref_span` → `StreamingChromRefFetcher::fetch` does `Vec::with_capacity(len)` per produced chunk; at `tvpc=256` that's one heap alloc + free per chunk on the producer. The `ChromRefFetcher` trait already exposes `fetch_into(start, len, &mut dst)`, so the producer could own one `ref_scratch: Vec<u8>` and `clear()`+`fetch_into` per chunk. The `read_uppercased_bases` self-time is the byte-strip/uppercase copy (intrinsic); the *allocation* is the removable part.
- **Measurement plan:** DHAT at N=50/T=1 before/after — producer-thread allocator self-time drops, wall at T≥2 improves, **zero VCF diff** (same `fetch_into` math). Threshold: any reproducible T=8 wall improvement with zero diff.
- **Complexity cost:** modest signature churn — `fetch_ref`'s closure (`pipeline.rs:286`) and `fetch_ref_span` change from returning/owning `Vec<u8>` to filling a producer-owned `&mut Vec<u8>`; `RefSpan` borrows or holds a producer slice. No `unsafe`. Honest caveat: the per-chunk `Vec` may be cheap enough at 256-variant granularity that the win is in the noise — hence the DHAT gate.

**L5: [src/var_calling/per_position_merger.rs:38-46](../../../src/var_calling/per_position_merger.rs#L38-L46),[306-339](../../../src/var_calling/per_position_merger.rs#L306-L339) + [types.rs:37-45](../../../src/var_calling/types.rs#L37-L45) — `per_sample: Vec<Option<PileupRecord>>` is an N-wide AoS the merge + grouper scan one field at a time; the `None` slots are dead cache lines**
- **Confidence:** Medium
- **Hot-path evidence:** pattern-match for the cache angle — the profile names EM and decode, not the per-position merge directly. At N=50 every emitted `PerPositionPileups` is a 50-slot vector; `variant_grouping::has_variant_observation` / `max_ref_span` scan all 50 slots per position reading only `r.alleles.len()` / `r.ref_span()`, chasing each `Option<PileupRecord>` into its heap `alleles` `Vec`, most `None` at a typical partially-covered variable position.
- **Pattern matched:** AoS-vs-SoA by access pattern — the grouper's seed-time scan reads exactly one derived field per slot but must pointer-chase to get there.
- **Mechanism:** projecting the grouper's working set (per-slot `is_variant: bool` + `ref_span: u32`) into a small parallel array refreshed once per emitted position would collapse the seed/extend scan to a contiguous sweep.
- **Measurement plan:** `sample` (or cachegrind) over `profile_cohort_e2e` at N=50; look at miss rate / self-time in `VariantGrouper::next` + `has_variant_observation`. Threshold: measurable miss-rate or self-time drop. Byte-identity: `variant_grouping.rs` tests + the cohort streaming-vs-reference equality test pass.
- **Complexity cost:** a derived cache array recomputed per `PerPositionPileups` + the sync invariant. Modest, no `unsafe`. Borderline `hot_loops`/restructuring (the merge path is not in the quoted leaves), hence Likely not Hot-path. Possibly better addressed in the producer by emitting columnar fold inputs rather than rebuilding `Vec<Option<PileupRecord>>`.

**L6: [src/var_calling/posterior_engine.rs:2761](../../../src/var_calling/posterior_engine.rs#L2761) — per-cell `is_finite()` branch inside the `e_step_simd` scalar tail loop**
- **Confidence:** Medium
- **Hot-path evidence:** the `e_step_simd` tail runs the scalar body for the `n_samples % 4` leftover (2 of every 50 at N=50); the per-cell `if !posterior.is_finite() { return Err(...) }` is a data-dependent branch inside the otherwise straight-line `exp` loop.
- **Pattern matched:** a tight numeric loop with an `if` that could be hoisted — accumulate an `is_finite` AND across the row, check once after the loop, store unconditionally.
- **Mechanism:** the per-element early-return prevents a clean straight-line `exp` + store; the error path is genuinely cold (a non-finite posterior is a malformed-record bug), so it belongs behind a `#[cold]` re-scan that locates the bad `(sample, genotype)` only when the fast AND fails.
- **Measurement plan:** `profile_posterior_engine` self-time on `e_step_simd` tail + scalar `e_step`. Worth it mainly **bundled with L2** (the E-step reshape). Byte-identity safe (reorders the *check*, not the math; preserves the exact error variant).
- **Complexity cost:** a hoisted `all_finite` accumulator + a `#[cold]` slow-path re-scan. Low.

**L7: [src/var_calling/cohort_integration.rs:516-520](../../../src/var_calling/cohort_integration.rs#L516-L520),[954-962](../../../src/var_calling/cohort_integration.rs#L954-L962) — per-row `kept_all.binary_search(&p)` builds the keep mask in O(n log m); a merge-walk against the ascending `kept_all` is O(n+m)**
- **Confidence:** Medium (Low hot-path; flagged by the concurrency agent as a cross-category route)
- **Hot-path evidence:** pattern-match only — inside `compact_samples` (rayon, per sample); not isolated as a leaf.
- **Pattern matched:** binary-search-in-a-loop against an already-sorted target where a linear merge-walk suffices.
- **Mechanism:** `kept_all` is ascending and the rows are walked in order, so a two-pointer merge replaces `n` binary searches with one `O(n+m)` sweep.
- **Measurement plan:** fold into the H2 / compaction measurement — DHAT/`sample` on `compact_samples` self-time. Threshold ≥ 5 % of `compact_samples` self-time.
- **Complexity cost:** a two-pointer rewrite of the mask build; small. Byte-identical (same mask).

### Speculative

**S1: [src/psp/reader.rs:1483](../../../src/psp/reader.rs#L1483) / [src/psp/block.rs:709](../../../src/psp/block.rs#L709),[721](../../../src/psp/block.rs#L721) — zstd content checksum (XXH64) is verified on every column frame and is not skippable from the in-scope reader; compression level 9 is tuned for write-side ratio**
- **Confidence:** Low
- **Hot-path evidence:** `XXH_INLINE_XXH64_update 2202` inside the ~11.2k zstd-decode bucket — pure verification over already-validated, self-produced `.psp` bytes (the block index already carries an XXH3-64 checksum, `reader.rs:186`).
- **Pattern matched:** redundant integrity work on trusted compressed input.
- **Mechanism:** the writer sets `include_checksum(true)` (`block.rs:721`); the reader verifies a trailing XXH64 per frame. But `zstd` 0.13's `bulk::Decompressor` exposes **no** `ZSTD_d_forceIgnoreChecksum` — skipping verification on the read side is unreachable without dropping to `zstd-safe`'s lower-level `DCtx`. The only API-reachable lever is the **write side** (`include_checksum(false)`, or a lower level), which **changes the `.psp` bytes on disk** — a wire-format/robustness change owned by the out-of-scope `pileup`/`psp-to-pileup` subcommands. Skipping the checksum does **not** change VCF output, but trades intermediate-file bit-rot/crash detection.
- **Measurement plan:** two write-side experiments, each gated on a full byte-identical VCF diff: (1) re-emit the corpus with `include_checksum(false)`, re-profile — does the XXH64 line disappear and wall drop? (2) sweep `ZSTD_COMPRESSION_LEVEL` down (6, 3) — decode self-time vs `.psp` size vs wall.
- **Complexity cost:** a constant flip, but a `.psp` **format-policy** decision (lose an integrity guard on a multi-GB on-disk intermediate). The read-side win is not reachable from in-scope code at all (API gap), so the lever lives entirely in the out-of-scope writer subcommands.

**S2: [src/var_calling/posterior_engine.rs:3171-3180](../../../src/var_calling/posterior_engine.rs#L3171-L3180) — the QUAL convolution FIR does not vectorize on apple-m1 despite a doc claiming "autovectorises to FMA" (site is now cold; QUAL is exempt)**
- **Confidence:** Low
- **Hot-path evidence:** `cargo asm` of `convolve_ac_linear` shows **zero** NEON FMA and zero vector-width loads (the axpy lowered to scalar, likely blocked by slice aliasing across `&mut` reborrows or the variable trip count). But the prior review's H1 took the QUAL convolution cold — `compute_qual_via_exact_af` is 869 samples and the convolution is a fraction.
- **Pattern matched:** confirm-autovec-with-asm — it isn't happening, cause is plausibly aliasing/variable-trip (fixable), but the payoff is small because the site is cold.
- **Mechanism:** scalar axpy leaves FMA units idle on the O(n²·ploidy) fold; restoring vectorization (assert non-alias, `chunks_exact` a fixed live window) would FMA it — only matters if QUAL becomes hot again at large N.
- **Measurement plan:** `profile_cohort_e2e` at N=200 under `sample`; act only if `compute_qual_via_exact_af`/`convolve_ac_linear` re-enter the top leaves. **Also fix or remove the stale "autovectorises to FMA" doc claim** (`:3042`/`:3107`/`:3167`) — it is currently false on apple-m1.
- **Complexity cost:** low (loop-shape change). **QUAL-only — float noise accepted**, so this is the one finding free to diverge numerically.

**S3: [src/var_calling/cohort_integration.rs:516-520](../../../src/var_calling/cohort_integration.rs#L516-L520),[862-870](../../../src/var_calling/cohort_integration.rs#L862-L870),[954-962](../../../src/var_calling/cohort_integration.rs#L954-L962) — per-chunk/per-segment `Vec<bool>` masks + `kept_all`/`variable` allocations built fresh each call**
- **Confidence:** Low
- **Hot-path evidence:** pattern-match only — not in the DHAT top attributions; `find_cut` is 8 samples. These are O(positions-in-fold) bool/`u32` vecs dwarfed by the columnar payload.
- **Pattern matched:** allocations outside hot loops — `keep`/`kept_all`/`variable` allocate per chunk (and `keep` per straddling segment). The fold double-buffers and `is_kept` are already `.clear()`-reused (the project's scratch pattern is applied to the fold).
- **Mechanism:** the remaining vecs are smaller, and some live inside the `compact_samples` `par_iter_mut` closure, so hoisting them onto `self` forces per-rayon-task scratch (`thread_local!` or a per-task struct) — trading churn for structure.
- **Measurement plan:** only worth it if a post-H2 DHAT re-attribution promotes them. Threshold: churn drop ≥ ~5 %.
- **Complexity cost:** per-rayon-task scratch ownership. Not worth it without a measurement promoting them — likely the 3–8 % "channel-bound" tier already judged low-value at N=50.

**S4: [src/var_calling/cohort_integration.rs:465-468](../../../src/var_calling/cohort_integration.rs#L465-L468) — `BufferedSegment::Ready(Box<ReadySegment>)` is a read-side pointer chase in the light-column scans**
- **Confidence:** Low
- **Hot-path evidence:** pattern-match only; not in the profile leaves. `coverage`/`watermark`/`rebuild_fold`/`drop_buffers_below` call `positions()`, dereferencing the `Box` for `Ready`.
- **Pattern matched:** pointer chase — but the boxing motive (keep the transient `Pending` variant small so the `Vec` stays dense) is correct; un-boxing would bloat every `Pending` slot.
- **Mechanism:** the scans only need `positions().first()/.last()`; storing `(first_pos, last_pos)` inline on `BufferedSegment` (both variants) would answer without dereferencing either `positions` `Vec`.
- **Measurement plan:** only if a profile attributes self-time to the buffer scans (current evidence does not). Defer.
- **Complexity cost:** two `u32` fields duplicated onto the enum + the sync invariant. Low but unmotivated by current evidence.

### Note

- **[src/psp/block.rs:451](../../../src/psp/block.rs#L451) `decode_list_column_csr` is codegen-clean** despite being the #2 leaf (3303 samples): `cargo asm --simplify` shows 0 `panic_bounds_check`/`slice_index_fail`/alloc in the hot body (the up-front `reserve` + DoS-guard let LLVM elide per-element bounds checks; the `to_string()` allocs are cold error closures). Its weight is the intrinsic chain-id `u64` payload — the lever is **not decoding it** (H2), not a loop fix.
- **[src/psp/varint.rs:82](../../../src/psp/varint.rs#L82) `decode_u64_leb128`** is already optimally split (`#[inline]` single-byte fast path + `#[cold] #[inline(never)]` multi-byte tail). No finding.
- **[src/var_calling/per_position_merger.rs:306](../../../src/var_calling/per_position_merger.rs#L306) `vec![None; n_samples()]`** per emitted position is genuine merger output (moved into `PerPositionPileups`), byte-identity-locked. A future redesign could slab/arena the slot-vectors; not actionable now without reshaping the merge contract.
- **[src/psp/reader.rs:1617](../../../src/psp/reader.rs#L1617) `read_compressed_blob`** allocates a fresh `Vec` per retained column — **intentional** (the blob must outlive the decode call until `set_variable_rows`). Flagged so it is not "optimized" into the per-column scratch, breaking the retain-then-inflate contract.
- **`buildFSETable`/HUF table reconstruction (529 + ~700 samples)** is intrinsic per-frame entropy-table rebuild under one-frame-per-column framing — the decompressor *context* is reused (confirmed), but tables are per-frame by design. Not a context-reuse miss; the only amortization is fewer/larger frames (a wire-format redesign, out of scope).
- **`em_posterior_calc::call_records` (3361, top leaf)** self-time is the EM driver + `merge_compacted_samples` record reconstruction + per-record `.clone()`s out of scratch at `posterior_engine.rs:2366-2370` — the clones route to `allocations` (subsumed by H2/L3); no hot_loops defect in the driver loop.

## 6. Out-of-scope observations

- **The maintained profiling driver profiles the wrong shape.** [examples/profile_cohort_e2e.rs:161](../../../../examples/profile_cohort_e2e.rs#L161) hard-codes `target_variants_per_chunk: 0` (legacy single-pull) and its doc (`:12-13`,`:42`) references the deleted `cohort_e2e_perf` bench + a stale `target-container` path. Fix the default to `256` (or add the flag the sibling `dhat_var_calling.rs` already has) and correct the doc when the new bench lands. Folded into the measurement plan (§3 items 1–2).
- **`.psp` write-side knobs** (`include_checksum`, compression level) that could cut decode-side cost live in the `pileup`/`psp-to-pileup` subcommands (out of scope) — see S1. Any change there is a format-policy decision for the format owners, gated on a byte-identical VCF.

## 7. What's already good

- **Persistent zstd `Decompressor` + reused scratch** — one `zstd::bulk::Decompressor<'static>` per reader, threaded through every column of every block ([src/psp/reader.rs:599](../../../src/psp/reader.rs#L599)); the DCtx workspace is not rebuilt per block. (Prior L1, confirmed still applied.)
- **BufReader seek short-circuit** — `seek_to_offset` ([src/psp/reader.rs:1851](../../../src/psp/reader.rs#L1851)) uses `SeekFrom::Current(delta)` with an early-return at target, never `SeekFrom::Start`, so the 64 KiB buffer survives contiguous block transitions (matching the tiny `read` self-time). (Prior PSP-reader H2, confirmed applied.)
- **`MergedRecord` flat row-major tables** ([src/var_calling/per_group_merger.rs:313-354](../../../src/var_calling/per_group_merger.rs#L313-L354)) collapse `n_samples + 1` small per-matrix heap allocations into one flat `Vec` with `*_row` accessors — the right shape for an N=50 cohort; do not regress to nested vecs.
- **Sound `thread::scope` topology with no hot-path locks** ([src/var_calling/pipeline.rs:234-259](../../../src/var_calling/pipeline.rs#L234-L259)) — bounded crossbeam channels (no unbounded-buffer OOM), `&caller` shared immutably (`Sync`, stateless `call_chunk`), channel handles dropped to close cleanly, first-error-wins join; no `Mutex`/`RwLock`/atomic anywhere in scope (the contention is scheduler oversubscription, H1, not locks).
- **VCF write is correctly buffered** ([src/vcf/sink.rs:72](../../../src/vcf/sink.rs#L72)) — 64 KiB `BufWriter`, explicit flush with `IntoInnerError` propagation, single end-of-stage fsync (file + parent dir), no per-record syscall.
- **The SIMD EM E-step genuinely vectorizes** — `cargo asm` confirms `wide::f64x4` ln/exp lowers to 260 NEON `*.2d` ops (the remaining L2 finding is loop *shape*, not a vectorization failure).

### Author response convention
Address each finding by its identifier (H1, L2, S1, …) with one of: `applied in <commit>` / `experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. The "experiment shows no gain" path is expected and welcome — that is what the measurement plan is for.
