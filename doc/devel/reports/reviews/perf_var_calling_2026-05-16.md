# Performance Review: var_calling
**Date:** 2026-05-16
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** all of `src/var_calling/` (Stages 4–5 of the multi-stage SNP caller) + the criterion harness at `benches/var_calling_perf.rs`
**Verdict:** Apply the listed wins
**Hot-path evidence:** five samply sampling profiles + fresh criterion baseline + one allocator-switch A/B experiment

---

## 1. Scope and constraints

- **What was reviewed:** the variant-calling pipeline, Stages 4 and 5:
  - [src/var_calling/mod.rs](../../../src/var_calling/mod.rs) (module index)
  - [src/var_calling/per_position_merger.rs](../../../src/var_calling/per_position_merger.rs) — Stage 4 k-way per-position merger
  - [src/var_calling/variant_grouping.rs](../../../src/var_calling/variant_grouping.rs) — Stage 4 overlap bundler
  - [src/var_calling/per_group_merger.rs](../../../src/var_calling/per_group_merger.rs) — Stage 5 allele unification + closed-form likelihood (rayon-parallel)
  - [benches/var_calling_perf.rs](../../../benches/var_calling_perf.rs) — the criterion harness driving every profile in this review
- **Reviewed against:** commit `387424c` on branch `main`. Build is `[profile.bench]` (inherits release + `lto = "fat"` + `codegen-units = 1` + `panic = "abort"` + `target-cpu = x86-64-v3`, with `debug = true`).
- **Performance intent and targets:** Stage 4 is a streaming serial pipeline (k-way per-position merge → overlap grouping). Stage 5 is the per-group merger, rayon-parallelised across groups. Target deployment: multi-core Linux (the dev host has 16 hyperthreads). Workload: WGS cohorts of N=10–1000 samples × 50–300 Mbp chromosomes. Throughput goal: practical cohort calls complete in minutes-to-hours, not days. Stage 5 dominates wall-clock per the criterion baseline (~22 K groups/s biallelic vs ~5.5 M positions/s grouper).
- **Correctness contract:** every per-sample × per-genotype log-likelihood must remain a finite log; numerical regressions are unacceptable. Every code-level finding below preserves bit-exact output.
- **Hot-path evidence available:** all evidence files live under [tmp/perf_review_2026-05-16_var_calling/](../../../tmp/perf_review_2026-05-16_var_calling/):
  - Fresh criterion baseline ([baseline_criterion.txt](../../../tmp/perf_review_2026-05-16_var_calling/baseline_criterion.txt)):
    ```
    var_calling_merger/dense_64_samples_200000_positions     904.11 ms   221 K elem/s
    var_calling_merger/sparse_64_samples_200000_positions    565.23 ms   354 K elem/s
    var_calling_grouper/snp_dense_16_samples_200000_pos       38.09 ms   5.25 M elem/s
    var_calling_grouper/overlap_extension_16_samples_5000_groups   3.86 ms  1.30 M elem/s
    var_calling_per_group_merger/biallelic_snp_64_samples_10000_groups        409.52 ms   24.4 K elem/s
    var_calling_per_group_merger/compound_all_anchored_64_samples_10000_groups   794.53 ms   12.6 K elem/s
    var_calling_per_group_merger/compound_half_chain_broken_64_samples_10000_groups   805.15 ms 12.4 K elem/s
    var_calling_per_group_merger/biallelic_snp_1_sample_10000_groups            29.63 ms  338 K elem/s
    ```
  - Five samply sampling profiles (3 Stage 5 + Stage 4 merger + grouper). Self-time and inclusive-time rollups in [all_symbolicated.txt](../../../tmp/perf_review_2026-05-16_var_calling/all_symbolicated.txt); the synthesised diagnosis (with libc symbols resolved via the system `libc6-dbg`) in [hot_paths_summary.md](../../../tmp/perf_review_2026-05-16_var_calling/hot_paths_summary.md).
  - **Headline diagnosis:** ~75 % of Stage 5 CPU is in glibc malloc — 30.8 % `futex_wait` + 23.3 % `__lll_lock_wake_private` + 6.6 % `_int_free_chunk` self-time. Call-chain attribution ([callers.py](../../../tmp/perf_review_2026-05-16_var_calling/callers.py)) says **93.4 %** of the futex_wait stacks bubble up from `drop_in_place<OverlappingVariantGroup>` — 16 rayon workers cross-thread free-ing the nested `Vec<Option<PileupRecord>>` → `Vec<AlleleObservation>` → `Vec<u8>` + `Vec<ChainId>` chain on the glibc arena lock.
  - **Allocator A/B:** mimalloc tried via the existing `alloc-mimalloc` feature ([mimalloc_full_biallelic.txt](../../../tmp/perf_review_2026-05-16_var_calling/mimalloc_full_biallelic.txt)). Result: **+47.7 % slower** on `biallelic_snp_64_samples_10000_groups` (604.90 ms vs 409.52 ms glibc); 3–9 % slower on compound paths. Switching allocators is not the fix here.
  - **Methodology caveat:** the Stage 4 merger and grouper profiles are **fixture-dominated** — only 13–21 % of their inclusive time is the function under test; the rest is `iter_batched` setup rebuilding 64×200 000 fixtures every iteration. All Stage 4 code-level findings below are tagged Likely / Speculative at best until the bench shape is fixed.
- **Deliberately out of scope:** Stage 6 (posterior engine, not implemented yet); Stage 3 (DUST filter, not implemented yet); legacy gVCF-merger code (`src/gvcf_parser.rs` and siblings — out by PROJECT_STATUS); the I/O code paths in `psp_reader` / `psp_writer` / CRAM reader (each has its own perf review).
- **Categories dispatched:** methodology (always; large bench/build-config surface here); allocations (75 % of CPU in malloc/free); data_layout (`Vec<Vec<…>>` nesting throughout Stage 5); concurrency (rayon par_iter, untuned batch size, cross-thread free pattern); hot_loops (per-genotype × per-sample inner loops in Stage 5). **Skipped:** io_and_syscalls — `src/var_calling/` performs no I/O directly; the `SharedRefFetcher` trait is pluggable and the benches mock it in memory.

---

## 2. Verdict

**Apply the listed wins.**

The dominant cost is well-evidenced and has a clear mechanism: cross-thread free of nested `Vec<u8>` allocations on the glibc arena lock. The mimalloc A/B experiment ruled out the easy build-config fix; the remaining levers are reducing per-group allocation count (Hot-path findings H1–H4 below) and tuning the rayon split granularity (H5). The methodology fixes (H6–H7) unblock measurement of the rest.

Each Hot-path finding ships independently and is measurable on the existing criterion harness (with the bench fixes from H6 / H7 for the Stage 4 ones). Run them as separate PRs — one hypothesis per measurement — and re-profile after each lands.

---

## 3. Measurement plan

The measurements break into three waves, ordered so each wave unblocks the next.

**Wave 1 — code-level changes already measurable on the current bench.** All three are Stage 5 internal changes that the existing `var_calling_per_group_merger/*` benches exercise. Each lands as one PR with a single hypothesis.

  1. Run `cargo bench --bench var_calling_perf -- var_calling_per_group_merger --save-baseline pre_fix_2026-05-16` once on `main`. This is the cross-commit baseline every Wave-1 PR compares against (the in-tree [baseline_criterion.txt](../../../tmp/perf_review_2026-05-16_var_calling/baseline_criterion.txt) is the same run captured for the review).
  2. For each Wave-1 PR: rebuild the bench, run `cargo bench --bench var_calling_perf -- var_calling_per_group_merger --baseline pre_fix_2026-05-16`. **Reject** any PR that regresses any of the three Stage 5 shapes (biallelic, compound_all_anchored, compound_half_chain_broken) by more than criterion's noise floor (~1–2 %) at p < 0.05. **Merge** any PR that shows ≥ 3 % improvement on at least one shape with no regression elsewhere. Per the methodology, do *not* bundle multiple Wave-1 changes into one PR — pick the order from H1 → H2 → H3 (described in §5).
  3. After every Wave-1 PR lands, re-capture a samply profile of `biallelic_snp_64_samples_10000_groups` (15 s `--profile-time`) and confirm the inclusive share on `drop_in_place<OverlappingVariantGroup>` is decreasing. The reproduction recipe is in [hot_paths_summary.md](../../../tmp/perf_review_2026-05-16_var_calling/hot_paths_summary.md).

**Wave 2 — bench-shape and config-tuning experiments.** Sequenced after Wave 1 because they are infrastructure for the Stage 4 and Stage 5 follow-ups.

  4. **H6 fixture-rebuild fix** — rewrite the Stage 4 merger and grouper benches to build their per-sample-streams *once* outside `b.iter`. Acceptance: a fresh samply profile of `var_calling_merger/dense_64_samples_200000_positions` shows `PerPositionMerger::next` ≥ 50 % inclusive (currently 21 %). Once this lands, all Stage 4 Likely/Speculative findings (L5, L7, S2–S4) become rankable against a clean profile.
  5. **H5 batch-size sweep** — add a parametric bench varying `PerGroupMergerConfig::batch_size ∈ {32, 128, 512, 2048}` for each of the three Stage 5 shapes. Capture medians, pick the lowest-cost batch_size that does not regress any shape, and reset `DEFAULT_BATCH_SIZE` (currently a self-declared "placeholder"). Save the result as a criterion baseline named `batch_sweep_baseline`.
  6. **H7 cohort-size sweep** — add bench variants at `N_SAMPLES ∈ {10, 64, 256, 1024}` (reduce `N_GROUPS` proportionally so wall time stays bounded) for `biallelic_snp` and `compound_all_anchored`. The output is wall-time per element across N; super-linear scaling at N=1024 vs N=64 names the next code-level priority.

**Wave 3 — Speculative experiments, only if Wave-1+2 doesn't close the gap.**

  7. PGO (S5) — only worth it once Stage 5 is back to CPU-bound rather than malloc-bound; do not run earlier.
  8. Streaming refill pipeline (S6, replaces the `.collect::<Vec<_>>()` barrier) — first add a per-batch min/max latency probe so the gap to close is measurable; only escalate to the rewrite if the gap is ≥ 20 % of median per-group time.
  9. `RAYON_NUM_THREADS` + CPU pinning (S7) — `taskset -c 0-15` + `RAYON_NUM_THREADS=16` vs unpinned; only document a non-default policy if it shows ≥ 10 % wall-time improvement.

---

## 4. Build / toolchain configuration

[Cargo.toml](../../../Cargo.toml) is already aggressive: `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`, `opt-level = 3` implicit, and `[profile.bench]` inherits release with `debug = true` for samply. [`.cargo/config.toml`](../../../.cargo/config.toml) pins `target-cpu = "x86-64-v3"`. `rust-toolchain.toml` pins 1.95 for reproducibility. The build-config defaults are not leaving wins on the floor.

Two configuration findings remain:

- **L1 [build] — `alloc-mimalloc` feature should be downgraded from "recommended A/B" to "documented negative result; off by default"** ([Cargo.toml:74-80](../../../Cargo.toml#L74-L80)). The current comment frames mimalloc as a generic perf knob, but the A/B run on this exact workload was +47.7 % slower; the next reader will retry it unless the negative result is in the comment. Documentation change only; the feature stays available for re-runs when allocation count drops. Full text proposed in [methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Likely §1.

- **S5 [build] — PGO is unmeasured.** All prerequisites for "everything else has been done" are set, so PGO is the next mechanical knob. Defer until after Wave 1 lands — PGO over a workload where 70 % of CPU is futex-wait is wasted work; once Stage 5 is CPU-bound rather than malloc-bound, run the `cargo pgo build` / `cargo pgo bench` / `cargo pgo optimize build` cycle and re-measure. Threshold: ≥ 3 % median improvement on the Stage 5 benches.

The defensive [build] note at [Cargo.toml:34](../../../Cargo.toml#L34) (release uses `debug = "line-tables-only"`, bench overrides to `debug = true`) is *already correct* — methodology sub-agent's `[Speculative]` filing was a guard against a future reviewer proposing to strip the bench debug info. No action.

---

## 5. Code-level findings

Severity codes are assigned here. Each finding has the full body in its sub-agent file; the synthesis below is the orchestrator's summary plus the order-of-application reasoning.

### Hot-path

#### H1: `src/var_calling/per_group_merger.rs:716-764` — Replace `BTreeMap<Vec<u8>, usize>` byte→allele index with `ahash::AHashMap`, and drop the second `Vec<u8>` key clone

- **Confidence:** High
- **Hot-path evidence:** profile attribution (1.0–1.7 % self-time at `project_scalars:1138` in [all_symbolicated.txt](../../../tmp/perf_review_2026-05-16_var_calling/all_symbolicated.txt)); every `BTreeMap` insert clones the projected `Vec<u8>` key ([per_group_merger.rs:743](../../../src/var_calling/per_group_merger.rs#L743), [:798](../../../src/var_calling/per_group_merger.rs#L798)) — those `Vec<u8>` clones live as long as the group and feed straight into the 71 % `drop_in_place<OverlappingVariantGroup>` inclusive cost.
- **Why first:** smallest diff (one type swap + one `use` + an optional key-borrow restructure), no API change, `ahash` is already a project dep. Even the mechanical first step (`BTreeMap<Vec<u8>, usize>` → `AHashMap<Vec<u8>, usize>`) removes the `O(log n)` tree-walk-per-insert and the SipHash overhead. The key-borrow restructure (which removes the *second* `Vec<u8>` clone entirely) is a follow-on commit if step 1 is on the profile.
- **Fix:** full diff in [allocations.md](../../../tmp/perf_review_2026-05-16_var_calling/allocations.md) Hot-path §1.

#### H2: `src/var_calling/per_group_merger.rs:1131-1133, :1100, :1183, :1421, :1464, :724-749` — Flatten per-group `Vec<Vec<T>>` matrices to row-major `Vec<T>` with a stride

- **Confidence:** High
- **Hot-path evidence:** `project_scalars:1138` directly named at 1.0 % self-time as "the `let mut scalars: Vec<Vec<AlleleSupportStats>>` allocation"; `compute_log_likelihoods:1490` at 1.9 % covers the same construction shape for `log_likelihoods` (`Vec<Vec<f64>>`). For 64 samples × 10 000 groups, every Stage 5 group pays `n_samples + 1` heap allocations per matrix × ~5 matrices (`scalars`, `log_likelihoods`, `chain_anchor_flags`, `enforce_max_alleles` pool, `compute_log_likelihoods` outer Vec). Drop runs on the worker thread and feeds the 71 % `drop_in_place<OverlappingVariantGroup>` figure.
- **Why second:** removes the single largest *count* of per-group allocations, but touches the `pub MergedRecord` API on `scalars`, `chain_anchor_flags`, `log_likelihoods`. The blast radius is currently small — no Stage 6 consumer exists, only tests + benches — so this is the right moment to ship the API change. Add `pub fn log_likelihood(sample, genotype) -> f64`, `pub fn scalars_row(sample) -> &[AlleleSupportStats]`, etc., so future Stage 6 work uses accessors rather than indexing the flat buffer directly. The internal-only `UnifiedAllele.per_sample_sources: Vec<Vec<(usize, usize)>>` ships in the same PR (it's `pub(crate)`-only; no API exposure) — flatten to either one `Vec<(u8, u32, u32)>` per allele (step a), or CSR `(offsets, items)` across the whole UnifiedAlleleSet (step b, only if step a is still on the profile).
- **Fix:** sketches in [allocations.md](../../../tmp/perf_review_2026-05-16_var_calling/allocations.md) Hot-path §2 and [data_layout.md](../../../tmp/perf_review_2026-05-16_var_calling/data_layout.md) Hot-path §1–§3.

#### H3: `src/var_calling/per_group_merger.rs:1620-1633` — Replace `chain_broken_log_likelihood`'s per-constituent `Vec<u32>` / `Vec<f64>` / `Vec<u32>` allocations with fixed-size `[T; MAX_BITMASK_ALLELES]` stack arrays

- **Confidence:** High
- **Hot-path evidence:** `compute_log_likelihoods:1481` at 1.4–3.4 % self-time on the compound profiles ([all_symbolicated.txt](../../../tmp/perf_review_2026-05-16_var_calling/all_symbolicated.txt)). The call is the inner body of a `(sample × genotype × constituent)` triple loop — three `Vec` allocations per cell, dropped on the worker thread. The chain-broken bench shape is the workload's most allocator-bound Stage 5 shape (805 ms vs 410 ms biallelic).
- **Why third:** Localised change inside one function, bit-exact (no FP-order change), directly analogous to the existing `p_counts: [u32; MAX_BITMASK_ALLELES]` pattern at [per_group_merger.rs:1546](../../../src/var_calling/per_group_merger.rs#L1546) — keeps the module's convention consistent. The most direct rewrite reads `rec.alleles[a].support.num_obs / q_sum` *inline* from the per-position record and only retains `per_pos_counts` as scratch; that removes all three allocations rather than just one. Cuts compound_half_broken's per-genotype allocation count to one stack-resident array.
- **Fix:** full replacement snippet in [hot_loops.md](../../../tmp/perf_review_2026-05-16_var_calling/hot_loops.md) Likely §3 / [allocations.md](../../../tmp/perf_review_2026-05-16_var_calling/allocations.md) Hot-path §3 (the two sub-agents converged on the same site).

#### H4: `src/var_calling/per_group_merger.rs:299-346, :1463` — Cache the `genotype_order(ploidy, n_alleles)` table on the merger; do not rebuild `Vec<Vec<u8>>` every group

- **Confidence:** High
- **Hot-path evidence:** `genotype_order` allocates ~`n_genotypes + 1` heap objects per group (3 for biallelic, 21 for default-cap multiallelic). For 10 000 groups × 22 allocations = 220 K allocations per bench iter — on the allocator-bound critical path. `ploidy` is fixed by the config; `n_alleles` is bounded by `max_alleles = 6` (default), so the realistic envelope is a handful of distinct tables.
- **Why fourth:** Smaller per-group allocation share than H1–H3, but the simplest fix in the list — build a `HashMap<n_alleles, Arc<Vec<Vec<u8>>>>` lazily on `PerGroupMerger` and look up by `n_alleles` per `process_group`. Bit-exact (same sort, same enumeration). Pairs naturally with H2 (both touch `compute_log_likelihoods`).
- **Fix:** snippet in [hot_loops.md](../../../tmp/perf_review_2026-05-16_var_calling/hot_loops.md) Likely §4.

#### H5: `src/var_calling/per_group_merger.rs:48-55` — `DEFAULT_BATCH_SIZE = 32` is a self-declared placeholder; sweep and re-set

- **Confidence:** High
- **Hot-path evidence:** the source comment is the evidence: *"**Placeholder** — chosen as a reasonable starting point for cache-friendly batches; no comparative benchmark has selected this specific value yet."* With 16 worker threads × 32-group batches, each worker handles ~2 groups before the [refill barrier](../../../src/var_calling/per_group_merger.rs#L459) (`.collect()`); 88.8 % of total samples are inclusive in `rayon_core::WorkerThread::wait_until_cold` (workers parked between batches). Whatever the optimum is, 32 is not measured to be it.
- **Why fifth:** Pure bench addition, no source change to the merger logic — the `PerGroupMergerConfig::batch_size` knob already exists. Ship the sweep, pick the winner, update the constant. Pairs with the Likely "streaming pipeline" finding (L8) as the diagnostic step: if the sweep flattens out at low batch sizes, the barrier is the bottleneck; if it improves linearly with batch size, the barrier is fine and per-batch overhead dominates.
- **Fix:** loop sketch in [concurrency.md](../../../tmp/perf_review_2026-05-16_var_calling/concurrency.md) Hot-path §1 and [methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Hot-path §2 (the two sub-agents converged on the same site).

#### H6: `benches/var_calling_perf.rs:151, :172, :295, :321` — Stage 4 benches rebuild the entire fixture inside `iter_batched` setup; profiles are unusable for Stage 4 code findings

- **Confidence:** High
- **Hot-path evidence:** `profile_merger_dense.json.gz` shows 67 % inclusive time in `c_with_alloca` / `iter_batched` and only 21 % in `PerPositionMerger::next`; `profile_grouper_snp_dense.json.gz` is 77 % vs 13 %. The fixture-rebuild dominates the sampling profile of the function under test, per the methodology rule *"`iter_batched` setup is sampled under `--profile-time`"*.
- **Why sixth:** Unblocks every Stage 4 code-level finding. Until this is fixed, no merger or grouper rewrite can be defended with profile evidence. The right shape is to build the per-sample-records master *once* outside `b.iter`, then per-iter clone an iterator view over it; sketches in [methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Hot-path §1.
- **Fix:** sketch in [methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Hot-path §1. Acceptance: `PerPositionMerger::next` inclusive ≥ 50 % in the new profile.

#### H7: `benches/var_calling_perf.rs:488` — No cohort-size sweep; Stage 5 / Stage 4 perf at N=10, 256, 1024 is unmeasured

- **Confidence:** High
- **Hot-path evidence:** every Stage 5 bench hardcodes `N_SAMPLES = 64`; the grouper hardcodes 16. The performance intent calls for N=10–1000. The deferred-review memo on 2026-05-16 also flagged this. Multiple sites (`UnifiedAllele.per_sample_sources`, `UnifiedAllele.chain_anchor_counts`, `project_scalars` matrix) scale with `n_samples`; whether the headline "70 % CPU in allocator" stays flat or grows super-linearly with N is currently unknown.
- **Why seventh:** Pure bench addition. Names the next code-level priority *after* H1–H4 land (per-sample data layout vs n_alleles-bound paths). Pair with a follow-up profile capture so the cohort-size N=1024 hot sites are known.
- **Fix:** `bench_with_input` sketch in [methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Hot-path §3.

### Likely

#### L1: `Cargo.toml:74-80` — `[build] alloc-mimalloc` feature comment should record the negative result

- **Confidence:** High
- **Evidence:** `mimalloc_full_biallelic.txt` shows +47.7 % slower on the dominant workload. Documentation-only change so a future reviewer (human or Claude) doesn't retry it without context.
- **Why Likely (not Hot-path):** doc fix, not a code-level perf win. Detailed wording in [methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Likely §1.

#### L2: `src/var_calling/per_group_merger.rs:1357-1373` — Replace `std::collections::HashMap` in `build_source_index` with pre-sized `ahash::AHashMap`

- **Confidence:** Medium
- **Evidence:** `hashbrown::RawTable::reserve_rehash` at 0.8 % self-time + `core::hash::BuildHasher::hash_one` at 0.5 % both attribute back to this single std-HashMap site ([all_symbolicated.txt](../../../tmp/perf_review_2026-05-16_var_calling/all_symbolicated.txt)). The map is built per group; keys are integer triples (no adversarial input).
- **Why Likely:** A 1.3 % cumulative self-time gain is real but small relative to H1–H3. Trivial change. Diff in [allocations.md](../../../tmp/perf_review_2026-05-16_var_calling/allocations.md) Likely §1.

#### L3: `src/var_calling/per_group_merger.rs:1471` — Precompute a per-sample `u64` chain-broken-compound mask; the per-cell `find_map` is recomputable

- **Confidence:** Medium
- **Evidence:** `compute_log_likelihoods:1471` at 0.2–0.3 % self-time across all three Stage 5 profiles; structurally a `Vec<UnifiedAllele>` indirection + nested `Vec<Vec<bool>>` lookup per `(sample, genotype, allele-in-genotype)` cell. Hot-loop sub-agent confirms a bit-exact rewrite using `g_bits & broken_compound_masks[sample_idx]`.
- **Why Likely:** small self-time share, but a clean cache-friendly rewrite with no algorithmic change. Lands as a small follow-on to H2 (which already flattens `chain_anchor_flags`). Code in [hot_loops.md](../../../tmp/perf_review_2026-05-16_var_calling/hot_loops.md) Likely §1.

#### L4: `src/var_calling/per_group_merger.rs:1547, :1565` — Rewrite `standard_log_likelihood`'s loops with `bits.trailing_zeros()` bit-iteration; preserves FP order

- **Confidence:** Medium
- **Evidence:** `standard_log_likelihood` reached via line 1490 at 1.9 % self-time. The data-dependent `if (g_bits >> a) & 1 == 0 { continue; }` blocks autovec; iterating set / unset bits via `trailing_zeros + bits &= bits - 1` produces straight-line code.
- **Why Likely:** Tiny gain expected because the loop body is short and the allocator dominates; the change is justified by codegen cleanliness more than measured wall time. Validate via `cargo asm` *and* a bench run; merge only if biallelic moves by ≥ 2 %. Snippet in [hot_loops.md](../../../tmp/perf_review_2026-05-16_var_calling/hot_loops.md) Likely §2.

#### L5: `src/var_calling/per_group_merger.rs:655, :1215, :1302, :1427` — `UnifiedAllele.chain_anchor_counts: Vec<u32>` → `Option<Vec<u32>>` (lazy allocation, compound-only)

- **Confidence:** Medium
- **Evidence:** pattern-match against the per-UnifiedAllele allocation shape called out at the 1.0 % `project_scalars:1138` hot site. Field is *only ever written for compound alleles*; for biallelic SNPs every allele pays an `n_samples`-wide zero vector that is never read. All three readers already gate on `if !allele.is_compound` so `None` ⇒ skip is the natural shape.
- **Why Likely:** Internal-only (`UnifiedAllele` is `pub(crate)`), trivial diff, makes biallelic groups strictly cheaper. Subsumed by H2 if H2's CSR variant goes all the way; if H2 only does row-major flattening, L5 is the cheaper sibling fix. Detail in [data_layout.md](../../../tmp/perf_review_2026-05-16_var_calling/data_layout.md) Likely §5.

#### L6: `benches/var_calling_perf.rs` — No multiallelic (`n_alleles ≥ 5`) bench; the inner-loop scaling against the documented cap of 6 is unmeasured

- **Confidence:** Medium
- **Evidence:** all Stage 5 benches sit at `n_alleles = 2` (biallelic) or 3 (compound). `genotype_order(2, 6) = 21` vs `(2, 2) = 3` — 7× more inner-loop work per sample. Whether `standard_log_likelihood` scales as expected at the documented cap is unverified.
- **Why Likely:** Pure bench addition. Names whether multiallelic workloads need their own attention or stay below the per-sample/cohort scaling concerns. Sketch in [methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Likely §2.

#### L7: `benches/var_calling_perf.rs:155-161, :303-307, :512-517` — Bench bodies have no internal-state assertions

- **Confidence:** Medium
- **Evidence:** every bench shape follows `for item in iter { let r = item.unwrap(); black_box(&r); count += 1; } black_box(count)` without `assert_eq!(count, EXPECTED)`. `Throughput::Elements(N_GROUPS)` claims per-group timing, but a future refactor that silently increases the drop rate (e.g. in `enforce_max_alleles` or REF-only group filtering) would show up as a faster bench.
- **Why Likely:** One-line addition per bench body. Per the methodology rule *"Benches that align with the API's internal state need an in-bench verification assertion"*. Expected counts: `N_GROUPS` for the per_group_merger biallelic+compound shapes; `N_POSITIONS` for the dense merger; `N_GROUPS_EXT` for grouper overlap_extension; `N_POSITIONS_SNP / snp_every_n = 1000` for grouper snp_dense. Detail in [methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Likely §3.

#### L8: `src/var_calling/per_position_merger.rs:306` — `vec![None; self.n_samples()]` per emitted position

- **Confidence:** Low (downgraded from Likely)
- **Evidence:** Stage 4 profile is fixture-dominated (depends on H6 first). For 64 samples × 200 000 positions per iter, this is 12.8 M Vec-allocations per bench iter — if Stage 4 ever shows up in an end-to-end profile, this is the candidate.
- **Why Likely (capped):** the *only* allocation-free version of this is a lending-iterator API change (`Iterator` cannot lend through `&mut self`), which is significant API surgery for a stage that currently isn't the bottleneck. The simpler "scratch Vec on the merger + `mem::take`" rewrite trades the `vec![None; n]` for a `mem::take + Vec::with_capacity(n) + resize(None)` — same allocation count, no win. Defer until H6 unblocks a clean Stage 4 profile.

### Speculative

#### S1: `src/var_calling/per_position_merger.rs:45` — `PerPositionPileups.per_sample: Vec<Option<PileupRecord>>` sparse representation for low-coverage cohorts

- **Confidence:** Low
- **Evidence:** Pattern-match — Stage 4 profile is fixture-dominated; the 71 % `drop_in_place<OverlappingVariantGroup>` cost is dominated by the *inner* allocations (`Vec<u8>`, `Vec<ChainId>`) more than the outer `Vec<Option<PileupRecord>>` slab.
- **Why Speculative:** big API change (sparse `Vec<(u32 sample_idx, PileupRecord)>`, breaks every `.len()` reader). Only a second-order win after the inner SmallVec-or-equivalent work; defer. [data_layout.md](../../../tmp/perf_review_2026-05-16_var_calling/data_layout.md) Likely §6.

#### S2: `src/per_sample_pileup/pileup/mod.rs:379-385` (cross-cutting) — `AlleleObservation.seq` and `chain_ids` SmallVec candidates

- **Confidence:** Low
- **Evidence:** Inclusive 71 % `drop_in_place<OverlappingVariantGroup>` bottoms out at *"each `Vec<AlleleObservation>` → `Vec<u8>` seq + `Vec<ChainId>` chain_ids"*. The 1-bp biallelic SNP case is the dominant workload; every allele's `seq` is one byte, which fits trivially inline.
- **Why Speculative + out of scope:** the change touches Stage 1 / Stage 2 / `.psp` reader / every test fixture — crate-wide blast radius. Logged here because the *cost* is consumed in Stage 5's drop cascade; the *fix* lives in pileup/mod.rs. **Route to a separate cross-cutting review.**

#### S3: `src/var_calling/per_group_merger.rs:459-477` — Replace `.collect::<Vec<_>>()` barrier in `refill` with a streaming pipeline

- **Confidence:** Low
- **Evidence:** 88.8 % inclusive in `wait_until_cold` is consistent with workers idle between batches. Compound paths are ~2× slower per group than biallelic — in mixed workloads the tail-latency-per-batch is the structural cap on throughput.
- **Why Speculative:** depends on H5 (batch-size sweep) outcomes; the order-preservation + error short-circuit semantics in [the current `refill` comment](../../../src/var_calling/per_group_merger.rs#L451-L458) are load-bearing and the streaming rewrite is substantial. Two-step plan in [concurrency.md](../../../tmp/perf_review_2026-05-16_var_calling/concurrency.md) Likely §1: add a per-batch min/max latency probe first; only escalate to a rewrite if the gap is ≥ 20 % of median per-group time.

#### S4: No `RAYON_NUM_THREADS` / CPU-pinning policy documented

- **Confidence:** Low
- **Evidence:** Pattern-match. The 88.8 % `wait_until_cold` figure is consistent with a default-sized pool with no pinning. Plausible interaction with the cross-thread free pattern that drives the allocator-contention diagnosis. Detail in [concurrency.md](../../../tmp/perf_review_2026-05-16_var_calling/concurrency.md) Speculative §4.
- **Why Speculative:** `taskset -c 0-15` + `RAYON_NUM_THREADS={8,16}` experiment first; documentation-only fix if it pays out. Threshold: ≥ 10 % improvement on any of the three Stage 5 shapes.

#### S5: `Cargo.toml` PGO

- **Confidence:** Low
- See §4 above. Defer until Wave 1 lands.

#### S6: `src/var_calling/per_group_merger.rs:459` — Per-`refill` `Vec<Result<…>>` allocation

- **Confidence:** Low
- **Why Speculative:** the result-vector allocation is one heap allocation per refill (32 entries at current batch size), dwarfed by the per-`MergedRecord` inner Vec<Vec<_>> drop cost (which H2 attacks). Re-evaluate only after H5 settles `batch_size` — at `bs = 2048` the result vector is 2048 entries and the picture may shift. Detail in [concurrency.md](../../../tmp/perf_review_2026-05-16_var_calling/concurrency.md) Speculative §2.

#### S7: `src/var_calling/per_group_merger.rs:359-364` — `Arc<dyn RefSeqFetcher + Send + Sync>` vtable indirection

- **Confidence:** Low
- **Evidence:** `ref_fetcher.fetch(...)` is called once per group; not visible in any of the three Stage 5 profiles. In production the underlying CRAM/FASTA I/O dwarfs vtable cost.
- **Why Speculative:** record-keeping only — if a future profile (with a real CRAM-backed fetcher) shows `fetch` self-time > 1 %, reopen.

#### S8: `src/var_calling/per_position_merger.rs:1` — Stage 4 serial throughput as feed-rate ceiling for Stage 5

- **Confidence:** Low
- **Why Speculative:** depends on H6 (fixture-dominated profile fix) *and* a new microbench isolating `self.upstream.next()` cost. If the pull-cost × `batch_size` ≥ ~10 % of parallel-batch wall-time, this becomes a real finding. Logged so the team knows to ask the question post-H6. Detail in [concurrency.md](../../../tmp/perf_review_2026-05-16_var_calling/concurrency.md) Speculative §3.

#### S9: `src/var_calling/per_group_merger.rs:459-477` — Arena/pool the per-`refill` batch + worker scratch

- **Confidence:** Low
- **Why Speculative:** the structural answer to the 70 % allocator cost, but the lifetime threading through `UnifiedAllele`, `UnifiedAlleleSet`, `ScalarProjection`, `CompoundCandidate`, `ChainProposal` is significant and the output `MergedRecord` still needs a heap copy on emit. Detail in [allocations.md](../../../tmp/perf_review_2026-05-16_var_calling/allocations.md) Speculative §1.

### Note

- **Positive sighting** — `Arc::clone(&self.ref_fetcher)` is hoisted *outside* the `into_par_iter().map(...)` closure at [per_group_merger.rs:450](../../../src/var_calling/per_group_merger.rs#L450), so the atomic increment happens once per refill, not per group. Guard the invariant — a future edit that moves `Arc::clone` inside the worker closure would regress this without test failure. ([concurrency.md](../../../tmp/perf_review_2026-05-16_var_calling/concurrency.md) Note §1.)
- **Stage 4 has no concurrency primitives** — `grep` confirms zero `rayon`/`Arc`/`Mutex`/`Atomic`/`channel`/`spawn` references in `per_position_merger.rs` or `variant_grouping.rs`. Stage 4 is serial by design (linear-scan k-way merge + serial grouper). No finding; confirms the orchestrator's expectation.
- **`xlogy` and `ln_factorial` are correct as-is.** `xlogy` is two ops + a branch; LLVM inlines it within the TU (and `lto = "fat"` makes the across-TU question moot). `ln_factorial` is table-backed for `n < 1024`; WGS depth per-allele is O(10)–O(100), so the iterative fallback never trips on real data. ([hot_loops.md](../../../tmp/perf_review_2026-05-16_var_calling/hot_loops.md) Note §1, §2.)
- **`AlleleSupportStats` field layout is already minimal.** `#[repr(Rust)]` reorders for minimum padding; explicit field ordering would be necessary only with `#[repr(C)]`. ([data_layout.md](../../../tmp/perf_review_2026-05-16_var_calling/data_layout.md) Speculative §10.)
- **No false sharing on `SharedRefFetcher`.** The trait method takes `&self`; the only shared mutable state in Stage 5 is the rayon scheduler's internal queues, which crossbeam-deque already pads correctly. ([data_layout.md](../../../tmp/perf_review_2026-05-16_var_calling/data_layout.md) Speculative §11.)
- **`sample_size = 10`** in the criterion harness is at the documented floor. Acceptable given the per-sample wall-time budget; cross-commit comparisons in the 5–10 % range need a revert experiment per the methodology rule. ([methodology.md](../../../tmp/perf_review_2026-05-16_var_calling/methodology.md) Note §1.)

---

## 6. Out-of-scope observations

- **`AlleleObservation.seq: Vec<u8>` and `chain_ids: Vec<ChainId>` are the structural origin of the 71 % `drop_in_place<OverlappingVariantGroup>` cost.** The fix (SmallVec / inline storage) lives in [src/per_sample_pileup/pileup/mod.rs:379-385](../../../src/per_sample_pileup/pileup/mod.rs#L379-L385), which is out of scope for this perf review. Once the Wave 1 changes here (H1–H4) reduce the per-group allocation count, route the residual drop cost to a cross-cutting review that covers Stage 1 / Stage 2 / `.psp` reader / the var_calling consumer in one pass.
- **`OverlappingVariantGroup` Drop cost as a Stage 4 ↔ Stage 5 seam redesign target.** The arena / bumpalo-backed group is the structural answer if the per-group merger needs another 2× after Wave 1+2. Filed S9 / [data_layout.md](../../../tmp/perf_review_2026-05-16_var_calling/data_layout.md) Speculative §7 as the future redesign target.
- **`detect_compound_candidates` keying its outer `BTreeMap` on `Vec<(usize, usize)>`** allocates a `Vec` per chain proposal as the key ([per_group_merger.rs:884](../../../src/var_calling/per_group_merger.rs#L884)). Similar pattern to H1's `BTreeMap<Vec<u8>, _>`; the same `AHashMap` swap applies but the call frequency is lower (one map per group vs O(n_alleles) inserts in H1). Follow-up after H1 if it shows up.
- **`build_chain_proposals` builds a `BTreeMap<ChainId, Vec<CompoundConstituent>>` per sample inside a per-sample loop** ([per_group_merger.rs:928](../../../src/var_calling/per_group_merger.rs#L928)). For biallelic SNP workloads with no compounds it constructs N empty trees. Cross-category; same `AHashMap` line of attack as H1 / L2 / the previous bullet.
- **`enforce_max_alleles` removes prunable alleles with `Vec::remove(i)` in a loop** ([per_group_merger.rs:1093](../../../src/var_calling/per_group_merger.rs#L1093)). Known anti-pattern (`O(n^2)` worst case) but `n_alleles ≤ MAX_BITMASK_ALLELES = 64` and the default cap is 6, so the bench never triggers it. Filed for future awareness if `max_alleles` is ever raised significantly.

---

## 7. What's already good

- **`process_group` is a pure function** taking `(group, &dyn RefSeqFetcher, &PerGroupMergerConfig)` and returning a `Result` — no shared mutable state, no `Mutex`, no atomic counters. The rayon parallelisation at [per_group_merger.rs:459-462](../../../src/var_calling/per_group_merger.rs#L459) inherits that cleanliness for free. The headline allocator-contention finding has nothing to do with how this code uses rayon; it's about the data shape being thrown into it.
- **`Arc::clone(&self.ref_fetcher)` is hoisted out of the `par_iter` closure** ([per_group_merger.rs:450](../../../src/var_calling/per_group_merger.rs#L450)) — one atomic increment per refill, not per group. The natural anti-pattern was avoided.
- **`standard_log_likelihood`'s bitmask scratch is already stack-resident** ([per_group_merger.rs:1546](../../../src/var_calling/per_group_merger.rs#L1546)) — `p_counts: [u32; MAX_BITMASK_ALLELES]` is the convention H3 propagates to `chain_broken_log_likelihood`. `ln_factorial` is table-backed via `LazyLock` ([per_group_merger.rs:1697](../../../src/var_calling/per_group_merger.rs#L1697)). Both are already what the perf review would recommend.

---

### Author response convention

Address each finding by its identifier (H1, L1, S1, …) with one of: `applied in <commit>` / `experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. The "experiment shows no gain" path is welcomed — it's what the measurement plan is for.
