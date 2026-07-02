# Performance Review: paralog-filter
**Date:** 2026-07-02
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the hidden-paralog filter's two scoring passes (calibrate S4 + write S5) and the pure per-locus kernel `score_locus_for_paralogy`
**Verdict:** Apply the listed wins
**Hot-path evidence:** wall-time measurement only (no sampling profile of the scoring path yet)

---

## 1. Scope and constraints

- **What was reviewed:** the paralog filter's post-main-pass work — the two spill-scoring passes and the marginal-likelihood kernel they call once per locus, twice each.
- **Reviewed against:** branch `tomato2-paralog-filter` as-provided (last commit `31a964d`).
- **Targets / input sizes:** 59-sample tomato2 cohort, 376 k scored loci, 58 usable samples/locus. The filter is **on by default**. Target: recover most of the 16× the main caller pass already gets, so the on-by-default path does not carry an +11.6× wall tax. Target hardware: dev = 6 fast cores (macOS Apple container), prod = 32 cores (Linux). **Hard constraint:** the surviving VCF must stay **byte-identical**, and the calibrate/write passes must produce **bit-identical** likelihood ratios (the FDR cut is applied to the histogram it was built from).
- **Hot-path evidence available:** wall measurement from [paralog_t1_t2_2026-07-02.md](../implementations/paralog_t1_t2_2026-07-02.md): filter on 272.2 s vs off 23.4 s (**+248.8 s, +11.6×**); peak RSS flat (+8 %); spill 3.5 GB read twice (~7 GB, ~15 s). ~660 µs/locus/pass. **No sampling profile and no criterion bench exist for this path** — so per the rubric, code-*internal* candidates are capped at Likely; the pass-level regression is directly measured.
- **In-scope files:**
  - [src/paralog/locus_score.rs](../../../../src/paralog/locus_score.rs) — the pure kernel
  - [src/paralog/mod.rs](../../../../src/paralog/mod.rs) — grid/param defaults
  - [src/var_calling/paralog_filter/calibrate.rs](../../../../src/var_calling/paralog_filter/calibrate.rs)
  - [src/var_calling/paralog_filter/write_pass.rs](../../../../src/var_calling/paralog_filter/write_pass.rs)
  - [src/var_calling/paralog_filter/spill.rs](../../../../src/var_calling/paralog_filter/spill.rs)
  - [src/var_calling/pipeline.rs](../../../../src/var_calling/pipeline.rs) (call site ~L755–830)
- **Deliberately out of scope:** the Stage-1 pileup path, the main caller EM, DUST — only the paralog filter's added work.
- **Categories dispatched:** methodology (always); concurrency (the passes are serial while the main pass is 16-way); hot_loops (the kernel is a transcendental-heavy numeric loop); allocations (per-locus `Vec`s); io_and_syscalls (spill read twice, per-survivor FASTA fetch).

## 2. Verdict

**Apply the listed wins.** The regression is not a memory or I/O problem — peak RSS is flat and the double spill read is ~15 s of the +249 s. The entire regression is the pure per-locus scorer running **single-threaded, twice**, over 376 k loci while 15 cores sit idle (the main pass already finished). Two levers are well-evidenced:

1. **Parallelize both scoring passes** (H1) — the scorer is pure and per-locus-independent; the histogram fold is order-independent; the write pass can score in parallel and re-serialize before the writer. This is a structural certainty independent of which kernel line is hottest, and directly divides the +249 s by the core count.
2. **Memoize the locus-invariant log-prior tables** (L1) — H1 recomputes the same 200×58 Wright genotype log-priors and H2 the same 40×58 carrier log-probs *at every locus* (~26 billion redundant `ln` in H1 alone), a pure bit-identical hoist.

Add the missing criterion bench (§3) as the gate for all of them, and — because parallelism introduces a new RAM cost against a deliberately RAM-flat subsystem — sweep RSS-vs-batch-size. A sampling profile is still worth collecting to *order* the kernel-internal micro-opts (L1/L2/L3), but it is not a prerequisite for the parallelism win.

## 3. Measurement plan

In the order that unblocks the rest:

1. **Add `benches/paralog_scoring_perf.rs`** (`harness = false`, register in `Cargo.toml`), two groups:
   - **Kernel microbench** — `score_locus_for_paralogy` over a fixed 58-sample synthetic locus (also N ∈ {20, 200} to bracket), fixtures built *once outside* `b.iter`, `black_box` in and out. Isolates the H1 (200-pt) + H2 (40×7) inner-loop cost.
   - **Streaming bench** — write M synthetic `ParalogSpillRecord` + sibling `WindowSpillRecord` to a `Cursor<Vec<u8>>` via the existing `test_support::{write_spill, write_window_spill}`, then time `calibrate(...)`. Exercises decode → join → score → fold so a kernel win that decode/join masks is caught. Add an in-bench assertion that M loci folded (fixture-wiring guard).
   - Metric: criterion median with `--save-baseline before` / `--baseline before`. Act threshold: ≥5 % on the streaming group, two clean measurements per side.
2. **Collect one `cpu-clock` sampling profile** on the Linux dev box (perf targets the VM kernel on macOS Apple container — see CLAUDE.md): `./scripts/dev.sh cargo flamegraph -e cpu-clock --bench paralog_scoring_perf -- --bench calibrate_streaming` (or `perf record -e cpu-clock -g`). Splits the serial time between decode/join, H1's 200-pt loop, and H2's `log_add_exp` loop — this *orders* L1/L2/L3, it does not gate H1.
3. **Thread sweep for the parallelism PRs** — 59-sample tomato2, filter on, `--threads {1,2,4,6}` (dev) and {…,32} (prod), best-of-3 median. Timer around `calibrate(...)` (pipeline.rs:787) and `run_write_pass(...)` (pipeline.rs:815) separately. **Correctness oracle (non-negotiable): `π`/`lr_threshold` bit-identical and the output VCF byte-identical** (`cmp`) vs the serial run.
4. **Batch-size / RSS sweep** — the parallel batch holds N owned `PosteriorRecord`s + joined windows at once (new RAM cost vs the current one-record-at-a-time design). Sweep N ∈ {1k, 8k, 64k}; pick the smallest N that saturates the pool; record peak RSS at that N (must stay within the memory budget — the spill exists precisely to keep RAM flat).
5. **DHAT for the allocation findings** — run the streaming bench under the existing `dhat-heap` feature; read allocation counts attributed to `enumerate_carrier_configs` / `h1_log_likelihood` / `h2_log_likelihood` before/after. Allocation hypotheses are confirmed by count, not by wall time.

## 4. Build / toolchain configuration

**No change.** The release profile is already fully tuned: [Cargo.toml:45-49](../../../../Cargo.toml#L45-L49) sets `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`; `.cargo/config.toml` sets per-target `target-cpu` (`x86-64-v3` Linux, `apple-m1` macOS); the toolchain is pinned (1.95); and `[profile.bench]` / `[profile.profiling]` carry debug info for profiling. The regression is a parallelism/algorithm problem, not dropped compiler speed. (PGO is theoretically available for this repeated lab workload but is a separate, later experiment that should not precede the parallelism win.)

## 5. Code-level findings

### Hot-path

**H1: [src/var_calling/pipeline.rs:784-827](../../../../src/var_calling/pipeline.rs#L784-L827) (+ calibrate.rs:278, write_pass.rs:153) — Both scoring passes run single-threaded while the whole 16-way pool sits idle**

`calibrate` ([calibrate.rs:278](../../../../src/var_calling/paralog_filter/calibrate.rs#L278)) and `run_write_pass` ([write_pass.rs:153](../../../../src/var_calling/paralog_filter/write_pass.rs#L153)) are serial `while let Some(record) = …next_record()` loops, invoked **sequentially, after** the main caller pass has finished — so the added ~249 s runs on one core while 5 (dev) / 31 (prod) are idle. The per-locus scorer `score_locus_for_paralogy` is documented pure and per-locus-independent. This is the dominant term in the +11.6× regression and the highest-value lever.

The fix is a **serial-read+join → rayon-score-batch → serial-reduce/emit** pipeline (the project's own "split data-shaping from math" pattern), because two pieces must stay single-threaded:
- **`WindowJoin` is a stateful forward cursor** (monotone `(chrom_id, tile)`, `debug_assert!`ed) — it cannot be shared across workers. Decode + join on one thread, producing self-contained `(record, window_gc, window_depths)` tuples; fan only the *scoring* out.
- **Calibrate reduce is free of ordering concerns:** the histogram fold is integer bin-count increments — commutative and associative — so per-worker histograms summed (or an in-order serial fold of the batch's `Option<f64>`) give **bit-identical** bins regardless of scoring order → identical `π`/cut.
- **Write pass must re-serialize before emit:** `check_reference_consistency` uses a monotonic-forward `StreamingChromRefFetcher` and `writer.handle` needs a gapless `chunk_order`. So `par_iter` computes each locus's `flagged: bool`, then a serial phase walks the batch **in index order** feeding survivors to the fetcher + writer — byte-identity holds because the survivor set and their order are unchanged; only the pure flag decision moved off-thread.

`obs_buf` must become per-worker scratch (`map_init`), not the single shared `&mut Vec`.

- **Measurement plan:** §3 steps 1, 3, 4. Ship when the scoring portion scales toward 1/N at T=6/32 with `π`/threshold/VCF identical. If a pass turns out spill-read-bound rather than scorer-bound (the ~15 s I/O figure says it won't), refute and stop.
- **Complexity cost:** a batching layer + per-worker scratch + (calibrate) per-worker histogram merge / (write) in-order re-serialization. New peak-RSS cost = N records held at once — bound via N (see H1↔memory tension below). rayon already in-tree; crate is `#![forbid(unsafe_code)]` so no unsafe. Largest change in the review; its own PR, its own before/after.

### Likely

**L1: [src/paralog/locus_score.rs:232-234](../../../../src/paralog/locus_score.rs#L232-L234) and [:296-301](../../../../src/paralog/locus_score.rs#L296-L301) — Precompute the two locus-invariant log-prior tables once, not 376 k×2 times**

Verified against the caller: `inbreeding` (F_s) and `single_copy_depth_sd` (σ₀) are cohort-length slices built once in the prepass and passed **unchanged** to every locus; `params` (both grids) is constant. Therefore:
- **H1:** `wright_genotype_log_priors(p_i, F_s)` is a pure function of `(p_i, F_s)` with a fixed p-grid and fixed F_s — the full **200×58** table of `(log_homref, log_het, log_homalt)` is identical at every locus, yet rebuilt 376 k×2 times. The inner loop calls it ~11,600×/locus, 3 `.ln()` each → ~34,800 `ln`/locus, **~26 billion `ln` across the run** — the bulk of H1's transcendental work.
- **H2:** `log_noncarrier[s]` / `log_carrier[s]` depend only on `(q_j, F_s)` — a **40×58** table, ~4,640 `ln`/locus removable. (The dominant H2 cost, ~16,240 `log_add_exp`/locus at [:307](../../../../src/paralog/locus_score.rs#L307), depends on per-locus AD/coverage and stays.)

Ship both as one signature revision: build a `PrecomputedTables` struct once in the driver and thread it in. Index by the **original cohort index** (add a `cohort_index` field to `UsableSample`), because `usable` is a variable-length per-locus subset.

- **Measurement plan:** kernel microbench before/after; `cargo asm` on `h1_log_likelihood`/`h2_log_likelihood` to confirm the inner `.ln()` sites vanish. **Bit-identity gate (exact float equality, not tolerance)** over a grid of synthetic loci — the suggested build reuses `sfs_grid_point`/`linspace_point`/`PROB_FLOOR` verbatim so every `.ln()` argument is bit-identical; calibrate↔write identity is automatic (same function, same table).
- **Complexity cost:** one precompute struct + a `cohort_index` field + one argument threaded through `score_joined_locus`/`score_spilled_locus`. No unsafe, no dependency. **Correctness-adjacent** (touches the scored value's provenance) — hence the exact-equality gate.

**L2: [src/paralog/locus_score.rs:133](../../../../src/paralog/locus_score.rs#L133) / [:319-345](../../../../src/paralog/locus_score.rs#L319-L345) — `enumerate_carrier_configs` rebuilds a locus-invariant `Vec` (+ 7×`sqrt`/`ln`) on every call**

The configs depend only on `params.carrier_copy_numbers`, constant for the run, but a fresh `Vec<CarrierConfig>` (7 elems, each `sqrt` + 2 `ln`) is allocated per locus inside the serial loop. Hoist to once-per-pass and pass `&[CarrierConfig]` in (requires `CarrierConfig` → `pub(crate)`; the `params.is_empty()` guard moves to the caller). Removes 376 k×2 allocations + the transcendental recompute; folds naturally into the L1 signature change and the H1 parallel workers (all borrow one shared slice).
- **Measurement plan:** kernel microbench (expect the largest relative effect at small N); DHAT allocation-count. Keep if ≥5 % at N=20 or a clear DHAT share. **Complexity cost:** minor signature churn, bit-identical.

**L3: [src/paralog/locus_score.rs:146,204,263-264,286-287](../../../../src/paralog/locus_score.rs#L146) — Six per-locus scratch `Vec`s can be reused like the existing `obs_buf`**

`usable`, `allele_by_genotype`, `noncarrier_base`, `carrier_branch` (also zero-fills 406 f64 overwritten immediately), `log_noncarrier`, `log_carrier` — all freshly allocated per locus (~6 allocs+frees × 752 k). The project already applies the load/use/clear/reload pattern to `obs_buf` one layer up; push it into the kernel via a `ScorerScratch` struct (per-worker once parallelized — sequence **after** H1 so the plumbing isn't reworked twice). RSS stays flat (scratch is `O(samples×configs)`, same order as the transient peak).
- **Measurement plan:** this is an *allocation* hypothesis — confirm with DHAT (count), not wall time. Merge on a meaningful allocation drop with non-negative wall delta. **Complexity cost:** one struct + `&mut` arg; `h1`'s `.collect()` becomes a manual clear+push. Bit-identical (same values, reused storage). May show no wall gain if math dominates — then file as "correct but not worth it."

### Speculative

**S1: [src/var_calling/pipeline.rs:784-828](../../../../src/var_calling/pipeline.rs#L784-L828) — Score once in calibrate, spill the per-locus LR (one f64) to disk, read it back in the write pass**

Removes the entire second scoring sweep (~half the ~249 s if symmetric). The documented recompute-not-cache invariant ([write_pass.rs:15-18](../../../../src/var_calling/paralog_filter/write_pass.rs#L15-L18)) forbids holding LRs **in RAM**; a **disk** LR spill (376 k × 8 B ≈ 3 MB, negligible vs the 7 GB spill) keeps RAM flat and sidesteps it — arguably *stronger* determinism (the write pass applies the cut to the exact stored value). **But:** (a) it overturns a stated design decision → owner sign-off; (b) it is *complementary* to, not a replacement for, H1 (parallelism divides by cores; caching halves the work — do both for ~2N×); (c) the write pass still reads the full record spill to feed the writer, so this saves the *recompute*, not the second read; (d) more code than the batched `par_iter` (new spill stream + join + round-trip proptest, including the "unscored → kept" case). **Do H1 first; consider S1 only after, with owner sign-off.**

**S2: [src/var_calling/paralog_filter/write_pass.rs:88](../../../../src/var_calling/paralog_filter/write_pass.rs#L88) — Per-survivor 1-byte reference fetch allocates a fresh `Vec`**

`StreamingChromRefFetcher::fetch` allocates a `Vec` for a single byte, ~356 k times. **Not** a per-call syscall (1 MiB sliding buffer refills ~once/MiB on a monotone walk) — the cost is allocator churn, not I/O. Trivial fix: `fetch_into` already exists; hoist a `Vec::with_capacity(1)` into `run_write_pass` and thread it through `check_reference_consistency`. File only if a profile shows it a non-trivial share (the report's framing says it won't be).

**S3: [src/paralog/locus_score.rs:302-311](../../../../src/paralog/locus_score.rs#L302-L311) — H2 `log_add_exp` accumulator autovectorization**

The densest H2 site (~16 k calls/locus), but `log_add_exp` is branchy and calls `exp`/`ln_1p`, which LLVM does not lane-vectorize without a vector-math lib. **Do not** reach for `std::simd`/`wide`; the transcendentals dominate and won't vectorize. No action unless a profile + `cargo asm` show a compiler limitation that isn't the `exp`/`ln_1p` calls — and L1/L2 remove far more work with less risk.

### Note

- **[src/paralog/mod.rs:65,70,83](../../../../src/paralog/mod.rs#L70)** — the grid sizes (SFS 200, carrier-freq 40, 7 configs) are the per-locus flop-count knob, but they are pinned to the validated tomato2 prototype; cutting them changes the science, not just the speed. Out of a perf reviewer's remit — flagged only so a science owner knows the direct cost knob.
- **[src/var_calling/paralog_filter/spill.rs](../../../../src/var_calling/paralog_filter/spill.rs)** — 64 KiB `BufReader`/`BufWriter` on all four spill streams is a sound default (~7 records/syscall at 9 KB/record); explicit `flush()?` in `finish()` (not a swallowed drop-flush); no `fsync` (correct — the file is ephemeral). No change without a syscall trace showing read/write time dominates.
- **[src/var_calling/paralog_filter/spill.rs:534-599](../../../../src/var_calling/paralog_filter/spill.rs#L534-L599)** — the ~9 KB/locus record carries the full genotype-posterior matrix; the **write pass needs it verbatim** to feed `VcfWriter` (byte-identity), so the shared record spill cannot be trimmed. Only S1 removes the calibrate side's dependence.

## 6. Out-of-scope observations

- **H1 ↔ memory tension:** the batched `par_iter` holds N owned records + joined windows at once, against a subsystem whose whole premise is RAM-flatness ([spill.rs:30](../../../../src/var_calling/paralog_filter/spill.rs#L30)). N is the RAM/scaling knob — gate the H1 merge on the RSS-at-chosen-N staying within the memory budget (§3 step 4). Consider wiring N / worker count behind the existing `PVC_CALLER_THREADS`-style env knobs for a consistent story with the main pass.
- **[spill.rs:628-657,773-780](../../../../src/var_calling/paralog_filter/spill.rs#L628-L657)** — `decode_record` walks the large `posteriors`/`gq_phred`/`allele_frequencies` f64 vectors with a per-field `f64::from_le_bytes`; a `chunks_exact(8).map(f64::from_le_bytes)` bulk decode would autovectorize. Cold relative to scoring — pursue only if a profile shows decode dominates a read pass.
- **The footgun from T1** (already noted in PROJECT_STATUS): running the on-by-default filter against summary-less `.psp` silently produces an *unfiltered* callset — a loud warn/error when few samples fit is a correctness follow-up, not a perf item.

## 7. What's already good

- **Memory-flatness by construction** — the spill holds nothing genome-wide; the reader decodes one record into a reused byte buffer ([spill.rs:30](../../../../src/var_calling/paralog_filter/spill.rs#L30)); the calibrate fold is a few-KB histogram. Peak RSS stayed flat (+8 %) across 376 k loci, exactly as designed.
- **The scorer is genuinely pure and deterministic** ([locus_score.rs:127](../../../../src/paralog/locus_score.rs#L127)) — which is what makes the H1 parallelization safe and the L1 memoization a bit-identical hoist. The `LogSumExp` online accumulator avoids materializing the p/q vectors.
- **The caller already reuses `obs_buf`** ([calibrate.rs:277](../../../../src/var_calling/paralog_filter/calibrate.rs#L277)) — the load/use/clear/reload pattern L3 just extends one layer down.

---

## Author response

Worked in the agreed order: bench first, then the non-parallel levers, then (pending) parallelisation.

- **Bench (§3.1) — added.** `benches/paralog_scoring_perf.rs` (kernel microbench over `score_locus_for_paralogy`, `N ∈ {20, 58, 200}`, fixtures + precompute built outside `b.iter`). Baseline `before` (pre-optimisation): N=20 **80.7 µs**, N=58 **245 µs**, N=200 **861 µs**. The N=58 cost × 376 k × 2 ≈ 184 s confirms the kernel is the bulk of the +249 s. Streaming bench over `calibrate` still deferred to the parallelisation phase (needs `pub(crate)` spill surface).
- **L1 + L2 — applied together** (one coherent "hoist all per-pass-invariant data" change). New `ParalogScorePrecompute` (public struct, private fields, `new(params, inbreeding)`) built once per pass in `calibrate` + `run_write_pass` and threaded through `score_joined_locus`/`score_spilled_locus` into `score_locus_for_paralogy`; holds the Wright genotype log-prior table (H1, `sfs_points × cohort`), the carrier/non-carrier log-prob table (H2, `carrier_freq_points × cohort`, also removed the `log_noncarrier`/`log_carrier` per-locus scratch), the carrier configs (L2), and the ε / SFS-weight constants. `UsableSample` gained `cohort_index` to key the tables by cohort position. **Result: −31 % on the kernel at every N** (N=58 245 → 167 µs; p < 0.05). **Bit-identical:** pure hoist (same `.ln()` args, same accumulation order); 1512 lib + 12 cohort-integration tests pass, incl. the parity anchor, the two-pass E2E test, and the byte-identity/determinism guards. Pipeline API unchanged (precompute is internal to the two passes).
- **L3 (per-locus scratch `Vec`s) — measured with DHAT; closed with a one-line presize, full scratch not applied.** Added `examples/dhat_paralog.rs` (scores 50 k synthetic loci, precompute built once). DHAT showed **8 allocations/locus** (~9.8 KB), of which **~4 were `usable: Vec::new()` regrowth** reallocations, plus one each for `allele_by_genotype` (collect), `noncarrier_base` (`with_capacity`), `carrier_branch` (`vec![0.0; …]`). Peak live heap was **11 blocks / 326 KB** (RAM flat, as designed). Applied the free win — `Vec::with_capacity(cohort_size)` on `usable` → **8 → 4 allocs/locus**, peak unchanged; bench N=58 167 → 165 µs. The remaining 4 small warm-allocator blocks are ~0.2 % of the 165 µs/locus (transcendental-bound: the residual cost is the per-locus `log_add_exp` `exp`/`ln_1p` in the H1/H2 LSE loops, which is neither memoisable nor allocation-related), so the full scratch-reuse struct threaded through the public scorer signature is **not worth it**. **Non-parallel levers exhausted: L1+L2+presize = −33 %** (245 → 165 µs at N=58).
- **H1 (parallelise) — pending**, the next phase per the agreed plan.

---

*Per-category audit trail: `tmp/perf_review_2026-07-02_paralog-filter/{methodology,concurrency,hot_loops,allocations,io_and_syscalls}.md`.*
