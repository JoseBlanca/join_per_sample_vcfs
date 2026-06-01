# Performance Review: var_calling
**Date:** 2026-06-01
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** Full `src/var_calling/` cohort pipeline (Stages 3–6 + chunk driver)
**Verdict:** Apply the listed wins
**Hot-path evidence:** sampling profile (macOS `sample` + `samply`) + criterion benches (added this run) + DHAT allocation profile + `cargo-show-asm`

---

## 1. Scope and constraints

- **Reviewed:** the whole `src/var_calling/` subtree — DUST filter, per-position merger, variant grouping, per-group merger, posterior engine (+ `backends`/`interp`/`shape`), contamination estimator, and the chunk driver (driver/loader/columns/column_span_reader/partition/worker + the three columnar kernels).
- **Reviewed against:** commit `9f73feb` on `main` (post the `from_psp` merge).
- **Targets / intent:** the project trades **RAM for sample-count scaling** — scale one cohort run to thousands of samples on a single ~64 GB host. Binding metrics: **wall-time scaling in N** and **peak RSS / allocator churn** as N grows. Hot path = the consume side at high N. Bench host: Apple Silicon arm64; production also runs Linux x86-64.
- **Hot-path evidence available:** YES, three converging streams (this is an *Apply-the-wins* review, not *Profile-first*):
  - **Two new committed benches** in [benches/var_calling_perf.rs](../../../../benches/var_calling_perf.rs): `var_calling_chunk_loader` (produce) and `var_calling_run_window` (consume = 3 kernels + per-group EM), swept over N ∈ {8, 64, 200}. The M10 gap from the 2026-06-01 code review is now closed.
  - **CPU sampling profile** (macOS `sample` over `run_window` N=200) — symbolicated self-time tally.
  - **DHAT allocation profile** (real tomato cohort N=18) + `cargo asm` on the hot math.
- **In-scope files:** all of `src/var_calling/` at `9f73feb`.
- **Deliberately out of scope:** `src/psp/` (PSP decode — separate perf reviews), `src/vcf/` (`build_samples`/`write_record`), `src/fasta/` (ref fetcher — reviewed 2026-05-23), `src/bam/`, `src/pop_var_caller/`. Read for context only.
- **Categories dispatched:** methodology (always), allocations (DHAT shows 51 M blocks), data_layout (columnar CSR + parallel state), concurrency (producer/worker threads + rayon), hot_loops (the EM/QUAL math), io_and_syscalls (PSP/REF reads on the produce side).

## 2. Verdict

**Apply the listed wins.** Three Hot-path candidates are well-evidenced (matching profile/DHAT site, plausible mechanism, contained fix): **H1** the QUAL-convolution scalar `log_sum_exp_2` (the single biggest CPU self-time, and the mechanism behind the N^1.7 consume scaling), **H2/H3** the two `unify_alleles` allocation sites that together account for ~23 M of the 51 M allocated blocks. Apply them one at a time with the new benches + a DHAT re-run as gates. Everything below Hot-path is a measurement-gated experiment.

A reframing worth stating up front: the profile's 51%-self-time `log_sum_exp_2` is **not** in the EM E-step (which takes a homogeneous-fixation fast path under the default `fixation_index = 0.0` config). It is in the once-per-record **QUAL** allele-count convolution `compute_qual_via_exact_af`. The fix target is therefore that convolution, not the EM.

## 3. Measurement plan

In priority order (each gates the finding of the same number):

1. **H1 — QUAL convolution.** Baseline already captured: `var_calling_run_window/snp_dense_{64,200}` = 127 ms / 967 ms. After replacing the inner pairwise `log_sum_exp_2` fold with one `log_sum_exp_slice` over the ≤`ploidy+1` terms, re-run; expect a measurable drop at N=200 (the site is ~51% of self-time). Gate: >5% wall improvement at N=200 with byte-identical VCF on the cohort fixture. Confirm vectorization with `cargo asm` on `compute_qual_via_exact_af`.
2. **H2/H3 — unify_alleles allocation churn.** Re-run DHAT (`cargo run --release --example dhat_var_calling --features dhat-heap -- …`); gate: total blocks drop from 51 M materially (H2 alone is ~15.7 M blocks, H3 ~7.3 M). Cross-check `var_calling_run_window` wall (allocator share rises with N).
3. **L1 — allocator A/B.** Build the *production* binary with vs without a mimalloc `#[global_allocator]`; run `profile_cohort_e2e` thread-scaling at N=200; gate: wall delta > run-to-run noise. (One change, isolated.)
4. **L5/L9/L11 — thread topology.** `profile_cohort_e2e` thread-scaling sweep T=1/2/4/8/16 at N=200, capturing **both** wall and peak RSS; test (a) a separate bounded rayon pool for the loader, (b) `par_iter().with_min_len(…)` on the loader, (c) a `BlockQueue` cap sweep. Lock-wait via `samply`.
5. **Bench hardening (L10):** add a sparse-variant tier (~1/200 density) to exercise the cohort variant-filter drop path, and an in-bench `assert!(partition.n_groups() > 0 && !output.is_empty())` so the N^1.7 curve can't silently measure an empty partition.

## 4. Build / toolchain configuration

`Cargo.toml` is already well-tuned for release/bench: `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`, `[profile.bench]` inherits release with full `debug = true` (so the quoted profiles are release+symbols — valid). `.cargo/config.toml` pins `target-cpu` floors (`x86-64-v3` / `apple-m1`), `rust-toolchain.toml` pins `1.95`. No build-config finding except:

- **L1 (allocator):** the **production binary uses the system allocator** — no `#[global_allocator]` in `src/` or `main.rs`; `alloc-mimalloc` is declared only in bench/example shims. With 15.5 GB / 51 M blocks allocated and allocator wall-share rising 11→38% with N (per the 2026-05-27 scaling measurement), mimalloc is a plausible one-line win. Treat as a build-config experiment (section 3 #3), separate from the code-level allocation fixes (H2/H3) — do not bundle.

## 5. Code-level findings

### Hot-path

**H1: [posterior_engine.rs:3037-3047](../../../../src/var_calling/posterior_engine.rs#L3037-L3047) — QUAL allele-count convolution does a scalar pairwise `log_sum_exp_2` fold in an O(n_samples²·ploidy²) loop**
**Category:** hot_loops. **Confidence:** High (profile + scaling + asm + mechanism all converge).
`compute_qual_via_exact_af` computes the site QUAL via a DP convolution over `P(AC = k)`: `for s in 0..n_samples { for k in 0..=(s+1)*ploidy { for c in 0..=ploidy { acc = log_sum_exp_2(math, acc, prev + ll) } } }`. The innermost [posterior_engine.rs:3043](../../../../src/var_calling/posterior_engine.rs#L3043) `log_sum_exp_2` is **~51% of `run_window` self-time** at N=200 (`sample` tally: 9462 of 18452; samply put 94.6% of leaf samples in one address inside it). It is `#[inline] fn … = m + math.ln(math.exp(a-m) + math.exp(b-m))` — `cargo asm` confirms it emits **out-of-line and fully scalar** (scalar `fmaxnm`, `bl exp`), so every accumulation pays 2×exp + 1×ln **plus** a call prologue. The O(n_samples²) loop is the mechanism behind the measured N^1.7 scaling of `run_window`. Per the 2026-05-18 posterior perf review this is the previously-predicted `log_sum_exp` hot spot (its H3).
**Fix:** rewrite the inner `c`-fold to gather the ≤`ploidy+1` terms `prev + ll` into a small stack buffer and call `log_sum_exp_slice` once — strictly fewer exp/ln (one max + one ln per `k` instead of per `c`), drops the per-call `bl` overhead, and is the prerequisite for vectorizing the `k`-loop (L4). **Complexity:** low — local rewrite of one function; no API/type change; byte-identity must be checked (log-sum-exp reassociation is not bit-identical, so confirm the VCF QUAL field is unchanged within the writer's float formatting, or gate on the existing fixture diff). **Measurement:** §3 #1.

**H2: [kernels/unify_alleles.rs:866](../../../../src/var_calling/kernels/unify_alleles.rs#L866) / [:871](../../../../src/var_calling/kernels/unify_alleles.rs#L871) — `build_chain_proposals_columnar` re-mallocs a per-chain `Vec` every sample×group despite the scratch hoist**
**Category:** allocations. **Confidence:** High (DHAT #1 site, bytes and blocks).
The DHAT #1 site: a `BTreeMap<u64, Vec<CompoundConstituent>>` `or_default()` + `push`, **1.97 GB / 5.35 M blocks (:866) + 664 MB / 10.4 M blocks (:871) ≈ 2.6 GB, 15.7 M blocks**. The 2026-05-29 review's M3 hoisted the *map* onto `UnifyAllelesScratch`, but `BTreeMap::clear()` **drops the value `Vec`s**, so each group re-allocates a fresh per-chain `Vec`. **Fix:** make the value a `SmallVec<[CompoundConstituent; 2]>` (constituents are usually 2 — the `< 2` skip), or replace the per-chain `Vec` with indices into a single reused pool/arena cleared (not freed) per group. Also consider `AHashMap` over `BTreeMap` for the chain key (no ordering need; matches the sibling `byte_index`). **Complexity:** low-medium (SmallVec is a drop-in; arena is more invasive). **Measurement:** §3 #2 (DHAT block count).

**H3: [kernels/unify_alleles.rs:590](../../../../src/var_calling/kernels/unify_alleles.rs#L590) / [:632](../../../../src/var_calling/kernels/unify_alleles.rs#L632) — `project_per_position_into_scratch` clones an owned `Vec<u8>` allele key per distinct allele per group**
**Category:** allocations. **Confidence:** High (DHAT #3/#4 by blocks, ~7.35 M blocks).
`byte_index: AHashMap<Vec<u8>, usize>` allocates a fresh owned-byte key per distinct allele, all dropped at each group's `clear()`. **Fix:** pool the retired key buffers, or move to a byte-arena keyed by `Range`. For the dominant biallelic-SNP case the key is a single base — a small-inline key type would erase it entirely. **Complexity:** medium (touches the dedup structure). **Measurement:** §3 #2.

### Likely

**L1: production system allocator vs mimalloc** — see §4. *methodology.* Build-config A/B on the production binary.

**L2: [worker.rs:453](../../../../src/var_calling/worker.rs#L453) — genotype table rebuilt per group** (`compute_log_likelihoods_columnar` is passed `genotype_tables: &[]`, so it rebuilds the `Vec<Vec<u8>>` via `genotype_order`, surfacing as `collect_non_decreasing` (per_group_merger.rs:554, 1.28 M blocks). The row-shape `PerGroupMerger` already caches this behind `Arc`. **Fix:** cache the per-`(ploidy, n_alleles)` table on `ColumnarPipelineScratch`. *allocations.* **Complexity:** low.

**L3: [posterior_engine.rs:2386](../../../../src/var_calling/posterior_engine.rs#L2386) — `validate_record_shape` builds a `Vec<Vec<u8>>` per record just to read `.len()`.** Replace with the closed-form multiset coefficient `C(n_alleles + ploidy − 1, ploidy)`. *allocations.* **Complexity:** low.

**L4: [posterior_engine.rs:3037-3047](../../../../src/var_calling/posterior_engine.rs#L3037-L3047) — vectorize the QUAL convolution `k`-loop** with a `log_sum_exp_slice_x4`-style 4-lane reduction (the SIMD machinery already exists for the E-step but isn't applied to the dominant QUAL loop). *hot_loops.* Sequenced **after H1** (one change per measurement). **Complexity:** medium (lane masking on ragged `k`-ranges).

**L5: [driver.rs:399](../../../../src/var_calling/driver.rs#L399) + loader `par_iter_mut` — thread oversubscription.** `n_workers` is the full core count spawned as dedicated OS threads, while the producer's `fill_block` fans out over the *global rayon pool* (also full cores); mid-decode, ~2× cores of runnable threads contend for the EM/kernel working sets. **Fix:** give the loader a separate, smaller rayon pool. *concurrency.* **Complexity:** medium; **must preserve** the producer-decode non-starvation invariant ([driver.rs:462-466](../../../../src/var_calling/driver.rs#L462-L466)) — route any change through a correctness check. **Measurement:** §3 #4.

**L6: [driver.rs:1014](../../../../src/var_calling/driver.rs#L1014) — DUST closure copies the fetcher's borrowed `&[u8]` into a fresh `Vec` per sub-span** (`fetcher.fetch(...).map(<[u8]>::to_vec)`), the bulk of the 546 MB DHAT DUST line (~1 MB × ~540 sub-spans). **Fix:** pass a borrow, or hoist one reusable `Vec` + `fetch_into`. *io_and_syscalls.* **Complexity:** low (closure signature). (Note: the DUST mask already streams in 1 MB sub-spans — not a whole-chromosome buffer.)

**L7: [driver.rs:1279](../../../../src/var_calling/driver.rs#L1279) — `ManualEvictChromRefFetcher::for_contig` re-opens the FASTA and re-parses the entire `.fai` per covered interval** (open + full-`.fai`-parse + O(contigs) scan scales with interval count, not chromosome count — costly on fragmented genomes). **Fix:** parse `.fai` once, share the immutable `ContigFai`, one fd per worker. *io_and_syscalls.* **Complexity:** medium (touches the fetcher ctor surface in `src/fasta/`, partly out of scope — coordinate). **Measurement:** `dtruss` open/read count on a real run.

**L8: [kernels/project_scalars.rs:49](../../../../src/var_calling/kernels/project_scalars.rs#L49) — `ProjectedScalarsColumns.scalars: Vec<AlleleSupportStats>` is AoS where the LH walk reads only 2 of 7 fields.** `standard_log_likelihood_columnar` ([compute_log_likelihoods.rs:260-278](../../../../src/var_calling/kernels/compute_log_likelihoods.rs#L260-L278)) touches `q_sum` + `num_obs` (~12 of ~40 bytes/cell); the five bias scalars are read only by out-of-scope `src/vcf/` code. SoA-splitting the two LH-hot fields would ~triple alleles-per-cache-line and unlock autovec. *data_layout.* **Confidence:** Medium — no profile line names that loop directly (its leaf time is downstream in `log_sum_exp_2`). **Measurement:** L1/L2 dcache-miss attribution (Instruments/`samply`) **before** committing; gate on a miss-rate drop.

**L9: [loader.rs:467-487](../../../../src/var_calling/loader.rs#L467-L487) — loader per-position `par_iter` has no `with_min_len`**, so cheap per-element binary searches incur task overhead / tail imbalance. **Fix:** `with_min_len(…)`. *concurrency.* **Complexity:** trivial; cheap experiment. **Measurement:** `var_calling_chunk_loader` bench.

**L10: bench hardening** — `bench_run_window`/`bench_chunk_loader` use a **100%-variant** synthetic workload (REF+ALT at every position), never exercising the cohort variant-filter pure-REF drop fast path (the dominant real-tomato cost); and `bench_run_window` lacks an assertion that the partition/output is non-empty (an empty partition would time as plausible noise). **Fix:** add a sparse tier (~1/200) and the assertion. *methodology.* **Complexity:** low.

**L11: [driver.rs:487](../../../../src/var_calling/driver.rs#L487) — `BlockQueue` cap = `n_workers*2` is unmeasured**; it trades peak RSS (the scaling thesis's other axis) for worker occupancy. **Fix:** sweep the depth while recording RSS. *concurrency.* **Measurement:** §3 #4.

### Speculative

- **S1:** [posterior_engine.rs:3057-3068](../../../../src/var_calling/posterior_engine.rs#L3057-L3068) — `ln_gamma` recomputed per-`k` in the QUAL Beta-Binomial prior loop (~1.6% self-time); a cumulative `ln(j!)` table cached per run by `max_k` removes it. *hot_loops.* Act only after H1/L4.
- **S2:** [kernels/unify_alleles.rs:799](../../../../src/var_calling/kernels/unify_alleles.rs#L799) — the per-group local `candidates` `BTreeMap` in `detect_compound_candidates_columnar`; not on the named DHAT sites — revisit only if it surfaces after H2/H3. *allocations.*
- **S3:** `DRIVER_PSP_BUFFER_BYTES` (64 KiB) could be swept larger, but loader scaling is linear-in-N (decode-bound, not `read(2)`-bound) and N×buffer RAM cuts against the thesis — likely leave as-is. *io_and_syscalls.*
- **S4:** `sample_size = 10` + a 5% `REGRESSION THRESHOLD` comment is a mismatch (5% is inside run-to-run noise at that sample size); cross-commit deltas need save-baseline pairs or the revert-experiment protocol. *methodology.*

### Note

- **No CI perf gate:** `ci.yml` runs `--lib --tests` and excludes benches, so perf regressions land silently on `main`. (Real infra cost; runner-noise caveats make a hard gate tricky — consider a manual baseline-compare job.)
- **DUST mask buffers** (dust_filter.rs:603 / driver.rs:1014) are 546 MB but only **541 blocks** — they touch the peak-RSS axis, not the churn story; L6 addresses the copy.
- **`build_posterior_record_columnar` clones** (worker.rs:539-541, 487 MB / 266 k blocks) are **legitimate owned output** moved into `PosteriorRecord` — not removable scratch; do not "fix".
- **False sharing: none.** Parallel workers hold private `WorkerScratch`/`WindowRunStats`, rolled up serially by the collector; the only `Arc`-shared atomics (`LhCapStats`) belong to the row-shape `PerGroupMerger` test oracle, not `run_window`.

## 6. Out-of-scope observations

- **`src/vcf/` `build_samples` / `write_record`** (record_encode.rs:534-536, writer.rs:226) — ~1.72 M blocks each in the DHAT run, on the cohort run path. Follow-up: a perf pass on the VCF encoder's per-record allocation.
- **`src/psp/block.rs` `decode_*_column`** (~92 MB ×2) — PSP decode allocation; covered by the PSP perf reviews, re-confirm there.
- **`src/fasta/` fetcher ctor** — L7's root cause (`for_contig` re-parsing `.fai`) lives in `fasta/fetcher.rs`; coordinate the fix with that module's owner.

## 7. What's already good

- **The columnar scratch-reuse design holds on the consume side** — `ColumnarPipelineScratch`/`WorkerScratch` retain capacity across groups/chunks; the DHAT churn is concentrated in two specific structures (BTreeMap value-Vecs, byte keys), not pervasive per-group allocation.
- **The E-step already has a homogeneous-fixation fast path** ([posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs)) that skips the per-genotype `log_sum_exp` work entirely when `fixation_index = 0` (the default) — which is why the hot spot is the QUAL convolution, not the EM.
- **Build config is already tuned** (`lto = "fat"`, `codegen-units = 1`, `panic = "abort"`, pinned `target-cpu` floors) — the methodology's highest-leverage knobs are set; only the production allocator (L1) is left on the table.

### Author responses

- **H1 — applied (2026-06-01, `f920a19`).** Replaced the QUAL convolution's pairwise `log_sum_exp_2` fold with a single `log_sum_exp_slice` over the gathered ≤`ploidy+1` finite terms (stack buffer). **Measured** (`var_calling_run_window` vs the `h1-before` criterion baseline, container): N=8 −0% (noise), **N=64 −16.0%**, **N=200 −21.3%** (p<0.05) — the win grows with N, matching the O(N²) mechanism. **Not byte-identical**: log-sum-exp reassociation shifts QUAL by ≤1 f32 ULP (max rel. delta 1.4e-6) in **71/5062 records (1.4%)** on the tomato fixture — QUAL column only, **no** AF/AC/GT/GQ/AD/FILTER changes, record count unchanged. PM accepted the trade and will rebaseline the out-of-tree byte-identity oracle. 409 lib + 11 cohort-integration tests pass; fmt/clippy clean. (Did **not** include SIMD — that is the separate L4.)
- **H2 — applied (2026-06-01, `65e8c4d`).** `build_chain_proposals_columnar`'s `BTreeMap<ChainId, Vec<CompoundConstituent>>` → `AHashMap<ChainId, SmallVec<[CompoundConstituent; 2]>>` (`BTreeMap::clear()` drops its nodes + value Vecs; `HashMap::clear()` retains buckets, SmallVec keeps the common ≤2-constituent compound inline). **Byte-identity-safe** (caller folds order-independently) — verified md5-identical to post-H1. **Measured** (DHAT, real cohort N=18): total **51.4M → 35.2M blocks (−31.5%)**; the #1 site eliminated. `run_window` synthetic wall unchanged (+0.7%, p=0.22) — an **allocator-pressure** win, not a single-thread wall win.
- **H3 — applied (2026-06-01, `c42f23e`).** `byte_index` dedup map key `Vec<u8>` → `SmallVec<[u8; 16]>` (SNP/short-indel keys stay inline; no per-allele heap key). **Byte-identical** to post-H2. **Measured** (DHAT): **35.2M → 27.86M blocks (−20.8%)**; combined H2+H3 **−45.8%** churn; the remaining top alloc sites are all in out-of-scope `src/vcf/`.
- **Peak-RSS check (H1–H3) — resolved, no regression.** DHAT t-gmax rose 816 → 1224 MB at N=18, which raised a RAM concern. **Verified against real peak RSS** (`/usr/bin/time -l`, profile_cohort_e2e N=50 threads=4, same input, HEAD `e55437b` vs pre-H1 `9f73feb`): **2.303 GB (pre-H1) → 2.286 GB (HEAD), −0.7 % — flat within noise.** The DHAT t-gmax difference was a **timing artifact** of the instrumented parallel run, not a real RSS increase: the `BlockQueue` is bounded (`cap = 2·n_workers`), so the producer cannot run far ahead however fast consume gets. So H1–H3 are clean wins with **no peak-RSS cost**, and **M13** (bound the result/recycle channels) is *not* needed to offset them — it stays a standalone latent-robustness item, untriggered here.

- **L4 — closed (won't do): target went cold after H1.** L4 proposed SIMD-vectorizing the QUAL convolution `k`-loop. A re-profile of `run_window` N=200 at HEAD `e55437b` (macOS `sample`) shows the convolution (posterior_engine.rs:303x) now holds **0 sampled frames** — H1 reduced it from the #1 self-time site to cold. The hot region has moved to the **EM E-step** (`e_step_simd`, posterior_engine.rs:2704/2723 — `log_sum_exp_slice_x4` + `exp_x4`), which is **already SIMD-vectorized**. So L4 would optimize cold code (predictably no gain); further E-step wins would be algorithmic (fewer iterations / SQUAREM / `log_indep` reuse), not a simple SIMD pass. Profiling before implementing avoided the wasted effort — exactly the skill's profile-first principle.

### Author response convention

Address each finding by id (`H1`, `L5`, …) with `applied in <commit>` / `experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. The "no gain" path is expected and welcome — that is what the §3 measurement plan is for. Apply H1/H2/H3 one at a time, each gated on its own measurement.

### Commands to re-verify

- Baseline benches (container): `./scripts/dev.sh cargo bench --bench var_calling_perf -- "var_calling_chunk_loader|var_calling_run_window"`
- CPU profile (host, native build): `cargo bench --bench var_calling_perf --no-run` then `sample <pid> 20 -file out.txt` against `… --bench --profile-time 30 'var_calling_run_window/snp_dense_200'`
- Allocation profile (host): `cargo run --release --example dhat_var_calling --features dhat-heap -- --psp-dir benchmarks/tomato1/results/ours/cohort/psp --n-samples 18 --reference ~/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa --output tmp/x.vcf --threads 4`
- Codegen check (host): `cargo asm --lib --simplify "pop_var_caller::var_calling::posterior_engine::compute_qual_via_exact_af"`

Per-category audit trail (incl. the parse scripts, the `sample` call tree, the samply `.json.gz`, and `dhat-heap.json`): `tmp/perf_review_2026-06-01_var_calling/`.
