# Performance Review: psp_to_vcf
**Date:** 2026-05-20
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** `.psp` → cohort-VCF pipeline (Stages 3–6, `pop_var_caller var-calling`)
**Verdict:** Apply the listed wins
**Hot-path evidence:** `perf record -F 997 -g --call-graph=dwarf` on a real tomato (SL4.0) cohort, T=1, 11.93 s wall, 14 K P-core samples; thread-scaling sweep at T=1/2/4/16; raw artefacts in `tmp/perf_review_2026-05-20_psp_to_vcf/` (perf.data, samply.json.gz, per-category findings)

---

## 1. Scope and constraints

- **Reviewed:** the cohort `.psp → VCF` pipeline (Stages 3–6) as wired by [src/pop_var_caller/var_calling.rs](../../../src/pop_var_caller/var_calling.rs) and [src/pop_var_caller/cohort_driver.rs](../../../src/pop_var_caller/cohort_driver.rs). Stages: DUST filter, k-way per-position merger, variant grouper, per-group merger, posterior engine, VCF writer. Bench infrastructure ([benches/cohort_e2e_perf.rs](../../../benches/cohort_e2e_perf.rs), [examples/profile_cohort_e2e.rs](../../../examples/profile_cohort_e2e.rs)) reviewed under `methodology`.
- **Reviewed against:** branch `main`, last commit `d4252da` (posterior_engine threshold relaxation 1e-4 → 1e-3, applied during the baseline pass that preceded this review).
- **Throughput / latency targets, expected input sizes, target hardware:**
  - **Production target:** x86_64 server (per memory `project_deployment_targets`); secondary Apple Silicon M5 (aarch64/NEON, dev-only).
  - **Real input:** tomato (SL4.0) `SRR7279725.small.psp` — 1 942 095 records across 13 chromosomes, 11.5 MB compressed; replicated N times to synthesise a cohort.
  - **Today's wall time:** 11.93 s at T=1 / N=10 on a 16-core hybrid CPU (Intel P+E cores) under `powersave` governor. Pipeline is **essentially single-threaded** — see Hot-path evidence below.
  - **Goal:** reduce end-to-end wall time on real data. Headline opportunity is parallelism (16 cores ~94 % idle today), with a substantial single-thread component in DUST and the per-record allocation churn.
- **Hot-path evidence available:**
  - Sampling profile (perf record, DWARF call-graph, 14 K P-core samples, 11.93 s).
  - Thread-scaling sweep on real data (T=1=12.7s, T=2=11.5s, T=4=12.0s, T=16=13.6s).
  - Prior per-stage reviews ([perf_var_calling_2026-05-16.md](perf_var_calling_2026-05-16.md), [perf_posterior_engine_2026-05-18.md](perf_posterior_engine_2026-05-18.md), [perf_psp_reader_2026-05-13.md](perf_psp_reader_2026-05-13.md), [perf_psp_writer_2026-05-13.md](perf_psp_writer_2026-05-13.md), [perf_baq_2026-05-12.md](perf_baq_2026-05-12.md)).
  - End-to-end criterion bench `cohort_e2e_perf` (synthetic, [benches/cohort_e2e_perf.rs](../../../benches/cohort_e2e_perf.rs)).
- **In-scope files:**
  - [src/var_calling/dust_filter.rs](../../../src/var_calling/dust_filter.rs)
  - [src/var_calling/per_position_merger.rs](../../../src/var_calling/per_position_merger.rs)
  - [src/per_sample_pileup/psp/reader.rs](../../../src/per_sample_pileup/psp/reader.rs)
  - [src/per_sample_pileup/psp/block.rs](../../../src/per_sample_pileup/psp/block.rs)
  - [src/per_sample_pileup/ref_fetcher.rs](../../../src/per_sample_pileup/ref_fetcher.rs)
  - [src/per_sample_pileup/pileup/mod.rs](../../../src/per_sample_pileup/pileup/mod.rs)
  - [src/pop_var_caller/cohort_driver.rs](../../../src/pop_var_caller/cohort_driver.rs)
  - [src/pop_var_caller/var_calling.rs](../../../src/pop_var_caller/var_calling.rs)
  - [src/pop_var_caller/common.rs](../../../src/pop_var_caller/common.rs)
  - [Cargo.toml](../../../Cargo.toml)
  - [benches/cohort_e2e_perf.rs](../../../benches/cohort_e2e_perf.rs)
  - [benches/var_calling_perf.rs](../../../benches/var_calling_perf.rs)
  - [examples/profile_cohort_e2e.rs](../../../examples/profile_cohort_e2e.rs)
- **Deliberately out of scope:**
  - **Per-group merger algorithmic core** + **posterior engine inner math** — both account for <2 % combined on real data despite the synthetic bench's >50 %; recently reviewed in [perf_var_calling_2026-05-16.md](perf_var_calling_2026-05-16.md) and [perf_posterior_engine_2026-05-18.md](perf_posterior_engine_2026-05-18.md). The synthetic-vs-real cost-distribution divergence is itself a finding (see L10).
  - **Stage 1 walker, BAQ, CRAM input** — different command path, not part of `var-calling`.
  - **VCF writer internals** — not in the top hot symbols; the [io_and_syscalls sub-agent](../../../tmp/perf_review_2026-05-20_psp_to_vcf/io_and_syscalls.md) confirmed the writer's flush/sync/rename discipline is fine.
- **Categories dispatched** (all 6): methodology (build + bench hygiene; DUST has no isolated bench), allocations (~21 % allocator self-time), data_layout (per-record memcpy + Vec-of-Vec), concurrency (sequential pull-iterator chain), hot_loops (33 % in DUST inner loop), io_and_syscalls (PSP block decode + M5 verify).

## 2. Verdict

**Apply the listed wins.** Multiple candidates are well-evidenced (DUST inner-loop `Vec::insert`, PSP reader BufReader-defeating seeks, M5 verify whole-contig allocation, per-chromosome parallelism, per-position merger Vec-alloc-per-emit). Methodology gaps (no isolated DUST bench, no allocator A/B) should land **before or alongside** the code changes they would validate — otherwise the listed code findings cannot be measured cleanly.

The single biggest lever by far is **per-chromosome parallelism (H1)**: the profile + thread sweep together show the pipeline is structurally single-threaded, the workload has 13 independent chromosomes, and the host has 16 cores. Theoretical max ~13× from this one change. Every other hot finding either fits underneath that lever (each per-chrom worker still runs DUST + per-record alloc) or is a serial-time reduction that compounds with it.

## 3. Measurement plan

In the order they unblock other findings:

1. **Add the DUST criterion bench** (H7). Synthetic `pure_ref_<N>_positions` + `mixed_<N>_positions` workloads in [benches/var_calling_perf.rs](../../../benches/var_calling_perf.rs), mirroring the SNP-cadence shape in `cohort_e2e_perf.rs`. Gates every subsequent DUST code change (H2, H3, L1, L2). Threshold: median time on `pure_ref_<N>` improves by ≥ 5 % at within-run CI ≤ 2 %, confirmed by revert.
2. **Flip `[profile.release] debug = "line-tables-only"` → `debug = true`** (L8) OR add a dedicated `[profile.bench-with-debug]`. Re-record `perf record --call-graph=dwarf`; the existing profile may have collapsed inlined frames under `lto = "fat"`. Diff the new flamegraph at the `DustFilter::next` subtree against `tmp/perf_review_2026-05-20_psp_to_vcf/perf.data`.
3. **Run the allocator A/B** (L9). Two back-to-back runs of `cohort_e2e_core/scaling_samples` at `N=64`: one with system glibc (`--save-baseline glibc`), one with `--features alloc-mimalloc --baseline glibc`. Threshold: ≥ 5 % wall-time improvement → flip `default = ["alloc-mimalloc"]`. Independent of any code change.
4. **Apply H2/H3** (DUST `find_perfect` `Vec::insert` + `VecDeque` indexing) and bench against (1). Single-PR change.
5. **Apply H4** (PSP `SeekFrom::Start` → `seek_relative`). Independent. Measure: `strace -c -e file,read,lseek` before/after for syscall counts; `samply` for the `RecordsIter::next` + `_raw_spin_lock` cluster.
6. **Apply H5** (M5 verify streaming). Independent. `md5::compress::compress` should drop out of the top symbols. Threshold: ≥ 0.5 s wall reduction.
7. **Apply H6** (per-position merger `vec![None; n_samples]` per emit). Pairs with the lending-iterator pivot in L3.
8. **Apply H1** (per-chromosome parallelism). Headline lever. Bench against `cohort_e2e_full/scaling_threads` over T=1, 2, 4, 8, 13, 16 — expected shape is roughly linear up to min(13, T), then flat. Once it lands, L1 (drop per-group `par_iter`) becomes mandatory (nested rayon is wasteful) and L5 (`SyncRefFetcher` pre-warm) becomes the next ceiling.
9. **Likely tier follow-ups** in any order based on benches: L2–L7, L10–L12.

## 4. Build / toolchain configuration

- **L8** (Likely): `[profile.release] debug = "line-tables-only"` in [Cargo.toml](../../../Cargo.toml). With `lto = "fat"` + `codegen-units = 1`, inlined frames may collapse under DWARF unwinding — `samply` / `perf` cannot then attribute time to inlined helpers of `DustFilter::next`. Fix: `debug = true` (full) or a dedicated `[profile.bench-with-debug]` profile.
- **L9** (Likely): No `alloc-mimalloc` A/B against real-data workload despite ~21 % allocator self-time. The feature exists ([Cargo.toml:89](../../../Cargo.toml#L89)) and the bench shims wire it up; the experiment is two `cargo bench` runs. If positive, `default = ["alloc-mimalloc"]`.
- **S1** (Speculative): No PGO baseline. The workload is the canonical PGO-friendly shape (same hot loops, same input shape). Defer until after H1/L9 land — *one change per measurement*.

## 5. Code-level findings

Grouped by severity. Within each, ordered by file path. Severity codes are dense and stable (referenced from PROJECT_STATUS Open lists and from any subsequent fixes-applied report).

### Hot-path

- **H1**: [src/pop_var_caller/var_calling.rs:230-367](../../../src/pop_var_caller/var_calling.rs#L230) + [src/pop_var_caller/cohort_driver.rs:99-131](../../../src/pop_var_caller/cohort_driver.rs#L99) — **[Hot-path]** Single sequential pull-iterator chain has no parallelism seam; the 13-chromosome workload + 16-core host is the natural decomposition
  - **Confidence:** High
  - **Hot-path evidence:** thread-scaling sweep on real tomato `.psp` cohort (N=10):
    ```
    T   Median wall
    1   12.72 s
    2   11.54 s
    4   12.01 s
    16  13.60 s (slower than T=1)
    ```
    `perf` shows the chain runs on the caller thread (no parallel work above `Iterator::next` between stages); `rayon::iter::plumbing::bridge_producer_consumer::helper` at 0.82 % is the per-group merger's intra-batch `par_iter` overhead, with <2 % useful work to amortise it against.
  - **Pattern matched:** *"Rayon's default splitting assumes roughly uniform per-item cost"* (negated — there is no `par_iter` over the chain itself, so adding threads cannot reduce wall time).
  - **Mechanism:** `run_var_calling` builds one `PerPositionMerger` over N readers, wraps it in `DustFilter → VariantGrouper → PerGroupMerger → PosteriorEngine → CohortVcfWriter`. Whole chain runs single-threaded. The writer's monotonicity constraint is **per-contig** (verified at [vcf_writer/writer.rs:71-73, 199-210](../../../src/var_calling/vcf_writer/writer.rs#L71)) — **not cross-contig** — so the 13 chromosomes can be processed end-to-end in parallel, with per-chrom outputs concatenated in contig-table order at the writer.
  - **Measurement plan:** Factor `drive_cohort_pipeline`'s DUST→…→write body into a `process_chromosome(chrom_id, per_chrom_readers) → Vec<PosteriorRecord>` helper; build one merger per chromosome over per-chrom sub-iterators of each `PspReader`; drive via `chromosomes.par_iter().map(process_chromosome).collect()`; walk the resulting `Vec<Vec<PosteriorRecord>>` in contig-table order and `writer.write_record` each. Bench: `cohort_e2e_full/scaling_threads` at T=1, 2, 4, 8, 13, 16 on real tomato. Threshold: ≥ 4× wall-time reduction at T=8 on the 13-chromosome input — below 2× the per-chrom-buffer memory cost is not justified.
  - **Complexity cost:** Either (a) per-chrom `Vec<PosteriorRecord>` buffers (memory ∝ records-per-largest-chrom; for SL4.0ch01 ~90 Mbp this is well below 100 MB at realistic SNP densities) or (b) 13 bounded crossbeam-channels + a dedicated writer thread draining channel 0 to exhaustion, then channel 1, etc. PSP reader needs `records_for_chrom(chrom_id)` (the `.psp` format has per-chromosome block boundaries — single seek per chrom). Per-chrom open multiplies open-file count by `13 × N_samples` (130 fds at N=10); within default ulimit but a `RLIMIT_NOFILE` bump is prudent in the binary's startup. The MD5 pre-warm in [common.rs:183-208](../../../src/pop_var_caller/common.rs#L183) already serialises every contig into `SyncRefFetcher`'s cache once, so reference is not a per-chrom contention surface (see L5 for what happens next).
  - **Suggested experiment / fix:** See [`concurrency.md` finding H1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/concurrency.md) for the numbered 5-step spike.

- **H2**: [src/var_calling/dust_filter.rs:353](../../../src/var_calling/dust_filter.rs#L353) — **[Hot-path]** `self.perf.insert(j, …)` inside `find_perfect` is an O(perf.len() − j) memmove per insert on the dominant hot path
  - **Confidence:** High
  - **Hot-path evidence:** call-graph at `tmp/perf_review_2026-05-20_psp_to_vcf/perf_callgraph_top3pct.txt`:
    ```
    33.10%  [.] DustFilter::next
         |--20.20%--find_perfect (inlined)
    ```
    Per-symbol report names `__memmove_avx_unaligned_erms` at 3.36 % — exactly the routine `Vec::insert` lowers to for the shift step.
  - **Pattern matched:** *"Bounds checks: avoid by structure"* / "`Vec::insert` in the middle of a vector is the canonical O(n) hot-loop wart".
  - **Mechanism:** `self.perf` is maintained sorted descending by `start`. `find_perfect` walks `i` right-to-left and inserts in descending-start order; on low-complexity tracts (homopolymer / dinuc / repeat) `perf` accumulates many candidates per window before `save_masked_regions` drains them, and every insert in the middle pays a `ptr::copy` to shift the tail. Flipping to ascending-start order + `Vec::push` + a sweep index removes the memmove entirely; `save_masked_regions` reads `perf.first()` / `swap_remove(0)` instead of `last()`.
  - **Measurement plan:** Swap layout (cheapest: ascending + sweep index). Re-run `samply` / `perf record` against the same real-data workload at T=1; expect `find_perfect` self-time to drop and `__memmove_avx_unaligned_erms` to fall out of the top 10. Threshold: ≥ 2 pp drop in DUST self-time. Gated by H7 (DUST bench).
  - **Complexity cost:** Ascending order requires flipping callers in `save_masked_regions` and the sweep-index discipline. ~30 lines, no `unsafe`, golden tests stay identical.
  - **Suggested experiment / fix:** See [`hot_loops.md` H1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/hot_loops.md) for the three considered designs and the recommended ascending+sweep-index variant.

- **H3**: [src/var_calling/dust_filter.rs:300-364](../../../src/var_calling/dust_filter.rs#L300) — **[Hot-path]** `find_perfect` indexes `VecDeque` (`self.window[i as usize]`) inside its inner loop; each read pays a bounds check + head-relative wrap modulo capacity
  - **Confidence:** High
  - **Hot-path evidence:** same call-graph attribution as H2 — `find_perfect` is 20.20 % of total. Inner body at lines 322-364 reads `self.window[i as usize]` once per iteration.
  - **Pattern matched:** *"Bounds checks: avoid by structure, not by `unsafe`"* + *"Autovectorization needs the compiler's confidence"*.
  - **Mechanism:** `VecDeque::index` performs `(self.head + idx) % cap` + length check per access — two extra branches per element. The loop body is data-dependent (`r += c[t]; c[t] += 1`) so it cannot vectorise even on a flat slice, but the per-element wrap and bounds check are pure overhead. Snapshotting the deque into a contiguous stack buffer (`SmallVec<[u8; 64]>` — typical window cap = 62 entries) at the top of `find_perfect` removes the wrap.
  - **Measurement plan:** Smallest viable change: snapshot via `let w: SmallVec<[u8; 64]> = self.window.iter().copied().collect();` and index `w[i]` instead. Re-bench. If win, do the proper fixed-size ring rewrite of `window`.
  - **Complexity cost:** SmallVec dep (if not already present); 30-40 lines to migrate `window` to a fixed-size ring. No `unsafe`.
  - **Suggested experiment / fix:** [`hot_loops.md` H2 sketch](../../../tmp/perf_review_2026-05-20_psp_to_vcf/hot_loops.md).

- **H4**: [src/per_sample_pileup/psp/reader.rs:585](../../../src/per_sample_pileup/psp/reader.rs#L585) + [src/per_sample_pileup/psp/reader.rs:803](../../../src/per_sample_pileup/psp/reader.rs#L803) — **[Hot-path]** Per-block `SeekFrom::Start` discards the 64 KiB `BufReader` prefetch buffer
  - **Confidence:** High
  - **Hot-path evidence:** profile names `RecordsIter::next` at 4.23 % + `_raw_spin_lock` at 1.41 % (likely `read(2)` kernel contention). Reader's own docstring acknowledges the per-open cost ("five `seek` + four `read_exact` calls just to open"); the per-block pattern repeats the same anti-pattern. `std::io::BufReader::seek` discards its internal buffer on every `SeekFrom::Start` (`std/io/buffered/bufreader.rs:519-524`, `result = self.inner.seek(pos)?; ...; self.discard_buffer();`).
  - **Pattern matched:** *"`File::read` without buffering is a per-call syscall"* — here the buffering exists but the seek pattern silently invalidates it.
  - **Mechanism:** `load_next_block` calls `seek(SeekFrom::Start(entry.block_offset))` on every block; `read_block_header` reads up to 4 KiB into a scratch `Vec`, decodes a typically-tiny header, then `seek(SeekFrom::Start(header_start + consumed))` to rewind past the over-read — both throw away the 4 KiB that the kernel just read.
  - **Measurement plan:** `strace -c -e file,read,lseek` on real data before vs after; compare `read` + `lseek` syscall counts and `% time`. Threshold to merge: `read` count drops by ~`n_blocks × n_samples` and wall drops ≥ 3 %.
  - **Complexity cost:** Small. (a) Replace the post-header rewind `SeekFrom::Start(header_start + consumed)` with `seek_relative(-((buf.len() - consumed) as i64))` — keeps BufReader's prefetch. (b) Gate the pre-block seek with "skip if already at this offset" via a tracked `expected_pos: u64`. No `unsafe`. Optionally: switch the header read to `BufRead::fill_buf` + `consume` (cleanest; removes the over-read entirely).
  - **Suggested experiment / fix:** [`io_and_syscalls.md` H1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/io_and_syscalls.md) gives the three variants in increasing-touch order.

- **H5**: [src/pop_var_caller/common.rs:183-208](../../../src/pop_var_caller/common.rs#L183) + [src/per_sample_pileup/ref_fetcher.rs:144-203](../../../src/per_sample_pileup/ref_fetcher.rs#L144) — **[Hot-path]** M5 verify materialises every contig as a fresh `Vec<u8>` and uppercases byte-by-byte before MD5
  - **Confidence:** High
  - **Hot-path evidence:** `md5::compress::compress` at 8.77 % of total (≈ 1.05 s of the 11.93 s wall at T=1). Tomato SL4.0 is 13 contigs totalling 770 MB — each contig allocates and uppercases its full bytes once. Peak RSS spikes by `max(contig_length)` ≈ 91 MB for the longest tomato chromosome.
  - **Pattern matched:** *"Reading a whole file into memory is the fastest path — when the file fits"* (negated — the materialisation buys nothing the MD5 needs; MD5 is incremental).
  - **Mechanism:** `verify_fasta_matches_psp_chromosomes` calls `fetcher.fetch(chrom_id, 1, contig.length)` per contig; `fetch_from_repository` allocates a brand-new `Vec<u8>` of `contig.length` and walks `.iter().map(|b| b.to_ascii_uppercase()).collect()` (the scalar branchy form does not autovectorise); `Md5::digest(&bytes)` then consumes the freshly allocated Vec.
  - **Measurement plan:** Stream the FASTA in fixed 64 KiB windows, feeding `Md5::update(chunk_upper)` from a reused buffer. `time` the binary before/after; expect `md5::compress::compress` to drop out of the top symbols and wall to drop ≥ 0.5 s on this fixture.
  - **Complexity cost:** Either (a) keep the `RefSeqFetcher` trait shape but loop over 64 KiB chunks (the trait still allocates per call, so this trades one 91 MB alloc for ~1400 × 64 KiB allocs per chromosome — probably worse; defer pending fetcher-trait changes), or (b) bypass `RefSeqFetcher` for the M5 verify and stream the FASTA on disk via `noodles_fasta::io::Reader` (peak RSS is one chunk regardless of contig size). (b) also avoids warming the noodles cache for contigs the cohort may not visit. The cache pre-warm comment in [common.rs:178-181](../../../src/pop_var_caller/common.rs#L178) is honest about the coupling.
  - **Suggested experiment / fix:** [`io_and_syscalls.md` H2 sketch](../../../tmp/perf_review_2026-05-20_psp_to_vcf/io_and_syscalls.md).

- **H6**: [src/var_calling/per_position_merger.rs:306](../../../src/var_calling/per_position_merger.rs#L306) — **[Hot-path]** Per-emission `vec![None; n_samples]` allocates a fresh `Vec<Option<PileupRecord>>` every output position
  - **Confidence:** High
  - **Hot-path evidence:** `PerPositionMerger::next` at 1.26 % self-time + the function emits one item per `(chrom_id, pos)` over 1.94 M records — ~19 M allocations from this site alone at N=10, contributing materially to the ~21 % allocator-family aggregate (`_int_malloc` 4.73 + `malloc` 4.49 + `_int_free_chunk` 3.55 + `cfree` 2.97 + `malloc_consolidate` 1.31 + `unlink_chunk` 0.95). Originally flagged as L8 in [perf_var_calling_2026-05-16.md](perf_var_calling_2026-05-16.md), re-confirmed in current code.
  - **Pattern matched:** *"Allocations belong outside hot loops"* + *"Pre-size containers"* (the size is fixed for the merger's lifetime).
  - **Mechanism:** A fresh `Vec` of length `n_samples` is built on every `next()`; the iterator currently moves `per_sample` into the emitted `PerPositionPileups`. Fix is a lending-iterator pivot (the merger owns a `per_sample_scratch: Vec<Option<PileupRecord>>` reused across emissions; consumer borrows for the duration of one position). This pairs with the same pivot for `DustFilter::next` (the upstream consumer) so the lend chain is continuous.
  - **Measurement plan:** Bench `cohort_e2e_core/scaling_samples` at N=64 before/after. Expect allocator-family aggregate to drop measurably and `PerPositionMerger::next` self-time to fall.
  - **Complexity cost:** API surface change. Two viable redesigns: (a) `'a [Option<PileupRecord>]` borrow tied to `&mut self` — breaks the `Iterator` trait, forces every downstream stage (DustFilter, grouper, …) to be re-typed; (b) `next_into(&mut buf)` companion method keeping `Iterator` as a thin adapter for tests. (b) is the lower-blast-radius path.
  - **Suggested experiment / fix:** [`allocations.md` H1 sketch](../../../tmp/perf_review_2026-05-20_psp_to_vcf/allocations.md).

- **H7**: [benches/var_calling_perf.rs](../../../benches/var_calling_perf.rs) — **[Hot-path]** No isolated criterion bench for `DustFilter::next` despite its being 33.10 % of real-data self-time
  - **Confidence:** High
  - **Hot-path evidence:** the single largest self-time symbol in the 14K-sample P-core profile, larger than the next four entries combined. Existing `var_calling_perf.rs` has four bench groups (merger, grouper, per-group-merger, posterior-engine) but no DUST group. `cohort_e2e_perf.rs` exercises DUST inside the end-to-end pipeline, but a DUST-only ≥ 10 % win there is diluted to ~3 % wall — comparable to the run-to-run noise the file's docstring describes (`feedback_bench_cpu_governor`).
  - **Pattern matched:** Methodology — *"No optimization without a benchmark"*.
  - **Mechanism:** Without an isolated bench, any code-level recommendation against `DustFilter::next` (H2, H3, L1, L2) lands without a baseline that can validate it before merge.
  - **Measurement plan:** Add `var_calling_dust_filter` group with `pure_ref_<N>_positions` + `mixed_<N>_positions` workloads. Threshold: median time on `pure_ref_<N>` improves by ≥ 5 % at within-run CI ≤ 2 %, confirmed by revert.
  - **Complexity cost:** One bench group + two fixture builders, no new dependencies, ~80 lines mirroring `bench_per_position_merger`.
  - **Suggested experiment / fix:** [`methodology.md` H1 sketch](../../../tmp/perf_review_2026-05-20_psp_to_vcf/methodology.md).

### Likely

- **L1**: [src/var_calling/per_group_merger.rs:594-608](../../../src/var_calling/per_group_merger.rs#L594) — **[Likely]** Per-group `par_iter` over `DEFAULT_BATCH_SIZE` groups is net cost — <1 % useful work + visible 0.82 % rayon overhead
  - **Confidence:** High (on the swap being safe; Likely overall because the gain is small in isolation)
  - **Hot-path evidence:** `rayon::iter::plumbing::bridge_producer_consumer::helper` at 0.82 % + per-group + posterior < 2 % combined.
  - **Pattern matched:** *"`rayon::par_iter` over short work is slower than the serial version"*.
  - **Mechanism:** Swap `batch.into_par_iter()` → `batch.into_iter()` (one line). Drops two `Arc::clone` per refill. Independent of H1, but **mandatory once H1 lands** (nested rayon under per-chrom outer parallelism is wasteful).
  - **Measurement plan:** Bench before/after at T=1 + T=16. Expect T=1 to improve by the rayon-bridge overhead (~0.5–1 %); T=16 to improve more (no longer kicking the global pool for 32-element jobs).
  - **Complexity cost:** Negative — removes a `use rayon::prelude::*` if not used elsewhere in the module + two `Arc::clone` lines.
  - **Suggested experiment / fix:** [`concurrency.md` L1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/concurrency.md).

- **L2**: [src/per_sample_pileup/pileup/mod.rs:377-385](../../../src/per_sample_pileup/pileup/mod.rs#L377) — **[Likely]** `AlleleObservation` carries `Vec<u8>` + `Vec<ChainId>` for what real-data shows are typically 1-byte alleles + short chain lists — every record materialisation does 2 heap allocations regardless of payload size
  - **Confidence:** High
  - **Hot-path evidence:** `_int_malloc` 4.73 + `malloc` 4.49 + `_int_free_chunk` 3.55 + `cfree` 2.97 = 15.74 % of cycles in malloc/free; at 1.94 M records × 10 samples × ~2 alleles per record this site contributes tens of millions of small allocations. Allele dump from real data (provided in the earlier merger-bug diagnosis) shows alleles like `A`, `T`, `AGA` — 1–4 bytes typical.
  - **Pattern matched:** *"Pointer chases fragment the cache. If the elements are uniform and small, store inline"*.
  - **Mechanism:** Swap `Vec<u8>` for `SmallVec<[u8; 15]>` and `Vec<ChainId>` for `SmallVec<[ChainId; 2]>` (`ChainId = u64`). Inlines the 1- and 2-allele common case. `mem::take` of inline storage is a byte copy of the struct, not a heap-vector header swap.
  - **Measurement plan:** Gated by H7. Bench before/after; expect `_int_malloc` cluster to drop, `RecordsIter::next` to drop, `__memmove_avx_unaligned_erms` to drop.
  - **Complexity cost:** New `smallvec` (or `tinyvec`) dependency. `AlleleObservation` grows from ~72 bytes to ~88–96 bytes — still within one cache-line spatial-prefetcher pair on x86_64. `#[non_exhaustive]` + `AlleleObservation::new` insulate external callers from the type change; internal callsites that do `record.alleles.iter()` etc. continue to compile (SmallVec derefs to slice).
  - **Suggested experiment / fix:** [`allocations.md` H2](../../../tmp/perf_review_2026-05-20_psp_to_vcf/allocations.md) + [`data_layout.md` L1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/data_layout.md).

- **L3**: [src/per_sample_pileup/psp/block.rs:380](../../../src/per_sample_pileup/psp/block.rs#L380) + [block.rs:412](../../../src/per_sample_pileup/psp/block.rs#L412) + [block.rs:521](../../../src/per_sample_pileup/psp/block.rs#L521) — **[Likely]** `decode_list_column` / `decode_bytes_split` produce `Vec<Vec<T>>`: one outer allocation per block + one inner `Vec::with_capacity(k)` per entry; writer's `encode_list_column_csr` already does the symmetric flat-layout fast path
  - **Confidence:** High
  - **Hot-path evidence:** `decode_list_column` 2.57 % + `decode_varint_column` 0.51 %; each block of the chain-ids column does `expected_count` inner allocations.
  - **Pattern matched:** *"Allocations belong outside hot loops"* + *"Array-of-Structs vs Struct-of-Arrays is a function of access pattern"*.
  - **Mechanism:** Replace `Vec<Vec<T>>` with a CSR pair `(data: Vec<T>, offsets: Vec<u32>)` — exactly the writer's existing layout. Decoder writes into one growable `data` buffer + one `offsets` buffer; `materialise_next_record` reads slices instead of moving owned inner Vecs. Outer + per-entry allocations collapse to two per block. For the LE-fixed-width fast path (`ChainId = u64`), a `decode_list_column_pod` using `bytemuck::cast_slice` is a single memcpy per row (mirror of `decode_scalar_column_pod`).
  - **Measurement plan:** Gated by L9 + H7. Bench `decode_list_column` directly + `cohort_e2e_perf`.
  - **Complexity cost:** New `decode_list_column_csr` (~25 lines) + `decode_list_column_pod` (~25 lines) + change to `DecodedBlock::allele_chain_ids` / `allele_seqs` field types + consumer in `materialise_next_record`. Symmetry with the writer's existing CSR keeps the design memorable.
  - **Suggested experiment / fix:** [`allocations.md` L2/L3](../../../tmp/perf_review_2026-05-20_psp_to_vcf/allocations.md) + [`hot_loops.md` L3](../../../tmp/perf_review_2026-05-20_psp_to_vcf/hot_loops.md) + [`data_layout.md` L3](../../../tmp/perf_review_2026-05-20_psp_to_vcf/data_layout.md).

- **L4**: [src/per_sample_pileup/ref_fetcher.rs:199-202](../../../src/per_sample_pileup/ref_fetcher.rs#L199) — **[Likely]** `fetch_from_repository` allocates a fresh `Vec<u8>` per call and walks `iter().map(to_ascii_uppercase).collect()` per byte — blocks vectorisation, no pre-sized capacity (log-N reallocations per chromosome load)
  - **Confidence:** Medium
  - **Hot-path evidence:** `fetch_from_repository` 0.69 % direct; called from DUST per chrom-change (≤13 times on tomato) and from per-group merger per group (tens of thousands of times). `__memmove_avx_unaligned_erms` 3.36 % is consistent with the slice copy.
  - **Pattern matched:** *"Pre-size containers when the size is known"* + branchless uppercase via `slice::make_ascii_uppercase` (stdlib autovectorises).
  - **Mechanism:** Minimal local fix: `let mut out = bytes[start..end].to_vec(); out.make_ascii_uppercase(); Ok(out)`. Bigger fix (L5): pre-warm into per-chrom `Arc<Vec<u8>>` cache, pre-uppercase at warm-time, drop the noodles `Repository` dependency for runtime fetches.
  - **Measurement plan:** Apply the minimal fix first; bench. If residual still meaningful, pair with L5.
  - **Complexity cost:** Trivial 2-line edit.
  - **Suggested experiment / fix:** [`hot_loops.md` L4](../../../tmp/perf_review_2026-05-20_psp_to_vcf/hot_loops.md) + [`io_and_syscalls.md` L1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/io_and_syscalls.md).

- **L5**: [src/per_sample_pileup/ref_fetcher.rs:109-138](../../../src/per_sample_pileup/ref_fetcher.rs#L109) — **[Likely]** `SyncRefFetcher` forwards through `noodles_fasta::Repository`, which uses `RwLock<HashMap>` per fetch — contention ceiling under per-chrom parallel workers (H1)
  - **Confidence:** High (on the lock confirmed by reading `~/.cargo/registry/src/.../noodles-fasta-0.61.0/src/repository.rs`; Likely on the contention impact because no parallel data exists yet)
  - **Hot-path evidence:** pattern-match on read-lock acquisition per fetch. Cache lookup is `O(1)` in the `HashMap` but gated by `RwLock::read()`; under 16 rayon workers each fetching the same contig's bytes thousands of times per second the read-lock cache line bounces. Today this is silent (single reader), but H1 surfaces it immediately.
  - **Pattern matched:** *"`RwLock` is not a free upgrade from `Mutex`"* + *"`Arc::clone` is one atomic increment, not free"*.
  - **Mechanism:** Stop forwarding to `Repository::get` at runtime. At construction, pre-warm into `bases: Vec<Arc<Vec<u8>>>` indexed by `chrom_id` (the contig table already declares every chrom the run will touch). `fetch` becomes `Ok(self.bases[chrom_id].as_ref()[start..end].to_vec())` — zero locking, zero uppercase work (uppercased at warm-time). Memory shape is unchanged (the existing M5 verify already pre-warms every contig).
  - **Measurement plan:** (1) Instrument current `SyncRefFetcher::fetch` with a `CachePadded<AtomicU64>` call-counter, run H1's per-chrom parallel variant at T=8, confirm read-lock is the suspect. (2) Apply the pre-warm; re-bench. Threshold: under H1, the parallel speedup ceiling lifts measurably (close to ~13× theoretical with the pre-warm vs. lower with the lock).
  - **Complexity cost:** Drops `noodles_fasta::Repository` runtime dependency (still used at startup for the one-shot read). Loses the repo's grow-on-demand cache — but the `.psp` workflow already requires every contig in the header up front, so the capability is unused. M5 verify (common.rs) is then re-shaped: the MD5 of cached bytes can be taken inside `SyncRefFetcher::new` and surfaced via `verify_md5(...)` — the "pre-warm as a side-effect of M5" coupling cleans up. Net: ~50 lines + an init-time MD5 sweep instead of a fetch-time loop.
  - **Suggested experiment / fix:** [`concurrency.md` Likely-2 sketch](../../../tmp/perf_review_2026-05-20_psp_to_vcf/concurrency.md).

- **L6**: [src/per_sample_pileup/psp/reader.rs:148-225](../../../src/per_sample_pileup/psp/reader.rs#L148) — **[Likely]** Open-sequence issues 5 seeks + 4 reads per `.psp`; at N=10 that's 50 seeks + 40 reads at startup; whole-file `mmap` sidesteps both this and H4
  - **Confidence:** Medium
  - **Hot-path evidence:** docstring explicitly lists the open-time pattern. `_raw_spin_lock` 1.41 % plausibly includes a fraction.
  - **Pattern matched:** *"`File::read` without buffering is a per-call syscall"* — buffering is present (64 KiB), but the open-sequence intentionally seeks to disparate file regions (file end for trailer, then `trailer.index_offset`, then file start for header), defeating prefetch.
  - **Mechanism:** Behind a CLI flag, try `PspReader::new(Cursor::new(&mmap[..]))` via `memmap2::Mmap`. Whole reader becomes no-syscall, page-fault-driven. Cross-cuts H4 (with `Cursor<&[u8]>`, `seek` is a single integer assignment with no buffer to discard).
  - **Measurement plan:** Add `--mmap-psp` to the example binary; bench wall-time on real data.
  - **Complexity cost:** `memmap2` dependency; SIGBUS handling on truncation (rare); page-cache pressure (several GB of address space at N=10, negligible on x86_64 server). Feature-flag the path so non-mmap consumers (tests, embedded) keep working.
  - **Suggested experiment / fix:** [`io_and_syscalls.md` L2](../../../tmp/perf_review_2026-05-20_psp_to_vcf/io_and_syscalls.md).

- **L7**: [src/var_calling/per_position_merger.rs:152](../../../src/var_calling/per_position_merger.rs#L152) — **[Likely]** `heads: Vec<Option<PileupRecord>>` has the k-way scan loop pull a 32-byte `PileupRecord` per slot just to read 8 bytes of `(chrom_id, pos)` key
  - **Confidence:** Medium
  - **Hot-path evidence:** `PerPositionMerger::next` 1.26 %; the head-scan reads only `(chrom_id, pos)` per head per emission. At N=10 the scan visits 320 bytes for 80 bytes of payload — wastes most of the cache line.
  - **Pattern matched:** *"Project hot fields into a cache struct (data-oriented design)"*.
  - **Mechanism:** Add a parallel `head_keys: Vec<u64>` (packed `(chrom_id, pos)` with `u64::MAX` sentinel) refreshed whenever `heads[i]` is updated. Scan becomes a pure tight loop over 8 bytes per slot, branchless comparison, autovectorisable on x86_64 / aarch64 NEON. Gain scales with cohort N — entry point at N=10, magnified at N=100+.
  - **Measurement plan:** Bench at N=10, 64, 256 before/after. Threshold: drop in `PerPositionMerger::next` self-time at higher N.
  - **Complexity cost:** One extra `Vec<u64>` field + two-field invariant (`head_keys[i]` ↔ `heads[i]`'s key) with a `debug_assert` helper. Touched at every `heads` mutation site.
  - **Suggested experiment / fix:** [`data_layout.md` Likely-4 sketch](../../../tmp/perf_review_2026-05-20_psp_to_vcf/data_layout.md).

- **L8** (build): [Cargo.toml:38](../../../Cargo.toml#L38) — **[Likely]** `[profile.release] debug = "line-tables-only"` may collapse inlined frames under DWARF unwinding given `lto = "fat"` + `codegen-units = 1`
  - **Confidence:** High
  - **Hot-path evidence:** the existing perf.data was useful at the leaf-symbol level, but a code-level finding that needed to attribute time to an inlined helper of `DustFilter::next` would not resolve it.
  - **Pattern matched:** Methodology — *"Profiles must come from release builds with debug info"*.
  - **Measurement plan:** Flip to `debug = true` or add `[profile.bench-with-debug] inherits = "release", debug = true`; re-record `perf record --call-graph=dwarf`; diff the `DustFilter::next` subtree's children.
  - **Complexity cost:** Larger release binary on disk (debug sections); no runtime cost. If the dedicated-profile route is taken, update `examples/profile_cohort_e2e.rs:23-46` build instructions.
  - **Suggested experiment / fix:** [`methodology.md` Likely-1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/methodology.md).

- **L9** (build): [Cargo.toml:89](../../../Cargo.toml#L89) — **[Likely]** No `alloc-mimalloc` A/B run against real-data workload despite ~21 % allocator self-time + `_raw_spin_lock` kernel symbol
  - **Confidence:** High
  - **Hot-path evidence:** allocator family sums to ~21 % in the profile; the feature exists and is wired correctly in both bench shims but no comparable mimalloc baseline cited.
  - **Pattern matched:** Methodology — *"Allocator choice is a one-line change"*.
  - **Measurement plan:** `cargo bench --bench cohort_e2e_perf -- cohort_e2e_core/scaling_samples --save-baseline glibc`, then `cargo bench --bench cohort_e2e_perf --features alloc-mimalloc -- cohort_e2e_core/scaling_samples --baseline glibc`. Threshold: ≥ 5 % at the N=64 cell → `default = ["alloc-mimalloc"]`.
  - **Complexity cost:** None beyond the existing feature flag.
  - **Suggested experiment / fix:** [`methodology.md` Likely-2](../../../tmp/perf_review_2026-05-20_psp_to_vcf/methodology.md).

- **L10** (bench hygiene): [benches/var_calling_perf.rs (multiple sites)](../../../benches/var_calling_perf.rs) + [benches/cohort_e2e_perf.rs:609](../../../benches/cohort_e2e_perf.rs#L609) — **[Likely]** Drained-count is computed in every bench body but never asserted on; silent-failure risk under refactor
  - **Confidence:** High
  - **Hot-path evidence:** ~9 sites in `var_calling_perf.rs` and 1 in `cohort_e2e_perf.rs` follow the `let mut count = 0; … black_box(count)` pattern with no `assert_eq!` against the expected row count.
  - **Mechanism:** A future refactor that silently changes which rows reach the timed body keeps producing reasonable-looking numbers. The bench fixtures encode load-bearing claims (`build_dense_per_sample_streams` claims every position has a record in every sample, etc.) — assertions surface drift.
  - **Measurement plan:** Add `assert_eq!(count, expected)` at the bottom of each `b.iter_batched` body. Cheap.
  - **Complexity cost:** Trivial — one `assert_eq!` per bench function.
  - **Suggested experiment / fix:** [`methodology.md` Likely-3](../../../tmp/perf_review_2026-05-20_psp_to_vcf/methodology.md).

- **L11** (bench hygiene): [benches/cohort_e2e_perf.rs:558-571](../../../benches/cohort_e2e_perf.rs#L558) — **[Likely]** `cohort_e2e_core/scaling_threads` builds rayon pools inside the `for t in thread_counts` loop without explicit `drop`; pools coexist across iterations until end-of-scope
  - **Confidence:** Medium
  - **Mechanism:** At the `max_threads = 16` cell on a 16-core hybrid CPU, the prior t=1, t=2, t=4 pools are still in scope (threads sleep but hold stacks/TLS) when the t=16 timing runs. Adds noise to the first-sample reading.
  - **Measurement plan:** Add explicit `drop(pool)` at the bottom of the for-loop body.
  - **Complexity cost:** Trivial.
  - **Suggested experiment / fix:** [`methodology.md` Likely-4](../../../tmp/perf_review_2026-05-20_psp_to_vcf/methodology.md).

- **L12** (DUST scratch): [src/var_calling/dust_filter.rs:301-302](../../../src/var_calling/dust_filter.rs#L301) — **[Likely]** `let mut c = self.cv; let mut r = self.rv;` copies 256 bytes of state on every `find_perfect` entry
  - **Confidence:** Medium
  - **Hot-path evidence:** `find_perfect` is 20 % of total runtime. `cv` is `[u32; SD_WTOT]` = `[u32; 64]` = 256 bytes per call.
  - **Mechanism:** Delta-revert: track touched indices in a `SmallVec<[u8; 32]>`, mutate `self.cv` directly, decrement touched entries at end. Skips the 256-byte copy per call.
  - **Measurement plan:** Gated by H7. Re-bench `var_calling_dust_filter`. Threshold: ≥ 0.5 pp drop in `find_perfect` self-time.
  - **Complexity cost:** Adds one scratch field + a revert pass, ~15 lines.
  - **Suggested experiment / fix:** [`hot_loops.md` Likely-1 sketch](../../../tmp/perf_review_2026-05-20_psp_to_vcf/hot_loops.md).

### Speculative

- **S1** (build): No PGO baseline despite the workload's stable-critical-path shape. Defer until after H1 + L9 land — *one change per measurement*. ([methodology.md Speculative-2](../../../tmp/perf_review_2026-05-20_psp_to_vcf/methodology.md))

- **S2** (per-fetch borrowed return): [src/per_sample_pileup/ref_fetcher.rs](../../../src/per_sample_pileup/ref_fetcher.rs) — change `RefSeqFetcher::fetch` to return `&[u8]` (lifetime tied to `&self`) or `Arc<Vec<u8>>` instead of `Vec<u8>`. Removes per-call slice allocation. Wide ripple through DustFilter / PerGroupMerger / test mocks. Defer; the L4/L5 changes capture most of the win without the API churn. ([io_and_syscalls.md Likely-1 complexity-cost discussion](../../../tmp/perf_review_2026-05-20_psp_to_vcf/io_and_syscalls.md))

- **S3** (lending-iterator pivot for Stage 3+): The H6 fix is one half of a larger pivot — every iterator from `PerPositionMerger` through `DustFilter` to the grouper could emit a borrowed `&PerPositionPileups` from a single owned scratch. Stops being `std::iter::Iterator` (no GAT-free lending iterator in stable). Bigger touch than H6 alone. Cross-cuts H6 + L3. ([data_layout.md Hot-path-1 sketch](../../../tmp/perf_review_2026-05-20_psp_to_vcf/data_layout.md))

- **S4** (PSP read_block_header zero-fill): [src/per_sample_pileup/psp/reader.rs:773-823](../../../src/per_sample_pileup/psp/reader.rs#L773) — `Vec::resize(prev + chunk, 0)` zero-fills ~4 KiB per block just to feed `read(...)` into the spare capacity. Replaceable with `Vec::spare_capacity_mut` + `Read::read` into uninit — adds one `unsafe` block. Profile does not name `read_block_header`; defer. ([io_and_syscalls.md Speculative-1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/io_and_syscalls.md))

- **S5** (example wall_time hygiene): [examples/profile_cohort_e2e.rs:176](../../../examples/profile_cohort_e2e.rs#L176) — single `wall_time` brackets `run_var_calling`, includes 8.77 % startup MD5. Split into `wall_time_full` + `wall_time_drive` (mirror `cohort_e2e_perf.rs`'s `run_core_iters`). ([methodology.md Speculative-1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/methodology.md))

### Note

- **N1**: `AlleleSupportStats` field order is fine — `repr(Rust)` reorders to 24 bytes; existing analysis included for the record so future reviewers don't re-investigate. ([data_layout.md Speculative-1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/data_layout.md))
- **N2**: VCF write side is fine — `BufWriter::with_capacity(64 KiB, file)` + once-per-finish `flush` / `sync_all` / `rename` / `sync_parent_dir`. No write-side syscall churn. ([io_and_syscalls.md Note-1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/io_and_syscalls.md))
- **N3** (forward-looking): under H1 (per-chrom parallel), the per-thread `compressed_scratch` / `decompressed_scratch` and per-merger `heads` arrays become per-thread state. If shared on adjacent cache lines they will false-share. Wrap each thread's scratch in `crossbeam_utils::CachePadded<T>` (128-byte pad on x86_64 for the spatial prefetcher). ([data_layout.md Note-2](../../../tmp/perf_review_2026-05-20_psp_to_vcf/data_layout.md))
- **N4**: DUST's `LoadedChrom` cache is already the right shape for per-chrom workers — per-worker mask, zero cross-worker sharing, *better* hit rate than the current single-threaded chain (zero mask rebuilds under per-chrom partition; current code pays one rebuild per chrom transition). Confirms H1's design. ([concurrency.md Note-1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/concurrency.md))
- **N5** (cold-path clones): three `self.sample_names[sample_idx].clone()` sites on error paths in `per_position_merger.rs`. Cold by construction (errors). `Arc<str>` would eliminate; not worth the ripple. ([allocations.md Likely-4](../../../tmp/perf_review_2026-05-20_psp_to_vcf/allocations.md))
- **N6**: `cohort_e2e_full` cannot sweep threads in one bench invocation — `configure_rayon_pool` / `ThreadPoolBuilder::build_global` is once-per-process. Already documented in the bench file's own docstring. ([methodology.md Note-1](../../../tmp/perf_review_2026-05-20_psp_to_vcf/methodology.md))

## 6. Out-of-scope observations

- **Synthetic-vs-real cost-distribution divergence.** The `cohort_e2e_perf` bench shows per-group merger + posterior engine at >50 % of synthetic cost; real-data profile shows them at <2 %. The synthetic SNP-cadence-of-50 with 30-depth REFs over-weights the per-group + posterior work and under-weights DUST + allocations. Worth adding a `real_data_shape` fixture to the bench (still synthetic, but with the SNP-cadence / depth distribution measured from a real `.psp`). Filed as a follow-up; not blocking.
- **Posterior engine `DidNotConverge` long-term fix.** Commit `d4252da` relaxed the convergence threshold from 1e-4 to 1e-3 to unblock real-data runs. A small fraction of records may still hit the iteration cap with `last_delta > 1e-3`; the right long-term fix is emit-with-flag (emit the record with `EmDiagnostics { converged: false, ... }` and let a downstream filter handle it). Out of scope for this perf review; flagged as a correctness follow-up.
- **Stage 1 walker, BAQ, CRAM input.** Not part of `var-calling`; if a future `pipeline` superseder runs them in one process the `configure_rayon_pool` silent-no-op (Likely concurrency-3 in the sub-agent finding) becomes a footgun.
- **VCF writer per-record `genotype_tables.entry(...).or_insert_with(...).clone()`** at [vcf_writer/writer.rs:212-216](../../../src/var_calling/vcf_writer/writer.rs#L212) — clones a `Vec<Vec<u8>>` on every record. Sub-1 % today; would matter under H1 if per-chrom workers each maintain their own writer fragment. Defer until after H1.

## 7. What's already good

- **PSP block decoder reuses scratch across blocks.** `RecordsIter` correctly clears + reuses `compressed_scratch` / `decompressed_scratch` (L1/L2 from the prior PSP reader review; [reader.rs:1056-1057](../../../src/per_sample_pileup/psp/reader.rs#L1056) and [:563-567](../../../src/per_sample_pileup/psp/reader.rs#L563)). The remaining allocation pressure is in per-allele Vec-of-Vec (L3), not in per-block scratch.
- **VCF writer's tmp-rename + sync_all + sync_parent_dir discipline is honest** ([vcf_writer/sink.rs:69-73](../../../src/var_calling/vcf_writer/sink.rs#L69)). One flush + one sync + one rename + one parent-dir sync, all at finish-time. No write-side syscall churn.
- **DUST's `LoadedChrom` per-chrom cache is the right shape for per-chrom workers** ([dust_filter.rs:585-589, 656-683](../../../src/var_calling/dust_filter.rs#L585)). The single-chrom cache that today rebuilds on chrom transition becomes a zero-miss per-worker cache under H1 — confirms the parallel decomposition.
