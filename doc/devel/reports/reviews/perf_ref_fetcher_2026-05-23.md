# Performance Review: ref_fetcher (post-cohort-migration)
**Date:** 2026-05-23
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** `src/per_sample_pileup/ref_fetcher.rs` + 3 hot-path consumers (DUST, PerGroupMerger, BAQ) after the Step-2 `ChromRefFetcher` migration
**Verdict:** **Apply the listed wins** — one structural fix (H1) is profile-named, has a contained mechanism (Mutex → RefCell), and converges from four independent categories. **Gated on H4** (benches don't compile — must be repaired before any code-level fix can be measured).
**Hot-path evidence:** `perf record --call-graph=dwarf,16384 -F 997` against `examples/profile_cohort_e2e` and `pop_var_caller pileup`, on real tomato (SL4.0) data with 4 host P-cores. 203,245 samples (cohort) + 30,620 samples (Stage 1). Verbatim output in [tmp/perf_review_2026-05-23_ref_fetcher/profile_summary.txt](../../../../tmp/perf_review_2026-05-23_ref_fetcher/profile_summary.txt) and `perf_callgraph_full.txt`.

---

## 1. Scope and constraints

- **What was reviewed:** the reference-FASTA fetcher module + the three hot-path consumers that were rewired to it in the Step-2 migration this morning.
- **Reviewed against:** branch `main` at HEAD; modified-but-uncommitted Step-2 changes in working tree (DustFilter / PerGroupMerger / drive_cohort_pipeline / var_calling_from_bam migrated to `ChromRefFetcher`; `SingleChromLegacyAdapter` deleted).
- **Throughput / latency targets, expected input sizes, target hardware:** tomato SL4.0 reference (795 MB, 13 chroms, contig lengths up to ~85 Mbp). Cohort N=10. Production target: 27.7 s end-to-end on a 13-thread x86_64 server (matches the 2026-05-20 cohort-per-chromosome-parallel result). Today's local 4-thread (i7-1260P P-cores) baseline: **58.85 s wall**.
- **Hot-path evidence available:** sampling profile (perf with dwarf call-graph). Both `cpu_atom` (E-cores) and `cpu_core` (P-cores) events were captured; the report quotes the higher of the two as the primary signal where they agree directionally.
- **In-scope files:**
  - [src/per_sample_pileup/ref_fetcher.rs](../../../../src/per_sample_pileup/ref_fetcher.rs) — fetcher module (2,190 LOC)
  - [src/var_calling/dust_filter.rs](../../../../src/var_calling/dust_filter.rs) — DUST consumer (streams `iter_bases`)
  - [src/var_calling/per_group_merger.rs](../../../../src/var_calling/per_group_merger.rs) — per-group `fetch` consumer; defines `SharedRefFetcher`
  - [src/per_sample_pileup/baq/engine.rs](../../../../src/per_sample_pileup/baq/engine.rs) — BAQ consumer (`ManualEvictChromRefFetcher`)
  - [src/pop_var_caller/cohort_driver.rs](../../../../src/pop_var_caller/cohort_driver.rs) — per-worker fetcher construction
  - [src/pop_var_caller/var_calling_from_bam.rs](../../../../src/pop_var_caller/var_calling_from_bam.rs) — same shape
  - [benches/](../../../../benches/) — measurement surface
  - [Cargo.toml](../../../../Cargo.toml)
- **Deliberately out of scope:**
  - `WalkerLegacyAdapter` and `pileup_to_psp` — deferred to a later migration round per user direction.
  - DUST sdust algorithm internals — covered by 2026-05-20 H2/H3 in [perf_psp_to_vcf_2026-05-20.md](perf_psp_to_vcf_2026-05-20.md).
  - BAQ probaln optimization — covered by [perf_baq_2026-05-12.md](perf_baq_2026-05-12.md).
- **Categories dispatched (all six):**
  - `methodology` — bench/profile/Cargo.toml audit (always).
  - `allocations` — `Vec<u8>` returns, `format!` adapters, `Box<dyn Iterator>`.
  - `data_layout` — `StreamState` layout + per-worker isolation.
  - `concurrency` — `Mutex<StreamState>` per-base lock (profile-named site).
  - `hot_loops` — `ChromRefBaseIter::next` + DUST `Map` adapter dispatch.
  - `io_and_syscalls` — raw `File::read` chunked reads, `.fai` re-parse.

## 2. Verdict

**Apply the listed wins** — H1 (drop the `Mutex<StreamState>`) is profile-named at 4.02 % of cohort wall time on the leaf of a clear DUST call-tree, has a one-paragraph mechanism (`Sync` is not load-bearing for any cohort caller; per-worker ownership is an existing invariant), and the smallest-diff fix is `Mutex<StreamState>` → `RefCell<StreamState>` + relax `SharedRefFetcher` to `Arc<dyn ChromRefFetcher + Send>`. Three of the other Hot-path findings (H2 trait `iter_bases` shape, H3 `fetch` returns `Vec<u8>`, H4 buffer-snapshot in iter) are entangled with H1 and naturally co-land.

**Gate**: every measurement-bearing finding in this review is **blocked by H4** — the cohort-side benches don't compile against the current API (`SyncRefFetcher` references, `CohortPipelineParams.fetcher`/`chromosomes` field-drift, `drive_cohort_pipeline` arg count drift). Repair the benches first; the rest follows.

## 3. Measurement plan

The benchmarks and profiles to add or run, in the order they unblock other findings:

1. **Repair the four broken benches** (H4). Verbatim from `./scripts/dev.sh cargo build --benches --release`:
   ```
   error[E0432]: unresolved import `pop_var_caller::per_sample_pileup::ref_fetcher::SyncRefFetcher`
   error[E0560]: struct `…cohort_driver::CohortPipelineParams` has no field named `fetcher`
   error[E0560]: struct `…cohort_driver::CohortPipelineParams` has no field named `chromosomes`
   error[E0061]: this function takes 7 arguments but 5 arguments were supplied
   error: could not compile `pop_var_caller` (bench "cohort_e2e_perf") due to 9 previous errors
   ```
   Fixed in one PR: `cohort_e2e_perf`, `var_calling_perf`, `baq_perf`, `pileup_walker_scaling`. Recipe in the methodology sub-report. **No code-level fix in this review can be merged with confidence until this compiles** (the Bench CPU governor caveat in memory says wall-clock numbers drift ~2× across sessions).
2. **Per-fetcher microbench** (`benches/ref_fetcher_perf.rs`) — `streaming_iter_bases_full_contig`, `streaming_fetch_per_group_window`, `manual_evict_fetch_then_evict`, `refill_warm_cache`. Localises the H1/H2/H3 wins below the cohort-e2e noise floor. Threshold per finding: ≥ 10 % on the microbench *and* ≥ 5 % on `cohort_e2e_perf::e2e_core/cohort_scaling_threads` before merging.
3. **H1 measurement** — apply the `Mutex → RefCell` patch, rerun `perf record --call-graph=dwarf` on `profile_cohort_e2e --threads 4`. Targets:
   - The `drop_in_place<MutexGuard<StreamState>>` line disappears from the `ChromRefBaseIter::next` call-tree.
   - `ChromRefBaseIter::next` self-time drops from 10.96 % / 7.46 % to ≤ 3 % (cpu_atom) / ≤ 2 % (cpu_core).
   - Wall time at `--threads 4` drops by at least ~1.5 s (the 4.02 % unlock floor) — directional only; back-to-back paired runs needed because of the `powersave` governor caveat.
4. **H2 + H3 follow-on** — once H1 lands, the GAT-monomorphic `iter_bases` (H2) and `fetch -> &[u8]` (H3) are the next levers, in that order. Each gets its own paired bench against the post-H1 baseline.
5. **mimalloc A/B** (L7) — once benches build, run `cargo bench --bench cohort_e2e_perf -- e2e_full/cohort_scaling_threads` paired with and without `--features alloc-mimalloc` on the 13-thread server (production target hardware, not the laptop — the laptop's `powersave` governor + hybrid P/E-core dispatch hides effects of this size). Threshold to add a production-binary `#[global_allocator]`: ≥ 5 % wall improvement, sign-confirmed by the revert experiment.

## 4. Build / toolchain configuration

What's already right (do not touch):

- `[profile.release] lto = "fat"`, `codegen-units = 1`, `panic = "abort"` — set.
- `[profile.bench] inherits = "release"` `debug = true` — set.
- `.cargo/config.toml` `rustflags = ["-C", "target-cpu=x86-64-v3"]` on `target.'cfg(all(target_arch = "x86_64", target_os = "linux"))'` — AVX2 + FMA + BMI2 enabled.
- `alloc-mimalloc` feature wired; `dhat-heap` feature wired.
- 1.95.0 toolchain pin implied by bench output (`info: syncing channel updates for 1.95-x86_64-unknown-linux-gnu`).

Two build-config changes recommended (both also called out in §5 as L6 + L7):

- **Add `[profile.release-with-debug]`** that `inherits = "release"` and sets `debug = true`. Document the recipe (`cargo build --profile release-with-debug --bin pop_var_caller`) — current production profile uses `debug = "line-tables-only"`, which strips the inlined-frame chains a future reviewer needs to reproduce the call-graph in this report.
- **Wire `#[global_allocator] mimalloc::MiMalloc` behind `alloc-mimalloc`** in `src/main.rs`. The feature exists but the production binary doesn't pick it up; glibc allocator symbols sum to ~14.5 % cpu_atom / ~23.2 % cpu_core in the cohort profile. Gate the merge on the L7 A/B measurement.

## 5. Code-level findings

### Hot-path

#### H1: [src/per_sample_pileup/ref_fetcher.rs:117](../../../../src/per_sample_pileup/ref_fetcher.rs#L117) — Drop `Sync` from `SharedRefFetcher`, replace `Mutex<StreamState>` with `RefCell<StreamState>`

- **Confidence:** High
- **Hot-path evidence:** verbatim from [profile_summary.txt](../../../../tmp/perf_review_2026-05-23_ref_fetcher/profile_summary.txt):
  ```
   10.96%   7.46%  ChromRefBaseIter::next
   --9.41%--ChromRefBaseIter::next
            |
            |--4.02%--swap (inlined)
            |        unlock (inlined)
            |        drop<StreamState> (inlined)
            |        drop_in_place<MutexGuard<StreamState>>
            |        ChromRefBaseIter::next  ← called from DUST sdust_mask_streaming
  ```
  Concurrence from four categories: methodology cross-routed it; allocations entangled it with the `Vec<u8>` return (H3); data_layout filed the equivalent "snapshot fields into iter"; hot_loops cross-routed it as the underlying cause of the vtable dispatch cost.
- **Pattern matched:** "Lock acquisition does not belong inside a hot loop" + "Atomic primitives that exist only to satisfy `Sync` on a single-threaded-by-construction site."
- **Mechanism:** `ChromRefBaseIter::next` ([ref_fetcher.rs:770-805](../../../../src/per_sample_pileup/ref_fetcher.rs#L770-L805)) takes `self.fetcher.inner.lock()` on every base — `sdust_mask_streaming` drives this iter once per base across the full contig (~85 Mbp for tomato ch00, ~800 Mbp/worker per N=10 cohort). The mutex exists only to satisfy `Sync` on `SharedRefFetcher = Arc<dyn ChromRefFetcher + Send + Sync>` ([per_group_merger.rs:485](../../../../src/var_calling/per_group_merger.rs#L485)), but three independent facts make the `Sync` requirement non-load-bearing:
  1. Per-chrom workers each construct their own fetcher ([cohort_driver.rs:476](../../../../src/pop_var_caller/cohort_driver.rs#L476), [var_calling_from_bam.rs:425](../../../../src/pop_var_caller/var_calling_from_bam.rs#L425)); the `Arc` never crosses a thread boundary — the comment at [ref_fetcher.rs:113-117](../../../../src/per_sample_pileup/ref_fetcher.rs#L113-L117) already states "each per-chrom worker owns its fetcher outright, so contention is zero by construction."
  2. `process_group` takes `&(dyn ChromRefFetcher + Send + Sync)` ([per_group_merger.rs:664](../../../../src/var_calling/per_group_merger.rs#L664)) and only calls `fetch(&self, …)` — the trait method is `&self`, not `&mut self`, so the call-site can stay identical under a `+ Send` bound.
  3. In-tree precedent: `ManualEvictChromRefFetcher` ([ref_fetcher.rs:1073](../../../../src/per_sample_pileup/ref_fetcher.rs#L1073)) is `Send`-only, takes `&mut self` on `fetch`, and BAQ consumes it per-worker via rayon `map_init` ([baq/engine.rs:111](../../../../src/per_sample_pileup/baq/engine.rs#L111)). The relaxed bound is already the in-tree pattern for the sibling fetcher.
- **Measurement plan:** Step 3 of §3. Threshold: `drop_in_place<MutexGuard<StreamState>>` disappears from the call-tree; cohort wall time at `--threads 4` drops by ≥ 1.5 s back-to-back.
- **Complexity cost:** Three coordinated edits:
  - `Mutex<StreamState>` → `RefCell<StreamState>` ([ref_fetcher.rs:117](../../../../src/per_sample_pileup/ref_fetcher.rs#L117)); the trait method stays `&self` because `RefCell::borrow_mut` is interior-mutability.
  - `SharedRefFetcher = Arc<dyn ChromRefFetcher + Send + Sync>` → `Arc<dyn ChromRefFetcher + Send>` ([per_group_merger.rs:485](../../../../src/var_calling/per_group_merger.rs#L485)).
  - `process_group`'s parameter type ([per_group_merger.rs:664](../../../../src/var_calling/per_group_merger.rs#L664)) drops `+ Sync`.
  No new dependency, no `unsafe`, `WalkerLegacyAdapter` unaffected (it owns its own `Mutex<Option<…>>` and is explicitly out of scope).
- **Suggested experiment / fix:**
  ```rust
  // ref_fetcher.rs:113-117
  pub struct StreamingChromRefFetcher {
      chrom_id: u32,
      contig_name: String,
      fai: ContigFai,
      // No Mutex — workers own their fetcher; SharedRefFetcher
      // no longer requires Sync. RefCell::borrow_mut is a
      // single non-atomic increment+check.
      inner: std::cell::RefCell<StreamState>,
  }

  // ChromRefBaseIter::next — ref_fetcher.rs:770
  fn next(&mut self) -> Option<Self::Item> {
      if self.done || self.next_base > self.fetcher.fai.length { return None; }
      let mut state = self.fetcher.inner.borrow_mut();
      // rest unchanged: same refill check, same indexing.
  }

  // ChromRefBaseIter::Drop — ref_fetcher.rs:808
  impl Drop for ChromRefBaseIter<'_> {
      fn drop(&mut self) {
          let mut state = self.fetcher.inner.borrow_mut();
          state.buf.clear();
          state.buf_start_base = 0;
      }
  }

  // per_group_merger.rs:485
  pub type SharedRefFetcher = Arc<dyn ChromRefFetcher + Send>;
  // per_group_merger.rs:664 (process_group signature)
  ref_fetcher: &(dyn ChromRefFetcher + Send),
  ```

#### H2: [src/var_calling/dust_filter.rs:728-730](../../../../src/var_calling/dust_filter.rs#L728-L730) + [ref_fetcher.rs:589-594](../../../../src/per_sample_pileup/ref_fetcher.rs#L589-L594) — `iter_bases` returns `Box<dyn Iterator>`, defeating devirtualization per base

- **Confidence:** High (profile-named via the underlying `ChromRefBaseIter::next` symbol at 10.96 % / 7.46 %)
- **Hot-path evidence:** the same call-tree quoted in H1. The chain inside `sdust_mask_streaming` is `Map<…>::next` → `Box<dyn Iterator>::next` (vtable) → `ChromRefBaseIter::next`. `DustFilter` is *already* generic over `F: ChromRefFetcher` ([dust_filter.rs:642-646](../../../../src/var_calling/dust_filter.rs#L642-L646)), but the trait's `iter_bases` returns `Box<dyn Iterator>` and re-boxes the type information.
- **Pattern matched:** "Static dispatch in hot loops" — the trait could expose a GAT-monomorphic iter type.
- **Mechanism:** Per byte: vtable call through `Box<dyn Iterator>::next` blocks LLVM from inlining `ChromRefBaseIter::next` into `sdust_mask_streaming`'s state machine. A GAT-typed `BaseIter<'a>` lets the concrete `ChromRefBaseIter` flow into the call site so the loop body inlines and (combined with H1) reduces to plain pointer arithmetic on a borrowed slice.
- **Measurement plan:** Two-step (per the hot_loops sub-report). Step A: change the trait to use a GAT. Step B: drop the `.map(|res| res.map_err(io::Error::other))` adapter by widening `sdust_mask_streaming`'s signature to accept `Result<u8, E: Into<io::Error>>` (or `ChromRefFetchError` directly). Targets: `<Box<dyn Iterator>>::next` disappears from the call-tree; cohort wall improves by ≥ 1 s on top of H1.
- **Complexity cost:** A GAT (`type BaseIter<'a>: Iterator<…> where Self: 'a;`) — mechanical update for every impl, but it does make the trait *non-object-safe* unless `iter_bases` is moved off the trait or guarded by a separate method on `dyn ChromRefFetcher` callers. Since `DustFilter` is already generic, this is the easier side; the cohort driver's `Arc<dyn ChromRefFetcher + Send>` hands the fetcher through the merger and DUST without ever calling `iter_bases` on the trait object — `iter_bases` is called only through the generic `F` parameter in `DustFilter`. So the GAT is feasible by keeping `iter_bases` off the dyn-safe call surface (or by splitting the trait into a dyn-safe core + an `IterableChromRefFetcher` subtrait that adds the GAT method).
- **Suggested experiment / fix:**
  ```rust
  // ref_fetcher.rs (sketch — keep the dyn-safe surface)
  pub trait ChromRefFetcher {
      fn length(&self) -> u32;
      fn fetch(&self, start_1based: u32, length: u32) -> Result<Vec<u8>, ChromRefFetchError>;
      // iter_bases stays here only on the dyn-safe surface returning Box.
      fn iter_bases<'a>(&'a self) -> Result<Box<dyn Iterator<Item = Result<u8, ChromRefFetchError>> + 'a>, ChromRefFetchError>;
  }
  // New subtrait for the generic-caller side (DUST):
  pub trait ChromRefFetcherTyped: ChromRefFetcher {
      type BaseIter<'a>: Iterator<Item = Result<u8, ChromRefFetchError>> where Self: 'a;
      fn iter_bases_typed(&self) -> Result<Self::BaseIter<'_>, ChromRefFetchError>;
  }
  ```
  Then `DustFilter<I, F: ChromRefFetcherTyped>` and `ensure_mask_loaded` use `iter_bases_typed()`. The vtable evaporates from the per-byte path; the dyn-safe `iter_bases` survives for any future caller that doesn't know the concrete type.

#### H3: [src/per_sample_pileup/ref_fetcher.rs:368,882](../../../../src/per_sample_pileup/ref_fetcher.rs#L368) — `StreamingChromRefFetcher::fetch` returns `Vec<u8>`; PerGroupMerger consumes it read-only and drops it

- **Confidence:** Medium (no symbol-attributed share, but profile shows allocator symbols at ~14.5 % cpu_atom / ~23.2 % cpu_core combined; the sibling `ManualEvictChromRefFetcher::fetch` already returns `&[u8]`)
- **Hot-path evidence:** verbatim from [profile_summary.txt](../../../../tmp/perf_review_2026-05-23_ref_fetcher/profile_summary.txt):
  ```
  3.87%/4.57%  _int_free_chunk
  3.01%/5.40%  malloc
  2.95%/6.60%  _int_malloc
  2.43%/3.02%  cfree
  1.53%/1.81%  PerGroupMerger::next
  ```
  `process_group` ([per_group_merger.rs:704](../../../../src/var_calling/per_group_merger.rs#L704)) does `let ref_seq = ref_fetcher.fetch(start, span)?;` then passes `&ref_seq` to `unify_alleles` ([:744](../../../../src/var_calling/per_group_merger.rs#L744)), `project_per_position_alleles` ([:860](../../../../src/var_calling/per_group_merger.rs#L860)), and `admit_compound_candidates` ([:861](../../../../src/var_calling/per_group_merger.rs#L861)) — never needs ownership. The fetcher allocates `state.buf[start_idx..end_idx].to_vec()` ([ref_fetcher.rs:882](../../../../src/per_sample_pileup/ref_fetcher.rs#L882)) per call.
- **Pattern matched:** "`Vec<T>` parameter where `&[T]` works is a finding."
- **Mechanism:** One heap alloc + one heap free per variant group, on the cohort PerGroupMerger path. With `--var-group-max-span` default 10 KB and 13 chroms × N=10 × thousands of groups, this is the largest allocator-traffic site directly attributable to the fetcher. Removing it means returning a borrow into `state.buf`, which is feasible **only after H1** because borrowing across a `MutexGuard` is not possible (the slice would outlive the guard). H3 lands as a free side-effect once H1's `RefCell` is in place: `state.buf` becomes directly borrowable from `&self`-bound `RefCell::borrow()`.
- **Measurement plan:** Per-fetcher microbench `streaming_fetch_per_group_window` from §3.2, paired with the cohort-e2e bench. Threshold: ≥ 1 % wall improvement on the cohort-e2e bench attributable to reduced allocator traffic. Validate via DHAT: `cargo run --example dhat_baq --features dhat-heap` style harness around `process_group`.
- **Complexity cost:** Trait signature change: `fn fetch(&self, …) -> Result<&[u8], …>` ties the result to a borrow on `&self`. With `RefCell` (H1), the borrow is the `RefCell::Ref` lifetime — the caller can keep it alive across the use sites because `process_group` is straight-line. The legacy `RefSeqFetcher::fetch` (still used by `WalkerLegacyAdapter` and BAQ via Stage 1) is unaffected — different trait.
- **Suggested experiment / fix:** Sketch in the allocations sub-report. Land in the same PR as H1 because the lifetime threading is the same refactor.

#### H4: [benches/](../../../../benches/) — Cohort-side benches don't compile; **the profile-named hot path is unbenchable today**

- **Confidence:** High
- **Hot-path evidence:** verbatim from `./scripts/dev.sh cargo build --benches --release`:
  ```
  error[E0432]: unresolved import `pop_var_caller::per_sample_pileup::ref_fetcher::SyncRefFetcher`
  error[E0560]: struct `…cohort_driver::CohortPipelineParams` has no field named `fetcher`
  error[E0560]: struct `…cohort_driver::CohortPipelineParams` has no field named `chromosomes`
  error[E0061]: this function takes 7 arguments but 5 arguments were supplied
  error: could not compile `pop_var_caller` (bench "cohort_e2e_perf") due to 9 previous errors
  error: could not compile `pop_var_caller` (bench "var_calling_perf") due to 1 previous error
  ```
  Same shape in `benches/baq_perf.rs` (PreparedRead field drift) and `benches/pileup_walker_scaling.rs`. Only `psp_reader_perf`, `psp_writer_perf`, `freebayes_bookkeeping` still build.
- **Pattern matched:** "No optimization without a benchmark" — combined with the `MEMORY.md` "Bench CPU governor caveat" note that wall-clock numbers drift ~2× across sessions. Without a working bench, every H1/H2/H3 measurement above relies on production-binary wall time alone — too noisy to detect 5 % effects.
- **Mechanism:** Benches reference the pre-Step-2 API: `SyncRefFetcher` no longer exists; `CohortPipelineParams.fetcher`/`.chromosomes` are gone (the current shape is in [var_calling_from_bam.rs:402-412](../../../../src/pop_var_caller/var_calling_from_bam.rs#L402-L412)); `drive_cohort_pipeline` takes 7 args now (`chrom_id, merger, params, fetcher, &vcf_out, metadata, writer_cfg` — see [cohort_driver.rs:200](../../../../src/pop_var_caller/cohort_driver.rs#L200)). The repair is mechanical — same recipe as the Step-2 migration that broke them.
- **Measurement plan:** Repair, then `cargo bench --bench cohort_e2e_perf -- e2e_core/cohort_scaling_threads --save-baseline ref-fetcher-baseline`. Sanity check: the bench wall time should land within ±15 % of the production-driver number on the same fixture at `--threads 4`.
- **Complexity cost:** Mechanical bench-fixture port. No new deps, no `unsafe`, no API surface change.
- **Suggested experiment / fix:** Fix `cohort_e2e_perf`, `var_calling_perf`, `baq_perf`, `pileup_walker_scaling` as one PR. Confirm via `cargo build --benches --release` in CI.

### Likely

#### L1: [src/per_sample_pileup/ref_fetcher.rs:447](../../../../src/per_sample_pileup/ref_fetcher.rs#L447) — `refill` reads raw 64 KiB chunks from `File` without `BufReader`

- **Confidence:** Medium
- **Hot-path evidence:** `0.94 % / 0.72 %  ref_fetcher::refill` in the cohort profile. Refill is ~once per 1 MB of bases; each refill emits 16 `read(2)` syscalls + 1 `seek(2)`.
- **Pattern matched:** "`File::read` without buffering are a per-call syscall."
- **Mechanism:** `Source::read_chunk` ([ref_fetcher.rs:166](../../../../src/per_sample_pileup/ref_fetcher.rs#L166)) goes straight through `File::read`. A `BufReader::with_capacity(64 KiB, file)` keeps the same on-disk byte volume but lets glibc/page-cache prefetch fill larger user-space buffers. Win is mostly the syscall trampoline cost; mirror the [cohort_driver.rs:488](../../../../src/pop_var_caller/cohort_driver.rs#L488) PSP-reader pattern (`BufReader::with_capacity(DEFAULT_BUFFERED_IO_CAPACITY, file)`).
- **Measurement plan:** Re-run `samply` on `profile_cohort_e2e`; threshold — `refill` symbol drops by at least the syscall-trampoline cost (~µs per call × refills-per-run); total wall not worse.
- **Complexity cost:** `Source::File(File)` → `Source::File(BufReader<File>)` + `seek_to`/`read_chunk` delegation. `BufReader` implements `Seek`; the buffer is invalidated on seek, which is what refill always does anyway.
- **Suggested experiment / fix:** Diff in the io_and_syscalls sub-report. Same change for `ManualEvictChromRefFetcher`'s `file` field (L8 below — parity).

#### L2: [src/per_sample_pileup/ref_fetcher.rs:461-469](../../../../src/per_sample_pileup/ref_fetcher.rs#L461-L469) + [:1390-1398](../../../../src/per_sample_pileup/ref_fetcher.rs#L1390-L1398) — Per-byte `to_ascii_uppercase()` + newline filter blocks autovec in `refill` / `read_uppercased_bases`

- **Confidence:** Medium
- **Hot-path evidence:** refill at 0.94 % / 0.72 %; cold-but-rises-if-H2-Option-2-lands.
- **Pattern matched:** "Autovectorization needs the compiler's confidence" — three data-dependent branches per byte (cap, newline, uppercase) + `Vec::push` defeats the LLVM autovec pass.
- **Mechanism:** Restructure as `extend(read_buf[..n].iter().filter(...).take(...).copied())` followed by `buf[before..].make_ascii_uppercase()` — the in-place `make_ascii_uppercase` is autovec-friendly (lowers to `vpor`/`vpand`-class instructions under `target-cpu=x86-64-v3`).
- **Measurement plan:** Confirm with `cargo asm pop_var_caller::per_sample_pileup::ref_fetcher::refill --rust` that the current loop is *not* vectorised; if it is, downgrade to Note. Threshold: refill drops below 0.3 % in the post-H1 profile.
- **Complexity cost:** Small. No new deps, no `unsafe`.
- **Suggested experiment / fix:** Diff in the hot_loops sub-report.

#### L3: [src/var_calling/dust_filter.rs:728-730](../../../../src/var_calling/dust_filter.rs#L728-L730) — `io::Error::other(format!("{e}"))` allocates a `String` per error wrap

- **Confidence:** Low (success-path bytes don't allocate; only the `Err` arm does)
- **Hot-path evidence:** pattern-match only.
- **Pattern matched:** `format!` whose output is immediately fed into `Error::other`.
- **Mechanism:** `ChromRefFetchError: std::error::Error`, so `io::Error::other(e)` works directly — no `format!`, preserves `source()`.
- **Measurement plan:** None — correctness-of-allocation-shape fix.
- **Complexity cost:** One-line change.
- **Suggested experiment / fix:** `bases.map(|res| res.map_err(io::Error::other))`. (Subsumed by H2's deeper fix that drops the `.map` adapter entirely.)

#### L4: [src/per_sample_pileup/ref_fetcher.rs:1353-1360](../../../../src/per_sample_pileup/ref_fetcher.rs#L1353-L1360) — `ManualEvictChromRefFetcher::prepend_backward` allocates a scratch `Vec` then `splice`s it

- **Confidence:** Low (BAQ profile shows the fetcher cold at 0.02 %; flagged for completeness)
- **Hot-path evidence:** pattern-match only.
- **Pattern matched:** "`Vec::with_capacity(n)` for a one-shot scratch + `Vec::splice(..0, …)`."
- **Mechanism:** Reserve once on `self.buf`, `copy_within` existing content right by `extra_bases`, then read straight into the prefix gap — eliminates the scratch.
- **Measurement plan:** Site is cold; defer until a backward-jump-heavy consumer exists.
- **Complexity cost:** Minimal. `read_uppercased_bases` currently appends; teaching it to write at an in-place offset is the only awkward bit.
- **Suggested experiment / fix:** Defer.

#### L5: [src/per_sample_pileup/ref_fetcher.rs:813](../../../../src/per_sample_pileup/ref_fetcher.rs#L813) — `Drop for ChromRefBaseIter` re-locks the mutex

- **Confidence:** Medium
- **Hot-path evidence:** pattern-match only; once per `iter_bases` call (~13× per cohort run × per-chrom-worker).
- **Pattern matched:** Lock pair at the cold edge of the same H1 anti-pattern.
- **Mechanism:** Folds into H1 for free — under `RefCell`, the Drop body is one borrow_mut + clear, no atomics.
- **Measurement plan:** Folded into H1.
- **Complexity cost:** Zero on top of H1.
- **Suggested experiment / fix:** Falls out of H1's `RefCell` switch.

#### L6: [Cargo.toml:35-39](../../../../Cargo.toml#L35-L39) — Add `[profile.release-with-debug]` for reproducible profiles

- **Confidence:** High
- **Hot-path evidence:** `[profile.release] debug = "line-tables-only"`; this report's call-graph required `debug = true`. Without a documented profile, the next reviewer cannot reproduce the inlined-frame chain from a clean checkout.
- **Pattern matched:** "Profiles must come from release builds with debug info."
- **Mechanism:** Reproducibility insurance; no wall-time threshold.
- **Measurement plan:** N/A.
- **Complexity cost:** One new profile block + one doc line.
- **Suggested experiment / fix:**
  ```toml
  [profile.release-with-debug]
  inherits = "release"
  debug = true
  ```

#### L7: [src/main.rs](../../../../src/main.rs) — Wire `mimalloc::MiMalloc` as `#[global_allocator]` behind `alloc-mimalloc`

- **Confidence:** Medium
- **Hot-path evidence:** allocator symbols sum to ~14.5 % cpu_atom / ~23.2 % cpu_core in the cohort profile. Feature is wired for benches; production binary doesn't pick it up.
- **Pattern matched:** "Allocator choice is a one-line change." Heavy allocator traffic + no A/B = Likely.
- **Mechanism:** Per-thread caches in mimalloc remove the inter-arena lock + fragmentation overhead that `_int_malloc` / `_int_free` show.
- **Measurement plan:** Step 5 of §3.
- **Complexity cost:** One `#[cfg_attr(feature = "alloc-mimalloc", global_allocator)]` block in `main.rs`. No new deps (mimalloc already in manifest).
- **Suggested experiment / fix:**
  ```rust
  // src/main.rs (top)
  #[cfg(feature = "alloc-mimalloc")]
  #[global_allocator]
  static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
  ```

#### L8: [benches/](../../../../benches/) — Add per-fetcher microbench `benches/ref_fetcher_perf.rs`

- **Confidence:** High
- **Hot-path evidence:** no bench exercises `ChromRefBaseIter::next`, `StreamingChromRefFetcher::fetch`, or `ManualEvictChromRefFetcher::fetch` directly. The cohort_e2e bench (once H4-repaired) measures the full pipeline; localising H1/H2/H3 to the fetcher needs an isolated bench.
- **Pattern matched:** "No optimization without a benchmark" — fetch + iter_bases are the right primitives to bench in isolation.
- **Mechanism:** Four `criterion` functions: `streaming_iter_bases_full_contig`, `streaming_fetch_per_group_window`, `manual_evict_fetch_then_evict`, `refill_warm_cache`.
- **Measurement plan:** See methodology sub-report for the bench skeleton.
- **Complexity cost:** ~300 LOC, one new file, no new deps.
- **Suggested experiment / fix:** Sketch in the methodology sub-report.

#### L9: [src/per_sample_pileup/ref_fetcher.rs:148](../../../../src/per_sample_pileup/ref_fetcher.rs#L148) — `Source` enum is sized for the test-only `Memory(Cursor<Vec<u8>>)` variant

- **Confidence:** Medium
- **Hot-path evidence:** pattern-match only; refill is 0.94 %.
- **Pattern matched:** Enum sized for the larger variant + discriminant + alignment (~40 B) on every `StreamState`, just to support a `#[cfg(test)]` code path.
- **Mechanism:** In production each `StreamState` pays ~32 bytes of layout for the test cursor. Cold today; promotes if the H1 fix re-densifies the per-byte working set into the cache line.
- **Measurement plan:** Defer until H1 lands; re-measure refill.
- **Complexity cost:** Either (a) gate `Memory` behind `#[cfg(test)]` (with a separate test fetcher constructor), or (b) replace the enum with `Box<dyn Read + Seek>` (vtable lookup on refill — likely worse).
- **Suggested experiment / fix:** Defer.

#### L10: [src/per_sample_pileup/ref_fetcher.rs:1371](../../../../src/per_sample_pileup/ref_fetcher.rs#L1371) — `read_uppercased_bases` reads raw 64 KiB chunks from `&mut File`

- **Confidence:** Low
- **Hot-path evidence:** 0.07 % cpu_core in Stage 1 BAQ profile (cold).
- **Pattern matched:** Same as L1 — `File::read` without `BufReader`.
- **Mechanism:** Same fix as L1 — `&mut BufReader<File>` parameter.
- **Measurement plan:** Combined with L1.
- **Complexity cost:** Same as L1.
- **Suggested experiment / fix:** Same diff shape as L1; pair the two fixes in one commit.

### Speculative

#### S1: [Cargo.toml](../../../../Cargo.toml) — PGO on the cohort var-calling binary

- **Confidence:** Low
- **Hot-path evidence:** pattern-match only. Methodology rule: PGO is on the table when LTO + codegen-units=1 are already set; both are.
- **Pattern matched:** Stable critical paths + repeatable training corpus = PGO candidate.
- **Mechanism:** Biases inlining + basic-block layout toward hot branches in a training run.
- **Measurement plan:** Defer until allocator A/B (L7) has landed; methodology says "one change per measurement."
- **Complexity cost:** Multi-stage build, representative training corpus, re-train when workload shape changes.
- **Suggested experiment / fix:** Defer.

#### S2: [src/per_sample_pileup/ref_fetcher.rs:670-733](../../../../src/per_sample_pileup/ref_fetcher.rs#L670-L733) — `.fai` re-parsed per fetcher construction (13× per cohort run)

- **Confidence:** Low
- **Hot-path evidence:** below noise floor (~13 ms / 58,850 ms = 0.02 %).
- **Pattern matched:** Repeated work that could be shared as `Arc<fai::Index>`.
- **Mechanism:** Shared parsed index passed by reference from `run_var_calling` down to `process_one_chromosome`.
- **Measurement plan:** Add a one-off `Instant::now()` around `fai::fs::read`; threshold to act on: > 50 ms total.
- **Complexity cost:** One new constructor variant + one new parameter on `process_one_chromosome`.
- **Suggested experiment / fix:** Defer.

#### S3: [src/per_sample_pileup/ref_fetcher.rs:174](../../../../src/per_sample_pileup/ref_fetcher.rs#L174) — `StreamState` field order under `repr(Rust)`

- **Confidence:** Low
- **Hot-path evidence:** pattern-match only.
- **Pattern matched:** Hot triplet `(buf.ptr, buf.len, buf_start_base)` should be cache-adjacent.
- **Mechanism:** Subsumed by H1's iter-local snapshot. Standalone, no action.
- **Measurement plan:** `dbg!(mem::size_of::<StreamState>(), mem::offset_of!(StreamState, buf), …)` if curiosity strikes.
- **Complexity cost:** None for the diagnostic.
- **Suggested experiment / fix:** Subsumed by H1.

#### S4: [rust-toolchain.toml](../../../../rust-toolchain.toml) — verify channel pin (likely already in place)

- **Confidence:** Low
- **Hot-path evidence:** bench-build output mentions `1.95-x86_64-unknown-linux-gnu`; verifies indirectly.
- **Pattern matched:** Compiler-version-drift hygiene.
- **Mechanism:** N/A if pin already exists.
- **Measurement plan:** `ls rust-toolchain.toml && cat rust-toolchain.toml`.
- **Complexity cost:** One file if missing.
- **Suggested experiment / fix:** Verify-then-no-op.

### Note

- **Good build config already in place** — `lto = "fat"`, `codegen-units = 1`, `panic = "abort"`, `target-cpu=x86-64-v3`, `[profile.bench] debug = true`, `dhat-heap`/`alloc-mimalloc` feature wiring.
- **Good test-fixture pattern** — `benches/cohort_e2e_perf.rs:184` uses `TempDir` for a real on-disk FASTA fixture; mirror this in `benches/ref_fetcher_perf.rs` (L8).
- **Good pre-sizing in `unify_alleles`** ([per_group_merger.rs:893-894](../../../../src/var_calling/per_group_merger.rs#L893-L894)) — `Vec<UnifiedAllele>` and `AHashMap<Vec<u8>, usize>` already pre-sized to the computed upper bound.
- **BAQ-side fetcher is cold on this fixture** — `ManualEvictChromRefFetcher::fetch` 0.02 % / 0.02 %, `read_uppercased_bases` 0.07 %, `refill` (manual-evict path) 0.01 %. Findings against the BAQ-side `ManualEvictChromRefFetcher` cap at Note (or Likely for ergonomic-parity changes like L10).
- **`ContigFai` pre-cached at construction** — refill reads only cached fields, no re-parse per call. Good shape.
- **Per-worker fetcher isolation** — each `StreamingChromRefFetcher` is its own allocation; refcount stays at 1–2 per worker; never crosses worker boundaries. No false-sharing risk; `CachePadded` rule does not apply.

## 6. Out-of-scope observations

- **PSP reader (zstd decode + RecordsIter) at 3.61–6.50 % wall** — outside this review's scope but visible in the same profile. Tracked under the 2026-05-20 review's L3 (`PSP CSR decoder`).
- **Allocator symbols charged to non-fetcher consumers** — `malloc`/`free` family sum ~14–18 % combined, only a fraction of which is the fetcher (`StreamingChromRefFetcher::fetch -> Vec`). The remainder is PSP decode + per-group merger internal buffers + posterior engine. L7 (mimalloc) helps all of them.
- **`md5::compress::compress` at 0.93 % / 1.02 %** — FASTA MD5 verify (now streaming-windowed per commit `839af2e`); the work is the MD5 itself, not the I/O. No action.

## 7. What's already good

- **Per-worker fetcher ownership invariant** is documented at the struct ([ref_fetcher.rs:113-117](../../../../src/per_sample_pileup/ref_fetcher.rs#L113-L117)) — H1's mechanism rests on this comment being load-bearing, and it is.
- **`ManualEvictChromRefFetcher` already returns `&[u8]`** from `fetch` ([ref_fetcher.rs:1188-1240](../../../../src/per_sample_pileup/ref_fetcher.rs#L1188-L1240)) — the in-tree precedent for H3.
- **DUST `ensure_mask_loaded` streams via `iter_bases`** ([dust_filter.rs:703-751](../../../../src/var_calling/dust_filter.rs#L703-L751)) instead of materialising the whole contig — the right call shape; the per-byte vtable cost (H2) is the *only* thing in the way of this design fully landing.

### Author response convention

Address each finding by its identifier (H1, H2, …, L1, L2, …, S1, …) with one of:
- `applied in <commit>`
- `experiment shows no gain — closing`
- `disputed because …`
- `deferred to <issue>`
- `won't fix because …`

The "experiment shows no gain" path is expected and welcome — that's what the measurement plan is for.
