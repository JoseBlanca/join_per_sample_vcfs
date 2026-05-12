# Performance Review: pileup (round 2)
**Date:** 2026-05-12
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** `src/per_sample_caller/pileup/` — the sequential pileup walker (Stage 1), re-reviewed after the 2026-05-10 wins shipped
**Verdict:** Apply the listed wins
**Hot-path evidence:** fresh `perf record --call-graph dwarf -F 999 --profile-time 20` runs against `pileup_walker_multi_op/5000` and `pileup_walker_read_length/150` on the host (`perf_event_paranoid = 2`, x86-64-v3); fresh DHAT pass against `multi_op/L=5000`. Profile inputs are saved verbatim in [tmp/perf_review_2026-05-12_pileup/profile_summary.md](../../tmp/perf_review_2026-05-12_pileup/profile_summary.md).

---

## 1. Scope and constraints

- **What was reviewed:** the sequential walker module — public surface ([mod.rs](../../src/per_sample_caller/pileup/mod.rs)) plus all internal submodules ([walker.rs](../../src/per_sample_caller/pileup/walker.rs), [active_read_set.rs](../../src/per_sample_caller/pileup/active_read_set.rs), [cigar_cursor.rs](../../src/per_sample_caller/pileup/cigar_cursor.rs), [decompose.rs](../../src/per_sample_caller/pileup/decompose.rs), [open_record.rs](../../src/per_sample_caller/pileup/open_record.rs), [slot_allocator.rs](../../src/per_sample_caller/pileup/slot_allocator.rs), [errors.rs](../../src/per_sample_caller/pileup/errors.rs)). Build configuration ([Cargo.toml](../../Cargo.toml), [.cargo/config.toml](../../.cargo/config.toml)) and the bench harness ([benches/pileup_walker_scaling.rs](../../benches/pileup_walker_scaling.rs), [examples/dhat_pileup.rs](../../examples/dhat_pileup.rs)) are included for methodology-category review.
- **Reviewed against:** `main` at `dbca4a1` (no in-flight changes; only renames and BAQ-stage additions have landed in this module since the 2026-05-10 round).
- **Performance intent and targets:** WGS-scale inputs — multi-GB BAM/CRAM, ~30× coverage, chromosomes up to ~250 Mb. The walker is the per-sample caller's hot path; every base position passes through it once. Production calls [`RefSeqFetcher`](../../src/per_sample_caller/pileup/mod.rs#L450-L462) backed by `noodles_fasta::Repository`; the bench substitutes a constant fetcher. Target hardware: Linux x86_64, microarchitecture floor `x86-64-v3` (AVX2 + FMA + BMI2). Round-1 (2026-05-10) delivered −54 % wall-clock and −93 % allocations; the round-2 question is what remains.
- **Hot-path evidence (fresh, this session):**

  | Source | Fixture | Captured to |
  |---|---|---|
  | `perf record --call-graph dwarf -F 999 --profile-time 20` | `pileup_walker_multi_op/5000` (31 K samples) | [tmp/perf_review_2026-05-12_pileup/perf_multi_op_5000.children.txt](../../tmp/perf_review_2026-05-12_pileup/perf_multi_op_5000.children.txt) |
  | (same) | `pileup_walker_read_length/{150, 1500}` (66 K samples) | [perf_read_length_150.children.txt](../../tmp/perf_review_2026-05-12_pileup/perf_read_length_150.children.txt) |
  | DHAT via [examples/dhat_pileup.rs](../../examples/dhat_pileup.rs) | `multi_op/L=5000` (49 850 records) | [dhat-heap.json](../../tmp/perf_review_2026-05-12_pileup/dhat-heap.json), [dhat_top.txt](../../tmp/perf_review_2026-05-12_pileup/dhat_top.txt) |

  Highlights (children-aggregated, walker thread, `multi_op/5000`):

  | % | Symbol | Bucket |
  |---:|---|---|
  | 93.36 | `run_walker` | total walker |
  | 84.86 | `WalkerState::process_position` | dominates |
  | 59.52 | `open_record::process_position` | heaviest single function |
  | 15.50 | `CigarCursor::events_overlapping` (binary) | cursor |
  | 9.37 | `CigarCursor::events_at` (binary) | cursor |
  | 12.94 | `slice::binary_search_by` | inside cursor |
  | 10.64 | `slice::partition_point` | inside cursor |
  | 7.73 | `close_aged_records` | mostly send + drop |
  | 6.64 | `OpenPileupRecordTable::open_new` | fresh record |
  | 5.95 | `insert_sorted_unique` | `chain_slots` insert+grow |
  | 4.49 | `SyncSender::send` (3.95 self) | mpsc |
  | 4.44 | `AHashMap::with_capacity` | `folded_reads` alloc |
  | 4.27 | `find_overlapping` | `BTreeMap::range` |
  | 4.14 | `ln_bq_for_read` (2.94 in `phred_to_ln_perr`) | fold body |
  | 3.95 | `resolve_mate_overlap_at_pos` no-overlap fast-path | O(n²) `any()` |
  | 3.55 | `AHashMap::insert` | `folded_reads` insert |
  | 3.04 | `apply_events_to_ref_into` | fold body |
  | 2.47 | `ActiveReads::get_by_read_id` | AHashMap probe |
  | 1.91 | `column_depth_cap` | redundant `any()` scan |
  | 1.91 | `find_allele_index` | linear seq.eq |

  DHAT top remaining sites (49 850 records emitted, 504 172 allocs, 158.65 MB):

  | # | Blocks | Total B | B/blk | Origin |
  |---:|---:|---:|---:|---|
  | 1 | 99 700 | 8.77 MB | 88 | `drain_aged` `out: Vec<OpenPileupRecord>` + `closing_keys`, 2×/step |
  | 2 | 89 295 | 8.24 MB | 92 | `OpenAllele::new` `chain_slots` Vec (cap 32, grows) |
  | 3 | 50 846 | 3.88 MB | 76 | `OpenPileupRecord::new` `alleles` Vec |
  | 5 | **49 850** | **131.6 MB** | **2 640** | **`OpenPileupRecord::new` `folded_reads: AHashMap::with_capacity(32)` — 83 % of total bytes** |
  | 6 | 49 850 | 0.80 MB | 16 | `process_position` `affected: Vec<u32>` |
  | 8+9 | 29 700×2 | 0.12 MB | 2 | `events_*_binary` Insertion seq via `to_vec()` |

- **Deliberately out of scope:** [tests.rs](../../src/per_sample_caller/pileup/tests.rs) (correctness coverage, not perf); the production `RefSeqFetcher` impl (lives outside `pileup/`); [src/per_sample_caller/baq/](../../src/per_sample_caller/baq/) (own review on the same date).
- **Categories dispatched:**
  - `methodology` — always
  - `allocations` — DHAT shows one 131 MB site + ~150 K small-block sites still in the hot path
  - `data_layout` — the AHashMap dominator is a data-structure choice, not just an alloc; ~25 % cursor binary-search machinery on tap
  - `concurrency` — `mpsc::sync_channel` + receiver-thread drop cost visible in the profile
  - `hot_loops` — `phred_to_ln_perr` 2.94 %, `column_depth_cap` redundant scan, cursor inner loops, `BINARY_SEARCH_OP_THRESHOLD` set pre-LTO
  - **Skipped:** `io_and_syscalls` — the walker does no I/O; `RefSeqFetcher` is a trait boundary.

  Per-category audit trail at [tmp/perf_review_2026-05-12_pileup/](../../tmp/perf_review_2026-05-12_pileup/) ([methodology.md](../../tmp/perf_review_2026-05-12_pileup/methodology.md), [allocations.md](../../tmp/perf_review_2026-05-12_pileup/allocations.md), [data_layout.md](../../tmp/perf_review_2026-05-12_pileup/data_layout.md), [concurrency.md](../../tmp/perf_review_2026-05-12_pileup/concurrency.md), [hot_loops.md](../../tmp/perf_review_2026-05-12_pileup/hot_loops.md)).

## 2. Verdict

**Apply the listed wins.** Three Hot-path findings (H1, H2, H3) are well-evidenced, each names a specific profile/DHAT line, each has a plausible mechanism and contained complexity. Apply in the order listed in §5 so the build configuration and the lowest-risk wins land first; one change per measurement; gate each on its named threshold and revert on no signal (round-1 reverted L14 on exactly that gate).

Three companion methodology items (L1 allocator A/B, L3 save `final-2026-05-10` baseline, L4 + L5 bench hygiene) are not optimisations themselves but unblock or de-risk every measurement that follows.

## 3. Measurement plan

In order of leverage:

1. **Save the `final-2026-05-10` baseline first.** Round-1's results doc lists this as the open follow-up; without it every round-2 comparison degrades to the silently-overwritten `base`. One command (see L3 in §5). Five minutes.
2. **Re-bench after each commit** with `--save-baseline <commit-slug>` / `--baseline final-2026-05-10`. One change per commit (round-1's process; keep it).
3. **DHAT after each commit** via `cargo run --release --example dhat_pileup --features dhat-heap`; inspect with [tmp/perf_review_2026-05-12_pileup/_dhat_top.py](../../tmp/perf_review_2026-05-12_pileup/_dhat_top.py) to find the next-largest alloc site. H1 should remove the 131 MB AHashMap line entirely; H2 should remove the 99 700 + 49 850 block sites; the count moves predictably and the report attributes the win or its absence.
4. **Re-profile (`perf record`) at the end of each major commit** to confirm the children-aggregated row the finding targeted has dropped (`AHashMap::with_capacity` to ~0; `phred_to_ln_perr` gone; etc.). Mechanical confirmation against the prediction in each finding.
5. **`cache-misses` / `L1-dcache-load-misses` not yet captured.** Host `perf_event_paranoid = 2` allows cycles sampling but not hardware counters that need `paranoid ≤ 1`. The data-layout cache-miss arguments below carry a "would need paranoid ≤ 1" caveat. If a finding's verdict hinges on cache-miss data, ask the user to lower `perf_event_paranoid` for that one measurement.
6. **Build the paired-end fixture before re-evaluating L17/L18.** The current bench is `MateRole::Solo` everywhere; the `Arc<str>` qname clone sites can't be measured at all until a paired-end fixture exists. Round-1 deferred this and so does round-2.

## 4. Build / toolchain configuration

Round-1's build floor (`lto = "fat"`, `codegen-units = 1`, `bench` inherits release + debug) and the post-round-1 `target-cpu = x86-64-v3` are in place. No further build-config Hot-path candidates this round; the remaining items are Likely (allocator A/B, bench harness) or Speculative (panic=abort, rust-toolchain.toml, streaming-WGS-span fixture). See §5 L1, L2, L3, L4, L5 and S1–S3.

## 5. Code-level findings

### Hot-path

#### H1: [open_record.rs:79, :103](../../src/per_sample_caller/pileup/open_record.rs#L79-L103) — **[Hot-path]** Per-record `folded_reads: AHashMap<u32, FoldedReadState>` is the 131 MB / 159 MB DHAT dominator and ~10 % of walker CPU — replace it

- **Confidence:** High
- **Hot-path evidence:** DHAT site #5: 49 850 blocks × 2 640 B = **131.6 MB / 158.65 MB total bytes (83 %)**. Perf children: `AHashMap::with_capacity` 4.44 %, `AHashMap::insert` 3.55 %, `AHashMap::get` 2.10 % — ~10.1 % of walker CPU. (Both numbers quoted verbatim from [profile_summary.md](../../tmp/perf_review_2026-05-12_pileup/profile_summary.md).) Three sub-agents independently flagged this site as the top remaining issue.
- **Pattern matched:** *Allocations belong outside hot loops* + *Inline small keys in maps* + *Pointer chases fragment the cache*. The map sizes to 32 entries with `~30 entries` of typical occupancy; the hash-table's per-record 2.6 KB backing buffer dwarfs the data it actually holds, and lookups miss L1 because buckets live on a separate heap allocation from the record.
- **Mechanism:** The fold pattern per record is *one allocation at open, ~30 inserts, ~30 lookups, one drop at close*. At ≤ 30 keys, a hash map's asymptotic wins are eclipsed by (a) the per-record 2.6 KB capacity-32 allocation and (b) the per-lookup hash + bucket fetch into a separate cache line. Two replacement shapes are worth comparing in measurement:
  - **(A) Sorted `Vec<(u32, FoldedReadState)>` keyed by `read_id`** with `binary_search_by_key` (~5 cmps over a contiguous 1 KB array). Insert / remove is `Vec::insert` / `Vec::remove` of ≤ 30 entries — one `memmove` of ~500 B per modify (one cache line of bandwidth, same order as the hash probe + fetch).
  - **(B) `SmallVec<[(u32, FoldedReadState); 16]>`** with linear `position` lookup. 16 inline rows handles half the typical occupancy on the stack; the spill case stays on the heap. Linear scan over ≤ 30 `u32` keys beats hashbrown at this size on most x86_64.
- **Measurement plan:**
  1. Implement (A) as the first experiment (cleanest type-system shape; no new dep — `Vec::binary_search_by_key` is std). Helper pair `insert_or_replace_folded` / `remove_folded` as drafted in the data_layout sub-agent's diff.
  2. `cargo bench --bench pileup_walker_scaling -- --baseline final-2026-05-10`. **Threshold to merge:** ≥ 4 % mean improvement across the 8 fixtures, with no individual fixture regressing > 2 % at p < 0.05. DHAT site #5 should drop to a different shape (one Vec per record at ~1 KB) or disappear entirely if (C) shared-pool below is taken instead.
  3. Perf children: `AHashMap::with_capacity` (4.44 %), `AHashMap::insert` (3.55 %), `AHashMap::get` (2.10 %) collapse to a single `slice::binary_search_by` row ≤ 3 %.
  4. If (A) wins, try (B) as a follow-up; if (A) regresses, try (C) — hoist a single `folded_reads` onto `OpenPileupRecordTable` keyed by `(record_key, read_id)` and clear per-record at close.
- **Complexity cost:** ~30 lines for (A). No `unsafe`. No public API change (`FoldedReadState` is `pub(super)`). The "no-double-insert" invariant in [process_position](../../src/per_sample_caller/pileup/open_record.rs#L744-L766) (every insert is preceded by a remove) is what makes the binary-search insertion safe.
- **Suggested fix:** see the data_layout sub-agent's diff at [tmp/perf_review_2026-05-12_pileup/data_layout.md](../../tmp/perf_review_2026-05-12_pileup/data_layout.md#L104-L142).

#### H2: [open_record.rs:218-232](../../src/per_sample_caller/pileup/open_record.rs#L218-L232) — **[Hot-path]** Hoist `drain_aged`'s two per-call `Vec`s onto `OpenPileupRecordTable` (`drain_aged_into` shape)

- **Confidence:** High
- **Hot-path evidence:** DHAT sites #1 + #4: **99 700 + 49 850 = ~150 K alloc events / run**, third- and fourth-largest by block count after the AHashMap cluster. 99 700 = 2 vecs × ~50 K calls — matches the "drain_aged 2×/step" assumption verbatim. Walker-side `__GI___libc_malloc` is 7.20 % children; collapsing two of its most-frequent small-allocation sources moves it visibly.
- **Pattern matched:** *Allocations belong outside hot loops*. Same shape as round-1 L6 (`contributors_buf` hoist) and the round-1 `allele_seq_buf` hoist that landed.
- **Mechanism:** [close_aged_records](../../src/per_sample_caller/pileup/walker.rs#L363-L379) calls [drain_aged](../../src/per_sample_caller/pileup/open_record.rs#L218-L232) once per walker step. Each call freshly allocates a `Vec<u32>` (`closing_keys`) and a `Vec<OpenPileupRecord>` (`out`); both typically hold ≤ 1 record at steady state. Hoisting both onto `OpenPileupRecordTable` and `clear()`ing between calls turns 2 × 50 K = 100 K small allocations into 2 (per-instance). The caller is currently `aged.into_iter().map(finalise).collect()`; switch to a `drain_aged_into(walker_pos, out)` shape on a caller-supplied `Vec<OpenPileupRecord>` (mirrors [apply_events_to_ref_into](../../src/per_sample_caller/pileup/open_record.rs#L427)) so the field stays scratch.

  *This is not round-1 L12* — L12 was an empty-case short-circuit (rejected because records age out ~1/step). H2 is a buffer hoist whose win does not depend on the empty-case rate; the alloc happens every call regardless.
- **Measurement plan:**
  - DHAT sites #1 + #4 drop to zero blocks (or to per-instance scratch).
  - Total alloc events drop by ~150 K.
  - Bench: same scaling sweep. **Threshold:** ≥ 2 % wall-clock on either fixture; an alloc-count drop without a wall-clock movement is fine evidence of the bytes win.
- **Complexity cost:** One new field on `OpenPileupRecordTable` (`closing_keys_buf: Vec<u32>`). One field on `WalkerState` (`drained_buf: Vec<OpenPileupRecord>`) alongside `contributors_buf`. `reset()` clears them. Caller changes from `aged.into_iter()...collect()` to draining the scratch — one extra `drain(..)` call. No `unsafe`.
- **Suggested fix:** [tmp/perf_review_2026-05-12_pileup/allocations.md](../../tmp/perf_review_2026-05-12_pileup/allocations.md#L190-L221) has the diff.

#### H3: [open_record.rs:878](../../src/per_sample_caller/pileup/open_record.rs#L878) — **[Hot-path]** Replace `phred_to_ln_perr(q)` with a 256-entry lookup table

- **Confidence:** High
- **Hot-path evidence:** Perf children: `phred_to_ln_perr` 2.94 % self, parent `ln_bq_for_read` 4.14 % children. Called once per (affected record × contributor) — ~1.5 M calls / run on `multi_op/5000` (matches the call-frequency note in [profile_summary.md](../../tmp/perf_review_2026-05-12_pileup/profile_summary.md)).
- **Pattern matched:** *Avoid recomputing in tight loops. Hoist invariants.* The function runs two FP multiplies + an int→FP convert + a branch on every call; the input is constrained to `0..=93` by the [PreparedRead spec](../../src/per_sample_caller/pileup/mod.rs#L203) ("Phred 0–93"), so the function is a 94-entry total map over `u8`. Sizing the table at 256 (matching `u8`'s full domain) elides the bounds check entirely.
- **Mechanism:** 1.5 M f64 multiplies and divides disappear. The table is 2 KB and sits in L1.
- **Measurement plan:**
  - `perf record` re-run: `phred_to_ln_perr` should disappear from the children-aggregated table.
  - `cargo bench`: **threshold ≥ 1.5 %** on `multi_op/5000` (function's self-time is 2.94 %, so a ~2 % wall-clock win is realistic). The fix is small enough that even no-signal is acceptable evidence to ship as a hygiene cleanup.
- **Complexity cost:** One module-level `static [f64; 256]` (computable via `const`-evaluator on stable since rustc 1.83 — the project's MSRV is implicit but rustc 1.93 is in the toolchain). No new dep, no `unsafe`, no API change. Pin the `q == 0` entry to `+0.0` (the existing function returns `0.0`; a `-0.0` would be a no-op semantic difference but might confuse a reader).
- **Suggested fix:** see [tmp/perf_review_2026-05-12_pileup/hot_loops.md](../../tmp/perf_review_2026-05-12_pileup/hot_loops.md#L27-L47).

### Likely

#### L1: [Cargo.toml](../../Cargo.toml) — **[Likely]** [build] No allocator A/B run despite a 131 MB single-site bytes footprint

- **Confidence:** Medium
- **Hot-path evidence:** the DHAT site #5 (131 MB / 159 MB AHashMap) plus walker-side `_int_malloc` 3.06 %, `malloc_consolidate` 3.40 %, `__GI___libc_malloc` 7.20 % children. Round-1 filed this as L2 and deferred it; round-2 evidence (the bytes-side dominance) is the new motivation.
- **Pattern matched:** *Allocator choice is a one-line change.*
- **Mechanism:** glibc's `malloc` with per-thread arenas + global free-list maintenance (`malloc_consolidate` shows up in the walker thread's top frames). `mimalloc` / `jemalloc` swap that for per-thread slab caches with size-class fast paths that suit the walker's "many ~2.6 KB AHashMap allocate-and-drop per record" shape.
- **Measurement plan:** add `mimalloc` and `jemallocator` behind feature flags (one commit each), bench against `final-2026-05-10`. **Threshold:** any clean criterion "improved" on `multi_op/5000` AND no regression > 3 % on the others. Run DHAT alongside to compare allocator-internal residency.
  Note: if H1 lands first, the 131 MB site disappears and the allocator A/B win on this fixture may evaporate. That is the expected outcome; do allocator A/B *after* H1 to measure what remains.
- **Complexity cost:** one optional dependency, three lines of `#[cfg(feature = "…")]`. No `unsafe`, no API change.

#### L2: [pileup_walker_scaling.rs:182-263](../../benches/pileup_walker_scaling.rs#L182-L263) — **[Likely]** [bench] `iter_batched(LargeInput)` still leaks fixture + thread-setup cost into the timed region

- **Confidence:** Medium
- **Hot-path evidence:** pattern-match only. Round-1 results §8 records this as open. Each timed iteration includes `mpsc::sync_channel` alloc, `thread::spawn`, `collector.join()` — none of which is walker work.
- **Pattern matched:** *Microbenchmarks lie about cache behavior* + *No optimization without a benchmark*. The methodology checklist's `iter_batched` rule: setup belongs in the setup closure, not the measured one.
- **Mechanism:** move channel + collector spawn into the setup closure passed to `iter_batched`; keep the walker `run(...)` call in the measured closure. Use `BatchSize::PerIteration` instead of `LargeInput` so per-iteration setup is paid once per measured iteration, not amortised across a batch that includes setup time.
- **Measurement plan:** re-run all eight fixtures; **threshold** is within-run CI shrinks and the L=150 ↔ L=5000 spread narrows further (any remaining gap is real walker scaling). Once the noise floor drops below round-1's ±5–9 %, the smaller deferred items (L11, L14, L17) become measurable.
- **Complexity cost:** ~20 bench lines; no walker change. Diff sketch in [methodology.md](../../tmp/perf_review_2026-05-12_pileup/methodology.md#L171-L195).

#### L3: [target-container/criterion/pileup_walker_*/](../../target-container/criterion/) — **[Likely]** [bench] `final-2026-05-10` saved baseline never created

- **Confidence:** High
- **Hot-path evidence:** verbatim from round-1 results §8. `find` on the criterion tree confirms most-recent saved baseline is `post-l7`; no `final-*` directory exists.
- **Pattern matched:** *Benchmark before and after, on the same machine, in the same configuration. Quote both numbers.*
- **Measurement plan:** one command — `./scripts/dev.sh cargo bench --bench pileup_walker_scaling -- --save-baseline final-2026-05-10`. Snapshots the post-round-1 main state so every round-2 commit can compare against a stable anchor. Run before L1 and before H1.
- **Complexity cost:** none. One bench run, ~5 minutes.

#### L4: [pileup_walker_scaling.rs:1-12, :50-52, :201-215, :232-239](../../benches/pileup_walker_scaling.rs#L1-L12) — **[Likely]** [bench] Bench narrative still stale (re-files round-1 H1)

- **Confidence:** High
- **Hot-path evidence:** round-1 H1 verbatim. Re-reading at the listed line ranges shows the pre-lazy-CIGAR-refactor prose ("eager CIGAR decomposition", "per-active-read clone cost", "every walker step the read is alive sees a Match event", "ceiling the lazy-CIGAR plan would raise") is still present. The lazy-CIGAR refactor shipped before round-1.
- **Pattern matched:** a benchmark whose documented hypothesis contradicts what it measures cannot serve as a regression guard.
- **Suggested fix:** three edits — replace the file docstring with a regression-guard prologue; strike per-active-read clone language at `:50-52` and `:201-205`; drop the "ceiling the plan would raise" at `:212-215`. Round-1 supplied verbatim replacement text in [tmp/perf_review_2026-05-10_pileup/methodology.md](../../tmp/perf_review_2026-05-10_pileup/methodology.md); round-2 [methodology.md](../../tmp/perf_review_2026-05-12_pileup/methodology.md#L291-L310) re-summarises.
- **Complexity cost:** documentation-only; zero code change.

#### L5: [pileup_walker_scaling.rs:267-277](../../benches/pileup_walker_scaling.rs#L267-L277) — **[Likely]** [bench] `config()` is dead-code shadowed by per-group calls (re-files round-1 L3)

- **Confidence:** High
- **Hot-path evidence:** `config()` returns `Criterion::default().sample_size(10).measurement_time(3 s)` and is wired into `criterion_group!`, but both `bench_walker_read_length` and `bench_walker_multi_op` override with per-group `sample_size(30).measurement_time(10 s)`. Runtime behaviour is correct; the file as written claims a configuration it does not use.
- **Pattern matched:** *Build configuration is part of the program* — dead bench configuration is a footgun for the next contributor. The mismatch gets copy-pasted into new benches.
- **Suggested fix:** move the `sample_size(30)` / `measurement_time(10s)` into `config()` and drop the per-group calls (or delete `config()` and use a plain `criterion_group!(benches, …)`).
- **Complexity cost:** negative — three to ten lines deleted.

#### L6: [walker.rs:372](../../src/per_sample_caller/pileup/walker.rs#L372) — **[Likely]** Per-record `SyncSender::send` shape is unbatched; ~8 % walker-CPU aggregate

- **Confidence:** Medium
- **Hot-path evidence:** `SyncSender::send` 4.49 % + `mpmc::Channel::write` 2.86 % + `SyncWaker::notify` 2.84 % + `Context::with` 2.08 % ≈ ~8 % aggregate. Per-record send fires ~49 850 times / run.
- **Pattern matched:** per-call overhead of a sync primitive in a tight hot loop; batching N records amortises the waker + TLS work.
- **Mechanism:** Buffer records into a `Vec<PileupRecord>` (capacity ~8) inside [close_aged_records](../../src/per_sample_caller/pileup/walker.rs#L363) / [flush_chromosome](../../src/per_sample_caller/pileup/walker.rs#L426) and send the whole batch through a `SyncSender<Vec<PileupRecord>>`. Most steps drain ≤ 1 record; the win comes from the multi-record-step case (the wide-deletion-unblocks-several-narrow-records case the [lifecycle-marks doc-comment](../../src/per_sample_caller/pileup/mod.rs#L400-L417) describes).
- **Measurement plan:** feature-gate a batched variant, bench, re-profile to confirm `SyncSender + waker + Context` children-aggregate drops below ~2 %. **Threshold:** ≥ 3 % wall-clock on either fixture.
- **Complexity cost:** the channel type becomes part of the public Stage-2 contract. Stage 2 does not exist yet; pick the batched shape now or carry the per-record shape into the encoder API. Picking now is cheaper than retrofitting. ~20 LOC walker change + matching bench harness change.

#### L7: [pileup_walker_scaling.rs:183](../../benches/pileup_walker_scaling.rs#L183) — **[Likely]** [bench] Two-thread collector hides receiver-side drop cost on a hybrid CPU; production may inherit it

- **Confidence:** Medium
- **Hot-path evidence:** receiver-thread profile (cpu_atom): 44 % `drop_in_place<PileupRecord>` + 40 % `drop_in_place<Vec<AlleleObservation>>` + 36 % `__libc_free` + 34 % `drop_in_place<[AlleleObservation]>` + 20 % `cfree`.
- **Pattern matched:** CPU-bound work on a worker thread that doesn't show up in your wall-clock because it ran on a different core. The bench's collector is *hiding* the drop cost on an E-core; single-threaded production pays the same work on the walker's critical path.
- **Measurement plan:** zero-code-change first — run the bench under `taskset --cpu-list 0` to force walker + collector onto one P-core. Compare wall-clock to the un-pinned bench. **Threshold to act:** if pinned ≥ 10 % slower than unpinned, the gap is what production inherits and the structural fix (batched send + drop on the producer side, or arena-backed `PileupRecord` allele buffers) is justified. If pinned < 5 % slower, the malloc impl handles cross-thread free fine and no action is needed.
- **Complexity cost:** zero for the measurement. Downstream fix (if needed) routes through allocations / data_layout — most of the alloc-cut findings in §5 reduce the receiver-side drop count by the same amount they reduce walker-side alloc count.

#### L8: [open_record.rs:609](../../src/per_sample_caller/pileup/open_record.rs#L609) — **[Likely]** Hoist `affected: Vec<u32>` onto `WalkerState`

- **Confidence:** Medium
- **Hot-path evidence:** DHAT site #6: 49 850 blocks × 16 B = 0.80 MB. One alloc per walker step. Volume cleanup more than CPU; `process_position` already dominates at 60 % and the per-call 16-byte alloc is small inside that.
- **Mechanism:** mirrors round-1 L6 — one new field on `WalkerState` (or a parameter), `clear()` per step. The `.contains()` dedupe at `:640` can become `sort_unstable + dedup` at the end with no behaviour change at this size.
- **Measurement plan:** DHAT site #6 → 0; alloc count down by ~50 K. **Threshold:** ≥ 1 % bench delta or wash.
- **Complexity cost:** one field, one signature change on `open_record::process_position` (the existing destructure-borrow trick used for `allele_seq_buf` extends to this).

#### L9: [open_record.rs:50-55](../../src/per_sample_caller/pileup/open_record.rs#L50-L55) — **[Likely]** `chain_slots` cap-32 mismatches the two allele populations (REF wants ≥ 30, ALT wants 1–2)

- **Confidence:** Medium
- **Hot-path evidence:** DHAT site #2: 89 295 blocks × 92 B = 8.24 MB. Perf children: `insert_sorted_unique` 5.95 %, `Vec::insert` 3.42 %, `RawVec::grow_one (chain_slots)` 2.18 %, `Vec::push` 5.57 %.
- **Pattern matched:** *Pre-size containers when the size is known or bounded* — but the pre-size of 32 is right for one allele class (REF) and ~16× over-provisioned for another (ALT).
- **Mechanism:** Split into `OpenAllele::new_ref(seq)` (`with_capacity(32)`) and `OpenAllele::new_alt(seq)` (`Vec::new()` — cap 0). REF buckets see 30 inserts and fit; ALT buckets see 1–2 inserts and amortise to one allocation. Alternative: `SmallVec<[SlotId; 4]>` for `chain_slots` — inlines the ALT case, REF spills heap at slot 5 — but changes the public type via `AlleleObservation.chain_slots`. Try the split-cap first.
- **Measurement plan:** DHAT site #2 bytes ≈ 8.24 MB → ~3-4 MB. `insert_sorted_unique` 5.95 % should also fall once the ALT path stops growing through 0→4. **Threshold:** ≥ 3 % walker wall-clock on `multi_op/5000`.
- **Complexity cost:** one new constructor pair. Public API unchanged.

#### L10: [cigar_cursor.rs:298, :432, :541, :645](../../src/per_sample_caller/pileup/cigar_cursor.rs#L298) — **[Likely]** `ReadEvent::Insertion::seq` `to_vec()` per insertion allocates ~60 K small blocks

- **Confidence:** Medium
- **Hot-path evidence:** DHAT sites #8 + #9: 29 700 + 29 700 = 59 400 blocks × 2 B avg = 120 KB across four cursor sites (linear + binary, ×2 query shapes).
- **Pattern matched:** *Tiny, often-empty collections are SmallVec candidates*. Round-1 L9 already did this for `EventsAt` and `EventsOverlapping`; the inserted-seq bytes inside `ReadEvent::Insertion` are the last `Vec` in the cursor's return surface.
- **Mechanism:** swap `seq: Vec<u8>` for `SmallVec<[u8; 4]>` (covers typical 1–2 bp inserts inline) or a bespoke `InsertSeq` newtype with inline-8 + boxed spill (smaller for the >4 bp tail).
- **Measurement plan:** DHAT sites #8/#9 → 0 blocks. **Threshold:** ≥ 2 % on `multi_op/5000` (binary-mode cursor; many insertions / read).
- **Complexity cost:** one type change, ripples through the test oracle in [decompose.rs](../../src/per_sample_caller/pileup/decompose.rs) and a few unit tests. `apply_events_to_ref_into` already takes a `&[u8]`, unchanged.

#### L11: [active_read_set.rs:37](../../src/per_sample_caller/pileup/active_read_set.rs#L37) — **[Likely]** `by_read_id: AHashMap<u32, usize>` secondary index → sorted parallel `Vec<(u32, u32)>`

- **Confidence:** Medium
- **Hot-path evidence:** `ActiveReads::get_by_read_id` 2.47 % children — called once per (affected record × contributor) inside the inner fold, ~1.5 M times / run.
- **Pattern matched:** *Inline small keys in maps* — same shape as H1 on the active-set side. ~30 entries at steady state; small enough that `binary_search_by_key` over a contiguous 240 B beats hashbrown's hash + bucket fetch.
- **Mechanism:** parallel `Vec<(u32, u32)>` sorted by `read_id`. **Note:** the read_id stream is monotonically increasing (`next_read_id += 1` per admit), so the insert path is `O(1)` push, not `binary_search + shift`. `swap_remove`'s "moved-entry index fixup" stays — `binary_search_by_key + assignment` (the moved entry's `read_id` doesn't change, so sortedness is preserved).

  *Not the same as round-1 L13.* L13 wanted to project the *read body* (which the per-position scan ends up dereferencing anyway via `events_at`). L11 projects only the secondary index — no read-buffer touch on this path.
- **Measurement plan:** `get_by_read_id` row drops from 2.47 % to ≤ 0.5 %. **Threshold:** ≥ 1 % at p < 0.05 across at least 5 of 8 fixtures.
- **Complexity cost:** same helper shape as H1. `admit`, `expire_passed`, `flush_all` swap to the new helpers. No `unsafe`.

#### L12: [walker.rs:519-523](../../src/per_sample_caller/pileup/walker.rs#L519-L523) — **[Likely]** `resolve_mate_overlap_at_pos` O(n²) duplicate-check fast-path costs ~4 % on the no-overlap edge

- **Confidence:** Medium
- **Hot-path evidence:** `resolve_mate_overlap_at_pos` 3.95 %, with 3.94 % self in the nested `any()`. At n ≈ 30 contributors × 50 K steps ≈ 22 M cmp+branch pairs / run.
- **Pattern matched:** O(n²) where O(n log n) suffices. The comment at [walker.rs:516-518](../../src/per_sample_caller/pileup/walker.rs#L516-L518) explicitly weighs this against the AHashMap-allocation alternative; a sort-based check sidesteps both.
- **Mechanism:** collect the `chain_slot_id`s into a stack-resident `SmallVec<[SlotId; 64]>`, `sort_unstable`, scan `windows(2).all(|w| w[0] != w[1])`. ~30 entries: sort is ~5 cmps amortised per element, no allocation. Alternative: a 4 KB `[u64; 64]` bitmap keyed by slot id (`MAX_ACTIVE_SLOTS = 4096`) hoisted onto `WalkerState`, cleared via a tracked-touched-slots list. The bitmap is bigger memory but the per-step work is purely O(n).
- **Measurement plan:** **threshold ≥ 2 % wall-clock**. If the bench is flat, the existing comment was right; revert.
- **Complexity cost:** ~30 LOC; one new field on `WalkerState` for the reusable buffer. No `unsafe`. Diff sketch at [hot_loops.md](../../tmp/perf_review_2026-05-12_pileup/hot_loops.md#L72-L87).

#### L13: [walker.rs:759-770](../../src/per_sample_caller/pileup/walker.rs#L759-L770) — **[Likely]** Fuse `column_depth_cap`'s any-indel scan into the existing contributors_buf build

- **Confidence:** Medium
- **Hot-path evidence:** `column_depth_cap` 1.91 % children — ~50 K steps × 30 contributors × ≤ 2 events ≈ 3 M `matches!` evaluations.
- **Pattern matched:** *Avoid recomputing in tight loops* — the function re-traverses `events_at_pos` on every contributor running the same `matches!(e, Insertion | Deletion)` test that `pair_has_indel` and `process_position`'s step-3 already touch.
- **Mechanism:** Track `any_indel_at_walker_pos: bool` while building `contributors_buf` in [WalkerState::process_position](../../src/per_sample_caller/pileup/walker.rs#L284-L315); pass it in instead of recomputing. Corner case: `resolve_mate_overlap_at_pos` may `swap_remove` indel losers, so the precomputed flag is "any indel pre-resolution" — the [spec](../../ia/specs/pileup_walker.md) says the cap is applied after resolve, but the bench shows the post-resolve indel-presence equals the pre-resolve set on every column (indel-overlap removes both losers, leaving 0 or ≥1 indels — never flipping the column to all-SNP). Worth pinning in a test before merging.
- **Measurement plan:** `column_depth_cap` row → ~0. **Threshold:** ≥ 1 % wall-clock on `multi_op/5000`.
- **Complexity cost:** one `bool` flag, one signature change on `column_depth_cap`. ≤ 20 LOC.

#### L14: [cigar_cursor.rs:106](../../src/per_sample_caller/pileup/cigar_cursor.rs#L106) — **[Likely]** Re-tune `BINARY_SEARCH_OP_THRESHOLD` now that `lto = "fat"` + `target-cpu = v3` are in effect

- **Confidence:** Medium
- **Hot-path evidence:** at L=5000 binary mode is 24 % of walker (15.50 + 9.05 + 12.94 + 10.64 inside cursor + binary_search_by + partition_point); at L=150 linear mode is ~20 %. The threshold (16) was set pre-LTO-fat and pre-target-cpu-v3 — two compiler effects that move the threshold in opposite directions.
- **Pattern matched:** *Autovectorization needs the compiler's confidence* + *Build configuration is part of the program*. With fat LTO, `emit_event_for_op_overlapping` may inline across the function boundary the manual-inlining comment at [cigar_cursor.rs:217-232](../../src/per_sample_caller/pileup/cigar_cursor.rs#L217-L232) was avoiding — binary mode benefits disproportionately, threshold should drop. With `target-cpu = v3`, the linear walk's branch ladder gets better cmov coverage — linear benefits disproportionately, threshold should rise. A sweep is the only way to tell.
- **Measurement plan:** sweep `BINARY_SEARCH_OP_THRESHOLD ∈ {8, 12, 16, 24, 32, 48}`, bench each. **Threshold to merge:** pick the value that wins on both `multi_op/5000` and `read_length/150` by ≥ 1 % on at least one and doesn't regress the other.
- **Complexity cost:** zero structural — a constant tweak with an updated rationale comment.

#### L15: [cigar_cursor.rs:247-330, :506-572, :359, :629-674](../../src/per_sample_caller/pileup/cigar_cursor.rs#L247-L330) — **[Likely]** Hoist `assert_eq!(offsets.len(), cigar.len() + 1)` on the four cursor index loops to elide bounds checks

- **Confidence:** Low
- **Hot-path evidence:** `events_overlapping` + `events_at` (binary + linear) together are ~25 % of walker CPU. The inner `for i in 0..n_ops` loops index `self.offsets[i]` and `read.cigar[i]`; rustc may or may not have proved both bounds equivalent via `n_ops = read.cigar.len()`.
- **Pattern matched:** *A leading `assert!(idx < slice.len())` (or `assert_eq!(slice.len(), N)`) hoists the bounds check out of the loop body.*
- **Mechanism:** one `assert_eq!` per function lets LLVM elide the per-iteration cmp/branch. The cursor invariant `offsets.len() == cigar.len() + 1` is true by construction but the compiler can't always see it.
- **Measurement plan:** `cargo asm --release` on the four cursor functions before / after; count `panic` / `bounds_check_failed` labels. **Threshold to merge:** ≥ 1 instruction per iteration eliminated AND ≥ 0.5 % bench delta. If the codegen already has zero panic labels, close as no-op.
- **Complexity cost:** four asserts + matching `let` bindings. Zero `unsafe`. Asserts run once per call (outside the loop).

#### L16: [cigar_cursor.rs:262-286, :398-419](../../src/per_sample_caller/pileup/cigar_cursor.rs#L262-L286) — **[Likely]** Confirm or refute L16 autovec on the Match-emit `for k in 0..span` body with `cargo asm`

- **Confidence:** Low
- **Hot-path evidence:** Match-emit body lives inside the 15.5 % / 9.4 % `events_overlapping` / `events_at` cost centres. Round-1 deferred this experiment because `cargo-show-asm` wasn't in the dev container; it now is.
- **Pattern matched:** *Autovectorization needs the compiler's confidence.* The per-base loop has two early `continue`s (read-N skip + adaptor-boundary skip) plus paired `bq_baq` + `seq` indexed loads; full autovec is unlikely but residual bounds checks may still be elidable.
- **Measurement plan:** `cargo asm` on `emit_event_for_op_overlapping` + `events_overlapping_linear` + `events_at_linear`. Look for panic labels in the inner `for k in 0..span` block. If present, add `assert_eq!(read.seq.len(), read.bq_baq.len())` at function entry (the invariant is already enforced by [PreparedRead::length](../../src/per_sample_caller/pileup/mod.rs#L307) at admit time but the cursor doesn't know that statically). If absent, close the round-1 L16 finding with a comment.
- **Complexity cost:** measurement only. Fix (if needed) is two LOC.

#### L17: [mod.rs:215](../../src/per_sample_caller/pileup/mod.rs#L215) — **[Likely]** `adaptor_boundary: Option<u32>` → `Option<NonZeroU32>` collapses 4 B of padding

- **Confidence:** Medium
- **Hot-path evidence:** `base_in_adaptor` is called on every Match emit inside the cursor (`events_at` + `events_overlapping`), order of (active reads × walker positions) ≈ hundreds of millions per run. Field lives inside the 108-byte `PreparedRead` and is touched on the hot loop.
- **Pattern matched:** *Niche optimization is real.* `Option<u32>` is 8 B (4 discriminant + 4 payload); `Option<NonZeroU32>` is 4 B with the same null semantic. 1-based positions start at 1; `0` is always invalid — niche is free.
- **Measurement plan:** print `size_of::<PreparedRead>()` before / after (expect 4 B drop). Bench: standalone impact likely below noise; file as cleanup, do not block on its measurement gate.
- **Complexity cost:** producer side wraps via `NonZeroU32::new(...)`; consumer side calls `.get()`. One line each.

### Speculative

- **S1: [Cargo.toml](../../Cargo.toml)** — `rust-toolchain.toml` is still absent (round-1 S1). Pin only if measured numbers stop reproducing across rustc versions; the round-2 numbers were taken on rustc 1.93.
- **S2: [Cargo.toml](../../Cargo.toml)** — `panic = "abort"` and PGO unexplored. `panic = "abort"` is a one-line change but removes unwind landing pads — audit `lib.rs` for `catch_unwind` first. PGO defer until the lab has a fixed representative WGS input in CI.
- **S3: [pileup_walker_scaling.rs](../../benches/pileup_walker_scaling.rs)** — Streaming-WGS-span fixture missing. The bench's span = 50 000 bp fits in L2; production spans are 5000× larger. Add a `bench_walker_span_streaming` group that sweeps span at fixed L=150, coverage=30 after L2 lands.
- **S4: [mod.rs:57](../../src/per_sample_caller/pileup/mod.rs#L57)** — `DEFAULT_OUTPUT_CHANNEL_CAPACITY = 64` has no measurement backing; the walker-only bench cannot exercise back-pressure. Defer until Stage 2 lands; the right answer is "instrument both ends and sweep capacity once a real consumer exists".
- **S5: [mod.rs:237](../../src/per_sample_caller/pileup/mod.rs#L237)** — `Arc<str>` qname clones (round-1 L17/L18) are still dormant on the solo-read fixture. Build the paired-end fixture first; only then file a fresh finding if the atomics show up.
- **S6: [active_read_set.rs:18](../../src/per_sample_caller/pileup/active_read_set.rs#L18)** — *narrowed* re-evaluation of round-1 L13: project only `(read_id, chain_slot_id, alignment_end)` into a parallel `Vec<ActiveReadHeader>` to fast-skip `expire_passed`. Differs from L13 because this projection does *not* touch read buffers — the `events_at` objection that killed L13 does not apply. Still Speculative until a microbench isolates the loop.
- **S7: [cigar_cursor.rs:67](../../src/per_sample_caller/pileup/cigar_cursor.rs#L67)** — `EventsOverlapping = SmallVec<[ReadEvent; 4]>` inline cap of 4. Halve to 2? At `ReadEvent` = 32 B, the inline buffer is 128 B (two cache lines); reducing the cap pays one heap alloc on the rare 3+ event return. Measure. Note: see also L19 below — the same enum width is a confounder.
- **S8: [mod.rs:456-462](../../src/per_sample_caller/pileup/mod.rs#L456-L462)** — `RefSeqFetcher::fetch -> Vec<u8>` allocates per call. Add `fetch_into(&self, …, out: &mut Vec<u8>)` with a default impl. Worth doing as a trait-shape adjustment now (low complexity) but the win on the bench is 49 850 × 1-byte allocs (negligible); revisit once the production noodles-backed fetcher's profile share is known.
- **S9: [decompose.rs:15](../../src/per_sample_caller/pileup/decompose.rs#L15)** — Box / inline `ReadEvent::Insertion.seq` so the enum shrinks 32 → 24 B; `EventsOverlapping` inline buffer drops 136 → 104 B. Low confidence given round-1 L14 (the prior cache-locality fusion that regressed 7.6 %). Pattern says measure; don't ship blind.

### Note (no action)

- **Receiver-thread drop cost** (44 % `drop_in_place<PileupRecord>` + 40 % `Vec<AlleleObservation>` + 36 % `__libc_free`) — most of the alloc-cut findings above (H1, H2, L9, L10) reduce this in lockstep with walker-side allocs. Not separately actionable; credit the receiver-side win when picking which finding to ship first.
- **`find_allele_index`** ([open_record.rs:576-578](../../src/per_sample_caller/pileup/open_record.rs#L576-L578)) — 1.91 % linear `seq.eq` scan over ≤ 4 alleles; the iterator chain compiles to a tight loop. A hash map at this size would lose to the linear scan. No change.
- **`slot_refcount: Vec<u8>` of 4 K entries** ([slot_allocator.rs:70](../../src/per_sample_caller/pileup/slot_allocator.rs#L70)) — 4 KB hot region. Indexed access pattern is right; only `active_count()` walks it (round-1 L15, cleanup-only). Filed so a future review doesn't propose packing as a bitmap (the bytes are written, not just read).
- **Channel is bounded** ([walker.rs:6](../../src/per_sample_caller/pileup/walker.rs#L6)) — positive sighting; `mpsc::sync_channel(N)`, not the unbounded `mpsc::channel()`. Rule passes.

## 6. Out-of-scope observations

- **[find_overlapping](../../src/per_sample_caller/pileup/open_record.rs#L250-L280)** does `BTreeMap::range(lo..=event_start).rev()` per event per step at 4.27 % children (3.16 % `BTreeMap::range` self). `BTreeMap` has known cache-unfriendliness vs a small sorted Vec at typical record counts (1–4). Cross-category note from `hot_loops`; investigate after H1 + H2 land. Filing as out-of-scope because the typical 1-candidate case makes the alternative win small relative to the AHashMap fixes above.
- **Cursor manual inlining at [cigar_cursor.rs:217-232](../../src/per_sample_caller/pileup/cigar_cursor.rs#L217-L232)** — the comment notes the per-op helper was extracted under a 15-20 % regression that LTO might now subsume. L14 (threshold re-tune) is the experiment that would also reveal whether the manual inlining is still needed; if the threshold sweep shows binary mode catching up to linear, consider whether the linear-mode duplication can collapse back to the shared helper. Pattern-match only.
- **[walker.rs:534](../../src/per_sample_caller/pileup/walker.rs#L534)** `AHashMap<SlotId, Vec<usize>>` in `resolve_mate_overlap_at_pos`'s slow path is dormant on the solo-read bench. If L12 ships, the fast-path bypasses it; if the paired-end fixture (S5) shows the slow path is hot, hoist the slow-path containers onto `WalkerState`.

## 7. What's already good

- **Round-1 results held intact through a month of renames.** DHAT post-experiment numbers (504 K allocs, 159 MB) match the saved `post-l7` baseline exactly. The reset semantics on [OpenPileupRecordTable::reset](../../src/per_sample_caller/pileup/open_record.rs#L195-L203) preserve the perf-hoisted `allele_seq_buf` across chromosome boundaries — round-1's Mi11 fix is doing what it was meant to.
- **L7 fast-path is doing its job** — [walker.rs:519-526](../../src/per_sample_caller/pileup/walker.rs#L519-L526) early-returns on no-overlap, and the dormant AHashMap allocation at `:534` does not appear in the round-2 perf trace. That's exactly the "watch this path on paired-end" gate round-1 set.
- **Cursor's adaptive linear/binary mode dispatch** at [cigar_cursor.rs:90-126](../../src/per_sample_caller/pileup/cigar_cursor.rs#L90-L126) works: at L=150 (5 ops) the linear path takes ~12 % combined, at L=5000 (199 ops) the binary path takes ~25 % combined. Both are textbook hot-loop work; the open question is whether the threshold (set pre-LTO) is still the right crossover (L14).

---

## Per-category audit trail

Sub-agent output is preserved at [tmp/perf_review_2026-05-12_pileup/](../../tmp/perf_review_2026-05-12_pileup/):

- [methodology.md](../../tmp/perf_review_2026-05-12_pileup/methodology.md) — L1, L2, L3, L4, L5, S1, S2, S3
- [allocations.md](../../tmp/perf_review_2026-05-12_pileup/allocations.md) — H1 (alternative shape), H2, L8, L9, L10, S8
- [data_layout.md](../../tmp/perf_review_2026-05-12_pileup/data_layout.md) — H1 (canonical shape), L11, L17, S6, S7, S9
- [concurrency.md](../../tmp/perf_review_2026-05-12_pileup/concurrency.md) — L6, L7, S4, S5
- [hot_loops.md](../../tmp/perf_review_2026-05-12_pileup/hot_loops.md) — H3, L12, L13, L14, L15, L16

Profile inputs:
- [profile_summary.md](../../tmp/perf_review_2026-05-12_pileup/profile_summary.md) — synthesised summary, the verbatim source the sub-agents quoted from
- [perf_multi_op_5000.children.txt](../../tmp/perf_review_2026-05-12_pileup/perf_multi_op_5000.children.txt), [perf_multi_op_5000.flat.txt](../../tmp/perf_review_2026-05-12_pileup/perf_multi_op_5000.flat.txt) — raw `perf report` output, L=5000 multi-op
- [perf_read_length_150.children.txt](../../tmp/perf_review_2026-05-12_pileup/perf_read_length_150.children.txt) — raw `perf report`, L=150/L=1500 read-length
- [dhat-heap.json](../../tmp/perf_review_2026-05-12_pileup/dhat-heap.json), [dhat_top.txt](../../tmp/perf_review_2026-05-12_pileup/dhat_top.txt) — DHAT data + top-by-block listing

## Round-1 status (carry-forward)

The round-1 audit trail at [perf_pileup_2026-05-10_results.md](perf_pileup_2026-05-10_results.md) records 7 changes applied, 1 reverted, 5 deferred. Round-2 maps the deferred items as follows:

| Round-1 finding | Round-2 disposition |
|---|---|
| L12 (drain_aged empty-case short-circuit) | Not refiled — record age-out rate makes the short-circuit moot. H2 (buffer hoist) addresses the alloc cost via a different mechanism. |
| L13 (ActiveRead DoD projection) | Not refiled in original shape. Narrower variant filed as S6 (`expire_passed` only). The active-set side of the AHashMap pattern is filed as L11. |
| L14 (OpRow fusion) | Stays reverted. The +7.6 % regression result is referenced in data_layout's [L19 caveat](../../tmp/perf_review_2026-05-12_pileup/data_layout.md). |
| L15 (`active_count` incremental counter) | Still cleanup-grade; not refiled. |
| L16 (Match-emit autovec) | Filed as L16 in round-2 — `cargo-show-asm` is now in the dev container; one-shot measurement decides. |
| L17 / L18 (`Arc<str>` qname clones) | Still dormant on solo-read fixture; build paired-end fixture first (S5). |
| L19 (`#[inline]` on small helpers) | Moot under `lto = "fat"`; not refiled. |

## Author response convention

Address each finding by its identifier with one of:
- `applied in <commit>`
- `experiment shows no gain — closing`
- `disputed because …`
- `deferred to <issue>`
- `won't fix because …`

The "experiment shows no gain" outcome is welcome and expected — round-1 reverted L14 on exactly that gate.
