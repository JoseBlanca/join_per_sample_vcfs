# Performance Review: bed-regions (`--regions`)
**Date:** 2026-06-16
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the BED `--regions` feature path (per-sample pileup Stage 1 + cohort var-calling restriction) as it exists after the segment-read-fetcher merge
**Verdict:** Profile first → **measured (2026-06-16)**: per-region overhead is negligible (≤2%); BAQ (`probaln_glocal`) dominates ~70% of pileup self-time. The `--regions` findings are **not** the lever; L1+L3 applied as perf-neutral cleanups. See §8.
**Hot-path evidence:** at review time none; **now** a fragmentation sweep + sampling profile on the tomato1 CRAM cohort (§8)

---

## 1. Scope and constraints

- **What was reviewed:** the `--regions` code path as it exists *now*, after the
  segment-read-fetcher merge — read fetch, the per-region pileup loop, the k-way
  segment merge, the per-column write clamp, and the cohort interval restriction.
- **Reviewed against:** branch `bed-regions-review` @ `3b843a6` (off `main`).
- **Targets / inputs / hardware:** cohorts of tens of samples; BEDs from tens to
  thousands of regions; BAM and CRAM input. Dev = macOS arm64 (~6 fast cores);
  production = Linux aarch64, up to 32 cores. The feature trades CPU for memory
  (the *superseded* query-path impl measured −71% peak RSS vs whole-file
  streaming at ~unchanged wall on a 1000-region GRCh38 HG002 run); the bar is
  "per-region fixed overhead stays negligible vs the per-position pileup+BAQ math,
  and the RSS win is not regressed."
- **Hot-path evidence available:** **none for the current code.** The only prior
  perf report ([bed_regions_perf_and_byte_identity_2026-06-03.md](../../bed_regions_perf_and_byte_identity_2026-06-03.md))
  measured the now-deleted `AlignmentMergedReader::query` path. No benchmark
  exercises `--regions`; no sampling profile of the `SegmentMergedReads` path
  exists. **All code-level findings below therefore cap at Likely per the rubric.**
- **In-scope files:**
  [src/regions.rs](../../../../src/regions.rs),
  [src/pop_var_caller/cli.rs](../../../../src/pop_var_caller/cli.rs),
  [src/pop_var_caller/stage1_pipeline.rs](../../../../src/pop_var_caller/stage1_pipeline.rs),
  [src/bam/segment_reader.rs](../../../../src/bam/segment_reader.rs),
  [src/bam/segment_merge.rs](../../../../src/bam/segment_merge.rs),
  [src/pileup/per_sample/pileup_to_psp.rs](../../../../src/pileup/per_sample/pileup_to_psp.rs),
  [src/var_calling/pipeline.rs](../../../../src/var_calling/pipeline.rs).
- **Deliberately out of scope:** the SNP posterior engine math, the SSR caller,
  the BAQ HMM internals, and noodles decode internals (flagged where they force
  whole-container decode, but not reviewed).
- **Categories dispatched:** methodology (always); allocations (per-region churn);
  io_and_syscalls (CRAM index/decode — the headline suspect); concurrency (serial
  region loop, reader-pool Mutex); hot_loops (per-record merge/filter); data_layout
  (per-record merge struct size).

> **Stale status note (acted on in §8):** `PROJECT_STATUS.md` records the
> segment reader as "no consumer wired yet" and "retrofit the SNP `--regions`
> path" as open (#5). The code disproves this:
> [cli.rs:414](../../../../src/pop_var_caller/cli.rs#L414) constructs
> `SegmentMergedReads`, and the comment at
> [cli.rs:374-375](../../../../src/pop_var_caller/cli.rs#L374-L375) states it
> replaces the old per-region `AlignmentMergedReader::query`. The retrofit has
> landed; the status file is behind.

## 2. Verdict

**Profile first.** There is not enough hot-path evidence to recommend code
changes. The feature's entire stated concern — per-region fixed overhead at scale
(fragmented BED, CRAM-sparse, many-contig) — sits on the timed path of *no*
committed benchmark, and the current read path has never been sampling-profiled.
The first deliverable is the measurement infrastructure in §3, not code edits.

Two caveats temper the "first" — both contained, both safe to act on *without* a
profile because the mechanism is a removed allocation, not a speculative re-design:

- **L1 (qname double-clone)** is a true per-read inner-loop allocation that three
  independent categories flagged; collapsing the double clone is low-risk.
- **L3 (dead per-region `Stage1Outputs` clones)** is pure waste the caller never
  reads.

Everything else — and especially the CRAM and concurrency findings — is gated on
the benchmark + profile.

## 3. Measurement plan

The primary deliverable. In dependency order:

1. **Add `benches/regions_pileup_perf.rs`** (criterion, `harness = false`,
   registered in `Cargo.toml`). Drive `run_pileup` (or, to isolate the seam,
   `with_stage1_chain` over a `SegmentMergedReads`) on a **span-count axis at a
   fixed total covered footprint**, so the only variable is per-region overhead:
   one synthetic contig covered by `{1, 10, 100, 1000}` equal-area BED spans, the
   1-span case as baseline. Cross with **format `{BAM, CRAM}`** (CRAM is the
   suspected worst case) and a **many-contig** axis (regions scattered across many
   contigs, to exercise the FASTA contig-transition reload). Build the fixture
   once outside `b.iter`; assert `RunSummary.records_emitted` is constant across
   span counts at fixed footprint (a fixture/region-resolution regression then
   trips the bench). Reuse the `build_cram`/`build_fasta` helpers already used by
   the segment_reader tests.
   - **Answers:** does per-region wall grow super-linearly with span count at
     fixed footprint? Threshold to act on any code finding: yes.
2. **Sampling-profile the heaviest combo** against a `--profile profiling` build
   (the `[profile.profiling]` profile already exists). `samply` on the macOS host
   binary, or `cargo flamegraph -c "record -e cpu-clock ..."` on the Linux dev box.
   - **Answers:** the inclusive-time split across `SegmentMergedReads::new`,
     `CramSegmentReads::refill` (`.crai` walk + container decode), the
     `with_stage1_chain` chain rebuild, and `repository.clear`/FASTA reload.
     Threshold: any site ≥ ~10% inclusive time on the heaviest combo is an
     actionable candidate; this is what promotes L2/L4/L5 above Likely.
3. **(CRAM-specific) cluster sweep** — same fragmented-CRAM fixture, regions
   spread one-per-container vs packed into few containers, at equal region count.
   If clustered is not cheaper, the per-region container re-decode (L4) is real and
   a last-container cache should close the gap.
4. **(concurrency) starvation probe** — instrument `produce_packets` send-blocked
   vs worker recv-empty counters (sketch in the concurrency category file) on a
   thousand-tiny-region BED. recv-empty ≫ send-blocked ⇒ the per-region pipeline is
   underfed and region-level parallelism (L7) is worth its complexity; otherwise
   it is not.

## 4. Build / toolchain configuration

**Sound — no change before code-level work.** `[profile.release]` =
`lto = "fat"`, `codegen-units = 1`, `panic = "abort"`,
`debug = "line-tables-only"`; a dedicated `[profile.profiling]`
(`inherits = "release"`, `lto = false`, `codegen-units = 16`, `debug = true`)
exists for symbolicated stacks; `.cargo/config.toml` pins `target-cpu` to the two
named deployments (`x86-64-v3`, `apple-m1`), not `native`; the toolchain is pinned
(`1.95`); an `alloc-mimalloc` A/B feature is wired. The only follow-on build
experiment is to re-run the new bench under `--features alloc-mimalloc` once it
exists (gated on the bench, not a fix).

## 5. Code-level findings

All cap at **Likely** (no profile names any site). Ordered by confidence, then
convergence across categories.

### L1 — [src/bam/segment_merge.rs:232-233,260-261](../../../../src/bam/segment_merge.rs#L232-L261) — qname cloned twice per merged read (the per-read inner loop)
- **Confidence:** Medium. Flagged independently by **allocations**, **hot_loops**, and **data_layout** — the strongest convergence in the review.
- **Hot-path evidence:** pattern-match only. `SegmentMergedReads::next` runs once per surviving read across every region — the true inner loop of the feature.
- **Mechanism:** line 233 does `read.qname.clone()` to read the head's keys while only a `&` peek is held; line 261 clones the same `qname` *again* into `new_key: ReadFingerprint`. Two `Vec<u8>` heap allocations + byte copies per emitted read. The second clone is paid unconditionally even on the common single-input case, where cross-file dedup is impossible.
- **Measurement plan:** micro-bench `SegmentMergedReads` draining N reads with ~20-byte qnames, single-file and 4-file, current vs one-clone variant; metric = allocs/read (dhat) + ns/read. `cargo asm` to confirm two `alloc` calls collapse to one. Merge on any non-noise reduction.
- **Complexity cost:** Low–moderate. Collapse to one clone by building `new_key` from the already-cloned local `qname` (move on the push path); additionally short-circuit the `current_locus_fingerprints` push/scan when `streams.len() == 1`. Keep the existing dedup tests as the guard. (L5's cache-key projection would further remove the line-233 clone for the in-order case.)

### L2 — [src/bam/segment_reader.rs:784-817](../../../../src/bam/segment_reader.rs#L784-L817) — CRAM `.crai` walk restarts from the file head every region
- **Confidence:** Medium.
- **Hot-path evidence:** pattern-match only. `CramSegmentReads` is seeded `next_index_record: 0` and `refill` linearly scans the flat `crai::Index = Vec<Record>` (one entry per container), `continue`-skipping all earlier-contig/earlier-container records — on every region, for every input file.
- **Mechanism:** a region on contig *k* at offset *o* re-scans every index record for contigs `0..k` plus every container before *o*, each call. A fragmented BED on a late contig turns an O(regions) walk into O(regions × index_head). This is exactly the deferred PROJECT_STATUS "binary-search the crai head" item.
- **Measurement plan:** multi-contig CRAM with reads on a late contig + fragmented BED (N = 100/1000/5000); profile the `refill` index-scan share as N and target-contig lateness grow. Confirm a `partition_point` seed flattens it.
- **Complexity cost:** Low. One `partition_point` over `index.as_slice()` to seed the cursor at the first record on the target contig that can overlap `segment.start`; keep the per-record overlap re-check as the correctness backstop; validate emitted-read byte-identity against the linear walk. Care needed at contig boundaries (land on the target contig's first record).

### L3 — [src/pop_var_caller/stage1_pipeline.rs:261-263](../../../../src/pop_var_caller/stage1_pipeline.rs#L261-L263) — per-region `ContigList::clone()` + `sample_name.to_string()` the caller never reads
- **Confidence:** Medium.
- **Hot-path evidence:** pattern-match only. `with_stage1_chain` runs once per region; each return builds `Stage1Outputs` with `contigs.clone()` (deep-clones the whole `Vec<ContigEntry>`, one heap `String` per contig) and `sample_name.to_string()`. The cli.rs caller ([cli.rs:435-449](../../../../src/pop_var_caller/cli.rs#L435-L449)) reads only `result`, `run_summary`, `stashed_upstream_error` — never `contigs`/`sample_name`.
- **Mechanism:** ~`n_regions × (n_contigs + 1)` pure-waste allocations for a whole-genome reference. Dropping the two fields removes them; the caller already owns its own copies.
- **Measurement plan:** dhat on a ~1000-region multi-contig run; confirm the `ContigList`/`String` clone frames vanish.
- **Complexity cost:** Low. Remove two fields from `Stage1Outputs`; verify no other workspace caller reads them (one tracked struct).

### L4 — [src/bam/segment_reader.rs:821-877](../../../../src/bam/segment_reader.rs#L821-L877) — CRAM decodes whole containers for sparse/clustered regions
- **Confidence:** Medium (mechanism), but the fix is non-trivial.
- **Hot-path evidence:** pattern-match only. Per overlapping container, `refill` decodes every slice/record then discards those outside the segment — intrinsic to CRAM granularity (noodles forces it). The project-side lever is the deferred "CRAM container cache": consecutive regions in the sorted `RegionSet` that fall in the same container re-seek and re-decode it once per region, because each `get_reads_from_segment` builds a fresh stateless `CramSegmentReads`.
- **Measurement plan:** the §3.3 cluster sweep (one-per-container vs packed). Merge a last-container cache if the clustered case shows >~10% wall in repeated decode of the same container offset.
- **Complexity cost:** Non-trivial. A per-handle `(last_container_offset, Vec<MappedRead>)` cache re-introduces mutable state into the deliberately-stateless pooled reader and must not double-count `filter_counts`. File as a measured follow-up, not a quick win.

### L5 — [src/bam/segment_merge.rs:190-201](../../../../src/bam/segment_merge.rs#L190-L201) — `argmin_head` dereferences a >100-byte `MappedRead` per peek to read only `(ref_id, pos)`
- **Confidence:** Medium.
- **Hot-path evidence:** pattern-match only. `argmin_head` runs once per emitted read over all streams (`O(n_files)` peeks); each peek pointer-walks into a `MappedRead` (four `Vec`s = 96 B + scalars, straddles ≥2 cache lines) but reads only two fields.
- **Mechanism:** the trailing `Vec` headers/scalars are pulled into cache for nothing on every comparison. Caching a 16-byte `head_key: Option<(usize,u64)>` in `PeekableSegmentStream`, refreshed in `peek()`, lets argmin iterate compact keys that pack many-per-line. Synergy: it also removes L1's line-233 clone for the in-order case (the order check then needs only the key).
- **Measurement plan:** restructure + micro-bench the merge on a 4–8-file cohort over a 1000-region BED (wall/merged-read); on macOS, Instruments Time Profiler to confirm the argmin frame's miss share drops. Merge only at `n_files >= 4` with a measurable improvement. (No HW cache counters on this machine — Instruments is the substitute.)
- **Complexity cost:** Low–moderate. One field + a "key mirrors head" invariant `peek()`/`next()` must maintain.

### L6 — [src/bam/segment_merge.rs:135-152](../../../../src/bam/segment_merge.rs#L135-L152) — `SegmentMergedReads::new` copies every input path per region
- **Confidence:** Low.
- **Hot-path evidence:** pattern-match only. Called once per region; line 143 does `file.path().to_path_buf()` for every input file — `n_regions × n_files` heap path copies whose contents never change. `paths` is read only to format error messages.
- **Mechanism:** the struct is already `'a`-parameterised over the borrowed `&'a [AlignmentFile]`, so `paths: Vec<&'a Path>` (or reading `files[idx].path()` at the error sites) removes the copies. Only material at the fragmented-BED × multi-input intersection.
- **Measurement plan:** dhat on a fragmented-BED multi-file run; look for `to_path_buf` allocations scaling with `n_regions × n_files`.
- **Complexity cost:** Low. One borrow; no `'static`, no dependency.

### L7 — [src/pop_var_caller/cli.rs:401-450](../../../../src/pop_var_caller/cli.rs#L401-L450) — pileup per-region loop is serial; region-level parallelism is available but trap-laden
- **Confidence:** Medium. **Measure-first; do nothing without evidence.**
- **Hot-path evidence:** pattern-match only. The loop is serial, but parallelism already lives *inside* each region: `with_stage1_chain` spawns a producer + `n_workers` threads at `--threads >= STAGED_MIN_THREADS`. So the cores are already saturated per region.
- **Mechanism:** regions are independent at the read level (the pool is `Sync`, proven by the `*_parallel_segments_match_sequential` tests), but two things force serialism: (1) all regions write one `&mut PspWriter` that must emit in coordinate order; (2) each region already fans out to `n_threads` workers, so a naive `par_iter` over the loop body yields `W × n_threads` threads — textbook two-pool oversubscription (the codebase guards against exactly this at [cli.rs:278](../../../../src/pop_var_caller/cli.rs#L278)). The *only* regime where outer region-parallelism helps is **many tiny regions** that underfeed the per-region pipeline.
- **Measurement plan:** the §3.4 starvation probe. Only if worker recv-empty ≫ send-blocked on a thousand-tiny-region BED is there parallelism to recover.
- **Complexity cost:** High. A correct fix caps total threads (`W` regions × `n_threads/W` workers) *and* adds a region→writer coordinate-ordered merge; must be gated on `.psp` byte-identity. Not worth attempting before the probe.

### Speculative

- **S1 — [src/bam/segment_reader.rs:426](../../../../src/bam/segment_reader.rs#L426)** — `pre_decode_config()` rebuilds a small `Copy` filter config per record (loop-invariant). LLVM may already sink/fold it. Confirm with `cargo asm`; cache on the iterator only if not elided.
- **S2 — [src/bam/segment_merge.rs:112,265-268](../../../../src/bam/segment_merge.rs#L265-L268)** — per-locus dedup is a linear scan whose `ReadFingerprint == ` compares the `Vec<u8>` qname first; the buffer is near-empty in sorted data (cleared per locus). If ever shown hot (deep same-base pileup, many files), reorder the fingerprint fields integers-first to short-circuit before the byte compare, or split AoS→SoA.
- **S3 — [src/bam/cram_input.rs:34-44](../../../../src/bam/cram_input.rs#L34-L44)** — the CRAM reader is an unbuffered `Reader<File>` (no `BufReader`/`mmap`). Softened because each container is one bulk read (not per-record). Gate on a `strace -e pread64,lseek -c` count before pursuing; `memmap2` adds `unsafe` + SIGBUS caveats not worth it absent evidence. (The BAM path is correctly bgzf-block-buffered — *not* a finding.)
- **S4 — [src/bam/segment_reader.rs:509,517](../../../../src/bam/segment_reader.rs#L509-L517)** — the per-file reader-pool `Mutex` is contention-free in the current serial wiring (lock held only for `Vec` pop/push, never across I/O) but would become a hot spot if L7's region-parallelism lands (`W` regions × `n_files` hammering the same per-file locks). If measured contended, shard the pool or use `crossbeam::queue::SegQueue` — **not** `parking_lot` (the section is short/uncontended, where `std::sync::Mutex` wins).

## 6. Out-of-scope observations

- **[src/var_calling/pipeline.rs:643-716](../../../../src/var_calling/pipeline.rs#L643-L716)** — the cohort side runs one DUST mask + REF-fetcher build per scheduled interval, and `--regions` multiplies the interval count. Whether per-interval fixed cost scales badly with a fragmented BED is worth a separate look (the cohort var-calling path, not the pileup path reviewed here). Surfaced by the methodology cross-category note.
- **[src/bam/segment_reader.rs:614](../../../../src/bam/segment_reader.rs#L614)** — the BAM `next()` loop builds a fresh `RecordBuf::default()` per record rather than reusing a scratch buffer. This is a general read-path allocation pattern (broader than `--regions`); fold into the scratch-buffer workstream, not this review.

## 7. What's already good

- **Reader-pool lock discipline** ([segment_reader.rs:509-518](../../../../src/bam/segment_reader.rs#L509-L518)) — the `Mutex<Vec<Handle>>` is locked only for pop/push, never held across the index query, seek, decode, or iteration. No lock-across-I/O bug.
- **FASTA repository is one-load-per-contig, not per-region** ([cli.rs:401-408](../../../../src/pop_var_caller/cli.rs#L401-L408)) — the sorted `RegionSet` + clear-on-contig-transition means each contig's sequence loads at most once and stays resident for all its regions. This *refutes* the brief's "reload per region" concern; the design is correct.
- **BAM I/O is correctly buffered** — the path rides `bgzf::io::Reader`'s internal block buffer and avoids the outer-`BufReader` trap bgzf forbids; static `enum { Bam, Cram }` dispatch on the per-read path (no `dyn`/vtable).
- **The rayon-pool sizing guard** ([cli.rs:278](../../../../src/pop_var_caller/cli.rs#L278)) deliberately skips a global rayon pool on the staged path to avoid `n` idle threads beside the pipeline's own `n` workers — the existing defense against the oversubscription L7 would risk.

### Author response convention
Address each finding by id (L1, L2, …, S1, …) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. The measurement plan exists so "no gain — closing" is a
welcome outcome.

## 8. Empirical results (2026-06-16) — the measurement was run

Measured on the **tomato1** benchmark (26→63-sample S. lycopersicum CRAM cohort,
CRAMs pre-sliced to 80×100 kb regions ≈ 8 Mb of SL4.0), single CRAM
`SRR7279481.p1`, `--threads 4`, `--profile profiling` host build (arm64).

### 8.1 Fragmentation sweep — per-region overhead (hyperfine, warmup 1, 3 runs)

Same 8 Mb covered bases, split into more regions to isolate per-region fixed cost:

| BED | regions | coverage | wall (mean ± σ) | vs x1 |
|---|---|---|---|---|
| x1   | 80   | 8 Mb | 13.455 s ± 0.067 | — |
| x10  | 800  | 8 Mb | 13.641 s ± 0.034 | +1.4% |
| x100 | 8000 | 8 Mb | 13.705 s ± 0.042 | **+1.9%** |

**100× more regions over identical coverage costs ~1.9% wall.** Per-region fixed
overhead is negligible on this workload. BED generator + BEDs:
`tmp/regions_bench/` (gen_frag_bed.py, regions_x{1,10,100}.bed).

### 8.2 Sampling profile (macOS `sample`, 12 s, single CRAM) — flat-leaf ranking

- `baq::probaln::probaln_glocal` — **7473 self-samples** (the dominant leaf, ~70%+)
- `pileup::walker::open_record::process_position` — ~839+ (the next chunk)
- `baq::scratch::ProbalnScratch::resize_for` — ~263
- `probaln` appears in **8064** stack lines; `segment_merge`/`SegmentMergedReads`
  in **151** (~53× fewer). `qname`/`ReadFingerprint`/`clone` did not register as
  named self-time nodes.

**The entire `--regions` read/merge path — where every finding below lives — is
~1–2% of runtime. BAQ dominates.** Profile: `tmp/regions_bench/sample_x1.txt`.

### 8.3 Per-finding disposition

- **L1 (qname double-clone) — `applied` (perf-neutral).** Collapsed the second
  clone to a move ([segment_merge.rs:259-264](../../../../src/bam/segment_merge.rs#L259-L264)).
  Re-measure: x1 13.80 s / x100 13.84 s — within noise of baseline, as expected
  for a ~1% path. Kept because it is a genuine (if invisible) dead allocation and
  the dedup tests still guard it; **not** a speedup.
- **L3 (dead `Stage1Outputs` clones) — `applied` (perf-neutral).** Dropped the
  unread `sample_name`/`contigs` fields
  ([stage1_pipeline.rs](../../../../src/pop_var_caller/stage1_pipeline.rs)); the
  one caller never read them. Cleanup, not a speedup.
- **L2 / L4 (CRAM `.crai` head-walk, container re-decode) — `experiment N/A on
  this data; untestable here`.** The tomato CRAMs are pre-sliced (tiny `.crai`,
  few containers), so the head-walk / re-decode scenarios these target don't
  reproduce. Plausible only on a **full, un-sliced CRAM + fragmented BED**, which
  this benchmark cannot provide. Left open as deferred (already mirrored by
  PROJECT_STATUS "binary-search the crai head" / "container cache").
- **L5 / L6 (merge cache-key projection, per-region path copies) — `experiment
  shows no gain — closing`.** The merge path is ~1–2% of runtime; projecting a
  16-byte key or borrowing paths cannot move wall. Filed as awareness only.
- **L7 (region-parallelism) — `won't fix (measured)`.** The fragmentation sweep
  shows the serial region loop costs ~2% even at 8000 regions; the per-region
  pipeline already saturates the cores. No starvation to recover; the
  oversubscription risk is not worth ~0 gain.
- **S1–S4 — closing.** All on the ~1–2% read/merge path; no measurable lever.

### 8.4 Where the real lever is (out of `--regions` scope)

BAQ (`probaln_glocal`, `ProbalnScratch::resize_for`) is ~70% of pileup self-time
— shared across **all** pileup, not `--regions`-specific, and already a
heavily-worked area ([ia/reviews/perf_baq_2026-05-12.md](../../../../ia/reviews/perf_baq_2026-05-12.md)).
Any pileup speedup lives there, not in the region machinery. Not pursued in this
`--regions`-scoped review.
