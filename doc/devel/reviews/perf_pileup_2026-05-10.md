# Performance Review: pileup
**Date:** 2026-05-10
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** `src/per_sample_caller/pileup/` — the sequential pileup walker (Stage 1)
**Verdict:** Run experiments
**Hot-path evidence:** saved criterion baselines for `pileup_walker_read_length` and `pileup_walker_multi_op` (`target-container/criterion/`); no flamegraph / heap profile / cache-miss profile available

---

## 1. Scope and constraints

- **What was reviewed:** the sequential walker module — public surface (`mod.rs`) plus all internal submodules (`walker.rs`, `active_set.rs`, `cigar_cursor.rs`, `decompose.rs`, `open_record.rs`, `slot_allocator.rs`, `errors.rs`).
- **Reviewed against:** `main` at `543d538` (crate state at session start; the `mod ` content of pileup has not changed since that commit).
- **Performance intent and targets:** WGS-scale inputs — multi-GB BAM/CRAM, ~30× coverage, chromosomes up to ~250 Mb. Walker is the per-sample caller's hot path; every base position passes through it once. Production calls `RefBaseFetcher` backed by `noodles_fasta::Repository`; the bench substitutes a constant fetcher.
- **Hot-path evidence (saved criterion baselines, span = 50_000, coverage = 30):**

  | Group | L=150 | L=500 | L=1500 | L=5000 |
  |---|---|---|---|---|
  | `pileup_walker_read_length` (pure-Match) | 415.6 ms | 446.4 ms | 440.9 ms | 459.2 ms |
  | `pileup_walker_multi_op` (periodic-insertion) | 468.3 ms | 494.0 ms | 491.0 ms | 531.4 ms |

  Roughly flat across a 33× sweep of L at fixed (span, coverage). Between-run noise floor previously measured at ±5–9 % in [pileup_clone_audit_2026-05-09.md](pileup_clone_audit_2026-05-09.md).

- **In-scope files:** [pileup/mod.rs](../../src/per_sample_caller/pileup/mod.rs), [walker.rs](../../src/per_sample_caller/pileup/walker.rs), [active_set.rs](../../src/per_sample_caller/pileup/active_set.rs), [cigar_cursor.rs](../../src/per_sample_caller/pileup/cigar_cursor.rs), [decompose.rs](../../src/per_sample_caller/pileup/decompose.rs), [open_record.rs](../../src/per_sample_caller/pileup/open_record.rs), [slot_allocator.rs](../../src/per_sample_caller/pileup/slot_allocator.rs), [errors.rs](../../src/per_sample_caller/pileup/errors.rs); plus `Cargo.toml` and [pileup_walker_scaling.rs](../../benches/pileup_walker_scaling.rs).
- **Deliberately out of scope:** `tests.rs` (correctness coverage, not perf); the production `RefBaseFetcher` impl (lives outside `pileup/`); other per-sample-caller modules.
- **Categories dispatched:**
  - `methodology` — always
  - `allocations` — heavy clone/Vec construction visible in inner loops
  - `data_layout` — wide hot-loop structs, possible cache pressure
  - `concurrency` — `Arc<str>` qname, `mpsc::sync_channel` to encoder
  - `hot_loops` — tight per-base / per-event loops in cursor and `apply_events_to_ref`
  - **Skipped:** `io_and_syscalls` — the walker itself does no I/O; `RefBaseFetcher` is a trait boundary

## 2. Verdict

**Run experiments.** One Hot-path finding is well-evidenced and is documentation-only (H1 below). Several Likely findings cluster at a small set of sites and have concrete measurement plans. **Apply the build-config changes (section 4) first** — they are cheap, low-risk, and may dwarf any code-level tweak. Code-level findings need a flamegraph and DHAT pass (section 3) to elevate above Likely; the prior clone audit's two reverted experiments failed precisely because the bench's between-run noise (±5–9 %) exceeded the candidate effects, and a profile would have decided each in one run.

## 3. Measurement plan

In order of leverage (each step unblocks the next):

1. **Add `[profile.release] debug = "line-tables-only"`** (or `[profile.bench] debug = true` inheriting from release) so subsequent profiles are readable.
2. **Flamegraph the walker.** `cargo flamegraph --bench pileup_walker_scaling` in the dev container. Top of the graph answers "is the time inside `events_overlapping_linear`'s op loop, `apply_events_to_ref`'s ref-cursor walk, the `BTreeMap` operations on `OpenPileupRecordTable`, or per-step `Vec`/`AHashMap` construction?". This is the prerequisite for elevating any code-level finding to Hot-path.
3. **Heap profile.** `dhat-rs` against `pileup_walker_multi_op/L=5000`. Confirms or refutes the per-step allocation hypotheses (L6–L12, L17–L18) in one run instead of round-tripping through criterion's noise floor.
4. **Cache-miss check.** `perf stat -e cache-references,cache-misses,L1-dcache-load-misses,LLC-load-misses` on the multi-op fixture. Pre-condition for prioritising data-layout findings (L13, L14) over allocation findings.
5. **Re-baseline criterion** after each build-config change in section 4, using `--save-baseline` / `--baseline`. One change per measurement.

Tools needed in the dev container, not yet installed: `cargo-flamegraph` (or `samply`), `dhat-rs` (crate dependency, not system tool), `perf`, optionally `cachegrind`. **Tell the user to add these to the dev container before step 2.**

## 4. Build / toolchain configuration

These are cheap and high-leverage. Apply each as its own commit so the bench delta is attributable.

### L1: `Cargo.toml` — **[Likely]** [build] No `[profile.release]` overrides

`lto`, `codegen-units`, `panic`, and `debug` are all at defaults. The cursor module already has a *manual* cross-function inlining workaround at [cigar_cursor.rs:213-215](../../src/per_sample_caller/pileup/cigar_cursor.rs#L213-L215) and [cigar_cursor.rs:343-347](../../src/per_sample_caller/pileup/cigar_cursor.rs#L343-L347) ("rustc apparently doesn't fully fuse the helper into this hot path") — that workaround is precisely the symptom LTO addresses automatically. Suggested:

```toml
[profile.release]
lto = "fat"
codegen-units = 1
panic = "abort"
debug = "line-tables-only"

[profile.bench]
inherits = "release"
debug = true
```

**Measurement plan:** save baseline → change `Cargo.toml` → re-bench → compare with `--baseline`. **Threshold:** any clean criterion "improved" verdict on `pileup_walker_read_length/L=150` (smallest reads = most setup overhead per useful work).
**Complexity cost:** six lines. Build time goes up; for a lab pipeline this is the right trade.

### L2: `Cargo.toml` — **[Likely]** [build] No allocator A/B

Walker per-step shape produces real allocation traffic (see L6–L12 below). `mimalloc` / `jemallocator` typically deliver size-class fast paths that fit this shape. Run the experiment, decide:

```toml
[features]
mimalloc = ["dep:mimalloc"]

[dependencies]
mimalloc = { version = "0.1", optional = true }
```
```rust
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
```

**Threshold:** clean criterion "improved" on `L=150` is the bar; below noise → revert.

### L3: [pileup_walker_scaling.rs:269-273](../../benches/pileup_walker_scaling.rs#L269-L273) — **[Likely]** [bench] Two `Criterion` configs disagree; `config()` is dead code

`group.sample_size(30)` / `group.measurement_time(10 s)` (set per-group) override the file-level `config()` at `:269-273` (`sample_size(10)`, `measurement_time(3 s)`). The unused config is a footgun for the next contributor. Either delete `config()` and use a plain `criterion_group!(benches, …)`, or align the two. Pick one.
**Cost:** negative — three lines deleted.

### L4: [pileup_walker_scaling.rs:184-228](../../benches/pileup_walker_scaling.rs#L184-L228) — **[Likely]** [bench] `iter_batched(LargeInput)` leaks fixture-build cost

`run_walker` allocates the channel, spawns the collector, and joins — all inside the timed region — and the input builder allocates ~3 MB of `PreparedRead` buffers per iteration at L=150. `BatchSize::PerIteration` with channel/thread setup in the setup closure measures only walker work. The clone audit already flagged this. **Threshold:** within-run CI shrinks below ~3 % *or* the L=150–L=5000 baselines move closer together (revealing real per-call walker scaling currently hidden behind setup noise).
**Cost:** ~20 bench lines; no walker change.

### L5: `Cargo.toml`, `benches/` — **[Likely]** [profile] No flamegraph / `samply` / `dhat` of the walker exists

Further code-level perf work is blind without one. The clone audit's two "Verdict — revert" sections failed because the bench cannot resolve the candidate effects through its noise floor; a flamegraph would have answered each in one run. Section 3 lists the exact commands. **Cost:** tooling install only.

### Speculative build/bench

- **S1: No `rust-toolchain.toml`.** Reproducibility insurance; cursor's manual-inlining workaround depends on rustc inliner decisions which drift between versions.
- **S2: No `target-cpu`.** The lab pipeline runs on dedicated hardware. `RUSTFLAGS="-C target-cpu=native"` (or `x86-64-v3`) would unlock AVX2 in the BAQ-min window and per-base scans. **Caveat:** binaries are not portable past the build CPU.
- **S3: One-fixture-per-bench is a microbenchmark.** Bench data fits in L2 (~50 kb × 30); production WGS spans are 5000× larger. Cache-warm-on-every-iter behaviour does not extrapolate. Suggested: a `bench_walker_span_streaming` group that varies span at fixed L=150, coverage=30.

## 5. Code-level findings

### Hot-path

#### H1: [pileup_walker_scaling.rs:1-12, 198-267](../../benches/pileup_walker_scaling.rs#L1-L12) — **[Hot-path]** [bench] Bench narrative is stale; the data validates the lazy-CIGAR refactor that already shipped

- **Confidence:** High
- **Hot-path evidence:** the bench output above + the bench's documented hypothesis ("flat or sub-linear scaling would mean the clone is not the bottleneck and the plan should be revisited") + source-tree state.
- **Pattern matched:** *No optimization without a benchmark — and no benchmark whose narrative contradicts its own data.*
- **Mechanism:** the lazy-CIGAR refactor in [pileup_lazy_cigar.md](../feature_implementation_plans/pileup_lazy_cigar.md) **already shipped** (commits `7088bac` / `401fe94` / `5886d0c` / `3a6dcc2`; closure recorded in [pileup_samtools_comparison_2026-05-07.md:172-194](pileup_samtools_comparison_2026-05-07.md)). `decompose` is now `#[cfg(test)]`, `ActiveRead.events` is gone, `ReadContribution.full_window_events` is gone, and `CigarCursor::events_at` / `events_overlapping` deliver the lazy lookup. Flat scaling across L is exactly what the cursor was supposed to deliver — it **validates** the refactor; it does not falsify it. The bench's prologue still describes the pre-refactor implementation. A future reader running the bench, seeing the flat data, and reading the prologue will conclude "the refactor is unnecessary" and either revert the cursor or repeat this investigation.
- **Measurement plan:** none required to confirm the staleness; the source-tree state plus the closure note in the samtools comparison are conclusive. Three documentation actions:
  1. Rewrite the file-level docstring at [pileup_walker_scaling.rs:1-12](../../benches/pileup_walker_scaling.rs#L1-L12) as a regression guard ("a >15 % slope between L=150 and L=5000 means the cursor's per-step cost has started scaling with L again — investigate").
  2. Update the function-level comments at `:50-52`, `:203-207` (remove "per-active-read clone cost", "clones a larger event vector at each walker step", "every walker step the read is alive sees a Match event").
  3. Move [pileup_lazy_cigar.md](../feature_implementation_plans/pileup_lazy_cigar.md) to `ia/reports/implementations/` (or add a status header) so it does not read as forward-looking.
- **Complexity cost:** documentation-only; zero code changes.
- **Suggested fix:** see methodology category file in `tmp/perf_review_2026-05-10_pileup/methodology.md` for verbatim replacement text.

### Likely

> Per-step / per-fold allocation cluster — file together; the flamegraph and DHAT pass in section 3 will tell which of these dominate. Each individually is below the bench noise floor; collectively they may not be.

#### L6: [walker.rs:231](../../src/per_sample_caller/pileup/walker.rs#L231) — **[Likely]** Per-step `Vec<ReadContribution>` reallocated every walker_pos
`Vec::new()` inside `process_position`, called ~50 000× per bench run (~3 G× per WGS run). Hoist onto `WalkerState` and `clear()` between steps; truncate-at-`column_depth_cap` already keeps the high-water short. **Measurement:** DHAT against `pileup_walker_multi_op/L=5000`. **Cost:** one extra `WalkerState` field.

#### L7: [walker.rs:432](../../src/per_sample_caller/pileup/walker.rs#L432) — **[Likely]** `resolve_mate_overlap_at_pos` rebuilds an `AHashMap` and two `Vec`s per walker step
Three fresh containers (`AHashMap<SlotId, Vec<usize>>`, `to_remove: Vec<usize>`, `bq_updates: Vec<…>`) on every call, even on the no-overlap common case. Either hoist + `clear()`, or add an early-exit fast-path. With chain slots bounded at `MAX_ACTIVE_SLOTS = 4096` a fixed-size `Vec<SmallVec<[usize; 2]>>` indexed by slot id avoids the map entirely. **Measurement:** DHAT; expect `AHashMap::new` + `Vec::new` sites to drop out.

#### L8: [open_record.rs:504](../../src/per_sample_caller/pileup/open_record.rs#L504) — **[Likely]** Per-step `Vec<u32>` for `affected` keys
Steady-state size 1–2. Cheapest experiment: `SmallVec<[u32; 4]>` local. **Cost:** one type change, no API impact.

#### L9: [cigar_cursor.rs:217](../../src/per_sample_caller/pileup/cigar_cursor.rs#L217), [:313](../../src/per_sample_caller/pileup/cigar_cursor.rs#L313), [:464](../../src/per_sample_caller/pileup/cigar_cursor.rs#L464), [:547](../../src/per_sample_caller/pileup/cigar_cursor.rs#L547) — **[Likely]** `events_at` / `events_overlapping` return owned `Vec<ReadEvent>` per call
The lazy-CIGAR plan in [pileup_lazy_cigar.md](../feature_implementation_plans/pileup_lazy_cigar.md) explicitly specified `SmallVec<[ReadEvent; 2]>` here ("a position can carry at most one Match plus one indel anchored there"); the shipped code uses plain `Vec`. Plan intent did not fully land. ~1.5 M `events_at` calls per bench run, similar for `events_overlapping`. **Suggested:** `SmallVec<[ReadEvent; 2]>` for `events_at`, `SmallVec<[ReadEvent; 4]>` for `events_overlapping`. Slice deref means call sites compile unchanged.

#### L10: [open_record.rs:362](../../src/per_sample_caller/pileup/open_record.rs#L362) — **[Likely]** `apply_events_to_ref` returns a fresh `Vec<u8>` per (record, contributor) per fold
Pre-sized (good) but constructed per call, then often dropped after the equality check at [open_record.rs:467](../../src/per_sample_caller/pileup/open_record.rs#L467) finds the existing REF/SNP allele. Pure waste in the steady-state SNP/REF column. **Suggested:** split into `apply_events_to_ref_into(&mut Vec<u8>, …)` plus a borrowed `find_allele_index(&[OpenAllele], &[u8]) -> Option<usize>`; only own the bytes when no match.

#### L11: [open_record.rs:362-442](../../src/per_sample_caller/pileup/open_record.rs#L362-L442) — **[Likely]** `apply_events_to_ref` byte-by-byte gap loop where `extend_from_slice` would memcpy
Distinct from L10 — same site, different angle. The inner `while ref_cursor < offset { … out.push(ref_seq[i]) … }` advances one byte at a time with two predicates per iteration. Once `ref_cursor >= consumed_until` the predicate is constant for the run; split into "skip until `consumed_until`" + `extend_from_slice(&ref_seq[start..end])` and the run becomes a single memcpy. Same shape applies to the tail loop. **Cost:** range bookkeeping; no `unsafe`. Verify with `cargo asm` that the emitted code calls `__memcpy`.

#### L12: [open_record.rs:142](../../src/per_sample_caller/pileup/open_record.rs#L142) — **[Likely]** `drain_aged` allocates two `Vec`s per walker step even when no record is aging
The empty-aged case is the steady state. One-line fast-path on `BTreeMap::iter().next().is_none_or(|(_, rec)| rec.footprint_end_exclusive() > walker_pos)` avoids the allocations.

#### L13: [active_set.rs:18](../../src/per_sample_caller/pileup/active_set.rs#L18) — **[Likely]** `ActiveRead` ~160 bytes; per-position scan touches ~28 bytes of it
The loop in `walker::process_position` reads `read_id`, `chain_slot_id`, `mq_log_err`, `is_reverse_strand`, `alignment_start`, `is_first_mate` plus the `cursor` — ~28 bytes scattered across all three cache lines occupied by each entry. `Vec<ActiveRead>` with entries inline → two extra cache lines fetched per iteration that the per-position scan never reads. DoD projection: parallel `Vec<ActiveScanRow>` (~24 bytes/entry) maintained alongside `reads`; dereference the full `ActiveRead` only inside the cursor call. **Cost:** second Vec to keep in sync (admit / `swap_remove`). **Threshold:** `perf stat -e cache-references,cache-misses,L1-dcache-load-misses` should drop measurably; criterion ≥ 5 % on multi-op L=5000 is the merge bar.

#### L14: [cigar_cursor.rs:113](../../src/per_sample_caller/pileup/cigar_cursor.rs#L113) — **[Likely]** `OpOffset` table and `CigarOp` cigar live in separate Vecs but are read in lockstep
Every cursor op-visit reads `self.offsets[i]` *and* `read.cigar[i]`. Fusing into `Vec<OpRow { ref_pos, read_pos, op }>` (16 bytes; 4/cache line) halves inner-loop fetch count. Cursor already owns the offsets table. **Cost:** small duplication of `read.cigar` data; cleaner than restructuring `PreparedRead`.

#### L15: [slot_allocator.rs:347](../../src/per_sample_caller/pileup/slot_allocator.rs#L347) — **[Likely]** `active_count()` linear-scans 4 KB on every admit
*Three sub-agents independently flagged this.* Each `set_refcount` call walks the full `slot_refcount: [u8; 4096]` to compute the high-water mark. Maintain `active_slots: u32` incrementally on 0 ↔ non-zero transitions. **Cost:** one field, one invariant; trivial. Below per-base hot path (admits are per-read), but the scan is gratuitous.

#### L16: [cigar_cursor.rs:240-260](../../src/per_sample_caller/pileup/cigar_cursor.rs#L240-L260), [:371-388](../../src/per_sample_caller/pileup/cigar_cursor.rs#L371-L388), [:485-494](../../src/per_sample_caller/pileup/cigar_cursor.rs#L485-L494), [:580-589](../../src/per_sample_caller/pileup/cigar_cursor.rs#L580-L589) — **[Likely]** Per-base Match-emit indexes `seq[i]` / `bq_baq[i]` inside `for k in 0..span`
`read_off = off.read_pos + (ref_pos - op_lo)` is monotonic and within bounds, but the compiler cannot prove it without a hoisted `assert!` or a slice + `iter().zip()` rewrite. Two separately-bounds-checked indexings per iteration also block straight-line autovectorization. **Suggested:** lift to `let bases = &read.seq[lo..hi]; let quals = &read.bq_baq[lo..hi]; for (k, (&b, &q)) in bases.iter().zip(quals.iter()).enumerate() { … }`. **Verify:** `cargo asm` shows the per-iteration `cmp+ja` is gone after.

#### L17: [slot_allocator.rs:350](../../src/per_sample_caller/pileup/slot_allocator.rs#L350) — **[Likely]** `evict_stale_pending` clones `Arc<str>` per stale entry
Collect-then-remove pattern forces an `Arc::clone` (atomic increment) per stale entry. `HashMap::retain` (or `extract_if`) avoids the intermediate. **Dormant on the existing solo-read benches** (`has_mate: false` everywhere); fires on paired-end inputs. **Measurement plan:** add a paired-end fixture before merging.

#### L18: [active_set.rs:174](../../src/per_sample_caller/pileup/active_set.rs#L174) — **[Likely]** `expire_passed` clones `Arc<str>` per expiring read
`let qname = self.reads[i].read.qname.clone();` exists only to release the partner ref by name. Borrow-split or a consolidated `ActiveSet::expire_one` removes the clone. Same dormancy as L17.

#### L19: [open_record.rs:681-699](../../src/per_sample_caller/pileup/open_record.rs#L681-L699) — **[Likely]** Inner-fold helpers not `#[inline]`
`add_contribution`, `subtract_contribution`, `phred_to_ln_perr`, `ln_bq_for_read`, `insert_sorted_unique`, `event_kind_rank`, `zero_event_bq` — all small, all called inside the per-(record, contributor) fold. Fine today (binary build, intra-crate); without LTO (L1 not yet applied) and if the walker is later linked as a library, calls are not inlined. **Cost:** seven attribute additions, no semantic change. **Conditional on:** L1 not landing or library extraction happening.

### Speculative

- **S4:** [walker.rs:312](../../src/per_sample_caller/pileup/walker.rs#L312), [:373](../../src/per_sample_caller/pileup/walker.rs#L373) — `aged.into_iter().map(finalise).collect()` then iterate-to-send. Per-emission, not per-step, so likely below noise; hoist `emit_buf` if a future profile shows it.
- **S5:** [slot_allocator.rs:353](../../src/per_sample_caller/pileup/slot_allocator.rs#L353) — `evict_stale_pending`'s outer `Vec<(Arc, SlotId)>`. Dormant on solo-read benches.
- **S6:** [open_record.rs:504-535](../../src/per_sample_caller/pileup/open_record.rs#L504-L535) — `affected.contains(&key)` inside the step-3 loop is small-N quadratic. Drop the guard; sort + dedup once at end.
- **S7:** [walker.rs:344-352](../../src/per_sample_caller/pileup/walker.rs#L344-L352), [active_set.rs:101-108](../../src/per_sample_caller/pileup/active_set.rs#L101-L108), [slot_allocator.rs:215-222](../../src/per_sample_caller/pileup/slot_allocator.rs#L215-L222) — `format!`-built `WalkerError::Internal` arms inline in fast-path functions. Hide behind `#[cold] #[inline(never)] fn …` constructors. Layout-only win; no throughput change expected on the happy path.
- **S8:** [mod.rs:128](../../src/per_sample_caller/pileup/mod.rs#L128) — `PreparedRead` carries 3× `Vec` headers + `Arc<str>` inline (~120 bytes); `Box<[T]>` / `Arc<[T]>` for the buffers shrinks it. Defer until after L13 lands.
- **S9:** [mod.rs:54](../../src/per_sample_caller/pileup/mod.rs#L54) — `DEFAULT_OUTPUT_CHANNEL_CAPACITY = 64` has no measurement backing; the bench's trivial collector cannot exercise back-pressure shape. Expose on `WalkerConfig` and revisit once Stage 2 is wired.

### Note (no action)

- `WalkerError` `format!` / `to_string()` on error construction — cold path; closed in the prior clone audit.
- `ReadEvent` enum width (~40 bytes due to inline `Vec<u8>` for `Insertion::seq`) — fine for typical Match-dominant workload.
- `AHashMap` already used everywhere a hash map is needed — the *Hot loops & codegen* checklist's FxHash recommendation does not apply.
- [cigar_cursor.rs:69](../../src/per_sample_caller/pileup/cigar_cursor.rs#L69) `base_in_adaptor` already `#[inline(always)]`.
- `AlleleObservation`, `PileupRecord`, `FiveScalars` widths checked, no change recommended.
- Record-at-a-time `tx.send` shape — fine today; flag only if Stage 2 wants batched dispatch.

## 6. Out-of-scope observations

- **`RefBaseFetcher::fetch` returns `Vec<u8>` per call** ([mod.rs:269-275](../../src/per_sample_caller/pileup/mod.rs#L269-L275)). The trait shape forces a per-call allocation regardless of the production impl over `noodles_fasta::Repository`. Walker calls `fetch` per `widen` and per `open_new` ([open_record.rs:260](../../src/per_sample_caller/pileup/open_record.rs#L260), [:322](../../src/per_sample_caller/pileup/open_record.rs#L322)) — once per affected record per step. Suggested follow-up: add `fn fetch_into(&self, …, out: &mut Vec<u8>)` so the walker can reuse a buffer. Bench currently hides this with `ConstFasta`; production walks may have a different allocation profile.
- **[open_record.rs:209-215](../../src/per_sample_caller/pileup/open_record.rs#L209-L215) `find_overlapping`** does `BTreeMap::range(lo..=event_start).rev()` per event per step. `BTreeMap` has known cache-unfriendliness vs. a small sorted `Vec` at typical record counts (1–4); investigate after the higher-leverage findings land.
- **[open_record.rs:466-473](../../src/per_sample_caller/pileup/open_record.rs#L466-L473) `find_or_create_allele_index`** does linear `Vec<u8>::eq` over every allele candidate. Small-N today; matters at over-cap regions.

## 7. What's already good

- **Per-call `Vec::with_capacity` is used where the size is known** (e.g. [pileup_walker_scaling.rs:58](../../benches/pileup_walker_scaling.rs#L58), [cigar_cursor.rs](../../src/per_sample_caller/pileup/cigar_cursor.rs) cursor construction) — preserves the pattern this review's allocation findings target at the next layer.
- **`AHashMap` is used consistently** throughout the module where the default `HashMap` would be HashDoS-overhead — exactly the recommendation in the *hot loops & codegen* checklist, applied as a project default.
- **The `CigarCursor`'s linear/binary-search dispatch on `BINARY_SEARCH_OP_THRESHOLD`** ([cigar_cursor.rs:90](../../src/per_sample_caller/pileup/cigar_cursor.rs#L90)) and the "early-break" optimization in op-loops are textbook hot-loop work — both already shipped (referenced in the bench's docstring as already-done).

---

## Per-category audit trail

Sub-agent output is preserved at `tmp/perf_review_2026-05-10_pileup/`:

- [methodology.md](../../tmp/perf_review_2026-05-10_pileup/methodology.md) — H1 plus bench/build findings
- [allocations.md](../../tmp/perf_review_2026-05-10_pileup/allocations.md) — L6–L12 plus speculative
- [data_layout.md](../../tmp/perf_review_2026-05-10_pileup/data_layout.md) — L13, L14
- [concurrency.md](../../tmp/perf_review_2026-05-10_pileup/concurrency.md) — L17, L18, S9
- [hot_loops.md](../../tmp/perf_review_2026-05-10_pileup/hot_loops.md) — L11, L15, L16, L19 plus speculative

## Tools requested for the dev container

To run the section 3 measurement plan: `cargo-flamegraph` (or `samply`), `perf` with hardware-counter access (or document the equivalent for the dev container's kernel), optionally `valgrind`/`cachegrind`. `dhat-rs` is a crate, no system install needed. **The user offered to install these on request.**
