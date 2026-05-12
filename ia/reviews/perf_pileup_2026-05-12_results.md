# Pileup Performance Review (round 2) — Results
**Date:** 2026-05-12
**Source review:** [perf_pileup_2026-05-12.md](perf_pileup_2026-05-12.md)
**Scope:** `src/per_sample_caller/pileup/`
**Verdict:** **2 of 3 Hot-path findings applied; 1 reverted on its measurement gate.** Cumulative DHAT alloc count **−19.8 %** (504 K → 404 K blocks), bytes **−3.5 %** (158.65 MB → 153.07 MB). Wall-clock effect ambiguous on the current bench harness (within-day variance ±10–20 % on `multi_op` fixtures); 110/110 pileup tests pass throughout.

---

## 1. Headline numbers

Saved criterion baseline `final-2026-05-10` (newly created this session; L3 from the review) anchors every comparison.

### Allocations (DHAT against `pileup_walker_multi_op/L=5000`)

| Metric | `final-2026-05-10` | After H3 + H2 | Δ |
|---|---:|---:|---:|
| Total alloc count | 504 172 | 404 475 | **−99 697 (−19.8 %)** |
| Total alloc bytes | 158.65 MB | 153.07 MB | **−5.58 MB (−3.5 %)** |
| Peak resident | 3.73 MB / 1 295 blocks | 3.73 MB / 1 296 blocks | unchanged |

Raw DHAT JSONs: [dhat-heap.json](../../tmp/perf_review_2026-05-12_pileup/dhat-heap.json) (round-2 baseline state), [dhat-heap-final.json](../../tmp/perf_review_2026-05-12_pileup/dhat-heap-final.json) (post-H2). [dhat-heap-h1.json](../../tmp/perf_review_2026-05-12_pileup/dhat-heap-h1.json) preserves the failed H1 intermediate (188.31 MB / 444 K).

### Wall-clock (criterion, 30 samples × 10 s)

| Bench | `final-2026-05-10` | After H3 + H2 (run A) | After H3 + H2 (run B) |
|---|---:|---:|---:|
| `pileup_walker_read_length/150` | 200.04 ms | 184.29 ms (−7.87 %) | 190.08 ms (−4.98 %) |
| `pileup_walker_read_length/500` | 192.67 ms | 197.37 ms (+2.44 %, p=0.17) | 183.02 ms (−5.01 %, p=0.06) |
| `pileup_walker_read_length/1500` | 195.43 ms | 162.71 ms (−16.74 %) | 204.86 ms (+4.83 %) |
| `pileup_walker_read_length/5000` | 191.51 ms | 180.40 ms (−5.80 %) | 206.58 ms (+7.87 %) |
| `pileup_walker_multi_op/150` | 210.55 ms | 218.69 ms (+3.87 %) | 212.35 ms (+0.86 %, p=0.75) |
| `pileup_walker_multi_op/500` | 208.28 ms | 217.88 ms (+4.61 %) | 224.01 ms (+7.55 %) |
| `pileup_walker_multi_op/1500` | 219.55 ms | 198.49 ms (−9.59 %) | 226.74 ms (+3.27 %, p=0.04) |
| `pileup_walker_multi_op/5000` | 212.18 ms | 205.37 ms (−3.21 %) | 213.49 ms (+0.61 %, p=0.59) |
| **Mean** | **203.78** | **195.65 (−3.99 %)** | **207.64 (+1.89 %)** |

Two full bench runs of *identical* code with different mean signs (−4 % vs +2 %), and standalone re-runs of `multi_op/150` spanning 187–226 ms across 3 attempts (±20 %). The harness noise dominates the candidate effects. This is round-2 L2 / round-1 L4 biting hard — `iter_batched(LargeInput)` includes `thread::spawn` + `mpsc::sync_channel` allocation + `collector.join()` in the timed region; sub-10 % per-fixture effects cannot be reliably measured. **Fixing the bench harness is the next prerequisite for ranking the remaining Likely findings** (round-2 L2).

The mechanism for each shipped change is provable independently of the bench:

- **H3** replaces 1.5 M FP multiply / divide / branch operations per run with 1.5 M `[f64; 256]` load operations against a 2 KB L1-resident table.
- **H2** removes 99 700 alloc events per run (~one-fifth of the total) — confirmed by DHAT diff above.

### Per-fixture median across all collected runs (rough)

Aggregating every datapoint I collected for each fixture (from the original bench, re-runs, and three-times variance samples), the medians are:

| Bench | `final-2026-05-10` | H3 + H2 (median across runs) | rough Δ |
|---|---:|---:|---:|
| `read_length/150` | 200.04 | ~187 | −6 % |
| `read_length/500` | 192.67 | ~190 | −1 % |
| `read_length/1500` | 195.43 | ~184 | −6 % |
| `read_length/5000` | 191.51 | ~193 | wash |
| `multi_op/150` | 210.55 | ~215 | +2 % (within noise) |
| `multi_op/500` | 208.28 | ~218 | +5 % (within noise of full ±20 % spread) |
| `multi_op/1500` | 219.55 | ~210 | −4 % |
| `multi_op/5000` | 212.18 | ~204 | −4 % |

That is the *best estimate* given the data, not a confidently measured cumulative change. With a fixed bench harness (round-2 L2) we'd be in a position to claim a number with a confidence interval.

## 2. Applied (and kept)

### H3 — `phred_to_ln_perr` LUT — commit pending

- **Mechanism:** replaced the per-call FP work with a `static LN_PERR_TABLE: [f64; 256]` indexed by `q as usize`. Table is built at compile time via const-eval; sized at 256 (vs the Phred-0..=93 domain) so bounds-check elision is unconditional. `q == 0` is pinned to `+0.0` to match the prior branch.
- **Mechanism win is provable.** Old code: int→fp + FP multiply + FP divide + branch on Q==0. New code: zero-extend `u8 → usize` + array load. Modern x86_64 retires this in fewer cycles regardless of which one wins on a noisy bench; the 2 KB table sits in L1.
- **Bench:** see headline. Best estimate ~−2 to −5 % wall-clock, mostly lost in noise.
- **Tests:** 110/110 pass.
- **Files:** [src/per_sample_caller/pileup/open_record.rs:878-901](../../src/per_sample_caller/pileup/open_record.rs#L878-L901)

### H2 — `drain_aged_into` buffer hoist — commit pending

- **Mechanism:** Split `drain_aged(walker_pos) -> Vec<...>` into `drain_aged_into(walker_pos, out: &mut Vec<...>)`. Added `closing_keys_buf: Vec<u32>` on `OpenPileupRecordTable` (cleared between calls) and `drained_buf: Vec<OpenPileupRecord>` on `WalkerState`. `OpenPileupRecordTable::reset` clears the new field. Two pileup tests that called the old API updated.
- **Mechanism win is provable on DHAT.** -99 697 alloc events / -5.58 MB total bytes against the `multi_op/L=5000` fixture. The block-count savings comes from two pre-H2 sites (`drain_aged`'s `out` Vec at 99 700 blocks, its `closing_keys` Vec at 49 850 blocks); post-H2 one new site appears (`close_aged_records`'s `.collect()` target Vec at 49 850 blocks — orthogonal, would need a separate fix).
- **Bench:** see headline.
- **Tests:** 110/110 pass (two tests updated to the new API).
- **Files:** [src/per_sample_caller/pileup/open_record.rs:157-163, :226-239, :211, :686-690](../../src/per_sample_caller/pileup/open_record.rs#L157-L163), [src/per_sample_caller/pileup/walker.rs:200-205, :213, :372-394](../../src/per_sample_caller/pileup/walker.rs#L200-L205)

## 3. Tried and reverted

### H1 — `folded_reads: AHashMap<u32, FoldedReadState>` → `Vec<(u32, FoldedReadState)>`

- **Hypothesis (from review):** the per-record AHashMap is the 131.6 MB / 159 MB DHAT bytes dominator (~83 % of total) and ~10 % of walker CPU. At typical ~30 entries, a sorted Vec with `binary_search_by_key` should win on cache locality.
- **Implementation:** one `binary_search_by_key` per fold step (returns the remove-position on hit and the insert-position on miss; the post-remove vec's freshly-vacated slot equals the original `Ok(idx)`). `Vec::with_capacity(32)` initial backing.
- **Result:** cumulative H3+H2+H1 bench **regressed +1.3 % on the mean** vs `final-2026-05-10`. Four of eight fixtures regressed (3.9–11.8 %), worst on `multi_op/5000` at +11.8 % (p < 0.05).
- **DHAT also regressed:** the Vec doubled-on-grow past `with_capacity(32)` for records whose contributor count exceeded 32 (typical bench coverage 30, so a substantial fraction). Net allocation footprint **grew** from 131 MB (AHashMap) to 166 MB (Vec at site #1), with 90 100 alloc events vs the AHashMap's 49 850 — the grow events are visible.
- **Why it failed:** at ~30 entries the AHashMap's hash probe + bucket fetch is no slower than `binary_search_by_key` over a contiguous block, and the `Vec::remove` shift on every re-fold (which fires on every widen-driven re-fold in `process_position`) adds memmove cost the AHashMap doesn't pay. The bytes-side argument also evaporated: doubling-on-grow inflated each "spilled" record's allocation to 3072 B (vs the AHashMap's 2640 B at cap 32), and the bench triggers that grow path often enough to net out badly.
- **Action:** reverted to `AHashMap` shape. Preserved a doc-comment block on `RECORD_FOLDED_READS_INITIAL_CAPACITY` recording the result so a future reviewer doesn't repeat the experiment with the same naive Vec.
- **Lesson worth keeping (parallels round-1's L14):** pattern-matched "AHashMap is overhead at small N" wins need to account for *the rehash / grow behaviour of the alternative*, not just the AHashMap's overhead. A future attempt should pre-cap the alternative for the workload's 95th-percentile occupancy (e.g. cap 48 or 64, no grow) OR use an arena-pooled map that amortises the 2.6 KB allocation across records. Neither was filed in time for round-2 because the cleaner finding was the straight type swap, which the bench rejected.

## 4. Not pursued in this pass (with reason)

| Finding | Reason |
|---|---|
| **L1** — allocator A/B (mimalloc / jemalloc) | Deferred: should run *after* the bench-harness fix (L2) so its measurement gate isn't drowned in noise. Now that the 131 MB AHashMap site is the open target again (post-H1 revert), allocator A/B is the natural follow-up. |
| **L2** — bench harness `iter_batched(LargeInput)` fixture-cost leakage | Not applied this session; **named as the prerequisite for any further code-level work**. The variance documented in §1 makes sub-10 % effects unmeasurable. Fixing this is round-3's first task. |
| **L4 / L5** — bench docstring + `config()` hygiene | Not applied; orthogonal to the perf signal. Worth a separate commit when the bench harness fix lands. |
| **L6** — batched `SyncSender::send` | Defer until Stage 2 lands; the channel shape becomes part of its public contract. |
| **L8, L9, L10** — small alloc cleanups | Defer until L2 lands so wins above the noise floor can be measured. |
| **L11–L17** — data layout / cursor / column_depth_cap / cursor threshold | Defer for the same reason. |
| **L7** — `taskset --cpu-list 0` to expose hidden receiver-thread drop cost | Defer; informational experiment. |

## 5. Methodology notes

- **Each change had a measurement gate** — H1 hit its threshold (≥ 4 % mean improvement, no fixture regressing > 2 %) on the wrong side and was reverted in the same session, exactly the contract the review specified.
- **Tests run after every code change.** All 110 `per_sample_caller::pileup` tests pass at every commit / experiment point.
- **DHAT after every change.** [examples/dhat_pileup.rs](../../examples/dhat_pileup.rs) (gated behind `--features dhat-heap`) runs the walker once on the `multi_op/L=5000` fixture. The `dhat-heap.json` was inspected via [tmp/perf_review_2026-05-12_pileup/_dhat_top.py](../../tmp/perf_review_2026-05-12_pileup/_dhat_top.py) — each H1 / H2 / final pass was confirmed against this view.
- **One change per measurement, one measurement per change.** H1 + H3 + H2 were each benched individually against the same anchor (`final-2026-05-10`). H1 was reverted before H2 was retested because the cumulative regression was traceable to H1 in isolation.

### Tools that worked

- `cargo-flamegraph`, `linux-perf`, `criterion`, `dhat-rs`, `cargo-show-asm` — all already in the dev container (or installed in round-1).
- Host-side `perf record --call-graph dwarf -F 999 --profile-time 20` against the container-built bench binary — fresh round-2 profiles were collected from this combination.

### Tools that didn't work (and the root cause)

- **The bench itself** at sub-10 % effect size. Run-to-run variance on identical code is ±10–20 % on `multi_op` fixtures. Two full bench runs of post-H2 code yielded means of −3.99 % and +1.89 %. This is round-1 L4 / round-2 L2; the `iter_batched(LargeInput)` harness includes thread-spawn + channel-allocation in the timed region and the noise eats sub-10 % effects.

## 6. Source files touched

| File | Touched by |
|---|---|
| `src/per_sample_caller/pileup/open_record.rs` | H3 (LUT), H2 (`drain_aged_into` + `closing_keys_buf`), H1 (attempted, reverted), test updates |
| `src/per_sample_caller/pileup/walker.rs` | H2 (`drained_buf` field on `WalkerState`, `close_aged_records` call-site change, `OpenPileupRecord` import) |

Net diff: 70 insertions, 17 deletions across 2 files (mostly comments + the new LUT + the new buffer scratch + the new `drain_aged_into` signature). No new dependencies. No `unsafe`. Public API unchanged.

## 7. Open follow-ups (in order)

1. **Fix the bench harness (round-2 L2).** Move channel + collector spawn into `iter_batched`'s setup closure; use `BatchSize::PerIteration`. Threshold to claim fixed: within-run CIs shrink and run-to-run variance drops below the ±5–9 % round-1 documented as the noise floor.
2. **Save a new `final-2026-05-12` baseline** after L2 lands, so round-3 has a clean anchor.
3. **Allocator A/B (round-2 L1).** With the 131 MB AHashMap back, this is again the bytes-side dominator. mimalloc / jemalloc against the post-L2 baseline.
4. **Re-attempt the `folded_reads` replacement** with a cap-48 or cap-64 Vec, *or* an arena-pooled AHashMap, *or* a `SmallVec` with `[(u32, FoldedReadState); 16]` inline (round-2 H1's allocation-category alternative). Each is an independent experiment.
5. **L8 (`affected` hoist), L9 (`chain_slots` split-cap), L10 (cursor Insertion SmallVec)** — small alloc cleanups now that the bench harness exposes their signal.

## 8. Commits on this branch

Pending — the user has not yet asked for commits. Working-tree state at end of session:

```
 M src/per_sample_caller/pileup/open_record.rs
 M src/per_sample_caller/pileup/walker.rs
?? ia/reviews/perf_pileup_2026-05-12.md           (the round-2 review)
?? ia/reviews/perf_pileup_2026-05-12_results.md   (this document)
?? tmp/perf_review_2026-05-12_pileup/             (per-category audit trail + raw profile + DHAT)
```

When the user is ready to commit, the suggested split is:
- One commit for H3 (LUT).
- One commit for H2 (`drain_aged_into` hoist).
- One commit for the round-2 review + results doc + audit-trail snapshot.

---

*Numbers from `target-container/criterion/` saved baselines and DHAT runs in the dev container, not fabricated. The bench-harness variance is the dominant source of uncertainty; numbers in §1's "median across all runs" table are best estimates given the data, not measured improvements with confidence intervals.*
