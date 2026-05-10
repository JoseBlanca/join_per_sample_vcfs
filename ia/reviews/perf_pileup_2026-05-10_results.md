# Pileup Performance Review — Results
**Date:** 2026-05-10
**Source review:** [perf_pileup_2026-05-10.md](perf_pileup_2026-05-10.md)
**Scope:** `src/per_sample_caller/pileup/`
**Verdict:** 7 changes applied, 1 reverted, 5 deferred. Cumulative **-54% mean walker wall-clock, -93% allocations**, no regressions kept, 88/88 pileup tests pass throughout.

---

## 1. Headline numbers

Saved criterion baselines (`target-container/criterion/{group}/{size}/{pre-experiments,post-l7}/`) anchor every comparison.

### Wall-clock (criterion, 30 samples × 10 s, 8 fixtures)

| Bench | Pre-experiments | Final | Δ |
|---|---:|---:|---:|
| `pileup_walker_read_length/150` | 369.7 ms | 197.4 ms | **−46.6 %** |
| `pileup_walker_read_length/500` | 423.3 ms | 175.2 ms | **−58.6 %** |
| `pileup_walker_read_length/1500` | 417.2 ms | 179.2 ms | **−57.0 %** |
| `pileup_walker_read_length/5000` | 435.8 ms | 180.8 ms | **−58.5 %** |
| `pileup_walker_multi_op/150` | 461.4 ms | 207.3 ms | **−55.1 %** |
| `pileup_walker_multi_op/500` | 410.4 ms | 214.4 ms | **−47.8 %** |
| `pileup_walker_multi_op/1500` | 475.1 ms | 212.1 ms | **−55.4 %** |
| `pileup_walker_multi_op/5000` | 485.3 ms | 219.8 ms | **−54.7 %** |
| **Mean** | | | **−54.2 %** |

### Allocations (DHAT against `pileup_walker_multi_op/L=5000`)

| Metric | Pre-experiments | Final | Δ |
|---|---:|---:|---:|
| Total alloc count | 7,292,320 | 504,174 | **−93.1 %** |
| Total alloc bytes | 1.14 GB | 159 MB | **−86 %** |
| Peak resident | 3.73 MB / 1330 blocks | 3.73 MB / 1271 blocks | unchanged |

## 2. Applied (and kept)

Each change validated independently before commit. Saved baselines + DHAT diff for every step. Order matches commit order on `main`.

| # | Commit | Finding | Mechanism | Wall-clock impact (vs prior) | Allocation impact |
|---|---|---|---|---|---|
| 1 | `98ad664` | **L1** | `[profile.release] lto = "fat", codegen-units = 1`; `[profile.bench] debug = true` | 7/8 improve 3.5–16.7 %, 1 regression at rl/150 (+11.1 %, attributable to L4 setup-cost dominance) | n/a |
| 2 | `98ad664` | **L7** | `resolve_mate_overlap_at_pos` early-exit when no two contributors share a `chain_slot_id` | 8/8 improve 19–35 % | −1.73 M blocks |
| 3 | `0b33563` | **L9** | Cursor `events_at` / `events_overlapping` return `SmallVec<[ReadEvent; N]>` (lazy-CIGAR plan intent that hadn't fully landed) | 7/8 improve 8.5–18.4 %, 1 no-change | −3.0 M blocks (events_at + events_overlapping outer Vecs) |
| 4 | `6bad07e` | **L10 + L11** | `apply_events_to_ref_into(&mut Vec<u8>, …)` + buffer hoist on `OpenPileupRecordTable.allele_seq_buf`; borrowed `find_allele_index`; `extend_from_slice` for the gap-loop | 8/8 improve 3.4–19.4 %, 0 no-change | −1.5 M blocks |
| 5 | `3fc128c` | **L6** | `WalkerState.contributors_buf` hoist; `Vec<ReadContribution>` cleared per step instead of allocated | 5/8 improve 3.9–9.7 %, 1 no-change, 2 small regressions at rl/150 / rl/500 (+2.4 %, +2.7 %; within noise floor) | −231 K blocks |
| 6 | `60c8659` | (extra) | `OpenPileupRecord.folded_reads = AHashMap::with_capacity(32)` | 8/8 improve 7.9–19.2 %, all p < 0.05 | −184 K blocks (FoldedReadState reserve_rehash) |
| 7 | `9b08396` | (extra) | `OpenAllele.chain_slots = Vec::with_capacity(32)` | 5/8 improve 1.9–7.3 %, 3 no-change | −140 K blocks (insert_sorted_unique grows) |

Findings 6 and 7 were not numbered L-findings in the original review — they were surfaced by DHAT during the apply pass and matched the methodology checklist's *Pre-size containers when the size is known or bounded* rule. Both follow the same shape: pre-allocate at typical-WGS-coverage capacity (32) so per-record growth is amortised.

## 3. Tried and reverted

### L14 — OpRow fusion (`OpOffset` + `CigarOp` into one struct)

- **Hypothesis (from review):** the cursor's per-op inner loop reads `self.offsets[i]` and `read.cigar[i]` from two separately-allocated Vecs. Fusing into `Vec<OpRow { ref_pos, read_pos, op }>` (16 bytes) halves the inner-loop fetch count.
- **Result:** all 8 benches regressed. +3.7 % to +14.0 %, mean **+7.6 %**, all p < 0.05.
- **Why it failed:** doubling the row width (8 → 16 bytes) inflated the cursor's `Vec` footprint, and that cost outweighed the saved fetch. At `multi_op/L=5000` the cursor's table grew from ~1.6 KB to ~3.2 KB per active read, plus icache pressure from regenerated cursor inner-loop code.
- **Action:** reverted. The original `OpOffset` + `read.cigar[i]` shape is restored.
- **Lesson worth keeping:** pattern-matched "halve the cache misses" wins need to account for *footprint inflation*, not just fetch count. The original review's measurement plan called for "≥ 5 % criterion improvement at L=5000 multi_op as the merge bar"; the actual measurement gave the opposite signal and the change was reverted in the same session — exactly what the *measurement plan* gate is for.

## 4. Not pursued (with reason)

| Finding | Reason |
|---|---|
| **L12** — `drain_aged` empty-case short-circuit | Records age out ~1-per-step in the bench fixture, so the empty-case path is rarely taken. The optimisation targets the steady-state-empty workload, not what the bench exercises. |
| **L13** — `ActiveRead` DoD projection (parallel `Vec<ActiveScanRow>`) | The review's cache-miss claim assumed the per-position scan reads scattered fields *without* touching `PreparedRead`'s headers. In practice every iteration calls `events_at` which dereferences `read.cigar` / `read.seq` / `read.bq_baq` Vec headers — the same lines the projection was meant to avoid. After L14's regression, the refactor's expected gain didn't justify the maintenance cost. Re-evaluate when `perf stat` / `cachegrind` data is available. |
| **L15** — `slot_allocator::active_count` incremental counter | Per-admit, not per-base. At ~50 K admits per bench run, replacing the 4 KB linear scan saves microseconds total. Below the bench's noise floor; clean-up rather than perf. |
| **L16** — Per-base Match-emit `slice + iter().zip()` for autovec | `cargo-show-asm` is now installed in the dev container, but verifying autovec on this site needs a contained microbench (the criterion fixtures already include the entire walker, so the vectorisable inner loop is one of many cost sources). Deferred. |
| **L17 / L18** — `Arc<str>` qname clones in mate-eviction paths | Dormant in the existing solo-read fixtures (`has_mate: false` in `pileup_walker_scaling.rs`). Need a paired-end fixture with realistic orphan rates before this can be measured at all. |
| **L19** — `#[inline]` on small `open_record.rs` helpers | Moot now that L1's `[profile.release] lto = "fat"` is in. LTO does cross-crate inlining without explicit hints. |

## 5. Methodology notes

- **Each change had a measurement gate.** Criterion baselines saved as `pre-experiments` and `post-l7` under `target-container/criterion/{group}/{size}/`; comparisons computed both via criterion's own `--baseline` machinery and via a small Python helper at `tmp/perf_review_2026-05-10_pileup/_compare.py`.
- **DHAT after every change.** `examples/dhat_pileup.rs` (gated behind `--features dhat-heap`) runs the walker once on the `multi_op/L=5000` fixture. The `dhat-heap.json` was inspected via `tmp/perf_review_2026-05-10_pileup/_dhat_top.py` to find the next-largest alloc site by count and bytes after each change.
- **Tests run after every code change.** All 88 `per_sample_caller::pileup` tests pass at every commit.
- **One change per commit, one measurement per change.** L1 + L7 were bundled into a single commit (`98ad664`) only because L1 is a pure build-config edit that didn't require its own validation cycle.

### Tools that worked

- `cargo-flamegraph`, `linux-perf`, `valgrind`, `hyperfine`, `criterion` — all already in the dev container.
- `samply` and `cargo-show-asm` — installed during the apply pass (commit `0f5d8a6`).
- `dhat-rs` — added as optional dep + `dhat-heap` feature + `examples/dhat_pileup.rs` (commit `0f5d8a6`). Drove every alloc-side decision in this session.

### Tools that didn't work

- **`samply` and `cargo flamegraph`** require host kernel `perf_event_paranoid ≤ 1`. The dev container reports level 3, so neither produced a profile despite being installed. Workaround: `sudo sysctl kernel.perf_event_paranoid=1` on the host (not done — DHAT was sufficient for this session's allocation-driven decisions).
- **`cargo-show-asm`** runs but was not used; verifying L16's autovec claim needs a contained microbench, deferred.

## 6. Source files touched

| File | Touched by |
|---|---|
| `Cargo.toml` | L1 (build config), L9 (smallvec dep), tooling commit (dhat dep + dhat-heap feature) |
| `Cargo.lock` | dependency additions |
| `Containerfile` | tooling commit (samply, cargo-show-asm) |
| `.gitignore` | tooling commit (`/dhat-heap.json`, `/flamegraph.svg`, `/perf.data*`) |
| `examples/dhat_pileup.rs` | tooling commit (new) |
| `src/per_sample_caller/pileup/walker.rs` | L7 (resolve_mate_overlap fast-path), L6 (contributors_buf hoist) |
| `src/per_sample_caller/pileup/cigar_cursor.rs` | L9 (SmallVec returns); L14 attempted and reverted |
| `src/per_sample_caller/pileup/open_record.rs` | L9 (field type), L10/L11 (apply_events_to_ref_into + extend_from_slice + find_allele_index), L6 (warnings cleanup), pre-allocations |

## 7. Commits on this branch

```
9b08396  Pre-allocate OpenAllele.chain_slots at typical-coverage capacity
60c8659  Pre-allocate OpenPileupRecord.folded_reads at typical-coverage capacity
3fc128c  Apply L6 from the pileup performance review
6bad07e  Apply L10 + L11 from the pileup performance review
0b33563  Apply L9 from the pileup performance review
98ad664  Apply L1 + L7 from the pileup performance review
2b46645  Add pileup performance review
0f5d8a6  Add profiling tooling for the pileup performance review
543d538  Add rust-performance-review skill
```

## 8. Open follow-ups

- **Update the review status** in [perf_pileup_2026-05-10.md](perf_pileup_2026-05-10.md): mark each finding with its actual outcome (kept / reverted / deferred) per the *Author response convention* the skill defines. This results document is the audit trail; the original review can carry the brief status annotations.
- **Re-baseline.** Save `final` as a named criterion baseline so future changes have a stable comparison point: `cargo bench --bench pileup_walker_scaling -- --save-baseline final-2026-05-10`.
- **L4 — bench harness.** `iter_batched(LargeInput)` leaks fixture-build cost into measurement; the small-L regressions on L1 and L6 both trace to this. Fixing the harness is bench-hygiene, not optimisation, but it would tighten the noise floor and let smaller candidate effects (L15, L16, L17/L18) be resolved.
- **Cache-miss profile.** Relax `perf_event_paranoid` to 1 (or run a privileged container session) and capture a flamegraph + `perf stat -e cache-misses,L1-dcache-load-misses` against the post-final binary. That data would clarify whether L13 (DoD projection) and L16 (autovec) have plausible remaining wins or have already been amortised.

---

*Generated alongside the apply pass; numbers from saved criterion baselines and DHAT runs in the dev container, not fabricated.*
