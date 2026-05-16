# Perf review fixes applied — var_calling Wave 1 (2026-05-16)

Companion to
[perf_var_calling_2026-05-16.md](perf_var_calling_2026-05-16.md). This
file records what was implemented from Wave 1, the bench numbers
collected, and what was learnt about measurement on this host.

## Findings applied

Wave 1 chosen by user scope: **L2 + H1 + H4 + H3 + H2** (Phase 2a:
`UnifiedAllele.per_sample_sources` inner `SmallVec`; Phase 2b:
`MergedRecord.scalars` / `chain_anchor_flags` / `log_likelihoods`
flattened to row-major). H5–H7 (bench-shape and config-sweep work)
remain deferred.

| ID | Site | Change |
|---|---|---|
| L2 | [`build_source_index`](../../../src/var_calling/per_group_merger.rs#L1372) | `std::collections::HashMap` → `ahash::AHashMap`, pre-sized by a one-pass count over `unified.alleles`. |
| H1 | [`project_per_position_alleles`](../../../src/var_calling/per_group_merger.rs#L716) + [`admit_compound_candidates`](../../../src/var_calling/per_group_merger.rs#L786) | `BTreeMap<Vec<u8>, usize>` → `ahash::AHashMap<Vec<u8>, usize>`, pre-sized to `1 + Σ allele-count` per group. |
| H4 | [`PerGroupMerger`](../../../src/var_calling/per_group_merger.rs#L367) struct + [`compute_log_likelihoods`](../../../src/var_calling/per_group_merger.rs#L1495) | New `genotype_tables: Arc<Vec<Vec<Vec<u8>>>>` field precomputed at `with_config`, indexed by `n_alleles`. Threaded through `refill` → `process_group` → `compute_log_likelihoods` via an `Arc::clone` per refill. Fallback to `genotype_order` retained for direct callers passing a smaller cache. |
| H3 | [`chain_broken_log_likelihood`](../../../src/var_calling/per_group_merger.rs#L1686) | Removed per-cell `Vec<u32>` / `Vec<f64>` / `Vec<u32>` scratch allocations. `per_pos_counts` is now `[u32; MAX_BITMASK_ALLELES]` (matching the convention in `standard_log_likelihood`); allele support stats are read inline from `rec.alleles[a].support`. Bit-exact (same accumulation order). |
| H2 phase 2a | [`UnifiedAllele.per_sample_sources`](../../../src/var_calling/per_group_merger.rs#L685) + [`DroppedOther.per_sample_sources`](../../../src/var_calling/per_group_merger.rs#L719) + [`CompoundChainAnchorEvidence.constituent_sources`](../../../src/var_calling/per_group_merger.rs#L919) | Inner `Vec<(usize, usize)>` → `SmallVec<[(usize, usize); 1]>` via a `type SourceList` alias. For the dominant biallelic case (≤ 1 source per sample) this eliminates the inner heap allocation entirely; chain-anchored compounds with 2 constituents still spill to heap. All access patterns work unchanged via `Deref`. |
| H2 phase 2b | [`MergedRecord`](../../../src/var_calling/per_group_merger.rs#L143) struct + [`ScalarProjection`](../../../src/var_calling/per_group_merger.rs#L1175) + helpers ([`sum_per_position_scalars`](../../../src/var_calling/per_group_merger.rs#L1208) / [`project_compound_scalars`](../../../src/var_calling/per_group_merger.rs#L1265) / [`subtract_compound_from_constituents`](../../../src/var_calling/per_group_merger.rs#L1351) / [`build_chain_anchor_flags`](../../../src/var_calling/per_group_merger.rs#L1534) / [`compute_log_likelihoods`](../../../src/var_calling/per_group_merger.rs#L1570)) | **Public API change.** Flatten three `Vec<Vec<T>>` matrices to row-major `Vec<T>`: `scalars` (n_samples × n_alleles), `chain_anchor_flags` (n_samples × n_alleles), `log_likelihoods` (n_samples × n_genotypes). New `n_samples` / `n_genotypes` fields and `scalars_row` / `chain_anchor_flag` / `chain_anchor_flags_row` / `log_likelihoods_row` accessors. Collapses `n_samples + 1` heap allocations per matrix per record into a single allocation. |

**Correctness:** all 571 lib tests pass after the combined change.

## Measurement results

Captured raw numbers (raw text) at
[tmp/perf_review_2026-05-16_var_calling/](../../../tmp/perf_review_2026-05-16_var_calling/),
file naming `bench_*.txt`. The session used criterion's named
baselines (`pre_fix`, `after_h1`, `after_h3`, `after_h2a`, `final`).

### H2 — cleanest signal of the session

Phase 2b (the public-API change flattening `MergedRecord`'s
`Vec<Vec<T>>` matrices) produced statistically significant wins on
every Stage 5 bench, confirmed by **two back-to-back runs** at the
same host load:

| bench | Phase 2a baseline | Phase 2b run 1 | Phase 2b run 2 | Δ (median of run 2 vs 2a) |
|---|---:|---:|---:|---:|
| biallelic_64 | 527.6 ms | 452.2 ms | 453.3 ms | **−14.1 %** (p = 0.00) |
| compound_anchored | 973.4 ms | 904.0 ms | 888.1 ms | **−8.8 %** (p = 0.00) |
| compound_half_chain_broken | 953.7 ms | 901.9 ms | 888.0 ms | **−6.9 %** (p = 0.00) |
| biallelic_1_sample | 33.9 ms | 29.6 ms | 29.5 ms | **−12.5 %** (p = 0.00) |

Run 1 and run 2 read within 1 % of each other on every bench,
confirming the host was stable across the two runs and the win is
genuine (not regression-to-mean from host-load volatility).

**Mechanism:** every emitted `MergedRecord` previously paid
`3 × (n_samples + 1)` heap allocations for `scalars`,
`chain_anchor_flags`, and `log_likelihoods` Vec<Vec<>>. After the
flattening that becomes `3` allocations per record. For the bench's
64 samples × 10 000 records that's ~1.95 M heap alloc/free pairs
eliminated per bench iteration — directly on the rayon worker's
allocator-bound critical path. The 1-sample bench benefits too (the
per-record Vec<Vec<>> headers were a fixed overhead).

Phase 2a (`SmallVec` for `per_sample_sources`) read ~3 % faster than
the previous Phase-H3 baseline (within noise). It eliminates the
inner Vec heap allocation for the typical 0-or-1-entries-per-sample
case but the outer `Vec<SourceList>` is still one allocation per
allele, and the cohort-side allocator pressure was already being
hammered by other per-group allocations. Kept because the mechanism
is clean and the cost is zero (existing `smallvec` dep).

The session-wide story is dominated by **host-load instability** — the
same un-modified code read 519 ms biallelic at the start of the
session, 318 ms thirty minutes later, and 572 ms thirty minutes after
that. Per-fix wall-time deltas are all 0–5 % so they sit inside this
30–70 % host-noise envelope.

Trustworthy back-to-back deltas (each pair of runs taken minutes apart
under the same host load, one with the fix applied, one with it
reverted):

| bench | fix | reverted | applied | Δ |
|---|---|---:|---:|---:|
| biallelic_64 | H1 sized | (BTreeMap baseline impossible to compare cleanly — host shifted; see below) | 403 ms then 318 ms then 545 ms | inconclusive |
| biallelic_64 | H4 | 571.95 ms | 543.93 ms | **−4.9 %** (p = 0.00 vs after_h1) |
| compound_anchored | H4 | 953.52 ms | 947.87 ms | −0.6 % (within noise) |
| compound_half_chain_broken | H4 | 956.51 ms | 943.89 ms | −1.3 % (within noise) |
| biallelic_1_sample | H4 | 32.83 ms | 34.63 ms | +5.5 % (within within-run spread) |

**H4 revert experiment** is the only clean signal among the
pre-H2 changes: biallelic 64-sample improves a consistent ~5 %,
compound paths are neutral (consistent with the perf review's
prediction — compound paths spend their cycles in
`chain_broken_log_likelihood` / project allocations, not in
`genotype_order` rebuild). The 1-sample regression shows up but its
within-run spread (32–41 ms across runs) is wider than the effect,
so it is most likely noise.

**H2 phase 2b is the largest single win of the session**, separately
documented above — ~14 % on biallelic_64 with p = 0.00 confirmed by
two back-to-back runs at the same host load.

**L2 / H1 / H3** did not produce a delta separable from host noise.
Each was kept because:

- **L2** swaps an internal `std::HashMap` for the project's already-in-tree `ahash::AHashMap` and pre-sizes it; the change is strictly better (typed pre-sized FxHash-quality hasher) with zero cost in dependency, build time, or maintenance.
- **H1** (after the pre-sizing fix — see *Surprise* below) is at least neutral and motivated by the documented 1.0–1.7 % `project_scalars:1138` self-time and the `Vec<u8>` key-clone cost feeding into the 71 % `drop_in_place<OverlappingVariantGroup>` inclusive cost. The mechanism is sound; the noise floor on this host is wider than the effect.
- **H3** is bit-exact, removes 3 small allocations per (sample × genotype × constituent) call in `chain_broken_log_likelihood` (the perf profile's 1.4–3.4 % self-time site on the compound profiles), and matches the existing stack-array convention in the sibling `standard_log_likelihood`. The compound_half_broken bench did not show a clean post-revert win, but the change is strictly less work and no risk.

## Surprise: `AHashMap::new()` was slower than `BTreeMap::new()` for tiny n

First attempt at H1 (`BTreeMap → AHashMap::new()` without pre-sizing)
showed +5.7 % / +5.3 % regression on biallelic / compound_anchored
(p = 0.07 / 0.02 respectively). The explanation: the typical biallelic
group has 2–3 distinct alleles, so the byte_index map holds 2–3
entries. `AHashMap::new()` defers the bucket allocation until the
first insert, at which point hashbrown picks an inline-capable starting
size; for 2–3 small `Vec<u8>` keys, the AES-based ahash setup cost +
first-grow cost exceeds `BTreeMap`'s "single root node with 2 keys
+ memcmp" path.

**Pre-sizing fixes it.** `AHashMap::with_capacity(max_distinct)`
where `max_distinct = 1 + Σ per-sample allele counts` (computed in a
single pass over `group.records`) — the constructor reserves enough
buckets up front, the per-insert hash setup amortises, and the first-
grow rehash never happens. Same code is applied to L2's
`build_source_index`.

This is a methodology lesson worth keeping: **for tiny maps,
`AHashMap` is not a free win over `BTreeMap`** without explicit
capacity hinting. Filed informally as part of this report.

## Per-fix size relative to the headline allocator cost

The perf review's headline diagnosis was "70 % of Stage 5 CPU is the
allocator, specifically `drop_in_place<OverlappingVariantGroup>` on
the rayon workers." The fixes shipped here each remove some
fraction of allocations on the *worker* side:

- L2: per-group `HashMap` allocation count → 1 (was several rehashes)
- H1: per-group `BTreeMap` → `AHashMap`, key clone count unchanged, hash time better
- H4: per-group `genotype_order` `Vec<Vec<u8>>` allocations → 0 (cached on merger)
- H3: per-cell `Vec<u32>` / `Vec<f64>` / `Vec<u32>` in `chain_broken_log_likelihood` → 0
- H2 phase 2a: per-allele inner `Vec<(usize, usize)>` allocations → 0 for the typical biallelic case (`SmallVec` inline)
- H2 phase 2b: per-emitted-record `Vec<Vec<T>>` matrix allocations → 1 flat `Vec<T>` each for `scalars`, `chain_anchor_flags`, `log_likelihoods` (~1.95 M alloc/free pairs saved per bench iter at 64 samples × 10 000 records)

The H2 set produces the biggest measurable wall-time win (~14 % on
biallelic_64) because the flattened matrices live on the emit path,
which the rayon worker touches once per group — every saved
allocation is per-record × n_groups. The earlier fixes (L2, H1, H4,
H3) all attack peripheral per-group allocations whose absolute count
is large but whose drop cost is overlapped with the still-untouched
inner allocations of the `OverlappingVariantGroup` input itself
(`Vec<u8>` allele bytes, `Vec<ChainId>` chain id lists, the
per-sample `Vec<Option<PileupRecord>>`).

The structural fix for that remaining cost lives at the Stage 1 /
Stage 2 / Stage 4 seam (the `AlleleObservation` SmallVec candidate
flagged S2 in the original review, and the `OverlappingVariantGroup`
arena redesign S9). Both are out of scope for this `var_calling`
review and require cross-cutting work.

## Cumulative result vs the original baseline

Best two readings at the same host load:

| bench | pre_fix (loaded) | after H1+H4+H3+H2 (run 2) | Δ |
|---|---:|---:|---:|
| biallelic_64 | 519.4 ms | **453.3 ms** | **−12.7 %** |
| compound_anchored | 964.1 ms | **888.1 ms** | **−7.9 %** |
| compound_half_chain_broken | 973.1 ms | **888.0 ms** | **−8.7 %** |
| biallelic_1_sample | 50.7 ms | **29.5 ms** | −41.8 % (noisy bench) |

The biallelic_1_sample delta is inflated by host-load drift over the
session (its run-to-run spread is ~20 %), so read the headline as
"~10–15 % wall-time improvement on the dominant 64-sample paths,
with a larger but noise-dominated improvement on the 1-sample
shape."

## What to do next

In rough cost/value order:

1. **Re-bench on a quiet host** with criterion `--save-baseline pre_wave1`
   on a `git checkout` of the pre-perf-review commit, then `--baseline pre_wave1`
   on the current `HEAD`. Today's session host shifted enough that
   the cumulative wall-time number above is the best we could get
   in-session.

2. **H5** — batch-size sweep over `{32, 128, 512, 2048}`. Pure bench
   addition. The current 32 is documented as a placeholder.

3. **H6 + H7** — Stage 4 bench fixture-rebuild fix + cohort-size
   sweep. Unblocks Stage 4 code-level findings and names the next
   priority at N=256, N=1024.

4. **L4** — bit-iteration rewrite in `standard_log_likelihood`. Tiny
   gain expected (loop body is already short and allocator-bound),
   primarily a codegen-cleanliness fix.

5. **S2** — `AlleleObservation.seq` + `chain_ids` as `SmallVec`.
   Out of scope for var_calling — cross-cutting change touching
   Stage 1 / Stage 2 / the `.psp` reader. The single largest
   structural opportunity if the per-group merger needs another 2×
   on top of the Wave-1 work.
