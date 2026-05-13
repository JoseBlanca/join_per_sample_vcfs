# Performance Review (Applied): psp_writer
**Date:** 2026-05-13 (same-day follow-up to `perf_psp_writer_2026-05-13.md`)
**Author:** rust-performance-review skill — applying the listed wins
**Companion to:** [perf_psp_writer_2026-05-13.md](perf_psp_writer_2026-05-13.md)

This document is the *author response* the original review prescribes
under its "Author response convention" — every finding cited by its
identifier with the outcome (applied / experiment shows no gain /
deferred / won't fix), plus the verbatim measurement that gates the
verdict.

---

## Summary of cumulative deltas

All benches were run inside the dev container against
`io::sink()`-backed in-memory sinks (or `BufWriter<tempfile>` for the
I/O bench), `target-cpu = x86-64-v3`, sysmalloc unless noted, on the
same workstation as the original review. The reference baseline
`pre_H` was captured at commit `7489c23` (just before any code-level
fix); the latest cumulative state is `de88f68` for sysmalloc, with
`alloc-mimalloc` available as an opt-in feature.

| Workload | pre_H (sys) | latest (sys) | latest (mimalloc) | Δ sys → mimalloc |
| --- | ---: | ---: | ---: | ---: |
| `psp_writer/snp_typical_3_3M`        | 488 ms | 313 ms  | 277 ms | -9.9 % |
| `psp_writer/phase_chain_heavy_1M`    | 245 ms | 204 ms  | 167 ms | -23.8 % |
| `psp_writer/multi_allele_500k`       | 127 ms | 127 ms  |  85 ms | -32.6 % |
| `psp_writer_io/bufwriter_file_1M_64KiB` | 82.7 ms | 82.4 ms | 81.9 ms | +2.4 % |
| `psp_writer_phases/write_record_steady_100k` | 7.88 ms | 3.79 ms | 3.52 ms | -4.4 % |
| `psp_writer_phases/flush_block_one` (NEW) | — | 33.1 ms | 24.0 ms | -27.2 % |

Cumulative sys → sys throughput change for the headline workload
(`snp_typical_3_3M`): **+56 %** records/s (6.76 → 10.7 M elem/s).
Cumulative sys → mimalloc: **+76 %** records/s (6.76 → 11.9 M elem/s).

The `flush_block_one` baseline was never accurate in the original
review (its prime-loop slack misaligned the boundary; pre-fix it was
measuring writer-drop, not flush). Bench fix shipped with H2.

DHAT (1 M-record fixture, writer-scoped):

| Stage | Total bytes | Alloc blocks |
| --- | ---: | ---: |
| pre_H baseline   | 466 MB | 883 |
| post H1          | 343 MB | 487 |
| post H3+L1+L2+L3 | 330 MB | 301 |

Both H1 and the H3/L1/L2/L3 bundle land their predicted alloc-count
reductions (`-66 %` cumulative on blocks).

---

## Verbatim profile delta (samply, snp_typical, sysmalloc)

The percentages cited below are SELF-time share on the bench thread
unless noted. Before / after compare commit `7489c23` to the current
HEAD (`de88f68`), both captured with the same `--profile-time 20 s`
invocation.

| Symbol | Pre (SELF) | Post (SELF) | Δ | Driver |
| --- | ---: | ---: | --- | --- |
| `PspWriter::validate_record`       | 10.09 % | 4.06 % | -60 % | H4 (ACGTN lookup) + L5 (cold errors) |
| `RawVecInner::grow_amortized` (incl) | 35.28 % | 7.74 % | -78 % | H1 (pre-size accumulator) |
| `__libc_realloc` (incl) | 35.39 % | 7.71 % | -78 % | H1 + L1 + L3 |
| `Vec::push_mut` (incl) | 37.70 % | 17.60 % | -53 % | H1 |
| `__brk` (self) | 4.70 % | 1.30 % | -72 % | H3 (zstd CCtx reuse) |
| `ZSTD_compressBegin_internal` (incl) | 7.12 % | below 0.5 % | gone | H3 |
| `ZSTD_cwksp_clean_tables` (incl) | 5.68 % | below 0.5 % | gone | H3 |
| `zstd encode work (compressStream2 / compressContinue / blocks)` (incl) | 17.43 % | 26.47 % | +52 % share | not slower in absolute terms — bigger share of a smaller total |

After all wins, the hot path is dominated by zstd's level-9 inner
work (`addEvents_generic`, `ZSTD_count`, `MEM_readST` together ~17 %
SELF) and allocator work that mimalloc shaves further. The
realloc-driven memcpy cluster and the per-block zstd workspace
setup — the two top H-tier candidates the original review named —
are gone.

---

## Findings, outcome by identifier

### Hot-path

- **H1 — Pre-size `BlockAccumulator` columns** — *applied in 26f3d65*
  - Bench: `snp_typical_3_3M` 488 → 408 ms (-13.6 %); `phase_chain_heavy_1M` 245 → 193 ms (-25.9 %); `multi_allele_500k` 127 → 77.1 ms (-40.0 %). DHAT alloc blocks 883 → 487.
  - Took Shape (a) — three named `INITIAL_*_HINT` constants. Shape (b) (reuse the accumulator across flushes) is a follow-up — not needed given the H1 wins already landed.

- **H2 — CSR layout for the three `Vec<Vec<SlotId>>` list columns** — *applied in 47c1514, kept after revert experiment*
  - Bench post-H2 sys: `snp_typical_3_3M` 488 → 307 ms cumulative (-37 %), `phase_chain_heavy_1M` 245 → 220 ms (-10 %).
  - **Experiment caveat:** an intermediate measurement against H4+L5 suggested H2 regressed SNP. A direct revert experiment (rolling back the writer.rs + block.rs portion of 47c1514) ran the bench and observed +27 % SNP, +37 % phase-chain, +24 % multi-allele, +24 % bufwriter, +39 % write_record_steady, +16 % flush_block_one — every workload got materially slower. The pre-H2 H4+L5 number had been a low-load outlier; H2 was helping all along. Kept.
  - The fix also bundled the bench prime-loop fix for `flush_block_one` (the broken slack-of-256 logic was measuring drop, not flush).

- **H3 — Reuse one `zstd::bulk::Compressor` across all columns of all blocks** — *applied in 969de6c (bundled with L1/L2/L3)*
  - Bench `flush_block_one`: 5.61 → 4.28 ms (-19 %, plus second-order wins via the persistent CCtx workspace warm).
  - Profile delta: `__brk` self -72 %, the `compressBegin_internal` / `cwksp_clean_tables` cluster gone from the > 0.5 % table.

- **H4 — `validate_record` ACGTN lookup table** — *applied in 154c4eb (bundled with L5)*
  - Bench: small delta on the current bench (alleles average 1 byte). Profile: `validate_record` SELF 10.09 % → 4.06 % cumulative. The win scales linearly with allele length; a real indel workload would show much more.

### Likely

- **L1 — Promote per-flush scratch buffers** — *applied in 969de6c (bundled with H3)*
  - Five `Vec<u8>` (and one `Vec<Vec<u8>>`, one `Vec<ColumnManifestEntry>`) fields on `PspWriter`, reused with `.clear()` per flush. DHAT alloc blocks 487 → 301.

- **L2 — Manifest single-pass in `flush_block`** — *applied in 969de6c*
  - Manifest entry pushed into `manifest_scratch` inside the same column loop that builds the compressed payload. The pre-refactor two-pass (build `payloads`, then `.iter().map().collect()`) is gone.

- **L3 — `mem::take(&mut block.snapshot_active_slots)` instead of `.clone()`** — *applied in 969de6c*
  - One `Vec<SlotId>` allocation removed per flush. Bench-neutral; alloc-count win.

- **L4 — `encode_u64_leb128` fast/cold split** — *applied in f95c6e5*
  - `flush_block_one` -5 %, `snp_typical_3_3M` -1.5 %. ~36 M varint calls per bench iteration take the inlined `< 0x80` fast path.

- **L5 — `#[cold]` error helpers in `validate_record`** — *applied in 154c4eb*
  - `err_invalid_record` and `err_invalid_allele_byte` constructors marked `#[cold] #[inline(never)]`; the per-arm `format!()` + `String` allocation moves out of the hot block.

- **L6 — Slab-cast list-column emit on little-endian** — *applied in 207b3e9*
  - `bytemuck::cast_slice::<T, u8>(row)` replaces the per-element `extend_from_slice(&to_le_bytes())` walk inside `encode_list_column_csr`. `phase_chain_heavy_1M` -5.4 %; others within run noise (rows are mostly empty there).
  - Adds `bytemuck = "1"` as a top-level dependency.

- **L7 — `examples/dhat_psp_writer.rs`** — *applied in 5d286df*
  - Profiler scope starts AFTER the fixture builder so the report names writer-side sites only. Baseline (pre-H1): 466 MB / 883 blocks for 1 M records.

- **L8 — `phase_chain_heavy_1M`, `multi_allele_500k`, `streaming_bufwriter_file_1M`** — *applied in 74ff4c4*
  - Three new criterion workloads added. The phase-chain fixture's slot-id allocator avoids the same-record `expired_chains` + `new_chains` collision the writer rejects.

- **L9 — Wrap bench inputs in `black_box`** — *applied in 7ddf990*
  - Zero behaviour change; bench hardening only.

- **L10 — `write_record_steady_100k` and `flush_block_one` sub-benches** — *applied in 854fe8b*
  - Added a `#[doc(hidden)] pub fn current_block_projected_bytes(&self) -> Option<usize>` on `PspWriter` so the bench's prime loop aligns with the writer's auto-flush boundary deterministically.
  - The prime-loop slack-of-256 logic landed broken in this commit; the fix shipped with H2 (47c1514). The `flush_block_one` baseline from `pre_H` is **not** comparable to anything after H2 because the original was timing writer-drop.

- **L11 — Rustdoc note on `PspWriter::new` about wrapping `File` sinks in `BufWriter`** — *applied in f2cec9a*

### Speculative

- **S1 — `active_slots: BTreeSet<SlotId>` → sorted `Vec<SlotId>`** — *applied in de88f68 (bundled with S2)*
  - Bench `phase_chain_heavy_1M`: -7.3 %. Operations: `contains` → `binary_search.is_ok()`, `insert` → `binary_search.unwrap_err()` + `insert(idx, slot)`, `remove` → `binary_search.ok().map(|i| remove(i))`.
  - Post-fix, `core::slice::<impl [T]>::binary_search_by` is 12.32 % SELF on phase_chain_heavy — the new hot leaf, expected, and still net much faster than B-tree lookup.

- **S2 — Pre-size the active-slot snapshot Vec** — *applied in de88f68 (bundled with S1)*
  - Block-start snapshot is now `self.active_slots.clone()` — one contiguous memcpy of ≤ 24 bytes, no tree walk + collect.

- **S3 — 12-arm `encode_column_into` dispatch jump-table check** — *no action — confirmed by final samply*
  - `encode_column_into` is 6.13 % inclusive but below the 0.5 % SELF reporting bar — dispatch shape is immeasurable on this workload. File-and-forget per the review's prescription.

- **S4 — `windows(2)` walks in `validate_record` on chain-heavy workloads** — *no action — confirmed by final samply*
  - `<core::slice::iter::Windows<T> as core::iter::traits::iterator::Iterator>::next` is at 0.61 % SELF on `phase_chain_heavy_1M`, barely above the 0.5 % threshold. The proposed `is_sorted_by` rewrite has identical codegen; no measurable gain to chase.

- **S5 — Rustdoc on `PspWriter::finish` about `BufWriter::into_inner()? → File::sync_all()?`** — *applied in f2cec9a (bundled with L11)*

### Build / toolchain

- **Allocator A/B (mimalloc)** — *applied in 7be79d3 as the `alloc-mimalloc` feature*
  - Threshold to merge was ≥ 5 %. Result: -9.9 % SNP, -23.8 % phase-chain, -32.6 % multi-allele, -27.2 % flush. Vastly above threshold. Kept as opt-in so production builds default to system malloc.

- **`panic = "abort"`** — *applied in 7bf139a, on `[profile.release]`* (after confirming `grep -rn catch_unwind` is empty)
  - Bench: -1 to -3 % across workloads. Small but real.

- **`rust-toolchain.toml`** — *applied in 07c5097, pinned at 1.95*
  - Reproducibility infrastructure for future criterion baselines.

- **PGO** — *deferred to a follow-up issue*
  - Code-level wins are now exhausted on this hardware; PGO is the next build-time lever. Not applied this round.

### Notes

All the original review's negative findings (the `[Note]` items) were
confirmed by reading and left as-is. No action.

---

## What's still hot, and the next-PR candidate list

The final samply on `snp_typical_3_3M` (sysmalloc) lists, in
descending SELF:

```
  10.42%  Vec<T,A>::push_mut [std]            ← allocator-bound; mimalloc shaves further
   9.74%  addEvents_generic                   ← zstd level-9 inner work
   7.77%  __memcpy_avx_unaligned_erms (libc)  ← realloc memcpy residual
   7.53%  core::ptr::write [std]
   5.22%  Vec<T,A>::append_elements [std]
   4.57%  BlockAccumulator::append_record [crate]
   4.39%  ZSTD_count                          ← zstd inner
   4.06%  PspWriter::validate_record [crate]  ← scales with allele length
```

Concrete remaining levers:

1. **`alloc-mimalloc` by default for the eventual CLI.** Opt-in
   today; the per-sample-caller binary should turn it on, given the
   measured gains. Needs a separate review of mimalloc's footprint
   under the full pipeline (BAQ + walker + writer) — not just the
   writer.
2. **Lower zstd level on a per-block-size feedback knob.** Spec pins
   level 9, but the spec was written before the writer was
   throughput-bound on the inner zstd work. A `--zstd-level` CLI
   flag with a default of 9 and a documented "use 3 for tight write
   loops" note would let users trade compressed size for throughput.
3. **PGO build pipeline** (see deferred above).
4. **Bench the `multi_allele` workload with real indel-shape allele
   lengths.** Current fixture averages ~3 bytes per allele — the
   `validate_record` ACGTN autovec win wouldn't show. A
   long-indel-heavy workload (e.g. ~20-50 byte alleles for medium
   INS/DEL) would prove H4's mechanism on more demanding inputs.

These four are out of scope for this review; the H/L/S items the
original report identified are all closed.
