# D2/D3 — paralog summary parity + cost (tomato2)

**Date:** 2026-06-29, branch `tomato2-paralog-filter`.
Empirical validation for the `.psp` per-sample summaries
([plan](../../implementation_plans/paralog_psp_summaries.md) milestones
D2/D3). Companion to the architecture doc
[hidden_paralog_psp_integration.md](../../architecture/hidden_paralog_psp_integration.md).

## D2 — coverage parity vs the Python prototype

**Method.** The coverage histogram is a pure function of the `.psp` body
(Premise 1), so the Rust `CoverageByGcAccumulator` was re-derived from the
existing tomato2 `.psp` files (no CRAM / reference re-run needed) via the
new `examples/dump_sample_summary.rs` (which falls back to re-deriving from
the record stream when a file carries no stored section). The result was
compared against the prototype's
[`window_cov.w500.parquet`](../../../../benchmarks/tomato2/results/window_cov.w500.parquet)
(built by `assemble_window_table.py` from the same `.psp` files).

Two definitions had to agree, and do:
- the prototype's `mean_depth = sum_depth / n_cov` is **covered-bases mean
  depth** — exactly the architecture's Premise-2C choice;
- the prototype's `win_start ≡ 1 (mod 500)` matches the accumulator's
  1-based tiling `tile = (pos − 1) / 500`.

**Result** (Rust re-derive vs prototype, per sample):

| sample | Rust n_tiles | proto windows | Rust mean depth | proto mean depth | Δ |
|---|---|---|---|---|---|
| SRR5079859 | 63 266 | 63 266 | 6.054 | 6.0582 | 0.07 % |
| SRR5079876 | 63 166 | 63 166 | 4.712 | 4.7095 | 0.05 % |
| SRR5080001 | 63 248 | 63 248 | 6.308 | 6.3072 | 0.01 % |

**Window counts are exact** across all three samples; the mean-depth
deltas are ≤ 0.07 %, accounted for entirely by the stored histogram's
0.5×-wide depth bins (the Rust figure is the bin-midpoint-weighted mean,
the prototype is the exact per-window mean). The depths also match the
spec's documented "~6× per-sample" tomato2 coverage. **Coverage parity
holds** — the introgression-safe core of the filter is faithfully
reproduced.

*Het note.* There is no original-prototype table for the het classifier
to compare against: the three-way binomial-LR classification (confident
het / hom-alt / ambiguous) is the amended design (spec §3, 2026-06-29),
not the prototype's. Its correctness is covered by the unit tests
(depth-awareness: 5/6@6× ambiguous, 59/60@60× hom-alt; exact-margin
boundaries). On SRR5079859 it yields 11 138 variant sites (3 111 het /
6 770 hom-alt / 1 257 ambiguous), Hobs(confident) = 0.315 — a plausible
selfer-ish tomato het level, with the ambiguous fraction (~11 %)
reflecting the ~6× depth.

## D3 — cost

Re-deriving the **full** summary (coverage + het) over a real tomato2
`.psp` (~31.6 M covered positions = 63 266 × 500 bp), including reading and
zstd-decoding every block, took **1.25 s single-threaded** per sample
(`time` over `dump_sample_summary`). The accumulator work is O(1) per
record — one coverage observation per position plus one binomial classify
per variant site — and is a small fraction of that 1.25 s (dominated by
block decode). In the live pileup it rides the walker's per-record
processing, which is far more expensive per record, so the C2 overhead is
in the noise, matching the architecture's expectation.

**Owed:** an exact pileup wall-time / RSS before-after on real CRAM was
not run here — the tomato2 reference is out-of-tree (`$HOME/genomes`,
unreachable from the build container) and container-built Linux binaries
do not run on the macOS host, so a clean host benchmark is a follow-up
(`benchmarks/tomato2/src/build_psp.sh` with a host release build). The
re-derive measurement above is the available evidence that the
accumulators are cheap.

## Artefact

`examples/dump_sample_summary.rs` — prints a `.psp`'s stored summary, or
re-derives it from the body for an older file. Doubles as the D1 consumer
demonstration (`PspReader::metadata()` + `SampleSummary::from_toml_bytes`)
and the D2 parity tool.
