# SSR fast-path investigation — should clean reads skip the pair-HMM?

**Date:** 2026-06-16
**Branch:** `ssr-pileup-review`
**Question:** other STR tools (HipSTR/GangSTR/ExpansionHunter) have a fast exact-count
path for clean reads and only realign ambiguous ones. The realign-everything
design parked our `count_pure_tiling` fast path as "try + measure". Now that the
pair-HMM is ~2× faster (H1/H2/P1), is a fast path worth wiring?

**Verdict: No — do not build a separate fast path.** The investigation instead
surfaced a simpler, mechanism-preserving win (the candidate window is ~2× wider
than needed) and a re-prioritization (on real data the realignment is often *not*
the end-to-end bottleneck — the fetch/catalog-walk is).

---

## Method (real-data fixture + instrumentation)

Built the missing end-to-end fixture (the review's open **L9**), on real tomato data:

- `trf-mod` compiled from the vendored `TRF-mod/` (`make -f compile.mak CC=clang`).
- A real ch01 SSR catalog: `ssr-catalog --reference <ch01.fa> --trf-mod-path TRF-mod/trf-mod` (~74 s; ch01 extracted with `samtools faidx`).
- Ran `ssr-pileup` on a real CRAM (`benchmarks/tomato1/crams/SRR7279482.p1.bench.cram`, full tomato reference, ch01 catalog).

Instrumented `analyze_read` (throwaway env-gated probe, since reverted) to tally,
per spanning read: the HMM's pruned profile (argmax length + #surviving lengths),
and whether several cheap fast-path predicates would have agreed.

## Data (1243 spanning reads, real CRAM)

```
HMM single-length result:                94.0%   (most reads → one confident length)

predicate A  region == a rung exactly:   20.9% recall, 98.8% precise
predicate C  isolated tract is a tiling:  30.2% recall, ~79% precise (naive isolation)
predicate D  observed_count == HMM argmax: 84.7% — BUT see below

min window to keep ALL surviving mass (|surviving − observed_count|):
   ≤2: 87.0%   3–4: 8.6%   5–6: 3.1%   >6: 1.3%
```

End-to-end wall on the real CRAM, varying the candidate window (CLI `--window`):

```
window=10 (21 rungs):  97.6 s   output 1,569,679 B
window=6  (13 rungs):  93.8 s   output 1,569,631 B   (−4% wall, output −0.003%)
window=4  (9 rungs):   86.0 s   output 1,569,503 B   (−12% wall, output −0.01%)
```

## What the data says

1. **The HMM resolves 94% of reads to a single length** — so most of the DP's
   precision is "wasted" on reads whose answer is trivial. Big headroom, in
   principle.

2. **But a *safe* fast path has low recall.** The only high-precision predicate
   (A: the read's window exactly equals a rung, ~99% precise) fires on just ~21%
   of reads — and is still not byte-identical. Skipping ~21% of an *already
   cheap* (post-P1) HMM is a modest win for a byte-identity break.

3. **The high-recall predicate is unsafe — and it confirms why realign-everything
   exists.** `observed_count` (the triage longest-run, free) matches the HMM for
   85%, but **~11% of reads are HMM-single with a *different* length than
   `observed_count`** — i.e. the HMM *corrected* the naive count. A fast path
   trusting `observed_count` would emit the wrong length on exactly those reads:
   the mapper/pre-probe mis-count that realign-everything was designed to fix.
   (`count_pure_tiling`-style isolation is also the hard part — predicate C's
   naive tract isolation mis-fired ~18% of the time, because isolating the tract
   *is* the alignment problem the HMM solves.)

4. **The real lever is the over-wide window, not a fast path.** The default
   `window = 10` is a documented calibration placeholder; the data shows
   **98.7% of reads keep identical surviving mass at window=6, 95.6% at window=4.**
   Tightening the window is the same speedup the fast path promised, but
   mechanism-preserving (the HMM still realigns and still corrects), a single
   parameter, and nearly output-identical (sizes above). It is still a mild
   approximation — the ~1.3% of reads needing window>6 are the large-correction
   tail — so it wants a genotype-concordance check before the default changes.

5. **Realignment is often not the end-to-end bottleneck.** On this real workload
   the window barely moved the wall (−4% at window=6) because realignment is a
   small fraction of pileup time — the catalog walk + per-locus indexed fetch +
   CRAM decode dominate. (The fixture's ch01-wide catalog over region-limited
   reads over-weights empty-locus fetch, but it directly corroborates the review's
   open **L5/L7** findings: the fetch path is where the end-to-end time goes.)

## Window concordance — measured, and the default lowered 10 → 6 (APPLIED)

Record-level diff of the `.ssr.psp` at each window vs window=10 (10425 loci,
1243 spanning reads; `ssr_psp_concordance` ignored test):

```
window  loci w/ any read diff   reads: profile differs   ARGMAX (called length) differs
  8      4 (0.04%)               10 (0.80%)               3 (0.24%)
  6      6 (0.06%)               16 (1.29%)               6 (0.48%)
  4      12 (0.12%)              55 (4.42%)               33 (2.66%)   <- elbow: don't go here
```

Realignment-only cost scales ~linearly in rung count even post-P1
(`ssr_realign/window`: 7 rungs 38 ms, 21 rungs 125 ms, 31 rungs 207 ms), so
`window=6` (13 rungs) is **~40% cheaper realignment** than `window=10` (21 rungs).
The elbow is below 6 (4→2.66% argmax change); 6 keeps 99.5% of read calls and
99.94% of loci identical, and ±6 brackets far more stutter than the ±1–2 units
real STRs show.

**Applied:** `DEFAULT_WINDOW` 10 → 6 ([driver.rs](../../../../src/ssr/pileup/driver.rs)).

## Recommendation

- **Don't wire a standalone fast path** (above).
- **Default window lowered 10 → 6** — done, validated on real data.
- **Prioritize the fetch path (L7) for end-to-end wins** — on real data it, not
  the realignment, is the wall. The realignment wins (H1/H2/P1 + the window
  tightening) pay most on high-coverage data where every locus has reads.

## Fetch path — profiled (the real end-to-end lever)

Single-threaded `sample` profile of the pileup over the real fixture
(`tmp/ssr_fix/fetch_sample.txt`). Excluding the rayon/idle scaffolding, the work
self-time is **per-locus-query CRAM decoding**, all under `CramSegmentReads::next`:

```
noodles_cram Slice::decode_blocks + Block::decode + Slice::records   ~18k samples  (CRAM slice decode)
md5::compress::compress                                               7478 samples  (~29%) reference-MD5 validation
alloc/realloc/free                                                    ~5k samples
```

The MD5 is noodles validating each slice's reference subsequence against the
slice-header `M5` — **per slice decode, hardcoded, no disable flag** (verified in
noodles-cram 0.93 `container/slice.rs:359`; the fasta `Repository` caches the
*sequence* but not the MD5). So the MD5 cannot be skipped; it can only be paid
fewer times by **decoding each slice fewer times**.

**Root cause = redundant decode.** One indexed query per locus; adjacent catalog
loci fall in the same CRAM slice and re-decode (and re-MD5) it. Coalescing a run
of adjacent loci into one query — decode the slice once, demux reads to each
locus — removes the redundancy (review finding **L7**), cutting *both* the decode
and the MD5.

**The catch (byte-identity):** the per-locus QC scalars `n_filtered` /
`mapped_reads` are the reader's filter-drop tallies *for that locus's window*
([fetch_reads.rs](../../../../src/ssr/pileup/fetch_reads.rs#L179-L190)). A single
union query returns one aggregate drop count for the whole group — the drops
can't be split back per locus. So byte-identical coalescing requires **moving the
read filter out of the reader and into the SSR fetch layer**: query the union
span *unfiltered* (primary reads only), then per locus apply window-overlap +
the MAPQ/dup/qc/length filter + reach-gate + reservoir. `yielded`/reservoir are
trivially recomputable; the filter relocation is the real work and must
byte-match the reader's `classify_segment_record` quality/dup logic. (Decode
caching in the shared pooled reader is the alternative, but it's cross-cutting
with the SNP path and has the same per-window-attribution issue.)

**Status:** profiled and designed; not yet implemented. It is a contained-but-
non-trivial change (new union-fetch + SSR-side filter, restructure the driver's
par_iter unit from locus → locus-group, byte-identity gate vs the per-locus path).
Recommended as the next focused step. On this sparse fixture most loci are empty
(catalog spans 90 Mb, reads only in ~900 kb of bench regions), which over-weights
empty-locus queries; the win is largest on dense whole-genome runs where every
locus decodes.

## Reproducing the fixture

`tmp/` is gitignored, so the fixture is local-only. To rebuild:
`TRF-mod` → `make -f compile.mak CC=clang`; `samtools faidx <ref> SL4.0ch01 > tmp/ssr_fix/ch01.fa`;
`ssr-catalog --reference tmp/ssr_fix/ch01.fa --output tmp/ssr_fix/ch01.ssr.catalog --trf-mod-path TRF-mod/trf-mod`;
`ssr-pileup <cram> --reference <full-ref> --catalog tmp/ssr_fix/ch01.ssr.catalog --output <out>`.
