# bed-regions: performance vs `main` + the whole-genome byte-identity investigation

**Date:** 2026-06-03
**Branch:** `bed-regions` (phases 1–5; see
[doc/devel/implementation_plans/bed_regions.md](../implementation_plans/bed_regions.md))
**Question:** does routing the pileup through the new region-driven path
(one indexed `query` per contig, replacing the whole-file streaming
reader) cost performance or change output vs `main`?

## TL;DR

- **Performance: a wash on time, a big win on memory.** Whole-genome
  pileup on the GIAB HG002 bottle CRAM over GRCh38 (the **worst case**:
  2580 contigs): wall time **+1.1 %** (`--no-baq` isolation) / **373 s
  vs 372 s** (full BAQ run); **peak RSS −71 %** (993 MB vs 3431 MB).
- **Output: not byte-identical, but the branch is *more correct*.** The
  `.psp` differs by **2972 records / 5.33 M (0.056 %)**, all at contig
  ends. Root cause: a **latent walker bug in `main`** that drops the
  last reads' tail columns at every chromosome transition. The
  per-region pileup processes each contig to end-of-input and so emits
  those columns correctly.
- **Accuracy: identical.** Against the GIAB truth within the benchmark
  BED, `main` and branch produce the **exact same** TP/FP/FN
  (SNP F1 0.9037, indel F1 0.2068). The branch's +10 calls are all in
  the read-footprint overhang **outside** the evaluation BED, so they
  don't affect the score.

**Conclusion:** the branch is safe to ship — same accuracy, same speed,
71 % less memory — and it recovers coverage `main` silently drops.

## Setup

Two host (macOS) release binaries built from the same tree state:
`main` (`tmp/pvc-main`) and `bed-regions` (`tmp/pvc-branch`). Input:
`benchmarks/human_genome_bottle/crams/HG002_reads_selected_1000_rg.cram`
(599 MB, pre-restricted to the 1000-interval GIAB BED) against
`GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna`
(**2580 contigs**). `main` streams the whole CRAM through one walker;
the branch issues one indexed `query` per contig and runs a fresh
walker per contig.

## Performance

Whole-genome `pileup --no-baq` (isolates the region/contig machinery;
BAQ is identical in both):

| metric | main (streaming) | branch (per-region) | Δ |
|---|---|---|---|
| wall time | 199.3 s | 201.4 s | **+1.1 %** |
| peak RSS | 3431 MB | 993 MB | **−71 %** |
| `.psp` size | 215,991,727 B | 216,134,796 B | +0.066 % |

Full pipeline with BAQ (the real benchmark): pileup **372 s** (main) vs
**373 s** (branch); var-calling 9 s vs 10 s.

The feared per-contig overhead (2580 queries, most over empty contigs)
is negligible. The memory win is the per-region pipeline holding one
contig's walker state instead of streaming the whole genome.

## Byte-identity investigation

The whole-genome `.psp` was **not** byte-identical. Run summaries showed
**identical reads admitted** (10,585,124) and identical filter counts,
but the walker emitted **+2972 records** on the branch
(`records_emitted` 5,324,306 vs 5,327,278; `record_widen_events` +27).

**Step 1 — rule out the reader.** A diagnostic
([examples/diag_reader_compare.rs](../../../examples/diag_reader_compare.rs))
dumped the post-filter read sequence (qname/flag/pos/mapq/cigar/
adaptor_boundary) from `AlignmentMergedReader::new` (streaming, filtered
to chr1) and from `AlignmentMergedReader::query("chr1")` (indexed):

```
streaming reads: 1092464, indexed reads: 1092464
total positions compared: 1092464, differing: 0
READ SEQUENCES IDENTICAL → divergence is in the walker, not the reader
```

The indexed `query` and streaming `new` return the **exact same reads
in the same order**. The divergence is downstream, in the walker.

**Step 2 — locate the walker divergence.** Streaming the record streams
through `psp-to-pileup` and diffing positions showed the first
divergence at the **tail of chr1's covered data**:

```
main  : chr1 …248275730            → chr2 786910      (jumps)
branch: chr1 …248275730 248275731 … 248275877 → chr2 786910
```

`main` drops `chr1:248275731–248275877` — **~147 positions ≈ one read
length** — and jumps straight to chr2. The 2972 total differences ÷
~147 ≈ ~20 contig boundaries: **one tail-drop per contig boundary,
systematic.**

**Step 3 — root cause.** In
[`fill_pending`](../../../src/pileup/walker/driver.rs) the moment a read
on a *new* contig is peeked, the `chrom_transition` branch calls
`flush_chromosome_into`, which finalizes *open* records but never
advances `walker_pos` through (and emits columns for) the positions
still covered by the current contig's last reads. The **end-of-input**
path has no early flush, so it steps through those tail positions and
emits them. `main` hits the transition path at every boundary (dropping
tails) but the end-of-input path only for the **last** contig — so it
inconsistently emits the last contig's tail and drops every other's.
That inconsistency is the tell: it's an oversight, not a design choice.

The per-region pileup gives **every** contig its own walker that ends
via end-of-input, so it emits every tail column. The branch is
therefore **more correct than `main`**; the 0.056 % difference is
`main` losing coverage at ~20 contig ends.

(The buggy chrom-transition path is now **dormant** on this branch:
pileup is single-contig per query, so the walker never transitions
contigs. The fix is moot once this branch replaces `main`.)

## Accuracy (HG002 bottle vs GIAB truth)

Full pipeline (BAQ on), scored with
`benchmarks/lib/compare_to_truth.sh` (bcftools allele concordance within
the 1000-interval BED):

| class | TP | FP | FN | precision | recall | F1 |
|---|---|---|---|---|---|---|
| SNP — **main & branch** | 6389 | 974 | 387 | 0.8677 | 0.9429 | **0.9037** |
| indel — **main & branch** | 139 | 0 | 1066 | 1.0000 | 0.1154 | **0.2068** |

**Identical.** The branch made **+10 variant calls** (a strict superset
of `main` — `main` makes zero calls the branch lacks), all at contig
ends (chr9/10/11/17/21, each within ~1 Mb of the contig's end). All
**10 fall outside the evaluation BED** (`bedtools intersect`: 0 inside)
— they're in the read-footprint overhang just past each interval edge,
so they don't change the score. The low indel recall is a pre-existing
property of the caller, identical on both.

## Reproduction

```sh
# Host release binaries (main via a worktree):
cargo build --release --bin pop_var_caller            # branch → tmp/pvc-branch
git worktree add tmp/main-wt main && (cd tmp/main-wt && cargo build --release)

CRAM=benchmarks/human_genome_bottle/crams/HG002_reads_selected_1000_rg.cram
REF=$HOME/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

# Perf isolation (no-baq), per binary:
/usr/bin/time -l <bin> pileup --reference $REF --output out.psp --no-baq $CRAM

# Reader isolation:
cargo build --release --example diag_reader_compare
./target/release/examples/diag_reader_compare --reference $REF --contig chr1 $CRAM

# Full accuracy pipeline + truth comparison, per binary:
POP_VAR_CALLER_BIN=<bin> OUT_ROOT=<out> \
  bash benchmarks/lib/run_ours.sh benchmarks/human_genome_bottle/bench.config.sh single
bash benchmarks/lib/compare_to_truth.sh benchmarks/human_genome_bottle/bench.config.sh <out>/ours/single_*.vcf
```

All run artifacts (binaries, worktree, per-run outputs) live under the
gitignored `tmp/`.
