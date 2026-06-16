# SNP `--regions` retrofit onto the segment read fetcher (implementation plan)

**Status:** plan, 2026-06-16. Increment **#5** of the segment-read-fetcher
work (the standalone primitive — increments #1–#3 — is merged to `main`;
the SSR Stage-1 wiring, #4, is independent and lives on its own branch).
This increment replaces the SNP pileup's per-region reader so that a
fragmented `--regions` BED stops paying the per-region file-reopen tax.

**Goal.** Make the Stage-1 `pileup` path open each input alignment file
**once** and reuse a pooled, re-seekable reader across all regions, instead
of re-opening the file (and re-walking the index) on every region — while
keeping the `.psp` output **byte-identical**. The mechanism already exists:
[`AlignmentFile`](../../../src/bam/segment_reader.rs)
(`BamFile`/`CramFile`) pools re-seekable readers. This increment routes the
pileup's per-region read fetch through it and retires the per-region
re-open path.

Cross-refs: the primitive itself
([segment_read_fetcher.md](segment_read_fetcher.md)); the reader-substitution
analysis that motivated this increment (this plan's §1–§2 distil it); the
`--regions` perf tax (`project_bed_regions_pileup_time_tax`); the
`bed_regions` feature that made the pileup region-driven
([bed_regions.md](bed_regions.md)).

---

## 1. What the pileup uses today

The production `pileup` path uses exactly one reader:
[`AlignmentMergedReader::query`](../../../src/bam/alignment_input.rs) (the
indexed path; `::new`/whole-file streaming is test-only since `bed_regions`).

- `run_pileup` ([cli.rs](../../../src/pop_var_caller/cli.rs)) builds a
  `RegionSet`, then loops **serially** over regions calling
  `AlignmentMergedReader::query(...)` **once per region**.
- `query` is backed by the per-format scanners
  **`OwnedIndexed{Bam,Cram}Records`**
  ([bam_input.rs](../../../src/bam/bam_input.rs),
  [cram_input.rs](../../../src/bam/cram_input.rs)). Each per-region call
  **re-opens the input file(s)** (`build_from_path` + re-read header) and
  re-walks the index. Headers, indexes, and the FASTA `Repository` are
  cached across regions; the file handles and the index walk are not.
- That per-region re-open is the `--regions` tax (~14% on indexed query;
  concentrated on **many small regions on one contig** — fragmented BEDs;
  negligible for whole-genome, where the region set is ~one span per contig).
- The region loop is **serial**; parallelism lives *inside* each region's
  BAQ/walker chain. So a pooled reader is borrowed one-at-a-time on the main
  thread and the pool holds ~one reader per input after warmup.

## 2. The contract `query` provides (must be preserved)

`query` is not just a fetcher — it is a merge/dedup/order/filter/count
assembly, and the pileup depends on every part. A replacement must preserve:

1. **Coordinate-sorted output** `(ref_id, pos)` — the walker is a forward
   sweep ([walker/driver.rs](../../../src/pileup/walker/driver.rs)); regressing
   input breaks it. **Load-bearing.**
2. **Multi-file k-way merge** — one sample may have several input files
   (e.g. per-lane BAMs); `query` merges them into one sorted stream.
   **Load-bearing** (multi-file is a supported input).
3. **Cross-file duplicate detection** — `ReadFingerprint`; **errors** with
   `DuplicateReadAcrossFiles` on a collision (not a silent drop).
4. **Within-file out-of-order detection + fuse-on-error** — errors on a
   coordinate regression, then fuses (`FusedIterator`).
5. **Cheap filter cascade** — `classify_pre_decode` (flags / MAPQ / dup) +
   `min_read_length`. The expensive filters (F1 mismatch, G2 bad-CIGAR, F3
   indel-norm) are **not** here; they already run downstream in
   `read_processor`. **Load-bearing for byte-identity.**
6. **`FilterCounts`** — accumulated and merged into the run summary
   ([stage1_pipeline.rs](../../../src/pop_var_caller/stage1_pipeline.rs)).
   **Reported to the user.**
7. **`ref_id` = canonical `ContigList` index** on `MappedRead` (and
   `mate_ref_id`). **Load-bearing.**

The segment reader already satisfies 1 (within one file/region), 5 (same
`classify_pre_decode` + `min_read_length` via `SegmentReadFilter`), and 7.
It does **not** provide 2, 3, 4, or 6 — those live in the merge layer.

## 3. Target architecture (the merge moves to the `MappedRead` level)

`AlignmentFile` replaces only the inner *"fetch the reads for one file +
region"* layer. The merge/dedup/order/count machinery stays — but it moves
from operating on `RecordBuf` streams (today) to operating on the
`MappedRead` streams the segment reader emits. This is the convergence the
primitive plan foresaw ("retire `OwnedIndexed{Bam,Cram}Records`").

```
PileupInputs (built once)
  └── one pooled AlignmentFile per input file  (open-once, re-seek per region)

per region (serial loop, unchanged):
  for each input:  AlignmentFile::get_reads_from_segment(chrom, start, end)
                       → MappedReadsInSegment  (cheap-filtered, coordinate-sorted, MappedRead)
  k-way merge the per-input MappedRead streams  ← NEW: merge at MappedRead level
     · argmin by (ref_id, pos)        (coordinate order)
     · cross-file dedup on (qname, flag, ref_id, pos)   → DuplicateReadAcrossFiles
     · within-file out-of-order check → OutOfOrderRead; fuse on first error
     · accumulate FilterCounts
  → the existing BAQ/walker chain (unchanged)
```

Why `MappedRead`-level merge (and not "lift only the pool under the existing
`RecordBuf` merge"):
- The dedup keys (qname, flag, ref_id, pos) and the order key (ref_id, pos)
  are all already `MappedRead` fields — no `RecordBuf` needed.
- The cheap filter already runs inside the fetcher, so `refill_heads`'
  `classify_pre_decode` pass disappears rather than being duplicated.
- It lets the old `RecordBuf`-streaming `query` + `OwnedIndexed*Records` be
  deleted outright instead of kept as a parallel shape.

The trade: the merge is re-homed, so the subtle emit-time semantics must be
reproduced exactly (see §5).

## 4. Increment sequence (each its own commit + tests)

1. **`FilterCounts` in the fetcher.** Teach `AlignmentFile` /
   `MappedReadsInSegment` to tally the cheap-filter drops it currently
   discards (per `FilterBucket`: duplicate / low-mapq / supplementary /
   secondary / unmapped / qc-fail / too-short), exposed via an accessor.
   Smallest, lowest-risk, independently testable step; closes the one
   functional gap (contract item 6). No consumer change yet.
2. **`MappedRead` k-way merge over `AlignmentFile`.** Build the merge layer
   (argmin order + cross-file dedup + within-file out-of-order + fuse + count
   reconciliation), consuming one pooled `AlignmentFile` per input. Land it
   behind the existing `AlignmentMergedReader::query` *call signature* (or a
   sibling the driver can switch to) so the swap in step 3 is a one-line
   change. Gate on a **byte-identity test** (§5).
3. **Flip `run_pileup` + retire the old path.** Hold pooled `AlignmentFile`s
   in `PileupInputs`; per region, drive the new merge. Then delete
   `OwnedIndexed{Bam,Cram}Records` and the now-dead `query` internals.
   Measure the `--regions` tax delta on a fragmented BED (§6).

Single-file fast path (fold into step 2): when a sample has exactly one
input, skip the k-way merge and drain one `MappedReadsInSegment` directly —
but **re-add the within-file out-of-order check** the segment reader does not
do today (it trusts the coordinate sort), so a corrupt unsorted single file
still errors instead of silently mis-piling.

## 5. Byte-identity (the gate)

The `.psp` output is contractually byte-identical and verified out-of-tree;
this swap must not change the surviving read set or its order. Risk areas to
pin with tests **before flipping the default**:

- **Dedup-vs-filter emit ordering.** The current merge has a deliberate
  subtlety: a read later dropped by G2/F1 (downstream) still occupies the
  at-locus dedup buffer, and that is documented as harmless because cross-file
  duplicates share the same verdict
  ([alignment_input.rs](../../../src/bam/alignment_input.rs) — the
  `record_buf_to_mapped_read` call-site comment). The new merge must
  reproduce *when* a read is counted, deduped, and emitted, not just the final
  set.
- **`FilterCounts` parity.** Step 1's tallies must equal the old per-region
  `filter_counts()` sums exactly (same buckets, same order of application).
- **Multi-file merge parity.** Same argmin tie-breaking and same
  `DuplicateReadAcrossFiles` / `OutOfOrderRead` triggers.

Add an in-tree oracle test: drive the **old** path and the **new** path over
the same multi-file, fragmented-BED fixture and assert identical
`Vec<MappedRead>` (+ identical `FilterCounts`). This also closes the
standing "byte-identity unguarded by a test" review debt.

## 6. Performance measurement

The win is the eliminated per-region re-open, so measure where it bites:

- **Fragmented BED** (many small regions on one contig) — the headline case;
  expect the bulk of the ~14% `--regions` tax to recover.
- **Whole-genome** (region set ≈ one span per contig) — expect ~neutral
  (few re-opens to begin with); confirm no regression.
- Real-data: a human bottle CRAM and a tomato cohort BAM/CRAM, `--regions`
  vs the no-regions baseline, before/after. Record both wall time and that
  the `.psp` is byte-identical.

(No `benches/` harness reaches this path yet; a criterion microbench over a
synthetic many-region fixture would localise the open-cost win below the
end-to-end noise — optional.)

## 7. Tests

- **Byte-identity oracle** — old `query` vs new merge, multi-file +
  fragmented BED, identical reads + `FilterCounts` (§5).
- **`FilterCounts` in the fetcher** — each bucket counted (a low-mapq read, a
  duplicate-flag read, a too-short read, …) for both BAM and CRAM.
- **Multi-file merge** — two inputs for one sample, interleaved coordinates,
  yields the merged sorted stream; a cross-file duplicate errors
  (`DuplicateReadAcrossFiles`).
- **Within-file out-of-order** — a single-file path with a coordinate
  regression still errors + fuses (the re-added check).
- **Pool reuse across regions** — N regions on one contig open each file
  **once** (assert via the pool's resting size / an open counter); the
  contig-transition seek path works (region on contig A then contig B reuses
  the same pooled reader, re-seeked).
- **Per-region clamp interaction** — a read spanning a region boundary is
  still written exactly once (the walker clamp in `pileup_to_psp` is
  unchanged, but confirm the overlap-yield from the fetcher feeds it the same
  reads).

## 8. Deferred (explicitly not this increment)

- **Binary-search the `.crai`/index head** — start each per-region fetch by
  jumping to the first relevant container/chunk instead of scanning the
  in-memory index from 0. Measure-first; only worth it if a profile shows the
  in-memory head-scan matters at fragmented-BED scale (the index is parsed
  once into RAM, so this is comparisons, not I/O — see the segment-reader
  review discussion).
- **CRAM container cache** (decode-once across adjacent same-container
  regions) — `segment_read_fetcher.md` §8.
- **Carry the index cursor forward across regions** — since regions are
  processed in sorted order, a per-file cursor could avoid even the
  binary-search; sidesteps the back-up/`Option`-ref-id edge cases. Only if
  the head-scan proves to matter.

## 9. Open items to pin while coding

1. Whether the new merge lands as a re-implemented `AlignmentMergedReader::query`
   body or a new sibling type the driver switches to (favour the latter for a
   clean A/B during the byte-identity gate, then collapse).
2. Where the pooled `AlignmentFile`s live — `PileupInputs` is the natural
   home (built once at startup, one per input file); confirm its lifetime
   spans the whole region loop.
3. Exact reconciliation point for `FilterCounts` between the fetcher's cheap
   tallies and the downstream `read_processor` G2/F1 tallies in the run
   summary (today both feed `stage1_pipeline`'s merge).
4. Whether to keep the old `query` path behind a flag for one release as the
   byte-identity oracle, or delete it once the in-tree oracle test lands.
