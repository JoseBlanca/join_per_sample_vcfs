# SSR Stage-1 fetch — per-thread `CramReader` with a slice cache (decode-once)

**Status:** settled design, 2026-06-17 (worked through in discussion). A focused
optimization of the [`ssr-pileup`](ssr_pileup.md) Stage-1 fetch path (the §8 read
pipeline). It does **not** change the evidence model, the catalog, the
`.ssr.psp` format, or the realignment — only *how reads are pulled from the
CRAM/BAM for each locus*. Goal: remove the dominant end-to-end cost (redundant
CRAM slice decoding) while keeping the per-locus output **byte-identical**.

Grounded in measurement: the fetch path was profiled on a real tomato catalog +
CRAM in
[ssr_fastpath_investigation_2026-06-16.md](../reports/reviews/ssr_fastpath_investigation_2026-06-16.md).

---

## 1. The problem (measured)

`ssr-pileup` issues **one indexed query per locus**
([fetch_reads.rs](../../../src/ssr/pileup/fetch_reads.rs) `fetch_locus_reads`):
for locus `L` it asks the reader for reads overlapping `L`'s ~150 bp window
(tract ± flank), reservoir-caps them, and analyzes.

A CRAM file stores reads in **slices** — a slice is a fixed batch of ~10 000
reads, stored columnar + compressed + reference-differential. A slice is the
**atomic unit of decode**: to read *any* read out you must decompress the whole
slice's blocks, reconstruct all ~10 000 records, and (in noodles) MD5-hash the
slice's reference subsequence. At ~30× coverage a slice spans roughly
`10 000 × 150 / 30 ≈ 50 kb` of genome.

SSR loci are small and dense (~1 per few kb), so **~10–25 loci fall inside one
slice**. The per-locus loop therefore decodes the *same* slice once per locus:

```
locus L1 (slice #42):  decode #42 in full → keep ~30 reads, discard ~9 970
locus L2 (slice #42):  decode #42 AGAIN   → keep ~30 reads, discard ~9 970
 ...                    (≈ 10–25 loci, all re-decoding #42) ...
```

noodles caches the reference *sequence bytes* but **not** the decoded slice and
**not** the MD5 result (verified in noodles-cram 0.93 `container/slice.rs:359`),
so every query re-does the full decompress + re-hash.

**Profile (single-thread, real fixture), work self-time excluding rayon idle:**

```
noodles_cram slice decode (decode_blocks + Block::decode + Slice::records)  ~18k samples
md5::compress  (per-slice reference-MD5 validation, hardcoded, no opt-out)   7.5k samples (~29%)
```

The MD5 cannot be disabled; both costs scale with the *number of slice decodes*.
The "useful" work (the ~30 reads each locus keeps) is a tiny fraction of each
decode. **The fix is to decode each slice once and serve every locus in it.**

## 2. The design in one picture

```
                         per worker thread (own file handle + own cache)
 loci (sorted) ──►  ┌──────────────────────────────────────────────────────┐
                    │ CramReader { file, index, repo, filter, slice_cache } │
   for each L:      │                                                       │
     fetch(L.win) ──┼─► slices overlapping L.win                            │
                    │      each slice:  cache hit → reuse                    │
                    │                   cache miss → decode + insert (FIFO)  │
                    │      classify reads in slice vs L.win:                 │
                    │         overlap? filter? → count + (copy out passers)  │
                    │   → (reads, FilterCounts)                              │ ─► reservoir → analyze
                    └──────────────────────────────────────────────────────┘
```

Two pieces:

1. **A per-thread `CramReader`** (§3) — owns the open file + a small FIFO cache of
   decoded slices, and serves per-locus windows through the existing read filter.
   Lock-free (one per worker, like `LocusScratch`).
2. **An optional `covered_regions()`** (§5) — exposes the slices' genomic spans so
   the orchestrator can *place* loci of the same slice on the same thread,
   maximizing cache hits and ensuring different threads decode different slices.
   A correctness-independent hint, added as a second step.

## 3. `CramReader` — the per-thread reader + slice cache

Drop-in for today's `AlignmentFile::get_reads_from_segment` + `filter_counts()`:

```text
CramReader::new(path, header, index, repository, filter, max_cached_slices = 3)
    // owns its own file handle; shares header / .crai index / fasta Repository via Arc

fn fetch_mapped_reads(&mut self, chrom, start, end) -> (Vec<MappedRead>, FilterCounts) {
    let mut reads  = Vec::new();
    let mut counts = FilterCounts::default();
    for slice in self.slices_overlapping(chrom, start, end) {   // from the .crai
        let decoded = self.slice_cache.get_or_decode(slice);    // FIFO, decode on miss
        for rec in decoded.records() {
            if !overlaps(rec, start, end) { continue; }
            match classify_segment_record(rec, &self.filter) {  // SHARED with the reader
                Keep(mapped)   => reads.push(mapped),           // copy out (never move)
                Drop(category) => counts.bump(category),        // count EVERY drop
            }
        }
    }
    (reads, counts)
}
```

Design decisions:

- **Per-thread, lock-free.** One `CramReader` per rayon worker (via `map_init`,
  exactly like `LocusScratch`). No shared cache → no lock on the per-locus hot
  path, no cache-stampede coordination. (Why not a shared cache: §7.1.)
- **The cache stores decoded *records*, not finished `MappedRead`s.** To count
  *every* filtered read, the classifier must run over all of a slice's reads —
  including ones that get dropped, some of which (unmapped-but-placed) cannot
  become a `MappedRead`. So the cache holds the decoded records (enough to
  classify and to build a `MappedRead` for passers); `MappedRead`s are
  materialized only for the reads actually returned.
- **The reader counts everything; the caller selects.** `fetch_mapped_reads`
  returns a *complete* `FilterCounts` (all categories, identical to today's
  reader), reusing the shared `classify_segment_record`
  ([segment_reader.rs](../../../src/bam/segment_reader.rs)) so categorization —
  including the precedence when a read matches several drop reasons — is identical
  *by construction*. `qc_counts` then uses only the four quality categories it
  cares about; that selection is the caller's job, not the reader's.
- **Serve by copy, never move.** Adjacent loci windows can overlap (window =
  tract ± flank; the catalog only guarantees ~flank separation), so one read can
  belong to two consecutive loci. Returned passers are *copied* out; the slice
  stays in the cache until FIFO eviction. "Serve" (copy) and "evict" (drop on
  cache overflow) are decoupled.
- **FIFO cache, cap = `max_cached_slices` (default 3).** A window overlaps ≤ 2
  slices (1 normally, 2 at a boundary), +1 for forward jitter. Because loci are
  processed in sorted order, FIFO ≈ LRU (the oldest slice is the furthest behind).
  **The cap is purely a perf knob:** a single `fetch` whose window overlaps more
  slices than the cap still works — it decodes, extracts, and lets older slices
  evict mid-call; the cap only governs *cross-call* reuse, never within-call
  correctness.

### Byte-identity (the gate)

A direct per-locus query for window `W` returns reads overlapping `W` from the
slices overlapping `W`, classified by `classify_segment_record`, with the drop
counts for `W`. `CramReader::fetch_mapped_reads(W)` iterates the *same* slices
(same `.crai` selection — §4), the same records, the same classifier, counts the
same drops, and copies out the same passers in the same order. The reservoir's
determinism depends only on reads arriving in `(ref_id, pos)` → file → record
order (arch §8.4), which the slice-ordered, record-ordered iteration reproduces.
∴ `(reads, FilterCounts)` are byte-identical to the per-locus path.

This is the **gate test**: for a real catalog + CRAM, assert
`CramReader::fetch_mapped_reads(L.window) == get_reads_from_segment(L.window)` for
every locus, and assert the end-to-end `.ssr.psp` is identical to the current
path's.

## 4. Slice selection must match noodles

`slices_overlapping(chrom, start, end)` must return exactly the slices noodles'
own `query()` would touch — the trap is boundary reads: a read overlapping
`[start,end]` can *start* in an earlier slice (long read / left overhang) and is
physically stored there. The `.crai` slice span covers its reads' full extents,
so an overlap test on the span includes that slice — but to avoid drifting from
noodles, the implementation reuses noodles' index-query slice selection (list the
overlapping slice offsets) rather than hand-rolling overlap math, then decodes
each selected slice (with caching) by offset.

*(Implementation note: confirm noodles-cram 0.93 exposes (a) listing the slice/
container offsets overlapping a region from the `crai::Index`, and (b) decoding a
slice/container by offset — the building blocks for interposing the cache. If a
clean seam isn't public, drive the container reader directly per the Query
iterator's internals.)*

## 5. `covered_regions()` + orchestrator slice-grouping (step 2, optional)

The `.crai` index is a flat list of slice entries carrying
`(reference_sequence_id, alignment_start, alignment_span, …)` — i.e. every
slice's genomic span, already loaded with the index.

```text
CramReader::covered_regions() -> Vec<(chrom, start, end)>   // None/empty for BAM
```

The orchestrator uses it to **place** loci, not for correctness:

- Assign each locus to a "home" slice — the last slice whose start ≤ the locus's
  window start (slices can overlap slightly, so a single rule makes the partition
  well-defined). Consecutive loci sharing a home slice form a **group**.
- `par_iter` over **groups** (an indivisible item per group). Each group runs on
  one worker, whose per-thread cache then decodes the group's slice once and
  serves all its loci from cache. Different groups → different slices → different
  threads → parallel decode of distinct slices.

This shrinks the residual cross-thread re-decode from "one boundary slice per
worker-pair" (what plain count-based `par_iter` splitting causes when a split
lands mid-slice) to just the handful of loci whose window literally straddles a
slice edge (assigned to their home group; processing them touches the neighbor
slice — an accepted small extra decode).

**It is a correctness-independent hint.** The cache (§3) gives byte-identical
output under *any* distribution; grouping only changes how many times slices are
decoded. So it ships *after* the cache and is A/B-measurable on its own.

**Trade to watch:** slice-aligned groups are *variable-sized* tasks, vs the even
ranges count-based splitting gives. A few deep slices → heavy groups → a
load-imbalance tail; rayon work-stealing + the per-locus reservoir cap mitigate
it. For a decode-bound workload, eliminating decode redundancy should beat
perfectly even load — but *measure* rather than assume (the incremental build
makes this a direct A/B).

## 6. Format handling

The optimization targets **CRAM** (heavy reference-differential decode + the
hardcoded per-slice MD5). `CramReader` and the cache are CRAM-specific.

- **CRAM:** `CramReader` (§3); slice-grouping (§5) when desired.
- **BAM:** `.csi`/`.bai` are bin-based (no slice concept) and BAM decode is far
  lighter (no reference-differential, no per-slice MD5). `covered_regions()` is
  empty, so the orchestrator falls back to count-based splitting, and the BAM
  path keeps today's per-segment reader (or a fixed-span buffer variant later if
  a profile justifies it). BAM is not where the measured cost is.

## 7. Alternatives considered (and rejected)

### 7.1 One shared, thread-safe `CramReader` (locked cache)

*Idea:* a single shared reader so each slice is decoded once **globally**.
*Rejected* because the benefit is small for this workload and the cost is real:

- **Tiny dedup opportunity.** Work is sorted + range-partitioned, so threads sit
  in *different* genome regions and rarely want the same slice — only at the seam
  between adjacent workers' ranges (~one boundary slice per worker-pair, a few %
  of decodes). Per-thread caches already capture the within-thread redundancy,
  which is the bulk.
- **The decode-under-lock dilemma.** Hold the lock during decompression and you
  serialize the bottleneck (worse than today). Don't, and two threads that miss
  the same slice both decode it (defeating the purpose) unless you add a
  single-flight/cache-stampede protocol (per-slice in-progress latches + waiters)
  — real, subtle concurrency code with its own contention.
- **No memory saving.** To avoid threads-at-different-regions evicting each other,
  the shared cache must hold ~`2 × threads` slices — the same total as per-thread
  caches.
- Plus Arc'd slices for safe copy-out and a shared-handle pool for parallel I/O.

Slice-grouping (§5) achieves "each slice decoded ~once, in parallel" by
*placement* instead of *locking* — same benefit, no shared mutable state.

### 7.2 Single fetch producer thread + parallel realignment

*Rejected:* decode **is** the bottleneck; a single producer serializes it and
caps throughput at one-thread decode speed. Decode must stay parallel.

### 7.3 Slice-grouping as a correctness requirement

*Rejected in favor of the cache:* making the per-locus QC counts depend on
explicit grouping (one union query per group, demux) breaks the per-locus
filter-drop attribution (one aggregate count for the group, unsplittable). The
cache keeps the per-locus query interface intact, so QC stays byte-identical, and
grouping becomes a pure placement optimization on top.

## 8. Build order

1. **`CramReader`** (§3): per-thread reader + FIFO slice cache, `fetch_mapped_reads`
   returning `(Vec<MappedRead>, FilterCounts)` via the shared classifier, cache of
   decoded records. Wire it into the SSR fetch (`fetch_locus_reads`) per worker
   via `map_init`, keeping the current `AlignmentFile` path behind a flag as the
   byte-identity oracle. **Gate test** (§3). This alone kills the within-thread
   redundancy — the bulk of the win — with no change to the driver's parallel
   structure.
2. **`covered_regions()` + slice-grouping** (§5): the orchestrator groups loci by
   home slice and `par_iter`s over groups. A/B vs step 1.
3. **Measure:** slice-decode count vs loci-with-reads (redundancy removed),
   re-profile (CRAM-decode + MD5 self-time drop), end-to-end wall on the fixture
   and a denser run.

## 9. Open questions

- **Group granularity (step 2):** one slice per group, or merge small adjacent
  slices to cut `par_iter` scheduling overhead?
- **Cache lifetime:** per-job (re-created by `map_init` on each rayon job, like
  `LocusScratch` today) vs persisted per-worker across jobs/batches. Per-job is
  simpler; persisting raises hit rate at batch seams. Start per-job.
- **noodles seam:** whether 0.93 exposes a clean "list overlapping slice offsets"
  + "decode slice by offset", or we drive the container reader internals (§4).
- **Where it lives:** `CramReader` in `ssr/pileup` first (SSR-specific payoff;
  the SNP `--regions` path already amortizes decode over a region). Promote to
  `bam/segment_reader` only if the SNP path later wants it.
