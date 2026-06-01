# `var-calling-from-bam` — unify on the chunk architecture

Implementation plan for re-pointing the **direct path** (`var-calling-from-bam`,
CRAM/BAM → VCF, single sample, no intermediate `.psp`) at the same
chunk-based columnar architecture the **`.psp` path** (`var-calling`) now
uses. Today the two subcommands run on two different engines:

- `var-calling` (`.psp` in) → `src/var_calling/from_psp/` —
  within-chromosome chunk loop, columnar, parallel windows
  ([`drive_cohort_chunked`](../../../src/var_calling/from_psp/driver.rs)).
- `var-calling-from-bam` (BAM/CRAM in) → `src/var_calling/from_bam/` —
  the **old streaming per-position pipeline**
  ([`drive_cohort_pipeline`](../../../src/var_calling/from_bam/pipeline.rs)),
  one rayon worker per chromosome, whole-chromosome record
  materialisation.

The goal is **one architecture, with as little direct-path-specific code
as possible.** The direct path exists chiefly as a like-for-like
comparison point against FreeBayes (single sample, no `.psp` round-trip);
in normal use the caller goes through `.psp`. So the bar here is
*architectural unity and correctness*, not a memory or wall-time headline.

## Spec / supporting documents

- Pipeline architecture spec:
  [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md).
- The chunk rewrite this plan piggybacks on:
  [cohort_within_chromosome_parallel.md](./cohort_within_chromosome_parallel.md)
  (+ the streaming-columnar produce and DUST-worker-pool follow-ups).
- The current direct-path structure being replaced:
  [var_calling_from_bam_per_chromosome.md](./var_calling_from_bam_per_chromosome.md).
- Module reorg context: `from_psp/` + `from_bam/` split, shared stages
  stay top-level in `src/var_calling/`.

## Design decisions settled before drafting

These were agreed in a design discussion and are the premises of the plan:

1. **Converge at the columnar chunk, not at row-shape records.** The
   production `.psp` path never materialises row-shape `PileupRecord`s —
   it decodes PSP columns straight into `SampleColumns` (the memory
   headline of the streaming-columnar rewrite). Forcing it back through a
   row-record seam to look symmetric with the direct path would regress
   the cohort path. So the shared boundary is a **materialised columnar
   block**, and only the direct path goes record → column.
2. **The direct path is always single-sample.** Multi-sample work goes
   through `.psp`. `n_samples == 1` collapses the loader's cohort
   machinery (cross-sample position union, the `min` watermark) to a
   degenerate case — this is *why* "little extra code" is achievable.
3. **The walker stream is the in-memory equivalent of the un-filtered
   PSP record stream.** A `.psp` file is a serialised Stage-1 walker run.
   So `BAM → walker → records` yields the same records `BAM → walker →
   .psp → read back` yields, pre-filter — which is what lets the direct
   path reuse the `.psp` path's variant filter verbatim.
4. **Memory is explicitly a non-goal for the direct path.** Holding a
   chunk (or even a whole chromosome) of single-sample records resident
   is fine. The only constraint is *not* paying that cost on the cohort
   path.

## The convergence seam: `ReadyBlock`

The `.psp` driver already factors cleanly into a **producer** and a
**source-agnostic consumer**, joined by one owned value:

```
ReadyBlock {
    seq_idx,                 // genomic production order
    chunk:    MaterialisedChunk,   // the loaded columns for [range.start, safe_end)
    partition: WindowPartition,    // the variant-group partition
    pre_fetched_ref_bytes: Vec<Vec<u8>>,  // one REF buffer per group
}
```

- **Producer** — [`BlockIterator`](../../../src/var_calling/from_psp/driver.rs)
  (`next_block` / `produce_block`). PSP-coupled: it owns
  `ColumnSpanReader` sources, discovers covered intervals from the PSP
  **block index**, drives the `StreamingBlockLoader` watermark fold,
  runs the `DustAheadPool`, partitions, and prefetches REF bytes.
- **Consumer** — [`drive_blocks_serial`](../../../src/var_calling/from_psp/driver.rs#L471)
  / `drive_blocks_parallel` → [`process_block`](../../../src/var_calling/from_psp/driver.rs#L1412)
  (→ `run_window`) → [`emit_or_drop`](../../../src/var_calling/from_psp/driver.rs#L753).
  **Already source-agnostic** — it only ever touches a `ReadyBlock`, the
  per-stage configs, the writer, and the counters.

So the unification reduces to: **give the direct path its own producer of
`ReadyBlock`s, and feed the existing consumer.** Everything from
`ReadyBlock` onward — partition math, per-group merger, posterior EM,
the post-EM downstream filters, VCF emit, the parallel block scheduler —
is shared unchanged.

The direct-path producer is *simpler* than the PSP one, not more complex,
because (per decision 2) it has no cross-sample fold and no block index:

| Producer step | PSP producer | Direct producer |
|---|---|---|
| Source | `ColumnSpanReader` over seekable, block-indexed `.psp` | walker `Result<PileupRecord, _>` stream, single sample |
| Interval discovery | union block indices → covered intervals | whole chromosome is one interval `[1, len]` |
| Span fold / filter | `StreamingBlockLoader::fill_block` (N-source watermark) | `load_chunk_from_iters` over the one walker iterator |
| Boundary safety | carryover + partition safe-gap | **identical** — carryover + partition safe-gap |
| DUST mask | `DustAheadPool` (parallel ahead-of-time) | `sdust_mask_for_span` inline (single sample, low volume) |
| Partition + REF prefetch | `partition_window` + `prefetch_window_ref_bytes` | **identical** |

Note the direct producer reuses `load_chunk_from_iters`
([loader.rs:263](../../../src/var_calling/from_psp/loader.rs#L263)), which
is *already retained in-tree as the streaming loader's
equivalence-test oracle* — there is an existing test contract asserting it
produces the same chunk `StreamingBlockLoader::fill_block` does. So the
direct path's variant filter is guaranteed identical to the cohort path's
**by a contract that is already maintained**, with no new filter code.

## Architecture after unification

```
              ┌─────────────────────────────────────────────┐
  .psp in ───▶│ PSP producer (BlockIterator)                 │──┐
              │  ColumnSpanReader · covered intervals ·      │  │
              │  StreamingBlockLoader · DustAheadPool        │  │ ReadyBlock
              └─────────────────────────────────────────────┘  │
                                                                ▼
              ┌─────────────────────────────────────────────┐  shared consumer
  BAM/CRAM ──▶│ direct producer (WalkerBlockProducer)        │  drive_blocks_{serial,parallel}
              │  Stage 1 (reader→BAQ→walker) ·               │──┘   → process_block (run_window)
              │  load_chunk_from_iters · sdust_mask ·        │      → emit_or_drop → CohortVcfWriter
              │  partition_window · prefetch_ref_bytes       │
              └─────────────────────────────────────────────┘
```

## Phases

The phasing mirrors the cohort rewrite: a pure refactor first
(byte-identical by construction), then the new producer (byte-identical
vs the old direct path at T=1), then the parallelism / retirement
decisions.

### Phase A — extract the source-agnostic block consumer (pure refactor)

No behaviour change; `var-calling` output byte-identical before/after.

1. Define a `BlockProducer` trait capturing what the consumer needs from
   a producer:
   ```rust
   trait BlockProducer {
       fn next_block(&mut self) -> Option<Result<ReadyBlock, ChunkDriverError>>;
       fn recycle(&mut self, block: ReadyBlock);
       fn chunk_counters(&self) -> (u64, u64); // chunks_loaded, chunk_variants_total
   }
   ```
   `BlockIterator` implements it (it already has these methods/fields).
2. Make `drive_blocks_serial` generic over `P: BlockProducer` instead of
   `BlockIterator<R>` directly. `drive_blocks_parallel` stays
   `BlockIterator`-specific **for now** (its producer thread, free-list
   recycling, and `Read + Seek + Send` bound are PSP-shaped; the direct
   path will start serial — see Phase C).
3. **Placement (decided): common machinery moves up to `src/var_calling/`,
   non-common pieces stay in `from_psp/` and `from_bam/`.** The
   source-agnostic pieces (`ReadyBlock`, `process_block`, `emit_or_drop`,
   `record_fails_mapq_diff_t`, `roll_window_stats`, the `BlockProducer`
   trait, `drive_blocks_serial`, `ChunkDriverParams` /
   `ChunkDriverStats` / `ChunkDriverError`) plus the already-shared
   columnar infrastructure (`columns`, `loader`, `partition`, `worker`,
   `kernels`) hoist out of `from_psp/` into a neutral top-level module
   under `src/var_calling/` (proposed: `src/var_calling/chunk_pipeline/`,
   a sibling of `from_psp/` and `from_bam/` — keeping the multi-file core
   grouped, alongside the already-flat shared stages `dust_filter.rs`,
   `per_group_merger.rs`, … in `src/var_calling/`). After the move:
   - `from_psp/` keeps only the PSP producer: `ColumnSpanReader`,
     `covered_intervals_for_chrom`, `DustAheadPool`, `BlockIterator`,
     `drive_cohort_chunked` + the `var-calling` CLI glue.
   - `from_bam/` keeps only the walker producer (Phase B) + the
     `var-calling-from-bam` CLI glue.
   Done as a **pure move** (no logic edits) so the refactor stays
   trivially reviewable and byte-identity is obvious. (Flat-files-vs-
   `chunk_pipeline/`-subdir is a cosmetic call settled at implementation
   time; either satisfies "common machinery in `src/var_calling/`".)
4. Validation: `var-calling` byte-identical (drop `^##`, md5) on the
   3-tomato fixture at serial + 8 threads; full lib test suite green.

### Phase B — direct-path walker producer (`WalkerBlockProducer`), serial, byte-identical

This is the bulk of the new code, and it is small.

1. **Stage 1 → record stream.** Reuse the existing per-chromosome Stage-1
   wiring (`AlignmentMergedReader::query` → BAQ → walker → error-shedding
   adapter), which `process_one_chromosome_from_bam` already builds
   ([from_bam/driver.rs](../../../src/var_calling/from_bam/driver.rs#L679)).
   It currently `collect()`s into `Vec<PileupRecord>`; keep that for
   Phase B (whole-chromosome materialisation, single sample — fine per
   decision 4). The `Vec<PileupRecord>::into_iter().map(Ok)` shape
   satisfies `load_chunk_from_iters`'s `IntoIterator<Item = Result<…>>`
   + `Send` bound (the `VecColumnSource` test fixture proves this works).
2. **`WalkerBlockProducer: BlockProducer` — one block per chromosome.**

   *Design correction (2026-06-01, found while coding Phase A):*
   `load_chunk_from_iters` is the **batch** loader — it *consumes* a
   carryover prefix but never *emits* one, and it compacts the whole
   requested range in a single call. Group-aware cutting + carry-forward
   of an open group across blocks is exclusively `StreamingBlockLoader`'s
   job, and that needs a `SpanColumnSource` (random-access peek), which a
   forward-only walker stream does not provide. So the batch loader
   **cannot** do chunk-at-a-time with carryover; the natural unit is one
   `load_chunk_from_iters` call over the **whole chromosome**, yielding
   **one `ReadyBlock` per chromosome**.

   This is byte-identical to the `.psp` path's many-small-blocks output:
   `partition_window` cuts the block into variant groups at
   `max_group_span` safe gaps, and a group is determined by variant
   positions + gaps, *independent of block boundaries*. The `.psp` path's
   safe-gap cutting guarantees it never splits a group across blocks, so
   both paths see the same natural groups → same `run_window` output. (The
   in-tree `var_calling_byte_identical_across_target_variants_per_chunk`
   test already proves the `.psp` path's block-boundary independence.) It
   also matches today's direct-path memory profile (whole-chromosome
   materialisation) — fine per decision 4.

   `WalkerBlockProducer` holds: the chromosome list + a cursor over it,
   a way to obtain each chromosome's single-sample records (Stage 1 —
   B.2), the per-chrom `StreamingChromRefFetcher`, a `ChunkLoadScratch`,
   a `PartitionScratch`, a DUST-mask buffer, and the per-stage configs.
   `next_block` advances to the next data-bearing chromosome and:
   1. `load_chunk_from_iters` once over `[1, chrom_length]` (single
      sample iterator; empty carryover) → the chromosome's `MaterialisedChunk`,
      applying the **shared** variant filter.
   2. `sdust_mask_for_span` over the block's `[range.start, safe_end)`
      span (unless `--no-complexity-filter`). The per-position DUST
      verdict is span-independent past the `MIN_DUST_HALO` barrier, so a
      single span call matches the `.psp` path's per-interval mask for the
      same positions. (For very long contigs, sub-span chunking like
      `dust_mask_for_interval` bounds the REF buffer — a memory refinement,
      not a correctness one; deferrable.)
   3. `partition_window` → `WindowPartition`; `prefetch_window_ref_bytes`
      per group.
   4. Assemble + return the `ReadyBlock`.
   `recycle` pushes spent buffers onto a free-list; `chunks_loaded` /
   `chunk_variants_total` return the running totals (1 chunk/chrom).
   - *Phase C note:* true within-chromosome chunking (multiple blocks per
     chromosome, the prerequisite for within-chromosome parallelism) needs
     a **walker-backed `SpanColumnSource`** feeding `StreamingBlockLoader`
     — a forward buffer with `peek_next_span`/`read_span` — not the batch
     loader. That's deferred with the parallelism work.
3. **Driver.** New `drive_from_bam_chunked` (sibling of
   `drive_cohort_chunked`): per chromosome, build a `WalkerBlockProducer`,
   call the shared `drive_blocks_serial` into one process-wide
   `CohortVcfWriter`. **No per-chrom fragments, no concat** — the chunk
   driver writes the whole VCF through one writer, exactly like
   `drive_cohort_chunked`. (This retires the
   `TempDir`/fragment/`concat_fragments` machinery for this path.)
4. **CLI.** `run_var_calling_from_bam` keeps steps 0–6 (index pre-flight,
   rayon pool, config build, contig harvest, header/index load, metadata)
   and swaps the step-9 `rayon::par_iter` over
   `process_one_chromosome_from_bam` for a call to
   `drive_from_bam_chunked`. `ChunkDriverStats` already mirrors
   `CohortDriveStats`, so the run-summary code is reused.
5. **Validation (the hard gate).** Two checks, serial:
   - **One-time regression guard:** before the old streaming driver is
     deleted (Phase D), capture a checked-in **golden VCF** from the
     *current* `var-calling-from-bam` on the from-bam integration
     fixtures, and assert the new path reproduces it byte-for-byte. This
     is the only window in which the old path can speak, so capture it
     here.
   - **Cross-path differential (the maintained oracle):** assert
     `var-calling-from-bam X.bam` == `pileup X.bam → X.psp; var-calling
     X.psp` byte-for-byte. Since the old direct path is being retired
     (decision below), the two **new-architecture** paths become each
     other's oracle going forward. This is contract (2) in
     "Byte-identity strategy" and depends on the PSP-losslessness spike.

### Phase C — within-chromosome parallelism (target; deferred)

**Decided direction: move to within-chromosome parallelism** (matching the
`.psp` path), but **deferred** — Phases A/B land serial first, and the
parallelism work is a later phase taken with a wall-time measurement in
hand.

The old direct path parallelised **across chromosomes** (one rayon worker
per contig). The shared chunk consumer parallelises **within a
chromosome** (`drive_blocks_parallel`, parallel windows) — the better
model for a single sample on a genome with one dominant chromosome. The
work is to generalise `drive_blocks_parallel` over `BlockProducer + Send`
(its producer-thread + free-list shape is currently `BlockIterator`-
specific) so the walker producer feeds it. Byte-identity must hold vs the
serial path — `drive_blocks_parallel` already proves seq-idx-ordered emit
is identical to serial, so this is a soundness/`Send` exercise, not an
output-changing one.

### Phase D — retire the old direct path / cleanup

**Decided: delete the production-dead code** (the saving is small in LoC,
but it leaves exactly one architecture — the point of this work). Once
Phase B's golden + cross-path differential are green:

1. **Delete `from_bam/pipeline.rs`'s `drive_cohort_pipeline`** (the old
   streaming driver) and `process_one_chromosome` if it has no remaining
   caller. The Phase-B golden VCF (checked in) replaces its
   regression-guard role; the maintained oracle is the cross-path
   differential between the two new-architecture paths.
2. Remove now-dead direct-path machinery: `process_one_chromosome_from_bam`,
   the fragment-path/`TempDir`/`concat_fragments` wiring in
   `run_var_calling_from_bam`, `all_chromosomes_to_contig_list` if unused.
3. Audit `from_bam/pipeline.rs`'s public surface (`CohortPipelineParams`,
   `CohortDriveStats`, the `DEFAULT_*` / `MAPQ_*` constants): keep only
   what the new path or the CLI still references; the rest goes with the
   driver.
4. `from_bam/mod.rs` doc + `PROJECT_STATUS.md` updates.

## Byte-identity strategy

The old direct path is being **retired**, so it cannot be a permanent
oracle. The validation therefore rests on:

1. **One-time golden, captured from the old path before deletion**
   (Phase B regression guard). Same caller, same input, same VCF body.
   The strongest guarantee is maximal reuse — only the producer differs,
   and even the producer reuses the filter via `load_chunk_from_iters`.
   This is a checked-in fixture, not a live code path.
2. **Cross-path differential — the maintained oracle.**
   `var-calling-from-bam X.bam` == `pileup X.bam → X.psp; var-calling
   X.psp`, byte-for-byte. With both callers on the new architecture and
   only their producers differing, this is the durable correctness check
   after Phase D.

**Prerequisite spike (gates contract 2) — DONE, PASSED (2026-06-01).**
Contract (2) holds **iff the PSP writer/reader round-trip is lossless for
every field the chunk pipeline consumes** — `pos`, alleles, `num_obs`,
`q_sum`, `fwd`, `placed_left`, `placed_start`, `mapq_sum`, `mapq_sum_sq`,
`ref_span` (= `alleles[0].seq.len()`), chain ids. A throwaway test built
an adversarial `Vec<PileupRecord>` (all seven scalars populated, chain
ids, a `num_obs == 0` bucket, deletion + insertion alleles, awkward
non-decimal `q_sum` f64s incl. `-1e-300`/`-0.0`, two chromosomes, a 64-byte
block target forcing multiple flushes so chain ids + position deltas cross
block boundaries), wrote it via `PspWriter`, read it back via
`PspReader::records()`, and asserted full-`PileupRecord` `PartialEq`
equality (`assert_eq!(got, want)`) — **passed**. Root cause it can't
regress: `q_sum`/`mapq_sum_sq` etc. are stored as raw fixed-width LE bytes
(`WireScalar`, `block.rs`), bit-exact by construction; `f64` already has
`to_bits()` round-trip tests. So the cross-path differential is a sound
oracle and no codec work is needed. (The spike was throwaway; a permanent
full-record round-trip test should be added as part of the Phase B test
scaffolding.)

New tests:
- A `BlockProducer`-level differential test: same single-sample records
  fed through (a) `WalkerBlockProducer` and (b) a `.psp` written from the
  same records → `BlockIterator`; assert identical `ReadyBlock` sequence
  (chunk columns, partition, REF bytes) — pins the seam directly.
- An integration test running both `drive_from_bam_chunked` and the old
  `drive_cohort_pipeline` on a fixture with ≥1 MNP, ≥1 deletion straddling
  a chunk boundary, ≥1 LH-cap site, ≥1 hom-REF group, ≥1 below-`min_alt_obs`
  site, ≥1 below-`qual_phred` site, ≥1 MAPQ-diff-drop site; assert VCF
  bodies byte-equal and counter sets field-equal.
- A deletion-straddling-the-chunk-boundary regression (exercises the
  carryover path with `chunk_genomic_span` set small enough to force a
  mid-group cut).

## Risks / things to watch

- **Forward-only source vs span re-slicing.** Phase B re-slices a
  materialised `Vec` per chunk, which sidesteps the forward-only problem
  by materialising up front. True streaming (Phase C/later) needs a
  persistent cursor; the `load_chunk_from_iters`-consumes-by-value shape
  is the thing to revisit there.
- **DUST equivalence.** The chunk path slices a per-interval DUST mask;
  the old streaming direct path runs `DustFilter` inline. Both call
  `sdust_mask_for_span` under the hood, but the *halo/boundary* handling
  (`MIN_DUST_HALO`, the W-barrier rule) must match — this is the known
  cross-path DUST-unification deferred item. Verify on a low-complexity
  fixture, not just random sequence.
- **Empty chromosomes.** A contig with zero walker records must produce
  zero blocks and contribute only its header presence — `next_block`
  returning `None` immediately for that chrom, exactly as the PSP
  producer handles an empty block index.
- **`min_alt_obs_per_sample` is pre-EM.** It runs inside `run_window`
  (merger→EM boundary), not in `emit_or_drop`. Already shared, but note
  it when reconciling counters against the old path.

## Decisions (resolved) and remaining spike

The four design questions are **resolved**:

1. **Module placement** — common machinery moves up to `src/var_calling/`
   (proposed `src/var_calling/chunk_pipeline/`); only the PSP producer
   stays in `from_psp/` and only the walker producer in `from_bam/`.
2. **Byte-identity** — the old direct path is retired, so the maintained
   oracle is the **cross-path differential** between the two new-
   architecture paths (+ a one-time golden captured from the old path
   before deletion).
3. **Old streaming driver** — **deleted** in Phase D (the saving is
   small, but it leaves one architecture).
4. **Parallelism** — target **within-chromosome** (matching the `.psp`
   path), **deferred** to a later phase after a serial Phase B lands.

**The prerequisite spike is done and passed (2026-06-01)** — the PSP
round-trip is lossless for every field the chunk pipeline consumes (see
"Byte-identity strategy"), so the cross-path differential oracle is sound
and there are **no open items** blocking Phase A. A permanent full-record
round-trip test folds into the Phase B test scaffolding.

## Out of scope

- Multi-sample direct calling (stays `.psp`-only by design).
- Any change to the `.psp` path's behaviour or output.
- Memory optimisation of the direct path (explicit non-goal — decision 4).
- `estimate-contamination` (the `.psp`-only contamination path is
  unaffected).
