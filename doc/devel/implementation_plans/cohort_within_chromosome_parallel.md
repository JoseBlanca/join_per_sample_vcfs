# Cohort `var-calling` — within-chromosome parallelism via chunked materialisation

Implementation plan for the architectural follow-up to the
[2026-05-27 scaling measurement](../reports/reviews/scaling_measurement_2026-05-27.md):
break the per-chromosome parallelism ceiling and decouple memory from
worker count by switching from streaming iterators to
chunk-at-a-time random-access materialisation, with workers
parallelised over windows of a chunk.

Today's reality, post-H1 per-chromosome parallelism and post-fcef495
512 KiB block defaults, measured at T=4 on a 2 Mbp tomato fixture
replicated synthetically:

| N (cohort) | wall (s) | peak RSS | dominant CPU buckets (inclusive) |
|---:|---:|---:|---|
| 50   | 12.8  | 1.3 GB  | per_group_merger 67 %, dust_filter 62 %, per_position_merger 16 % |
| 200  | 33.9  | 3.8 GB  | per_group_merger 64 %, dust_filter 56 %, per_position_merger 39 %, allocator 30 % |
| 1000 | 185.5 | 16.7 GB | per_group_merger 58 %, dust_filter 48 %, per_position_merger 45 %, psp_reader 44 %, allocator 38 %, posterior_engine 11 % |

Two structural problems that this rewrite addresses:

1. **Memory ceiling.** Peak RSS grows linearly at ~16 MB/sample on
   synthetic data (~8 MB/sample on real). At N=5000 that's 40–80 GB
   on a single host — fragile or impossible. The dominant term is
   `T_chrom × N × decoded_block_size` because every per-chrom worker
   holds N PSP readers alive for the duration of its chromosome.
2. **Per-chromosome parallelism ceiling.** H1 (May 2026) realised
   3.85× wall reduction at T=13 on the multi-chrom tomato fixture, but
   the cliff is hard at T=13 because there are 13 chromosomes. ch00
   carries ~13× the per-chrom record load on tomato, so even at T=13
   the run wall is bounded by ch00's serial walk. Adding more workers
   does nothing.

Switching to chunk-based within-chromosome parallelism solves both:
parallelism scales with worker count rather than chromosome count
(no more ch00 cliff), and the memory footprint is decided by chunk
size × N rather than `T_chrom × N`, so memory becomes orthogonal to
both parallelism and chromosome count.

## Spec / supporting documents

- **Scaling measurement (the motivating data):**
  [scaling_measurement_2026-05-27.md](../reports/reviews/scaling_measurement_2026-05-27.md).
  Memory linearity, per-stage CPU shares, per-stage heap shares
  across N=50/200/1000 establish the architectural cliff this plan
  addresses.
- **Per-chromosome parallelism (the predecessor):**
  [cohort_per_chromosome_parallel.md](cohort_per_chromosome_parallel.md)
  + [cohort_per_chromosome_parallel_2026-05-20.md](../reports/implementations/cohort_per_chromosome_parallel_2026-05-20.md).
  This plan replaces H1's per-chrom rayon decomposition with a
  within-chrom decomposition; the per-chrom layer goes away. H1's
  contig-table ordering invariant survives (chromosomes are
  processed sequentially in contig order). The bgzf-aware concat
  module (`src/var_calling/vcf_writer/concat.rs`) is *not* used by
  the new var-calling path — assembly of worker outputs is a
  var_calling-internal iterator chain, not a vcf_writer concern.
  The concat module stays in place because `var-calling-from-bam`
  still uses it for its per-chrom fragments.
- **Pipeline architecture spec:**
  [calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md).
  §"Stage 3 — low-complexity filter" through §"Stage 6 — posterior
  engine" remain unchanged at the algorithm level; what changes is
  the *driver* — how records get fed through them.
- **PSP reader random-access primitive:**
  [src/per_sample_pileup/psp/reader.rs](../../src/per_sample_pileup/psp/reader.rs)
  `region_records(chrom_id, start, end)` already exists and returns
  an iterator over records in the genomic range, using the per-block
  offsets table for seek. This is the load-bearing primitive that
  makes chunk-based materialisation cheap.
- **Per-position merger (gets replaced):**
  [src/var_calling/per_position_merger.rs](../../src/var_calling/per_position_merger.rs).
  The k-way merge over N PSP iterators becomes the chunk loader's
  responsibility; the streaming merger as currently written is no
  longer needed once the rewrite lands.
- **Bench harness for the A/B comparison:**
  [benchmarks/tomato1/scripts/perf_scaling_synthetic.py](../../benchmarks/tomato1/scripts/perf_scaling_synthetic.py)
  (sweep + samply + dhat) and Pair D of
  [benchmarks/tomato1/scripts/perf_dashboard.py](../../benchmarks/tomato1/scripts/perf_dashboard.py).

## Why chunk-based materialisation, not streaming-pipeline parallelism

Two alternatives were considered and rejected on the discussion path
that led here. Recording so the choice is not re-litigated mid-implementation:

- **Producer-consumer over the current streaming pipeline.** One
  per-chrom producer runs the existing `PerPositionMerger +
  VariantGrouper`, feeds completed groups into a bounded queue,
  T worker threads consume groups in parallel. *Rejected because*
  the producer is itself the bottleneck at large N — `PerPositionMerger::next`
  is O(N) per emitted position (a linear scan over N stream heads
  + an N-element allocation), so consumer-side parallelism is
  bounded by the single-threaded producer rate. At N=1000 the
  per_position_merger inclusive share is 45 %, almost the same as
  per_group_merger's 58 %. Parallelising consumers alone leaves
  half the wall untouched.
- **Iterator-level windowed parallelism (workers each open their
  own iterator chain over their window).** *Rejected because* the
  current grouper/merger are pull-iterator-based — they can't look
  back across positions. Workers can't independently scan a window
  starting mid-record-stream without missing context from groups
  that began before the window. Making the iterator chain
  random-access-capable requires materialising records anyway,
  which is what this plan does directly.

The remaining design is: materialise a genomic chunk of records for
all N samples into RAM (variant positions only), run a sequential
pre-pass that partitions the chunk into autonomous windows whose
boundaries fall in `max_group_span`-wide safe gaps, spawn parallel
workers (one per window). Each worker runs **the full pipeline**
on its window — per_position_merger + variant_grouper +
per_group_merger + posterior_engine (EM) — and emits its final
variant records into a per-window output buffer. The driver streams
window outputs to `vcf_writer` in order (window 0 of chunk 0 first,
then window 1, …) as soon as each window's records become available
— no per-chromosome accumulation buffer, no global EM pass at the
end. The EM lives inside the worker because each variant group's EM
is independent of every other group's (priors are per-group, not
cross-group).

## Scope

**In scope:**

- New chunk-based driver in `src/var_calling/cohort_block.rs`
  replacing the per-chrom rayon decomposition with chunk-level
  loading + within-chunk worker decomposition. (Distinct file from
  the existing `src/pop_var_caller/cohort_driver.rs`, which stays
  as the CLI-side wiring layer.)
- New `MaterialisedChunk` data structure (columnar per-sample
  storage) holding N samples' records for a genomic range, with
  random-access retrieval by `(sample_idx, position)`.
- **Cohort-wide variant-position filter** in the chunk loader:
  drop every position where no sample carries any non-reference
  allele. Chunks hold only variant positions; non-variant
  positions never reach the worker pipeline. ~30–100× reduction
  in per-chunk record count and memory on real cohorts.
- Sequential **pre-pass** (`fix_boundaries`) after each chunk's
  load + filter: chooses the chunk's `safe_end` and partitions
  `[chunk.range.start, safe_end)` into windows whose boundaries
  all fall in `max_group_span`-wide safe gaps. Records past
  `safe_end` carry over to the next chunk's load.
- Workers process windows in parallel. **Fully autonomous** — by
  construction nothing crosses a window boundary, so each worker
  emits every group it sees. No lookback, no canonical-position
  rule, no inter-worker coordination.
- Pipeline: chunk K+1 loads in background while chunk K is being
  processed. Bounded by a 2-chunk steady state.
- Assemble window outputs into an iterator of records inside
  `var_calling` — chained per-chunk, then per-chromosome — and hand
  *that single iterator* to `vcf_writer`. **vcf_writer has no
  knowledge of chunks, workers, windows, or parallelism**; it
  receives an iterator of post-EM records in genomic order and
  writes one cohort VCF top-to-bottom. The bgzf fragment-concat
  module is not used by this path.
- File-descriptor budget management: at large N the design opens
  many per-sample `region_records` iterators. Plan calls for a
  `RLIMIT_NOFILE` raise at startup; **fail-fast with a clear
  "raise your ulimit to N" error if the raise can't satisfy the
  required budget** — no software fallback.
- Removing the per-chrom-level `rayon::par_iter` from H1; replaced
  by the chunk-loop's parallel windows. H1's contig-table ordering
  invariant survives (chromosomes processed in declared order
  inside `run_var_calling`).

**Out of scope (this plan):**

- **DUST → BED cache.** Cheap independent win; deferred to a separate
  follow-up. The per-chunk worker still runs DUST today; the rewrite
  doesn't depend on or block the cache.
- **SQUAREM acceleration in the posterior engine.** Algorithmic
  change to the EM; orthogonal to the architectural rewrite. Could
  follow Phase D if the new architecture's perf surface motivates it.
- **rayon-over-records in the EM.** The EM stays single-threaded
  in this plan. It's a global pass after the parallel chunk phase
  ends. Folding it into the chunk-parallel phase changes the
  M-step reduction's correctness story (cross-record dependency
  on sample-level priors) and is a separate workstream.
- **PSP reader H1 (CSR collapse of ragged `Vec<Vec<u8>>`).** The
  reader's per-block decoded representation is what the chunk loader
  consumes; if H1 lands first, the chunk loader benefits transparently,
  but neither side blocks the other.

## Design

### Terminology

Three units of work appear in this design — explicit definitions so
the rest of the document reads cleanly:

| Term | Meaning | Memory unit |
|---|---|---|
| **PSP block** | On-disk unit inside `.psp` files (post-`fcef495` default = 512 KiB compressed; decoded into a `DecodedBlock`) | File-format detail |
| **chunk** | The cohort-level unit *loaded into memory* — one genomic range across all N samples. Spans 1–few PSP blocks per sample | Loader / per-iteration memory budget |
| **window** | A sub-range of a chunk processed by **one parallel worker** | Parallelism granularity |

Workers process windows; chunks contain windows; chunks are
materialised by reading from PSP blocks. The boundary between
adjacent chunks (and between adjacent windows inside a chunk) is
not a fixed slice of the genome — it's chosen by the pre-pass
boundary-placement routine described below so that no variant
group can span it.

### Materialised chunk — columnar layout

A chunk covers a genomic range `[G_start, G_end)` on one chromosome,
where `G_start` is the previous chunk's `safe_end` and `G_end` is the
chunk's target end (which the pre-pass may revise down to a `safe_end`
≤ `G_end`).

**The chunk holds only records at variant positions** — the chunk
loader applies a cohort-wide filter immediately after reading and
drops every position where no sample carries any non-reference
allele (see §"Variant-position filter" below). The columnar arrays
in `SampleColumns` are therefore sparse along the genome but dense
where they have entries, with the same position list across all
samples (each sample may have a record at a variant position or
have no record there, but only variant positions are represented).

**Storage is columnar per sample from Phase A.** Two reasons:

1. **It matches PSP's on-disk layout.** PSP blocks are already
   column-oriented (fixed-width per-record scalars + CSR-style
   ragged columns for per-allele data — the writer's
   `encode_list_column_csr`). Loading is a column-by-column copy
   from a decoded PSP block into the chunk's columns; no
   intermediate per-record struct synthesised on the read path.
2. **It keeps the SIMD path open.** The per_group_merger's
   `compute_log_likelihoods` walks N samples × G genotypes; with
   per-sample-then-per-column storage, a future SIMD pass can
   gather four samples' stats at a time without a layout
   refactor. (Whether the SIMD pass actually pays back is a
   Phase D question, but Phase A's job is to not foreclose it.)

```rust
pub struct MaterialisedChunk {
    pub chrom_id: u32,
    pub range: Range<u32>,         // [G_start, G_end) — initial target range
    pub safe_end: u32,             // set by fix_boundaries; <= range.end
    pub windows: Vec<Range<u32>>,  // set by fix_boundaries; partitions [G_start, safe_end)
    pub per_sample: Vec<SampleColumns>,
}

/// Columnar per-sample chunk storage. Records are sorted by position.
/// Per-allele data uses CSR (offsets + flat values) so a record's
/// `k`-th allele lives at index `allele_offsets[i] + k` in the flat
/// per-allele arrays.
pub struct SampleColumns {
    // Per-record (fixed-width):
    pub positions:           Vec<u32>,    // length = n_records
    pub ref_allele_id:       Vec<u8>,     // length = n_records
    // Per-record CSR offsets into the flat per-allele arrays
    // (length = n_records + 1, last entry = total allele count).
    pub allele_offsets:      Vec<u32>,
    // Flat per-(record, allele) fixed-width columns:
    pub allele_num_obs:      Vec<u32>,
    pub allele_q_sum:        Vec<f64>,
    pub chain_anchor_flags:  BitVec,
    pub chain_ids:           Vec<ChainId>,
    // Variable-length allele sequences (nested CSR):
    pub allele_seq_offsets:  Vec<u32>,    // length = total alleles + 1
    pub allele_seq_bytes:    Vec<u8>,
    // (plus whatever else `PileupRecord` carries today)
}
```

Position lookup at a target `pos` is `binary_search` on `positions`;
on a hit, the record's data is gathered via `allele_offsets[i]..allele_offsets[i+1]`
across the flat per-allele arrays. Consumers that today take
`&PileupRecord` get a borrowed-view type — `PileupRecordRef<'a>`
with `&[u8]` / `&[ChainId]` slices — synthesised at access time;
this preserves the existing call shapes through per_position_merger,
variant_grouper, and per_group_merger without forcing them to
become columnar-aware in Phase A. (Phase D may push the columnar
view deeper into the merger algorithms if perf review surfaces it.)

### Chunk loader

```rust
pub fn load_chunk(
    psp_paths: &[Path],
    ref_fetcher: &ChromRefFetcher,
    chrom_id: u32,
    range: Range<u32>,
    carryover: Vec<SampleColumns>,   // already filtered to variant positions
) -> Result<MaterialisedChunk, ChunkLoadError>
```

Three logical steps, one function:

1. **Raw load (per-sample columnar copy).** For each sample, open
   `region_records(chrom_id, range.start, range.end)` and append
   decoded block columns directly into a temporary `SampleColumns`
   — no intermediate `PileupRecord` rows synthesised on the read
   path. The carryover from the previous chunk's pre-pass
   (already variant-position-filtered) is *prepended* to each
   sample's columns before the new records, preserving
   sorted-position order. PSP's existing per-block offsets table
   makes this O(records-in-window) disk-side, with one or two
   blocks per sample typically loaded.
2. **Filter (cohort-wide variant-position selection).** A single
   N-way merge over per-sample position arrays produces the set
   of distinct positions present in any sample's raw columns.
   For each such position, check whether **any** sample's record
   carries any non-reference allele (predicate detailed in
   §"Variant-position filter"). Positions failing the check are
   dropped from the chunk entirely.
3. **Compact.** For each sample, copy records at the surviving
   variant positions into the chunk's final `SampleColumns`.
   Discard the raw load buffers. Steady-state memory after this
   step is the post-filter footprint (~30–100× smaller than the
   raw footprint on real cohorts where variant density is ≪ 1 %).

Peak load memory is ~2× the post-filter footprint, briefly, during
the compact step.

The Phase A loader is single-threaded across samples (sequential
loop over `psp_paths`). Phase B parallelises sample loading inside
`load_chunk` via rayon if it shows up as a bottleneck.

### Variant-position filter

A position `P` is **variant** iff at least one sample's record at
`P` has any allele whose sequence differs from the reference base
at `P`.

Default (safe) implementation: for each sample's record at `P`,
walk its allele list (CSR-indexed); if any allele's `seq_bytes` slice
≠ the single ref-base byte from `ChromRefFetcher`, the record is
variant for this sample, and `P` is variant for the cohort.

```
for each position P in the chunk's position-union timeline:
    ref_base = ref_fetcher.base_at(chrom_id, P)
    is_variant = false
    for each sample s with a record at P:
        for each allele k in record's CSR range:
            if allele_seq_bytes[k_offsets] != ref_base:
                is_variant = true; break
        if is_variant: break
    if not is_variant: drop P
```

If the PSP record schema carries an explicit "this allele is the
reference" marker (e.g. a `ref_allele_index` field per record),
the inner loop simplifies to "is there an allele index ≠
`ref_allele_index` with `num_obs > 0`?" — a fast path worth using
if available. The walker's actual schema decides this at
implementation time; the safe byte-level check is the fallback.

The filter cannot be per-sample-streaming: a sample's record at a
position must be retained whenever **any other** sample carries a
variant at that position, even if the sample itself is homref
there. The merger needs the homref evidence to compute joint
likelihoods. So the filter is post-load, with all N samples'
position data already in memory.

Edge cases:

- **Sample missing at a variant position.** If sample S has no
  record at position P but some other sample T carries variant at
  P, P is kept; sample S's columns simply have no entry at P
  (the merger sees `None` for sample S at P, as today).
- **All samples carry the same alt.** A monomorphic alt (cohort
  is fixed-alt) is still a variant position by this predicate —
  alleles differ from the reference base, so `is_variant = true`.
  This matters for population calling: fixed alts are legitimate
  calls.
- **DUST low-complexity sites.** Currently filtered downstream by
  `dust_filter`. The variant filter doesn't subsume DUST (a
  variant low-complexity site would survive the variant filter
  but be dropped by DUST). DUST stays in the worker pipeline,
  applied to the much smaller variant-only stream.

### Pre-pass: safe-boundary placement

After loading a chunk and before spawning parallel workers, a
sequential pre-pass walks the chunk's records once to:

1. Choose `safe_end ≤ range.end` so no variant group can span the
   chunk's right boundary into the next chunk.
2. Partition `[range.start, safe_end)` into T windows whose
   boundaries also fall in safe gaps, so workers can run on their
   windows fully autonomously (no scan-vs-emit-range distinction,
   no canonical-position rule, no inter-worker coordination).
3. Split off records with `position ≥ safe_end` as carryover for
   the next chunk's loader.

A boundary position `B` is **safe** iff:

- No two records straddling `B` are within `max_group_span` of each
  other (the grouper would join them into one group). Equivalently:
  there's a gap of width > `max_group_span` at `B`.
- No allele with reference span > 1 starts before `B` and ends
  after `B` (a deletion / MNP / complex allele rooted in the left
  side crossing into the right side). For typical
  `max_allele_ref_span ≤ max_group_span` this is implied by the
  gap-width check; we test both for robustness against long
  deletions.

Pseudocode for the pre-pass, columnar-friendly:

```
fix_boundaries(chunk, target_window_count = T):
    # one columnar pass per sample to precompute allele extents
    for each sample s:
        for each record i in s:
            allele_end_max[s][i] = positions[s][i]
                + max_ref_span(alleles_at_record_i)

    # build the *union* of (record-position, allele_end_max) timeline
    # across all samples (an N-way merge in sorted-by-position order
    # of all record positions, with their allele_end_max)
    timeline = merge_sorted(samples)

    # scan timeline backwards from range.end to find chunk safe_end:
    safe_end = largest position <= range.end such that
        - the gap to the next record below is > max_group_span, AND
        - no allele below has allele_end_max > safe_end

    # split records past safe_end into carryover (per-sample slices)
    carryover = take_records_with(position >= safe_end)

    # partition [range.start, safe_end) into T windows by walking
    # the timeline forward and placing T-1 boundaries in safe gaps
    # nearest to evenly-spaced positions
    windows = partition_into_safe_gaps(timeline, T)

    return (safe_end, carryover, windows)
```

The pre-pass is cheap: one columnar scan per sample to compute
`allele_end_max`, one N-way merge over per-sample position arrays
(which are already sorted), a single backward scan for `safe_end`,
and a forward placement of T-1 window boundaries. No record-level
decode, no allele-data touch beyond `ref_span`.

After the pre-pass, the chunk's `windows` field holds T (or fewer)
disjoint ranges that completely tile `[range.start, safe_end)`,
each of which sits inside a `max_group_span`-wide gap on both
sides. Workers process windows in parallel; each worker runs the
**full pipeline** — per_position_merger + variant_grouper +
per_group_merger + posterior_engine (EM) — on its window's records
and emits final variant records (post-EM, VCF-ready) into a
per-window output buffer. No inter-worker coordination, no
canonical-position rule, no per-chrom EM post-pass.

The same mechanism handles both chunk boundaries (with carryover
to the next chunk) and window boundaries (no carryover — adjacent
windows are both in the same chunk and stay in memory). One
routine, one correctness argument, applied at both levels.

#### Edge cases

- **No safe gap exists inside the loaded chunk** (the entire chunk
  is a single dense run of records < `max_group_span` apart).
  Fallback: extend the chunk's load range by another step and
  retry. Hard cap at 4× nominal chunk size — past that, fail loudly
  with a "pathological input" error rather than silently producing
  a multi-GB chunk.
- **Fewer safe gaps than T workers.** Some workers idle on this
  chunk; the next chunk picks up the slack. Not a correctness
  issue, just a transient throughput hit. Chunk-size is the lever
  if this happens often.
- **Chromosome start / end.** No carryover at chromosome start
  (`chunk.range.start = chrom.start`). Chromosome end is always a
  safe boundary by construction (no records beyond it), so the
  last chunk's `safe_end = chrom.end` and carryover is empty.

### Window count vs. chunk size

The chunk size and worker count are independent tuning knobs:

- **Chunk size (`CHUNK_GENOMIC_SPAN`)** controls memory:
  `N × chunk_records × sizeof(record-columns)` per in-flight chunk.
  Larger chunks → fewer chunk-transition overheads, more memory.
  Default target: chunk fits in 1–2 PSP blocks per sample at the
  median record density, so ≈ 100 Kbp on tomato.
- **Worker count (`COHORT_WORKERS`)** controls compute parallelism:
  T workers consume windows of the current chunk in parallel.
  Pre-pass picks ≤ T windows depending on how many safe gaps the
  chunk contains.

Pipeline: while workers process chunk K, the chunk loader is
preparing chunk K+1 in a background thread. Bounded queue of 2
chunks → steady-state memory is `2 × N × chunk_records × sizeof(record-columns)`,
independent of T.

### Grouped-variants buffer — columnar layout

The per_group_merger's output (`PerGroupMerger::process_group`)
is the bridge between the grouper and the posterior engine. In
the current design each emit is a row-shaped record carrying the
group's geometry, unified alleles, per-sample log-likelihoods, and
per-sample scalars / chain-anchor flags. The chunk-parallel
rewrite replaces the iterator emit with a **columnar batch
accumulator** local to each worker, mirroring the `SampleColumns`
shape:

```rust
/// Per-window output of per_group_merger; the input to the
/// worker's per-group EM. Lives only inside one worker — never
/// crosses thread boundaries.
pub struct GroupedVariantsBatch {
    // Per-group (fixed-width):
    pub chrom_ids:        Vec<u32>,
    pub starts:           Vec<u32>,
    pub ends:             Vec<u32>,
    // Per-group CSR into per-allele data (length = n_groups + 1):
    pub allele_offsets:   Vec<u32>,
    // Per-(group, allele) flat columns:
    pub allele_seq_offsets: Vec<u32>,  // nested CSR for sequences
    pub allele_seq_bytes:   Vec<u8>,
    pub allele_is_ref:      BitVec,
    // Per-(group, sample) flat (CSR over groups via allele_offsets * n_samples
    // and a separate genotype CSR for variable n_genotypes per group):
    pub genotype_offsets:   Vec<u32>,
    pub log_likelihoods:    Vec<f64>,  // length = sum(n_samples * n_genotypes per group)
    // Per-(group, sample, allele) for the chain-anchor signals:
    pub chain_anchor_flags: BitVec,
    pub scalars:            Vec<AlleleSupportStats>,
}
```

Why columnar here too:

- **Symmetric with `SampleColumns`** — the worker's input and
  internal-intermediate use the same parallel-array discipline,
  no shape conversion mid-worker.
- **SIMD-ready for the EM.** The posterior engine already does
  SIMD-across-samples per group; with the batch stored
  column-wise, a Phase D pass could go further (SIMD across
  groups for the same genotype shape, batched
  `compute_mixture_log_likelihoods_simd` calls, etc.) without a
  layout refactor.
- **No per-emit allocation.** The merger appends to existing
  column tails; no per-group `Vec::new` for the log-likelihoods
  matrix or the allele sequences.

The batch is consumed by the worker's local posterior engine call
in one shot (or in genotype-shape-aligned sub-batches). Output is
final variant records pushed into the worker's output buffer.

### Worker pipeline

Each worker, given a window `[W_start, W_end)` of a
`MaterialisedChunk`:

1. Runs per_position_merger + variant_grouper over the window's
   columnar slices (using `PileupRecordRef<'a>` borrowed views).
2. Each emitted group is appended to a local `GroupedVariantsBatch`
   via per_group_merger's `process_group`.
3. The local `GroupedVariantsBatch` is fed to the per-group
   posterior engine (EM), which writes final variant records
   into the worker's per-window output buffer (`Vec<FinalVariantRecord>`
   or similar — concrete shape decided at implementation time).
4. Worker returns its output buffer to the driver.

Workers are fully autonomous — no shared state, no shared queues,
no inter-worker synchronisation during their compute. The EM lives
inside the worker because **each variant group's EM converges
independently of every other group's**: the EM's M-step aggregates
across samples within a single group, not across groups. There is
no global prior estimation phase that would force a chrom-level
EM post-pass.

### Driver and streaming emit

The chunk-loop driver runs sequentially per chromosome:

```rust
for chrom in contig_table_order(...) {
    let mut carryover = empty;
    while chunk_start < chrom.end {
        let chunk = load_chunk(..., carryover);   // load + filter + compact
        fix_boundaries(&mut chunk);               // sets safe_end + windows
        carryover = chunk.split_carryover();

        // T workers in parallel, output collected by window index:
        let per_window_outputs: Vec<Vec<FinalVariantRecord>> =
            chunk.windows.par_iter()
                .map(|w| process_window(&chunk, w))
                .collect();

        // Stream each window's records to vcf_writer in window-index order:
        for window_out in per_window_outputs {
            for record in window_out {
                vcf_writer.write(record)?;
            }
        }
        chunk_start = chunk.safe_end;
    }
}
```

Within a chunk, `par_iter().collect()` returns outputs in window-index
order, which is genomic order (windows tile `[range.start, safe_end)`
ascending). Across chunks within a chromosome, chunks tile the
chromosome in ascending order, so concatenation preserves position
order. Across chromosomes, the outer loop iterates in contig-table
order. **No per-chromosome buffer**, **no per-cohort buffer** — the
worker output is the only intermediate, and it flows to the writer
as soon as the chunk finishes.

`vcf_writer`'s contract stays as today's: *"consume an iterator of
records, emit a cohort VCF."* Sub-window / chunk / chromosome
structure is entirely a var_calling internal detail; the writer
would not change if the parallel decomposition were redesigned
again later.

### Chromosome boundaries

Chromosomes remain independent. Chunks never cross chromosome
boundaries — the chunk-loop resets carryover to empty at each
new chromosome (the grouper's `max_group_span` is intra-chromosome
only). The driver above handles this naturally by iterating
chromosomes one at a time.

### File descriptor budget

At large N, per-chunk loading opens N `region_records` iterators
(one per sample per chunk). Across the 2-chunk pipeline that's 2N
open file descriptors at steady state. `RLIMIT_NOFILE` default on
Linux is 1024; at N=512 we hit the cap.

The driver calls `setrlimit(RLIMIT_NOFILE, …)` at startup,
requesting `max(soft, 4N)` (capped at the hard limit). **If the
raise fails — meaning the soft limit can't be raised high enough
to cover 4N fds — the driver fails fast at startup with a clear
error message naming the required limit and pointing the user at
`ulimit -n` / `/etc/security/limits.conf` / their cluster's
equivalent.** We do not implement an fd-pool fallback. Rationale:
the OS already does this job well, the user always knows the
target N at submission time, and a software workaround would
materially complicate the loader (per-sample re-open on each
chunk transition, locking, error handling around partial
re-opens) for a problem that's trivially solved by configuring
the environment correctly. The error message tells the operator
exactly what to set; that's the right division of labour.

## Phasing

Phased so we can bail out at the cheapest possible point if anything
goes sideways. Each phase ends with a measurable, mergeable result.

### Phase A — single-threaded chunk loader, columnar from day one

Goal: byte-identical VCFs at T=1.

- New `MaterialisedChunk` with the **columnar `SampleColumns`
  layout** described above — committed from Phase A, not deferred.
- New `load_chunk` (single-threaded sample loop) doing the
  load → filter → compact sequence: raw column copy from PSP
  blocks, cohort-wide variant-position filter against the
  reference (via `ChromRefFetcher`), compact to variant-only
  `SampleColumns`. Accepts and prepends the previous chunk's
  carryover (already variant-filtered from the prior chunk's run).
- New `fix_boundaries` pre-pass on the **post-filter** chunk:
  computes per-record `allele_end_max`, scans the N-way-merged
  variant-position timeline, picks the chunk's `safe_end`, splits
  carryover, and (at T=1) produces a single window
  `[range.start, safe_end)`. Safe gaps are easier to find than on
  the pre-filter timeline because non-variant positions are gone.
- `PileupRecordRef<'a>` borrowed-view type so existing
  per_position_merger + variant_grouper + per_group_merger
  consumers see the same call shapes (`&PileupRecord` →
  `PileupRecordRef<'a>`). The merger algorithms themselves are
  unchanged in this phase.
- New `GroupedVariantsBatch` columnar accumulator (parallel
  arrays for group geometry, alleles, log-likelihoods, etc.).
  per_group_merger's `process_group` appends into this batch
  instead of emitting row records.
- New chunk-loop driver that processes chunks sequentially, one
  per call. At T=1 the single window equals
  `[range.start, safe_end)`; the worker runs the *full* pipeline
  (merger + grouper + per_group_merger + posterior_engine EM) on
  it and emits final variant records into a per-window output
  buffer.
- Driver streams per-window output buffers to `vcf_writer` in
  window-index order as each chunk completes (no per-chrom
  buffer, no per-cohort buffer). `vcf_writer`'s record-iterator
  contract is unchanged from main.
- All integration tests pass byte-identical against `main`.
- Wall + peak RSS measured at T=1 on N=50 / N=200 / N=1000 via
  `perf_scaling_synthetic.py`. Wall is expected to be similar to
  or slightly worse than main at T=1 (load + filter + pre-pass
  overhead with no parallelism gain yet, partially offset by the
  downstream stages seeing only variant positions); peak RSS
  should already be substantially lower because (a) we're holding
  one chunk × N samples instead of one whole-chromosome × N
  samples per per-chrom worker, and (b) the variant-position
  filter cuts the per-chunk record count by a large factor on
  real data.

Exit criterion: byte-identical VCFs. If correctness fails here,
everything stops and we re-examine.

### Phase B — parallel within-chunk workers

Goal: byte-identical VCFs at T=4. Memory at par with Phase A.

- Extend `fix_boundaries` to place T-1 internal window boundaries
  in safe gaps (forward scan placing boundaries near
  evenly-spaced positions). Output: `windows: Vec<Range<u32>>` of
  T (or fewer, if the chunk has < T safe gaps).
- T workers run the **full pipeline** (merger + grouper +
  per_group_merger + posterior_engine EM) on their respective
  windows in parallel via rayon. **Fully autonomous** — by
  construction no group can span a window boundary, so each
  worker emits every group it sees. Each worker returns its own
  output buffer of final variant records. No lookback, no
  canonical-position rule, no inter-worker coordination.
- Driver streams per-window output buffers to `vcf_writer` in
  window-index order as each chunk completes; per-chunk streams
  concatenate naturally in genomic order across chunks within a
  chromosome.
- Integration tests pass byte-identical (deterministic scheduling
  on the test fixture is required — `rayon::ThreadPool::install`
  with a fixed per-chunk work order). Window assignment is
  deterministic for a given chunk because `fix_boundaries` is a
  pure function of the chunk's records.
- Wall + peak RSS measured vs Phase A and vs main at T=4 on N=50
  / N=200 / N=1000. This is the headline A/B point.

Exit criterion: byte-identical VCFs and a clear-enough perf signal
to decide whether the design is delivering. No hard rule on a
specific number — decision is data-driven once we see the
comparison.

### Phase C — pipelined chunk loading

Goal: wall reduction vs Phase B by overlapping I/O and compute.

- Background-thread chunk loader: while workers are processing
  chunk K, the loader is preparing chunk K+1.
- Bounded queue (depth 2) of `MaterialisedChunk` — natural
  back-pressure if compute is faster than I/O or vice versa.
- No correctness change. Same byte-identical guarantee.
- Measure wall vs Phase B. Expectation: ~10–20 % wall reduction
  from hiding chunk-load latency behind compute. If the chunk
  loader is already much faster than the per-chunk compute,
  this phase delivers little — fine, ship it for hygiene.

### Phase D — performance review pass

Goal: deep performance analysis of the new architecture using the
project's [`rust-performance-review` skill](../skills/rust-performance-review/),
including SIMD opportunities flagged by the review.

The columnar `SampleColumns` layout was committed in Phase A, so
SIMD candidates surfaced here can be applied directly without a
layout refactor:

- Per-sample-batch SIMD in `compute_log_likelihoods` (gather four
  samples' per-(record, allele) scalars into `f64x4` lanes;
  evaluate `standard_log_likelihood` lanewise; scatter back).
- Cache-locality measurement on the per_position_merger pass over
  the columnar storage — confirm the layout pays back in L2 hit
  rate / instructions-retired vs. an artificial row-wise control.
- Whatever else the per-category sub-agents (allocations,
  data_layout, concurrency, hot_loops, io_and_syscalls,
  methodology) surface on the new shape.

Apply the same review methodology used for psp_to_vcf 2026-05-20:
baseline profile via `perf record` on real tomato data + the
`cohort_e2e_perf` bench, dispatch sub-agents per category, produce
a report at `doc/devel/reports/reviews/perf_within_chrom_parallel_<date>.md`.

Phase D is contingent on Phases A–C landing successfully and on the
A/B perf comparison; it's the "what comes next" pass, not a gate
on merge.

## Correctness validation — byte-identical VCFs

Non-negotiable. The rewrite reorders work but must produce identical
VCF records (same chrom/pos/ref/alt/qual/info/format/genotype) as
`main` for every fixture in `tests/`. Validation strategy:

- **Test fixture parity.** All existing integration tests in
  [tests/cohort_cli_integration.rs](../../tests/cohort_cli_integration.rs)
  pass under the rewrite (same input → identical output VCF).
- **Determinism harness.** Window assignment is a pure function of
  the chunk's records (`fix_boundaries` is deterministic), and
  per-window output is collected into a position-indexed buffer
  rather than a shared queue, so emit order is fixed by window
  index regardless of worker completion order. The integration-test
  driver additionally pins `rayon::ThreadPool::install` per chunk
  for belt-and-braces (mirrors the H1 determinism pattern from
  commit `0b1e958`, applied at the window level here).
- **Real-data parity.** Before merge: run `cohort_e2e_perf` and
  `examples/profile_cohort_e2e` against tomato N=10 / N=200 /
  N=1000 on both main and the branch, diff the output VCFs
  field-by-field. Diff must be empty.
- **`bcftools view` smoke.** Manual smoke against a real cohort
  VCF for sanity. Documented but not gating.

## Sequencing

1. **Plan revision on `main`.** This document is on `main` for
   review and is the implementation spec the feature branch will
   reference. Updates to scope / design land on `main` first, with
   the branch rebasing.
2. **Cut feature branch `cohort-within-chromosome-parallel` from
   current `main`.** No work yet on the branch.
3. **Phase A on branch** → impl report committed →
   integration tests pass byte-identical → measure (T=1).
4. **Phase B on branch** → integration tests pass → measure (T=4).
   Headline A/B point: read off wall + peak RSS at N=50/200/1000
   on `main` vs branch.
5. **Decision point.** Look at the data; no commitment to a
   specific decision rule beforehand. Either Phase C/D follow,
   or we pause and revisit the plan.
6. **Phases C, D** as decided.
7. **Merge to `main`** when satisfied. PROJECT_STATUS update with
   final numbers and any deferred follow-ups.

## Open work / non-goals

- **DUST → BED cache** — deferred. Trivial, cheap; can land on
  `main` independently any time.
- **SQUAREM EM acceleration** — deferred. May fall out of Phase D
  if the new architecture's perf review surfaces it as worthwhile.
- **rayon-over-records in the EM** — *subsumed*. The EM now
  runs inside each per-window worker, so per-group EM work is
  already parallel across windows / chunks. Explicit
  rayon-over-records inside a single worker's EM is no longer
  a relevant lever.
- **PSP reader H1 (CSR ragged-column collapse)** — independent of
  this plan; lands on whichever side ships first.
- **`var-calling-from-bam` parallelism** — this rewrite is for the
  `var-calling` (PSP → VCF) path. The from-bam path got its own
  per-chrom rayon in 2026-05-24; revisiting it for chunk-level
  parallelism is a follow-up if/when it matters on real workloads.
- **Sub-cohort batching as a memory fallback** — was on the
  scaling-report recommendation list as a safety net if the
  architecture-rewrite memory targets don't land. If Phase B's
  data shows memory hasn't dropped enough for N=5000 to fit, the
  batching fallback comes back on the table.

## Risks

- **Chunk-size tuning is empirical.** Too small → chunk-transition
  overhead (open N fds, decompress N blocks, materialise N record
  arrays, run filter + pre-pass) dominates; too large → memory
  hit during the raw-load step (before the filter compacts) and
  load imbalance across workers. Plan: ship with a default that
  fits ~1 PSP block per sample at the median density, expose a
  hidden CLI knob for re-tuning, measure on tomato. The
  post-filter footprint is the *steady-state* memory; the
  raw-load step transiently holds ~2× that during the compact
  step.
- **Variant-filter predicate against the PSP record schema.** The
  safe predicate (per-allele byte comparison against the
  reference base) is straightforward but allocates one
  `ChromRefFetcher` lookup per chunk position. If PSP records
  carry an explicit "this is the reference allele" marker, the
  predicate degenerates to "any allele with index ≠
  ref_allele_index has num_obs > 0" — much cheaper. Phase A
  implementation should check the actual `PileupRecord` schema
  and pick the fast path if it exists. If the predicate is
  expensive enough to bottleneck the chunk loader, the
  per-record check is the lever; the higher-level cohort-wide
  filter design doesn't change.
- **Pre-pass can fail to find a safe gap.** Pathological input
  with records every < `max_group_span` bp end-to-end for the
  whole chunk leaves no place to put a boundary. Mitigation:
  dynamically extend the chunk's load range and retry, with a
  hard cap at 4× nominal chunk size; past the cap, fail loudly
  rather than producing a multi-GB chunk silently. Real data has
  inter-record gaps every few hundred bp at most, so this is a
  defensive backstop rather than a primary concern.
- **Determinism for tests.** Window assignment is a pure function
  of the chunk's records, so `fix_boundaries` is deterministic.
  Worker scheduling order can still affect output if workers
  push to a shared queue rather than position-indexed slots; the
  integration test harness pins a `rayon::ThreadPool::install`
  per chunk and uses position-indexed per-window output buffers
  so the concat order is fixed by window index, not by completion
  order.
- **File-descriptor budget.** Default `RLIMIT_NOFILE` of 1024 is
  exceeded at N=512 (worst case 2N = 1024). Plan calls for
  `setrlimit` at driver startup; on failure we fail-fast at
  startup with a clear message naming the required limit. No
  software fallback (the OS does this job well; complicating the
  loader to work around an unraised ulimit isn't worth the
  surface area). Risk reduces to "user gets a clear error and
  has to set their environment up", which is the standard
  contract for any tool with a large fd footprint.
- **`max_group_span` sizes the safe-gap threshold.** Currently
  `DEFAULT_MAX_VARIANT_GROUP_SPAN` (~1 kb). Larger
  `max_group_span` means fewer safe gaps in the chunk, which
  means fewer windows per chunk (less parallelism per chunk).
  Tested today's value works; if a flag bumps it significantly,
  chunk size may need a corresponding bump.
- **Per-group EM independence is load-bearing for the
  worker-runs-EM design.** Each variant group's EM converges
  using only that group's per-sample log-likelihoods; there is no
  cross-group prior estimation. If a future change introduces a
  cross-group dependency (e.g., genome-wide allele-frequency
  estimation feeding back into per-record priors), the EM would
  have to move back out of the worker into a global pass over
  the per-chrom output buffer, with the corresponding memory
  cost. Pre-implementation spike: confirm against the current
  `posterior_engine` that the M-step is purely intra-group, and
  bake a regression-test into the EM module so this invariant
  doesn't drift silently.

## Forward references

These were considered during plan drafting and are not part of this
plan, recorded for context:

- **Sub-cohort batching fallback** (`--max-cohort N` flag chunking
  the cohort into independent groups). If memory targets don't
  land, this is a recommended safety net per
  [scaling_measurement_2026-05-27.md §5](../reports/reviews/scaling_measurement_2026-05-27.md).
- **Cross-group EM if the algorithm ever changes.** Today's EM
  is per-group, which is why this plan can fit the EM inside
  each worker and stream emits to the writer without any
  per-chrom buffer. If a future design adds a genome-wide
  allele-frequency prior or a similar cross-group dependency,
  the per-chrom buffer (records-per-chrom × `sizeof(MergedRecord)`,
  sub-GB on tomato at N=1000) reappears and the EM has to come
  back out of the worker.
- **mimalloc** — the scaling report's allocator-pressure findings
  (38 % inclusive at N=1000) suggest mimalloc may help, but it
  affects main and branch equally. Independent workstream.

## Estimated effort

- Phase A: **5–7 days** (columnar `SampleColumns` + CSR per-allele
  layout + `PileupRecordRef<'a>` view type + `load_chunk` + chunk
  loop driver + tests; algorithm-side per_position_merger /
  variant_grouper / per_group_merger reused with the borrowed-view
  signature change). Larger than a row-wise Phase A would have been
  because columnar is committed up front, but smaller than the
  combined row-wise + later-refactor path the previous draft proposed.
- Phase B: **3–5 days** (extend `fix_boundaries` to multi-window
  partitioning, deterministic per-window output buffers, byte-identity
  validation under T=4).
- Phase C: **2–3 days** (background loader thread, bounded queue).
- Phase D: **3–5 days** (perf review skill run + applying any
  high-value SIMD findings, with the columnar layout already in
  place from Phase A).
- A/B comparison + integration testing + report writeup: **2–3 days**.

Total: **2–4 weeks** of focused work. The columnar-from-Phase-A
decision saves roughly 3–6 days vs. the row-then-refactor sequence
the previous draft had.
