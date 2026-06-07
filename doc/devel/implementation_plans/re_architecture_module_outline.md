# Appendix — module & type outline (new architecture)

Companion to
[re_architecture_streaming_pipeline.md](re_architecture_streaming_pipeline.md)
(read §2 there first). This is the **bird's-eye view**: the modules, the
structs they hold, and the main functions — no algorithmic pseudocode. Its
jobs are to (1) **settle naming**, (2) expose **parts not yet detailed
enough** (collected in the last section), and (3) act as the map during
implementation. The ordered, gated phasing is the separate
[execution plan](re_architecture_execution_plan.md).

**Legend.** `[REUSE]` exists today, kept as-is · `[RENAME]` exists, renamed ·
`[REVISE]` exists, reshaped · `[REVIVE]` existed pre-columnar, restore from
git · `[NEW]` does not exist · `[OPEN]` undecided.

**Build note.** This is built in a parallel **`var_calling_new::`** package
and renamed to `var_calling::` at the swap (plan §7). Module paths below show
the **target** `var_calling::` names; read them as `var_calling_new::` during
the build. The byte-identity-sensitive numeric kernels are **copied verbatim**
from the old package (not rewritten); the structure is built from scratch.

---

## 0. Naming decisions (settle these first)

| concept | proposed name | note |
|---|---|---|
| psp **format** block reader | `psp::PspReader<R>` | `[REUSE]` — **already exists**; do **not** reuse this name for the per-sample reader |
| per-sample chunk reader | **`SamplePspReader`** | `[RENAME]` of today's `ColumnSpanReader` — avoids the `PspReader` collision above. **[OPEN]** vs renaming the format reader instead |
| one sample's chunk | **`SamplePspChunk`** | `[NEW]` |
| one variable cohort position (all samples) | **`CohortPileupRecord`** | `[REVIVE]` (cf. `PerPositionPileups`) — the producer's output unit |
| a work-chunk of those records | **`PileupCohortChunk`** | `[NEW]` — list of `CohortPileupRecord`s + `RefSpan`; **not** the old columnar `MaterialisedChunk` |
| grouped overlapping records | **`OverlappingPileupRecords`** | `[REVIVE]` record form (see `OverlappingVariantGroup`) |
| final output record | **`Variant`** | `[RENAME]` of `PosteriorRecord` |
| per-chunk genomic-order counter | **`chunk_order`** | `[RENAME]` of `seq_idx` (field) |
| section 1 (produce chunks) | **`CohortChunkIntegrator`** | `[NEW]` (wraps today's `BlockIterator`) |
| section 2 (group + call) | **`VariantCaller`** | `[NEW]` — groups *and* genotypes |
| section 3 (write) | **`VcfWriter`** | `[REUSE]`-ish |

---

## A. Per-sample reading — module `var_calling::sample_reader`

*(today: `column_span_reader.rs`, backed by `psp::reader`)*

**`SamplePspReader<R: Read + Seek>`** `[RENAME ColumnSpanReader]`
One per sample. Created from a `psp::PspReader<R>` + the region/bed segments.
No dust (§2.3). `!Send`.
- `fn new(reader, regions) -> Self`
- `fn next_chunk(&mut self, span_end) -> Result<SamplePspChunk>` — the unit it hands to the producer
- `fn peek_next_span(&self) -> Option<u32>` `[REUSE]` — the **next psp segment
  (block) boundary** (`last_pos`): how far this reader serves without a new
  decode
- `fn reset(&mut self, …)` `[REUSE]`

**Segment-aligned reads (keep — it's a perf invariant).** The natural read
unit is the psp **segment (block)**. The producer takes `min(peek_next_span)`
across readers as the watermark, so each `read_span` costs **at most one new
block decode per sample** ([column_span_reader.rs:102-123](../../src/var_calling/column_span_reader.rs#L102)).
Because the per-sample `.psp` files are segment-aligned, the watermark lands
on a shared block boundary and every reader pulls **one whole aligned
segment** in lockstep — all columns for those records co-located, sequential
I/O. This **composes with column-selective decode**: from each segment, decode
the **light** columns to decide which positions vary, then the **heavy**
columns for just those positions (building the records), then **drop the
segment** — decoded once, no retained compressed bytes. A `PileupCohortChunk`
(work-chunk) spans the records of several segments (accumulated to
`target_variants`), but every disk read is one aligned segment per sample.

**`SamplePspChunk`** `[NEW]`
One sample's compressed columns for **one psp segment** (the natural read
unit — see "Segment-aligned reads" above), plus its column manifest. **Which
columns load when is encoded in the methods, not a tag type:** `new()` decodes
the light fold columns eagerly; each getter decodes its heavy column(s)
lazily, for the variable rows only. Each getter **moves the result out** and
frees the underlying compressed bytes (fetched at most once).
- fields: the segment's compressed column store + manifest, the cached light
  fold slices (positions / per-position non-ref obs / per-position ref-span)
- **light accessors** (decoded + cached at `new`, read by the cohort fold;
  no owned struct handed out): `fn positions(&self) -> &[u32]`,
  `fn nonref_obs(&self) -> &[u32]`, `fn ref_spans(&self) -> &[u32]`
- **typed getters** (no generic `Column` view — both ends statically typed):
  `fn take_seq(&mut self, keep: &[bool]) -> PerAlleleSeq`,
  `fn take_chain_ids(…) -> PerAlleleChainIds`,
  `fn take_fixed(…) -> PerAlleleFixed` (whole scalar SoA —
  `num_obs`/`q_sum`/`fwd`/… — together; merger + EM need all of it). The
  producer's per-position merge assembles these per-sample pieces into
  `CohortPileupRecord`s.
- **[DEFERRED]** low-memory mode (re-read from disk instead of holding the
  segment compressed) — only if measurement vs `main` calls for it (gap #9)

> **[RESOLVED #11]** With the record-streaming producer (§2.2 / §B) there is
> **no mid-segment-cut problem**: a `SamplePspChunk` is **one segment**,
> decoded **once** into per-position records; the work-chunk boundary is a
> safe gap *between whole records*, so no compressed segment is ever split or
> re-decoded.

**No `ColumnTag` type.** The §4.1 "decode-some, skip-the-rest" lever rides on
the psp format's **existing** per-block manifest: each `ColumnManifestEntry`
already carries a stable `tag: u16` + compressed length
([psp::block](../../src/psp/block.rs#L829)), so a method decodes the tags it
wants and seeks past the others. The var_calling side needs no parallel enum.

**No `LightColumns` type.** The light fold inputs — per position:
`positions` (from `DeltaPos`), non-ref obs (from `AlleleObsCount`, summed over
non-ref alleles), and the REF allele's `ref_span`/reach (from
**`AlleleSeqLen`**) — are decoded + cached inside `SamplePspChunk` at `new`
and read by the fold through the light accessors above. No owned struct is
handed across; the fold reads each chunk in place (cf. today's
`fold_sample(&SampleColumns)` reading via `position_at` / `non_ref_obs_sum_at`
/ `ref_span_at`). This works *because* the format stores allele **lengths**
(`AlleleSeqLen`) in a column separate from the sequence **bytes**
(`AlleleSeq`), so `ref_span` needs no heavy decode. (`AlleleObsCount`/`num_obs`
is also part of `PerAlleleFixed` (§C), so it could be retained from the fold
rather than re-decoded at materialise — minor.)

The take-getters reuse the existing destination sub-structs (§C):
`PerAlleleSeq` (CSR `offsets`+`bytes`), `PerAlleleChainIds` (CSR),
`PerAlleleFixed` (scalar SoA).

---

## B. Cohort producer (section 1) — module `var_calling::cohort_integration`

*(today: `driver.rs` `BlockIterator`, `two_pass.rs`, `loader.rs`)*

**`CohortChunkIntegrator`** `[NEW]` (wraps today's `BlockIterator`)
Owns the N `SamplePspReader`s, the `DustAheadPool`, and the REF fetcher.
**Streams `CohortPileupRecord`s, sliced into chunks at safe gaps.** Runs §2.2
steps 1–5.
- `fn produce_chunk(&mut self, out: &mut PileupCohortChunk) -> Result<u32>`
  (returns variable-record count; `0` = exhausted) — cf. `produce_block`
- `fn run(self, tx: Sender<PileupCohortChunk>)` — the section's thread body;
  stamps each chunk with a monotonic `chunk_order` (genomic order), the ticket
  the writer reorders on (§E)

**`CohortPerPositionMerge`** `[REVIVE per_position_merger + REVISE two_pass]`
The cohort join: walks all samples' light columns by position to find the
variable positions (AC / `min_alt_obs`), then builds one `CohortPileupRecord`
per variable position from the readers' `take_*` getters. This is the old
`per_position_merger` (`PerPositionPileups`) brought back as the producer's
core — it's the fan-in barrier.
- `fn next_variable_record(&mut self) -> Option<CohortPileupRecord>` (or a
  batched form) — emits per-position cohort records in genomic order
- AC keep rule `[REUSE derive_is_kept]`; obs/ref-span fold across samples
  `[REUSE fold_sample logic]`

**Free functions** `[REUSE]`
- `find_block_cut(...) -> u32` — safe-gap chunk boundary (from record reach)
- apply dust: drop dust-masked positions (after the obs/keep decision)
- `fetch_chunk_ref_span(...)` `[REVISE prefetch_window_ref_bytes]` — one
  contiguous REF buffer per chunk (replaces per-group `Vec<Vec<u8>>`)

**Control flow** (record-streaming). To produce one chunk:
1. **lockstep read** all readers' next segment up to the **watermark** =
   `min` over readers of `peek_next_span()` (segment boundary), decoding only
   the **light** columns;
2. **merge across samples by position** → variable positions (AC /
   `min_alt_obs`), then **apply dust** (drop masked positions — the caller's
   grouper bridges the gaps they leave);
3. **build a `CohortPileupRecord`** per variable position (heavy columns via
   `take_*`, variable rows only); each segment is decoded **once**;
4. **accumulate records; cut at a safe gap** (`find_block_cut`) — the chunk
   boundary falls *between whole records*, so nothing is split;
5. **fetch the chunk's REF span**; ship the chunk (records + `RefSpan`).

**Memory — the key property (must be preserved).** The cohort-wide in-memory
footprint is **only the consolidated records: variable positions, AC /
`min_alt_obs` already applied**. Non-variant positions (~96 %) are dropped at
step 2 and **never built or held cohort-wide**. So RAM ∝ *variable* positions
× samples × in-flight chunks — not all positions. This preserves today's
AC-pushdown saving; never hold a full-coverage cohort structure.

**Scope of change vs today** (≈ moderate, net-simpler):
- `[REUSE]` the watermark/segment loop, the AC keep rule, the safe-gap cut
  (today's `fill_block` / `rebuild_position_union_and_is_kept` /
  `find_block_cut`); `[REVIVE]` the per-position merge (`per_position_merger`)
  as the record builder;
- **leaves** the integrator: `partition_window` (grouping → `VariantCaller`)
  and per-group REF prefetch (→ one `fetch_chunk_ref_span`) — *simpler*;
- the real new work: column-selective decode (light to decide variable, heavy
  only for the ~4 % that vary), and emitting **records** instead of a columnar
  block. The allocation churn of many small records is the thing to measure.

---

## C. Shared data types — module `var_calling::types` (or co-located)

**`CohortPileupRecord`** `[REVIVE — cf. `PerPositionPileups`]`
One **variable cohort position** with every sample's data at it. The unit the
producer emits and the grouper consumes.
- per position: `chrom_id`, `pos`, and per sample (those with a record there)
  its alleles + obs/stats. Carries enough for grouping (reach/`ref_span`) and
  for the per-group merge + EM downstream.

**`PileupCohortChunk`** `[NEW — replaces columnar `MaterialisedChunk`]`
The producer→caller **work-unit**: a list of `CohortPileupRecord`s for one
safe-gap-bounded span, **not a columnar block**.
- fields: `chunk_order: u64`, `records: Vec<CohortPileupRecord>`,
  `ref_span: RefSpan`
- "chunk" = one ordered piece of parallel work; its contents are records, so
  the caller groups them without first splitting a columnar block.

**`SampleColumns`** `[REUSE — within `SamplePspChunk`]` — one sample's columnar
slice (CSR ragged columns + parallel scalar columns); the decode target inside
a `SamplePspChunk`, consumed by the per-position merge.

**`RefSpan`** `[NEW]` — the chunk's contiguous REF bytes + genomic start,
**owned** (it crosses the producer→caller queue with the chunk; can ride the
block-buffer recycling so it isn't a fresh alloc per chunk).
- `{ genomic_start: u32, bytes: Vec<u8> }` — covers `[genomic_start,
  genomic_start + bytes.len())`
- `fn slice(&self, group_start: u32, group_end: u32) -> &[u8]` — offset
  arithmetic; per-group REF span is exactly `[start, end]`, **no padding**
- **Fetched only on the producer** (the single thread): one
  **monotonic-forward** fetch per chunk of `[first_kept_pos, max_reach]`,
  `max_reach = max(pos + ref_span − 1)` over kept positions (`ref_span` is a
  light column, already in hand), clamped to `chrom_length`. *Required*, not
  just tidy: `StreamingChromRefFetcher` panics on non-monotonic access
  ([fetcher.rs:299](../../src/fasta/fetcher.rs#L299)), and the parallel
  callers can't fetch monotonically. Replaces the per-group `Vec<Vec<u8>>`.
- **No halo / no cross-chunk REF**: the safe-gap cut makes groups
  self-contained (`max_reach < safe_end`), so the chunk's own span covers
  every group's slice.

**`OverlappingPileupRecords`** `[REVIVE]` — record-based overlapping
group(s). **`VariantCaller`-local**, never queued.
- **[RESOLVED #3]** the row shape is **in-tree, not a git relic**: it's
  `variant_grouping::OverlappingVariantGroup` (+ `per_group_merger::MergedRecord`),
  currently kept `#[cfg(test)]` as the oracle the columnar path is validated
  against. We promote it to production.

**`Variant`** `[RENAME PosteriorRecord]` — final emitted record (GT/GQ/AD/AF/
AC/FILTER/QUAL).

---

## D. VariantCaller (section 2) — modules `var_calling::pileup_overlaps` + `var_calling::em_posterior_calc`

Section 2 is **two modules that work together** — grouping then calling —
plus the thin worker that composes them. Split because they're distinct
concerns (build overlapping record groups vs. run the math on them), even
though one worker drives both back-to-back.

**`VariantCaller`** `[NEW]` (wraps today's `process_block`/`run_window`) —
the composer. W parallel workers; `PileupCohortChunk` → `Vec<Variant>`; owns
the per-worker scratch. Two steps: **(1) generate** the
`OverlappingPileupRecords` (`pileup_overlaps`), then **(2) process** each
record → `Variant` (`em_posterior_calc`); slices REF per group from the
chunk's `RefSpan`. (Home: co-located with `em_posterior_calc`, the stage that
emits the output.)
- `fn call_chunk(&mut self, chunk: &PileupCohortChunk) -> Result<Vec<Variant>>`
- `fn run(self, rx, tx)` — pop chunks, push one `CalledChunk` per chunk
  (tagged with `chunk_order`), **even when empty** (see §E / gap #7)

**Record-based grouping; SoA + SIMD EM kept.** Grouping (D.1) goes
record-based — its row builder `build_overlapping_variant_group` is in-tree
as the `#[cfg(test)]` oracle, so we *promote* it (byte-identical by
construction, and grouping isn't perf-critical). The EM (D.2) keeps its
per-record **SoA** layout + SIMD `ln`/`exp` backend — it is **not** reverted
to the AoS scalar oracle, so no SIMD loss. "Record-based" is about the
grouping and the iteration unit (already per-record), not the EM's inner
numeric layout.

### D.1 `var_calling::pileup_overlaps` — grouping
*(today: the `#[cfg(test)]` `build_overlapping_variant_group` in `worker.rs`,
+ `variant_grouping.rs`)*
the chunk's `CohortPileupRecord`s → `OverlappingPileupRecords` (§C):
walk the records, joining ones whose reach overlaps into groups (bridging the
gaps left by dropped dust positions). `[REVIVE]` the row builder to
production.
- `variant_grouping` — overlapping-group construction (`GrouperConfig`,
  `OverlappingVariantGroup`)
- **[RESOLVED #4]** `per_position_merger` is **revived — in the *producer***
  (§B), where it builds the `CohortPileupRecord`s (the per-sample→per-position
  cohort join). It is *not* here in the caller: the caller already receives
  per-position records and only groups them.

### D.2 `var_calling::em_posterior_calc` — per-group merge + per-record EM
*(today: `per_group_merger.rs`, `posterior_engine/`, `kernels/`, `worker.rs`
— `[REUSE]`)*
`OverlappingPileupRecords` → `Vec<Variant>`, **one record (group) at a time**.
- **per-record EM** — each group runs an independent EM to convergence; the
  iteration count is per-record, so batching records column-wise would waste
  iterations on already-converged ones. `run_em_columnar` is *already*
  per-record.
- **keep SoA + SIMD** — the speed is the lane-of-4 `ln`/`exp` backend
  (`InterpUnivariateSimdMath`), and its lanes fill from the group's **SoA**
  likelihood cells (samples × genotypes). So the EM's *inner* buffers stay
  columnar: `compute_log_likelihoods_columnar` + `run_em_columnar` are
  `[REUSE]`, **not** reverted to the AoS scalar oracle (`run_em_for_record`
  stays a test oracle).
- `per_group_merger` — the per-group merged alleles feeding the EM; row-shape
  `MergedRecord` is fine *as long as* it hands **SoA** likelihood buffers to
  the EM (else keep the columnar merge).
- **[RESOLVED #5]** per-record iteration (already so) + SoA-within-record for
  SIMD → **no SIMD loss**. "Record-based" is the grouping (D.1) and the
  iteration unit, not the EM numeric layout.

---

## E. Writer (section 3) — module `var_calling::vcf_writer`

**`VcfWriter`** `[REUSE/REVISE]`
`Variant`s → VCF, **single consumer**, reorders by `chunk_order`.
- `fn run(self, rx: Receiver<CalledChunk>)`
- holds the `BTreeMap<chunk_order, …>` reorder buffer `[REUSE]`, draining on
  the *next expected* `chunk_order`
- **gapless invariant**: every chunk yields exactly one `CalledChunk` (empty
  ones included), or the drain stalls; add a debug-assert that the consumed
  `chunk_order` sequence is contiguous, so a dropped chunk fails loudly in
  tests instead of silently deadlocking
- **[OPEN]** where the reorder buffer lives — in the writer, or a shim stage

**`CalledChunk`** `[RENAME BlockResult]` — `{ chunk_order: u64, Vec<Variant>,
stats }`, the caller→writer payload. `records` may be **empty** (a chunk that
produced no `Variant`s still ships, to keep `chunk_order` gapless); an empty
`Vec` doesn't allocate, and the writer treats it uniformly.

---

## F. DUST helper — module `var_calling::dust`

*(today: `dust_filter.rs`)* — fine as-is (§2.3).
- `DustAheadPool` `[REUSE]` — ordered look-ahead mask pool
- `IntervalDustMask` `[REUSE]` — per-interval mask, owned by the producer

---

## G. Driver / wiring — module `var_calling::pipeline`

*(today: `driver.rs` `drive_blocks_parallel`)*
Spawns the three sections, owns the two bounded `crossbeam-channel`
hand-offs, joins, propagates first-error.
- `run_var_calling(args) -> Result<RunStats>` `[REUSE entry]`
- `fn drive_pipeline(producer, callers, writer) -> Result<…>` `[REVISE]`
- queues: `producer → caller` (MPMC, `PileupCohortChunk` = records + `RefSpan`),
  `caller → writer` (single-consumer, `CalledChunk`)
- queue admission: simple **count cap** (crossbeam's built-in bound);
  byte-budget deferred (gap #8)

---

## H. Out of scope / orthogonal

- `contamination_estimation` `[REUSE]` — separate analysis, not on the
  restructure path
- `psp::*` format layer `[REUSE]` — only `PspReader<R>` is touched (the
  per-sample reader builds on it); column-selective decode (§4.1) extends
  `decode_block_payload`

---

## I. Parts still not detailed enough (the gaps to close)

1. ~~`PspReader` name collision~~ **[RESOLVED]** — per-sample reader is
   `SamplePspReader`; the format reader keeps `psp::PspReader`. (§0)
2. ~~`Column` / `ColumnId` shape~~ **[RESOLVED]** — typed getters
   (`take_seq`/`take_chain_ids`/`take_fixed`) returning the existing
   sub-structs; no cross-boundary `Column` type, and no `ColumnTag` enum —
   eager-at-`new` vs lazy-on-getter is encoded in the methods, and the
   skip-the-rest lever uses the psp manifest's existing per-column tags. (§A)
3. ~~`OverlappingPileupRecords` fields~~ **[RESOLVED]** — it's the in-tree
   `OverlappingVariantGroup` (+ `MergedRecord`), the `#[cfg(test)]` oracle;
   promote, not recover-from-git. (§C/§D)
4. ~~`per_position_merger`'s fate~~ **[RESOLVED]** — **revived in the
   producer** (§B): it builds the `CohortPileupRecord`s (the per-sample→
   per-position cohort join). Not in the caller. (Flipped from the earlier
   "gone" once we moved to record-streaming.)
5. ~~EM internals: columnar vs record~~ **[RESOLVED]** — per-record iteration
   (unchanged) **with SoA-within-record + SIMD `ln`/`exp` kept**
   (`run_em_columnar`). 'Record-based' is the grouping, not the EM numerics —
   **no SIMD loss**. (§D.2)
6. ~~`RefSpan` API~~ **[RESOLVED]** — producer-side single monotonic-forward
   fetch of `[first_kept_pos, max_reach]` (forced by the forward-only
   fetcher) → owned `RefSpan` shipped with the chunk → caller slices by
   offset. No halo: safe-gap cut makes groups self-contained. (§C)
7. ~~`seq_idx` ownership~~ **[RESOLVED]** — renamed **`chunk_order`**;
   producer stamps it (monotonic, genomic order), caller carries it through
   untouched, writer reorders via `BTreeMap`. Gapless contract: one
   `CalledChunk` per chunk (empty included) + a contiguity debug-assert in
   the writer. (§E)
8. ~~Queue admission~~ **[DECIDED]** — simple **count cap** now (crossbeam's
   built-in bound). Byte-budget + `--memory-budget` is a deferred, additive
   bonus, to revisit with a profile (§4.2). (§G)
9. **Low-memory mode** — **[DEFERRED — measurement-gated]** don't build it
   until the new architecture's RAM is measured against `main`. The count-cap
   + hold-compressed-segment design may already be enough; only if the
   N×threads peak is a problem do we add a `SamplePspChunk` re-read mode. (§A)
10. ~~Module homes~~ **[DECIDED]** — **flat** module structure (no subtree
    unless it later reads clearer). Built in a parallel **`var_calling_new::`**
    package and renamed to `var_calling::` at the swap (build strategy, plan
    §7). The module paths in this appendix are the **target** (`var_calling::`)
    names; during the build they live under `var_calling_new::`.
11. ~~`SamplePspChunk` granularity & mid-segment cut~~ **[RESOLVED — dissolved
    by record-streaming]** a `SamplePspChunk` is **one psp segment**, decoded
    **once** into per-position `CohortPileupRecord`s; the work-chunk cut falls
    at a safe gap *between whole records*, so no compressed segment is split or
    re-decoded. (§A / §2.2)
