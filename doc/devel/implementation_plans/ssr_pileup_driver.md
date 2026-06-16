# SSR Stage 1 — read I/O driver + Stage-1 driver (implementation plan)

**Status:** refreshed 2026-06-16 after the fetcher landed. **Increment #1
(the fetcher) is BUILT** (commit `696bf7e`) — this doc now records what shipped
and re-scopes the remaining work to the **driver loop (#2)** and **CLI (#3)**.
Grounded in the actual seams: the pooled `segment_reader::AlignmentFile`,
`MappedRead`, `triage::reaches_locus`/`analyze_read`/`aggregate`,
`PspWriter<SsrKind>`.

> **What changed from the original sketch.** The fetcher was sketched against
> `AlignmentMergedReader::query` (per-call re-open). It instead landed on the
> **pooled indexed-segment reader** (`AlignmentFile::get_reads_from_segment`) —
> the primitive built for exactly the ~10⁶-tiny-loci case — and the admission
> gate was settled as **spanning-capable-only** (not "admit flanking/in-repeat
> too"). Both are reflected below.

Cross-refs: architecture §8 (worker topology); `fetch_reads.rs` module doc;
[`segment_read_fetcher.md`](segment_read_fetcher.md) (the reader primitive);
[`ssr_pileup.md`](ssr_pileup.md) (overall Stage-1 plan);
[`ssr_stage1_remaining.md`](ssr_stage1_remaining.md) §2.

---

## 0. The one big architectural decision — confirm first

**SSR Stage 1 does NOT reuse the SNP pileup walker / BAQ / `with_stage1_chain`.**
It is **locus-oriented**, not position-streaming. The seams make this clean:

- The pooled `AlignmentFile` yields `MappedRead` (qname, pos, mapq, cigar,
  uppercased seq, raw qual, …), applying only the **cheap** `SegmentReadFilter`
  (flags / MAPQ / length). It deliberately does **not** apply the reference-
  dependent F1/F2/F3 filters (mismatch fraction, bad-CIGAR) or indel left-
  alignment — SSR realigns every spanning read with the pair-HMM, so it distrusts
  the mapper's CIGAR those filters lean on (see `SegmentReadFilter` docs).
- `triage_read(&MappedRead, &Locus)` and `analyze_read(&MappedRead, &Locus, …)`
  consume that `MappedRead` directly.
- The catalog's `Locus` **is** `crate::ssr::types::Locus` — the exact type the
  analysis takes. No locus conversion.

So the driver iterates raw `MappedRead`s from `get_reads_from_segment` straight
into triage/analyze — **no BAQ, no `PreparedRead`, no `PileupWalker`**. (BAQ is a
SNP false-positive-near-indels mitigation; SSR's pair-HMM models the indels
itself, and HipSTR/GangSTR use raw base qualities. If BAQ ever proves wanted, it
slots in as a per-read transform later — not a v1 concern.)

This follows the project's "build the worker natively over the data, don't
adapter-shim the old streaming pipeline" rule. **Confirmed** — the shipped
fetcher (#1) consumes raw `MappedRead`s with no walker/BAQ.

---

## 1. File layout

```
src/ssr/pileup/
├── fetch_reads.rs   # DONE: Reservoir + locus_seed + reaches_locus consumer;
│                    #   fetch_locus_reads -> LocusReads (per-locus query + gate + cap)
├── driver.rs        # NEW: the Stage-1 run loop — catalog walk -> fetch -> analyze
│                    #   -> aggregate -> name->chrom_id adapt -> write_locus
└── mod.rs           # re-exports
src/pop_var_caller/
└── ssr_pileup.rs    # NEW: SsrPileupArgs (clap) + run_ssr_pileup -> driver::run
```

Increment boundaries (each its own commit + tests):
1. ~~**Fetcher** (`fetch_reads.rs`)~~ — **DONE** (`696bf7e`): `fetch_locus_reads`
   over `&[AlignmentFile]` → reach gate → reservoir → `LocusReads`.
2. **Driver loop** (`driver.rs`): build `Vec<AlignmentFile>` → catalog walk →
   fetch → analyze → aggregate → adapt → write; header build; name→chrom_id.
   **Structure all per-locus work as `process_locus(locus, …, &mut scratch) ->
   record`** so the loop is a one-liner and the later parallel step is just
   `for` → `par_iter().map_init(…)` (see §4).
3. **CLI** (`ssr_pileup.rs`): `ssr-pileup` subcommand wired into `PopVarCallerCommand`.

**We build #2 next.**

---

## 2. Inputs & shared handles — build the `AlignmentFile`s once

The fetcher takes `&[AlignmentFile]` (the pooled per-file readers), so the driver
constructs them once at startup from the SNP loader's `PileupInputs`:

```rust
// once, at startup:
let inputs = load_pileup_inputs(&alignment_files, &reference, build_index)?; // PileupInputs
//   -> sample_name, contigs: ContigList, headers: [Arc<Header>], indexes: [AlignmentIndex]
let repo   = build_fasta_repository(&reference)?;   // load_pileup_inputs does NOT return it; CRAM needs it

let files: Vec<AlignmentFile> = inputs.indexes.into_iter().zip(inputs.headers).enumerate()
    .map(|(i, (index, header))| AlignmentFile::from_input(
        alignment_files[i].clone(), header, index, Some(repo.clone()),
        SegmentReadFilter::default(),   // mapq/dup/qc-fail/length — the SNP-shared cheap filter
        /* source_file_index = */ i,
    ))
    .collect::<Result<_, _>>()?;
```

`ssr-pileup` accepts `&[PathBuf]` alignment files (1+), mirroring `pileup`'s
multi-file input — same `load_pileup_inputs`, same cross-file `@SQ`/sample checks.
The pooled reader applies the cheap `SegmentReadFilter` (mapq / dup / qc-fail /
length) before yielding, so reads reaching the gate are already quality-filtered,
and the drops are counted (`LocusReads.filtered` → `n_filtered`).

---

## 3. The fetcher (increment 1) — `fetch_reads.rs` — **BUILT**

As shipped in `696bf7e`:

```rust
/// One locus's reads + the fetch-pass tallies the aggregator can't recover.
pub(crate) struct LocusReads {
    pub(crate) reads: Vec<MappedRead>,  // reservoir-sampled admitted reads (<= cap)
    pub(crate) yielded: u64,            // reader-yielded reads at the locus (post cheap filter)
    pub(crate) filtered: FilterCounts,  // the reader's flag/MAPQ/length drops -> n_filtered
}

pub(crate) fn fetch_locus_reads(
    files: &[AlignmentFile],   // pooled per-file readers, built once (§2)
    locus: &Locus,
    cap: usize,                // MAX_READS_PER_LOCUS
) -> Result<LocusReads, AlignmentInputError>;
```

Mechanics (as built):
1. **Window = the locus's embedded `ref_bytes` window** (tract ± the catalog's
   flank margin), converted 0-based half-open → **1-based inclusive**:
   `seg_start = ref_bytes_start + 1`, `seg_end = ref_bytes_start + ref_bytes.len()`.
   No separate `flank_bp` param and no contig-end clamp needed — `ref_bytes` is
   already clamped at contig ends by the catalog builder. *(Off-by-one pinned by
   the segment_reader's 1-based-inclusive tests + the fetcher tests.)*
2. For each file: `file.get_reads_from_segment(locus.chrom(), seg_start, seg_end)?`
   — the **pooled** reader (no per-call file re-open).
3. One per-locus `Reservoir::new(cap, locus_seed(locus.chrom(), locus.start()))`,
   shared across the sample's files (order-independent → per-file concatenation).
4. For each yielded `read`: `yielded += 1`; **admission gate**
   `if reaches_locus(&read, locus) { reservoir.offer(read) }`. The gate admits a
   read only if it is **spanning-capable** — footprint brackets both tract ends,
   *or* it is soft-clipped (always admitted; the clip may carry a long allele).
   Reach-discarded reads get no reservoir slot (settled: the cap budget is spent
   on spanning evidence; we do **not** reservoir flanking/in-repeat reads).
5. After draining each file, `filtered.merge(reads.filter_counts())` — the reader
   owns the cheap filter, so it owns the `n_filtered` source.
6. `LocusReads { reads: reservoir.into_held(), yielded, filtered }`.

**Determinism:** the reservoir is seeded only by `locus_seed(chrom, start)` and
fed in fixed (file, coordinate) order, so the sample is identical across runs and
(later) thread counts.

> **Consequence of spanning-only admission (accepted v1 limitation).** Because
> flanking/in-repeat reads are *not* admitted, `aggregate`'s `n_flanking`/`n_frr`
> (derived from worker outcomes) count only the **soft-clipped** non-spanning
> reads that slip through the gate — they are not a full flanking/in-repeat
> census. Accepted per the "don't reservoir the discarded reads" decision; the
> spanning evidence and the depth cap are what matter for genotyping.

---

## 4. The driver loop (increment 2) — `driver.rs`

```rust
pub(crate) struct SsrPileupConfig {
    pub alignment_files: Vec<PathBuf>,
    pub reference: PathBuf,
    pub catalog: PathBuf,
    pub output: PathBuf,
    pub merged_reader: AlignmentMergedReaderConfig,  // mapq/dup filters
    pub window: u16,                                  // analyze_read rung window (NOT flank_bp)
    // ... block layout, threads (later) ...
}

pub(crate) fn run(cfg: &SsrPileupConfig, date: String) -> Result<(), SsrPileupError>;
```

### The per-locus unit — `process_locus` (the parallelization seam)

All per-locus work is one self-contained function over **shared `&` inputs** plus
a `&mut LocusScratch`, returning the container-ready record. The driver loop
becomes a one-liner, so the future parallel step is "turn the `for` into a `map`"
— nothing inside changes.

```rust
/// Reusable per-locus scratch (the alloc-churn win). One per worker thread.
pub(crate) struct LocusScratch {
    cands: Vec<CandidateAllele>,
    hmm: PairHmmScratch,
}
impl LocusScratch { fn new() -> Self { /* … */ } }

/// Everything for one locus: fetch → analyze → aggregate → name→id adapt.
/// Pure over its shared inputs (`files`/`model`/`name_to_id` are `&`, and
/// `AlignmentFile` is `Sync`), so it is safe to call concurrently for different
/// loci. Returns the chrom_id-keyed record ready for `write_locus`.
pub(crate) fn process_locus(
    files: &[AlignmentFile],
    locus: &Locus,
    params: &LocusParams,                 // { window, cap } — the calibration scalars
    model: &HmmModel,
    name_to_id: &HashMap<&str, u32>,
    scratch: &mut LocusScratch,
) -> Result<psp::registry_ssr::SsrLocusRecord, SsrPileupError> {
    let fetched = fetch_locus_reads(files, locus, params.cap)?;
    let outcomes: Vec<_> = fetched.reads.iter()
        .map(|r| analyze_read(r, locus, params.window, model,
                              &mut scratch.cands, &mut scratch.hmm))
        .collect();
    let record_named = aggregate(locus, &outcomes, qc_counts(&fetched)); // chrom NAME
    to_container_record(record_named, name_to_id)                         // chrom_id
}
```

`run()` shape — the loop is now trivial:

```text
inputs = load_pileup_inputs(files, reference); repo = build_fasta_repository(reference)
files  = build Vec<AlignmentFile> from inputs + repo (see §2)
catalog = CatalogReader::new(File::open(cfg.catalog))          // header carries flank_bp, reference_md5
header  = build_ssr_writer_header(sample, &inputs.contigs, ...) // ChromosomeEntry md5/length, params
writer  = PspWriter::<_, SsrKind>::new_ssr(File::create(output), header)
name_to_id = { contig.name -> index } from inputs.contigs
model = HmmModel::default(); params = LocusParams { window, cap }
let mut scratch = LocusScratch::new()
for locus in catalog.read_locus():                             // sorted (chrom, start)
    let locus = locus?
    writer.write_locus(&process_locus(&files, &locus, &params, &model, &name_to_id, &mut scratch)?)?
finish: writer.finish()? + sync; (atomic temp->final rename like the SNP path)
```

The fetcher owns the per-file FASTA repository (handed to each `AlignmentFile` in
§2). The SNP path's per-contig `repo.clear()` was a streaming-memory tactic; for
v1's per-locus access it is **not** needed — revisit only if profiling shows the
repository cache growing unbounded.

### Parallelizing later (why `process_locus` makes it cheap)

`process_locus` takes only shared `&` inputs (and `AlignmentFile` is `Sync` — the
reader pool was built for concurrent per-segment calls) plus thread-owned
`&mut LocusScratch`. So the next step is essentially the `for` → `map`:

```rust
let records: Vec<_> = loci.par_iter()
    .map_init(LocusScratch::new, |s, locus| process_locus(&files, locus, &params, &model, &name_to_id, s))
    .collect::<Result<_, _>>()?;
for r in records { writer.write_locus(&r)?; }   // collect preserves catalog order
```

`map_init` gives each worker its own reused scratch (keeps the alloc-churn win).
Two caveats this simple form carries, to address in the §8.4 step — not now:
- **Ordering:** `write_locus` wants sorted `(chrom,start)`; `collect()` preserves
  input order, so order holds — but it **materializes all records** (~10⁶). The
  arch §8.4 bounded-queue + ordered-collector is the memory-bounded replacement.
- **Determinism:** unaffected — each locus's reservoir is seeded by
  `locus_seed(chrom,start)`, independent of thread/scheduling.

### `LocusReads` → `QcCounts` (§7.3 — **decided**)

All three columns describe **independent primary alignments** at the locus:
duplicates are auditable as filtered but excluded from coverage; secondary /
supplementary (non-independent) and unmapped (not at the locus) are excluded
everywhere.

```rust
fn qc_counts(f: &LocusReads) -> QcCounts {
    let d = &f.filtered;
    // Primary reads dropped for quality (not duplicates, not non-primary).
    let primary_quality_drops = d.qc_fail + d.low_mapq + d.too_short;
    QcCounts {
        // Usable primary reads considered (passed the cheap filter + overlapped;
        // uncapped — pre reach-gate / reservoir).
        depth:        f.yielded as u32,
        // Primary reads the gate dropped: quality + duplicates.
        n_filtered:  (primary_quality_drops + d.duplicate) as u32,
        // Dup-free primary coverage denominator (excludes dup / secondary /
        // supplementary / unmapped).
        mapped_reads: (f.yielded + primary_quality_drops) as u32,
    }
}
```

Relationships: `mapped_reads = depth + (qc_fail + low_mapq + too_short)` so
`mapped_reads ≥ depth` (coverage ≥ usable depth); `n_filtered` adds duplicates on
top, so it is *not* simply `mapped_reads − depth` — duplicates show as filtered
without counting as coverage. (No `FilterCounts::total()` needed; the buckets are
named explicitly. `secondary`/`supplementary`/`unmapped` are intentionally
unused; `high_mismatch`/`bad_cigar`/`baq` are always 0 from `segment_reader`.)

### Name → chrom_id adapter (the §1 boundary)

```rust
/// In-memory ssr::pileup::SsrLocusRecord (chrom NAME) -> container
/// psp::registry_ssr::SsrLocusRecord (chrom_id). The only field that changes.
fn to_container_record(
    r: ssr::pileup::locus_record::SsrLocusRecord,
    name_to_id: &HashMap<&str, u32>,
) -> Result<psp::registry_ssr::SsrLocusRecord, SsrPileupError> {
    let chrom_id = *name_to_id.get(r.chrom.as_ref())
        .ok_or(SsrPileupError::ContigNotInReference { name: r.chrom.into() })?;
    Ok(SsrLocusRecord { chrom_id, start: r.start, end: r.end, depth: r.depth,
        n_flanking: r.n_flanking, n_frr: r.n_frr, n_filtered: r.n_filtered,
        mapped_reads: r.mapped_reads, spanning: r.spanning })
}
```
(This resolves the *two `SsrLocusRecord` types* review item: the adapter is the
SNP-symmetric name→id boundary; the two types stay distinct by design.)

### Header build

Mirror `cli.rs::build_writer_header`: `ChromosomeEntry { name, length, md5 }` per
contig (md5 from `inputs.contigs` or recomputed), `subcommand: "ssr-pileup"`,
`input_crams`/`input_fasta` basenames, `parameters` (flank_bp, window, mapq,
catalog md5, reservoir cap). **Bind the catalog**: record the catalog's
`reference_md5` / path in `parameters` so a `.ssr.psp` is tied to its catalog.

---

## 5. The CLI (increment 3) — `ssr_pileup.rs`

`SsrPileupArgs`: `--reference`, `--catalog`, `--output`, `--sample` (or derive),
the alignment files (positional), `--trf-mod-path` n/a, mapq/dup filter knobs
(reuse `shared_args::Stage1Args` subset), `--window`, block-layout knobs.
`run_ssr_pileup` maps to `SsrPileupConfig` + stamps `date`, calls `driver::run`.
Wire `SsrPileup(SsrPileupArgs)` into `PopVarCallerCommand` + `main.rs` (exactly
like `ssr-catalog`).

---

## 6. Per-locus query vs per-contig sweep (perf decision)

**Chosen for v1: per-locus query** (arch §8 / `fetch_reads` doc) — one
`get_reads_from_segment` per locus. Pros: simple; natural per-locus reservoir +
seed + (later) per-locus parallel tasks. The pooled reader already removes the
per-call **file re-open + header re-parse** (the cost that motivated
`segment_reader`); what remains is one **index seek** per locus (a catalog can
be ~10⁵–10⁶ loci) — cf. the SNP `--regions` ~14% seek tax.

**Fallback if too slow (measure):** per-contig stream — `query(contig, None)`
once, sweep reads against the sorted loci. The postprocess **bundle drop**
guarantees loci are isolated by ≥ `bundle_threshold ≥ flank_bp`, so **each read
overlaps at most one locus window** → a clean sweep-line with no multi-locus
fan-out. Keep this in reserve; do not build it speculatively.

This is a measurement item, not a blocker — v1 ships per-locus, we profile on a
real catalog + BAM.

---

## 7. Decisions

**Resolved (in the fetcher, #1):**

1. ~~Window off-by-one~~ — **resolved**: the segment is the `ref_bytes` window,
   0-based half-open → 1-based inclusive (`+1` start, inclusive end); no contig
   clamp needed (the catalog already clamps `ref_bytes`). Covered by tests.
2. ~~Admission gate rule~~ — **resolved**: `reaches_locus` admits **spanning-capable
   reads only** (footprint brackets both ends) plus **always** any soft-clipped
   read. Flanking/in-repeat reads are *not* reservoir-admitted (so `n_flanking` /
   `n_frr` are a clipped-only subset — accepted; see §3).

3. ~~`QcCounts` semantics~~ — **decided** (see §4 `qc_counts`): all three columns
   count **independent primary alignments**. `depth = yielded`;
   `n_filtered = qc_fail + low_mapq + too_short + duplicate`;
   `mapped_reads = yielded + qc_fail + low_mapq + too_short` (dup-free coverage).
   Secondary / supplementary / unmapped excluded everywhere.

**Still open (pin while building the driver, #2):**

4. **`window` (rung) default** for `analyze_read` — a Stage-1 calibration param,
   separate from the catalog's flank.
5. **Adapter field mapping** — does the container `registry_ssr::SsrLocusRecord`
   store `n_spanning`, or derive it from `spanning.len()`? Verify when wiring
   `to_container_record` (§4).

---

## 8. Tests (the safety net, per increment)

- ~~**Fetcher**~~ — **DONE** (`696bf7e`): an inline indexed-BAM fixture asserts
  the window query + spanning gate (a spanning read admitted, a flanking read
  reach-skipped, a low-MAPQ read reader-filtered → `filtered.low_mapq`), and
  multi-file concatenation. The `reaches_locus` cases (spanning / one-end /
  buried / soft-clipped) are covered in `triage.rs`; the reservoir cap + seed
  determinism in the existing `Reservoir` tests. *Gap to add with the driver:* a
  deep locus (> cap) exercised through `fetch_locus_reads` (cap held, `yielded`
  reports the full count).
- **Driver** — end-to-end on a tiny synthetic reference + catalog + BAM with one
  clean spanning allele: `.ssr.psp` round-trips (read back via
  `records_of::<SsrKind>()`) to the expected `SsrLocusRecord`
  (chrom_id, coords, n_spanning, one profile). Reuses the container round-trip.
- **Adapter** — `to_container_record` maps name→id, errors on an unknown contig.
- **CLI** — `ssr-pileup --help` parses; the args→config mapping; a full run on the
  synthetic fixture produces a non-empty `.ssr.psp`.

Determinism gate (later, with the pool): byte-identical `.ssr.psp` across
`--threads ∈ {1, N}`, including the reservoir subsample.

---

## 9. What this plan deliberately defers

- The **fetcher-thread / bounded-queue / worker-pool / ordered-collector**
  topology (§8.4) — single-threaded first.
- **Off-ladder** alleles + columns (not generated yet).
- The **measured `count_repeats` fast path** (realign-everything is the default;
  the shortcut is a later, benchmarked option).
