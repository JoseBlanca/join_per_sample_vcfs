# SSR Stage 2 — `ssr-call` reading & merge (implementation plan)

**Status:** draft, 2026-06-21, branch `ssr-cohort`. **Phase 1 of Stage 2** (`ssr-call`):
read the N per-sample `.ssr.psp` files, k-way-merge by catalog locus into the
`CohortLocus` work-item, and feed it through a producer→queue→worker→writer topology.
The EM/statistics (Phases 2–3) are **out of scope** here — this plan stops at *a
deterministic, catalog-ordered stream of `CohortLocus` reaching a worker stub*.

**Specs & architecture this implements:**
- [spec §4.1](../specs/ssr_cohort_mark2.md) — reading & orchestration *intent*
  (settled 2026-06-19). **Where the spec's topology prose predates the architecture
  refinement (it still says "batches of K", "dedicated EM pool", `fetch`), the
  architecture doc wins** — those were resolved this session.
- [architecture: ssr_call_reading.md](../architecture/ssr_call_reading.md) — the
  settled module/struct shape (the authoritative design for this plan).

**Companion plans (not yet written):** parameters pre-pass and genotyping/EM —
the consumers of this stream (arch [parameters](../architecture/ssr_call_parameters.md) /
[genotyping](../architecture/ssr_call_genotyping.md)).

**Build philosophy:** incremental — each phase compiles, tests, and is reviewable on
its own; pause between phases. The reading layer is built **synchronous-first**
(inline decode), and the shared decode pool + prefetch (arch §3/§5) is a **later,
profiling-gated** phase, not designed into the first cut.

---

## 0. What exists to build on (verified surfaces)

| surface | path | what we use |
|---|---|---|
| `PspReader<R>` | [src/psp/reader.rs](../../../src/psp/reader.rs) | `new()`, `header()`, `block_index() -> &[BlockIndexEntry]`, `into_column_blocks() -> BlockColumnReader`, `records_of::<SsrKind>()` |
| `BlockIndexEntry` | [src/psp/index.rs](../../../src/psp/index.rs) | `{ chrom_id: u32, first_pos: u32, last_pos: u32 (inclusive), n_records: u32, block_offset: u64 }`; overlap = `first_pos ≤ end && last_pos ≥ start` |
| `BlockColumnReader<R>` / `BlockColumns<'a>` | [src/psp/reader.rs](../../../src/psp/reader.rs) | `peek_block()`, `load_current()`, `columns()`, `seek_to(chrom,start,end)`, `advance()` — explicit one-block control |
| `SsrKind` / `SsrDecoder` / `SsrLocusRecord` | [src/psp/registry_ssr.rs](../../../src/psp/registry_ssr.rs) | `SsrDecoder::decode_block(...)` + `next_record() -> Option<Result<SsrLocusRecord>>`; record = `{ chrom_id, start (1-based), end (excl), depth, n_filtered, mapped_reads, n_low_quality, n_border_off_end, observed: Vec<(Box<[u8]>, u32)> }` |
| `SsrLocusObs` / `SsrQc` | [src/ssr/pileup/locus_tally.rs](../../../src/ssr/pileup/locus_tally.rs) | the Stage-1 in-memory mirror `SampleEvidence` shadows (chrom-**name**, 0-based) |
| `CatalogReader<R>` / `Locus` / `Motif` | [src/ssr/catalog/io.rs](../../../src/ssr/catalog/io.rs), [src/ssr/types.rs](../../../src/ssr/types.rs) | `read_locus() -> Option<Result<Locus>>` (sequential, file-order = coord-order), `header()` (md5); `Locus` = chrom-**name**, **0-based half-open**, `motif`, `ref_bytes`/`ref_tract()`/flanks |
| SNP cohort topology precedent | [src/var_calling/pipeline.rs](../../../src/var_calling/pipeline.rs), [vcf_writer.rs](../../../src/var_calling/vcf_writer.rs) | `crossbeam_channel::bounded`, `rayon::ThreadPool`, writer reorder via `BTreeMap<u64, _>` by seq |
| CLI wiring | [src/pop_var_caller/cli.rs](../../../src/pop_var_caller/cli.rs), [src/main.rs](../../../src/main.rs) | add `SsrCall(SsrCallArgs)` to `PopVarCallerCommand` + dispatch (mirrors `SsrPileup`) |

**Does not exist yet:** any SSR *cohort* reader/merger; a catalog coordinate index
(only sequential `read_locus()` — fine, we drive off it, Q-R2); a "shared decode pool"
abstraction (Phase 5 builds it).

---

## 1. Two reconciliations to settle before coding (correctness gotchas)

These are not in the arch doc because they're code-level; both are load-bearing.

1. **Coordinate frames disagree.** Catalog `Locus` is **0-based half-open**, chrom by
   **name**; the container `SsrLocusRecord`/`BlockIndexEntry` is **1-based**, chrom by
   **id** (per-file dictionary). The cohort-match key `LocusId` (arch §2: `(chrom,
   start, end)`) must be defined in **one** frame and every comparison converted to it.
   **Decision:** `LocusId` carries a **global chrom id** + the **catalog's 0-based
   half-open** `(start, end)`; each cursor converts its file's records into that frame
   (per-file `chrom_id` → name via the psp header dict → global id; `+1`/`−1` offset).
   The conversion lives at the cursor's decode boundary so the merger only ever sees
   the catalog frame.

2. **`observed` vs `seq_counts`.** The container decodes to `SsrLocusRecord.observed`;
   the arch doc renamed the cohort field to `SampleEvidence.seq_counts`. **Decision
   (default):** keep the container/Stage-1 name `observed` untouched; the cursor's
   decode adapter maps `record.observed → SampleEvidence.seq_counts`. (Renaming
   Stage-1 too is a separate, optional cleanup — out of scope here.)

> Resolve both in Phase 0's `types.rs` doc-comment so later phases inherit one frame.

---

## 2. Phase plan (dependency order)

### Phase 0 — scaffolding: module, CLI slot, core types
- New module `src/ssr/cohort/` (arch §9): `mod.rs`, `types.rs`.
- `types.rs`: `LocusId` (global chrom id + 0-based half-open `start,end`),
  `CohortLocus` (sparse SoA — `present: Vec<u32>` + `samples: Vec<SampleEvidence>`,
  arch §2 / Q-R1), `SampleEvidence { seq_counts: Vec<(Box<[u8]>, u32)>, qc: SsrQc }`.
- CLI: `SsrCallArgs` (inputs: N `.ssr.psp` paths + `.ssr.catalog`; `--threads`;
  output VCF path; `--queue-depth` knob, Q-R3), register `SsrCall` in
  `PopVarCallerCommand`, dispatch `run_ssr_call` in `main.rs` (stub returning
  `Ok` for now). Mirrors `SsrPileup`.
- **Tests:** type round-trips; CLI parses. **Exit:** `cargo build` + `ssr-call --help`.

### Phase 1 — `SampleEvidenceCursor` (synchronous, inline decode)
The heart. Per sample (arch §3). **No pool yet** — decode inline in `advance()`.
- **Open:** `PspReader::new` per file; verify catalog md5 against the `CatalogReader`
  header (hard error on mismatch, spec §4.1); load `block_index()`; build the
  per-file `chrom_id ↔ global` map from the psp header dict.
- **State:** block index, a `BlockColumnReader` (or direct `SsrDecoder` over blocks),
  within-block position, `held: Option<(LocusId, SampleEvidence)>`,
  `last_query: Option<LocusId>`.
- **`evidence_at(q: LocusId) -> Option<SampleEvidence>`** — the settled contract
  (arch §3): `debug_assert!(last_query < q)` monotonic guard; then match `q` vs
  `held`: `==` → take + `advance()` + `Some`; `<` → `None` (Absent); `>` → `panic`
  ("merger skipped a stored locus"). Exhausted (`held == None`) → `None`.
- **`advance()`** — move to next stored record, converting frames (§1); decode it into
  `held`; cross blocks via the block index when needed (drop current, decode covering
  block — inline for now). No more records → `held = None`. Called from `new` (preload)
  and after each hit.
- **Decode adapter:** `SsrLocusRecord → SampleEvidence` (`observed → seq_counts`,
  QC scalars → `SsrQc`), with coordinate conversion.
- **Tests (synthetic `.ssr.psp` fixtures — mirror the Stage-1 fixture helpers):**
  ascending hits; Absent on a skipped locus (`q < held`); exhaustion → permanent
  `None`; rewind → panic; skip (`q > held`) → panic; multi-block file crosses
  boundaries correctly; catalog-md5 mismatch → error.

### Phase 2 — catalog-driven k-way merger (one `CohortLocus` at a time)
- `merge.rs`: drive off `CatalogReader::read_locus()` (catalog-driven, Q-R2); for each
  catalog `Locus`, build its `LocusId`, call `evidence_at` on **every** cursor, gather
  the `Some`s into a `CohortLocus` (record `present[k]` = cursor index), attach the
  catalog frame (`motif`, `ref_tract`); **≥1 present → emit; all absent → drop**
  (sparse-omit). Tag each emitted `CohortLocus` with a **monotonic locus seq**.
- **One at a time** onto the consumer — **no K-batching** (arch §4; deferred).
- **Tests:** dense cohort (every sample present) → every catalog locus emitted in
  order; sparse fixture → all-absent loci omitted, seq numbers monotonic; determinism
  (same stream regardless of how many cursors).

### Phase 3 — driver: producer + bounded queue + worker stub + writer
- `driver.rs` (arch §5): producer thread runs the Phase-2 merge, pushes one
  `CohortLocus` per item onto a `crossbeam_channel::bounded(queue_depth)`; a worker
  pool (`rayon::ThreadPool`) pops items and runs a **stub** (Phases 2–3 EM lands
  later) that just passes the locus through with its seq; a writer thread reorders by
  **locus seq** (`BTreeMap<u64, _>`, mirroring `vcf_writer.rs`) and emits in catalog
  order.
- This makes the reading layer **end-to-end testable** before the EM exists (e.g. the
  stub emits a TSV/debug dump of the merged stream).
- **Tests:** output order == catalog order across `--threads 1..K` (determinism, arch
  §7 — the stub must be order-pure); back-pressure (small queue depth) doesn't
  deadlock or drop loci.

### Phase 4 — two-pass restart (re-read)
- The pre-pass consumes the merge stream, then Phase 3 consumes it again (arch §8,
  Q-R5: **re-read, not cache**). Make the merge layer **cheap to restart**: a
  `restart()` that re-seeks every cursor from its block index (re-`new`/re-`advance`),
  no global state to rebuild. Expose the merge as a restartable iterator the driver
  can run twice.
- **Tests:** two full passes over the same inputs yield identical streams; restart
  allocates no new file handles beyond reopening (or reuses them).

### Phase 5 — shared decode pool + prefetch (PROFILING-GATED, arch §3/§5, Q-R4)
**Build only when profiling shows the producer is decode-bound** (interacts with
Q-R6). Converts inline decode → offloaded:
- `decode_pool.rs`: one shared pool across all cursors; a decode task is a **stateless
  pure transform** (block's compressed bytes + column set → decoded `SsrLocusRecord`s),
  result returned on a `oneshot`/crossbeam receiver (the "future").
- Cursor gains `current` / `next` `BlockSlot { Decoded | Pending(Receiver) | Empty }`;
  **self-double-buffers** — on entering `current`, submit `next`'s decode; on crossing,
  `recv()` (usually ready), promote, re-prefetch. Blocking `recv` = rare degraded path.
- Pool is **shared with the EM workers, decode-priority** (arch §5 — prefetch makes
  decode latency-insensitive, so no rigid split; dodges the mis-split penalty).
- **Decision deferred to build time:** how the decode task reads the block's compressed
  bytes — a per-worker file handle (reopen / `Arc<Mmap>`) vs the producer pre-reading
  raw bytes and shipping them. Pick against the `BlockColumnReader`/`SsrDecoder` I/O
  shape then.
- **Memory:** `current + next` = N × 2 blocks resident (vs N × 1); tunable by prefetch
  depth (arch §7).
- **Tests:** byte-identical stream vs Phase-1 inline decode (decode placement must not
  change content, only timing); determinism across threads preserved.

---

## 3. Testing & invariants (cross-phase)

- **Determinism is the headline invariant** (arch §7): the `CohortLocus` stream and the
  VCF order are a pure function of inputs + catalog order, identical across `--threads`.
  Every phase that touches threading re-asserts it.
- **Synthetic fixtures**, Bucket-1 style (mirror the Stage-1 `.ssr.psp` test helpers):
  cover dense, sparse, single-sample, multi-block-per-sample, and exhausted-cursor
  cases. Use them for the cursor panics (rewind/skip) too.
- **Stage-1 writer invariant we depend on** (arch §3, spec §4.1): a locus lives wholly
  inside one block — never split across a boundary. Add a **debug assertion / fixture
  check** so a future Stage-1 change that violates it fails loudly here.
- **Lockstep memory** (arch §7): assert resident decompressed blocks ≈ N (Phase 1–4) or
  N × 2 (Phase 5), via a counting test, not just prose.

---

## 4. Open decisions carried in (from arch §10), and their build status

| arch Q | status for this plan |
|---|---|
| **Q-R1** sparse SoA | settled — `present: Vec<u32>` + `samples` (Phase 0) |
| **Q-R2** catalog-driven | settled — drive off `CatalogReader` (Phase 2) |
| **Q-R3** queue depth | a CLI knob (`--queue-depth`), default measured (Phase 3) |
| **Q-R4** shared decode-priority pool + prefetch | design settled; **build = Phase 5, profiling-gated** |
| **Q-R5** re-read vs cache | settled — re-read (Phase 4) |
| **Q-R6** genomically-aligned Stage-1 blocks | **out of scope** — a Stage-1-writer follow-up; this plan only adds the in-block invariant check |
| §1 coord-frame + `observed`→`seq_counts` | **new code-level decisions** — settle in Phase 0 `types.rs` |

---

## 5. Out of scope (explicit)

- The EM, candidate assembly, HipSTR likelihood, `F` loop, VCF *content/genotypes*
  (Phases 2–3 — companion plans). Phase 3's worker is a **stub**.
- Parameter pre-pass *logic* (parameters doc); this plan only provides the
  **restartable merge stream** it will consume (Phase 4).
- Genomically-aligned Stage-1 blocks (Q-R6) and any Stage-1 rename of `observed`.
