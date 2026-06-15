# `.psp` container generalization — implementation sketch

**Status:** sketch for alignment, 2026-06-15. The *build plan* for architecture
§10 ([ssr_genotyping_architecture.md](../architecture/ssr_genotyping_architecture.md)
§10 is the architecture — written in-arch because it refactors existing shared
code). Goal: make [`src/psp/`](../../src/psp/) host **two schemas** (`snp`, `ssr`)
over one generic core, so Stage-1 `ssr-pileup` can write `.ssr.psp`. Files +
the `PspSchema` shape + the incremental step order — not full pseudocode.

**Regression gate, every step:** the SNP caller's end-to-end tests pass. *Not*
byte-identity to the old `.psp` (pre-alpha, no on-disk back-compat — arch §3.2).

Grounding (confirmed by reading the module, 2026-06-15):
- **Already generic — untouched:** [block.rs](../../src/psp/block.rs) (column
  codecs, `BlockHeader`, `ColumnManifestEntry`, zstd), [varint.rs](../../src/psp/varint.rs),
  [trailer.rs](../../src/psp/trailer.rs), and the columnar **read** layer
  (`BlockColumnReader` / `BlockColumns` in [reader.rs](../../src/psp/reader.rs)).
- **The SNP coupling, concentrated:** [registry.rs](../../src/psp/registry.rs)
  (`V1_0_COLUMNS` + the `ColumnKey` enum); [writer.rs](../../src/psp/writer.rs)
  (`BlockAccumulator` of SNP column `Vec`s + `append_record(&PileupRecord)`;
  `encode_column_into` = exhaustive `ColumnKey` dispatch; `encode_and_compress_columns`
  walks `V1_0_COLUMNS`; `WriterScratch.compressed` sized to `V1_0_COLUMNS.len()`;
  flush reads `block.delta_pos.len()` / `allele_seq_len.len()` for n_records /
  n_total_alleles); [reader.rs](../../src/psp/reader.rs) (two `V1_0_COLUMNS` walks +
  the `ColumnKey` decode match → `PileupRecord`).
- **One semantic tweak — [index.rs](../../src/psp/index.rs):** `BlockIndexEntry`
  already has `{chrom_id, first_pos, last_pos, block_offset}`; the change is what
  fills `last_pos` (SSR: `max(record.end)`) + an interval overlap test.
- **[header.rs](../../src/psp/header.rs):** already self-describing (`[[column]]`
  array); needs a `kind` tag + cross-check against the `kind`-selected registry.

---

## 1. The schema abstraction (the §10.7-Q1 fork — *decide at step 1*)

A schema supplies the four things the generic core can't know (arch §10.2):
the **registry table**, a **record type**, **record→columns** (encode), and
**columns→record** (decode). Sketch as a `PspSchema` trait:

```rust
/// A `.psp` schema family (snp | ssr): what the generic container core needs
/// beyond the schema-agnostic block / index / header machinery.
pub(crate) trait PspSchema {
    type Record;
    type Block: BlockAccumulator<Record = Self::Record>;

    const KIND: &'static str;                       // header tag "snp" | "ssr" (§10.3)
    fn columns() -> &'static [ColumnDef];           // the registry table

    /// Encode one column from a filled block to its uncompressed wire bytes —
    /// the per-schema exhaustive `ColumnKey` dispatch (today's encode_column_into).
    fn encode_column(def: &ColumnDef, block: &Self::Block, out: &mut Vec<u8>)
        -> Result<(), PspWriteError>;

    /// Materialise one record from decoded columns. Kept per-schema (each schema
    /// has its own closed decode); the shared substrate is the columnar reader
    /// beneath it (arch §10.2) — Stage-2 SSR rides that columnar layer, not this.
    fn record_from_columns(/* BlockColumns view + row index */) -> Self::Record;
}

/// The per-block accumulator a schema fills as records arrive: its column
/// buffers + the structural metadata the generic flush and block index need.
trait BlockAccumulator {
    type Record;
    fn new(chrom_id: u32, first_pos: u32) -> Self;
    fn append(&mut self, record: &Self::Record);
    fn last_pos(&self) -> u32;        // SNP: last record start; SSR: max(record.end) (§10.5)
    fn n_records(&self) -> u32;
    fn n_entries(&self) -> u32;       // per-allele total (SNP); the CSR entry total (SSR)
    fn projected_bytes(&self) -> usize;
    // chrom_id / first_pos are owned by the generic writer (it opens the block).
}
```

`PspWriter<W, S: PspSchema>` and the reader's typed iterator become generic over
`S`, with `SnpSchema` as the **only** impl in steps 1–3.

> **The fork to resolve in code (§10.7 Q1):** if the generics stay clean, keep
> the `PspSchema` trait. If they get unwieldy — `'static`-bound creep (house-style
> smell), or the decode side fighting the borrow of `BlockColumns` — fall back to
> **two concrete builders** (`SnpPspWriter` / `SsrPspWriter`) sharing the generic
> mechanics as **free functions** (flush, manifest, zstd, header). Decide from the
> *real* step-1 ergonomics, not up front.

**Keeping call sites stable:** give `S` a default (`PspWriter<W, S = SnpSchema>`)
and keep the public constructors (`new` / `new_with_block_target` /
`new_with_block_layout`) on the `S = SnpSchema` impl, so the ~19 existing
`PspWriter::new(...)` call sites + tests are untouched; `write_record(&S::Record)`
is `&PileupRecord` for them. (If the default-param inference fights us, a
`type SnpPspWriter<W> = PspWriter<W, SnpSchema>` alias is the fallback.)

---

## 2. Build sequence (arch §10.6 — each step green before the next)

1. **Parameterize the writer on a schema, behaviour-preserving (SNP-only).**
   Introduce `PspSchema` + `BlockAccumulator`; move the SNP `BlockAccumulator` /
   `append_record` / `encode_column_into` / `V1_0_COLUMNS` behind a `SnpSchema`
   impl; make `PspWriter<W, S>` generic with `S = SnpSchema` default + constructors
   pinned to SNP. **This is where the trait-vs-builders fork resolves.** SNP e2e
   green. *(Split: writer first (1a), then the reader's typed path (1b) — the two
   biggest files, each landed green.)*
2. **Add the `kind` tag** (§10.3). SNP writes `kind="snp"`; the reader selects the
   registry by `kind` and `cross_check_against_registry` validates against it
   (not the hardcoded `V1_0_COLUMNS`). SNP e2e green.
3. **Generalize the block-index overlap to intervals** (§10.5): `last_pos` =
   block max-`end`; region-overlap test against `[start, end)` (SNP = degenerate
   point, `end = start + 1`). SNP e2e green.
4. **Add `registry_ssr` + `SsrLocusRecord`** (§10.4) — new code, zero SNP impact
   (next section).
5. **Round-trip test** a synthetic `.ssr.psp`: write → columnar read → typed read.
   The first Bucket-1 container hook for Stage 1.

Steps 1–3 are the "extract from the live SNP consumer" work; 4–5 add the second
consumer. Nothing here needs the SSR math.

---

## 3. The SSR schema (step 4 — `registry_ssr` + `SsrLocusRecord`)

`registry.rs` → **`registry_snp.rs`** (rename via `git mv`; `V1_0_COLUMNS` + SNP
`ColumnKey`). New **`registry_ssr.rs`**: a `ColumnDef` table for the §4.3 columns
+ an `SsrLocusRecord` + its encode/decode. **No new codecs** — every SSR column
lands on a codec `block.rs` already has (arch §10.4):

- scalars (`start`, `end`, `depth`, `n_spanning`, …) → varint/scalar columns;
- the per-read `Qᵣ` CSR (`amb_read_offsets` / `amb_lengths` / `amb_logliks`) and
  any `offl_*` → the existing **CSR ragged-list** codecs;
- `chrom` / `offl_seqs` strings → the existing **bytes/dict** codecs.

> **Reconcile with the all-CSR storage decision (2026-06-15):** Stage 1's
> [`SsrLocusRecord`](../../src/ssr/pileup/locus_record.rs) dropped the histogram +
> weight and stores one renormalized `Qᵣ` profile per spanning read. So the SSR
> column table is the **`amb_*` CSR + QC scalars**, *not* §4.3's `hist_*` (and
> `n_flank_indel` is dropped). The container's `SsrLocusRecord` mapping flattens
> the in-memory `Vec<Vec<(u16, f32)>>` profiles into the CSR columns; `offl_*`
> stays empty until off-ladder wiring lands. (Amend spec §4.3 to match — the same
> follow-up the all-CSR + realign-everything decisions already owe.)

`SsrLocusRecord` placement: the in-memory one lives in `ssr/pileup/locus_record.rs`
today; the container's record↔columns mapping lives in `registry_ssr.rs`. Whether
they're one type (the container consumes the pileup record) or the writer adapts
the pileup record into a container record is a step-4 call (lean: one type, moved
to a shared peer if the back-reference rule demands it).

---

## 4. Files touched (flat layout — arch §10.5)

| file | change |
|---|---|
| `block.rs`, `varint.rs`, `trailer.rs` | none |
| `index.rs` | `last_pos` = max-`end`; interval overlap (point = degenerate) |
| `header.rs` | add `kind`; cross-check against the `kind`-selected registry |
| `writer.rs`, `reader.rs` | parameterize on `S: PspSchema`; SNP becomes one impl |
| `registry.rs` → `registry_snp.rs` | `git mv`; `V1_0_COLUMNS` + SNP `ColumnKey` |
| `registry_ssr.rs` (new) | SSR `ColumnDef` table + `SsrLocusRecord` mapping |

Directory stays flat; `core/`/`snp/`/`ssr/` submodules deferred to the §3.2
crowding trigger.

---

## 5. Deliberately left vague / deferred

- **Trait vs two builders** — resolved at step 1 from real ergonomics (§1).
- **Exact `SsrLocusRecord` type identity** (one shared type vs a container-side
  record the writer adapts to) — step-4 call.
- **`offl_*` columns** — schema slots defined, populated only when off-ladder
  candidate generation lands (empty for now).
- **The Stage-1 driver + fetcher I/O** that *uses* the writer — separate pass;
  this plan stops at "the writer/reader can round-trip an `SsrLocusRecord`."
- **Spec §4.3 amendment** (all-CSR, drop `hist_*`/`n_flank_indel`) — a doc
  follow-up owed jointly with the realign-everything + all-CSR decisions.
