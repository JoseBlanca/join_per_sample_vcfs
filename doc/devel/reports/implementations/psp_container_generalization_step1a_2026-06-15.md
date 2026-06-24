# `.psp` container generalization — step 1a (writer on a `PspKind` trait)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [psp_container_generalization.md](../../../doc/devel/implementation_plans/psp_container_generalization.md)
(architecture §10 of [ssr_genotyping_architecture.md](../../../doc/devel/architecture/ssr_genotyping_architecture.md))

## Goal

Step 1a of the container generalization: parameterize the `.psp` **writer**
on a `PspKind` trait so a second schema (`ssr`) can be added later without
touching the block/index/header/zstd machinery. SNP is the only `PspKind`
impl; the change is **behaviour-preserving** (regression gate = the SNP
end-to-end tests, not byte-identity — pre-alpha, no on-disk back-compat).
The decode side (reader's typed path) is deliberately left to step 1b.

## What changed

### New: [src/psp/kind.rs](../../../src/psp/kind.rs)
Two `pub` traits:

- **`PspKind`** — `type Record`, `type Block: BlockAccumulator<Record = …>`,
  `const KIND: &'static str` (header tag, consumed in step 2), `fn columns()
  -> &'static [ColumnDef]`, and `fn encode_column(def, &block, out)`. No
  `record_from_columns` yet — that is the step-1b decode half.
- **`BlockAccumulator`** — `new_block`, `append`, `chrom_id`, `first_pos`,
  `last_pos`, `n_records`, `n_entries`, `projected_bytes`.

### [src/psp/writer.rs](../../../src/psp/writer.rs)
- Private struct `BlockAccumulator` → `pub struct SnpBlock`; its inherent
  `new` → trait `new_block`, `append_record` → trait `append`; added the
  accessors (`chrom_id`/`first_pos`/`last_pos` return the fields,
  `n_records` = `delta_pos.len()`, `n_entries` = `allele_seq_len.len()`,
  `projected_bytes` returns the field). `impl BlockAccumulator for SnpBlock`.
- New `pub struct SnpKind; impl PspKind for SnpKind` — `columns()` =
  `V1_0_COLUMNS`, `KIND = "snp"`, `encode_column` delegates to the renamed
  `encode_snp_column_into` (the old `encode_column_into`, still folding in
  the M5 `predict_uncompressed_len` self-check).
- `PspWriter<W: Write>` → `PspWriter<W: Write, S: PspKind = SnpKind>`, with
  `block: Option<S::Block>`. The default keeps every `PspWriter<W>` /
  `PspWriter::new(...)` call site (~19, plus tests) unchanged.
- The three flush-phase free functions (`encode_and_compress_columns`,
  `assemble_block_header`, `emit_block_to_sink`) are now generic over `S`,
  driven by `S::columns()` / `S::encode_column` / the `BlockAccumulator`
  accessors. `WriterScratch::new(n_columns)` sizes its per-column buffers
  to `S::columns().len()` (the constructors pass `SnpKind::columns().len()`).

### [src/psp/mod.rs](../../../src/psp/mod.rs)
- `pub mod kind;`
- `pub use registry::ColumnDef;` — `ColumnDef` appears in the `pub`
  `PspKind` signatures, so it must be reachable at `pub` visibility. Done
  via a single re-export rather than widening the whole `pub(crate)`
  registry module.

## Fork resolution (arch §10.7 Q1 — trait vs two builders)

**Kept the `PspKind` trait.** The generics stayed clean (no `'static`
creep, no decode-borrow fights — the decode side isn't in this step). The
one nuance: rather than making `write_record` itself generic, the split is

- **record ingest + validation** (`write_record`, `validate_record`,
  `apply_record_to_block`, the `should_flush` decision) stays on the
  concrete `impl<W: Write> PspWriter<W, SnpKind>` — it is intrinsically
  `PileupRecord`-shaped (ACGTN/`q_sum`/chain-id rules, `(chrom_id, pos)`
  extraction), and the listed `BlockAccumulator` surface deliberately has
  no validate/coordinate methods;
- **the column flush/encode machinery** (`finish`, `flush_block`, the three
  phase functions) is generic over `S` — this is the actual refactor that
  stops hardcoding `V1_0_COLUMNS`.

SSR will add its own inherent ingest path (`write_locus`-style) on
`impl<W> PspWriter<W, SsrKind>` and reuse the same generic flush machinery.
This is the trait + shared-generic-functions hybrid; it satisfies both the
"keep `PspKind`" lean and the minimal trait surface.

**Visibility note:** because `psp` is `pub mod` and `PspWriter` is reachable
by external bench/example crates, the trait machinery it bounds on must also
be `pub` — so `PspKind`, `BlockAccumulator`, `SnpKind`, `SnpBlock`, and the
re-exported `ColumnDef` are all `pub`. (`SnpBlock`'s fields stay private.)
The compiler's `private_bounds` / `E0446` enforced this; it is not optional.

## Deviation from the plan sketch

The plan's `BlockAccumulator` sketch commented that `chrom_id` / `first_pos`
are "owned by the generic writer." Implemented instead as **accessors on the
trait**: the accumulator already holds them (SNP needs `first_pos` for delta
encoding), so exposing them avoids duplicating that state into the writer as
a parallel `Option`. Net simpler; documented in `kind.rs`.

## Gates (all in the dev container)

- `cargo build --lib` — clean
- `cargo clippy --all-targets -- -D warnings` — clean
- `cargo fmt --check` — clean
- `cargo test --lib` — **1133 passed, 0 failed, 1 ignored**
- SNP e2e: `psp_to_pileup_integration` (4), `pileup_cli_integration` (10),
  `thread_budget_integration` (1) — all green.

## Next (remaining §10.6 steps)

- **1b** — reader's typed path on the schema (the `record_from_columns`
  half + registry-by-kind selection for the columnar reader).
- **2** — `kind` header tag (`SnpKind::KIND` is already wired; the header +
  `cross_check_against_registry` consume it).
- **3** — interval-overlap block index (`last_pos` = block max-`end`).
- **4** — `registry_ssr` + `SsrLocusRecord` (all-CSR: `amb_*` CSR + QC
  scalars, no `hist_*`).
- **5** — synthetic `.ssr.psp` round-trip test.
