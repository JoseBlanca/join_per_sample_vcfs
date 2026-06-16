# `.psp` container generalization — step 4 (`registry_ssr` + SSR schema) + step 5 (round-trip)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [psp_container_generalization.md](../../../doc/devel/implementation_plans/psp_container_generalization.md)
(architecture §10.4 of [ssr_genotyping_architecture.md](../../../doc/devel/architecture/ssr_genotyping_architecture.md))
**Builds on:** steps [1a](psp_container_generalization_step1a_2026-06-15.md) /
[1b](psp_container_generalization_step1b_2026-06-15.md) /
[2](psp_container_generalization_step2_2026-06-15.md) /
[3](psp_container_generalization_step3_2026-06-15.md)

## Goal

Add the **second `PspKind`** — the `ssr` schema — and a round-trip test
(step 5 folded in: the new serialization is untestable without write+read).
The SNP path stays byte-identical.

## The mapping (architecture §10.4, reconciled with all-CSR §11)

Stage-1's per-locus evidence (`Vec<Vec<(u16, f32)>>` profiles) fits the SNP
2-level container shape: **locus = record**, **spanning-read profile = entry**.
`n-spanning` is the per-record grouping count (exactly as `n-alleles` groups
alleles), and each profile carries two parallel per-profile CSR ragged lists —
`amb-lengths` (u16) + `amb-logliks` (f32). No `hist_*` (all-CSR); off-ladder
columns deferred. Every column lands on an existing `block.rs` codec.

`SsrLocusRecord` (container form) is **`chrom_id`-keyed** (mirroring how
`PileupRecord` relates to the block), not the Stage-1 chrom-name-keyed record —
the future Stage-1 driver adapts name↔id at the boundary, as for SNP.

## What changed

### New: [src/psp/registry_ssr.rs](../../../src/psp/registry_ssr.rs)
Self-contained SSR schema: `SSR_KIND`, `SSR_COLUMNS` (10 columns),
`SsrColumnKey` (+ `tag`/`from_tag`), `SsrLocusRecord`, `SsrKind: PspKind`,
`SsrBlock: BlockAccumulator`, `SsrDecoder: BlockDecoder`, and the round-trip
test.

### `ColumnDef` made schema-agnostic (the unblocker)
`ColumnDef` carried `key: ColumnKey` (the **SNP** dispatch enum), but
`SSR_COLUMNS` must also be `&[ColumnDef]` (it flows through
`build_header_bytes_for` / `columns_for_kind`). So the SNP-specific `key` field
blocked a second schema. **Removed `key` from `ColumnDef`**; each schema now
maps its own `tag → key` via a private `from_tag` (`ColumnKey::from_tag` for
SNP, `SsrColumnKey::from_tag` for SSR). The exhaustive `match Self` in each
key's `tag()` keeps the M4 guarantee. SNP encode/decode dispatch via
`ColumnKey::from_tag(def.tag)` — **byte-identical output** (verified).

### Block-header invariant relaxed (a SNP leak)
`validate_block_header_invariants` enforced `n_total_alleles >= n_records` —
a SNP semantic ("every record has ≥1 allele") that an SSR locus with **zero
spanning profiles** legitimately violates. Removed from the generic core; SNP
integrity is still enforced at `write_record` (zero-allele rejection) + the
read-side per-allele column-count checks. The `AllelesLessThanRecords` variant
is retained (pub API) but no longer raised.

### Writer ([writer.rs](../../../src/psp/writer.rs))
- Factored a schema-generic private `PspWriter::<W, S>::open` (header build via
  `S::KIND`/`S::columns()`, scratch sized to `S::columns().len()`); SNP's
  `new_with_block_layout` delegates to it (byte-identical), and new SSR
  constructors `new_ssr` / `new_ssr_with_block_layout` reuse it.
- `write_locus` + `validate_locus` on `impl<W> PspWriter<W, SsrKind>` (mirror
  of `write_record`; windows on locus `start`, block `last_pos` = max `end`).

### Reader ([reader.rs](../../../src/psp/reader.rs))
- `records_of::<S>()` — schema-typed sequential iterator (`records()` is the
  SnpKind case).
- New `pub(crate) read_and_inflate_column` (shared read+decompress, used by
  `SsrDecoder`); `read_compressed_blob` made `pub(crate)` (unknown-column skip).

### Wiring
- [registry.rs](../../../src/psp/registry.rs): `columns_for_kind` adds the
  `"ssr"` arm; `KNOWN_KINDS = "snp, ssr"`.
- [mod.rs](../../../src/psp/mod.rs): `pub mod registry_ssr`.

## Verification

- **SNP byte-identity:** a deterministic SNP `.psp` produced at step 3
  (`5d9134c`, via worktree) vs step 4 — `cmp` **identical**. The `ColumnDef`
  refactor + invariant removal are implementation-only for SNP.
- **SSR round-trip** (`ssr_round_trips_through_the_container`): write a
  multi-block `.ssr.psp` (multiple loci over windows, a 0-profile locus, a
  bimodal profile, a chromosome boundary) → reopen → `header().kind == "ssr"`,
  ≥4 blocks, and `records_of::<SsrKind>()` yields records **equal** to the
  input.

## Gates (dev container)

- `cargo build --lib` / `clippy --all-targets -- -D warnings` / `fmt --check` — clean
- `cargo test --lib` — **1138 passed, 0 failed, 1 ignored**
- SNP e2e: `psp_to_pileup_integration`, `pileup_cli_integration`,
  `cohort_cli_integration` — green

## Deviations from the plan (documented)

- **Two `SsrLocusRecord` types** (Stage-1 chrom-name-keyed vs container
  chrom_id-keyed), not the plan's "lean one type" — the container needs
  `chrom_id` (block-keyed, SNP-symmetric); the driver adapts. Documented on the
  container type.
- **`ColumnDef.key` removed** + **block invariant relaxed** — neither was
  anticipated as a touch point, but both are SNP leaks that hosting a second
  schema necessarily exposes; both verified SNP-byte-safe.

## Next

The §10.6 refactor is **complete** (steps 1a–5). Remaining SSR Stage-1 work
(separate from the container): the `fetch_reads` I/O driver, the Stage-1 driver
that produces `SsrLocusRecord`s and adapts them to the container's chrom_id
form, and the `ssr-pileup` CLI. Off-ladder columns + the measured fast path
remain deferred.
