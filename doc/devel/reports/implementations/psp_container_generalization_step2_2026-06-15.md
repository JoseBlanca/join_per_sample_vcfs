# `.psp` container generalization — step 2 (the `kind` header tag)

**Date:** 2026-06-15
**Branch:** `ssr-architecture`
**Plan:** [psp_container_generalization.md](../../../doc/devel/implementation_plans/psp_container_generalization.md)
(architecture §10.3 of [ssr_genotyping_architecture.md](../../../doc/devel/architecture/ssr_genotyping_architecture.md))
**Builds on:** [1a](psp_container_generalization_step1a_2026-06-15.md),
[1b](psp_container_generalization_step1b_2026-06-15.md)

## Goal

Add the `kind` schema-family tag to the header (§10.3). The writer emits
`kind = S::KIND`; the reader selects the column registry by `kind` and
`cross_check_against_registry` validates the `[[column]]` array against *that*
registry rather than the hardcoded `V1_0_COLUMNS`.

## Format change — byte-identity intentionally broken

Unlike 1a/1b (pure refactors, byte-identical output), **step 2 changes the
produced `.psp`**: the TOML header gains a `kind = "snp"` line. This is a
deliberate format addition (pre-alpha, no on-disk back-compat — arch §3.2), so
the gate is **round-trip + SNP e2e green + back-compat**, not byte-identity. The
change is surgical: the only header delta is the one `kind` line; pre-tag files
(no `kind` key) still read as the SNP default.

## What changed

### [src/psp/registry.rs](../../../src/psp/registry.rs)
- `SNP_KIND = "snp"` (single source of truth for the tag), `KNOWN_KINDS` (for
  diagnostics), and `columns_for_kind(kind) -> Option<&'static [ColumnDef]>` —
  the kind→registry selector (`"snp"` → `V1_0_COLUMNS`; `None` otherwise).

### [src/psp/errors.rs](../../../src/psp/errors.rs)
- New `PspReadError::UnknownKind { kind, known }`.

### [src/psp/header.rs](../../../src/psp/header.rs)
- `WireHeader` gains `kind: String` (after `format-version`), with
  `#[serde(default = "default_kind")]` → `"snp"` so pre-tag files parse.
- `ParsedHeader` gains `pub kind: String`.
- Build path parameterized: `build_header_toml_for` / `build_header_bytes_for`
  take `(kind, columns)`; the no-arg `build_header_toml` / `build_header_bytes`
  stay as SNP-default wrappers (the ~20 test call sites are untouched).
  `wire_from_writer_header` now emits the passed `kind` + `columns`.
- Parse path: selects the registry via `columns_for_kind(&wire.kind)` →
  `UnknownKind` on miss; threads it into `cross_check_against_registry(columns,
  registry)` (no longer hardcoded to `V1_0_COLUMNS`); records `kind` on
  `ParsedHeader`.

### [src/psp/writer.rs](../../../src/psp/writer.rs)
- `SnpKind::KIND = registry::SNP_KIND` (ties the write-side tag to the
  registry's single source of truth).
- The constructor now calls `build_header_bytes_for(&header, SnpKind::KIND,
  SnpKind::columns())` — schema-faithful, so an `SsrKind` constructor (step 4+)
  drops in by passing its own `KIND`/`columns()`.

## Tests added (header.rs)

- `header_carries_kind_snp_and_round_trips` — body declares `kind = "snp"`,
  round-trips into `ParsedHeader::kind`.
- `unknown_kind_is_rejected` — a mutated `kind` → `UnknownKind`.
- `missing_kind_defaults_to_snp` — a header with the `kind` line stripped parses
  as the SNP default (pre-tag back-compat).

## Gates (dev container)

- `cargo build --lib` / `clippy --all-targets -- -D warnings` / `fmt --check` — clean
- `cargo test --lib` — **1136 passed, 0 failed, 1 ignored** (+3 new)
- SNP e2e: `psp_to_pileup_integration`, `pileup_cli_integration`,
  `cohort_cli_integration` — green

## Next (remaining §10.6 steps)

- **3** — generalize the block-index overlap to intervals (`last_pos` = block
  max-`end`; SNP = degenerate point).
- **4** — `registry_ssr` + `SsrLocusRecord` (all-CSR: `amb_*` CSR + QC scalars).
- **5** — synthetic `.ssr.psp` round-trip test.
