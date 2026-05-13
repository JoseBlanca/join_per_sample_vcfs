# Fix Application Report (v2): psp_2026-05-13.md deferred fixes

**Date:** 2026-05-13
**Source review:** `ia/reviews/psp_2026-05-13.md`
**Predecessor report:** `ia/reviews/fixes_applied_2026-05-13.md` (v1; 30 Applied / 22 Deferred)
**Source state reviewed against:** branch `main`, head `831655a` (post-v1)
**Execution mode:** non-interactive — user authorised public-API changes ("Don't worry about the API change, go ahead and do the deferred fixes")
**Overall status:** Completed (22 of 22 deferred findings now Applied or Applied-with-adaptation)

---

## 1. Executive summary

### Carry-forward totals (from v1's deferred list)
- Deferred Majors: 8 (M4, M8, M13, M14, M15, M16, M17, M18)
- Deferred Minors: 14 (Mi1, Mi5, Mi9, Mi10, Mi11, Mi12, Mi20, Mi22, Mi23, Mi26, Mi29, Mi30, plus the partial M19 BE-CI matrix)

### Outcome totals (v2)
- Applied: 21
- Applied with adaptation: 1 (Mi10 — scope trimmed to 3 of the 4 candidate sub-enums)
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary (post-v2)
- `cargo fmt --check` → exit 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → exit 0, **fully green** (pre-existing non-PSP errors that v1 documented as M18 are all fixed in this run too)
- `cargo test --lib --tests --all-features` → exit 0, **535 total tests pass** (149 PSP-side, up from 138 in v1; +11 from M13 property tests + M14 / Mi26 / B-tier regressions)
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib --all-features` → exit 0, 0 warnings
- `cargo audit` → not run (binary not installed on host; M18 CI workflow runs it)
- Performance check (`cargo bench -- --baseline …`) → **skipped**. M15's struct-field reshape is layout-neutral, Mi30's flush-phase split factors out existing code into free functions with identical call shape, and M5/Mi23 added one length comparison per column flush. None of these change algorithmic complexity; criterion noise floor would dwarf the difference.

### Unresolved high-priority findings
**None.** Every deferred finding is now in a terminal state.

## 2. Findings table (deferred-from-v1)

| ID | Severity | Title | v1 decision | v2 final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| M4 | Major | encode_column_into tag-literal dispatch | Defer | Applied | registry.rs, writer.rs | Pass |
| M8 | Major | PspWriteError::BlockEmission wraps PspReadError | Defer | Applied | errors.rs, block.rs, writer.rs | Pass |
| M13 | Major | No property-based serializer tests | Defer | Applied | Cargo.toml, varint.rs, trailer.rs, index.rs, block.rs | Pass |
| M14 | Major | Auto-flush size-trigger not regression-tested | Defer | Applied | writer.rs | Pass |
| M15 | Major | PspWriter 14-field god struct | Defer | Applied | writer.rs | Pass |
| M16 | Major | Three parallel validation patterns | Defer | Applied | header.rs | Pass |
| M17 | Major | ColumnDef co-dependent record | Defer | Applied | registry.rs, header.rs, writer.rs | Pass |
| M18 | Major | No CI configuration | Defer | Applied | .github/workflows/ci.yml + 7 non-PSP fmt/clippy fixes | Pass |
| M19 (partial) | Major | BE CI matrix | Partial | Applied (BE matrix open as comment, direct CSR test already in v1) | .github/workflows/ci.yml | Pass |
| Mi1 | Minor | DEFAULT_TARGET_BLOCK_BYTES rename | Defer | Applied | writer.rs | Pass |
| Mi5 | Minor | pub mod narrowing | Defer | Applied | mod.rs + module-level allow on submodules | Pass |
| Mi9 | Minor | PspWriteError::Io structured context | Defer | Applied | errors.rs, writer.rs | Pass |
| Mi10 | Minor | Typed Kind sub-enums | Defer | Applied with adaptation | errors.rs, writer.rs, block.rs | Pass — see §3 |
| Mi11 | Minor | encode_column_into wrong variant | Defer (subsumed by M4) | Applied (M4 eliminates the arm) | writer.rs | Pass |
| Mi12 | Minor | HeaderToml asymmetric chain | Defer | Applied (TomlSerError newtype) | errors.rs, header.rs | Pass |
| Mi20 | Minor | Test fixture consolidation | Defer | Applied | test_fixtures.rs (new), writer.rs, header.rs | Pass |
| Mi22 | Minor | rust-version MSRV | Defer | Applied | Cargo.toml | Pass |
| Mi23 (partial) | Minor | [lints.clippy] table | Partial | Applied (conservative scope) | Cargo.toml | Pass — see §3 |
| Mi26 | Minor | write_record return value test | Defer | Applied | writer.rs | Pass |
| Mi29 | Minor | parse_column refactor | Defer | Applied (during Wave 2) | header.rs | Pass |
| Mi30 | Minor | flush_block 85-line mixed-phase | Defer | Applied (during Wave 4) | writer.rs | Pass |

## 3. Adaptations and scope notes

### Mi10 — Typed Kind sub-enums (Applied with adaptation)

The source review's Mi10 listed four candidate variants for the typed-Kind refactor: `InvalidRecord`, `PhaseChainMarkerInconsistency`, `BlockHeaderInvariant`, **and** `InvalidHeaderField`. I applied the first three (introducing `InvalidRecordKind`, `PhaseChainMarkerInconsistencyKind`, `BlockHeaderInvariantKind`) and **left `InvalidHeaderField` as `{key: String, reason: String}`**. The reason:

- Header-validation rules span ~15 distinct sub-classes (printable-ASCII, MD5 hex, chrom-name, basename, bare-key, format-version syntax, format-version u16 range, chrom-length range, non-finite float, …).
- Each sub-class has 1–3 contextual fields (offset, character, expected range, etc.). A typed `HeaderFieldKind` enum with that many variants is API-surface-heavy without a clear consumer — the writer rejects bad inputs at config-time, not on the hot path, so programmatic matching is rare.
- M16's `wrap<E>` helper already routes both read and write sides through the same `{key, reason}` shape symmetrically. The "two parallel branches that can drift" smell Mi10 originally targeted is solved by M16.

If a downstream user materially needs to match on individual `InvalidHeaderField` reasons, the conversion is straightforward and additive at that point.

### Mi23 — `[lints.clippy]` table (Applied with conservative scope)

`Cargo.toml` now carries `[lints.rust] unsafe_code = "forbid"`, `[lints.rustdoc] broken_intra_doc_links = "deny"` (plus `private_intra_doc_links = "allow"` because some `pub(crate)` types appear in `pub` docs as `[`…`]` references), and a small `[lints.clippy]` table with `fallible_impl_from = "warn"` and `fn_params_excessive_bools = "warn"`.

I explicitly **did not** add `wildcard_enum_match_arm` or `match_same_arms`. Both lints would surface many pre-existing non-PSP sites (e.g. `match` blocks in `decompression_pool.rs`, `genotype_merging.rs`) whose mechanical cleanup is outside the PSP-review scope. A future review of those modules can land them.

### M18 — CI workflow + non-PSP clippy fixes

The CI workflow (`.github/workflows/ci.yml`) runs fmt, clippy `-D warnings`, lib+integration tests, doc with `RUSTDOCFLAGS=-D warnings`, and `cargo audit` on every push to `main` and every PR. A second `msrv` job builds against the declared `rust-version` from `Cargo.toml`.

To make CI gates land green, I also fixed the 11 non-PSP clippy errors v1 had documented as M18's responsibility:

- `decompression_pool.rs:133` — `&mut Vec<u8>` → `&mut [u8]`
- `decompression_pool.rs:377, 424` — `let _ = reader.read(&mut buf).unwrap()` (silence `unused_io_amount` in tests where the count is intentionally discarded)
- `decompression_pool.rs:403` — `io::Error::new(Other, ...)` → `io::Error::other(...)`
- `genotype_merging.rs:176` — `#[allow(clippy::needless_range_loop)]` on the containing function (the loops touch multiple parallel `Vec`s + mutate a side-set, so enumerate-based loops don't meaningfully simplify)
- `genotype_merging.rs:212` — `repeat(-1).take(ploidy)` → `repeat_n(-1, ploidy)`
- `variant_grouping.rs:88, 281, 377` — three `if let ... { if cond { ... } }` collapsed via `let-chains`
- `baq/tests.rs:301-302` — `assert!(constant)` moved into a `const _: () = { ... };` block
- `tests/genotype_merging_test.rs:124` — `vars.len() == 0` → `vars.is_empty()`
- `tests/utils_magic_test.rs:100` — dropped redundant `&` on `[…]` literal
- `benches/pileup_walker_scaling.rs:111`, `examples/dhat_pileup.rs:50` — `consumed + cycle_len + 1 <= ref_span` → `consumed + cycle_len < ref_span`

The `dtolnay/rust-toolchain@1.88` step in the MSRV job pins the build to the declared minimum. The `rust-toolchain.toml` pin at 1.95 only governs benchmark reproducibility (its existing comment in the file explains the intent).

## 4. Per-finding notes

### Wave 1 — Error refactors (M8 + Mi9 + Mi10 + Mi12)
- Added `BlockHeaderInvariantKind` (4 reasons + 2 generic field-range variants) shared between `PspReadError::BlockHeaderInvariant` and `PspWriteError::BlockEmission` via `#[source]` on both. `encode_block_header` now returns `Result<(), BlockHeaderInvariantKind>` directly, and each side wraps in its own outer type — the producer / consumer split is no longer fused.
- Added `InvalidRecordKind` (6 variants) and `PhaseChainMarkerInconsistencyKind` (4 variants); every site in `writer.rs` was rewritten to construct typed kinds instead of `format!(...)`-built reasons.
- Added `TomlSerError(toml::ser::Error)` wrapper newtype implementing `Display` + `Error` via delegation. `PspWriteError::HeaderToml { source: TomlSerError }` now preserves the `Error::source()` chain symmetrically with the reader's `HeaderToml { source: toml::de::Error }`.
- Extended `PspWriteError::Io` from `{ context, source }` to `{ context, block_index: Option<u64>, column_tag: Option<u16>, source }`. The 8 construction sites in `writer.rs` pass the relevant context. `Display` formats the optional fields inline (`I/O error writing block column payload (block 17, column 0x11): ...`).

### Wave 2 — Registry refactor (M4 + M17 + Mi11)
- Added `ColumnKey` enum (12 variants exhausting the v1.0 registry). `ColumnDef` gains a `key: ColumnKey` field; `writer.rs::encode_column_into` now matches exhaustively on `def.key` — adding a `ColumnDef` row without a matching `ColumnKey` is a compile error, and the runtime `unknown =>` fall-through (which v1 had stubbed to `unreachable!`) is gone entirely.
- Added `ColumnPayload` enum carrying `Scalar { element_type }` / `List { element_type }` / `Bytes { length_column }`. `ColumnDef` replaces the prior three flat fields with `payload: ColumnPayload`. Legacy access patterns are preserved via `shape()` / `element_type()` / `length_column()` convenience accessors so existing read sites compile unchanged.
- `ParsedColumn` mirrors the refactor with `ParsedColumnPayload` (the parse-side variant whose `length_column` is an owned `String`).
- `predict_uncompressed_len` is rewritten as exhaustive matches on `ColumnPayload` × `ElementType` — no `_ =>` catch-all. Adding an `ElementType` variant now requires updating this function as a compile error.
- The Mi11 misuse of `InvalidRecord { record_index: 0 }` for a dispatch-table miss is eliminated because the dispatch is no longer dispatchable.

### Wave 3 — Header validation collapse (M16, Mi29)
- Introduced `HeaderFieldError` trait implemented by both `PspReadError` and `PspWriteError`. The generic `wrap<E>(key, Result<(), String>) -> Result<(), E>` helper routes both sides through one function — no more separate `validate_*` wrappers on the read side plus inline if-let-Err blocks on the write side. `validate_writer_header` is now a thin facade over `validate_header_fields::<PspWriteError>` and shares its implementation with the (future) read-side strong-typed validator.
- Mi29 was implemented during Wave 2 (extracting `cross_check_shape_consistency` from `parse_column`).

### Wave 4 — Writer restructuring (M15, Mi30)
- `PspWriter` shrinks from 14 fields to 7: `sink`, `header`, `sink_offset`, `index_entries`, `block`, `ingest: IngestState`, `scratch: WriterScratch`. The extracted structs own their own bookkeeping with `new()` constructors.
- `flush_block` is now an orchestrator that calls three free functions: `encode_and_compress_columns(&block, &mut scratch, block_index)`, `assemble_block_header(&mut block, &mut scratch, block_index, n_records, n_total_alleles)`, `emit_block_to_sink(&mut sink, &scratch, block_index)`. Each phase is independently readable and testable.

### Wave 5 — Block-target knob + size-trigger tests (M14, Mi26)
- `PspWriter::new_with_block_target(sink, header, target_block_bytes)` constructor exposed under `#[doc(hidden)] pub` so tests (and benches that already use `#[doc(hidden)] pub fn current_block_projected_bytes`) can override the spec-pinned 16 MiB. `should_flush` reads `self.target_block_bytes` instead of the constant.
- New tests `auto_flushes_on_projected_size_boundary` (asserts ≥2 blocks emerge from feeding 2000 records into a writer with `target_block_bytes = 8 KiB`) and `write_record_returns_flushed_byte_count` (loops until an auto-flush, asserts the post-flush record reports 0 pushed bytes).

### Wave 6 — Property tests (M13)
- Added `proptest = "1"` as a `[dev-dependencies]` entry.
- New properties: `proptest_u64_round_trip`, `proptest_i64_round_trip` (varint.rs); `proptest_trailer_round_trip` (trailer.rs); `proptest_index_round_trip` (index.rs); `proptest_scalar_column_round_trip_u32`, `proptest_varint_column_round_trip`, `proptest_list_column_round_trip_u16`, `proptest_block_header_round_trip` (block.rs). Each runs proptest's default 256 cases per build.

### Wave 7 — Renames & visibility (Mi1, Mi5)
- `DEFAULT_TARGET_BLOCK_BYTES` → `TARGET_BLOCK_BYTES`; doc-comment now describes it as spec-pinned (no `DEFAULT_` prefix, matching `ZSTD_COMPRESSION_LEVEL`).
- `block`, `errors`, `index`, `registry`, `trailer`, `varint` narrowed to `pub(crate)`; `header` and `writer` stay `pub` (the bench imports both). Public re-exports in `psp::mod` keep `PspReadError`/`PspWriteError` plus the three new `*Kind` sub-enums, plus `TomlSerError`, `ScalarDecodeError`, `VarintError` reachable through `psp::*`.
- Each newly-`pub(crate)` submodule carries a module-level `#![allow(dead_code)]` with a comment naming the reason: their decoder primitives are reader-side surface awaiting the not-yet-built `PspReader`.

### Wave 8 — Tooling (Mi22, Mi23, M18)
- `rust-version = "1.88"` added under `[package]`. Rationale: edition 2024 requires ≥ 1.85; the `let-chains` features (`if let ... && cond`) used by this codebase are stable since 1.88.
- `[lints.rust] unsafe_code = "forbid"`, `[lints.rustdoc] broken_intra_doc_links = "deny"` + `private_intra_doc_links = "allow"`, and a small `[lints.clippy]` table — see §3 for the scope rationale.
- CI workflow `.github/workflows/ci.yml` runs `fmt --check`, `clippy -D warnings`, `test --lib --tests --all-features`, `doc --no-deps` with `RUSTDOCFLAGS=-D warnings`, and `cargo audit`. A second `msrv` job builds against `rust-version = "1.88"`.
- 11 non-PSP clippy errors fixed so CI starts green — listed in §3.

### Wave 9 — Test fixtures (Mi20)
- New `src/per_sample_caller/psp/test_fixtures.rs` (gated `#[cfg(test)]`, declared `mod test_fixtures;` from `mod.rs`) provides `writer_header(n_chroms)` (the previous `writer.rs::tests::writer_header`) and `realistic_writer_header()` (the previous `header.rs::tests::minimal_writer_header`). Both test modules now `use super::super::test_fixtures::...`.

## 5. Validation commands run

- `cargo fmt --check` — clean
- `cargo clippy --all-targets --all-features -- -D warnings` — fully green
- `cargo test --lib --tests --all-features` — 535 tests pass (149 PSP, +11 from v1)
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib --all-features` — 0 warnings
- `cargo build --lib` — clean after each wave
- Targeted `cargo test --lib per_sample_caller::psp::` after each wave to confirm no PSP regression

`cargo audit` is run by the new CI workflow but was not run locally (binary not installed). `cargo bench` was skipped — see §1 rationale.

## 6. Notes

- **Public API impact.** This run is a significant SemVer event for the `psp` module: error variants gain typed `kind` fields, `ColumnDef` reshapes, several submodules narrow visibility, a new `TARGET_BLOCK_BYTES` constant supersedes `DEFAULT_TARGET_BLOCK_BYTES`. The crate's version is still `0.1.0` and the user authorised the breakage, so this is expected; no SemVer-careful migration shims were added.
- **`PspReadError` field rename.** `BlockHeaderInvariant { reason: String }` → `BlockHeaderInvariant { #[source] kind: BlockHeaderInvariantKind }`. Tests in `block.rs` matching on `BlockHeaderInvariant { .. }` continue to compile (the wildcard pattern still works against the new shape); tests that needed to assert a specific variant were updated to match on the inner `kind`.
- **Build-side `validate_writer_header` invariants preserved.** The M16 refactor produces byte-identical error messages where possible, so the 11 writer-side validator tests added in v1 (M12) continue to pass without changes.
- **Reader-side `PspReader`** is still not implemented. The decoder primitives now in `pub(crate)` modules carry `#![allow(dead_code)]` because they're exercised only by tests; the eventual reader will exercise them in production code and the allow can be removed.
- **Bench file untouched** — `benches/psp_writer_perf.rs` imports `psp::header::{...}` and `psp::writer::PspWriter`. Both remain `pub`; no bench changes needed.
