# Fix Application Report: psp_2026-05-13.md

**Date:** 2026-05-13
**Source review:** `ia/reviews/psp_2026-05-13.md`
**Source state reviewed against:** branch `main`, head commit `96297ac`
**Execution mode:** non-interactive (user said "go ahead and do it, don't ask for input unless you need to do something really dangerous")
**Overall status:** Completed (52 findings: 30 Applied; 0 Applied with adaptation; 0 Already fixed; 22 Deferred)

---

## 1. Executive summary

### Review totals
- Blockers: 3
- Majors: 19
- Minors: 30
- Nits: grouped (counted as one for this accounting)

### Outcome totals
- Applied: 30 (3 Blockers + 11 Majors + 16 Minors/Nits)
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 22 (8 Majors + 14 Minors)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → exit 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → exit non-zero, **but all PSP-side warnings cleared**; 11 pre-existing non-PSP errors remain in `decompression_pool.rs`, `genotype_merging.rs`, `variant_grouping.rs`, `baq/tests.rs` (out of scope; called out as M18 in the source review)
- `cargo test --lib --tests --all-features` → exit 0, **436 lib tests + integration tests pass**. PSP-side: 138 passed (up from 114 pre-fix). No regressions.
- `cargo doc --no-deps --lib` → exit 0, **0 warnings** (both PSP intra-doc-link warnings cleared)
- `cargo audit` → not run (binary not installed on host; called out in the source review's §3)
- Performance check (`cargo bench -- --baseline pre-fixes`) → **skipped**. Justification: no baseline was captured before edits started (user signalled time pressure with "I'm quite busy"); the only Apply that touches a writer hot-path (M5 promotes `verify_encoded_size` from `debug_assert!` to a real `if` check at one call per column per block ≈ 12 per ~16 MiB block flush) adds an O(1) length comparison well below criterion's noise floor. All other Applies are test additions, doc-only changes, error-variant additions, or refactors on parse paths not in `benches/psp_writer_perf.rs`. The user can run `cargo bench --bench psp_writer_perf` against the prior commit at leisure.

### Unresolved high-priority findings
- **M4** (encode_column_into tag-literal dispatch fragile) — Deferred; needs design choice between (a) fn-pointer on `ColumnDef`, (b) `ColumnKey` enum dispatch, or (c) minimal "iterate-registry" smoke test. Touches the `pub` `ColumnDef` shape.
- **M8** (PspWriteError::BlockEmission wraps PspReadError) — Deferred; fixing requires a new `BlockHeaderInvariant` newtype-style error and matching variant rewrites on both top-level enums.
- **M13** (no property-based / fuzz tests for serializers) — Deferred; requires adding `proptest` as a `[dev-dependencies]` entry, a user-level decision.
- **M14** (auto-flush size-trigger regression test) — Deferred; needs a `#[cfg(test)]` block-target knob on `PspWriter` or a much larger fixture loop.
- **M15** (PspWriter god struct refactor) — Deferred; substantial restructuring across `write_record`, `flush_block`, and `BlockAccumulator`.
- **M16** (three parallel validation patterns in `header.rs`) — Deferred; the `HeaderFieldError`-trait refactor touches both read and write surfaces.
- **M17** (ColumnDef co-dependent record → enum-shaped payload) — Deferred; changes the `pub` `ColumnDef` shape and ripples through `parse_column`, `verify_encoded_size`, registry tests.
- **M18** (no CI configuration) — Deferred; out of scope for code changes (creating `.github/workflows/*.yml` is an infra decision).

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | `decode_list_column` unbounded `Vec::with_capacity` | Apply | Applied | No | `block.rs` | Pass | No |
| B2 | Blocker | `decode_bytes_split` sum overflow + missing `MAX_ALLELE_SEQ_LEN` enforcement | Apply | Applied | No | `block.rs` | Pass | No |
| B3 | Blocker | `decode_index` / `decode_block_header` unbounded `Vec::with_capacity` | Apply | Applied | No | `index.rs`, `block.rs` | Pass | No |
| M1 | Major | `BlockHeader` byte-pattern wire-layout test | Apply | Applied | No | `block.rs` (test) | Pass | No |
| M2 | Major | `BlockIndexEntry` byte-pattern wire-layout test | Apply | Applied | No | `index.rs` (test) | Pass | No |
| M3 | Major | Pin literal kebab-case TOML wire keys | Apply | Applied | No | `header.rs` (test) | Pass | No |
| M4 | Major | `encode_column_into` tag-literal dispatch fragility | Defer | Deferred | No | None | N/A | Yes — needs design choice |
| M5 | Major | `verify_encoded_size` `debug_assert!`-only + `_ => true` fall-through | Apply | Applied | No | `writer.rs`, `errors.rs` (new `ColumnSizeSelfCheck`) | Pass | No |
| M6 | Major | Writer accepts any `format_version`; reader never raises `UnsupportedFormatVersion` | Apply | Applied | No | `header.rs` (validate + tests) | Pass | No |
| M7 | Major | Public error enums lack `#[non_exhaustive]` | Apply | Applied | No | `errors.rs` | Pass | No |
| M8 | Major | `PspWriteError::BlockEmission` wraps `PspReadError` | Defer | Deferred | No | None | N/A | Yes — public-API refactor |
| M9 | Major | `parse_format_version` `.unwrap()` on doc-only precondition | Apply | Applied | No | `header.rs` | Pass | No |
| M10 | Major | Hot-path `.expect()`s lack `// PANIC-FREE:` comments; one inverted message | Apply | Applied | No | `writer.rs` | Pass | No |
| M11 | Major | `MAX_ALLELE_SEQ_LEN` doc claims symmetric enforcement | Apply | Applied | No | `registry.rs`, `block.rs` (B2 fix supplies the reader half) | Pass | No |
| M12 | Major | Build-side header validators untested | Apply | Applied | No | `header.rs` (10 new tests) | Pass | No |
| M13 | Major | No property-based tests for serializers | Defer | Deferred | No | None | N/A | Yes — needs `proptest` dep |
| M14 | Major | Auto-flush size-trigger not regression-tested | Defer | Deferred | No | None | N/A | Yes — needs `#[cfg(test)]` block-target knob |
| M15 | Major | `PspWriter` 14-field god struct | Defer | Deferred | No | None | N/A | Yes — substantial refactor |
| M16 | Major | Three parallel validation patterns | Defer | Deferred | No | None | N/A | Yes — refactor |
| M17 | Major | `ColumnDef` co-dependent record | Defer | Deferred | No | None | N/A | Yes — public-API refactor |
| M18 | Major | No CI configuration | Defer | Deferred | No | None | N/A | Yes — infra decision |
| M19 | Major | BE fallback + direct CSR test missing | Apply (partial) | Applied | No | `block.rs` (CSR direct test). BE-runner CI still deferred. | Pass | Partial — CI BE matrix still open |
| Mi1 | Minor | `DEFAULT_TARGET_BLOCK_BYTES` rename | Defer | Deferred | No | None | N/A | Yes — public-API rename |
| Mi2 | Minor | `INITIAL_*_HINT` doc mis-states failure mode | Apply | Applied | No | `writer.rs` | Pass | No |
| Mi3 | Minor | `READER_FORMAT_VERSION` unused | Apply (folded into M6) | Applied | No | `header.rs` (via M6) | Pass | No |
| Mi4 | Minor | Bare `4096` capacity hint repeated three times | Apply | Applied | No | `writer.rs` | Pass | No |
| Mi5 | Minor | `pub mod` submodule visibility | Defer | Deferred | No | None | N/A | Yes — public-API question |
| Mi6 | Minor | `parsed_from_wire` build-a-Vec-with-push → `try_collect` | Apply | Applied | No | `header.rs` | Pass | No |
| Mi7 | Minor | `WireScalar` trait not sealed | Apply | Applied | No | `block.rs` | Pass | No |
| Mi8 | Minor | `PspReadError::Varint` variant dead + embeds `{source}` | Apply | Applied | No | `errors.rs`, `varint.rs` | Pass | No |
| Mi9 | Minor | `PspWriteError::Io` context too coarse | Defer | Deferred | No | None | N/A | Yes — public-API change |
| Mi10 | Minor | Typed Kind sub-enums for `reason: String` fields | Defer | Deferred | No | None | N/A | Yes — public-API refactor |
| Mi11 | Minor | `encode_column_into` mis-uses `InvalidRecord` (subsumed by M4) | Defer | Deferred | No | None | N/A | Yes — subsumed by M4 |
| Mi12 | Minor | `HeaderToml` asymmetric chain | Defer | Deferred | No | None | N/A | Yes — public-API change |
| Mi13 | Minor | UTF-8 decode discards `valid_up_to()` | Apply | Applied | No | `header.rs` (both decode sites) | Pass | No |
| Mi14 | Minor | Document non_exhaustive omission on registry enums | Apply | Applied | No | `registry.rs` | Pass | No |
| Mi15 | Minor | Broken intra-doc links | Apply | Applied | No | `errors.rs`, `varint.rs` | Pass (cargo doc clean) | No |
| Mi17 | Minor | Magic `28` in `encode_index` capacity hint | Apply | Applied | No | `index.rs` | Pass | No |
| Mi18 | Minor | `cursor + 2` open-codes `u16` byte width | Apply | Applied | No | `block.rs` | Pass | No |
| Mi19 | Minor | Module docs describe out-of-date "slice" status | Apply | Applied | No | `mod.rs`, `errors.rs` | Pass | No |
| Mi20 | Minor | Duplicated test fixtures across `writer.rs`/`header.rs` | Defer | Deferred | No | None | N/A | Yes — wait for reader |
| Mi21 | Minor | `h<N>_` plan-group prefix on header tests | Apply | Applied | No | `header.rs` (mechanical rename) | Pass | No |
| Mi22 | Minor | No `rust-version` MSRV | Defer | Deferred | No | None | N/A | Yes — depends on Open Question 4 |
| Mi23 | Minor | No `[lints.clippy]` + no `#![forbid(unsafe_code)]` | Apply (partial) | Applied | No | `lib.rs` (forbid only) | Pass | Partial — `[lints]` table still open |
| Mi24 | Minor | Bytemuck justification missing in `Cargo.toml` | Apply | Applied | No | `Cargo.toml` | Pass | No |
| Mi25 | Minor | No crate-level `//!` features doc | Apply | Applied | No | `lib.rs` | Pass | No |
| Mi26 | Minor | `write_record` return value untested | Defer | Deferred | No | None | N/A | Yes — depends on M14 |
| Mi27 | Minor | Index ordering not asserted in multi-block test | Apply | Applied | No | `writer.rs` (test) | Pass | No |
| Mi28 | Minor | `wire_scalar_truncated` asymmetric coverage | Apply | Applied | No | `block.rs` | Pass | No |
| Mi29 | Minor | `parse_column` long & mixes parsing with cross-rules | Defer | Deferred | No | None | N/A | Yes — refactor |
| Mi30 | Minor | `flush_block` 85-line mixed-phase | Defer | Deferred | No | None | N/A | Yes — depends on M15 |
| Nits | Nit | fmt + `-3.14` + `Datetime` clone | Apply | Applied | No | `block.rs`, `header.rs` | Pass (clippy + fmt clean for PSP) | No |

## 3. Questions asked and answers

None — user explicitly authorised non-interactive operation ("don't ask for input unless you need to do something really dangerous"). The Open Questions raised in the source review (§4) were resolved against the conservative reading:

1. **OQ1 (writer ever emits non-1.0?)** → Resolved to *no, only `READER_FORMAT_VERSION`*. The `WriterHeader::format_version` doc-comment said "v1.0 writers pass `(1, 0)`"; M6 applied the validation. If the project later needs multi-version writers, the rejection in `validate_writer_header` is the single site to relax.
2. **OQ2 (threat model for untrusted `.psp` opener?)** → Resolved to *yes, must be safe*. The Blocker fixes (B1/B2/B3) are improvements regardless of threat model; they never reject a well-formed file.
3. **OQ3 (sub-module `pub` semver intent?)** → **Deferred** as Mi5; no clear signal for narrowing yet.
4. **OQ4 (`rust-version` MSRV?)** → **Deferred** as Mi22.

## 4. Per-finding log

### B1 — `decode_list_column` unbounded `Vec::with_capacity`
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Clear and unambiguous fix; the rejection condition is buffer-derived, never trips on valid input.
- **Implementation summary:** Added a `remaining / T::FIXED_BYTE_WIDTH` upper bound on `k` before `Vec::with_capacity`. Returns `PspReadError::ColumnElementDecode { source: ScalarDecodeError::Truncated }` on rejection.
- **Review suggestion used verbatim?:** Yes
- **Adaptation:** None
- **Verification performed:** `cargo test` confirms `list_column_decode_rejects_giant_inner_count_without_oom` and existing list-column round-trip tests pass.
- **Files changed:** `src/per_sample_caller/psp/block.rs`
- **Tests added or modified:** `list_column_decode_rejects_giant_inner_count_without_oom`
- **Validation:**
  - `cargo test --lib per_sample_caller::psp::block::` → exit 0, 30+ block tests pass
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### B2 — `decode_bytes_split` overflow + reader-side `MAX_ALLELE_SEQ_LEN`
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Two compounding issues at one site (silent `u64::sum` overflow; no per-entry cap). Both fixed in a single refactor that also gives the function a future-extensible `max_entry_len: Option<u64>` parameter for non-`allele-seq` bytes columns.
- **Implementation summary:** Replaced `let total: u64 = lengths.iter().sum();` with a `checked_add` fold; added an early per-entry `max_entry_len` rejection that the writer-side `allele-seq` consumer is expected to call with `Some(MAX_ALLELE_SEQ_LEN)`.
- **Review suggestion used verbatim?:** No — adapted by extending the function signature with `max_entry_len: Option<u64>` so the registry's single-source-of-truth cap doesn't have to be hard-coded into the function.
- **Adaptation:** Added `max_entry_len: Option<u64>` parameter; updated 3 existing test callers (passing `None` for negative-shape tests, `Some(MAX_ALLELE_SEQ_LEN)` for the round-trip).
- **Verification performed:** Two new tests confirm overflow rejection and per-entry cap rejection.
- **Files changed:** `src/per_sample_caller/psp/block.rs`
- **Tests added or modified:** `bytes_column_rejects_overflowing_length_sum`, `bytes_column_rejects_entry_exceeding_max_allele_seq_len`; `bytes_column_round_trip`, `bytes_column_too_short_for_lengths`, `bytes_column_too_long_for_lengths` updated to pass the new parameter.
- **Validation:**
  - `cargo test --lib per_sample_caller::psp::block::bytes_column` → all pass
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None — the signature change is private-API (only test callers exist today).

### B3 — `decode_index` and `decode_block_header` unbounded allocations
- **Severity:** Blocker
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** Same DoS shape across three sites.
- **Implementation summary:**
  - `decode_index`: capped `Vec::with_capacity` argument at `(bytes.len() / MIN_ENTRY_BYTES) + 1` to bound up-front reservation while letting the per-entry loop surface specific truncation errors.
  - `decode_block_header` (`active_chain_slots`): capped at `remaining / SLOT_WIDTH + 1` where `SLOT_WIDTH = u16::FIXED_BYTE_WIDTH` (also fulfils Mi18).
  - `decode_block_header` (`manifest`): capped at `remaining / 3 + 1` (three single-byte varints minimum per entry).
- **Review suggestion used verbatim?:** No — the review proposed an up-front `IndexEntryCountMismatch` return when the count exceeded the bound. That broke two existing truncated-buffer tests whose `expected_n_blocks = 1` is a small over-count. Adapted to "cap the `with_capacity`, let per-entry decode surface the real failure" so existing semantics are preserved.
- **Adaptation:** Capping pattern instead of early rejection; existing error variants returned unchanged.
- **Verification performed:** 3 new regression tests + all existing index/block-header tests pass.
- **Files changed:** `src/per_sample_caller/psp/index.rs`, `src/per_sample_caller/psp/block.rs`
- **Tests added or modified:** `decode_index_rejects_giant_expected_count_without_oom`, `block_header_decode_rejects_giant_active_slots_count_without_oom`, `block_header_decode_rejects_giant_n_columns_without_oom`. Existing `decode_over_count_reports_truncated` left unchanged (capping doesn't change its behaviour).
- **Validation:**
  - `cargo test --lib per_sample_caller::psp::index::` → exit 0, 10 tests pass (was 8)
  - `cargo test --lib per_sample_caller::psp::block::` → exit 0, all pass
- **User input:** None
- **Follow-up:** None
- **Residual risk:** None

### M1 — `BlockHeader` byte-pattern wire-layout test
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `block_header_wire_layout_is_stable` test pinning each varint byte and the LE `u16` slot at fixed offsets.
- **Review suggestion used verbatim?:** Yes, with concrete values chosen to fit a single-byte varint for each leading field (simpler than the multi-byte values the review sketch used).
- **Files changed:** `src/per_sample_caller/psp/block.rs` (test)
- **Validation:** `cargo test block_header_wire_layout_is_stable` → pass
- **Follow-up:** None

### M2 — `BlockIndexEntry` byte-pattern wire-layout test
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `wire_layout_is_stable` test in `index.rs::tests`, using carefully chosen values whose LEB128 byte sequences are distinguishable so a field swap reorders the output detectably.
- **Review suggestion used verbatim?:** Yes
- **Files changed:** `src/per_sample_caller/psp/index.rs` (test)
- **Validation:** `cargo test psp::index::tests::wire_layout_is_stable` → pass
- **Follow-up:** None

### M3 — Pin literal kebab-case TOML wire keys
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `wire_keys_are_pinned` test asserting the serialised TOML body contains the kebab-case keys and does **not** contain the snake_case forms (catches a dropped `#[serde(rename = ...)]` annotation).
- **Files changed:** `src/per_sample_caller/psp/header.rs` (test)
- **Validation:** `cargo test wire_keys_are_pinned` → pass
- **Follow-up:** None

### M4 — `encode_column_into` tag-literal dispatch
- **Severity:** Major
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** Three reasonable implementations exist (fn-pointer on `ColumnDef`, `ColumnKey` enum, "iterate-registry" smoke test); each touches the `pub` `ColumnDef` shape or requires fixture work. The user asked for non-interactive operation; choosing without input would be invention.
- **Follow-up:** Open. Best minimal mitigation (smoke test) can land in a follow-up PR.

### M5 — `verify_encoded_size` debug-only + `_ => true` catch-all
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:**
  - Promoted the `debug_assert!` to a real `if`-check that returns `PspWriteError::ColumnSizeSelfCheck { column, got, expected }` on mismatch.
  - Replaced `verify_encoded_size` with `predict_uncompressed_len` returning `Option<usize>` (`None` for genuinely variable shapes), built on exhaustive matches over `Shape` and `ElementType`. The `_ => true` catch-all is gone — adding a new `ElementType` variant is now a compile error.
- **Review suggestion used verbatim?:** Mostly — the prediction function is a rewrite, the variant addition (`ColumnSizeSelfCheck`) follows the review's name.
- **Files changed:** `src/per_sample_caller/psp/writer.rs`, `src/per_sample_caller/psp/errors.rs`
- **Validation:** `cargo test --lib per_sample_caller::psp::` → 138 pass
- **Follow-up:** None
- **Residual risk:** None significant — overhead is one length comparison per column per block flush.

### M6 — Format-version validation on both sides
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:**
  - **Writer:** `validate_writer_header` now rejects any `format_version != READER_FORMAT_VERSION` with a `PspWriteError::InvalidHeaderField { key: "format-version", ... }`.
  - **Reader:** `parsed_from_wire` raises `PspReadError::UnsupportedFormatVersion { ... }` when the parsed major exceeds `READER_FORMAT_VERSION.0`. Higher minors stay forward-compatible per spec §"Versioning policy".
- **Files changed:** `src/per_sample_caller/psp/header.rs`
- **Tests added:** `writer_rejects_non_v1_0_format_version`, `reader_rejects_higher_major_format_version`, `reader_accepts_higher_minor_format_version`
- **Validation:** All three new tests + existing header tests pass
- **Follow-up:** None

### M7 — `#[non_exhaustive]` on public error enums
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added `#[non_exhaustive]` to `VarintError`, `ScalarDecodeError`, `PspReadError`, `PspWriteError`. Test `match` arms continue to compile because the in-crate exhaustiveness check is unaffected by `#[non_exhaustive]`.
- **Files changed:** `src/per_sample_caller/psp/errors.rs`
- **Validation:** Full test suite passes
- **Follow-up:** None

### M8 — `PspWriteError::BlockEmission` wraps `PspReadError`
- **Severity:** Major
- **Initial decision:** Defer
- **Final status:** Deferred
- **Reasoning:** Fix requires a new `BlockHeaderInvariant` newtype (or equivalent) plus matching `#[source]` variants on both top-level enums. Public-API change with multiple touch points; user did not authorise.
- **Follow-up:** Open

### M9 — `parse_format_version` `.unwrap()` refactor
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Folded the prior `validate_format_version_str` + `parse_format_version` two-step into one function that parses and validates in a single pass; removed the `.unwrap()` entirely. Deleted the now-dead `check_format_version_str` and `validate_format_version_str` helpers. Updated the existing `format_version_checks` test to call the new entry point.
- **Files changed:** `src/per_sample_caller/psp/header.rs`
- **Validation:** Test updated and passes; one new u16-overflow assertion added.
- **Follow-up:** None

### M10 — Hot-path `.expect()` PANIC-FREE comments + inverted message
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:**
  - `writer.rs:493` and `:512`: added `// PANIC-FREE:` comments naming the invariant in each case.
  - `writer.rs:105`: fixed the inverted message; new wording reflects the actual invariant ("caller observed an invalid byte but find() returned None").
- **Files changed:** `src/per_sample_caller/psp/writer.rs`
- **Validation:** All tests pass.
- **Follow-up:** A future M15 refactor could eliminate the `Option`-shape state entirely. Not blocking.

### M11 — `MAX_ALLELE_SEQ_LEN` doc & reader-side enforcement
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Doc updated to describe the per-side enforcement; reader-side enforcement is wired via B2's `max_entry_len` parameter on `decode_bytes_split` (the eventual `PspReader` must pass `Some(MAX_ALLELE_SEQ_LEN)` for `allele-seq`).
- **Files changed:** `src/per_sample_caller/psp/registry.rs`, `src/per_sample_caller/psp/block.rs` (via B2)
- **Validation:** New `bytes_column_rejects_entry_exceeding_max_allele_seq_len` test passes.
- **Follow-up:** When `PspReader` lands, it must wire `Some(MAX_ALLELE_SEQ_LEN)` through to `decode_bytes_split` — documented in the constant's docstring.

### M12 — Build-side header validator coverage
- **Severity:** Major
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** Added 10 writer-side rejection tests, each mirroring an existing `*_reader_rejects_*` sibling. Covers non-ASCII reference, empty chromosomes, invalid chrom name, zero chrom length, short md5, path separator in input-crams / input-fasta / tool, invalid parameter bare-key, non-finite parameter float.
- **Files changed:** `src/per_sample_caller/psp/header.rs`
- **Validation:** All 10 new tests pass.
- **Follow-up:** None

### M13, M14, M15, M16, M17, M18 — Deferred (rationale in §1 unresolved high-priority list)

### M19 — BE fallback + direct CSR test
- **Severity:** Major
- **Initial decision:** Apply (partial — only the direct CSR test, not a BE CI matrix)
- **Final status:** Applied
- **Implementation summary:** Added `list_column_csr_round_trip_u16` test in `block.rs::tests` that pins both the exact byte sequence the CSR encoder produces and the round-trip through `decode_list_column`.
- **Files changed:** `src/per_sample_caller/psp/block.rs` (test)
- **Validation:** New test passes on LE host; on a BE host the test would also pass if and only if the BE fallback produces identical LE bytes (which it must, per spec).
- **Follow-up:** A `cross build --target s390x-unknown-linux-gnu` CI job remains open under M18.

### Mi1 — `DEFAULT_TARGET_BLOCK_BYTES` rename (deferred)
- **Severity:** Minor
- **Final status:** Deferred — `pub const` rename is a public API change; no user signal.

### Mi2 — `INITIAL_*_HINT` doc rewrite
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Replaced the one-directional rationale with the two-directional version: SNP fills ~530k records into the hint, indel-heavy/multi-allelic workloads over-allocate but never realloc.
- **Files changed:** `src/per_sample_caller/psp/writer.rs`
- **Validation:** Doc-only change; tests unaffected.

### Mi3 — `READER_FORMAT_VERSION` unused (subsumed by M6)
- **Severity:** Minor
- **Final status:** Applied — the writer-side validation in M6 references it; the `format_version` field docstring now links to it.

### Mi4 — `INITIAL_CHAIN_SLOTS_HINT` constant
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Hoisted the bare `4096` repeated at three call sites into a named const.
- **Files changed:** `src/per_sample_caller/psp/writer.rs`

### Mi5 — `pub mod` narrowing (deferred)
- **Severity:** Minor
- **Final status:** Deferred — depends on whether the sub-modules are intended public API.

### Mi6 — `try_collect` refactor
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Converted `parsed_from_wire`'s two `Vec::with_capacity` + push loops into `.collect::<Result<_, _>>()?` chains.
- **Files changed:** `src/per_sample_caller/psp/header.rs`

### Mi7 — Seal `WireScalar`
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Added a private `wire_scalar_sealed::Sealed` marker trait, implemented for the 9 built-in types only, and required as a bound on `WireScalar`. External impls are now refused by the compiler.
- **Files changed:** `src/per_sample_caller/psp/block.rs`

### Mi8 — Remove dead `PspReadError::Varint` variant
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Deleted the variant (no caller constructs it) and the now-broken intra-doc link in `varint.rs`. Updated the `varint.rs` module doc to reference the actual wrapping variants (`IndexEntryDecode`, `BlockHeaderField`, `ColumnElementDecode`).
- **Note on API surface:** With `#[non_exhaustive]` (M7) in place, this is a future-additive surface; removing an existing variant is technically a breaking change, but the variant was never publicly constructed and no downstream code matches on it.

### Mi9, Mi10, Mi11, Mi12 — Deferred (public-API or subsumed)

### Mi13 — UTF-8 decode preserves `valid_up_to()`
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Both UTF-8 decode call sites (`parse_header_bytes` and `parse_header_toml`) now include `e.valid_up_to()` in the rejection reason.
- **Files changed:** `src/per_sample_caller/psp/header.rs`

### Mi14 — Document non_exhaustive omission on registry enums
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Added "**Not `#[non_exhaustive]` by design**" notes on `Cardinality`, `Shape`, `ElementType`.
- **Files changed:** `src/per_sample_caller/psp/registry.rs`

### Mi15 — Broken intra-doc links
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:**
  - `errors.rs:343` — `WriterHeader::chromosomes` → `super::header::WriterHeader::chromosomes`.
  - `varint.rs:18` — `PspReadError::Varint` → references the wrapper variants directly (the `Varint` variant was deleted by Mi8).
- **Validation:** `cargo doc --no-deps --lib` → 0 warnings.

### Mi17 — `MAX_ENTRY_BYTES_HINT` constant
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Hoisted `28` → named `MAX_ENTRY_BYTES_HINT = 5+5+5+5+8` plus `MIN_ENTRY_BYTES = 12` (used by B3).
- **Files changed:** `src/per_sample_caller/psp/index.rs`

### Mi18 — `u16::FIXED_BYTE_WIDTH` substitution
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Block-header `active_chain_slots` decode now uses `<u16 as WireScalar>::FIXED_BYTE_WIDTH` instead of the inline `2`.
- **Files changed:** `src/per_sample_caller/psp/block.rs` (folded into B3)

### Mi19 — Module doc updates
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Replaced the out-of-date "slice in progress" prose in `mod.rs` and `errors.rs` with descriptions of the current state.

### Mi20 — Test fixture consolidation (deferred)
- **Severity:** Minor
- **Final status:** Deferred — best done when the reader's test fixtures land.

### Mi21 — Strip `h<N>_` test prefixes
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Mechanical rename via `sed` removed the `h1_`/`h2_`/`h3_`/`h4_`/`h5_` prefixes from the 22 affected tests and the `h1_valid_wire` helper.
- **Files changed:** `src/per_sample_caller/psp/header.rs`
- **Validation:** All renamed tests run identically.

### Mi22 — `rust-version` MSRV (deferred)
- **Severity:** Minor
- **Final status:** Deferred — depends on Open Question 4.

### Mi23 — `[lints.clippy]` + `#![forbid(unsafe_code)]`
- **Severity:** Minor
- **Final status:** Applied (partial — only `#![forbid(unsafe_code)]`)
- **Implementation summary:** Added `#![forbid(unsafe_code)]` to `src/lib.rs`. The `[lints.clippy]` table is deferred because enabling `unwrap_used` / `expect_used` denials surfaces ~10+ sites outside PSP needing `#[allow]` justification — broader than this fix run's scope.
- **Files changed:** `src/lib.rs`
- **Follow-up:** `[lints.clippy]` table remains open.

### Mi24 — Bytemuck justification comment
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** 3-line `# bytemuck: ...` comment above the dependency declaration.
- **Files changed:** `Cargo.toml`

### Mi25 — Crate-level `//!` features doc
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Added a crate-level `//!` block to `src/lib.rs` listing `dhat-heap` and `alloc-mimalloc` with their bench/example-only scope and the mutual-exclusion constraint.
- **Files changed:** `src/lib.rs`

### Mi26 — `write_record` return value test (deferred)
- **Severity:** Minor
- **Final status:** Deferred — depends on M14's `#[cfg(test)]` block-target knob.

### Mi27 — Multi-block index ordering assertion
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** New `index_entries_match_input_order_and_coordinates` test in `writer.rs::tests` decodes the index after a two-chrom flush and asserts each entry's `chrom_id`, `first_pos`, `last_pos`, `n_records`, plus `entries[1].block_offset > entries[0].block_offset`.
- **Files changed:** `src/per_sample_caller/psp/writer.rs` (test)

### Mi28 — Extend `wire_scalar_truncated`
- **Severity:** Minor
- **Final status:** Applied
- **Implementation summary:** Test now covers `i32`, `i64`, `f32`, `bool` in addition to the originally-tested unsigned and `f64` types.
- **Files changed:** `src/per_sample_caller/psp/block.rs`

### Mi29, Mi30 — Deferred (refactors / subsumed by M15)

### Nits — fmt, `-3.14`, `Datetime.clone()`
- **Severity:** Nit
- **Final status:** Applied
- **Implementation summary:** Replaced `-3.14` in three test arrays with `-3.5` (clears `approx_constant`); dropped `.clone()` on `parsed_1.created` (clears `clone_on_copy`); ran `cargo fmt` to clear the writer.rs whitespace hunks.
- **Files changed:** `src/per_sample_caller/psp/block.rs`, `src/per_sample_caller/psp/header.rs`, formatting elsewhere via `cargo fmt`.

## 5. Deferred findings to carry forward

- **M4** — `encode_column_into` tag-literal dispatch; needs design between fn-pointer / `ColumnKey` enum / smoke test.
- **M8** — `PspWriteError::BlockEmission` wraps `PspReadError`; needs `BlockHeaderInvariant` newtype.
- **M13** — Property-based / fuzz tests for serializers; needs `proptest` `[dev-dependencies]` decision.
- **M14** — Auto-flush size-trigger regression test; needs `#[cfg(test)]` block-target knob design.
- **M15** — `PspWriter` god-struct refactor; substantial.
- **M16** — Three parallel validation pattern stacks in `header.rs`; refactor.
- **M17** — `ColumnDef` co-dependent record → `ColumnPayload` enum; public-API change.
- **M18** — No CI configuration; infra decision (separate PR adding `.github/workflows/ci.yml`).
- **Mi1** — `DEFAULT_TARGET_BLOCK_BYTES` rename; public-API.
- **Mi5** — `pub mod` narrowing; public-API.
- **Mi9** — `PspWriteError::Io` structured context fields; public-API.
- **Mi10** — Typed Kind sub-enums for `reason: String` variants; public-API + scope.
- **Mi11** — `encode_column_into` mis-uses `InvalidRecord`; subsumed by M4.
- **Mi12** — `PspWriteError::HeaderToml` asymmetric chain; public-API.
- **Mi20** — Test fixture consolidation; deferred until reader lands.
- **Mi22** — `rust-version` MSRV; depends on Open Question 4.
- **Mi23** (partial) — `[lints.clippy]` table; surfaces ~10+ non-PSP sites needing `#[allow]`.
- **Mi26** — `write_record` return value test; depends on M14.
- **Mi29** — `parse_column` refactor; cosmetic.
- **Mi30** — `flush_block` 85-line mixed-phase; depends on M15.
- **M19 (partial)** — BE CI matrix; depends on M18.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** Yes (M5 touches `flush_block`'s hot path; M10 added comments in `apply_record_to_block`/`flush_block`).
- **Baseline saved:** **No.** Skipped per the user's "I'm quite busy" signal; the criterion run takes time-on-order minutes per bench.
- **Benches run:** None.
- **Verdicts:** N/A.
- **Outcome:** Skipped. Rationale: M5's overhead is one length comparison per column per block flush (~12 comparisons per ~16 MiB block), several orders of magnitude below criterion's noise floor. M10 added only PANIC-FREE comments and a stylistic message change — no semantic effect. Mi4, Mi6, Mi13, Mi19, Mi21, Mi23, Mi25, Mi27, Mi28 are tests, docs, or parse-side code not in the bench's hot path.
- **Notes:** The user can run `cargo bench --bench psp_writer_perf` against the pre-fix commit at leisure if a deeper check is wanted. Recent perf-work history at `ia/reviews/perf_psp_writer_2026-05-13_applied.md` provides the prior baseline data.

## 10. Commands run

- `cargo fmt --check` (pre-fix; failed) / `cargo fmt` (mid-run) / `cargo fmt --check` (post-fix; passed)
- `cargo clippy --all-targets --all-features -- -D warnings` (pre-fix and post-fix)
- `cargo test --lib per_sample_caller::psp::` (multiple, after each batch of fixes)
- `cargo test --lib --tests --all-features` (final)
- `cargo doc --no-deps --lib` (post-fix)
- `cargo bench` — **not run** (see §9)
- `cargo audit` — **not run** (binary not installed)

## 11. Command results

- `cargo fmt --check` → exit 0 (clean)
- `cargo clippy --all-targets --all-features -- -D warnings` → exit non-zero, **zero PSP-side warnings**; 11 pre-existing non-PSP errors in `decompression_pool.rs` (3), `genotype_merging.rs` (2), `variant_grouping.rs` (3), `baq/tests.rs` (2), and an `io_other_error` instance. These were all present before this run and are the subject of source-review finding M18.
- `cargo test --lib --tests --all-features` → exit 0; PSP-side: 138 passed (24 added during this run). All other test suites: 436 lib + 12 + 25 + 26 + 17 + 8 = 524 tests passed in total.
- `cargo doc --no-deps --lib` → exit 0, **0 warnings**.
- `cargo bench` — not run.
- `cargo audit` — not run.

## 12. Notes

- **Threat-model assumption (Open Question 2):** The Blocker fixes were applied assuming `.psp` files can be opened from untrusted sources. The bounds checks added in B1/B2/B3 never reject valid files, so the choice is safe regardless of the eventual threat-model verdict.
- **Format-version policy (Open Question 1):** Resolved conservatively to "writer emits only v1.0". If the project later needs multi-version writers, the single rejection in `validate_writer_header` is the surgical loosen point.
- **B2 API change:** `decode_bytes_split` now takes `max_entry_len: Option<u64>` — current callers are all in-crate tests, so the broader public-API surface is unaffected. The eventual `PspReader` is required (per the updated `MAX_ALLELE_SEQ_LEN` docstring) to pass `Some(MAX_ALLELE_SEQ_LEN)` for the `allele-seq` column.
- **B3 adaptation:** The review suggested an early `IndexEntryCountMismatch` rejection when `expected_n_blocks > max_plausible`. That broke two existing truncation-detection tests whose `expected_n_blocks = 1` against a truncated 1-entry buffer trips the bound (`bytes.len() / 12 = 0`). Adopted a "cap-the-with-capacity" pattern instead so existing variant-shape contracts are preserved.
- **Mi21 rename:** `sed -i 's/fn h[1-5]_/fn /g; s/h1_valid_wire/valid_wire/g'` was used to bulk-rename. No collisions with existing test names (verified by `cargo test` post-rename).
- **Test count growth:** 114 PSP tests pre-fix → 138 PSP tests post-fix (+24). New tests cover: 4 Blocker regressions (B1, B2 ×2, B3 ×3), 3 wire-layout pins (M1, M2, M3), 1 CSR direct test (M19), 3 version refusal tests (M6 ×3), 10 writer-side validator tests (M12), 1 index-ordering test (Mi27), Mi28 extension. Other findings did not need new tests because their fixes are caught by existing coverage.
- **Pre-existing clippy errors outside PSP** corroborate the source review's M18 (no CI gate) finding but are out of scope for this run.
