# Fix Application Report: psp_reader_2026-05-13.md

**Date:** 2026-05-13
**Source review:** `ia/reviews/psp_reader_2026-05-13.md`
**Source state reviewed against:** commit `2d7f211` on `main`
**Execution mode:** non-interactive (autonomous per the project rule
"go ahead, don't ask unless something is really dangerous")
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 3 (B1, B2, B3)
- Majors: 16 (M1-M16)
- Minors: 19 (Mi1-Mi19)
- Nits: ~10

### Outcome totals
- Applied: 28 (3 B + 12 M + 13 Mi)
- Applied with adaptation: 1 (M1 — restructured to lift the unwrap into the caller via `match`, not `// PANIC-FREE:` comment alone)
- Already fixed: 0
- Deferred: 9 (M14 fuzz infra, M15 7-step extraction refactor, Mi10 decode_one_column extraction, Mi11 round_trip helper, Mi12 open_reader helper, Mi14 record/allele helpers, Mi18 u32→usize counters, plus most Nits)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → exit 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → exit 0, clean
- `cargo test --lib --tests --all-features` → exit 0, **562 passed / 0 failed** (474 lib + 88 integration)
- `cargo test --doc --all-features per_sample_caller::psp::reader` → exit 0, 2 doctests pass (new `PspReader::new` + `region_records` examples)
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib` → exit 0, clean
- `cargo audit` → not run (project policy)
- **Performance check** → not triggered (no `benches/psp_reader*.rs` exists; recommendation to add one is captured in M13's deferred-follow-up note)

### Unresolved high-priority findings

None. All 3 Blockers and 14 of 16 Majors landed. M14 (fuzz infra)
and M15 (open-sequence extraction refactor) are deferred to a
follow-up review cycle — neither is a correctness fix.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | Per-block manifest required-tag coverage missing | Apply | Applied | reader.rs, errors.rs | Pass; new test `manifest_missing_required_column_for_every_required_tag` |
| B2 | Blocker | Untrusted `compressed_len` drives 4 GiB allocation | Apply | Applied | reader.rs | Pass; new test `decode_one_column_rejects_oversized_compressed_len` |
| B3 | Blocker | `sum(n_alleles) == n_total_alleles` unenforced | Apply | Applied | reader.rs, errors.rs | Pass; new test `n_alleles_sum_mismatch_against_n_total_alleles` |
| M1 | Major | `materialise_next_record` `.expect` without PANIC-FREE | Apply | Applied with adaptation | reader.rs | Pass; lifted to `match self.cur_block.as_mut() { Some(b) => b, None => return Err(Io { … }) }` instead of an `expect` + comment |
| M2 | Major | `try_into().unwrap()` on `[u8; 12][4..12]` slice | Apply | Applied | reader.rs | Pass |
| M3 | Major | `VarintError::Overflow` wrong source for "header exceeds cap" | Apply | Applied | reader.rs, errors.rs | Pass; new `PspReadError::BlockHeaderExceedsCap` |
| M4 | Major | `field: "header truncated mid-decode"` is a description | Apply | Applied | reader.rs, errors.rs | Pass; new `PspReadError::BlockHeaderTruncated` |
| M5 | Major | `BlockIndexOffsetInvalid` with `min: 0` factually wrong | Apply | Applied | reader.rs | Pass; range + monotonic checks now share one half-open interval after the header parse |
| M6 | Major | `NonFiniteFloat` provably unreached today; F5 deferred | Apply | Applied | reader.rs | Pass; column name now sourced from registry. F5 fixture-patch test still deferred (requires `replace_column_payload_f64` helper) |
| M7 | Major | Spec check #10 first-block-of-chrom snapshot empty | Apply | Applied | reader.rs, errors.rs | Pass; new `BlockHeaderInvariantKind::NonEmptySnapshotAtChromStart` + test `region_mode_rejects_non_empty_snapshot_at_chrom_start` |
| M8 | Major | Reader uses `BTreeSet<SlotId>` instead of sorted `Vec` | Apply | Applied | reader.rs | Pass; mirrors writer's S1+S2 pattern (sorted Vec + binary_search) |
| M9 | Major | Three `#[allow(dead_code)]` fields are genuinely dead | Apply | Applied | reader.rs | Pass; `header_end_offset`, `n_total_alleles`, `active_chain_slots_snapshot` removed |
| M10 | Major | `DecodedColumn::Unknown` swallows new keys silently | Apply | Applied | reader.rs | Pass; `decode_one_column` returns `Option<DecodedColumn>`, no `Unknown` variant |
| M11 | Major | F5 finite-float sweep hard-coded to AlleleQSumLog | Apply | Applied | reader.rs, registry.rs | Pass; new `ColumnDef::finite_constraint: bool`, set true on `allele-q-sum-log` |
| M12 | Major | `BlockHeader::active_chain_slots: Vec<u16>` not tied to SlotId | Apply | Applied | block.rs | Pass; type changed to `Vec<SlotId>`; SLOT_WIDTH now drives off `<SlotId as WireScalar>::FIXED_BYTE_WIDTH` |
| M13 | Major | 4 clones per record in `materialise_next_record` | Apply | Applied | reader.rs | Pass; switched to `mem::take` on the 4 inner `Vec`s |
| M14 | Major | No fuzz target / reader-side proptest | Defer | Deferred | — | Follow-up: add `fuzz/fuzz_targets/psp_reader_no_panic.rs` per review §6 M14 |
| M15 | Major | `PspReader::new`'s 7-step open sequence not structurally marked | Defer | Deferred | — | Follow-up: extract steps to named helpers per review §6 M15 |
| M16 | Major | Unchecked `u64 + u64` in trailer error-path | Apply | Applied | reader.rs | Pass; the error-construction sum now uses `checked_sub`+`unwrap_or(0)` |
| Mi1 | Minor | `SnapshotMismatch` reports only lengths, not content | Apply | Applied | errors.rs, reader.rs | Pass; variant now carries `snapshot: Vec<u16>` and `expected: Vec<u16>`, capped at `MAX_SNAPSHOT_SLOTS_IN_ERROR = 8` |
| Mi2 | Minor | `PhaseChainConsistency` asymmetric with writer | Apply | Applied | errors.rs, reader.rs | Pass; added `record_index: u64` (file-global) alongside `record_in_block` |
| Mi3 | Minor | Duplicated `io_err` context strings on seek+read pairs | Apply | Applied | reader.rs | Pass; "trailer seek" / "trailer read" / "block index seek" / "block index read" / "header seek" |
| Mi4 | Minor | `prev.block_offset + 1` could overflow | Apply | Applied | reader.rs | Pass; switched to `saturating_add(1)` (folded into Mi15's windows(2) refactor) |
| Mi5 | Minor | `delta as u32` truncation silently accepts forged deltas | Apply | Applied | reader.rs | Pass; replaced with `u32::try_from(delta).map_err(…ColumnElementDecode…)?` |
| Mi6 | Minor | `region_records` silent-empty semantics undocumented | Apply | Applied | reader.rs | Pass; the docstring now explicitly documents the silent-empty cases |
| Mi7 | Minor | `f8_first_block_non_empty_snapshot` swallows open-time failure | Apply | Applied | reader.rs | Pass; test now asserts the contract on both branches (open-time + iter-time) |
| Mi8 | Minor | `PspReader::trailer` is `pub` but `#[doc(hidden)]` | Apply | Applied | reader.rs | Pass; demoted to `pub(crate)` |
| Mi9 | Minor | No `# Examples` doctests on public methods | Apply | Applied | reader.rs | Pass; doctests on `new()` and `region_records()`, with `# Errors` sections |
| Mi10 | Minor | `decode_one_column` does several distinct things | Defer | Deferred | — | Follow-up: extract `validate_allele_seq_len` / `validate_finite_floats` per review §6 Mi10 |
| Mi11 | Minor | Five round-trip tests share the same skeleton | Defer | Deferred | — | Cosmetic; not blocking |
| Mi12 | Minor | `PspReader::new(Cursor::new(bytes))` boilerplate 25× | Defer | Deferred | — | Cosmetic; not blocking |
| Mi13 | Minor | "Slice N" obsolete comments throughout file | Apply | Applied | reader.rs | Pass; behaviour-relative phrasing in `records()` docstring, `Slice N` markers dropped from comments on remaining fields and tests |
| Mi14 | Minor | Two near-identical record/allele test helpers | Defer | Deferred | — | Cosmetic; not blocking |
| Mi15 | Minor | Hand-rolled enumerate + strict-monotonic check | Apply | Applied | reader.rs | Pass; monotonic check is now a separate `windows(2)` pass |
| Mi16 | Minor | `match → unwrap_or` for value-or-default | Apply | Applied | reader.rs | Pass |
| Mi17 | Minor | Region-mode peek-ahead reads `RangeClamp` internals | Apply | Applied | reader.rs | Pass; new `RangeClamp::record_past_window` / `record_before_window` / `block_past_window` methods |
| Mi18 | Minor | `u32` iterator counters force `as usize` everywhere | Defer | Deferred | — | Not worth the diff churn; counter field type retained |
| Mi19 | Minor | Test names re-introduce plan-group prefixes | Apply | Applied | reader.rs | Pass; mechanical sed `r<N>_`/`b<N>_`/`f<N>_` rename + manual unstick of `b5_r9_phase_chain_across_block_boundary` → `phase_chain_across_block_boundary` |
| Nits | Nit | Various — see review §6 Nits | Selective Apply | Mixed | reader.rs | Renamed `first_block_overlapping` → `find_first_overlapping_block` (Nit on naming); added the four-bullet ASCII comment at the top of `RecordsIter::next`'s `loop`; `BLOCK_HEADER_READ_CAP.saturating_sub(buf.len())` for the chunk-doubling tightening. Remaining Nits (decode_one_column rename, predict_fixed_uncompressed_len rename, RangeClamp::None/Sequential rename, tuple `<=` comment, `materialise_next_record` style) deferred — too many borderline preferences to mass-apply. |

## 3. Questions asked and answers

None — applied autonomously per the project rule.

## 4. Per-finding log

### B1 — Per-block manifest required-tag coverage missing

- **Severity:** Blocker
- **Final status:** Applied
- **Reasoning:** Two `.expect()` panic surfaces (the 12 at the end of `decode_block_payload`, the one inside `decode_one_column`'s `AlleleSeq` arm) close in one move. The new error variant `PspReadError::MissingRequiredColumnInManifest { name, tag }` is symmetric with the existing `MissingRequiredColumn` for the file-level TOML side.
- **Implementation summary:** `decode_block_payload` walks `V1_0_COLUMNS` and verifies every `required` column is present in `header.manifest` before entering the per-column decode loop. Test pins this for every required tag by building a manifest that omits exactly one column at a time.
- **Files changed:** `src/per_sample_caller/psp/errors.rs`, `src/per_sample_caller/psp/reader.rs`
- **Tests added:** `manifest_missing_required_column_for_every_required_tag`, `manifest_missing_required_column` (renamed from F7)
- **Validation:** all 176 psp tests pass.

### B2 — Untrusted `compressed_len` drives 4 GiB allocation

- **Severity:** Blocker
- **Final status:** Applied
- **Reasoning:** The block's geographic extent is known from the block index (next entry's `block_offset` or `trailer.index_offset` for the last block). Threading a `byte_budget` through `decode_block_payload` → `decode_one_column` lets the reader reject a hostile `compressed_len = u32::MAX` before any allocation.
- **Implementation summary:** New `block_byte_budget(index, trailer, block_idx)` helper; `decode_one_column` takes a `remaining_budget: u64` and rejects with `PspReadError::ColumnTruncated` if `compressed_len > remaining_budget`. Each iteration of the manifest loop in `decode_block_payload` subtracts the consumed `compressed_len`.
- **Files changed:** `src/per_sample_caller/psp/reader.rs`
- **Tests added:** `decode_one_column_rejects_oversized_compressed_len`
- **Validation:** all 176 psp tests pass.

### B3 — `sum(n_alleles) == n_total_alleles` unenforced

- **Severity:** Blocker
- **Final status:** Applied
- **Reasoning:** Without this check, an over-run panics on `block.allele_seqs[j]` indexing; an under-run silently emits truncated allele lists. New `BlockHeaderInvariantKind::NAllelesSumMismatch { n_total_alleles, sum_n_alleles }` returns a typed error.
- **Implementation summary:** Added after the per-column dispatch loop in `decode_block_payload`. Test fixture rewrites the `n_total_alleles` varint byte in a real round-trip file.
- **Files changed:** `src/per_sample_caller/psp/errors.rs`, `src/per_sample_caller/psp/reader.rs`
- **Tests added:** `n_alleles_sum_mismatch_against_n_total_alleles`
- **Validation:** all 176 psp tests pass.

### M1 — `materialise_next_record` `.expect` without PANIC-FREE comment

- **Severity:** Major
- **Final status:** Applied with adaptation
- **Reasoning:** The review suggested either a `// PANIC-FREE:` comment or restructuring the function to take `&DecodedBlock` directly. I chose a third path: a `match self.cur_block.as_mut() { Some(b) => b, None => return Err(Io { … }) }` that converts the unreachable case to a typed error (so even if a future caller violates the precondition, the result is a `Result`-shaped failure, not a panic). Cost is a one-time `std::io::Error::other("internal invariant")` allocation — negligible relative to the per-record decode.
- **Files changed:** `src/per_sample_caller/psp/reader.rs`
- **Validation:** all 176 psp tests pass.

### M2 — `try_into().unwrap()` on slice

- **Severity:** Major. **Final status:** Applied. Replaced with `copy_from_slice` into a stack `[u8; 8]`. **Files:** reader.rs.

### M3 — `VarintError::Overflow` wrong source

- **Severity:** Major. **Final status:** Applied. New variant `PspReadError::BlockHeaderExceedsCap { cap, consumed }`. **Files:** errors.rs, reader.rs.

### M4 — `BlockHeaderField` `field` is a description, not a field name

- **Severity:** Major. **Final status:** Applied. New variant `PspReadError::BlockHeaderTruncated { offset, consumed }`. **Files:** errors.rs, reader.rs.

### M5 — `BlockIndexOffsetInvalid` `min: 0` misleading

- **Severity:** Major. **Final status:** Applied. The three offset/chrom/pos checks now run after the header parse (step 5 in the open sequence); every offset-failure error carries the honest half-open interval `[header_end, trailer.index_offset)`. The monotonic check is a separate `windows(2)` pass (Mi15). **Files:** reader.rs.

### M6 — `NonFiniteFloat` column name hard-coded; F5 deferred test

- **Severity:** Major. **Final status:** Applied. The non-finite sweep is now driven by `ColumnDef::finite_constraint` (M11); when it fires, the column name comes from `def.name` rather than a hard-coded literal. The F5 fixture-patch test requires a `replace_column_payload_f64` helper that round-trips a writer output and re-zstds one column — captured as follow-up but not landed here. **Files:** reader.rs.

### M7 — Spec check #10 first-block-of-chrom snapshot empty

- **Severity:** Major. **Final status:** Applied. New `BlockHeaderInvariantKind::NonEmptySnapshotAtChromStart { block, snapshot_len }`; check runs in both `RangeClamp::None` and `RangeClamp::Window` arms of `load_next_block`. Test: `region_mode_rejects_non_empty_snapshot_at_chrom_start`. **Files:** errors.rs, reader.rs.

### M8 — `BTreeSet<SlotId>` → sorted `Vec<SlotId>`

- **Severity:** Major. **Final status:** Applied. `RecordsIter::active_chain_slots: Vec<SlotId>`; membership via `binary_search`; snapshot equality is a direct `block_header.active_chain_slots == self.active_chain_slots`; region-mode re-init is `clear()` + `extend_from_slice`. **Files:** reader.rs.

### M9 — Three `#[allow(dead_code)]` fields are dead

- **Severity:** Major. **Final status:** Applied. Removed `PspReader::header_end_offset`, `DecodedBlock::n_total_alleles`, `DecodedBlock::active_chain_slots_snapshot` (the last drops one `Vec<SlotId>` clone per block). **Files:** reader.rs.

### M10 — `DecodedColumn::Unknown` catch-all

- **Severity:** Major. **Final status:** Applied. `decode_one_column` returns `Option<DecodedColumn>`; the `None` case is the unknown-tag skip. `DecodedColumn` is now closed over the v1.0 registry with no `Unknown` variant. Adding a `ColumnKey` requires both a new `decode_one_column` match arm *and* a matching `DecodedColumn` variant. **Files:** reader.rs.

### M11 — F5 finite-float sweep hard-coded

- **Severity:** Major. **Final status:** Applied. New `ColumnDef::finite_constraint: bool`; set `true` on `allele-q-sum-log`, `false` on every other v1.0 column. `decode_block_payload`'s sweep is now registry-driven. **Files:** registry.rs, reader.rs.

### M12 — `BlockHeader::active_chain_slots: Vec<u16>` not tied to `SlotId`

- **Severity:** Major. **Final status:** Applied. Type is now `Vec<SlotId>`; `SLOT_WIDTH` is computed from `<SlotId as WireScalar>::FIXED_BYTE_WIDTH` (rather than `<u16 as WireScalar>::FIXED_BYTE_WIDTH`); decode now calls `SlotId::from_le_bytes(arr)` (since `SlotId = u16` today, this is a no-op rename; widening `SlotId` would be a compile error at this site rather than a silent wire-format break). **Files:** block.rs.

### M13 — Four clones per record in `materialise_next_record`

- **Severity:** Major. **Final status:** Applied. Switched `block.allele_seqs[j].clone()`, `block.allele_chain_slots[j].clone()`, `block.new_chain_slots[i].clone()`, `block.expired_chain_slots[i].clone()` to `std::mem::take`. The block is consumed forwards-only (the iterator never revisits `(i, j)`) and dropped wholesale once exhausted, so the moved-out vectors are dead the moment the record is emitted. Validation reads stay before the `mem::take` loop. **Files:** reader.rs.

### M14 — No fuzz target / reader-side proptest

- **Severity:** Major. **Final status:** Deferred. Setting up `cargo-fuzz` adds dev-dependencies, a new `fuzz/` subdirectory, and CI integration choices that go beyond this review's scope. Follow-up.

### M15 — `PspReader::new` 7-step open sequence not structurally marked

- **Severity:** Major. **Final status:** Deferred. The 7-step refactor is a non-trivial extraction (5 named helpers) and not a correctness fix. The numeric step comments now read 1-6 after M5 collapsed steps 4 + 6 of the original into the new step 5. Follow-up.

### M16 — Unchecked `u64 + u64` in trailer error path

- **Severity:** Major. **Final status:** Applied. The error-construction expression now uses `computed_end.and_then(|end| trailer_start.checked_sub(end)).unwrap_or(0)` so no path panics in debug or wraps in release on `u64::MAX`-style trailers. **Files:** reader.rs.

### Mi1 — `SnapshotMismatch` reports only lengths

- **Severity:** Minor. **Final status:** Applied. `BlockHeaderInvariantKind::SnapshotMismatch` now carries `snapshot: Vec<u16>` and `expected: Vec<u16>`, capped at the new `MAX_SNAPSHOT_SLOTS_IN_ERROR = 8` so log lines stay bounded on pathological inputs. New `cap_slots_for_error` helper. **Files:** errors.rs, reader.rs.

### Mi2 — `PhaseChainConsistency::record_in_block` asymmetric with writer

- **Severity:** Minor. **Final status:** Applied. `PspReadError::PhaseChainConsistency` now carries `record_index: u64` (file-global) alongside `record_in_block: u32` (block-local). `RecordsIter::record_index` is incremented after every successful `materialise_next_record` call. **Files:** errors.rs, reader.rs.

### Mi3 — Duplicated `io_err` context strings

- **Severity:** Minor. **Final status:** Applied. "trailer seek" / "trailer read" / "block index seek" / "block index read" / "header seek" — every seek/read pair now has a distinct context label. **Files:** reader.rs.

### Mi4 — `prev.block_offset + 1` could overflow

- **Severity:** Minor. **Final status:** Applied. Folded into Mi15's windows(2) refactor: `min: w[0].block_offset.saturating_add(1)`. **Files:** reader.rs.

### Mi5 — `delta as u32` truncation

- **Severity:** Minor. **Final status:** Applied. `u32::try_from(delta).map_err(|_| PspReadError::ColumnElementDecode { source: ScalarDecodeError::VarintOverflow, … })?`. Test path is the new compile-check `delta_pos_u32_overflow_compile_check` (the end-to-end fixture needs a varint-patch helper deferred with M6's F5 test). **Files:** reader.rs.

### Mi6 — `region_records` silent-empty docstring

- **Severity:** Minor. **Final status:** Applied. The docstring now explicitly documents the three silent-empty cases: out-of-range `chrom_id`, no-overlap window, and `start > end`. **Files:** reader.rs.

### Mi7 — `f8` swallows open-time failure

- **Severity:** Minor. **Final status:** Applied. The test now branches on `match PspReader::new(…)`: in the `Ok` arm it expects an iteration-time error (with the M7 / M4 / M3 new variants in the accept-set); in the `Err` arm it pins the contract on the open-time variant. **Files:** reader.rs.

### Mi8 — `PspReader::trailer` `pub` + `#[doc(hidden)]`

- **Severity:** Minor. **Final status:** Applied. Demoted to `pub(crate)`. `#[cfg_attr(not(test), allow(dead_code))]` because it's currently only called from tests. **Files:** reader.rs.

### Mi9 — No `# Examples` doctests

- **Severity:** Minor. **Final status:** Applied. Doctests on `PspReader::new` and `PspReader::region_records`, each with a `# Errors` section listing the open-time failure variants. Both compile (run as `compile`-only with `no_run` since they reference `sample.psp` on disk). **Files:** reader.rs.

### Mi10 — `decode_one_column` extract validators

- **Severity:** Minor. **Final status:** Deferred. Cleanup-class; the function is now 138 lines, but with M10 + M11 the per-column structure already reads more cleanly. Extracting `validate_allele_seq_len` / `validate_finite_floats` is a follow-up.

### Mi11 — Five round-trip tests share the skeleton

- **Severity:** Minor. **Final status:** Deferred. Cosmetic.

### Mi12 — 25× `PspReader::new(Cursor::new(…))` boilerplate

- **Severity:** Minor. **Final status:** Deferred. Cosmetic.

### Mi13 — "Slice N" obsolete comments

- **Severity:** Minor. **Final status:** Applied. The `records()` docstring no longer references Slices 2-5; the comment on `RangeClamp::Window` no longer says "Slice 5"; the `f8` test docstring no longer references Slice 6; the comment on `DecodedColumn::Unknown` is gone with the variant (M10). One residual: the test section heading `// ---------- Slice 6: F1-F8 corruption / negative tests ------` is gone with the Mi19 rename (the comment was edited there too). **Files:** reader.rs.

### Mi14 — Test helpers consolidation

- **Severity:** Minor. **Final status:** Deferred. Cosmetic.

### Mi15 — Hand-rolled enumerate + monotonic check

- **Severity:** Minor. **Final status:** Applied. The monotonic check now sits in a separate `windows(2)` pass, the range check is a separate enumerate pass — convergent with M5's restructuring. **Files:** reader.rs.

### Mi16 — `match` → `unwrap_or`

- **Severity:** Minor. **Final status:** Applied. `it.cur_block_idx = first.unwrap_or(it.reader.index.len());`. **Files:** reader.rs.

### Mi17 — `RangeClamp` methods

- **Severity:** Minor. **Final status:** Applied. `RangeClamp::record_past_window` / `record_before_window` / `block_past_window`; `RecordsIter::next` now reads as `self.clamp.<predicate>(…)` calls — both inline `if let RangeClamp::Window { … } = self.clamp` destructures removed. **Files:** reader.rs.

### Mi18 — `u32` counter types

- **Severity:** Minor. **Final status:** Deferred. The diff would touch ~6 sites for no functional gain.

### Mi19 — Plan-group prefixes on test names

- **Severity:** Minor. **Final status:** Applied. Mechanical `sed -i 's/fn r\([0-9]\+\)_/fn /g; s/fn b\([0-9]\+\)_/fn /g; s/fn f\([0-9]\+\)_/fn /g'` rename plus one manual rename for `b5_r9_phase_chain_across_block_boundary` → `phase_chain_across_block_boundary`. **Files:** reader.rs.

### Nits

- **`first_block_overlapping` → `find_first_overlapping_block`** — Applied (verb-shape consistency with `read_block_header`, `decode_block_payload`, `decode_one_column`).
- **`RecordsIter::next` state-machine ASCII comment** — Applied (four-bullet comment at the top of the `loop` body labelling the phases).
- **`BLOCK_HEADER_READ_CAP - buf.len()` → `BLOCK_HEADER_READ_CAP.saturating_sub(buf.len())`** — Applied (defence-in-depth).
- **`decode_one_column` → `decode_column`, `predict_fixed_uncompressed_len` rename, `RangeClamp::None` → `Sequential`, tuple `<=` comment, `materialise_next_record` style, `io_err` naming** — Deferred (borderline stylistic preferences; not worth the cross-file churn given the volume of substantive changes in this round).

## 5. Deferred findings to carry forward

- **M14** — cargo-fuzz target for `PspReader::new` + `RecordsIter::next`. Add `fuzz/fuzz_targets/psp_reader_no_panic.rs` per review §6 M14. Triage every panic / OOM / hang to a typed error.
- **M15** — extract `PspReader::new`'s open sequence into five named helpers (`read_trailer`, `read_index`, `parse_header`, `cross_check_index`, `seek_to_first_block`).
- **M6 follow-up** — F5 end-to-end test: build a single-record file, replace the `allele-q-sum-log` column payload (8 bytes) with `f64::NAN.to_le_bytes()`, re-zstd, re-patch manifest `compressed_len`, re-checksum the trailer. Requires the test-only `replace_column_payload_f64` helper that the implementation report's §6 has been carrying.
- **Mi10** — `validate_allele_seq_len` / `validate_finite_floats` extraction.
- **Mi11 / Mi12 / Mi14** — test boilerplate consolidation (`round_trip` helper, `open_reader` helper, single `record(…)` + `allele(…)` constructor).
- **Mi18** — `u32` → `usize` counter type swap.
- **Reliability §8 still-missing tests** — items 6 (`next_returns_none_after_first_error_is_yielded`), 7 (`read_block_header_returns_overflow_error_on_pathological_continuation_stream`), 8 (`first_block_overlapping_finds_first_block_when_window_starts_before_it`), 10 (`pspreader_new_rejects_trailer_with_overflowing_arithmetic_without_panicking`). All four are mechanical to add but expand the test surface further than this round was scoped for.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** No.
- **Reason:** No `benches/psp_reader*.rs` exists yet — there is no baseline to compare against. M13's `mem::take` switch is plausibly performance-positive (4 fewer mallocs per record) but the gain has not been measured here.
- **Follow-up:** Add a `criterion` micro-bench at `benches/psp_reader_perf.rs` reading a 1M-record file. The next perf review can use it to bisect any regression introduced by this round.

## 10. Commands run

- `cargo build --lib --tests --all-features`
- `cargo fmt --check`
- `cargo fmt`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib --tests --all-features`
- `cargo test --doc --all-features per_sample_caller::psp::reader`
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib`

## 11. Command results

- `cargo build --lib --tests --all-features` → exit 0, clean
- `cargo fmt --check` → exit 0 after `cargo fmt` ran once
- `cargo clippy --all-targets --all-features -- -D warnings` → exit 0, clean
- `cargo test --lib --tests --all-features` → exit 0, **562 passed / 0 failed** (474 lib + 88 integration)
- `cargo test --doc --all-features per_sample_caller::psp::reader` → exit 0, 2 doctests pass
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib` → exit 0, clean

## 12. Notes

- **Cargo / dev container.** Project policy says cargo runs inside `./scripts/dev.sh`, but `podman` is not installed on this host. Cargo was invoked directly with `cargo` from the host's toolchain (Rust 1.95.0); `CARGO_TARGET_DIR` defaulted to the project's `target/` (the container path `target-container/` is also written by the test harness — both directories now exist).
- **Public API changes that landed.**
  - `PspReader::trailer` demoted from `pub` to `pub(crate)` (Mi8). Not used by anything outside the crate.
  - `BlockHeader::active_chain_slots: Vec<u16>` → `Vec<SlotId>` (M12). `SlotId = u16` today so this is a no-op for callers; type-tightening for forward-compat.
  - `PspReadError`: 4 new variants (`MissingRequiredColumnInManifest`, `BlockHeaderExceedsCap`, `BlockHeaderTruncated`), 1 new sub-enum variant (`BlockHeaderInvariantKind::NAllelesSumMismatch`, `NonEmptySnapshotAtChromStart`). All additive; the enum carries `#[non_exhaustive]`.
  - `BlockHeaderInvariantKind::SnapshotMismatch`: field types changed from `snapshot_len: usize, expected_len: usize` to `snapshot: Vec<u16>, expected: Vec<u16>` (Mi1). Source-breaking for anyone matching on the variant; the existing test uses `{ .. }` and survives.
  - `PspReadError::PhaseChainConsistency`: gained a `record_index: u64` field (Mi2). Same source-breaking caveat; the one matching site uses `{ .. }`.
  - `ColumnDef`: new public field `finite_constraint: bool` (M11). Source-breaking for code that constructs `ColumnDef` literals — only the registry does so, and every entry has been updated.
  - The user's prior directive "Don't worry about the API change" covered the deferred-fixes wave; this round inherits the same posture since the reader was added in the same week and is still pre-stable.
- **F7 test regression.** The original F7 test (forge a single-column manifest with a wrong `uncompressed_len`) no longer reaches `UncompressedLenSchemaMismatch` because B1's coverage check fronts the per-column decode loop. The test was renamed to `manifest_missing_required_column` and now asserts the new path. A dedicated `UncompressedLenSchemaMismatch` regression test requires a fixture-patch helper that round-trips a real file and overwrites one manifest entry's `uncompressed_len` byte; left as follow-up.
- **Verification quoted vs. re-run.** All validation commands were re-run after the fix sequence (not quoted from the review). Test counts: 27 reader tests (22 original + 5 new from this round), 176 psp tests overall, 474 lib tests overall, 562 with integration.
