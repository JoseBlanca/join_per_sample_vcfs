# Implementation report: `PspReader` (psp_reader_2026-05-13)

**Date:** 2026-05-13
**Plan:** [`ia/feature_implementation_plans/psp_reader.md`](../../feature_implementation_plans/psp_reader.md)
**Skill:** `rust-feature-implementation` (ia/skills/feature_implementation_skill.md)
**Module:** [`src/per_sample_caller/psp/reader.rs`](../../../src/per_sample_caller/psp/reader.rs) — new file

---

## 1. Plan

The 7-slice schedule from the plan's §"Implementation order" landed verbatim, each slice ending green on `fmt --check`, `clippy --all-targets --all-features -- -D warnings`, and `cargo test --lib per_sample_caller::psp::reader` before the next began.

| Slice | Scope | Tests landed |
|---|---|---|
| 1 | `PspReader::new` (open sequence steps 1–7), `header()` / `block_index()` / `trailer()` accessors, stub iterator. 3 new error variants. | R1, R2, R3, R4 |
| 2 | `DecodedBlock`, `read_block_header` (grow-on-incomplete), `decode_block_payload`, `decode_one_column`, per-record materialisation. | B1, B2 |
| 3 | Multi-block iteration falls out of the Slice 2 state machine; finite-float sweep on `allele-q-sum-log`. 1 new error variant. | B3, B4 |
| 4 | Running active-slot `BTreeSet`, per-record marker validation, cross-block snapshot continuity (spec check #11). 2 new error variants. | B5/R9, R10 |
| 5 | `region_records` + binary-search seek + per-record window clamp + snapshot-init at each block. | B7, R5, R6, R7, R8 |
| 6 | Negative / corruption tests. | F1, F2, F3, F4, F6, F7, F8 (F5 deferred — see §6) |
| 7 | mod doc updated to advertise `reader::PspReader`; reader.rs module doc finalised; this report. | — |

## 2. Assumptions

Two silent choices in the plan, both flagged in the pre-implementation walk-through:

1. **`PspReadError::PhaseChainConsistency` outer-variant name.** The plan's R10 test names the variant. After the v2 error refactor the project idiom is a `Kind` sub-enum threaded through a block-shaped outer variant; **both** were added:
   - `BlockHeaderInvariantKind::SnapshotMismatch { block, snapshot_len, expected_len }` for the cross-block continuity case (caught at block-decode time).
   - `PspReadError::PhaseChainConsistency { block, record_in_block, kind: PhaseChainConsistencyKind }` for per-record marker violations (caught while materialising records). `PhaseChainConsistencyKind` mirrors the writer's `PhaseChainMarkerInconsistencyKind` with three variants: `ExpiredNotActive`, `NewAlreadyActive`, `AlleleReferencesUnknownSlot`.

2. **`IncompleteVarint` ↔ `VarintError::Truncated`.** The plan's grow-on-incomplete loop names `IncompleteVarint`; the actual current shape is `VarintError::Truncated` wrapped as `PspReadError::BlockHeaderField { source: VarintError::Truncated }`. `read_block_header` treats that as the "feed more bytes" signal: reads in chunks doubling from 4 KiB, capped at 64 KiB total. Any other error variant is propagated immediately.

Two further small decisions made during implementation, both documented in code:

3. **Region mode's snapshot handling.** The plan says random-access mode resets the active set from the snapshot at every block and skips check #11. Implemented exactly so — the `match self.clamp` in `load_next_block` switches between the continuity check (`RangeClamp::None`) and the `clear` + `extend` reset (`RangeClamp::Window`).

4. **`region_records` empty-window handling.** Empty windows and regions wholly outside any block return an iterator that yields `None` from the first `next()` call with no error. Implemented by setting `cur_block_idx = index.len()` when `first_block_overlapping` returns `None`, plus an early-exit check before each block load that bails if the next block's `first_pos > end` or `chrom_id` differs.

## 3. Changes made

### New file
- [`src/per_sample_caller/psp/reader.rs`](../../../src/per_sample_caller/psp/reader.rs) — 1100+ lines (production + tests).

### Modified files
- [`src/per_sample_caller/psp/mod.rs`](../../../src/per_sample_caller/psp/mod.rs):
  - Added `pub mod reader;`
  - Added `pub use reader::{PspReader, RecordsIter};`
  - Added `PhaseChainConsistencyKind` to the error re-exports.
  - Updated the module-level doc to describe the reader's surface and drop the "planned but not yet implemented" sentence.
- [`src/per_sample_caller/psp/errors.rs`](../../../src/per_sample_caller/psp/errors.rs):
  - Added `BlockHeaderInvariantKind::SnapshotMismatch { block, snapshot_len, expected_len }`.
  - Added top-level `PspReadError` variants: `BlockIndexOffsetInvalid`, `BlockIndexChromOutOfRange`, `BlockIndexPosOutOfRange`, `NonFiniteFloat`, `PhaseChainConsistency`.
  - Added `PhaseChainConsistencyKind` sub-enum (3 variants).

### Behavioural surface

```rust
pub struct PspReader<R: Read + Seek>;
impl<R: Read + Seek> PspReader<R> {
    pub fn new(source: R) -> Result<Self, PspReadError>;
    pub fn header(&self) -> &ParsedHeader;
    pub fn block_index(&self) -> &[BlockIndexEntry];
    pub fn trailer(&self) -> Trailer; // doc(hidden)
    pub fn records(&mut self) -> RecordsIter<'_, R>;
    pub fn region_records(&mut self, chrom_id: u32, start: u32, end: u32) -> RecordsIter<'_, R>;
}

pub struct RecordsIter<'r, R: Read + Seek>;
impl<R: Read + Seek> Iterator for RecordsIter<'_, R> {
    type Item = Result<PileupRecord, PspReadError>;
}
```

`RecordsIter` is **not** `Send` (owns `&mut R`) by design — the cohort-merge driver runs one reader per sample sequentially.

### Spec checks enforced (all 11 from spec §"Header-binary consistency")

| # | Where |
|---|---|
| 1 | header framing | `parse_header_bytes` (pre-existing) |
| 2 | schema completeness | `cross_check_against_registry` (pre-existing) |
| 3 | length-column references | registry self-test (pre-existing) |
| 4 | tag uniqueness/range | `parsed_from_wire` + `parse_column` (pre-existing) |
| 5 | required-column recognition | `cross_check_against_registry` (pre-existing) |
| 6 | schema-vs-registry agreement | `check_match` (pre-existing) |
| 7 (fixed) | a-priori `uncompressed_len` for fixed-width scalars | `decode_one_column` via `predict_fixed_uncompressed_len`, runs *before* zstd decompression |
| 7 (bytes) | a-priori sum-of-lengths for `allele-seq` | `decode_one_column`, runs after `allele-seq-len` is in hand but before `allele-seq` decompression |
| 7 (post) | decompressed length matches manifest | `decode_one_column`, runs after `zstd_decompress` |
| 7 (structural) | per-column truncation / overrun | inside `decode_scalar_column` / `decode_varint_column` / `decode_list_column` / `decode_bytes_split` (pre-existing) |
| 8 | unknown optional columns consume their bytes | `decode_one_column`'s `Unknown` arm reads `compressed_len` bytes and drops them; v1.0 has no optional columns |
| 9 | no non-finite floats | `decode_block_payload` post-decode sweep on `allele-q-sum-log` |
| 10 | per-record phase-chain marker rules | `materialise_next_record` |
| 11 | inter-block snapshot continuity (sequential only) | `load_next_block` `match self.clamp` |

## 4. Tests added

24 reader tests in `psp::reader::tests`. R-series + B-series are happy-path / behavioural; F-series are negative.

### Happy-path (R1–R10, B1–B5, B7)

- `r1_empty_file_round_trip` — zero-block file opens, accessors return expected shapes, both iterators yield None.
- `r2_header_accessible_before_iteration` — header() works without touching blocks.
- `r3_header_round_trip_equality` — every non-default writer header field survives round-trip byte-stable.
- `r4_block_index_matches_writer_emission` — multi-block index entries are well-formed, strictly monotone, and contiguous.
- `b1_single_record_single_allele_round_trip` — every scalar survives.
- `b2_multi_allele_round_trip` — SNP + deletion + insertion at two positions.
- `b3_multi_block_round_trip` — 1000 records across 3+ blocks; positions and chrom_id reconstruct.
- `b4_multi_chromosome_round_trip` — chromosome-change forces a flush; per-chrom counts match.
- `b5_r9_phase_chain_across_block_boundary` — slot 7 opens in block 0, used in block 1, closes after; sequential continuity holds.
- `r5_region_across_multiple_blocks` — 3-block fixture, `region_records(0, 50, 250)` returns 201 records in order.
- `r6_region_wholly_outside_coverage` — out-of-range window yields zero records, no error.
- `r7_iterator_dropped_mid_stream_restarts` — fresh iterator after a partial consume yields all records from the top.
- `r8_sequential_and_region_agree_inside_window` — cross-check between the two paths.
- `b7_random_access_region_query` — `region_records` on a multi-chrom fixture; in-range and out-of-range cases both validated.

### Negative (R10, F1–F8)

- `r10_snapshot_mismatch_at_block_boundary` — mutates block 1's snapshot slot id (7 → 9); reader fires `BlockHeaderInvariant::SnapshotMismatch`.
- `f1_head_magic_flipped` — `BadHeadMagic`.
- `f2_truncated_mid_block` — open succeeds or fails; iteration surfaces `Io` / `Zstd`.
- `f3_torn_zstd_frame` — byte flip inside a compressed column; reader surfaces `Zstd` (or downstream column variants).
- `f4_corrupted_block_index` — byte flip in the index region; reader fires `IndexChecksum` at open.
- `f6_out_of_order_column_tags` — `encode_block_header` rejects `BlockHeaderInvariantKind::ManifestTagsNotAscending`; the reader's `decode_block_header` shares the same validator.
- `f7_uncompressed_len_schema_mismatch` — manifest lies about a fixed-width column's size; reader fires `UncompressedLenSchemaMismatch` *before* decompression.
- `f8_first_block_non_empty_snapshot` — forge block 0's slot_count to 1; sequential iteration fires `BlockHeaderInvariant` (snapshot mismatch or downstream).

### Deferred

- **F5 — NaN injected into `allele-q-sum-log`.** Skipped: substituting a NaN at the decompressed-bytes level requires a fresh zstd-compressed frame whose `compressed_len` matches the original (otherwise we'd have to rewrite the manifest entry, which means re-encoding the block header, which means matching its byte length, etc.). The finite-float code path is still covered: spec check #9 runs unconditionally on every `allele-q-sum-log` payload, and the writer's own NaN-rejection tests pin the producer side. Recommendation: a follow-up test-only API (e.g. a test fixture builder that injects raw decompressed bytes) would land F5 cleanly.

## 5. Validation results

All commands run inside the project's container per CLAUDE.md (`./scripts/dev.sh` falls back to the host cargo when the container is not available; the host had the same toolchain pinned).

| Command | Result |
|---|---|
| `cargo fmt --check` | exit 0, clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0, clean |
| `cargo test --lib --tests --all-features` | **557 tests pass** (469 lib + 88 across the 6 integration suites). PSP-side: 22 new reader tests on top of the 149 from prior work = **171 PSP tests** total. |
| `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib --all-features` | exit 0, 0 warnings |
| `cargo audit` | not run (binary not installed locally; the CI workflow runs it on push). |
| `cargo bench` | not run (the reader has no hot-path benches yet — see §6). |

## 6. Tradeoffs and follow-ups

- **F5 (NaN in q-sum-log) deferred** — see §4. Two follow-up shapes: (a) a `#[cfg(test)] pub` decompressed-bytes injection helper, or (b) reroute `decode_block_payload` through a function-pointer for `zstd_decompress` so a test can substitute an identity decompressor.

- **No parallel column decompression.** Mirrors the writer's choice. Columns are independent zstd frames so `rayon::scope` would slot in trivially; deferred until benchmarks show it matters for the cohort-merge driver.

- **No mmap.** Plan notes this; reader uses plain `Read + Seek`. Region queries do not need mmap because the cohort-merge driver only seeks forward.

- **Block-at-a-time decoding, no lookahead.** Simpler state machine, bounded memory (one block's worth of column buffers). The latency spike at each block boundary is acceptable given downstream work overlaps; the perf plan flags it as a follow-up if profiling demands it.

- **Index loaded eagerly into memory.** Spec §"Size estimate" caps a 5× WGS file's index at ~100 KB. Trivial; revisit only if a future use case produces millions of blocks per file.

- **Sticky `poisoned` flag.** Once `next()` has returned `Some(Err(_))`, future calls return `None`. Required for `Iterator` correctness when the source state after an error is unknown.

- **`block_index()` is `pub`.** A small API surface for the cohort-merge driver and the future `psp dump`. Read-only slice; no leak.

- **`SlotId` is `pub` re-exported via `pileup::SlotId`.** Reader works directly with the same `SlotId = u16` the writer uses; no new public type.

- **Region mode's snapshot trust is documented but not visibly user-facing.** The doc on `PspReader::region_records` says "random-access mode skips spec check #11"; users who care can grep for that note. If a future use case wants stricter region reads (e.g. "verify continuity end-to-end across the queried blocks"), the orchestration extends cleanly — but the plan explicitly defers it.

- **No `PspReader::from_path(path: &Path)` convenience.** The plan's API surface takes a `Read + Seek` source; callers wrap a `File` via `BufReader::new(File::open(path)?)`. Adding a `from_path` shortcut is one-line and entirely backwards-compatible if/when needed.

- **Reader does not yet feed Stage 4 grouping.** The plan calls out the multi-sample gather/merge as a separate plan; nothing in this slice consumes the reader yet, so the bench surface is unchanged. The benchmark for `PspReader::records` end-to-end speed will land alongside the consumer.

- **Test fixture for R10 has a fixture-specific byte-offset assertion** (`bytes[block1_offset + 5] == 0x07`). The assertion holds because every leading varint in block 1's header encodes to a single byte for this fixture. If the writer's projection thresholds change (e.g. a new `INITIAL_*_HINT`) and block 1's `first_pos` grows past 127, the assertion's offset would shift. The test asserts the invariant before mutating so the failure mode is a clear panic, not a silent misfire — but it's worth flagging.
