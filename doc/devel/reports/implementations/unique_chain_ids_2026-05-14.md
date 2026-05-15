# Phase chains: recycled `u16` → unique `u64` chain ids

Implementation date: 2026-05-14. Plan:
[ia/feature_implementation_plans/unique_chain_ids.md](../../feature_implementation_plans/unique_chain_ids.md).

## Plan

Switch the pipeline's phase-chain identifier mechanism from
recycled `u16` slot ids with per-record `new_chains` /
`expired_chains` lifecycle markers to **unique-per-`.psp`-file
`u64` chain ids** with no recycling. Per the plan, the goal was
correctness-by-construction (the integration finding pinned by the
previous commit could not exist under the new design) and the
removal of a category of complexity that lived across three
modules (walker slot_allocator, psp writer's active-set tracking,
psp reader's active-view tracking, plus the block format's
snapshot field).

## Assumptions / silent choices

- **`u64` over `u32`.** Plan picked `u64` for absolute future-
  proofing; implemented as `pub type SlotId = u64`. Wire encoding
  is fixed-width LE `u64` (8 bytes per id; zstd absorbs leading
  zeros on small ids).
- **Wire encoding for `allele_chain_slots`: fixed `u64`, not
  `varint`.** Plan proposed varint; implementation chose fixed
  `u64` to reuse the existing `encode_list_column<T: WireScalar>`
  / `decode_list_column<T>` paths without adding a varint-list
  codec. Net effect after zstd: probably indistinguishable from
  varint for the dominant cases; the simpler codec wins.
- **Overflow guard:** `allocate_for_read` does a `checked_add(1)`
  on the per-file monotonic counter and surfaces
  `WalkerError::ChainIdSpaceExhausted { chrom_id, pos }` on
  overflow. Test-only constructor
  `SlotAllocator::with_next_id_for_testing(start: u64)` seeds the
  counter near `u64::MAX` so two tests pin both the off-by-one
  boundary and the exhaustion behaviour without literally
  allocating 2⁶⁴ chains.
- **`SlotExhausted` renamed to `ActiveReadsExhausted`.** Mechanical
  rename; the variant still trips on the
  `DEFAULT_MAX_ACTIVE_SLOTS = 4096` cap, but the name is now
  honest about what's exhausted (concurrent active reads, not a
  bounded id space).
- **`active_count` is bumped once per actual admitted read**, not
  once per anticipated pair. The previous "anticipate the partner"
  refcount bump is gone; the orphan-cleanup path
  (`release_pending_partner_ref_if_present`) only clears the
  `pending_mates` entry now. The first mate's own
  `release_slot` decrements `active_count` when it exits.
- **`PileupRecord` shape change.** `new_chains` and
  `expired_chains` removed. `PileupRecord::new` constructor
  signature drops the two arguments — examples and benches
  updated to match.
- **The pre-existing `paired_mates_share_chain_slot_id` test
  was a false positive** that relied on slot-id recycling to
  produce a matching `[0]` literal. Replaced with two tests that
  pin the real semantics: overlapping mates share, gap-past-end
  mates don't (with a `pending_mates`-lifetime note flagging it as
  a follow-up, since today the entry is tied to active-set lifetime
  rather than to the `mate_lookup_window`).
- **`format_version` stays `(1, 0)`.** No file has shipped; this
  is pre-1.0 design churn. The block-format change is incompatible
  with pre-change writers' output, but no such files exist in tree
  or on disk anywhere.

## Changes made

### Code

- **`src/per_sample_caller/pileup/slot_allocator.rs`** — rewrote
  from scratch as a monotonic `u64` counter plus `pending_mates`
  AHashMap. Removed `free` / `pending_free` / `slot_refcount` /
  `new_marks` / `expired_marks`. Added overflow check on every
  increment. Added test-only `with_next_id_for_testing(start)`
  constructor.
- **`src/per_sample_caller/pileup/mod.rs`** — `SlotId = u64`;
  `PileupRecord` loses `new_chains` and `expired_chains` fields;
  `PileupRecord::new` loses the two arguments;
  `WalkerConfig::max_active_slots`'s doc-comment updated; the
  whole "lifecycle markers" doc block on `PileupRecord` removed.
- **`src/per_sample_caller/pileup/walker.rs`** — dropped the
  `drain_lifecycle_marks` + stamp logic from both
  `close_aged_records_into` and `flush_chromosome_into`. Module-
  level doc comment updated.
- **`src/per_sample_caller/pileup/open_record.rs`** —
  `OpenPileupRecord::finalise` no longer populates the gone
  `new_chains` / `expired_chains` fields.
- **`src/per_sample_caller/pileup/active_read_set.rs`** — comment
  updated; the `release_slot` / `release_pending_partner_ref_if_present`
  call shape unchanged (the API stays for the caller's
  convenience, but the implementation simplified).
- **`src/per_sample_caller/pileup/errors.rs`** — renamed
  `SlotExhausted` → `ActiveReadsExhausted`; added
  `ChainIdSpaceExhausted { chrom_id, pos }`.
- **`src/per_sample_caller/psp/registry.rs`** — dropped
  `NewChainSlots` and `ExpiredChainSlots` from the `ColumnKey`
  enum and from the `V1_0_COLUMNS` table; updated
  `AlleleChainSlots`'s element-type from `U16` to `U64`.
- **`src/per_sample_caller/psp/block.rs`** — dropped
  `active_chain_slots` field from `BlockHeader`, the
  encode/decode passes, and the `ActiveSlotsNotAscending`
  validation. `SlotId` import dropped.
- **`src/per_sample_caller/psp/writer.rs`** — dropped
  `IngestState.active_slots`,
  `apply_record_to_block`'s 50-line marker-validation block, the
  chromosome-boundary "active set must be empty" check, the
  `BlockAccumulator.snapshot_active_slots` field, the
  `new_chain_slots` / `expired_chain_slots` `ListColumn`s, the
  `BlockAccumulator::new` snapshot parameter, and the
  `is_first_block_on_chrom` helper.
- **`src/per_sample_caller/psp/reader.rs`** — dropped
  `DecodedBlock.new_chain_slots` / `expired_chain_slots` fields,
  the `RecordsIter.active_chain_slots` running view, the block-
  start snapshot decode, both spec-check #10 / #11 paths in
  `load_next_block`, the lifecycle-marker apply-and-validate
  block in `materialise_next_record`, the
  `cap_slots_for_error` helper, and the
  `MAX_SNAPSHOT_SLOTS_IN_ERROR` import.
- **`src/per_sample_caller/psp/errors.rs`** — dropped
  `PspWriteError::PhaseChainMarkerInconsistency`,
  `PspReadError::PhaseChainConsistency`,
  `PhaseChainMarkerInconsistencyKind`,
  `PhaseChainConsistencyKind`,
  `BlockHeaderInvariantKind::ActiveSlotsNotAscending`,
  `BlockHeaderInvariantKind::SnapshotMismatch`,
  `BlockHeaderInvariantKind::NonEmptySnapshotAtChromStart`, and
  the `MAX_SNAPSHOT_SLOTS_IN_ERROR` constant.
- **`src/per_sample_caller/psp/mod.rs`** — public re-exports
  shrunk: dropped `PhaseChainConsistencyKind` and
  `PhaseChainMarkerInconsistencyKind`.
- **`src/per_sample_caller/pileup_to_psp.rs`** — flipped the
  previous regression-pin test
  (`writer_rejects_walker_output_when_chain_expires_…`) into the
  real roundtrip-parity test
  (`drive_pileup_to_psp_roundtrips_records_through_in_memory_sink`).
  Module-level doc comments updated.

### Examples and benches

- **`examples/dhat_psp_writer.rs`** — `PileupRecord::new` call
  drops the two empty `Vec::new()` args.
- **`benches/psp_writer_perf.rs`**, **`benches/psp_reader_perf.rs`** —
  `build_phase_chain_heavy_records` rewritten to mint unique `u64`
  ids via a local monotonic counter; the previous "active set"
  fixture (256-id sliding window with marker-aware insertion)
  collapsed to a simpler "sliding window of active ids" structure.
  Other call sites' `PileupRecord::new` calls updated.

### Tests

- **`src/per_sample_caller/pileup/slot_allocator.rs`** —
  rewrote the test block. Kept tests for: solo allocation,
  monotonic increasing ids, mate-pair sharing, no-recycling,
  active-cap exhaustion, stale-pending eviction, reset behaviour,
  high-water warning. Added two overflow-guard tests
  (at `u64::MAX` and at `u64::MAX - 1`).
- **`src/per_sample_caller/pileup/tests.rs`** —
  dropped `lifecycle_markers_appear_on_emitted_records` and
  `orphan_first_mate_emits_balanced_lifecycle_marks`; added
  `chain_ids_are_unique_and_monotonically_allocated`,
  `chain_ids_persist_across_chromosome_boundaries`,
  `paired_mates_share_a_single_chain_id`,
  `paired_mates_with_overlapping_positions_share_chain_id`, and
  `paired_mates_with_gap_past_first_mate_get_distinct_chain_ids`.
  Replaced the old `paired_mates_share_chain_slot_id` test (which
  relied on slot recycling to fake a match).
- **`src/per_sample_caller/pileup/active_read_set.rs`** —
  `flush_all_releases_every_slot` renamed to
  `flush_all_drops_every_active_read`; dropped the
  `drain_lifecycle_marks` assertions that no longer apply.
- **`src/per_sample_caller/psp/block.rs`** — dropped 3 tests
  pinned to `active_chain_slots` invariants (unsorted, duplicate,
  giant-count). Updated the wire-layout pin test and the property
  test to drop the snapshot field.
- **`src/per_sample_caller/psp/writer.rs`** — dropped
  `rejects_phase_chain_marker_inconsistency_expired_not_active`,
  `rejects_phase_chain_marker_inconsistency_double_open`,
  `rejects_allele_chain_slot_not_in_active_set`,
  `active_slot_snapshot_carries_across_blocks`, and
  `rejects_chain_active_at_chromosome_boundary`. The `record()`
  helper signature dropped the two trailing `Vec<u16>` args; the
  `allele()` helper switched its `chain_slots` parameter from
  `&[u16]` to `&[u64]`. All `record(...)` callers shrunk from 5
  args to 3 via a Python-rewrite pass.
- **`src/per_sample_caller/psp/reader.rs`** — dropped
  `snapshot_mismatch_at_block_boundary`,
  `first_block_non_empty_snapshot`, and
  `region_mode_rejects_non_empty_snapshot_at_chrom_start`. Renamed
  `phase_chain_across_block_boundary` →
  `chain_id_round_trips_across_block_boundary` and simplified.
  `record_open_slot` / `record_use_slot` / `record_close_slot`
  test-builder triple collapsed to a single
  `record_with_chain_id` helper.
- **`src/per_sample_caller/psp/header.rs`** —
  `reader_rejects_chain_slots_per_allele_to_per_record` renamed
  and re-pointed at `n-alleles` (since `new-chain-slots` is gone).

### Specs

- **`ia/specs/phase_chain.md`** — rewrote §6 ("The chain-id
  encoding, with a worked example") to describe the unique-`u64`
  design; updated the worked example to drop `new_chains` /
  `expired_chains` rows; rewrote §8's "How Stage 5 walks the
  records" (no active-view bookkeeping); updated the side-by-side
  comparison table; replaced the recap's "two design decisions"
  paragraph with the one decision that remains.
- **`ia/specs/per_sample_pileup_format.md`** — dropped the
  `new-chain-slots` / `expired-chain-slots` `[[column]]` entries
  from the header example; the `allele-chain-slots` entry now
  declares `element-type = "u64"`; dropped the
  `active_chain_slots_at_block_start` field from the block-layout
  pseudocode; rewrote consistency check #10 (well-formedness
  only); deleted check #11; rewrote §"Phase-chain state across
  blocks" to a single paragraph plus an explanatory note; deleted
  the snapshot invariants from §"Block invariants"; rewrote
  Q-FC8 and Q-PL4 with `**Superseded 2026-05-14**` markers; the
  column-tag cheat-sheet now lists only `0x22` for chain ids with
  the LE-`u64` payload description; the prose enumeration of
  v1.0 features updated to "per-allele phase chain identifiers —
  unique-per-file `u64`s, no recycling, no lifecycle markers".
- **`ia/specs/calling_pipeline_architecture.md`** — narrative
  already uses "phase chain identifier" abstractly; no
  implementation-level changes needed.

## Validation

All inside the dev container:

- `cargo fmt --check` — clean.
- `cargo clippy --all-targets --all-features -- -D warnings` —
  clean.
- `cargo test --tests --lib --examples --all-features` — **all
  green**:
  - 466 unit tests in `src/`,
  - 25 + 26 + 17 + 8 + 12 = 88 integration tests in `tests/`,
  - 0 doc/example tests.

Notable green tests pinning the new design:

- `pileup::slot_allocator::allocator_errors_when_chain_id_counter_would_overflow`
- `pileup::slot_allocator::allocator_succeeds_at_chain_id_counter_one_below_max`
- `pileup::tests::chain_ids_are_unique_and_monotonically_allocated`
- `pileup::tests::chain_ids_persist_across_chromosome_boundaries`
- `pileup::tests::paired_mates_with_overlapping_positions_share_chain_id`
- `pileup::tests::paired_mates_with_gap_past_first_mate_get_distinct_chain_ids`
- `psp::reader::tests::chain_id_round_trips_across_block_boundary`
- `pileup_to_psp::tests::drive_pileup_to_psp_roundtrips_records_through_in_memory_sink`
  — the test deferred at the seam-implementation commit, now
  passing.

## Tradeoffs and follow-ups

- **File-size cost.** Plan estimated +10–20% larger `.psp` files.
  Switching `allele_chain_slots` to fixed LE `u64` (rather than
  varint) probably leaves us at the higher end of that range
  uncompressed; after zstd, the gap shrinks substantially because
  the high-order zero bytes compress trivially. One real-data
  benchmark on a representative input would confirm.
- **Pre-existing `pending_mates`-lifetime wart**, surfaced by the
  unique-id refactor: mate pairs whose second mate arrives *after*
  the first mate has expired from the active set get distinct
  chain ids today, even though they should share. The old
  recycling design hid this with an accidental id collision; the
  new design exposes the truth. Pinned by
  `paired_mates_with_gap_past_first_mate_get_distinct_chain_ids`.
  Fix is to tie the `pending_mates` entry lifetime to
  `mate_lookup_window`, not to active-set residence — out of scope
  for this refactor.
- **The `release_slot` / `release_pending_partner_ref_if_present`
  API surface is preserved** for the active-read-set caller's
  convenience, even though the semantics under unique ids are
  thinner ("decrement active_count by 1" / "clear pending_mates
  entry, count eviction"). A future cleanup could simplify the
  active-read-set side to call simpler primitives. Not blocking.

## File touch list

Code:

- `src/per_sample_caller/pileup/slot_allocator.rs` (rewrite)
- `src/per_sample_caller/pileup/mod.rs`
- `src/per_sample_caller/pileup/walker.rs`
- `src/per_sample_caller/pileup/open_record.rs`
- `src/per_sample_caller/pileup/active_read_set.rs`
- `src/per_sample_caller/pileup/errors.rs`
- `src/per_sample_caller/pileup/tests.rs`
- `src/per_sample_caller/psp/registry.rs`
- `src/per_sample_caller/psp/block.rs`
- `src/per_sample_caller/psp/writer.rs`
- `src/per_sample_caller/psp/reader.rs`
- `src/per_sample_caller/psp/errors.rs`
- `src/per_sample_caller/psp/mod.rs`
- `src/per_sample_caller/psp/header.rs`
- `src/per_sample_caller/pileup_to_psp.rs`
- `examples/dhat_psp_writer.rs`
- `benches/psp_writer_perf.rs`
- `benches/psp_reader_perf.rs`

Specs:

- `ia/specs/phase_chain.md`
- `ia/specs/per_sample_pileup_format.md`
