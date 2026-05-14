# Phase chains: switch from recycled `u16` slot ids to unique `u64` chain ids

Proposal date: 2026-05-14.

## Domain intent

Replace the current "recycled slot id + per-record lifecycle markers"
mechanism with a far simpler one: every read (or read-pair) gets a
**unique `u64` identifier** that is never reused within the lifetime
of a `.psp` file. The pre-existing identifier-recycling rationale
(keep `SlotId` small to save bytes) is being explicitly retired in
favour of correctness-by-construction and a substantial reduction in
cross-module complexity.

The downstream semantic — *"two allele observations sharing a chain
id come from the same physical molecule in this sample"* — is
unchanged. Stage 5's compound-haplotype check (chain-set
intersection) is also unchanged in shape, only simpler in
implementation (no active-set bookkeeping required).

## Why now

This plan is the resolution chosen for the integration finding pinned
in [src/per_sample_caller/pileup_to_psp.rs](../../src/per_sample_caller/pileup_to_psp.rs)
(`writer_rejects_walker_output_when_chain_expires_on_same_record_as_final_reference`)
and described in
[ia/reports/implementations/pileup_to_psp_seam_2026-05-14.md](../reports/implementations/pileup_to_psp_seam_2026-05-14.md).
That finding is one symptom of a broader category: the walker stamps
`expired_chains` on a record whose alleles still reference the
expiring slot, and the writer's "apply expirations before validating
alleles" rule rejects it.

The recycling mechanism that generates the conflict exists *only* to
keep ids small. Once we accept slightly larger per-record chain
bytes, the entire mechanism — and the bug class it enables — goes
away.

The previous rejection in
[ia/specs/per_sample_pileup_format.md](../specs/per_sample_pileup_format.md)
Q-FC8 / Q-PL4 of "never recycle slot ids across a file" was on the
basis that the `u16` namespace would overflow on real workloads.
With `u64`, that objection no longer applies.

## Out of scope

- Cohort-wide chain identifier reconciliation (no such thing exists;
  chains are per-sample). Different `.psp` files can reuse the same
  `u64` values internally and Stage 5 cares only about within-sample
  intersections.
- Re-tuning zstd compression level or block size. Existing settings
  are kept.
- File-format version bumps. The format has not shipped; this is
  pre-1.0 design churn and stays at `format_version = (1, 0)`.
- Performance optimization for the larger chain columns. Take the
  hit (~10–20% larger psp files in raw byte estimates; less after
  zstd) as the price of complexity removal.

## What changes, in one sentence per layer

- **Walker (pileup).** `SlotId` widens from `u16` to `u64`. The slot
  allocator becomes a monotonic counter — no free list, no
  refcount-driven release, no mark queues. Records no longer carry
  `new_chains` / `expired_chains`.
- **`.psp` block format.** Two columns (`new_chain_slots`,
  `expired_chain_slots`) are removed. `allele_chain_slots` switches
  from a list of `u16` to a list of `varint`. The block header's
  `active_chain_slots_at_block_start` snapshot is removed.
- **Writer.** All active-set tracking, marker validation, and
  per-block active-slot snapshotting is deleted (`active_slots`
  field, the `PhaseChainMarkerInconsistency*` validations, the
  `snapshot_active_slots` plumbing).
- **Reader.** Active-view bookkeeping during iteration is deleted;
  random-access region reads are simpler.
- **Stage 5 (future code).** Will not need an active-set view at
  all; the chain-intersection check becomes a direct `Vec<u64> ∩
  Vec<u64>` with no per-record state to maintain.
- **Specs.** `phase_chain.md`, `per_sample_pileup_format.md`, and
  `calling_pipeline_architecture.md` lose their recycling/lifecycle
  sections.

## Assumptions and silent choices

- **`SlotId = u64`** (not `u32`). Memory cost is 4 extra bytes per id
  in memory; disk cost is unchanged because varint encodes by
  magnitude, not type width. `u64` removes any future-proofing
  concern (covers up to ~10¹⁹ chains per file; any imaginable
  workload sits well under 10¹⁰).
- **Wire encoding for chain ids: `varint` (LEB128).** The format
  already uses LEB128 for other unsigned fields and the registry has
  `"varint"` as an `element-type`. Typical chain ids on real
  workloads will compress to 3–5 bytes per id; small ids (early in
  the file) take 1–2 bytes. Better than committing to a fixed `u64`
  on disk.
- **The chain id starts at `0` per file** and increments
  monotonically. The walker resets nothing across chromosomes — a
  chain id never repeats within one `.psp`.
- **`SlotExhausted` error variant is repurposed, not removed.** It
  currently fires when the active-read cap (`MAX_ACTIVE_SLOTS =
  4096`) is reached. That cap remains as a bound on
  *concurrent active reads*, not a cap on the id space. The variant
  is renamed to `ActiveReadsExhausted` (or similar — see "Open
  questions") with the same trigger condition.
- **`u64` chain id space is checked, not just assumed.** Although
  `2^64` chains per `.psp` file is astronomically beyond any
  imaginable workload (~10¹⁹ reads), the counter must still error
  cleanly on overflow rather than silently wrapping back to 0 and
  producing colliding ids that the reader cannot disambiguate. A
  hard check is cheap (one `checked_add` per allocation) and
  surfaces a clear failure rather than a silent corruption. See
  "Overflow guard" below.
- **Pair handling is unchanged.** The slot allocator's
  `pending_mates` map already returns the first mate's id when the
  second mate arrives; the new monotonic counter is consulted only
  for the first mate of a pair. One chain per molecule remains
  invariant.

## Overflow guard

The chain id counter is `u64`, so the absolute upper bound is `2^64 - 1`
≈ 1.8 × 10¹⁹ allocations per `.psp` file. That ceiling is unreachable
in practice — the worst plausible per-file workload sits ten orders
of magnitude below it — but silent wrap-around would be catastrophic:
two distinct molecules would carry the same chain id, and Stage 5's
intersection check would conflate them with no way to detect the
collision after the fact.

The allocator therefore performs a checked increment on every fresh
allocation:

```rust
let id = self.next_id;
self.next_id = self
    .next_id
    .checked_add(1)
    .ok_or(WalkerError::ChainIdSpaceExhausted { /* locus context */ })?;
```

`WalkerError::ChainIdSpaceExhausted` is a new variant in
`pileup/errors.rs`, carrying the locus where the overflow occurred:

```rust
#[error(
    "phase-chain id space exhausted at chrom_id={chrom_id} pos={pos}: \
     this psp file has reached 2^64 unique read identifiers, the per-file \
     limit imposed by the u64 chain id encoding"
)]
ChainIdSpaceExhausted { chrom_id: u32, pos: u32 },
```

The error is a hard fail — the walker stops and the partial `.psp` is
discarded, same shape as the other hard-fail variants. There is no
recovery path; the only mitigations are splitting the input into
smaller per-sample shards (which the pipeline already supports — one
`.psp` per BAM is the unit of work) or, if this ever became a real
constraint, widening to `u128`.

**Why we don't just `assert!` / `expect("u64 cannot overflow")`.** A
panic on this path would crash the walker with a backtrace and no
operator-actionable message. Returning a typed error keeps the failure
inside the `WalkerError` channel, surfaces a message that names the
limitation, and lets calling code log a sensible "this BAM is too
large for the current encoding; split it" hint.

**Test strategy for the overflow path.**

- Add a `#[cfg(test)]` constructor on the allocator —
  `SlotAllocator::with_next_id_for_testing(start: u64)` — that lets
  a test seed `next_id` to `u64::MAX`. Production code never sees this
  constructor (`#[cfg(test)] pub(crate)` keeps it invisible to the
  public API).
- New test
  `pileup/tests.rs::walker_errors_when_chain_id_counter_would_overflow`:
  seed the allocator to `u64::MAX`, run the walker against an input
  read, expect the iterator's first yielded item to be
  `Err(WalkerError::ChainIdSpaceExhausted { … })`.
- New test
  `pileup/tests.rs::walker_succeeds_at_chain_id_counter_one_below_max`:
  seed to `u64::MAX - 1`, run a single-read input, confirm the read
  is admitted (returns `u64::MAX` as its id) and the walker completes
  cleanly. Pins the off-by-one boundary so a future change to
  `checked_add`'s semantics doesn't silently relax the check.

## Module-by-module change list

### `src/per_sample_caller/pileup/slot_allocator.rs`

- `pub type SlotId = u16` → `pub type SlotId = u64`.
- Replace the free-list-based allocator with a struct holding
  `next_id: u64` (monotonic counter), `pending_mates: AHashMap<Arc<str>, …>`
  (kept; pair logic unchanged), plus counters
  (`slot_allocations`, `slot_high_water`, `mate_lookup_evictions`).
- Remove fields: `free: Vec<SlotId>`, `pending_free: Vec<SlotId>`,
  `refcount: AHashMap<SlotId, u8>`, `new_marks: Vec<SlotId>`,
  `expired_marks: Vec<SlotId>`, `max_active_slots`.
- Remove methods: `release_slot`, `drain_lifecycle_marks`,
  `set_refcount` and the refcount-aware path.
- `allocate_slot` becomes: if the qname is in `pending_mates`,
  return the stored id; else read `next_id` as the new id, then
  do `next_id = next_id.checked_add(1).ok_or(...)?` to bump it —
  returning `WalkerError::ChainIdSpaceExhausted { chrom_id, pos }`
  on overflow (see "Overflow guard" above).
- `slot_high_water` is tracked off the active-read-set's `len()`,
  not off an allocator-internal active set (which no longer exists).
  Sourced from the walker's `active_reads` instead.
- `SlotAllocatorCounters` simplifies — `slot_allocations` is still
  meaningful (number of fresh `next_id` increments), `slot_high_water`
  moves to a walker-level computation, `mate_lookup_evictions` stays.
- Add `#[cfg(test)] pub(crate) fn with_next_id_for_testing(start: u64) -> Self`
  to let the overflow-guard tests seed `next_id` near `u64::MAX`.

### `src/per_sample_caller/pileup/active_read_set.rs`

- `expire_passed` and `flush_all` no longer call `release_slot`;
  they just drop the `ActiveRead` entries. Slot ids the dropped
  reads carried are simply unused going forward (uniqueness ensures
  no aliasing).

### `src/per_sample_caller/pileup/mod.rs`

- `PileupRecord` loses two fields: `new_chains: Vec<SlotId>` and
  `expired_chains: Vec<SlotId>`. The struct's doc-comment loses the
  whole §"Lifecycle markers" block.
- `PileupRecord::new` constructor signature changes (no
  `new_chains` / `expired_chains` arguments). External callers
  (examples, benches, doc tests) are updated.

### `src/per_sample_caller/pileup/open_record.rs`

- `OpenAllele.chain_slots: Vec<SlotId>` is unchanged in shape but
  now holds `u64`s.
- `finalise()` unchanged (no longer copies lifecycle marks because
  they don't exist).

### `src/per_sample_caller/pileup/walker.rs`

- `close_aged_records_into`: drop the `drain_lifecycle_marks` call
  and the `out.get_mut(batch_start)` stamp.
- `flush_chromosome_into`: same — drop the
  `drain_lifecycle_marks` + stamp pair.
- The `expire_passed_reads` / `close_aged_records_into` / `advance`
  step ordering inside `fill_pending` is preserved (no longer
  needed for marker semantics, but the existing ordering still
  matches the closure rule).

### `src/per_sample_caller/pileup/errors.rs`

- Rename `WalkerError::SlotExhausted` → `WalkerError::ActiveReadsExhausted`
  (mechanical; same trigger and meaning, but the new name is
  honest about what's exhausted now). The cap is still
  `DEFAULT_MAX_ACTIVE_SLOTS` / `WalkerConfig::max_active_slots`;
  consider renaming those too for consistency.
- Add `WalkerError::ChainIdSpaceExhausted { chrom_id: u32, pos: u32 }`
  for the `u64` counter overflow path. See "Overflow guard" above
  for the message text and rationale.

### `src/per_sample_caller/psp/writer.rs`

- Remove the `IngestState::active_slots: Vec<SlotId>` field and all
  uses.
- Remove `apply_record_to_block`'s entire marker-validation and
  marker-application block (lines ~528–560).
- Remove the per-chrom-boundary "active_slots must be empty" check.
- Remove `BlockAccumulator::snapshot_active_slots` and the cloning
  of `active_slots` at block open (writer.rs ~308–325).
- Remove the `new_chain_slots` and `expired_chain_slots` list
  columns from `BlockAccumulator` and `append_record`.
- `BlockHeader` emission no longer carries
  `active_chain_slots_at_block_start`.

### `src/per_sample_caller/psp/reader.rs`

- Remove the active-slot running view from
  `RecordsIter`; drop the snapshot-init logic at block start.
- Sequential and random-access paths converge — both just decode
  the records without bootstrapping any state.
- Remove block-boundary continuity checks tied to the snapshot
  (spec checks #10–11 in `per_sample_pileup_format.md`).

### `src/per_sample_caller/psp/block.rs`

- Drop the `new_chain_slots` and `expired_chain_slots` `ListColumn`
  encode/decode paths.
- Switch `allele_chain_slots` from u16-payload to varint-payload.
  (The `ListColumn` infrastructure already supports varint payloads
  for other columns; this is reusing existing machinery.)
- Update the `BlockHeader` codec to drop the
  `active_chain_slots_at_block_start` field.

### `src/per_sample_caller/psp/registry.rs`

- Drop the `(0x20, "new-chain-slots")` and
  `(0x21, "expired-chain-slots")` entries from `V1_0_COLUMNS`.
- Update `(0x22, "allele-chain-slots")`'s `element-type` to
  `"varint"`.

### `src/per_sample_caller/psp/errors.rs`

- Drop `PspWriteError::PhaseChainMarkerInconsistency` and the
  `PhaseChainMarkerInconsistencyKind` enum entirely.
- Drop `PspReadError::PhaseChainConsistency` (block-level
  consistency check #10) and `PhaseChainConsistencyKind`.
- Drop the inter-block continuity error variant (check #11).

### `src/per_sample_caller/pileup_to_psp.rs`

- The regression test
  `writer_rejects_walker_output_when_chain_expires_on_same_record_as_final_reference`
  is rewritten into a real roundtrip:
  - Drive a walker → writer → bytes → reader pipeline.
  - Assert per-record equality against an independent walker run.
  - Assert that every emitted `PileupRecord` has
    `alleles[*].chain_slots` containing monotonically-allocated unique
    `u64` ids.
- The `drive_pileup_to_psp_passes_records_through_unmodified` test is
  kept as-is.

### Examples and benches

- `examples/dhat_pileup.rs`, `examples/dhat_psp_writer.rs`: any
  `PileupRecord::new(..., new_chains, expired_chains, ...)` callers
  drop the two arguments.
- `benches/pileup_walker_scaling.rs`, `benches/psp_writer_perf.rs`,
  `benches/psp_reader_perf.rs`: same. Bench fixtures simplify.

## Spec / documentation changes

- `ia/specs/phase_chain.md`:
  - §6 "The slot encoding, with a worked example" rewritten — the
    encoding becomes "each read gets a unique `u64`; downstream
    intersects sets". The worked example loses `new_chains` /
    `expired_chains` lines from every per-position record.
  - §8 "What Stage 5 does with the chains" simplifies: no active
    view to maintain.
- `ia/specs/per_sample_pileup_format.md`:
  - §"Block header" loses `active_chain_slots_at_block_start`.
  - §"Column registry" loses `new_chain_slots` and
    `expired_chain_slots`; `allele_chain_slots` element-type
    becomes `varint`.
  - Consistency checks #10 (phase-chain active-set within block)
    and #11 (inter-block phase-chain continuity) are removed.
  - Q-FC8 marked superseded; Q-PL4 marked superseded with a pointer
    to this plan.
- `ia/specs/calling_pipeline_architecture.md`:
  - §"Phase chain identifiers" subsection (under Stage 1) updated
    to reflect the unique-id design.
  - Any Stage 5 references to active-set bookkeeping simplified.

## Test strategy

### Tests removed (no longer applicable)

- `pileup/tests.rs::lifecycle_markers_appear_on_emitted_records` and
  any sibling tests that assert on `new_chains` / `expired_chains`
  totals.
- `psp/writer.rs::tests` cases that assert on
  `PhaseChainMarkerInconsistency` variants.
- `psp/reader.rs::tests` cases that assert on the active-slot
  snapshot or on consistency checks #10–11.
- The regression-pin test in `pileup_to_psp.rs` (replaced by a real
  roundtrip — see above).

### Tests added

- `pileup/tests.rs::chain_ids_are_unique_and_monotonically_allocated`
  — assert every chain id observed across emitted records is
  unique, and that ids appear in non-decreasing order of first
  reference position.
- `pileup/tests.rs::paired_mates_share_a_single_chain_id` — pin
  the pair-collapsing semantic at the new boundary.
- `pileup/tests.rs::walker_errors_when_chain_id_counter_would_overflow`
  — seed the allocator to `u64::MAX` via the test-only constructor,
  run the walker against an input read, expect the first yielded
  item to be `Err(WalkerError::ChainIdSpaceExhausted { … })`.
- `pileup/tests.rs::walker_succeeds_at_chain_id_counter_one_below_max`
  — seed to `u64::MAX - 1`, single-read input, confirm the read is
  admitted with id `u64::MAX` and the walker completes cleanly.
  Pins the off-by-one boundary so a future change to the
  `checked_add` semantics doesn't silently relax the check.
- `pileup_to_psp.rs::drive_pileup_to_psp_roundtrips_records_through_in_memory_sink`
  — the test deferred by the previous plan; now achievable.
- `psp/writer.rs::tests` updated: hand-built records no longer carry
  lifecycle marks; existing record-validation tests stay.
- `psp/reader.rs::tests` updated: drop snapshot/continuity checks;
  add a direct sequential-decode test that confirms
  `allele.chain_slots` round-trips losslessly with varint
  encoding.

### Tests preserved as-is

- All non-chain pileup walker tests (SNP folding, deletion
  anchoring, REF widening, mate-overlap, etc.).
- All non-chain psp writer / reader tests (header, trailer, block
  index, allele sequence encoding, scalar columns, region queries).
- The integration tests in `tests/`.

## Migration sequencing

A single cohesive commit makes most sense — the format is changing
in a way that doesn't have a meaningful intermediate state. Suggested
ordering within the commit:

1. **Spec first.** Update `phase_chain.md`,
   `per_sample_pileup_format.md`, and `calling_pipeline_architecture.md`
   to reflect the new design. Reviewers can validate intent before
   reading code.
2. **`SlotId` widening + walker.** `pileup/mod.rs` (type), `slot_allocator.rs`
   (counter), `walker.rs` (drop stamp logic), `errors.rs` (rename).
3. **PileupRecord shape.** Drop the two fields; update
   `PileupRecord::new` callers (examples, benches, doc tests).
4. **psp block format.** `registry.rs`, `block.rs`, `header.rs`
   (codec), `writer.rs` (drop bookkeeping), `reader.rs` (drop
   bookkeeping), `errors.rs` (drop variants).
5. **Tests.** Remove deprecated tests, add the new ones, flip the
   regression-pin to a real roundtrip.
6. **`cargo fmt` / `cargo clippy` / `cargo test` inside the
   container.**

## Validation

Inside the dev container:

- `cargo fmt --check` clean.
- `cargo clippy --all-targets --all-features -- -D warnings` clean.
- `cargo test --tests --lib --all-features` all green, including
  the resurrected roundtrip parity test.
- `cargo check --benches --all-features` clean.

A small one-off sanity script (or just an extra logging line in the
roundtrip test, gated behind a custom feature) can compare the
post-change `.psp` byte size against a saved pre-change baseline on
the same input — useful for confirming the "~10-20% larger" estimate
empirically, but not blocking.

## Risks

- **Active-read cap interaction.** The current `MAX_ACTIVE_SLOTS =
  4096` cap is a memory bound on the active read set, not on the
  slot id space. The rename clarifies this but does not change the
  cap. Pathologically deep regions (depth > 4096) still error; that's
  unchanged.
- **Mate-pair correctness.** The `pending_mates` mechanism is the
  only place where the new monotonic-counter allocator behaves
  non-trivially. The new `paired_mates_share_a_single_chain_id` test
  pins it; manual review of the allocator's mate-merge path is
  warranted before commit.
- **Reader breakage on stale files.** Any `.psp` files generated by
  pre-change writers can no longer be read after this change. No
  such files exist in tree (verified by `cargo test`); production
  hasn't shipped. If anyone has stale local files, they will need to
  regenerate them — the format version stays at `(1, 0)` but the
  block schema differs in a way the new reader won't recognise.
- **File-size estimate uncertainty.** "10–20% larger" is a back-of-
  envelope figure. Real measurement on a 30× WGS sample (or a
  representative slice) could come in lower (zstd is good at
  redundancy) or higher (if chain density turns out higher than
  expected). Not blocking, but worth a one-shot measurement before
  committing the migration.

## Out-of-scope follow-ups

- A `chain_id_first_seen_chrom_id` column or similar metadata if a
  future stage wants to detect chains spanning unexpected
  chromosomes. Skip until needed.
- Eviction strategy if the active-read set ever needs to grow
  beyond `MAX_ACTIVE_SLOTS`. The current behaviour (error out) is
  fine for now; tuning happens when real data demands it.
- Anywhere downstream that wants `chain_slots` to imply ordering
  *across* records (e.g., "chain X was allocated before chain Y" →
  X < Y always). The monotonic counter gives this for free, but no
  consumer relies on it yet; not a documented contract.

## File touch list

Code:

- `src/per_sample_caller/pileup/slot_allocator.rs` (heavy)
- `src/per_sample_caller/pileup/active_read_set.rs` (light)
- `src/per_sample_caller/pileup/mod.rs` (medium)
- `src/per_sample_caller/pileup/open_record.rs` (light)
- `src/per_sample_caller/pileup/walker.rs` (light)
- `src/per_sample_caller/pileup/errors.rs` (rename + cleanup)
- `src/per_sample_caller/pileup/tests.rs` (medium — remove/add tests)
- `src/per_sample_caller/psp/writer.rs` (heavy)
- `src/per_sample_caller/psp/reader.rs` (medium)
- `src/per_sample_caller/psp/block.rs` (medium)
- `src/per_sample_caller/psp/header.rs` (light — drop snapshot field)
- `src/per_sample_caller/psp/registry.rs` (light)
- `src/per_sample_caller/psp/errors.rs` (cleanup)
- `src/per_sample_caller/pileup_to_psp.rs` (flip regression test to real roundtrip)
- `examples/dhat_pileup.rs`, `examples/dhat_psp_writer.rs`
- `benches/pileup_walker_scaling.rs`, `benches/psp_writer_perf.rs`,
  `benches/psp_reader_perf.rs`

Specs:

- `ia/specs/phase_chain.md`
- `ia/specs/per_sample_pileup_format.md`
- `ia/specs/calling_pipeline_architecture.md`
