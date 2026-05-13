# `.psp` reader — `PspReader` implementation

Implementation plan for the **reader half** of the per-sample pileup
format. Covers the missing `src/per_sample_caller/psp/reader.rs`
file: opening a `.psp` from a `Read + Seek` source, validating it,
and yielding `PileupRecord`s.

This is a focused delta on top of the original combined plan
[per_sample_pileup_writer_reader.md](per_sample_pileup_writer_reader.md),
which already specifies the reader's public API, validation order,
and round-trip tests. That plan is still the byte-format reference;
this one only documents the work remaining and the orchestration
decisions that the original plan deferred to the implementer.

Out of scope:

- The multi-sample cohort gather/merge that joins N `PspReader`s
  into a per-position cohort stream. Owns its own plan; the Stage 3
  DUST filter and Stage 4 grouping consume that, not the reader
  directly.
- `psp head` / `psp dump` CLI utilities.
- Parallel column decompression, mmap, streaming zstd — all deferred
  in the original plan.
- Writer-side work — already shipped.

## Cross-references

- [per_sample_pileup_format.md](../specs/per_sample_pileup_format.md) —
  byte layout; in particular §"Header-binary consistency: required
  reader checks" (the 11-item checklist this plan maps to reader.rs)
  and §"Block layout".
- [per_sample_pileup_writer_reader.md](per_sample_pileup_writer_reader.md) —
  the combined writer/reader plan. §"Reader public API" defines the
  signatures; §"Reader-side validation order" defines the order of
  operations; §"Block layout — read path" describes the per-block
  decode. This plan refines and extends.
- [calling_pipeline_architecture.md §"Stage 3"](../specs/calling_pipeline_architecture.md) —
  consumer context (DUST filter sits between `.psp` reading and
  Stage 4 grouping).
- [design_principles.md](../specs/design_principles.md) §3 — "Errors
  must not pass silently". Every check in the spec's checklist is
  load-bearing.

## What is already in place

The decoder primitives the reader will compose are all implemented
and unit-tested. The job of `reader.rs` is to **orchestrate** them,
not to add new wire-level codecs.

| File | What it already does |
|---|---|
| [psp/varint.rs](../../src/per_sample_caller/psp/varint.rs) | LEB128 + zig-zag svarint encode/decode, edge-case tests. |
| [psp/trailer.rs](../../src/per_sample_caller/psp/trailer.rs) | `decode_trailer(&[u8; 32]) -> Result<Trailer, PspReadError>`. Magic + arithmetic checks inside. |
| [psp/index.rs](../../src/per_sample_caller/psp/index.rs) | `decode_index(bytes, expected_n_blocks) -> Result<Vec<BlockIndexEntry>, PspReadError>`. Checksum verification + monotonicity inside. |
| [psp/header.rs](../../src/per_sample_caller/psp/header.rs) | `parse_header_bytes(&[u8]) -> Result<(ParsedHeader, usize), PspReadError>` runs magic, length-prefix range, TOML parse, per-field rules, schema-vs-registry agreement (checks 1–6 of the spec checklist), and sentinel cross-check. Also `parse_header_toml(body)` for the body-only path. The `ParsedHeader` / `ParsedChromosome` / `ParsedColumn` types are public and final. |
| [psp/block.rs](../../src/per_sample_caller/psp/block.rs) | `decode_block_header(&[u8]) -> Result<(BlockHeader, usize), …>`, `decode_scalar_column`, `decode_varint_column`, `decode_list_column`, `decode_bytes_split`, `zstd_decompress`. Block-header invariants (`n_records ≥ 1`, manifest tag-ascending, snapshot ascending, etc.) are checked inside `decode_block_header`. |
| [psp/registry.rs](../../src/per_sample_caller/psp/registry.rs) | `V1_0_COLUMNS` source of truth. |
| [psp/errors.rs](../../src/per_sample_caller/psp/errors.rs) | `PspReadError` with every variant the spec checklist names. The reader does not need to add variants. |
| [psp/writer.rs](../../src/per_sample_caller/psp/writer.rs) | The `PspWriter` half — used by round-trip tests. |

What is missing:

- [psp/reader.rs](../../src/per_sample_caller/psp/reader.rs) — does not
  yet exist. This plan describes its contents.
- The `pub mod reader;` line and re-exports in
  [psp/mod.rs](../../src/per_sample_caller/psp/mod.rs).
- The reader-side tests in Groups B (round-trip) and F
  (corruption) listed in the original plan — they become executable
  once `reader.rs` lands.

## Module additions

`src/per_sample_caller/psp/mod.rs`:

- Add `pub mod reader;`.
- Add `pub use reader::{PspReader, RecordsIter};` (the parsed-header
  types are already re-exported from `header`).
- Drop the "The reader-side public API (`PspReader`) is planned but
  not yet implemented" sentence in the module-level doc.

`src/per_sample_caller/psp/reader.rs` — new file. Sections below.

## Public surface (recap, with deltas)

The original plan's signatures stand; one delta from what the
writer actually built:

```rust
// src/per_sample_caller/psp/reader.rs
use std::io::{Read, Seek};

use super::header::ParsedHeader;
use super::index::BlockIndexEntry;
use super::errors::PspReadError;
use crate::per_sample_caller::pileup::PileupRecord;

pub struct PspReader<R: Read + Seek> { /* private */ }

impl<R: Read + Seek> PspReader<R> {
    pub fn new(source: R) -> Result<Self, PspReadError>;
    pub fn header(&self) -> &ParsedHeader;
    pub fn block_index(&self) -> &[BlockIndexEntry];
    pub fn records(&mut self) -> RecordsIter<'_, R>;
    pub fn region_records(
        &mut self,
        chrom_id: u32,
        start: u32,
        end: u32,
    ) -> RecordsIter<'_, R>;
}

pub struct RecordsIter<'r, R: Read + Seek> { /* private */ }

impl<'r, R: Read + Seek> Iterator for RecordsIter<'r, R> {
    type Item = Result<PileupRecord, PspReadError>;
    fn next(&mut self) -> Option<Self::Item>;
}
```

Deltas from the original plan:

- The `ParsedHeader` type lives in `header.rs` and uses
  `toml::value::Datetime` for the `created` field, not
  `chrono::DateTime`. The original plan suggested chrono; the writer
  ended up using `toml::Datetime` to keep the dependency surface
  small. Reader inherits.
- Add `block_index()` accessor. Useful for tests, `psp dump`, and
  the eventual cohort merge driver that wants to know per-sample
  block layout before scheduling reads. Free to expose; the index
  is already in memory.

`RecordsIter` is **not** `Send`. It borrows the source mutably, so
only one iterator lives at a time per `PspReader`. The cohort
merge driver runs one `PspReader` per sample and pulls from each
sequentially; per-sample threading is the merge driver's
concern, not this module's.

## Open sequence (`PspReader::new`)

Order chosen so each step gates the next on data that is already
known to be well-formed. Every step is a hard fail on any error.

1. **Seek to end; size the file.** `source.seek(SeekFrom::End(0))`
   → total length `L`. If `L < HEADER_FRAMING_BYTES + TRAILER_BYTES`
   (i.e. cannot contain even an empty file's framing), reject with
   `PspReadError::FileTooShort` or `Io { context: "trailer" }` after
   the trailer read fails. Single bound check up front avoids two
   separate fails-after-IO error paths.
2. **Read trailer.** `source.seek(SeekFrom::End(-32))`, read 32 bytes,
   pass to `trailer::decode_trailer`. This validates `PSPE` magic and
   the arithmetic `index_offset + index_byte_length ≤ L - 32`.
3. **Read block index.** `source.seek(SeekFrom::Start(trailer.index_offset))`,
   read `trailer.index_byte_length` bytes, pass to
   `index::decode_index(bytes, trailer.n_blocks)`. This verifies the
   XXH3-64-trunc32 checksum, decodes exactly `n_blocks` entries, and
   enforces coordinate monotonicity (chrom_id non-decreasing;
   first_pos strictly increasing within a chromosome).
4. **Sanity-check block offsets.** Walk the decoded entries once: every
   `block_offset` must be `≥ HEADER_FRAMING_BYTES + 1` (after the
   header), strictly increasing across entries, and the last entry's
   `block_offset` must be `< trailer.index_offset` (blocks end before
   the index starts). Mismatch → `PspReadError::BlockIndexOffsetInvalid`
   (variant name follows existing convention; add if missing — see
   "Error variants" below).
5. **Read and parse header.** `source.seek(SeekFrom::Start(0))`,
   read the first 12 bytes (magic + length prefix), validate magic,
   range-check `toml_body_length`, allocate the body buffer
   (capped at `MAX_HEADER_BODY_BYTES`), `read_exact` body + sentinel,
   hand the whole slice (12 + body + 17 bytes) to
   `header::parse_header_bytes`. That call runs checks 1–6 of the
   spec checklist plus the sentinel cross-check.
6. **Bind chromosomes to index entries.** Walk the block index a
   second time: every `chrom_id` must be `< header.chromosomes.len()`,
   every `first_pos` and `last_pos` must lie in
   `[1, chromosomes[chrom_id].length]`. Mismatch →
   `PspReadError::BlockIndexChromOutOfRange` or
   `BlockIndexPosOutOfRange`.
7. **Stash state.** Construct `PspReader { source, header, index,
   ... }`. Position the source cursor at end-of-header (offset 12 +
   body_len + 17 = `header_end_offset`, returned by
   `parse_header_bytes`'s second return value). This positions the
   cursor for an immediate sequential `records()` call without
   needing another seek.

Tail-first is correct: an empty file (no blocks) still has a
well-formed header and trailer; the trailer gives us `n_blocks = 0`,
the index is zero bytes, the open path returns successfully, and
`records()` immediately yields `None`.

## Sequential iteration (`records()`)

`RecordsIter` is a state machine driven by `Iterator::next`. State:

```rust
struct RecordsIter<'r, R: Read + Seek> {
    reader: &'r mut PspReader<R>,
    // Which entry of the block index we are currently inside; None
    // before the first block has been loaded, Some(i) while
    // iterating block i, == n_blocks when exhausted.
    cur_block: BlockCursor,
    // Decoded payload of the current block (None if cur_block has
    // not yet been pulled from disk). Owning all per-column buffers
    // for the block; freed when the block is exhausted.
    cur_payload: Option<DecodedBlock>,
    // Index of the next record to materialise inside cur_payload.
    next_record_in_block: u32,
    // Running active phase-chain set, carried across blocks for the
    // sequential-mode continuity check (spec check #11). BTreeSet
    // is fine: peak size is bounded by max active slots (the writer
    // caps this; well under 1k in practice).
    active_chain_slots: BTreeSet<SlotId>,
    // Region clamp; set by `records()` to "no clamp", by
    // `region_records()` to a chrom+window. See next section.
    clamp: Clamp,
    // Sticky error flag: once `next()` has yielded an Err, future
    // calls must return None. Bookkeeping for the Iterator contract.
    poisoned: bool,
}
```

`next()`:

1. If `poisoned`, return `None`.
2. If `cur_payload.is_none()`:
   - Advance `cur_block` to the next block index entry consistent
     with `clamp`. If none remains, return `None`.
   - `source.seek(SeekFrom::Start(entry.block_offset))`.
   - Decode the block header: read enough bytes to feed
     `block::decode_block_header`. **Bytes-to-read is unknown
     up-front** because the header is a varint stream. Strategy:
     read in a small initial chunk (start with 256 bytes, doubling
     up to a hard cap of e.g. 64 KiB), call `decode_block_header`;
     if it returns `IncompleteVarint`, read more and retry. Helper
     `read_varint_framed_section` lives in `reader.rs`.
   - Validate snapshot vs. running set (spec check #11) — see
     "Phase-chain bookkeeping" below.
   - Walk the column manifest, decoding columns in tag-ascending
     order. For each column:
       - Apply the a-priori `uncompressed_len` prediction where
         the schema allows (spec check #7 / case 1: fixed-width
         scalars). For `bytes` columns, defer the prediction until
         after the `length-column` payload has been decompressed
         (case 2).
       - `read_exact(manifest.compressed_len)`, `zstd_decompress`,
         confirm the result length equals `manifest.uncompressed_len`.
       - Dispatch on `(cardinality, shape, element_type)` to the
         matching `decode_*` helper in `block.rs`. Store the result
         in a column-typed slot inside `DecodedBlock`.
   - After all columns decoded, run finite-float check on every f32
     and f64 column (spec check #9).
   - Set `next_record_in_block = 0`.
3. Materialise one record:
   - `pos = last_pos + delta_pos[next_record_in_block]` (where
     `last_pos = first_pos` at index 0, else the previous record's
     pos).
   - Read `new_chain_slots[i]`, `expired_chain_slots[i]` for this
     record. Apply to `active_chain_slots`: `expired` then `new`,
     in that order (mirrors the writer; the spec says the running
     set after record `i` is `(prev ∪ new_i) \ expired_i`).
   - Validate the phase-chain consistency: `expired_i ⊆ prev_active`,
     `new_i ∩ prev_active = ∅`, every slot in any allele's
     `chain_slots` appears in the post-application active set (spec
     check #10).
   - Slice into per-allele column buffers using the running allele
     offset (cumulative sum over `n_alleles`).
   - Construct one `PileupRecord` + its `AlleleObservation`s.
   - If `clamp` is set and the record falls outside the window,
     decide: skip (continue to next record) if before window;
     return `None` if past window; otherwise yield.
   - Increment `next_record_in_block`. If it equals
     `cur_payload.n_records`, drop `cur_payload` (free memory) and
     leave `cur_block` pointing at the just-finished block. The
     next `next()` call advances `cur_block`.
4. On any error, set `poisoned = true`, return `Some(Err(_))`.

## Region iteration (`region_records`)

`region_records(chrom_id, start, end)` returns a `RecordsIter` with:

- `clamp = Clamp::Window { chrom_id, start, end }`.
- `cur_block` positioned at the first index entry whose
  `chrom_id == chrom_id && last_pos >= start && first_pos <= end`.
  Use binary search on `block_index` (`chrom_id`, then `first_pos`).
- `active_chain_slots` initialised from the entry block's
  `active_chain_slots_at_block_start` snapshot — **not** from a
  running set carried from earlier blocks (the spec requires this
  for random-access reads; check #11 is not performed).
- Iteration stops when the next block's `first_pos > end` or
  `chrom_id != requested chrom_id`.

Inside each in-range block, records before `start` are silently
skipped (they cost a column slice but no allocation beyond the
buffers already in `DecodedBlock`); records strictly after `end`
terminate the iterator.

Empty regions and regions wholly outside any block both return
zero records with no error.

## Phase-chain bookkeeping (reader side)

The reader maintains `active_chain_slots: BTreeSet<SlotId>` as the
running active set. The slot id type is
[`SlotId`](../../src/per_sample_caller/pileup/slot_allocator.rs#L21)
= `u16` (re-exported via `pileup::SlotId`).

Per spec check #10, at each record:

- The pre-application active set must contain every slot in
  `expired_chain_slots` for this record.
- The pre-application active set must **not** contain any slot in
  `new_chain_slots` for this record.
- After applying `expired` then `new`, every slot id mentioned in
  any allele's `chain_slots` must be a member of the active set.

On a block boundary in sequential mode (check #11): the snapshot in
the just-decoded block header must equal the running set carried
forward from the previous block. The first block of a chromosome
(including block 0 of the file) has an empty snapshot, and the
running set must also be empty at that point — the writer drains
all active slots before crossing a chromosome boundary, and the
reader cross-checks.

In random-access mode the running set is reset from the snapshot at
every block; the inter-block check is structurally impossible.

## Error variants

`PspReadError` already covers every check the spec lists. The reader
may need a couple of new variants for the block-offset cross-checks
in steps 4 and 6 of the open sequence:

- `BlockIndexOffsetInvalid { block, offset, expected_range }` — for
  out-of-bounds or non-monotonic `block_offset` values discovered
  after the index has decoded.
- `BlockIndexChromOutOfRange { block, chrom_id, n_chroms }`.
- `BlockIndexPosOutOfRange { block, chrom_id, pos, chromosome_length }`.

Each new variant gets a one-line addition to `errors.rs` with a
spec-aligned message. The negative tests in Group F that need them
ship in the same commit.

All I/O failures wrap into `PspReadError::Io { context, source }`.
Distinct `context` strings per I/O site (`"trailer"`, `"block index"`,
`"header magic"`, `"header body"`, `"block header"`, `"block payload
column {name}"`) so a corrupt-file user can tell which section
failed without running a debugger.

## Mapping of spec checks to enforcement site

| Spec check | Where it fires |
|---|---|
| 1. Header well-formedness | `header::parse_header_bytes` (already enforced) |
| 2. Schema completeness | `header::parse_header_bytes` |
| 3. `length-column` references | `header::parse_header_bytes` |
| 4. Tag uniqueness/range | `header::parse_header_bytes` |
| 5. Required-column recognition | `header::parse_header_bytes` |
| 6. Schema-vs-registry agreement | `header::parse_header_bytes` |
| 7. Per-block manifest agreement, fixed-width prediction | `reader.rs` block-decode loop, before each column's zstd read |
| 7. Per-block manifest agreement, bytes prediction (post-length-column) | same loop, deferred until length-column decoded |
| 7. Decompressed size matches manifest | same loop, after each `zstd_decompress` |
| 7. Structural correctness (list/varint/bytes overrun-or-underrun) | inside the `decode_*` helpers (already enforced); reader checks return value |
| 8. No surprises in skipped columns | `reader.rs` walks manifest entries; skipping any column would consume exactly `compressed_len` bytes anyway. v1.0 has no optional columns; this is a structural invariant. |
| 9. No non-finite floats | `reader.rs` per-column post-decompress sweep over f32/f64 columns |
| 10. Phase-chain active-set consistency (per record) | `reader.rs` record-materialise loop |
| 11. Inter-block continuity (sequential only) | `reader.rs` block-boundary handler |

The decoder primitives already enforce most invariants that can be
checked locally. The reader.rs orchestrator is responsible for
checks that span multiple primitives (snapshot vs. running set;
manifest vs. block header; finite-float sweeps; cross-block
continuity).

## Tests

The original plan lists Group B (B1–B7) round-trip and Group F
(F1–F8) corruption tests. Those land as part of `reader.rs` and
are unblocked the moment the reader is wired up. This section
lists only **reader-specific additions** — assertions the round-trip
tests don't cover.

### Group R — reader-only behaviour

In `psp/tests.rs` (sibling to the writer's tests; created the same
slice the reader lands in).

#### R1 — empty file round-trip

Write zero records; finish. `PspReader::new` succeeds; `header()`
matches; `block_index().is_empty()`; `records()` yields `None`
immediately; `region_records(0, 1, 100_000)` also yields `None`.

#### R2 — `header()` accessible before iteration

Open a real file; call `header()`; do not call `records()`; drop.
Reading the header alone does not require any block reads.

#### R3 — header round-trip equality

Build a `WriterHeader` with non-default values for every field
(sample name with full 64-char width, two chromosomes with random
MD5s, three input CRAMs, four parameters of mixed types). Round-trip
via writer → bytes → reader; assert the resulting `ParsedHeader`'s
public fields equal the input. This is the strongest single
guarantee the header layer can offer.

#### R4 — `block_index()` matches writer emission

Round-trip a file with 5 blocks deliberately sized (target block
size = 16 KiB, 200 small records). Compare `block_index()` entry by
entry against what the writer emitted via its public side (compare
counts and ranges; offset equality is implicit in monotonicity).

#### R5 — region across multiple blocks

3 blocks on chrom 0 covering positions 1–100, 101–200, 201–300.
`region_records(0, 50, 250)` returns records with `50 ≤ pos ≤ 250`,
in order, across all three blocks.

#### R6 — region wholly outside coverage

Same fixture. `region_records(0, 1_000_000, 2_000_000)` yields no
records, no error.

#### R7 — iterator dropped mid-stream; new iterator restarts

Open a file with 100 records. Take a `records()` iterator; consume
20; drop. Take a fresh `records()` iterator; assert it yields all
100 from the top. This documents the reader's "iterators do not
share progress with the source" contract.

#### R8 — sequential and region paths produce same records inside a region

Iterate the whole file; collect records inside `[start, end]` on
`chrom_id`. Independently call `region_records(chrom_id, start,
end)`. Assert equality. Cross-validates the two paths.

#### R9 — phase-chain continuity across blocks (positive)

A fixture where slot 7 spans the boundary between block 0 and
block 1. Sequential `records()` succeeds. Bonus: collect the
`PileupRecord`s and assert slot 7 appears in both blocks' relevant
allele `chain_slots`.

#### R10 — phase-chain snapshot vs. running set mismatch

Take the R9 fixture; mutate block 1's
`active_chain_slots_at_block_start` snapshot byte to drop slot 7.
Sequential `records()` returns a
`PhaseChainConsistency { block: 1, … }` error. `region_records(0,
…)` *starting at block 1* succeeds (trusts the snapshot) — different
contract, different test.

### Groups B and F (already specified)

B1–B7 and F1–F8 in the original plan execute end-to-end after
`reader.rs` exists. No changes; they become live tests in the same
commit.

## Implementation order

Suggested slicing (each step lands green):

1. **`reader.rs` skeleton + open path.** Add the file with
   `PspReader::new`, `header()`, `block_index()`, and stubbed
   `records()` / `region_records()` returning an empty iterator.
   Implements steps 1–7 of the open sequence. Test: R1, R2, R3, R4
   pass; round-trip header equality (the strongest invariant
   testable without block decode).

2. **Block decode loop.** Implement single-block forward iteration:
   pull one block, decode its header, decode each column, materialise
   records. No region clamp, no phase-chain check yet. Test: B1, B2.

3. **Multi-block sequential iteration.** Add the `cur_block`
   advance, block-boundary handling, finite-float sweep. Test: B3
   (multi-block) and B4 (multi-chrom).

4. **Phase-chain bookkeeping.** Add the running active set, per-record
   checks, cross-block continuity. Test: B5 (positive), R9, R10.

5. **`region_records`.** Add region clamp, block-index binary
   search, snapshot-only active set init. Test: B7, R5, R6, R8.

6. **Negative tests.** F1–F8: forge or corrupt fixtures, assert
   error variants. Adds the new `BlockIndex*` error variants
   identified in §"Error variants" as their tests demand them.

7. **Mod export + doc.** `pub mod reader;`, re-export `PspReader`
   and `RecordsIter`; update mod.rs's status sentence.

Each step ends with `cargo fmt --check`, `cargo clippy --all-targets
-- -D warnings`, and `cargo test per_sample_caller::psp` (targeted)
all green. Run inside the container per CLAUDE.md.

## Tradeoffs and follow-ups

- **Block-at-a-time decoding, no lookahead.** Simpler state machine
  and bounded memory (one block's worth of column buffers). Cost: a
  small latency spike at each block boundary while the next block
  decompresses. Acceptable; the cohort merge driver will overlap
  decode with downstream work anyway. Streaming or lookahead is a
  follow-up if profiling demands it.

- **Index loaded eagerly into memory.** Spec §"Size estimate" caps
  this at ~100 KB on a 5× WGS. Trivial. Revisit only if a future
  use case produces millions of blocks per file.

- **`RecordsIter` is not `Send`.** It borrows `&mut R`. The cohort
  merge driver runs one reader per sample; if it wants to decode
  multiple samples on different threads it needs one reader per
  thread, which is the natural shape anyway (one `Read + Seek`
  source per sample). Forcing `Send` would mean unsafe stitching of
  the source between threads; not worth it.

- **No parallel column decompression.** Mirrors the writer's
  decision. Columns are independent zstd frames, so `rayon::scope`
  inside the block-decode loop is trivial to add; defer until
  benchmarks show it matters.

- **No mmap.** Same reasoning. Stage 3+ exercises random-access
  region queries only via the cohort merge's seek-forward pattern,
  which `pread`-style reads handle fine.

- **Sticky `poisoned` flag.** Required for `Iterator` correctness:
  once `next()` has returned `Some(Err(_))`, future calls must
  return `None`. The alternative (let the iterator try to recover)
  would mean leaving the source in an unknown state, which the
  spec's "errors must not pass silently" rule rules out.

- **`block_index()` is `pub`.** A small API surface for the cohort
  merge driver and `psp dump`. Read-only slice; no leak.

- **Reading the block header is a "grow until parse succeeds"
  loop.** Varint-framed headers do not have a length prefix at the
  outer level; the standard idiom for streaming varint sections is
  read-some, try-decode, grow-on-incomplete. The reader's helper
  caps the grow at a sensible upper bound (64 KiB is far more than
  any sane block header) so a malformed file cannot induce unbounded
  read.
