# Code Review: psp_reader
**Date:** 2026-05-13
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** `src/per_sample_caller/psp/reader.rs` (new file, commit `2d7f211`), plus the new error variants in `errors.rs` and re-exports in `mod.rs`
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the new `PspReader` implementation that shipped in commit `2d7f211`.
- **Reviewed against:** commit `2d7f211` on `main`.
- **In-scope files:**
  - `src/per_sample_caller/psp/reader.rs` (1864 lines including tests)
  - `src/per_sample_caller/psp/errors.rs` — new variants (`BlockHeaderInvariantKind::SnapshotMismatch`, `PspReadError::BlockIndexOffsetInvalid`, `BlockIndexChromOutOfRange`, `BlockIndexPosOutOfRange`, `NonFiniteFloat`, `PhaseChainConsistency`, plus the sub-enum `PhaseChainConsistencyKind`)
  - `src/per_sample_caller/psp/mod.rs` — re-exports of `PspReader`, `RecordsIter`, `PhaseChainConsistencyKind`
- **Deliberately out of scope:**
  - The writer (`writer.rs`) and the decoder primitives it composes (`block.rs`, `index.rs`, `trailer.rs`, `header.rs`, `varint.rs`) — all reviewed previously and recently re-audited (see `ia/reviews/psp_2026-05-13.md`). Referenced only where the reader's expectations cross those boundaries.
  - The CLI (`src/main.rs`) and the rest of the pipeline.
- **Categories dispatched (7 in parallel):**
  - `reliability` — always-on; the reader is the consumer-side parser of a potentially-untrusted on-disk format.
  - `errors` — always-on; six new error variants and 16 `.expect`/`.unwrap` sites in non-test code.
  - `naming` — always-on; new public types `PspReader`, `RecordsIter`, plus private `DecodedBlock`, `RangeClamp`, helpers.
  - `idiomatic` — always-on; ~620 lines of new non-test Rust.
  - `refactor_safety` — always-on; new layer between the decoder primitives and the consumer.
  - `smells` — always-on.
  - `extras` — applies (parser/validator on **untrusted input**, **stable-on-disk format** decode half, **public-ish library surface**, hot-path candidate feeding Stage 4 grouping).
- **Categories not dispatched:**
  - `defaults` — no new defaults or configuration.
  - `unsafe_concurrency` — `#![forbid(unsafe_code)]` at crate level; no new `Arc`/`Mutex`/atomics/threads.
  - `tooling` — `Cargo.toml` unchanged.

## 2. Verdict

**Request-changes.** Three Blockers and twelve Majors. The Blockers are the kind of issue a reader-side review exists to catch: hand-crafted bytes can drive a `panic!` (twice) or a multi-GiB allocation, instead of returning a typed `PspReadError`. The Majors include refactor-safety holes (silent column-drop, `SlotId`-widening would break the wire format invisibly), error-typing regressions (one variant misused, one error message factually wrong), and a hot-path allocation pattern symmetric to the one the writer just optimised away. None of the Blockers requires structural redesign — each has a localised fix sketched below.

## 3. Execution status

| Command | Result |
|---|---|
| `cargo fmt --check` | exit 0 (clean) — re-quoted from the implementation report |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 (clean) |
| `cargo test --lib per_sample_caller::psp::reader` | 22 passed, 0 failed, 0 ignored |
| `cargo doc --no-deps --lib` | succeeded, 0 warnings |
| `cargo audit` | not run (project-policy: skipped on this review) |

These commands were not re-run during synthesis; the per-category subagents accepted the orchestrator's quoted output. Findings labelled "Needs verification": **0** (every finding cites a file and line that was read).

## 4. Open questions and assumptions

1. **Q (B1, B2, M2):** Is shipping a typed `PspReadError::MissingRequiredColumnInManifest` variant in this slice acceptable, or does the v1.0 schema commitment require the column-coverage check to be silent-panic-free without introducing a new variant? — recommendation is to add the variant; this is consistent with the existing `MissingRequiredColumn` for the TOML side.
2. **Q (M10, M16):** The writer's S1+S2 perf review (`perf_psp_writer_2026-05-13.md`) changed the active-chain-slot store from `BTreeSet` to a sorted `Vec`. Is the reader expected to mirror that decision now, or wait for a reader-side perf bench? — recommendation is to mirror, *and* add a `criterion` micro-bench so the next perf review has a baseline.
3. **A (B2):** Assumes per-block byte budget = `next_block_offset - end_of_block_header`. The block index records `block_offset` only; the end of the block header is `source.stream_position()` after `read_block_header` returns. Both are available without re-reading from disk.
4. **A (M14):** Assumes `SlotId` will eventually widen (the writer-side review flagged this earlier as a forward-compat concern); if `SlotId` is frozen at `u16` for v1.0+v1.x, the fix collapses to a one-line type alias change with no wire impact.

## 5. Top 3 priorities

1. **B1** — `decode_block_payload`'s 12 `.expect()`s + `decode_one_column`'s `allele_seq_len.expect()` panic on hand-crafted input. Add a per-block manifest required-tag coverage check (the symmetric layer to `cross_check_against_registry`) so the panic surface is gone before the dispatch loop runs. **One fix closes two sites.**
2. **B2** — `compressed_len` drives `read_exact` allocation up to 4 GiB per column with no file-derived bound. Thread the block's byte budget into `decode_one_column` and reject any `compressed_len` exceeding the remaining budget before allocating. **Single attack-surface concession we'd regret leaving open.**
3. **B3** — Cross-column `sum(n_alleles[i]) == n_total_alleles` is unenforced; a single varint flip in `n_alleles` panics with index-out-of-bounds in `materialise_next_record`. Verify after `decode_block_payload`; emit a typed `BlockHeaderInvariantKind::NAllelesSumMismatch`.

## 6. Findings

### Blockers

#### B1: src/per_sample_caller/psp/reader.rs:823-846, 966 — Per-block manifest is not coverage-checked against the v1.0 required-column set; missing required tags drive `.expect` panics

**Confidence:** High.
**Categories:** refactor_safety, errors, smells (12 repeated `.expect`s), extras (cross-category observation on the same shape).

**Problem.** `decode_block_payload` walks `header.manifest`, dispatches each entry through `decode_one_column`, and ends with twelve `.expect("<column> required by v1.0 schema")` calls (lines 829-846). The comment at lines 819-822 says these are "structurally unreachable on a schema-valid file" because `cross_check_against_registry` enforces required-column presence. That justification is **wrong**: `cross_check_against_registry` (header.rs:679-708) operates on the file-level TOML `[[column]]` array, not on the per-block manifest. The per-block manifest is decoded independently in `block.rs:599 decode_block_header`, and only four invariants are validated (`n_records >= 1`, `n_total_alleles >= n_records`, `active_chain_slots strictly ascending`, `manifest tags strictly ascending`). There is no check that the manifest contains every required tag.

A hand-crafted file with a per-block manifest that omits, say, tag `0x01` (delta-pos) passes header parse, passes `decode_block_header`, walks the manifest in `decode_block_payload` without ever assigning `delta_pos`, and **panics** at line 829.

The same anti-pattern fires a second time at line 966 — `decode_one_column`'s `ColumnKey::AlleleSeq` arm `.expect`s that `allele_seq_len` is `Some`. The justification ("manifest tags strictly ascending; 0x03 < 0x04") proves "if both are present, lens comes first", **not** "both are present". A manifest carrying only tag 0x04 (allele-seq) passes every existing block-header invariant and panics here.

**Why it matters.** A panic on user-controlled bytes is the textbook reader failure mode: it crashes the host process rather than returning a typed `PspReadError`. Downstream code calling `PspReader::new(...)?` and iterating cannot catch a panic in a `Result` chain. The F1-F8 negative-test surface targets exactly this kind of crafted-input robustness, yet the most central invariant (required-tag coverage) is unenforced.

**Suggested fix.** Add a required-tag coverage check at the top of `decode_block_payload` (or, equivalently, at the end of `decode_block_header`). New error variant `PspReadError::MissingRequiredColumnInManifest { name: String, tag: u16 }` parallel to the existing `MissingRequiredColumn` for the TOML side.

```rust
fn decode_block_payload<R: Read>(
    source: &mut R,
    header: &BlockHeader,
) -> Result<DecodedBlock, PspReadError> {
    // Coverage check: every v1.0 required tag must be present in
    // the per-block manifest. Symmetric with the header-TOML check
    // in cross_check_against_registry; the two are independent
    // layers and both need enforcing.
    for def in super::registry::V1_0_COLUMNS {
        if def.required && !header.manifest.iter().any(|e| e.tag == def.tag) {
            return Err(PspReadError::MissingRequiredColumnInManifest {
                name: def.name.to_string(),
                tag: def.tag,
            });
        }
    }
    // ...rest unchanged. The `.expect("…")`s now hold structurally because the
    // loop above has verified each tag's presence, and `decode_one_column`'s
    // `AlleleSeq` arm `expect` holds because `allele-seq-len` is required.
}
```

Test: hand-build a `BlockHeader` whose manifest omits one required tag and assert the new error variant fires. After this fix the 12 `.expect()`s plus the `allele_seq_len.expect()` are genuinely structurally unreachable — and the inline prose can be retagged `// PANIC-FREE:` per `ia/skills/code_review/errors.md` (resolves M1, M2, M3 as a side effect).

---

#### B2: src/per_sample_caller/psp/reader.rs:889, 926 — Untrusted `compressed_len` drives `read_exact` allocation up to 4 GiB with no file-derived bound

**Confidence:** High.
**Categories:** extras (parser-security trigger).

**Problem.** `decode_one_column` allocates `compressed = vec![0u8; entry.compressed_len as usize]` (line 926) and the unknown-column-skip path `sink = vec![0u8; entry.compressed_len as usize]` (line 889) directly from `ColumnManifestEntry.compressed_len`. `compressed_len` is a `u32` decoded from the block header (`block.rs:637`), so on a 64-bit host the upper bound is `u32::MAX` ≈ 4 GiB **per column**. A block with 12 columns can pre-reserve ~48 GiB before any byte is read. The trailer cross-check at line 127 bounds the *index* region against file size, and `read_block_header` caps the *block header* read at 64 KiB (line 56), but no equivalent bound applies to a column's compressed payload. The block's byte extent (this `block_offset` → next `block_offset` or `trailer.index_offset`) is known up front and could cap the allowable `compressed_len`, but is not consulted.

**Why it matters.** A single 30-byte forged manifest entry can drive a multi-GiB allocation that either OOM-kills a pipeline worker or stalls it. The prior reviews hardened the decoder primitives (`decode_list_column`, `decode_bytes_split`, `decode_block_header`) against exactly this shape; the reader integration layer regressed it.

**Suggested fix.** Compute the block's byte budget once per block (next-block-offset minus end-of-block-header) and thread it through `decode_block_payload` → `decode_one_column`, tracking remaining bytes as columns consume them.

```rust
// In load_next_block, after read_block_header succeeds:
let block_end_offset = self
    .reader
    .index
    .get(self.cur_block_idx + 1)
    .map(|e| e.block_offset)
    .unwrap_or(self.reader.trailer.index_offset);
let header_end_pos = self
    .reader
    .source
    .stream_position()
    .map_err(io_err("post-block-header pos"))?;
let column_budget = block_end_offset
    .checked_sub(header_end_pos)
    .ok_or_else(|| /* malformed block index */)?;

// In decode_one_column, before line 926:
if entry.compressed_len as u64 > remaining_budget {
    return Err(PspReadError::ColumnTruncated {
        column: column_name.to_string(),
        decoded: 0,
        expected: entry.compressed_len as usize,
    });
}
```

Test: forge a manifest with `compressed_len = u32::MAX` against a small block region; assert a specific typed error variant fires before allocation. Pair with a `cargo-fuzz` target (see M17).

---

#### B3: src/per_sample_caller/psp/reader.rs:574-605 — Per-allele cross-column inconsistency can panic on corrupt / forged input

**Confidence:** High.
**Categories:** reliability.

**Problem.** `materialise_next_record` computes `allele_end = self.next_allele_in_block as usize + n_alleles[i] as usize` and then indexes `block.allele_seqs[j]`, `block.allele_obs_count[j]`, etc. for `j in allele_start..allele_end`. Per-allele columns have length `header.n_total_alleles` (enforced by `decode_scalar_column` / `decode_list_column`'s `expected_count` argument). Block-header validation enforces `n_total_alleles >= n_records` (`BlockHeaderInvariantKind::AllelesLessThanRecords`) but **does not enforce `sum(n_alleles[i]) == n_total_alleles`**. A varint flip in the `n_alleles` column can make the sum exceed `n_total_alleles`, causing `block.allele_seqs[j]` to panic with index-out-of-bounds rather than producing a typed `PspReadError`. The reverse case (sum < n_total_alleles) is silent: the reader emits valid-looking records and ignores trailing per-allele entries, so the consumer (Stage 4 grouping) ingests truncated allele lists.

The spec's §"Per-block manifest agreement" rule says "Any over- or under-run is a hard error" — this is exactly such a mismatch, but at the *cross-column* level rather than within a single column.

**Why it matters.** A panic on adversarial input is a DoS vector at the reader; a silent under-run on forged input feeds wrong allele data downstream. Spec §10's "errors must not pass silently" principle is violated for both directions. `f3_torn_zstd_frame` (line 1719) admits `Zstd | UncompressedLenMismatch | ColumnTruncated | ColumnTrailingBytes | ColumnElementDecode` — a flip that lands inside the `n_alleles` column and survives the zstd frame check currently surfaces as a panic, not as any of those.

**Suggested fix.** In `decode_block_payload`, after all columns are decoded and before constructing `DecodedBlock`, assert `sum(n_alleles) == n_total_alleles` and return a new `BlockHeaderInvariantKind::NAllelesSumMismatch { n_total_alleles, sum_n_alleles }`:

```rust
let n_alleles_ref = n_alleles
    .as_ref()
    .expect("n-alleles column required by v1.0 schema");
let n_alleles_sum: u64 = n_alleles_ref.iter().sum();
if n_alleles_sum != header.n_total_alleles as u64 {
    return Err(PspReadError::BlockHeaderInvariant {
        kind: BlockHeaderInvariantKind::NAllelesSumMismatch {
            n_total_alleles: header.n_total_alleles,
            sum_n_alleles: n_alleles_sum,
        },
    });
}
```

Pair with the missing test `materialise_next_record_returns_typed_error_when_n_alleles_sum_exceeds_n_total_alleles` (see §8).

### Major

#### M1: src/per_sample_caller/psp/reader.rs:535 — `materialise_next_record` `.expect()` without a `// PANIC-FREE:` invariant comment

**Confidence:** High.
**Categories:** errors.

**Problem.** `self.cur_block.as_ref().expect("materialise_next_record requires a loaded block")` carries a message describing the precondition but no `// PANIC-FREE:` comment naming the invariant that gates the call. The single caller (`<RecordsIter as Iterator>::next` at line 631) gates this behind `if let Some(block) = &self.cur_block && self.next_record_in_block < block.n_records`, but that invariant lives implicitly in the caller. A future refactor that adds a second call site loses the guard silently.

**Suggested fix.** Tag the prose:

```rust
// PANIC-FREE: the only caller (`<Self as Iterator>::next`) gates this call
// behind `if let Some(block) = &self.cur_block && ...`, so `cur_block` is
// always `Some` here.
let block = self
    .cur_block
    .as_ref()
    .expect("materialise_next_record requires a loaded block");
```

Better still: restructure the function to take `&DecodedBlock` directly, lifting the unwrap into the caller and turning the invariant into a type-level fact.

---

#### M2: src/per_sample_caller/psp/reader.rs:190 — `try_into().unwrap()` on a `[u8; 12][4..12]` slice; replace with the panic-free form

**Confidence:** High.
**Categories:** errors.

**Problem.** `let body_len = u64::from_le_bytes(head_prefix[4..12].try_into().unwrap());`. The conversion is infallible (slicing exactly 8 bytes out of a 12-byte array always yields an 8-byte slice), but the rule demands a `// PANIC-FREE:` comment for non-test `unwrap`. The cleaner idiom removes the `unwrap` entirely.

**Suggested fix.**

```rust
let mut buf = [0u8; 8];
buf.copy_from_slice(&head_prefix[4..12]);
let body_len = u64::from_le_bytes(buf);
```

---

#### M3: src/per_sample_caller/psp/reader.rs:738-743 — `VarintError::Overflow` is the wrong `#[source]` for "header exceeds 64 KiB read cap"

**Confidence:** High.
**Categories:** errors.

**Problem.** When the block-header grow loop hits the `BLOCK_HEADER_READ_CAP` ceiling without `decode_block_header` succeeding, the function returns `PspReadError::BlockHeaderField { field: "header exceeds read cap", source: VarintError::Overflow }`. `VarintError::Overflow` is documented as "More than 10 LEB128 continuation bytes were consumed", which is a per-varint failure mode. The rendered log line — `"block header field header exceeds read cap: varint overflow: continuation bytes exceeded the 10-byte cap"` — is doubly misleading: a buffer-cap failure at a different layer is presented as a varint bug.

**Suggested fix.** Add a dedicated variant under `PspReadError`:

```rust
/// The block-header grow-on-incomplete loop hit its read cap
/// without `decode_block_header` succeeding.
#[error("block header exceeds {cap}-byte read cap (read {consumed} bytes without a successful decode)")]
BlockHeaderExceedsCap { cap: usize, consumed: usize },
```

Return it at the cap-hit site. Pair with the missing test `read_block_header_returns_overflow_error_on_pathological_continuation_stream` (see §8).

---

#### M4: src/per_sample_caller/psp/reader.rs:712-715 — `field: "header truncated mid-decode"` is a description, not a field name; doubled prose in the rendered error

**Confidence:** High.
**Categories:** errors.

**Problem.** The return `Err(PspReadError::BlockHeaderField { field: "header truncated mid-decode", source: VarintError::Truncated })` produces `"block header field header truncated mid-decode: varint truncated: continuation bit set on the final available byte"` — doubly redundant and grammatically broken. Every other `BlockHeaderField` instance uses an actual field name (`"chrom_id"`, `"n_records"`, etc.).

**Suggested fix.** Capture the source position into a dedicated variant:

```rust
/// EOF hit while pulling more bytes into the block-header buffer.
#[error("block header truncated at offset {offset} after reading {consumed} bytes")]
BlockHeaderTruncated { offset: u64, consumed: usize },
```

The existing `header_start` offset (line 700) is exactly the diagnostic data the operator needs.

---

#### M5: src/per_sample_caller/psp/reader.rs:160-167 — `BlockIndexOffsetInvalid` with `min: 0` is factually wrong

**Confidence:** High.
**Categories:** errors.

**Problem.** Step 4's "every entry's offset must be strictly less than `trailer.index_offset`" check reports `min: 0`. The variant Display is `"block index entry {block}: block_offset {offset} outside expected range [{min}, {max})"`. An operator reading `"block_offset 1234 outside expected range [0, 5678)"` reasons "1234 >= 0 and < 5678, so this error is impossible" — the message contradicts the failure. The real lower bound is `header_end_offset` (enforced at step 6, lines 209-218), but step 4 fires before the header is parsed.

**Suggested fix.** Delay the upper-bound check until step 6 so the lower bound is always honest, and report `[header_end_offset, trailer.index_offset)` at every offset failure:

```rust
// Move all three offset checks to step 6, after the header is parsed.
// The error variant then always carries the correct
// [header_end_offset, trailer.index_offset) interval.
```

Alternatively split into three distinct variants — one per failure mode (offset >= index_offset, offset < header_end_offset, non-increasing) — and stop carrying a `min` that is sometimes lying.

---

#### M6: src/per_sample_caller/psp/errors.rs:494-503, reader.rs:798-805 — `NonFiniteFloat` variant is provably unreached today; F5 deferred test risk

**Confidence:** High.
**Categories:** errors.

**Problem.** The variant fires at lines 798-805 inside `decode_block_payload`'s `AlleleQSumLog` arm. The writer enforces finite-only q-sum at produce time, so the reader can only fire this variant on a hand-crafted or torn file. No test constructs such a file (F5 is deferred per the implementation report §6). Concrete drift risks:

- The loop returns the first non-finite value but `entry` is the array index into the per-allele column, not the `(record, allele)` pair the writer-side `InvalidRecordKind::NonFiniteQSum` carries.
- The column name string `"allele-q-sum-log"` is hard-coded (line 801) instead of using `column_name` from `def.name` already in scope.

**Suggested fix.** Two parts:

1. Replace `"allele-q-sum-log".to_string()` with `column_name.to_string()` so a registry rename doesn't break the message silently.
2. Land the F5 test (see §8 missing test #9). Requires the test-only `replace_column_payload_f64` helper from the implementation report's §6.

---

#### M7: src/per_sample_caller/psp/reader.rs:489-513 — Spec check #10's "snapshot empty at first block of chromosome" is not enforced in random-access mode

**Confidence:** High.
**Categories:** reliability.

**Problem.** Spec §"Header-binary consistency" check #10 says, unconditionally ("checked at each block"): "The snapshot is empty whenever the block is the first block of its chromosome." In `load_next_block`:

- `RangeClamp::None` (sequential) runs the carry-in vs snapshot equality check (check #11). For a chrom-change boundary, this *indirectly* enforces "snapshot empty at first block of chromosome" — but only because the writer flushes on chrom change. A forged file whose previous-block phase-chain markers leave a non-empty active set, and whose new-chrom block carries that same non-empty snapshot, passes check #11 but violates check #10.
- `RangeClamp::Window` (region) does `clear() + extend(snapshot)` with no first-block-of-chrom check. A forged or genuinely buggy file is silently trusted.

**Suggested fix.** Run an explicit "first block of chromosome → snapshot must be empty" check inside `load_next_block`, in both `RangeClamp` arms, before applying / comparing against the carry-in set:

```rust
let is_first_of_chrom = self.cur_block_idx == 0
    || self.reader.index[self.cur_block_idx].chrom_id
        != self.reader.index[self.cur_block_idx - 1].chrom_id;
if is_first_of_chrom && !block_header.active_chain_slots.is_empty() {
    return Err(PspReadError::BlockHeaderInvariant {
        kind: BlockHeaderInvariantKind::NonEmptySnapshotAtChromStart {
            block: self.cur_block_idx,
            snapshot_len: block_header.active_chain_slots.len(),
        },
    });
}
```

Test: `load_next_block_rejects_non_empty_snapshot_at_chrom_change_in_region_mode` (see §8).

---

#### M8: src/per_sample_caller/psp/reader.rs:441 — Reader uses `BTreeSet<SlotId>` instead of the writer's sorted `Vec<SlotId>`; documented divergence

**Confidence:** High.
**Categories:** idiomatic, refactor_safety (cross-category from perf).

**Problem.** `RecordsIter::active_chain_slots` is a `BTreeSet<SlotId>` (line 441), with set-equality against the snapshot at lines 493-503 and per-record `contains`/`insert`/`remove` at lines 540-565. The writer's symmetric `IngestState::active_slots` is a sorted `Vec<SlotId>` with `binary_search`, justified inline at `writer.rs:130-141` ("S1+S2"): the typical active count is ~8, so one contiguous allocation beats a tree, and the block-start snapshot is a memcpy. The reader's access pattern is identical; the divergence is undocumented. Snapshot-equality at line 493 currently builds the snapshot as a `BTreeSet` just to compare — a sorted Vec skips that intermediate.

**Suggested fix.** Mirror the writer:

```rust
// reader.rs:441
active_chain_slots: Vec<SlotId>,

// reader.rs:493-503 — snapshot is already a sorted Vec<u16> in
// BlockHeader::active_chain_slots; direct equality:
if block_header.active_chain_slots != self.active_chain_slots {
    return Err(PspReadError::BlockHeaderInvariant {
        kind: BlockHeaderInvariantKind::SnapshotMismatch { ... },
    });
}

// reader.rs:509-511 — region mode re-init
self.active_chain_slots.clear();
self.active_chain_slots
    .extend_from_slice(&block_header.active_chain_slots);

// reader.rs:540-565 — membership / mutation via binary_search.
```

Add the same "S1+S2 rationale lives in writer.rs:135-140" comment so the tree doesn't get re-introduced.

---

#### M9: src/per_sample_caller/psp/reader.rs:76-77, 373-374, 377-378 — Three `#[allow(dead_code)]` fields are genuinely dead, not "consumed by a later slice"

**Confidence:** High.
**Categories:** smells, idiomatic, refactor_safety (convergent — three subagents).

**Problem.** Three struct fields carry `#[allow(dead_code)]` with comments deferring justification to a future slice that has already shipped:

- `PspReader::header_end_offset` (line 77) — "Slice 5's region init may also use it". Slice 5 has shipped; `region_records` does not read this field. The local `header_end_offset` binding inside `new` does the validation work and is copy-stored on the struct without ever being read again.
- `DecodedBlock::n_total_alleles` (line 374) — "sanity-checked against the manifest in Slice 3+". Slice 3 shipped; no read site exists for the struct field. The function-local `n_total_alleles` parameter at line 900 does the manifest cross-check.
- `DecodedBlock::active_chain_slots_snapshot` (line 378) — "consumed by Slice 4's continuity check." Slice 4 shipped; the check reads `block_header.active_chain_slots` directly *before* `DecodedBlock` is built. The field is set (with a `.clone()`!) but never observed.

**Why it matters.** Each `#[allow(dead_code)]` is a deferred decision that has now expired. The `active_chain_slots_snapshot` clone is dead allocation work per block. The "REVIEW WHEN SLICE N LANDS" pattern loses its signalling value once half are stale.

**Suggested fix.** Delete all three fields, their `#[allow(...)]` attributes, the field assignments in the `DecodedBlock { ... }` literal at lines 827-828, and the `header_end_offset:` field on the `Self { ... }` literal at line 262. Keep the local `header_end_offset` binding in `new` for its validation uses.

---

#### M10: src/per_sample_caller/psp/reader.rs:786-815, 867-868 — `decode_block_payload`'s match dispatch is structurally **not** exhaustive over `ColumnKey`; `DecodedColumn::Unknown` is a silent-drop escape hatch

**Confidence:** High.
**Categories:** refactor_safety.

**Problem.** The implementer's stated design goal (echoed in this review's prompt) is: "adding a `ColumnKey::Foo` variant forces the user to also add `DecodedColumn::Foo` (compile error there too)." Tracing the flow: `decode_one_column` (line 944) has an exhaustive match on `def.key: ColumnKey`. Adding `ColumnKey::Foo` forces a new arm there. But that arm must produce some `DecodedColumn` value, and **nothing forces the implementer to add `DecodedColumn::Foo`** — they can route a new key through `DecodedColumn::Unknown` (the "drop bytes silently" arm) and the code compiles. The match at lines 786-815 accepts `Unknown` and does nothing, so a new column is invisible to `DecodedBlock`.

The writer's `encode_column_into` is hardened (every `ColumnKey` arm produces concrete bytes). The reader has a load-bearing variant that swallows any new key the user forgets to wire up.

**Suggested fix.** Remove the `DecodedColumn::Unknown` catch-all. Move the unknown-tag branch in `decode_one_column` (lines 886-895) into a returning `Option<DecodedColumn>` (`None` = unknown), so the dispatch loop becomes `Some(col) => match col { /* twelve variants, no Unknown */ }` plus `None => /* skip */`. Then `DecodedColumn` is closed over the registry: adding a `ColumnKey` adds an exhaustive-match miss at three sites.

```rust
/// Returns None for unknown optional columns (bytes already consumed);
/// Some(decoded) for any known v1.0 column.
fn decode_one_column<R: Read>(...) -> Result<Option<DecodedColumn>, PspReadError> {
    let Some(def) = lookup_by_tag(entry.tag) else {
        let mut sink = vec![0u8; entry.compressed_len as usize];
        source.read_exact(&mut sink).map_err(io_err("unknown optional column payload"))?;
        return Ok(None);
    };
    // ...remove DecodedColumn::Unknown entirely
}
```

---

#### M11: src/per_sample_caller/psp/reader.rs:792-808 — F5 finite-float sweep is hard-coded to `AlleleQSumLog`; a future F32/F64 column needs the reviewer to remember to add a parallel sweep

**Confidence:** High.
**Categories:** refactor_safety.

**Problem.** The spec-check-#9 enforcement at lines 798-805 fires only for `DecodedColumn::AlleleQSumLog`. The decision of *which* columns need the constraint lives nowhere structural — no `ColumnDef` field, no `ColumnKey::has_finite_constraint()` method. Adding a new F-typed column means adding a `ColumnKey` variant (which `decode_one_column`'s match flags) and a `DecodedColumn` variant (currently optional — see M10), but nothing flags that the finite sweep is missing.

**Suggested fix.** Add a `finite_constraint: bool` field to `ColumnDef`. Default it to `false` for non-float columns; set it to `true` for `allele-q-sum-log` and future F-typed columns. In `decode_block_payload`, after the per-column dispatch, run a generic sweep based on the registry flag.

---

#### M12: src/per_sample_caller/psp/block.rs:554, reader.rs:441 — No compile-time link between `BlockHeader::active_chain_slots: Vec<u16>` and `SlotId`; widening `SlotId` silently breaks the wire format

**Confidence:** High.
**Categories:** refactor_safety.

**Problem.** `SlotId = u16` is defined once in `pileup/slot_allocator.rs:21`. The reader's running set uses `BTreeSet<SlotId>` (line 441). `BlockHeader::active_chain_slots` is hard-coded to `Vec<u16>` (block.rs:554) — not `Vec<SlotId>`. The active-slot collect at reader.rs:493-494 does `block_header.active_chain_slots.iter().copied().collect::<BTreeSet<SlotId>>()`, which silently works because `SlotId` *is* `u16` today. If `SlotId` is ever widened to `u32`, the writer at block.rs:584-586 (`slot.to_le_bytes()` on `u16` → 2 bytes) silently breaks the file format: the writer emits 4-byte ids, the reader parses `u16`s.

**Suggested fix.** Replace `Vec<u16>` with `Vec<SlotId>` in `BlockHeader::active_chain_slots` (block.rs:554) and the surrounding slot-handling paths. Make the snapshot-decode loop at block.rs:611-624 generic over `<SlotId as WireScalar>::FIXED_BYTE_WIDTH` so widening forces every assumption-of-`u16` site to flag.

---

#### M13: src/per_sample_caller/psp/reader.rs:578-613 — `materialise_next_record` clones four `Vec`s per record where `mem::take` would work

**Confidence:** Medium (no reader benchmark yet; perf claim extrapolated from the per-record arithmetic and the writer's recently-landed S1+S2 work).
**Categories:** idiomatic, extras (hot-path trigger), smells (cross-category observation).

**Problem.** Four `.clone()` calls per record:
- `block.allele_seqs[j].clone()` (line 595) — `Vec<u8>` per allele
- `block.allele_chain_slots[j].clone()` (line 603) — `Vec<u16>` per allele
- `block.new_chain_slots[i].clone()` (line 610) — `Vec<u16>` per record
- `block.expired_chain_slots[i].clone()` (line 611) — `Vec<u16>` per record

On a 1-allele 1M-record stream that's ~4M `malloc`s the consumer immediately frees the originals. The block is consumed forwards-only (`next_record_in_block` increments at line 616, the block is dropped wholesale at line 660, no retry path revisits index `i` or `j`), so each slot is read exactly once before the whole `DecodedBlock` drops — `std::mem::take` yields the same `PileupRecord` shape with zero per-record allocation.

The writer just landed a Major+Hot-path round (`ia/reviews/perf_psp_writer_2026-05-13.md`) collapsing `Vec<Vec<SlotId>>` precisely to kill per-record inner-`Vec` mallocs. Re-allocating the same shape on the read side is the symmetric regression.

**Suggested fix.** Take `&mut DecodedBlock` and `mem::take` the four slots. Validation loops at lines 540-557, 580-592 must run *before* the `mem::take`s (the existing structure already does this).

```rust
let block = self.cur_block.as_mut().expect("...");
let i = self.next_record_in_block as usize;
// ... validation against active_chain_slots (unchanged) ...

let mut alleles = Vec::with_capacity(n_alleles_here);
for j in allele_start..allele_end {
    alleles.push(AlleleObservation {
        seq: std::mem::take(&mut block.allele_seqs[j]),
        support: AlleleSupportStats { /* scalars unchanged */ },
        chain_slots: std::mem::take(&mut block.allele_chain_slots[j]),
    });
}
let record = PileupRecord {
    chrom_id: block.chrom_id, pos,
    new_chains: std::mem::take(&mut block.new_chain_slots[i]),
    expired_chains: std::mem::take(&mut block.expired_chain_slots[i]),
    alleles,
};
```

Pair with a `criterion` micro-bench at `benches/psp_reader_perf.rs` reading a 1M-record file, with `// REGRESSION THRESHOLD: 5%`.

---

#### M14: src/per_sample_caller/psp/reader.rs (no fuzz target, no reader-side proptest) — Reader's "never panic / OOM / hang on malformed bytes" contract is not property-tested

**Confidence:** High.
**Categories:** extras.

**Problem.** The writer ships proptests (`varint.rs:271`, `block.rs:1448-1530`, `index.rs:463`, `trailer.rs:178`) but they cover encode→decode round-trips, not "arbitrary file bytes handed to `PspReader::new` either error or succeed, never panic". The reader covers 8 corruption fixtures (F1-F8, R10) — each manually constructed against the small fixture. Each fixture hits only the corruption shape the implementer thought of. The two top findings in the prior writer review (`decode_list_column`'s `u64::MAX` allocation, `decode_bytes_split`'s wrap-on-sum) would not have been caught by the round-trip proptests as they stood; both are first-line "obvious" once spotted, and both have unit tests *after*. A reader-side fuzz target catches the next round-of-obvious before it ships.

**Suggested fix.** Add a `fuzz/` directory with a single target. Run periodically (a few CPU-hours on CI cron) and triage every panic / OOM / hang to a typed error.

```rust
// fuzz/fuzz_targets/psp_reader_no_panic.rs
#![no_main]
use libfuzzer_sys::fuzz_target;
use std::io::Cursor;
use merge_per_sample_vcfs::per_sample_caller::psp::PspReader;
fuzz_target!(|data: &[u8]| {
    if let Ok(mut reader) = PspReader::new(Cursor::new(data)) {
        for r in reader.records().take(10_000) {
            let _ = r;
        }
    }
});
```

Pair with a `tests/corpus_psp_reader/` directory committing every interesting input the fuzzer finds.

---

#### M15: src/per_sample_caller/psp/reader.rs:91-264 — `PspReader::new`'s 7-step open sequence has no structural marker forcing a new step to be added at the right place

**Confidence:** Medium.
**Categories:** refactor_safety.

**Problem.** Steps numbered in a comment (1, 2, 3, 4, 5, 6, 7). Adding a v1.1 header-annex parse step means a developer thinks to add it *and* puts it after step 5 and before step 7. Nothing in the type system flags that step 6's cross-check should now also cover the annex, or that step 7's seek target should account for it. `header_end_offset` doubles as the seek target on line 254 and the manifest lower-bound on line 209 — adding an annex would shift `header_end_offset`'s semantics and a refactorer would need to remember to update both call sites.

**Suggested fix.** Extract each step into a named helper (`fn read_trailer`, `fn read_index`, `fn parse_header`, `fn cross_check_index`, `fn seek_to_first_block`). The open sequence then reads as a literal list of typed calls.

---

#### M16: src/per_sample_caller/psp/reader.rs:130-132 — Unchecked `u64 + u64` in the trailer error-path arithmetic

**Confidence:** High.
**Categories:** errors (Minor), extras (Major), refactor_safety (Major), reliability (Minor) — convergent across four subagents.

**Problem.** Line 127 correctly does `trailer.index_offset.checked_add(trailer.index_byte_length)` for the equality check. Line 131 (inside the error-construction branch entered *exactly when* `checked_add` returned `None`) recomputes the same expression without `checked_add`:

```rust
trailing_bytes: trailer_start
    .saturating_sub(trailer.index_offset + trailer.index_byte_length)
    as usize,
```

In the `None`-overflow case the raw `+` panics in debug (debug-build CI panics on a hostile file with `index_offset = u64::MAX - 5`, `index_byte_length = 10`) and silently wraps in release (producing a misleading `trailing_bytes` count). The blast radius is the error message, not the control flow, but the project's own CI would trip on a hand-crafted input that reaches this branch.

**Suggested fix.**

```rust
let computed_end = trailer.index_offset.checked_add(trailer.index_byte_length);
if computed_end != Some(trailer_start) {
    let trailing_bytes = computed_end
        .and_then(|end| trailer_start.checked_sub(end))
        .unwrap_or(0) as usize;
    return Err(PspReadError::IndexTrailingBytes { trailing_bytes });
}
```

Test: `pspreader_new_rejects_trailer_with_overflowing_arithmetic_without_panicking` (see §8).

### Minor

- **Mi1: src/per_sample_caller/psp/reader.rs:495-503 — `SnapshotMismatch` reports only set lengths, not content delta.** The actual comparison is `BTreeSet` equality (content-aware) but the error carries `snapshot_len` and `expected_len` only. A pair of `[7]` vs `[9]` produces `snapshot_len: 1, expected_len: 1` — useless for diagnosis. Extend the variant to carry the two sets (capped at ~8 slots for log safety). (Category: refactor_safety.)

- **Mi2: src/per_sample_caller/psp/errors.rs:509-514 — `PhaseChainConsistency::record_in_block` asymmetric with the writer's `PhaseChainMarkerInconsistency::record_index`.** Writer side: `record_index: u64` (file-global). Reader side: `record_in_block: u32` (block-local). Cross-validation code wanting "the offending record" must branch on which side produced the error. Recommendation: add a derived `record_index: u64` (writer-stream-global) by carrying a running record count in `RecordsIter`. (Category: errors.)

- **Mi3: src/per_sample_caller/psp/reader.rs:115, 119, 140, 144 — Duplicated `io_err` context strings across seek+read pairs.** `io_err("trailer")` is used for both the seek (line 115) and the read (line 119); same for "block index" (lines 140, 144). Use distinct strings ("trailer seek" / "trailer read" / "block index seek" / "block index read"). (Category: errors.)

- **Mi4: src/per_sample_caller/psp/reader.rs:175 — `prev.block_offset + 1` could overflow when `trailer.index_offset == u64::MAX`.** Replace with `prev.block_offset.saturating_add(1)`. (Category: errors.)

- **Mi5: src/per_sample_caller/psp/reader.rs:567-572 — `delta as u32` truncation silently accepts forged 64-bit deltas.** `delta_pos` is decoded as `u64`; the cast to `u32` silently truncates a varint with top 32 bits set; `saturating_add` then caps `pos` at `u32::MAX`. Reject `delta > u32::MAX` at decode time via `u32::try_from`, returning `PspReadError::ColumnElementDecode { source: ScalarDecodeError::VarintOverflow, ... }`. (Categories: reliability, errors — convergent.)

- **Mi6: src/per_sample_caller/psp/reader.rs:308-321 — `region_records` silently treats out-of-range `chrom_id` and `start > end` as empty windows with no diagnostic.** Either document the silent-empty semantics in the docstring, or return `Result<RecordsIter, PspReadError>` and reject. Cheapest fix is to document. (Category: reliability.)

- **Mi7: src/per_sample_caller/psp/reader.rs:1851-1862 — `f8_first_block_non_empty_snapshot` swallows the open-time failure path with no assertion.** The early-return on `Err` is unconditional; a regression that makes open silently accept the forged file would pass the test. Track whether open succeeded; assert at end-of-test that the file was rejected at either open or iteration. (Category: reliability.)

- **Mi8: src/per_sample_caller/psp/reader.rs:285-287 — `PspReader::trailer` is `pub` but `#[doc(hidden)]`; the visibility is a type smell.** Hidden-but-public means the symbol is part of the crate's stable API for SemVer but invisible on docs.rs. Demote to `pub(crate)` — the test + future `psp dump` both live in the crate. (Category: extras.)

- **Mi9: src/per_sample_caller/psp/reader.rs:97-322 — No `# Examples` doctest on any public method; `PspReader::new` has no `# Errors` section.** Add a doctest on `PspReader::new` and `region_records`. Document the open-time failure variants (`BadHeadMagic`, `BadTrailerMagic`, `IndexChecksum`, `BadHeaderLength`, `UnsupportedFormatVersion`, `BlockIndexOffsetInvalid`, `BlockIndexChromOutOfRange`, `BlockIndexPosOutOfRange`). Consider `#![deny(missing_docs)]` once the round is complete. (Category: extras.)

- **Mi10: src/per_sample_caller/psp/reader.rs:877-1014 — `decode_one_column` does several distinct things across ~138 lines, 4 phases.** Hoist the `AlleleSeqLen` per-entry cap into `validate_allele_seq_len(&[u64], &str) -> Result<(), PspReadError>` (matches the finite-float sweep's shape, gathering "post-decode per-column validator" into one section). Optionally fold the spec-check-#7 `AlleleSeq` special-case into a renamed `predict_uncompressed_len`. (Category: smells.)

- **Mi11: src/per_sample_caller/psp/reader.rs:1229-1444 — Five round-trip tests (B1-B5) share the same setup → write → finish → read → collect skeleton.** Centralise in `fn round_trip(header, records, block_target) -> Vec<PileupRecord>` helper. (Category: smells.)

- **Mi12: src/per_sample_caller/psp/reader.rs (test module) — `PspReader::new(Cursor::new(bytes))` boilerplate appears 25 times.** Add `open_reader(bytes)` / `try_open_reader(bytes)` test-module helpers. (Category: smells.)

- **Mi13: src/per_sample_caller/psp/reader.rs:75, 291-293, 373, 376, 440, 867, 1663 — Comments reference development "Slice N" milestones that have all shipped.** Replace plan-relative phrasing with behaviour-relative phrasing. The dead fields at lines 75, 373, 376 drop with M9; the others rewrite in place. (Category: smells.)

- **Mi14: src/per_sample_caller/psp/reader.rs:1064-1224 — Two near-identical record/allele test helpers (`one_allele_record`, `allele`) plus three single-record builders (`record_open_slot`/`use`/`close`) that differ only in 2-3 fields.** Collapse to one general `record(...)`, one `allele(...)`, two trivial wrappers. (Category: smells.)

- **Mi15: src/per_sample_caller/psp/reader.rs:160-180 — `for (i, entry) in index.iter().enumerate()` with hand-rolled `if i > 0 { &index[i - 1] }` strict-monotonic check.** The "block_offset in range" check needs `i` for the error; the "strictly increasing" check is the textbook `index.windows(2)` shape. Split into two passes. (Category: idiomatic.)

- **Mi16: src/per_sample_caller/psp/reader.rs:315-320 — `match first { Some(idx) => ..., None => ... }` for "set field to value-or-default" is `unwrap_or` in disguise.** Replace with `it.cur_block_idx = first.unwrap_or(it.reader.index.len());`. (Category: idiomatic.)

- **Mi17: src/per_sample_caller/psp/reader.rs:642-654, 667-678 — Region-mode peek-ahead in `next()` reads `RangeClamp`'s internals; a method on `RangeClamp` would localise the logic.** Add `record_past_window` / `record_before_window` / `block_past_window` methods on `RangeClamp`. Three call sites do the same overlap check today. (Category: idiomatic.)

- **Mi18: src/per_sample_caller/psp/reader.rs:427, 431 — `u32` iterator counter fields force `as usize` at every index site (three sites in `materialise_next_record`).** Store as `usize`; convert to `u32` only at the error-variant construction site. (Category: idiomatic.)

- **Mi19: src/per_sample_caller/psp/reader.rs:1088-1862 — 22 of 22 reader tests carry plan-row prefixes `r1_`/`b1_`/`f1_`; regression on writer review Mi21 (mechanical `sed` rename of `h<N>_`).** This is the same pattern `psp_writer`'s review explicitly removed; `b5_r9_phase_chain_across_block_boundary` even shows the cost in microcosm — a hyphenated double-prefix encodes that one test satisfies two plan rows. Apply the same mechanical rename:

  ```bash
  sed -i 's/fn r\([0-9]\+\)_/fn /g; s/fn b\([0-9]\+\)_/fn /g; s/fn f\([0-9]\+\)_/fn /g' \
      src/per_sample_caller/psp/reader.rs
  ```

  Resolve the one collision by hand: `b5_r9_phase_chain_across_block_boundary` → `phase_chain_across_block_boundary`. Keep the plan-row reference inside the doc-comment above each test. (Category: naming.)

### Nits

- `reader.rs:877` — `decode_one_column`'s "one" mirrors the prose contrast with "every column" but `decode_column` is unambiguous. Matches `index.rs:131 decode_one_entry` either way — rename both or leave both.
- `reader.rs:327` — `first_block_overlapping` is a noun phrase; every other helper in the file is verb-shaped (`read_block_header`, `decode_block_payload`, `decode_one_column`). `find_first_overlapping_block` matches.
- `reader.rs:1021` — `predict_fixed_uncompressed_len` qualifies the wrong noun: "fixed" modifies the column, not the length. `predict_fixed_width_uncompressed_len` or `fixed_column_uncompressed_len`.
- `reader.rs:400-409` — `RangeClamp::None` shadows `Option::None` in the call site (`RecordsIter::new(self, RangeClamp::None)`); type annotation disambiguates. Alternative: `RangeClamp::Sequential` / `RangeClamp::Window` names the *mode*. Stylistic.
- `reader.rs:823-847` — The 12 repeated `.expect("... required by v1.0 schema")` calls are a duplication trigger at three-or-more, but the expand-once-per-column shape is clearer than a helper today. Will become structurally unreachable after B1; revisit if v1.x adds another required column. (Category: smells, downgraded.)
- `reader.rs:625-688` — `RecordsIter::next`'s `loop { ... }` is a four-phase state machine and reads cleanly as-is; decomposing into helpers would force shared `&mut self` recreating the same logic across method boundaries. Adding a four-bullet ASCII comment at the top of the loop body labelling the phases is the only suggestion. (Category: smells, low.)
- `reader.rs:738-744` — `BLOCK_HEADER_READ_CAP - buf.len()` is safe by construction (line 738 returns when `buf.len() >= cap`); belt-and-braces `BLOCK_HEADER_READ_CAP.saturating_sub(buf.len())` if defence-in-depth is desired.
- `reader.rs:344` — The tuple `<=` ordering is lexicographic and correct (verified by case analysis), but a one-line comment explicitly stating "lexicographic ordering required by the block-index sort key on `(chrom_id, first_pos)`" would help.
- `reader.rs:535` — Style nit: could be `let block = self.cur_block.as_ref()?;` returning `Option`. The panic-on-invariant is intentional; keep.
- `reader.rs:1043` — `io_err` is a closure factory; `io` is acronym-shaped (treated as a unit per Rust convention) and `err` reads as the obvious shorthand. Used 17 times in a 1900-line file; the abbreviation pays back. Acceptable.

## 7. Out of scope observations

- `block.rs:554` — `BlockHeader::active_chain_slots: Vec<u16>` (M12) is in `block.rs` (out of scope for this review) but is the wire-format side of the reader's `SlotId`-widening hazard. Fix lives in `block.rs`; the reader-side mirror (`reader.rs:441`) is in scope here.
- Several spec-check sites listed in `extras.md`'s Nits section (#1-#6 land in `header::parse_header_bytes` at `header.rs:489, 684, 701, 713`; manifest tag-ascending in `block.rs::validate_block_header_invariants`) — out of scope, included for cross-reference only. All 11 spec checks land somewhere reachable from `PspReader::new` + `RecordsIter::next`.
- The previous writer-side perf review's S1+S2 finding ([perf_psp_writer_2026-05-13.md](perf_psp_writer_2026-05-13.md)) — applied on the writer; reader-side mirror is M8 above.

## 8. Missing tests to add now

Sourced from the `reliability` subagent's challenge-tests pass. Per-test signatures use the `function_returns_expected_on_condition` form.

### 1. `materialise_next_record_returns_typed_error_when_n_alleles_sum_exceeds_n_total_alleles`

- **Input class:** Block whose decoded `n_alleles` column sums to more than `header.n_total_alleles`.
- **Bug caught:** B3 — index-out-of-bounds panic instead of typed error.
- **Body:** Build a 2-record fixture (`n_alleles=[1,1]`), then flip the second `n_alleles` varint byte from 0x01 to 0x03 (requires a test-only `replace_column_payload` helper that re-zstds and re-checksums the trailer). Assert `BlockHeaderInvariant { kind: NAllelesSumMismatch { .. } } | ColumnTruncated { .. }` (never a panic).

### 2. `load_next_block_rejects_non_empty_snapshot_at_chrom_change_in_region_mode`

- **Input class:** Multi-chromosome file with block N (first block of chrom 1) forged to declare a non-empty `active_chain_slots` snapshot.
- **Bug caught:** M7 — random-access mode silently trusts a non-empty cross-chromosome snapshot.
- **Body:** 10 records on chrom 0, 10 on chrom 1, locate block 1's snapshot slot_count byte (`block1_offset + 4` for the canonical fixture), patch from 0 to 1, splice `[0x07, 0x00]`. `region_records(1, 1, 10)` must yield `BlockHeaderInvariant { kind: NonEmptySnapshotAtChromStart { .. } }`.

### 3. `region_records_yields_exactly_one_record_for_start_equals_end_window`

- **Input class:** `region_records(chrom_id, p, p)` with exactly one record at `p`.
- **Bug caught:** Off-by-one in the inclusive-end window clamp (fencepost regression: `>=` or `<=` would return 0 or 2 records).
- **Body:** Three-block fixture, call `region_records(0, 150, 150)`, assert `got.len() == 1 && got[0].pos == 150`.

### 4. `region_records_yields_empty_iter_when_start_exceeds_end`

- **Input class:** `start > end` (caller bug).
- **Bug caught:** Drift from the silent-empty contract (regression to panic or infinite loop).
- **Body:** Three-block fixture, `region_records(0, 200, 100)`, assert `got.is_empty()`.

### 5. `region_records_yields_empty_iter_for_chrom_id_past_table_end`

- **Input class:** `chrom_id` greater than `header.chromosomes.len()`.
- **Bug caught:** A future change that uses `chrom_id` to index without bounds-checking would panic.
- **Body:** `region_records(5, 1, 1000)` on a 1-chrom fixture; assert empty.

### 6. `next_returns_none_after_first_error_is_yielded`

- **Input class:** Iterator that has already yielded `Some(Err(_))` once.
- **Bug caught:** Regression in the `poisoned` flag — a refactor that removes `self.poisoned = true` could let records flow after an error.
- **Body:** Reuse R10's snapshot-mutation fixture; drive `next()` to completion, count: exactly one `Err`, zero items after.

### 7. `read_block_header_returns_overflow_error_on_pathological_continuation_stream`

- **Input class:** Block-header stream of 70 KiB of `0x80` bytes (continuation set, payload 0).
- **Bug caught:** Regression in `BLOCK_HEADER_READ_CAP` check (flipping `>=` to `>`) or in the chunk-doubling arithmetic.
- **Body:** Direct call to `read_block_header` on `Cursor::new(vec![0x80u8; 70 * 1024])`; assert `BlockHeaderField { source: VarintError::Overflow }` or `Truncated` — or, with M3 applied, `BlockHeaderExceedsCap`.

### 8. `first_block_overlapping_finds_first_block_when_window_starts_before_it`

- **Input class:** Window `[0, 50]` against a fixture whose first block covers positions 1..=100.
- **Bug caught:** Off-by-one in `partition_point` candidate-window arithmetic (`pivot.saturating_sub(1)..(pivot + 1).min(index.len())`) missing the first block when `pivot == 0`.
- **Body:** Assert `region_records(0, 0, 50)` yields exactly 50 records (positions 1..=50).

### 9. `decode_block_payload_rejects_nan_in_allele_q_sum_log`

- **Input class:** F5 deferred case — block with NaN in `allele-q-sum-log`.
- **Bug caught:** M6 — semantic drift in `NonFiniteFloat`'s `entry` index, plus regression in the finite-float sweep.
- **Body:** Build a single-record single-allele file, replace the `allele-q-sum-log` column payload (8 bytes for one `f64`) with `f64::NAN.to_le_bytes()`, re-zstd, re-patch manifest `compressed_len`, re-checksum the trailer. Assert `PspReadError::NonFiniteFloat { column: "allele-q-sum-log", entry: 0, value: nan }`. Requires the test-only `replace_column_payload_f64` helper that the implementation report's §6 already recommends.

### 10. `pspreader_new_rejects_trailer_with_overflowing_arithmetic_without_panicking`

- **Input class:** Trailer whose `index_offset + index_byte_length` overflows `u64::MAX`.
- **Bug caught:** M16 — debug-build panic on the unchecked `+` in the error-construction expression.
- **Body:** Build a minimal `.psp` file, overwrite the trailer's `index_offset` to `u64::MAX` and `index_byte_length` to 1, leaving the magic intact. Assert `Err(IndexTrailingBytes { .. } | Io { .. } | BadTrailerMagic { .. })` without panicking.

## 9. What's good

- **Spec-check coverage is in real call paths, not in headers.** All 11 spec checks land somewhere reachable from `PspReader::new` + `RecordsIter::next` (verified by reading each site; see `extras.md` Nits). Check #11 in particular is implemented with the right asymmetry — sequential mode runs the equality check, region mode snapshot-inits.
- **`#[non_exhaustive]` is consistent across the six new variants** — inherited from `PspReadError`'s enum-level attribute at `errors.rs:229`; the new sub-enum `PhaseChainConsistencyKind` carries its own `#[non_exhaustive]` at `errors.rs:183`.
- **Block-header grow loop bounds the budget at 64 KiB** (`BLOCK_HEADER_READ_CAP`, line 56) with chunk-doubling capped at the remainder. Right shape for "varint streams have no delimiter, must grow on incomplete".
- **`partition_point` + candidate-window scan for `first_block_overlapping`** is the appropriate idiom for the half-open lookup; `binary_search_by` would not be cleaner because the predicate is not the overlap predicate. The comment block at lines 333-343 carries the weight of the trace.
- **`poisoned` flag is set at every error-producing path** in `materialise_next_record` and `load_next_block`, so once an iterator yields `Err`, subsequent `next()` calls return `None`. Sticky-poison is the right contract; only the test that pins it is missing (see §8 #6).

## 10. Commands to re-verify

Run inside the container (`./scripts/dev.sh ...`):

- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib --tests --all-features per_sample_caller::psp::reader`
- `cargo test --lib --tests --all-features` (full suite, for any cross-module breakage)
- `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --lib`

After the fix-application pass introduces a `criterion` bench (M13) and a `cargo-fuzz` target (M14):

- `cargo bench --bench psp_reader_perf` (informational; baseline for future regressions)
- `cargo fuzz run psp_reader_no_panic -- -max_total_time=600` (10-minute smoke run; longer on CI cron)

---

### Author response convention

Address each finding by its identifier (B1-B3, M1-M16, Mi1-Mi19) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer open questions §4.1-§4.4 first.
