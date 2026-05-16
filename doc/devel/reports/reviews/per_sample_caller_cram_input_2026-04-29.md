# Code Review: per_sample_caller / cram_input

**Date:** 2026-04-29
**Reviewer:** Claude (Opus 4.7, applying `ia/skills/code_review_skill.md`)
**Module:** `src/per_sample_caller/{mod, cram_input, errors, record_specs, cram_files}.rs`
**Status:** Approve-with-changes

---

## 1. Scope

- Reviewed: the diff in commit `94bbad4` ("implement CRAM input slice").
- In-scope files:
  - [src/per_sample_caller/mod.rs](../../src/per_sample_caller/mod.rs)
  - [src/per_sample_caller/cram_input.rs](../../src/per_sample_caller/cram_input.rs)
  - [src/per_sample_caller/errors.rs](../../src/per_sample_caller/errors.rs)
  - [src/per_sample_caller/record_specs.rs](../../src/per_sample_caller/record_specs.rs) (test fixture)
  - [src/per_sample_caller/cram_files.rs](../../src/per_sample_caller/cram_files.rs) (test fixture)
  - [Cargo.toml](../../Cargo.toml) and [src/lib.rs](../../src/lib.rs) (one line each)
- Deliberately out of scope:
  - Existing modules (`buffered_peekable.rs`, `gvcf_parser.rs`, `genotype_merging.rs`, etc.) — not touched by this commit.
  - `ia/reports/implementations/per_sample_caller_cram_input_2026-04-29.md` — narrative artefact, not code.
  - `CLAUDE.md` and `.claude/settings.json` — operational, gitignored, not code under review.

## 2. Verdict

**Approve-with-changes.** The slice meets its contract — coordinate-sorted N-CRAM merge with header validation and the per-read filter cascade — and the test mass is appropriate. There are five **Major** findings that should be fixed before this code is consumed by downstream slices: numeric truncation on positions/lengths, error variants reused outside their domain, and a silent fallback path in MD5 parsing. None of them are Blockers, but each is the kind of bug that surfaces only on inputs the test suite does not currently exercise.

## 3. Execution Status

| Command | Result |
|---|---|
| `cargo fmt --check` (host) | exit 0, no diff |
| `cargo clippy --all-targets --all-features -- -D warnings` (host) | clean for `per_sample_caller`; pre-existing warnings/errors in `decompression_pool.rs`, `genotype_merging.rs`, `vcf_writer.rs` remain (out of scope) |
| `./scripts/dev.sh cargo test --lib per_sample_caller` (container) | 22 passed, 0 failed |
| `./scripts/dev.sh cargo test --lib` (container) | 83 passed, 0 failed |
| `cargo doc --no-deps` | not run |
| `cargo audit` | not run (not installed in container) |

## 4. Open Questions and Assumptions

1. **MAPQ-missing semantics.** *Resolved 2026-04-29 by Jose:* reads with MAPQ unavailable (SAM 0xFF / noodles `None`) are treated as MAPQ 0 and therefore rejected under any non-zero `min_mapq`. The current code already does this via `unwrap_or(0)`. No behavioural change required; Mi6 becomes an action item to (a) add a doc-comment on `DEFAULT_MIN_MAPQ` recording the decision and rationale, (b) add a regression test pinning the behaviour. See Mi6 below.
2. **Fail-vs-tolerate on malformed CRAM `M5`.** *Resolved 2026-04-29 by Jose:* a malformed `M5` is a hard error. Per the project's "errors must not pass silently" principle, a present-but-invalid `M5` tag must surface a typed error rather than collapse to `None`. The fix proposed in M3 (add `CramInputError::MalformedMd5`, make `decode_md5_hex` fallible, thread the result up through `extract_header`) is the agreed implementation. See M3 below.
3. **Iterator behaviour after `Some(Err)`.** *Resolved 2026-04-29 by Jose:* fuse on first error. After `next()` returns `Some(Err)`, every subsequent call returns `None`. Implement `std::iter::FusedIterator`. Rationale: matches the spec's "halts Stage 1" language, keeps the iterator API safe under any caller pattern, and silently fixes the latent "OutOfOrderRead returned in a loop" defect (the offending head is currently not consumed before the error returns, so a naive retry sees the same error). See Mi3 below for the concrete action.
4. **Maximum supported contig length.** *Resolved 2026-04-29 by Jose:* widen positions and contig lengths to `u64` everywhere they cross the API surface. Rationale: plant and amphibian genomes have contigs > 4 Gb; `u64` is the minimum-friction way to support them and matches noodles' own internal width (`Position::get() -> usize`, `fai::Record::length() -> u64`). The cost of switching is small and isolated to this slice's surface area. See M1 below for the concrete change list.

## 5. Top 3 Priorities

1. **M1 — silent integer truncation** at multiple call sites converting `usize`/`u64` to `u32` via `as`. Largest correctness risk in the slice.
2. **M2 — `MultipleSampleNames` misused** for within-file SM disagreements (`path_a == path_b`), producing a misleading error message.
3. **M3 — silent `decode_md5_hex` tolerance** of malformed `M5` tags violates the project's "errors must not pass silently" principle.

## 6. Findings

### Major

#### M1: [src/per_sample_caller/cram_input.rs:329, 470, 987, 1045, 1056](../../src/per_sample_caller/cram_input.rs#L329) — silent `as u32` truncation on position and length values *(decision recorded 2026-04-29: widen to u64)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Widen positions and contig lengths to `u64`. Rationale: plant and amphibian genomes routinely exceed `u32::MAX` (~4.29 Gb); `u64` matches noodles' own internal width (`Position::get() -> usize`, `fai::Record::length() -> u64`); the change is local to this slice's surface. No fallible conversion is needed — `usize → u64` is lossless on every supported target.
- **Problem:**
  Five places convert from a wider integer to `u32` with the bare `as` cast:
  ```rust
  // src/per_sample_caller/cram_input.rs:329
  let length: u32 = usize::from(ref_seq_map.length()) as u32;
  // src/per_sample_caller/cram_input.rs:470
  let fai_length = fai_record.length() as u32;
  // src/per_sample_caller/cram_input.rs:987
  let pos = rb.alignment_start()?.get() as u32;
  // src/per_sample_caller/cram_input.rs:1045
  .get() as u32;
  // src/per_sample_caller/cram_input.rs:1056
  let mate_pos = rb.mate_alignment_start().map(|p| p.get() as u32);
  ```
  `Position::get()` returns `usize` (64-bit on the project's targets); `fai_record.length()` returns `u64`; `ref_seq_map.length()` is `NonZero<usize>`. If any of these exceeds `u32::MAX` (~4.29 Gb) the cast silently wraps, then validation against the FASTA — done in `u32` — passes against a wrong number, or a position is silently re-anchored to a wrong locus.

- **Why it matters:**
  Plant and amphibian genomes routinely have chromosomes near or above 4 Gb (lily, lungfish, etc.). For human genomes the values are within range, but a CRAM wrongly produced from a synthetic / concatenated reference would be silently corrupted instead of flagged. Per the project's design principles (`ia/specs/design_principles.md` §3), every out-of-range value is a hard error — silent wrap is the opposite.

- **Action — type changes (production code):**
  | File | Field | Before | After |
  |---|---|---|---|
  | `cram_input.rs` | `MappedRead::pos` | `u32` | `u64` |
  | `cram_input.rs` | `MappedRead::mate_pos` | `Option<u32>` | `Option<u64>` |
  | `cram_input.rs` | `ContigEntry::length` | `u32` | `u64` |
  | `cram_input.rs` | `WindowKey::pos` | `u32` | `u64` |
  | `cram_input.rs` | `HeadKey::pos` | `u32` | `u64` |
  | `cram_input.rs` | `PerFileOrder` (`Option<(usize, u32)>`) | `u32` | `u64` |
  | `cram_input.rs` | `window_anchor: Option<(usize, u32)>` | `u32` | `u64` |
  | `cram_input.rs` | `prev_per_file: Vec<Option<(usize, u32)>>` | `u32` | `u64` |
  | `errors.rs` | `OutOfOrderRead.{prev_pos, this_pos}` | already `u64` (delete `encode_order_key` per M5; widen the new structured `prev_pos` / `this_pos` introduced there to `u64` from the start) | `u64` |
  | `errors.rs` | `DuplicateReadAcrossFiles.pos` | `u32` | `u64` |

- **Action — call-site changes:**
  ```rust
  // cram_input.rs:329 — usize is already ≤ u64 on every supported target
  let length: u64 = usize::from(ref_seq_map.length()) as u64;

  // cram_input.rs:470 — fai_record.length() is already u64
  let fai_length: u64 = fai_record.length();

  // cram_input.rs:987 — Position::get() returns usize
  let pos = rb.alignment_start()?.get() as u64;

  // cram_input.rs:1045 (record_buf_to_mapped_read)
  let pos = rb.alignment_start().ok_or_else(/* … */)?.get() as u64;

  // cram_input.rs:1056
  let mate_pos = rb.mate_alignment_start().map(|p| p.get() as u64);
  ```
  Use of `as u64` is acceptable here because `usize → u64` is lossless on every supported target (32-bit and 64-bit). No fallible conversion or new error variant is required.

- **Action — test fixtures (these will recompile straightforwardly once the production types change):**
  | File | Field | Before | After |
  |---|---|---|---|
  | `record_specs.rs` | `RecordSpec::pos` | `u32` | `u64` |
  | `record_specs.rs` | `RecordSpec::mate_pos` | `Option<u32>` | `Option<u64>` |
  | `cram_files.rs` | `ContigSpec::length` | `u32` | `u64` |
  | `cram_files.rs` | `HeaderOverrides::length_overrides: Vec<(String, u32)>` | `u32` | `u64` |

  In `record_specs::record_spec`, change the `Position::new(spec.pos as usize)` calls to use `usize::try_from(spec.pos).ok().and_then(noodles_core::Position::new)` so a u64 position larger than `usize::MAX` (only possible on 32-bit) returns `None` instead of silently truncating. The `noodles_fasta::io::indexed_reader` API on 32-bit will reject such inputs anyway, but the conversion site should be explicit.

- **Action — assertions in existing tests:**
  Every numeric position literal in tests (e.g. `assert_eq!(out[0].pos, 100);`) becomes a `u64` literal automatically by type inference. No textual change required for the literals themselves; they are already untyped integer literals.

- **Action — regression test:**
  Add a fixture-driven test that pins the wider range:
  ```rust
  #[test]
  fn position_above_u32_max_round_trips_through_mapped_read() {
      let big_pos: u64 = (u32::MAX as u64) + 100;
      let mut rec = pass_record("R", 0, big_pos);
      // The record's CIGAR/SEQ length must still be > DEFAULT_MIN_READ_LENGTH
      // and the contig length in default_contigs() must be widened too.
      let cram = open_cram_from_records("a.cram", vec![record_spec(rec)]);
      let mut contigs = default_contigs();
      contigs.entries[0].length = big_pos + 1_000;
      let reader = CramMergedReader::from_open_crams(
          vec![cram], contigs, "sample".into(),
          CramMergedReaderConfig::default(),
      ).expect("reader");
      let mut iter = reader;
      let read = iter.next().expect("ok").expect("ok");
      assert_eq!(read.pos, big_pos);
  }
  ```
  Note: this test does *not* go through the noodles writer (Group A only), so we are not constrained by what a real CRAM writer accepts.

#### M2: [src/per_sample_caller/cram_input.rs:401-407](../../src/per_sample_caller/cram_input.rs#L401-L407) — `MultipleSampleNames` reused for within-file SM disagreement, leading to a self-referencing error *(decision recorded 2026-04-29: add dedicated variant)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Approved as proposed — add `CramInputError::MultipleSampleNamesInFile` for the within-file case and use it in `extract_single_sample_name`. The cross-file `MultipleSampleNames` variant stays as-is for its actual domain.
- **Assumptions:** None.
- **Problem:**
  When a single CRAM has multiple `@RG` entries with disagreeing `SM` values, `extract_single_sample_name` returns:
  ```rust
  return Err(CramInputError::MultipleSampleNames {
      path_a: path.to_path_buf(),
      sm_a: existing.clone(),
      path_b: path.to_path_buf(),  // same path
      sm_b: sm,
  });
  ```
  The `MultipleSampleNames` variant is documented (and intended in the spec) for the *cross-file* case. Here `path_a == path_b`. The user-facing message will say something like *"multiple sample names across CRAMs: '/a.cram' has SM:'foo', '/a.cram' has SM:'bar'"* — which is grammatically wrong and obscures the actual bug class.

- **Why it matters:**
  Two distinct failure modes (within-file SM disagreement vs cross-file SM disagreement) are surfaced through the same variant, which makes them indistinguishable to downstream code that wants to react to one but not the other, and produces a confusing error message at the CLI.

- **Suggested fix:**
  Add a dedicated variant and use it here:
  ```rust
  // errors.rs
  #[error(
      "multiple sample names within '{path}': @RG '{rg_a}' has SM:'{sm_a}', @RG '{rg_b}' has SM:'{sm_b}'"
  )]
  MultipleSampleNamesInFile {
      path: PathBuf,
      rg_a: String,
      sm_a: String,
      rg_b: String,
      sm_b: String,
  },
  ```
  Track each `rg_id` alongside its `sm` value in `extract_single_sample_name` and emit the new variant when they disagree.

#### M3: [src/per_sample_caller/cram_input.rs:361-372](../../src/per_sample_caller/cram_input.rs#L361-L372) — `decode_md5_hex` silently tolerates malformed `M5` *(decision recorded 2026-04-29: hard error)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** A present-but-malformed `M5` tag is a hard error. The fix below is approved as the implementation direction.
- **Assumptions:** None — the project's design principle 3 ("errors must not pass silently") is explicit on this case.
- **Problem:**
  ```rust
  fn decode_md5_hex(hex: &[u8]) -> Option<[u8; 16]> {
      if hex.len() != 32 {
          return None;
      }
      // …
  }
  ```
  If the CRAM `@SQ` carries an `M5` tag but it is not exactly 32 hex chars, `md5_from_reference_sequence` returns `None`, and the contig is treated as if it had no MD5 at all. `ContigEntry::eq`'s wildcard rule then matches it against any other MD5 silently. The spec explicitly says MD5 mismatches are a hard error and never a warning; corrupt-and-ignored is strictly worse.

- **Why it matters:**
  A subtly wrong CRAM (lowercase hex, extra space, half-byte truncated by an upstream tool) is accepted as if its `@SQ` line was clean. We then trust it for the downstream BAQ window and allele extraction.

- **Suggested fix:**
  Promote `decode_md5_hex` to fallible and surface a typed error:
  ```rust
  // errors.rs
  #[error("malformed @SQ M5 in '{path}' for contig '{contig}': {detail}")]
  MalformedMd5 { path: PathBuf, contig: String, detail: String },
  ```
  ```rust
  // cram_input.rs — replace md5_from_reference_sequence with
  fn md5_from_reference_sequence(
      path: &Path,
      contig_name: &str,
      ref_seq_map: &sam::header::record::value::Map<…>,
  ) -> Result<Option<[u8; 16]>, CramInputError> {
      use noodles_sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;
      let Some(raw) = ref_seq_map.other_fields().get(&MD5_CHECKSUM) else {
          return Ok(None);
      };
      decode_md5_hex(raw.as_ref()).map(Some).ok_or_else(|| {
          CramInputError::MalformedMd5 {
              path: path.to_path_buf(),
              contig: contig_name.into(),
              detail: format!(
                  "expected 32 hex chars, got {} bytes",
                  raw.as_ref().len()
              ),
          }
      })
  }
  ```
  Thread the resulting `Result` up through `extract_header`.

#### M4: [src/per_sample_caller/cram_input.rs:874-902](../../src/per_sample_caller/cram_input.rs#L874-L902) — `peek_head_keys` reports `MalformedRecord` for legitimate kept-unmapped reads *(decision recorded 2026-04-29: always drop unmapped, remove the toggle)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Always drop unmapped reads. Remove `drop_unmapped` from `CramMergedReaderConfig`. Rationale: an unmapped read contributes no allele evidence to any position; keeping it has no downstream consumer in this pipeline. The toggle was speculative flexibility that turned out to be a footgun.
- **Problem:**
  When `drop_unmapped = false`, an unmapped record (no `reference_sequence_id`, no `alignment_start`) survives `refill_heads`. `peek_head_keys` then sees a record where `head_key()` returns `None` and reports:
  ```rust
  Err(CramInputError::MalformedRecord {
      …,
      source: io::Error::new(
          io::ErrorKind::InvalidData,
          "record has no reference_sequence_id or alignment_start",
      ),
  })
  ```
  The current test suite avoids this path because every "unmapped" record in `a10` is built by `pass_record_for_b`-shaped helpers that *do* set `ref_id` and `pos` despite the `0x4` flag. In a real CRAM, an unmapped record routinely has `reference_sequence_id == None`.

- **Why it matters:**
  Real CRAMs aligned by bwa/STAR produce many unmapped records. With `drop_unmapped` removed and unmapped always filtered, this defect cannot fire — the `0x4` flag is enough to drop the record before `head_key()` is ever consulted.

- **Action — config and code changes:**
  ```rust
  // 1. cram_input.rs — drop the field from the struct.
  pub struct CramMergedReaderConfig {
      pub min_mapq: Option<u8>,
      pub min_read_length: Option<u32>,
      // pub drop_unmapped: bool,  ← remove
      pub drop_secondary: bool,
      pub drop_supplementary: bool,
      pub drop_qc_fail: bool,
      pub drop_duplicate: bool,
  }

  // 2. Update Default impl.
  impl Default for CramMergedReaderConfig {
      fn default() -> Self {
          Self {
              min_mapq: Some(DEFAULT_MIN_MAPQ),
              min_read_length: Some(DEFAULT_MIN_READ_LENGTH),
              drop_secondary: true,
              drop_supplementary: true,
              drop_qc_fail: true,
              drop_duplicate: true,
          }
      }
  }

  // 3. classify_pre_decode — unmapped is now unconditional.
  // Slot it in at the same hit-rate-ordered position (still after
  // duplicate/MAPQ/supplementary/secondary). The compile-time
  // guarantee "unmapped is always dropped" lets us remove the
  // defensive arm in peek_head_keys (next bullet).
  if (flag & FLAG_UNMAPPED) != 0 {
      return PreDecodeOutcome::Drop(FilterBucket::Unmapped);
  }

  // 4. peek_head_keys — the "head_key returned None" arm becomes
  // unreachable in the normal path because every record reaching
  // peek_head_keys has been filtered by classify_pre_decode and
  // therefore has a non-None reference_sequence_id and
  // alignment_start. Either:
  //   (a) keep the defensive arm but log/comment it as "should be
  //       unreachable; surfaces malformed CRAMs from upstream tools",
  //       OR
  //   (b) replace it with debug_assert! and an io::Error pointing at
  //       a malformed input.
  // Recommend (a) — defence in depth against malformed CRAM bytes.
  ```
- **Action — test changes:**
  - `a9_each_flag_drop_one_at_a_time` ([cram_input.rs:1438](../../src/per_sample_caller/cram_input.rs#L1438)) — remove the `unmapped` row from `cases` (it has no toggle to flip any more); the "all defaults" sub-assertion still expects `counts.unmapped == 1` and is correct unchanged.
  - `a10_all_flag_drops_disabled_passes_everything` ([cram_input.rs:1532](../../src/per_sample_caller/cram_input.rs#L1532)) — must be updated. With `drop_unmapped` gone, the unmapped record in `six_flagged_records()` is dropped unconditionally:
    ```rust
    // Now expects 5 records (not 6), and counts.unmapped == 1.
    assert_eq!(out.len(), 5);
    let mut expected = FilterCounts::default();
    expected.unmapped = 1;
    assert_eq!(counts, expected);
    ```
    Rename the test for accuracy (e.g. `all_optional_flag_drops_disabled_keeps_everything_except_unmapped`), or restructure `six_flagged_records()` to a `five_optional_flag_records()` helper that returns the five records *whose drop is still configurable* and keep `a10`'s assertion as-is.
- **Action — regression test pinning the new contract:**
  ```rust
  #[test]
  fn truly_unmapped_record_with_no_ref_id_is_filtered_not_errored() {
      // Build a record with flag=FLAG_UNMAPPED, ref_id=None,
      // alignment_start=None — i.e. a real-CRAM unmapped record.
      let unmapped_realistic = RecordSpec {
          qname: "U".into(),
          flag: FLAG_UNMAPPED,
          ref_id: 0,           // record_spec sets it; the test's point
          pos: 0,              // is that mapping_quality / pos / etc.
          mapq: 0,             // are not what the filter relies on.
          cigar_ops: vec![],
          seq: vec![],
          qual: vec![],
          mate_ref_id: None,
          mate_pos: None,
      };
      // (record_spec already maps pos=0 to "no alignment_start" via
      // its `if spec.pos > 0` guard, and ref_id=0 with flag=0x4 plus
      // empty SEQ/CIGAR is the closest the helper comes to a
      // real-CRAM unmapped record. If a tighter fixture is needed,
      // extend record_spec to take Option<usize> for ref_id.)
      let cram = open_cram_from_records("a.cram", vec![record_spec(unmapped_realistic)]);
      let mut reader = CramMergedReader::from_open_crams(
          vec![cram],
          default_contigs(),
          "sample".into(),
          CramMergedReaderConfig::default(),
      ).expect("reader");
      assert!(reader.next().is_none(), "unmapped record must be silently dropped");
      assert_eq!(reader.filter_counts().unmapped, 1);
  }
  ```
  If `record_spec` does not currently allow constructing a record with `ref_id == None`, extend it as part of this fix so the test exercises the actual real-CRAM shape rather than a synthetic stand-in.

#### M5: [src/per_sample_caller/cram_input.rs:751-757](../../src/per_sample_caller/cram_input.rs#L751-L757) — `OutOfOrderRead` carries packed `(ref_id, pos)` keys in fields named `prev_pos` and `this_pos` *(decision recorded 2026-04-29: approved as proposed)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Approved as proposed — replace the packed-u64 `prev_pos`/`this_pos` with structured `prev_ref_id, prev_pos, this_ref_id, this_pos`, drop `encode_order_key`. The `prev_pos`/`this_pos` fields are already `u64` per M1's widening decision, so M5's new shape lands naturally on top.
- **Assumptions:** None.
- **Problem:**
  ```rust
  return Some(Err(CramInputError::OutOfOrderRead {
      path: self.paths[chosen_idx].clone(),
      qname: …,
      prev_pos: encode_order_key(prev_ref, prev_pos),
      this_pos: encode_order_key(head.ref_id, head.pos),
  }));
  ```
  with `encode_order_key` producing `(ref_id << 32) | pos`. The variant declares the fields as `prev_pos: u64` and `this_pos: u64`, but the value carries *both* a chromosome index and a position in one number. The error message renders something like *"out-of-order read … at 4294967396 regresses from 4294967396"* — uninterpretable, and silently identical when the two records are on different chromosomes whose positions happen to differ by less than `2^32`.

- **Why it matters:**
  The error is one of the spec-mandated hard errors and is the user's only signal that a CRAM is mis-sorted. A garbled message slows diagnosis.

- **Suggested fix:**
  Replace the packed key with the actual `(ref_id, pos)` pair:
  ```rust
  // errors.rs
  #[error(
      "out-of-order read in '{path}': QNAME '{qname}' at \
       (ref_id={this_ref_id}, pos={this_pos}) regresses from \
       (ref_id={prev_ref_id}, pos={prev_pos})"
  )]
  OutOfOrderRead {
      path: PathBuf,
      qname: String,
      prev_ref_id: usize,
      prev_pos: u64,
      this_ref_id: usize,
      this_pos: u64,
  },
  ```
  and remove `encode_order_key`. Update `a4_out_of_order_within_a_single_stream` to assert on the structured fields rather than the packed value.

### Minor

#### Mi1: [src/per_sample_caller/cram_input.rs:543-551](../../src/per_sample_caller/cram_input.rs#L543-L551) — empty-input check returns `CramInputError::Io` *(decision recorded 2026-04-29: approved as proposed)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Approved as proposed — add `CramInputError::NoInputs` and use it in `CramMergedReader::new` instead of the synthesised `Io` error.
- **Problem:** `crams.is_empty()` is a programming error from the caller, not an I/O failure. Reusing the `Io` variant with an empty `PathBuf` and a synthesized `io::Error` is misleading; a reader of the error sees "I/O error on '': at least one CRAM input is required" with no path that exists.
- **Why it matters:** Error variants are part of the typed-error contract — code that matches on `CramInputError::Io { path, .. }` would catch this case alongside genuine I/O failures.
- **Suggested fix:** Add a dedicated `NoInputs` variant or simply make the constructor reject empty input via type (e.g. require at least one path through a `&[PathBuf; N]` of compile-time length, or accept `(PathBuf, Vec<PathBuf>)`). At minimum:
  ```rust
  #[error("at least one CRAM input is required")]
  NoInputs,
  ```

#### Mi2: [src/per_sample_caller/cram_input.rs:621, 632](../../src/per_sample_caller/cram_input.rs#L621-L632) — defensive `.unwrap_or_default()` on `canonical_path` masks an invariant *(decision recorded 2026-04-29: approved as proposed)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Approved — replace the two `.unwrap_or_default()` calls with `.expect("canonical_path is set when canonical_contigs is Some")`, matching the style at lines 657-660. The fixing agent may also pull `canonical_contigs`, `canonical_sample`, and `canonical_path` into a single struct that's only `Option<>` until the first iteration completes — that is a stronger restructure but optional; the simple `.expect` is sufficient.
- **Problem:**
  ```rust
  reference_path: canonical_path.clone().unwrap_or_default(),
  ```
  By construction, `canonical_path` is `Some` whenever `canonical_contigs` or `canonical_sample` is `Some`, because they are all set in the same `None` branch of the iteration. Lines 657-660 already use `.expect(...)` for the same invariant. The two error paths above silently substitute `PathBuf::new()` if the invariant breaks instead of panicking, which would hide the bug.
- **Suggested fix:** Replace with the same `expect` text:
  ```rust
  reference_path: canonical_path
      .clone()
      .expect("canonical_path is set when canonical_contigs is Some"),
  ```
  Or refactor so `canonical_contigs`, `canonical_sample`, `canonical_path` are pulled from the first iteration into a single struct, eliminating the option entirely from second-iteration onwards.

#### Mi3: [src/per_sample_caller/cram_input.rs:725-823](../../src/per_sample_caller/cram_input.rs#L725-L823) — iterator behaviour after `Some(Err)` is unspecified and can loop the same error *(decision recorded 2026-04-29: fuse on first error)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Option A — fuse on first error. After `next()` returns `Some(Err)`, every subsequent call returns `None`. Implement `std::iter::FusedIterator`. Rationale: matches the spec's "halts Stage 1" language, makes the iterator safe under any caller pattern (including `.collect::<Result<Vec<_>, _>>()` and `for r in reader { … }` without explicit break-on-error), and silently fixes the latent infinite-loop defect described below.
- **Problem:** When `next()` returns `Some(Err(OutOfOrderRead{...}))` or `Some(Err(DuplicateReadAcrossFiles{...}))`, the offending head is *not* consumed — the error is returned before the `peekers[chosen_idx].next()` call. A caller that naively retries (or uses `for r in reader { … }` without breaking on `Err`) gets the same error in a tight loop. Even for variants where the head *is* consumed (`MalformedRecord`), a re-driven iterator's behaviour is undefined and brittle.
- **Why it matters:** The Rust iterator convention does not enforce fusing on error; the spec's "halts Stage 1" language does. Without code-level enforcement, the contract relies on caller discipline that is easy to break.
- **Action — implement fuse-on-error:**
  ```rust
  // 1. Add a field to CramMergedReader (cram_input.rs).
  pub struct CramMergedReader {
      // … existing fields …
      /// Set to true the first time `next()` yields `Some(Err)`. Once
      /// set, every subsequent call returns `None`. See
      /// `ia/reviews/per_sample_caller_cram_input_2026-04-29.md` (Mi3).
      fused: bool,
  }

  // 2. Initialise in both constructors.
  fused: false,

  // 3. Guard at the top of next() and flip on every Err return.
  fn next(&mut self) -> Option<Self::Item> {
      if self.fused {
          return None;
      }
      // … existing body, but every `return Some(Err(e))` becomes:
      self.fused = true;
      return Some(Err(e));
      // (or factor the assignment into a small helper)
  }

  // 4. Mark the iterator as fused for combinator correctness.
  impl std::iter::FusedIterator for CramMergedReader {}
  ```
  Helper to keep the body tidy:
  ```rust
  impl CramMergedReader {
      fn fail(&mut self, e: CramInputError) -> Option<<Self as Iterator>::Item> {
          self.fused = true;
          Some(Err(e))
      }
  }
  // then: return self.fail(CramInputError::OutOfOrderRead { … });
  ```
- **Action — pin behaviour with a regression test:**
  ```rust
  #[test]
  fn iterator_returns_none_after_first_error() {
      // Same fixture as a4: pos 200 then pos 100 in a single stream.
      let cram_a = open_cram_from_records(
          "a.cram",
          vec![
              record_spec(pass_record("R1", 0, 200)),
              record_spec(pass_record("R2", 0, 100)),
          ],
      );
      let mut reader = CramMergedReader::from_open_crams(
          vec![cram_a],
          default_contigs(),
          "sample".into(),
          CramMergedReaderConfig::default(),
      )
      .expect("reader");
      assert!(reader.next().expect("ok item").is_ok());      // pos=200
      assert!(reader.next().expect("err item").is_err());    // OutOfOrderRead
      assert!(reader.next().is_none(), "iterator must fuse after Err");
      assert!(reader.next().is_none(), "fuse is sticky");
  }
  ```
- **Action — update the doc comment on `impl Iterator`:**
  ```rust
  /// Iterates over surviving reads in coordinate order.
  ///
  /// **Fuse-on-error semantics.** Once `next()` returns `Some(Err(_))`,
  /// every subsequent call returns `None`. This matches the spec's
  /// "halts Stage 1" requirement (`ia/specs/per_sample_caller.md`
  /// §"Errors") and lets callers use the iterator with any consumer
  /// (`for`, `collect`, `try_fold`) without separate "stop on first
  /// error" bookkeeping.
  impl Iterator for CramMergedReader { … }
  ```

#### Mi4: [src/per_sample_caller/cram_input.rs:840-902](../../src/per_sample_caller/cram_input.rs#L840-L902) — `MalformedRecord` errors with empty `qname` context *(decision recorded 2026-04-29: approved as proposed)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Approved — change `MalformedRecord.qname` to `Option<String>` and adjust the `#[error]` format string to omit the qname clause when `None`. The fixing agent should pick whichever of the two equivalent shapes (`Option<String>` or "empty string means absent" with a custom `Display`) is most idiomatic; the user-visible behaviour — no more `qname=''` clauses — is what matters.
- **Problem:** Both `refill_heads` (l. 848-852) and the `peek_head_keys` Err arm (l. 897-901) emit `qname: String::new()`. The error message will say *"malformed record in '/foo.cram' (qname=''): …"* — the empty qname does not localise the failure within the file.
- **Suggested fix:** When `peek()` returns `Err`, we have lost the record but we can at least try to read the previous accepted qname for that file (via `prev_per_file[idx]`'s position) or simply omit the qname field from the message when empty:
  ```rust
  #[error("malformed record in '{path}'{qname_part}: {source}")]
  MalformedRecord {
      path: PathBuf,
      qname: String,
      #[source]
      source: std::io::Error,
  }
  // with a Display impl that omits the qname=… clause when empty,
  // or change the error variant to take Option<String>.
  ```
  Lower-effort: store the previous accepted qname per file alongside `prev_per_file` and reuse it as the locator.

#### Mi5: [src/per_sample_caller/cram_input.rs:806-813](../../src/per_sample_caller/cram_input.rs#L806-L813) — `min_read_length` filter pays for full record decode before rejecting *(decision recorded 2026-04-29: approved as proposed)*

- **Confidence:** Medium
- **Decision (2026-04-29, Jose):** Approved — read SEQ length from the `RecordBuf` and apply `min_read_length` *before* allocating the `MappedRead` (clone of qname/CIGAR/seq/qual). Keep this check after the order and dedup checks so the duplicate-detection window stays consistent.
- **Problem:** `record_buf_to_mapped_read` (l. 1027-1070) clones `qname`, walks and converts the CIGAR vector, uppercases SEQ, copies QUAL, and only then returns. The post-decode min-read-length check at l. 808-813 then drops the result. For datasets where many reads are short (RNA-seq UMI-tagged reads, ancient DNA), this is wasted work.
- **Suggested fix:** Read SEQ length from the `RecordBuf` before allocating the `MappedRead`:
  ```rust
  if let Some(min) = self.config.min_read_length
      && (record.sequence().as_ref().len() as u32) < min
  {
      self.filter_counts.too_short += 1;
      let _ = self.peekers[chosen_idx].next();
      continue;
  }
  let mapped = match record_buf_to_mapped_read(&record, chosen_idx) { … };
  ```
  Note: the order check, dup check, and consume-the-head still happen first, so the duplicate window stays consistent.

#### Mi6: [src/per_sample_caller/cram_input.rs:1046](../../src/per_sample_caller/cram_input.rs#L1046) — MAPQ-missing collapses to 0 *(decision recorded; documentation + regression test required)*

- **Confidence:** High
- **Decision (2026-04-29, Jose):** Reads with MAPQ unavailable (SAM 0xFF / noodles `mapping_quality()` returning `None`) are treated as MAPQ 0 and rejected under any non-zero `min_mapq`. Rationale: at the project's low coverage targets (2-10×), an "unknown" placement is not trustworthy enough to admit; matches the bcftools/freebayes convention.
- **Current code:** Already correct — `rb.mapping_quality().map(u8::from).unwrap_or(0)` at both [cram_input.rs:1046](../../src/per_sample_caller/cram_input.rs#L1046) (`record_buf_to_mapped_read`) and [cram_input.rs:952](../../src/per_sample_caller/cram_input.rs#L952) (`classify_pre_decode`).
- **Action 1 — document the decision so a future refactor doesn't silently flip it:**
  ```rust
  /// Reads with MAPQ strictly below this are dropped. Matches bcftools'
  /// default and the spec recommendation in
  /// `ia/specs/per_sample_caller.md` §"Read filters".
  ///
  /// Reads with MAPQ unavailable (SAM 0xFF / noodles `mapping_quality()`
  /// returning `None`) are treated as MAPQ 0 and therefore rejected
  /// under any non-zero minimum. Decision recorded 2026-04-29: at the
  /// project's 2-10× coverage targets an "unknown" placement is not
  /// trustworthy enough to admit; matches the bcftools/freebayes
  /// convention. See `ia/reviews/per_sample_caller_cram_input_2026-04-29.md`
  /// (Mi6) for the full rationale.
  pub const DEFAULT_MIN_MAPQ: u8 = 20;
  ```
- **Action 2 — pin the behaviour with a regression test.** A naive refactor that swaps `unwrap_or(0)` for `unwrap_or(u8::MAX)` would let unknown-MAPQ reads through silently; the test below would fail loudly.
  ```rust
  #[test]
  fn missing_mapq_is_treated_as_zero_and_filtered_by_default() {
      // record_spec maps mapq=0 to "no MappingQuality set" because of
      // the `if spec.mapq > 0` guard, so this exercises the
      // mapping_quality()==None path.
      let mut rec = pass_record("R", 0, 100);
      rec.mapq = 0;
      let cram = open_cram_from_records("a.cram", vec![record_spec(rec)]);
      let reader = CramMergedReader::from_open_crams(
          vec![cram],
          default_contigs(),
          "sample".into(),
          CramMergedReaderConfig::default(),
      )
      .expect("reader");
      let (out, counts, err) = run_to_completion(reader);
      assert!(err.is_none());
      assert!(out.is_empty(), "MAPQ-missing read should be dropped");
      assert_eq!(counts.low_mapq, 1);
  }
  ```

### Nits

*(decision recorded 2026-04-29: approved as proposed except where noted)*

- **`OpenCram::path_for_errors`** — verbose. `OpenCram::path` is unambiguous given the surrounding code. No behavioural change. **Approved.**
- **`window: VecDeque<WindowEntry>`** ([cram_input.rs:530](../../src/per_sample_caller/cram_input.rs#L530)) — only `push_back`, linear `find`, and `clear()` are used. `Vec<WindowEntry>` is simpler and has the same complexity. Replace. **Approved.**
- **WindowKey clone** ([cram_input.rs:765](../../src/per_sample_caller/cram_input.rs#L765)) — `head.qname.clone()` could be a move if `head` is no longer needed afterwards. Restructure or accept the small clone (qname is short). **Approved — fixing agent's discretion which way.**
- **`encode_order_key`** ([cram_input.rs:1001-1005](../../src/per_sample_caller/cram_input.rs#L1001-L1005)) — gone after M5 fix. **Approved (subsumed by M5).**
- **Implementation file size** — `cram_input.rs` is 2117 lines; ~1100 are tests. Splitting `head_key` / `record_buf_to_mapped_read` / `cigar_to_ops` into a `conversion.rs` and `extract_header` / `validate_fasta_agreement` into a `header_validation.rs` would each be ~150 lines. **Won't fix now (2026-04-29, Jose):** keep `cram_input.rs` as a single file for this slice. The split may make sense once later slices (BAQ, pileup walker) reveal which functions actually want to be reused; doing it pre-emptively is premature. The fixing agent must NOT split this file as part of the review fixes.
- **Unused conversion `Container::default()`** at [cram_input.rs:648](../../src/per_sample_caller/cram_input.rs#L648) — fine, but adding `Container::new()` would be more idiomatic if noodles offered it. Pure aesthetics. **Approved (no concrete change required; revisit on noodles upgrade).**

## 7. Out of Scope Observations

These are pre-existing issues in untouched files. They are noted because they were observed while running validation commands; **none block this PR**.

- `src/decompression_pool.rs:403` — `io::Error::new(io::ErrorKind::Other, "read failed")` flagged by `clippy::io_other_error`. Suggest `io::Error::other(...)` in a follow-up.
- `src/decompression_pool.rs:424` — `reader.read(&mut buf).unwrap()` is flagged by `clippy::unused_io_amount`. The test happens to use a `Cursor` that reads to completion, so it works in practice; should be `read_exact` or handle short reads.
- Multiple modules carry `clippy::repeat_take`, `clippy::collapsible_if`, `clippy::needless_range_loop` warnings. Aggregate them into a one-shot lint-cleanup PR.

## 8. Missing Tests to Add Now

Each test addresses a specific bug class above. Each name follows
`function_returns_expected_on_condition`.

### Against `decode_md5_hex` / `extract_header` (M3)

```rust
#[test]
fn extract_header_rejects_malformed_md5() {
    // Build a sam::Header whose @SQ M5 is "not_hex" (7 bytes, not 32).
    // Assert extract_header returns CramInputError::MalformedMd5.
}

#[test]
fn extract_header_rejects_short_md5_hex() {
    // 30 hex chars instead of 32. Same expectation.
}
```

### Against `extract_single_sample_name` (M2)

```rust
#[test]
fn within_file_sm_disagreement_uses_dedicated_variant() {
    // Two @RG entries in the same file: SM:foo, SM:bar.
    // Assert CramInputError::MultipleSampleNamesInFile, and that
    // path_a == path_b is NOT reported (the variant is gone).
}
```

### Against `OutOfOrderRead` (M5)

```rust
#[test]
fn out_of_order_error_reports_structured_ref_id_and_pos() {
    // Build the same fixture as a4 but assert the structured fields
    // prev_ref_id, prev_pos, this_ref_id, this_pos directly.
}
```

### Against the merge iterator (Mi3)

```rust
#[test]
fn iterator_returns_none_after_first_error() {
    // Drive the same fixture as a4 (out-of-order). After the first
    // Err, every subsequent next() returns None.
}
```

### Against `drop_unmapped: false` and unmapped real records (M4)

```rust
#[test]
fn truly_unmapped_record_with_no_ref_id_does_not_error() {
    // RecordSpec with flag=FLAG_UNMAPPED AND no reference_sequence_id /
    // no alignment_start. With drop_unmapped=false:
    // - Either yielded as a MappedRead (if option 2 is taken)
    // - Or filtered with FilterCounts.unmapped += 1 (if option 1 is taken)
    // Either way, NO MalformedRecord.
}
```

### Against numeric truncation (M1)

If the project picks the "widen to u64" route:

```rust
#[test]
fn contig_length_above_u32_max_is_preserved() {
    // FASTA with a synthetic contig of length u32::MAX as u64 + 1.
    // Assert ContigList.entries[0].length is the correct u64.
}
```

If it picks the "checked u32" route:

```rust
#[test]
fn contig_length_above_u32_max_returns_typed_error() {
    // FASTA same as above. Assert CramInputError::ContigTooLong.
}
```

### CIGAR fidelity beyond `Match`

```rust
#[test]
fn cigar_to_ops_roundtrips_every_kind() {
    // Build a RecordSpec carrying a CIGAR with M, I, D, N, S, H, P, =, X.
    // Assert MappedRead.cigar matches our CigarOp enum entries 1:1.
}
```

## 9. What's Good

- **Two-constructor split** ([cram_input.rs:537-706](../../src/per_sample_caller/cram_input.rs#L537-L706)) cleanly separates I/O from logic; tests exercise the merge with zero filesystem traffic. This is a model the rest of the pipeline should reuse.
- **`OwnedCramRecords`** ([cram_input.rs:247-299](../../src/per_sample_caller/cram_input.rs#L247-L299)) sidesteps a real lifetime trap (noodles `Records` borrows reader + header) without unsafe — by reproducing the slice loop on owned state. Worth pointing future readers at.
- **`ContigEntry::PartialEq`** ([cram_input.rs:97-109](../../src/per_sample_caller/cram_input.rs#L97-L109)) treats absent `M5` as a wildcard rather than as "different from any value". Compact, well-tested, and the right semantics for cross-file CRAM comparison.
- **Hit-rate ordered filter cascade** ([cram_input.rs:936-975](../../src/per_sample_caller/cram_input.rs#L936-L975)) is documented in code with the per-filter expected hit rate, making the order auditable. Good template for similar filter pipelines.
- **`record_specs.rs` builders** ([record_specs.rs:23-76](../../src/per_sample_caller/record_specs.rs#L23-L76)) keep Group A tests under ten lines apiece. Reusable verbatim by every later slice that needs synthetic `RecordBuf`s.

## 10. Commands to Re-verify

- `cargo fmt --check` — must remain clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — must remain clean for `per_sample_caller`.
- `./scripts/dev.sh cargo test --lib per_sample_caller` — should grow from 22 to ~28 passing tests after the M-finding fixes (one new test per fix, plus M1/M2/M3/M5 each retire / replace one assertion).
- `./scripts/dev.sh cargo test --lib` — still 83+ passing total, no regressions.

## 11. Author Response Template

Address each finding by its code (M1, Mi3, etc.) with one of:

- `fixed in <commit>` — link to follow-up commit.
- `disputed because …` — argue the finding away with evidence.
- `deferred to <issue/plan>` — point at where the work is tracked.
- `won't fix because …` — justify keeping current behaviour.

Resolve Open Questions 1-4 before responding to individual findings, since several findings are downstream of those answers.
