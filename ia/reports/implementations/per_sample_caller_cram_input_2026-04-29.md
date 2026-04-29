# Per-sample caller — CRAM input slice — implementation report

**Date:** 2026-04-29.
**Plan:** [`ia/feature_implementation_plans/per_sample_caller_cram_input.md`](../../feature_implementation_plans/per_sample_caller_cram_input.md).
**Spec:** [`ia/specs/per_sample_caller.md`](../../specs/per_sample_caller.md) §"Inputs", §"Multi-CRAM ingestion", §"Read filters", §"Errors".
**Architecture:** [`ia/specs/calling_pipeline_architecture.md`](../../specs/calling_pipeline_architecture.md).
**Design principles:** [`ia/specs/design_principles.md`](../../specs/design_principles.md).

## Summary

First slice of Stage 1 of the multi-sample calling pipeline:
turning N coordinate-sorted CRAMs + an indexed reference FASTA into a
single coordinate-sorted stream of `MappedRead`s, with the cheap
per-read filter cascade applied along the way. No BAQ, no pileup
walker, no allele extraction, no `.psf` writing — those are later
slices.

The slice ships `CramMergedReader` with two constructors that match
the **side-effects-at-the-edges** principle: `new` opens files and
validates headers; `from_open_crams` is the I/O-free merge/filter core
tests drive directly with injected `Vec<RecordBuf>` streams.

## Code added

```
src/per_sample_caller/
    mod.rs                    pub mod re-exports + (test) helper modules
    cram_input.rs             header validation, peek-and-scan merge,
                              filter cascade, MappedRead + CigarOp,
                              CramMergedReader{,Config}, FilterCounts,
                              ContigList, OwnedCramRecords (private)
    errors.rs                 CramInputError (thiserror enum)
    record_specs.rs           Group A test fixture (no I/O)
    cram_files.rs             Group B test fixture (real CRAM/FASTA)
```

[`src/lib.rs`](../../../src/lib.rs) gets one new line: `pub mod per_sample_caller;`.

`Cargo.toml` gains `noodles-cram = 0.92`, `noodles-sam = 0.84`,
`noodles-fasta = 0.60`, `noodles-core = 0.19`, `anyhow = 1.0` (kept
for future Stage 1 orchestration in `mod.rs`), and `bstr = 1` as a
dev dependency for fixture record assembly.

## Highlights vs. plan

- **Owned record iteration over noodles.** `noodles_cram::io::reader::Records`
  borrows from `Reader` *and* `Header` — it cannot be packed into a
  `Box<dyn Iterator + 'static + Send>`. Implemented `OwnedCramRecords`
  (private to `cram_input.rs`) which replicates the noodles records
  loop manually: it owns the `Reader<File>`, the `sam::Header`, and a
  cloned `fasta::Repository`, drives `read_container` →
  `Container::compression_header` → `Container::slices` →
  `Slice::decode_blocks` → `Slice::records` → `RecordBuf::try_from_alignment_record`,
  and yields owned `RecordBuf`s. The merge then wraps it in
  `BufferedPeekable` exactly as the plan specifies.
- **Filter cascade as a free function.** Originally a method on
  `CramMergedReader`; the borrow checker rejected `self.classify_pre_decode(rb)`
  while a `peek()` borrow was live on `self.peekers`. Refactor to
  free function `classify_pre_decode(&CramMergedReaderConfig, &RecordBuf)`
  (config is `Copy`) — same behaviour, no false-positive borrow
  conflict, and it makes the cascade trivially unit-testable in
  isolation if the need arises later.
- **Trust .fai length, accept the spec's MD5 trust.** The CRAM `@SQ`
  validates against the canonical `ContigList` built from the *first*
  CRAM seen; FASTA agreement is `(name, length)` only, against the
  `.fai`. We do not recompute MD5 from FASTA bytes (per the user
  decision recorded in the plan). MD5 disagreement *between* CRAMs
  is still surfaced.
- **`ContigEntry` MD5 wildcard equality.** Custom `PartialEq` so that
  a CRAM that omits `M5` is not in conflict with a CRAM that carries
  one — only `Some(a) != Some(b)` is a mismatch. Tested in P1.
- **Pre-flight validation is strict.** Every failure mode listed in
  the spec's §"Errors" maps to a concrete `CramInputError` variant
  carrying the offending file path / contig / qname / position. No
  `Other(String)` catch-all.
- **Hit-rate-ordered cascade.** Implemented exactly per the plan:
  duplicate → MAPQ → supplementary → secondary → unmapped → qc_fail,
  then post-decode `min_read_length`. Each filter gated by its
  `drop_*` config bool; `min_mapq` and `min_read_length` are
  `Option`s where `None` means "no minimum". Verified by tests A7,
  A8, A9, A10, A11, A12.
- **No magic numbers.** `DEFAULT_MIN_MAPQ`, `DEFAULT_MIN_READ_LENGTH`,
  `PER_PEEKER_BUFFER_SIZE`, plus the SAM/BAM flag bits (`FLAG_UNMAPPED`,
  etc.) live as named constants with doc comments. Tests assert the
  documented defaults so a silent change is loudly caught.

## Departures from plan

- **Empty-records writes for header-mismatch tests.** The plan called
  for synthetic CRAMs with records throughout. In practice
  `noodles_cram::io::Writer` cross-checks each record's `@SQ` name
  against the FASTA repository, so a CRAM whose `@SQ` is renamed
  fails at write time, not at our reader's header check. To exercise
  the cross-CRAM `ContigListMismatch` case for *name*/*length*/*md5*
  disagreements we write CRAMs with empty record vectors — the @SQ
  list still propagates through the file header and our reader's
  pre-flight validation runs against it. Documented in the comments
  on `b5_contig_list_mismatch_across_crams`. Rejected the alternative
  of bypassing the noodles writer with hand-crafted CRAM bytes —
  far more brittle than building empty-but-valid CRAMs.
- **B5 sub-case 1 reframed as order disagreement.** The plan suggests
  testing "name disagreement" with renamed contigs. We instead build
  a 2-contig FASTA and write two CRAMs with the contigs listed in
  *opposite order*; both match the FASTA, so the FASTA-agreement
  check passes and the cross-CRAM check fires with detail
  `name disagreement at index 0`. This is a closer model of a real
  failure mode: lanes/shards typically don't rename contigs but can
  emit them in different orders.
- **No separate `read_filter.rs` module.** Per the plan and Tradeoff
  in the plan, the filter cascade lives in `cram_input.rs`. The spec
  module-layout sketch's `read_filter.rs` line is intentionally
  superseded.

## Tests added

22 tests total, all passing:

### Pure-type (1)

- `p1_contig_list_md5_wildcard_equality` — `ContigEntry` `PartialEq`
  treats a `None` MD5 as a wildcard against `Some`.

### Group A — via `from_open_crams` (14)

- `a1_single_stream_pass_through` — three records, fields preserved.
- `a2_multi_stream_merge_order` — interleaved positions across two
  streams, sources alternate correctly.
- `a3_tiebreaker_on_equal_coordinates` — file-index breaks ties.
- `a4_out_of_order_within_a_single_stream` — `OutOfOrderRead` on
  regression.
- `a5_duplicate_read_across_streams` — `DuplicateReadAcrossFiles` on
  same `(qname, flag, ref_id, pos)` in two files.
- `a6_duplicate_window_clears_on_advance` — same QNAME at different
  positions is not a duplicate.
- `a7_min_mapq_filter` + asserts `DEFAULT_MIN_MAPQ == 20`.
- `a8_min_mapq_none_disables` — `Option::None` keeps everything.
- `a9_each_flag_drop_one_at_a_time` — six records, one per flag bit
  plus a clean control; iterates one drop-bool active at a time and
  with all defaults; verifies `FilterCounts` increments and which
  records survive.
- `a10_all_flag_drops_disabled_passes_everything`.
- `a11_min_read_length_drops_short_and_empty` + asserts
  `DEFAULT_MIN_READ_LENGTH == 30`; `None` disables the filter.
- `a12_filter_precedence_is_hit_rate_ordered` — three combo records
  pin the cascade order at duplicate → MAPQ → unmapped without
  enumerating every link.
- `a13_empty_stream`.
- `a14_mixed_empty_and_non_empty`.

### Group B — via `new` (7)

- `b1_header_parsing_on_known_good_cram` — `sample_name() == "s1"`,
  `contigs()` matches.
- `b2_cram_4x_rejected_at_header_time` — manually-written 26-byte
  file definition with major=4 → `UnsupportedCramVersion{major:4}`.
- `b3_non_coordinate_sort_rejected` — `@HD SO:queryname` →
  `NotCoordinateSorted`.
- `b4_sample_tag_handling` — single SM:foo, missing SM, two CRAMs
  with disagreeing SMs.
- `b5_contig_list_mismatch_across_crams` — name/order, length, md5
  disagreements.
- `b6_fasta_agreement` — happy path, missing `.fai`, length mismatch.
- `b7_end_to_end_smoke` — three records via `new` produce three
  `MappedRead`s with the expected fields.

## Validation results

Run inside the project's container (`./scripts/dev.sh`) and on the host
where appropriate:

| Command | Outcome |
|---|---|
| `cargo build` (container) | clean |
| `cargo test --lib per_sample_caller` (container) | 22 / 22 passing |
| `cargo test --lib` (container, full lib) | 83 / 83 passing (no regressions in pre-existing tests) |
| `cargo fmt --check` (host) | clean after running `cargo fmt` |
| `cargo clippy --all-targets --all-features -- -D warnings` (host) | clean for `per_sample_caller`; pre-existing warnings/errors in `decompression_pool.rs` and `genotype_merging.rs` remain — out of scope for this slice |

`cargo fmt` and `cargo clippy` are not installed inside the project's
container image and were therefore run on the host. The container has
`cargo` + `rustc` only, which is enough to build and test.

The pre-existing benchmark `cargo bench gvcf_perf` requires a real
VCF on the host (`/home/jose/analyses/g2psol/source_data/TS.vcf.gz`)
and panics if absent — unchanged by this slice.

## Tradeoffs / follow-ups (deferred per plan)

- **Single-threaded decoders.** Plan calls for one decoder thread
  per CRAM eventually. The `Iterator` API is stable for that swap.
- **Per-record allocations.** Each `MappedRead` clones `qname`, `seq`,
  `qual`, and the CIGAR vector. Acceptable as the cost of decoupling
  from noodles; profile before optimising.
- **MD5 trust.** `.fai` length is trusted; CRAM `@SQ M5` is
  trusted. We do not recompute MD5 from the FASTA bytes.
- **`--region` / `.crai`-driven seeking.** Out of scope; the iterator
  API stays stable when added.
- **Stage 1 orchestrator.** `per_sample_caller/mod.rs` is currently
  a pub-mod-only file. The CLI entry, anyhow-context-wrapped error
  flow, and downstream slice wiring (`baq.rs`, `pileup_walker.rs`,
  `phase_chain.rs`, `psf_writer.rs`) are all separate slices.

## Files touched

- New: `src/per_sample_caller/mod.rs`, `cram_input.rs`, `errors.rs`,
  `record_specs.rs`, `cram_files.rs`.
- Modified: `src/lib.rs` (+1 line), `Cargo.toml` (deps), `.gitignore`
  (gitignore `CLAUDE.md` and `.claude/`).
- Operational scaffolding (collateral, not functional): `CLAUDE.md`,
  `.claude/settings.json` — gitignored. These let the assistant work
  autonomously in this project (host read-only allowlist, container
  for build/test, Explore-agent recommendation for dependency API
  lookups).
