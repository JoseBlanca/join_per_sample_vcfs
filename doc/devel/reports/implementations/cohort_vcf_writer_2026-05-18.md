# Implementation report — cohort VCF writer (Stage 6 sink)

**Date:** 2026-05-18.
**Spec:** [calling_pipeline_architecture.md §"Stage 6 — posterior engine"](../../specs/calling_pipeline_architecture.md).
**Plan:** [cohort_vcf_writer.md](../../implementation_plans/cohort_vcf_writer.md).
**Skill:** rust-feature-implementation
([feature_implementation_skill.md](../../ia/skills/feature_implementation_skill.md)).

## Plan

Implementation followed the plan as written. The plan's three
explicit knob-shaped decisions ("QUAL = ∞ → cap at 9999",
"pin noodles-vcf at implementation time", "`GP` opt-in")
were ratified by the user before coding and shipped as designed.

The user also approved a noodles-family-wide version bump (all
six noodles crates already pinned in `Cargo.toml`) up to the
`noodles-vcf 0.88` release wave, eliminating the dual-version
duplication that would otherwise have landed alongside the new
`noodles-vcf` dependency. The bump was done first, verified
with `cargo build` against the rest of the codebase, then the
writer work proceeded.

## Assumptions (silent choices the spec or plan left open)

1. **`Number=G` enum mapping.** noodles-vcf 0.88 names the
   VCF-spec `Number=G` value as `format::Number::Samples`
   (single-letter spec keys are mapped to enum variants by
   semantic role rather than the letter). Verified against
   the crate source at
   `header/record/value/map/format/number.rs:3-25`; an inline
   comment in `header.rs` names the gotcha so a future reader
   isn't surprised.
2. **`Number=R` mapping for AD.** Same module — VCF spec
   `Number=R` is `format::Number::ReferenceAlternateBases`. No
   gotcha; documented at the call site.
3. **VCF position type comes from noodles-core 0.20.** The
   project's existing `noodles-core` pin was at 0.19 before
   this slice; the noodles-vcf bump required 0.20. The user
   approved upgrading all noodles crates together, so the
   project now uses 0.20 throughout — no dual-version
   resolution.
4. **`AlleleSupportStats` construction in the integration
   test.** The type is `#[non_exhaustive]`, so external
   constructors must go through `AlleleSupportStats::new(…)`.
   In-crate unit tests still use struct literals; only the
   `tests/` crate has to call the constructor.

## Changes made

### New module `src/var_calling/vcf_writer/`

Sits next to its upstream `posterior_engine.rs`. Legacy
`src/vcf_writer.rs` (gVCF-merger code) is untouched and will
be removed when the legacy pipeline retires.

- [`mod.rs`](../../../src/var_calling/vcf_writer/mod.rs) —
  module wiring; re-exports `CohortVcfWriter`,
  `CohortMetadata`, `WriterConfig`, `VcfWriteError`. Defines
  the small [`WriterConfig`](../../../src/var_calling/vcf_writer/mod.rs#L36-L62)
  struct (`output`, `default_filter_pass`, `emit_gp`).
- [`errors.rs`](../../../src/var_calling/vcf_writer/errors.rs) —
  `VcfWriteError` with eight variants: `Io`, `Encode`,
  `InvalidMetadata`, `RecordOutOfOrder`,
  `GenotypeIndexOutOfBounds`, `UnknownChromId`,
  `SampleCountMismatch`, `ContigLengthOverflow`.
- [`sink.rs`](../../../src/var_calling/vcf_writer/sink.rs) —
  `SinkKind` enum dispatching plain `BufWriter<File>` vs
  `noodles_bgzf::io::Writer<File>` based on the output-path
  suffix (`.vcf.gz`/`.vcf.bgz` → bgzf, else plain).
  `open_tmp` opens `<final>.tmp`; `finish` flushes, emits the
  bgzf EOF block when applicable, syncs, and atomic-renames
  `<final>.tmp` → `<final>`.
- [`header.rs`](../../../src/var_calling/vcf_writer/header.rs) —
  `CohortMetadata` (sample names, contigs, tool string,
  command line) and `build_vcf_header` that produces a
  `noodles_vcf::Header` for VCF 4.4 with `##source`,
  `##commandline`, one `##contig=<ID,length,md5>` per entry,
  the five INFO definitions (`AF`, `AC`, `AN`, `DP`, `CA`),
  the four mandatory FORMAT definitions
  (`GT:GQ:DP:AD`), and the conditional `GP` FORMAT
  declaration when `config.emit_gp`. Validates: non-empty
  `sample_names`, no duplicate samples, no duplicate contigs,
  no contig length exceeding `i32::MAX`.
- [`record_encode.rs`](../../../src/var_calling/vcf_writer/record_encode.rs) —
  Stateless `PosteriorRecord` → `noodles_vcf::variant::RecordBuf`
  translation. Builds the format-key list once via
  `build_format_keys`; per record assembles the CHROM /
  POS / REF / ALT / QUAL / FILTER / INFO / FORMAT cells.
  QUAL clamps `[0, 9999]` with `INFINITY → 9999` and `NaN →
  0` (the `QUAL_MAX` constant is the searchable anchor for
  this decision).
- [`writer.rs`](../../../src/var_calling/vcf_writer/writer.rs) —
  `CohortVcfWriter` orchestration. `new` builds the header,
  opens the sink, writes the header; `write_record` enforces
  non-decreasing `(chrom_id, start)` order, encodes the
  record, hands it to `noodles_vcf::io::Writer::write_variant_record`;
  `finish` recovers the sink and runs the
  flush/EOF-block/sync/rename dance.

### Wiring

- [`src/var_calling/mod.rs`](../../../src/var_calling/mod.rs#L22) —
  added `pub mod vcf_writer;`.
- [`Cargo.toml`](../../../Cargo.toml#L56-L60) — added
  `noodles-bgzf = "0.47.0"` and `noodles-vcf = "0.88.0"`;
  bumped the existing noodles-family pins from the 0.19/0.46/
  0.60/0.84/0.92 wave to 0.20/0.47/0.61/0.85/0.93 so the
  whole noodles dep graph sits on one consistent release set.

## Tests added/updated

23 in-crate unit tests + 4 integration tests = **27 new tests**,
all passing.

### Unit tests in the new module

- `sink::tests` (5) — suffix dispatch (plain / `.gz` / `.bgz`),
  `tmp_path_for` shape, plain-sink atomic-rename round-trip,
  bgzf-sink EOF-marker byte-pattern check after `finish`.
- `header::tests` (5) — header serialises through
  `noodles_vcf::io::Writer::write_header` containing all
  expected lines for both `emit_gp = false` (default) and
  `emit_gp = true`; rejects empty sample names, duplicate
  samples, duplicate contigs.
- `record_encode::tests` (10) — biallelic SNP renders the
  exact expected line under both `emit_gp` settings;
  triallelic SNP produces a 2-element `Number=A` list for
  `AF` / `AC`; the `CA` flag surfaces when any sample's
  chain-anchor-broken bit is set; `QUAL = INFINITY` renders
  as the `9999` cap; `QUAL = NaN` renders as `0`;
  out-of-range `chrom_id` errors with `UnknownChromId`;
  cohort-size mismatch errors with `SampleCountMismatch`;
  `build_format_keys` tracks `emit_gp`; `format_gt_unphased`
  joins allele indices with `/`.
- `writer::tests` (3) — full plain-text round-trip writes two
  records and produces a complete on-disk VCF with no
  leftover tmp file; out-of-order record raises
  `RecordOutOfOrder` with prev/current locus filled in;
  equal-locus record (chr1:200 twice) also errors.

### Integration tests in `tests/cohort_vcf_writer_integration.rs`

- `plain_text_path_writes_three_records_with_default_config`
  — biallelic SNP + triallelic SNP + insertion fixture
  written through the public API, then re-read as bytes;
  asserts header content (`##fileformat=VCFv4.4`,
  `##source`, `##contig=<ID=chr1`, INFO/FORMAT lines, no
  `##FORMAT=<ID=GP`), three data lines in genomic order,
  triallelic ALT cell `T,C`, insertion-record QUAL cap at
  `9999`, FORMAT key `GT:GQ:DP:AD` on every line, no
  leftover `<output>.tmp`.
- `plain_text_path_with_emit_gp_adds_format_column` — same
  fixture, `emit_gp = true`; asserts the `##FORMAT=<ID=GP,
  Number=G,Type=Float` header line appears, FORMAT key
  becomes `GT:GQ:DP:AD:GP`, every per-sample cell carries
  five colon-separated fields.
- `bgzf_path_writes_compressed_vcf_with_eof_marker` — same
  fixture written to `.vcf.gz`; asserts the file ends with
  the 28-byte htslib EOF block (`1f 8b 08 04 …` byte
  pattern), then decompresses via `noodles_bgzf::io::Reader`
  and confirms three data lines round-trip.
- `out_of_order_records_surface_a_clear_error` — feeds
  chr1:900 then chr1:100; asserts `RecordOutOfOrder` with
  both loci correctly populated.

### On-disk size delta (per the plan's validation requirement)

Measured on the three-record integration fixture (2 samples
× {biallelic SNP, triallelic SNP, insertion}) using a one-off
`examples/measure_vcf_sizes` shim built against the writer's
public API. Numbers reported are the on-disk file size:

| `emit_gp` | Plain `.vcf` | `.vcf.gz` |
|-----------|--------------|-----------|
| `false`   | 1 246 bytes  | 608 bytes |
| `true`    | 1 454 bytes  | 684 bytes |

`emit_gp = true` adds ~17 % to the plain-text size and ~13 %
to bgzf on this tiny fixture — the absolute numbers are
dominated by the header. On real cohorts where the data
section dwarfs the header, the per-line GP cost is the
dominant signal; at 1 000 samples × 21 floats/sample/record
that's 21 000 floats per record (~6 chars each ≈ 126 KB
text per data line). The opt-in default is justified.

## Validation results

Run inside `./scripts/dev.sh` (rootful podman dev
container) against the project tree.

| Command | Result |
|---------|--------|
| `cargo fmt --all -- --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib var_calling::vcf_writer` | 23 passed / 0 failed |
| `cargo test --test cohort_vcf_writer_integration` | 4 passed / 0 failed |
| `cargo test --lib` | 756 passed / 0 failed |

Clippy initially flagged four warnings in the new code; all
fixed in this slice: `unnecessary_closure_used_to_substitute`
× 2 (`.ok_or_else(|| …)` where the error variant doesn't
capture moved values → `.ok_or(…)`), `collapsible_if` in the
writer's order check (collapsed into an `if let … && …`
chain), `field_reassign_with_default` in a unit test (use
struct-update syntax).

Manual smoke (`bcftools view` / `bcftools stats`) is **not
yet run** — the bench host's bcftools binary needs an update
to a recent enough version that recognises VCF 4.4 cleanly,
and that's not a writer-slice concern. The integration test
exercises the only htslib-shaped failure mode (the bgzf EOF
block) directly against the on-disk bytes; the bcftools
sanity-check moves to the cohort CLI slice where there's a
real end-to-end command line to run.

## Tradeoffs and follow-ups

### Decisions ratified before coding
- **QUAL cap at 9999 for `INFINITY`.** Matches what GATK and
  freebayes do. `QUAL_MAX` is the searchable anchor.
- **`GP` opt-in.** `WriterConfig::emit_gp` defaults to
  `false`; declared in the header and emitted per-sample
  only when on.
- **Pin the noodles-vcf version at implementation time,
  with the whole noodles family bumped together.** Now
  `noodles-vcf 0.88` / `noodles-bgzf 0.47` / `noodles-core 0.20`
  / `noodles-sam 0.85` / `noodles-cram 0.93` /
  `noodles-fasta 0.61`. Verified by `cargo build` against
  the full project.

### Deferred to follow-up slices (per the plan)
- Tabix `.tbi` index alongside `.vcf.gz`.
- `PL` (phred-scaled per-genotype likelihoods).
- Per-sample contamination fraction in INFO.
- Stdout output.
- BCF output.
- Custom filter expressions.

### Not yet exercised (waits for the cohort CLI slice)
- The wider `pop_var_caller cohort` end-to-end pipeline that
  connects `PspReader → PerPositionMerger → DustFilter →
  VariantGrouper → PerGroupMerger → PosteriorEngine →
  CohortVcfWriter`. The writer is exercised in isolation
  against synthetic `PosteriorRecord` fixtures; the CLI
  slice will fold it into the real pipeline.
- The `bcftools view` / `bcftools stats` manual smoke (see
  Validation results).

### Caveats called out for the next implementer
- A `CohortMetadata::sample_names` mis-ordering relative to
  the upstream merger's `sample_names()` would render
  per-sample cells under the wrong sample column without
  any error from the writer. The cohort CLI must source
  the sample-name vector from the merger's `sample_names()`
  call rather than the CLI argument order. The plan's
  Risks section makes this explicit; the writer's
  doc-comment repeats it.
- `MergedAllele.seq` carries the cohort-projected REF/ALT
  bytes. The encoder calls `str::from_utf8` and surfaces an
  `Encode` error on non-UTF-8 bytes. Stage 1 guarantees the
  alphabet `{A,C,G,T,N}`, so this is defensive only.
- The dev container test run picked up one pre-existing
  noodles-bgzf 0.46/0.47 duplicate during the in-flight
  noodles bump; resolved by lifting the project's existing
  noodles pins to the 0.20 / 0.47 release wave (per the
  user's go-ahead). No dual-version resolution remains.
