# `pop_var_caller psp-to-pileup` — implementation report

Date: 2026-05-15.
Plan: [psp_to_pileup.md](../../feature_implementation_plans/psp_to_pileup.md).
Slice: small read-side utility that streams a `.psp` artefact as
samtools-mpileup-style text plus a trailing column carrying PSP's
per-allele aggregates.

## Plan

Followed the plan as written, with no scope changes during
implementation. A new subcommand `pop_var_caller psp-to-pileup`
opens an input `.psp` with `PspReader`, maps each record to a
7-column tab-delimited line (CHROM, POS, REF, DEPTH, BASES,
QUALS, ALLELE_DETAILS), and streams to a file or stdout. Optional
`--region chrom[:start-end]` clamps via
`PspReader::region_records`; `--show-chain-ids` appends a chain-id
list per allele in the trailing column.

## Assumptions

Silent choices made where the plan or the spec left details open:

- **`'!'` placeholder for column 6.** Plan said "`'!'`" or
  "single sentinel". Used `'!'` (ASCII 33, Phred 0) because every
  mpileup parser treats it as "missing quality" and it cannot be
  confused with `'*'`, which already has a meaning on column 5
  (deleted base from upstream).
- **Empty `--show-chain-ids` field renders as a trailing
  colon.** When a `--show-chain-ids` run hits an allele whose
  `chain_slots` list is empty, the chain field is emitted as an
  empty token after the colon (`A:5:3:0:0:-1.234:`). Documented
  in `emit_line_chain_ids_empty_list_when_no_chains`. The
  alternative of omitting the colon would break the rule "column
  shape is constant when `--show-chain-ids` is on", which would
  bite downstream parsers.
- **Multi-event deletion alleles are approximated.** When
  `allele.seq.len() < ref_span`, the deleted ref bases are taken
  as the suffix of `REF.seq` past the allele's length. This is
  exact for single-event deletions at offset 0 (the common case)
  and an approximation for multi-event alleles. The plan's
  "MNP encoding loses per-position direction" risk note covers
  this; no extra metadata is recorded.
- **`PathBuf::as_os_str() == "-"`** is how the orchestrator
  distinguishes "stdout" from "file path". This matches what other
  bioinformatics CLIs do; any filesystem with a literal `-` file
  must use `./-` to disambiguate.

## Changes made

- New: [src/pop_var_caller/psp_to_pileup.rs](../../../src/pop_var_caller/psp_to_pileup.rs)
  — clap `PspToPileupArgs`, `PspToPileupError`,
  `run_psp_to_pileup` orchestrator, `parse_region`, `emit_line`
  encoder, indel-marker formatters, and 18 unit tests
  (~430 lines).
- New: [tests/psp_to_pileup_integration.rs](../../../tests/psp_to_pileup_integration.rs)
  — 4 integration tests building a `.psp` via `PspWriter`
  directly (no CRAM round-trip), running the utility against a
  tempfile sink, and asserting the emitted text columns
  (~230 lines).
- Modified: [src/pop_var_caller/mod.rs](../../../src/pop_var_caller/mod.rs)
  — added `pub mod psp_to_pileup;` and re-exports of `PspToPileupArgs`,
  `PspToPileupError`, `run_psp_to_pileup`.
- Modified: [src/pop_var_caller/cli.rs](../../../src/pop_var_caller/cli.rs)
  — added the `PspToPileup(PspToPileupArgs)` variant to
  `PopVarCallerCommand`.
- Modified: [src/main.rs](../../../src/main.rs) — imported
  `run_psp_to_pileup`, added the dispatch arm for the new
  subcommand.

Unmodified (no changes needed to `per_sample_caller/*`,
`Cargo.toml`, or other binary infrastructure — every needed
piece was already public on the PSP reader API).

## Tests added/updated

**Unit (in [psp_to_pileup.rs](../../../src/pop_var_caller/psp_to_pileup.rs)):**

- `parse_region_chrom_only` / `parse_region_chrom_start` /
  `parse_region_chrom_start_end` — region parser happy paths.
- `parse_region_rejects_empty` /
  `parse_region_rejects_missing_chrom` /
  `parse_region_rejects_nonnumeric_start` /
  `parse_region_rejects_nonnumeric_end` /
  `parse_region_rejects_end_before_start` /
  `parse_region_rejects_zero_start` — region parser errors.
- `emit_line_ref_only_all_forward` /
  `emit_line_ref_only_mixed_strand` — REF-only column 5 encoding
  (`.` × fwd + `,` × rev).
- `emit_line_snp_alt_mixed_strand` — SNP ALT emits
  uppercase/lowercase ALT base per strand, REF block precedes
  ALT block in PSP allele order.
- `emit_line_insertion_emits_plus_marker_per_read` — insertion
  emits `<anchor>+N<inserted>` per supporting read, lowercase
  for reverse strand.
- `emit_line_deletion_emits_minus_marker_per_read` — deletion
  emits `<anchor>-N<deleted-ref-bases>` per supporting read,
  lowercase for reverse strand.
- `emit_line_chain_ids_field_off_by_default` /
  `emit_line_chain_ids_field_on_when_flag_set` /
  `emit_line_chain_ids_empty_list_when_no_chains` — `--show-chain-ids`
  toggles the trailing chain field; empty list → trailing colon.
- `allele_detail_includes_placed_counters` — `placed_left` /
  `placed_start` surface through to column 7 in the expected
  position.

**Integration (in [tests/psp_to_pileup_integration.rs](../../../tests/psp_to_pileup_integration.rs)):**

- `streams_every_record_with_seven_columns_per_line` — writes a
  3-record `.psp` (REF-only SNP, insertion, deletion), runs the
  full subcommand against a tempfile sink, asserts exact column
  content for every line.
- `region_flag_clamps_to_one_record` — `--region chr1:20-20`
  filters out records at positions 10 and 30, leaves position 20.
- `unknown_chromosome_in_region_is_a_hard_error` — passing a
  chromosome name not in the input's `[[chromosome]]` array
  returns `PspToPileupError::UnknownChromosome`.
- `show_chain_ids_appends_chain_field` — chain ids appear in the
  trailing position with the documented `;` separator.

## Validation results

All commands run inside the project container via `./scripts/dev.sh`:

- `cargo build --all-targets` — clean.
- `cargo test --lib psp_to_pileup` — **18 passed; 0 failed**.
- `cargo test --test psp_to_pileup_integration` — **4 passed; 0
  failed**.
- `cargo test --lib --tests` — **592 passed; 0 failed**
  across all test targets in the workspace (498 lib + 94 across
  the seven integration suites).
- `cargo clippy --all-targets --all-features -- -D warnings` —
  clean.
- `cargo fmt` — applied once after first writing; rustfmt
  reformatted two long assert macros and one long
  `AlleleObservation::new` literal; no semantic changes.

`cargo test --all-targets --all-features` also runs the benches,
and `gvcf_perf` panics during setup because its fixture
(`/home/jose/analyses/g2psol/source_data/TS.vcf.gz`) does not
exist in this environment. Unrelated to this slice; pre-existing.

Manual smoke test of the CLI surface:

```
$ ./target-container/debug/pop_var_caller psp-to-pileup --help
Stream a .psp as samtools-mpileup-style text plus a trailing column with per-allele aggregates PSP carries

Usage: pop_var_caller psp-to-pileup [OPTIONS] --input <INPUT>

Options:
      --input <INPUT>    Input .psp file
      --output <OUTPUT>  Output path. Use `-` (the default) to stream to stdout [default: -]
      --region <REGION>  Restrict output to a region: …
      --show-chain-ids   Append the per-allele `chain_slots` list (semicolon-separated) …
  -h, --help             Print help
```

## Tradeoffs and follow-ups

Deferred (no work needed now, but on the radar):

- `--bed FILE` for multi-region filtering. Mechanical loop over
  the existing region path; raise when the one-region path is in
  use.
- `--per-allele-row` long-format mode (one line per
  AlleleObservation rather than one line per position).
  Mechanically trivial on top of `emit_line`.
- Compressed output (`.tsv.gz`). User can pipe through `gzip` /
  `bgzip`; not worth a flag.

Non-goals reaffirmed during implementation:

- No FASTA input. PSP's `alleles[0] = REF` invariant gives us
  the REF base for free, including the deleted-ref bases needed
  for the `-N<refbases>` encoding.
- No per-base BQ reconstruction. PSP aggregates BQ into
  `q_sum`; column 6 is filled with `'!'` so consumers cannot
  mistake it for a per-read quality string.
- No atomic tmp+rename for the output. This is a text-dump
  utility; partial output on failure matches `samtools mpileup` /
  `bcftools view` / `awk` convention. Rerun overwrites.

Two implementation-time micro-decisions worth flagging for
future readers:

1. **Region tuple shape.** The orchestrator stores the parsed
   region as `(chrom_id, start, end)` rather than the structured
   `RegionSpec`, dropping the chromosome name once we have the id.
   The id is the only thing `PspReader::region_records` needs, and
   resolving the name up-front means the lookup happens at the CLI
   boundary, where the error message can name the missing
   chromosome.
2. **Drain closure.** The whole-file and region paths share a
   `drain` closure that runs the `RecordsIter` to exhaustion and
   writes lines to the sink. The duplication-avoidance is small
   but the alternative — two near-identical `for` loops — would
   drift on future edits.
