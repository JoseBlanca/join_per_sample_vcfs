# ng reference_info — Milestone A (types + the cheap `.fai` reader)

*Implementation report, 2026-07-19. Plan-driven (Milestone A, steps A1–A4), one commit
per step. Design authority: spec [`../../ng/spec/reference_info.md`](../../ng/spec/reference_info.md)
and arch [`../../ng/arch/reference_info.md`](../../ng/arch/reference_info.md); build order
[`../../ng/impl_plan/reference_info.md`](../../ng/impl_plan/reference_info.md). All code in
[`src/ng/reference_info.rs`](../../../../src/ng/reference_info.rs) (+ one `pub mod` line in
[`src/ng/mod.rs`](../../../../src/ng/mod.rs)).*

## Plan

Stand up the reference-info reader's foundations: the data types, the sibling-path helper,
the committed `samtools` oracle, and the cheap `.fai` reader — the `Fasta` streaming pass
(the heart) is Milestone B. Types first within the milestone, one implement → review →
apply → commit loop per step.

## Assumptions / recorded deviations

- **`ReferenceInfoError` variant fields pinned when coding** (arch §6 licenses this). The
  *set* of eight variants is from arch §1.3; the fields are chosen to serve B/C/E
  (`MalformedFasta { path, contig, byte_offset, detail }`, `FastaFaiMismatch { fasta, fai,
  contig, field, detail }`, etc.). Field-guard failures in the `Fai` arm reuse `FaiRead`
  with a synthesised `InvalidData` error — the shape `RawChromReader` already wraps
  `ContigFai::validate` in — rather than adding a variant.
- **A4 reads the `.fai` twice** (a raw pre-scan for the FASTQ six-column shape, then
  `noodles fai::fs::read`). Inherent: noodles' `splitn(5, '\t')` folds a sixth column into
  `line_width` and rejects it as a generic parse error, so the `FastqIndex` diagnosis
  (spec §3.8) cannot be made after the parse. Negligible cost on a <1 KB file; kept over a
  more complex single-read variant.
- **The `Fasta` arm is a `todo!()`** until Milestone B — the intended incremental state
  (no consumer calls `read_reference_info` yet).

## Changes made (by step)

- **A1** (`29e1399`) — the data types: `ContigInfo` (name/length/geometry/md5 in one
  place), `ReferenceInfo` (whole-reference md5 + `Vec<ContigInfo>` in file order) with the
  transitional `contig_list()` projection, `ReferenceSource` (`Fai` | `Fasta { fasta, fai }`),
  and the `#[non_exhaustive]` `ReferenceInfoError` (eight variants, per-variant docs). Leaf
  imports only (`crate::fasta::{ContigEntry, ContigList}`).
- **A2** (`3557e63`) — `sibling_fai_path`: `<fasta>` + `.fai`, no I/O (five lines copied from
  `with_fai_extension`). Sibling *discovery* stays the caller's policy (spec §3.6).
- **A3** (`88d2095`) — test infra + committed `samtools` oracle: `write_wrapped_fasta` (a
  width-configurable FASTA writer, since `ref_seq.rs::build_fasta` is one-line-per-contig);
  the golden reference's `LN`/`M5`/`.fai` from `samtools 1.16.1` (dev container) as
  constants + the committed `synthetic_ref.fa.fai`; and the `faidx.5` man-page worked
  example (LF + CR-LF), verbatim from vendored `htslib/faidx.5`.
- **A4** (`8c8d087`) — `read_reference_info`, the `Fai` arm: parse via `noodles fai::fs::read`
  into `ContigInfo`s (`md5: None`); reject a FASTQ index (`FastqIndex`), a duplicate name
  (`DuplicateContigName`, T2), and a `.fai` failing the field guards (`FaiRead`, T3); a
  missing `.fai` errors (T4).

## Tests added

13 tests in the module's `#[cfg(test)]`:

- **A1**: `contig_list` projects name/length/md5 in order; drops geometry only;
  `MalformedFasta` Display names the contig (T3).
- **A2**: `sibling_fai_path` appends `.fai`, replaces no extension, touches nothing.
- **A3**: the committed golden `.fai` equals the constants; `write_wrapped_fasta`'s geometry
  locates the first and last base of a wrapped contig; the `faidx.5` worked-example tuples
  agree with the constructed bytes under LF and CR-LF (guards B2's vector), incl. the T7
  first-word-is-the-name check.
- **A4**: the golden `.fai` parses to the expected `ContigInfo`s with `md5: None`; the three
  rejections (dup name, FASTQ index, both field guards) each error with the named variant;
  a missing `.fai` errors.

## Validation

Every step validated in the dev container (`./scripts/dev.sh`): `cargo fmt --check` clean;
`cargo clippy --lib --tests --all-features -D warnings` clean; `cargo test` green
(full lib suite **1904 passed** after A4). Cross-checked on the **host** (cargo 1.95.0):
the module and its 13 tests build and pass natively.

**Pre-existing, unrelated:** `cargo clippy --all-targets --all-features -D warnings` has one
failure in `examples/ssr_psp_seqdump.rs:41` (`unnecessary_sort_by`) — it fires on the clean
tree (a toolchain 1.95 patch bump; a warning normally, an error only under `-D warnings`),
is outside `src/ng`, and was left untouched per the production freeze. The library and its
tests are clean under `-D warnings`.

## Tradeoffs / follow-ups (Milestone B onward)

- The `Fasta` streaming pass (geometry + MD5, the heart) — Milestone B, where the committed
  M5/`.cat`/worked-example oracles get consumed.
- `write_fai` (C), the single-flight cache (D), the two entry points (E).
- The MD5 constants (`GoldenContig.md5_hex`) are committed now but only *sanity-checked* in
  A (32 hex chars); B2 is where they become the byte-parity oracle.
