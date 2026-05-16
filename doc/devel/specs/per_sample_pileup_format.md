# `.psp` — Per-Sample Pileup: byte-layout specification

**Status:** Design closed, 2026-05-13. Ready for implementation; v1.0
to be ratified by the first round-trip test. Replaces the 2026-04-21
PSP draft.

This document specifies the **on-disk byte layout** of the `.psp`
(per-sample pileup) artefact — the file Stage 1 of the calling pipeline
writes per sample, and that Stages 3–6 read back without ever
revisiting the CRAM. It is the byte-level companion to the architecture
document's Stage 2 section. The artefact's high-level shape:

- per-position records (no variant / non-variant distinction);
- five per-allele scalars (observation count, `Σ max(ln_BQ, ln_MQ)`,
  forward-strand count, placed-left count, placed-start count);
- per-allele phase chain identifiers — unique-per-file `u64`s, no
  recycling, no lifecycle markers (see Q-FC8);
- columnar storage within zstd-compressed blocks;
- a mandatory tail-mounted **block index**;
- a plain-text **TOML file header** framed by a `PSP\n` magic, an
  8-byte little-endian length prefix, and a `---END-HEADER---\n`
  sentinel — any tool with `head` / `cat` / `grep` or a TOML parser
  can inspect file metadata without our
  software installed;
- a **VCF-style binary schema** carried in the TOML header as a
  `[[column]]` array — the binary body is self-describing;
- `(major, minor)` versioning, with the rule that minor bumps add
  only optional content.

What is fixed by other specs is referenced here, not re-derived. What
this document fixes for the first time: exact byte widths, varint
encoding, the column-tag registry, the per-record / per-allele wire
shape inside a block, and the file trailer. The reasoning behind
each major decision lives in the "Design decisions log" appendix at
the end.

## Cross-references

- [calling_pipeline_architecture.md](calling_pipeline_architecture.md) §"Stage 2 — per-sample file contract" — high-level record shape, rationale for scalars and phase chains, schema-evolution policy.
- [per_sample_pileup.md](per_sample_pileup.md) — Stage 1, the producer. Defines what the writer puts into every column.
- [freebayes_posterior_gt_probs.md](freebayes_posterior_gt_probs.md) — semantics of the five scalars.

## Encoding conventions

The file has two regions with different encoding rules: a plain-text
**header** that can be inspected with `head` / `cat` / `grep` / any
TOML parser, followed by a **binary body** of zstd-compressed
columnar blocks plus index plus trailer.

### Header (plain-text TOML)

| Convention | Choice |
|---|---|
| Format | [TOML v1.0.0](https://toml.io/en/v1.0.0). |
| Character set | 7-bit US-ASCII for every value in v1.0. UTF-8 is reserved for description-type fields added in future minor bumps (mirroring SAM v1.6's "ASCII by default, UTF-8 only in specific named fields"). |
| Line terminator | LF only (`\n`, `0x0A`). No CR, no CRLF. |
| String quoting | Basic strings (`"..."`). Literal strings (`'...'`) are accepted on read but writers emit basic strings for consistency. |
| Numbers | Decimal integers and TOML's offset date-time for timestamps. |

#### Per-field character set

Per-field rules follow SAM v1.6 §1.2.1 wherever the field has a SAM
analogue. SAM's global rule is "7-bit US-ASCII except in specific
named fields"; we adopt the same posture.

| TOML key | Allowed character set | SAM analogue |
|---|---|---|
| `format-version` | `[0-9.]+` matching `[0-9]+\.[0-9]+` | n/a |
| `sample` | ASCII printable (`[ -~]+`) | `@RG SM` |
| `reference` | ASCII printable (`[ -~]+`) | filenames are ASCII in our domain |
| `created` | TOML offset date-time (all ASCII) | n/a |
| `chromosome.name` | SAM `@SQ SN` regex: `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*` | `@SQ SN` |
| `chromosome.length` | positive integer in `[1, 2^31 − 1]` | `@SQ LN` |
| `chromosome.md5` | exactly 32 lowercase hex characters | `@SQ M5` (SAM §1.3.2) |
| `writer.tool`, `writer.version`, `writer.subcommand` | ASCII printable (`[ -~]+`); no `/`, `\` | n/a (`@PG PN`/`VN` analogue) |
| `writer.input-crams[*]`, `writer.input-fasta` | ASCII printable filename (`[ -~]+`), no `/` or `\` — recorded as **basename only** to avoid leaking host paths (see §"File header — provenance and path stripping") | n/a |
| `writer.parameters.<key>` | TOML bare key for `<key>` (`[A-Za-z0-9_-]+`); value type and range defined per parameter in the producer's config struct | n/a |
| `column.name` | lowercase kebab-case ASCII (`[a-z][a-z0-9-]*`) | n/a — matches the column-tag registry |
| `column.cardinality` | enum: `per-record` \| `per-allele` | n/a |
| `column.shape` | enum: `scalar` \| `list` \| `bytes` | n/a |
| `column.element-type` | enum: `u8`, `u16`, `u32`, `u64`, `i32`, `i64`, `f32`, `f64`, `varint`, `svarint`, `bool` | n/a |
| `column.length-column` | the `name` of another column (lowercase kebab-case ASCII) | n/a |
| `column.required` | bool | n/a |
| `column.description` | ASCII printable; canonical text reproduced verbatim from the spec for the matching `tag` | n/a |
| (future) `description`, `notes` | UTF-8 allowed | `@SQ DS`, `@PG DS`, `@CO` |

The producer validates every value against its rule before writing;
the consumer validates on parse. Non-conforming values are a hard
error (`design_principles.md` §"Errors must not pass silently").

### Body (binary)

Everything from the byte after the header sentinel through the last
byte of the trailer is binary, with these conventions:

| Convention | Choice |
|---|---|
| Endianness | Little-endian for every fixed-width integer and float. |
| Integer widths | `u8`, `u16`, `u32`, `u64` denote unsigned 1/2/4/8-byte little-endian integers; same with `i…` for signed. |
| Float widths | `f32` is IEEE-754 binary32 little-endian; `f64` is binary64. |
| Variable-length unsigned ints (`varint`) | Unsigned LEB128: 7 data bits per byte, MSB set on continuation. Max width capped at 10 bytes (covers `u64`). |
| Variable-length signed ints (`svarint`) | Zig-zag encoded into LEB128 (`(n << 1) ^ (n >> 63)` then unsigned LEB128). Only used where signed deltas appear. |
| Allele sequence | ASCII bytes over `{A, C, G, T, N}`, uppercase. zstd absorbs the run-of-`A`/`C`/`G`/`T` redundancy; a packed encoding is not pursued. |
| Boolean | Single `u8`: `0` for false, `1` for true; other values are an error. |

Files use the file extension `.psp` (per-sample pileup). The format is
not gzip-wrapped at the outer level; zstd compression happens only
within block column payloads (§Block layout).

## File-level structure

```
+-----------------+
| file header     |   (uncompressed, prefix-decodable)
+-----------------+
| block 0         |   ┐
+-----------------+   │ zero or more blocks,
| block 1         |   │ written in genomic order
+-----------------+   │ (chrom_id non-decreasing,
| ...             |   │ pos non-decreasing within chrom).
+-----------------+   │ A block never crosses a
| block N-1       |   │ chromosome boundary.
+-----------------+   ┘
| block index     |   one entry per block; coordinate-keyed
+-----------------+
| file trailer    |   fixed 32 bytes; locates the index, detects
+-----------------+   truncation
```

- A `.psp` may legally contain zero blocks (sample with no covered
  positions, e.g. an empty CRAM). In that case the file is
  `header + empty index + trailer`; the index carries zero entries
  and the trailer still locates it.
- Blocks are written sequentially in genomic order. Their byte
  offsets are recorded once, in the tail-mounted **block index**
  (§"Block index"). A reader that needs region-keyed access reads
  the trailer, then the index, then seeks to the relevant blocks.
- The tail-index pattern is what modern columnar binary formats
  use (Parquet, ORC, Arrow IPC, ZIP's end-of-central-directory).
  It keeps the writer single-pass — blocks are written as they
  finalise, the index accumulates in memory and is dumped before
  the trailer — and the index is atomic with the file content
  rather than a sidecar that can drift.

## File header

The file header is plain-text TOML, framed by a small fixed-byte
magic at the start, an 8-byte little-endian `u64` length prefix
giving the exact byte count of the TOML body, and a fixed sentinel
line at the end. The length prefix is authoritative for locating
the body boundary; the sentinel is kept as a structural cross-check
and so that `head file.psp` ends cleanly on a recognisable line.
Goal: every key piece of metadata about the file (sample, reference,
contigs, writer version, creation time) is readable with standard
Unix tools, with no `.psp`-specific software installed:

```
$ head file.psp
PSP
# .psp v1.0 header
format-version = "1.0"
sample         = "NA12878"
reference      = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
...
```

### Layout

```
+-------------------------+
| PSP\n                   | 4 bytes — magic. Printable ASCII; lets
+-------------------------+   `file` / libmagic identify a `.psp`
                              from offset 0.
| toml_body_length        | 8 bytes — u64 little-endian. Authoritative
+-------------------------+   byte count of the TOML body that follows.
                              Decoded directly from the bytes; no TOML
                              parser is invoked to find the body's
                              extent.
| <TOML body>             | exactly `toml_body_length` bytes of valid
|                         |   TOML v1.0.0, restricted per §"Per-field
|                         |   character set". A writer SHOULD make the
|                         |   last byte a `\n` so the sentinel that
|                         |   follows starts on its own line for clean
|                         |   `head` output, but the byte boundary is
|                         |   set by `toml_body_length`, not by
|                         |   scanning for a newline.
+-------------------------+
| ---END-HEADER---\n      | 17 bytes — sentinel line. Redundant with
+-------------------------+   the length prefix for locating the body
                              boundary, but verified by the reader as
                              a structural cross-check (see "Why this
                              framing" below). The next byte is the
                              first byte of block 0.
```

Constraints:

- The magic is exactly the four bytes `'P'`, `'S'`, `'P'`, `'\n'`.
- `toml_body_length` is in the range `[1, 1048547]` — i.e. at least
  one byte of TOML, and at most 1 MiB minus the framing overhead
  (4 magic + 8 length + 17 sentinel = 29). A reader rejects any
  value outside this range before allocating buffers.
- Total header size (magic through sentinel inclusive) is therefore
  bounded by a hard limit of **1 MiB**.
- The 17 bytes immediately after the TOML body must equal the exact
  byte sequence `---END-HEADER---\n`. Mismatch is a hard error
  ("length prefix and sentinel disagree; writer is buggy or file
  is corrupt").
- Per-field rules (§"Per-field character set") forbid `\n` in
  every v1.0 field, so a well-formed TOML body never accidentally
  contains the sentinel pattern. This is a writer-side invariant;
  the reader does not need to enforce it (the length prefix is
  authoritative).

Why this framing:

- **`PSP\n` as fixed magic** lets `file` / libmagic identify a
  `.psp` from offset 0, and lets a binary reader bail out on
  non-`.psp` input before invoking a TOML parser.
- **Length prefix as the authoritative body boundary.** A
  sentinel-only framing (an earlier draft of this spec) finds the
  body end by scanning for `\n---END-HEADER---\n`. That works only
  if the sentinel byte sequence does not appear inside any TOML
  value — which the per-field rules guarantee for well-formed v1.0
  files but cannot guarantee against producer bugs (TOML v1.0.0
  multi-line basic / literal strings `"""..."""` / `'''...'''`
  legally contain `\n`, and a writer that fails to validate would
  silently produce a file where the scan ends in the wrong place).
  The length prefix moves the body boundary from "wherever the
  reader finds the pattern" to "wherever the writer says the body
  ends," which is what the format actually means. A bug in the
  writer surfaces as a sentinel-mismatch hard error, not as a
  successfully-parsed file with a truncated header.
- **Sentinel kept as a structural cross-check.** Two reasons not
  to delete it: (a) `head file.psp` still ends on a clearly
  recognisable line, preserving the "any tool with head / cat /
  grep can inspect" property the format wants; (b) a length
  prefix that gets corrupted to a smaller value would otherwise
  silently truncate the body — keeping the sentinel turns that
  into a hard mismatch.
- **TOML rather than a custom key-value grammar** so we inherit a
  single canonical spec, mature Rust support (`toml` crate, used
  by Cargo itself), native dates, comments, and existing query
  tooling (`taplo`, `dasel`, `tomlq`) — at the cost of one small
  runtime dependency.
- **Writer cost: none beyond what was already paid.** The writer
  must validate every TOML field before emitting the file
  (§"Per-field character set"); that validation runs on an
  in-memory string. Once validated, that string's byte length is
  known and goes into the length prefix. No back-patching, no
  two-pass write.

### Reader protocol

1. `pread` the first 4 bytes; verify they are `PSP\n`. Non-match →
   not a `.psp` file, abort with a clear error.
2. `pread` the next 8 bytes; decode as a little-endian `u64` →
   `toml_body_length`. Verify the value falls in
   `[1, 1048547]`. Out of range → hard error before any allocation
   (`toml_body_length` is attacker-controllable if a malicious file
   is ever fed to a reader, so the bound is checked first to keep
   the next `pread` from allocating an arbitrary buffer).
3. `pread` `toml_body_length` bytes starting at file offset 12.
   Parse them with a TOML v1.0.0 parser. Any parse failure → hard
   error.
4. Validate every TOML field against the per-field rules
   (§"Per-field character set"). Any rule violation is a hard
   error; the file is rejected.
5. `pread` the next 17 bytes (at file offset `12 + toml_body_length`)
   and verify they exactly equal `---END-HEADER---\n`. Mismatch is
   a hard error ("length prefix and sentinel disagree"); never a
   warning, never a recovery.
6. Block 0 starts at file offset `12 + toml_body_length + 17`
   (i.e. the byte after the sentinel's trailing `\n`).

### TOML schema (v1.0)

Required top-level keys:

| Key | TOML type | Notes |
|---|---|---|
| `format-version` | string | `"MAJOR.MINOR"`. v1.0 writers emit `"1.0"`. |
| `sample` | string | Sample name from the input CRAMs' `@RG SM`, validated identical across all input CRAMs (Stage 1 §"Multi-CRAM ingestion"). |
| `reference` | string | Free-form, typically the FASTA's basename. Diagnostic. |
| `created` | offset date-time | When the file was finalised. RFC 3339 / ISO 8601, e.g. `2026-05-13T10:00:00Z`. |

Required array-of-tables `[[chromosome]]`, in `@SQ` order, **with at
least one entry**. The zero-based index into this array is the
`chrom_id` referenced by blocks (Stage 1 validates that this order
is identical across all input CRAMs and identical to the FASTA). A
file whose `[[chromosome]]` array is empty (zero entries) is rejected
at header parse — there is no such thing as a `.psp` over a reference
with no contigs; the merged CRAM input must always carry an `@SQ`
header, and a producer that somehow tries to write a zero-contig
file has a bug. (This is distinct from the zero-*block* case at
§"File-level structure", which is legal and represents a sample
with no covered positions: a non-empty `[[chromosome]]` list plus
zero data blocks.)

| Key | TOML type | Notes |
|---|---|---|
| `name` | string | Contig name, matching SAM `@SQ SN`. |
| `length` | integer | Contig length in base pairs (≥ 1). |
| `md5` | string | 32-character lowercase hex MD5 of the uppercase reference sequence for this contig (i.e. the contig's bases, uppercased, with no whitespace or non-base characters — exactly the SAM `@SQ M5` semantics). Stage 3 compares this against the FASTA it is given **per contig** and refuses to proceed on any mismatch — never silently. Per-contig (rather than whole-file) MD5 means re-wrapping a FASTA, regenerating its `.fai` sidecar, or re-downloading from a different mirror does **not** cause spurious rejections, but any actual change in a contig's sequence does, with a clear per-contig error. |

Required table `[writer]` — provenance for what produced the file.
See §"File header — provenance and path stripping" below for the
full schema and the rationale for the design (no literal `argv`,
basenames only, every effective parameter recorded).

Required array-of-tables `[[column]]` — the **binary schema**: a
declaration of every column the file's binary body carries, with
its tag, cardinality, shape, type, and a canonical one-line
description. See §"File header — binary schema" below.

**Unknown top-level keys.** A reader on a v1.0 file MUST NOT reject
the file because an unrecognised top-level key appears; it skips
the key. This is the mechanism by which future minor bumps can
add optional file-scope content (see §"File-scope extensions
(deferred to v1.x)" and §"Versioning policy").

### File header — provenance and path stripping

The `[writer]` table records exactly what produced the file, in a
form that supports reproducibility on any host without leaking the
producer's directory structure or username.

Required keys directly under `[writer]`:

| Key | TOML type | Notes |
|---|---|---|
| `tool` | string | Executable name, e.g. `"join_per_sample_vcfs"`. No path component. |
| `version` | string | Tool version, e.g. `"0.3.0"`. Format is `tool`-defined; conventionally `MAJOR.MINOR.PATCH`. |
| `subcommand` | string | The CLI subcommand actually invoked, e.g. `"per-sample"`. Single token, no path component. |
| `input-crams` | array of strings | The CRAM files Stage 1 ingested, recorded as **basenames only** (no directory component, no `/` or `\`). One entry per input CRAM, in the order passed on the CLI. |
| `input-fasta` | string | The reference FASTA, recorded as **basename only**. |

Required sub-table `[writer.parameters]` — one key per CLI-exposed
knob, value = the effective value the writer ran with. Defaults and
user-overrides are recorded identically; a re-run that explicitly
sets a default to its default is a no-op, so the distinction adds
no value here. Examples (one row per knob the producer's config
structs expose):

```toml
[writer.parameters]
min-mapq               = 30
max-snp-column-depth   = 8000
max-indel-column-depth = 250
drop-qc-fail           = true
drop-duplicate         = true
min-read-length        = 30
# ... etc.
```

The parameter set the writer emits matches its CLI surface
exactly: every knob exposed on the CLI appears here, no more, no
less. A new knob added in a v1.x writer becomes a new key in the
parameters sub-table; readers that don't know it ignore it (per
the unknown-key rule).

**Why no literal `command-line` / `argv`.** A literal `argv` array
necessarily carries full filesystem paths, which leak host
directory structure (`/home/jblanca/projects/cohort_2026/data/...`)
and the running user's name. The structured form above carries
strictly more *useful* information for re-running and strictly
less *incidental* information about the producer's host. A
reader that wants to re-run on its own host substitutes its own
paths to the recorded basenames:

```
$ join_per_sample_vcfs per-sample \
      --reference <your-path>/GRCh38.fa \
      --min-mapq 30 \
      --max-snp-column-depth 8000 \
      ... \
      <your-path>/sample_lane1.cram
```

**Why basenames, not `~`-anonymised paths.** Replacing `$HOME`
with `~` would still leak project structure
(`~/projects/cohort_2026/...`). Basename-only forecloses both
classes of leak with a single rule and is sufficient for
re-running.

**What's not in `[writer]`.** No `cwd`, no environment variables,
no hostname, no UID/GID. The format records what the writer *did*,
not where it lived.

### File header — binary schema

The `[[column]]` array of tables in the header declares every column
the binary body carries: its tag (which is what the per-block manifest
references), its scope (per-record vs per-allele), its on-the-wire
shape, its element type, whether it is required, and a one-line
description. The intent is **VCF-style self-description**: a tool
that has only this file and the closed type vocabulary below can
decode the binary body without consulting the spec.

The `[[column]]` array doubles as the file-level "columns used"
manifest: it lists every column the file's binary body contains,
in full schema detail rather than as a bare list of names.

#### `[[column]]` entry schema

Each `[[column]]` table carries the following keys:

| Key | TOML type | Notes |
|---|---|---|
| `tag` | integer | The numeric column tag the per-block manifest references. Conventionally written in hex (`0x11`); decimal is also legal TOML. Must be unique within the file and fall in a range valid for the file's `format-version` (see §"Reserved column-tag ranges"). |
| `name` | string | Lowercase kebab-case, e.g. `"allele-q-sum-log"`. Stable across versions. |
| `cardinality` | string | Closed enum: `"per-record"` (one entry per record, count = `n_records`) or `"per-allele"` (one entry per allele, count = `n_total_alleles`). |
| `shape` | string | Closed enum: `"scalar"`, `"list"`, or `"bytes"`. |
| `element-type` | string | Closed enum: `"u8"`, `"u16"`, `"u32"`, `"u64"`, `"i32"`, `"i64"`, `"f32"`, `"f64"`, `"varint"`, `"svarint"`, `"bool"`. Required when `shape != "bytes"`; absent when `shape == "bytes"`. |
| `length-column` | string | Only present (and required) when `shape == "bytes"`. Names another column with the same `cardinality` whose values give per-entry byte lengths. |
| `required` | bool | If true, a reader that does not recognise this `name`/`tag` pair must reject the file. If false, the reader may skip the column. |
| `description` | string | One-line human-readable text. Canonical for the column's `tag` at a given `format-version`; the writer reproduces the spec text verbatim. |

The TOML allows any extra keys per entry; a reader must skip
unknown keys (per the unknown-key rule). This is the hatch for
v1.x columns that need additional metadata (e.g. a per-column
checksum policy, a unit-of-measure annotation).

#### Encoding rules per (`cardinality`, `shape`)

The column's on-the-wire encoding is fully determined by its
`cardinality` and `shape`:

| (`cardinality`, `shape`) | Encoding inside the column's payload |
|---|---|
| `(per-record, scalar)` | `n_records` `element-type` values, packed end to end. |
| `(per-allele, scalar)` | `n_total_alleles` `element-type` values, packed end to end. |
| `(per-record, list)` | `n_records` entries, each: `varint` count `k`, then `k` `element-type` values. |
| `(per-allele, list)` | `n_total_alleles` entries, each: `varint` count `k`, then `k` `element-type` values. |
| `(per-record, bytes)` | flat byte stream; chunk by walking the `length-column` (which must be `(per-record, scalar)` with `element-type = "varint"`). |
| `(per-allele, bytes)` | flat byte stream; chunk by walking the `length-column` (which must be `(per-allele, scalar)` with `element-type = "varint"`). |

A reader that knows the closed type vocabulary and the closed
`(cardinality, shape)` matrix can decode any column the header
declares, even ones it has never seen before.

#### v1.0 binary schema

A v1.0 file carries one `[[column]]` entry per row of §"Required
columns in v1.0", in the same order. The writer emits the
description text verbatim from the spec. The v1.0 schema is
reproduced inline in §"Example v1.0 header" below.

#### Header-binary consistency: required reader checks

> **Note to whoever implements the reader.**
>
> A `.psp` file is **only valid** when the header's binary schema and
> the binary body's actual layout agree on every detail. If a reader
> discovers a disagreement, it must **abort with a clear error**.
> Never silently reconcile. Never fall back to "what the spec says
> for this format-version". Never patch up the data. The whole
> point of the in-file schema is to make the file self-describing;
> a writer that contradicts itself has produced a corrupt or buggy
> artefact, and continuing to read would silently corrupt downstream
> analyses.
>
> This is the project's "errors must not pass silently" principle
> (`design_principles.md` §General principle 3) applied at the
> sharpest edge of the format.

Concrete checks the reader must perform at file-open (and per-block
where noted):

1. **Header well-formedness.** TOML parses cleanly; every required
   key/table is present; every value matches its per-field rule
   (§"Per-field character set").
2. **Schema completeness.** Every `[[column]]` entry has the keys
   its `(cardinality, shape)` requires (`element-type` for non-`bytes`
   shapes; `length-column` for `bytes` shapes; both are mutually
   exclusive).
3. **`length-column` references.** For every `shape = "bytes"`
   column, its `length-column` names another column in the schema,
   that column has the same `cardinality`, `shape = "scalar"`, and
   `element-type = "varint"`. Otherwise reject.
4. **Tag uniqueness and range.** Tags are unique within the file;
   each tag falls in a range valid for the file's `format-version`
   (no minor-bump tags in a v1.0 file, no major-bump tags in a
   v1.x file, no `0x00` sentinel).
5. **Required-column recognition.** Every `[[column]]` with
   `required = true` has a `name`/`tag` the reader recognises.
   Otherwise abort, naming the column.
6. **Schema-vs-registry agreement.** For every recognised column,
   the in-file schema's `cardinality`, `shape`, `element-type`,
   `length-column`, and `required` flags match what the reader's
   own column-tag registry says for that `name`/`tag` at the file's
   `format-version`. Any disagreement → abort, naming the column
   and the disagreeing field.
7. **Per-block manifest agreement** (checked at each block):
    - Every tag in the per-block manifest is declared by a
      `[[column]]` entry in the header.
    - The manifest's column ordering satisfies §"Column ordering
      within a block" (strictly ascending tags, no duplicates).
    - **`uncompressed_len` a-priori prediction, where the schema
      allows it.** The reader cross-checks `uncompressed_len`
      *before* decompression for columns whose byte size is
      determined by the block-level counts and the column's
      schema alone — i.e. independent of the column's actual
      values. Two cases:
        1. **Fixed-width scalars.** `(per-record, scalar)` and
           `(per-allele, scalar)` columns whose `element-type` is
           one of `u8`, `u16`, `u32`, `u64`, `i32`, `i64`, `f32`,
           `f64`, `bool`. Predicted size:
           `n × sizeof(element-type)`, where `n` is `n_records` for
           per-record columns and `n_total_alleles` for per-allele
           columns.
        2. **`bytes` columns, once their `length-column` has been
           decompressed.** `(per-record, bytes)` and `(per-allele,
           bytes)` columns. Predicted size: the sum of the
           `length-column`'s varint entries. (Because this needs
           the `length-column` payload, it is checked *after* the
           length-column has been decompressed; the column ordering
           rule guarantees the length-column appears first since
           every `bytes` column's `length-column` has a numerically
           lower tag in the v1.0 registry, but a reader implementing
           this rule generically should treat it as a deferred
           check rather than relying on that ordering coincidence.)
      A mismatch in any of these is a hard error: the manifest and
      the schema disagree about a column whose size is knowable in
      advance.
    - **`uncompressed_len` is not a-priori predictable for the
      remaining shapes.** Variable-width and list shapes carry data
      whose byte size depends on the values themselves and cannot
      be predicted from the manifest alone. The reader does not
      attempt a prediction check for them. Concretely:
        - `(per-record, scalar)` and `(per-allele, scalar)` with
          `element-type = "varint"` or `"svarint"`;
        - every `(per-record, list)` and `(per-allele, list)`
          column, regardless of `element-type`, because the
          per-entry leading varint count makes the per-entry
          size value-dependent.
      For these columns, the manifest's `uncompressed_len` is what
      the writer wrote; the only pre-decompression check is that
      it is non-zero whenever there is at least one entry of the
      appropriate cardinality. Structural correctness is checked
      during decode (see below).
    - **Decompressed size matches the manifest, always.** After
      zstd decompression of any column payload, the produced byte
      count must equal the manifest's `uncompressed_len`. This
      applies uniformly to every column shape and catches torn or
      truncated zstd frames whose internal checksum happens to
      pass.
    - **Structural correctness during decode** (for shapes where
      it can fire):
        - **List columns.** Walking the decompressed payload, each
          per-entry varint count `k` must be followed by `k`
          element-type values that fit within the remaining
          buffer; after processing `n_records` (per-record list)
          or `n_total_alleles` (per-allele list) entries, the
          cursor must land exactly on the buffer's end.
        - **Varint / svarint scalar columns.** After `n_records`
          (per-record) or `n_total_alleles` (per-allele) varint
          decodings, the cursor must land exactly on the buffer's
          end. No partial varint at the tail.
        - **`bytes` columns.** Already covered above: the sum of
          the length-column entries must equal the decompressed
          size; concatenating the sliced byte ranges must consume
          the buffer fully.
      Any over- or under-run is a hard error, naming the column,
      the block, and the entry index at which the mismatch was
      detected.
8. **No surprises in skipped columns.** A reader skipping a
   `required = false` column it does not understand may not
   silently consume bytes beyond the per-block-manifest's
   `compressed_len` for that column.
9. **No non-finite floats.** Every IEEE 754 float column
   (`element-type = "f32"` or `"f64"`) is checked entry-by-entry
   after decompression: NaN, +∞, and −∞ are forbidden. A
   non-finite value is a producer bug — the per-allele scalars are
   finite sums by construction — and a reader that encountered one
   and proceeded would silently corrupt downstream likelihoods.
   The error names the column, the block, and the offending entry
   index.
10. **Phase-chain identifier ordering** (per record). Every
    `allele_chain_ids` list is strictly ascending with no
    duplicates. The reader rejects any record that violates this
    well-formedness check, naming the block, the record index, and
    the allele index.

(Checks #10–11 in earlier drafts of this document — the
"phase-chain active-set consistency" and "inter-block phase-chain
continuity" rules — are gone along with the recycled-chain-id
design. The unique-per-file `u64` chain ids defined in
`phase_chain.md` §6 make active-set bookkeeping unnecessary on
both the writer and reader sides.)

Failing #10 is a hard error.

### Example v1.0 header

The 8-byte `toml_body_length` prefix that sits between `PSP\n` and
the TOML body is binary and is omitted from this rendering for
readability; on disk, those eight bytes are a little-endian `u64`
giving the exact byte count of the TOML between the prefix and the
sentinel.

```
PSP
# .psp v1.0 header
format-version = "1.0"
sample         = "NA12878"
reference      = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
created        = 2026-05-13T10:00:00Z

[[chromosome]]
name   = "chr1"
length = 248956422
md5    = "6aef897c3d6ff0c78aff06ac189178dd"

[[chromosome]]
name   = "chr2"
length = 242193529
md5    = "f98db672eb0993dcfdabafe2a882905c"

# ... more chromosomes ...

[[chromosome]]
name   = "chrM"
length = 16569
md5    = "c68f52674c9fb33aef52dcf399755519"

[writer]
tool        = "join_per_sample_vcfs"
version     = "0.3.0"
subcommand  = "per-sample"
input-crams = ["sample_lane1.cram", "sample_lane2.cram"]
input-fasta = "GRCh38_full_analysis_set_plus_decoy_hla.fa"

[writer.parameters]
min-mapq               = 30
max-snp-column-depth   = 8000
max-indel-column-depth = 250
drop-qc-fail           = true
drop-duplicate         = true
min-read-length        = 30

# ---- Binary schema: declarations of every column in the body ----

[[column]]
tag          = 0x01
name         = "delta-pos"
cardinality  = "per-record"
shape        = "scalar"
element-type = "varint"
required     = true
description  = "Distance to the previous record's position. First record in a block has value 0 and uses the block header's first-pos."

[[column]]
tag          = 0x02
name         = "n-alleles"
cardinality  = "per-record"
shape        = "scalar"
element-type = "varint"
required     = true
description  = "Number of alleles in this record (>= 1). Sum across records equals n_total_alleles."

[[column]]
tag          = 0x03
name         = "allele-seq-len"
cardinality  = "per-allele"
shape        = "scalar"
element-type = "varint"
required     = true
description  = "Byte length of each allele's sequence. Drives chunking of allele-seq."

[[column]]
tag           = 0x04
name          = "allele-seq"
cardinality   = "per-allele"
shape         = "bytes"
length-column = "allele-seq-len"
required      = true
description   = "Allele sequence bytes (uppercase ASCII over {A,C,G,T,N}). REF is the first allele in each record."

[[column]]
tag          = 0x10
name         = "allele-obs-count"
cardinality  = "per-allele"
shape        = "scalar"
element-type = "u32"
required     = true
description  = "Observation count: reads supporting this allele."

[[column]]
tag          = 0x11
name         = "allele-q-sum-log"
cardinality  = "per-allele"
shape        = "scalar"
element-type = "f64"
required     = true
description  = "Sum over supporting reads of max(ln_BQ_BAQ, ln_MQ). The freebayes per-read combination of base and mapping quality."

[[column]]
tag          = 0x12
name         = "allele-fwd-count"
cardinality  = "per-allele"
shape        = "scalar"
element-type = "u32"
required     = true
description  = "Forward-strand count: reads on the forward strand supporting this allele."

[[column]]
tag          = 0x13
name         = "allele-placed-left-count"
cardinality  = "per-allele"
shape        = "scalar"
element-type = "u32"
required     = true
description  = "Reads whose 5' end is to the left of the record's position (freebayes' placedLeft)."

[[column]]
tag          = 0x14
name         = "allele-placed-start-count"
cardinality  = "per-allele"
shape        = "scalar"
element-type = "u32"
required     = true
description  = "Reads whose 5' end is the record's position (freebayes' placedStart)."

[[column]]
tag          = 0x22
name         = "allele-chain-ids"
cardinality  = "per-allele"
shape        = "list"
element-type = "u64"
required     = true
description  = "Phase-chain identifiers contributing to each allele observation. Ascending fixed-width little-endian u64s. Identifiers are unique within the .psp file and never recycled; two observations sharing an identifier came from the same read or read-pair in this sample. Tags 0x20 and 0x21 (the old new-chain-ids / expired-chain-ids lifecycle markers) are reserved-unused under the unique-id design — see Q-FC8."
---END-HEADER---
```

Comments are TOML-legal and survive the header roundtrip; writers
may emit them or not, readers ignore them.

### Header size in practice

The chromosome list is the dominant contributor. Per-contig TOML
overhead is ~90–100 bytes (name + length + 32-char `md5`):

- ~200 contigs (typical clean reference) → ~20 KB.
- ~500 contigs (with decoys / alts) → ~50 KB.
- ~3000 contigs (some plant references with unplaced scaffolds)
  → ~300 KB — still well under the 1 MiB cap.

Trivial against multi-GB block bodies.

## Block layout

A block holds the records of a single contiguous run of positions on
a single chromosome, encoded columnar and zstd-compressed.

```
block:
  block_header:                       (uncompressed)
    chrom_id:           varint
    first_pos:          varint        // 1-based reference position of the
                                       // first record in the block
    n_records:          varint        // number of per-position records in
                                       // the block
    n_total_alleles:    varint        // sum of n_alleles over all records;
                                       // lets readers size per-allele
                                       // columns up front
    n_columns:          varint        // number of column-manifest entries
                                       // that follow. May be less than the
                                       // header's [[column]] count if the
                                       // block legitimately omits an
                                       // optional column it has no data for.
    per column (n_columns entries):
      column_tag:       varint        // selects from the header's [[column]]
                                       // array. Required-ness, cardinality,
                                       // shape, and element-type are
                                       // resolved against the header
                                       // schema; the per-block manifest
                                       // does not redeclare them.
      compressed_len:   varint        // byte length of the column's
                                       // zstd-compressed payload
      uncompressed_len: varint        // expected size after decompression.
                                       // The reader always checks that the
                                       // actual decompressed payload is
                                       // exactly this many bytes (torn-frame
                                       // catcher). For column shapes whose
                                       // size is determined by the schema
                                       // and the block counts alone
                                       // (fixed-width scalars; bytes columns
                                       // once their length-column is known)
                                       // the reader additionally cross-checks
                                       // this value against the schema-
                                       // predicted size, before decompression.
                                       // For variable-width / list shapes
                                       // no a-priori prediction is possible
                                       // and the structural correctness of
                                       // the payload is checked during
                                       // decode instead. Full rules in
                                       // §"Header-binary consistency".
  column_payloads:                    (concatenated, each independently
                                       zstd-compressed)
    payload_0:          compressed_len[0] bytes
    payload_1:          compressed_len[1] bytes
    ...
    payload_{n_columns-1}: compressed_len[n_columns-1] bytes
```

The per-block manifest carries no redundant schema fields. The
header's `[[column]]` array is the single source of truth for what a
tag means; the per-block manifest only adds what genuinely varies
block-to-block: which columns this block actually wrote, and how
big their payloads are.

### Phase-chain state across blocks

Phase chains routinely straddle block boundaries: a chain that
begins in block N and is still active when the writer cuts the
block continues into block N+1. Under the unique-`u64`-id design
no special handling is required at the block boundary — chain ids
are unique within the entire `.psp` file, so the same id can
appear in any block and refer to the same molecule it referred to
in any earlier block. The block header therefore carries no
active-chain snapshot; a reader landing on block N+1 via the tail
index needs no carried-in state to materialise records.

(Earlier drafts of this document specified an
`active_chain_ids_at_block_start` snapshot field plus inter-block
continuity checks. Both have been retired along with the recycled-
chain-id design; see Q-FC8.)

### Why per-column zstd, not whole-block zstd

This is the form the architecture's §"Per-position record layout
and compression" already commits to. Recording the consequence for
the byte layout here: a reader that does not recognise a non-required
column tag can skip its `compressed_len` bytes and continue, without
having to decompress data it would only discard anyway. Per-column
framing also lets readers parallelise decompression by column when
the work-stealing tradeoff is favourable (rayon over columns,
trivially independent).

### Block sizing

**Target uncompressed block size:** ~16 MiB (writer hardcoded; not
exposed on the CLI in v1). The writer accumulates records until
either:

- the projected uncompressed column-payload total reaches the
  target, **or**
- a chromosome boundary is reached (blocks never cross chromosomes).

The target is a writer-side knob, not part of the on-disk contract.
A reader handles any block size from 1 record up to the varint
maximum (`2^64 − 1` records, in practice limited by memory).

**Block invariants (must hold for every block).** A reader rejects
the file if any block in the body violates these:

- `n_records ≥ 1` — empty blocks are forbidden. A writer that has
  no records to emit (e.g. a chromosome with zero coverage) emits
  no block for that chromosome, not a zero-record one.
- `n_total_alleles ≥ n_records` — every record has at least one
  allele (the REF allele, always present per architecture doc
  §"Allele and record conventions").
- All records in a block share the same `chrom_id` (the block
  header's `chrom_id`), and their positions are strictly increasing
  in genomic order.
- **Per-record `delta_pos` invariant.** The first record in a block
  has `delta_pos = 0` (and its position is the block header's
  `first_pos`); every subsequent record has `delta_pos ≥ 1`. A
  `delta_pos = 0` on any record other than the first, or a
  `delta_pos` value that would push the position past the contig's
  declared `length`, is a hard error. This makes the "strictly
  increasing positions" rule above checkable column-locally,
  without having to materialise positions first.
- Every per-allele `allele_chain_ids` list is strictly ascending
  and pairwise distinct.

Rationale for 16 MiB rather than the smaller (~4 MiB) figure earlier
drafts considered: the artefact is expected to be large — WGS at
the project's target coverage already produces multi-GB files —
so any single sample run cares far more about compression ratio
than about per-block decompression latency. Larger blocks give
zstd more context (better ratio) and have negligible decompression
overhead on modern CPUs. A future benchmark on real data may
change the default; that's a writer change, not a format change.

### Compression

zstd is mandatory for every column payload in v1. The compression
level is a writer-internal choice; the on-disk format does not
record it and a reader decompresses any level transparently.

The v1 writer uses **zstd level 9**, hardcoded. It is not exposed
on the CLI in v1: the artefact is written once per sample and read
many times, the level-vs-ratio trade-off is well within zstd's
useful range at 9, and exposing it as a knob is more configuration
surface than a standard user benefits from. Should benchmarking
on real data later justify a different default, that is a
writer change (and the new tool version is recorded in
`[writer]`), not a format change.

**Dictionary:** not used in v1. Revisit only if observed compression
ratio at small block sizes is poor.

**Per-block checksum:** rely on zstd's built-in per-frame checksum.
Each column payload is its own zstd frame and so carries its own
checksum; corruption inside a payload fails the frame integrity
check at decompression time. No separate block-level CRC is added.

## Logical record content and column shapes

Each block expresses a sequence of N per-position records, each with
some number of alleles. The logical content of one record:

- a position (encoded implicitly via `delta_pos` + the block's
  `first_pos`);
- `n_alleles`: at least 1;
- per allele: a sequence (the literal allele) and five scalars;
- phase chain lifecycle markers (`new_chains`, `expired_chains`);
- per allele: the list of phase chain ids that contributed
  this observation.

Stored columnar: each of these items becomes one column entry per
record (per-record items) or one column entry per allele
(per-allele items). The block's `n_records` sizes the per-record
columns; `n_total_alleles` sizes the per-allele columns; an
`n_alleles` column maps records onto their allele runs.

The column-tag registry below names each column, gives its
per-entry type, and declares whether it's required in v1.

### Required columns in v1.0

Every v1.0 block carries exactly these required columns, in this
order. Optional v1.0 columns are introduced under §"Optional columns
defined in v1.0" — currently none.

This table **is** the v1.0 column-tag registry. It is also what a
v1.0 writer emits verbatim into the file header's `[[column]]` array
(§"File header — binary schema"). The spec and the in-file schema are
two views of the same data; readers verify they agree
(§"Header-binary consistency: required reader checks").

| Tag (hex) | Name | Per-record / per-allele | Type per entry | Notes |
|---|---|---|---|---|
| `0x01` | `delta_pos` | per-record (N entries) | `varint` | Distance to the previous record's position. For the first record in a block, the value is `0` and the position is `first_pos` (from the block header). |
| `0x02` | `n_alleles` | per-record (N entries) | `varint` | At least `1` per record. The sum equals the block's `n_total_alleles`. |
| `0x03` | `allele_seq_len` | per-allele (M entries, M = `n_total_alleles`) | `varint` | Length in bytes of each allele's sequence. **Minimum `1`** — no allele is zero-length (SNP = 1, deletion ALT = 1 anchor base, insertion ALT = 1+ bytes). **Maximum `10_000`** — a hard sanity bound on a single allele's sequence length. Real biological alleles are far shorter (the longest realistic insertion in the project's domain is well under 1 kb); a value over the cap is a producer bug, and a reader that proceeded would either OOM on a 4 GiB varint or silently consume the rest of the block as allele bytes. A reader rejects any record carrying `allele_seq_len < 1` or `allele_seq_len > 10_000`, naming the block, the record, and the allele index. Within a record, the first allele is REF (architecture doc §"Allele and record conventions"). |
| `0x04` | `allele_seq` | per-allele bytes | concatenation | The allele sequences laid end to end, decoded by walking `allele_seq_len`. Bytes are uppercase ASCII over `{A, C, G, T, N}`. |
| `0x10` | `allele_obs_count` | per-allele (M entries) | `u32` | Observation count scalar (reads supporting this allele). `u32` covers the architecture's per-column depth caps with headroom. |
| `0x11` | `allele_q_sum_log` | per-allele (M entries) | `f64` | `Σ max(ln_BQ_BAQ, ln_MQ)` over supporting reads. Stored in ln-units. Matches the implementation's `AlleleSupportStats::q_sum: f64` ([pileup/mod.rs:357](../../src/per_sample_pileup/pileup/mod.rs#L357)); chosen to preserve the precision the accumulator uses, with no truncation at the I/O boundary. |
| `0x12` | `allele_fwd_count` | per-allele (M entries) | `u32` | Forward-strand count. |
| `0x13` | `allele_placed_left_count` | per-allele (M entries) | `u32` | Reads whose 5′ end is to the left of the record's position (freebayes' `placedLeft`). |
| `0x14` | `allele_placed_start_count` | per-allele (M entries) | `u32` | Reads whose 5′ end *is* the record's position (`placedStart`). |
| `0x22` | `allele_chain_ids` | per-allele list | varint count + `u64` LE ids | Per allele: a `varint` count `k`, then `k` little-endian `u64` chain ids in ascending order. Ids are unique within the `.psp` file and never recycled; two observations sharing an id came from the same read or read-pair in this sample. Matches the implementation's `ChainId = u64` ([pileup/chain_id_allocator.rs](../../src/per_sample_pileup/pileup/chain_id_allocator.rs)). Tags `0x20` and `0x21` (the old `new_chain_ids` / `expired_chain_ids` lifecycle markers) are reserved-unused under the unique-id design — see Q-FC8. |

A reader walks the columns in tandem; the per-record columns give
the position and allele-count structure, and the per-allele columns
are sliced by `n_alleles` to attribute entries to records.

### Optional columns defined in v1.0

None. The optional-column hatch exists for additions that v1.x will
introduce as minor bumps without breaking v1.0 readers; the
architecture document anticipates one already (an exact-mixture
quality scalar for contamination handling). When v1.1 ships such a
column, this section lists it with a stable tag in the range
reserved below.

### Reserved column-tag ranges

| Range | Use |
|---|---|
| `0x00` | Reserved sentinel; never a valid tag. |
| `0x01` – `0x0F` | Per-record structural columns (position, allele structure). |
| `0x10` – `0x1F` | Per-allele scalar columns. |
| `0x20` – `0x2F` | Phase-chain columns. |
| `0x30` – `0x7F` | Reserved for future per-position / per-allele columns introduced as minor bumps. |
| `0x80` – `0xFE` | Reserved for future major-version-only columns. |
| `0xFF` and `≥ 0x100` | Available; encoder/decoder must handle varint widths. |

Tags inside reserved-for-future ranges that a v1.0 file carries are
an error: a v1.0 writer never emits them, and a reader encountering
one in a file claiming `version_major = 1` rejects the file.

### Column ordering within a block

**Strict requirement: per-block manifest entries are listed in
ascending `column_tag` order, with no duplicates.** Writers sort
at emit time; readers reject any block whose manifest is
out-of-order or contains duplicates with a clear error.

The writer-side cost is zero (sorting a small list of tags) and
the wins are concrete:

- byte-deterministic output across writers makes `psp dump` and
  hex-level diffs reproducible;
- the reader gets a cheap O(n) check that catches a class of
  writer bugs (forgot to sort, accidentally emitted a column
  twice) without parsing payloads;
- file-format-level consistency stays in the same posture as the
  header/binary consistency checks (§"Header-binary consistency").

## File-scope extensions (deferred to v1.x)

v1.0 defines no file-scope extension mechanism beyond the binary
schema and the writer/provenance tables. The architecture
document anticipates eventually adding file-scope extras (e.g. a
per-sample contamination estimate, a notes field) but none exist
yet, so the v1.0 spec carries no dedicated table for them.

When the first such extension is needed, it is introduced in a
v1.x minor bump as either:

- a new top-level TOML key or table (e.g. `contamination-estimate = 0.012`
  or `[notes]`), under the rule below, or
- a new `[[column]]` if it is per-record/per-allele data that belongs
  in the binary body.

**Minor-bump rule (load-bearing).** Any content a v1.x minor bump
adds to the header MUST be **optional** — a v1.0 reader, which
applies the "unknown top-level keys are silently skipped" rule
(§"TOML schema"), must remain able to consume the file. Required
additions necessitate a major bump. The same rule applies to
`[[column]]` entries: minor bumps add `required = false` columns
only; promoting an optional column to required is a major bump.

This rule is what keeps the cohort re-callable property
(architecture doc constraint 5) without a one-time legacy shim: a
v1.0 writer's `.psp` files are forever consumable by future minor-
bump readers, and minor-bump writers' files remain consumable by
v1.0 readers.

## Block index

A coordinate-keyed index of every block, written immediately after
the last data block and before the trailer. The index is
**mandatory** in v1: every v1 writer emits it and every v1 reader
can rely on its presence.

```
block_index:
  per block (n_blocks entries, in the order the blocks were
             written — which is genomic order):
    chrom_id:         varint
    first_pos:        varint   // 1-based reference position of the
                                // block's first record
    last_pos:         varint   // 1-based reference position of the
                                // block's last record (inclusive).
                                // first_pos == last_pos for single-
                                // record blocks.
    n_records:        varint
    block_offset:     u64      // absolute file offset of the block's
                                // first byte (the block_header's first
                                // byte, not the column payloads)
```

Notes:

- Entries are listed in file order. A reader doing a region query
  binary-searches by `(chrom_id, first_pos, last_pos)`.
- `block_offset` is the **absolute** offset from the start of the
  file. We do not encode it as a delta from the previous entry: at
  ~16–24 bytes per entry and a few hundred KB of index even on a
  ×5 WGS, the size win is not worth the encode/decode complexity.
- `n_records` is redundant with reading the block header itself but
  cheap to keep; it lets a reader budget memory before decompressing
  and helps `psp head` summarise the file without scanning blocks.
- The index is **uncompressed**. Compressing it would save ~50% on a
  structure that is already a fraction of a percent of file size,
  and would force a reader to instantiate zstd before answering "what
  blocks cover chrN?". Not worth it.
- A `.psp` with zero data blocks still has a zero-entry index; the
  trailer's `index_offset` points to the byte immediately after the
  header (where blocks would have started).

### Why a tail index, not a head index or a sidecar

- **Streaming writer.** A head index would force two-pass write or
  back-patched offsets; tail-mounted means the writer accumulates
  index entries in memory and dumps them once at the end. Memory
  cost is bounded by the block count (a few thousand entries on a
  WGS-scale file).
- **Atomic with the file.** Sidecar indexes (`.bai`, `.crai`, `.tbi`)
  can drift relative to the data file when one is regenerated and
  the other isn't, or get lost in transit. A tail index travels
  with the bytes it indexes.
- **Standard pattern.** Parquet, ORC, Arrow IPC File, and ZIP all
  mount their metadata at the file tail and locate it via a fixed
  trailer. The shape below is the same.

### Size estimate

Roughly 16–24 bytes per entry (varints + one `u64`). At the writer's
default target block size (16 MiB uncompressed):

- 5× human exome (~100 MB compressed) → ~6 blocks → ~140 B index.
- 5× human WGS (~50 GB compressed) → ~3 k blocks → ~60 KB index.

Negligible against the data body.

## File trailer

A fixed 32-byte trailer at the very end of the file. Its job is
threefold: (a) locate the block index, (b) detect truncation /
incomplete writes that would otherwise read as a shorter but
valid-looking file, (c) detect silent corruption of the index
region itself (the one region whose payload is not already covered
by a zstd frame checksum).

| Offset (from trailer start) | Field | Type | Notes |
|---|---|---|---|
| 0 | `index_offset` | `u64` | Absolute file offset of the block index's first byte. |
| 8 | `index_byte_length` | `u64` | Length of the block index in bytes. `index_offset + index_byte_length` must equal the trailer's own start offset. |
| 16 | `n_blocks` | `u64` | Total number of blocks written. Must equal the number of entries the reader decoded from the index. |
| 24 | `index_checksum` | `u32` | XXH3-64 hash of the `index_byte_length` index bytes, truncated to its low 32 bits and stored little-endian. Same hash and same truncation that zstd uses for its frame content checksum (RFC 8878 §3.1.1), so the same XXH3 implementation already pulled in by the zstd dependency computes both. Covers the index region — the only region in the file that is not already covered by a zstd frame's built-in checksum. For an empty index (zero blocks) the value is the XXH3-64 hash of the empty byte sequence, truncated. |
| 28 | `trailer_magic` | 4 B | ASCII bytes `'P'`, `'S'`, `'P'`, `'E'` (for "PSP End"). Different from the head magic so truncation that copies the head pattern would not pass the tail check. Placed last so the EOF-aligned `pread` of the last four bytes is enough to reject a non-`.psp` or truncated file. |

A reader's tail-read protocol:

1. `pread` the last 32 bytes; validate `trailer_magic`.
2. Seek to `index_offset`; read `index_byte_length` bytes.
3. Recompute the XXH3-64 of the read bytes, truncate to its low 32
   bits, and compare to `index_checksum`. Mismatch is a hard
   corruption error — no fall-back, no partial reads. The
   architecture's principle is that errors do not pass silently
   (`design_principles.md` §3); the index is the one file region
   whose silent corruption could redirect every subsequent block
   read to wrong bytes, so it is checked before any block offset
   from it is trusted.
4. Decode the index; assert the decoded entry count equals
   `n_blocks`.
5. Region queries binary-search the in-memory index by
   `(chrom_id, first_pos, last_pos)`; sequential reads ignore the
   index entirely and stream from block 0.

A reader that reaches EOF without seeing the trailer, whose trailer
arithmetic disagrees with `n_blocks`, or whose recomputed
`index_checksum` disagrees with the stored value, reports a hard
truncation / corruption error. There is no auto-recovery.

## Versioning policy (current draft)

The version pair lives in the TOML header as
`format-version = "MAJOR.MINOR"` (string).

**Ordering is numeric, not lexical.** A reader parses the string by
splitting on the single `.`, converting each half to a non-negative
integer with no leading-zero allowed (except for the literal `0`),
and comparing those integers as numbers. So `"1.10"` is greater than
`"1.9"`, even though lexically it would sort the other way. Each
integer half fits in a `u16`; values outside `[0, 65535]` are a hard
error. The per-field rule at §"Per-field character set" already
restricts the surface syntax (`[0-9]+\.[0-9]+`); this paragraph is
what gives those characters their order.

Reader behaviour:

- **Same major, same or lower minor.** Reader fully understands the
  file. Unknown header keys are skipped (TOML schema is open by
  convention here); block-level column dispatch goes by tag; required
  unknown tags are an error (forbidden in a properly-issued minor
  bump), optional unknown tags are skipped.
- **Same major, higher minor than the reader supports.** Reader
  consumes the file by tag dispatch. Unknown header keys are
  skipped. Required column tags the reader does not know abort
  with a clear error — but minor bumps are not allowed to add
  required content (§"File-scope extensions (deferred to v1.x)"),
  so a properly-issued v1.x file is always consumable by a v1.0
  reader. Optional unknowns are silently skipped.
- **Higher major.** Reader aborts immediately with "this file was
  written by a newer format major version; please update". No
  attempt to parse.
- **Lower major (downgrade).** Reader aborts immediately.

A writer always emits the *lowest* `(major, minor)` that covers the
content it wrote. A v1.1 writer that happens to emit no v1.1-only
columns produces a v1.0 file. This keeps the producer side honest
about which readers can consume the artefact.

## Non-goals

- **No partial-file recovery.** A truncated file is an error, not
  a recoverable state.
- **No BCF/VCF interoperability at the byte level.** Optional
  exporters live in their own module and do not constrain the on-disk
  format (architecture doc §"Why custom binary rather than BCF").
- **No backward compatibility for the v0 PSP draft** that this
  document supersedes; nothing was ever written in the PSP
  layout, so there is no legacy to read.

## Appendix — Design decisions log

This appendix preserves the reasoning behind each major design
decision so future readers can understand *why* the format looks
the way it does, not just *what* it requires. Each entry was once
an open question during the design discussion; the resolution is
recorded here. Two groups: **forward-compatibility** (the shape of
the versioning / extensibility story) and **parameters & layout
details** (defaults and field widths).

### Forward-compatibility

**~~Q-IDX — should the format carry a random-access index?~~**
*Resolved 2026-05-13: mandatory tail-mounted block index (option 2),
specified in §"Block index".*

**~~Q-HDR — should the file header be plain text rather than binary?~~**
*Resolved 2026-05-13: yes — TOML body framed by `PSP\n` magic and
`---END-HEADER---\n` sentinel, ASCII per SAM v1.6 §1.2.1 conventions
for v1.0 values. Specified in §"File header" and §"Encoding
conventions / Header (plain-text TOML)". The framing was later
revised in Q-HDR2 to add an authoritative byte-count length prefix
between the magic and the TOML body.*

**~~Q-HDR2 — how is the end of the TOML body located: by sentinel
scan, by length prefix, or both?~~**
*Resolved 2026-05-13: an 8-byte little-endian `u64` length prefix
between `PSP\n` and the TOML body is authoritative for the body
boundary; the `---END-HEADER---\n` sentinel is kept as a structural
cross-check that the reader verifies immediately after consuming
the body. Specified in §"File header" / "Layout" and "Reader
protocol".*

*The earlier sentinel-only framing (Q-HDR original resolution)
relies on the byte sequence `\n---END-HEADER---\n` not appearing
inside any TOML value. The per-field rules forbid `\n` in every
v1.0 field, so well-formed files are safe — but TOML v1.0.0
multi-line basic strings (`"""..."""`) and literal strings
(`'''...'''`) legally contain `\n`, and a producer that fails the
field-rule validation could silently emit a file whose sentinel
scan ends in the wrong place. The "errors must not pass silently"
posture (`design_principles.md` §3) doesn't tolerate a framing
that depends on the writer being bug-free to find the writer's
bug; the length prefix moves the body boundary from "wherever the
reader finds the pattern" to "wherever the writer says the body
ends," and a mismatch with the sentinel becomes a hard error
rather than a silent prefix-eat.*

*Alternatives considered: (a) drop the sentinel entirely once the
length prefix is in place — rejected because `head file.psp` would
no longer end on a recognisable line and a corrupted length prefix
could silently truncate the body without any cross-check firing;
(b) use an ASCII-decimal length prefix on its own line for
`head`-cleanliness — rejected because it puts framing data
syntactically inside the TOML region (a TOML parser would parse
the framing line as a key), and the binary u64 form is shorter,
unambiguous about leading zeros, signedness, and base, and trivial
to read with a single fixed-offset `pread`. The 8 bytes of binary
between `PSP\n` and the TOML body are mostly null (the body length
fits in 3 bytes for any realistic header) and terminals handle the
nulls cleanly; `psp head` is the user-facing inspection tool, raw
`head` is the fallback. Writer cost is zero beyond what was
already paid: the writer assembles the TOML body in memory to
validate it (§"Per-field character set"), and the validated
buffer's length is the prefix value — no back-patching, no
two-pass write.*

**~~Q-FC1 — versioning granularity: file-level only, or also
per-column?~~**
*Resolved 2026-05-13: file-level only. The `format-version` in
the TOML header is the single version axis. Per-column versioning
adds machinery for a problem that can be addressed more simply by
minting a new column tag when a column's semantics change (the
Q-FC3 "Option 2" path). Knock-on: drop option 3 from Q-FC3.*

**~~Q-FC2 — should the file header carry a "columns used" manifest
in addition to per-block manifests?~~**
*Resolved 2026-05-13: yes, and the manifest carries the full schema
(VCF-style), not just a list of names. Specified in §"File header
— binary schema". The header's `[[column]]` array is the file's
single source of truth for column semantics; the per-block manifest
references it by tag and carries only the per-block-varying fields
(tag, compressed_len, uncompressed_len).*

**~~Q-FC3 — how is a "semantic change to an existing column" expressed?~~**
*Deferred 2026-05-13: no v1.0 commitment. The case ("we changed the
meaning of an existing column's bytes without changing their type")
is hypothetical and may never arise. When/if it does, the policy
will be decided then; the choice is between a major bump and
minting a new column tag (the option 3 "per-column version flag"
path is closed by the Q-FC1 resolution).*

**~~Q-FC4 — should minor bumps be allowed to promote an
optional tag to required?~~**
*Resolved 2026-05-13 (alongside the capability-list YAGNI cleanup):
no. Minor bumps add only optional content; promoting optional to
required is a major bump. Stated as the load-bearing "minor-bump
rule" in §"File-scope extensions (deferred to v1.x)".*

**~~Q-FC5 — should a higher-minor reader be required to roundtrip a
lower-minor file losslessly?~~**
*Deferred 2026-05-13: no v1.0 commitment. No repacking / relocation
tool exists or is planned; when one does, the requirement is easy
to retrofit (it's a property of the tool, not the format).*

**~~Q-FC6 — capability list: dependencies between capabilities?~~**
*Moot 2026-05-13: v1.0 has no capability list (YAGNI; see §"File-
scope extensions (deferred to v1.x)"). If a future minor bump
introduces inter-extension dependencies, the question reopens then.*

**~~Q-FC7 — what reserves do we keep for "we ran out of room"?~~**
*YAGNI 2026-05-13. The column-tag-range reserves already in §"Reserved
column-tag ranges" are kept; no new reserved-flag bits are pre-baked.
The per-block manifest's `flags` byte was dropped entirely when the
binary schema moved to the header (Q-FC2 resolution), so the question
about reserving bits in it is moot. Varint widths grow indefinitely
as needed. We'll add reserves when a concrete need appears.*

**~~Q-FC8 — how is phase-chain state preserved across block boundaries
for random access?~~**
*Originally resolved 2026-05-13 with a per-block
`active_chain_ids_at_block_start` snapshot, on the assumption that
chain ids needed recycling to keep the namespace small. **Superseded
2026-05-14**: switched to **unique-per-`.psp`-file `u64` chain ids**
(alternative (b) above, originally rejected on `u16` namespace
grounds; with `u64` the namespace concern is gone). Under the new
design there is no active-set state to preserve across blocks —
chain ids are unique throughout the file — so the snapshot field
is gone from the block header and inter-block continuity checks
have been retired. The full rationale and migration is in
[ia/feature_implementation_plans/unique_chain_ids.md](../feature_implementation_plans/unique_chain_ids.md).*

### Parameters & layout details

**~~Q-PL1 — default zstd compression level.~~**
*Resolved 2026-05-13: writer hardcodes zstd level 9. Not CLI-
exposed in v1 (too much detail for the standard user). The format
does not record the level; readers decompress any level. A future
benchmark on real data may change the writer's default — that is a
writer change, not a format change. See §"Compression".*

**~~Q-PL2 — default target uncompressed block size.~~**
*Resolved 2026-05-13: writer hardcodes 16 MiB. Files are expected
to be large, so ratio matters more than per-block latency; bigger
blocks compress better. Not CLI-exposed in v1; not recorded in the
format. A future benchmark on real data may revise the default —
that is a writer change, not a format change. See §"Block sizing".*

**~~Q-PL3 — `f32` vs `f64` for `allele_q_sum_log`.~~**
*Resolved 2026-05-13: `f64`, matching the implementation
([pileup/mod.rs:357](../../src/per_sample_pileup/pileup/mod.rs#L357),
`AlleleSupportStats::q_sum: f64`). Storing at the accumulator's
native precision avoids a truncation step at the I/O boundary; the
extra 4 B per allele is negligible against the rest of the per-
allele payload.*

**~~Q-PL4 — `u8` vs `u16` for phase-chain ids.~~**
*Originally resolved 2026-05-13: `u16` over `u8` for headroom under
recycling. **Superseded 2026-05-14**: phase-chain ids are now
**unique-per-`.psp`-file `u64`** (`pub type ChainId = u64` in
[pileup/chain_id_allocator.rs](../../src/per_sample_pileup/pileup/chain_id_allocator.rs)),
with no recycling. The on-disk encoding is fixed-width LE `u64`
under tag `0x22` (`allele-chain-ids`); zstd absorbs the leading
zero bytes on small ids. Overflow at `u64::MAX + 1` surfaces as
`WalkerError::ChainIdSpaceExhausted` rather than silent wrap; the
ceiling (~1.8 × 10¹⁹ chains per file) is unreachable on any
realistic workload. See
[ia/feature_implementation_plans/unique_chain_ids.md](../feature_implementation_plans/unique_chain_ids.md).*

**~~Q-PL5 — reserve a v1.x optional column for 2-bit packed allele
sequences?~~**
*Resolved 2026-05-13: no. Trust zstd. The size win over uppercase-
ASCII bytes is expected to be minimal once zstd has absorbed the
A/C/G/T redundancy across adjacent records, and a packed encoding
would add a parallel writer/reader path for negligible gain.*

**~~Q-PL6 — writer constraint on column ordering within a block.~~**
*Resolved 2026-05-13: strict requirement. Per-block manifest entries
are tag-ascending with no duplicates; writers sort at emit time,
readers reject violations. See §"Column ordering within a block".*

**~~Q-PL7 — trailer shape.~~**
*Resolved alongside Q-IDX: trailer is now 32 B carrying
`index_offset`, `index_byte_length`, `n_blocks`, `index_checksum`,
and `trailer_magic`. (The original draft had a 4 B reserved padding
slot in place of `index_checksum`; that slot was used up by Q-PL11.)
See §"File trailer".*

**~~Q-PL8 — per-block checksum.~~**
*Resolved 2026-05-13: zstd's built-in per-frame checksum is
sufficient; no extra block-level CRC. See §"Compression".*

**~~Q-PL9 — float byte-layout details: NaN handling.~~**
*Resolved 2026-05-13: hard error. NaN (or ±∞) in any float column
is a producer bug; readers reject the file with a clear message
naming the column, the block, and the offending entry. Applies to
every IEEE 754 float column the format ever carries, not just
`allele-q-sum-log`. Consistent with `design_principles.md`
§"Errors must not pass silently".*

**~~Q-PL11 — integrity coverage for the block-index region.~~**
*Resolved 2026-05-13: the trailer carries a 32-bit `index_checksum`
(XXH3-64 of the index bytes, truncated to its low 32 bits — the
same hash and truncation zstd uses for its frame content checksum).
Reason: zstd frame checksums cover every column payload (§"Compression"),
the trailer's `trailer_magic` + `n_blocks` arithmetic covers EOF
truncation, and the file header is plain TOML whose corruption
surfaces as a parse failure. The block index was the one file region
left uncovered — bytes that are uncompressed varints plus absolute
`u64` block offsets, sitting between the data and the trailer. A
single flipped bit in an offset would silently redirect a region
query to the wrong block, where decompression succeeds (the wrong
column payloads still pass their own zstd checksums) and the reader
returns wrong genotypes downstream. That is precisely the
"errors-must-not-pass-silently" failure the project bans
(`design_principles.md` §3), so the gap has to be closed.
Alternatives considered: (a) CRC32C — equivalent error-detection
strength but a second hash family to bring into the codebase
(rejected: zstd already ships XXH3); (b) full 64-bit checksum
(rejected: 32-bit is the ecosystem convention — zstd, BGZF,
Parquet, ZIP all use 32 bits for region-scale checksums, the
collision probability against accidental corruption is 1 in 2³²,
and the saved 4 bytes are not the relevant cost). The 4 bytes
this checksum occupies came from the trailer's previous
`trailer_padding` field — exactly what that field's documented
forward-compat purpose was reserved for ("a tiny tail-side feature,
e.g. file-scoped checksum"), used up now rather than later. Future
tail-side features past v1.0 will need their own room.*

**~~Q-PL10 — reference identity: per-file or per-contig hash?~~**
*Resolved 2026-05-13: per-contig `md5` on each `[[chromosome]]`
entry, with SAM `@SQ M5` semantics (MD5 of the uppercase reference
sequence for that contig). An earlier draft used a top-level
`reference-md5` of the literal FASTA file bytes; that turned out
to reject the same reference re-wrapped, regenerated, or
re-downloaded from a different mirror — i.e. it fired on FASTA-file
identity rather than reference-sequence identity, which is the
opposite of what the consistency check is meant to enforce. The
per-contig `@SQ M5` form has been the GA4GH / refget / samtools
canonical reference hash since long before this project existed and
is what every other tool in the pipeline (CRAM in particular)
already keys reference identity on. Adopting it makes Stage 3's
mismatch error specific (the offending contig, by name) rather
than whole-file, and aligns this format with the broader
ecosystem at no extra cost beyond ~50 bytes per contig in the
TOML header. See §"File header" / `[[chromosome]]` `md5` key.*
