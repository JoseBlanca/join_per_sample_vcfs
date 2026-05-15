# `pop_var_caller psp-to-pileup` — `.psp` → text pileup dump

Implementation plan for a small read-side utility: stream a `.psp`
artefact through the PSP reader and emit a human-readable, mostly
`samtools mpileup`-compatible text dump with extra trailing columns
for PSP-specific per-allele aggregates. Intended for inspection,
debugging, comparison against `samtools mpileup` runs on the same
input, and quick eyeballing of regions.

Spec/context:
[per_sample_pileup_format.md](../specs/per_sample_pileup_format.md)
(the `.psp` byte layout), and the existing `pop_var_caller pileup`
subcommand
([cli.rs](../../src/pop_var_caller/cli.rs)) which produces the
`.psp` artefacts this utility consumes.

## Why this slice exists

Stage 1 produces `.psp` files but they're binary and only readable
via Rust calls to `PspReader`. There is no command-line way to
look at what a `.psp` actually contains. This utility closes that
gap with a text-mode dump similar enough to `samtools mpileup`
that existing bioinformatics tooling (awk pipelines, comparison
scripts, eyeballing in `less`) works on the output.

Read-side surface that already exists and we just orchestrate:

- [`PspReader::new`](../../src/per_sample_caller/psp/reader.rs#L132)
  — opens a `.psp`, validates header/trailer/index.
- [`PspReader::header`](../../src/per_sample_caller/psp/reader.rs#L304)
  — returns the parsed `[[chromosome]]` table.
- [`PspReader::records`](../../src/per_sample_caller/psp/reader.rs#L328)
  — sequential iterator over every `PileupRecord`.
- [`PspReader::region_records`](../../src/per_sample_caller/psp/reader.rs)
  — coordinate-clamped iterator (used by `--region`).
- `PileupRecord` carries `alleles[0] = REF` invariant
  ([pileup/mod.rs:394](../../src/per_sample_caller/pileup/mod.rs#L394))
  — REF base is derivable from the record, no FASTA needed.

## Scope

In scope:

- A new subcommand `pop_var_caller psp-to-pileup` that takes a
  `.psp` and writes text records to a file or stdout. **Stdout
  output is the natural default here** (text is line-streamable,
  unlike the random-access `.psp` writer path).
- A 7-column tab-delimited line format: the 6 classical
  `samtools mpileup` columns, plus a trailing column with per-allele
  PSP aggregates encoded in a parseable form. See "Output format"
  below for the column contract.
- Indel encoding matches `samtools mpileup`: insertions emit
  `<anchor>+N<bases>` per supporting read; deletions emit
  `<anchor>-N<refbases>` per supporting read. Reference bases for
  the deletion are reconstructed from `alleles[0].seq` — no FASTA
  required.
- A `--region chrom[:start-end]` flag that routes through
  `PspReader::region_records` for fast spatial filtering.
- Toggles to include or omit the optional trailing PSP-specific
  column(s): `--show-chain-ids` flips on per-allele
  `chain_slots` dumping (off by default — chain ids are noisy and
  most consumers don't need them).
- Output is **text streamed line-by-line** — no atomic
  tmp+rename, no random access on read or write. A killed run
  leaves a partially written file by design (matches every
  text-dump tool's behaviour, including `samtools mpileup`).

Deferred to a follow-up slice:

- `--bed FILE` for many-region filtering (region list). The
  one-region path lands first; a list is mechanically the same
  loop over `region_records`.
- `--per-allele-row` mode (a true long-format companion: one row
  per AlleleObservation rather than one row per position). Useful
  for analytics; the per-position dump is the common case.

Not planned at all:

- **Per-base quality reconstruction in column 6.** PSP aggregates
  per-base BQ into per-allele `q_sum`; reconstructing per-read
  qualities would be a lie. Column 6 is emitted with a fixed
  placeholder character (see "BQ column" below) so mpileup
  parsers that *demand* six columns still work, but the user
  cannot mistake the value for real per-base BQ.
- **A `--reference FASTA` flag.** PSP's REF invariant
  (`alleles[0]` is always REF) means we can produce every column
  without consulting a FASTA. Adding `--reference` would invite
  unverifiable cross-checks against a FASTA that may have changed.
- **`bgzf` or `.tsv.gz` compressed output.** The whole point of
  this utility is "human-readable text". Compression is the
  user's job (pipe through `gzip`/`bgzip`).
- **Cross-sample joining.** Strictly one input `.psp`, one
  output stream. Multi-sample work is Stage 2+.

## Module layout

```
src/
  main.rs                         — extend the Subcommand enum
  pop_var_caller/
    mod.rs                        — pub mod psp_to_pileup
    cli.rs                        — unchanged
    cli/
      error_bridge.rs             — unchanged
    psp_to_pileup.rs              — NEW: args, error, run, encoder
```

The utility is small enough to fit in a single file. Encoder
helpers are private functions; only the args struct, error type,
and `run` function are public.

## Output format

One line per `PileupRecord` (i.e. per covered reference position
with non-zero coverage). Tab-delimited. Columns:

| # | Name              | Source                                                                   |
|---|-------------------|--------------------------------------------------------------------------|
| 1 | `CHROM`           | `ParsedHeader.chromosomes[record.chrom_id].name`                         |
| 2 | `POS`             | `record.pos` (1-based, as PSP stores it)                                 |
| 3 | `REF`             | First base of `alleles[0].seq` (the anchor REF base)                     |
| 4 | `DEPTH`           | `Σ allele.support.num_obs` across all alleles in the record              |
| 5 | `BASES`           | samtools-mpileup-style read-bases string (see "Read-bases encoding")     |
| 6 | `QUALS`           | `'!'` × `DEPTH` — placeholder; PSP does not retain per-base BQ           |
| 7 | `ALLELE_DETAILS`  | per-allele aggregates (see "ALLELE_DETAILS encoding")                    |

### Read-bases encoding (column 5)

For each allele in `record.alleles`, in PSP order (REF first), append:

- **REF (`alleles[0]`, with `seq` exactly one base):**
  `'.'` × `fwd` followed by `','` × `(num_obs - fwd)`.

- **REF column on a multi-base REF span (`ref_span > 1`):**
  Same as the single-base case for the anchor base — `.` / `,`.
  The deleted ref bases live on the deleting allele's `-N`
  marker, not on the REF column. PSP's REF-on-DEL-anchor entry
  represents reads that did *not* delete, so they get `.`/`,`.

- **SNP ALT (`allele.seq.len() == 1` and != REF base):**
  uppercase `allele.seq[0]` × `fwd` + lowercase × `(num_obs - fwd)`.

- **MNP ALT (`allele.seq.len() == ref_span > 1` and != REF):**
  Same as SNP but using `allele.seq[0]` (the anchor mismatch).
  PSP does not retain per-base direction within a multi-base
  allele; we emit only the anchor character. Documented in the
  `--help` text.

- **Insertion-bearing ALT (`allele.seq.len() > ref_span`):**
  For each supporting read (`num_obs` total, of which `fwd` are
  forward):
  - Anchor character: `.` if forward, `,` if reverse (the anchor
    base matches REF).
  - Followed immediately by `+N<inserted-bases>` where
    `N = inserted_bases.len()` and `inserted-bases =
    allele.seq[ref_span..]`. Uppercase for forward reads,
    lowercase for reverse — samtools' convention.

- **Deletion-bearing ALT (`allele.seq.len() < ref_span`):**
  For each supporting read:
  - Anchor character (`.` / `,`).
  - Followed by `-N<deleted-ref-bases>` where the deleted
    bases are `REF.seq[allele.seq.len()..ref_span]`. Uppercase
    for forward reads, lowercase for reverse.

The ordering between REF-block characters and ALT-block characters
in column 5 is deterministic (PSP allele order). It is **not** the
samtools-mpileup convention of interleaving by read order — PSP
does not retain read order. The `--help` text documents this so
users diffing against `samtools mpileup` know to sort each
column 5 character-wise before comparison.

### `ALLELE_DETAILS` encoding (column 7)

Comma-separated list of allele records, in PSP order (REF first).
Each allele record is

```
<seq>:<num_obs>:<fwd>:<placed_left>:<placed_start>:<q_sum>
```

- `seq`: the allele bases (ASCII over `{A,C,G,T,N}`), or the
  literal `*` if the allele's seq is empty (pure deletion to
  zero bases).
- `q_sum`: formatted with `{:.3}` (three decimal places) — enough
  for downstream sorting/clustering, tight enough that the column
  doesn't bloat on synthetic samples.

A line for a clean ref-only SNP position thus looks like:

```
chr1<TAB>1234<TAB>A<TAB>30<TAB>...........,,,,,,,,,,<TAB>!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!<TAB>A:30:20:5:3:-12.345
```

`--show-chain-ids` appends one more colon-delimited field per
allele record:

```
<seq>:<num_obs>:<fwd>:<placed_left>:<placed_start>:<q_sum>:<chain_ids>
```

with `<chain_ids>` as a `;`-separated list of `u64`s (semicolon
because comma is the allele-record separator). When the list is
empty, the field is empty (`...:`).

### BQ column placeholder

Column 6 is filled with `'!'` (ASCII 33 = Phred 0). Rationale:

- `'!'` is the standard Phred+33 encoding of quality 0 and is
  recognised by every mpileup parser as "missing quality".
- It is **not** the same character as the column-5 deletion-event
  marker `'*'` (which appears in canonical mpileup when a base
  at this column was deleted by an earlier event). Keeping them
  visually distinct avoids confusing column 5 with column 6.
- A single sentinel — rather than a per-allele computed
  approximation from `q_sum / num_obs` — keeps the lossy
  aggregation explicit. Anyone reading column 6 of this output
  is immediately told "PSP doesn't have this; ignore it".

## CLI surface (clap-derive)

Subcommand added to the existing `Subcommand` enum in
[src/main.rs](../../src/main.rs):

```rust
#[derive(clap::Subcommand)]
enum Subcommand {
    /// Stage 1: run BAQ + pileup over one sample's CRAMs, emit a .psp.
    Pileup(PileupArgs),

    /// Stream a .psp as samtools-mpileup-style text (plus PSP per-allele aggregates).
    PspToPileup(PspToPileupArgs),
}
```

Args struct (lives in `src/pop_var_caller/psp_to_pileup.rs`):

```rust
#[derive(clap::Args)]
pub struct PspToPileupArgs {
    /// Input .psp file.
    #[arg(long)]
    pub input: PathBuf,

    /// Output path. Use `-` (or omit) to stream to stdout.
    #[arg(long, default_value = "-")]
    pub output: PathBuf,

    /// Restrict output to a single region: `chrom` or
    /// `chrom:start-end` (1-based, inclusive on both ends).
    /// Without this flag the whole .psp is streamed.
    #[arg(long)]
    pub region: Option<String>,

    /// Append chain_slot ids to each allele record in column 7
    /// (semicolon-separated, after q_sum).
    #[arg(long)]
    pub show_chain_ids: bool,
}
```

This is the entire surface — five flags including positional
inputs. No `--threads` (single-threaded by design), no
`--reference` (PSP carries REF), no `--format` (output format is
fixed and documented).

## `run_psp_to_pileup` — the orchestrator

```rust
pub fn run_psp_to_pileup(args: &PspToPileupArgs) -> Result<(), PspToPileupError>;
```

Behaviour:

1. Open `PspReader::new(BufReader::with_capacity(64 * 1024,
   File::open(&args.input)?))`. Extract the parsed header so we
   can map `chrom_id → name`.
2. Open the sink:
   - `args.output == "-"` → `Box<dyn Write> = Box::new(io::stdout().lock())`.
   - Otherwise → `Box::new(BufWriter::new(File::create(&args.output)?))`.
   (Buffered because we emit many small writes per line.)
3. Build the chromosome-id → name lookup table from
   `reader.header().chromosomes`.
4. If `--region` is set, parse it (a small free-function helper),
   resolve the chrom name to a `chrom_id` via the lookup, and
   pull records from `reader.region_records(chrom_id, start, end)`.
   Otherwise pull from `reader.records()`.
5. For each `Result<PileupRecord, PspReadError>` yielded:
   - On error: short-circuit, surface as `PspToPileupError::Read(e)`.
   - On Ok: format the seven columns into a `String` via
     `emit_line(&mut buf, record, chrom_name, args.show_chain_ids)?`
     and write `buf` to the sink, followed by `\n`. Reuse the
     buffer across iterations (one `String::clear` per record).
6. On end-of-iterator: `sink.flush()?` and return `Ok(())`.

The sink is *not* synced to disk — text output is line-streamed
and the buffer flush is sufficient. If a downstream consumer reads
from a pipe, it sees every flushed line; if a file is the sink,
the OS commits as usual on close.

### Region parsing

```rust
fn parse_region(s: &str) -> Result<RegionSpec, PspToPileupError>;

struct RegionSpec {
    chrom: String,
    start: Option<u32>,  // 1-based, inclusive; None ⇒ 1
    end: Option<u32>,    // 1-based, inclusive; None ⇒ u32::MAX
}
```

Accepted forms:
- `chr1` — whole chromosome.
- `chr1:1000` — from position 1000 to chromosome end.
- `chr1:1000-2000` — from 1000 to 2000 inclusive.

Anything else (multiple colons, non-numeric range, end < start) is
a hard error with the parsed-input context in the message.

### Error type

```rust
#[derive(thiserror::Error, Debug)]
pub enum PspToPileupError {
    #[error("PSP read: {0}")]   Read(#[from] PspReadError),
    #[error("io: {0}")]         Io(#[from] io::Error),
    #[error("region: {0}")]     Region(String),
    #[error("chromosome '{0}' not found in input .psp")]
    UnknownChromosome(String),
}
```

## Test strategy

### Unit tests (in `psp_to_pileup.rs`)

- `parse_region` happy paths (`chr1`, `chr1:1000`, `chr1:1000-2000`).
- `parse_region` rejects: empty string, multiple colons, non-numeric
  range, end < start, position 0.
- `emit_line` SNP-only, all-forward: assert column 5 is `.` × n.
- `emit_line` SNP-only, mixed strand: assert column 5 is
  `.` × fwd + `,` × (obs - fwd).
- `emit_line` with a SNP ALT: column 5 contains the right mix of
  REF chars and uppercase/lowercase ALT chars.
- `emit_line` with an insertion ALT: column 5 contains
  `.+N<bases>` / `,+N<bases>` per supporting read.
- `emit_line` with a deletion ALT: column 5 contains
  `.-N<refbases>` / `,-N<refbases>` per supporting read; the
  deleted ref bases come from `REF.seq[allele.seq.len()..]`.
- `emit_line` column 6 is `'!'` × `DEPTH`.
- `emit_line` column 7 with `show_chain_ids = false` omits the
  chain field.
- `emit_line` column 7 with `show_chain_ids = true` appends
  `;`-separated chain ids per allele record.

All `emit_line` tests build `PileupRecord` literals via the
existing `PileupRecord::new` / `AlleleObservation::new` /
`AlleleSupportStats::new` constructors — no in-memory `.psp`
needed at this layer.

### Integration test (`tests/psp_to_pileup_integration.rs`)

Two cases. Both build a small in-memory `.psp` via the existing
`PspWriter` + a few hand-built `PileupRecord`s (no CRAM round-trip
needed; this is a read-side tool and the input contract is
`PspReader`-shaped):

1. **Full file, no region.** Write a 3-record `.psp` (one
   SNP-only, one with an insertion, one with a deletion), call
   `run_psp_to_pileup` with `args.output` pointing at a tempfile,
   read it back, assert exact line-by-line content.

2. **Region clamp.** Same `.psp`, this time
   `args.region = Some("chr1:2-2".into())`, assert only the second
   record's line is emitted.

A third case using the end-to-end binary (spawning the actual
`pop_var_caller psp-to-pileup` from a unit test via
`std::process::Command::new(env!("CARGO_BIN_EXE_pop_var_caller"))`)
is **deferred**; the in-process integration test above already
exercises the same code path.

### Manual smoke

After landing:

```
$ cargo run --release --bin pop_var_caller -- psp-to-pileup \
    --input sample.psp \
    | head -20
```

against a real `.psp` produced by `pop_var_caller pileup`. Expected:
20 well-formed lines, each with 7 tab-separated columns; the
ALLELE_DETAILS column parses as comma-separated colon-tuples.

## Risks and trade-offs

- **Column 5 ordering is not samtools-faithful.** Samtools'
  column 5 interleaves bases by read order in the BAM; PSP
  doesn't retain read order. We emit per-allele runs (all REF
  chars first, then ALTs in PSP order). Tools that diff against
  `samtools mpileup` need to sort each column 5 character-wise
  before comparing. Documented in `--help` and in the `ALLELE_DETAILS`
  description.
- **MNP encoding loses per-position direction.** A multi-base
  ALT (`allele.seq.len() == ref_span > 1`) emits only the anchor
  mismatch character on column 5. PSP itself does not retain
  per-position direction inside a multi-base allele, so any
  attempt at "expand MNP into N adjacent columns" would
  manufacture per-position data we do not have. The
  ALLELE_DETAILS column carries the full `seq`.
- **`'!'` for column 6.** Honest about PSP's lossiness; some
  third-party parsers may treat `'!'` as "real low quality" and
  filter it. Help text spells out the choice.
- **No atomic write.** Text dump utilities never atomic-write.
  A killed process leaves a partial file; rerun overwrites it.
  Matches `samtools mpileup`, `bcftools view`, `awk`, etc.
- **Stdout default.** Unlike Stage 1's `.psp` output (always a
  file, never stdout), this tool defaults to stdout because text
  is naturally pipe-consumable. Documented in `--help`.

## Validation

The slice is done when:

1. `cargo test --test psp_to_pileup_integration` passes.
2. Full test suite passes (`cargo test`).
3. Manual smoke run produces a 7-column tab-delimited file whose
   column count is constant and whose `ALLELE_DETAILS` column
   parses round-trip in a quick Python one-liner.
4. `cargo clippy --workspace --all-targets -- -D warnings` clean.
5. `cargo fmt -- --check` clean.
6. An implementation report lands under
   [ia/reports/implementations/](../reports/implementations/)
   noting any column-format adjustments made during
   implementation.

## File touch list

New:

- `src/pop_var_caller/psp_to_pileup.rs` — `PspToPileupArgs`,
  `PspToPileupError`, `run_psp_to_pileup`, `emit_line`,
  `parse_region`, and the unit tests (~350 lines).
- `tests/psp_to_pileup_integration.rs` — two-case integration
  test (~120 lines).
- `ia/reports/implementations/psp_to_pileup_<date>.md` —
  implementation report once landed.

Modified:

- `src/pop_var_caller/mod.rs` — add
  `pub mod psp_to_pileup; pub use psp_to_pileup::{run_psp_to_pileup,
  PspToPileupArgs, PspToPileupError};`.
- `src/main.rs` — extend the `Subcommand` enum with the
  `PspToPileup(PspToPileupArgs)` variant and route to
  `run_psp_to_pileup` in the dispatch match.

Unmodified:

- All `per_sample_caller` modules. The reader's existing public
  API (`PspReader::new`, `records`, `region_records`, `header`)
  is sufficient.

## Open questions for future slices

These are recorded so they aren't forgotten, but they do not block
this slice.

- A `--per-allele-row` flag flipping the output from
  one-line-per-position to one-line-per-allele. Mechanically
  trivial on top of `emit_line`; mostly a question of UX.
- A `--bed FILE` flag for many-region filtering. Adds nothing
  surprising over the existing region path; loop `read_records`
  with successive `region_records` calls, skip records already
  emitted.
- A future paired tool `pileup-to-psp` (PSP → samtools-mpileup-style
  bytes) is **not** something we expect to want — PSP is the
  source of truth, mpileup the lossy projection. Listed here so
  the asymmetry is on record.
