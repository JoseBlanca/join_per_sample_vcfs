# Cohort VCF writer (Stage 6 sink)

Implementation plan for the final stage of the pipeline: a streaming
writer that consumes
[`PosteriorRecord`](../../src/var_calling/posterior_engine.rs)
items from the posterior engine and emits a multi-sample VCF (plain
text or bgzipped) to disk.

Spec sections:
[calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)
§"Stage 6 — posterior engine" (record shape, posteriors, QUAL/GQ
semantics) and §"Stage 5 — per-group processing" step 3
(chain-anchor-broken compound calls → `CA` flag).

Background:
[freebayes_posterior_gt_probs.md](../specs/freebayes_posterior_gt_probs.md)
for QUAL / FORMAT conventions, and
[per_sample_pileup_format.md](../specs/per_sample_pileup_format.md)
§"Chromosome table" for the contig metadata shape (name, length,
md5) that flows forward into VCF `##contig=` lines.

This is glue + format encoding, not new algorithms. Every quantity
the writer needs is already on `PosteriorRecord`; the work is
mapping it into the right VCF cells and producing valid VCF 4.4
bytes.

## Scope

In scope:

- A `CohortVcfWriter` that streams `PosteriorRecord` →
  one VCF data line per record, plus a header built from cohort
  metadata at construction time.
- Output sink dispatch: a `.vcf` path writes plain text; a
  `.vcf.gz` path writes **bgzf** (not plain gzip) so the file is
  htslib-compatible and indexable by `tabix` in a follow-up.
- VCF version: **4.4** (matches the field semantics we need: `AD`
  as `Number=R`, `GP` as `Number=G`, `<*>` allele not required).
- Header construction from a `CohortMetadata` input (sample names,
  contig table with `name`/`length`/`md5`, project tool/version
  string for `##source`, command-line string for
  `##commandline=...` provenance).
- Per-record encoding from `PosteriorRecord`:
  - CHROM / POS / ID (`.`) / REF / ALT / QUAL / FILTER (`PASS` by
    default) / INFO / FORMAT.
  - INFO: `AF`, `AC`, `AN`, `DP`, plus `CA` (flag) when any
    sample in the record has a chain-anchor-broken call.
  - FORMAT per sample: `GT:GQ:DP:AD` by default; `GT:GQ:DP:AD:GP`
    when `WriterConfig::emit_gp` is enabled. `GP` is opt-in because
    its size grows with allele count (~21 floats/sample at
    ploidy=2, n_alleles=6) and the typical consumer doesn't read
    it. See Risks.
- Atomic write: bytes go to `<output>.tmp`, `fs::rename` to the
  final path on `finish()`. (Same convention adopted in the Stage
  1 `pileup` CLI plan for `.psp`.)
- Backwards-emit safety: writer enforces that records arrive in
  non-decreasing `(chrom_id, pos)` order and surfaces an error on
  regression. (The upstream posterior engine already preserves
  Stage 4's order, but the writer is the last line of defence.)
- Unit tests on header round-trip and FORMAT-cell formatting.
  Integration test on a hand-built `Vec<PosteriorRecord>` fixture
  exercising every FORMAT field, the `CA` flag, and the
  multi-allelic case (≥ 3 alleles).

Deferred to a follow-up slice (intended to land later, not in this
cut):

- **Tabix `.tbi` index alongside `.vcf.gz`.** The bgzf format
  carries the virtual offsets a tabix builder needs; we just don't
  build the index in v1. Easy to slot in via `noodles-tabix` later.
- **A `--region` reader-side feature.** This is a writer; regions
  are imposed upstream (CLI cohort flag), not here.
- **`PL` (phred-scaled per-genotype likelihoods).**
  `PosteriorRecord` does not currently carry the pre-EM
  log-likelihoods (only the posteriors). Adding `PL` means either
  forwarding `MergedRecord.log_likelihoods` through the posterior
  engine into a new `PosteriorRecord.log_likelihoods: Vec<f64>`
  field, or recomputing them on the fly. Out of scope here;
  separate slice.
- **Per-sample contamination fraction in INFO.** The Stage 6
  contamination side-pass writes `ContaminationEstimates` to its
  own JSON sidecar; the spec does not require it in the VCF.
  Revisit if a downstream consumer needs it.
- **A custom `--include`/`--exclude` filter expression DSL.** Out
  of scope; the writer emits `PASS` on every record and a future
  filter slice can post-process.

Not planned at all:

- **Streaming write to a pipe via `--output -`.** The bgzf path
  needs explicit EOF block emission (the 28-byte empty-bgzf-block
  trailer); a half-closed pipe loses it silently. Plain-text
  stdout writes are not blocked by the format, but supporting
  stdout doubles the path-handling surface for one use case.
  Not in v1; revisit if asked.
- **BCF output.** The internal pipeline artefacts are `.psp`; the
  final consumer-facing format is VCF text. BCF would be a
  separate sink module if ever needed. Not on the roadmap.
- **Multi-allelic decomposition / left-alignment / normalisation.**
  Records emitted by Stage 5 are already normalised (anchor-base
  VCF convention, `alleles[0] == REF`). Re-normalising in the
  writer would duplicate Stage 5 logic and could mask Stage 5
  bugs. The writer is a pass-through encoder.

## Dependencies (Cargo.toml)

Add:

```toml
[dependencies]
noodles-vcf = "0.86"     # version pinned to match the noodles-cram = 0.92 era
noodles-bgzf = "0.41"    # for .vcf.gz output; same era as noodles-cram
```

Exact minor versions are pinned at implementation time after a
`cargo update -p noodles-vcf --dry-run` check against the existing
`noodles-cram = 0.92.0` / `noodles-core = 0.19.0` pins — the
noodles family co-releases across crates but version numbers are
per-crate, so the implementer confirms the matching pair. If a
version conflict surfaces, prefer bumping the existing pins
upward in a separate prep commit rather than downgrading the new
crates.

## Module layout

```
src/
  var_calling/
    vcf_writer/
      mod.rs              — NEW: pub use of public surface; submodule wiring
      writer.rs           — NEW: CohortVcfWriter, write_record, finish
      header.rs           — NEW: build_vcf_header(CohortMetadata) → noodles_vcf::Header
      record_encode.rs    — NEW: PosteriorRecord → noodles_vcf::Record translation
      sink.rs             — NEW: SinkKind enum dispatch (plain / bgzf) + tmp/rename
      errors.rs           — NEW: VcfWriteError
    mod.rs                — MODIFIED: pub mod vcf_writer;
```

The submodule lives under `var_calling/` so it sits next to its
upstream (`posterior_engine.rs`) rather than at the crate root
alongside the legacy `vcf_writer.rs`. The two names do not collide
(legacy: `crate::vcf_writer`, new: `crate::var_calling::vcf_writer`).

`record_encode.rs` is broken out because the
`PosteriorRecord → noodles_vcf::Record` mapping is the meat of the
plan and is worth unit-testing in isolation.

`sink.rs` is broken out because the plain-vs-bgzf dispatch has its
own well-defined contract (atomic rename, EOF-block-on-finish for
bgzf) and stays simple if it isn't tangled with VCF formatting.

## Public API

### `CohortMetadata`

```rust
/// Everything the writer needs at construction time. Composed by
/// the CLI from the upstream merger / engine handles.
pub struct CohortMetadata {
    /// Sample names in the cohort, in the order that
    /// `PosteriorRecord.posteriors` / `.scalars` rows are laid
    /// out. The writer emits a `#CHROM ... SAMPLE_0 SAMPLE_1 ...`
    /// header line in this order.
    pub sample_names: Vec<String>,
    /// Contig table — name, length, md5. Sourced from the
    /// `PerPositionMerger`'s `chromosomes()` slice (which is in
    /// turn sourced from the `.psp` headers).
    pub contigs: Vec<ParsedChromosome>,
    /// String to put in `##source=`. Conventionally
    /// `"pop_var_caller <version>"`.
    pub tool_string: String,
    /// String to put in `##commandline=`. The CLI builds this from
    /// `std::env::args()`; library users may pass an empty string.
    pub command_line: String,
}
```

### `WriterConfig`

```rust
#[derive(Debug, Clone)]
pub struct WriterConfig {
    /// Output path. Suffix selects the sink kind:
    ///   `.vcf.gz` / `.vcf.bgz` → bgzf
    ///   anything else          → plain text
    pub output: PathBuf,
    /// Always emit `PASS` in the FILTER column for now. Reserved
    /// for a future filter slice; v1 has no filter expressions.
    pub default_filter_pass: bool,
    /// When true, emit the per-sample `GP` (genotype posteriors)
    /// FORMAT field; the header gains the `##FORMAT=<ID=GP,...>`
    /// declaration and every data line carries a `GP` cell.
    /// Off by default — `GP` is `Number=G` so its per-cell size
    /// grows as `(ploidy + n_alleles - 1) choose ploidy` (21 floats
    /// per sample at ploidy=2, n_alleles=6), and most downstream
    /// consumers don't read it.
    pub emit_gp: bool,
}

impl Default for WriterConfig {
    fn default() -> Self {
        Self {
            output: PathBuf::new(), // caller must set
            default_filter_pass: true,
            emit_gp: false,
        }
    }
}
```

`WriterConfig` is intentionally narrow today; structured so future
flags (`--no-info-dp`, `--strict-monomorphic-as-ref`, etc.) slot
in without changing the public constructor signature.

### `CohortVcfWriter`

```rust
pub struct CohortVcfWriter { /* private */ }

impl CohortVcfWriter {
    /// Open the sink at `<output>.tmp`, write the header, and stand
    /// the writer up ready for `write_record` calls.
    pub fn new(
        metadata: CohortMetadata,
        config: WriterConfig,
    ) -> Result<Self, VcfWriteError>;

    /// Encode and write one record. Enforces non-decreasing
    /// `(chrom_id, pos)` order against the prior record.
    pub fn write_record(
        &mut self,
        record: &PosteriorRecord,
    ) -> Result<(), VcfWriteError>;

    /// Flush, emit the bgzf EOF block (when applicable), close the
    /// file, sync, atomic-rename `<output>.tmp` → `<output>`.
    /// Consumes self so a forgotten `finish()` shows as a missing
    /// output file rather than a silently truncated one.
    pub fn finish(self) -> Result<(), VcfWriteError>;
}
```

Forgetting `finish()` deliberately leaves the `<output>.tmp` file
on disk and no `<output>`: this is loud failure by design. The
type does not implement `Drop` with a side-effect.

### Error type

```rust
#[derive(thiserror::Error, Debug)]
pub enum VcfWriteError {
    #[error("io: {0}")] Io(#[from] io::Error),
    #[error("vcf encode: {0}")] Encode(#[from] noodles_vcf::Error),
    #[error("invalid metadata: {0}")] InvalidMetadata(String),
    #[error("record out of order: \
             record at {chrom_id}:{pos} <= previous {prev_chrom_id}:{prev_pos}")]
    RecordOutOfOrder {
        chrom_id: u32,
        pos: u32,
        prev_chrom_id: u32,
        prev_pos: u32,
    },
    #[error("record at {chrom_id}:{pos}: sample {sample_idx} \
             best_genotype index {got} out of bounds for n_genotypes {n_genotypes}")]
    GenotypeIndexOutOfBounds {
        chrom_id: u32,
        pos: u32,
        sample_idx: usize,
        got: usize,
        n_genotypes: usize,
    },
    #[error("contig '{name}' length {length} exceeds VCF i32::MAX")]
    ContigLengthOverflow { name: String, length: u32 },
}
```

`noodles_vcf::Error` is wrapped transparently (the noodles encoder
surfaces TOML/string-validation errors as its own type). `InvalidMetadata`
covers the few constructor-time checks the writer enforces (e.g. empty
`sample_names`).

## Algorithms / step-by-step

### Sink dispatch (`sink.rs`)

```rust
enum SinkKind {
    Plain(BufWriter<File>),
    Bgzf(noodles_bgzf::Writer<File>),
}
```

- Suffix selection in `SinkKind::new(path)`:
  `.vcf.gz` or `.vcf.bgz` → `Bgzf`, anything else → `Plain`. Any
  other suffix is accepted as plain text (the writer does not
  refuse `.txt` etc.; the CLI may enforce a stricter rule).
- All bytes routed through `Write`. The bgzf writer flushes blocks
  on its own schedule; the plain writer uses a 64 KiB `BufWriter`.
- `SinkKind::finish(self) -> io::Result<File>`:
  - `Plain`: `BufWriter::into_inner()?` returns the underlying
    `File`.
  - `Bgzf`: `noodles_bgzf::Writer::finish()` writes the empty-block
    EOF marker (`htslib`-required), then returns the inner `File`.
- `CohortVcfWriter::finish` then calls `File::sync_all` and
  `fs::rename(<output>.tmp, <output>)` — the standard atomic-write
  closing dance, same as the Stage 1 CLI does for `.psp`.

### Header construction (`header.rs`)

`build_vcf_header(metadata: &CohortMetadata) -> Result<noodles_vcf::Header, VcfWriteError>`

Sections in fixed order:

1. `##fileformat=VCFv4.4`.
2. `##source=<tool_string>` (e.g.
   `pop_var_caller 0.x.y`).
3. `##commandline=<command_line>` when non-empty.
4. One `##contig=<ID=<name>,length=<len>,md5=<hex>>` per entry in
   `metadata.contigs`. Length is `u32` upstream; the spec's
   `##contig=length` is a signed 32-bit integer per VCF 4.4 — we
   error on overflow (`ContigLengthOverflow`) rather than silently
   truncate, although in practice this is unreachable for any
   real chromosome length.
5. INFO definitions (numbers / types pinned to VCF 4.4 spec):
   - `AF` `Number=A Type=Float` "Allele frequency from EM"
   - `AC` `Number=A Type=Integer` "Allele count in called genotypes"
   - `AN` `Number=1 Type=Integer` "Total number of called alleles"
   - `DP` `Number=1 Type=Integer` "Total depth (sum across samples)"
   - `CA` `Number=0 Type=Flag` "At least one sample is a
     chain-anchor-broken compound call (Stage 5 fallback)"
6. FORMAT definitions (always emitted):
   - `GT` `Number=1 Type=String` "Genotype"
   - `GQ` `Number=1 Type=Integer` "Genotype quality (Phred)"
   - `DP` `Number=1 Type=Integer` "Read depth"
   - `AD` `Number=R Type=Integer` "Allelic depths for REF + ALT"

   Conditionally emitted (only when `config.emit_gp`):
   - `GP` `Number=G Type=Float` "Genotype posterior probabilities"

   Omitting the `##FORMAT=<ID=GP,...>` header when no data line
   carries `GP` keeps strict parsers (bcftools `--no-update`,
   pysam) from warning about an unused declaration.
7. `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t<SAMPLE_0>\t<SAMPLE_1>\t...`

`build_vcf_header` validates: `metadata.sample_names` is non-empty;
no duplicate sample names; no duplicate contig names.

### Record encoding (`record_encode.rs`)

```rust
fn encode(
    record: &PosteriorRecord,
    contigs: &[ParsedChromosome],
    config: &WriterConfig,
) -> Result<noodles_vcf::Record, VcfWriteError>;
```

Field-by-field mapping from `PosteriorRecord`:

| VCF column | Source on `PosteriorRecord` |
|---|---|
| `CHROM`  | `contigs[record.locus.chrom_id as usize].name` |
| `POS`    | `record.locus.pos` (already 1-based) |
| `ID`     | `.` (no IDs assigned in v1) |
| `REF`    | `record.alleles[0].seq` (REF invariant per `MergedAllele`) |
| `ALT`    | `record.alleles[1..]` joined with `,` |
| `QUAL`   | `record.qual_phred` formatted as integer-or-decimal; `INFINITY` → emit the writer's max (`9999`) as the `.psp` writer does — see Risks |
| `FILTER` | `PASS` when `config.default_filter_pass`, else `.` |

INFO assembly (one buffer per record, no per-allele allocation):

- `AF` = `record.allele_frequencies[1..]` (ALT only — REF is
  implied) formatted as `,`-separated floats with `format!("{:.6}", x)`
  precision (6 dp matches what bcftools writes).
- `AC` per ALT = count of `record.best_genotype[s]` slots that
  resolve to that ALT, summed across samples. Resolution: each
  `best_genotype[s]` is an index into `genotype_order(ploidy,
  n_alleles)` — the writer looks up the per-allele counts from
  the cached genotype table.
- `AN` = sum of (ploidy) across samples whose `best_genotype` is
  not the "no-call" sentinel (v1: every sample is called).
- `DP` = sum over samples of (sum over alleles of
  `record.scalars_row(s)[a].num_obs`).
- `CA` flag = present iff any element of
  `record.chain_anchor_flags` is `true`.

FORMAT-per-sample assembly, in `GT:GQ:DP:AD` order (with `:GP`
appended when `config.emit_gp`):

- `GT`: decode `record.best_genotype[s]` against
  `genotype_order(record.ploidy, record.alleles.len())` to get a
  multiset of allele indices; format as
  `<a>/<b>...` (unphased — Stage 5 does not currently emit phase).
  Per spec, sort the multiset ascending (`0/1` not `1/0`).
- `GQ`: `record.gq_phred[s].round().clamp(0, 99) as u8`. (Engine
  already clamps to `max_gq_phred`; this is defensive.)
- `DP`: per-sample sum of `record.scalars_row(s)[a].num_obs`.
- `AD`: `record.scalars_row(s)[a].num_obs` for each `a in
  0..alleles.len()`, `,`-joined.
- `GP` (only when `config.emit_gp`): `record.posteriors_row(s)`
  floats, `,`-joined with the same `{:.6}` precision as `AF`.

The FORMAT-key string (`GT:GQ:DP:AD` vs `GT:GQ:DP:AD:GP`) is
constant per writer instance — built once in
`CohortVcfWriter::new` from `config.emit_gp` and reused on every
record. The per-sample cells are assembled by appending colons
between the four (or five) values, so the GP path is one extra
`join` call per sample, never any per-record branching cost.

The genotype table — `genotype_order(ploidy, n_alleles)` — is
already constructible from
[`var_calling::per_group_merger::genotype_order`](../../src/var_calling/per_group_merger.rs#L355).
The writer caches one per `(ploidy, n_alleles)` pair it has seen,
keyed by `n_alleles` since `ploidy` is fixed per record. Avoids
rebuilding the table once per record.

### Order check

`CohortVcfWriter` keeps `last_locus: Option<RecordLocus>`. Each
`write_record` compares against it and errors with
`RecordOutOfOrder` on regression. Cost: one tuple compare per
record.

## Test strategy

### Unit tests (in `record_encode.rs`)

- Encode a 2-sample biallelic SNP (`A` → `T`) with
  `best_genotype = [0/0, 0/1]`, `posteriors` rows summing to 1,
  scalars filled in. Assert the rendered FORMAT cells byte-for-byte.
  Default config (`emit_gp = false`): FORMAT key string is
  `GT:GQ:DP:AD` and the per-sample cells carry four colon-separated
  values, no trailing `GP`.
- Same fixture with `emit_gp = true`: FORMAT key string becomes
  `GT:GQ:DP:AD:GP` and each per-sample cell gains a trailing
  `GP` field with the right `Number=G` arity.
- Encode a 1-sample triallelic SNP (`A`, `T`, `C`) and check that
  `AF`, `AC`, `AN` produce `Number=A`-shaped lists of length 2.
- Encode a record with `chain_anchor_flags[s][a] = true` for some
  `(s, a)` and assert the INFO column carries `CA` as a flag.
- Encode a record with `qual_phred = f64::INFINITY` and assert it
  renders as the writer's max-QUAL cap (`9999`, per Risks).
- `genotype_order` decoding round-trip: `best_genotype` values from
  `0..n_genotypes` decode through the cached table to a
  monotonically-ascending allele-index multiset that the formatter
  joins with `/`.

### Unit tests (in `header.rs`)

- Header round-trips through `noodles_vcf::Header::to_string()` ↔
  `parse` for a 3-sample, 2-contig metadata input with the
  default `emit_gp = false`. Assert the resulting `Header` has no
  `GP` FORMAT entry.
- Same input with `emit_gp = true`. Assert the resulting `Header`
  carries a `GP` FORMAT entry with `Number=G, Type=Float`.
- Empty `sample_names` → `InvalidMetadata`.
- Duplicate sample name → `InvalidMetadata`.
- Duplicate contig name → `InvalidMetadata`.

### Unit tests (in `sink.rs`)

- `.vcf` path → `SinkKind::Plain`.
- `.vcf.gz` path → `SinkKind::Bgzf`.
- `.vcf.bgz` path → `SinkKind::Bgzf`.
- Bgzf sink: after `finish()` the on-disk file ends with the 28-byte
  empty-bgzf EOF marker (assert the byte pattern directly — htslib
  refuses files missing it).
- Plain sink: `finish()` leaves no trailing nulls and the file
  size matches the bytes written.
- Atomic rename: after `finish()`, `<output>.tmp` is gone and
  `<output>` exists with the expected bytes.

### Integration test (`tests/cohort_vcf_writer_integration.rs`)

Build a hand-crafted `Vec<PosteriorRecord>` (3 records: a
biallelic SNP, a triallelic SNP, an indel) and drive the sinks
through `CohortVcfWriter`:

1. **Plain text path, default config.** Write to
   `tests/tmp/out.vcf` with `emit_gp = false`. Re-open with
   `noodles_vcf::Reader`, assert: header parses; record count = 3;
   per-record `chrom`, `pos`, `ref`, `alt` round-trip; per-sample
   `GT`, `GQ`, `DP`, `AD` round-trip; the FORMAT key string is
   `GT:GQ:DP:AD` (no `GP`); the header has no `##FORMAT=<ID=GP,…>`
   line.
2. **Plain text path, `emit_gp = true`.** Same fixture, second
   config. Assert the FORMAT key string is `GT:GQ:DP:AD:GP`; each
   per-sample row carries a `GP` cell whose float-list length
   matches `genotype_order(ploidy, n_alleles).len()`; per-sample
   `GP` rows sum to 1 within tolerance after round-trip parsing.
3. **Bgzf path.** Same fixture as (1), write to
   `tests/tmp/out.vcf.gz`. Open with `noodles_bgzf::Reader` +
   `noodles_vcf::Reader` and repeat (1)'s assertions. Also assert
   the on-disk file's tail matches the htslib EOF block
   byte-pattern.
4. **Out-of-order error.** Write record at `chr1:200`, then
   record at `chr1:100`. Assert `RecordOutOfOrder` with the right
   chrom/pos fields.

### Manual smoke

Before declaring the slice done, run a downstream smoke through
the public tools:

```
$ bcftools view tests/tmp/out.vcf.gz | head
$ bcftools stats tests/tmp/out.vcf.gz
```

`bcftools view` should print the header cleanly and the data
lines. `bcftools stats` should report a record count matching the
fixture. Document the commands in the implementation report.

## Risks and trade-offs

- **`QUAL = INFINITY` from the engine.** Per
  [posterior_engine.rs:434](../../src/var_calling/posterior_engine.rs#L434),
  `qual_phred` is `f64::INFINITY` when every sample is certainly
  variant. VCF has no sentinel for "infinite quality"; bcftools
  treats huge values as numbers but some downstream tools cap or
  reject. **Decision: cap at `9999`** (same cap freebayes and
  GATK use). The cap is hardcoded in the encoder for v1; a CLI
  knob is not justified. The implementation report should note
  the cap so a future user grepping for "9999" in a VCF can find
  the explanation.
- **`Number=A`-style fields and the multi-allelic case.** `AF`,
  `AC` per VCF 4.4 are `Number=A` (one entry per ALT). With
  `n_alleles = 1` (REF-only, no ALTs) the writer must emit `.` for
  these fields, not an empty list. The encoder special-cases
  `record.alleles.len() == 1`: such records would mean Stage 5
  emitted a non-variant site. Today the per-group merger drops
  REF-only groups upstream (see
  [per_group_merger.rs:432](../../src/var_calling/per_group_merger.rs#L432)),
  so this is defensive — if it ever fires, we emit a record with
  empty `ALT` (`.`) and `.` for the `Number=A` info fields. A test
  for this branch is included in the unit tests.
- **Sample-name order coupling.** `PosteriorRecord.posteriors`
  rows are indexed by `sample_idx` matching the upstream
  `PerPositionMerger`'s sample order. The writer trusts the
  `CohortMetadata.sample_names` it was constructed with to be in
  that same order. The CLI plan must source `sample_names` from
  `merger.sample_names()` rather than the user's CLI argument
  order, even when they look identical — single source of truth.
  Worth a doc-comment on `CohortMetadata` and a CLI-side test in
  the cohort-CLI slice.
- **`noodles-vcf` version drift.** The existing `noodles-cram =
  0.92` / `noodles-sam = 0.84` / `noodles-fasta = 0.60` pins are
  internally consistent (same release wave). Pulling
  `noodles-vcf` requires picking a matching pair; the noodles
  family co-releases but the per-crate version numbers diverge.
  Implementer runs `cargo update --dry-run` before pinning. A
  version mismatch fails to compile cleanly (noodles' shared
  `noodles-core` is the bottleneck) — there's no silent runtime
  surprise.
- **Atomic-rename predictable-tmp-name.** Same convention used by
  the Stage 1 `.psp` writer
  ([pop_var_caller_pileup_cli.md "Risks"](pop_var_caller_pileup_cli.md))
  — predictable `<output>.tmp` is left behind on crash, safe to
  delete, overwritten by re-run. Worth mentioning once in the
  `--help` text the CLI slice writes.
- **Mid-write `BrokenPipe` on plain stdout.** Not in scope (no
  stdout sink), but worth flagging: if v2 adds stdout, the
  writer needs the same `broken_pipe` graceful-shutdown the
  legacy `src/vcf_writer.rs` has. The bgzf path additionally
  has to suppress its EOF-block on `BrokenPipe`.
- **`GP` size on large allele sets — opt-in.** `Number=G` is
  `(ploidy + n_alleles - 1) choose ploidy` — for ploidy=2 and
  n_alleles=6 (the per-group-merger's max), that's 21 floats per
  sample per record. At 1000 samples that's a meaningful per-line
  size for most consumers who never read `GP`. v1 defaults
  `WriterConfig::emit_gp = false`; the cohort CLI exposes a
  `--emit-gp` flag that users opt into when they actually want
  the posterior table on disk. Toggle is honoured both in the
  header (no `##FORMAT=<ID=GP,…>` declaration when off) and in
  the per-sample cells (no trailing `:<gp_list>`), so the on-disk
  file stays consistent and bcftools doesn't warn about a
  declared-but-unused FORMAT.

## Validation

The slice is done when:

1. The unit tests in `record_encode.rs`, `header.rs`, and
   `sink.rs` pass.
2. The integration test (`cargo test --test
   cohort_vcf_writer_integration`) passes.
3. The full test suite passes (`cargo test` inside the container).
4. The manual smoke (`bcftools view` and `bcftools stats`)
   reports a clean header and a record count matching the
   fixture.
5. `cargo clippy --workspace --all-targets -- -D warnings` is clean.
6. `cargo fmt -- --check` is clean.
7. An implementation report is added under
   [doc/devel/reports/implementations/](../reports/implementations/)
   with the standard date-stamped filename, summarising what
   landed, what was deferred (tabix index, `PL`, contamination
   in INFO, stdout support), the on-disk size delta
   `emit_gp = true` vs `false` on the integration fixture
   (a single concrete number is enough; this is the only place
   the opt-in's payoff gets quantified before the CLI slice),
   and the bcftools commands that passed in the smoke.

## File touch list

New:

- `src/var_calling/vcf_writer/mod.rs` — module wiring; re-exports
  `CohortVcfWriter`, `CohortMetadata`, `WriterConfig`,
  `VcfWriteError` (~20 lines).
- `src/var_calling/vcf_writer/writer.rs` — `CohortVcfWriter`,
  `new` / `write_record` / `finish`, order check, genotype-table
  cache (~250 lines).
- `src/var_calling/vcf_writer/header.rs` — `build_vcf_header`,
  metadata validation, unit tests (~200 lines).
- `src/var_calling/vcf_writer/record_encode.rs` —
  `PosteriorRecord → noodles_vcf::Record` mapping, unit tests
  (~350 lines — the FORMAT/INFO cell formatting is the largest
  surface).
- `src/var_calling/vcf_writer/sink.rs` — `SinkKind`, suffix
  dispatch, atomic rename, unit tests (~150 lines).
- `src/var_calling/vcf_writer/errors.rs` — `VcfWriteError`
  (~40 lines).
- `tests/cohort_vcf_writer_integration.rs` — three-case
  integration test (~250 lines).
- `doc/devel/reports/implementations/cohort_vcf_writer_<date>.md`
  — implementation report once the slice lands.

Modified:

- `Cargo.toml` — add `noodles-vcf` and `noodles-bgzf`
  dependencies with pins matching the existing noodles family.
- `src/var_calling/mod.rs` — add `pub mod vcf_writer;`.

Unmodified (verified — the existing modules expose everything the
writer needs):

- `src/var_calling/posterior_engine.rs` — `PosteriorRecord`,
  `RecordLocus`, `posteriors_row`, `scalars_row`,
  `chain_anchor_flags_row` cover every field the writer reads.
- `src/var_calling/per_group_merger.rs` — `genotype_order` is the
  one public symbol the writer pulls from outside its own module.
- `src/per_sample_pileup/psp/header.rs` — `ParsedChromosome` is
  re-exported through `var_calling::per_position_merger` and
  used in `CohortMetadata` verbatim.
- Legacy `src/vcf_writer.rs` — untouched. Will be deleted later
  when the legacy gVCF-merger code path is retired (separate
  cleanup pass).

## Open questions for future slices

These are recorded so they aren't forgotten, but they do not block
this slice.

- Tabix `.tbi` index: when (and if) cohort VCFs grow large enough
  that random-access matters, add a `noodles-tabix`-based index
  builder fed by virtual offsets from the bgzf writer.
- `PL` FORMAT field: add `MergedRecord.log_likelihoods` forwarding
  through the posterior engine into `PosteriorRecord`, then a
  separate writer slice surfaces it. Decide whether `PL` is
  default-on (matches GATK) or opt-in (smaller files).
- Per-sample contamination in INFO: a `CONTAM` INFO field
  carrying the `c_s` estimate per sample is straightforward once
  the cohort CLI threads `ContaminationEstimates` into the
  writer's record path. Out of scope here; raise when the CLI
  slice lands.
- VCF 4.4 vs 4.3 vs 4.2: v1 emits 4.4. If a major consumer
  (downstream lab pipeline) refuses 4.4 we can downgrade the
  `##fileformat` line with no other change — the FORMAT/INFO
  shapes are 4.2-compatible. Decide on real-user feedback.
