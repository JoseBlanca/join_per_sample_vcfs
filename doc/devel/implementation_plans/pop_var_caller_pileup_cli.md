# `pop_var_caller pileup` — Stage 1 CLI

Implementation plan for the first user-visible binary of the multi-sample
SNP caller: a `pop_var_caller pileup` subcommand that consumes one
sample's CRAMs and a reference FASTA, runs the Stage 1 pipeline
(CRAM input → BAQ → pileup walker → `.psp` writer), and emits a
single `.psp` artefact.

Spec: [per_sample_caller.md](../specs/per_sample_caller.md) §"CLI",
§"Errors", §"Parallelism".
Architecture context:
[calling_pipeline_architecture.md](../specs/calling_pipeline_architecture.md)
§"Stage 1 — per-sample caller".

This is glue, not new algorithms. Every stage already exists:

- CRAM input → `MappedRead` stream:
  [src/per_sample_caller/cram_input.rs:681](../../src/per_sample_caller/cram_input.rs#L681)
  (`CramMergedReader::new`).
- BAQ stage (`MappedRead` → `PreparedRead`, rayon-parallel per chunk):
  [src/per_sample_caller/baq/stream.rs:111](../../src/per_sample_caller/baq/stream.rs#L111)
  (`BaqStream::new`).
- Pileup walker (`PreparedRead` → `PileupRecord`):
  [src/per_sample_caller/pileup/walker.rs:38](../../src/per_sample_caller/pileup/walker.rs#L38)
  (`pileup::run`).
- `.psp` writer:
  [src/per_sample_caller/psp/writer.rs:226](../../src/per_sample_caller/psp/writer.rs#L226)
  (`PspWriter::new` + `write_record` + `finish`).
- Walker → writer seam:
  [src/per_sample_caller/pileup_to_psp.rs:50](../../src/per_sample_caller/pileup_to_psp.rs#L50)
  (`drive_pileup_to_psp`).
- Memory-bounded `RefSeqFetcher`:
  [src/per_sample_caller/ref_fetcher.rs:44](../../src/per_sample_caller/ref_fetcher.rs#L44)
  (`ChromBoundaryRefFetcher::new`).

What this plan adds is the orchestration: a new binary entry point, a
clap-based subcommand layer, the `MappedRead → PreparedRead` error
bridge between `BaqStream` (`Result<…, CramInputError>`) and the
walker (`Item = PreparedRead`), a `WriterHeader` factory, a top-level
`run_pileup` driver function, and an end-to-end integration test.

## Scope

In scope:

- A new `[[bin]]` target `pop_var_caller` with a clap-derive
  subcommand surface: `pop_var_caller pileup [opts] <crams…>`.
- A `pileup` orchestration module under
  [src/per_sample_caller/](../../src/per_sample_caller/) that
  constructs the four stages, threads them together, and writes one
  `.psp` file.
- **Full configurable surface on the CLI**: every tunable field of
  `CramMergedReaderConfig`, `BaqConfig`, `BaqStream` chunk size, and
  `WalkerConfig` is exposed as a `--flag` with the module's
  default. Spec-pinned values (the PSP writer's 16 MiB block
  target) are **not** exposed — see "Not planned at all".
- `--no-baq` to skip the BAQ HMM entirely. The pipeline still
  produces `PreparedRead`s — `bq_baq` is just a copy of the raw
  `qual` from the CRAM record. Needs one small helper in
  [baq/engine.rs](../../src/per_sample_caller/baq/engine.rs) (see
  "Wiring change for `--no-baq`" below).
- Hard-error enforcement that every CRAM's `@SQ M5` is present
  (the `WriterHeader` requires a 32-hex-char md5 per chromosome and
  the spec demands it).
- A run-summary stderr log: `FilterCounts` (from the merged
  reader) + `BaqSkipCounts` (from the BAQ stream) +
  `RunSummary` (from the walker), printed once at end-of-run.
- Threading: a single `--threads` flag sized into the global rayon
  pool (BAQ already uses rayon via `par_drain`); no extra
  pipelining work this slice.
- An integration test that drives a tiny CRAM (built from a fixture
  FASTA + synthetic SAM at test setup time, using the existing
  [cram_files.rs](../../src/per_sample_caller/cram_files.rs)
  helpers) through the full pipeline and reads the resulting `.psp`
  back with `PspReader` for record-level parity. A second test
  case exercises the `--no-baq` path.

Deferred to a follow-up slice (intended to land later, not in this
cut):

- `--region`. Adds a CRAI dependency, region parsing, and a second
  noodles read path (`query` vs `records`). Easy to slot in later;
  no API of the modules below depends on this slice being
  region-aware.
- Computing missing per-contig MD5 from FASTA. Today: hard-error
  with a clear message. Easy to relax later.

Not planned at all (deliberately not part of this product
direction — the items below are listed only to head off "you forgot
to mention…" questions on review):

- **Cross-sample parallelism inside one process.** Pileups for N
  samples are produced by running N independent
  `pop_var_caller pileup` processes in parallel. No in-process
  multi-sample driver, ever.
- **A subcommand wrapping the legacy gVCF-merger code paths.** The
  legacy `merge_alleles_and_genotypes` pipeline is not exposed
  from the new binary. It stays as library code only until it is
  deleted in a separate cleanup pass.
- **A CLI flag for the PSP writer's block-size target.** The
  16 MiB target is spec-pinned (the doc-comment on
  [`new_with_block_target`](../../src/per_sample_caller/psp/writer.rs#L235)
  says explicitly: not for production use; the spec pins it). A
  CLI flag would invite users to violate that contract.
- **Writing the `.psp` to stdout.** The `.psp` format is
  semi-random-access (the trailer points back into the body, the
  index sits at the tail); a pipe is not a useful sink because no
  downstream consumer can seek into it. Output is always a regular
  file. **Deviates from the spec sketch in
  [per_sample_caller.md](../specs/per_sample_caller.md) §CLI**
  ("`--output -` writes to stdout") — the spec sketch should be
  updated to match in the same cleanup pass, but the deviation is
  not a blocker for this slice.

## Dependencies (Cargo.toml)

Add:

```toml
[dependencies]
clap = { version = "4", features = ["derive"] }

[[bin]]
name = "pop_var_caller"
path = "src/main.rs"
```

The `[package].name` stays `merge_per_sample_vcfs` so the existing
library re-imports under that crate name continue to work. The
binary's published name is `pop_var_caller`, independent of the
crate name. Rename the crate later as a separate cleanup.

The existing `merge_alleles_and_genotypes` pipeline stays compiled
as library code (it's `pub` from `src/lib.rs`) but the new
`src/main.rs` no longer invokes it. Removing the legacy code is a
separate deletion pass, not in this plan.

## Module layout

```
src/
  main.rs                              — replaced: clap dispatch, no orchestration logic
  lib.rs                               — add `pub mod pop_var_caller;`
  pop_var_caller/
    mod.rs                             — NEW: pub mod cli; re-export run_pileup
    cli.rs                             — NEW: PileupArgs, PileupCliError, WriterHeader
                                                factory, run_pileup, stderr summary
    cli/
      error_bridge.rs                  — NEW: BaqStream → walker error-shedding adapter
  per_sample_caller/                   — unchanged
    cram_input.rs, baq/, pileup/, psp/, pileup_to_psp.rs, ref_fetcher.rs
```

Separation of concerns: `per_sample_caller` is the algorithm-side
library slice (header validation, BAQ, walker, `.psp` writer/reader).
The new `pop_var_caller` module is the binary's organising namespace —
CLI struct, subcommand dispatch, orchestration. Future subcommands
(e.g. a multi-sample call) slot in as siblings of `cli.rs` here
without touching the per-sample pipeline modules.

The `error_bridge` is broken out only because the adapter has a
non-trivial invariant worth unit-testing on its own; it could fold
back into `cli.rs` if it ever shrinks.

## Public API

### CLI surface (clap-derive)

```rust
#[derive(clap::Parser)]
#[command(name = "pop_var_caller", version)]
struct Cli {
    #[command(subcommand)]
    cmd: Subcommand,
}

#[derive(clap::Subcommand)]
enum Subcommand {
    /// Stage 1: run BAQ + pileup over one sample's CRAMs, emit a .psp.
    Pileup(PileupArgs),
}

#[derive(clap::Args)]
struct PileupArgs {
    // ===== COMMON FLAGS (visible in `-h`) =====================

    // ----- I/O ------------------------------------------------
    /// Reference FASTA. A sibling `.fai` is required.
    #[arg(long)]
    reference: PathBuf,

    /// Output .psp path. Must be a regular file — stdout is not
    /// supported (the `.psp` format needs random access on read,
    /// so a pipe is not a useful sink).
    #[arg(long)]
    output: PathBuf,

    /// Worker threads for the BAQ stage and any other rayon work.
    /// Default: all logical cores.
    #[arg(long)]
    threads: Option<usize>,

    /// One or more coordinate-sorted CRAMs for the sample.
    #[arg(required = true)]
    crams: Vec<PathBuf>,

    // ----- common filters / toggles ---------------------------
    /// Drop reads with MAPQ < N. `0` admits everything.
    /// Default: 20.
    #[arg(long, default_value_t = DEFAULT_MIN_MAPQ)]
    min_mapq: u8,

    /// Skip the BAQ HMM. `bq_baq` becomes a copy of the raw CRAM
    /// `QUAL`. The walker, writer, and run-summary path are
    /// unchanged; only BAQ is bypassed.
    #[arg(long)]
    no_baq: bool,

    // ===== ADVANCED FLAGS (hidden from `-h`, shown in `--help`) =====
    // Marked with `hide_short_help = true`. `--help` shows them
    // organised under `next_help_heading` sections so the long
    // help isn't a flat 13-flag list.

    // ----- advanced CRAM input filters -----------------------
    #[command(next_help_heading = "Advanced — CRAM input filters")]

    /// Drop reads with decoded SEQ length < N. `0` admits any length.
    /// Default: 30.
    #[arg(long, hide_short_help = true, default_value_t = DEFAULT_MIN_READ_LENGTH)]
    min_read_length: u32,

    /// Keep reads with the QC-fail flag (0x200) set.
    /// Default: drop them.
    #[arg(long, hide_short_help = true)]
    keep_qc_fail: bool,

    /// Keep reads with the duplicate flag (0x400) set.
    /// Default: drop them.
    #[arg(long, hide_short_help = true)]
    keep_duplicates: bool,

    /// Drop reads whose M-op mismatch fraction exceeds X. Must be
    /// in [0.0, 1.0]. Pass `0.0` to disable the filter entirely.
    /// Anything negative, > 1.0, NaN, or ±∞ is a hard error.
    /// Default: 0.10.
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MAX_READ_MISMATCH_FRACTION,
          value_parser = parse_mismatch_fraction)]
    max_read_mismatch_fraction: f32,

    /// BQ floor below which a mismatch does not count toward
    /// `--max-read-mismatch-fraction`. `0` makes every mismatch
    /// count. Default: 10.
    #[arg(long, hide_short_help = true, default_value_t = DEFAULT_MISMATCH_BQ_FLOOR)]
    mismatch_bq_floor: u8,

    // ----- advanced BAQ HMM tuning ---------------------------
    #[command(next_help_heading = "Advanced — BAQ HMM")]

    /// BAQ gap-open probability. samtools/htslib default: 1e-3.
    #[arg(long, hide_short_help = true,
          default_value_t = SAMTOOLS_ILLUMINA_GAP_OPEN_PROB)]
    baq_gap_open_prob: f32,

    /// BAQ gap-extension probability. samtools/htslib default: 0.1.
    #[arg(long, hide_short_help = true,
          default_value_t = SAMTOOLS_ILLUMINA_GAP_EXTEND_PROB)]
    baq_gap_extend_prob: f32,

    /// BAQ band half-width. samtools/htslib default: 7. The engine
    /// auto-widens per-read on long indels — this is the floor.
    #[arg(long, hide_short_help = true,
          default_value_t = SAMTOOLS_ILLUMINA_BAND_HALF_WIDTH)]
    baq_band_half_width: i32,

    /// BAQ batch size: reads processed in parallel per rayon
    /// chunk. Larger = better throughput, more peak memory.
    /// Default: 1024.
    #[arg(long, hide_short_help = true, default_value_t = DEFAULT_BAQ_CHUNK_SIZE)]
    baq_chunk_size: usize,

    // ----- advanced walker caps ------------------------------
    #[command(next_help_heading = "Advanced — Pileup walker")]

    /// Max contributors folded at a pure-SNP/REF column.
    /// Default: 8000.
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MAX_SNP_COLUMN_DEPTH)]
    max_snp_column_depth: u32,

    /// Max contributors folded at a column carrying any indel
    /// observation. Default: 250.
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MAX_INDEL_COLUMN_DEPTH)]
    max_indel_column_depth: u32,

    /// Hard cap on per-record reference span. Default: 5000.
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MAX_RECORD_SPAN)]
    max_record_span: u32,

    /// How far past a first mate the walker keeps a pending-mates
    /// entry before evicting and treating the first mate as solo.
    /// Default: 10000.
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MATE_LOOKUP_WINDOW)]
    mate_lookup_window: u32,

    /// Hard cap on concurrently-active reads (defensive bound;
    /// exceeding it is a hard error). Default: 4096.
    #[arg(long, hide_short_help = true,
          default_value_t = DEFAULT_MAX_ACTIVE_SLOTS)]
    max_active_reads: u32,
}
```

**Tiered help.** Common flags (`--reference`, `--output`,
`--threads`, `--min-mapq`, `--no-baq`, plus the positional CRAMs)
show up in the short `-h` output. The remaining 12 flags are
tagged `hide_short_help = true` so they only appear under
`--help`, organised by `next_help_heading` into three
sub-sections: "Advanced — CRAM input filters", "Advanced — BAQ
HMM", and "Advanced — Pileup walker". This mirrors `cargo`, `rg`,
and `git`'s convention — typing `-h` gets you a tight overview,
typing `--help` gets you the full power-user surface. No second
source of truth (no `--config` file), every effective value still
records itself into the `.psp`'s `WriterProvenance.parameters` on
output.

Conventions used above:

- **All defaults are pulled from the module's existing `pub const
  DEFAULT_*` symbols** so the source of truth stays in the module.
  Any future tweak to a module default is automatically picked up
  by the CLI; the `--help` text shows the live value via
  clap's `default_value_t`.
- `Option<u8>` / `Option<u32>` "disable" knobs in
  `CramMergedReaderConfig` (`min_mapq`, `min_read_length`) are
  surfaced as plain `u8` / `u32` on the CLI with `0` as the
  "disable" sentinel — the orchestrator translates back to
  `None`. This keeps clap's help text tidy and matches what users
  expect from samtools/bcftools.
- `max_read_mismatch_fraction: Option<f32>` is surfaced as plain
  `f32`. The CLI accepts values in `[0.0, 1.0]`. `0.0` translates
  to `None` (filter disabled); a value in `(0.0, 1.0]` translates
  to `Some(x)`. Anything outside that range — negative, > 1.0,
  NaN, or ±∞ — is **rejected at parse time** by a small
  `parse_mismatch_fraction` `value_parser` (clap surfaces the
  error before `run_pileup` is even called). The user gets a
  clear message at the CLI boundary rather than a silent
  reinterpretation of an odd input. Sketch:

  ```rust
  fn parse_mismatch_fraction(s: &str) -> Result<f32, String> {
      let v: f32 = s.parse().map_err(|e| format!("not a number: {e}"))?;
      if !v.is_finite() {
          return Err(format!("must be finite, got `{s}`"));
      }
      if !(0.0..=1.0).contains(&v) {
          return Err(format!("must be in [0.0, 1.0], got `{s}`"));
      }
      Ok(v)
  }
  ```
- The two default-on `drop_*` booleans become opt-out `keep_*`
  flags so the absence of the flag matches the default.
- The `max_active_slots` field is renamed on the CLI to
  `--max-active-reads`, matching its current doc-comment ("under
  the unique-`u64`-chain-id design it bounds concurrent active
  reads, not a slot id namespace") — internal field stays as-is.
- `--region` and any future flags can be added without
  rearrangement; the struct stays append-only.

### `run_pileup` — the orchestrator

In `src/pop_var_caller/cli.rs`:

```rust
pub fn run_pileup(args: &PileupArgs) -> Result<(), PileupCliError>;
```

Behaviour:

1. Validate inputs cheaply (FASTA exists, `.fai` exists — the
   merged reader's own `new` already enforces this, but doing it up
   front gives a better error message before the CRAM opens).
2. Translate CLI values into module configs:
   - `CramMergedReaderConfig` — set every field from the
     corresponding CLI flag. The translations are:
     - `min_mapq`: `0` → `None`, `n` → `Some(n)`.
     - `min_read_length`: `0` → `None`, `n` → `Some(n)`.
     - `max_read_mismatch_fraction`: `0.0` → `None`, `x` →
       `Some(x)`. The CLI's `value_parser` already rejected
       out-of-range values, so no validation needed here.
     - Booleans are `drop_qc_fail = !args.keep_qc_fail`,
       similarly for duplicates.
   - `BaqConfig` — set every field from the `--baq-*` flags.
   - `WalkerConfig` — set every field from the corresponding flag.
3. Size rayon's global pool from `--threads` (early-return error if
   the pool is already configured by an earlier call; the binary
   only ever calls this once per process).
4. Open `CramMergedReader::new(&crams, &reference, cram_cfg)?`.
   Extract `sample_name()` and `contigs()` for the writer header.
5. Build `WriterHeader` via `build_writer_header(...)` (below).
   This step **hard-errors** if any contig has no `@SQ M5`.
6. Construct the production `RefSeqFetcher`:
   `ChromBoundaryRefFetcher::new(&reference, contigs.clone())?`.
7. Build the `PreparedRead` source. Branch on `args.no_baq`:
   - `false` → `BaqStream::new(reader, baq_cfg, &fetcher,
     args.baq_chunk_size)`.
   - `true` → a `map`-adapter over the merged reader that calls
     `baq::prepare_passthrough` per record. See "Wiring change
     for `--no-baq`".
8. Wrap that stream in the error-shedding adapter (next
   subsection) so its `Item` becomes `PreparedRead` and the first
   `CramInputError` is held aside for end-of-run propagation.
9. Construct `PileupWalker` via `pileup::run(prepared_reads,
   &fetcher, &walker_cfg)`.
10. Open the output sink as `BufWriter<File>` pointing at
    `<output>.tmp` (the rename to the final `<output>` happens on
    success — see Risks).
11. Construct `PspWriter::new(sink, header)?`.
12. Call `drive_pileup_to_psp(walker, writer)?`. Collect
    `(_sink, RunSummary)`.
13. Check the error bridge for a stashed upstream error — if
    present, return it (the walker would have already stopped
    pulling, so this is the right place to surface it).
14. Pull `FilterCounts` from the reader, and (BAQ-on) the
    `BaqSkipCounts` from the BAQ stream. Format the three counter
    blocks to stderr. When `--no-baq`, the BAQ block prints
    `baq: disabled`.
15. `BufWriter::into_inner()?` to recover the underlying `File`,
    then `File::sync_all()`, then atomic
    `fs::rename(<output>.tmp, <output>)`.

### Wiring change for `--no-baq`

The walker takes `Iterator<Item = PreparedRead>`. With BAQ on, that
iterator is a `BaqStream` (after error-shedding). With BAQ off, we
need a stream that produces `PreparedRead`s from `MappedRead`s
without running the HMM.

The engine module already has the private function
[`mapped_to_prepared`](../../src/per_sample_caller/baq/engine.rs#L381)
that builds a `PreparedRead` from a `MappedRead` plus a `bq_baq`
buffer. The no-BAQ path is exactly that, with `bq_baq` =
`read.qual.clone()` (or moved). One small change to the BAQ module:

```rust
// src/per_sample_caller/baq/engine.rs — promote one private fn

pub fn prepare_passthrough(read: MappedRead, chrom_id: u32) -> PreparedRead {
    let bq_baq = read.qual.clone();
    mapped_to_prepared(read, chrom_id, bq_baq)
}
```

Re-exported from `baq/mod.rs`. No new types, no new errors; the
helper just exposes the work the engine already does in
[`process`](../../src/per_sample_caller/baq/engine.rs) for the
case where the HMM step is bypassed.

The orchestrator branches on `args.no_baq`:

```rust
let prepared: Box<dyn Iterator<Item = Result<PreparedRead, CramInputError>>> =
    if args.no_baq {
        Box::new(reader.map(|r| {
            r.map(|read| {
                let chrom_id = u32::try_from(read.ref_id).expect("ref_id fits u32");
                prepare_passthrough(read, chrom_id)
            })
        }))
    } else {
        Box::new(BaqStream::new(reader, baq_cfg, &fetcher, args.baq_chunk_size))
    };
```

A `Box<dyn Iterator>` adds one indirection per read; at the
per-read budget the pipeline already pays (BAQ HMM, walker
bookkeeping, writer encoding), this is invisible. If profiling
ever flags it, an enum dispatch is a one-day refactor.

The BAQ-disabled path bypasses `BaqStream`, so its
`BaqSkipCounts` is the default zero — the stderr summary shows
`baq: disabled` instead of the counts block in that case.

### Error bridge: `BaqStream` → walker

`BaqStream` yields `Result<PreparedRead, CramInputError>`. The walker
takes `Iterator<Item = PreparedRead>` — no `Result`. The bridge
adapts shape *and* preserves error context:

```rust
// src/pop_var_caller/cli/error_bridge.rs

pub struct ErrorSheddingAdapter<I>
where
    I: Iterator<Item = Result<PreparedRead, CramInputError>>,
{
    inner: I,
    stashed: Rc<RefCell<Option<CramInputError>>>,
}

impl<I> ErrorSheddingAdapter<I> {
    pub fn new(inner: I) -> Self;
    /// Cloned handle the caller keeps so it can inspect (and take)
    /// the stashed error after the walker exhausts.
    pub fn error_handle(&self) -> ErrorHandle;
}

pub struct ErrorHandle(Rc<RefCell<Option<CramInputError>>>);

impl ErrorHandle {
    pub fn take(&self) -> Option<CramInputError>;
}

impl<I> Iterator for ErrorSheddingAdapter<I> {
    type Item = PreparedRead;
    fn next(&mut self) -> Option<PreparedRead> {
        match self.inner.next()? {
            Ok(r) => Some(r),
            Err(e) => { *self.stashed.borrow_mut() = Some(e); None }
        }
    }
}
```

The adapter is single-threaded; `Rc<RefCell<…>>` is fine. The walker
sees the upstream's first error as a clean end-of-stream and runs its
own flush logic (closing the final chromosome and emitting any
pending records). After `drive_pileup_to_psp` returns, the
orchestrator calls `ErrorHandle::take()` once — if it yields
`Some(e)`, that becomes the `PileupCliError::CramInput(e)` returned
from `run_pileup`.

Notes on this design:

- It deliberately runs the walker to natural completion on what it
  did get. The records already produced before the upstream error are
  not retained on disk because we surface the error before reporting
  success; the `.psp` file's `flush` is not called. (See "Risks /
  partial output" below.) The walker still completes cleanly so its
  internal state can produce a meaningful `RunSummary` for the
  error log.
- An alternative — propagate the error through the walker via a
  `Result` type parameter — would require touching every walker
  call site. The adapter localises the impedance mismatch to one
  place and the walker API stays unchanged.

### `WriterHeader` factory

```rust
fn build_writer_header(
    sample: &str,
    fasta_path: &Path,
    contigs: &ContigList,
    cram_paths: &[PathBuf],
    args: &PileupArgs,
) -> Result<WriterHeader, PileupCliError>;
```

Responsibilities:

- Set `format_version = (1, 0)`.
- Convert each `ContigEntry` to a `ChromosomeEntry`:
  - `name`: `String` from `ContigEntry.name`.
  - `length`: `u32::try_from(entry.length)` (error if > `u32::MAX` —
    won't happen in practice for SAM contigs).
  - `md5`: `entry.md5.ok_or(MissingMd5 { contig })?` — formatted as
    lowercase hex via `format!("{:032x}", u128::from_be_bytes(...))`
    or equivalent. **This is the hard-error point for missing `@SQ M5`.**
- Set `created = chrono::Utc::now()` formatted as a TOML `Datetime`
  (the header type uses `toml::value::Datetime`, see
  [psp/header.rs:95](../../src/per_sample_caller/psp/header.rs#L95);
  no new dep required — `toml` is already in `Cargo.toml`).
- Build `WriterProvenance`:
  - `tool = "pop_var_caller"`
  - `version = env!("CARGO_PKG_VERSION")`
  - `subcommand = "pileup"`
  - `input_crams`: basenames only (the writer header doc-comment is
    explicit about this — no directory components).
  - `input_fasta`: basename of `fasta_path`.
  - `parameters`: **every effective CLI knob**, including ones left
    at their default — the on-disk file is a self-describing
    record of what produced it, so a re-run from the same `.psp`
    can reproduce. Suggested key set (BTreeMap, deterministic
    order):
    - `min_mapq` (Integer), `min_read_length` (Integer),
      `drop_qc_fail` (Boolean), `drop_duplicate` (Boolean),
      `max_read_mismatch_fraction` (Float; recorded as the
      effective value the user supplied — `0.0` when the filter
      is disabled, the actual fraction otherwise — matching the
      CLI's "0.0 disables" convention), `mismatch_bq_floor`
      (Integer).
    - `baq_enabled` (Boolean), `baq_gap_open_prob` (Float),
      `baq_gap_extend_prob` (Float), `baq_band_half_width`
      (Integer), `baq_chunk_size` (Integer).
    - `max_snp_column_depth` (Integer),
      `max_indel_column_depth` (Integer),
      `max_record_span` (Integer), `mate_lookup_window`
      (Integer), `max_active_reads` (Integer).
    - `threads` (Integer — record the *effective* count after
      rayon was sized, not `None`).
  - `parameters` does **not** record the input/output paths or the
    sample name — those live in dedicated header fields.

The PSP header parser will then accept the file we write — its own
test fixture
[psp/test_fixtures.rs:46](../../src/per_sample_caller/psp/test_fixtures.rs#L46)
shows the shape this matches.

### Error type

```rust
#[derive(thiserror::Error, Debug)]
pub enum PileupCliError {
    #[error("CRAM input: {0}")] CramInput(#[from] CramInputError),
    #[error("PSP writer: {0}")] Psp(#[from] PspWriteError),
    #[error("pipeline: {0}")]  Pipeline(#[from] PileupToPspError),
    #[error("contig '{contig}' has no @SQ M5 in any input CRAM")]
    MissingMd5 { contig: String },
    #[error("rayon thread pool already initialised — refusing to override")]
    RayonAlreadyConfigured,
    #[error("io: {0}")] Io(#[from] io::Error),
}
```

This is the only crate-public error type added; the writer and
walker errors flow through it transparently.

## Algorithms / step-by-step

This is glue, not algorithms. The step list in
[Public API → `run_pileup`](#run_pileup--the-orchestrator) above is
authoritative. Two non-obvious bits worth flagging:

- **Order of error checks.** The orchestrator checks the error
  bridge *after* `drive_pileup_to_psp` returns `Ok(...)`. If the
  bridge has a stashed error, that error wins even though the seam
  itself succeeded. (The walker stopped naturally because the
  adapter shed an error from the upstream; from the walker's
  perspective there was no error.) This is the only place the
  error-then-clean-flush ordering matters.
- **Output sink lifetime.** `drive_pileup_to_psp` returns the
  underlying `W` (here a `BufWriter<File>`) after `finish()`. We
  then call `BufWriter::into_inner()` → `File::sync_all()` →
  `fs::rename`.

## Test strategy

### Unit tests (in `pop_var_caller/cli.rs`)

- `build_writer_header` happy path on a `ContigList` with M5 set.
- `build_writer_header` errors with `MissingMd5 { contig }` when a
  `ContigEntry.md5` is `None`.
- `build_writer_header` strips directory components from
  `input_crams` and `input_fasta`.
- `build_writer_header` formats MD5 as exactly 32 lowercase hex
  chars (`check_md5_hex` from the writer side then accepts it).

### Unit tests (in `pop_var_caller/cli/error_bridge.rs`)

- Bridge passes through `Ok` items unchanged.
- Bridge stashes the first `Err` and then returns `None` forever.
- `ErrorHandle::take` returns the stashed error exactly once;
  subsequent `take` calls return `None`.
- Bridge composes with the walker without panicking when the
  underlying CRAM stream errors mid-stream (synthetic
  `Vec<Result<PreparedRead, …>>` as the upstream).

### Integration test (`tests/pileup_cli_integration.rs`)

Three cases. All build their CRAMs in-test using
[per_sample_caller/cram_files.rs](../../src/per_sample_caller/cram_files.rs)
on a small FASTA fixture (use the existing
[tests/data/baq](../../tests/data/baq) layout or a new fixture
under `tests/data/pileup`):

1. **Happy path, default config.** Two-CRAM input (same sample SM
   tag, identical `@SQ` list with M5), one short contig, ~50
   reads. Run `run_pileup` programmatically (i.e. call the library
   entry — the binary itself is not the test surface). Re-open the
   `.psp` with `PspReader` and assert:
   - The chromosome list matches what came out of `contigs()`.
   - `RunSummary.records_emitted` matches the per-position record
     count in the `.psp`.
   - `WriterProvenance.subcommand == "pileup"`.
   - `WriterProvenance.input_crams` are basenames, not paths.
   - Every CLI-exposed knob appears as a key in
     `WriterProvenance.parameters`, with the default values.

2. **`--no-baq` path.** Same fixture as (1); set `args.no_baq =
   true`. Assert: the `.psp` still round-trips through `PspReader`;
   `WriterProvenance.parameters["baq_enabled"] == Boolean(false)`;
   per-record `bq_baq` in the emitted records matches the raw
   `QUAL` from the CRAM input (no values capped, no BAQ HMM
   smoothing applied).

3. **Missing-M5 hard-error.** Build one CRAM whose `@SQ` lacks
   `M5` (the `cram_files.rs` helper has an `md5_overrides` field;
   bypass it). Call `run_pileup`; assert
   `PileupCliError::MissingMd5` with the expected contig name.

The integration test does *not* exercise BAQ correctness — that is
already covered by
[src/per_sample_caller/baq/tests.rs](../../src/per_sample_caller/baq/tests.rs).
It is only there to lock the glue.

### Manual smoke

Before declaring the slice done, run

```
$ cargo run --release --bin pop_var_caller -- pileup \
    --reference reference.fa \
    --output /tmp/out.psp \
    sample.cram
```

inside the container against a sample CRAM and the bundled FASTA;
confirm:
- the binary runs to completion, exit 0;
- the stderr summary block prints `reads_admitted`,
  `records_emitted`, and the `FilterCounts` / `BaqSkipCounts`
  fields;
- the `.psp` opens cleanly with a follow-up `PspReader`-based
  smoke (a one-off `cargo run --example` script, no separate test).

## Risks and trade-offs

- **Partial `.psp` on upstream error — handled via write-tmp-then-
  rename.** The error bridge lets the walker drain a partially-read
  upstream cleanly. If we then bail before calling
  `BufWriter::flush` / `File::sync_all`, the partial bytes are
  typically lost to the OS buffer — but a process kill or power
  failure between the seam returning and our error check could in
  principle leave a partial file. Worse, a truncated `.psp` whose
  TOML header still parses but whose body ends mid-block would look
  "present" to Stage 2 and to pipeline tooling
  (`Path.exists()`-style checks). We handle this by writing to
  `<output>.psp.tmp` and only `fs::rename`-ing into place after
  `finish()` + `flush` + `sync_all` succeed. POSIX `rename(2)` is
  atomic within a filesystem, so the path either holds the
  previous version or a complete new version, never a half-write.

  **Convention note for future readers.** This is standard practice
  in systems software (editors, `git`, `dpkg`, snakemake/nextflow
  output publishing, Python's `Path.replace`) but a deliberate
  departure from common bioinformatics-CLI norms — samtools,
  bcftools, bwa, picard, GATK, and freebayes all write straight
  to the output path and leave half-baked files behind on failure.
  We adopt the safer convention because (a) it costs one
  `fs::rename` line, (b) `.psp` is a structured artefact whose
  partial form is silently parseable, and (c) it composes cleanly
  with workflow managers that key off output-file existence.
  One edge case worth knowing: a killed process leaves
  `<output>.psp.tmp` behind. The name is predictable, the file is
  safe to delete, and a re-run overwrites it. The `--help` text
  mentions this.
- **`@SQ M5`-required is strict.** Some old or hand-rolled CRAMs
  omit M5. The error message must spell out the fix
  ("run `samtools quickcheck` / regenerate the CRAM with an
  M5-aware tool, or wait for the
  `--compute-missing-md5` follow-up"). Listed in
  [Out of scope](#scope).
- **Thread-pool init.** Calling `rayon::ThreadPoolBuilder::build_global`
  twice in a process panics. The CLI does it once; library callers of
  `run_pileup` (the integration test) must not pre-configure rayon.
  We handle this by attempting the init and converting the error to
  `RayonAlreadyConfigured`; in the test we never set it. Worth
  documenting on `run_pileup`.
- **No `.psp` block-size flag.** The writer flushes blocks at
  16 MiB by default
  ([psp/writer.rs:236](../../src/per_sample_caller/psp/writer.rs#L236)
  carries the override-aware constructor); we do not expose it
  (the override is `#[doc(hidden)]` and the spec pins the target).
- **Parameter surface size — managed via tiered help.** The
  `pileup` subcommand exposes ~17 flags. The CLI tiers them: 5
  common flags (plus the positional CRAMs) appear in `-h`; the
  other 12 are tagged `hide_short_help = true` and grouped under
  `next_help_heading` sections for `--help`. Defaults are pulled
  from each module's `pub const DEFAULT_*` so the live help text
  always shows the current defaults, no manual sync.

  If recurring overrides become a pattern — e.g. you find yourself
  passing the same five advanced flags on every invocation — the
  signal to revisit is a `--config foo.toml` option layered on top
  of the existing CLI (CLI overrides config, config overrides
  defaults). Not in this slice; raise it then.
- **`max_read_mismatch_fraction` input validation.** Module
  config holds `Option<f32>`; the CLI flattens to `f32`. Decided
  convention:
  - `0.0` → `None` (filter disabled).
  - `x ∈ (0.0, 1.0]` → `Some(x)`.
  - Anything else (negative, > 1.0, NaN, ±∞) → hard error at
    parse time via a clap `value_parser`. A user passing such a
    value is doing something odd and deserves a clear message
    at the CLI boundary, not a silent reinterpretation.

  `WriterProvenance.parameters["max_read_mismatch_fraction"]`
  records the effective f32 the user supplied — `0.0` when
  disabled, the real fraction otherwise. The provenance reads
  back identically to whatever the user originally typed, no
  magic sentinel needed.

## Validation

The slice is done when:

1. The integration test passes (`cargo test --test pileup_cli_integration`).
2. The full test suite passes (`cargo test`).
3. The manual smoke run produces a `.psp` that round-trips through
   `PspReader` and reports a non-zero `records_emitted`.
4. `cargo clippy --workspace --all-targets -- -D warnings` is clean.
5. `cargo fmt -- --check` is clean.
6. An implementation report is added under
   [ia/reports/implementations/](../reports/implementations/) using
   the standard date-stamped filename, summarising:
   what landed, what was deferred (`--region`, FASTA-based MD5
   compute), and the smoke-test command line that worked.

## File touch list

New:

- `src/main.rs` — rewritten as the clap dispatch shell, ~30 lines.
- `src/pop_var_caller/mod.rs` — `pub mod cli; pub use
  cli::{run_pileup, PileupArgs, PileupCliError};` (~10 lines).
- `src/pop_var_caller/cli.rs` — `PileupArgs`, `PileupCliError`,
  `build_writer_header`, `run_pileup`, parameter-to-config
  translation, and the summary-formatting helpers (~400 lines —
  the parameter surface is wide).
- `src/pop_var_caller/cli/error_bridge.rs` —
  `ErrorSheddingAdapter`, `ErrorHandle`, unit tests (~120 lines).
- `tests/pileup_cli_integration.rs` — three-case integration test
  (~200 lines).
- `ia/reports/implementations/pop_var_caller_pileup_cli_<date>.md` —
  implementation report once the slice lands.

Modified:

- `Cargo.toml` — add `clap` dependency, add `[[bin]] name =
  "pop_var_caller"` section.
- `src/lib.rs` — add `pub mod pop_var_caller;`.
- `src/per_sample_caller/baq/engine.rs` /
  `src/per_sample_caller/baq/mod.rs` — promote one helper
  (`prepare_passthrough`) to `pub` for the `--no-baq` path.
  Pure code-motion: no behaviour change to BAQ-on.

Unmodified (verified — the existing modules expose everything the
glue needs):

- `src/per_sample_caller/{cram_input.rs, pileup/*, psp/*,
  pileup_to_psp.rs, ref_fetcher.rs}`. The walker, writer, and
  seam are touched only through their existing public API.
- `src/per_sample_caller/mod.rs` — the per-sample slice does not
  gain a `cli` submodule; the CLI lives under `pop_var_caller/`.

## Open questions for future slices

These are recorded so they aren't forgotten, but they do not block
this slice.

- Does `--region` need to be ergonomic when the sample's `.psp` is
  going to be regenerated from scratch anyway? Stage-2 callers
  might prefer "whole-CRAM only" as the per-sample contract.
- A `pop_var_caller join` (or similar) subcommand for Stage 3+
  will share the same clap surface. Keep the subcommand enum
  open-ended.
