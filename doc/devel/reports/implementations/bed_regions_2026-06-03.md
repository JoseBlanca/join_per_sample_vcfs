# Analysis regions from a BED file — implementation report (2026-06-03)

Implementation report for
[bed_regions.md](../../implementation_plans/bed_regions.md): add
`--regions <bed>` to both `pileup` and `var-calling`, restricting
analysis to the listed regions and *seeking* to them via the on-disk
indexes instead of streaming the whole genome.

Landed on branch `bed-regions` (8 commits, `866f230`…`efe1bf9`) and
merged to `main` (`6ecd22c`, `--no-ff`).

## Domain intent

The pipeline analyzed the whole genome with no way to (a) call only a
few regions of interest, or (b) skip known-pathological regions. This
feature lets the user pass a BED file to either stage; the rest of the
genome is ignored, and because the reference (`.fai`) and the
alignment/`.psp` files are indexed, the pipeline seeks to the regions
rather than walking the genome.

The organising principle (settled with the PM up front): **there is one
code path, and it always operates on a region set.** When the user
supplies a BED, that is the set; when they do not, the set is one
full-length span per contig. "Whole genome" is not a special case — it
is the region set whose every span covers an entire contig. This
deleted the whole-file streaming branch rather than adding a parallel
one.

## Plan (as executed)

Five incremental phases, each fmt/clippy/doc-clean with the full suite
green before the next:

- **Phase 1** (`866f230`) — the `regions` module: `RegionSet` + BED
  parser, pure logic, no callers.
- **Phase 2** (`9a161b6`) — `AlignmentMergedReader::query` gains a
  sub-contig range (`ContigInterval`).
- **Phase 3.1** (`bb79910`) — `load_pileup_inputs`: index pre-flight +
  per-input handle/metadata loader for the indexed pileup path.
- **Phase 3** (`1b99ef0`) — region-driven pileup core + CLI wiring;
  retire the streaming path.
- **Phase 4** (`664f04a`) — var-calling region restriction in the
  cohort driver.
- **Phase 5** (`21dac5b`) — command line + analysis-regions provenance
  in the `.psp` header.

## Design

**`regions` module** ([src/regions.rs](../../../../src/regions.rs)). A
top-level peer (consumed by both pipeline stages, so it depends on no
stage's contig type). `Region` is a 1-based **inclusive** span on one
contig; `RegionSet` is a sorted, overlap-merged list, built either from
a BED (`from_bed_reader` / `from_bed_path`) or as one full-length span
per contig (`whole_contigs`). The neutral `ContigBounds { name, length }`
input lets each caller adapt its own contig list. `regions_for(chrom_id)`
returns the per-contig slice (binary search).

The internal 1-based-inclusive convention was **verified, not assumed**,
to match what `PspReader::region_records` and the pileup walker already
use; the BED parser converts the 0-based half-open BED span
`[s, e)` → `[s+1, e]` at the single point where BED text is read.

**Pileup** ([src/pop_var_caller/cli.rs](../../../../src/pop_var_caller/cli.rs)).
`run_pileup` now: loads metadata + index handles
([`load_pileup_inputs`](../../../../src/bam/alignment_input.rs), pre-flight
build-or-error), builds the region set (BED or whole-contig default),
then iterates it — per region it opens an indexed reader
([`AlignmentMergedReader::query`](../../../../src/bam/alignment_input.rs)
with a `ContigInterval`), runs the BAQ→walker chain over the overlapping
reads ([`with_stage1_chain`](../../../../src/pop_var_caller/stage1_pipeline.rs),
extracted from the old streaming helper so it works over any opened
reader), and writes only the in-region columns
([`drive_region_into_writer`](../../../../src/pileup/per_sample/pileup_to_psp.rs))
to one shared `PspWriter`. Reads straddling a region edge are fetched
(BAQ keeps its flanking context) but their out-of-range columns are
clamped off; regions are disjoint, so each column is written once. The
whole-file streaming path is gone — pileup now always uses the index.
Per-region counters are totalled via new `merge` helpers on
`FilterCounts` / `RunSummary` / `BaqSkipCounts`.

**var-calling** ([src/var_calling/driver.rs](../../../../src/var_calling/driver.rs)).
The cohort driver already works off per-chromosome "covered intervals"
that both the block producer and the DUST-ahead pool walk in lockstep.
The restriction plugs in at exactly that point: `build_dust_plans`
intersects each chromosome's covered intervals with the region set
(`restrict_intervals_to_regions`, a two-pointer sweep), so block
production, DUST, grouping, and VCF emission all flow from the narrowed
intervals. Threaded as `Option<&RegionSet>`; when `None` (no `--regions`)
the intersection is skipped entirely, so the whole-genome path runs the
**identical** code and is byte-identical by construction.

**Provenance** ([src/psp/header.rs](../../../../src/psp/header.rs)).
`WriterProvenance`/`ParsedWriter` gain `command_line` (the full argv,
always recorded; serialised as the TOML `command-line` key with
`#[serde(default)]` so older `.psp` files still parse). With `--regions`,
the header also records `regions_bed` (basename) and `regions_count`.

## Assumptions / decisions

No silent assumptions — the ones the plan left open were settled with
the PM and are recorded here:

- **BED handling.** Overhang past the contig end is clamped to the
  length; a start past the contig end, `end <= start`, an unknown
  contig, and malformed lines are typed `BedError`s carrying line
  numbers. Overlapping **and directly-adjacent** spans are merged (exact
  interval union — identical base selection, fewer seeks); genuine gaps
  are preserved. `track`/`browser`/`#` lines and extra columns are
  ignored; tokenising is whitespace-lenient.
- **Index policy: error by default, opt-in build.** A missing
  `.crai`/`.csi`/`.bai` hard-errors naming the path; `--build-map-file-index`
  builds it in place. The reference `.fai` is required and never built.
- **Behaviour change: pileup now always requires an index** (even
  whole-genome), and mixed-format / unknown-extension inputs now surface
  via the index pre-flight rather than `new()`. Deliberate, consistent
  with the one-path principle.
- **Whole-genome output is *not* byte-identical to `main` — and the
  branch is more correct.** Investigated thoroughly (separate report:
  [bed_regions_perf_and_byte_identity_2026-06-03.md](../bed_regions_perf_and_byte_identity_2026-06-03.md)).
  The reader returns byte-identical reads; the difference is a **latent
  walker bug in `main`** that drops the last reads' tail columns at every
  chromosome transition. The per-region pileup ends each contig via
  end-of-input and emits them, recovering ~one read length of coverage
  per contig boundary. The recovered positions sit outside the GIAB
  evaluation BED, so accuracy is unchanged.

## Changes made

- **New:** [src/regions.rs](../../../../src/regions.rs) (RegionSet + BED
  parser + `BedError`); `regions` registered in
  [src/lib.rs](../../../../src/lib.rs).
- **Reader:** [src/bam/alignment_input.rs](../../../../src/bam/alignment_input.rs)
  — `ContigInterval` + sub-contig range on `query`; `overlaps_record`;
  `PileupInputs` + `load_pileup_inputs`; `FilterCounts::merge`.
  [src/bam/cram_input.rs](../../../../src/bam/cram_input.rs) /
  [bam_input.rs](../../../../src/bam/bam_input.rs) — per-region narrowing +
  overlap filter in the indexed CRAM/BAM iterators.
  [src/bam/errors.rs](../../../../src/bam/errors.rs) — `AlignmentInputError::AlignmentIndex`
  bridge.
- **Pileup:** [src/pop_var_caller/cli.rs](../../../../src/pop_var_caller/cli.rs)
  (region loop, `--regions` + `--build-map-file-index`, provenance);
  [stage1_pipeline.rs](../../../../src/pop_var_caller/stage1_pipeline.rs)
  (`with_stage1_chain`); [pileup_to_psp.rs](../../../../src/pileup/per_sample/pileup_to_psp.rs)
  (`drive_region_into_writer`); `RunSummary::merge`
  ([walker/driver.rs](../../../../src/pileup/walker/driver.rs)),
  `BaqSkipCounts::merge` ([baq_stream.rs](../../../../src/pileup/per_sample/baq_stream.rs)).
- **var-calling:** [src/var_calling/driver.rs](../../../../src/var_calling/driver.rs)
  (`restrict_intervals_to_regions`, threaded `Option<&RegionSet>`);
  [src/pop_var_caller/var_calling.rs](../../../../src/pop_var_caller/var_calling.rs)
  (`--regions`, region-set build, `BedError` bridge).
- **Provenance:** [src/psp/header.rs](../../../../src/psp/header.rs)
  (`command_line` on the wire + parsed/writer types).
- **Diagnostic:** [examples/diag_reader_compare.rs](../../../../examples/diag_reader_compare.rs)
  (used to prove reader equivalence).

## Tests added / updated

- **regions** (unit, 23+3): BED coordinate conversion, sort/overlap/
  adjacent merge, gap preservation, no cross-contig merge, line skipping,
  clamping, every `BedError` path with line numbers; `regions_for`.
- **reader** (unit, 7): CRAM/BAM sub-region selection, start-edge
  spanning, whole-contig≡unbounded, empty range; `load_pileup_inputs`
  build-missing-index + missing-index-error.
- **driver** (unit, 5): `restrict_intervals_to_regions` clip/split/
  no-overlap/empty.
- **walker seam** (unit, 1): `drive_region_into_writer` clamp + shared
  writer.
- **integration:** pileup `--regions` (BED [5,9] → exactly those
  positions = whole-genome clamped) + header provenance assertions;
  var-calling `--regions` (cohort SNPs at 11 and 30; BED [11,13] →
  only pos 11); existing fixtures opt into `build_map_file_index`; the
  two error-path tests updated for the pre-flight-sourced variants.

## Validation results

All in the dev container unless noted.

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features
  -- -D warnings` clean; `cargo doc --no-deps --lib` clean.
- `cargo test --lib`: 1101 pass (1 ignored); all integration suites
  pass (cohort 12, pileup 7, others unchanged).
- **Benchmark (host, GIAB HG002 bottle vs truth)** — branch vs `main`:
  accuracy **identical** (SNP F1 0.9037, indel F1 0.2068); pileup time
  373 s vs 372 s; peak RSS **−71 %** (993 MB vs 3431 MB, `--no-baq`
  isolation). Details in the companion report.

## Tradeoffs and follow-ups

- **Reference fetcher per region.** The pileup rebuilds a streaming
  `MultiChromStreamingRefFetcher` per region, which re-streams from the
  contig start. Harmless for the whole-genome default (one per contig)
  but wasteful for many sub-contig regions on a large contig; an
  indexed `.fai`-seek fetch would fix it. Deferred.
- **`var-calling` `--regions` performance** not yet benchmarked on a
  sparse BED (the headline benchmark exercised whole-genome).
- **Latent walker bug in `main`** (chrom-transition tail-drop) is now
  *dormant* on the merged code: pileup is single-contig per query, so
  the walker never transitions contigs. A general fix (drain the tail
  before the transition flush) would make any future multi-contig walker
  use correct; tracked but not required.
- **Indel recall** (0.1154) is a pre-existing caller property, identical
  before and after this change.
