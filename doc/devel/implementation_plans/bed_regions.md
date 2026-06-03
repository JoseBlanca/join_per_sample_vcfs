# Analysis regions from a BED file

Implementation plan for restricting the calling pipeline to a set of
genomic regions supplied as a BED file. Motivated by the
[TODO note](../TODO.txt): *"let's add regions to analyze using a bed
file. For instance to do only some selected regions or to avoid
pathological regions."*

**Problem.** Today both pipeline stages process the whole genome.
`pileup` stream-walks every input alignment file end to end
([`AlignmentMergedReader::new`](../../src/pop_var_caller/stage1_pipeline.rs#L107)),
and `var-calling` iterates every chromosome of the cohort `.psp` files.
There is no way to (a) analyze only a handful of regions of interest, or
(b) skip known-pathological regions — and when only a few regions are
wanted, walking the entire genome is pure waste even though both the
reference (`.fai`) and the alignment files (`.crai`/`.csi`/`.bai`) are
indexed and could be seeked directly.

**Fix.** Let the user pass `--regions <bed>` to either stage. Restrict
all work to those regions, and — because the inputs are indexed — *seek*
to each region rather than walking the genome.

## Design principle: there is always a region set

The simplifying decision (agreed up front): **we do not add a
region-restricted code path alongside the whole-genome one. There is
only one path, and it always operates on a region set.** When the user
supplies a BED, that is the set. When they do not, we synthesize one
full-length interval per contig from the contig list. "Whole genome" is
not a special case — it is the region set whose every interval spans an
entire contig.

This deletes the whole-file streaming branch rather than adding a
parallel one, and it maps cleanly onto machinery that already exists:

- [`AlignmentMergedReader::query`](../../src/bam/alignment_input.rs#L792)
  already does indexed, **per-contig** random access (k-way merge across
  a sample's CRAMs, coordinate-sorted). A default full-contig region is
  exactly `query(contig)`; a BED sub-region is `query(contig, start,
  end)` — the only extension needed is a position range.
- The `.psp` reader already exposes region random access:
  [`region_records(chrom_id, start, end)`](../../src/psp/reader.rs#L383)
  and [`BlockColumnReader::seek_to`](../../src/psp/reader.rs#L954), both
  binary-searching the internal block index. A full-chromosome region is
  just `start = 0, end = len`. `var-calling` already drives these in
  chunks via [`column_span_reader`](../../src/var_calling/column_span_reader.rs#L93).
- [`preflight_alignment_indexes`](../../src/bam/index_preflight.rs#L199)
  already implements the index policy we want: detect a missing
  `.crai`/`.csi`/`.bai` and either **build it in place** or
  **hard-error**, gated by a `build_if_missing` flag. The `.fai` is
  already required (typed `MissingFastaIndex`).

So the bulk of this feature is *glue*, not new seek machinery.

## Decisions (settled)

- **Both stages, independent flags.** `--regions <bed>` on both `pileup`
  and `var-calling`. Each stage builds its own region set from its own
  contig source (`.fai` for pileup, PSP header for var-calling). A user
  can restrict at pileup time (smaller `.psp`), or call a region subset
  from a whole-genome `.psp`, or both.
- **Index policy: error by default, opt-in build.** A missing
  `.crai`/`.csi`/`.bai` or `.fai` hard-errors with a message naming the
  expected path; building happens only behind an explicit flag (reuse
  the existing `--build-map-file-index` pattern). We never write files
  next to user inputs unasked.
- **Interval normalization: sort + merge overlapping/adjacent only.** No
  gap-coalescing threshold — emit exactly the requested positions. (A
  merge-gap knob can be added later if fragmented BEDs prove a seek
  bottleneck.)

## Consequence to accept

Routing pileup through the indexed path means **pileup now always
requires an alignment index + `.fai`**, even for a plain whole-genome
run that today streams without one. This is consistent with the "one
way" principle and the error-by-default policy; the index build is a
one-time, ~one-pass cost. It is a deliberate behavior change, not a
surprise.

## What already exists vs. what is new

| Piece | Status |
|---|---|
| Indexed per-contig alignment read (`query`) | exists, **no prod caller** |
| Sub-contig range on `query` | **new** (small extension) |
| Index build-or-error preflight | exists |
| `.psp` region read (`region_records`/`seek_to`) | exists, used by var-calling |
| BED parse + `RegionSet` type | **new** |
| `--regions` CLI flag (both stages) | **new** |
| Pileup region loop + retire streaming `new()` | **new** |
| var-calling ∩ region set | **new** |
| Region provenance in PSP header | **new** |

## Edge cases / correctness

- **Coordinate convention.** BED is 0-based half-open; our internal and
  VCF coords are 1-based. Convert at exactly one boundary (the BED
  parser), and unit-test it, or every region edge will be off by one.
- **Reads spanning a region edge & BAQ context.** We must *fetch* every
  read that overlaps a region (noodles query does this), but only *emit*
  pileup columns at positions inside it. BAQ's HMM needs the read's full
  span — fine, since the whole read is fetched. The clamp is on emitted
  positions, not on fetched reads.
- **Variant/indel straddling a boundary.** Inclusion rule: a variant is
  in-region iff its **anchor (start) position** falls in a region.
  Matches the column-emission clamp; conventional.
- **Contig-name consistency.** BED contig names must validate against
  the contig list; an unknown contig is a typed error (`query` already
  rejects unknown contigs).
- **Region-restricted `.psp` is partial.** It has no data outside its
  regions. Record the regions in the PSP header so the file
  self-describes, and confirm var-calling never reads "no data" as
  "hom-ref" (it should not today).
- **Whole-genome equivalence.** Whole-genome-via-per-contig-`query` must
  produce a **byte-identical** `.psp` to today's streaming path (it
  should — same coordinate-sorted reads). Verify both output equality
  and that there is no wall-time regression vs streaming.

## Phasing (incremental — pause between phases)

**Phase 1 — `RegionSet` type + BED parser.** A shared module (e.g.
`src/regions/`) with a genomic-interval type and a `RegionSet` (sorted,
overlap-merged, indexed by contig). A BED parser (0-based half-open →
internal) and a constructor that synthesizes the full-contig default
set from a contig list. Validate contig names against a provided list.
Pure logic, fully unit-tested; no CLI or reader wiring yet.

**Phase 2 — sub-contig range on `query`.** Extend
[`AlignmentMergedReader::query`](../../src/bam/alignment_input.rs#L792)
to accept an optional `[start, end)` within the contig (noodles supports
a region range natively). Whole-contig default = the existing behavior.
Unit-test the range path against the existing multi-contig fixtures.

**Phase 3 — pileup onto the region loop.** Wire `--regions` into the
pileup CLI; run index preflight (error by default, opt-in build); build
the region set (BED or full-contig default); iterate it via `query`,
clamping emitted columns to each interval. Retire the streaming `new()`
path. Verify byte-identical whole-genome output + no perf regression.

**Phase 4 — var-calling ∩ region set.** Wire `--regions` into the
var-calling CLI; intersect the existing chunk/region iteration with the
region set (clamp chunk edges to interval boundaries). Default =
full-chromosome set from the PSP header.

**Phase 5 — provenance + benchmarks.** Record the analyzed regions in
the PSP header parameters. Add a small region-restricted benchmark and
confirm the seek path is faster than whole-genome for a sparse BED.

## Testing

- Unit: BED coordinate conversion, sort/merge, full-contig default,
  unknown-contig rejection.
- Integration: a small multi-contig fixture with a BED selecting a
  sub-region of one contig; assert pileup emits columns only in-region
  and var-calling calls only in-region; assert no-`--regions` output is
  byte-identical to the pre-feature baseline.
