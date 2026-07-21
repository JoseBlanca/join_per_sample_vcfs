# ng read input — implementation plan

**Status:** draft, 2026-07-20. The build order for **step 1's input edge**: `src/ng/read/input/` —
opening and validating one alignment file, serving a region as an ordered filtered stream, and
merging a sample's k files into one. Design is settled in two spec/arch pairs:
[`spec/alignment_file.md`](../spec/alignment_file.md) + [`arch/alignment_file.md`](../arch/alignment_file.md)
(one file) and [`spec/sample_reads.md`](../spec/sample_reads.md) +
[`arch/sample_reads.md`](../arch/sample_reads.md) (the sample), under the shared arch docs
([step interfaces](../arch/ng_step_interfaces.md), [module layout](../arch/module_layout.md)).
**One plan covers both**, because the sample layer is a thin cap on the file layer and shares its
module and its fixtures.

This roadmap turns that design into build order; it is **not** a place for new design. Both specs
record every open question as resolved.

---

## Scope

**In:** `src/ng/read/input/` (`mod.rs`, `open_bam.rs`, `region_query.rs`, `merge.rs`); the
`GenomePosition` addition to `types.rs`; two small extensions to `read/filtering.rs` (a probe-free
`ReadFilter` constructor and buffer hand-off); `AlignmentFile` with its validate-on-open gate, reader
pool, BAM and CRAM region queries and order guard; `check_assembly`; `SampleReads` with the cross-file
sample-name check, the k-way merge, and the merge-free single-file arm.

**Out (later plans / owners):**

- **Wiring into a CLI or driver** — no `pop_var_caller_exp` subcommand consumes this yet; the
  `type-regions` CLI plan ([`typed_regions_cli.md`](typed_regions_cli.md)) owns command surfaces.
- **The CRAM decoded-container cache** — deferred to a measured pass (spec `alignment_file.md` §8).
- **Batched-locus forward sweep** — same measured pass (spec `alignment_file.md` §8).
- **Cross-sample cohort assembly** — whatever drives the cohort holds N `SampleReads`
  (spec `sample_reads.md` §8).
- **A generic k-way merger** — explicitly not built; reopened only at the cohort layer
  (spec `sample_reads.md` §5).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The gate before the reader.** `AlignmentFile::open` establishes the invariant every later step
  leans on (`ref_id == ContigId`). Build and test the validation before anything reads a record, so
  no step downstream is ever written against an unvalidated handle.
- **Simplest impl first, as the oracle for the next.** BAM before CRAM: BAM's one-record-at-a-time
  read is the simpler container, and once it works it becomes CRAM's parity oracle (T8) — the same
  reads written both ways must produce the same ordered stream.
- **Isolate the silent failures.** Two steps here fail *quietly* — a region query that misses reads
  at a chunk edge, and an order guard that lets a regression through — and both produce wrong
  genotypes rather than a crash. Each lands as its **own commit with its oracle green**, so a
  `git bisect` can find it. Marked **own commit** below.
- **Verify against ground truth, not self-consistency.** The region query's oracle is the *existing
  whole-file* `BamRecordSource`: scan the entire file linearly, keep what overlaps the region, and
  the indexed query must return exactly that. It is an independent implementation of the same
  question, which is what makes it worth more than any assertion the new code could make about
  itself.
- **Reuse over rewrite.** `ReadFilter`, `NoodlesRawRecord`, `load_alignment_index`,
  `ContigList::first_disagreement`, `ReferenceInfo::contig_list` are called as-is; production's
  `segment_reader` is the **model**, not a dependency (both O3 decisions resolved to rebuild).
- **Incremental, with pauses.** One milestone, then stop for review.
- **Ungated / container builds.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); a native host build
  at completion.

## Preconditions (already in place — confirm before A1)

- **Step 1 (read filtering) is complete**: `ReadFilter` ([`filtering.rs:581`](../../../../src/ng/read/filtering.rs)),
  the `RecordSource`/`RawRecord` seam (`:310`, `:285`), `NoodlesRawRecord` (`:344`),
  `ReadFilterConfig`/`Counts`/`Error` (`:47`, `:117`, `:550`), and the whole-file
  `BamRecordSource`/`CramRecordSource` (`:379`, `:429`) — the last two are this plan's **oracle**,
  not its dependency.
- **`reference_info` is complete**: `ReferenceInfo` + `contig_list()`
  ([`reference_info.rs:71`, `:91`](../../../../src/ng/reference_info.rs)), and
  `VerificationHandle::join` (`:931`) for the deferred assembly check.
- **`RefSeq` is complete**: `RawRefSeq` ([`ref_seq.rs:180`](../../../../src/ng/ref_seq.rs)) with
  `InMemoryRefSeq` (`:240`) as the test reference.
- **Production helpers are `pub`**: `load_alignment_index` / `preflight_alignment_indexes`
  ([`index_preflight.rs:120`, `:199`](../../../../src/bam/index_preflight.rs)),
  `ContigList::first_disagreement` ([`fasta/mod.rs:69`](../../../../src/fasta/mod.rs), `pub(crate)`).
  **Note:** `extract_header` / `extract_single_sample_name`
  ([`alignment_input.rs:292`, `:378`](../../../../src/bam/alignment_input.rs)) are **module-private**
  — ng reads the `sam::Header` itself (arch `alignment_file.md` §5). Do not plan on calling them.

---

## The steps

### Milestone A — vocabulary, scaffold, and the two `filtering.rs` extensions

**A1. `GenomePosition` in `types.rs`.**  ✅
`pub struct GenomePosition { pub contig: ContigId, pub position: Position }` with the standard
derives; **`Ord` derived so field order is genome order**. Doc comment says what the value *is* (one
base, genome-wide) and that a bare `Position` does not identify a base. Unit test: sorting a shuffled
vector yields contig-major, position-minor order. A shared-vocabulary addition, so it lands alone.
*Source:* arch `sample_reads.md` §1.1, `ng_step_interfaces.md` §1.

**A2. Scaffold `read/input/`.**  ✅
`src/ng/read/input/{mod.rs, open_bam.rs, region_query.rs, merge.rs}` with `#[cfg(test)]` blocks;
`pub mod input;` in `read/mod.rs`. Declare `AlignmentFileError` and `IngestError` with their variants
and doc comments, no logic. *Depends:* A1. *Source:* arch `alignment_file.md` §Module home + §2,
arch `sample_reads.md` §2.

**A3. Extend `read/filtering.rs`: probe-free construction + buffer hand-off.**  ✅
Two additions, both ng-owned (not production edits): (a) a constructor that **skips** `new`'s
per-contig resolve probe — that probe is O(contigs) reference fetches, and building a `ReadFilter`
per region query would pay it ~10⁶ times, while the open gate proves something strictly stronger;
(b) construct-from / release-to parts so the record buffer and `ref_buf` can be handed in from the
pool and taken back, keeping a per-query filter allocation-free. `ReadFilter::new` is unchanged for
whole-file callers. Unit test: the probe-free constructor yields identical output to `new` on a fake
source. *Depends:* A2. *Source:* arch `alignment_file.md` §5 (decision 2), §1.2.

> **Checkpoint A:** `GenomePosition` sorts in genome order; the module tree compiles; `ReadFilter`
> has both constructors with identical behaviour. Pause for review.

### Milestone B — the open gate (validation before any read)

**B1. Header extraction, pure.**  ✅
Read `@HD SO`, the `@SQ` list (name + length + `M5`), and the `@RG SM` set off a noodles
`sam::Header` — ng's own, since production's are private. Returns the pieces the gate needs; no file
I/O, so it unit-tests against hand-built headers. Tests: missing/other `SO`; two distinct `SM`s;
`@SQ` with and without `M5`. *Depends:* A2. *Source:* spec `alignment_file.md` §3.1, arch §5.

**B2. `AlignmentFile::open` — the gate.**  ✅
Compose B1 with `load_alignment_index` and the `@SQ`↔reference comparison via
`ContigList::first_disagreement` against `ReferenceInfo::contig_list()`; capture `sq_md5s` by
`ContigId`; store the parsed index. Fail-fast, in this order: `SO` → `@SQ` → index → `SM`. Tests
**T1** (permutation caught — mutation-verify against a resolves-only check, which passes it), **T2**
(name/length/count mismatch each naming the right field and index; MD5 mismatch caught only when the
`ReferenceInfo` carries digests), **T3** (`SO` wrong/missing), **T12a** (two `SM`s in one file).
*Depends:* B1. *Source:* spec `alignment_file.md` §3.1, arch §3.

> **Checkpoint B:** a file is either validated or an error; the permutation hole reference_info.md §1
> named is closed and proven closed by mutation. No read has been read yet. Pause for review.

### Milestone C — serving a region (the silent-failure zone)

**C1. The reader pool.**  ✅
`ReaderHandle` (reader + `record_buf` + `ref_buf`) and the `Mutex<Vec<ReaderHandle>>` on
`AlignmentFile`; borrow-or-open, return on `Drop`. The index and header stay outside the handle,
shared — never re-parsed. Test: N sequential borrows open the file once; a returned handle is reused.
*Depends:* B2, A3. *Source:* spec `alignment_file.md` §3.3, arch §1.2.

**C2. `BamRegionSource`.**  ✅ **Own commit — do not bundle.**
Index query → chunks (`BinningIndex::query` on the already-parsed index), seek the borrowed handle
per chunk, fill the reused buffer, drop non-overlapping records **uncounted**, early-stop past `end`.
**Oracle (T5):** for a set of regions over a fixture BAM, the indexed query must return *exactly*
what a whole-file `BamRecordSource` scan filtered to the same region returns — same reads, same
order. This is the step whose failure is silent (a missed chunk edge is a wrong genotype, not a
crash), so it lands alone with that oracle green. *Depends:* C1. *Source:* spec §3.2/§3.3, arch §4.

**C3. `OrderVerified`.**  ✅ **Own commit — do not bundle.**
The adapter comparing each read's `GenomePosition` against the previous, hard-erroring
`OutOfOrderRead { previous, current }` on a **strict** decrease. State lives in the iterator, not the
handle. Tests **T4a** (planted position regression; mutation-verify that removing the check lets it
through), **T4b** (contig-order regression), **T4c** (equal positions are legal), **T4d** (querying
region B then region A is not a regression). The second silent failure: a guard that never fires
looks exactly like a guard that works. *Depends:* A1, C2. *Source:* spec §3.2, arch §4.

**C4. `reads_in_region(&self)` — the composed chain.**  ✅
Source → `ReadFilter` (probe-free, buffers from the handle) → `OrderVerified` → `RegionReads`,
returning the handle on `Drop` including on the error path; `counts()`. Tests **T9** (a #8
high-mismatch read is dropped — proves the full filter is composed in, not the cheap subset),
**T10** (truncated file → one `Err` then `None`), **T13** (many queries parse the index once — via a
counting wrapper). **Pin here:** the `counts()` mechanism under `&self` (per-handle tallies summed, or
atomics) — arch §7 left it as an impl-time confirmation. *Depends:* C2, C3. *Source:* spec §3.2, arch §4.

**C5. `CramRegionSource`.**  ✅
`.crai` walk with a **carried cursor or binary search** — not production's rescan-from-entry-0, which
is O(n) per locus ([`segment_reader.rs:1043`](../../../../src/bam/segment_reader.rs)) — container
decode, drop non-overlapping, container-level early stop. **Oracle (T8):** the same reads written as
BAM and as CRAM produce the same ordered `MappedRead` sequence through `reads_in_region`. *Depends:*
C4. *Source:* spec §3.2/§3.3, arch §4.

> **Checkpoint C:** region queries work for BAM and CRAM, agree with the linear-scan oracle and with
> each other, run the full step-1 filter, and reject an unsorted file. The index is parsed once.
> Pause for review.

### Milestone D — the deferred assembly check

**D1. `check_assembly`.**  ✅
The pure comparison of captured `@SQ M5` tags against a verified `ReferenceInfo`'s per-contig
digests; a contig is compared only when **both** sides carry one; returns
`AssemblyCheck { compared, total }`. Missing `M5` is **never** an error and never a warning. Test
**T2b**: mismatch → `AssemblyMismatch` naming contig and both digests; no tags → `compared == 0`,
`Ok`; partial → `compared` counts exactly the tagged contigs. Pure, so no fixture needed.
*Depends:* B2. *Source:* spec `alignment_file.md` §3.1, arch §3.

> **Checkpoint D:** wrong-assembly detection works and stays out of the startup path. **The file
> layer is complete.** Pause for review.

### Milestone E — the sample layer

**E1. `SampleReads::open` + the sample-name agreement.**  ✅
Open k files through `AlignmentFile::open`, then check they all name one sample; wrap any per-file
error with its `source_file_index`. Test **T12b** (two files, different `SM` → `SampleNameMismatch`
before any read). *Depends:* B2. *Source:* spec `sample_reads.md` §3.1, arch §3.

**E2. `MergedRegionReads` — the argmin merge.**  ✅ **Own commit — do not bundle.**
Heads as `Option<MappedRead>` with **keys held beside them** in a parallel array; linear argmin, ties
to the lowest file index; emit with `Option::take` and refill only the winning slot. The
same-file-twice check runs **only on a tie**, comparing `flag` before `qname`. Tests **T6**
(interleaving, correct `source_file_index`, tie-break to lower index, run-to-run identical output),
**T7** (same file twice → error at the *first* collision; same-position-different-`qname` reads both
survive — the case the cheap-first comparison order must not get wrong). A mis-ordered merge is the
third silent failure here. *Depends:* A1, C4. *Source:* spec §3.2, arch §4.

**E3. `SampleRegionReads` — the two arms.**  ✅
The `Single | Merged` enum (never `Box<dyn Iterator>`), `reads_in_region(&self)` building one chain
per file and a merge only when k > 1, plus `counts()` (per file, never summed) and `sq_md5s()`.
Test **T11**: a one-file sample yields exactly what the same file yields as one of two inputs with
the other empty — asserted through the same type, so the arms are provably indistinguishable.
**Pin here:** whether the per-file reference accessor needs `R: Clone` or `&R` — arch §7's other
impl-time confirmation. *Depends:* E1, E2. *Source:* spec §3.4, arch §3.

**E4. The merge's per-read budget bench.**  ✅
**T14:** a synthetic k-file stream asserting zero per-read allocations and no `MappedRead` clone in
the merge loop (counting allocator or `dhat`). Worth its own step because this failure is invisible
to every correctness test — a stray `clone()` costs throughput without changing an output byte.
*Depends:* E2. *Source:* spec §3.2 (the budget), arch §4.

> **Checkpoint E:** a sample's reads stream in coordinate order from one file or from several, with
> per-file counts and the cross-file guards. **Read input is complete.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | `GenomePosition` sorts contig-major; both `ReadFilter` constructors agree on a fake source |
| B | T1 (permutation, **mutation-verified** against a resolves-only check), T2, T3, T12a — all at `open`, before any read |
| C | **T5: indexed query ≡ whole-file linear scan filtered to the region** (independent oracle); T8 BAM≡CRAM; T4a–d order guard (T4a mutation-verified); T9 full filter; T10 fused; T13 index parsed once |
| D | T2b — the three `check_assembly` outcomes, pure, no fixture |
| E | T12b cross-file `SM`; T6/T7 merge order, tie-break, determinism and the same-file-twice error; T11 single-arm ≡ merged-arm; **T14 per-read budget bench** |

Fixtures throughout: tiny in-memory BAM/CRAM via noodles writers + an in-crate CSI/CRAI build
(production's `segment_reader`/`segment_merge` test modules and the `index_preflight` tests are the
templates), over a small multi-contig reference indexed by `reference_info`.

## Out of scope (next plans)

- **CRAM decoded-container cache** and **batched-locus forward sweep** — one measured performance
  pass over this module, once there is a real workload to measure (spec `alignment_file.md` §8).
- **Consuming this module** — the pileup walker and the STR tract extractor pull from
  `SampleReads`; each lands with its own step's plan.
- **Cross-sample cohort assembly** — holds N `SampleReads`; also where the generic-merger question
  reopens (spec `sample_reads.md` §5, §8).
- **Parallelism** — the pool and the `&self` signatures exist so this needs no redesign, but no
  threading lands here.
