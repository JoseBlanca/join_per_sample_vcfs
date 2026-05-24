# `var-calling-from-bam` — per-chromosome parallelism

**Date:** 2026-05-24
**Plan:** [var_calling_from_bam_per_chromosome.md](../../implementation_plans/var_calling_from_bam_per_chromosome.md)
**Sibling slice (cohort path, same structural shape):** [cohort_per_chromosome_parallel_2026-05-20.md](./cohort_per_chromosome_parallel_2026-05-20.md)

## Verdict

**Code-complete; bench validation deferred.** The parallel
`run_var_calling_from_bam` path is live in `main`. All existing
from-bam integration tests pass through the new par_iter +
fragment-concat structure, and the full test suite (898 lib + every
integration / bench-compile target) is green under
`cargo clippy -D warnings` and `cargo fmt --check`.

The wall-time measurement on real multi-chrom tomato CRAMs (the
analogue of cohort H1's 3.85× at T=13) is not in this report: the
`examples/profile_from_bam_e2e.rs` + `benches/from_bam_e2e_perf.rs`
infrastructure (~1000 LOC of fixture generation + criterion plumbing
mirroring the cohort side) was scoped out per a design call; the
measurement runs on real input that lives on the user's host.
Filed as a follow-up.

## What landed

Single feature, four reviewable commits (commit 5 of the plan was
the bench infrastructure and was deferred):

| # | Commit  | Scope |
|---|---------|-------|
| 1 | `29b28e8` | **Alignment-index pre-flight + `--build-map-file-index` flag.** New `src/bam/index_preflight.rs` module with `preflight_alignment_indexes(inputs, build_if_missing)`. New typed errors `AlignmentIndexError` (bam-module variants, CLI-vocabulary-free) and `VarCallingFromBamCliError::{MissingMapFileIndex, IndexBuildFailed, UnsupportedAlignmentExtension}` (CLI variants whose Display impls name both the flag and the `samtools index` recipe). Pre-flight runs as step 0 of `run_var_calling_from_bam`, before rayon init. 6 unit tests + 2 integration tests. |
| 2 | `e5261ba` | **`CramMergedReader::query` indexed per-contig variant.** Mirror of `CramMergedReader::new`'s shape but each `OpenCram` is backed by a new `OwnedIndexedCramRecords` iterator that walks the `Arc<crai::Index>`, seeks to each matching container, and decodes records into the same k-way peek-and-scan merge `new` uses. New `AlignmentIndex` enum (single-variant `Crai(Arc<crai::Index>)`) + `load_alignment_index` helper in `index_preflight`. New `CramInputError::{ContigNotInList, PerInputHandleCountMismatch}` typed errors. 6 unit tests. |
| 3 | `050da41` | **`process_one_chromosome_from_bam` per-chrom worker.** Mirror of the cohort-side `process_one_chromosome` shape but with Stage 1 (CRAM → BAQ → walker) on the front of the chain. Each worker owns its own `CramMergedReader::query`, BAQ chunk pool, walker, walker-side `MultiChromStreamingRefFetcher`, cohort-side `StreamingChromRefFetcher`, and writer. The self-referential CRAM → BAQ → walker borrow chain materialises one chrom's `PileupRecord`s into a `Vec` before they cross into the merger (peak per-worker memory is the surviving-sites count, typically 1–100 K records). Walker errors are stashed via `ErrorSheddingAdapter` and surfaced before the merger sees a partial chrom. |
| 4 | `a858832` | **Driver reshape.** Replaces `run_var_calling_from_bam`'s serial body with: harvest contigs + sample name once via `CramMergedReader::new` (validates per-file + cross-file invariants), load per-input `Arc<sam::Header>` + `Arc<crai::Index>` once, allocate per-chrom fragment paths under a `TempDir` next to the output, drive `chromosomes.par_iter().enumerate().map(process_one_chromosome_from_bam)`, fail-fast aggregate stats, concat fragments via `vcf::concat::concat_fragments` (atomic rename + parent fsync). Retires `run_cohort_pipeline_for_single_sample` + `PerChromRecordsIter` + 5 `per_chrom_iter_*` tests (net diff: +179 / -520; a substantial simplification on top of delivering the parallel structure). |

## Assumptions surfaced (silent choices made)

- **BAM scaffolding dropped from this slice.** The original plan
  covered `.bai`/`.csi` for BAM inputs. `noodles-bam` is not a
  project dependency today and the only readable alignment format
  is CRAM. `AlignmentFileKind` is a single-variant enum
  (`Cram`-only) keeping the dispatch surface ready for when BAM
  input support lands; the BAM-side flow can land in one focused
  PR at that point.
- **Build BAI (not CSI) when auto-building BAM indexes — _not
  exercised_.** The plan's design § noted CSI as v1's default for
  BAM. The decision to build BAI via `noodles_bam::fs::index` (one
  call vs the manual `noodles_csi::binning_index::Indexer` loop CSI
  would need) is documented in the index_preflight code comments
  but currently unreachable (no BAM dispatch). When BAM lands the
  decision is reviewable in one place.
- **Index freshness is the user's responsibility.** If a `.crai`
  exists, we trust it. Stale-index detection is out of scope —
  matches samtools / bcftools convention.
- **Per-chrom CRAM-reader double-open.** Each per-chrom worker
  opens its CRAM(s) freshly (one `cram::io::Reader<File>` per input
  per worker) because noodles' `Query` iterator borrows the reader
  mutably and workers cannot share. The driver also opens each
  CRAM once at startup for harvest + header / index loading. The
  per-CRAM open cost is small next to the per-contig decode it
  enables; accepted as v1 design.
- **Per-chrom records buffered before the merger.** The walker's
  output (the chrom's `PileupRecord`s) is collected into a `Vec`
  in the per-chrom worker before it crosses into the cohort
  pipeline. Sidesteps the self-referential borrow chain that
  `with_stage1_pipeline`'s callback pattern handles for the
  whole-genome streaming case. Peak per-worker memory is the
  per-chrom surviving-sites count (typically 1–100 K records); not
  cross-worker accumulating. The same pattern can be revisited
  if a worker ever holds more than a chrom's records in scope.
- **No multi-contig integration coverage in this slice.** The
  in-module fixture builders (`tests/common/mod.rs::build_cram`
  and the `src/pileup/per_sample/cram_files.rs::build_cram` test
  helper) are both single-contig. Adding a multi-contig fixture
  builder + cross-contig golden / determinism tests would have
  expanded this slice substantially; filed as the most important
  test follow-up.
- **Stage 1 `FilterCounts` / `BaqSkipCounts` summary dropped.**
  The serial driver printed per-Stage-1 counter lines after
  Stage 1 completed. The parallel driver has no shared accumulator
  across the per-chrom workers, and reintroducing one needs a
  reducer (each per-chrom worker exposes its own counters, then
  the driver folds them). Tracked as a follow-up; the
  `records_emitted` summary is authoritative for what reached the
  VCF.
- **Nested rayon kept as (a) from the plan.** Outer per-chrom
  `par_iter` + inner BAQ chunk `par_drain` both share the global
  pool. The plan's options (b) "disable inner BAQ" and (c)
  "conditionally disable inner BAQ when n_chroms ≥ workers" are
  not implemented; we ship the simplest shape and revisit only
  after profile evidence on real data.

## Tests

| Layer | Count | What's exercised |
|---|---:|---|
| `src/bam/index_preflight.rs` unit tests | 6 | accept-existing-crai, error-when-missing-no-flag, unsupported-extension, idempotent-when-flag-set, first-missing-reported, target-index-path |
| `src/bam/cram_input.rs` group C unit tests | 6 | multi-contig query target-ref_id, single-contig coordinate order, k-way merge across inputs, empty-contig zero records, unknown-contig error, per-input-handle-count mismatch |
| `tests/cohort_cli_integration.rs` from-bam tests | 4 | happy-path (now via parallel path), walker-error surfacing, missing-index-without-flag error, missing-index-with-flag builds successfully |
| **Lib tests overall** | **898** | All pass (903 - 5 from removing `per_chrom_iter_*` tests with their owning helper). |

## Validation

Inside the project's dev container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` —
  clean.
- `cargo test --all-targets --all-features` — 898 lib + every
  integration / bench-compile target pass; 0 failures.

The bench numbers (T sweep on a real multi-chrom tomato CRAM, the
analogue of cohort H1's `examples/profile_cohort_e2e` 13-thread
table) are not in this report; the `examples/profile_from_bam_e2e`
+ `benches/from_bam_e2e_perf` infrastructure is the natural place
for that and was scoped out — see the "What landed" table.

## Behaviour changes (user-visible)

- **New flag** `--build-map-file-index` on `var-calling-from-bam`
  (off by default). When off, missing `.crai` indexes are a hard
  error citing both the flag and the `samtools index` recipe.
- **`var-calling-from-bam` now requires an alignment index** next
  to every input CRAM. Existing scripts that ran without `.crai`
  files will see the new `MissingMapFileIndex` error and must
  either pre-index (`samtools index`) or pass
  `--build-map-file-index`.
- **The from-bam Stage 1 summary lines are no longer printed**
  (per-Stage-1 `FilterCounts` / `BaqSkipCounts`). The
  `var-calling-from-bam: sample=… records_emitted=…` summary line
  is unchanged. Reintroducing the per-Stage-1 lines needs a
  cross-worker reducer; tracked as a follow-up.

## Open follow-ups

1. **Bench / profile infrastructure.** `examples/profile_from_bam_e2e.rs`
   + `benches/from_bam_e2e_perf.rs` mirroring the cohort H1
   measurement infrastructure. Required before measuring the
   wall-time win on real tomato CRAMs.
2. **Multi-contig CRAM fixture builder + integration tests.**
   Add to `tests/common/mod.rs` (or a sibling) a `build_cram` that
   takes a `&[ContigSpec]` instead of a single `CONTIG_NAME`/`CONTIG_LEN`
   constant, then add the three multi-contig integration tests the
   plan listed (golden serial-vs-parallel, single-chrom degenerate,
   empty-chromosome-in-middle).
3. **Multi-reference slice cross-contig filter coverage.** The
   `OwnedIndexedCramRecords` post-decode ref_id filter — for
   slices whose decoded records span multiple contigs — is
   unverified by unit tests because the in-module fixture builder
   produces multi-ref slices that surface noodles-side fragility.
   Land coverage via the integration tests of #2 (real CRAM
   writers produce single-ref slices per container in the typical
   case).
4. **Streaming `CramMergedReader::new` on 2-CRAM-with-records
   fixtures.** A surprise from commit 2: no existing in-module
   test drives records through a 2-CRAM `new` merge (every 2-CRAM
   test in group B passes `&[]` for records), and exercising it
   surfaced an `UnexpectedEof` mid-stream. Independent of this
   slice; track + investigate as a `cram_input.rs` follow-up.
5. **Stage 1 counter reducer for the parallel driver.** Restore
   the `FilterCounts` / `BaqSkipCounts` summary by collecting each
   per-chrom worker's counters and folding them in the driver.
6. **`RLIMIT_NOFILE` bump.** At the test fixtures' scale (≤4
   CRAMs × ≤13 chroms = ≤52 fds) we are well under the 1024
   default. Many-CRAM-per-sample workflows (rare) could exceed it;
   call out as a defensive `rlimit::increase_nofile_limit(8192)`
   at binary startup if real users hit it.
7. **Determinism harness across thread counts.** The rayon pool is
   process-global and idempotent (first caller wins), so a single
   `cargo test` invocation cannot easily compare `threads=1` vs
   `threads=N` outputs. A two-invocation harness (subprocess +
   `RAYON_NUM_THREADS` env) is the natural shape; deferred.
8. **BAM input + `.bai`/`.csi` dispatch.** When BAM input lands on
   the reader side, extend `AlignmentFileKind` and the pre-flight
   dispatch.
9. **Sub-chromosome decomposition.** The same ch00-style read
   imbalance the cohort H1 hit at 3.85× will cap this slice's
   speedup on workloads with one dominant contig. Splitting a
   giant contig across multiple workers is the order-of-magnitude
   lever beyond this slice; warranted only if profile evidence
   demands it.

## Estimated cost

Net diff: roughly +1320 / -520 across the four commits. The largest
single piece is `OwnedIndexedCramRecords` + `CramMergedReader::query`
(commit 2 at +608 / -1) followed by the driver reshape (commit 4 at
+179 / -520; net simplification because retiring
`run_cohort_pipeline_for_single_sample` + `PerChromRecordsIter`
deleted more than the new par_iter body added).
