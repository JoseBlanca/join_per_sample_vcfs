# Project status

<!--
ABOUT-PARAGRAPH-START — do not edit this paragraph.
Skills and agents are instructed to leave it untouched.
-->
> **About this project.** Multi-sample SNP caller:
> per-sample pileup → `.psp` artefact → DUST filter → variant grouping →
> per-group merger → posterior engine. The authoritative design document is
> [doc/devel/specs/calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md);
> read it before anything else. Companion design context is in
> [doc/devel/specs/design_principles.md](doc/devel/specs/design_principles.md).
> All work below is graded against that spec. For the AI assistant's
> per-skill instructions on reading and
> updating this file, see [doc/devel/ia/skills/](doc/devel/ia/skills/) —
> every skill defines a "Project status protocol" section.
<!-- ABOUT-PARAGRAPH-END -->

> **Current focus.** _Maintained by skills (last-completed) and the human
> project manager (next-task)._
>
> - **Last completed task (2026-06-17):** **ssr-pileup Mark-2 review fixes applied**
>   (branch `ssr-pileup-mark2`, [fixes_applied_2026-06-17_v2.md](doc/devel/reports/reviews/fixes_applied_2026-06-17_v2.md)).
>   Applied the Mark-2 code review: **all 3 Blockers** + **9 of 12 Majors** + 9 Minors. **B1** doc gate restored
>   (`cargo doc` green); **B2** length-inconsistent reads (empty-`QUAL`/over-consuming CIGAR) dropped+counted at
>   the fetch boundary (`LocusReads.malformed`→`n_filtered`) + `extract_region` `r_start` clamp — a test confirmed
>   noodles decodes `*` QUAL to an empty buffer, the exact pre-fix panic; **B3** added the missing byte-identity
>   tests (multi-chunk e2e + cap-bites e2e + reservoir-subset); **M2** the Stage-0 empty-flank gate (`finish_locus`);
>   **M1** split the catch-all `Io` into 5 operation-named variants; **M3** `reach_min_flank_bp` recorded in the
>   `.ssr.psp` header; **M5** byte-identity scoped to one target/toolchain; **M6/M7/M8/M9/M10/M11** + Mi1/Mi3/Mi4/Mi5/
>   Mi6/Mi7(local)/Mi-docs. **Deferred:** **M4** (HMM-model provenance — schema call), **M12** (`DpState` enum
>   refactor — now guarded by the B3 tests), Mi2/Mi8/Mi9/Mi10/Mi11/Mi12 + Nits. Gates: fmt/clippy `-D warnings`/doc
>   clean, **1129 lib + integration pass** (the only `--all-targets` failure is the pre-existing `psp_writer_perf`
>   bench panic). Audit trail `tmp/review_2026-06-17_ssr-pileup-mark2/`.
> - **Prior task (2026-06-17):** **ssr-pileup Mark-2 code review**
>   (branch `ssr-pileup-mark2`, [ssr_pileup_mark2_2026-06-17.md](doc/devel/reports/reviews/ssr_pileup_mark2_2026-06-17.md),
>   11 categories). First *correctness* review of the Mark-2 rebuild (`src/ssr/pileup/`). **Verdict:
>   Request-changes** — 3 Blockers, 12 Major, 14 Minor + Nits. Gates: fmt/clippy `-D warnings` clean,
>   40 ssr::pileup lib tests pass — but **`cargo doc` FAILS** (**B1**: fetch_reads.rs:17 links the deleted
>   `super::locus_record::aggregate`; `locus_record`→`locus_tally`, `aggregate`→`tally`). **B2**: `process_locus`
>   slices `read.qual` by a `seq.len()`-derived range with no length-consistency guard → an empty-`QUAL` record
>   (legal SAM; the SNP BAQ engine guards it at baq_engine.rs:117, the SSR path doesn't) panics the whole run;
>   an over-consuming CIGAR also panics via the unclamped `r_start`. **B3**: the byte-identity-for-any-thread-count
>   contract is untested on the only paths that can break it (reservoir eviction never bites — 1 read/locus,
>   cap 1000; `par_chunks` never splits — 3 loci < MIN_FETCH_CHUNK 64). Determinism *mechanism* verified sound
>   by inspection; diff matches the stated empirical-candidate intent. Majors: M1 `Io(#[from])` collapses 5 I/O
>   sites; M3/M4 MIN_FLANK_BP + HMM model constants shape output but aren't in the `.ssr.psp` header; M6
>   `..FilterCounts::default()` test masks new buckets; M8 `ref_to_read` indel branches untested. Audit trail
>   `tmp/review_2026-06-17_ssr-pileup-mark2/`. **(suggested follow-up: apply fixes — B1 first, then B2/B3;
>   answer Q1 [does Stage-0 emit empty-flank loci → M2 severity] and Q4 [drop vs error for malformed records].)**
> - **Prior task (2026-06-17):** **SSR Stage 1 (`ssr-pileup`) rebuilt as Mark-2**
>   (branch `ssr-pileup-mark2`, [ssr_pileup_mark2_2026-06-17.md](ia/reports/implementations/ssr_pileup_mark2_2026-06-17.md)).
>   Replaced the Mark-1 reference-anchored **rung** model with the **empirical-candidate** model: candidate
>   alleles are observed sequences, the reference is only a coordinate frame, no on/off-ladder, no Stage-1
>   likelihood. Per read: a per-Q Viterbi+traceback **delimits** the repeat region → first-quartile **quality
>   gate** (Phred 15) → tally **observed sequence → count** per locus. New `registry_ssr` schema
>   (per-observation `obs-count`/`obs-seq-len`/`obs-seq` Bytes; no version bump, pre-alpha). Built in 7 steps;
>   the old `src/ssr_mark1/` tree + its bench/example were deleted at the cutover. fmt/clippy `-D warnings`
>   clean; **1116 lib + integration + doctests, 0 failed**, incl. end-to-end catalog→BAM→`.ssr.psp` +
>   thread-count determinism. Design: [ssr_ladder_model.md](doc/devel/architecture/ssr_ladder_model.md),
>   [ssr_pileup_mark2.md](doc/devel/architecture/ssr_pileup_mark2.md). **Latest review:**
>   [ssr_pileup_mark2_2026-06-17.md](doc/devel/reports/reviews/ssr_pileup_mark2_2026-06-17.md) (Request-changes);
>   fixes [fixes_applied_2026-06-17_v2.md](doc/devel/reports/reviews/fixes_applied_2026-06-17_v2.md) (3 Blockers + 9 Majors applied).
>   Open (deferred from the fix run): **M4** HMM-model provenance tag; **M12** `DpState` enum refactor; Mi2/Mi8/Mi9/Mi10/Mi11/Mi12;
>   `WriterProvenance.input_crams` crate-wide rename; allele-diverse cap-bites + `BorderOffEnd` e2e tests.
>   Then: calibrate Q1/cap on real data; a Mark-2 bench; spec §4 amendment; Stage 2 (`ssr-call`).

---

## Pipeline stages

Stage descriptions are one-line reminders; the spec is authoritative.

### Stage 1 — per-sample pileup (BAM → `.psp`)

Stage 1 reads each BAM/CRAM once per sample and writes one `.psp` artefact.

#### Analysis regions from a BED file (`--regions`)
- **Status:** fixes-applied (2026-06-16, on branch `bed-regions-review`; re-review + fixes after the segment-reader retrofit); feature merged to `main` (`6ecd22c`), fix branch not yet merged
- **Plan:** [bed_regions.md](doc/devel/implementation_plans/bed_regions.md)
- **Impl report:** [bed_regions_2026-06-03.md](doc/devel/reports/implementations/bed_regions_2026-06-03.md)
- **Perf + byte-identity investigation:** [bed_regions_perf_and_byte_identity_2026-06-03.md](doc/devel/reports/bed_regions_perf_and_byte_identity_2026-06-03.md)
- **Cross-cutting:** Stage 1 (pileup) + var-calling. The pipeline always operates on a `RegionSet` — the `--regions` BED, or one full-length span per contig (whole genome). The whole-file pileup streaming path is retired; pileup now always seeks via the index.
- **Code:** [src/regions.rs](src/regions.rs) (`RegionSet` + BED parser); region-driven `run_pileup` ([src/pop_var_caller/cli.rs](src/pop_var_caller/cli.rs)) now opens pooled per-file readers once and k-way-merges per region via `SegmentMergedReads` ([src/bam/segment_merge.rs](src/bam/segment_merge.rs) + [src/bam/segment_reader.rs](src/bam/segment_reader.rs)) — the old per-region `AlignmentMergedReader::query` is **gone** (the segment-read-fetcher retrofit landed; cf. cli.rs:374-375); var-calling restriction in [src/var_calling/pipeline.rs](src/var_calling/pipeline.rs) (`restrict_intervals_to_regions`); command-line + regions provenance in [src/psp/header.rs](src/psp/header.rs).
- **Latest code review:** [regions_2026-06-16.md](doc/devel/reports/reviews/regions_2026-06-16.md) — **Approve-with-changes** (re-review of the current segment-reader path; 10 categories): 0 Blockers, 6 Major, 12 Minor + Nits. Top: M4 worker `JoinHandle`s dropped (a worker panic could silently truncate the read stream — only `thread::scope`'s re-raise saves it); M1 silent empty output on a non-selecting BED; M3 + M2 = prior M1/M3 still unfixed. Both this branch's cleanup commits verified behavior-identical. Audit trail `tmp/review_2026-06-16_regions/`.
- **Latest fixes-applied:** [fixes_applied_2026-06-16_v2.md](doc/devel/reports/reviews/fixes_applied_2026-06-16_v2.md) — 17 Applied + 2 Applied-with-adaptation, 4 Deferred, 1 Disputed (Mi11). All gates green (1165 lib + integration tests, clippy `-D warnings`, doc). Notably resolves the two carried-over 2026-06-03 Majors (counter-fold field-safety, `BedError::Io` context). Q1=empty-BED→hard error (`BedError::NoRegions`); Q2=`MappedReadSource` kept as a documented DI seam.
- **Prior review:** [bed_regions_2026-06-03.md](doc/devel/reports/reviews/bed_regions_2026-06-03.md) — Approve-with-changes on the now-superseded `query` impl; 0 Blockers, 6 Major, 9 Minor. HG002-bottle accuracy identical to `main` (SNP F1 0.9037), peak RSS −71%.
- **Latest perf review:** [perf_bed-regions_2026-06-16.md](doc/devel/reports/reviews/perf_bed-regions_2026-06-16.md) — Profile-first → **measured** on tomato1 (§8): a fragmentation sweep (80→8000 regions over identical 8 Mb = **+1.9% wall**) + a sampling profile (`probaln_glocal` ~70%; whole `--regions` read/merge path ~1–2%) **refute the per-region/read-path findings as a perf lever.** L1 (qname double-clone) + L3 (dead `Stage1Outputs` clones) applied as **perf-neutral cleanups**; rest closed/deferred. Real lever is BAQ (out of scope). Audit trail `tmp/perf_review_2026-06-16_bed-regions/`.
- **Open (from the 2026-06-03 review):**
  - **M1** — the three `merge()` counter-folds (`RunSummary`/`FilterCounts`/`BaqSkipCounts`) aren't field-addition-safe (exhaustive-destructure + 3 tests).
  - **M2** — `ContigInterval { start: 0, .. }` is constructible and aborts `query_interval`'s `Position::new(..).expect()` (checked constructor / drop `pub` fields).
  - **M3** — `BedError::Io` drops path/line and flattens its `io::Error` into `Display` (no `#[source]`).
  - **M4** — FASTA-repository construction duplicated ×3; cross-file validation copied between `new` and `load_pileup_inputs`.
  - **M5** — `--build-map-file-index` (index management) folded into this PR (split or document the coupling).
  - **M6** — whole-contig range collapse applied silently with no runtime trace.
  - **Minors:** per-region repeated work (FASTA repo rebuild + clones, undercuts "sparse BED is cheap"); pileup-vs-var-calling provenance asymmetry; empty-BED → silent empty output; `build_map_file_index` third-name / `stashed_upstream` noun; `#[from] AlignmentIndexError` collapses two origins; `#[allow(too_many_arguments)]` without justification. Plus the §8 missing tests (merge folds, `command-line` serde-default round-trip, adversarial BED coordinates, empty-BED pin).
  - **Efficiency follow-up:** indexed `.fai`-seek fetcher to replace the per-region streaming `MultiChromStreamingRefFetcher` rebuild (many sub-contig regions on a large contig).
- **Code-review fixes (2026-06-16, see [fixes_applied_2026-06-16_v2.md](doc/devel/reports/reviews/fixes_applied_2026-06-16_v2.md)):**
  - **Fixed:** M1 (empty BED → `BedError::NoRegions`), M2 + M3 (the carried-over 2026-06-03 Majors: `BedError::Io` context + counter-fold field-safety), M4 (explicit worker join + doc), M6 (cross-thread determinism test), and Mi1–Mi10 + nits. Mi11 disputed (kept as a documented DI seam).
  - **Still open — M5** (deferred): staged/inline producer interval-walk duplication (`pipeline.rs:431-595`) — regression-prone; extract a shared `drive_schedule` helper as its own byte-identity-gated effort.
  - **Still open — Mi12** (deferred): wire `cargo audit` + `deny.toml` (cargo-audit not installed; no new deps this branch).
  - **Still open — partial Mi4:** the per-line overhang-*clamp* warning (the mode-announce half landed); minor follow-ups Nit-bedlines / Nit-cache.
- **Resolved by the 2026-06-16 perf review (measured — see report §8):**
  - **Applied (perf-neutral cleanups, gates green):** L1 qname double-clone → move ([segment_merge.rs](src/bam/segment_merge.rs)); L3 dropped the unread `Stage1Outputs` `sample_name`/`contigs` fields ([stage1_pipeline.rs](src/pop_var_caller/stage1_pipeline.rs)). Neither moves wall (~1% paths) — kept as cleanups, not speedups.
  - **Closed — no gain:** L5/L6 (merge path ~1–2% of runtime); L7 region-parallelism (serial loop ~2% even at 8000 regions; cores already saturated).
  - **Deferred — untestable on pre-sliced CRAMs:** L2 `.crai` head binary-search, L4 container cache — need a full un-sliced CRAM + fragmented-BED workload (tomato1 CRAMs are pre-sliced; tiny `.crai`). Still mirrored by the segment-read-fetcher block's deferred items.
  - **Real lever (out of this feature's scope):** BAQ (`probaln_glocal` ~70% of pileup self-time) — see [perf_baq_2026-05-12.md](ia/reviews/perf_baq_2026-05-12.md).

#### Alignment-file input (CRAM + BAM)
- **Status:** fixes-applied (BAM slice); shipped (CRAM)
- **Plans:**
  - CRAM slice: [per_sample_caller_cram_input.md](doc/devel/implementation_plans/per_sample_caller_cram_input.md)
  - BAM slice: [bam_input_support.md](doc/devel/implementation_plans/bam_input_support.md)
- **Impl reports:**
  - CRAM slice: [per_sample_caller_cram_input_2026-04-29.md](doc/devel/reports/implementations/per_sample_caller_cram_input_2026-04-29.md)
  - BAM slice (2026-05-24): [bam_input_support_2026-05-24.md](doc/devel/reports/implementations/bam_input_support_2026-05-24.md)
- **Code:** [src/bam/](src/bam/) — `alignment_input.rs` (merge + filter + header validators, format-agnostic), `cram_input.rs` + `bam_input.rs` (per-format owned record-stream decoders), `index_preflight.rs` (CRAI / CSI / BAI detection + build), `errors.rs`.
- **Latest reviews:**
  - BAM slice (2026-05-24): [bam_input_support_2026-05-24.md](doc/devel/reports/reviews/bam_input_support_2026-05-24.md) — Approve-with-changes: 0 Blockers, 19 Major (M1–M19), 22 Minor (Mi1–Mi22), grouped Nits. Architecture sound (no `unsafe`, no shared mutable state, `'static + Send` iterators correct by construction); Major findings cluster around (a) rename-debt — CRAM vocabulary still in error displays + variant name `PileupCliError::CramInput` mislabels every BAM error; (b) missing `#[non_exhaustive]` on `AlignmentInputError` / `PileupCliError` / `BamIndex`; (c) `load_per_input_headers` reimplements opener work AND skips the CRAM-version gate that the per-format helper has — CRAM 4.x silently passes header-load; (d) `.csi`/`.bai` policy + CSI depth defaults live in function bodies, not named constants; (e) reliability test gaps for the new `AlignmentIndexFormatMismatch`, `UnsupportedExtension`, indexed-BAM chunk-walk, and `load_alignment_index` branches. Six missing-test specs in §8 of the report.
- **Latest fixes-applied (BAM slice):** [fixes_applied_2026-05-24.md](doc/devel/reports/reviews/fixes_applied_2026-05-24.md) — **Completed**: every Major (M1–M19) Applied; 13 of 22 Minors Applied; 7 Minors Deferred (Mi5 / Mi7 / Mi12 / Mi13 / Mi15 / Mi20 / Mi21) + M17's redesign half. 10 fix commits `9fc1df0` → `7a8569b`. **Total +10 new tests** (1 cram_input regression + 4 index_preflight `load_alignment_index` triplet + 1 alignment_input AlignmentIndexFormatMismatch + 3 bam_input [CSI arm, multi-chunk, truncated-BAM err] + 1 pileup_cli UnsupportedExtension + 1 cohort_cli MixedFormat lock-step). 924 lib + 45 integration tests pass; cargo fmt / clippy / test all clean.
- **Latest reviews / fixes (CRAM):** [per_sample_caller_cram_input_2026-04-29.md](doc/devel/reports/reviews/per_sample_caller_cram_input_2026-04-29.md), [fixes_applied_2026-05-01.md](doc/devel/reports/implementations/fixes_applied_2026-05-01.md)
- **Open (BAM slice — after the 2026-05-24 fix run):**
  - Wall-time validation on real multi-chrom BAMs (analogue of cohort H1's 3.85× at T=13 on tomato CRAMs). Picked up alongside the `examples/profile_from_bam_e2e.rs` + `benches/from_bam_e2e_perf.rs` infrastructure already deferred from the prior task.
  - **Deferred Minors** (from the fix run): **Mi5** (`process_one_chromosome_from_bam` rename — subcommand-name coupling), **Mi7** (`input_crams` PSP field rename — reach beyond `src/bam/`), **Mi12 + M17 redesign half** (lift `MixedAlignmentFileFormats` into a shared sub-enum), **Mi13** (`From<AlignmentIndexError>` 4-arm passthrough — design call), **Mi15** (`crate::bam` → `crate::alignment` module rename), **Mi20** (`OwnedIndexedBamRecords::next` phase extraction — cosmetic), **Mi21** (per-record `RecordBuf::default()` allocation — folds into parallelisation-tuning workstream), **Mi9 partial** (`tests/common::build_csi` copy stays until a test-support feature flag is introduced).
  - Lift the no-mixing restriction (CRAM + BAM in one invocation) if a real workload appears. Pre-flight gate is the only place the restriction lives; merge core is already format-agnostic.

#### Segment read fetcher (shared indexed-segment read source)
- **Status:** fixes-applied (2026-06-16; increments #1–#3; standalone primitive, no consumer wired yet)
- **Plan:** [segment_read_fetcher.md](doc/devel/implementation_plans/segment_read_fetcher.md)
- **Impl report:** [segment_read_fetcher_2026-06-16.md](doc/devel/reports/implementations/segment_read_fetcher_2026-06-16.md)
- **Code:** [src/bam/segment_reader.rs](src/bam/segment_reader.rs) — `AlignmentFile` (`BamFile`/`CramFile`), `from_input`, `get_reads_from_segment`, pooled re-seekable readers (`Mutex<Vec<Handle>>`), `MappedReadsInSegment` per-call iterators, narrow `SegmentReadFilter` config, shared `classify_segment_record` per-record cascade. Reuses `classify_pre_decode` / `record_buf_to_mapped_read` / `query_interval` / `ContigInterval::overlaps_record` (visibility lifted to `pub(super)`); two new typed errors (`MissingCramReference`, `InvalidSegment`).
- **Tests:** 21 in-module (overlap + 1-based-inclusive edges, read spanning two segments yielded by both, pool open-once / resting size, poisoned-pool recovery, `par_iter` determinism on a shared `&AlignmentFile`, empty segment, CRAM container early-stop, target-contig + MAPQ filters, both format/index mismatch arms, error paths; BAM + CRAM, `.crai` via `cram::fs::index`). All gates green (fmt / clippy `-D warnings` / doc; 1059 lib + integration tests).
- **Latest review:** [segment_reader_2026-06-16.md](doc/devel/reports/reviews/segment_reader_2026-06-16.md) — **Approve-with-changes**: 0 Blockers, 6 Major, 9 Minor + Nits. Concurrency design verified sound; happy path well tested both formats. Audit trail: `tmp/review_2026-06-16_segment-reader/`.
- **Latest fixes-applied:** [fixes_applied_2026-06-16.md](doc/devel/reports/reviews/fixes_applied_2026-06-16.md) — all 6 Majors **Applied** (M1 Mutex-poison recovery + test; M2 explicit format/index mismatch arms + test; M3 shared `classify_segment_record` chokepoint; M4 CRAM container early-stop + test; M5 narrow `SegmentReadFilter`; M6 report-location convention fixed in the skills on `main` `8c32b00` + report moved) plus Mi1 (handle-invariant docs) and Mi2 (`#[cfg_attr(not(test), allow(dead_code))]`). 1059 lib + integration tests pass; the only `--all-targets` failure is the pre-existing `psp_writer_perf` bench panic.
- **Open (deferred from the 2026-06-16 fix run):**
  - **Minors:** `get_reads_from_segment` `get_`-prefix / diverges from `query` (Mi3, naming — decide at #5); pool-never-shrinks doc note (Mi4); residual inline `Io` sites (Mi5, reader-scoped borrow constraint); cast-spelling + boundary test (Mi6); CRAM container-skip `span==0` test (Mi7, blocked on a multi-container fixture knob); `from_input` path↔header validation (Mi8, when #4 lands); test-only back-ref into `pileup::per_sample::cram_files` (Mi9). Plus the §8 missing tests not yet added (handle-leak-on-index-error, BAM error-fuse, min-length-filters-everything ×2, CRAM error paths, downstream-container completeness, `resolve_segment` boundary) — several blocked on the multi-container fixture limitation. Nits grouped for a cosmetic pass.
  - **#4** — wire the SSR Stage-1 `fetch_locus_reads` onto it (deferred: SSR module not on this branch).
  - **#5** — retrofit the SNP `--regions` path (replace per-BED-interval re-`query`), measuring the `--regions` tax improvement; converge the duplicated chunk/container-walk with `OwnedIndexed{Bam,Cram}Records` then. Remove the module-level dead-code allow when a consumer lands.
  - Deferred (plan §8): CRAM container cache; binary-search the crai head (measure-first, when #4 provides a real workload).

#### Pileup walker
- **Status:** shipped
- **Code:** [src/per_sample_pileup/pileup/](src/per_sample_pileup/pileup/)
- **Subfeature plans & impl reports:**
  - lazy CIGAR — [plan](doc/devel/implementation_plans/pileup_lazy_cigar.md), [impl](doc/devel/reports/implementations/pileup_lazy_cigar_2026-05-07.md)
  - fold cache — [impl](doc/devel/reports/implementations/pileup_fold_cache_2026-05-07.md)
  - freebayes-style bench — [plan](doc/devel/implementation_plans/pileup_freebayes_style_benchmark.md), [impl](doc/devel/reports/implementations/pileup_freebayes_bench_c_2026-05-08.md)
  - pull-iterator walker — [plan](doc/devel/implementation_plans/pileup_pull_iterator.md), [impl](doc/devel/reports/implementations/pileup_pull_iterator_2026-05-14.md)
  - unique chain ids — [plan](doc/devel/implementation_plans/unique_chain_ids.md), [impl](doc/devel/reports/implementations/unique_chain_ids_2026-05-14.md)
- **Latest reviews:** [pileup_2026-05-11.md](doc/devel/reports/reviews/pileup_2026-05-11.md), [perf_pileup_2026-05-12.md](doc/devel/reports/reviews/perf_pileup_2026-05-12.md)
- **Latest fixes-applied:** [fixes_applied_2026-05-13.md](doc/devel/reports/implementations/fixes_applied_2026-05-13.md)
- **Open:** none

#### BAQ
- **Status:** shipped
- **Plan:** [baq.md](doc/devel/implementation_plans/baq.md)
- **Impl report:** [baq_2026-05-12.md](doc/devel/reports/implementations/baq_2026-05-12.md)
- **Latest review / perf review / fixes:** [baq_2026-05-12.md](doc/devel/reports/reviews/baq_2026-05-12.md), [perf_baq_2026-05-12.md](doc/devel/reports/reviews/perf_baq_2026-05-12.md), [fixes_applied_2026-05-12.md](doc/devel/reports/implementations/fixes_applied_2026-05-12.md)
- **Open:** none

#### Indel normalization (left-alignment)
- **Status:** implemented + reviewed + restructured (branch `indel-normalization`, not yet merged to `main`)
- **Spec:** [calling_pipeline_architecture.md §"Indel normalization (left-alignment)"](doc/devel/specs/calling_pipeline_architecture.md) (commit `8d7dd17`)
- **Plan:** [indel_normalization.md](doc/devel/implementation_plans/indel_normalization.md)
- **Final approach (revised — option 2):** the code review surfaced that the
  BAM/CRAM input cascade **already** left-aligned every read's indels (the
  always-on "F3" pass `left_align_indels` in
  [src/bam/alignment_input.rs](src/bam/alignment_input.rs), pre-existing on
  `main`). Rather than add a second normalization in the BAQ prep stage, the
  feature now **replaces F3's single-forward-pass shifter with the GATK
  `leftAlignIndels` port** at F3's call site: F3's `left_align_indels` /
  `try_apply_indel_shift` (+ helpers, + its `f3_*` tests) are deleted, and
  the reader calls [`indel_norm::left_align_indels`](src/pileup/walker/indel_norm.rs)
  (GATK port: right-to-left, trim-first, dual ref+read check, **collision-merge**
  — the refinement F3 explicitly punted on). The reader's existing reference
  fetch (for the F1 mismatch filter) is reused; nothing else in the merged
  reader changes. Because the reader is upstream of everything, normalization
  is structurally mandatory — including under `--no-baq` (no BaqEngine
  wiring; the earlier prep-stage detour was reverted).
- **Code:** [src/pileup/walker/indel_norm.rs](src/pileup/walker/indel_norm.rs) (the GATK port + tests); call site at [src/bam/alignment_input.rs](src/bam/alignment_input.rs) ("F3").
- **Latest review:** [indel_normalization_f3_replacement_2026-05-29.md](doc/devel/reports/reviews/indel_normalization_f3_replacement_2026-05-29.md) — Approve-with-changes, of the F3-replacement diff. Confirmed the call-site integration is correct (`ref_seq[0]` is the read's first aligned base → `read_start=0`; reference fetch reused from the F1 filter; F3 fully deleted, port a faithful drop-in). Fixed during the review: a stale F3 comment falsely claiming "adjacent indels invariant" (the port merges them), a broken `[left_align_prepared]` doc link (a `cargo doc` regression this PR introduced), the missing `// UNREACHABLE:` comment, and the untested `left_align_indels` wrapper/debug-invariant (added `wrapper_*` tests). The [first review](doc/devel/reports/reviews/indel_normalization_2026-05-29.md) (BaqEngine detour) is superseded.
- **Open (from the F3-replacement review):**
  - **M2** — add a reader-level integration test (`AlignmentMergedReader` → `fetch_ref_for_read` → `left_align_indels`); the deleted `f3_*` unit tests have no integration equivalent, so a slice-offset/`read_start` wiring regression would be invisible.
  - **M4** — hot-path allocation regression: F3 shifted in place (zero-alloc); the port allocates 4–5 `Vec`s per indel read on the single-threaded reader loop. Add scratch reuse + a normalization criterion bench (no bench covers this path).
  - **M5** — move `indel_norm` to `src/bam/` (beside its sole consumer; removes a `bam → pileup` edge) — or a neutral peer if a second consumer appears.
  - **M6** (needs verification) — the port can merge colliding indels into an adjacent `D`/`I` pair that G2 `cigar_is_bad` (run earlier) doesn't re-check; confirm reachability + walker handling with a test.
  - **Minors:** read-consume guard emits no telemetry (Mi1); demote `left_align_cigar`/`LeftAlignResult` to private (Mi2); `indel_norm/tests.rs` dir convention (Mi3); plus carryovers `Range`→`IndexRange`, predicate-helper factoring, `build_cigar` splice→fold.
  - Deferred measurement: HG002/tomato indel-recall delta (the original acceptance signal).

#### Pileup → psp seam
- **Status:** shipped
- **Code:** [src/per_sample_pileup/pileup_to_psp.rs](src/per_sample_pileup/pileup_to_psp.rs)
- **Impl report:** [pileup_to_psp_seam_2026-05-14.md](doc/devel/reports/implementations/pileup_to_psp_seam_2026-05-14.md)
- **Open:** none

#### `pop_var_caller` CLI
- **Status:** Stage 1 CLI shipped (subcommands `pileup`, `psp-to-pileup`);
  cohort CLI **shipped** — every Deferred finding from the
  2026-05-19 review is Applied (Waves 1 – 5).
- **Plans:**
  - Stage 1 CLI (`pileup`, `psp-to-pileup`):
    [pop_var_caller_pileup_cli.md](doc/devel/implementation_plans/pop_var_caller_pileup_cli.md)
  - Cohort CLI (`var-calling`, `estimate-contamination`,
    `var-calling-from-bam`):
    [pop_var_caller_cohort_cli.md](doc/devel/implementation_plans/pop_var_caller_cohort_cli.md)
- **Impl reports:**
  - Cohort slice:
    [pop_var_caller_cohort_cli_2026-05-19.md](doc/devel/reports/implementations/pop_var_caller_cohort_cli_2026-05-19.md)
  - `var-calling-from-bam` per-chromosome parallelism (2026-05-24):
    [var_calling_from_bam_per_chromosome_2026-05-24.md](doc/devel/reports/implementations/var_calling_from_bam_per_chromosome_2026-05-24.md)
- **Latest reviews (cohort slice):**
  [ref_fetcher_2026-05-23.md](doc/devel/reports/reviews/ref_fetcher_2026-05-23.md)
  (code review of `src/per_sample_pileup/ref_fetcher.rs` — status
  *Request-changes*; 3 Blockers, 23 Major, 19 Minor. Production
  code paths are correct; gaps are in the public-trait contract
  surface. Highlights: B1 `.fai` `line_bases=0` panic on
  attacker-influenced input; B2 `legacy_io_error` drops
  `contig_name` from every legacy-walker diagnostic; B3 the
  `OutOfPattern` contract has no production-buffer-size regression
  test; convergent Majors: stale doc link at line 282 — 6
  categories, missing `#[non_exhaustive]` on `ChromRefFetchError`,
  `Source::Memory` cfg-gating, vestigial legacy surface on
  `StreamingChromRefFetcher` after the Step-2 migration, ~80 LOC
  constructor duplication across three sites.);
  [perf_ref_fetcher_2026-05-23.md](doc/devel/reports/reviews/perf_ref_fetcher_2026-05-23.md)
  (FASTA fetcher post-Step-2-migration — verdict: *Apply the listed
  wins*, gated on H4; 4 Hot-path, 10 Likely, 4 Speculative;
  headline: `ChromRefBaseIter::next` 10.96% wall with 4.02% on
  `MutexGuard<StreamState>` drop alone — H1 drops the `Sync`
  requirement, swaps `Mutex<StreamState>` for `RefCell<StreamState>`);
  [perf_psp_to_vcf_2026-05-20.md](doc/devel/reports/reviews/perf_psp_to_vcf_2026-05-20.md)
  (end-to-end perf review against real tomato data — verdict: *Apply
  the listed wins*; 7 Hot-path findings, 12 Likely, 5 Speculative;
  headline: pipeline is essentially single-threaded, DUST 33 %,
  allocations 21 %; H1 per-chromosome parallelism is the
  order-of-magnitude lever);
  [cohort_cli_2026-05-19.md](doc/devel/reports/reviews/cohort_cli_2026-05-19.md)
  (correctness — Request-changes: 0 Blockers, 14 Major (M1–M14), 23
  Minor + grouped Nits).
- **Latest fixes-applied (cohort slice):**
  Wave 5 of the deferred follow-up,
  [cohort_cli_2026-05-19_applied_wave5.md](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave5.md) —
  **Mi20**, **Mi23**, **M1/M2 follow-up**, **M5 follow-up**
  (FASTA → `.psp` MD5 enforcement is the only behaviour change
  across the 5-wave plan). 890 lib + 9 cohort integration
  tests pass. Earlier passes:
  [Wave 4](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave4.md)
  (M13);
  [Wave 3](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave3.md)
  (Mi8 / Mi19 / M9-followup / Mi18 / M10 / M11);
  [Wave 2](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave2.md)
  (M4 / Mi2 / Mi14 / Mi21);
  [Wave 1](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied_wave1.md)
  (M8 / Mi5 / Mi6 / Mi13);
  [original pass](doc/devel/reports/reviews/cohort_cli_2026-05-19_applied.md)
  (20 Applied + Mi7 Disputed-with-test).
  **Status: shipped** — every originally-Deferred finding is
  Applied; the 5-wave follow-up plan is complete.
- **Code:** [src/pop_var_caller/](src/pop_var_caller/)
- **Integration tests:**
  [tests/pileup_cli_integration.rs](tests/pileup_cli_integration.rs)
  (Stage 1) and
  [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs)
  (cohort subcommands).
- **Open (from 2026-05-23 ref-fetcher perf review):**
  - **H4 (gate)** — Repair the four broken benches before any other
    finding can be measured: `cohort_e2e_perf`, `var_calling_perf`,
    `baq_perf`, `pileup_walker_scaling`. API drift from the Step-2
    `ChromRefFetcher` migration (stale `SyncRefFetcher` import,
    `CohortPipelineParams.fetcher`/`.chromosomes` field drift,
    `drive_cohort_pipeline` arg count). Mechanical port; no new deps.
  - **H1** — `src/per_sample_pileup/ref_fetcher.rs:117` —
    `Mutex<StreamState>` → `RefCell<StreamState>` and relax
    `SharedRefFetcher = Arc<dyn ChromRefFetcher + Send + Sync>` to
    `Arc<dyn ChromRefFetcher + Send>`. Per-base mutex unlock is
    4.02 % of cohort wall time; the `Mutex` is uncontended by
    construction (per-worker ownership documented at L113-117).
  - **H2** — `iter_bases` returns `Box<dyn Iterator>`; GAT-monomorphic
    `type BaseIter<'a>` (or a `ChromRefFetcherTyped` subtrait that
    keeps the dyn-safe surface intact) removes per-byte vtable
    dispatch on the DUST mask scan.
  - **H3** — `StreamingChromRefFetcher::fetch` returns owned `Vec<u8>`
    where every PerGroupMerger consumer only borrows. Switch to
    `&[u8]` — entangled with H1's `RefCell` (lands together).
  - **L1 / L10** — Wrap `Source::File(File)` in `BufReader<File>` in
    both `StreamingChromRefFetcher` and `ManualEvictChromRefFetcher`.
  - **L2** — Autovec uppercase pass in `refill` /
    `read_uppercased_bases` — strip-newlines via `filter+take+extend`
    then `make_ascii_uppercase` on the destination slice.
  - **L3** — `src/var_calling/dust_filter.rs:728` —
    `io::Error::other(format!("{e}"))` → `io::Error::other(e)`
    (drops the `String` alloc, preserves `source()`).
  - **L5** — `ChromRefBaseIter::Drop` re-locks the mutex; folds
    into H1 for free.
  - **L6** — Add `[profile.release-with-debug]` to `Cargo.toml`
    that `inherits = "release"` + `debug = true`, so future
    profile captures reproduce the inlined-frame call-graph from
    a clean checkout.
  - **L7** — Wire `mimalloc::MiMalloc` as `#[global_allocator]`
    behind `alloc-mimalloc` in `src/main.rs`. Glibc allocator
    symbols sum to ~14.5 % cpu_atom / ~23.2 % cpu_core; gate
    merge on a paired A/B against the 13-thread server.
  - **L8** — Add `benches/ref_fetcher_perf.rs` with four
    `criterion` functions (`streaming_iter_bases_full_contig`,
    `streaming_fetch_per_group_window`, `manual_evict_fetch_then_evict`,
    `refill_warm_cache`) so H1/H2/H3 effects can be localised below
    the cohort_e2e noise floor.
  - **L4 / L9 / S1–S4** — see the full report.
- **Open (from 2026-05-20 perf review):**
  - **H1 — closed 2026-05-20** (per-chromosome parallelism shipped).
    Five-commit PR `309a5be` → `0b1e958`; impl report
    [cohort_per_chromosome_parallel_2026-05-20.md](doc/devel/reports/implementations/cohort_per_chromosome_parallel_2026-05-20.md).
    Realised **3.85× wall-time speedup at T=13** on the multi-chrom
    real-data fixture (106.6 s → 27.7 s). Workload imbalance
    (ch00 carries 1.4 M of 2.6 M total reads as the unplaced/decoy
    contig) gates the ceiling below the plan's predicted 6–10×
    range; L5 contention absorbs another ~23 %.
  - **L1 — closed 2026-05-20** (per-group merger inner `par_iter`
    removed, commit `309a5be`). Mandatory prep for H1 — nested
    rayon under per-chrom outer is wasteful and would have
    polluted the H1 measurement.
  - **H2 / H3** — DUST `find_perfect` inner-loop: replace
    `Vec::insert(j, …)` (O(n) memmove) and `VecDeque` indexing
    (per-element wrap + bounds check). Gated by **H7**.
  - **H4** — PSP reader: replace per-block `SeekFrom::Start` (which
    discards the 64 KiB BufReader) with `seek_relative` for the
    post-header rewind and "skip if already at offset" for the
    pre-block seek (or pull-style `fill_buf`/`consume`).
  - **H5** — M5 verify: stream the FASTA in 64 KiB windows
    feeding `Md5::update` instead of materialising whole contigs (up
    to 91 MB) and uppercasing byte-by-byte.
  - **H6** — `PerPositionMerger::next`'s `vec![None; n_samples]` per
    emit (~19M allocations / run at N=10): lending-iterator pivot,
    paired with the same fix in `DustFilter`.
  - **H7** — add an isolated `var_calling_dust_filter` criterion bench
    (gates H2 / H3 / L1 / L2 / L12 measurements).
  - **L2–L12** — see the full report; highlights: **L5 is now the
    next ceiling under H1** (`SyncRefFetcher` `RwLock<HashMap>::read()`
    on every fetch — pre-warm into `Vec<Arc<Vec<u8>>>` indexed by
    `chrom_id`, drop the noodles `Repository` runtime dep);
    L2 `AlleleObservation` `Vec` → `SmallVec`,
    L3 PSP CSR decoder, L4 `fetch_from_repository`
    `make_ascii_uppercase`, L8 `[profile.release] debug = "line-tables-only"`
    → `debug = true`, L9 `alloc-mimalloc` A/B against real-data
    workload, L10 missing drained-count assertions in benches.
  - **Per-chrom parallel follow-ups** (from the impl report): streaming
    concat (v2 — append finished fragments while slower chroms run);
    block-level bgzf concat (v2 — skip decompress+re-encode);
    sub-chromosome decomposition (push below
    `max(per-chrom-time)` on imbalanced workloads like tomato's ch00);
    `RLIMIT_NOFILE` bump (defaults to 1024; N=256 × n_chrom > 1024 fds);
    `var-calling-from-bam` + `estimate-contamination` parallelisation;
    `profile_cohort_e2e --em-convergence-threshold` knob; posterior
    `DidNotConverge` emit-with-flag long-term fix; multi-chrom
    integration-test fixture in `tests/common/mod.rs`.
- **Open (from cohort-slice review):**
  *None — the 16 originally-Deferred findings are all Applied.*
  - **Closed in Wave 5 (2026-05-19):** **Mi20**, **Mi23**,
    **M1+M2-followup**, **M5-followup**.
  - **Closed in Wave 4 (2026-05-19):** **M13**.
  - **Closed in Wave 3 (2026-05-19):** **Mi8**, **Mi19**,
    **M9-followup**, **Mi18**, **M10**, **M11**.
  - **Closed in Wave 2 (2026-05-19):** **M4**, **Mi2**, **Mi14**,
    **Mi21**.
  - **Closed in Wave 1 (2026-05-19):** **M8**, **Mi5**, **Mi6**,
    **Mi13**.
  - Selected deferred Nits — drop unused `#[from]` variants, vestigial
    `let _ = cfg;`, `Stage1RunSummary` renames, `DEFAULT_BATCH_ID`
    value-mismatch, etc.
  - `bcftools view` / `bcftools stats` manual smoke against real
    cohort data (pre-existing).

---

### Stage 2 — per-sample file (`.psp`) contract

Stage 2 is the on-disk artefact format that Stage 1 writes and Stages 3–6
consume. Not a runtime step — an interface.

#### `.psp` writer + reader
- **Status:** shipped
- **Spec:** [per_sample_pileup_format.md](doc/devel/specs/per_sample_pileup_format.md)
- **Plans:**
  - reader: [psp_reader.md](doc/devel/implementation_plans/psp_reader.md)
  - writer/reader bootstrap: [per_sample_pileup_writer_reader.md](doc/devel/implementation_plans/per_sample_pileup_writer_reader.md)
  - psp → pileup roundtrip: [psp_to_pileup.md](doc/devel/implementation_plans/psp_to_pileup.md)
- **Impl reports:** [psp_reader_2026-05-13.md](doc/devel/reports/implementations/psp_reader_2026-05-13.md), [psp_to_pileup_2026-05-15.md](doc/devel/reports/implementations/psp_to_pileup_2026-05-15.md)
- **Latest reviews:** [perf_psp_reader_2026-05-23.md](doc/devel/reports/reviews/perf_psp_reader_2026-05-23.md) (cohort-hot-path re-review against the post-H3 profile — verdict: *Apply the listed wins*; 3 Hot-path / 10 Likely / 4 Speculative; H1 CSR ragged-column collapse subsumes 87% of project-side allocator pressure, H2 = `SeekFrom::Current` for the 2026-05-20 H4 that was never applied, H3 = per-allele bounds-check hoist), [psp_2026-05-13.md](doc/devel/reports/reviews/psp_2026-05-13.md), [psp_reader_2026-05-13.md](doc/devel/reports/reviews/psp_reader_2026-05-13.md), [perf_psp_writer_2026-05-13.md](doc/devel/reports/reviews/perf_psp_writer_2026-05-13.md), [perf_psp_reader_2026-05-13.md](doc/devel/reports/reviews/perf_psp_reader_2026-05-13.md)
- **Latest fixes-applied:** [fixes_applied_psp_reader_2026-05-13.md](doc/devel/reports/implementations/fixes_applied_psp_reader_2026-05-13.md), [perf_psp_reader_2026-05-13_applied.md](doc/devel/reports/reviews/perf_psp_reader_2026-05-13_applied.md), [perf_psp_writer_2026-05-13_applied.md](doc/devel/reports/reviews/perf_psp_writer_2026-05-13_applied.md)
- **Open (from 2026-05-23 PSP reader perf review):**
  - **H1** — CSR collapse of `DecodedBlock`'s
    `allele_seqs: Vec<Vec<u8>>` and
    `allele_chain_ids: Vec<Vec<ChainId>>` into
    `(data: Vec<u8>, offsets: Vec<u32>)` + per-`RecordsIter`
    scratch reuse. Symmetric with the writer's existing
    `encode_list_column_csr`. Stage A keeps `AlleleObservation.seq`
    as `Vec<u8>` (one alloc per emit); Stage B (cross-cutting)
    switches it to `Box<[u8]>` or arena handle. Gate Stage B
    behind dhat data via the missing `examples/dhat_psp_reader.rs`
    (L7).
  - **H2** — `SeekFrom::Current(delta)` helper for the per-block
    seek at `reader.rs:587` and the post-block-header rewind at
    `reader.rs:807`. `SeekFrom::Start` discards `BufReader`'s
    64 KiB user-space buffer on every block transition even
    though blocks are written contiguously on disk; the 2026-05-20
    review's H4. This was never applied.
  - **H3 (gated by H1)** — leading `assert!(allele_end <=
    block.<col>.len())` for the 9 ragged columns at the top of
    `materialise_next_record`. Collapses the per-allele inner
    loop's bounds checks to one each outside the loop.
  - **L1** — `decode_list_column_pod<T: Pod>` LE-slab cast for
    `ChainId = u64` (mirror of the writer's
    `encode_list_column_csr`'s Pod fast path).
  - **L2** — drop `stream_position()` from `read_block_header` by
    threading `entry.block_offset` from the caller (saves one
    `lseek(2)` per block).
  - **L3** — `region_records` first-block seek (subsumed by H2).
  - **L4 (gated by H1)** — pack the 7 fixed-width per-allele
    columns into a single `#[repr(C)]` row for the
    `materialise_next_record` gather.
  - **L5** — `AlleleObservation` SmallVec (cross-cuts the
    2026-05-20 L2; coordinate with that finding).
  - **L6** — `benches/psp_reader_perf.rs` is `Cursor`-backed;
    add a file-backed `BufReader<File>` group so H2 / L1 / L2 /
    L3 have a microbench signal.
  - **L7** — `examples/dhat_psp_reader.rs` (deferred since
    2026-05-13 L9; the post-H3 profile justifies it as a
    prerequisite for H1's Stage B decision).
  - **L8** — bench bodies should `assert_eq!(count, expected)`
    so a silent-truncate refactor can't pass.
  - **S1** — sweep `DEFAULT_BUFFERED_IO_CAPACITY` (64 KiB
    today) after H2 lands.
  - Background-confirmed: L5 (varint fast/cold), L6 (LE-slab cast),
    L8 (BufReader doc) from 2026-05-13 are all still in place.

---

### Stage 3 — low-complexity (sdust) filter

Streaming per-position filter that computes sdust scores from the
reference and silently drops low-complexity records. No intermediate
mask file. Algorithm ported from `lh3/sdust` (vendored at `sdust/`,
gitignored).

- **Status:** fixes-applied
- **Spec section:** `## Stage 3 — low-complexity filter` in [calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md)
- **Plan:** [dust_filter.md](doc/devel/implementation_plans/dust_filter.md)
- **Code:** [src/var_calling/dust_filter.rs](src/var_calling/dust_filter.rs)
- **Tests:** 38 tests in the module (algorithmic core, config
  validation, iterator plumbing, golden vector against committed
  `lh3/sdust` outputs, seeded-random invariant sweep over ~40 K
  input combinations, half-open boundary, threshold-strictness
  boundary, `pos == 0` latch, exhaustion latch). `lh3/sdust` is
  not a build- or test-time dependency.
- **Impl report:** [dust_filter_2026-05-17.md](doc/devel/reports/implementations/dust_filter_2026-05-17.md)
- **Latest review:** [dust_filter_2026-05-17.md](doc/devel/reports/reviews/dust_filter_2026-05-17.md) — Request-changes: 0 Blockers, 8 Major, 24 Minor + 9 Nits.
- **Latest fixes-applied:** [dust_filter_2026-05-17_applied.md](doc/devel/reports/reviews/dust_filter_2026-05-17_applied.md) — 36 of 41 findings Applied (7 of 8 Majors + 22 of 24 Minors + 7 Nits); 1 Applied-with-adaptation (Mi21 partially blocked on parallel `posterior_engine` WIP); 2 Won't-fix per project preference (M4 out-of-order check is the merger's responsibility; Mi6 informational tracing conflicts with the "no logs" project preference); 2 Nits won't-fix.
- **Open:**
  - **Mi21** — full `cargo doc --no-deps --all-features` with `-D warnings` to be re-run once the parallel `posterior_engine` refactor compiles; intra-doc links in `dust_filter.rs` itself have been visually verified.
  - Criterion bench — deferred until the cohort CLI runs against
    real cohort data and we can measure DUST's share of real
    wall time.
  - **Closed 2026-05-19** (cohort CLI slice): CLI parser bindings
    (`--complexity-window`, `--complexity-threshold`,
    `--no-complexity-filter`) shipped with `var-calling`;
    end-to-end integration test exists in
    [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs).

---

### Stage 4 — grouping

Walks the filtered per-position stream and bundles overlapping / adjacent
candidate positions into `OverlappingVarGroup`s. Sequential, but emits
independent groups that downstream stages can process in parallel.

#### Multi-way per-position iterator
The k-way merge that produces the per-position stream consumed by the
grouper.
- **Status:** shipped
- **Plan:** [multi_way_per_position_iterator.md](doc/devel/implementation_plans/multi_way_per_position_iterator.md)
- **Impl report:** [multi_way_per_position_iterator_2026-05-15.md](doc/devel/reports/implementations/multi_way_per_position_iterator_2026-05-15.md)
- **Code:** [src/var_calling/per_position_merger.rs](src/var_calling/per_position_merger.rs)
- **Bench:** `var_calling_merger/*` in [benches/var_calling_perf.rs](benches/var_calling_perf.rs)
- **Latest reviews:** [per-position-merger_2026-05-15.md](doc/devel/reports/reviews/per-position-merger_2026-05-15.md), [perf_var_calling_2026-05-16.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16.md)
- **Latest fixes-applied:** bundled into [fixes_applied_2026-05-16.md](doc/devel/reports/implementations/fixes_applied_2026-05-16.md) (cohort run) + [fixes_applied_2026-05-16_v2.md](doc/devel/reports/implementations/fixes_applied_2026-05-16_v2.md) (Mi14/Mi15 OutOfOrder+ChromosomeMismatch land here).
- **Open:**
  - **H6** (perf review): the Stage 4 merger bench is fixture-dominated
    (67 % inclusive in `iter_batched` setup vs 21 % in
    `PerPositionMerger::next`); rebuild fixtures once outside `b.iter`
    before any Stage 4 code-level change can be defended.
  - **L8** (perf review, depends on H6): `vec![None; self.n_samples()]`
    per emission at [per_position_merger.rs:306](src/var_calling/per_position_merger.rs#L306)
    is the candidate Stage 4 allocation site to revisit once the bench
    is unblocked.

#### Variant grouper (Stage 4 bundler)
- **Status:** fixes-applied (2026-05-16); no separate implementation report
- **Plan:** [cohort_variant_grouping.md](doc/devel/implementation_plans/cohort_variant_grouping.md)
- **Code:** [src/var_calling/variant_grouping.rs](src/var_calling/variant_grouping.rs)
- **Bench:** `var_calling_grouper/*` in [benches/var_calling_perf.rs](benches/var_calling_perf.rs)
- **Latest reviews:** [cohort_2026-05-16.md](doc/devel/reports/reviews/cohort_2026-05-16.md), [perf_var_calling_2026-05-16.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16.md)
- **Latest fixes-applied:** [fixes_applied_2026-05-16.md](doc/devel/reports/implementations/fixes_applied_2026-05-16.md) + [fixes_applied_2026-05-16_v2.md](doc/devel/reports/implementations/fixes_applied_2026-05-16_v2.md) (v2 closes the entire v1 deferred list)
- **Open:**
  - Implementation report not yet written.
  - **H6** (perf review): the Stage 4 grouper bench is also
    fixture-dominated (77 % inclusive in `iter_batched` setup vs 13 %
    in `VariantGrouper::next`); rebuild fixtures once outside `b.iter`
    before any Stage 4 code-level change can be defended. Shared
    bench-shape fix with the merger above.

---

### Stage 5 — per-group processing

Allele unification + per-sample likelihood reconstruction + phase-chain
consistency on one `OverlappingVarGroup` at a time. Parallel across groups
via rayon.

#### Per-group merger (allele unification + likelihood)
- **Status:** fixes-applied (2026-05-16)
- **Plan:** [cohort_per_group_merger.md](doc/devel/implementation_plans/cohort_per_group_merger.md)
- **Code:** [src/var_calling/per_group_merger.rs](src/var_calling/per_group_merger.rs)
- **Bench:** `var_calling_per_group_merger/*` in [benches/var_calling_perf.rs](benches/var_calling_perf.rs).
  Slowest stage per element on the criterion baseline (~22 K groups/s
  biallelic, ~12 K groups/s compound) and the primary focus of the
  2026-05-16 perf review.
- **Impl report:** not yet saved.
- **Latest reviews:** [cohort_2026-05-16.md](doc/devel/reports/reviews/cohort_2026-05-16.md), [perf_var_calling_2026-05-16.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16.md)
- **Latest fixes-applied:** [fixes_applied_2026-05-16.md](doc/devel/reports/implementations/fixes_applied_2026-05-16.md) + [fixes_applied_2026-05-16_v2.md](doc/devel/reports/implementations/fixes_applied_2026-05-16_v2.md) (v2 closes the v1 deferred list) + [perf_var_calling_2026-05-16_applied.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16_applied.md) (Wave 1 of the perf review: L2 + H1 + H4 + H3 + H2 applied — H2 phase 2b is the headline ~14 % biallelic win, p = 0.00, confirmed by two back-to-back runs; all 571 tests pass)
- **Open:**
  - Implementation report for the Stage 5 merger has not been saved; the
    next `rust-feature-implementation` run for this feature should produce
    one and link it here.
  - Phase-chain integration tests for the likelihood calculation (see
    *Standing items* below).
  - Re-bench Wave 1 on a quieter host with a clean pre-perf-review
    checkout baseline to put an absolute number on the cumulative
    L2 + H1 + H4 + H3 + H2 effect. The in-session ~12–14 %
    biallelic_64 figure is reliable for H2 (back-to-back confirms);
    the absolute pre-vs-post number drifted with host load. See
    [perf_var_calling_2026-05-16_applied.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16_applied.md).
  - Remaining Hot-path findings to apply (see
    [perf_var_calling_2026-05-16.md](doc/devel/reports/reviews/perf_var_calling_2026-05-16.md)
    §5):
    - **H5** — sweep `DEFAULT_BATCH_SIZE` (currently a self-declared
      placeholder at [per_group_merger.rs:55](src/var_calling/per_group_merger.rs#L55))
      and reset to the measured optimum.
    - **H6** — Stage 4 bench fixture-rebuild fix (the per_position_merger
      / variant_grouper benches' iter_batched setup is fixture-dominated
      per the original perf review). Unblocks Stage 4 code-level findings.
    - **H7** — cohort-size sweep at N=10/64/256/1024 samples.
  - Wave 2 / 3 / Likely / Speculative findings tracked in the report.

#### Within-chromosome chunk-parallel rewrite (now flat under `src/var_calling/`; formerly `from_psp/`, originally `cohort_block/`)
- **Status:** **superseded (2026-06-05)** — the columnar
  `driver`/`worker`/`loader`/`columns`/`partition`/`two_pass`/`kernels`
  chain this block describes was deleted and replaced by the
  record-streaming pipeline; see the "Re-architected record-streaming
  pipeline" block below. (Historical status follows.) fixes-applied (2026-06-01 review — see "Applied" note below); reviewed (2026-06-01, full-but-prioritized subtree review); parallel block-consume shipped
  (`0d49cf8`, `51b5c63`); **streaming-columnar produce rewrite —
  Stages 1, 2 & 3 implemented**, plus the **DUST worker pool**
  (parallel DUST-ahead). Memory fix landed (Stage 2: N=26 peak RSS
  3963→550 MB); serial DUST moved off the critical path (Stage 3) then
  parallelised across the independent covered intervals (DUST pool). The
  branch now **beats `main` on both wall and RSS** at N=8 (4.1 s / 79 MB
  vs 6.0 s / 142 MB) and N=26 (7.8 s / 291 MB vs 8.9 s / 403 MB), 8
  threads. Optional follow-up: sub-span DUST jobs if a single giant
  interval ever dominates the floor.
- **Plans:**
  - Streaming-columnar produce (current): [cohort_produce_streaming_columnar.md](doc/devel/implementation_plans/cohort_produce_streaming_columnar.md)
  - DUST worker pool (implemented; parallel DUST-ahead — the measured wall floor after Stage 3): [cohort_dust_worker_pool.md](doc/devel/implementation_plans/cohort_dust_worker_pool.md)
  - Master: [cohort_within_chromosome_parallel.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel.md)
  - Phase A.2 column-native EM: [cohort_within_chromosome_parallel_phase_a2_em.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_a2_em.md)
  - Phase B prereq (variant-bounded chunks): [cohort_within_chromosome_parallel_phase_b1_variant_bounded_chunks.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b1_variant_bounded_chunks.md)
  - Phase B prereq (estimate-contamination migration): [cohort_within_chromosome_parallel_phase_b2_estimate_contamination.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b2_estimate_contamination.md)
  - Phase B parallel windows: [cohort_within_chromosome_parallel_phase_b_parallel_windows.md](doc/devel/implementation_plans/cohort_within_chromosome_parallel_phase_b_parallel_windows.md)
- **Impl reports:**
  - Phase A (2026-05-28): [cohort_within_chromosome_parallel_phase_a_2026-05-28.md](doc/devel/reports/implementations/cohort_within_chromosome_parallel_phase_a_2026-05-28.md)
  - Streaming produce Stage 1 — span-addressable columnar PSP reader (2026-05-31): [cohort_produce_streaming_columnar_stage1_2026-05-31.md](ia/reports/implementations/cohort_produce_streaming_columnar_stage1_2026-05-31.md)
  - Streaming produce Stage 2 — streaming fold+compact producer / memory fix (2026-05-31): [cohort_produce_streaming_columnar_stage2_2026-05-31.md](ia/reports/implementations/cohort_produce_streaming_columnar_stage2_2026-05-31.md)
  - Streaming produce Stage 3 — DUST-ahead queue / serial DUST off the critical path (2026-05-31): [cohort_produce_streaming_columnar_stage3_2026-05-31.md](ia/reports/implementations/cohort_produce_streaming_columnar_stage3_2026-05-31.md)
  - DUST worker pool — parallel DUST-ahead (2026-05-31): [cohort_dust_worker_pool_2026-05-31.md](ia/reports/implementations/cohort_dust_worker_pool_2026-05-31.md)
- **Code:** [src/var_calling/](src/var_calling/) — the chunk driver now lives flat in the `var_calling` module root: `driver.rs`, `loader.rs`, `worker.rs`, `columns.rs`, `column_span_reader.rs`, `partition.rs`, `test_helpers.rs`, plus `kernels/{mod, unify_alleles, project_scalars, compute_log_likelihoods}.rs`. (Originally `cohort_block/`, then the `from_psp/` submodule; the `from_psp/` layer was merged up on 2026-06-01 once the sibling `from_bam/` path was gone. The streaming-columnar rewrite earlier folded `pre_pass.rs` into the loader/partition path.)
- **Tests:** 88 unit tests in the module (per `cargo test --lib var_calling::cohort_block` at commit `36989d6`); 3 integration tests in [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs) (`var_calling_emits_deterministic_vcf_across_runs`, `var_calling_byte_identical_across_worker_windows_per_chunk`, `var_calling_byte_identical_across_target_variants_per_chunk`).
- **Latest fixes-applied:** [cohort_block_2026-05-29_applied.md](doc/devel/reports/reviews/cohort_block_2026-05-29_applied.md) — **Wave 1**: all 5 Blockers Applied (B1 / B2 / B3 / B4 / B5-deferred-to-Wave-2 per Q2) + M5 (bundled with B1) + M11 (bench/example unblocker) + M14 (carryover snapshot helper) + M17 (drop trailing `..`) + M18 (delete dead `chain_id_scratch`) + M19 (`debug_assert!` on sorted `masked_intervals`). 1 026/1 026 lib pass (+3); 21/21 cohort_cli integration pass; fmt clean; criterion baseline saved. 47 findings deferred to Waves 2–3 per Q4 ("apply all structural refactors now"). 2 Won't fix per Q1 (Mi8 / Mi13). Out-of-scope edits flagged in §12.
- **Latest review:** [var_calling_2026-06-01.md](doc/devel/reports/reviews/var_calling_2026-06-01.md) — **Request-changes** (full-but-prioritized review of the whole `src/var_calling/` subtree, 27 489 LoC, commit `3c9ebf2`; prioritized the post-reorg `from_psp/` integration seams, spot-checked the mid-May-reviewed stage internals). 0 Blockers, 13 Major, 10 Minor + grouped Nits. No correctness Blockers — `unsafe_concurrency` verified the parallel-worker soundness holds after the reorg (the `Send + !Sync` REF fetcher never crosses a thread boundary; `seq_idx`-ordered emit + additive stats back byte-identity). Drivers of the verdict: the `cargo doc` CI gate is **red** (M1 — two redundant explicit intra-doc links, `column_span_reader.rs:3` + `posterior_engine.rs:925`); the rewrite's hard byte-identity contract is **unguarded by any test** (M2 — serial-vs-parallel; M9 — no committed golden VCF); and the `from_bam` removal left dead code + dead doc links (M3 dead `ChunkDriverError::PerGroupMerger` variant+`From`; M4 dead `into_shared_ref_fetcher`; M7 three dangling `from_bam` intra-doc links that render as silent dead URLs). Other Major: M5 inert/misleading `chunk_genomic_span` knob (+ phantom `MAX_CHUNK_SPAN_GROWTH` doc ref); M6 startup log omits every per-stage default; M8 `from_psp` public surface should be `pub(crate)` + its facade is bypassed; M10 no perf bench/threshold on the chunk-driver hot path; M11 `BlockQueue` hand-rolled concurrency untested; M12 `DustAheadPool::shutdown` swallows worker panics; M13 unbounded `mpsc` result/recycle channels. Four open questions gate M5/M8/M2/M9 (chunk_genomic_span fate; from_psp public-vs-internal; byte-identity oracle now that the streaming driver is deleted; keep both serial+parallel drive paths?). Per-category audit trail at `tmp/review_2026-06-01_var_calling/`.
- **Latest perf review:** [perf_var_calling_2026-06-01.md](doc/devel/reports/reviews/perf_var_calling_2026-06-01.md) — **Apply the listed wins.** Closes review-finding **M10**: two committed criterion benches added to `benches/var_calling_perf.rs` (`var_calling_chunk_loader` = produce, `var_calling_run_window` = consume kernels + EM), swept N ∈ {8,64,200}. Evidence: bench (loader ~linear in N, **`run_window` ~N^1.7**), CPU sampling profile (macOS `sample`/`samply`), DHAT (real cohort N=18: 15.5 GB / 51.4 M blocks), `cargo-show-asm`. **3 Hot-path:** **H1** the QUAL allele-count convolution `compute_qual_via_exact_af` (posterior_engine.rs:3037-3047) does a *scalar pairwise* `log_sum_exp_2` in an O(n_samples²·ploidy²) loop = ~51 % of `run_window` self-time and the mechanism behind the N^1.7 scaling (the hot `log_sum_exp_2` is the QUAL convolution, **not** the EM, which takes the homogeneous-fixation fast path under default config); **H2/H3** the two `unify_alleles` allocation sites (`build_chain_proposals_columnar` BTreeMap value-Vecs ~15.7 M blocks — M3's scratch-hoist is defeated by `BTreeMap::clear` dropping the Vecs; `project_per_position_into_scratch` byte-key clones ~7.3 M blocks). 11 Likely (mimalloc on the production binary; per-group genotype-table cache; thread oversubscription producer-rayon vs worker-threads; DUST `to_vec` copies; `for_contig` re-parsing `.fai` per interval; AoS→SoA on the LH-hot scalars; bench hardening; …), 4 Speculative, Notes. Build config already tuned (`lto=fat`, `codegen-units=1`, `panic=abort`, pinned `target-cpu`); production allocator is the one build lever. Per-category audit trail + profiles at `tmp/perf_review_2026-06-01_var_calling/`.
- **Latest fixes-applied (2026-06-01, branch `var-calling-review-fixes`):** PM-directed subset of the review, in two phases. **Phase 1** (commits `d2a1397` review-fixes + `31ddf5a` doc cleanup): **M1** applied (+ all 8 in-scope broken/redundant intra-doc links, and a bonus crate-wide pass fixing 14 *pre-existing* broken links outside `var_calling` — `cargo doc --no-deps` is now green for the whole crate; the original review under-reported this, having tailed the doc output); **M3** + **M4** applied (dead `PerGroupMerger` variant/`From` + dead `into_shared_ref_fetcher` deleted); **M5** applied (inert `ChunkSizingParams.chunk_genomic_span` field deleted — `DEFAULT_CHUNK_GENOMIC_SPAN` kept; the contamination chunked-stream genuinely uses it); **M7** applied (dangling `from_bam` doc links repointed/de-linked). Plus the PM's **always-parallel** decision: `drive_blocks_serial` removed, `drive_cohort_chunked` always runs `drive_blocks_parallel` with `n_workers = current_num_threads().max(1)` — which **moots M2** (no two paths to keep byte-identical). **Phase 2** (commit `c068af4`): the `from_psp/` → `var_calling/` merge above. **Dispositions:** **M2** moot (serial path removed); **M9** won't-fix as written — byte-identity is verified out-of-tree (human-in-the-loop benchmark + branch-vs-`main` diff during development), not via a committed in-tree golden VCF. Verification (in container): `fmt --check`, `clippy --all-targets --all-features -D warnings`, `doc --no-deps` all clean; `cargo test --lib` 1059 pass; cohort integration 17 pass. **Phase 3** (commits `79a6856` + `78a357b` + `cac838f`): **M6** (startup log now dumps every per-stage config), **M11** (5 `BlockQueue` concurrency tests), **M12** (`DustAheadPool::shutdown` logs a panicking worker), **M13** (mpsc in-flight bound documented), **Mi1** (Display source double-render), **Mi7** (`fill_block` doc), **Mi8** (Relaxed-atomic ordering comment) applied. **M8** attempted then reverted — `pub(crate)` un-masks a dead / test-only-code layer (~15 items) that needs its own removal pass first; the rationale is recorded at the module surface. Final verification: `fmt --check` / `clippy --all-targets --all-features -D warnings` / `doc --no-deps` clean; `cargo test --lib` 1064 pass; cohort integration 17 pass. **Still open:** M8 (dead-code cleanup → then internalise), M10 hot-path benches **added** (2026-06-01 perf review — `var_calling_chunk_loader`/`var_calling_run_window`), so M10 is satisfied and the perf wins it surfaced (H1–H3) are now the follow-up, Minors Mi2/Mi3/Mi4/Mi5/Mi6/Mi9/Mi10, and the out-of-scope `benches/psp_writer_perf.rs:385` panic.
- **Prior review:** [cohort_block_2026-05-29.md](doc/devel/reports/reviews/cohort_block_2026-05-29.md) — **Request-changes**: 5 Blockers (B1 writer-tmp leak on driver-error path; B2 full-chrom `Vec<u8>` materialisation in `compute_dust_mask_for_chrom` defeats per-chunk memory contract; B3 `NoSafeGap` retry is a no-op when `target_variants_per_chunk > 0`; B4 `NAllelesExceedsBitmask` silently rewritten as `DegenerateLikelihood { usize::MAX, … }`; B5 missing cross-driver byte-identity oracle test + missing unit tests for `drive_cohort_chunked` / `drive_one_chrom_generic` / `load_and_run_chunk_with_retry` / `emit_or_drop` / `compute_dust_mask_for_chrom`), 32 Major, 26 Minor, grouped Nits. `unsafe_concurrency` returned `No findings.` — the parallel-section soundness is statically enforced by the `Send + !Sync` typedef on `SharedRefFetcher`. Four open questions for the author (stable-API intent on the new pub data structs; streaming `drive_cohort_pipeline` oracle's long-term fate; sentinel-vs-`NonZero` policy for `target_variants_per_chunk` / `target_window_count`; filter-order equivalence vs streaming pipeline) gate several Major findings. Per-category audit trail at `tmp/review_2026-05-29_cohort_block/`.
- **Open (from the 2026-06-01 review — see report §6/§8 for full text + fixes):**
  - **M1** — Drop the explicit `(path)` target on the two redundant intra-doc links (`column_span_reader.rs:3`, `posterior_engine.rs:925`) to turn the `cargo doc` gate green.
  - **M2** — Add a serial-vs-parallel byte-identity test (drive the same fixture through `drive_blocks_serial` at 1 thread and `drive_blocks_parallel` at N threads; assert VCF bytes + every `ChunkDriverStats` field equal).
  - **M3 / M4** — Delete the dead `ChunkDriverError::PerGroupMerger` variant + `From` impl + import (`driver.rs:296-323`, `:52`); delete the dead `into_shared_ref_fetcher` (`worker.rs:702-715`, `mod.rs:41`).
  - **M5 / M7** — Resolve the inert `chunk_genomic_span` knob (delete vs. mark reserved; drop the phantom `MAX_CHUNK_SPAN_GROWTH` doc ref) per Q1; repoint or de-link the three dangling `from_bam` intra-doc links (`driver.rs:23,150`, `per_group_merger.rs:25`).
  - **M6** — Extend the startup config log to dump every effective per-stage default (DUST / grouper / per-group / `PosteriorEngineConfig`).
  - **M8** — `pub(crate)` sweep over the `from_psp` submodules + re-exports (gated on Q2); settle `DEFAULT_CHUNK_GENOMIC_SPAN` on one import path.
  - **M9 / M10** — Commit a golden VCF + byte-identity test (gated on Q3); add a `from_psp` chunk-driver + columnar-kernel bench with `// REGRESSION THRESHOLD: N%`.
  - **M11 / M12 / M13** — Add `BlockQueue` concurrency tests (mirror the `DustAheadPool` suite); surface/log `DustAheadPool::shutdown` worker-panic joins; document or `sync_channel`-bound the unbounded `mpsc` result/recycle channels.
  - **Minors** — `Display` source double-render (Mi1); `DownstreamFilterParams`/`ChunkSizingParams` default-binding docs (Mi2); co-dependent mapq-filter fields (Mi3); duplicated `DRIVER_PSP_BUFFER_BYTES`/`MAPQ_FILTER_MIN_READS_PER_SIDE` consts (Mi4); reconcile `chunk` vs `block` and stale `window`/`chunk.windows` vocabulary (Mi5/Mi6); `fill_block` doc says `bool` but returns `u32` (Mi7); `Relaxed`-atomic ordering comment (Mi8); direct `find_block_cut` tests (Mi9); malformed-PSP decode tests (Mi10). Plus the out-of-scope `benches/psp_writer_perf.rs:385` panic (breaks `cargo test --all-targets`) and the `genotype_order` placement question.
- **Open (from the 2026-05-29 review):**
  - **B1** — Add `CohortVcfWriter::abort()` that takes the tmp path it actually used and removes it; call from the error branch. Add an integration test injecting a mid-loop error and asserting no leftover tmp on disk.
  - **B2** — Stream `fetcher.iter_bases()?` directly into `sdust_mask_streaming` (the helper already accepts `Iterator<Item = io::Result<u8>>`). Add a regression test using a synthetic fetcher that panics on `.collect::<Vec<_>>()`.
  - **B3** — Thread `attempt_span` into the loader's `max_span` so retries can actually load more data, or push the retry inside the loader as "extend until safe gap or chrom cap". Add a test for a chunk that satisfies `target_variants_per_chunk` on the first attempt but has no safe gap.
  - **B4** — Add a dedicated `PerGroupMergerError::NAllelesExceedsBitmask { n_alleles, chrom_id, group_start, group_end }` variant; drop the `let _ = n_alleles;`. Add a test triggering the path with `cfg.max_alleles > MAX_BITMASK_ALLELES`.
  - **B5** — Add one integration test in [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs) that runs both `drive_cohort_chunked` and the streaming `drive_cohort_pipeline` on a multi-position fixture (≥3 samples, ≥1 MNP, ≥1 LH-cap site, ≥1 hom-REF group, ≥1 below-`min_alt_obs` site, ≥1 below-`qual_phred` site, ≥1 above-`mapq_diff_t` site); assert VCF bodies byte-equal **and** field-by-field equal counter sets. If the streaming driver is scheduled for removal, capture a checked-in golden VCF first.
  - **M1** — Replace `#[derive(Default)]` on `SampleColumns` with a hand-written `impl Default for SampleColumns { fn default() -> Self { Self::empty() } }`.
  - **M2** — Rewrite `prefetch_window_ref_bytes` so the outer `Vec<Vec<u8>>` resizes-without-dropping and each inner `Vec<u8>` is cleared-in-place rather than freshly allocated.
  - **M3** — Move `detect_compound_candidates_columnar`'s two `BTreeMap`s into `UnifyAllelesScratch` for per-group reuse; add a permutation-invariance proptest.
  - **M4 / M5 / M21** — Reshape `ChunkDriverError`: drop `#[from]` on `Io` / `PspRead`; rename variants by operation (`OpenPsp`, `WriteVcf { chrom_id, start, end, … }`, `FetchRefBases`, …); drop `: {0}` interpolation; add chrom/range context to `WriteVcf`. Replace `let _ = std::fs::remove_file(...)` with a structured `tracing::warn!` event.
  - **M6** — Surface `u32_from_usize` as a typed `ChunkLoadError::CsrOffsetOverflow` error or use `try_into().expect(...)` to panic loudly with a named invariant rather than wrap silently.
  - **M7 / M9 / Mi20 / Mi21** — Either propagate `ZeroTargetWindowCount` from the pre-pass (dropping the `.max(1)`) or name the default with `pub const DEFAULT_WORKER_WINDOWS_PER_CHUNK` and emit `tracing::debug!` when applied. Same shape for `target_variants_per_chunk == 0` (name `TARGET_VARIANTS_DISABLED` or lift to `Option<NonZeroU32>`). Add startup `tracing::info!` listing every effective `ChunkDriverParams` value. Add `chunks_with_fewer_windows_than_requested: u64` counter to `ChunkDriverStats`.
  - **M8** — Either reject `max_span < initial_span` with a new `ChunkLoadError::MaxSpanBelowInitial` variant or delete the `effective_initial_span` no-op chain and document the invariant.
  - **M10 + Mi26** — Fix the 2 in-scope unresolved-link errors and 3 redundant-link warnings (`worker.rs:17`, `worker.rs:236`, `driver.rs:105`, `driver.rs:108`, `worker.rs:9`).
  - **M11** — Add the two new `VarCallingArgs` fields to `benches/cohort_e2e_perf.rs:286`, `examples/profile_cohort_e2e.rs:152`, `examples/dhat_var_calling.rs:121` with the legacy defaults (`target_variants_per_chunk: 0`, `worker_windows_per_chunk: 1`). Optional follow-up: add a `VarCallingArgs::for_profiling(...)` constructor.
  - **M12 / M13 / M14 / M15 / M16 / M17 / M18 / M19 / M20 / M22 / M23 / M24 / M25 / M26 / M27 / M28 / M29 / M30 / M31 / M32** — see the report's §6 Findings section and §8 Missing tests. Highlights: split `load_and_run_chunk_with_retry`'s 19-param body into three phase helpers; group `load_chunk_from_iters`'s span/variant knobs into `ChunkLoadExtent`; add `SampleColumns::clone_from_columns` and use it for both carryover snapshot / restore loops; split per-window counters out of `ColumnarPipelineScratch`; drop the trailing `..` from the `AlleleSupportStats` destructure at `columns.rs:116`; delete or wire `chain_id_scratch` (`#[allow(dead_code)]`); validate `masked_intervals` sorting in `partition_window`; pin filter order in `emit_or_drop` with five per-category unit tests; pin `enforce_max_alleles_columnar` tie-break against the row-shape kernel; clamp `safe_end` to `chrom_one_past_end` on the last chunk; rename `*_cfg` vs `*_config` to a single form crate-wide; rename `shared_ref_fetcher` to `into_shared_ref_fetcher`; add tests for `SampleCountMismatch` / `CarryoverLengthMismatch` (both loader and pre-pass); add a `par_iter_mut` vs sequential equivalence test.
  - **Mi-class** (~26 minors): `#[non_exhaustive]` on the new pub data structs (gated on Open Question 1); `pub mod` → `pub(crate) mod` for every submodule that has no out-of-crate consumer (gated on Open Question 1); rename `MaterialisedChunk::clear_data` → `clear`; rename `WorkerSlot.output_buf` / `WorkerSlot.scratch` to carry domain nouns; collapse `Arc::new(StreamingChromRefFetcher)` to a borrow; drop `chunk.windows.clone()`; take `&PosteriorEngineConfig` in `run_window`; remove double-clone in `push_allele_into_scratch`; demote `pub` items with no caller; merge `Ok(idx) | Err(idx)` arms; split `unify_alleles.rs`/`worker.rs`/`loader.rs` along their existing internal sub-step boundaries; move `build_overlapping_variant_group` out of `worker.rs` and into `test_helpers.rs`; convert the three `Vec<Vec<_>>` jagged arrays to CSR; consider `OneBasedPos` / `OneBasedRange` / `ChromId` newtypes; group `ChunkDriverParams` along stage boundaries; add `// REGRESSION THRESHOLD: N%` to `benches/cohort_e2e_perf.rs`; add `--ignored` should-panic regression for `u32_from_usize` overflow.
  - **Nits** — single mechanical pass to clear the 16 in-scope clippy errors (`single_range_in_vec_init` ×8, `type_complexity` ×4, `bool_assert_comparison` ×2, `doc_lazy_continuation` ×2) plus add per-call-site justification comments to the 14 `#[allow(...)]` annotations (12 `clippy::too_many_arguments` + `clippy::arc_with_non_send_sync` + `clippy::needless_range_loop`).

#### Re-architected record-streaming pipeline (replaces the chunk-parallel rewrite)
- **Status:** fixes-applied (2026-06-08) — the 2026-06-08 review applied on
  branch `var-calling-review-fixes` (6 commits: M1/M4/M6/Mi6 deletions+move; M3
  typed errors; Mi4 rename + M7/Q1 doc sweep + Minors; M5/M8 tests; **M2**
  `pub(crate)` demotion + the dead-code cleanup it unmasks — delete fully-dead
  `is_empty`/`clear`, `#[cfg(test)]`-gate the eager-decode oracle the prior M8
  named as a prerequisite). All gates green (fmt / clippy `-D warnings` / `doc`;
  1021 lib + integration tests). Open questions resolved with the PM: P7-swap
  doc debt resolved now (Q1), dead code deleted (Q2), cheap in-tree tests over
  the out-of-tree oracle (Q3), modules demoted to `pub(crate)` (Q4). Deferred:
  the harder M8 leftovers (`StalledCut` / `compact_samples` straddler /
  `dust_mask_for_interval` sub-span paths + a multi-thread e2e byte-identity
  test). Prior: reviewed (2026-06-08) —
  re-reviewed on `main` post-merge (subtree review excluding `posterior_engine`;
  verdict Approve-with-changes, 0 Blockers / 8 Major / 13 Minor; see Latest
  review). Prior: fixes-applied
  (2026-06-05) — all 8 Major review findings +
  2 Minors applied across 3 commits (`b4e767c` review, `8210f46`
  M1/M2/M5/Mi2, `ed141ff` M3/M4/M6/M8/Mi8). Shipped to production on branch
  `re-architect` (Phase 7 swap, `1d34f85`). The columnar driver/worker/
  loader/columns/partition/two_pass/kernels chain is **deleted**; the only
  cohort `.psp` → VCF path is now the three-component record-streaming
  topology (producer `CohortChunkIntegrator` → W `VariantCaller` callers →
  `VcfWriter`), wired with two bounded crossbeam channels inside a
  `std::thread::scope`. Goal was **clarity** (runtime stages map to named
  logical sections), with byte-identical calls vs `main` as the hard
  contract and memory/wall as a measured guardrail.
- **Spec / plan:**
  [re_architecture_streaming_pipeline.md](doc/devel/implementation_plans/re_architecture_streaming_pipeline.md)
  (architecture, constraints, build strategy),
  [re_architecture_module_outline.md](doc/devel/implementation_plans/re_architecture_module_outline.md)
  (module/type map),
  [re_architecture_execution_plan.md](doc/devel/implementation_plans/re_architecture_execution_plan.md)
  (the byte-identity-gated phases),
  [re_architecture_p6_measurement.md](doc/devel/implementation_plans/re_architecture_p6_measurement.md)
  (Phase 6 vs `main` at scale).
- **Code:** [src/var_calling/](src/var_calling/) — new modules `types.rs`,
  `sample_reader.rs`, `cohort_integration.rs`, `pileup_overlaps.rs`,
  `em_posterior_calc.rs`, `vcf_writer.rs`, `pipeline.rs`, rewritten
  `mod.rs`; the numeric kernels (`per_group_merger`, `posterior_engine`,
  `variant_grouping`, `per_position_merger`, `dust_filter`,
  `contamination_estimation`) are carried verbatim. CLI rewired in
  [src/pop_var_caller/var_calling.rs](src/pop_var_caller/var_calling.rs);
  contamination estimator ported to the record-based `PerPositionMerger`
  in [src/pop_var_caller/estimate_contamination.rs](src/pop_var_caller/estimate_contamination.rs).
- **Tests:** 996 lib + all cohort integration pass (in container, commit
  `ef93b67`); `fmt --check` + `clippy --all-targets --all-features -D warnings`
  clean.
- **Latest review:** [var_calling_2026-06-08.md](doc/devel/reports/reviews/var_calling_2026-06-08.md)
  — **Approve-with-changes** (orchestrator skill, 10 categories; `main` @
  `35a6b67`; subtree review excluding the PM-deferred `posterior_engine`). 0
  Blockers, 8 Major, 13 Minor + Nits. All gates green (fmt / clippy `-D warnings`
  / `doc`; 1052 lib + integration tests pass; only `--all-targets` failure is the
  pre-existing out-of-scope `psp_writer_perf` bench panic). `unsafe_concurrency`
  clean (deadlock-free channels, order-independent parallel fold,
  `VariantCaller: Sync`). Drivers: M1 dead `PileupCohortChunk`, M2 over-broad
  `pub` surface (extends prior M8), M3 `PipelineError::Config(String)`
  flattening, M4 caller→producer coupling (`merge_compacted_samples`), M5
  no worker-count byte-identity test, M6 dead/mislabeled `SamplePspChunk` API,
  M7 stale verbatim-provenance pointing at deleted modules, M8 untested
  release-guards/byte-identity helpers. Per-category audit trail at
  `tmp/review_2026-06-08_var-calling/`.
- **Previous review:** [re_architecture_pipeline_2026-06-05.md](doc/devel/reports/reviews/re_architecture_pipeline_2026-06-05.md)
  — **Request-changes** (orchestrator skill, all 11 categories; HEAD
  `ef93b67`). 0 Blockers, 8 Major, 16 Minor, grouped Nits. No correctness
  Blockers in emitted calls; `unsafe_concurrency` clean (no `unsafe`,
  sound channel-close ordering, order-independent parallel-decode/serial-fold);
  `refactor_safety` verified `derive_is_kept`/`find_cut`/`chunk_cuts`/
  `merge_block_ranges`/`emit_or_drop`/`passes_min_alt_obs`/
  `merge_group_with_ref` are line-by-line identical to `main`. Per-category
  audit trail at `tmp/review_2026-06-05_re-architect-pipeline/`.
- **Latest perf review:** [perf_var_calling_cohort_2026-06-06.md](doc/devel/reports/reviews/perf_var_calling_cohort_2026-06-06.md)
  — **Run experiments** (orchestrator skill, 6 categories; HEAD `37e02c2`).
  First CPU sampling profile collected at this commit (macOS `sample`, T=1 30 s
  + T=8 6 s, N=50 real tomato cohort, tvpc=256). **3 Hot-path:** **H1** the T≥2
  wall gap vs `main` is *scheduler oversubscription* — the producer's rayon
  compaction pool and the crossbeam caller pool are two independent pools each
  sized to `--threads` (producer thread **64 % blocked on rayon** at T=8;
  `swtch_pri`/`__ulock_wait2`/`mutexwait` elevated vs T=1; tracks the
  on-par-T=1 / +6/+12/+15 %-T=2/4/8 matrix); **H2** the chain-id dead-weight
  (REF-allele chain_ids = 65 % of in-flight payload, ~96.6 % of all chain_ids,
  never read — `per_group_merger.rs:1304` skips allele 0; verified by 3
  categories); **H3** `ln_factorial` not `#[inline]` (out-of-line `bl` on the
  likelihood triple-loop, `cargo asm`-confirmed; cheap byte-identical apply).
  7 Likely (producer serial-on-main floor; `e_step_simd` sample-axis gather;
  `records_all` per-allele clones; per-chunk REF-span `Vec`; `Vec<Option<>>`
  AoS merge scan; E-step `is_finite` branch hoist; `binary_search` keep-mask →
  merge-walk), 4 Speculative, Notes. **No cohort end-to-end criterion bench
  exists** (deleted with from-bam) — building it is measurement-plan item 1.
  Build config already tuned (`lto=fat`/`codegen-units=1`/`panic=abort`);
  allocator A/B is the one unstruck build lever. Per-category audit trail +
  profiles at `tmp/perf_review_2026-06-06_var-calling-cohort/`.
- **Applied (2026-06-05, commits `8210f46` + `ed141ff`):** **M1** (doc gate
  green: 7 dead intra-doc links repointed/removed + 5 redundant targets
  dropped; verified under `RUSTDOCFLAGS=-D warnings`), **M2** (restored
  `current_command_line()` for the `##commandline` header), **M3**
  (zero-allele `.psp` rejected at the decode boundary via new
  `BlockHeaderInvariantKind::ZeroAlleleRecord` + `validate_n_alleles_column`
  helper + 4 unit tests), **M4** (typed `ChromRefFetchError` boxed through
  `ProducerError::Ref` — new `pub type RefFetchError`), **M5** (dropped the
  crate-wide `#![allow(dead_code)]` + deleted the 3 zero-caller items
  `into_reader`/`n_positions`/`chunk_cuts`), **M6** (`ProducerError::StalledCut`
  + `WriterError::MissingChunks` replace the two release-load-bearing
  `debug_assert!`s), **M8** (named `QUEUE_DEPTH_PER_WORKER`, startup log of
  resolved `workers`/`queue_cap`/`target_variants_per_chunk`, rewrote stale
  `--threads`/`--target-variants-per-chunk` CLI docs), **Mi2** (removed the
  inert `--target-variants-per-chunk` from `estimate-contamination` + its
  vacuous test), **Mi8** (refreshed stale "Phase 4 / `!Send`" module docs).
  Verified in container: `fmt` / `clippy --all-targets -D warnings` / `doc -D
  warnings` clean; 1000 lib + cohort integration tests pass.
- **Open (from the 2026-06-08 code review — see report §6/§8 for full text + fixes):**
  - **M1** — delete the dead `pub` `PileupCohortChunk` (types.rs); retarget its 3 doc refs at `RawCohortChunk`.
  - **M2** — demote the 5 zero-external-consumer rebuilt modules (`types`/`sample_reader`/`cohort_integration`/`pileup_overlaps`/`em_posterior_calc`) to `pub(crate)`; consider re-exporting `WriterStats` from `pipeline` to demote `vcf_writer` too.
  - **M3** — replace `PipelineError::Config(String)` flattening with operation-named variants carrying each typed cause via `#[source]`.
  - **M4** — move `merge_compacted_samples` (+ `ok_record`/`KeptRecordIter`) from the producer module into the caller.
  - **M5** — add a `var_calling_byte_identical_across_worker_counts` integration test + a direct `VcfWriter` out-of-order reorder unit test.
  - **M6** — delete dead `append_kept`/`SamplePspChunk::chrom_id()`; `#[cfg(test)]`-gate or relabel the eager `from_block`/`take_*` oracle path; fix the misdirecting module doc.
  - **M7** — pin each "copied verbatim from `two_pass`/`driver`/`worker`/`loader`" provenance to the source commit hash; rewrite the migration-future-tense `mod.rs`/`types.rs`/kernel docs to present tense (P7 swap has happened); resolve the deferred `Variant`-alias decision; fix stale `doc/devel/...` doc links.
  - **M8** — add the targeted unit tests for `MissingChunks`, `StalledCut`, `emit_or_drop` ordering/counters, `overlapping_groups`, `restrict_intervals_to_regions`/`dust_mask_for_interval`, `rebuild_fold` reduce-order independence, and the `compact_samples` straddler.
  - **Minors** — Mi1 startup-log provenance tags; Mi2 `BUFFERED_IO_CAPACITY` const drift; Mi3 collapse `DownstreamFilters` mapq fields to an enum; Mi4 rename `em_posterior_calc`→`variant_caller`; Mi5 `CohortPileupRecord`/`PileupCohortChunk` near-anagram (resolved by M1); Mi6 `RefSpan::empty()` fictional-doc; Mi7 `CallStats`→`WriterStats` exhaustive-destructure; Mi8 `// PANIC-FREE:` comments on release `.expect()`s; Mi9 caller-panic→typed-error vs fatal-by-design; Mi10 `SamplePspReader::new(r,0,1,1)` placeholder comment; Mi11 `sd`/`so`/`cd`/`co` scratch names; Mi12 hot-path bench regression threshold; Mi13 `passes_min_alt_obs` layout-coupling.
  - **Out of scope (follow-ups):** `dust_filter::is_masked` release-path passthrough on an unloaded mask (verbatim kernel correctness); the pre-existing `psp_writer_perf.rs:386` bench panic; the deferred `posterior_engine` review.
- **Open (perf — from [perf_var_calling_cohort_2026-06-06.md](doc/devel/reports/reviews/perf_var_calling_cohort_2026-06-06.md), verdict Run experiments):**
  - **Bench gap (do first)** — add `benches/cohort_var_calling_perf.rs` (criterion,
    `harness=false`) sweeping N∈{1,8,50}×T∈{1,2,8} at tvpc=256; without it no
    code-level perf finding is rankable cross-commit. Also fix
    `examples/profile_cohort_e2e.rs:161` (`target_variants_per_chunk: 0` →
    `256`) so the maintained driver profiles the production shape.
  - **H1 (wall gap, highest priority)** — collapse the rayon-compaction ⟂
    crossbeam-caller two-pool oversubscription; sweep cap-rayon-below-`--threads`
    / serialize-`compact_samples` / move-compaction-onto-callers, gated on
    byte-identical calls + a re-`sample` showing the oversubscription frames fall.
  - **H2 (memory win)** — drop the never-read REF-allele chain_ids (step 1:
    `records_all`/`records_for` REF slot → empty `Vec`; step 2: CSR stores
    non-REF alleles only). DHAT-gated, zero-VCF-diff gate.
  - **H3 (cheap apply)** — `#[inline]` `ln_factorial` + `#[cold]` tail; pure codegen, byte-identical.
- **Open (deferred — lower-priority Minors + the test additions):**
  - **M7** — Add unit tests for the writer reorder buffer (permuted +
    buffered-future-chunk; the `MissingChunks` gap path is now testable),
    `emit_or_drop` per-gate ordering/counters, a multi-thread
    `read_samples`-vs-reference test, and the new `StalledCut` guard. (PM:
    byte-identity is verified out-of-tree vs the previous version + the GIAB
    benchmark, so no in-tree A/B oracle is needed.)
  - **Minors** — dropped `chunks_loaded`/`avg_variants_per_chunk` run
    summary (Mi1, needs a chunk counter in `WriterStats`); ploidy-vs-inverted-range
    error precedence test in `merge_group_with_ref` (Mi4); missing
    `// PANIC-FREE:` on the thread joins (Mi5) and `fetch_ref_span`
    `binary_search().expect()` (Mi6); `ProducerError::Merge` not `#[source]`
    (Mi7); const-via-comment drift (Mi9); two `0 ⇒ ?` target conventions
    (Mi10); `em_posterior_calc` named after a sub-step (Mi11);
    `CohortPileupRecord`/`PileupCohortChunk` near-anagram (Mi12);
    co-dependent mapq-filter fields (Mi13); pipeline → CLI `VarCallingArgs`
    back-reference (Mi14); stale Cargo.toml crossbeam comment path (Mi15);
    unguarded `max_reach - first + 1` (Mi16). Plus the out-of-scope
    `benches/psp_writer_perf.rs:386` panic (breaks `cargo test --all-targets`;
    not the CI gate, which uses `--lib --tests`).
  - **Open questions — resolved by the PM (2026-06-05):** the `##commandline`
    change was unintentional → fixed (M2); `--target-variants-per-chunk`
    should not stay a no-op → removed from `estimate-contamination` (Mi2);
    byte-identity is verified out-of-tree (vs the previous version + the GIAB
    benchmark), so no in-tree oracle is needed; dead code should be removed →
    the `#![allow(dead_code)]` and its dead surface are gone (M5).

---

### Stage 6 — posterior engine

EM over merged records → final multi-sample VCF.

#### Posterior engine
- **Status:** fixes-applied
- **Spec sections:** `## Stage 6 — posterior engine` in [calling_pipeline_architecture.md](doc/devel/specs/calling_pipeline_architecture.md); background in [freebayes_posterior_gt_probs.md](doc/devel/specs/freebayes_posterior_gt_probs.md) and [gatk_em_calculation.md](doc/devel/specs/gatk_em_calculation.md)
- **Plan:** [posterior_engine.md](doc/devel/implementation_plans/posterior_engine.md)
- **Code:** [src/var_calling/posterior_engine.rs](src/var_calling/posterior_engine.rs)
- **Tests:** unit tests in the module + [tests/posterior_engine_integration.rs](tests/posterior_engine_integration.rs)
- **Impl reports:** [posterior_engine_2026-05-16.md](doc/devel/reports/implementations/posterior_engine_2026-05-16.md); perf history: [perf 2026-05-17](doc/devel/reports/implementations/posterior_engine_perf_2026-05-17.md), [SIMD analysis 2026-05-18](doc/devel/reports/implementations/posterior_engine_simd_analysis_2026-05-18.md), [perf wave-1 2026-05-18](doc/devel/reports/implementations/posterior_engine_perf_2026-05-18.md), [post-H4 profile](doc/devel/reports/implementations/posterior_engine_post_h4_profile_2026-05-18.md), [perf wave-2 (H4 + RecordScratch) 2026-05-18](doc/devel/reports/implementations/posterior_engine_perf_2026-05-18_v2.md)
- **Latest reviews:** [posterior_engine_2026-05-16.md](doc/devel/reports/reviews/posterior_engine_2026-05-16.md) — Request-changes: 2 Blocker test-gap findings, 11 Major, 13 Minor; [perf_posterior_engine_2026-05-18.md](doc/devel/reports/reviews/perf_posterior_engine_2026-05-18.md) — Run-experiments: 4 Hot-path (H4 homogeneous-fixation hoist subsumes H1 SIMD `log_sum_exp_2_x4 -INF` short-circuit; H2 `log_indep` cross-batch reuse; H3 GATK-style natural-log `softplus_neg` table for `log_sum_exp_2_x4`), 11 Likely (per-engine `RecordScratch` lift for 13 per-record allocs, `log_likelihoods` batch-of-4 transpose, autovec-blocking branch in `accumulate_expected_counts`, `#[inline(always)]` / `#[cold]` audits, in-bench self-validation), 11 Speculative (SQUAREM, warm-start across records, adaptive genotype pruning à la Octopus, cross-record tensor-style batching at same-shape, `f64x8` AVX-512 backend, PGO, `x86-64-v4` floor, …). Built on the May 17/18 bench evidence; the deep web + local-source (`gatk/`, `freebayes/`, `bcftools/`) research that fed it surfaced that none of the three reference callers vectorise ln/exp (GATK's `JacobianLogTable` is the closest precedent for H3) — `wide::f64x4`-native ln/exp already leads the field.
- **Latest fixes-applied:** [posterior_engine_2026-05-16_applied.md](doc/devel/reports/reviews/posterior_engine_2026-05-16_applied.md) — both Blockers fixed (B1 `NonFinitePosterior` test, B2 trivial-record row-sum bug + tests); 10 of 11 Majors applied (M3 config validation deferred — needs engine-vs-CLI boundary decision); 8 of 13 Minors applied (Mi5 perf bench, Mi6 golden test, Mi12 `with_config -> Result`, Mi13 fixture consolidation deferred as standalone follow-ups). 622 lib + 7 integration tests pass, clippy-clean on the in-scope files.
- **Open:**
  - **Mi5** — `benches/posterior_engine_perf.rs` regression-threshold criterion bench.
  - **Mi6** — golden `tests/golden/posterior_engine/*` fixtures locking emitted `f64` values across the representative-record matrix.
  - **Mi13** — fixture consolidation across the unit / integration crate boundary.
  - **Closed 2026-05-19** (cohort CLI slice): **M3** (engine-side
    config validation for `F`, pseudocounts) and **Mi12**
    (`Config::new -> Result<Self, _>`) — both shipped as
    `PosteriorEngineConfig::new` with full range validation. The
    CLI surface mirrors the same ranges in
    [src/pop_var_caller/cli/parsers.rs](src/pop_var_caller/cli/parsers.rs).
  - **Perf — applied 2026-05-18 (closed):** H1 (subsumed by H4), H2 (`log_indep` cross-batch reuse), H4 (homogeneous-fixation hoist) — see commit `013b49f`. L5–L7 (per-engine `RecordScratch` lift) — see commit `9594533`. Full session report at [posterior_engine_perf_2026-05-18_v2.md](doc/devel/reports/implementations/posterior_engine_perf_2026-05-18_v2.md).
  - **Perf — refuted by hardware-counter profile (demoted to Note):** L8/L9 (`log_likelihoods` / `posteriors` batch-of-4 transpose). Post-RecordScratch IPC = 2.97, L1 hit 99.66 %, branch-miss 0.11 % — the engine is execution-bound, not memory-bound. Cache layout has no remaining room to help.
  - **Perf — demoted to Speculative:** H3 (GATK-style `softplus_neg` table). H4 removed most `log_sum_exp_2_x4` calls; `__ieee754_log_fma` is now only ~3 % of cycles.
  - **Perf — promoted to Likely:** S5 (AVX-512 `f64x8` backend). With IPC=2.97 single-thread saturated, wider lanes are one of the few remaining single-thread levers (~2× on Skylake-X+/Zen 4-5 hosts).
  - **Perf next big lever — rayon-over-records.** Order-of-magnitude lever (16–128× on a server). RecordScratch lift is structured for per-thread ownership: each worker gets its own scratch. Separate plan; deferred from the wave-2 session per user's call.
  - **Perf remaining small-but-cheap:** L3 (`#[inline(always)]` on `log_sum_exp_*_x4`, gated by `cargo asm` check), L10 (in-bench `debug_assert!` on `diagnostics.iterations`).
  - **Closed 2026-05-19** (cohort CLI slice): end-to-end
    PspReader → … → PosteriorEngine integration test lives in
    [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs).

#### Contamination-estimation side-pass
- **Status:** fixes-applied
- **Spec:** [contamination_estimation.md](doc/devel/specs/contamination_estimation.md)
- **Plan:** [contamination_estimation.md](doc/devel/implementation_plans/contamination_estimation.md)
- **Code:** [src/var_calling/contamination_estimation.rs](src/var_calling/contamination_estimation.rs);
  Stage 6 consumer-side wiring in
  [src/var_calling/posterior_engine.rs](src/var_calling/posterior_engine.rs)
- **Tests:** unit + 3 proptests in the module +
  [tests/contamination_estimation_integration.rs](tests/contamination_estimation_integration.rs)
- **Impl report:** [contamination_estimation_2026-05-17.md](doc/devel/reports/implementations/contamination_estimation_2026-05-17.md)
- **Latest review:** [contamination_estimation_2026-05-17.md](doc/devel/reports/reviews/contamination_estimation_2026-05-17.md) — Request-changes: 2 Blockers, 16 Major, 23 Minor + grouped Nits.
- **Latest fixes-applied:** [fixes_applied_2026-05-17.md](doc/devel/reports/reviews/fixes_applied_2026-05-17.md) — 33 Applied (both Blockers + 14 of 16 Majors + 17 of 23 Minors), 2 Applied-with-adaptation (M15 cohort_alleles per-position alloc instead of cross-position scratch; Mi11 INF-delta test inverted), 1 Already-fixed (Mi17 covered by M9), 6 Deferred (M13 / M14 / Mi9 / Mi19 / Mi21 / Mi22). 668 lib + 109 integration tests pass; clippy + fmt clean. Perf comparison: 6 of 8 benches unchanged or improved, 2 regressed (var_calling_merger/dense +1.99 %, var_calling_grouper/overlap_extension +11.31 %) — both kept per skill policy (all fixes are correctness/test-coverage grade), most likely code-layout / inlining drift from the added test fixtures in adjacent modules.
- **Open:**
  - **M13** — `..Default::default()` in 14 new test sites; mechanical fix as a focused PR before any new `ContaminationEstimationConfig` field lands. See fixes_applied §M13.
  - **M14** — `OnlineEmState` god-struct split (3 sub-structs: `EmRunningStats`, `EmFrozenParameters`, `StabilityTracker`); pair with Mi21 (`estimate_contamination` long flat state machine). See fixes_applied §M14.
  - **Mi9** — Test fixture duplication between in-module and integration tests; needs cross-crate shared test module design.
  - **Mi19** — `DEFAULT_*` constants lack spec citations; pure-doc churn across 5 sites.
  - **Mi21 / Mi22** — paired with M14 and M13 respectively.
  - Partition-parallel scans across chromosomes (spec §"Cost and
    parallelism") — add only if real cohort runs show wall time
    matters.
  - `--external-allele-frequencies` for substituting a reference-
    panel `q_b` (v2 follow-up).
  - `--contamination-estimates`-only short path that estimates
    `q_b` while freezing user-supplied `c_s`
    (`ContaminationEstimateSource::Mixed` slot exists but is unused
    in v1).
  - **Closed 2026-05-19** (cohort CLI slice): CLI parser bindings
    shipped with `pop_var_caller estimate-contamination`.
  - Compound-allele `q_b = 0` choice (Assumption 2 in the impl
    report) — calibration impact in real cohorts is open.
  - Two perf regressions surfaced by the fix-application run
    (`var_calling_merger/dense`, `var_calling_grouper/overlap_extension`)
    — likely code-layout effects from new code in adjacent modules;
    revisit if real-cohort wall time matters.

#### Cohort VCF writer (Stage 6 sink)
- **Status:** fixes-applied
- **Plan:** [cohort_vcf_writer.md](doc/devel/implementation_plans/cohort_vcf_writer.md)
- **Code:** [src/var_calling/vcf_writer/](src/var_calling/vcf_writer/)
- **Tests:** 41 unit tests in the module + 4 integration tests in
  [tests/cohort_vcf_writer_integration.rs](tests/cohort_vcf_writer_integration.rs).
- **Impl report:** [cohort_vcf_writer_2026-05-18.md](doc/devel/reports/implementations/cohort_vcf_writer_2026-05-18.md)
- **Latest review:** [cohort_vcf_writer_2026-05-18.md](doc/devel/reports/reviews/cohort_vcf_writer_2026-05-18.md)
- **Latest fixes-applied:** [cohort_vcf_writer_2026-05-18_applied.md](doc/devel/reports/reviews/cohort_vcf_writer_2026-05-18_applied.md) — 31 Applied (1 Blocker + 14 of 15 Majors + 16 of 19 Minors + Nits) + 1 Applied-with-adaptation (Mi16) + 4 Deferred.
- **Open:**
  - `bcftools view` / `bcftools stats` manual smoke against real
    cohort data — synthetic integration-test fixture is too tiny to
    be meaningful.
  - **Closed 2026-05-19** (cohort CLI slice): end-to-end exercise
    via `pop_var_caller var-calling` lives in
    [tests/cohort_cli_integration.rs](tests/cohort_cli_integration.rs).
  - Tabix `.tbi` index alongside `.vcf.gz` — out of v1 scope; add
    when random-access matters.
  - `PL` (phred-scaled likelihoods) FORMAT field — needs Stage 5 →
    `PosteriorRecord` forwarding of `log_likelihoods`.
  - Per-sample contamination fraction in INFO — wires in once the
    cohort CLI threads `ContaminationEstimates` into the writer.
  - **M15** — `benches/vcf_writer_perf.rs` (criterion bench at
    1000 samples × biallelic SNP / `emit_gp = true` /
    encode-only); deferred from the review fix pass.
  - **Mi11 + Mi12** — `WriterConfig` → `CohortVcfWriterConfig`
    + `tool_string` → `source_label`/`tool_name` public-API
    renames; pair into a coordinated naming pass before
    publishing the crate.
  - **Mi14** — three-way test-fixture deduplication via a
    `pub(crate) test_fixtures` module + `test-support` feature
    flag; deferred.
  - **M13 follow-up** — fault-injection test for `new()` tmp
    cleanup; lands when the sink-injection seam exists.
  - **Mi16 follow-up** — full `BGZF_EOF` deduplication
    (currently shared in-crate but the integration test keeps
    its own copy).

#### Posterior engine — approximate-LUT inner loop
- **Status:** not yet implemented (config flag wired only)
- **Plan section:** [posterior_engine.md §"Approximation via precomputed lookup tables"](doc/devel/implementation_plans/posterior_engine.md)
- **Open:** evaluation methodology pinned in the plan; needs the
  exact-math engine bench numbers before deciding which candidates
  to land.

---

## SSR/STR caller (independent pipeline)

A second, independent caller — microsatellite/STR **length** genotyping from
aligned reads — sharing only low-level alignment I/O and a few numerical kernels
with the SNP caller above, never its records or math. Mirrors the SNP shape:
`ssr-catalog → catalog → ssr-pileup → .ssr.psp → ssr-call → VCF`. Design:
[ssr_genotyping.md](doc/devel/specs/ssr_genotyping.md) (model + why) and the
architecture docs ([overall](doc/devel/architecture/ssr_genotyping_architecture.md),
[shared types](doc/devel/architecture/ssr_shared_types.md),
[Stage 0 catalog](doc/devel/architecture/ssr_catalog.md)). Decision E and the
type model are settled; built in data-flow order (types → Stage 0 → Stage 1/2).

### Stage 0 — `ssr-catalog` (reference → catalog)

#### Shared domain types (`src/ssr/types.rs`)
- **Status:** fixes-applied
- **Design:** [ssr_shared_types.md](doc/devel/architecture/ssr_shared_types.md)
- **Code:** [src/ssr/types.rs](src/ssr/types.rs), [src/ssr/mod.rs](src/ssr/mod.rs) (commit `74b5a2d`; no separate impl report — commit message carries context)
- **Tests:** 16 unit tests in the module (6 original + 10 boundary/error tests from the fix run)
- **Latest review:** [ssr_types_2026-06-12.md](doc/devel/reports/reviews/ssr_types_2026-06-12.md) — Request-changes (1 Blocker, 1 Major, 5 Minor + missing tests)
- **Latest fixes-applied:** [fixes_applied_2026-06-12.md](doc/devel/reports/reviews/fixes_applied_2026-06-12.md) — **Completed**: all 9 findings resolved (B1 + M1 + Mi1–Mi5 Applied; Nit1 kept; Nit2 already consistent). M1 → private `Locus` fields + validated `Locus::new`/`LocusError`; Mi5 → `pub(crate)` surface (with a temporary module `#![allow(dead_code)]` until Stage 0 consumers land). `cargo doc`/ssr clippy/fmt clean; 1039 lib tests pass.
- **Open:**
  - Remove the temporary `#![allow(dead_code)]` in `src/ssr/mod.rs` once `ssr-catalog` (Stage 0) wires these types up (tracked as Mi5 follow-up in the fix report).

#### Catalog builder (`src/ssr/catalog/`) — *planned*
- **Status:** planned
- **Plan:** [ssr_catalog.md](doc/devel/implementation_plans/ssr_catalog.md) (implementation sketch: files, structs, fn signatures)
- **Notes:** detector = lh3/TRF-mod (shell-out via temp files, no FFI), genome-wide; post-process drops compound/bundled loci (GangSTR-style, no split); worker-per-contig with an ordered collector; `--num-chroms-in-parallel` is a speed⇄RAM knob.

### Stage 1 — `ssr-pileup` (per-sample evidence extraction)

#### `ssr-pileup` stage
- **Status:** **implemented — rebuilt as Mark-2** (2026-06-17, branch `ssr-pileup-mark2`). The Mark-1 detail below (the reference-anchored rung model + its perf/review history) is **superseded**, kept for chronology.
- **Mark-2 rebuild (current):** [ssr_pileup_mark2_2026-06-17.md](ia/reports/implementations/ssr_pileup_mark2_2026-06-17.md). Empirical-candidate model — observed sequences, the reference is only a coordinate frame, no on/off-ladder, no Stage-1 likelihood. Design: [ssr_ladder_model.md](doc/devel/architecture/ssr_ladder_model.md) + [ssr_pileup_mark2.md](doc/devel/architecture/ssr_pileup_mark2.md); plan [ssr_pileup_mark2.md](doc/devel/implementation_plans/ssr_pileup_mark2.md). New `src/ssr/` (`types`/`catalog`/`pileup::{footprint,fetch_reads,alignment,locus_tally,driver}`) + rewritten `psp/registry_ssr.rs`; old `src/ssr_mark1/` + its bench/example deleted. fmt/clippy `-D warnings` clean; 1116 lib + integration + doctests green (end-to-end + thread-determinism). **Open:** calibrate `MIN_REGION_Q1` (15) / `MAX_READS_PER_LOCUS` on real data; a Mark-2 bench (the Mark-1 one was deleted); spec §4.2/§4.3/§5.1 amendment; then Stage 2 (`ssr-call`).
- **Latest fixes-applied:** [fixes_applied_2026-06-17.md](doc/devel/reports/reviews/fixes_applied_2026-06-17.md) — B1 (doc gate restored), M5 (shared `decode_cram_container` for refill + caching path; folds in M2 EOF-stop + Mi3/Mi6), M3 (ContainerCache eviction/dedup unit tests), Mi1/Mi4/Mi5/Mi7/Mi8/Mi9/Mi10 + doc nits. M1/M4 reworded-as-doc per the review discussion (no `.ssr.psp` baseline → H1/H2 ULP-difference is inert; cache cap is non-output-affecting). All gates green: fmt/clippy `-D warnings`/doc clean, 1171 lib (+4) + 45 integration tests. **Open:** M3 multi-container *fetch* test (needs a multi-container CRAM fixture), Mi2 (cfg-gate bench seam), Mi11 (bench baseline guard), Mi12 (driver error context — re-assess vs `par_chunks`).
- **Architecture:** [ssr_pileup.md](doc/devel/architecture/ssr_pileup.md) (every structural question decided)
- **Plan:** [ssr_pileup.md](doc/devel/implementation_plans/ssr_pileup.md) (implementation sketch: build order, modules, structs, fn signatures)
- **Build order:** (1) allele types in `types.rs` → (2) lift `normalize_alleles` to `src/norm_seqs/` → (3) container SSR schema → (4) stage modules (`count_repeats` → `pair_hmm` → `candidate_generation` → `triage` → `fetch_reads` → `mod`).
- **Latest code review:** [ssr_pileup_2026-06-17.md](doc/devel/reports/reviews/ssr_pileup_2026-06-17.md) (branch `ssr-pileup-review`, 11 categories) — **Request-changes**: 1 Blocker, 5 Major, 12 Minor + Nits. **B1** `cargo doc` gate broken (fetch_reads.rs:9 still links the removed `AlignmentFile` — verified failure). **M1** (the crux): the H1/H2 math rewrites (`ln_1p` fold + single-pass `ln_sum_exp3`) are ULP-different from `main` (probe: 4.74% of `ln_sum_exp2` evals differ in f64 bits) and the "byte-identical" claim is **asserted, never tested vs `main`** — the in-tree bit-identity test compares the new code against itself; run the shipped `ssr_psp_concordance` with window pinned to 10 to settle it. **M2** `CachingCramReader` `.crai` walk diverges from `refill` on `decoded==0` (continues vs stops); **M3** cache eviction/multi-container path untested (single-container gate fixture); **M4** `DEFAULT_MAX_CACHED_CONTAINERS` not in the `.ssr.psp` header (uninspectable); **M5** `.crai`-walk/decode duplicated with `refill` (the source of M2's drift). The shared-prefix DP (P1) is **provably bit-identical + well-gated**, and the `map_init` per-worker concurrency is sound. Audit trail `tmp/review_2026-06-17_ssr-pileup-review/`.
- **Latest perf review:** [perf_ssr_pileup_2026-06-16.md](doc/devel/reports/reviews/perf_ssr_pileup_2026-06-16.md) (branch `ssr-pileup-review`) — **Apply the listed wins.** First measurement of the realignment hot path: new criterion bench [benches/ssr_pileup_perf.rs](benches/ssr_pileup_perf.rs) + profiling driver [examples/profile_ssr_pileup.rs](examples/profile_ssr_pileup.rs), both driving a new `#[doc(hidden)]` [bench_harness](src/ssr/pileup/bench_harness.rs) seam. The `sample` profile pins **~74.6% of self-time on `exp`+`log`** in the pair-HMM `ln_sum_exp2/3`; cost is **perfectly linear in the `2·window+1` rung count** (13.3 ms/rung). Allocation is cold (~0.02%). Audit trail `tmp/perf_review_2026-06-16_ssr-pileup/`.
- **Impl reports:**
  - Task 1 — allele representation (`Allele`/`NormalizedSeq` + `to_sequence`/`repeat_count`) in [src/ssr/types.rs](src/ssr/types.rs): [ssr_pileup_task1_allele_types_2026-06-15.md](ia/reports/implementations/ssr_pileup_task1_allele_types_2026-06-15.md)
  - Task 2 — lift `normalize_alleles` (+ `IndexRange`) to the shared [src/norm_seqs.rs](src/norm_seqs.rs); SNP CIGAR path ([src/pileup/walker/indel_norm.rs](src/pileup/walker/indel_norm.rs)) now wraps the kernel. Behaviour-preserving; full lib suite green: [ssr_pileup_task2_norm_seqs_lift_2026-06-15.md](ia/reports/implementations/ssr_pileup_task2_norm_seqs_lift_2026-06-15.md)
  - Fast-path counter (reordered ahead of task 3) — `count_pure_tiling` core in [src/ssr/pileup/count_repeats.rs](src/ssr/pileup/count_repeats.rs); needs only `types`. Triage-typed wrapper deferred. [ssr_pileup_count_repeats_2026-06-15.md](ia/reports/implementations/ssr_pileup_count_repeats_2026-06-15.md)
  - Slow-path forward — `forward` + `HmmModel` + `PairHmmScratch` in [src/ssr/pileup/pair_hmm.rs](src/ssr/pileup/pair_hmm.rs); 3-state log-space pair-HMM, Dindel emission. Constant gap-open (homopolymer-indexed table deferred) + unbanded (banding deferred), both calibration per arch §14. [ssr_pileup_pair_hmm_2026-06-15.md](ia/reports/implementations/ssr_pileup_pair_hmm_2026-06-15.md)
  - Candidate generation Job 1 — `CandidateAllele` + `build_rungs` (on-ladder rungs `left_flank + motif×L + right_flank`) in [src/ssr/pileup/candidate_generation.rs](src/ssr/pileup/candidate_generation.rs); composes `Allele::to_sequence`. Job 2 (off-ladder normalization) deferred pending a contract decision. [ssr_pileup_candidate_generation_2026-06-15.md](ia/reports/implementations/ssr_pileup_candidate_generation_2026-06-15.md)
  - `score_candidates` in [src/ssr/pileup/pair_hmm.rs](src/ssr/pileup/pair_hmm.rs) — joins `build_rungs` + `forward` into the dense per-read `Qᵣ` (one `(allele, log-lik)` per candidate, raw scores; pruning/renorm is the aggregator's job). Small addition, context in commit; 2 added tests.
  - Candidate generation Job 2 — `build_offladder` + `normalize_offladder` (off-ladder candidates) in [src/ssr/pileup/candidate_generation.rs](src/ssr/pileup/candidate_generation.rs); contract A + B1, canonical form is the **verbatim full tract** (left-alignment is provably a no-op on a full-tract key; clean flanks). **`norm_seqs` not needed for B1** — task-2's second consumer doesn't materialize. `candidate_generation` now complete. [ssr_pileup_candidate_generation_job2_2026-06-15.md](ia/reports/implementations/ssr_pileup_candidate_generation_job2_2026-06-15.md)
- **Perf (review + applied wins, 2026-06-16):**
  - **Applied (≈−54% on realignment):** H1 fold `exp(0)` in `ln_sum_exp2`; H2 single-pass `ln_sum_exp3`; **P1 shared-prefix DP** in `score_candidates` (score the rungs' longest common prefix once, continue each tail from the saved seam). **P1 is provably bit-identical** to per-candidate `forward` (gated by `score_candidates_is_bit_identical_to_per_candidate_forward`). **H1/H2 are algebraically equivalent but not bit-identical to the prior form** (`ln_1p` is more accurate near zero; the single-pass reduction reorders the sum) — f64 scores differ within ULP, absorbed by the `f32` profile storage. Since there is **no released `.ssr.psp` baseline**, this is not an output-compat concern; the first produced catalog is the baseline (review M1, [fixes_applied_2026-06-17.md](doc/devel/reports/reviews/fixes_applied_2026-06-17.md)).
  - **Fast-path investigation (closed, 2026-06-16):** [ssr_fastpath_investigation_2026-06-16.md](doc/devel/reports/reviews/ssr_fastpath_investigation_2026-06-16.md) — measured on a real ch01 catalog + tomato CRAM (built the L9 fixture: vendored `trf-mod` + `ssr-catalog`). **Verdict: do NOT build a separate exact-match fast path** — safe predicates have ~21% recall + break byte-identity; the high-recall `observed_count` predicate would emit the wrong length on the ~11% of reads where the HMM *corrects* the pre-probe count (the realign-everything rationale, confirmed real). Instead: (a) **applied — `DEFAULT_WINDOW` lowered 10 → 6** (validated on real data: 0.48% of read-calls / 0.06% of loci change, ~40% cheaper realignment; concordance via the new `ssr_psp_concordance` ignored test); (b) on real data the **fetch/catalog-walk, not realignment, is the end-to-end wall** (confirms L5/L7) — the next lever.
  - **Open — P1 follow-up:** incremental motif-by-motif seam advance would also share the per-rung tract tail (a further ~20%), more complex; deferred.
  - **Fetch path (L7) — step 1 DONE (2026-06-17), ~17× faster, byte-identical.** Design: [ssr_pileup_read_buffer.md](doc/devel/architecture/ssr_pileup_read_buffer.md). Added a per-worker [`CachingCramReader`](src/bam/segment_reader.rs) (FIFO cache of decoded CRAM containers) + `WorkerReader` (CRAM=cache, BAM=pooled), wired into `fetch_locus_reads`/`process_locus` via `map_init`. Decodes each container once instead of once-per-locus. Gate: `caching_cram_reader_matches_per_call_path` (reader-level byte-identity) + end-to-end `.ssr.psp` diff on the real fixture = **0 record diffs**. **Measured (cap A/B + decode counter):** cap=0 (no cache) 8931 decodes / 85.3s → cap=3 341 decodes / 4.23s = **~26× fewer decodes, ~20× wall, byte-identical**. **Step 2 (slice-grouping) DROPPED — measured unnecessary:** cap=3 is already at the decode floor (cap=64 → 420, no improvement), so cross-thread boundary re-decodes are noise-level; the doc's prediction held. Decode redundancy eliminated; end-to-end factor is smaller on high-coverage runs (realignment, unchanged, dominates there).
  - **Open — L1 (experiment):** band the DP (arch §5.5 `PAIR_HMM_BAND_BP`); approximation, gate on genotype concordance. With P1 landed, the remaining `exp` (~48% self-time) is the natural target for L1 / approx-math (S1).
  - **Open — B1 (build):** add an `(aarch64, linux)` `target-cpu` floor to `.cargo/config.toml` (prod misses Neoverse SIMD; macOS bench flatters prod).
  - **Open — L9 (measurement gap):** no end-to-end bench — fetch path + driver batch loop unmeasured; build one from the in-process `stage1_fixture` (no `trf-mod`).
  - Design lever (route to SSR owners, not perf): the `window` candidate multiplier touches the realign-everything contract + Stage 2 support points.
- **Open:**
  - **Task 3 (container generalization, arch §10) deferred by choice** — a large multi-step refactor of the production `.psp` writer (58KB) + reader (136KB); gates the `.ssr.psp` writer/round-trip but not the stage's compute modules. Trait-vs-builders fork (§10.7 Q1) still to be decided at its step 1.
  - `OnLadder::to_sequence` is a clean tiling; imperfect-locus interruptions deferred to `candidate_generation.rs`. The SSR off-ladder adapter onto `norm_seqs` lands with `candidate_generation` (first consumer beyond the SNP path).
  - **Architecture revised (2026-06-15): realign-everything** — v1 realigns every spanning read (pair-HMM), dropping the CIGAR-trusting two-tier fast/slow gate; the direct-count fast path (`count_repeats`, built + parked) becomes a **measured optimization the user requires be tried**. Docs amended: arch `ssr_pileup.md` §2/§14, plan §4.
  - `triage` complete (realign-everything) — `find_longest_stretch` pre-probe + `read_footprint`/`brackets`/`extract_region`/`triage_read` in [src/ssr/pileup/triage.rs](src/ssr/pileup/triage.rs): coverage classification (footprint position + clips, no flank-byte match) → region extract → window centre. Read seam = `MappedRead`. Recovers soft-clipped long alleles via the clip-included pre-probe (test: mapper's 3 units → true 6). 22 tests. [ssr_pileup_triage_2026-06-15.md](ia/reports/implementations/ssr_pileup_triage_2026-06-15.md)
  - Per-read analysis wired — `analyze_read` + `ReadOutcome` in [src/ssr/pileup/read_analysis.rs](src/ssr/pileup/read_analysis.rs): composes `triage_read` → `build_rungs` → `score_candidates` → dense `Qᵣ` (Spanning) / Flanking / InRepeat. Reuses scratch buffers. Integration-tested incl. the realign-everything win (soft-clipped 6-unit allele the mapper called as 3 → ranks rung 6). 4 tests. Off-ladder candidate generation deferred (needs anchored observed-tract isolation; rare). Context in commit.
  - `locus_record` aggregation — `aggregate` + `SsrLocusRecord` + `QcCounts` in [src/ssr/pileup/locus_record.rs](src/ssr/pileup/locus_record.rs). **Storage model revised → all-CSR** (no histogram, no weight): every spanning read stored as one pruned + renormalized `Qᵣ` profile; `AMB_LL_DROP` pruning; renormalize provably lossless (`Z_r` cancels); derive `n_spanning`/`n_flanking`/`n_frr`, take `depth`/`n_filtered`/`mapped_reads`; drop vestigial `n_flank_indel`. In-memory (CSR flattening = deferred container's job). Diverges from spec §4.3 (amend). 6 tests. [ssr_pileup_locus_record_2026-06-15.md](ia/reports/implementations/ssr_pileup_locus_record_2026-06-15.md)
  - `fetch_reads` started — the per-locus depth cap: `Reservoir<T>` (Algorithm R) + deterministic `locus_seed(chrom, start)` (FNV-1a) + inline `SplitMix64` PRNG, in [src/ssr/pileup/fetch_reads.rs](src/ssr/pileup/fetch_reads.rs). Net-new (SNP has only a column cap); byte-identity-critical (deterministic seed + caller's fixed read order, §8.3/§8.4). `MAX_READS_PER_LOCUS=1000` placeholder. 6 tests. Context in commit.
  - **Container generalization** (the deferred §10 refactor — gates Stage-1 output): implementation sketch written, [psp_container_generalization.md](doc/devel/implementation_plans/psp_container_generalization.md) — `PspSchema` trait (or two-builders, decide at step 1), 5-step sequence (parameterize → `kind` tag → interval index → `registry_ssr`+`SsrLocusRecord` → round-trip), SNP e2e as the gate each step. Reconciles with all-CSR (SSR table = `amb_*` CSR + QC scalars, no `hist_*`). Next: step 1 (writer parameterization, SNP-only).
  - Stage modules still to build (need the container writer / real I/O): the `fetch_reads` catalog-walk + `query` driver (mirror SNP `run_pileup`: load handles once, share repo, clear per contig) + admission gate (reuse triage footprint) + bundles; the driver (`mod.rs run()`); the container schema/writer (deferred refactor — flattens `SsrLocusRecord` profiles → CSR). Single-threaded semantics first, then the fetcher-thread/pool (determinism-gated). Deferred: off-ladder wiring; the measured `count_repeats` fast-path shortcut.
  - `pair_hmm` follow-ups: homopolymer-indexed gap-open + banding (calibration/optimization, arch §14).

---

## Standing project-wide items

Not bound to a specific stage block; pick up whenever the active feature
list is clear.

- **BED-region skip.** Add a Stage 1 CLI flag that takes a BED file of
  regions to skip (avoid pathological regions). Source:
  [doc/devel/TODO.txt](doc/devel/TODO.txt).
- **Phase-chain integration tests.** Add integration tests asserting that
  phase chains are correctly carried through the Stage 5 likelihood
  calculation. Source: [doc/devel/TODO.txt](doc/devel/TODO.txt).
- **Parallel-optimization integration perf benches.** Build end-to-end
  criterion (or equivalent) benches that measure wall time for the two
  pipeline arms whose parallelism still needs tuning: CRAM → `.psp`
  (Stage 1, `pileup` / `var-calling-from-bam`) and `.psp` → cohort VCF
  (Stages 3–6, `var-calling`). These integration benches are the
  ground truth for the deferred parallelisation-tuning pass
  (rayon-over-records, `--per-group-batch-size`, per-group batch
  sizing) — micro-benches alone can't catch end-to-end scaling
  artefacts.
  - **`.psp` → cohort VCF arm — shipped 2026-05-19:**
    [benches/cohort_e2e_perf.rs](benches/cohort_e2e_perf.rs). Two
    bench-group families: `cohort_e2e_core/*` times
    [`drive_cohort_pipeline`](src/pop_var_caller/cohort_driver.rs) in
    isolation (PSP open + FASTA verify outside the timed region) with
    sub-groups `scaling_samples` (N ∈ {10, 64, 256}), `scaling_region`
    (L ∈ {1 000, 5 000, 20 000}), and `scaling_threads` (T ∈ {1, 2,
    4, max-cores} via per-bench local `rayon::ThreadPool` +
    `pool.install(...)`); `cohort_e2e_full/*` times
    [`run_var_calling`](src/pop_var_caller/var_calling.rs) end-to-end
    over `scaling_samples` + `scaling_region`. The full group cannot
    sweep thread count within one `cargo bench` invocation —
    `configure_rayon_pool` / `ThreadPoolBuilder::build_global` is
    once-per-process; run separate invocations under
    `RAYON_NUM_THREADS=N`. Bench-surface lift:
    [`drive_cohort_pipeline`](src/pop_var_caller/cohort_driver.rs) +
    `CohortPipelineParams` promoted from `pub(crate)` to
    `#[doc(hidden)] pub`. 15 bench variants pass `cargo bench --bench
    cohort_e2e_perf -- --test`.
  - **`.psp` → cohort VCF arm — perf review against real data
    2026-05-20:**
    [perf_psp_to_vcf_2026-05-20.md](doc/devel/reports/reviews/perf_psp_to_vcf_2026-05-20.md).
    Used the bench above + an `examples/profile_cohort_e2e.rs`
    one-off + `perf record` on real tomato (SL4.0)
    `SRR7279725.small.psp × N=10` to produce the headline
    diagnosis: pipeline is essentially single-threaded
    (T=1=12.7s ≈ T=16=13.6s); DUST 33 %, allocations 21 %, PSP
    decode 15 %; per-group merger + posterior together <2 %. Seven
    Hot-path + 12 Likely findings; **H1 per-chromosome parallelism**
    is the order-of-magnitude lever.
  - **`.psp` → cohort VCF arm — H1 per-chromosome parallelism shipped
    2026-05-20:**
    [cohort_per_chromosome_parallel_2026-05-20.md](doc/devel/reports/implementations/cohort_per_chromosome_parallel_2026-05-20.md).
    Realised **3.85× wall-time reduction at T=13** on the multi-chrom
    real-data fixture
    (`tmp/SRR7279727.multichrom.psp` — 2 Mbp from each of 13 SL4.0
    chroms via `samtools view --regions`, N=10 cohort: 106.6 s →
    27.7 s). Workload imbalance (ch00 unplaced/decoy reads at
    13× the median per-chrom count) gates the ceiling below the
    plan's predicted 6–10×; L5 contention is the next ceiling
    (acknowledged for follow-up). Includes L1 (per-group inner
    `par_iter` removed) and a new pure-Rust bgzf-aware concat
    module (`src/var_calling/vcf_writer/concat.rs`).
  - **CRAM → VCF arm (`var-calling-from-bam`) — REMOVED 2026-06-01.**
    The direct single-sample BAM/CRAM → VCF subcommand was deleted
    entirely (see the Current-focus "Last completed task"); the only
    route to a VCF is now `pileup` → `.psp` → `var-calling`. The
    history below is retained for the record. Per-chromosome
    parallelism shipped 2026-05-24:
    [var_calling_from_bam_per_chromosome_2026-05-24.md](doc/devel/reports/implementations/var_calling_from_bam_per_chromosome_2026-05-24.md);
    plan
    [var_calling_from_bam_per_chromosome.md](doc/devel/implementation_plans/var_calling_from_bam_per_chromosome.md).
    Four-commit PR (`29b28e8` → `e5261ba` → `050da41` → `a858832`):
    new `src/bam/index_preflight.rs` + `--build-map-file-index`
    flag (opt-in `.crai` auto-build; off by default with a
    samtools-pointing error), new
    `CramMergedReader::query` indexed per-contig variant +
    `OwnedIndexedCramRecords` iterator, new
    `process_one_chromosome_from_bam` per-chrom worker, and a
    `run_var_calling_from_bam` reshape that retires the serial
    `run_cohort_pipeline_for_single_sample` + `PerChromRecordsIter`
    helpers in favour of `rayon::par_iter` over chromosomes +
    `vcf::concat::concat_fragments` (no new file-format module —
    reuses the cohort H1's pure-Rust bgzf concat). Net diff
    +1320 / -520 across the four commits; 898 lib + every
    integration test pass; clippy + fmt clean. Wall-time
    validation on real multi-chrom tomato CRAMs (the analogue of
    cohort H1's 3.85× at T=13) is deferred: the
    `examples/profile_from_bam_e2e.rs` + `benches/from_bam_e2e_perf.rs`
    infrastructure (commit 5 of the plan) was scoped out.
  - **CRAM → `.psp` arm (`pileup`): not planned** (design call
    2026-05-24). The typical PSP workflow runs N independent samples
    and the orchestrator (Snakemake, Nextflow, GNU parallel)
    already provides the per-sample parallelism; a single `pileup`
    invocation processes exactly one sample and stays serial.
    Revisit only if a real workload appears where one sample's
    `.psp` build dominates a pipeline's wall time.

