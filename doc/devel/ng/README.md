# `ng` — next-generation caller docs

All documentation for the **ng** effort (the step-decomposed, benchmark-driven
algorithm lab — see [`spec/ng_proposal.md`](spec/ng_proposal.md)) lives here, split
by document kind:

- **`spec/`** — proposals and design specs (what to build and why).
  - [`ng_proposal.md`](spec/ng_proposal.md) — the plan: the step-decomposed caller
    taxonomy, the benchmark strategy, and the single-phase lab.
  - [`read_filtering.md`](spec/read_filtering.md) — step 1 (the whole-read keep/drop
    prelude) + the ng foundations it settles (skeleton, `types.rs` seed, conventions).
  - [`read_preparation.md`](spec/read_preparation.md) — step 2 (the per-read transform),
    the shared discipline; two path specs:
    [`read_preparation_generic.md`](spec/read_preparation_generic.md) (SNP/indel: left-align +
    BAQ → `PreparedRead`) and
    [`read_preparation_ssr.md`](spec/read_preparation_ssr.md) (STR: tract extraction →
    `SsrTractObs`).
  - [`typed_regions.md`](spec/typed_regions.md) — step 3 (the typed-region generator): walks
    the reference and cuts it into `TypedRegion`s (SsrSegment / SsrBundle / Generic / Satellite).
    First integrating spec — stands on `RefSeq`, the tandem-repeat scanner, and the STR catalog.
  - [`typed_regions_cli.md`](spec/typed_regions_cli.md) — the **`pop_var_caller_exp`** binary
    (ng's command surface, kept out of the production CLI) and its first subcommand,
    `type-regions`: step 3's walk driven from the command line, writing the genome's
    partition to a text file (kind + motif + repeat count per region).
  - [`ref_seq.md`](spec/ref_seq.md) — the `RefSeq` reference-sequence accessor
    (foundational infra: resident + streaming + in-memory impls). Read filtering #8 and
    the pileup depend on it.
  - [`reference_info.md`](spec/reference_info.md) — reading a reference's info (contig table +
    content digest + reconstructed index) from a FASTA (optionally checked against a `.fai`) or
    a `.fai` alone; with the MD5s, the fasta↔fai check, a caller-held cache, and a `.fai` writer.
    Builds the `ContigList` every `RefSeq` impl takes but none can build.
  - [`alignment_file.md`](spec/alignment_file.md) — **one** alignment file: open-and-validate
    (`SO`, `@SQ`↔reference, index, `@RG SM`) and serve a region as an ordered, filtered read
    stream. Closes the `@SQ` permutation hole and the sort-order check; step 1's input edge.
  - [`sample_reads.md`](spec/sample_reads.md) — **the sample**: k files (usually several
    experiments) merged into one coordinate-ordered stream, with the cross-file checks. Stands
    on `alignment_file.md`.
- **`arch/`** — architecture (the shared types and the interfaces implementations
  plug into).
  - [`ng_step_interfaces.md`](arch/ng_step_interfaces.md) — the common domain
    newtypes and one swappable trait per pipeline step.
  - [`module_layout.md`](arch/module_layout.md) — the `src/ng/` module tree: one
    folder per step (trait + impls + tests together), shared vocabulary, `bench/`.
  - [`read_filtering.md`](arch/read_filtering.md) — step 1's types & interfaces,
    distilled (the code-facing companion to the spec).
  - [`typed_regions.md`](arch/typed_regions.md) — step 3's types & interfaces (the
    typed-region generator); companion to `spec/typed_regions.md`.
  - [`typed_regions_cli.md`](arch/typed_regions_cli.md) — the `pop_var_caller_exp` /
    `type-regions` types & interfaces (the CLI that drives step 3); companion to
    `spec/typed_regions_cli.md`.
  - [`reference_info.md`](arch/reference_info.md) — the reference-info reader's types &
    interfaces (`ReferenceInfo`/`ContigInfo`, the cache, the writer, the background verify);
    companion to `spec/reference_info.md`. Foundational infra, not a step.
  - [`alignment_file.md`](arch/alignment_file.md) — `AlignmentFile` (the validated handle),
    the region-query `RecordSource` impls, the order guard, the reader pool;
    companion to `spec/alignment_file.md`.
  - [`sample_reads.md`](arch/sample_reads.md) — `SampleReads`, the argmin merge and its
    per-read budget; seeds the shared `GenomePosition`. Companion to `spec/sample_reads.md`.
- **`impl_plan/`** — step-by-step implementation plans (build order, not new design).
  - [`foundations.md`](impl_plan/foundations.md) — the first ng code: skeleton,
    `types.rs` seed, and the `RefSeq` accessor (three impls).
  - [`read_filtering.md`](impl_plan/read_filtering.md) — step 1: the `read/` module,
    the cascade, the `RecordSource`/`RawRecord` seam, the `ReadFilter` iterator.
  - [`read_input.md`](impl_plan/read_input.md) — step 1's input edge (`read/input/`): the
    validate-on-open gate, the BAM/CRAM region queries, the order guard, and the k-file
    merge. Covers both `alignment_file` and `sample_reads`.
  - [`typed_regions.md`](impl_plan/typed_regions.md) — step 3: the catalog rebase/knobs,
    the windowed substrate, and the `region_typing.rs` walk (resident → windowed).
  - [`typed_regions_cli.md`](impl_plan/typed_regions_cli.md) — the `pop_var_caller_exp` /
    `type-regions` build order: binary skeleton, `--min-copies` parser, output writer, and
    the `run_typed_regions` driver.
  - [`reference_info.md`](impl_plan/reference_info.md) — the reference-info reader: types →
    `.fai` reader → FASTA pass (the heart) → writer → cache → the two entry points.

This mirrors the repo-wide `doc/devel/{specs,architecture,implementation_plans}`
convention but scoped to ng, so the growing set of ng docs stays together.
