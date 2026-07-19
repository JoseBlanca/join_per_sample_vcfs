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
    the reference and cuts it into `TypedRegion`s (SsrLocus / SsrBundle / Generic / Satellite).
    First integrating spec — stands on `RefSeq`, the tandem-repeat scanner, and the STR catalog.
  - [`ref_seq.md`](spec/ref_seq.md) — the `RefSeq` reference-sequence accessor
    (foundational infra: resident + streaming + in-memory impls). Read filtering #8 and
    the pileup depend on it.
  - [`reference_info.md`](spec/reference_info.md) — reading a reference's info (contig table +
    content digest + reconstructed index) from a FASTA (optionally checked against a `.fai`) or
    a `.fai` alone; with the MD5s, the fasta↔fai check, a caller-held cache, and a `.fai` writer.
    Builds the `ContigList` every `RefSeq` impl takes but none can build.
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
  - [`reference_info.md`](arch/reference_info.md) — the reference-info reader's types &
    interfaces (`ReferenceInfo`/`ContigInfo`, the cache, the writer, the background verify);
    companion to `spec/reference_info.md`. Foundational infra, not a step.
- **`impl_plan/`** — step-by-step implementation plans (build order, not new design).
  - [`foundations.md`](impl_plan/foundations.md) — the first ng code: skeleton,
    `types.rs` seed, and the `RefSeq` accessor (three impls).
  - [`read_filtering.md`](impl_plan/read_filtering.md) — step 1: the `read/` module,
    the cascade, the `RecordSource`/`RawRecord` seam, the `ReadFilter` iterator.
  - [`typed_regions.md`](impl_plan/typed_regions.md) — step 3: the catalog rebase/knobs,
    the windowed substrate, and the `region_typing.rs` walk (resident → windowed).
  - [`reference_info.md`](impl_plan/reference_info.md) — the reference-info reader: types →
    `.fai` reader → FASTA pass (the heart) → writer → cache → the two entry points.

This mirrors the repo-wide `doc/devel/{specs,architecture,implementation_plans}`
convention but scoped to ng, so the growing set of ng docs stays together.
