# ng — module layout

*Status: architecture draft (2026-07-10), companion to
[`ng_step_interfaces.md`](ng_step_interfaces.md). That doc defines the shared **types
and step traits**; this one is their **physical home** — how the ng caller is laid out
as Rust modules, and the rules that keep the layout serving the lab (§3 of the
[spec](../spec/ng_proposal.md)) rather than fighting it.*

## Where ng lives

ng is a **new module inside `pop_var_caller`**, not a separate crate: `src/ng/`. That
keeps it able to **reuse the existing code directly** (the filters, `ReadLikelihoodModel`,
`GenotypeEmModel`, …; see the reconciliation table in `ng_step_interfaces.md` §6) and
means the eventual port-back of winning steps into the production engine is an
in-crate move, not a dependency dance.

## The tree

```
src/ng/
├── mod.rs            – module declarations + re-exports only (kept minimal)
├── types.rs          – the shared vocabulary, one file to start: the domain newtypes
│                       (Bp, SsrPeriod, LogProb, ids, ReadWeight…) plus the wire structs
│                       (GenomeRegion, LocusKind, Genotype, ModelParams, …). A deliberate
│                       *temporary* catch-all: splits into concept modules (units, locus,
│                       genotype, params) as clusters grow — see principle 3.
│
│   # one module per pipeline step — each owns its trait + its swappable impls + tests
├── read_admission/   – step 1  (a fixed prelude; wraps the existing bam/alignment_input filters)
├── read_prep/        – step 2  ReadPrep + impls (trust_mapper, reassemble, pair_hmm)
├── locus_router/     – step 3  LocusRouter + LocusSource + impls (catalog, active_region)
├── pre_pass/         – step 4  Caller, SampleSummarizer, CohortEstimator + impls
├── allele_candidates/ – step 6  CandidateGenerator + impls (rung_ladder [STR], assembly [generic])
├── likelihood/       – step 7  ReadLikelihood + impls (stutter models, pair-HMM)
├── genotype_prior/   – step 8  GenotypePrior + impls (flat, dirichlet, sfs, marginalized)
├── inference/        – step 9  GenotypeModel + impls (ml_grid, em_af, joint)
├── phasing/          – step 10 (physical / SNP-based)
├── locus_filter/     – step 11  the locus-level filters (add more as needed):
│                        · hidden_dup.rs  – ArtifactFilter (11a: paralog / hidden-duplication)
│                        · emission/      – EmissionModel (11b: heuristic, bic, freebayes)
├── allele_representation/ – step 12 AlleleRepresentation
├── quality/          – step 13 QualityModel
│
├── pipeline.rs       – the CallerRecipe + the driver that runs it end-to-end (single-phase)
└── bench/            – the standards harness: gold / silver / synthetic scoring
```

(Steps 5 — STR read-class / spanning — and the read-class machinery live inside
`locus_router/` or `allele_candidates/`; see *Open items*.)

## Organizing principles

**1. One module per pipeline step — trait, impls, and tests together.** This is the
load-bearing rule. Each step folder holds its trait *and* every competing implementation
*and* their tests, so a step's alternatives sit **side by side**. That is exactly what
the bake-off needs — "swap `allele_candidates::assembly` for `allele_candidates::rung_ladder`,
hold the rest, re-measure." It also satisfies the naming rule ([naming.md](../../../../ai/skills/rust-code-review/code_review/naming.md)):
modules are named for the **concept** they own (`allele_candidates`, `likelihood`, `genotype_prior`),
never for a layer or pattern (`models`, `services`, `common`, `utils`). *One step may host
sub-modules when the spec frames it as one step with several questions* — `locus_filter/`
(step 11) holds `hidden_dup` (artifact) and `emission` (emit decision) side by side, each
with its own trait and impls.

**2. STR-ness is not a separate subtree.** An STR candidate generator is just
`allele_candidates/rung_ladder.rs` sitting next to the generic `allele_candidates/assembly.rs`; the
router (step 3) decides which runs per locus. We do **not** split the pipeline into
`ssr/` vs `generic/` — that would scatter each step's variants across two trees and make
the per-step comparison awkward. STR-ness is a property of certain *implementations*, not
a top-level division. (STR domain *types* — `SsrMotif`, `SsrPeriod`, `SsrLocus` — still
carry the `Ssr` prefix and live in `units`/`locus` beside their generic peers.)

**3. Shared vocabulary starts in one `types.rs`, splits by concept as it grows.** The
cross-step types (the scalar newtypes + the wire structs) begin together in `types.rs`
rather than pre-split into many small files — the idiomatic Rust rhythm is to let a file
grow and extract a coherent chunk *when it tells you to*. `types.rs` is a common Rust
convention and an honest name for a mixed starting file. naming.md discourages a
*permanent* generic catch-all; we satisfy that by **splitting into concept modules** —
`units` (scalars), `locus`, `genotype`, `params` — as each cluster grows. Start minimal,
end concept-named.

**4. `bench/` is first-class.** The lab's whole purpose is *measuring* (spec §2), so the
standards harness — the representation-normalising truth comparison, the gold/silver/
synthetic scorers — is a real module the pipeline is built around, not a `tests/`
afterthought. A step's winner is decided here.

**5. Reuse over rewrite.** New modules, standing on existing code: `read_admission` wraps
the filters in [bam/alignment_input.rs](../../../../src/bam/alignment_input.rs),
`likelihood` builds on [ssr/cohort/read_model/](../../../../src/ssr/cohort/read_model/),
`inference` on [var_calling/posterior_engine.rs](../../../../src/var_calling/posterior_engine.rs).
The §6 reconciliation table is the map of what to reuse.

## Anatomy of a step module

Every step folder follows the same shape — the trait in `mod.rs`, one file per
implementation, tests beside the code:

```
allele_candidates/
├── mod.rs          – the CandidateGenerator trait + re-exports of the impls
├── rung_ladder.rs  – the STR repeat-length implementation  (+ #[cfg(test)] tests)
├── assembly.rs     – the generic local-assembly implementation (+ tests)
└── …               – further impls as the bake-off grows
```

New implementation of a step = a new file in that step's folder implementing the trait.
Nothing else in the tree changes — that locality is the point.

## How it assembles

- **`CallerRecipe`** (in `pipeline.rs`; defined in `ng_step_interfaces.md` §4) names one
  impl per step — it *is* one experiment ("freebayes candidates + HipSTR stutter
  likelihood + our cohort prior").
- **`pipeline.rs`** drives a recipe over its inputs, single-phase, in memory.
- **`bench/`** scores the pipeline's output against the standards and reports the frontier.

So: the step folders provide the *parts*, the recipe *selects* a set, the pipeline *runs*
it, and bench *judges* it. Swapping one part and re-running is the unit of work.

## Crate boundary and the port-back

ng stays a single-phase module inside `pop_var_caller` (spec §3): no `.psp` split, one
thread, reuse freely. The module tree here is the
*research* home; the production modules remain the *scaling* home.

## Naming to confirm

- `read_prep` (step 2) — "prep" is a mild abbreviation; alternatives `realignment` or
  `read_preparation`. The trait is `ReadPrep`, so `read_prep` matches it.
- `types.rs` (the one shared-types file) — a common Rust convention, honest for a mixed
  starting file; naming.md leans against a *permanent* generic module, so the plan is to
  split it into concept modules (`units`/`locus`/`genotype`/`params`) as it grows
  (principle 3). If you'd rather avoid `types.rs` entirely, start with `units` + concept
  files instead.

## Open items

- **Where step 5 (STR read-class / spanning) lives** — it is STR-only and feeds candidate
  generation; likely a submodule of `allele_candidates/` or `locus_router/`, not its own top-level step
  folder. Decide when the STR path is built.
- **Feature-gating.** If ng grows heavy, gate it behind a `cargo` feature so the production
  build need not compile the lab. Decide once there is code to gate.
- **`bench/` vs the existing `benchmarks/` tree.** `benchmarks/` holds data + scripts; the
  ng `bench/` module is in-crate scoring code. Keep the boundary explicit so they don't
  blur.
