# SSR/STR genotyping — architecture

**Status:** first draft, 2026-06-11 — a discussion starter, not a settled
design. Companion to the
[`ssr_genotyping.md`](../specs/ssr_genotyping.md) specification: the spec
says *what* we build and *why* (the model, the statistics, the trade-offs);
this document says *how the code is structured* to build it — crates,
modules, data flow, and what existing plumbing it reuses. Where the two
disagree on intent, the spec wins; where they disagree on code layout, this
document wins.

Follows [`design_principles.md`](../specs/design_principles.md). Modelled on
the SNP caller's
[`calling_pipeline_architecture.md`](../specs/calling_pipeline_architecture.md),
which is the closest existing analog (a two-stage extract→genotype pipeline
over the same `.psp` plumbing).

> **What this document is — and where detail lives.** This is the **overall
> architecture**: the component view, the load-bearing structural decisions
> (E), the reuse map, and the one refactor of *existing shared code* the SSR
> caller forces — the `.psp` container generalisation (§10, here because it
> touches the SNP path, not because it's an SSR module). That work is now
> **settled**, and this doc has reached its job.
>
> **Per-module detail does *not* keep growing this document.** Detail written
> far ahead of code rots — a type decided in prose loses to the first compile.
> So each SSR module is detailed **just-in-time in its own document /
> implementation plan, immediately before it is built** (design → build →
> next), in data-flow order. The first such module is the **shared types**
> (`src/ssr/types.rs`), in a dedicated document. This doc stays the stable
> overall reference the per-module docs point back to. **Open points are
> flagged inline** — proposals to argue with, not conclusions.

---

## 1. Scope of this document

In: the executable/crate decomposition, the module tree, the data-flow
wiring between stages, the map of reused vs. net-new code, and the
cross-cutting concerns (parallelism, memory, errors, determinism).

Out: the statistical model (spec §5), the on-disk byte formats (spec §3.2,
§4.3 — those are format specs, owned by the spec and by a future
`per_sample_ssr_format.md` mirroring
[`per_sample_pileup_format.md`](../specs/per_sample_pileup_format.md)), and
parameter defaults (spec §6). This document references those; it does not
restate them.

**Naming.** We say **SSR** (not the broader "TR") deliberately — our scope is
period ≤ 6 / short-read / pop-gen, the domain where SSR/microsatellite is the
precise and customary term; "TR" is the long-read/expansion umbrella we don't
target. **TRF** is the *detector tool* (Stage 0), never a feature word. Full
rationale in the spec glossary's *SSR vs STR vs TR* note.

---

## 2. The big picture

**Three** CLI subcommands over three persistent artifacts, exactly mirroring
the SNP path's `pileup → .psp → var-calling` shape (the simulator is *not* a
subcommand — it is crate/test code, §2.1):

```
                          ┌──────────────────────────────────────────────┐
  reference FASTA ───────►│  ssr-catalog   (Stage 0)                      │──► catalog
                          └──────────────────────────────────────────────┘    .ssr_catalog.bed.gz
                                                                                    │
                                       ┌────────────────────────────────────────────┘
                                       ▼
  BAM/CRAM (sample 1) ─┐    ┌──────────────────────────────────────────────┐
  BAM/CRAM (sample 2) ─┼───►│  ssr-pileup    (Stage 1, per sample, ∥)       │──► evidence
        ...            │    │  reads + catalog → per-read Qᵣ(L)              │    sampleN.ssr.psp
  BAM/CRAM (sample N) ─┘    └──────────────────────────────────────────────┘    (one per sample)
                                                                                    │
                                       ┌────────────────────────────────────────────┘
                                       ▼
                          ┌──────────────────────────────────────────────┐
                          │  ssr-call      (Stage 2, cohort, one process) │──► VCF
                          │  N evidence + catalog → joint EM → posteriors  │    (GangSTR-compatible)
                          └──────────────────────────────────────────────┘

  (test/dev only, no CLI) ssr simulator: genotypes + stutter → synthetic BAM
                          and/or evidence + truth table — a crate module the
                          test suite calls directly (§2.1)
```

The names mirror the SNP caller's `pileup` (per-sample evidence) and
`var-calling` (cohort calls) one-for-one — `ssr-pileup` gathers per-sample
read evidence, `ssr-call` makes the cohort-wide genotype calls — so the
two-callers parallel is legible at the command line.

Two framings from the spec drive every structural choice below:

- **Two independent callers, one BAM** (spec §1.1). The SSR pipeline shares
  only *low-level alignment I/O and the vendored sdust* with the SNP caller —
  not its records, not its math, not its VCF. Architecturally this means the
  reuse boundary is drawn at the **plumbing layer** (I/O, container,
  numerical kernels), never at the **caller layer** (records, mergers,
  variant semantics). Keeping that line clean is the single most important
  structural constraint here.
- **Two-stage extract→genotype** (spec §1.2): heavy per-sample work once,
  summarised to a columnar artifact, then a light cohort math stage. This is
  the same memory-for-scaling bet the SNP path makes, so the same plumbing
  (`.psp` blocks, `--block-window-bp` decode unit, coordinate-merge scan)
  applies — and is reused, not reinvented.

### 2.1 The simulator is crate/test code, not a subcommand

**Decided 2026-06-11.** The STR-aware simulator (spec §7) is **not** a shipped
CLI subcommand. It is a **library module in the crate** (`src/ssr/simulate/`)
that the test suite calls directly — unit, property, and statistical tests
construct synthetic evidence (or synthetic reads) and a truth table in-process,
then assert the caller recovers the injected genotypes/`π`/`θ`/`F`. It is
test/dev scaffolding, not part of the production pipeline a user runs, so it
does not belong on the production CLI.

- **Why a module, not a subcommand:** tests need a *programmatic* API (build a
  cohort, perturb the stutter model, read back the truth), which a CLI would
  only wrap awkwardly. Keeping it library-shaped also keeps the shipped CLI
  surface to exactly the three pipeline stages a user actually invokes.
- **Reachability:** it is ordinary (non-`#[cfg(test)]`) crate code so both
  in-crate `#[cfg(test)]` tests **and** `tests/` integration tests can use it;
  gating behind a `sim` cargo feature is an option if we want it excluded from
  release builds (decide when building it).
- **If a command-line generator is ever wanted** — e.g. to materialise a
  reusable benchmark BAM once and commit it — add a dev-only
  [`examples/`](../../examples/) binary or a hidden dev subcommand then, driven
  by a concrete need. Not now, and not on the production surface.
- **Anti-tautology boundary (spec §7) still holds:** the simulator's model
  definition stays separate from the caller's code regardless of where it
  lives — a module in the same crate is fine as long as it does not *import the
  caller's* `(u,d,ρ)` to generate the data it then tests the caller against.

---

## 3. Decision E — the load-bearing structural decision

The spec (§8, §11) named three coupled structural questions collectively **E**:
repo/crate placement, the shared `.psp` container, and CLI naming. They had to
be settled first because every module path below depends on them. **All three
are now decided (2026-06-11)** — recorded below with the rationale and options
considered, so the trail survives.

### 3.1 Repo / crate placement — **DECIDED: C (same crate, new module tree)**

**Settled 2026-06-11.** The SSR caller is built as a new `src/ssr/` module
tree with `ssr-*` subcommands on the existing binary — no new repo, no
workspace split. Rationale below; the §"escape hatch" records the cheap path
to a later split if compile times or independent release ever demand it.

**Options considered:**

| option | what it means | cost |
|---|---|---|
| **A. New repo** | SSR caller is its own project | re-vendoring or cross-repo-crate-splitting noodles, sdust, the `.psp` plumbing, posterior engine; two CI/release tracks; the "shared only raw alignments" boundary becomes a published API |
| **B. New crate(s), same workspace** | split this single crate into a workspace; SSR is a sibling crate | introduces a workspace (currently one crate); forces the shared plumbing into its own crate *now*, up front |
| **C. Same crate, new module tree + subcommands** (recommended) | `ssr-*` subcommands on the existing binary, code under a new top-level `src/ssr/` tree | none structural; reuse is direct `use` of existing modules; shared plumbing extracted *in place* as it's needed |

**Recommendation: C.** The reuse surface (§5) is broad and lives *inside*
this crate — the `.psp` container, the posterior engine's IBD prior, the BAQ
banded-forward kernel, the FASTA/BAM readers, sdust. A new repo (A) or even a
workspace split (B) turns every one of those into a cross-crate API that must
be designed and stabilised before the SSR caller can compile a line. Option C
lets the SSR caller `use crate::psp::…` today and lets the shared/private line
be drawn *empirically*, by extraction from two live consumers (§3.2) — which
is exactly what the spec §11 asks for ("extracted *from the two concrete
consumers* as SSR is built, not designed up front").

The existing binary is already a multi-subcommand dispatcher (`pileup`,
`var-calling`, `estimate-contamination`, … under
[`src/pop_var_caller/`](../../src/pop_var_caller/)), so adding `ssr-catalog`
/ `ssr-pileup` / `ssr-call` is the established pattern, not a new one.

> **Escape hatch:** if SSR-side compile times or a desire to ship the SSR
> caller independently later forces a split, C → B is a mechanical move
> (promote `src/ssr/` and the extracted `psp` container to crates) and is
> *easier* after the shared boundary has been found by extraction than before.
> So C is also the lowest-regret path to B.

### 3.2 The shared `.psp` container — **DECIDED: generic core + two specializations**

**Settled 2026-06-11.** Reuse [`src/psp/`](../../src/psp/) with a **generic
columnar-block core** and **two thin specializations** — `snp` and `ssr` —
riding on it. Not a fork (two copies), not a registry flag bolted onto a
SNP-shaped container, but a clean separation where the core knows *bytes and
blocks* and the specializations know *columns and records*.

**This is a conceptual three-layer split, not (yet) a directory tree.** The
layering below is a rule about *what each file is allowed to know*, enforced by
naming, not by nesting. The physical layout **starts flat** in `src/psp/` and
only grows submodules if it gets crowded — see "Physical layout" below.

**Pre-alpha — no backwards compatibility.** We are not constrained to keep the
existing on-disk `.psp` layout readable. The container version can be revved,
the SNP schema restructured, the header reshaped — whatever makes the
core/specialization split clean. So the refactor is a *redesign of the SNP
container into its generic form*, not a careful byte-preserving lift. The
regression gate is therefore **"the SNP caller still produces correct calls"**
(its end-to-end tests pass), **not** "the SNP `.psp` is byte-identical to the
old format" — old files do not need to read, and byte-identity is no longer a
constraint we owe anyone. (Byte-identity across thread counts / runs within
*one* format version is still a useful determinism property — that's separate
from format compatibility across versions.)

**The three conceptual layers (roles, not directories):**

- **Generic core** — header skeleton, block grid, tail index, trailer, varint,
  per-column zstd framing, column-selective decode, the region-query seek.
  Knows the *shapes* a column can take (`Scalar` / `List` CSR-ragged / `Bytes`
  dict) but **nothing about which columns exist or what they mean**. The schema
  is data it's handed, never baked in. The `kind` header field selects the
  specialization (spec §4.3). This is the bulk of today's files — `block`,
  `header`, `index`, `trailer`, `varint`, `reader`, `writer`, `errors`, `mod`.
- **`snp` specialization** — today's `V1_0_COLUMNS` registry plus the SNP
  record ↔ columns mapping, expressed against the generic core. This is the
  existing [registry.rs](../../src/psp/registry.rs) content (renamed to mark it
  schema-specific, below).
- **`ssr` specialization** — the spec §4.3 column table (`int32` scalars,
  `hist_*`/`amb_*`/`offl_*` CSR lists, `chrom`/`offl_seqs` dict strings) plus
  the SSR locus-record ↔ columns mapping. New (~2 files).

The registry already models exactly the three shapes both schemas need, so the
SSR schema is **a registry table + a record mapping against the generic core**,
never a new container. The win of keeping the core genuinely generic (over a
registry-as-parameter bolt-on) is that the *core has no SNP assumptions left in
it* — which is what makes the SSR side a clean addition rather than a pile of
`if kind == ssr` branches. That win is about the *code boundary*, and is
independent of whether the files sit flat or in submodules.

**Physical layout — flat first, submodules on a trigger.** `src/psp/` has 11
files today; ~9 are the generic core, 1 is the SNP schema (`registry.rs`), and
SSR adds ~2 — landing around **13–14 files**. That is navigable flat,
especially since [mod.rs](../../src/psp/mod.rs) already documents the layout. So:

- **Start flat.** Keep the generic-core files as they are; add the SSR registry
  + record mapping as new files beside them.
- **Name schema-specific files to mark the boundary** — e.g. rename
  `registry.rs` → `registry_snp.rs` and add `registry_ssr.rs` (and any
  `record_snp.rs` / `record_ssr.rs` mappings), so a flat `ls` still reads as
  "generic core + `*_snp` + `*_ssr`." The naming *is* the enforcement of the
  conceptual split at this size.
- **Carve `core/` `snp/` `ssr/` submodules only on a trigger** — if
  schema-specific files exceed ~3 per schema, or the directory pushes past ~18
  files. That move is purely mechanical (a path reshuffle, like the var_calling
  reorg) and is cheaper done once the real shape is known than guessed now.

This is also the more faithful reading of spec §11 ("extracted *from the two
concrete consumers* as SSR is built, not designed up front"): we commit to the
code boundary on day one and defer the directory commitment until file count
earns it.

**The one structural difference the core must absorb (spec §4.3):** the block
index keys on **interval extent** (`last_end`), because SSR records are
intervals `[start, end)` whereas SNP records are points. Since the core is now
genuinely generic, point-keying becomes the degenerate case of interval-keying
(`last_end == last_start + 1`), so the core carries the interval form and the
SNP specialization is just the narrow case — not two index implementations.
This lives in [index.rs](../../src/psp/index.rs) and the region-query path, and
is the **first concrete extraction task** (§8 roadmap item 2).

> **Two separate "how much up front" questions — don't conflate them.**
> (1) *Code-boundary depth:* make the core genuinely schema-agnostic now (not a
> registry flag on a SNP-shaped container) — **yes, do this up front**; the
> pre-alpha freedom means no compatibility tax for restructuring, and a clean
> core is cheaper than two schemas fighting a half-generalised container later.
> (2) *Directory structure:* nest `core/`/`snp/`/`ssr/` submodules now — **no,
> defer**; the file count doesn't earn it yet, and the conceptual split is
> enforced by naming regardless. The first is about correctness of the
> abstraction; the second is just where files sit.

### 3.3 CLI naming — **DECIDED**

**Settled 2026-06-11.** Three subcommands, named to mirror the SNP caller's
roles one-for-one (the simulator is not a subcommand, §2.1):

- **`ssr-catalog`** — reference FASTA → catalog of loci (Stage 0); embeds the
  local reference into the catalog (`ref_seq`), so it is the only stage that
  reads the FASTA for the SSR algorithm.
- **`ssr-pileup`** — BAM/CRAM + catalog → per-sample evidence `.ssr.psp`
  (Stage 1). The SSR analog of the SNP `pileup`. No reference for the SSR math
  (it reads `ref_seq` from the catalog); a reference is needed only to *decode
  CRAM* input, not BAM (spec §3.2 caveat).
- **`ssr-call`** — N evidence files + catalog → cohort VCF (Stage 2). The SSR
  analog of the SNP `var-calling`.

The `ssr-` prefix keeps them from colliding with the SNP subcommands and makes
the two-callers split legible at the command line. Shared flags
(`--reference`, `--regions`, `--threads`, `--block-window-bp`) reuse the
existing [shared_args.rs](../../src/pop_var_caller/cli/shared_args.rs) parsers.
(Earlier draft names `ssr-extract` / `ssr-genotype` were renamed for clarity —
"pileup/call" say what the stage *does* and match the SNP vocabulary.)

---

## 4. Module tree (proposed)

Under decision **C**, a new top-level `src/ssr/` tree, one sub-module per
stage plus shared SSR types, with the CLI subcommands wired in beside the
existing ones. **Sketch — names and boundaries are first-draft:**

```
src/
├── psp/                     # EXISTING — stays FLAT; conceptual 3-layer split by naming (§3.2)
│   ├── block,header,index,  #   GENERIC CORE (~9 files): knows shapes, NOT columns.
│   │   trailer,varint,       #     index.rs gains interval-extent keying (point = degenerate case)
│   │   reader,writer,        #
│   │   errors,mod            #
│   ├── registry_snp.rs      #   snp schema: V1_0_COLUMNS + SNP record↔columns (renamed from registry.rs)
│   └── registry_ssr.rs      #   NEW ssr schema: spec §4.3 columns + locus-record↔columns
│                            #   → carve core/ snp/ ssr/ submodules later IF it crowds (§3.2 trigger)
│
├── ssr/                     # NEW — the SSR caller, all stages
│   ├── mod.rs
│   ├── types.rs             #   ★ FOUNDATIONAL, built first — shared domain types (Locus,
│   │                        #   Motif, allele repr/key, candidate set). Doc: ssr_shared_types.md
│   ├── catalog/             # Stage 0: TRF wrapper + post-process + catalog I/O
│   │   ├── trf.rs           #   run/parse TRF — shell-out OR vendor/port, NO FFI (§6)
│   │   ├── postprocess.rs   #   period≤6, purity/score filter, merge, split compounds, mappability
│   │   └── format.rs        #   self-describing bgzip+tabix BED-like read/write
│   ├── pileup/              # Stage 1 (ssr-pileup): per-sample reads → .ssr.psp
│   │   ├── read_handling.rs #   pull/anchor reads, soft-clip recovery, spanning test (spec §4.1)
│   │   ├── fast_path.rs     #   flank-anchored exact motif count (spec §4.2 fast)
│   │   ├── pair_hmm.rs      #   NEW bespoke banded pair-HMM forward (spec §4.2 slow)
│   │   ├── ladder.rs        #   on-/off-ladder candidate construction + normalization
│   │   └── locus_record.rs  #   in-memory per-locus aggregation; persisted via psp ssr schema (§3.2)
│   ├── cohort/              # Stage 2 (ssr-call): cohort EM → VCF
│   │   ├── candidate_set.rs #   A_ℓ assembly: cohort-union peak+adjacent rule (spec §5.1)
│   │   ├── stutter.rs       #   covariate stutter kernel S_θ + M-step (spec §5.2)
│   │   ├── seed.rs          #   confident-homozygote pre-pass + per-cell cutoff (spec §5.4)
│   │   ├── em.rs            #   joint deconvolution EM loop (π, θ; latent G)
│   │   ├── prior.rs         #   thin adapter onto the SNP IBD-mixture prior (spec §5.3)
│   │   ├── base_measure.rs  #   Dirichlet base measure G₀ (spec §5.5)
│   │   └── vcf_out.rs       #   GangSTR-compatible VCF (spec §5.9)
│   └── simulate/            # test/dev only, NOT a subcommand (§2.1): evidence- + read/BAM-level
│
└── pop_var_caller/
    └── cli/                 # EXISTING — add ssr-catalog/ssr-pileup/ssr-call dispatch here
```

This is a placeholder partition, not a contract — each sub-tree is fixed in its
own per-module pass (§8), `types.rs` first. The intent it encodes:
**`types.rs` is the shared spine, built first**; **Stage 0/1/2 are separate
module families** (separate executables, separate failure modes); **the net-new
numerical machinery is named explicitly** (`pair_hmm`, `stutter`, `em`,
`candidate_set`, `seed`); and **the reuse adapters are thin and named**
(`prior.rs`, `locus_record.rs`).

---

## 5. Reuse map — what the SSR caller rides on

The reuse boundary is the plumbing layer (§2). Concretely:

| existing component | path | how SSR uses it |
|---|---|---|
| **`.psp` container** | [`src/psp/`](../../src/psp/) | the evidence file *is* this container — generic core shared, a new SSR schema (`registry_ssr`) added (§3.2). Direct, central reuse. |
| **IBD-mixture genotype prior** | [posterior_engine.rs](../../src/var_calling/posterior_engine.rs) | spec §5.3: reuse the multiallelic `F`-adjusted prior **verbatim** — it already enumerates non-decreasing ploidy-tuples with the `F·π_i + (1−F)·π_i^ploidy` homozygous term. SSR feeds it a repeat-allele set instead of base alleles. Adapter only (`cohort/prior.rs`). |
| **BAQ banded forward** | [probaln.rs](../../src/baq/probaln.rs), [scratch.rs](../../src/baq/scratch.rs) | spec §4.2: the slow-path pair-HMM **learns the banded-forward pattern** (and the scratch-buffer discipline) but does **not** couple to it — it's a *bespoke* 3-state forward. Pattern reuse, not code reuse; the scratch-buffer ethos transfers directly. |
| **Indel left-align kernel** | [indel_norm.rs](../../src/pileup/walker/indel_norm.rs) (`normalize_alleles`) | spec §4.2: off-ladder canonicalization (§4 of the types doc) **reuses the actual GATK-ported kernel**, not a re-implementation — it's already a representation-neutral `(seqs, bounds)` primitive. Lift it to a neutral shared module (two real users now); SNP wrapper + SSR adapter both call it. **Code reuse**, unlike BAQ. |
| **BAM/CRAM read I/O** | noodles via [`src/bam/`](../../src/bam/) | spec §4.1: pull reads overlapping a locus + flank. Shared low-level I/O — the sanctioned cross-caller reuse. |
| **FASTA reference reader** | [`src/fasta/`](../../src/fasta/) | **Stage 0 only** for the SSR algorithm — scans the reference to build the catalog and embed `ref_seq` (§3.2). Stages 1–2 read reference bases from the catalog, not the FASTA (CRAM decoding aside). Same reader the pileup walker shares. |
| **vendored sdust** | (vendored) | Stage 0 optional TRF prefilter (spec §3.1) — same masker the SNP DUST filter uses. |
| **VCF writer** | [`src/vcf/`](../../src/vcf/) | Stage 2 output; SSR writes its own records/headers (GangSTR-compatible) but reuses the writer plumbing. |
| **`--regions` / block index** | [regions.rs](../../src/regions.rs), [index.rs](../../src/psp/index.rs) | region-restricted `ssr-pileup`/`ssr-call`; the interval-keying change (§3.2) is the one delta. |
| **CLI shared args** | [shared_args.rs](../../src/pop_var_caller/cli/shared_args.rs) | `--reference`, `--threads`, `--block-window-bp`, `--regions`. |

**Net-new machinery with no analog** (must be built; these are where the
risk and the work are): the **banded pair-HMM forward** (`pileup/pair_hmm.rs`),
the **covariate stutter kernel + its EM M-step** (`cohort/stutter.rs`), the
**candidate-set peak+adjacent rule** (`cohort/candidate_set.rs`), the
**confident-homozygote seed pre-pass** (`cohort/seed.rs`), the
**integer-partition genotype enumeration over a repeat-allele set**
(ConSTRain-style, feeding the reused prior), the **TRF wrapper +
post-process** (`catalog/`), and the **STR-aware simulator** (`simulate/`,
test/dev only — §2.1).

The split is the headline: **Stage 2's *prior* is reused, but its
*likelihood* (stutter convolution) and *candidate assembly* and *EM* are all
new** — because the SNP and SSR generative models diverge precisely at the
likelihood (spec §1.1).

---

## 6. Data-flow walkthrough (architectural, not statistical)

Each stage as a wiring diagram: what it reads, the module that does the work,
what it writes. The math is the spec's; this is the plumbing.

### Stage 0 — `ssr-catalog`
- **In:** reference FASTA (+ optional sdust mask).
- **Work:** `catalog/trf.rs` runs/parses TRF over the genome (or sdust
  windows); `catalog/postprocess.rs` applies period≤6, purity/score filters,
  overlap merge, compound split, mappability drop.
- **Out:** one self-describing bgzip+tabix BED-like TSV (`catalog/format.rs`),
  `##`-header carrying reference md5 + TRF params (spec §3.2).
- **Shape:** batch, single process, embarrassingly parallel by contig.
- **TRF integration — no FFI.** TRF is an external C tool, not a crate. **FFI is
  ruled out.** The choice is **shell-out to a `trf` binary** (parse its output;
  a runtime/container dependency) **vs. vendor-and-port** the detection into
  Rust. Settled at the **top of the Stage 0 pass**, *informed by reading the
  vendored TRF / HipSTR / GangSTR sources* (how they invoke or reimplement
  detection) — not pre-decided here.

### Stage 1 — `ssr-pileup` (per sample, parallel processes)
- **In:** one BAM/CRAM + catalog (the catalog's `ref_seq` supplies all
  reference bases; a FASTA is needed only to decode CRAM, spec §3.2).
- **Work, per locus:** `pileup/read_handling.rs` pulls + anchors reads and
  runs soft-clip recovery / the spanning test (spec §4.1);
  `pileup/fast_path.rs` counts the clean majority; `pileup/pair_hmm.rs`
  realigns the ambiguous/soft-clip-recovered/off-ladder minority;
  `pileup/ladder.rs` builds the on-/off-ladder candidates and normalizes
  off-ladder sequences. Results aggregate into a per-locus record (confident
  reads → histogram; ambiguous → sparse CSR; off-ladder → `offl_*`).
- **Out:** `sampleN.ssr.psp` via `pileup/locus_record.rs` → the `psp` SSR
  schema (`registry_ssr`) + generic container writer. Sparse (no-coverage loci
  absent), block-gridded.
- **Shape:** one process per sample, run in parallel à la the SNP pileup
  (memory:
  [benchmark-each-caller-with-its-native-parallelism]); within a process,
  parallel over catalog loci / blocks. **All stored likelihoods stutter-free**
  (spec §4.2) — the architectural reason Stage 1 needs no cohort context and
  is perfectly parallel.

### Stage 2 — `ssr-call` (cohort, one process)
- **In:** N `.ssr.psp` files + catalog.
- **Work:** a **coordinate-merge synchronized scan** over the N files
  (block-aligned by the shared window grid, spec §4.3) yields, per locus, the
  cohort-aggregate evidence. Then per locus: `cohort/candidate_set.rs`
  assembles `A_ℓ` (peak+adjacent + off-ladder union); `cohort/seed.rs` seeds
  the kernel from confident homozygotes; `cohort/em.rs` runs the joint
  deconvolution EM, calling `cohort/stutter.rs` (kernel + convolution),
  `cohort/prior.rs` (reused IBD prior), and `cohort/base_measure.rs` (G₀);
  `cohort/vcf_out.rs` emits the call.
- **Out:** one GangSTR-compatible VCF (spec §5.9).
- **Shape:** one process; **parallel across loci** (each locus is an
  independent EM), with the stutter kernel pooled across loci within a
  covariate cell — so the cross-locus pooling is a reduce step, not a
  per-locus dependency. This is the one place the architecture must reconcile
  "loci are independent" with "the kernel pools across loci"; the seed
  pre-pass (a cohort-wide reduce *before* the per-locus EM) is how the spec
  resolves it (§5.4), and the module layout mirrors that: `seed.rs` runs
  first and globally, `em.rs` runs per-locus after.

### Simulator (`simulate/`, test/dev only — not a subcommand, §2.1)
- A crate module the test suite calls directly. Emits at two levels (spec §7):
  **evidence-level** synthetic `.ssr.psp` (to test Stage 2 before Stage 1
  exists — the critical-path order) and **read/BAM-level** synthetic BAM +
  truth table (to test the whole pipeline incl. the pair-HMM). Anti-tautology
  rule: its model definition is kept separate from the caller's code.

---

## 7. Cross-cutting concerns

- **Parallelism model — adapted to the SSR situation, *not* inherited from the
  SNP path.** The natural decomposition is different: Stage 0 by contig; Stage 1
  by sample (process) × by locus (within); Stage 2 by independent locus with a
  global seed/kernel reduce. That shape is **more embarrassingly parallel** than
  the SNP cohort path (loci are independent units; no per-position streaming
  pipeline), so the SNP thread-budget/pool machinery is **not assumed to carry
  over** — it was built for a different problem. The concrete parallelism for
  each stage is a **per-stage design decision settled in that stage's pass**
  (Stage 1 / Stage 2), against the real work shape, not pre-defaulted here. Open
  item, tracked in §9.
- **Memory.** The `--block-window-bp` decode unit is the SSR evidence's
  primary RSS lever exactly as on the SNP path
  ([block-window memory lever]) — it controls Stage 2's merge-scan working
  set. Scratch-buffer-over-fresh-allocation discipline
  ([scratch buffers]) applies to the pair-HMM (per-read forward matrices) and
  the EM (per-locus responsibility tables) — the two hot loops.
- **Self-describing artifacts.** Catalog `##` header, evidence TOML header
  (`kind="ssr"`, `reference_md5`, `catalog_md5`, ploidy, extraction params),
  VCF `##` headers — no sidecars (spec §2). md5 binding (reference + catalog)
  is checked at every stage boundary, reusing the SNP path's md5 convention
  ([common.rs](../../src/pop_var_caller/common.rs) `FastaVerify`).
- **Errors.** Typed per-stage error enums in the [errors.rs](../../src/psp/errors.rs)
  house style (writer-side vs reader-side split), bridged to the CLI via the
  existing [error_bridge.rs](../../src/pop_var_caller/cli/error_bridge.rs).
- **Determinism.** The SNP path treats byte-identical output as a
  regression gate; the SSR Stage 1 (stutter-free, no cohort context) should
  hold the same bar. Stage 2's EM is iterative — determinism there means
  fixed iteration order + a pinned convergence rule, not byte-identity across
  thread counts (a property to specify in the Stage 2 plan).

---

## 8. Roadmap

**Settled here (the overall architecture):**

- ✅ **Decision E** (§3) — placement (same crate, `src/ssr/`), the `.psp`
  container split (generic core + `snp`/`ssr` schemas, flat), CLI names
  (`ssr-catalog`/`ssr-pileup`/`ssr-call`; simulator not a subcommand).
- ✅ **Container generalisation + SSR schema** (§10) — the one refactor of
  *existing shared code*; detailed and ready to spin into a
  `per_sample_ssr_format.md` format spec + a container-refactor implementation
  plan. Build-ready independent of the SSR modules.

**From here — per-module, just-in-time (design → build → next), in data-flow
order.** Each module is detailed in **its own document / implementation plan**
immediately before it is built, *not* by growing this doc (see the header
note). Each carries Bucket-1 synthetic tests on the critical path (spec §7/§11):

0. **`types.rs` — shared domain model** *(first; foundational; drafted in
   [ssr_shared_types.md](ssr_shared_types.md)).* The allele representation
   (sequence-identity, on-/off-ladder hybrid, normalized off-ladder key),
   `Locus`, `Motif`, the candidate-set type. The spine every stage shares —
   settled before the stages that use it.
1. **Stage 0 — `ssr-catalog`.** TRF integration strategy (shell-out vs
   vendor/port, **no FFI** — §6, decided at this pass's start from the vendored
   tools), the sdust-prefilter decision harness, the post-process filter
   pipeline, catalog I/O types.
2. **Stage 1 — `ssr-pileup`.** The banded pair-HMM type design (bands, scratch,
   emission model), the fast/slow dispatch, the ladder/normalization types, the
   locus record → `ssr` columns mapping; **the stage's parallelism** (§7).
3. **Stage 2 — `ssr-call`.** The EM loop structure, the stutter-kernel
   representation + its two open modelling cores (spec §5.2 Option 1/2; §5.4
   `auto` cutoff), the candidate-set + genotype-enumeration types, the prior
   adapter, the VCF mapping; **the stage's parallelism** (§7).
4. **Simulator.** The two emission levels and the anti-tautology boundary.

Implementation plans land under
[`doc/devel/implementation_plans/`](../implementation_plans/).

---

## 9. Open questions — deferred to the pass that owns them

Not blockers for the overall architecture; each is settled inside the
data-flow pass named:

1. **Per-stage parallelism** (§7, *Stage 1 / Stage 2 passes*) — the concrete
   decomposition is adapted to the SSR work shape, **not** inherited from the
   SNP thread machinery. Includes whether a hypervariable cohort ever forces
   the Stage-2 merge-scan to chunk (as the SNP cohort path eventually needed).
2. **TRF integration** (§6, *Stage 0 pass*) — shell-out vs vendor/port (no
   FFI), decided from the vendored TRF/HipSTR/GangSTR sources.
3. **Where the pair-HMM lives** (*Stage 1 pass*) — `src/ssr/pileup/pair_hmm.rs`
   (SSR-private, recommended) vs. a shared `src/align/` next to `baq/` (if a
   third user ever appears). Kept SSR-private until a second user exists.
4. **Trait vs two builders for the container schema** (§10.7, *container
   build*) — one `PspSchema` trait vs two concrete writer/reader builders;
   settle against real generics ergonomics when implementing.

---

## 10. Container generalisation + SSR schema (detailed here, by exception)

*This detail lives **in** the architecture doc — unlike the per-module SSR
detail, which goes in its own documents (header note) — because it is a refactor
of **existing shared code** ([`src/psp/`](../../src/psp/)), not a new SSR module.
It is the **code architecture** of making `src/psp/` host two schemas — not the
byte format (a future `per_sample_ssr_format.md`, mirroring
[`per_sample_pileup_format.md`](../specs/per_sample_pileup_format.md)) and not
the step-by-step build (a future implementation plan). Grounded in a read of the
current module, 2026-06-11.*

### 10.1 Current state — what's already generic, what's SNP-coupled

A read of the module sorts its files into three buckets, which is what makes
the "generic core + two schemas" split (§3.2) cheap rather than a rewrite:

- **Already generic — reuse verbatim, no change:**
  [block.rs](../../src/psp/block.rs) is pure column codecs
  (`encode_scalar_column<T>`, `decode_list_column_csr<T>`, the bytes/CSR
  encoders, the zstd helpers, `BlockHeader` / `ColumnManifestEntry`) — it knows
  *shapes* (`Scalar` / `List` / `Bytes` / varint) and wire *types*
  (`WireScalar`), never which columns exist. [varint.rs](../../src/psp/varint.rs)
  and [trailer.rs](../../src/psp/trailer.rs) likewise. The **columnar read
  layer** (`BlockColumnReader` / `BlockColumns` in
  [reader.rs](../../src/psp/reader.rs)) is also schema-agnostic at its core —
  it hands back decoded columns.
- **Generic with one semantic tweak — [index.rs](../../src/psp/index.rs):**
  `BlockIndexEntry` *already* stores a per-block range `{ chrom_id, first_pos,
  last_pos, block_offset }`. So the interval change (§3.2) is **not** a new
  structure — it is (a) what fills `last_pos` (for SSR, `max(record.end)` over
  the block, not the last record's start) and (b) the region-overlap test
  comparing the query against record intervals `[start, end)`. Smaller than it
  first looked.
- **The actual SNP coupling — the work concentrates here:**
  [writer.rs](../../src/psp/writer.rs) iterates `V1_0_COLUMNS` directly and
  takes `&PileupRecord`; [reader.rs](../../src/psp/reader.rs)'s *typed-record*
  path decodes via an **exhaustive `ColumnKey` match** into `PileupRecord`;
  [registry.rs](../../src/psp/registry.rs) is `V1_0_COLUMNS` + the `ColumnKey`
  enum. These three are where "schema" is currently hardcoded.
- **Header — [header.rs](../../src/psp/header.rs):** the file is *already*
  self-describing — its `[[column]]` array **is** the binary schema, and the
  reader cross-checks it (`cross_check_against_registry`). What's missing is a
  **`kind` schema-family tag** (spec §4.3): today the cross-check validates
  against the *hardcoded* `V1_0_COLUMNS`; it must instead validate against the
  registry the `kind` selects.

### 10.2 The generalisation pattern — "schema as data"

The generic core owns everything in the first two buckets: blocking, the block
index, per-column zstd, the header skeleton + `[[column]]` self-description,
and the columnar read/write of a *given* set of `ColumnDef`s. A **schema**
(`snp` or `ssr`) supplies exactly four things:

1. its **registry table** — a `&'static [ColumnDef]` (today's `V1_0_COLUMNS` is
   the `snp` one);
2. a **record type** (`PileupRecord`; new `SsrLocusRecord`);
3. **record → column accumulators** (the encode half of `write_record`);
4. **columns → record** (the decode half of the typed-record iterator).

Captured as a `PspSchema` trait (or, if generics get unwieldy, two parallel
concrete writer/reader builders — decide when implementing, §10.6). The writer
becomes generic over the schema instead of closing over `V1_0_COLUMNS` +
`PileupRecord`; the reader's typed path likewise.

**Key leverage: the SSR cohort stage rides the *columnar* layer, not the typed
one.** Stage 2's coordinate-merge scan (§6) wants columns, not materialised
records — and `BlockColumnReader` / `BlockColumns` are already schema-agnostic.
So the per-record SSR materialisation (item 4) can stay **thin**, and the
heavy cross-sample path reuses the generic columnar reader directly. This also
means the exhaustive-`ColumnKey`-match decode does **not** need to be
generalised into one match over both schemas — each schema keeps its own typed
decode (a closed enum per schema), and the shared substrate is the columnar
layer beneath it.

### 10.3 The `kind` tag + self-describing validation

Add `kind` to the TOML header (`"snp"` | `"ssr"`, spec §4.3). On read, `kind`
selects the registry; `cross_check_against_registry` is parameterised on *that*
registry rather than the `V1_0_COLUMNS` constant. The `[[column]]` array stays
the authoritative binary schema (the file remains self-describing even for a
reader that doesn't know the `kind`); `kind` is the fast schema-family
discriminator and the thing a wrong-caller-wrong-file mismatch trips on early.

### 10.4 The SSR schema (`registry_ssr`)

A new `ColumnDef` table for the spec §4.3 columns + an `SsrLocusRecord` type
and its mapping. No new codecs are needed: every SSR column lands on a codec
`block.rs` already has —

- scalars (`depth`, `n_spanning`, `start`, `end`, …) → `encode_scalar_column` /
  varint columns;
- the on-ladder histogram and ambiguous CSR (`hist_*`, `amb_*`) and the
  off-ladder block (`offl_*`) → the existing **CSR ragged-list** codecs
  (`encode_list_column_csr` / `decode_list_column_csr`);
- `chrom` and `offl_seqs` strings → the existing **bytes/dict** codecs.

So `registry_ssr` is genuinely "a table + a record mapping," as §3.2 claimed —
the byte-level work is choosing tags and wiring the mapping, not new
compression or wire code.

### 10.5 Files touched (flat layout, §3.2)

| file | change |
|---|---|
| `block.rs`, `varint.rs`, `trailer.rs` | none (already generic) |
| `index.rs` | `last_pos` = block max-`end`; interval overlap test (point = degenerate) |
| `header.rs` | add `kind`; cross-check against the `kind`-selected registry |
| `writer.rs`, `reader.rs` | parameterise on schema (registry + record encode/decode); SNP becomes one impl |
| `registry.rs` → `registry_snp.rs` | rename; `V1_0_COLUMNS` + SNP `ColumnKey` (mark schema-specific) |
| `registry_ssr.rs` (new) | SSR `ColumnDef` table + `SsrLocusRecord` + mapping |

Directory stays flat; `core/`/`snp/`/`ssr/` submodules deferred to the §3.2
crowding trigger.

### 10.6 Refactor sequence (incremental; SNP e2e tests are the gate)

Per the project's step-by-step preference, each step is a green-tests
checkpoint, and the regression gate is **the SNP caller's end-to-end tests
passing — not byte-identity** to the old format (pre-alpha, §3.2):

1. **Parameterise on schema, behaviour-preserving.** Lift `writer`/`reader` off
   the `V1_0_COLUMNS` + `PileupRecord` constants onto a schema parameter with
   SNP as the *only* implementation. The risky, no-feature-change refactor — do
   it first and alone, SNP e2e green before moving on.
2. **Add the `kind` tag** (SNP writes `kind="snp"`; reader selects registry by
   `kind`). SNP e2e green.
3. **Generalise the index overlap to intervals** (SNP = degenerate point). SNP
   e2e green.
4. **Add `registry_ssr` + `SsrLocusRecord`** (new code, zero SNP impact).
5. **Round-trip test** a synthetic `.ssr.psp` through write → columnar read →
   typed read (the first Bucket-1 hook for Stage 1).

Steps 1–3 are the "extract from the live SNP consumer" work (spec §11); 4–5 add
the second consumer. Nothing here needs Stage 1's pair-HMM or any SSR math — so
this pass can land before the Stage-0/1 design is finalised, which is why it is
roadmap item 2.

### 10.7 Open container questions

1. **Trait vs two concrete builders** (§10.2) — one `PspSchema` trait with
   generic writer/reader, or two parallel concrete builders sharing free
   functions? Lean trait, but settle against the real generics ergonomics when
   implementing (watch for `'static`-bound creep, per house style).
2. **How thin can SSR's typed-record path be?** If Stage 2 consumes only the
   columnar layer (§10.2), `SsrLocusRecord` materialisation may be needed only
   for tests and Stage-1 write — worth confirming before building a full typed
   reader for it.
3. **`ColumnKey` per schema** — confirmed-leaning: each schema keeps its own
   closed decode enum; the generic core stops at the columnar layer. Flag only
   if a shared typed-decode abstraction turns out to pay for itself (doubtful).
