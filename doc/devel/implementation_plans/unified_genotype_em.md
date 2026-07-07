# Implementation plan: shared genotype-EM core (Phase 2)

**Spec:** [../specs/unified_genotype_em.md](../specs/unified_genotype_em.md) ·
**Architecture:** [../architecture/unified_genotype_em.md](../architecture/unified_genotype_em.md) ·
**Branch:** `em-convergence-criterion`

Phase 2 extracts a shared, hook-based EM loop from
`src/var_calling/posterior_engine.rs` with the SNP caller riding it
**byte-identically**. Each step below is one `implement → code-review →
apply-fixes → commit` cycle. Every step keeps the full test suite green and the
tomato1 cohort VCF byte-identical against the pre-Phase-2 binary.

## Invariants (hold at every step)

- **SNP byte-identical.** No math change. Guardrail: 1613-test suite green +
  tomato1 cohort VCF diff = 0 vs the pre-Phase-2 binary (psp set in
  `tmp/em_bench/`), and tomato1 timing unchanged.
- **Static dispatch only.** The loop is generic/monomorphised over
  `(M: MathBackend, Model: GenotypeEmModel)`; no `dyn`, so the SNP copy is the
  current machine code.
- **Minimal diff per step.** Mechanical extraction; no opportunistic cleanup.

## Steps

### 2.1 — Architecture doc + this plan (design only) — *DONE when committed*
Pin the trait surface, module layout, SNP↔SSR mapping. No code.

### 2.2 — Introduce `GenotypeEmModel` + `SnpModel`; route M-step + driver
- Define the `GenotypeEmModel` trait (`e_step`, `m_step`, `convergence_delta`) in
  `posterior_engine.rs` (moves to its own module in 2.4).
- Add `SnpModel` (zero-sized) implementing `m_step` (calls `m_step_p_hat` +
  `m_step_f_hat_compound`) and `convergence_delta` (the current cohort/single
  branch). `e_step` on the trait for now delegates straight to `dispatch_e_step`.
- Make `run_em_loop` generic over `Model`, call `model.m_step` +
  `model.convergence_delta`; keep the E-step call as `dispatch_e_step` directly
  (routed via the model only in 2.3).
- Test: existing suite; add a unit test asserting `SnpModel` round-trips a known
  record identically to the pre-refactor path (or rely on the existing
  value tests + regression test).

### 2.3 — Move the E-step behind the trait
- `SnpModel::e_step` delegates to `dispatch_e_step` (no math change).
- `run_em_loop` and the post-loop final E-step in `run_em_columnar` call
  `model.e_step` instead of `dispatch_e_step`.

### 2.4 / 2.5 — Module relocation — **DEFERRED to Phase 3** (2026-07-07)
Originally: move the trait + `run_em_loop` into `src/var_calling/genotype_em/`.
**Folded into Phase 3** instead. Rationale: the architecture doc's real module
home is crate-level `src/genotype_em/`, reached in Phase 3 alongside the
`RecordScratch` split. Doing an intermediate move into `var_calling/genotype_em/`
now would (a) move the code twice and (b) force `pub(crate)` bumps on several
`EmContext` / `RecordScratch` fields (`n_samples`, `p_hat`, `p_hat_next`),
`EmStepPhase`, `EmDiagnostics` that the Phase-3 split reworks anyway. The
architectural value of Phase 2 — the `GenotypeEmModel` trait boundary + a generic
`run_em_loop`, SNP byte-identical — is already delivered by 2.2–2.3. Phase 3 does
the relocation once, at the crate-level destination, with the type split that
makes the visibility boundary clean. (Matches the arch doc's "avoid premature
migration" principle.)

**Phase 2 is therefore complete at 2.3.**

## Deferred to Phase 3 (not in this plan)

- Relocate the trait + generic loop to crate-level `src/genotype_em/` (was 2.4).
- Crate-level hoist + splitting `RecordScratch` into shared-core vs per-model
  state.
- The `finalize` hook (SSR posterior-homozygosity vs SNP exact-AF QUAL).
- `SsrModel` (mode-centred `G₀` seed + stutter read-model as model state) and the
  SSR benchmark-gated call change.

## Validation per step

`./scripts/dev.sh cargo fmt --check`, `cargo clippy --all-targets`,
`cargo test`, plus the tomato1 byte-identity diff (`tmp/em_bench/compare_vcf.py`)
against the pre-Phase-2 binary for the extraction steps.
