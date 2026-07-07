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

## Phase 3 — give SSR the improved prior (lightweight; SNP untouched)

Design: architecture doc §7 (lightweight, decided 2026-07-07). Goal: SSR adopts
the improved marginalize**+LOO** DM prior, seeded by its own mode-centred `G₀`,
sharing `genetics::dirichlet_multinomial_log_priors` + a small convergence-
discipline helper. **The SNP E-step is not touched** (protects the SIMD gains).
Two slim loops stay separate.

- **3.2** — Build the SSR marginalize+LOO+`F` DM prior wrapper (mode-centred `G₀`
  seed, via the shared DM primitive), as **additive** SSR code with unit tests:
  high-concentration limit ≈ plug-in HWE; hand-computed small cases; the LOO
  exclusion. SNP untouched; SSR unchanged until the 3.4 toggle.
- **3.3** — Extract a small shared convergence-discipline helper (test the driver
  `expected_counts`, not the reported `π`); the new SSR prior path uses it. SNP
  `run_em_loop` adopts it only if byte-identical, else left as-is.
- **3.4** — Wire the new prior into `run_pi_em` behind an `EmCfg` toggle
  (default = current plug-in). SSR byte-identical at default (ssr_tomato1
  unchanged); new prior opt-in.
- **3.5** — Benchmark ssr_tomato1 vs HipSTR (plug-in vs marginalized,
  `ssr_vs_hipstr_dashboard.py`). Watch specifically for the `G₀`-as-concentration
  effect (too-diffuse / too-tight prior), not just overall concordance. Flip the
  default only if it holds/improves; else keep opt-in and record the result.

Decisions folded in (2026-07-07): lightweight loop-sharing (not a generic
`run_em`); marginalize **+ LOO** for SSR; `G₀`-as-DM-concentration is an explicit
validation target in 3.5; SNP E-step untouched for perf safety.

## Validation per step

`./scripts/dev.sh cargo fmt --check`, `cargo clippy --all-targets`, `cargo test`.
Phase 2 (2.2–2.3) used the tomato1 byte-identity diff vs the pre-Phase-2 binary.
Phase 3: 3.2–3.3 additive/SNP-untouched (unit tests); 3.4 ssr_tomato1
byte-identity at default; 3.5 ssr_tomato1-vs-HipSTR concordance.
