# Code Review: ssr_call_driver

**Date:** 2026-06-24
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** the `ssr-call` Stage-2 cohort-SSR caller driver-wiring effort (milestones G→H→I→J), diff `ba96344..HEAD` on branch `ssr-cohort` — turning the placeholder TSV dump into a real VCF.
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** PR diff `ba96344..HEAD` (24 commits, +1811/−216 across 8 files; the 6 substantive files below).
- **Reviewed against:** branch `ssr-cohort`, HEAD `8e9f296`; arch [ssr_call_driver.md](../../architecture/ssr_call_driver.md), plan [ssr_call_driver.md](../../implementation_plans/ssr_call_driver.md), spec [ssr_cohort_mark2.md §4.5](../../specs/ssr_cohort_mark2.md).
- **In-scope files:**
  - [src/ssr/cohort/driver.rs](../../../../src/ssr/cohort/driver.rs) — two-pass `run()`, `build_param_set`, `genotype_locus`, `write_genotyped_chunk`, `collect_burn_in_subset`, `check_unique_sample_names`, `SsrCallError`, integration tests.
  - [src/ssr/cohort/em.rs](../../../../src/ssr/cohort/em.rs) — per-locus genotype EM + `θ_locus`/level refit (`attribute_locus`, `refit_level_multiplier`, `compute_data_ll`, `run_pi_em`, `final_calls`, `EmCfg`).
  - [src/ssr/cohort/prepass.rs](../../../../src/ssr/cohort/prepass.rs) — `G₀` decay fit, variable-loci spread accumulation.
  - [src/ssr/cohort/param_estimation.rs](../../../../src/ssr/cohort/param_estimation.rs) — `G0FitCfg`, `AlleleSpreadAccum`, `SlipProfile::add_slip`, `ParamSet`, `FixedPointAccum`.
  - [src/ssr/cohort/merge.rs](../../../../src/ssr/cohort/merge.rs) — `CohortMerger::chromosomes()` / `sample_names()` + the chromosome-table field.
  - [src/ssr/cohort/vcf_out.rs](../../../../src/ssr/cohort/vcf_out.rs) — `write_vcf_header`, `format_vcf_record`, FP control, `FILTER_DESCRIPTIONS`.
- **Deliberately out of scope:** the milestone A–F genotyping algorithm internals (rung ladder, likelihood, pair-HMM, candidate set, `em_init`, stutter kernel, simulator); the reading/merge cursor layer (`reader.rs`, `types.rs`) reviewed 2026-06-21; the SNP path. `inbreeding.rs` (`run_cohort_em`) read only to confirm the burn-in determinism chain.
- **Categories dispatched (11, in parallel):** reliability (always); errors (always); naming (always); defaults (config + fallback values); idiomatic (always); refactor_safety (always); module_structure (multi-file); unsafe_concurrency (rayon `par_iter`/`ThreadPool`, byte-identity contract); smells (always); extras (stable output + hot path + untrusted `.ssr.psp` boundary + diff-matches-intent); tooling (crate).

## 2. Verdict

**Request-changes.** One Blocker: per-sample VCF columns are emitted in present-order, not cohort-order, so any locus a subset of samples covers (the normal real-cohort case) produces a data row with the wrong column count and genotypes printed under the wrong sample names — wrong results + invalid VCF on the core output path, masked by tests that give every sample evidence at every locus. The determinism / byte-identity machinery, the decision-E hard-error path, the ploidy panic, and the emit/drop policy are otherwise correct and well-reasoned; the remaining findings cluster around silent fallbacks that contradict the no-silent-default invariant the module elsewhere upholds loudly, and untrusted-header validation.

## 3. Execution status

Run in the dev container (`./scripts/dev.sh`, Apple `container`, rustc 1.95.0), quoted verbatim:

- `cargo fmt --check` — **exit 0** (no output).
- `cargo clippy --all-targets --all-features -- -D warnings` — **exit 0**, `Finished \`dev\` profile [unoptimized + debuginfo] target(s) in 11.25s`.
- `cargo test --all-targets --all-features` — **exit 101**, the **only** failure being the pre-existing baseline bench panic: `thread 'main' panicked at benches/psp_writer_perf.rs:386:60: index out of bounds: the len is 3300000 but the index is 3300000` → `error: test failed, to rerun pass \`--bench psp_writer_perf\``.
- `cargo test --lib --all-features` — **exit 0**, `test result: ok. 1280 passed; 0 failed; 2 ignored; 0 measured; 0 filtered out` (the 142 `ssr::cohort` tests all pass).
- `cargo doc --no-deps` — **exit 0**, generated cleanly.
- `cargo audit` — **not run** (cargo-audit not installed; this diff adds no new dependencies — `Cargo.toml`/`Cargo.lock` untouched — so no new supply-chain surface).

Findings labeled "Needs verification": 0. (B1 was independently re-verified by the orchestrator against `types.rs` and `reader.rs::evidence_at`.)

## 4. Open questions and assumptions

1. **Does Stage-1 `ssr-pileup` ever emit per-sample `.ssr.psp` files that omit a catalog locus the sample didn't cover (sparse), or is every file dense over the catalog?** Affects whether B1 manifests today vs latently. The merger's `evidence_at` returns `Ok(None)` (Absent) for a sample lacking a record, the sparse-SoA `present` vector exists precisely to carry partial coverage, and a real cohort over differing BAMs/regions will have partial-coverage loci — so B1 is reachable regardless; the answer only changes how soon. **B1 is the fix either way.** (B1)
2. **Can a `CohortLocus` ever reach `genotype_locus` with `chrom_id ≥ chrom_names.len()`?** The merger hard-errors `UnknownCatalogChrom` otherwise, so the `"?"` fallback is a "cannot happen" path — but it is silent if it ever does. (M1)
3. **Is a frozen `ParamSet` guaranteed dense over `0..n_samples` for `group_of_sample`/`error_per_sample_group`/`level_per_group` at every `sample_chemistry` call?** Decision E + `build_param_set` say yes; the three `unwrap_or` fallbacks contradict that contract if it is ever broken. (M2)
4. **Are contig/sample names from `.ssr.psp` headers guaranteed free of tabs/newlines/VCF-structured chars upstream, or must `ssr-call` validate them?** The `write_vcf_header` doc says the caller validates, but only sample-name *uniqueness* is checked. (M5)

## 5. Top 3 priorities

1. **B1** — [vcf_out.rs:271](../../../../src/ssr/cohort/vcf_out.rs#L271) `format_vcf_record` emits present-order sample columns; expand to dense cohort width keyed by `locus.present`, with `./.:.:.` for absent samples. Wrong genotypes + invalid VCF on partial-coverage loci.
2. **M2 / M1** — the `sample_chemistry` triple `unwrap_or` ([em.rs:153](../../../../src/ssr/cohort/em.rs#L153)) and the `"?"` contig fallback ([driver.rs:360](../../../../src/ssr/cohort/driver.rs#L360)) silently fabricate values where decision-E / the merger contract intend a hard failure — make them total (index/panic) or loud.
3. **M3** — `period_decay`'s `FALLBACK { p: 0.5 }` ([em.rs:142](../../../../src/ssr/cohort/em.rs#L142)) is a second copy of `G0FitCfg::fallback_p`, defeating the "single fallback source of truth" the driver's own backfill doc claims; share one named const.

## 6. Findings

### Blocker

**B1: `src/ssr/cohort/vcf_out.rs:271` — Per-sample VCF columns written in present-order, not cohort-order**
- **Categories:** extras (convergent: module_structure, errors cross-category note the same "output not validated" theme)
- **Confidence:** High
- **Problem:** `format_vcf_record` builds one sample field per element of `call.calls` (`for sample in &call.calls`, [vcf_out.rs:271-295](../../../../src/ssr/cohort/vcf_out.rs#L271-L295)). `LocusCall.calls` is "Per present sample (parallel to `locus.samples`)" ([em.rs:103-104](../../../../src/ssr/cohort/em.rs#L103-L104)); its length is the *present* count, and `samples[k]`/`calls[k]` belong to cohort sample `locus.present[k]` ([types.rs:92-96](../../../../src/ssr/cohort/types.rs#L92-L96)). The cursor returns `Ok(None)` (Absent) for any sample lacking a record at a locus ([reader.rs:113-114](../../../../src/ssr/cohort/reader.rs#L113-L114)), and the merger drops only *all-absent* loci ("sparse-omit", [merge.rs:253-255](../../../../src/ssr/cohort/merge.rs#L253-L255)). So a locus covered by a *subset* of the cohort is emitted with `present` shorter than `n_samples`. The header always writes `n_samples` columns (`sample_names()`, all inputs). Nothing between `run_locus_em_with` and `format_vcf_record` re-expands present-order calls to full cohort width. Result: (a) the data row has `present_count` columns, not `n_samples` — a column-count mismatch with `#CHROM`, an invalid VCF; (b) the present samples' `GT:GQ:REPCN` land under the wrong sample names.
- **Why it matters:** Any real cohort has loci not every sample covers (differing depth/regions); each is silently malformed and mis-genotyped — wrong results on the headline code path. It is untested because every driver test (`run_emits_…`, the byte-identity and chunk tests) builds cohorts via `write_cohort`, which zips a record for *every* sample at *every* locus, so `present` always equals the full cohort and the bug is masked.
- **Suggested fix:** Expand to cohort width before formatting, placing each present call by `locus.present[k]` and `SampleCall::no_call()` (`./.:.:.`) for absent samples. Pass `present` + `n_samples` into `format_vcf_record` (it currently has neither). Add a driver test with a locus a strict subset of samples covers, asserting the row has exactly `n_samples` columns and that an absent sample's column is `./.:.:.` in the right position.
  ```rust
  // dense per-cohort fields (n_samples = sample_names.len(), threaded into the formatter)
  let mut dense: Vec<Option<&SampleCall>> = vec![None; n_samples];
  for (k, c) in call.calls.iter().enumerate() {
      dense[locus.present[k] as usize] = Some(c);
  }
  for slot in &dense {
      match slot {
          Some(c) => { /* format GT:GQ:REPCN as today */ }
          None => sample_fields.push("./.:.:.".to_string()),
      }
  }
  ```
  Note: `is_variable`, `apply_fp_control`, and `site_qual` are internally consistent on present-order and need no change — only the emitted *row* must be dense.

### Major

**M1: `src/ssr/cohort/driver.rs:360` — `chrom_names.get(...).unwrap_or("?")` silently emits a `"?"` contig instead of failing**
- **Categories:** errors, defaults (convergent; extras cross-category)
- **Confidence:** Medium (reachability requires a broken merger invariant; the corruption-if-reached is certain)
- **Problem:** `genotype_locus` resolves the contig name as `chrom_names.get(locus.locus.chrom_id as usize).map(String::as_str).unwrap_or("?")` ([driver.rs:360-363](../../../../src/ssr/cohort/driver.rs#L360-L363)). On an out-of-range `chrom_id` it writes a record with a literal `"?"` CHROM — a contig declared in no `##contig` header line, a structurally-valid-but-wrong VCF row — with no error and no log. The merger guarantees every emitted locus indexes `chrom_names` ([merge.rs:212-216](../../../../src/ssr/cohort/merge.rs#L212-L216) hard-errors `UnknownCatalogChrom`), so this is a "cannot happen" fallback; but if it does, the failure is a corrupt file, not a stop.
- **Why it matters:** Contrast the loud `UnresolvedSamples` / `DuplicateSampleName` errors a few lines away — this fallback masks a merger/id invariant break with parseable-but-wrong output.
- **Suggested fix:** Make it total. Since `genotype_locus` returns `Option<String>`, either index directly so an out-of-range id panics with context, or pre-validate once that every emitted `chrom_id < chrom_names.len()`:
  ```rust
  let chrom = chrom_names.get(locus.locus.chrom_id as usize).map(String::as_str)
      .unwrap_or_else(|| panic!("locus chrom_id {} out of range (merger invariant broken)", locus.locus.chrom_id));
  ```

**M2: `src/ssr/cohort/em.rs:153` — `sample_chemistry` silently substitutes three magic-number defaults on a frozen-param lookup decision-E says cannot miss**
- **Categories:** defaults, errors, reliability (convergent)
- **Confidence:** High (the code is as described); reachability Medium (needs a broken invariant)
- **Problem:** `sample_chemistry` resolves group/ε/level entirely through fallbacks: `group_of_sample.get(global).copied().unwrap_or(SampleGroupId(0))`, `error_per_sample_group.get(...).map(|e| e.0).unwrap_or(0.01)`, and `level_per_group.get(...).copied().unwrap_or(StutterLevel { baseline: 0.05, slope: 0.0 })` ([em.rs:153-178](../../../../src/ssr/cohort/em.rs#L153-L178)). Decision E guarantees every present sample resolves to a group (the driver hard-errors `UnresolvedSamples` in `build_param_set`, [driver.rs:118-122](../../../../src/ssr/cohort/driver.rs#L118-L122)), and the per-group vectors are sized per group — so all three branches are "cannot happen." If any invariant breaks (length disagreement, a future index bug), the EM silently genotypes the sample under fabricated chemistry (group 0 / ε=0.01 / level baseline 0.05) — wrong likelihoods → wrong calls, no error, no log. This is exactly the silent-cohort-default that `UnresolvedSamples` exists to prevent, applied one layer deeper, and the three values are un-named magic numbers in a function body.
- **Why it matters:** A broken frozen-param invariant yields confidently-wrong genotypes instead of surfacing — a swallowed-failure class bug contradicting the loud `assert_eq!(ploidy, 2, …)` style two functions up.
- **Suggested fix:** Make the lookup total — index directly so an out-of-range access panics with a clear index, matching the no-silent-default invariant:
  ```rust
  let group = params.group_of_sample[global];                // decision-E: dense; panics loudly on a broken invariant
  let eps = params.error_per_sample_group[group.0 as usize].0;
  let level = level_per_group[group.0 as usize];
  ```
  If a soft path is genuinely wanted, route it through named `const`s with a doc on when it can fire and a `tracing::warn!` naming the sample. At minimum, document why a miss is impossible.

**M3: `src/ssr/cohort/em.rs:142` — `period_decay` `FALLBACK { p: 0.5 }` is a second source of truth for the `G₀` fallback decay**
- **Categories:** defaults, smells, module_structure (convergent)
- **Confidence:** High
- **Problem:** `period_decay` hardcodes `const FALLBACK: G0PseudocountDecay = G0PseudocountDecay { p: 0.5 };` ([em.rs:142](../../../../src/ssr/cohort/em.rs#L142)), while `G0FitCfg::dev_default()` independently sets `fallback_p: 0.5` ([param_estimation.rs:98-104](../../../../src/ssr/cohort/param_estimation.rs#L98-L104)). `build_param_set` backfills every *characterized* period from `g0_cfg.fallback_p` *precisely so* "`g0_cfg.fallback_p` is the single fallback source of truth" ([driver.rs:99-102](../../../../src/ssr/cohort/driver.rs#L99-L102)). But the EM layers its own, unlinked literal underneath for any period that slips through (an unobserved period). If `fallback_p` is recalibrated and this literal is not, they diverge silently — and only on the exact periods the driver's single-source claim says it eliminated.
- **Why it matters:** Recalibration of the documented fallback fails to propagate to the EM's own fallback, producing a `G₀` prior the operator cannot trace to any config field.
- **Suggested fix:** One named source referenced by both:
  ```rust
  // param_estimation.rs
  /// Coded G₀ decay for a period with no data-driven fit (and the EM's last-resort fallback).
  pub(crate) const DEFAULT_G0_FALLBACK_P: f64 = 0.5;
  // G0FitCfg::dev_default(): fallback_p: DEFAULT_G0_FALLBACK_P,
  // em.rs period_decay():    const FALLBACK = G0PseudocountDecay { p: DEFAULT_G0_FALLBACK_P };
  ```
  Or plumb `g0_cfg.fallback_p` down into `period_decay`.

**M4: `src/ssr/cohort/vcf_out.rs:193` and `src/ssr/cohort/em.rs:460-464` — uncommented `partial_cmp(...).unwrap()` argmax over `f64` on the parallel VCF-emit path**
- **Categories:** errors, idiomatic, smells, unsafe_concurrency (convergent — flagged by four reviewers)
- **Confidence:** Medium
- **Problem:** `site_qual` picks the modal allele with `.max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())` over `call.pi` ([vcf_out.rs:188-195](../../../../src/ssr/cohort/vcf_out.rs#L188-L195)), and `final_calls` picks the MAP genotype with `.max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap()).unwrap()` over `log_joint` ([em.rs:460-464](../../../../src/ssr/cohort/em.rs#L460-L464)). `partial_cmp` returns `None` (and the unwrap panics) on a NaN. Both run inside the chunk-parallel sweep (`write_genotyped_chunk` → `genotype_locus`), so a NaN panics a rayon worker and aborts the run with no context. The surrounding code does establish the no-NaN invariants (normalized π from `run_pi_em`; `lambda`-floored finite read likelihoods; `k ≥ 1` for a PASS locus), but neither call carries the `// PANIC-FREE:` marker the crate convention requires, and `total_cmp` would remove the unwrap entirely. The concurrency review confirmed this does not affect determinism (deterministic given finite inputs).
- **Why it matters:** A hidden-invariant break becomes an abrupt worker-thread panic on the production emit path instead of a typed error; the convention exists so the next editor cannot silently invalidate the precondition.
- **Suggested fix:** Use the total comparator and document the invariant:
  ```rust
  // PANIC-FREE: pi values are normalized finite frequencies from run_pi_em (Σ=1, no NaN).
  .max_by(|(_, a), (_, b)| a.total_cmp(b))
  // em.rs final_calls: log_joint is non-empty (a PASS locus has ≥1 genotype) and finite-or-−inf.
  .max_by(|(_, a), (_, b)| a.total_cmp(b)).expect("a PASS locus enumerates ≥1 genotype")
  ```

**M5: `src/ssr/cohort/vcf_out.rs:88` / `merge.rs:280` — no guard against tabs / control / VCF-reserved chars in contig and sample names from untrusted `.ssr.psp` headers**
- **Categories:** extras
- **Confidence:** High
- **Problem:** Contig names come from the first input's `.ssr.psp` chromosome table ([merge.rs:136](../../../../src/ssr/cohort/merge.rs#L136), via `chromosomes()`); sample names from path basenames ([merge.rs:287-292](../../../../src/ssr/cohort/merge.rs#L287-L292)). `write_vcf_header` writes both verbatim into `##contig=<ID={name},…>` and the tab-separated `#CHROM` columns ([vcf_out.rs:98-128](../../../../src/ssr/cohort/vcf_out.rs#L98-L128)). The function doc claims the *driver* validates "VCF-clean sample/contig names (no tabs)" ([vcf_out.rs:83-87](../../../../src/ssr/cohort/vcf_out.rs#L83-L87)), but the driver checks only sample-name *uniqueness* (`check_unique_sample_names`) — nothing checks for `\t`/`\n`/`\r`/`<`/`>`/`,`, and contig names are never validated at all. A `.ssr.psp` is untrusted input at the merge boundary: a tab in a sample basename shifts every later column; a newline in a contig name splits the meta line; a `,`/`>` breaks the structured `##contig`.
- **Why it matters:** Untrusted header content can silently produce a malformed/mis-parsed VCF or smuggle extra columns — a validator-boundary gap the doc asserts is handled but no caller fulfills.
- **Suggested fix:** Add a validation pass in `run` before `write_vcf_header` rejecting any contig/sample name containing `\t`/`\n`/`\r` (and VCF-structured chars for contigs) with a typed error, or escape per VCF 4.4 §1.4.7. Add malformed-input tests asserting the error variant for a tab-bearing contig and a tab-bearing sample basename.

**M6: `src/ssr/cohort/driver.rs:318` — `genotype_locus`'s filtered-locus *emit* branch has no test**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** `genotype_locus` implements the whole emit/drop/FILTER policy: PASS+monomorphic ⇒ `None` (drop), PASS+variable ⇒ `Some`, and any non-PASS admission ⇒ `Some` (kept with reason). The drop and PASS+variable branches are covered e2e by `run_emits_a_vcf_with_a_pass_variant_and_drops_a_monomorphic_locus`, but the **filtered-locus-emitted** branch (`candidates.admit != Admission::Pass` must still produce a record line) is never exercised through `run`/`genotype_locus`. `vcf_out.rs` tests `format_vcf_record` on a hand-built `LowDepth` call, but nothing pins that `genotype_locus` *emits* (rather than drops) a filtered locus. If the `Admission::Pass &&` guard ([driver.rs:356](../../../../src/ssr/cohort/driver.rs#L356)) were ever widened/inverted, filtered loci would be silently dropped and no test would fail — a silently-incomplete VCF (wrong output without a panic).
- **Why it matters:** "Filtered loci kept with reason, monomorphic dropped" is the driver's headline contract; the kept-with-reason half is untested end-to-end.
- **Suggested fix:** Add `run_emits_a_filtered_locus_with_its_reason` — a cohort with one resolvable locus plus one locus below the depth floor (`Admission::LowDepth`); assert the VCF contains a data record for the thin locus with `FILTER == "lowDepth"` and `./.:.:.` columns, and that the PASS variant still emits. (See §8.)

### Minor

**Mi1: `src/ssr/cohort/em.rs:507` / `prepass.rs:159` / `vcf_out.rs:134` — read-to-nearest-allele attribution duplicated across 3 in-scope sites (4 with the out-of-scope `inbreeding::reduce_level`) with no shared primitive.**
Categories: module_structure, smells (convergent). The "find the nearest parent allele by `|read_units − allele_units|`, bin faithful/slipped" kernel is re-derived in `attribute_locus`, `accumulate_locus`, and `allele_balance`. Confidence High that it is duplicated; the smells reviewer's "inconsistent tie-break" claim was checked and does **not** hold today — `min_by_key` (first-min) and `allele_balance`'s `<=` both resolve a tie toward the first/smaller allele, so they agree. The finding is maintainability/drift, not a current behavioral divergence: the modules all document a planned soft-split refinement that would have to touch every copy. Extract one helper near the `SlipProfile`/`add_bin` sinks (a small `attribution` peer, since callers span three files) returning `(delta, parent, count)` per read, with the tie-break explicit and shared.

**Mi2: `src/ssr/cohort/em.rs:341` — `compute_data_ll`'s 11-positional-argument signature, repeated verbatim at two call sites, warrants a `LocusModel` bundle.**
Categories: smells, tooling (cross-category). `compute_data_ll(locus, candidates, params, level_per_group, theta, level_mult, cand_units, genotypes, lambda, distinct, scratch)` is called at [em.rs:282](../../../../src/ssr/cohort/em.rs#L282) and [em.rs:310](../../../../src/ssr/cohort/em.rs#L310) with identical 11-arg lists. Several args (`cand_units`, `genotypes`, `distinct`) are once-computed locus shaping re-passed positionally; two adjacent `f64`s (`level_mult`, `lambda`) and two adjacent slices invite a silent transposition the compiler won't catch. Bundle the iteration-invariant locus shaping into a `LocusModel { candidates, cand_units, genotypes, distinct }` built once in `run_locus_em_with`, dropping `compute_data_ll` to ~6 args and removing the duplicated call — the project's "split data-shaping from math" preference.

**Mi3: `src/ssr/cohort/em.rs:355` — `compute_data_ll` reallocates the full per-sample likelihood structure every refit round.**
Category: idiomatic. Each call builds a fresh `Vec<Vec<f64>>` (`data_ll`), and per sample a fresh `obs_qr: Vec<(u32, Vec<f64>)>`, a `row: Vec<f64>` of length `k` per distinct obs, and a `sample_ll: Vec<f64>`. Only `theta`/`level_mult` change between rounds; the iteration *shape* is fixed per locus. On the hot per-locus path in the chunk-parallel sweep this is `O(rounds · samples · obs · k)` short-lived `Vec`s — the project's "prefer scratch over per-iteration fresh allocations" lever. Hoist `data_ll`/`obs_qr` into `run_locus_em_with` (reuse via `clear()`), or extend the already-threaded `LikelihoodScratch` to own them. Confirm with a real-cohort bench before committing.

**Mi4: `src/ssr/cohort/driver.rs:337` — `f_present` allocates a fresh `Vec<f64>` per locus inside the parallel sweep.**
Category: idiomatic. `genotype_locus` builds `f_present` per `par_iter` item. A `par_iter().map_init(Vec::new, |buf, locus| { buf.clear(); buf.extend(...); … })` lets each rayon worker reuse one buffer. Small (× locus count); verify with a bench before committing — at small `present` lengths it may be noise.

**Mi5: `src/ssr/cohort/prepass.rs:110,251` — `(u16, u64, u64)` / `(u16, u64)` tuple bins with linear-scan `find` accumulation (primitive obsession).**
Category: smells. `add_bin` reads `bin.1 += faithful; bin.2 += slipped;` — opaque positional fields whose meaning lives only in a doc comment, the same tuple shape repeated in `merge_sample_stats` and `fit_level`. The O(n²) scan is fine (distinct-length cardinality is tiny). Replace the tuple with `struct LengthBin { length: u16, faithful: u64, slipped: u64 }` so `bin.faithful += faithful` reads itself; keep the `Vec` + a one-line note that small cardinality makes the scan intentional.

**Mi6: `src/ssr/cohort/driver.rs:372`, `em.rs:234`, `em.rs:340` — three new bare `#[allow(clippy::too_many_arguments)]` carry no justification comment.**
Category: tooling. The project discipline (cf. the rationale on `result_large_err = "allow"` in `Cargo.toml`) is that per-site `#[allow]` records *why*. Add a one-liner at each (e.g. "the frozen params + dev-config structs are threaded explicitly to keep the per-locus call a pure function of its inputs — the byte-identity contract; arg count intentional"), or bundle the configs (see Mi2).

**Mi7: `src/ssr/cohort/em.rs:864-878` — `LocusSlipFit` test helpers use `..Default::default()`, silently absorbing the `profile` field.**
Category: refactor_safety. `let matched = LocusSlipFit { slipped: 10, expected_slipped: 10.0, ..Default::default() };` (and `stuttery`). If `LocusSlipFit` later gains a field feeding `refit_level_multiplier`, these literals keep compiling and test the new field at its `Default` — masking the change. Spell `profile: SlipProfile::default()` explicitly (zero behavior change, makes a future field a compile error).

**Mi8: `src/ssr/cohort/em.rs:281` (and :148/:303/:309/:365) — `level_mult` uses a non-standard abbreviation for "multiplier".**
Category: naming. A long-lived local in `run_locus_em_with`, a `compute_data_ll` parameter, and `new_level_mult` — `mult` is not a Rust-standard shorthand, and the spelled form ("per-locus rate multiplier", `refit_level_multiplier`) is used everywhere else for the same value. Rename `level_mult` → `level_multiplier` throughout em.rs (mechanical).

**Mi9: `src/ssr/cohort/driver.rs:175` — `FrozenParams.params` is a bare generic noun.**
Category: naming. Beside the well-named `f_per_sample`/`level_per_group`, `params: ParamSet` reads as "the params of the frozen-params." Its doc already calls it "Frozen chemistry"; rename the field to `chemistry` and update the three reads in `genotype_locus` + the constructor in `run`.

**Mi10: `src/ssr/cohort/driver.rs:545` — stale "the sweep is serial" comment in `run_is_byte_identical_across_thread_counts`.**
Categories: reliability, refactor_safety (cross-category). Post-Milestone-J the sweep is chunk-parallel (`write_genotyped_chunk` uses `par_iter`), as the module header and `chunk_parallel_sweep_orders_records_and_is_deterministic` state. The comment understates what the test now guards (it runs `threads = 4`). Reword to "the burn-in reduces and the chunk-parallel sweep's order-preserving `par_iter` collect are both byte-identical across threads."

**Mi11: `src/ssr/cohort/em.rs:44,218` — `EmCfg::inbreeding_f` is effectively test-only.**
Categories: defaults, smells (convergent). The production path (`genotype_locus` → `run_locus_em_with`) builds `f_present` from `frozen.f_per_sample` and never reads `cfg.inbreeding_f`; only the `run_locus_em` convenience wrapper (tests) uses it. A caller setting `EmCfg { inbreeding_f: 0.1, .. }` and calling `run_locus_em_with` directly gets it silently ignored — a latent misuse trap. Document on the field that it seeds only the `run_locus_em` wrapper, or drop it and pass `f: f64` to that wrapper.

**Mi12: `src/ssr/cohort/driver.rs:207,246` — `config.threads.max(1)` and `queue_depth == 0 ⇒ DEFAULT_SWEEP_CHUNK` coercions are unsurfaced.**
Category: defaults. `threads == 0` quietly becomes single-threaded with no field doc and no log; the `queue_depth == 0` sentinel is documented on the const but not on the `SsrCallConfig::queue_depth` field callers read, and neither logs the resolved value. Low impact (output is chunk-size-invariant, tested) — document the `0 ⇒ default` sentinels on the fields and emit a `debug`/`eprintln` recording the resolved thread count + chunk size, or reject `0` threads.

**Mi13: `src/ssr/cohort/vcf_out.rs:257-261` — `String::from_utf8_lossy` silently substitutes U+FFFD for a non-UTF8 REF/ALT tract.**
Category: extras. Allele bytes come from untrusted catalog/psp evidence; a non-ACGTN/invalid-UTF8 byte becomes a replacement char in the VCF rather than a loud failure. Unlikely (alleles should be ACGTN) but it is the untrusted-input case. Validate alleles are ASCII/ACGTN at the merge boundary (typed error), or use `std::str::from_utf8(...)?`; reserve `from_utf8_lossy` only with a documented upstream guarantee.

**Mi14: `src/ssr/cohort/vcf_out.rs:300-305` — a variable PASS locus whose `site_qual` clamps to `0.0` emits `QUAL=.` rather than a number.**
Category: extras. The emit/drop decision is correctly on `is_variable`, not QUAL (good). But `qual_field` prints `.` for `qual <= 0.0`. A locus emitted *because* variable but with proxy QUAL clamped to exactly 0 shows `QUAL=.` — valid VCF (spec §4.5 allows `.`), but downstream QUAL filters see "missing" instead of "lowest." Either always print numeric QUAL for emitted records, or document the `.`-means-zero convention and pin it with a test.

**Mi15: `src/ssr/cohort/em.rs:299` — the per-locus refit loop's max-rounds-exhausted (non-convergence) path is untested.**
Category: reliability. Tests cover `refit_max_rounds: 0` and converging cases; none pins that on non-convergence the returned `calls`/`posterior_hom` are the last-recomputed values (they are — each iteration recomputes at the end), not stale pre-loop values. A refactor moving the recompute would go uncaught. Add `run_locus_em_with_returns_calls_consistent_with_final_theta_on_non_convergence` (`refit_max_rounds: 1` on a stuttery locus; assert the call equals the `refit_max_rounds: 3` genotype at high depth).

**Mi16: `src/ssr/cohort/driver.rs:217` — `grouped.level_per_group.clone()` clones a `Vec` dead after the call.**
Category: idiomatic. Inside the `pool.install` closure, `build_param_set(&est, &grouped, …)` borrows `grouped`, then `run_cohort_em(…, grouped.level_per_group.clone(), …)` clones; `grouped` is unused afterward. Once-per-run, so low impact — restructure to move `level_per_group` after the `build_param_set` borrow ends, or accept with a `// once-per-run clone` note.

### Nits

- `src/ssr/cohort/driver.rs:64` — `SsrCallError::Write(#[from] std::io::Error)` collapses every write/create/flush site into one variant with no path/operation context; mirror `SsrMergeError::Open { path, source }` if richer I/O context is wanted (internal error, low priority).
- `src/ssr/cohort/em.rs:522-526` — `.expect("a non-empty call has ≥1 allele")` is genuinely unreachable (the `is_empty` `continue` two lines up) but lacks the `// PANIC-FREE:`/`// UNREACHABLE:` marker the crate convention uses.
- `src/ssr/cohort/param_estimation.rs:147` / `em.rs:562` — `SlipProfile::add_slip` (the `delta == 0` / `|Δ| > MAX_SLIP` boundaries) and `phred_gq` (the `1e-10` floor + 99 cap) have no direct unit test; both are only exercised transitively (see §8).
- `src/ssr/cohort/param_estimation.rs:228` — `FixedPointAccum::add` bounds finiteness with a `debug_assert!` only; the magnitude contract (`< i128::MAX / 2^40`) is documented but a release-mode overflow saturates silently. (Mostly outside the changed-line set; noted.)
- `src/ssr/cohort/vcf_out.rs:298` — the `format!` template mixes positional (`filter_text`, `sample_fields`) and named (`qual_field`) captures; all-named would be less order-coupled.
- `src/ssr/cohort/driver.rs:153` — `const PLOIDY: u8 = 2` + `assert_eq!(ploidy, 2, …)` thread a raw `u8` through every signature; accepted primitive obsession while v1 is diploid-only, the obvious first casualty when polyploid lands.
- `src/ssr/cohort/driver.rs:235` — `eprintln!` for the cohort warning rather than `tracing`/`log`; consistent with the crate-wide convention (no `tracing` dep), flagged only for completeness.

## 7. Out of scope observations

- `src/ssr/cohort/driver.rs:159-171` (`BURN_IN_MAX_LOCI` positional first-`cap` subset) — the decision-E resolution check runs on this subset, so positional selection could in principle reject a sample whose confident loci all fall past the cap. Acknowledged in-comment as a calibration follow-up (reading Q-R5); raised by reliability/concurrency/errors as a cross-category note. Not re-litigated per the review scope.
- `src/ssr/cohort/inbreeding.rs` (`run_cohort_em`, `reduce_f`, `reduce_level`) — read only to confirm the burn-in float reduces route through `FixedPointAccum`/integer `add_bin` (they do; byte-identity sound). Out of the changed-line set; no findings.
- The `dev_default()` constant *values* (`BURN_IN_MAX_LOCI`, `DEFAULT_SWEEP_CHUNK`, `theta_*`, `level_*`, `max_iters`, `LEVEL_MULT_MAX`, FP-control thresholds) are uncalibrated by design — the separate downstream calibration effort, not re-filed here.

## 8. Missing tests to add now

Grouped by function under test; `function_returns_expected_on_condition` naming.

**`format_vcf_record` / `run` (B1):**
- `run_emits_dense_sample_columns_for_a_partial_coverage_locus` — a cohort where a strict subset of samples covers a locus; assert the data row has exactly `n_samples` tab-separated sample columns and an absent sample's column is `./.:.:.` in its cohort position. Catches B1 (present-order columns).

**`genotype_locus` / `run` (M6, emit policy):**
- `run_emits_a_filtered_locus_with_its_reason` — one resolvable locus + one `LowDepth` locus; assert a data line with `FILTER == "lowDepth"`, `./.:.:.` columns, and that the PASS variant still emits. Catches a widened/inverted drop guard.
- `run_writes_header_only_for_an_all_monomorphic_cohort` — every locus hom-ref for every sample; assert a `#CHROM` line and zero non-`#` data lines. Catches chunk-flush / header-record interleaving with zero emitted records.
- `run_on_an_empty_psp_file_list_errors` — `psp_files: vec![]`; assert `Err(SsrCallError::Merge(SsrMergeError::NoInputs))` reached through `run`.

**`run_locus_em_with` (Mi15, refit termination + determinism):**
- `run_locus_em_with_returns_calls_consistent_with_final_theta_on_non_convergence` — `refit_max_rounds: 1` on a stuttery locus; assert the capped call equals the `refit_max_rounds: 3` genotype at high depth.
- `sweep_is_byte_identical_across_threads_when_loci_refit` — drive a shape-mismatched high-stutter cohort (one that actually triggers ≥1 refit round, exercising the `expected_slipped` float sum) through `driver::run` at `threads = 1` vs `4`; assert byte-identical VCFs. The current byte-identity tests use clean loci that may not exercise the refit float path.

**`SlipProfile::add_slip` (Nit):**
- `add_slip_drops_magnitudes_beyond_the_cap_and_ignores_zero` — `delta = 0`, `±(MAX_SLIP + 1)`; assert all bins stay zero, then `add_slip(MAX_SLIP, 3)` lands in `up[MAX_SLIP - 1]`. Catches an off-by-one in `(1..=MAX_SLIP)` / `magnitude - 1`.

**`phred_gq` (Nit):**
- `phred_gq_caps_at_99_and_floors_at_zero` — `assert_eq!(phred_gq(1.0), 99); assert_eq!(phred_gq(0.0), 0);` plus a mid value. Catches removal of the `1e-10` floor (→ `-inf`) or the 99 cap.

**`sample_chemistry` (M2, pin current behavior until made total):**
- `sample_chemistry_falls_back_to_coded_defaults_on_a_missing_group` — a deliberately short `group_of_sample`; assert the documented `(0.01, baseline 0.05)` fallback, so any change (e.g. to a panic per M2) is a conscious one.

## 9. What's good

- **The cross-thread byte-identity argument is genuinely sound and was traced end-to-end.** The prepass `fold/reduce` folds only integer accumulators ([prepass.rs:68-83](../../../../src/ssr/cohort/prepass.rs#L68-L83)); the one per-locus `f64` reduce (`attribute_locus::expected_slipped`) is order-fixed because it is consumed entirely within the single serial task that produces it ([em.rs:502-506](../../../../src/ssr/cohort/em.rs#L502-L506) documents exactly this); the sweep's `par_iter().collect::<Vec<Option<String>>>()` is an indexed collect over a slice, order-preserving regardless of thread count.
- **The `PrepassStats::merge` doc honestly distinguishes the determinism layer from the struct form** — "the struct itself is not a canonical form; determinism tests must compare `EstimatedParams`, not raw stats" ([prepass.rs:56-67](../../../../src/ssr/cohort/prepass.rs#L56-L67)) — and the tests obey it.
- **`FixedPointAccum` debug-asserts finiteness on `add`** ([param_estimation.rs:228-234](../../../../src/ssr/cohort/param_estimation.rs#L228-L234)) so a NaN/overflow cannot silently defeat the associativity guarantee, with bit-identical order/grouping-independence unit tests.
- **`filter_text` + `FILTER_DESCRIPTIONS` match `Admission` with no `_` arm** ([vcf_out.rs:46-69](../../../../src/ssr/cohort/vcf_out.rs#L46-L69)), so the header declarations and per-record FILTER values are one compiler-checked source of truth and a new verdict is caught at both.
- **The decision-E hard error names the offending sample indices and explains the degenerate-input causes** ([driver.rs:65-77](../../../../src/ssr/cohort/driver.rs#L65-L77)) — the no-silent-default invariant made loud and legible (the very contract M2's deeper `sample_chemistry` fallbacks should match).

## 10. Commands to re-verify

Re-run in the dev container (`./scripts/dev.sh bash -c '<cmd>'`):
- `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --lib --all-features` (1280 pass; the only `--all-targets` failure is the pre-existing `psp_writer_perf` bench panic)
- `cargo doc --no-deps`

New invocations the review introduces (after B1/M6 fixes land): the §8 tests, in particular `run_emits_dense_sample_columns_for_a_partial_coverage_locus` and `run_emits_a_filtered_locus_with_its_reason`.

### Author response convention

Address each finding by id (B1, M1…M6, Mi1…Mi16) with one of: `fixed in <commit>` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer the §4 open questions first — especially Q1, which determines how urgently B1 manifests (the fix is required either way).
