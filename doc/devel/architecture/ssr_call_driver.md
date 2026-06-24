# SSR Stage 2 — wiring the genotyper into the `ssr-call` driver (architecture sketch)

**Status:** settled, 2026-06-23, branch `ssr-cohort`. Decisions A–E resolved (§11) after
discussion; ready to implement. The genotyping algorithm (milestones A–F) and the
reading/merge spine are both implemented, reviewed, and tested; they already compose
end-to-end in one test
([`inbreeding.rs::full_pipeline_calls_and_emits_a_variant_vcf_line`](../../../src/ssr/cohort/inbreeding.rs)).
What is missing is the driver that runs that composition over a real cohort and writes a
**VCF** instead of the current Phase-1 TSV dump. **Headline (decision A):** *not* the
test's in-memory materialization — a **two-pass streaming** design (a bounded
pre-pass/burn-in that freezes the cross-locus-pooled parameters, then a single streaming
genotyping sweep of independent per-locus EMs), so memory stays bounded at cohort scale.

It is the fourth Stage-2 sketch, the *how-it-all-runs* companion to the three that
came before:

- [ssr_call_reading.md](ssr_call_reading.md) — the reader/merge spine (the
  `CohortMerger` this driver consumes), incl. §8's two-pass (re-read) decision.
- [ssr_call_parameters.md](ssr_call_parameters.md) — the pre-pass (`run_prepass_stats`
  → `estimate` → `group_samples`) that produces the chemistry/priors.
- [ssr_call_genotyping.md](ssr_call_genotyping.md) — candidate assembly, the
  likelihood, the per-locus EM, the `F` loop, and the VCF semantics (§6).

Where this doc and the spec [ssr_cohort_mark2.md](../specs/ssr_cohort_mark2.md)
disagree on *intent*, the spec wins; on *code layout*, this doc wins.

---

## 1. What this phase produces

Today `driver::run` opens the cohort, streams `CohortLocus`es, and writes a
catalog-ordered **TSV dump** (a placeholder, `--threads`/`--queue-depth` reserved).
After this task `ssr-call` produces a **valid VCF**: a header (contigs, FORMAT, INFO,
FILTER, the `#CHROM … <samples>` line, the apparent-`F_IS` warning) plus one data line
per emitted locus, genotypes correct at depth, PASS and filtered records present,
monomorphic loci dropped — *or* a hard error if the pre-pass cannot resolve confident
genotypes for every sample (§7, decision E).

Everything the body needs already exists as functions; the driver is **orchestration,
not new algorithms**. The one genuinely new piece of logic is `build_param_set`
(§3) — currently inlined in the test.

---

## 2. The composition

The *algorithm* is the one the working test
(`inbreeding.rs::full_pipeline_calls_and_emits_a_variant_vcf_line`) already composes
end-to-end — the same pieces. The driver runs them in the **two-pass streaming** shape
of decision A (§4), **not** in the test's in-memory `Vec<CohortLocus>` (materialization
is rejected — §4):

```
 PASS 1 — pre-pass / burn-in   (bounded: a stratified SUBSET of loci held in RAM, §4)
   run_prepass_stats ─► estimate ─► group_samples ─► fit G₀ (§9) ─► estimate F
        │
        └─ build_param_set ─► the FROZEN ParamSet                          ← NEW (§3)
              (ε, θ_period parent, per-group level line, G₀ decay, F_i — all frozen)

 PASS 2 — genotyping   (single streaming sweep over ALL loci; ≈ N×lockstep blocks resident)
   for each merged CohortLocus, on the frozen params:
     build_rungs ─► assemble_candidates ─► run_locus_em (π + genotypes;
        Step 2: + LOCAL θ_locus / rate refit, shrunk toward the frozen priors)
     ─► apply_fp_control
     ─► is_variable?  ── no & PASS ─► DROP (monomorphic)
                      └─ yes / filtered ─► site_qual ─► format_vcf_record ─► write
```

The driver wraps this in: open → pre-pass + freeze → write VCF header → stream-genotype-emit
→ flush. The header needs the cohort sample names and the contig table (§5); the body
recomputes the per-locus candidate assembly per streamed locus (cheap, pure — §6). Note
there is **no `run_cohort_em` whole-cohort outer loop** in the genotyping path: freezing
the pooled params (§4) collapses it to one streaming sweep of independent per-locus EMs.

---

## 3. `build_param_set` — the pre-pass → EM interface

The test hand-assembles a `ParamSet` from the pre-pass outputs. That assembly is the
real, missing public function:

```rust
fn build_param_set(
    est: &EstimatedParams,      // eps, shape_by_period, level_by_sample
    grouped: &GroupedParams,    // group_of_sample, eps/level/shape per group, n_groups
    n_samples: usize,
) -> ParamSet
```

Field-by-field it is mechanical, with two judgement calls:

| `ParamSet` field | source |
|---|---|
| `error_per_sample_group` | `grouped.eps_per_group` → `PerBaseError` |
| `stutter_shape_parent` | `est.shape_by_period` |
| `stutter_shape_by_cell` | `grouped.shape_by_group_period` |
| `level_seed` | `grouped.level_per_group` |
| `group_of_sample` | dense `Vec<SampleGroupId>` over `0..n_samples`; **hard error** if any sample is absent (no confident genotype) |
| `pseudocount_decay_per_loci_group` | `est.g0_by_period` — **now fit in the pre-pass** (decision B, settled 2026-06-23, §9) |
| `f0_seed` | `0.0` (the `F` loop estimates it) |

**Two things to write down (no silent assumptions):**

1. **`group_of_sample` — a sample with no confident genotype is a hard error
   (decision E, settled 2026-06-23).** `grouped.group_of_sample` is a
   `HashMap<u32, SampleGroupId>` keyed only on samples that produced confident
   genotypes in the pre-pass. The `ParamSet` wants a dense `Vec` of length
   `n_samples`. A sample absent from that map never resolved a confident genotype in
   the whole cohort — that is **weird** (a blank/degenerate sample, a mis-supplied
   file, a catalog/chemistry mismatch), so `build_param_set` **errors hard** naming the
   offending sample(s) rather than silently giving it cohort-default chemistry. If
   real-data runs show this is common, we investigate the cause and the right
   remedy then; until then, loud failure over a silent default. This makes
   `build_param_set` fallible (`Result<ParamSet, _>`; a new `SsrCallError` variant
   naming the unresolved samples), and it **subsumes the m2(a) cohort-wide fallback** —
   see §7.

2. **`G₀` decay (decision B — fit it).** `build_param_set` reads
   `est.g0_by_period` straight into `pseudocount_decay_per_loci_group`; the pre-pass now
   fits the per-period decay `p` from the cohort's confident-genotype allele spread (the
   coded `{ period 2 → p = 0.5 }` default disappears from the call path). The pre-pass
   always emits a `p` for every period present — a real fit where the data support it, a
   coded fallback for thin periods — so the driver synthesizes nothing. See §9.

---

## 4. Decision A — memory model (settled 2026-06-23): two-pass streaming, freeze the pooled params

**Materialization is rejected.** Holding every `CohortLocus` in RAM (the test's
in-memory `Vec`) breaks the project's cohort-scaling thesis at genome scale (~1M loci).
The settled design is **two streaming passes** — a bounded pre-pass/burn-in, then a
single streaming genotyping sweep — made coherent by one principle:

> **Freeze the *pooled* (cross-locus) parameters after the pre-pass; refine the *local*
> (per-locus) ones during genotyping.**

**Why genotyping looked multi-pass.** Two quantities are estimated by pooling across
loci, and the current `run_cohort_em` re-estimates them in an outer loop (`reduce_f`,
`reduce_level`) that re-reads every locus's reads each round (≈8 rounds): the
per-individual inbreeding `F`, and the per-group stutter **level**. The level enters the
*read likelihood* (recomputed inside `run_locus_em_with`), so re-fitting it needs the
raw reads again every round — the thing that would otherwise force materialization or
~8× re-decompression.

**Why it collapses to one pass.** `F` enters only the genotype *prior*, never the reads.
The stutter level/shape *parents* are chemistry properties, estimated far more precisely
by pooling than HipSTR's per-locus fits (research §6f). So we **freeze the pooled
quantities after the pre-pass** and let each locus adapt **locally**:

| quantity | scope | settled |
|---|---|---|
| `F` (per individual) | pooled | **frozen** after the pre-pass — estimated in the burn-in, then a fixed prior |
| stutter level line (per group) | pooled | **frozen** — the chemistry/length anchor |
| `θ_period` shape parent (per period) | pooled | **frozen** |
| `ε` (per group) | pooled | **frozen** (already) |
| `G₀` decay `p` (per period) | pooled | **frozen** — fit in the pre-pass (§9) |
| `θ_locus` shape (this locus) | **local** | **refined per-locus**, shrunk toward `θ_period` (Step 2) |
| this locus's stutter **rate** | **local** | **refined per-locus**, shrunk toward the group line (Step 2) |

Per-locus refinement uses only the locus's own reads (already in hand), so it adds **no**
cross-locus coupling: genotyping stays a single streaming pass, each locus a pure
function of (its reads + the frozen priors) — **no outer loop, no likelihood cache, no
re-read** — and byte-identical across threads by construction (no cross-locus reduce
whose order could matter).

**The burn-in (lives in the pre-pass).** The spec's burn-in does the iterating the outer
loop used to: it iterates the global estimates — `ε`, `θ_period`, the per-group level
line, `G₀`, and `F_i` — over a **representative stratified subset** of loci held in RAM
(bounded by the *subset*, not the cohort), settles, and **freezes** them. This is the
analog of HipSTR's "several rounds before genotyping," except we pool across a subset of
loci rather than within one locus. The `reduce_f` / `reduce_level` logic moves here;
`run_cohort_em`'s whole-cohort outer loop is **retired** for the genotyping path.

**Shrinkage is the anti-oscillation knob (the per-locus stutter refit).** A high-depth
locus's own reads dominate the shrinkage → `θ_locus` / rate move to the locus's MLE
(HipSTR-like adaptation; removes the frozen-param bias). A thin / ambiguous locus cannot
overcome the prior → it stays at the frozen group value (no oscillation). No hard depth
threshold — the shrinkage strength does it smoothly (`shape_shrink_strength`, reused from
the per-`(group, period)` fit). This is HipSTR's per-locus stutter **plus** a
cohort-pooled shrinkage anchor HipSTR lacks — strictly better at low coverage (research
§6f, line 782).

**Sequencing — two milestones, both streaming.**
- **Step 1 — ships the streaming VCF.** Freeze *all* params from the pre-pass (`ε`,
  `θ_period`, group level, `G₀`, `F`), **drop the outer loop**, stream loci →
  `run_locus_em` on frozen params per locus → emit. This matches the per-locus EM's
  current capability (it already runs on fixed params — `em.rs` freezes `ε`/`θ`/level
  today), so it is a *small* change, and it produces a real VCF for calibration. It
  carries a bounded frozen-param bias, flagged honestly.
- **Step 2 — the per-locus adaptation.** Build the deferred per-locus slip attribution →
  the `θ_locus` shape M-step **+** the shrunk per-locus rate refit (`em.rs:12` calls the
  `θ_locus` M-step "deferred until D wires the slip accumulators"). Local, so it drops
  into the same single-pass structure and removes the Step-1 bias where the data support
  it.

**Reading-layer fit.** The genotyping sweep is the consumer of
[ssr_call_reading.md §5](ssr_call_reading.md)'s producer → bounded-queue → worker →
writer topology; the pre-pass is a separate, earlier consumer of the same merger
(reading §8: *re-read, don't cache* — moot under materialize, live here). The pre-pass
reads a **bounded subset**; genotyping streams all loci at ≈ N × lockstep blocks
resident (`--block-window-bp` + queue depth are the RSS knobs).

> **J realization (settled 2026-06-24): chunk-parallel `par_iter`, not the channel
> pipeline.** Milestone J achieves the topology's guarantees — parallel, bounded,
> ordered, byte-identical — with far less machinery: the sweep accumulates a bounded
> chunk of loci, genotypes it on the `--threads` pool with an **order-preserving
> `par_iter`** (so `lines[i]` stays aligned to `chunk[i]`), and writes the chunk in
> catalog order; `config.queue_depth` is the chunk size. Because the chunk is index-
> ordered there is **no `seq`-reorder** to do (the reorder-by-seq writer is only needed
> when workers finish out of order, which they can't here). Trade-off vs the full
> producer/worker/writer channel pipeline: chunking reads a chunk *then* processes it, so
> it does not overlap the merger read with genotyping — a small serialization, dwarfed by
> the per-locus EM at any reasonable chunk size. The fully-overlapping channel pipeline
> is a **measure-first** follow-up, not built.

---

## 5. The VCF header + writer (decision C)

`format_vcf_record` ([vcf_out.rs](../../../src/ssr/cohort/vcf_out.rs)) already emits a
**plain-text** data line (`CHROM…FORMAT <samples>`, `GT:GQ:REPCN`, `PERIOD` INFO,
PASS/`notPeriodic`/`tooManyAlleles`/`lowDepth` FILTER). What is missing is the
**header** and the line that names the samples.

**Decision C — recommend a small dedicated SSR text-header writer in `vcf_out.rs`,
not the SNP `src/vcf/writer.rs`.** Rationale:

- `format_vcf_record` already produces text, not `noodles_vcf` records. The SNP writer
  (`build_vcf_header` → `noodles_vcf::Header`) declares **SNP** INFO/FORMAT
  (`GP`, `EMNoConv`, …) and serializes via noodles. Bolting the SSR text-record path
  onto a noodles header mixes two serialization styles for one file. A handful of
  `writeln!(out, "##…")` lines matching the existing text records is simpler and keeps
  the whole SSR VCF in one vocabulary.
- The header is small and fully SSR-specific: `##fileformat=VCFv4.4`; one `##contig`
  per chromosome; `##INFO=<ID=PERIOD,…>`; `##FORMAT` for `GT`/`GQ`/`REPCN`;
  `##FILTER` for `notPeriodic`/`tooManyAlleles`/`lowDepth` (PASS is implicit); two
  `##` warning comments (§7); then the `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
  <sample_1> … <sample_N>` line.

**Contig source — a correction.** The task brief says "contigs from `merger.chrom_names()`
+ **catalog** lengths." The catalog is a lengthless TSV
(`#chrom start end motif purity_fraction ref_seq_start ref_seq`) — it has **no contig
lengths**. The lengths (and md5) live in the **`.ssr.psp` headers**, as
`ParsedChromosome { name, length, md5 }`, which the merger already reads at
`from_parts`. So the header's contigs come from a **new merger accessor over the psp
chromosome table**, not the catalog:

```rust
// on CohortMerger, alongside chrom_names()
pub(crate) fn chromosomes(&self) -> Vec<ParsedChromosome>   // name+length+md5, global-id order
pub(crate) fn sample_names(&self) -> Vec<String>            // merger/catalog order (= `labels`)
```

This mirrors the SNP `CohortMetadata { sample_names, contigs: Vec<ParsedChromosome> }`
exactly (the SNP path already sources contigs from psp headers), so the SSR header
writer can take the same two inputs. `sample_names()` exposes the existing private
`labels` field (the `#CHROM` sample columns must be in merger order — the same
single-source-of-truth rule the SNP `CohortMetadata` doc-comments call out).

> **Open C-1.** `labels` currently holds the **path string** passed to `open`
> (`path.display().to_string()`), not a clean sample name. For VCF columns we want a
> sample *name* (basename without `.ssr.psp`, or a `--sample-names` override). Decide:
> derive a basename in the accessor, or thread real names through `open`. Basename
> derivation is the cheap default; flag it.

---

## 6. The emit loop + recomputed candidates

The per-locus EM (`run_locus_em`) returns a `LocusCall` but **not** the
`CandidateSet` it built internally. Emission needs the candidates again (for
`is_variable`, `site_qual`, `format_vcf_record`). Two options: thread the candidates out
of the EM, or **recompute** `build_rungs` + `assemble_candidates` for the streamed locus
(what the test does). Recompute is cheap (the candidate set is tiny and pure) and keeps
the EM's return type unchanged — **recommend recompute** right where each locus is
genotyped in the streaming sweep. If profiling later shows it matters, returning the
candidates is an easy localized change. (Flag: recompute must use the **same
`RungCfg`/`CandidateCfg`** the EM used, or the allele indices won't line up — pass one
config struct through both.)

**Emit policy (confirm against spec §4.5 / genotyping §6):**

```
for each streamed CohortLocus (catalog order, via the bounded-queue workers):
    rungs   = build_rungs(locus, rung_cfg)
    cands   = assemble_candidates(locus, &rungs, ploidy, cand_cfg)
    call    = run_locus_em(locus, &rungs, &cands, &FROZEN params, ploidy, ...)   // §4
    apply_fp_control(locus, &mut call, &fp)
    emit =  cands.admit != PASS              // filtered → always emit with its reason
         || is_variable(&call, &cands)       // PASS + polymorphic → emit
    if !emit: continue                        // PASS + monomorphic → DROP
    qual = site_qual(&call, &cands, &fp)
    writeln!(out, format_vcf_record(chrom_name, locus, &cands, &call, qual))
    // writer reorders by locus seq → catalog-ordered VCF (reading §5)
```

The rule is **emit-iff-variable for PASS loci, but never silently drop a filtered
locus** — a `notPeriodic`/`tooManyAlleles`/`lowDepth` site is emitted with its FILTER
reason regardless of variability (genotyping §6, spec §4.5). Output stays in catalog
order because the streaming writer reorders finished loci by their monotonic locus
sequence number (reading §5), exactly as the SNP writer does.

---

## 7. Warnings; and the m2(a) fallback is now a hard error

- **Apparent-`F_IS` warning.** `vcf_out::f_is_warning(&calls.f_per_sample, &fp)`
  returns `Some(msg)` when the cohort mean `F` is implausibly high. Emit it as a `##`
  header comment *and* to stderr. The label must say **apparent `F_IS`** (it folds in
  population structure / mapping artifacts, not pure inbreeding — genotyping §6 / spec
  §4.5). This is the only header warning that survives decision E.

- **m2(a) chemistry fallback — dropped; it is the degenerate case of decision E's hard
  error.** m2(a) was: when the pre-pass finds **no confident genotypes cohort-wide**,
  fall back to coded literature defaults and warn. But "no confident genotypes
  cohort-wide" means **every** sample is absent from `grouped.group_of_sample` — so the
  per-sample hard error (§3, decision E) fires on the first sample *before* any
  fallback is reachable. The whole-cohort and per-sample cases collapse into one
  behaviour: **estimate chemistry or fail loud, naming the unresolved samples.** No
  literature-default path, no silent fallback. (If real data shows confident genotypes
  are routinely scarce, we revisit — possibly reinstating a *warned* literature-default
  mode behind an explicit flag — but not in v1.) This removes the "chemistry not
  estimated → literature defaults" header warning the brief asked for.

---

## 8. Decision D — ploidy; and threads

- **Ploidy.** `run_locus_em` asserts diploid. The driver **hard-codes ploidy 2** and
  documents the non-diploid panic as the contract until polyploidy lands. Confirm (the
  task lists this as decision D); it matches every test and the spec's diploid v1.
- **Threads.** `config.threads` maps to the rayon pool the parallel stages run under
  (the pre-pass over the subset, and the per-locus genotyping workers of the streaming
  sweep) — install a `rayon::ThreadPoolBuilder::new().num_threads(config.threads)` scope,
  as the byte-identity tests do. Determinism is guaranteed: the pre-pass reduces are
  integer/fixed-point, and each genotyped locus is a pure independent function of its
  reads + the frozen priors (no cross-locus reduce in the sweep — §4), so `--threads`
  changes speed, never output.
- **`config.queue_depth` is now live** — it is the bounded-queue depth of the streaming
  genotyping sweep ([ssr_call_reading.md §5](ssr_call_reading.md)), the peak-resident-loci
  RSS knob alongside `--block-window-bp`. (It was "reserved" only under the rejected
  materialize-all model.)

---

## 9. Decision B — fit the `G₀` decay in the pre-pass (settled 2026-06-23)

`G₀` gives each candidate allele a prior pseudocount `p^|Δ|` that decays with its
repeat-unit distance `Δ` from the locus's modal allele (a small-N regularizer + a
false-positive guard against far-out spurious candidates; the `G0_FLOOR` keeps a far
*real* allele recoverable). `p` is **fit per period** from the cohort, **not** a coded
default. This is a small pre-pass extension that lands **before/with** the driver
wiring — not deferred driver glue.

**What we measure (germline spread, not stutter — the load-bearing distinction).** `p`
encodes how far *real, inherited alleles* sit from the locus mode across the cohort.
That is the spread of **called alleles** (one entry per chromosome — a confident `8/8`
contributes length 8 twice, a confident `6/10` contributes 6 and 10 once each), **never**
the spread of *reads* (read scatter is stutter, which the stutter model already owns —
fitting `p` off reads would re-measure stutter and make the prior fight it). The
pre-pass already gates confident genotypes (homozygotes ∪ separated hets) for ε / shape
/ level; `G₀` reuses the same calls, tabulating a quantity currently discarded.

**The accumulator (mirrors `SlipProfile`).** Inside `accumulate_locus` — which already
sees every present sample's confident genotype for the locus — tally chromosome counts
per allele length, take the modal length (deterministic tie-break on the smaller
length), and add each called allele's count-weighted folded distance `k = |length −
mode|` into a **per-period histogram** `C[k]`. Integer counts ⇒ order-independent
merge ⇒ byte-identical across thread counts, exactly like the slip profile.

**Only loci that vary feed the fit (settled).** A locus is included **iff it shows ≥2
distinct confident alleles cohort-wide**; monomorphic loci are skipped. Reason:
monomorphic loci contribute only `k = 0` and would drag `p → 0`, yielding a uselessly
tight prior that over-suppresses the rare real variants at the loci that *do* vary —
precisely where the prior matters most (weak data). So `p` means "*given* a locus
varies, how spread are its alleles." Included variable loci keep **all** their
allele-copies (the mode copies too — the mode/off-mode ratio is what `p` encodes).

**The fit (closed form).** With `K̄ = Σ_k k·C[k] / Σ_k C[k]` the mean distance-from-mode,
the MLE of the symmetric geometric `p^|Δ|` is

```
p = ( √(1 + K̄²) − 1 ) / K̄          # K̄→0 ⇒ p→0 (tight); K̄ large ⇒ p→1 (flat)
```

— always in `(0,1)`, monotone in `K̄`. Same estimator *family* as the stutter-shape
decay (`(mean−1)/mean`), adapted for a two-sided distribution that includes the mode
(every `k≥1` has two alleles, mode±k). Symmetric only: a directional `G₀` is a larger
model change — noted, not built.

**Thin-period guard (default minimum, agreed).** Below a minimum included-allele-copy
count for a period, **fall back to the coded default** `p` — so tetra/penta loci with
too few confident variable alleles keep today's behaviour rather than fitting noise.
Because the fallback still emits a `p`, the pre-pass yields a `p` for **every** period
present, and the driver synthesizes nothing.

**Conservative-prior safeguards (agreed caveats).** (a) Clamp the fitted `p` so it can
never go *steeper* than a floor — a relatedness/structure-compressed cohort (the Wahlund
caveat that also haunts `F`) makes real spread look artificially tight and would
otherwise over-tighten the prior; `G0_FLOOR` still protects far real alleles. (b) It is
**empirical Bayes** (the prior is set from the same data it later regularizes) — mild
and standard, since `p` is fit from the *easy* confident/variable loci and only gently
shapes the hard ones; documented, not a blocker.

**Wiring + one alignment check.** `EstimatedParams` gains `g0_by_period: HashMap<u8,
G0PseudocountDecay>`; `estimate` computes it (variable-loci filter, closed form, min-count
fallback, clamp); `build_param_set` reads it straight into
`pseudocount_decay_per_loci_group`. Define "distance from the mode" at fit time the same
way `g0_pseudocounts` defines it at call time (it uses `rungs.modal_length()`) so the
two notions of "the locus's most common allele" stay consistent.

**Tests.** Recovery (simulate germline alleles drawn around a mode with known spread →
fitted `p` ≈ the formula); thin-period falls back to the default; variable-only filter
(adding monomorphic loci does not move `p`); byte-identical across thread counts.

---

## 10. Proposed driver shape + module layout

Keep it in `driver.rs`; lift the pure helpers so they're testable:

```
src/ssr/cohort/
  driver.rs    # run(SsrCallConfig): open → pre-pass+freeze → header → stream-genotype-emit → flush
               #   build_param_set() -> Result (§3, pure, unit-tested; hard-errors on unresolved samples)
               #   emit decision + format      (§6, pure per-locus, unit-tested)
               #   apparent-F_IS warning        (§7)
  prepass.rs   # + the burn-in: iterate ε/θ_period/level/G₀/F over the SUBSET → freeze  (§4);
               #   + per-period G₀ accumulator + g0_by_period fit  (§9, decision B)
  sample_groups.rs # + estimate+freeze F_i in the burn-in  (§4; reduce_f logic moves here)
  em.rs        # Step 2: per-locus θ_locus M-step + shrunk per-locus rate refit (LOCAL)  (§4)
  inbreeding.rs# run_cohort_em's whole-cohort outer loop RETIRED for genotyping (§4);
               #   genotyping = stream loci → run_locus_em on frozen params → emit
  vcf_out.rs   # + write_vcf_header(out, &chromosomes, &sample_names, &warnings)  (§5)
  merge.rs     # + chromosomes() + sample_names() accessors  (§5)
```

`run` orchestrates; `build_param_set`, the emit decision, and the header writer are
pure and individually tested. The genotyping sweep calls the existing per-locus
`run_locus_em` on each streamed `CohortLocus` under the frozen `ParamSet` — there is no
whole-cohort outer loop in the sweep (§4). The driver test (definition of done) runs the
**real** `run()` over a simulated-cohort `.ssr.psp` set to a VCF and asserts header +
records (reuse `sim.rs` to write real psp inputs, or assemble via `test_support` as the
existing driver tests do).

`SsrCallConfig` keeps its fields; `threads` becomes live, `queue_depth` becomes the
streaming-sweep bounded-queue depth (§8). `SsrCallError` gains **one** variant:
`UnresolvedSamples` (or similar) carrying the names of samples the pre-pass produced no
confident genotype for (§3/§7, decision E) — the only new non-`io` failure. Emit
failures stay `io::Error`; the merge errors stay `SsrMergeError`.

---

## 11. Open items / decisions to resolve

| # | decision | resolution |
|---|---|---|
| **A** | memory model | **two-pass streaming** (settled 2026-06-23, §4): materialization rejected. Burn-in iterates ε/θ_period/level/G₀/`F` over a stratified **subset** → **freeze the pooled params** → **single streaming genotyping sweep** of independent per-locus EMs. Per-locus stutter (shape + rate) refined **locally**, shrunk toward the frozen priors. `F` frozen after the pre-pass. Step 1 freezes everything (ships the VCF); Step 2 adds the local per-locus stutter refit. |
| **B** | `G₀` decay: coded default vs fit in pre-pass | **fit it** (settled 2026-06-23, §9): per-period closed-form geometric MLE on confident germline allele spread, **variable loci only**, min-count fallback to the default per period, clamp against over-tightening. |
| **C** | VCF writer: dedicated SSR text header vs adapt SNP `writer.rs` | **dedicated text header** in `vcf_out.rs`; contigs from a new `merger.chromosomes()` over psp headers (**not** the catalog). |
| **C-1** | sample names: psp path string vs basename vs `--sample-names` | **basename** default; flag. |
| **D** | ploidy | **hard-code 2**, document the non-diploid panic as the contract. |
| **E** | sample with no confident genotype | **hard error** naming the sample(s) (settled 2026-06-23) — subsumes the m2(a) cohort-wide fallback, which is dropped. Revisit only if real data shows it is common. |
| — | candidates for emit: recompute vs thread out of the EM | **recompute** per streamed locus in the emit path (cheap, pure); share one cfg. |
| — | pre-pass subset | the burn-in runs on a **stratified subset** held in RAM (§4); selection strategy is the parameters doc's (Q-R5). |
| — | `queue_depth` | **live** — the streaming-sweep bounded-queue depth (§8). |

**Settled upstream, not re-opened here:** emit-iff-variable for PASS / filtered-always-
emitted (spec §4.5, genotyping §6); apparent-`F_IS` labelling (M5); byte-identity across
thread counts (pre-pass reduces + per-locus EM, already tested); the re-read (not cache)
decision for the two passes (reading §8).

**Out of scope (do not do here):** calibrating the `dev_default` constants (incl. the
`G₀` fallback `p` and clamp floor); the deferred refinements (exact-AF QUAL kernel,
beta-binomial, polyploidy). Note the **per-locus stutter rate/shape refit is now in
scope as Step 2** (§4), not deferred; the *pooled* in-genotyping level refit it replaces
is dropped (the level line is frozen). No fabricated calibration numbers.
