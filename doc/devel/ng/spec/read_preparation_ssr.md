# ng step 2 (STR path) — read preparation: tract extraction

*Status: design spec, 2026-07-14. The **STR path** of step 2. Inherits the shared contract and
discipline from the preamble [`read_preparation.md`](read_preparation.md) — read that first;
this spec covers only what is specific to the STR path. Grounded in the production `ssr/pileup`
Stage-1 pipeline ([src/ssr/pileup/](../../../../src/ssr/pileup/)). **No code yet.** Naming:
**STR** in prose, `ssr` in code — the module is `src/ng/read/` (the impl file `pair_hmm.rs`) and
the STR types carry the `Ssr` prefix.*

---

## 1. Scope — goals and non-goals (path-specific)

**Goal:** prepare an STR read by **extracting its repeat tract** — center a window on the tract,
delimit the read's repeat region against the reference with a per-read pair-HMM, and read out the
observed tract bytes. The output is the observed sequence between the tract borders.

**Non-goals** (beyond the shared preamble's): the delimiter is a **ruler, not a scorer** — it
finds the read's tract boundaries; it does **not** score candidate alleles. Scoring
(`Lr = P(read | allele)`) is step 7, a *separate* forward pair-HMM over the extracted bytes (§6).
Deduplicating identical tracts into per-sequence counts is the STR **gatherer's** job (the tract
tally), not preparation's (§2).

The output type is **`SsrTractObs`** and the consumer is the **tract tally** (preamble §2).

---

## 2. The transform — fetch → extract → delimit → gate

Ported from the production `ssr/pileup` Stage-1 pipeline. Per spanning read, in order:

1. **Footprint + extract.** Map the reference tract window onto the read's coordinates across any
   internal indels, centering a window on the tract; widen it when a long allele reaches past the
   window (long-allele recovery).
2. **Delimit** — a per-read **Viterbi (max-path) pair-HMM**. Align the extracted read slice
   against `left_flank + reference_tract + right_flank` with a **tract-aware affine gap** (gaps
   are cheap inside the tract, where length variation is expected, and dear in the flanks), and
   read the repeat region off the two flank-junction columns.
3. **Quality-gate → observation.** Classify the delimited read: its extracted tract bytes if it
   clears the base-quality gate, else a tallied non-observation (§4) → `None`.

**Jargon, once — delimitation.** *Delimiting* a read means finding where its repeat tract starts
and ends, so the bytes between can be read out as the observed allele. It is alignment used as a
*ruler*, not as a *score* — the distinct, later scoring HMM (step 7) is a forward sum, not a
Viterbi max-path.

---

## 3. The output type — `SsrTractObs`

```rust
/// The STR prepared read — the output of the STR ReadPrep impl (`type Prepared`). The observed
/// repeat-tract bytes, ready for the tract tally to deduplicate into per-sequence counts.
pub struct SsrTractObs {
    pub tract: Box<[u8]>,     // the extracted repeat tract (verbatim observed bytes)
    pub widened: bool,        // true if the window was widened for a long allele (long-allele recovery)
}
```

`SsrTractObs` is a reduced observation — a byte sequence — not a decomposable read: the STR path
has *already* localised the evidence to the tract, so there is no downstream per-position
decomposition (contrast the generic `PreparedReadNg`, which the pileup walker still decomposes).

The STR implementation:

```rust
pub struct SsrDelimitPrep { /* config: quality gate, widen policy */ }
impl ReadPrep for SsrDelimitPrep {
    type Prepared = SsrTractObs;
    fn prepare_read(&self, read: &MappedRead, window: &LocusWindow) -> Option<SsrTractObs> { /* footprint → delimit → gate */ }
}
```

The delimiter needs the **tract structure** (motif, borders, flank sequences) — the `LocusWindow`
on the STR path carries the routed STR locus (from the `ssr-catalog`), not just reference bases.

---

## 4. Reference dependency and the "no observation" reasons

**Reference.** The delimiter aligns against `left_flank + reference_tract + right_flank`, all
**canonical** bases (`RefSeq`), plus the motif/borders carried on the routed locus. No raw-byte
path is needed on the STR side.

**`None` reasons (tallied):** `BorderOffEnd` — the tract border falls off the end of the read
(the read spans only part of the tract); `LowQuality` — the extracted window fails the base-quality
gate; `WindowTruncated` — the reference/read window was truncated. Each is a tallied per-read
outcome, not a run error (preamble §5).

**Read selection (upstream of prep).** Only reads that *span* the tract (with adequate flanking)
are prepared; the spanning/flank gate is part of fetching reads per STR locus. Non-spanning reads
never reach `prepare_read` — they are not filtered here so much as *not selected* for this locus.
(Whether flank-only and partial-spanning reads inform length via an insert-size fragment model is
**step 5**, read-class determination — not step 2.)

---

## 5. Cross-cutting concerns

- **Performance / parallelism.** Pairwise-independence (preamble §3) makes STR prep parallel per
  read. The Viterbi delimiter is the cost; reuse a per-worker `ViterbiScratch` for its matrices
  rather than allocating per read (the project's scratch-buffer discipline — the production
  delimiter already threads a reusable scratch).
- **Determinism.** Extraction depends only on the read and the locus, so a read delimits to the
  same `SsrTractObs` regardless of thread interleaving.

---

## 6. Deferred, with a recommended home

- **Tract tally (dedup into per-sequence counts) → the STR gatherer.** Identical `SsrTractObs`
  tracts are aggregated into a `Vec<(tract_bytes, depth)>` — the STR-path counterpart of the
  pileup walker, producing `LocusEvidence`. Aggregation across reads is not a per-read
  preparation step.
- **Read likelihood `Lr = P(read | allele)` → step 7.** The scoring HMM (`ReadLikelihoodModel`,
  Model A = HipSTR) is a *separate* **forward-sum** pair-HMM over the tract bytes, cleanly split
  from this step's **Viterbi** delimiter. Keeping delimitation (ruler) apart from scoring (score)
  is what let the read-likelihood bake-off swap models independently (`ng_proposal.md`).
- **Read-class / spanning-with-fragment-model → step 5.** Reaching lengths beyond the read via an
  insert-size model is a distinct per-locus step, not read preparation.

**Alternatives to bench (recorded, not built now):** GangSTR-style **realign every read with
Striped Smith-Waterman** against a synthetic repeat reference — an alternative *delimiter*
producing an `SsrTractObs`-shaped output, so the tally consumes it unchanged. (HipSTR's
read↔haplotype HMM is a *scoring* model — it belongs to step 7's likelihood bake-off, not to this
delimitation step; noting it here only to place it correctly.)

---

## 7. Reuse over rewrite — the map to production

The parity oracle is the production extracted tract: a ported STR impl is correct when its
`SsrTractObs.tract` matches the production `ReadObs::Sequence` payload on a fixture.

| what | existing code | ng reuse |
|---|---|---|
| fetch spanning reads (per STR locus) | `fetch_locus_reads` ([ssr/pileup/fetch_reads.rs](../../../../src/ssr/pileup/fetch_reads.rs)) | the read-selection gate upstream of `prepare_read` |
| footprint + extract/center | `read_footprint` / `extract_region` / `widen_region` ([ssr/pileup/footprint.rs](../../../../src/ssr/pileup/footprint.rs)) | call directly |
| delimit (Viterbi) | `delimit_read` ([ssr/pileup/alignment.rs](../../../../src/ssr/pileup/alignment.rs)), `ViterbiScratch` | model for `SsrDelimitPrep` |
| the observation | `ReadObs::Sequence` ([ssr/pileup/locus_tally.rs](../../../../src/ssr/pileup/locus_tally.rs)) | model for `SsrTractObs` |
| reference | `RefSeq` (canonical flanks + tract) ([ref_seq.md](ref_seq.md)) | reuse as-is |

---

## 8. Open questions

- **`LocusWindow` / STR locus shape.** The delimiter needs motif, borders, and flank sequences
  from the routed locus; the exact carrier is co-owned with the router / STR-catalog spec.
- **Widening policy.** How far to widen the window for long alleles (long-allele recovery) is a
  tunable inherited from production; whether it stays a fixed policy or a config knob is open.
- **Where STR prep is *invoked* — leaning: per STR locus, by the STR locus processor**, which
  fetches spanning reads and calls `prepare_read` on each (the compose-not-subsume resolution,
  preamble §2). Confirm the call site when the STR gatherer is specced.
