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

**Second goal — keep the partially-covering reads.** A read whose alignment anchors one tract
border but runs off the read's end before the other **still carries evidence**: it cannot fix the
repeat count, but it proves a **lower bound** on it ("at least *k* units"). Production discards
these (they become a `BorderOffEnd` non-observation). Keeping them is strictly more evidence at
the same depth — which matters most exactly where we are weakest: lone-carrier hets at ~3
reads/plant. So the preparer emits a **censored** observation for them rather than dropping them
(§3).

**The prior art is GangSTR, not HipSTR.** Worth stating precisely, because an earlier draft of
this spec got it wrong and the mistake would have sent an implementer to the wrong source:

- **GangSTR** ([Mousavi et al. 2019]) is the model to follow. It has an explicit **bounding-read**
  class (`RC_BOUND`) whose likelihood (`FlankingClass`) consumes exactly this lower bound. §6
  records the two details that make it work.
- **HipSTR** does **not** do this. Its allele set comes only from fully-enclosing reads, and a
  read without a high-confidence seed outside every stutter block receives a *flat, all-zero
  likelihood row* — it is not censored evidence, it is **no** evidence. This is why HipSTR cannot
  call alleles longer than a read. Do not look here for partial-read handling.
- **freebayes** uses partially-covering reads, but by a **different mechanism**: it splits a
  partial's support `1/N` across the N candidate haplotypes it cannot distinguish (an *ambiguity*
  weight, independent of how much was covered) — not a length bound. It is also not STR-aware.
  Our `LowerBound` payload is sufficient to express freebayes' rule too (§3), but the rule we are
  implementing is GangSTR's.

> **This is the one place the STR path is *not* a port.** Everything else here reproduces
> production's Stage-1 pipeline; partial-read bounds are new behaviour with **no production
> oracle**, which scopes the parity test (§7) and ripples downstream (§6). Called out because a
> spec that quietly mixes "port" and "improvement" makes the parity test a lie.

**Non-goals** (beyond the shared preamble's): the delimiter is a **ruler, not a scorer** — it
finds the read's tract boundaries; it does **not** score candidate alleles. Scoring
(`Lr = P(read | allele)`) is step 7, a *separate* forward pair-HMM over the extracted bytes (§6).
Deduplicating observations into per-sequence counts is the STR **gatherer's** job (the tract
tally), not preparation's (§2).

**Alleles longer than a read are out of scope for this step — by construction, not oversight.**
Step 2 emits only **per-read facts**: what *this* read shows about the tract. GangSTR reaches
beyond read length with two classes that are *not* per-read facts and therefore cannot be
expressed here at all:

- **FRR** (fully-repetitive read — the read lies wholly inside the tract). Its datum is not a
  property of the read's own sequence: it is **the mate's bp offset from the tract**, scored
  against a sample-level **insert-size distribution** (a CDF-difference term), plus a Poisson
  **count** term needing locus coverage. An FRR read touches no border and observes only repeat
  bases, so in *any* per-read sequence-keyed schema it degenerates to the same constant for every
  FRR at the locus — all its signal lives in the mate.
- **Spanning pairs.** Same shape: the datum is an **insert size**, a pair-level quantity.

Both are **step 5** (read-class determination + fragment model), and this is the load-bearing
justification for the step-2/step-5 split: they are pair- and distribution-level, not per-read.
Step 2 therefore covers alleles up to roughly read length, plus correct one-sided constraints
beyond it. Reaching *past* read length is step 5's job, and it will need the mate and the insert
distribution that this step deliberately does not carry.

The output type is **`SsrTractObs`** and the consumer is the **tract tally** (preamble §2).

---

## 2. The transform — fetch → extract → delimit → gate

Ported from the production `ssr/pileup` Stage-1 pipeline, with one deliberate extension (step 3).
Per candidate read, in order:

1. **Footprint + extract.** Map the reference tract window onto the read's coordinates across any
   internal indels, centering a window on the tract; widen it when a long allele reaches past the
   window (long-allele recovery).
2. **Delimit** — a per-read **Viterbi (max-path) pair-HMM**. Align the extracted read slice
   against `left_flank + reference_tract + right_flank` with a **tract-aware affine gap** (gaps
   are cheap inside the tract, where length variation is expected, and dear in the flanks), and
   read the repeat region off the flank-junction columns.
3. **Classify by how many borders were anchored → observation.** This is where the extension
   lives, and note that the delimiter *already computes it* — the same alignment, read differently:
   - **both borders anchored** → the tract length is pinned: an **exact** observation (the
     verbatim tract bytes), exactly as production emits today;
   - **one border anchored, the read running off its end mid-tract** → the length is **censored
     below**: a **lower-bound** observation carrying the repeat the read did prove. Production
     calls this `BorderOffEnd` and drops it; we keep it (§1);
   - **no border anchored** (the read lies wholly inside the tract) → still `None` here. Such a
     read bounds the allele only via an insert-size fragment model, which is **step 5**, not a
     per-read fact the delimiter can state.
4. **Quality-gate.** A delimited read that fails the base-quality gate is a tallied
   non-observation (§4) → `None`, whichever class it fell in.

**Jargon, once — delimitation.** *Delimiting* a read means finding where its repeat tract starts
and ends, so the bytes between can be read out as the observed allele. It is alignment used as a
*ruler*, not as a *score* — the distinct, later scoring HMM (step 7) is a forward sum, not a
Viterbi max-path.

---

## 3. The output type — `SsrTractObs`

```rust
/// The STR prepared read — the STR ReadPreparer's `type Prepared`. What one read *shows* about
/// the tract: either its exact sequence, or a lower bound when the read ran out mid-tract.
pub enum SsrTractObs {
    /// Both borders anchored: the tract length is pinned. The verbatim observed bytes — what
    /// production emits today (`ReadObs::Sequence`).
    Exact {
        tract: Box<[u8]>,
        /// true if the window was widened for a long allele (long-allele recovery)
        widened: bool,
    },
    /// One border anchored, the read running off its end mid-tract: the length is **censored
    /// below**. The read cannot fix the repeat count — it proves the allele is *at least* this
    /// long. Kept as bytes, not a unit count, for the same reason `Exact` is: interior
    /// interruptions are sequence facts a count would erase.
    LowerBound {
        /// The repeat the read did prove, verbatim (a prefix or suffix of the true tract).
        observed: Box<[u8]>,
        /// Which border was anchored — i.e. which side `observed` is flush with.
        anchored: TractBorder,   // Left | Right
    },
}
```

`SsrTractObs` is a reduced observation — bytes, not a decomposable read: the STR path has
*already* localised the evidence to the tract, so there is no downstream per-position
decomposition (contrast the generic `PreparedRead`, which the pileup walker still decomposes).

**Why bytes and not a minimum repeat count.** The lower bound *is* "at least *k* units", and `k`
falls straight out of `observed.len() / period` — but the spec keeps the bytes and lets the
consumer derive `k`. Storing the count instead would throw away the sequence, and the
interrupted-repeat work established that same-length interior substitutions are real, distinct
alleles: a count cannot represent them, and a length-keyed ladder collapses them. Bytes keep
`LowerBound` on the same sequence-keyed footing as `Exact`. (GangSTR reduces each bounding read
to a single `int32_t` copy number precisely because it works in *count* space; we do not.)

**Why keep `anchored` when GangSTR does not.** GangSTR discards the side — its pre-flank and
post-flank reads both collapse into one bounding class, and its placement term integrates over
both orientations. That is sound **in count space**, where a prefix and a suffix of *k* units are
the same fact. It does not transfer to us: in **sequence space** a prefix and a suffix of the
tract are *different constraints* — an interrupted repeat's prefix is not its suffix, which is
the whole reason we key alleles by sequence. `anchored` is also exactly what a freebayes-style
prefix/suffix ambiguity match would need (§1). One field, two independent justifications; it
stays.

The STR implementation:

```rust
pub struct SsrDelimitPreparer<Canon: RefSeq> {
    canon_ref: Canon,             // canonical bases for the flank+tract+flank alignment
    scratch: ViterbiScratch,      // reused per-worker DP matrices, not per-read allocation
    /* config: quality gate, widen policy */
}

impl<Canon: RefSeq> ReadPreparer for SsrDelimitPreparer<Canon> {
    /// Unlike the generic path, the STR delimiter **does** need the locus: its motif, borders
    /// and flanks are what make the gap penalties tract-aware.
    type Locus = SsrLocus;
    type Prepared = SsrTractObs;
    fn prepare_read(&self, read: &MappedRead, locus: &SsrLocus) -> Option<SsrTractObs> {
        /* footprint → delimit → gate */
    }
}
```

The `SsrLocus` (from the `ssr-catalog`, via the router) carries the **tract structure** — motif,
borders, flank sequences. This is the one asymmetry with the generic path, whose `type Locus` is
`()`: production's `delimit_read(region_seq, region_qual, locus, model, scratch)` takes the locus
for exactly this reason, while `process_read` takes none. The reference accessor is a **field**,
as on the generic path (preamble §3) — no window is passed in.

---

## 4. Reference dependency and the "no observation" reasons

**Reference.** The delimiter aligns against `left_flank + reference_tract + right_flank`, all
**canonical** bases (`RefSeq`), plus the motif/borders carried on the routed locus. No raw-byte
path is needed on the STR side.

**`None` reasons (tallied):** `NoBorderAnchored` — the read lies wholly inside the tract, so the
delimiter can state no per-read fact (only a fragment model could, which is step 5);
`LowQuality` — the extracted window fails the base-quality gate; `WindowTruncated` — the
reference/read window was truncated. Each is a tallied per-read outcome, not a run error
(preamble §5).

**Production's `BorderOffEnd` is deliberately *not* a `None` reason here.** That is exactly the
one-border-anchored case, which we now keep as `SsrTractObs::LowerBound` (§1, §3). The counts
should still record how many observations were bounded rather than exact — a run whose STR
evidence is mostly lower bounds is telling you something (short reads vs. long tracts), and
"no silent caps" applies to *kept* reads too, not only dropped ones.

**Read selection (upstream of prep).** Reads are fetched per STR locus, and the fetch gate must
now admit the partially-covering reads it previously excluded — otherwise the `LowerBound` class
is unreachable and this extension is dead code. The gate that remains is *relevance* (does the
read touch this tract at all), not *spanning*. This is a change to the production fetch stage
(`fetch_locus_reads`), not just to the delimiter — flagged here because it is easy to miss: the
new observation class is worthless if selection still drops its inputs upstream.

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

- **Tract tally (dedup into per-sequence counts) → the STR gatherer.** `Exact` observations
  aggregate into a `Vec<(tract_bytes, depth)>` — the STR-path counterpart of the pileup walker,
  producing `SampleLocusObservations`. Aggregation across reads is not a per-read preparation step.
  **`LowerBound` observations do not dedup the same way** and the gatherer must say how it carries
  them: a bound of "≥ k units" is compatible with *many* alleles, so it cannot collapse into the
  exact-sequence count table without losing its meaning. Deciding that representation is the
  gatherer's job, but it is a **direct consequence of this step's §3** — flagged so it is not
  discovered late.
- **Read likelihood `Lr = P(read | allele)` → step 7 — and it must learn censoring.** The scoring
  HMM (`ReadLikelihoodModel`, Model A) is a *separate* **forward-sum** pair-HMM over the tract
  bytes, cleanly split from this step's **Viterbi** delimiter. But `Lr` for a `LowerBound` is
  **not** the exact-match likelihood. Model A scores `Lr(obs | cand)` assuming `obs` is a
  *complete* tract; feed it a prefix and it scores the bound as a **short allele** — actively
  misleading, not merely imprecise. This is the largest ripple of §1, and the step-7 spec must
  resolve it before `LowerBound` is consumed (§7). Three specifics from GangSTR's `FlankingClass`
  that a from-scratch attempt would plausibly miss:

  1. **The hard inequality.** `allele >= bound` or the term is zero (`NEG_INF`). One-sided
     constraint, nothing more.
  2. **The `1/allele` dilution — the load-bearing part.** When the constraint *holds*, the term is
     **not** a point mass: it is a placement probability times `1/allele`. A longer allele makes
     any *specific* partial observation less likely, so the likelihood is **not monotonically
     increasing in allele size**. Without this, "allele ≥ k" makes every allele ≥ k *equally*
     good — the bound would only permit, never discriminate, and would bias calls upward. A
     censored term that merely truncates is worse than useless: it is a systematic inflation.
  3. **The conservative discount.** GangSTR bounds at `nCopy − 1`, deliberately, because the
     terminal partial unit is unreliable ("flanking is always picked up +1"). And it
     **outlier-filters** bounds before use (drop bounds above ~3× the median) — a wild bound
     otherwise drags the ceiling arbitrarily.

  Implement from the **paper** (Mousavi et al. 2019), not from the vendored source: GangSTR is
  **GPL-3-or-later** and HipSTR **GPL-2**, so a transliteration would pull copyleft onto the
  caller — the same constraint the TRF-mod work ran into. (freebayes is MIT, if its 1/N ambiguity
  rule is ever wanted.)
- **Read-class / spanning + FRR with the fragment model → step 5.** Reaching lengths beyond the
  read is a distinct per-locus step, and §1 records *why* it cannot be step 2: those classes are
  pair- and distribution-level, not per-read. When step 5 is specced it will need, per GangSTR:
  the **mate's offset** from the tract (FRR) or the **insert size** (spanning); a sample-level
  **empirical insert-size distribution** (its FRR term is a CDF difference over it); and **locus
  coverage** for the Poisson FRR-count term that is the actual expansion signal. None of those
  are facts about the read this step prepares. The wholly-inside-the-tract read is `None` here
  for exactly that reason (§4).

**Alternatives to bench (recorded, not built now):** GangSTR-style **realign every read with
Striped Smith-Waterman** against a synthetic repeat reference — an alternative *delimiter*
producing an `SsrTractObs`-shaped output, so the tally consumes it unchanged. (HipSTR's
read↔haplotype HMM is a *scoring* model — it belongs to step 7's likelihood bake-off, not to this
delimitation step; noting it here only to place it correctly.)

---

## 7. Reuse over rewrite — the map to production

The parity oracle is the production extracted tract — but it **only covers the `Exact` class**. A
ported STR impl is correct when, on a fixture, every `SsrTractObs::Exact { tract, .. }` matches
the production `ReadObs::Sequence` payload byte for byte, **and** every read production called
`BorderOffEnd` now yields a `LowerBound` instead of vanishing (a count-level check: nothing
production kept may be lost, and the new class must be non-empty).

**`LowerBound` itself has no production oracle** — it is new behaviour (§1), so parity cannot
validate it. It needs a different standard: the `benchmarks/ssr_tomato1` harness and
HipSTR-concordance, which can answer the question that actually matters — *does keeping the
partial reads improve calls at low coverage, or just add noise?* Until step 7 handles censored
observations (§6), `LowerBound` should be **emitted and tallied but not yet consumed**, so this
step can land and be measured without a half-finished likelihood silently mis-scoring bounds as
short alleles.

| what | existing code | ng reuse |
|---|---|---|
| fetch reads (per STR locus) | `fetch_locus_reads` ([ssr/pileup/fetch_reads.rs](../../../../src/ssr/pileup/fetch_reads.rs)) | the read-selection gate upstream of `prepare_read` — **must be widened**: its spanning/flank gate currently drops the partially-covering reads `LowerBound` depends on (§4) |
| footprint + extract/center | `read_footprint` / `extract_region` / `widen_region` ([ssr/pileup/footprint.rs](../../../../src/ssr/pileup/footprint.rs)) | call directly |
| delimit (Viterbi) | `delimit_read` ([ssr/pileup/alignment.rs](../../../../src/ssr/pileup/alignment.rs)), `ViterbiScratch` | model for `SsrDelimitPreparer` |
| the observation — exact | `ReadObs::Sequence` ([ssr/pileup/locus_tally.rs](../../../../src/ssr/pileup/locus_tally.rs)) | model for `SsrTractObs::Exact` — the parity target (§7) |
| the observation — bounded | `ReadObs::BorderOffEnd` (a **drop** today, same file) | **not a port**: production discards it; we keep it as `SsrTractObs::LowerBound` (§1). No production oracle (§7) |
| reference | `RefSeq` (canonical flanks + tract) ([ref_seq.md](ref_seq.md)) | reuse as-is |

---

## 8. Open questions

- **`SsrLocus` shape.** The delimiter needs motif, borders, and flank sequences from the routed
  locus; the exact carrier (`type Locus = SsrLocus`) is co-owned with the router / STR-catalog
  spec (`LocusKind`). *Settled here:* it is the **locus** that is passed, not a reference window
  — the accessor is a field on the preparer (preamble §3); the earlier `LocusWindow` idea is
  dropped on both paths.
- **Widening policy.** How far to widen the window for long alleles (long-allele recovery) is a
  tunable inherited from production; whether it stays a fixed policy or a config knob is open.
- **Does keeping the partial reads actually pay? — open, and empirical.** The premise (§1) is that
  a lower bound is real evidence and that more evidence helps most exactly where we are weakest:
  lone-carrier hets at ~3 reads/plant. But a bound is *weak* evidence, and it is not free — it
  needs a censored likelihood (§6). It could plausibly buy little over the exact reads, or add
  noise. Measure on `benchmarks/ssr_tomato1` (silver recall + HipSTR concordance) with the class
  on vs. off before committing. Note the regime difference that makes this genuinely open:
  **GangSTR's bounding class exists to serve pathogenic expansions far beyond read length**,
  whereas our catalog is period 2–6 with alleles mostly *inside* a read — where the exact class
  already fires. The mechanism transfers; the *payoff* may not. **This question decides whether
  the extension stays.**
- **How much censoring does step 7 need — blocking `LowerBound`'s consumption.** `LowerBound`
  cannot be *used* until the likelihood models it (§6): the hard inequality, the `1/allele`
  dilution, the conservative discount. Until then it is emitted and tallied but not consumed.
  Whether our sequence-keyed frame can reuse GangSTR's count-space placement term unchanged, or
  needs a sequence-aware analogue, is step 7's call and gates switching this on.
- **Which border-anchoring cases are real in our data.** The design distinguishes one-border
  (`LowerBound`) from no-border (FRR → `None`, needing step 5's mate + insert distribution).
  Whether the no-border class is common enough at our read length and tract lengths to justify
  step 5's pair-level machinery is unmeasured — count both classes while measuring the above,
  since that count is what decides whether step 5's fragment model is worth building at all.
- **Where STR prep is *invoked* — leaning: per STR locus, by the STR locus processor**, which
  fetches spanning reads and calls `prepare_read` on each (the compose-not-subsume resolution,
  preamble §2). Confirm the call site when the STR gatherer is specced.
