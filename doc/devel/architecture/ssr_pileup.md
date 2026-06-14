# Stage 1 — `ssr-pileup` (per-sample evidence extraction)

**Status:** settled design, 2026-06-13 (worked through from an initial discussion
draft). The second *stage* pass, after the shared types
([ssr_shared_types.md](ssr_shared_types.md)) and Stage 0
([ssr_catalog.md](ssr_catalog.md)); follows the overall
[architecture](ssr_genotyping_architecture.md) (§8 module 2) and the spec
([ssr_genotyping.md](../specs/ssr_genotyping.md) §4). The spec settles *what* the
evidence is and *why* (the two-tier likelihood, the stutter-free invariant, the
columns); this designs *how the extractor is structured* — the read pipeline, the
pair-HMM, the parallelism, the writer. **Every structural question is decided
(§14)**, grounded where it matters in the vendored HipSTR/GangSTR + Dindel sources;
what remains is *data-validation of named constants*, not design — ready to spin
into an implementation plan.

This is also where two parameters the earlier passes deferred land: **`flank_bp`**
(catalog §5 / [ssr_catalog.md](ssr_catalog.md) §10.1 explicitly punted its *value*
to this stage) and the **Stage-1 parallelism shape** (architecture §7/§9 open item
1). Both are owned here.

---

## 1. What Stage 1 does

`ssr-pileup` takes **one sample's BAM/CRAM + the catalog** and emits **one
`.ssr.psp` evidence file** — a sparse, columnar, per-locus summary of what the
reads say about each locus's repeat length, **with stutter deliberately left
out** (it is learned and applied in Stage 2). It is the SSR analog of the SNP
`pileup`: heavy per-sample read work, done once, summarised to a columnar
artifact a light cohort stage then merges.

```
 BAM/CRAM ─┐                      ┌──────────────── per locus ─────────────────┐
 (1 sample)│   pull reads over    │  fast path: clean spanning reads → count L  │
 catalog ──┼──► [locus ± flank] ─►│  slow path: ambiguous/soft-clip → pair-HMM  │──► sampleN.ssr.psp
 (loci +   │   (triage)           │  ladder:    on-/off-ladder candidates       │    (sparse, block-gridded,
  ref_seq) ┘                      │  aggregate: histogram + CSR + offl_*         │     stutter-free likelihoods)
                                  └─────────────────────────────────────────────┘
```

Three framings from the spec drive every choice below:

- **Stutter-free (§4.2).** Stage 1 scores reads under *sequencing/alignment error
  only*. This is *the* architectural simplifier: with no stutter and no cross-read
  haplotype assembly, **every locus is an independent unit of work** — no cohort
  context, no per-position streaming state machine, no ordering dependency between
  loci. Stage 1 is far more embarrassingly parallel than the SNP pileup (§8).
- **Two tiers, two storage regimes (§4.2/§4.3).** A read produces *either* one
  confident length (collapses into a locus histogram, not stored per-read) *or* a
  short sparse `(length, log-lik)` distribution (stored per-read in CSR columns).
  The fast/slow code split and the histogram/CSR storage split are the same split.
- **The catalog is the only reference (§3.2).** Flanks, the motif, the reference
  tract — every reference base the SSR math needs comes from the catalog's
  embedded `ref_bytes` ([types.rs](../../src/ssr/types.rs) `Locus::ref_bytes`). A
  FASTA is opened **only to decode CRAM** (the noodles layer), never by the SSR
  algorithm.

---

## 2. The two-tier likelihood — the core algorithm

> **The capillary-trace intuition (read this first if the math below is opaque).**
> Genotyping a microsatellite on a capillary, you never see one clean peak — you
> see a **main peak at the true allele** plus smaller **stutter peaks** a repeat
> or two shorter, and your trained eye says "the small 11 peak is just stutter,
> the real allele is 12." That judgement is really **two steps**: first you *read
> the peaks off the trace* (record what's physically there), then you *interpret*
> them (decide which are stutter). Our two stages are exactly those two steps.
>
> - **Stage 1 (`ssr-pileup`) = reading the peaks.** It looks at one read at a time
>   and asks only *"what length molecule does this read look like it came from?"* —
>   accounting for **sequencing error** (a misread base makes a messy read a bit
>   unsure: "70% an 11, 30% a 12"; a clean read just says "12"). It records the
>   stutter peaks too — it simply does **not decide** which peaks are stutter. It
>   never uses any stutter *rate*. That is all "stutter-free" means.
> - **Stage 2 (`ssr-call`) = interpreting the trace.** It pools every read's
>   measurements across the cohort and decides the true alleles, using stutter
>   rates it **learns from the whole population** ("that 11 shoulder is slippage
>   off the 12 → genotype 12/12").
>
> Two kinds of noise, kept in separate stages: **stutter** is PCR slippage that
> happens *in the tube before sequencing* (it makes a physically-11 molecule from
> a true 12) — Stage 2's concern; **sequencing error** is the machine misreading
> bases — Stage 1's concern. Stage 1 measures the molecules it sees; Stage 2
> explains them as *true allele + stutter*.
>
> **Why split it?** Reading each read carefully is expensive and only has to be
> done **once**; the interpretation is **repeated many times** (Stage 2 tries many
> stutter rates to find the best cohort fit). Measure once, save, reuse — that is
> the whole reason for two stages. (The notation below makes this precise; the
> trace picture is the same thing in words.)

**The unit of work is one spanning read at one locus, and the question asked of
it is: which repeat-allele length(s) could this read have come from, and with what
likelihood?** The answer is `Qᵣ(L)` — the probability the read's observed bases
arose from a **molecule whose tract is `L` repeat units long**, under
sequencing/alignment error. `Qᵣ` does not attribute that length to stutter: a
length-`L` molecule may be a faithful copy of an `L`-unit allele *or* a stutter
product of a different allele, and `Qᵣ` is blind to which — that attribution is
Stage 2's stutter step. So `Qᵣ(L)` is not the whole story; it is **one factor** of
the full per-read likelihood `P(read | allele a, θ) = Σ_L Qᵣ(L)·S_θ(L | a)` (spec
§5.1) — the factor Stage 1 owns. That per-read function *is* this stage's
algorithm; everything here exists to compute it and aggregate it per locus.

Concretely, the algorithm **acts on two sequences**: the read's bases (including
soft-clipped ones, §3) and a set of **candidate haplotypes** reconstructed from
the catalog — `H_L = left_flank + (motif × L) + right_flank`, built from
`Locus::ref_bytes` + `motif` for each plausible repeat count `L`. Scoring the read
against `H_L` for a range of `L` yields `Qᵣ(L)`: a point mass when the read pins
one length confidently, a spread when several lengths are plausible.

**The plausible `L` are not a fixed global ladder — they are a per-read window
around the read's own observed length.** Each read first yields a *single*
observed count: the direct motif tally for a clean read, or the
**soft-clip-recovered** count for a clipped one (§3.2 — the realignment of the
clip establishes the true, possibly large, count). A *confident* read stops there
— its `Qᵣ` is a point mass at that count, stored as a single length (§11), no
window. An *ambiguous* read instead has its forward evaluated over the window
`L ∈ [count − W, count + W]` (`W = STUTTER_WINDOW_UNITS`, §5), because under
sequencing error alone several nearby lengths may genuinely fit; only the lengths
with real mass are then kept (§11). Either way a read is scored only against rungs
*near what it looks like*, never the locus's whole allele range — which is what
keeps the pair-HMM band narrow and the stored profile sparse.

> **"Stutter-free" means stutter-*parameter*-free, not stutter-unaware — the
> distinction matters.** Two things make Stage 1 look entangled with stutter, and
> both are real: the candidate window's width is named `STUTTER_WINDOW_UNITS`, and
> `Qᵣ(L)` is a *factor* of the stutter-aware likelihood `Σ_L Qᵣ(L)·S_θ(L | a)`. So
> Stage 1 genuinely computes *part* of the stutter calculation. What is absent
> from `Qᵣ` is any **stutter parameter** (`u, d, ρ`): each value conditions on a
> realized molecule length `L` and scores sequencing error only; the kernel `S_θ`
> and the sum over `L` are Stage 2's. The window then merely decides *which* `L`
> get such a value — sized so Stage 2's convolution has support points where
> stutter can reach. **The payoff of factoring it this way is θ-independence:**
> Stage 2's EM re-fits `θ` every iteration, and because `Qᵣ` does not depend on
> `θ`, the expensive pair-HMM runs **once** per read in Stage 1 and is reused
> across all EM iterations — the per-sample/cohort split (§8) rests on exactly
> this. (Spec §4.2 / §5.1.)

**The two tiers are not two biological classes of read — they are two ways to
*reach the same `Qᵣ`*, routed by how ambiguous the length read-out is.** The
classification the gate performs is therefore epistemic, not biological: *can the
length be counted off the read with confidence, or must it be inferred as a
distribution?* The end in both cases is identical — per-read length evidence for
Stage 2 — but the cheap, confident majority is handled by direct counting and the
ambiguous minority by full probabilistic realignment, because forcing a single
length onto an ambiguous read would fabricate certainty the data don't support
(precision-first, §1.2 of the spec). Restating the shape the code must take (spec
§4.2 is the contract):

| | **fast path** | **slow path** |
|---|---|---|
| **who** | read spans tract, both flanks cleanly aligned (≥ `MIN_FLANK_BP`, no soft-clip), tract is a **pure integer tiling**, boundary `Q ≥ MIN_BASE_QUAL` | everything else — interior indel/mismatch, impure tract, low boundary Q, **and every soft-clip-recovered read** |
| **work** | count motif copies in O(read len) | banded pair-HMM forward over `2W+1` candidate lengths |
| **output** | one confident on-ladder length `L*` (+ weight) | a sparse `Qᵣ(L)` distribution; optionally an off-ladder sequence |
| **stored as** | tally into the locus **histogram** (`hist_*`) — *not per-read* | per-read **CSR** rows (`amb_*` / `offl_amb_*`) |
| **share** | the bulk (~majority of reads) | the minority — but the **cost driver** (§4.2 "slow-path fraction") |

The gate between them is one predicate (`FAST_PATH_GATE`, spec §4.2): *both flanks
cleanly aligned ∧ pure tiling ∧ no interior sequencing indel ∧ boundary Q ok.* A
read that fails *any* clause falls to the slow path. A **soft-clipped read never
qualifies** (a clipped flank is not "cleanly aligned"), so it always takes the
slow path — which is exactly where its true (possibly large) length is *recovered*
from the clip (§3).

> **Recommendation:** build the gate as a single function returning
> `enum Tier { Fast(FastHit), Slow(SlowReason) }`, where `SlowReason` records
> *why* (soft-clip / interior-indel / impure / low-Q) — that reason is what we
> report in the per-stage diagnostics the spec demands ("the fast/slow read
> split, and the slow-path share", §4.2). Don't collapse it to a bool; the split
> is a headline benchmark number.

---

## 3. Read triage (`pileup/triage.rs`)

Given the reads already near a locus, **triage** them: anchor each to the flanks,
recover soft-clipped tract ends, and sort each into spanning / allele-too-long /
flanking / FRR / filtered. The subtle stage — everything downstream trusts this
sorting. (Note this is a *different* decision from the fast/slow gate of §2: that
gate decides *how to compute* a usable read's likelihood; triage decides *which
reads are usable evidence* at all.)

**3.1 Where the reads come from (not this module's job).** Reads overlapping
`[start − flank_bp, end + flank_bp]` arrive in a bundle from the **fetcher** (§8),
not from `triage.rs`. The existing
[`AlignmentMergedReader::query`](../../src/bam/alignment_input.rs) already yields
every read whose footprint overlaps a region, sharing one FASTA `Repository` for
CRAM decode — the sanctioned cross-caller I/O reuse
([reuse map](ssr_genotyping_architecture.md), §5). The fetcher also applies the
per-locus depth cap (§8.3) before handing the bundle off, so `triage.rs` just
takes "the (capped) reads near this locus" and classifies them.

**3.2 Anchor & soft-clip recovery — the hard part.** When a read carries an
allele **longer than the reference**, the aligner often cannot place the extra
units and **soft-clips** them, pushing a flank *into the clip*. Naive CIGAR
parsing of the aligned span then silently undercounts long alleles — the failure
mode that makes SSRs hard. So we **use the full read sequence including the clip**
and recover the length from its content, not from the CIGAR.

**The recipe — a cheap content pre-probe, then the targeted pair-HMM (decided,
matching both reference tools).** The two vendored STR callers agree on the shape,
and we follow it: use the full read (clip included), do a cheap `O(read length)`
scan of the read's own repeat content *first* to estimate the length and bound the
work, then realign only over that bounded set. **Neither tool blindly wide-band-
aligns**, so neither do we.

- **Pre-probe = count the motif in the read, GangSTR-style.** Scan the read
  (clip included) for the longest contiguous motif run + total motif copies —
  GangSTR's `find_longest_stretch`
  ([realignment.cpp:28-64](../../GangSTR/src/realignment.cpp)), an `O(read len)`
  pass that bounds its realignment range and filters junk. This is cheaper and
  more robust than aligning the clip to the flank (it needs no flank match to
  estimate the count). It gives the `count` that centres the pair-HMM window
  (`count ± W`, §5).
- **Plus a clean-flank check** against the catalog's embedded flank
  (`Locus::left_flank` / `right_flank`) — this is the **junk test and the
  spanning test** (§3.3): a real long allele shows *extra clean motif units + a
  clean flank surviving in the clip*; an adapter / low-quality end shows *no clean
  flank* and is rejected (and is the "allele ≥ read length" case if the tract runs
  off the end). HipSTR likewise requires a soft-clip to carry the flank before it
  trusts a clip-spanning read
  ([bam_processor.cpp:172-186](../../HipSTR/src/bam_processor.cpp)).
- **Then the pair-HMM (§5)** runs over the narrow `count ± W` window — the
  pre-probe is what lets the band stay narrow (`W`-sized), the whole point of
  banding. Recovery is a **slow-path** operation, not a CIGAR count.

> **Why not just re-align everything (HipSTR's choice)?** HipSTR re-aligns *every*
> read from scratch (Needleman-Wunsch in a 75 bp window,
> [AlignmentOps.cpp](../../HipSTR/src/SeqAlignment/AlignmentOps.cpp)) — it doesn't
> trust the CIGAR even for clean reads. We deliberately **don't**: our two-tier
> design trusts the alignment for the clean fast-path majority (§4) and pays the
> pre-probe + pair-HMM only on the soft-clipped / ambiguous minority. That is our
> efficiency edge, and it rests on the fast-path gate (§2) being conservative
> enough that a misclassified long allele is rare. (GangSTR similarly keeps a fast
> CIGAR-based estimate and only falls back to realignment when the CIGAR indels
> aren't a clean motif multiple.)

**3.3 Spanning is a property of recovered *content*, not the CIGAR (§4.1).** A
read is **spanning** iff its sequence carries ≥ `MIN_FLANK_BP` clean flank on
**both** sides of the full tract — whether those bases were *aligned* or
*recovered from a clip*. Classification:

- **spanning** (aligned flanks *or* clip-recovered clean flanks) → used for the
  length likelihood. A clip with a clean flank carries a **real long allele**.
- **allele ≥ read length** (clip yields *no* clean flank — tract runs off the
  read end) → **counted, not used** (the §1.4 scope boundary). Surfaces as
  `n_frr` / a flanking count.
- **flanking / in-repeat** (non-spanning, not FRR) → counted, not used in v1.

The counts feed the QC scalar columns (`depth`, `n_spanning`, `n_flanking`,
`n_frr`, `n_filtered`, `n_flank_indel`, `mapped_reads`); only the **spanning**
class flows into the likelihood.

> **MAPQ owns mappability (catalog §4 decision).** Low-MAPQ reads are filtered
> here (→ `n_filtered`); unmappable/paralogous loci self-suppress (low depth →
> Stage-2 no-call). No reference-side mappability track. The MAPQ threshold is a
> Stage-1 parameter (§10).

---

## 4. Fast path (`pileup/fast_path.rs`)

The cheap bulk. A read that passes `FAST_PATH_GATE`: walk the tract between the
two clean flank anchors, confirm it is a pure integer tiling of `Locus::motif`,
emit one on-ladder length `L*` (repeat units) + a base-quality weight. O(read
length), no DP. A fast-path read is **on-ladder by construction** (a pure tiling
is a clean rung).

These reads are **not written individually** — they tally into the locus
histogram (`hist_lengths` / `hist_counts` / `hist_weight`), per-read identity
discarded. This is the storage win: ~majority of reads cost one histogram bump,
not a stored profile.

> **Recommendation:** keep `fast_path.rs` allocation-free per read — it runs on
> the hot majority. The tract walk is a strided motif compare; the only output is
> `(units: u16, weight: f32)` pushed to a per-locus `HashMap<u16, (count,
> weight)>` (or a small sorted vec, since lengths cluster). Scratch-buffer
> discipline ([scratch buffers](../../src/baq/scratch.rs) ethos) applies.

---

## 5. Slow path — the banded pair-HMM (`pileup/pair_hmm.rs`)

**The net-new numerical machinery, and the stage's main risk** (architecture §5).
For every read failing the fast gate, compute a **forward (sum-over-alignments)**
score against each candidate haplotype
`H_L = left_flank + (motif × L) + right_flank` for `L ∈ [count − W, count + W]`
(`W = STUTTER_WINDOW_UNITS`), where `count` is the **recovered** repeat count
(§3.2) — so the window centres on the true length even for a soft-clipped long
allele, not on the reference.

The design below is **grounded in the two vendored STR callers** (read directly):
**HipSTR** does exactly this — a forward, base-quality-aware pair-HMM realigning
reads to STR haplotypes ([HipSTR/src/SeqAlignment/](../../HipSTR/src/SeqAlignment/),
`HapAligner`, `AlignmentModel`), itself a port of **Dindel** (Albers et al. 2011);
**GangSTR** takes the simpler route ([realignment.cpp](../../GangSTR/src/realignment.cpp))
— Smith-Waterman *Viterbi*, base-quality ignored, one realign per copy number.
Spec §4.2 chose forward + base-quality, so **HipSTR is our model**. The headline
finding: **our HMM is HipSTR's minus its hardest part** (next).

**5.1 What we keep from HipSTR, and the big thing we drop.** HipSTR's aligner is
complicated mostly because it **bakes stutter into the alignment** — its
`StutterAlignerClass` marginalizes over PCR-artifact (stutter) sizes *during* the
DP ([StutterAlignerClass.cpp](../../HipSTR/src/SeqAlignment/StutterAlignerClass.cpp)).
Our stage split puts stutter in Stage 2, so **we delete that entirely**: each
`H_L` is a plain fixed haplotype and the pair-HMM scores the read against it under
**sequencing error only**. We keep HipSTR's flank pair-HMM (states, emission,
transitions); we throw away the stutter machinery. This is the simplification the
§1 stage split buys us.

> **Trade-offs vs HipSTR's integrated DP — what the split does and doesn't cost.**
> Recorded here so it isn't re-litigated. **The factorization is exact, not a
> shortcut:** the generative chain `allele a → (stutter) → molecule length L →
> (seq error) → read` factors as `P(read | a, θ) = Σ_L S_θ(L|a)·Qᵣ(L)`, and `Qᵣ(L)`
> provably does **not** depend on the stutter parameters `θ` (θ governs how the
> molecule reached length `L`, not how a length-`L` molecule is read). So freezing
> the stutter-free `Qᵣ` in Stage 1 and convolving `θ` in Stage 2 loses nothing *at
> the `Qᵣ` level* — there is no "stutter should have informed the alignment"
> feedback being discarded, and on clean loci the split is exact. HipSTR computes
> the same sum, just inside its aligner. What the split **does** give up is four
> things — the first two approximation, the last two capability:
>
> 1. **Sparsification loss.** We store `Qᵣ` only on the `count ± W` window (§5.5)
>    and then prune to a sparse set (`AMB_LL_DROP`, §11). HipSTR *also* truncates
>    its stutter range (±6 units), so the window isn't the difference — the
>    **pruning is**: dropping an `L` with non-negligible `Qᵣ(L)·S_θ(L|a)` mass
>    loses it. Controllable (size `W`, loosen `AMB_LL_DROP`), but a real lossy step.
> 2. **Impure loci.** The factorization is exact only where "length `L`" is a
>    sufficient statistic — i.e. a length-`L` molecule has one sequence (clean
>    rungs, where `H_L` pins the interruption). At an **interrupted tract**, two
>    same-length molecules can differ by *where* a stutter unit landed; HipSTR's
>    `StutterAlignerClass` marginalizes that position inside the DP, whereas we
>    collapse to reconstructed rungs + per-read off-ladder candidates. Weaker
>    exactly at impure loci — already the stutter model's flagged weak regime (spec
>    §5.2).
> 3. **No cross-read/cross-sample assembly (the biggest give-up).** Stage 1 is
>    per-read, per-sample, no cohort context; the spec forbids assembly (§4.2). So
>    an allele present only in error-buried form across many reads — or **two
>    same-length alleles differing by a SNP inside the repeat** — which HipSTR's
>    haplotype assembly can surface, we capture only if a *single* read supports it
>    cleanly as off-ladder. A precision-first **sensitivity/recall** cost.
> 4. **No physical phasing / joint read modelling.** The in-memory design can phase
>    the STR against nearby SNPs and reason jointly over reads; we discard per-read
>    identity for the confident majority and excluded SNP-phasing by design (spec
>    §1.1).
>
> **Net:** for our scope — length genotyping, period ≤ 6, pop-gen, precision-first,
> **large cohorts** — these are small-and-tunable (1), in an already-weak regime
> (2), or deliberate precision/scope trades (3, 4). In exchange the split buys what
> HipSTR's in-memory re-alignment-inside-EM cannot: θ-independence (run the pair-HMM
> **once**, not per EM iteration), per-sample parallelism, a memory-bounded
> columnar intermediate, and cohort scaling. The split trades some sensitivity and
> some impure-locus accuracy for cohort scale and speed; it is **not** a shortcut
> on the common (clean-locus) case.

**5.2 Model — 3-state forward, log-space.** Match / Insertion / Deletion;
**forward (log-sum-exp), not Viterbi** — a true likelihood that reports honest
ambiguity (a cleanly-long soft-clipped read just yields a sharply-peaked `Qᵣ`).
HipSTR confirms forward is the right call here; GangSTR's Viterbi is the corner we
do *not* cut. Log-space with a rescale/`log_sum_exp` guard against underflow (BAQ's
`RESCALE_THRESHOLD` pattern, [probaln.rs](../../src/baq/probaln.rs)).

**5.3 Emission — Dindel base-quality model (already in our codebase).**
`match = log(1 − 10^(−Q/10))`, `mismatch = log(10^(−Q/10) / 3)`, from a **256-entry
per-Q lookup**. This is *exactly* what our BAQ engine already computes
([scratch.rs](../../src/baq/scratch.rs) `Q2P`; [probaln.rs](../../src/baq/probaln.rs)
`EM = 0.33333…` = 1/3) and what HipSTR uses (`BaseQuality.h`). We reuse the
*pattern* (a process-wide `LazyLock` lookup), not the BAQ code, per the
no-coupling decision (§5.6).

**5.4 Transitions — Dindel homopolymer-indexed affine gaps (adopt verbatim).**
HipSTR ports Dindel's parameters directly, and we take the same values:

- **M→indel (gap open) is indexed by homopolymer-run length**, not constant —
  Dindel's `{2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4,
  1.0e-3, 1.4e-3}` for run lengths 1–10, linear extrapolation above, capped at 15
  ([AlignmentModel.cpp](../../HipSTR/src/SeqAlignment/AlignmentModel.cpp)).
- **gap extend** fixed at `e^−1 ≈ 0.368` (`INS→INS` / `DEL→DEL`); `gap→M` =
  `1 − e^−1`.

These model **sequencing** indels — *including* the elevated indel rate inside
homopolymer runs, which is exactly the sequencing-error inflation a mononucleotide
SSR suffers (and is harmless for period ≥ 2, where the in-motif homopolymer
context is short). This is **not** allele stutter — that is Stage 2's `S_θ`. The
two slippage phenomena stay cleanly separated: homopolymer *sequencing* slippage
here, allele *stutter* there.

**5.5 The `2W+1` layout & band width — start independent + narrow, optimize to
shared-flank.** The two vendored tools bracket the choice: GangSTR realigns
**independently** per length; HipSTR **shares** the flank DP across all candidates
(align each flank once, reuse, recompute only the repeat block). The clarifying
point is that **band width follows from the layout**:

- **Independent forwards (recommended first).** Align the read against each
  *right-length* `H_L` on its own. The read then sits **near-diagonal**, so the
  band need only cover **sequencing-indel slack** — a handful of bp (BAQ's default
  bandwidth ≈ 7 is the precedent). The `W·period` candidate spread is absorbed by
  *enumerating* `H_L`, not by a wide band. Simple, obviously correct, and
  `2W+1 ≈ 5–13` forwards per slow read is cheap.
- **Shared-flank / nested lattice (the optimization).** Because our `H_L` differ
  *only* in the middle repeat count and share **byte-identical flanks** (even
  cleaner than HipSTR, which also juggles interior SNVs), the left-flank forward
  can be computed **once** and the candidates read off as "exit the repeat into the
  right flank at copy `L`." This is HipSTR's seed/shared-flank trick specialized to
  pure length variation; band then widens to `≈ W·period`. Adopt **only if
  profiling shows the slow path binds** (spec §4.2 says measure).

**5.6 Reuse the pattern, not the code (§4.2 / [reuse map](ssr_genotyping_architecture.md), §5).**
The BAQ engine is a *port of htslib `probaln_glocal`* tuned for byte-parity with
htslib — coupling to it would drag that parity contract into the SSR path for no
benefit. We borrow the **banded-forward loop shape**
([probaln.rs:48](../../src/baq/probaln.rs#L48)) and the **scratch discipline**
(below), but the math is bespoke and **simpler** than `probaln_glocal` (no backward
pass, no posterior decoding, no htslib parity literals). A clean ~200-line 3-state
forward.

**5.7 Scratch — a rolling two-row buffer (we can beat HipSTR here).** HipSTR keeps
**full** `read × hap` matrices because it does **traceback** (for its alignment
visualizer) and `new`/`delete`s them per read (it flags this as a perf gap). **We
need only the final forward sum `Qᵣ(L)`, never an alignment**, so we keep just the
**previous DP row** to compute the next — `O(band)` memory per state, not
`O(read × hap)`. One `PairHmmScratch` per worker, `resize_for(read_len, band)`
grow-and-keep (BAQ's [scratch.rs](../../src/baq/scratch.rs) pattern), **zero
per-read allocation** on the hot path.

**5.8 Off-ladder emission (§4.2).** If a read's best alignment is a definite
*non-rung* sequence (an in-frame count + a 1 bp indel, a partial unit, a variable
interruption), the forward *also* emits that **normalized** sequence as its own
candidate with its forward score — so Stage 2 treats it as a distinct allele.
Normalization is the shared indel-norm kernel (§7), not a re-implementation. The
off-ladder leg is the §3.2-of-the-types-doc canonical key that makes the
cross-sample union work.

**5.9 Flanks are clean by construction (§4.2).** Stage 0 dropped bundled/compound
loci, so every catalog locus is isolated (no detected repeat within `flank_bp`).
The realigner anchors on `Locus::left_flank` / `right_flank` with **no inner-flank
and no compound special case** — a real simplification we inherit from the catalog
decisions.

---

## 6. The ladder — candidate construction & off-ladder normalization (`pileup/ladder.rs`)

`ladder.rs` is the **translation layer between reads and alleles.** Up to this
point the stage works in raw read bytes and a repeat *count* ("this read looks
like 12 repeats"). Everything after it — the stored evidence, the cohort merge,
the VCF — works in a typed **allele**, of which there are exactly two kinds (the
two variants of the `Allele` type, [ssr_shared_types.md](ssr_shared_types.md) §2,
explained in full here so this section stands alone):

- **on-ladder** — a clean rung of the repeat ladder, identified by a **single
  number**, its repeat count (e.g. "the 12-unit allele"). The number is enough
  because the full sequence is rebuildable from the reference: it is the left
  flank, the motif tiled that many times, then the right flank.
- **off-ladder** — a sequence that is **not** a clean rung (12 repeats plus one
  stray base, a partial unit, an internal interruption). It has no rung number, so
  we must carry its literal sequence.

The module has two jobs, one per kind.

**Job 1 — build the candidate rungs as sequences (what the pair-HMM aligns
against).** The repeat ladder is just the set of plausible lengths — rung 10, rung
11, rung 12, … The pair-HMM (§5) cannot align a read against a *number*; it needs
an actual stretch of DNA. So `ladder.rs` writes each rung out in full:

```
rung L  =  left_flank  +  (motif × L)  +  right_flank
```

— e.g. rung 12 of a `CAG` repeat is `…flank… CAGCAGCAG…(12×)… flank…`. It builds
these from the catalog locus, which already carries the reference bases and the
motif (so no FASTA is opened), and only for the handful of rungs in the read's
window (`count ± W`, §5), not the whole ladder. The property that matters
downstream: because every rung is built from the **catalog** — one shared
reference for the whole cohort — **rung 12 is the byte-identical sequence in every
sample.** That is exactly what lets Stage 2 recognise "the 12-unit allele" as the
same allele across all samples by its number alone.

**Job 2 — give every off-ladder oddball one canonical spelling.** An off-ladder
allele is carried as its literal sequence, and that creates a subtle trap: the
*same* allele can be written more than one way. It is the **indel-representation
ambiguity familiar from VCFs** — an inserted base can be slid left or right and
still describe the same variant. If sample A spells the allele one way and sample
B another, the cohort step sees two different sequences and splits one real allele
into two phantom ones, scattering its evidence across both.

To prevent that, `ladder.rs` rewrites every off-ladder sequence into **one
agreed-upon canonical form** (left-aligned, minimally trimmed) before it is
stored. Then the rule is exact: **the same allele, in any two samples, always
produces identical bytes** — so the cohort union counts it once. The canonicaliser
is not new code: it is the very routine the SNP caller already uses to left-align
indels (a port of GATK / `bcftools norm`), reused here rather than re-implemented
(see §7) so the two callers can never drift apart on what "the same allele" means.

So the module is named for the *ladder* because it does both halves of the
ladder's bookkeeping: it constructs the rungs, and it decides — for everything
that lands *between* the rungs — the single form in which it will be remembered.

> **Dependency this stage forces on `types.rs`.** Today
> [types.rs](../../src/ssr/types.rs) carries only `Motif` and `Locus` (the
> catalog's output). Stage 1 is the first consumer of the **allele
> representation** the shared-types doc designed but did not yet build —
> `Allele` (`OnLadder`/`OffLadder`), `NormalizedSeq`, and the
> `to_sequence`/`repeat_count` methods. **Building these is the first task of the
> Stage-1 pass**, before `ladder.rs`. They are designed (types doc §2–§5); this
> is the build trigger.

---

## 7. Shared-code lift this stage forces — `normalize_alleles`

The off-ladder path (§5.3/§6) **reuses the SNP caller's actual indel-norm
kernel**, not a copy (types doc §4, decided). The mechanics land here:

- Lift `normalize_alleles` (+ its `Range`/helpers) out of
  [pileup/walker/indel_norm.rs](../../src/pileup/walker/indel_norm.rs) into a
  **shared, public, representation-neutral module** — it is already a
  CIGAR-agnostic `(seqs, bounds)` primitive with two real users now.
- The SNP CIGAR path keeps its wrapper and calls the kernel; **regression gate =
  SNP end-to-end tests pass** (it is reference-validated, GATK/freebayes
  cross-checked — must not diverge).
- The SSR off-ladder path is a **thin adapter** building `(seqs, bounds)` from
  `(off-ladder candidate, ref tract)` — *not* a CIGAR-faking shim.

> Same incremental rule as the container refactor (architecture §10.6): the lift
> is a behaviour-preserving move, landed and green *before* the SSR adapter is
> wired on top.

---

## 8. Data-flow & parallelism

Architecture §7/§9 deferred Stage-1 parallelism to *this* pass, explicitly **not**
inheriting the SNP thread machinery (built for a per-position streaming pipeline;
the SSR shape is different). The stutter-free / locus-independent property (§1)
makes loci embarrassingly parallel, so the only real questions are *how reads
reach loci* and *how the file access is kept sane*. The design below separates
**I/O from compute**: a single thread owns all BAM access, and a pool does the
CPU work from a queue of self-contained bundles.

### 8.1 Why indexed access, not a whole-BAM scan

Catalog loci are **sparse** — Stage 0's whole point is a small curated locus set,
covering a tiny fraction of the genome. Two facts follow, and they settle the
read-fetching strategy without a benchmark:

- **Reads are coordinate-ordered**, so for one locus you seek to the first read
  overlapping `[start − flank_bp, end + flank_bp]` and read forward until the
  first read past it — noodles' region query
  ([`AlignmentMergedReader::query`](../../src/bam/alignment_input.rs)) does exactly
  this, sharing one FASTA `Repository` for CRAM decode.
- **Loci are far apart**, so jumping to each via the index and decoding only its
  reads is far cheaper than streaming the whole BAM and discarding the megabases
  between loci. (A full scan would only win for a pathologically *dense* catalog,
  which ours is not.)

So: **indexed access, driven off the sorted catalog.** The one cost to manage is
seeking — addressed next by centralizing it.

### 8.2 The pipeline — one fetcher, a bounded queue, a worker pool, an ordered collector

```
 loci (sorted) ─►┌─ FETCHER (1 thread) ──────────────┐         ┌─ POOL (N workers) ──────────┐
                 │ forward-only walk of the catalog;  │ bounded │ per bundle (no BAM access):  │
 BAM/CRAM ──────►│ index-jump to each locus, read its │  queue  │  triage → fast/slow pair-HMM │ ─►┌ collector ┐
                 │ reads, build a (locus + reads)     │ ──────► │  → ladder → aggregate to     │   │ reorder   │─► sampleN
                 │ bundle, with a depth cap (§8.3)    │ (back-  │  SsrLocusRecord              │   │ by (chrom,│   .ssr.psp
                 └────────────────────────────────────┘ press.) └──────────────────────────────┘   │ start)    ┘
                                                                                                    └───────────┘
```

- **Single fetcher thread — this is what kills seek-thrash.** The back-and-forth
  that worries a per-locus design only happens if *many* threads seek
  independently. With **one** thread walking the catalog in coordinate order,
  every seek is **forward-only** (monotonic), and the workers never touch the file
  at all. The fetcher decodes only the covered fraction, builds a self-contained
  `(locus + reads)` bundle, and hands it off.
- **Compute is per-locus, not per-read.** A locus is the natural parallel unit
  because its outputs are **reductions** over its reads — the length histogram, the
  candidate set, the off-ladder union are all "combine all of this locus's reads."
  Keeping a locus whole on one worker makes that aggregation purely local (no
  locking); splitting a locus's reads across workers would force a locked merge of
  partial per-locus state. A locus is also a meaty enough unit (several reads, some
  slow-path pair-HMM) to amortize task overhead, where a single read is too
  fine-grained.
- **Bounded queue = back-pressure + memory bound.** If the fetcher outruns the
  pool, an *unbounded* queue would balloon RAM. A bounded queue blocks the fetcher
  when full — the same mechanism as the SNP pipeline's bounded channel. Bound it by
  **reads/bytes in flight, not bundle count**, since a deep locus is a large
  bundle. (The depth cap, §8.3, bounds each individual bundle; this bounds the
  number in flight.)
- **Ordered collector (determinism + the block grid).** Workers finish out of
  order, but the `.ssr.psp` block grid + tail index need records in
  `(chrom, start)` order — so a reorder buffer restores that order before the
  `registry_ssr` writer ([architecture](ssr_genotyping_architecture.md), §10).
  Records are small and sparse, so buffering is cheap. This is also what makes the
  output byte-identical across `--threads` (§8.4).
- **A read overlapping two nearby loci** is simply copied into both bundles — rare
  (loci far apart, reads short), so no special handling; just don't assume bundles
  partition the reads.
- **One free fetch optimization:** where several loci sit close together, the
  fetcher can cover them with **one query over the cluster's span** rather than one
  query per locus (fewer seeks, sequential decode). Where loci are far apart — the
  common case — a cluster is one locus and this changes nothing.

**Scope: one invocation = one sample (not the cohort).** `ssr-pileup` processes a
**single sample** and writes that sample's one `.ssr.psp`. It accepts **one or
more** BAM/CRAM inputs, but they must all be the **same sample** (lanes /
replicates / re-sequencing of that individual); the fetcher k-way coordinate-merges
them into one ordered read stream — the same merge the SNP `pileup` already does
([`AlignmentMergedReader`](../../src/bam/alignment_input.rs)). It writes the sample
name into the `.ssr.psp` header and rejects inputs whose read groups disagree.

Running **one process per sample across a cohort is the user's orchestration**, not
this subcommand's concern — exactly as for the SNP `pileup`
([benchmark each caller with its native parallelism]). So the cross-sample axis
does not appear in this design at all; `--threads` here is purely the
*within-sample* pool of the pipeline above. (The output RSS lever is
`--block-window-bp`, as on the SNP path, [block-window memory lever];
`locus_record.rs` does the per-locus aggregation → `SsrLocusRecord` — confident →
`hist_*`, ambiguous → `amb_*` CSR, off-ladder → `offl_*` — written via the generic
container + SSR schema, architecture §10.4.)

> **If the fetcher ever becomes the floor** (doubtful — it decodes only the
> covered fraction, while the pair-HMM pool does the heavy work), the escape hatch
> is to **shard the fetcher by contig**: several fetchers on *disjoint* regions,
> which still never seek over each other. Defer until measured.

### 8.3 The per-locus depth cap — reservoir sampling

A hypervariable, high-depth locus could make one bundle enormous (and one worker's
task fat). Cap each locus at `MAX_READS_PER_LOCUS`, mirroring the SNP path's
depth cap — but do it **in the fetcher's single pass**, since the total depth is
not known until the locus is fully read.

**Reservoir sampling (Algorithm R)** does exactly this. Keep the first `K =
MAX_READS_PER_LOCUS` reads; for the `i`-th read with `i > K`, keep it with
probability `K / i`, and if kept, drop one of the `K` currently held uniformly at
random. At the end every read that covered the locus has the **same** probability
`K / n` of being in the bundle (`n` = total depth) — an unbiased uniform sample,
computed in one streaming pass with **`O(K)` memory and no second read of the
BAM**, exactly as you'd want for the fetcher.

> **Determinism caveat (it interacts with §8.4).** Random subsampling must **not**
> make the output depend on timing. Seed the per-locus reservoir RNG
> **deterministically from the locus** (e.g. from `(chrom, start)`), not from
> wall-clock or thread id, so the *same* reads survive the cap on every run and at
> every `--threads`. With a per-locus seed the cap stays byte-identical. (The cap
> threshold and the seed scheme go in the `.ssr.psp` header's `extraction_params`,
> §10, so the subsample is reproducible and self-describing.)

### 8.4 Determinism (invariant)

Stage 1 is stutter-free with no cohort context, so it holds the **byte-identical
output regardless of `--threads`** bar — the same standard as the catalog
([ssr_catalog.md](ssr_catalog.md) §8.4) and the SNP pileup. Three things give it:
the **fetcher** reads loci in a fixed coordinate order; the **collector** emits in
`(chrom, start)` order regardless of which worker finished first; the **pair-HMM**
is deterministic (fixed iteration order, no parallel float reduction within a
read). The one subtlety is the depth cap — handled by the per-locus deterministic
seed (§8.3). A regression gate.

---

## 9. `flank_bp` — pinning the catalog parameter (owned here)

Two different "flank" sizes, easy to conflate:

- **`MIN_FLANK_BP`** — how much *clean flank a read must show* on each side to count
  as spanning (the per-read anchor requirement, §3.3).
- **`FLANK_BP`** — how much *reference flank the catalog embeds* around each tract
  (`ref_seq` = tract + `FLANK_BP` each side, [ssr_catalog.md](ssr_catalog.md) §5) —
  the reference context the pair-HMM aligns the read's flank against. The catalog
  stores it but punted the *value* here, because Stage 1 is the binding consumer.

`FLANK_BP` must be the larger: `FLANK_BP ≥ MIN_FLANK_BP + (pair-HMM band, §5)`, so
the realigner has reference context to anchor the flank against and the band can't
run off `ref_bytes`.

**Values — decided, converging on the two reference tools.** Both vendored callers
independently land on the same anchor and a ~30 bp embedded flank:

- **`MIN_FLANK_BP = 5`.** HipSTR's `MIN_FLANK = 5` (a read must extend ≥5 bp past
  the region each side, [bam_processor.cpp:287](../../HipSTR/src/bam_processor.cpp))
  and GangSTR's `min_match = 5` ("minimum matching basepairs on each end of an
  enclosing read", [options.cpp:88](../../GangSTR/src/options.cpp)) agree. (HipSTR
  *also* checks flank *uniqueness* over a 15 bp window; we don't need that — MAPQ
  owns mappability, §3.3.)
- **`FLANK_BP = 30`.** HipSTR carries ~30–35 bp of reference flank to build/align
  haplotypes (`MAX_REF_FLANK_LEN = 30` / `REF_FLANK_LEN = 35`,
  [seq_stutter_genotyper.h:153](../../HipSTR/src/seq_stutter_genotyper.h)); 30
  comfortably clears `5 + band(~7)` and is tiny in the catalog. (GangSTR's
  `realignment_flanklen = 100` is larger because it targets expansions/long reads —
  overkill for our spanning-only, period ≤ 6 scope.)

**The tension that makes `FLANK_BP` *not* free — the reason it was left open.**
`FLANK_BP` is coupled to the catalog's `bundle_threshold`: Stage 0 **drops** any
locus with another repeat within `bundle_threshold` bp (to guarantee clean unique
flanks), and we **set `bundle_threshold = FLANK_BP`** so the clean-flank guarantee
and the embedded extent coincide. The consequence:

> **A bigger `FLANK_BP` discards more catalog loci.** Every SSR with a neighbour
> within `FLANK_BP` is bundled out. So `FLANK_BP` trades *flank cleanliness /
> alignment quality* against *catalog completeness (recall)* — and in a
> repeat-dense plant genome, closely-spaced SSRs are common, so pushing `FLANK_BP`
> 30 → 100 would silently drop a real chunk of markers. **30 bp is the sweet
> spot:** enough to anchor + band, small enough not to bundle away too many
> neighbours.

**To measure, not design:** in the Stage-0 accuracy harness, the surviving-locus
count as a function of `FLANK_BP` (the completeness cost), and on the read side the
anchor pass-rate at 5 vs 8 vs 10 bp on real data — to *confirm* 30/5 rather than
assert them. Both stay named `const`s with a CLI override, recorded in the catalog
`##` header and the `.ssr.psp` `extraction_params`.

---

## 10. Parameters / constants (spec §4.2/§4.3, named here)

All become named `const`s (house style — no magic numbers), recorded in the
`.ssr.psp` TOML header's `extraction_params` so the evidence is self-describing
and Stage 2 can verify the cohort was extracted consistently:

| const | role | §ref |
|---|---|---|
| `MIN_FLANK_BP` | clean flank required each side to call a read spanning | §3.3 |
| `MIN_BASE_QUAL` | boundary base-quality floor for the fast gate | §2 |
| `FAST_PATH_GATE` | the bundled fast-path predicate | §2 |
| `STUTTER_WINDOW_UNITS` (`W`) | half-width of the candidate-`L` window | §5 |
| `AMB_LL_DROP` | sparse-profile cutoff (drop candidates > this below per-read max, then renormalize) | §11 |
| `FLANK_BP` | embedded-reference margin (≥ `MIN_FLANK_BP` + band; = `bundle_threshold`) | §9 |
| `MIN_MAPQ` | read mappability filter | §3.3 |
| `MIN_MOTIF_RUN_BY_PERIOD[]` | pre-probe trigger: min motif-run in the read to attempt slow-path recovery (period-indexed, GangSTR-style: ~5/4/3 for period 2/3/≥4) | §3.2 |
| `MAX_READS_PER_LOCUS` (`K`) | per-locus depth cap; excess reservoir-sampled (+ the per-locus seed scheme) | §8.3 |
| `PAIR_HMM_BAND_BP` | banded-forward half-width (sequencing-indel slack, ~7) | §5.5 |
| `DINDEL_GAP_OPEN[]` / `GAP_EXTEND_PROB` | Dindel homopolymer-indexed gap-open table + fixed extend (`e^−1`), adopted verbatim | §5.4 |

> **Naming note (spec §4.2):** `STUTTER_WINDOW_UNITS` sizes the candidate-`L`
> range *even though Stage 1 is stutter-free* — the window must be wide enough to
> cover plausible stutter excursions so **Stage 2** has support points to convolve
> the kernel against. It bounds *which lengths get a likelihood*, not a stutter
> computation here.

---

## 11. How `Qᵣ` is stored — the two tiers map onto two regimes (§4.3)

The output columns are settled in spec §4.3 and the container schema in
[architecture](ssr_genotyping_architecture.md) §10.4. Stage 1's job is to *fill*
them; the mapping the aggregator implements:

- **Fast-path reads → histogram.** Confident `L*` collapses into `hist_lengths`
  (ascending distinct `uint16` lengths) / `hist_counts` / `hist_weight`. Per-read
  identity discarded — the bulk, stored as a tally.
- **Slow-path reads → sparse CSR.** The forward is *evaluated* over `2W+1`
  lengths but **stored sparse**: drop candidates > `AMB_LL_DROP` below the
  per-read max, renormalize over the survivors (typically 1–3 lengths with real
  mass), write into `amb_read_offsets` (CSR prefix) / `amb_lengths` /
  `amb_logliks`. The only reads stored individually.
- **Off-ladder → `offl_*`, mirroring the on-ladder regimes.** Confident
  off-ladder reads tally into `offl_seqs` (dict) / `offl_counts` / `offl_weight`;
  a read ambiguous *between* an off-ladder sequence and ladder rungs carries its
  off-ladder leg in `offl_amb_*`. **Empty at the vast majority of loci** — cost
  paid only where real off-ladder signal exists.

All stored log-liks are **stutter-free** (the invariant). `units` lengths are
`uint16` (non-negative; the only signed length quantity, the ref-offset Δ, is
derived in Stage 2, never stored — types doc §2).

---

## 12. Testing (spec §7/§11 — Bucket-1 on the critical path)

Two levels, both via the crate-internal simulator
([architecture](ssr_genotyping_architecture.md) §2.1, not a subcommand):

- **Read/BAM-level → tests the whole stage incl. the pair-HMM.** Synthesize reads
  from known genotypes with injected sequencing error (and deliberately
  *soft-clipped long alleles*, the §3.2 failure mode), run `ssr-pileup`, assert
  the recovered `hist_*`/`amb_*`/`offl_*` match the injected lengths. This is the
  only way to exercise soft-clip recovery and the forward together.
- **Unit-level → the pair-HMM in isolation.** A known read × known haplotype with
  a hand-computable forward score; the fast/slow gate on boundary cases; ladder
  reconstruction round-trips (`to_sequence` ∘ count = identity); off-ladder
  normalization gives identical keys for the same allele built two ways.
- **Determinism gate (§8.4):** byte-identical `.ssr.psp` across `--threads`,
  *including* the depth-cap subsample (the per-locus seed, §8.3).

> The **evidence-level** simulator (synthetic `.ssr.psp` directly, no reads) is a
> *Stage-2* tool — it lets `ssr-call` be tested before `ssr-pileup` exists
> (architecture's critical-path order). Stage 1's tests are the *read*-level ones,
> because Stage 1 is precisely what turns reads into evidence.

---

## 13. Module layout (proposed — architecture §4 sketch, refined)

```
src/ssr/pileup/
├── mod.rs              # ssr-pileup driver: wire fetcher → bounded queue → pool → collector → write (§8)
├── fetch.rs            # single I/O thread: forward-only index walk, depth-cap reservoir, build bundles (§8.2/§8.3)
├── triage.rs           # anchor reads, soft-clip recovery, spanning classification (§3)
├── fast_path.rs        # flank-anchored exact motif count (§4)
├── pair_hmm.rs         # NEW bespoke banded 3-state forward + PairHmmScratch (§5)
├── ladder.rs           # on-/off-ladder candidate construction + normalization (§6)
└── locus_record.rs     # per-locus aggregation → SsrLocusRecord → registry_ssr (§8.2/§11)
```
(plus the `ssr-pileup` subcommand under
[pop_var_caller/cli/](../../src/pop_var_caller/cli/), reusing the shared
`--reference`/`--threads`/`--regions`/`--block-window-bp` parsers; `--reference`
is **CRAM-decode-only** here, spec §3.2.)

Net-new and shared-lift dependencies this stage pulls in, in build order:
1. **Extend `types.rs`** with `Allele`/`NormalizedSeq` + methods (§6 — first).
2. **Lift `normalize_alleles`** to a shared module (§7).
3. **Build the container SSR schema** (`registry_ssr` + `SsrLocusRecord`) if the
   container refactor (architecture §10) hasn't landed yet — it is a prerequisite
   for the writer and is independent of the SSR math.
4. Then the stage modules above, fast path → slow path → ladder → aggregation.

---

## 14. Decisions & residual validation

Every design question this pass raised is now **decided**; nothing structural is
left open. Recap, with the deferred items and the data-validation that remains.

**Inherited from upstream:** stutter-free Stage 1 (spec §4.2); two-tier fast/slow +
histogram/CSR storage (§2/§11); catalog is the only reference, FASTA for CRAM
decode only (§1); flanks clean by construction — no inner-flank case (§5.9); MAPQ
owns mappability (§3.3); share the SNP indel-norm kernel (§7).

**Settled in this pass:**

- **Parallelism (§8)** — indexed access off the sorted catalog; single fetcher
  (forward-only, no seek-thrash) → bounded queue → per-locus worker pool → ordered
  collector; per-locus depth cap by reservoir sampling with a deterministic
  per-locus seed (§8.3).
- **Pair-HMM (§5)** — HipSTR's flank model minus its stutter marginalization: a
  3-state forward in log-space; emission = the Dindel/BAQ per-Q model we already
  have; transitions = Dindel's homopolymer-indexed affine gaps (verbatim); layout =
  independent narrow-band forwards first, shared-flank lattice as a measured
  optimization; scratch = a rolling two-row buffer (no traceback ⇒ `O(band)`).
- **Soft-clip / long-allele recovery (§3.2)** — cheap content pre-probe (motif
  count in the read, GangSTR `find_longest_stretch`) + clean-flank check, then the
  targeted pair-HMM. We don't re-align every read (HipSTR's choice); the fast path
  trusts the alignment for the clean majority.
- **Flank sizes (§9)** — `MIN_FLANK_BP = 5` (HipSTR `MIN_FLANK` / GangSTR
  `min_match`) and `FLANK_BP = 30 = bundle_threshold` (HipSTR's embedded-flank
  size).
- **Where the pair-HMM lives** — `src/ssr/pileup/pair_hmm.rs`, **SSR-private**.
  Promote to a shared `src/align/` (next to `baq/`) only if a second user ever
  appears (architecture §9.3); kept private until then per the project's
  extract-on-second-use rule.

**Deferred by design (not v1):**

- **Read-pair merging (spec §1.4)** — deferred to a later version. It raises the
  spanning ceiling but is *not free* (combined-quality model + chimeric-merge
  scrutiny = the same slow path), so v1 is **spanning-only**; the merge is a clean
  later extension, not a v1 entanglement.
- **Stage-2 segdup/mappability revisit (spec §5.8)** — a Stage-2 concern, flagged
  there, not here.

**Residual validation — calibration, not design (owned by §12 / Stage-0 harness).**
These do not block building; they tune named `const`s already in place:

- the pair-HMM constants on simulated + real repeat-rich data — the homopolymer
  extrapolation past length 10, any STR-period adjustment to the `e^−1` gap-extend,
  the band-width value (§5);
- whether the single fetcher ever becomes the throughput floor (escape hatch:
  shard by contig, §8.2);
- the `FLANK_BP` locus-survival curve (completeness cost) and the read anchor
  pass-rate, to confirm `30 / 5` (§9);
- the shared-flank lattice optimization, built only if the slow path measurably
  binds (§5.5).
