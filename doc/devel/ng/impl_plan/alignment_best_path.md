# ng alignment — best-path aligners: implementation plan (1 of 3)

**Status:** draft, 2026-07-23. The build order for the **aligners that return one line-up**: the
`src/ng/alignment/` skeleton, the shared aligner types, the repeat-aware delimiter (the parity
anchor), its banding, the two-penalty comparison entry, and the general-purpose affine aligner.
Design is settled in [`../spec/alignment.md`](../spec/alignment.md) (§3, §4, §10) and
[`../arch/alignment.md`](../arch/alignment.md). This turns that design into build order; it is **not**
a place for new design.

The first of three: followed by [`alignment_marginal.md`](alignment_marginal.md) (plan 2) and
[`alignment_normalization.md`](alignment_normalization.md) (plan 3).

> **This plan goes first because it unblocks work that is stopped today.** Milestone B produces the
> `align_read` that [`locus_generation_ssr.md`](locus_generation_ssr.md) Milestone D has been waiting
> on — the only externally-blocking item in the whole module — and it is also the module's **only
> byte-parity oracle**, so building it first gives every later comparison a fixed reference point.
> Everything after B is comparison work; if the generator is the priority, B is the stopping point
> that matters.
>
> **One divergence from spec §10.3, recorded.** The spec interleaves algorithm 5 between 3 and 4, to
> land "reproduce production" as one unit. Grouping the plans by interface separates them. That is
> harmless here: algorithm 4 is compared against **algorithm 3**, not against 5, so Milestone D's
> comparison needs nothing from plan 2. The end-to-end production comparison spans this module and the
> genotyping in any case (spec §10.3).

---

## Scope

**In:** `src/ng/alignment/mod.rs` (the module skeleton every plan builds on); the `Alignment` type;
the `Emission` component and its two implementations; `RepeatSpan`, `RepeatGeometry`,
`RepeatContext`; `StutterModel`; the `BestPathAligner` trait; **algorithm 3** (the repeat-aware
delimiter, ported unbanded, then banded as a separate change); **algorithm 4** (whole-unit slips
priced apart); **algorithm 2** (the general-purpose affine aligner), gated.

**Out (later plans):**

- **`LogProb`, `MarginalAligner`, algorithms 5 and 6** → plan 2,
  [`alignment_marginal.md`](alignment_marginal.md).
- **`AlignmentNormalizer` and the three normalizers** → plan 3,
  [`alignment_normalization.md`](alignment_normalization.md). They consume the `Alignment` type this
  plan defines and nothing else from it.
- **The STR read preparer that calls algorithm 3**, and the tract-tally that consumes its output →
  [`locus_generation_ssr.md`](locus_generation_ssr.md) Milestone D and read preparation's plan.
- **What triggers re-alignment on the generic path** — unresolved, and it gates Milestone F
  (`../spec/read_preparation.md` §4).
- **The stutter parameters' estimation** — this plan consumes a `StutterModel`; fitting one is the
  genotyping's (spec §5.2).

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The parity anchor first.** Algorithm 3 is the only algorithm in this module with a byte-parity
  oracle; it is built first so every later comparison has a fixed reference point, and because the
  STR generator is blocked on it (spec §10.3).
- **Transcribe, then change — never both at once.** The port lands **unbanded**, matching production;
  banding is a separate, separately-proved commit (Milestone C). A behaviour change smuggled into a
  port destroys the parity that justifies the port.
- **Isolate the silently-wrong steps.** Three steps here produce a quietly-wrong *number* rather than
  a crash: the delimiter port (a wrong repeat length), the flank-anchor readout (a wrong observation
  class), and banding (silently-lost long alleles). Each is its own commit with its oracle green
  before and after.
- **Reuse over rewrite.** Algorithm 3 ports `delimit_read`'s recurrence, `HmmModel`'s constants and
  `ViterbiScratch`'s shape; the emission table is ported including its quality-zero floor (arch §5).
- **Build last what cannot yet fire.** Algorithm 2's calling mode has no trigger, so it is gated
  behind that decision rather than built speculatively (spec §10.3).
- **Incremental, with pauses.** One milestone, stop for review.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; native host build at the end.

## Preconditions (already in place)

- **`src/ng/types.rs` exists** (foundations A1). This plan adds nothing to it; `LogProb` lands in
  plan 2, where its first user is.
- **The design is settled**: spec §3 (the axes, banding), §4.1/§4.2 (both aligners, the costs, the
  determinism rules, the one-vs-two-penalty comparison), §5.2 (the stutter model); arch §2.1–§2.4,
  §3, §5.
- **The reuse targets and the parity oracle:** `delimit_read`
  ([alignment.rs:171](../../../../src/ssr/pileup/alignment.rs#L171)), `HmmModel` ([:80](../../../../src/ssr/pileup/alignment.rs#L80)),
  `ViterbiScratch` ([:111](../../../../src/ssr/pileup/alignment.rs#L111)), `Delimited`
  ([:141](../../../../src/ssr/pileup/alignment.rs#L141)), `EMISSION_LN` ([:59](../../../../src/ssr/pileup/alignment.rs#L59)),
  the gap constants ([:38](../../../../src/ssr/pileup/alignment.rs#L38), [:50](../../../../src/ssr/pileup/alignment.rs#L50));
  and for the stutter distribution, `stutter_pmf`
  ([hipstr.rs:168](../../../../src/ssr/cohort/read_model/hipstr.rs#L168)) with `MAX_SLIP`
  ([param_estimation.rs:21](../../../../src/ssr/cohort/param_estimation.rs#L21)).
- **A fixture** shallow enough to hand-check, carrying a clean repeat, a long allele, an interrupted
  repeat, and reads that run off each end.

---

## The steps

### Milestone A — module skeleton and the aligner types (no logic)

**✅ A0. Scaffold `src/ng/alignment/` and define `Alignment`.**
`mod.rs` declaring the module's files, wired into `ng/mod.rs`. A **folder**, not a file: it holds
competing implementations (arch §Module home). Plus `Alignment { reference_offset, ops: Vec<CigarOp> }`
reusing production's `CigarOp` rather than minting a parallel operation type — its users are
algorithm 2 (Milestone E, gated) and the normalizers in plan 3. No algorithms. *Source:* arch §Module
home, §2.1.

**✅ A1. `Emission` + its two implementations.**
The trait, plus the per-base-quality implementation (port `EMISSION_LN` **including its
quality-zero floor**, without which a Q0 base annihilates every path through it) and the flat-rate
one. Note the two have **no common source** for the inserted-base score: production's per-quality
path uses `ln(1/4)`, its flat path has no emission for it at all, so the flat value is a decision
recorded at implementation, not a port (arch §2.3). *Depends:* A0. *Source:* spec §4.2, arch
§2.3, §5.

**✅ A2. `RepeatSpan`, `RepeatGeometry`, `RepeatContext`, and the `BestPathAligner` trait.**
`RepeatSpan`'s four cases (both flanks, left only, right only, neither); the geometry with
**measured, never assumed-equal** flank lengths; the per-call context bundling geometry and stutter;
the trait with its `Scratch`, `Output` and `Context` associated types. No algorithms. *Depends:* A1.
*Source:* arch §2.1, §2.2, §3.

**✅ A3. `StutterModel` — the parameters and the distribution.**
The seven parameters and `probability(bp_diff, period)`, following spec §5.2 term by term. Built
**per locus**, from the locus's stutter shape and its *reference* allele length — never from a
candidate's length, which is the genotyping's dimension (arch §5). Unit tests: the distribution
reproduces the published formula term by term; **a test that fails if the geometric is inverted**
(the success-probability-versus-decay trap); a test pinning that `equal` is floored, so "the five
masses sum to one" is *not* asserted (arch §2.4). *Depends:* A2. *Source:* spec §5.2, arch §2.4.

> **Checkpoint A:** the types compile; the emission table reproduces production's values including
> the Q0 floor; the stutter distribution matches the published formula and its two traps are pinned
> by tests. Pause for review.

### Milestone B — algorithm 3, the parity anchor *(this unblocks the STR generator)*

**✅ B1. The two-regime matrix and traceback, unbanded. Own commit; do not bundle.**
Port `delimit_read`: match/insertion/deletion states, affine gaps, and the **tract-aware gap-open**
(stiff in the flanks, soft inside the repeat) keyed by which reference column a gap touches. Fill the
whole matrix, as production does — **unbanded is deliberate** and is what makes parity provable
(spec §3). Carry the two determinism rules: tie-break **match, then deletion, then insertion**, and a
junction gap belongs to the block on its **5′ side**. **Isolated because its failure is silent:** a
wrong recurrence is a wrong repeat length, not a panic. *Depends:* A3. *Source:* spec §4.2, arch §5.

**✅ B2. The `RepeatSpan` readout — the widening. Own commit; do not bundle.**
Walk the traceback and read the repeat off the two flank-junction columns, reporting **which flanks
anchored**. This is where ng goes beyond production: `Delimited::BorderOffEnd` is side-blind, and the
STR preparer cannot otherwise tell a measurement from a lower bound (arch §5). **Isolated because its
failure is silent:** a mis-assigned side is a wrong observation class, and the read still looks fine.
Tests: all four `RepeatSpan` cases reachable; a clean repeat measures exactly; a long allele does not
collapse to the reference; an interrupted repeat comes out verbatim. *Depends:* B1. *Source:* spec
§4.2, arch §2.1, §5.

**✅ B3. Byte-parity against production — the port anchor.**
On the shared fixture, every read whose production result is a `Delimited::Region` must measure the
**same bytes** here. Reads production called `BorderOffEnd` must now land in a one-flank or no-flank
case, and that reclassification is checked by count, since it has no production oracle. *Depends:*
B2. *Source:* spec §10.3, arch §Test & bench shape.

> **Checkpoint B:** algorithm 3 measures repeats byte-identically to production on the fixture, and
> reports which flanks anchored. **The STR locus generator is unblocked** — its Milestone D can
> proceed against this interface. Pause for review.

### Milestone C — banding, as its own change

**✅ C1. Band the delimiter. Own commit; do not bundle.**
Restrict the matrix to a per-read band: a **floor of |read length − reference length|** — forced,
because these are global alignments — plus headroom for a path that strays past the minimum and
returns. The headroom is a named constant carrying its reasoning, **never a literal**, and **never
derived from the slip cutoff**, which is a scoring parameter and would both be too narrow and let the
scoring model blind the ruler (spec §3, §9). **Isolated because its failure is silent:** a too-narrow
band loses long alleles without any error. Two tests, both required: the parity fixture's output is
**unchanged**, and a **long-allele case at the extreme** where a too-narrow band would show.
*Depends:* B3. *Source:* spec §3, §10.3, arch §4, §6.

> **Checkpoint C:** banding changes no output on the fixture and holds at the long-allele extreme;
> the headroom constant carries its derivation. Pause for review.

### Milestone D — algorithm 4, the comparison entry

**✅ D1. The two-penalty aligner.**
Algorithm 3's matrix, with whole-unit slips priced from the `StutterModel` (A3) instead of as a gap
of that many bases. Three things must be right or it fails for its implementation rather than its
premise: gaining and losing are **priced separately**; out-of-frame changes **keep a route**; and the
model reads the shared stutter parameters rather than a second copy. *Depends:* C1. *Source:* spec
§4.2, arch §5.

**✅ D2. The 3-versus-4 comparison on synthetic reads.**
Score both on simulated reads with known repeat lengths, measuring **calibration as well as
accuracy**. Run **period 1 as well as period 2 and above** — the two regimes exercise different parts
of the model and neither is skippable. At period 1, score indels **of the repeat's own base**
separately from indels of **a different base**: only the first is genuine stutter, the second is
mis-routed by the arithmetic in-frame test, and averaged together the two effects cancel. *Depends:*
D1. *Source:* spec §4.2, §10.3.

> **Checkpoint D:** both aligners run on synthetic truth; the comparison reports calibration and
> accuracy, split by period and — at period 1 — by whether the indel is the repeat's own base. Pause
> for review.

### Milestone E — algorithm 2, the general-purpose aligner *(gated)*

> **Gate.** Build this only once the generic path has decided **what marks a region as
> not-to-be-trusted** (`../spec/read_preparation.md` §4). Until then the mode that would call this
> aligner cannot fire, and building it earlier buys nothing (spec §10.3). If the gate is still open at
> Checkpoint D, **stop there**; this milestone moves to read preparation's plan.

**☐ E1. The affine best-path aligner.**
Standard match/mismatch/insertion/deletion under affine gaps, one gap model along the whole
reference, scoring with **per-base qualities** (A1). Returns an `Alignment`. Note this rules out the
wavefront, bit-vector and difference-recurrence cores — they work in fixed costs and cannot consume a
per-base quality (spec §4.1). *Depends:* D2, and the gate. *Source:* spec §4.1, arch §5.

> **Checkpoint E:** the affine aligner re-places reads on a fixture with known answers. **Plan 2 is
> complete.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the module compiles; emission values match production including the Q0 floor; the stutter distribution matches the published formula term by term, with the inverted-geometric trap and the `equal` floor pinned by tests |
| B | **byte-parity of the measured repeat against `delimit_read`** on the shared fixture; all four `RepeatSpan` cases reachable; long allele does not collapse; interrupted repeat verbatim |
| C | the parity fixture's output is **unchanged** by banding, **and** a long-allele case at the extreme |
| D | synthetic reads with known lengths; calibration and accuracy, split by period, and at period 1 split by whether the indel is the repeat's own base |
| E | re-placement on a fixture with known answers |

## Out of scope (next plans)

- **`LogProb`, `MarginalAligner`, algorithms 5 and 6** → plan 2,
  [`alignment_marginal.md`](alignment_marginal.md).
- **`AlignmentNormalizer` and the three normalizers** → plan 3,
  [`alignment_normalization.md`](alignment_normalization.md).
- **The STR read preparer and the tract tally** that consume algorithm 3 →
  [`locus_generation_ssr.md`](locus_generation_ssr.md) Milestone D.
- **Algorithm 2's calling mode** — blocked on the trigger decision (`../spec/read_preparation.md` §4).
