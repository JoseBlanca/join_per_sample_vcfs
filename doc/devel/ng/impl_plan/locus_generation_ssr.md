# ng STR locus generator — implementation plan

**Status:** draft, 2026-07-22. The build order for the **first `LocusGenerator`**: the STR
generator that turns one `SsrSegment` into one `SampleLocusObservations` — fetch the reads over a
tract, align each to read off its repeat, tally the answers. Design is settled in
[`../spec/locus_generation_ssr.md`](../spec/locus_generation_ssr.md) (spec) and
[`../arch/locus_generation_ssr.md`](../arch/locus_generation_ssr.md) (types & interfaces); it
inherits the locus type, the contract, and the error model from
[`locus_generation.md`](locus_generation.md) (its plan, the shared shape). This turns that design
into build order; it is **not** a place for new design — its two remaining open items are
*empirical* (spec §8), settled by measurement on the finished generator, not here.

It is a **port** of production [`src/ssr/pileup/`](../../../../src/ssr/pileup/), adapted at two
seams the spec names (the split of coordinates from bases, and the widened admission gate). The
north-star is byte parity with production's tract tally on the class production keeps.

> **Sequencing — this plan is gated on the ng STR aligner.** The per-read `align_read` (the STR
> path of step 2, `ReadPreparer`) is **specced but not built** — `read_preparation_ssr.md` has no
> code and no impl plan yet. Milestones **A–C need none of it** and can land now; Milestones **D
> and E consume `align_read` as an interface** and are blocked until that generator's own plan
> lands. This plan calls the aligner; it never reimplements it (spec §1 non-goal).

---

## Scope

**In:** the prerequisite `region_typing` rename (`flank_bp` → `bundle_threshold`);
`src/ng/locus_generation/ssr.rs`; `SsrGeneratorConfig` / `DEFAULT_SSR_MAX_READS_PER_LOCUS` /
`SsrGeneratorCounts`; `SsrLocus` + the margin fetch; the ported `Reservoir` + `locus_seed`; the
fetch → align → tally → output transform; `SsrGenerator: LocusGenerator<SsrSegment>`;
`examples/ng_ssr_loci_dump.rs` and the parity oracle.

**Out (later plans / other work):**

- **The ng STR aligner (`align_read`)** — the pair-HMM tract delimiter, step 2's STR path
  (`read_preparation_ssr.md`); its own plan. This generator calls it.
- **The censored likelihood** that makes partial observations usable in genotyping — step 7, and
  it gates their consumption (spec §7 there). This generator only *records* partials.
- **A second STR generator** (GangSTR-style realign-everything) — the bake-off sibling in
  `locus_generation/`, its own plan (spec §7).
- **Widening policy for long alleles** — the aligner's call (spec §7).
- **The empirical experiments** — do partials pay for genotyping; does `flank_bp` equal the bundle
  threshold (spec §8). Sweeps on the finished generator, scored on the standards, not build order.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The prerequisite rename first.** `region_typing` calls its bundle-clustering radius `flank_bp`;
  this generator needs that name for the flank it fetches. Free the name before writing against it
  (spec §7) — a behaviour-preserving pass of its own, so the generator is written against the freed
  name from the start.
- **Transcribe, then change — never both at once.** The port's two silent-failure seams — the
  reservoir **seed** and the margin **coordinates** — are wrong-answer-not-crash risks. Port each
  faithfully, pinned by production's own tests, before it gains anything new; the reservoir seed
  lands as **its own commit** so a `git bisect` can find it if the kept set moves.
- **Reuse over rewrite.** Call `SampleReads::reads_in_region`, the ported `Reservoir` + `locus_seed`,
  `RefSeq::fetch_into`, and `align_read` — and **do not** port `reaches_locus` (its spanning gate
  makes partial observations unreachable, spec §2). The tally is production's, *extended*.
- **Verify against ground truth.** The north-star is **byte parity** with production's
  `SsrLocusObs.observed` on a shared fixture — complete observations, cap disabled, on a fixture
  shallower than the cap (else parity cannot pass by construction, spec §4). Partial observations
  are **new behaviour with no oracle** — measured, not verified (spec §6).
- **Incremental, with pauses.** One milestone, stop for review, next.
- **Container builds.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); native host build at
  completion. **No `bench/`** yet — the second STR generator makes the bake-off; this one is
  measured by the parity oracle + the dump tool (spec §9).

## Preconditions

- **The shared shape is done** ([`locus_generation.md`](locus_generation.md)) —
  `SampleLocusObservations`, `ObservedSequence`, `ReadCoverage`, `LocusKind` / `SsrDetail`,
  `LocusGenerator<S>`, the dispatcher, `LocusGenerationError`.
- **Step 3 is done** — `SsrSegment`, `Motif`
  ([region_typing/segment_criteria.rs](../../../../src/ng/region_typing/segment_criteria.rs)),
  `RegionKind::SsrSegment`.
- **`SampleReads` + `RefSeq`** are in place (merged read ingestion + foundations).
- **The production oracle is readable:** `Reservoir` + `locus_seed` + `MAX_READS_PER_LOCUS` and
  their tests ([fetch_reads.rs](../../../../src/ssr/pileup/fetch_reads.rs)), `reaches_locus`
  ([footprint.rs:223](../../../../src/ssr/pileup/footprint.rs#L223)) — the gate **not** to port,
  `tally` + `SsrLocusObs` ([locus_tally.rs:77](../../../../src/ssr/pileup/locus_tally.rs#L77)),
  `delimit_read` ([alignment.rs:171](../../../../src/ssr/pileup/alignment.rs#L171)), `Locus`
  ([ssr/types.rs:136](../../../../src/ssr/types.rs#L136)); and `examples/ssr_psp_seqdump.rs` as the
  dump-tool exemplar. A shared fixture shallower than the cap for the parity run.
- **NOT yet in place — the blocker:** the ng STR aligner (`align_read`), specced
  `read_preparation_ssr.md`, unbuilt. Milestones **D, E** are gated on it (Sequencing, above).

---

## The steps

### Milestone A — the prerequisite rename (`region_typing`, behaviour-preserving)

**✅ A1. Rename `region_typing`'s `flank_bp` → `bundle_threshold`.**
The single `SsrSegmentCriteria` field (the collapsed bundle radius,
[segment_criteria.rs:560](../../../../src/ng/region_typing/segment_criteria.rs#L560)) + its ~88
call sites, the `region_typing` tests, the CLI arg, and `typed_regions.md` / `typed_regions_cli.md`.
**Own commit; a pure rename, no logic change** — the `region_typing` suite and `.cat` parity are
**byte-identical before and after**. Frees `flank_bp` for this generator's genuine flank use and
reverses the documented collapse (segment_criteria.rs:483). *Source:* spec §7.

> **Checkpoint A:** `bundle_threshold` is the region-typing name; every `region_typing` test and
> `.cat` parity green, output unchanged. Pause for review.

### Milestone B — config, counts, and `SsrLocus` + the margin fetch

**☐ B1. `SsrGeneratorConfig` / counts / the cap constant.**
`SsrGeneratorConfig { flank_bp: Bp, max_reads_per_locus: Option<u32> }`; `SsrGeneratorCounts` (reads
fetched/discarded, complete/partial obs, the no-observation reasons); `DEFAULT_SSR_MAX_READS_PER_LOCUS`
= 1000, **ng's own const, distinct from production's `MAX_READS_PER_LOCUS`** (spec §4); the
`flank_bp <= bundle_threshold` check at generator construction (the region-typing radius, passed in).
No logic beyond the check. *Depends:* A. *Source:* spec §4, arch §1. Verified: `Default`; the flank
check rejects `flank_bp > bundle_threshold`.

**☐ B2. `SsrLocus` + the margin fetch.**
`SsrLocus { segment, tract_with_margin_bases, margin_start }`; fetch the tract ± `flank_bp` through
the generator's own `RefSeq::fetch_into`, **clamped at contig ends** so the span may be shorter than
`2·flank_bp + tract`, **measuring each flank, never assuming** (spec §2). **Own step** — this is the
"adaptation, not a lift", where an off-by-one in the fetched span is a silent wrong-bases. *Depends:*
B1. *Source:* spec §2, arch §1, §5. Verified: the fetched span coordinates match `margin_start` +
length; each flank measured from the clamp, so a tract near a contig end yields **unequal** flanks.

> **Checkpoint B:** the STR config/counts compile; `SsrLocus` fetches its margin with each flank
> measured and clamped. Pause for review.

### Milestone C — the reservoir cap (a faithful port, with the seed trap)

**☐ C1. Port `Reservoir` + `locus_seed` — own commit, do not bundle.**
Transcribe production's reservoir sampler and its seed as ng's own, keyed to
`DEFAULT_SSR_MAX_READS_PER_LOCUS`. **The seed trap:** FNV-1a over the contig **name** bytes folded
with the **0-based** tract start — feed the name and `start - 1`, **not** the `ContigId` or the
1-based start, and assert it at the call site (spec §4). The offer order is fixed to `SampleReads`'
merge order. **Own commit** — a wrong seed is a silently-different kept set that fails the parity
oracle *looking like an aligner bug*, so it must be `git bisect`-able. *Depends:* B. *Source:* spec
§4, arch §5. **Oracle:** production's own reservoir tests
([fetch_reads.rs](../../../../src/ssr/pileup/fetch_reads.rs) —
`locus_seed_is_deterministic_and_distinguishes_loci` and the reservoir-uniformity tests) ported and
green, plus the seed derivation (name + `start - 1`) asserted directly.

> **Checkpoint C:** the reservoir is a byte-faithful port; its determinism and seed derivation are
> pinned by production's tests. Pause for review.

### Milestone D — the transform: fetch → align → tally → output *(gated on `align_read`)*

**☐ D1. Fetch + the relevance gate + the cap.**
Fetch reads over **the tract plus the margin** via `SampleReads::reads_in_region`, admit on
**relevance** (overlap with the query span, which `SampleReads` already applies) — **do not** port
`reaches_locus` — and offer survivors to the reservoir (C). *Depends:* B, C. *Source:* spec §2, arch
§5. Verified: the query span is tract+margin (not the tract); a margin-only read is admitted; no
spanning gate is applied, so partially-covering reads reach the tally.

**☐ D2. Align + tally.**
Call `align_read(&read, &ssr_locus)` per kept read; classify each result **complete** vs **partial**;
tally into `observed_sequences` with the dedup key **`(bases, read_coverage)`** (a complete and a
partial of identical bases stay separate rows); record each observation's `read_coverage` **before
the cap** (so the derived depth is the sample's, not the reservoir's); count the no-observation
reasons into `SsrGeneratorCounts`. The port of `tally`, **extended** with partials and the support
moments. *Depends:* D1, **the ng STR aligner**. *Source:* spec §2, §3, §4, arch §3, §5.

**☐ D3. Output assembly + `SsrGenerator: LocusGenerator<SsrSegment>`.**
Assemble the `SampleLocusObservations`: `region` = the tract coords, `reference_bases` = the tract
**only** (no flanks), `kind = Ssr(SsrDetail { motif, left_flank, right_flank })` split from the
fetched margin, the two read-drop counters; both tables sorted so the result is order-independent.
`begin_segment` clears a produced flag; the first `next_locus` runs the four steps and returns the
locus, the second returns `None` — **one locus per `SsrSegment`, including when no read covers it**
(empty tallies, zeroed counts). *Depends:* D2. *Source:* spec §2, §3, arch §1, §2. Verified: exactly
one locus per segment; a zero-coverage tract yields an empty-but-present locus.

> **Checkpoint D:** the STR generator produces one locus per tract, complete and partial
> observations tallied and tagged. Pause for review.

### Milestone E — the port anchor: dump tool + parity *(gated on D)*

**☐ E1. `examples/ng_ssr_loci_dump.rs` + a committed fixture.**
Following `examples/ssr_psp_seqdump.rs`: positional args, a `#`-prefixed `key=value` counts header,
a TSV column line, tab-separated rows (partial reads on their own rows). Asserted on a small
committed fixture: **one locus per `SsrSegment`** (incl uncovered), **every fetched read accounted
for** (complete / partial / capped / no-observation), **byte-identical across repeated runs**, and
**unchanged when the cap is raised above the fixture's deepest locus** (a cap below must change the
output). *Depends:* D. *Source:* spec §9.

**☐ E2. The parity oracle — the port anchor.**
Against production's `SsrLocusObs.observed` on the shared fixture: every **complete** observation
matches **byte for byte** (bases and count), **with the cap disabled**, on a fixture **shallower
than the cap** — the one precondition, without which parity cannot pass by construction (spec §4,
§6). And **partial observations exist**, which proves the relevance gate (D1) actually admitted the
partially-covering reads. *Depends:* E1, the production oracle. *Source:* spec §6, §9.

> **Checkpoint E:** complete-observation parity with production is green; partials are present and
> measured. **The STR generator is complete.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | `region_typing` suite + `.cat` parity **byte-identical before and after** the rename |
| B | `Default` + flank-check tests; the fetched span coordinates and **unequal clamped flanks** |
| C | **production's reservoir tests** ported and green; the seed derivation (name + `start - 1`) asserted |
| D | query span = tract+margin, relevance gate (margin-only read admitted, no spanning gate); one locus per segment incl zero-coverage |
| E | dump-tool fixture (accounting + determinism + cap-invariance); **byte parity vs production `SsrLocusObs.observed`** (complete class, cap disabled) + partials present |

## Out of scope (next plans)

- **The ng STR aligner (`align_read`)** — `read_preparation_ssr.md`; this plan's blocker and its
  own plan.
- **The censored likelihood** (step 7), **a second STR generator** (bake-off sibling), **widening
  policy** — spec §7.
- **The empirical sweeps** — partials-for-genotyping, `flank_bp` vs the bundle threshold (spec §8);
  scored on `benchmarks/ssr_tomato1`, not build order.
