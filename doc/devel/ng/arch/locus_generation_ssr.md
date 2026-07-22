# ng — the STR locus generator: types & interfaces

*Architecture draft (2026-07-22), the code-facing companion to
[`../spec/locus_generation_ssr.md`](../spec/locus_generation_ssr.md) — the first
[`LocusGenerator`](locus_generation.md), consuming `SsrSegment` and producing one locus per
tract. It **inherits** the locus type, the contract, the dispatch and the error model from
[`locus_generation.md`](locus_generation.md) (arch) — read that first; this doc adds only what is
STR-specific. Grounded in production [`src/ssr/pileup/`](../../../../src/ssr/pileup/). Naming:
**STR** in prose, `ssr` in code. The per-read operation is **`align_read`** — on the STR path it
is a pair-HMM alignment, not the left-align+BAQ "preparation" the generic path does;
`read_preparation_ssr.md` still specs it as `ReadPreparer::prepare_read`, a name to reconcile to
`align_read` for this path (flagged in §5, not fixed here). Signatures are illustrative; the
**contract** is the deliverable.*

## Module home

`src/ng/locus_generation/ssr.rs` — beside `locus_generation/mod.rs`. It **starts** as one file; if it
grows past a comfortable single-file size it splits into an `ssr/` **folder** (as the generic
generator is already `locus_generation/pileup/`) — the file-grows-into-a-folder rhythm
`module_layout.md` principle 3 sets, an internal change that touches nothing outside it. Either
shape is one `LocusGenerator` implementation; a second, *alternative* STR generator (deferred,
spec §7) sits beside it — the bake-off's side-by-side shape (principle 1).

## 1. The types

This generator fills the shared `SampleLocusObservations` (`locus_generation.md` §1); it defines
no locus type of its own. What is STR-local: the working input it builds, the `LocusKind::Ssr`
payload it mints, and its config/counts.

```rust
/// An STR locus ready to align against: the segment plus the reference bases the aligner
/// needs. The ng counterpart of production's `Locus`, **split** so the coordinates come
/// from region typing and the bases are fetched here (spec §2). This split is what makes
/// the port an adaptation, not a lift — the most likely place for it to go subtly wrong.
pub struct SsrLocus {
    pub segment: SsrSegment,
    /// The tract plus its query margin — canonical bases, clamped at contig ends, so this
    /// may be shorter than `2 * flank_bp + tract`. **Each flank must be measured, never
    /// assumed** (production threads the two lengths separately) (spec §2).
    pub tract_with_margin_bases: Box<[u8]>,
    /// 1-based position of `tract_with_margin_bases[0]`.
    pub margin_start: Position,
}

/// The `LocusKind::Ssr` payload this generator mints — the extras a microsatellite locus
/// carries beyond the shared fields. Part of the shared `LocusKind` (ng_step_interfaces.md
/// §2 / spec §3); its fields get their code home here, where the STR path fills them by
/// splitting the fetched `tract_with_margin_bases`.
pub struct SsrDetail {
    /// The repeat unit.
    pub motif: Motif,
    /// The reference flanks the read model aligns the tract against — left and right of the
    /// tract. Lengths are `.len()`; clamped at contig ends, so the two can differ (spec §3).
    pub left_flank: Box<[u8]>,
    pub right_flank: Box<[u8]>,
}

/// This generator's knobs — owned and taken at construction (`locus_generation.md` §5).
pub struct SsrGeneratorConfig {
    /// The flanks fetched either side of the tract — the aligner's anchor and the read
    /// query margin. `flank_bp <= bundle_threshold` must hold; a cross-config relation, so
    /// the generator's constructor checks it (no newtype invariant can) (§3; spec §4).
    pub flank_bp: Bp,
    /// Reads kept per locus, reservoir-sampled. `None` = no cap (spec §4).
    pub max_reads_per_locus: Option<u32>,
}

/// ng's own cap constant — **not** production's `MAX_READS_PER_LOCUS`, so the two can
/// diverge. Starts at 1000 (matching production) but is never-measured and soft: to be set
/// by experiment (spec §4).
pub const DEFAULT_SSR_MAX_READS_PER_LOCUS: u32 = 1000;

/// Run-level counts, alongside the shared `LocusCounts`. The locus records *that* reads
/// yielded nothing (`reads_without_observation`); **why** is this generator's to report,
/// because the reasons are specific to how it reads a tract (spec §4).
pub struct SsrGeneratorCounts {
    pub reads_fetched: u64,
    pub reads_discarded_by_cap: u64,
    pub observations_complete: u64,
    pub observations_partial: u64,
    /// Reads that reached `align_read` and yielded nothing, by reason
    /// (`read_preparation_ssr.md` §4).
    pub no_border_anchored: u64,
    pub low_quality: u64,
    pub window_truncated: u64,
}
```

## 2. The interface

The STR generator is one `LocusGenerator<SsrSegment>` — it holds its own accessors and turns one
tract into one locus.

```rust
/// The STR generator. Holds its own reference accessor, the aligner, and reusable scratch
/// (reservoir buffer, alignment matrices) as fields — the "a generator holds its own
/// accessors" convention (`locus_generation.md` §2; spec §2).
pub struct SsrGenerator<R> { /* reference: R, aligner, reservoir scratch, config, counts */ }

impl<R: RawRefSeq + EvictableRefSeq> LocusGenerator<SsrSegment> for SsrGenerator<R> {
    fn begin_segment(&mut self, region: GenomeRegion);
    fn next_locus(&mut self, segment: &SsrSegment, reads: &SampleReads)
        -> Result<Option<SampleLocusObservations>, LocusGenerationError>;
}
```

**Contract.** Exactly **one** `SampleLocusObservations` per `SsrSegment`, **including when no read
covers it** (empty tallies, zeroed counts — "we looked and saw nothing" ≠ "we never looked", spec
§2). A one-locus generator meets the streaming contract trivially: `begin_segment` clears an
"already produced" flag; the first `next_locus` runs the four steps below and returns the locus;
the second returns `None` (spec §2). The four steps happen inside that first call — a tract is not
divisible into partial results.

The four steps (spec §2): **(1)** build `SsrLocus` — fetch the flanks through the generator's own
`RefSeq` accessor. **(2)** fetch reads over **the tract plus the margin**, admitting on *relevance*
(overlap with the query span), **not** spanning (§3). **(3)** `align_read(&read, &ssr_locus)` per
kept read. **(4)** tally into the locus, both tables sorted so the result is independent of read
arrival order.

## 3. Decisions — decided (why in the spec)

- **`SsrLocus` splits the segment from the fetched bases** — no fork: production embedded
  `ref_bytes` on its `Locus`, ng's `SsrSegment` carries coordinates only, so the generator fetches
  (spec §2).
- **Admit on relevance, not spanning — do NOT port `reaches_locus`** — genuine fork: production's
  gate drops any non-soft-clipped read that fails to bracket both tract ends by 5 bp, which is
  exactly what a partial observation is made of, so porting it makes partials unreachable. The
  gate here is the overlap `SampleReads` already applies, and nothing more (spec §2).
- **`align_read` is invoked here, per locus, per read** — closes `read_preparation_ssr.md` §8.
  *Rejected:* aligning reads once in a stream upstream — the tract-aware aligner needs the locus to
  place its gap penalties, so it cannot run before the locus exists (spec §2).
- **One table tagged by `(bases, read_coverage)`, not a separate partials table** — a `Complete`
  and a `Partial` of the same bases are different evidence (exact vs. lower bound) and stay
  separate rows; two partials with the same bases/coverage merge (identical constraint). This
  answers `read_preparation_ssr.md` §6 (spec §3).
- **Store the observed sequence, not a repeat count** — two alleles of one length can differ by an
  interior substitution (an interrupted repeat), which a count cannot distinguish (spec §3).
- **Adopt the reservoir cap; ng's own constant; cap sits between fetch and align** — bounds the
  per-locus pair-HMM bill, on from the start so fixtures aren't rebaselined later. Two call-site
  traps the spec details: seed from the contig **name** + **0-based** start (`start - 1`, not the
  `ContigId` / 1-based ng speaks), and run parity with the cap **disabled** (it defeats the oracle
  at deep loci) (spec §4).
- **Support moments filled though unconsumed** — the generator holds the `MappedRead`s, so filling
  strand/BQ/MAPQ moments is free; nothing STR-side consumes them today and it is not established
  one will. **Soft** — the parity oracle checks only bytes and counts (spec §3).
- **No-observation *total* on the locus; reasons only at run level** — per-locus reason breakdowns
  are a genotyping input nothing consumes yet; if a model wants them they go back on the locus.
  **Soft, cheap to reverse.** (Border-off-end reads are kept as partials, not dropped — more than
  production, not less.) (spec §4).
- **`flank_bp <= bundle_threshold`, equal by default** — the flank must stay within the width
  region typing guarantees repeat-free; a wider fetch may hit a neighbouring repeat and lose the
  aligner's clean anchor (spec §4).

## 4. Open items

- `OPEN:` **Do partial observations pay *for genotyping*?** — empirical; *not* whether to record
  them (`num_obs_along_locus()` needs them regardless). Inherited from `read_preparation_ssr.md`
  §8; settled by `benchmarks/ssr_tomato1` silver recall + HipSTR concordance, partials fed to the
  likelihood vs. not (spec §8).
- `OPEN:` **Should `flank_bp` equal the bundle threshold?** — equal by default; whether a wider
  flank buys long-allele recovery worth the weakened anchor is unexamined. Settled by counting
  `window_truncated` against flank width on a fixture. Soft (spec §8).
- *Impl-time confirmations, not design items:*
  - **`prepare_read` → `align_read` name reconciliation** — the STR path calls the per-read op
    `align_read` (alignment, not preparation); `read_preparation_ssr.md` still names it
    `prepare_read`. **Flag, do not fix that spec** in this work.
  - **Prerequisite:** region typing's `flank_bp` is renamed **`bundle_threshold`** *before* this
    generator is coded (spec §7) — a ~88-site cargo-verified refactor of a built module, its own
    pass. This doc references `bundle_threshold` as the intended name; **the rename is not done
    here.**

## 5. Reconciliation with existing code

Every row read at the cited line (2026-07-22). Convergence, not new types.

| what | existing code | ng action |
|---|---|---|
| `SsrSegment` (input) | [region_typing/segment_criteria.rs:240](../../../../src/ng/region_typing/segment_criteria.rs#L240) | consume as-is — **coordinates, no bases** ([RegionKind note, mod.rs:171](../../../../src/ng/region_typing/mod.rs#L171)), which is why step 1 fetches |
| `Motif` | [segment_criteria.rs:148](../../../../src/ng/region_typing/segment_criteria.rs#L148) | reuse (ng's own port); the `SsrDetail.motif` value |
| `SsrLocus` ← production `Locus` | [ssr/types.rs:136](../../../../src/ssr/types.rs#L136) (embeds `ref_bytes`, `ref_bytes_start`) | **adapt** — split: coordinates from `SsrSegment`, bases fetched into `tract_with_margin_bases` |
| read fetch | `SampleReads::reads_in_region` [read/input/mod.rs:440](../../../../src/ng/read/input/mod.rs#L440) | reuse as-is; **do not** port `fetch_locus_reads`' spanning gate (read on `main`) |
| the spanning gate NOT to port | `reaches_locus` [ssr/pileup/footprint.rs:223](../../../../src/ssr/pileup/footprint.rs#L223) | **omit** — porting it makes partials unreachable (§3; spec §2) |
| depth cap | `Reservoir` [ssr/pileup/fetch_reads.rs:80](../../../../src/ssr/pileup/fetch_reads.rs#L80) + `locus_seed` [:57](../../../../src/ssr/pileup/fetch_reads.rs#L57) | **port** directly, with the seed trap (§3); `DEFAULT_SSR_MAX_READS_PER_LOCUS` is ng's own, **not** `MAX_READS_PER_LOCUS` [:30](../../../../src/ssr/pileup/fetch_reads.rs#L30) |
| margin fetch | `RefSeq::fetch_into` [ng/ref_seq.rs:142](../../../../src/ng/ref_seq.rs#L142) | reuse as-is; replaces `Locus.ref_bytes` |
| read alignment | `delimit_read` [ssr/pileup/alignment.rs:171](../../../../src/ssr/pileup/alignment.rs#L171) → `Delimited` [:141](../../../../src/ssr/pileup/alignment.rs#L141) | **call** via `align_read` (`read_preparation_ssr.md`), do not reimplement |
| tally | `tally` [ssr/pileup/locus_tally.rs:77](../../../../src/ssr/pileup/locus_tally.rs#L77) → `SsrLocusObs.observed` [:70](../../../../src/ssr/pileup/locus_tally.rs#L70) | **model** for `observed_sequences`; **extended** with partial observations and the support moments (§3) |
| border-off-end reads | `BorderOffEnd` [locus_tally.rs:41](../../../../src/ssr/pileup/locus_tally.rs#L41), counted+dropped [:91](../../../../src/ssr/pileup/locus_tally.rs#L91) | **keep** as partial observations, not drop — the one place this path is new |
| bundle radius (the `<=` bound) | `SsrSegmentCriteria.flank_bp` [segment_criteria.rs:560](../../../../src/ng/region_typing/segment_criteria.rs#L560) | reference as **`bundle_threshold`** (rename pending, §4) — `flank_bp <= bundle_threshold` |
| `SampleLocusObservations` / `ObservedSequence` / `SsrDetail` | [locus_generation.md](locus_generation.md) §1, ng_step_interfaces.md §2 | fill — the shared types |

## Test & bench shape

The definition of done is a dump tool, not "compiles": **`examples/ng_ssr_loci_dump.rs`**
(following `examples/ssr_psp_seqdump.rs`), a `#`-header of counts plus TSV rows, asserted on a
small committed fixture (spec §9). The regression anchors: one locus per `SsrSegment` including
uncovered; every fetched read accounted for; **complete observations match production's
`SsrLocusObs.observed` byte for byte, with the cap disabled** (the parity oracle and its one
precondition, spec §6); and partial observations *exist* (proving the relevance gate admitted the
partially-covering reads). Partial observations are new behaviour with **no oracle** — measured,
not verified.

**Licence.** The partial-covering-read treatment is GangSTR prior art — implement from Mousavi et
al. 2019, **the paper**, never the GPL-3 source (same rule for HipSTR/TRF-mod); freebayes (MIT)
and htslib are freely portable (spec §3).
