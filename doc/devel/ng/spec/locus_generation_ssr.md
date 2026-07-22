# ng — the STR locus generator: one tract segment → one locus

*Status: design spec, 2026-07-21. **No code yet.** The first
[`LocusGenerator`](locus_generation.md) — it consumes `SsrSegment` and produces one locus per tract.
Inherits the locus type, the contract, the dispatch and the error model from
[locus_generation.md](locus_generation.md); read that first. Drives the STR read aligner
([read_preparation_ssr.md](read_preparation_ssr.md)). Grounded in production
[src/ssr/pileup/](src/ssr/pileup/). Naming: **STR** in prose, `ssr` in code — and the per-read
operation is **`align_read`**: on the STR path it is a pair-HMM alignment, not the left-align+BAQ
"preparation" the generic path does. `read_preparation_ssr.md` still specs it as the `ReadPreparer`
trait's `prepare_read`; that name should be reconciled to `align_read` for the STR path.*

---

## 1. Scope — goals and non-goals

**What it does.** For one microsatellite tract blessed by the typed-region generator, gather what one
sample's reads say about it: fetch the reads that touch the tract, align each read to read off what
repeat it shows, and tally the answers into a locus.

**Goals.**

1. Produce exactly one `SampleLocusObservations` per `SsrSegment`, **including when no read covers it**.
2. Decide how lower-bound observations are carried — reads that partially cover the tract, proving it
   is *at least* that long but not how long — which `read_preparation_ssr.md` §6 explicitly leaves to
   whoever writes this — §3.
3. Close `read_preparation_ssr.md` §8's last open question: where the read alignment (`align_read`) is
   invoked — §2.
4. Reach parity with production's tract tally on the class production keeps, and say exactly where
   parity stops applying — §6.

**Non-goals.**

- **Consuming partial observations.** They are recorded and deliberately not fed to the likelihood
  until step 7 can score them — because a model that assumes a complete tract reads a partial one as
  a *short* allele (§3).
- **Genotyping.** No likelihood, no prior, no QUAL. This emits evidence.
- **The alignment itself.** Finding a read's tract borders belongs to the aligner
  (`read_preparation_ssr.md`); this generator calls it.

---

## 2. The transform — four steps per segment

**How it meets the streaming contract.** `LocusGenerator` yields loci one at a time
(`locus_generation.md` §4), which a one-locus generator satisfies trivially: `begin_segment` clears a
"already produced" flag, the first `next_locus` runs the four steps below and returns the locus, and
the second returns `None`. All four steps happen inside that first call — a tract is not divisible
into partial results, so there is nothing to stream within it.

**1. Build the locus the aligner needs.** `SsrSegment` carries motif, borders and purity but **no
bases** ([src/ng/region_typing/segment_criteria.rs:240](src/ng/region_typing/segment_criteria.rs#L240)),
while production's `Locus` embedded `ref_bytes` and every downstream function reads them off it
([src/ssr/types.rs:136](src/ssr/types.rs#L136)). So this generator fetches the flanks itself, through
its own `RefSeq` accessor — consistent with the "a generator holds its own accessors" rule
(`locus_generation.md` §4), but it makes the port an **adaptation, not a lift**, and it is the most
likely place for the port to go subtly wrong.

```rust
/// An STR locus ready to genotype: the segment plus the reference bases the aligner
/// aligns against. The ng counterpart of production's `Locus`, split so the coordinates
/// come from region typing and the bases are fetched here.
pub struct SsrLocus {
    pub segment: SsrSegment,
    /// The tract plus its query margin, canonical bases, clamped at contig ends — so this
    /// may be shorter than `2 * flank_bp + tract`, and **each flank must be
    /// measured, never assumed**. Production hit exactly this and threads the two flank
    /// lengths separately.
    pub tract_with_margin_bases: Box<[u8]>,
    /// 1-based position of `tract_with_margin_bases[0]`.
    pub margin_start: Position,
}
```

**2. Fetch the reads.** The query span is **the tract plus the margin — not the tract** — matching
production, which queries the full `ref_bytes` span
([src/ssr/pileup/fetch_reads.rs:175](src/ssr/pileup/fetch_reads.rs#L175)). A read reaching only into
the margin still carries a border anchor, and the aligner needs it.

**Jargon, once — the aligner.** The pair-HMM that `align_read` runs to find where a read's repeat
starts and ends (it *delimits* the tract), by aligning the read against flank + tract + flank
(`read_preparation_ssr.md` §2).

The admission gate is **relevance, not spanning**. Production's `reaches_locus`
([src/ssr/pileup/footprint.rs:223](src/ssr/pileup/footprint.rs#L223)) drops a non-soft-clipped read
that does not bracket both tract ends by 5 bp (soft-clipped reads it admits anyway) — and a read that
runs off mid-tract is exactly what a partial observation is made of. Porting the gate would make
partial observations unreachable (`read_preparation_ssr.md` §4). So the gate here is overlap with the
query span, which `SampleReads` already applies, and nothing further.

**3. Align each read** — `align_read(&read, &ssr_locus)` per kept read. **This closes
`read_preparation_ssr.md` §8: the alignment is invoked here, by this generator, per locus, per read.**
It confirms that spec's leaning. The alternative — aligning reads once in a stream upstream of any
locus — cannot work on the STR path, because the tract-aware aligner needs the locus to place its
gap penalties.

**4. Tally** the per-read observations into the locus (§3), sorting both tables so the result does
not depend on the order reads arrived.

**A tract with no reads still produces a locus** — empty tallies, zeroed counts. *"We looked here and
this sample showed nothing"* and *"we never looked here"* are different facts, and dropping the locus
destroys the difference. Zero coverage is a real observation; at low depth it is most of what a
sample has to say.

---

## 3. The output — what the reads showed about the tract

This generator fills the shared `SampleLocusObservations` (`locus_generation.md` §3); it defines no locus type of
its own. What follows is how a tract's reads map onto that type.

**Jargon, once — complete and partial observations.** A read showing **both** borders of the tract
pins its length: a *complete* observation. A read showing one border and then running off its own end
mid-tract proves the tract is *at least* that long and no more: a *partial* observation. Statistically
that is a **censored** observation — you know the value exceeds a threshold, not what it is, like a
trial patient still alive when the study ends. Production keeps only the complete ones, counting the
rest as `BorderOffEnd` and discarding them
([src/ssr/pileup/locus_tally.rs:91](src/ssr/pileup/locus_tally.rs#L91)); keeping them is the one place
this path is new rather than a port.

| `SampleLocusObservations` field | what the STR generator puts there |
|---|---|
| `region` | the tract's coordinates, from the `SsrSegment` |
| `reference_bases` | the reference tract — the REF sequence, the tract **only**, no flanks (those are in `SsrDetail`) |
| `observed_sequences` | one entry per distinct `(bases, read_coverage)`, sorted by bytes |
| `kind` | `LocusKind::Ssr(SsrDetail { motif, left_flank, right_flank })` — the flanks split out of the fetched `tract_with_margin_bases`, the read model's alignment context |
| `reads_without_observation` | reads that reached `align_read` and yielded nothing |
| `reads_discarded_by_cap` | reads the reservoir dropped (§4) |

Per-position depth is *not* a stored field — `num_obs_along_locus()` derives it from the observations'
`read_coverage`, complete reads counting at every tract position and partials only where they reached
(`locus_generation.md` §3). **But the read-coverage each observation carries must be counted before the
cap**, over the reads that reached the tract — otherwise the derived depth is shaped by the reservoir
rather than the sample, and the windowed statistics read a subsample as the truth. That is the one
place this generator's read-coverage bookkeeping differs from "reads that produced an observation".

**One table, tagged — not a separate table for partials.** The observations live in one list, each
carrying its `read_coverage`. Two reads with the **same bases but different coverage are different
evidence**: a `Complete` `ATATAT` says the allele *is* `ATATAT`, a `Partial` `ATATAT` says it is *at
least* `ATATAT`. So the dedup key is `(bases, read_coverage)`, not bases alone, and the two stay
separate rows — collapse them and the partial gets counted as an exact observation of `ATATAT`,
claiming exact-length evidence it never gave. Two *partials* with the same bases and coverage do merge
into one count: they are the identical constraint, so the count loses nothing.

`PartialLeft`/`PartialRight` records which border the read was flush with — in sequence space a prefix
and a suffix are different constraints (`read_preparation_ssr.md` §3) — and the `u16` records how far
it reached, which drives the derived depth. **This answers `read_preparation_ssr.md` §6**, which left
the carrying of partial observations to this spec.

**We store the observed sequence, not a repeat count.** Two alleles of the same length can differ by
an interior substitution — an interrupted repeat — and a count cannot tell them apart.

**Support: the shared type carries the moments, so STR fills what it already holds.** `ObservedSequence`
carries strand, base-quality and MAPQ moments beside the read count. The STR generator holds the
`MappedRead`s, so filling them is a free by-product — no second pass, and few distinct observations per
tract, so the bytes are negligible here (unlike the generic path's per-position scale). **Nothing
STR-side consumes them today** — production's SSR model gates quality per sample rather than carrying
it ([src/ssr/cohort/types.rs:60](src/ssr/cohort/types.rs#L60)) — and it is *not* established that a
downstream STR filter is the mechanism this caller wants; they are filled on the chance that it is,
not on a plan. **Soft:** unconsumed, and the parity oracle (§6) checks only bytes and counts. Read-position-bias
fields are not among them — those are anchor-relative and degenerate on a tract, so they are the
generic path's, added by the pileup generator (`locus_generation.md` §3).

**Recorded, not fed to a model.** Partial observations are recorded regardless — `num_obs_along_locus()`
needs them for per-position depth (§8). What this generator does *not* do is feed them to a genotyping
model: whether and how a likelihood uses a lower bound — the censored term — is the caller's to design
(step 7, `read_preparation_ssr.md` §6), not this step's. The shared type keeps them out of harm's way
meanwhile: `complete_observations()` is what a likelihood consumes, so a partial is never scored as if
it were complete (`locus_generation.md` §3).

**Licence:** the treatment of partially-covering reads is GangSTR's prior art — implement from Mousavi et al. 2019,
**the paper**, never the vendored source, which is GPL-3-or-later. Same rule for HipSTR (GPL-2) and
TRF-mod (AGPL). freebayes (MIT) and htslib are freely portable.

---

## 4. Depth cap, and where it sits

Production reservoir-samples at most 1000 reads per locus, seeded from the locus so the kept set is
independent of thread count ([src/ssr/pileup/fetch_reads.rs:57](src/ssr/pileup/fetch_reads.rs#L57)).
`SampleReads` does no capping and says so, and `locus_generation.md` §4 leaves it to the generator.
**Adopt the reservoir mechanism — decided; the value is ng's own, and experimental.** The reservoir
bounds the read buffer and, more importantly, the per-locus pair-HMM bill (§5); production's own note
calls its value a calibration placeholder, there only to stop one locus making a worker's task huge
([src/ssr/pileup/fetch_reads.rs:26](src/ssr/pileup/fetch_reads.rs#L26)). ng gives the cap a **separate
constant** from production's, starting at 1000 but free to diverge — this and the other caps are to be
set by experiment, not inherited. The mechanism is on from the start because it changes calls at deep
loci, and switching it on after fixtures exist would rebaseline them.

**The cap sits between fetch and align** — reads are reservoir-sampled as they are fetched, and only
survivors are aligned. That is production's order, and it is what makes the cap bound the pair-HMM
bill rather than just the buffer. The consequence to know: the kept reads include some that will
yield no observation, so effective depth is at most the cap, not equal to it.

> **Two traps, both cheap to avoid and expensive to debug.**
>
> **The seed.** It is FNV-1a over the contig **name** bytes, folded with the **0-based** tract start.
> ng speaks `ContigId` and 1-based positions. Seeding from the id, or from the 1-based start,
> silently produces a *different kept set* — the run stays deterministic and the parity test fails
> for a reason that looks like a aligner bug. Feed it the contig name and `start - 1`, and say so
> at the call site.
>
> **The cap defeats the parity oracle.** Production samples 1000 reads that already passed
> `reaches_locus`; ng samples from the strictly larger relevance-admitted population (§2). At any
> locus deeper than the cap the two kept sets differ, so the observations differ and §6's
> byte-for-byte parity **cannot pass by construction**. Run parity with the cap disabled, on a
> fixture shallower than production's cap. This is not a defect in either design — it is what
> widening the gate means.

Determinism also needs the offer order fixed: the reservoir must see reads in `SampleReads`' merge
order — coordinate order, ties to the lowest file index — which is already guaranteed.

**Config.** Per `locus_generation.md` §7, a generator owns its knobs and takes them at construction.

```rust
pub struct SsrGeneratorConfig {
    /// The flanks fetched either side of the tract — the aligner's anchor and the read
    /// query margin (§2).
    pub flank_bp: Bp,
    /// Reads kept per locus, reservoir-sampled. `None` = no cap.
    pub max_reads_per_locus: Option<u32>,
}
```

**Counts.** The locus records *that* reads yielded nothing (`reads_without_observation`); **why**
they did is this generator's to report, because the reasons are specific to how it reads a tract and
mean nothing to a pileup. Run-level, alongside the shared `LocusCounts`:

```rust
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

*One granularity difference from production.* Production carries the no-observation reasons **per
locus** into the cohort ([src/ssr/cohort/types.rs:39](src/ssr/cohort/types.rs#L39)); ng keeps the
per-locus *total* (`reads_without_observation` on the locus) but breaks it down by reason only at run
level. Per-locus counts answer "is *this* locus untrustworthy?", a genotyping input; if a model wants
that, the reasons go back on the locus as fields of their own. **Soft, and cheap to reverse.**
(Border-off-end reads are *not* lost this way — production counts and drops them, whereas ng keeps
them as partial observations on the locus, which is more, not less.)

**`flank_bp` — and its relation to the bundle threshold.** This is the flank fetched either side of
the tract, for the aligner's anchor and the read query (§2). It must not exceed the width the
typed-region generator guarantees clean: that generator drops a segment unless its flanks are
repeat-free out to `bundle_threshold`, the radius within which a neighbouring repeat clusters it into
a bundle ([src/ng/region_typing/segment_criteria.rs:560](src/ng/region_typing/segment_criteria.rs#L560)).
So **`flank_bp ≤ bundle_threshold`, equal by default** — fetch a wider flank and the query may hit a
neighbouring repeat, leaving the aligner's anchor no longer guaranteed repeat-free. (Region typing
currently calls that radius `flank_bp` too; it is being renamed `bundle_threshold` to free the name
for this genuinely-flank use — see below.)

`DEFAULT_SSR_MAX_READS_PER_LOCUS` is ng's own constant — **not** production's `MAX_READS_PER_LOCUS`
([src/ssr/pileup/fetch_reads.rs:30](src/ssr/pileup/fetch_reads.rs#L30)), so the two can diverge. It
starts at 1000, matching production, but that value is **never measured and soft**: to be set by
experiment.

---

## 5. Cross-cutting concerns

**Memory.** One locus at a time: its margin bases, up to `max_reads_per_locus` reads, and the tally.
The margin fetch drives `evict_before` behind the current locus, so the reference window stays
proportional to one locus plus the next read's forward reach.

**Determinism.** The two places it could leak are the reservoir seed (§4) and tally ordering, and
both are pinned: the seed derives only from the locus, and both tables are sorted before emission.

---

## 6. Reuse over rewrite — the map to production

| what | existing code | ng reuse |
|---|---|---|
| read fetch | `SampleReads::reads_in_region` | reuse as-is; **do not** port `fetch_locus_reads`' spanning gate (§2) |
| depth cap | [src/ssr/pileup/fetch_reads.rs:80](src/ssr/pileup/fetch_reads.rs#L80) | port `Reservoir` + `locus_seed` directly, with the seed trap (§4) |
| margin fetch | [src/ng/ref_seq.rs:142](src/ng/ref_seq.rs#L142) | reuse as-is; replaces `Locus.ref_bytes` |
| read alignment | `align_read` per `read_preparation_ssr.md` | call, do not reimplement |
| tally | [src/ssr/pileup/locus_tally.rs:77](src/ssr/pileup/locus_tally.rs#L77) | model for filling `observed_sequences`; **extended** with partial observations and the support moments (§3) |

**Parity oracle:** production's `SsrLocusObs.observed`, on a shared fixture — every **complete**
observation must match byte for byte, in bases and count. It covers only the complete class, only
with the cap disabled (§4), and not the support moments, which production does not compute. Partial
observations are new behaviour with no oracle: they are measured, not verified.

---

## 7. Deferred, with a recommended home

- **The censored likelihood** that makes partial observations usable — step 7, and it gates their
  consumption (§3).
- **A second STR generator.** The contract exists so alternatives can sit side by side
  (`locus_generation.md` §4); GangSTR-style realignment of every read against a synthetic repeat
  reference is the recorded candidate (`read_preparation_ssr.md` §6). **Home: `src/ng/locus/`**,
  beside this one.
- **Widening policy for long alleles** — inherited from production; whether it stays fixed or becomes
  a knob is the aligner's call (`read_preparation_ssr.md` §8).
- **Rename `region_typing`'s `flank_bp` → `bundle_threshold` — do this *before* coding this generator.**
  Region typing calls its bundle-clustering radius `flank_bp`; this generator uses `flank_bp` for the
  flanks it fetches, a different thing. Renaming the region-typing field to `bundle_threshold` (its
  production name, before ng collapsed the two,
  [segment_criteria.rs:483](src/ng/region_typing/segment_criteria.rs#L483)) frees the name and makes
  the relation `flank_bp ≤ bundle_threshold` legible. **A code change to a built module** (~88 sites,
  tests, the CLI arg, and `typed_regions.md` / `typed_regions_cli.md`), and it reverses that documented
  collapse — so it is its own cargo-verified pass, done first so this generator is written against the
  freed name. **Home: `src/ng/region_typing/`.**

---

## 8. Open questions

- **Do partial observations pay *for genotyping*? — open, empirical.** Not whether to record them:
  they are needed regardless, because `num_obs_along_locus()` counts each read over the positions it
  reached, so a read that runs off mid-tract still contributes depth there — drop it and the
  per-position depth is wrong. What is open is only whether **feeding them to the genotyping
  likelihood** improves calls. Inherited from `read_preparation_ssr.md` §8: the mechanism transfers
  from GangSTR, but the payoff may not — its partially-covering reads catch pathogenic expansions far
  beyond read length, while this catalog is period 2–6 with alleles mostly *inside* a read, where the
  complete observations already fire. *What would settle it:* `benchmarks/ssr_tomato1`, silver recall
  and HipSTR concordance, with partials fed to the likelihood versus not.
- **Should `flank_bp` equal the bundle threshold?** Equal by default, and the clean-flank guarantee
  only reaches that far (§4). Whether a wider flank buys long-allele recovery worth more than the
  weakened anchor guarantee is unexamined. *What would settle it:* count `window_truncated` outcomes
  against flank width on a fixture. Soft; safe to move.

---

## 9. Acceptance test

This generator emits no variant calls, so "done" must not decay into "compiles". The definition of
done is a dump tool over a fixture whose output is inspectable and whose counts are asserted.

**`examples/ng_ssr_loci_dump.rs`**, following `examples/ssr_psp_seqdump.rs`: positional arguments, a
`#`-prefixed `key=value` header carrying the counts, then a bare TSV column line and tab-separated
rows. `depth` sums the reads behind **complete** observations only — the ones that pinned the tract —
and partial reads are shown on their own rows, because conflating them is the mistake §3 exists to
prevent.

```
ng_ssr_loci_dump <reference.fa> <sample.bam> [region]

# ssr_loci=1284 zero_coverage=91 reads_capped=0 reads_without_observation=1662
# obs_complete=35102 obs_partial=1447
contig  start   end     motif  ref_tract   depth  read_coverage  observed      reads
chr1    10442   10461   AT     ATATATAT…   14     complete       ATATATATAT…   9
chr1    10442   10461   AT     ATATATAT…   14     partial:left   ATATAT        2
```

Asserted on a small committed fixture, so this is a regression test and not a demo:

1. **Every `SsrSegment` produces exactly one locus**, including uncovered ones — the emitted count
   equals the typed-region walk's `TypedRegionCounts::ssr_loci` (admitted segments, not bundled
   tracts). This makes the zero-coverage rule (§2) checkable.
2. **Every fetched read is accounted for**: it supports a complete observation, supports a partial
   one, was discarded by the cap, or is counted in `reads_without_observation`.
3. **Complete observations match production byte for byte** on the shared fixture, **with the cap
   disabled** — the parity oracle and its one precondition (§4, §6).
4. **Partial observations exist**, which is what proves the relevance gate (§2) actually admitted the
   partially-covering reads. If there are none, the gate was ported wrong.
5. **Output is byte-identical across repeated runs**, and unchanged when the cap is raised above the
   fixture's deepest locus. A cap *below* the depth must change the output — that is what a cap does,
   so it is not a determinism test.
