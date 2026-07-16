# ng step 3 — the typed-region generator

*Status: design spec (2026-07-15). **No code yet.** Third spec under
[`ng_proposal.md`](ng_proposal.md) (§1 step 3, *The locus stream*), and the first that **integrates**
rather than builds: it stands on `RefSeq` ([`ref_seq.md`](ref_seq.md)), the tandem-repeat scanner
([`ssr_repeat_scanner.md`](ssr_repeat_scanner.md)), and the STR catalog. It closes the three
questions the read-preparation specs deferred to "the router spec" — `LocusWindow`, the `SsrLocus`
shape, and `ReadPreparer::Locus` ([`read_preparation_ssr.md`](read_preparation_ssr.md) §8) — and
retires `locus_router/` from [`../arch/module_layout.md`](../arch/module_layout.md). Amends four
sibling docs; see §9. Code-facing companion: [`../arch/typed_regions.md`](../arch/typed_regions.md).*

*Naming: **STR** in prose, `ssr` in code. **Licensing:** implement from papers, never transliterate
— TRF-mod is AGPL-3, GangSTR GPL-3, HipSTR GPL-2. Nothing here needs their source.*

---

## 1. What it is

Walk the reference end to end — streaming, never loading it all into memory — and cut it into
consecutive pieces, each one a region plus **what the sequence there is**. Every base gets exactly
one label; the pieces tile what was asked for, gap-free. That partition is the design's spine (§2.3).

**This step types regions. It does not decide their fate.** The label is a statement about the
sequence — *this is a tandem array longer than a kilobase*, *this is an STR locus with this motif and
these borders*, *nothing more specific can be said here from the reference alone*. What a consumer
then does with each — genotype it, pile it up, mask it, skip it, call it anyway — is a decision
downstream, and this spec has no opinion on it. The one place that boundary is easy to lose is the
`SsrBundle`, which exists precisely so the decision *can* be taken later, with the evidence in hand,
rather than by deleting it here (§2.4).

### 1.1 Vocabulary

Five words here all meant "a range of DNA", which is why this was hard to read. One word per
concept, each naming what the value *is*:

| word | meaning | code |
|---|---|---|
| **region** | a **physical** piece of DNA: contig + 1-based inclusive range, no genetic claim | `GenomeRegion` |
| **typed region** | a region **+ what the sequence there is** | `TypedRegion { region, kind }` |
| **region kind** | that "what" — one of four | `RegionKind` |
| **locus** | a **genetic** object: has alleles, segregates, gets genotyped | `Locus` (the catalog's, rebased — §4) |
| **tract** | the repeat's own extent inside an STR locus | `Locus::ref_tract()` |
| **window** | the walk's memory unit. Never appears in output (§2.3) | `window_bp` |

**Locus (genetic) vs region (physical) is the split that explains the design**, not just a
convention. Exactly one of the four kinds is genetic:

| kind | is | why |
|---|---|---|
| `SsrLocus` | **locus** | the only kind the reference alone hands you as a genetic object |
| `SsrBundle` | region | real repeats, none with clean flanks, so no locus can be named (§2.4) |
| `Generic` | region | nothing more specific can be said from the reference alone |
| `Satellite` | region | a tandem array too long to be a microsatellite |

That is why the STR kind is special, in one line. It also vindicates `ng::tandem_repeat::Region`:
`Repeat`/`Satellite`/`Unique` are physical classifications, so `Region` was the right word there.

Every typed region is a region of DNA. it is what `ng_step_interfaces.md` §1 already reserved, and it consolidates
`regions::Region` and `bam::ContigInterval` as that doc intended.

**What the user asks for needs no word of its own** — it is a set of genome regions
(`GenomeRegions`, §2.5). The three sets in play are adjectives on one noun: the **requested**
regions, the **scanned** regions (grown by a margin), the **effective** regions (requested, grown to
hold a straddling locus or bundle).

**Jargon, once.** A **tandem repeat** is a stretch of DNA where a short unit repeats back to back
(`CAGCAGCAG…`); the unit is the **motif**, its length the **period**. An **STR** (short tandem
repeat; a.k.a. SSR / microsatellite) is one short enough to genotype from short reads — *where that
line falls is open and empirical* (§5.2, §10), not a definition. **Satellite DNA** is the same
structure at megabase scale. A locus's **flanks** are the unique sequence either side, which is what
lets a read be anchored to it. **Purity** is the fraction of a tract matching a perfect motif
tiling — a degree, not a flag. An **STR bundle** is a cluster of repeats packed so close that none
has clean flanks (§2.4).

### 1.2 Goals, non-goals

**Goals.** Produce the partition. **Make the typing policy a parameter**, so the experiments that
are ng's whole point can move it (§5, §5.2) — reusing the catalog's *implementation* rather than
rewriting it, and starting at its settings only so there is something to compare against (§8). Give
the STR path the locus it asked for — motif, borders, flanks — in a type it can already consume
(§4). Be a pure function of the reference, the regions, and the config.

**Non-goals**, deliberately: **classifying a region handed to us** (nothing is handed to this step —
it walks and produces); **finding loci in generic regions** (the pileup's job, from the data — §9);
**fetching or filtering reads** (this step never opens an alignment file); **data-driven STR
discovery** (deferred by `ng_proposal.md` §step-3's own decision rule); **a bake-off** (§6);
**parallelising the walk** (§7 — a rejection, not a deferral).

**A typed region carries coordinates, not sequence** — with one exception, and the reason is *size, not
need*. Every consumer needs reference bases; the pileup fetches its own. The difference is how much:
a generic region is unbounded, an STR locus's bases are capped at ~1.1 kb by the satellite cap and are
typically ~200 bytes. So `Generic`/`Satellite` carry only their region, `SsrLocus` embeds its tract+flanks,
and an **open generic run costs 8 bytes however many megabases it spans** — which is what makes
§2.1's emission rule affordable.

---

## 2. The walk

### 2.1 The algorithm

Per window, five steps. The first three are the catalog's *implementation*, driven by whatever
settings the config carries (§5) — running them at the catalog's own settings is what makes §8's
oracle a check on the machinery.

1. **Detect** — `find_tandem_repeats(bases, periods, params)` → raw, overlapping candidate
   intervals (coordinates, period, score). No policy.
2. **Clean** — the catalog's pre-filter: per-period copy floor, then period-multiple redundancy
   (one tract is re-detected at every multiple of its period). **Not optional** (§5).
3. **Admit** — `build_loci` → the STR loci *and* the tracts it set aside as bundle members (§5).
4. **Cap** — merge cleaned intervals into coverage runs; a run over `max_repeat_len` (1 kb) is
   satellite. Drop admitted loci *and bundles* inside one.
5. **Partition** — emit `SsrLocus` at each surviving tract, `SsrBundle` across each surviving
   cluster's hull, `Satellite` at each run, `Generic` across the rest.

**Emission lags decision.** A generic run is decided when a repeat ends it; a locus when the walk
has seen far enough past it to know its cluster closed and its coverage run stayed under 1 kb.
**Running out of window decides nothing** — the walk carries the open run across as many windows as
it takes. A generic region flushed at a window edge is a split run, and an indel spanning that join is
callable by neither half (§2.3). The window must never appear in the output.

**A satellite array is typed as one object, not searched for loci inside it.** Satellite runs are
excluded from admission's input, so nothing admitted can land inside one and the one-label-per-base
rule needs no tie-break. That is a **typing** claim — past `max_repeat_len` the tandem structure is
an array, not a microsatellite — and `max_repeat_len` is its parameter, not a constant of nature
(§10 asks whether 1 kb is right). Whether anyone calls a satellite region is their business.

### 2.2 A rejected repeat is generic territory, not a hole

The generic path is the **default**; the other three, STR, STR bundle, and satellite are exceptions carved out of it.

**Flanks reach across borders; ownership does not.** An STR locus embeds ±50 bp of flank, which lies
in generic territory. The partition is over **ownership** — who calls variants at a base — not over
what sequence a span reads. The STR locus emits only tract alleles; the neighbouring `Generic` spans
own and call the flank bases. After §2.4 this is true by definition: a locus is admitted only if no
other repeat lies within `flank_bp` of it.

### 2.3 The invariant — the acceptance test

Within one requested region, the typed regions are **contiguous** (`start == prev.end + 1`, 1-based
inclusive), **non-overlapping**, **complete** (union == the effective region, §2.5), and **maximal**
(no two consecutive share a kind). Across the walk, regions appear in genomic order, each in one
unbroken block, none revisited.

One property: **concatenating the regions reconstructs what was asked for, exactly.** Every way this
design fails shows up as a violation — a rejected repeat left as a hole breaks completeness; a flank
counted as ownership breaks non-overlap; a generic run flushed at a window edge breaks maximality.

**Maximality is a correctness requirement for `Generic`, not tidiness.** A generic region is territory
the pileup mints loci *inside*, so its reach is bounded by the region it was handed. Split a run at *p*
and an indel spanning *p* is callable by neither half — it never appears, and nothing fails. The
windowed walk makes this live: emitting `Generic [core_start, core_end)` per window is the obvious
implementation and would lose every indel across every 100 kb boundary. Two boundaries where
generic regions legitimately stop and must **not** merge: the end of a requested region, and a contig
transition.

**And the partition must not depend on `window_bp`** — the window is a memory knob. This is the
scanner's window-count-invariance gate, extended over admission, and §2.6 earns it.

Nothing overlaps: an admitted locus has clean flanks, so the next repeat either way is ≥ `flank_bp`
off; within a bundle the tracts are *closer* than that, which is why the cluster is one region. A
bundle's span is the hull of its tracts. Zero-length contigs never reach the walk (`GenomeRegions`
drops them).

### 2.4 STR bundles

An STR needs flanks — unique sequence either side to anchor reads against. When two or more repeats
sit so close together that no flank can be built around them, they are not STRs any more: they are
an **STR bundle**, and the region is typed as one.

Algorithm notes:

- **The threshold is `flank_bp`, and there was never a reason for two knobs.** They are two
  histories, not two designs. `bundle_threshold` is GangSTR's `THRESH=50` (`2_trim.sh`), ported —
  and GangSTR's panel build has **no flank concept at all**, so over there the number relates to
  nothing. `flank_bp` is ours: the margin we embed in `ref_bytes`. Both landed on 50 independently,
  and the documented `bundle_threshold >= flank_bp` invariant records that coincidence rather than
  resolving it. Collapse them — **the flank requirement is the primitive, bundle-ness is derived** —
  so when we explore what flank the delimiter actually needs (§10), the bundle definition follows
  for free. It also gives us a reason GangSTR never had: no source comment of theirs says why
  bundles are dropped at all (§10). Since this is the catalog's config, it moves the catalog (§9).
- **Membership is a local test:** a repeat is a bundle member iff another repeat lies within
  `flank_bp` on either side. The cluster falls out of that — there is no separate transitive rule to
  implement. It also selects exactly the records `drop_bundles` sets aside today, so §8's oracle is
  unaffected.
- **The bundle's region is the hull of its tracts** — first start to last end. The gaps between
  members are inside it: they are shorter than a flank, so nothing can be anchored in them either.
- `build_loci` must hand the members back rather than dropping them (§5); carrying an open cluster
  across a window boundary is §2.6.

### 2.5 What to walk — `GenomeRegions`, and why a BED must not change the answer

A user BED is not a special case: `regions.rs` already settled it — *"'Whole genome' is not a special
case — it is the region set whose every region covers an entire contig."*

**The problem, which production does not have.** We compute
during the walk, so a BED edge can cut a decision in half: a repeat at 999 is discarded if a
neighbour sits at 1,030, but with a BED of `[1, 1000]` we never see 1,030 and **admit a locus the
whole-genome run rejects**. Same reference, different calls, because of `--regions`.

**The property: whether a base is STR / bundle / generic / satellite must not depend on the BED.**
The BED chooses what you are *shown*, never what things *are*.

**Decision: scan wider than you emit.** Grow each region by `max_repeat_len`, re-coalesce, scan those,
emit only what overlaps the user's regions. Preprocessing on the region set — the walk needs no special
logic, because a grown region is just a longer continuous run. Two sets: **scan** (grown) and **emit**
(the user's). The residual — a cluster chaining past 1 kb without a 50 bp break — is dense-repeat
territory heading for `Satellite`; **§8 tests for it directly** rather than trusting the argument.

Rejected: scan only the user's region (the bug above); always scan whole contigs (exact, and the
fallback if the test fails on real data, but makes a 10 kb region pay for a 90 Mb chromosome).

**A STR or STR bundle region straddling an edge is emitted whole, and the region grows to hold it.**
That grown region is the **effective region** the invariant is stated over. It terminates in one
step: neither loci nor bundles overlap anything (§2.4), so at most one of them straddles each edge.
**Generic regions are clipped to the user's region.**

### 2.6 Windows have edges; the features don't respect them

You look at ~100 kb at a time. A repeat can sit across your edge, a satellite can be longer than your
whole window, and the code you call was written for a whole contig. Three problems, three fixes.

**1. Read more than you keep.**

```
    fetch ──►│◄─ 1 kb ─►│◄──────── core: 100 kb ────────►│◄─ 1 kb ─►│
                        └─ keep the repeats that START here ─┘
```

Look only at your own 100 kb and a repeat lying across the right edge gets measured as a fragment —
the detector decides where a repeat starts and stops from what it can see, so showing it half a
tract gets you half an answer. Fix: fetch the core **plus 1 kb either side** (`max_repeat_len`), run
detection on the whole slice, then keep only the repeats whose **start** lands in the core.

That gives you both properties at once: each repeat is found **exactly once** (by the one window
that owns its start) and **whole** (an STR is at most 1 kb by definition, so the extra 1 kb always
covers it — and anything longer is a satellite, handled in 3). It also covers admission for free: a
locus needs 50 bp of clean sequence each side to build its flanks, and 1 kb is a lot more than 50.

**This is `collect_windowed` — already written, already tested. Reuse it, don't rewrite it** (§6.1).

**2. Tell `build_loci` where the contig actually ends.**

It works out a locus's right flank as `(new_end + flank_bp).min(contig_seq.len())` and throws the
locus away if that clamped to nothing — that is how it detects "this tract is at the end of the
chromosome, there is nothing to anchor against". Hand it a 100 kb slice and `contig_seq.len()` is
100 kb, so **it believes the chromosome ends at your window edge**: every locus within 50 bp of
every boundary gets silently dropped, and a different set of them for every `window_bp` you pick.

Pass the real contig length in (§5). This is the easiest way to build this wrong and it makes no
noise at all when you do.

**3. Don't emit what you can't decide yet.**

Two questions genuinely can't be answered from the sequence in front of you, and no margin fixes
either — a bigger window just moves the edge. But **both answer themselves if you wait**, so the
rule is simply: keep the unfinished thing in a variable and carry on reading.

- *Is this repeat in a bundle?* You need to see the next repeat, and bundles chain — A–B–C–D each
  30 bp apart runs past any margin you choose. **But a bundle ends the moment a gap of ≥ 50 bp turns
  up.** So hold the open cluster, and resolve it when a far-enough repeat arrives (or the contig
  does). The data tells you when it's over; you don't have to guess a reach.
- *Is this coverage run a satellite?* It is once it passes 1 kb, and a satellite can be megabases.
  **But the answer only moves one way**: a run that reaches 1 kb is a satellite and can never go
  back to being anything else. So carry the open run — at worst you're holding 1 kb of undecided
  territory before the answer forces itself.

So at any moment the walk is holding three unfinished things: the open cluster, the open coverage
run, and the open generic run (§2.1). All three are just coordinates, which is why this costs
nothing (§7).

---

## 3. The types

Live in `src/ng/region_typing.rs`; `GenomeRegion` is shared vocabulary and lands in `ng::types`.

```rust
/// A genome region plus what the sequence there IS — the walk's output.
/// The span is a field, not a per-variant repeat: every typed region has one, structurally, and
/// the invariant (§2.3) reads it off directly. It is the one place ng's 1-based base is stated.
pub struct TypedRegion { pub region: GenomeRegion, pub kind: RegionKind }

/// Exactly one of the four is a genetic object; the other three are physical (§1.1).
/// `Generic` and `Satellite` carry nothing because they *are* just spans.
pub enum RegionKind {
    /// The catalog's `Locus`, used directly — motif, borders, purity, and the embedded
    /// flank+tract+flank bases. No ng wrapper: `Locus` is 1-based like everything else (§4), and
    /// `TypedRegion` already carries the region. It is `ReadPreparer::Locus`, closing
    /// `read_preparation_ssr.md` §8.
    SsrLocus(Locus),                               // Locus = ssr::types::Locus
    /// A cluster of repeats none of which has clean flanks (§2.4). Carries the tracts as
    /// coordinates — enough to see the structure (each interval has its period) without this step
    /// pre-deciding what it is for. The hull is the `TypedRegion`'s own span.
    SsrBundle { tracts: Box<[RepeatInterval]> },   // >= 2, coordinate-ordered
    Generic,
    Satellite,
}

/// Walks the reference: region by region in genomic order, within each in coordinate order, gap-free.
/// Holds one window's bases plus three unfinished coordinates (§2.6) — never a contig, let alone
/// the genome.
/// Concrete over the accessor: the walk needs raw bytes and eviction, which are impl
/// capabilities, not trait methods (§7). Owns its inputs so it can be moved onto a producer
/// thread (§7) — a borrowing walk would need scoped threads or an `Arc`, and nothing else wants
/// them (each worker needs its own accessor, §9).
pub struct TypedRegionIterator { /* … */ }

impl Iterator for TypedRegionIterator {
    /// Every window reads the reference, so a read can fail mid-walk — and a `None` meaning both
    /// "end of the walk" and "a window in chromosome 7 failed" would silently un-call the rest of
    /// the genome. Fused: `Some(Err(_))` once, then `None`.
    type Item = Result<TypedRegion, TypedRegionError>;
}

impl TypedRegionIterator {
    /// The only constructor. Infallible: `GenomeRegions` has already validated everything.
    pub fn over_regions(reference: WindowedRefSeq, spans: GenomeRegions,
                      config: TypedRegionConfig) -> Self;
    /// Running tally (§3.1) — readable mid-walk, complete once exhausted.
    pub fn counts(&self) -> &TypedRegionCounts;
}

/// A fatal, walk-level error. `#[non_exhaustive]`. One variant, because `GenomeRegions` owns
/// everything else that could go wrong (unknown contig, bad BED line, region past a contig end).
pub enum TypedRegionError { Reference(RefSeqError) }

/// A set of genome regions: sorted, non-overlapping, coalesced, clamped, genomic order.
/// Wraps production's `RegionSet`, which already parses BED, coalesces, clamps, converts BED's
/// 0-based text at the one boundary where BED is read, and drops zero-length contigs. The wrapper
/// owns the one u32 → u64 widening ng needs here (§4) rather than scattering it.
pub struct GenomeRegions { /* RegionSet */ }
impl GenomeRegions {
    pub fn whole_contigs(contigs: &[ContigBounds]) -> Self;                  // the default
    pub fn from_bed_path(bed: &Path, contigs: &[ContigBounds]) -> Result<Self, BedError>;
}
```

The contig's name and length come from the `ContigList` the accessor already holds (needs a
`contigs()` accessor added). The length is doubly load-bearing: it fixes the partition's right edge
**and** it is what admission must be *told* rather than infer from its slice (§2.6).

### 3.1 Config and counts

Mirrors `ReadFilterConfig`/`ReadFilterCounts` (`read_filtering.md` §2.4): `Option<T>` for absent,
defaults as named `pub const`s, `Default` = what the lab runs, no dormant knobs.

```rust
pub struct TypedRegionConfig {
    pub periods: PeriodRange,     // starts at 2..=6 — the catalog's setting, for
                                  //   comparability only. A parameter, not a policy (§5)
    pub scan: ScanParams,         // the scanner's defaults (2 / 7 / 2)
    pub max_repeat_len: Bp,       // the satellite cap AND the scan margin — one field because
                                  //   they must be the same number (§2.6)
    pub window_bp: Bp,            // the memory knob; default 100 kb. Must not change the output
    pub catalog: CatalogParams,   // admission's rules — ALL of them, once the hardcoded period
                                  //   scope and copy floors move in (§5). `flank_bp` (50) is now
                                  //   the bundle threshold too, so `bundle_threshold` goes (§2.4)
}

/// Running tally. "No silent caps": a base typed away from the STR path must be accounted for.
pub struct TypedRegionCounts {
    pub spans: u64,                        // requested regions walked
    pub ssr_loci: u64,
    pub ssr_bundles: u64,
    pub ssr_bundle_bp: u64,                // THE number for §10's bundle question — never before
                                           //   seen, because the answer was deleted uncounted
    pub generic: u64,
    pub satellites: u64,
    pub satellite_bp: u64,
    pub repeat_bp_with_no_locus: u64,      // repeat coverage that yielded no locus (§2.2)
    pub rejected_by_reason: RejectionCounts,  // { copy_floor, purity, compound, no_clean_trim,
                                           //   flank_clamped } — no `str_bundle`: that is a route
                                           //   now, not a rejection
}
```

`repeat_bp_with_no_locus` is in **bp, broken out by reason**, because a per-repeat count answers the
wrong question twice: admission trims every survivor, so a repeat that admits one locus and sheds
200 bp contributes nothing to a per-repeat counter; and one total cannot separate a purity rejection
from a copy-floor one, which is exactly the distinction §10 needs. These are the live-caller view of
the catalog's measured ~35% STR coverage gap.

---

## 4. Coordinates

**1-based inclusive, `u64`, ids stay `u32` — everywhere, with no seam.** ng decided 1-based
(`ng_step_interfaces.md` §5). The scanner and the catalog are 0-based half-open today, so **we change
them** (§9) rather than converting at every boundary forever: ng is early enough that consistency is
worth more than the legacy, and a straddle in the coordinate system is the kind of legacy that gets
more expensive every month it survives.

That decision pays for itself immediately. An earlier draft kept `Locus` 0-based and wrapped it in an
`ng::SsrLocus` to present a 1-based surface — a type whose *entire* purpose was to hide the
mismatch. Rebase `Locus` and the wrapper has nothing left to do: `RegionKind::SsrLocus(Locus)` uses
the catalog's type directly, `TypedRegion` already carries the region, and
`ReadPreparer::Locus = ssr::types::Locus` (`read_preparation_ssr.md` §8). **The workaround was
bigger than the fix.**

**`TypedRegion` carries its 1-based region as a struct field**, converted once at construction — one
span, one base, one place. Rejected: repeating it per variant with a `region()` accessor over four,
which makes "every typed region has a region" a convention rather than a fact of the type. (Evidence it
fails: an earlier draft, written against the accessor design, stated the invariant in 0-based
phrasing — in the property the spec calls its spine.)

**`u64` for coordinates and lengths** (owner's call, against this spec's recommendation, recorded so
the reasoning is not re-derived from the outcome). Reason: **simplicity, chosen over economy**. ng is
the lab; one width means no narrowing, no `as`, no checked conversion, no off-by-width bug — worth
more than four bytes on a struct we hold hundreds of. Memory is no argument: the walk holds ~102 kb
of *bases*, and bases are bytes. The rejected alternative was `u32` throughout, on the grounds that
every coordinate ng has *built* is `u32` and that ng addresses by `(ContigId, offset-in-contig)` —
never a genome-wide offset — so only one *contig* must fit in 32 bits and none nears 4.29 Gb. It
loses because it optimises for a port-back `ng_proposal.md` §3 explicitly makes a later problem.

**One finding makes it more than a wash:** `ref_seq.rs:87` already does
`u32::try_from(contig_len).unwrap_or(u32::MAX)` — a >4 Gb contig **silently clamps**, in built ng
code, today. `u64` deletes it rather than guarding it.

**The `u32` holdouts widen in this step, not later** — ng's own (`RefSeq::fetch_into` and its three
impls, the scanner's `RepeatInterval`/`RegionSpan`/`SegmentOptions`, `Bp`) **and the catalog's**
(`Locus`, `build_loci`, `CatalogParams::flank_bp`), which the rebase is touching anyway (§9).
Deferring any of it leaves ng mixed-width, which is the state the decision exists to end.

**That leaves exactly one conversion seam, not two**: `GenomeRegions` widening `RegionSet`'s `u32`
on the way in. Widening the catalog *removes* a seam rather than moving one — `build_loci`/`Locus`
stop converting entirely.

The line stops at `regions.rs` on purpose: it is *"a top-level peer consumed by both pipeline
stages"*, so widening it is a change to the whole production caller, not to ng. (It is also already
inconsistent with itself — `ContigBounds::length` is `u32` while `ContigEntry.length` is `u64` — so
a >4 Gb contig is unrepresentable to `RegionSet` whatever ng does.) The BAM readers are not in scope
either: this step never opens an alignment file (§1.2).

---

## 5. Admission — the policy is a parameter

**Which repeats become STR loci is the experiment, not a fixture.** ng exists to find the best way
to call SNPs, indels, and STRs; where the STR route beats the generic one is one of the things it is
here to measure (§5.2, §10). So this step does not *have* an admission policy — it takes one, and
the partition is a function of it.

**We reuse the catalog's implementation; we do not adopt its policy.** Those are different, and an
earlier draft of this spec collapsed them into a goal ("reuse the catalog's admission policy"),
which was backwards. `postprocess::build_loci`
([postprocess.rs:69](../../../../src/ssr/catalog/postprocess.rs)) is a working, tested
implementation of the whole rule set — period scope, score gate, compound-motif drop, bundle drop,
minimal trim, copy floor, purity floor, flank embed, contig-edge drop — and rewriting it would be
daft. **v1 starts at its settings for one reason only: comparability** (§8). The catalog is a
yardstick, not an authority.

**The rules must therefore be parameters — and half of them are not.** This is a real gap, and it
bites exactly where it hurts:

| rule | how it is set today |
|---|---|
| purity floor, score gate, flank bp, bundle radius | `CatalogParams` — real parameters |
| **period scope** (`MIN_PERIOD` = 2, `MAX_PERIOD` = 6) | **hardcoded `const`** |
| **copy floor per period** (`copy_number_floor`) | **hardcoded `const fn`** — and the pre-filter has a *second*, disagreeing table (§10) |

**The two dimensions §5.2 wants measured — period and length — are precisely the two that are not
knobs.** So the experiment cannot be run against this code as it stands. Moving them into
`CatalogParams` (and reconciling the two copy-floor tables) is a co-requisite of this step, not a
follow-up (§9): a config that cannot express the question is not a config.

**One rule we do re-decide.** The **bundle drop** we keep as a *selection* and reject as a
*disposal* (§2.4) — the same records are set aside, but handed back and routed rather than deleted.
That is a design decision, not a parameter, and it costs no comparability because the selection is
unchanged.

**Three constraints on how it is called.**

**(a) It must be generalised to a window, and must return what it set aside.** Not a behaviour
change:

```rust
// today
fn build_loci(recs: Vec<TrfRecord>, chrom: &str, contig_seq: &[u8], p: &CatalogParams) -> Vec<Locus>

// u64 throughout: the catalog widens with ng (§4), so nothing converts in this call.
fn build_loci(recs: Vec<RepeatInterval>, chrom: &str,
              bases: &[u8], bases_start: u64, contig_len: u64,
              p: &CatalogParams) -> Admitted

struct Admitted {
    loci: Vec<Locus>,                 // exactly what it returns today
    bundled: Vec<RepeatInterval>,     // NEW: the cluster members it silently drops today
}
```

`bases_start` + `contig_len` are the two facts `contig_seq` silently stood in for (§2.6).
**Whole-contig is the degenerate case**, so `catalog::run` is unchanged and byte-identical **by
construction** — that is what preserves the oracle.

**Returning the bundled tracts rather than re-deriving them is what makes this safe.** The obvious
alternative — run the flank test ourselves in the pre-filter, before `build_loci` — **breaks parity
on an ordering subtlety**: compound motifs are dropped *before* bundling, so testing earlier lets a
compound record take out a real neighbour. Silently.

**(b) A pre-filter must run first** — a validated finding, not a guess. `trf-mod` hands `build_loci`
a clean set; our scanner is deliberately permissive (`min_copies = 2`), so it also emits low-copy
noise and every period-multiple of every real tract. Fed straight in, that noise trips the bundle
drop — which runs *before* the copy floor — and cascades the real loci away. With the two cleanups
(`copy floor`, period-multiple redundancy) the scanner reproduces the golden catalog at 16/16.

**(c) `build_loci`'s door is not open yet.** It takes `Vec<TrfRecord>`, and the only bridge from a
scanner interval is `TrfRecord::for_test`, which is `#[cfg(test)]` — so non-test code **cannot call
it at all**. `ssr_repeat_scanner.md` §6 already settled the fix (delete `TrfRecord`, substitute
`RepeatInterval`). **This step is blocked on, or co-lands with, that substitution.**

**A trap in the reused defaults:** `CatalogParams::default()` sets `min_score: 0`, and the field's
comment says why — *"leaves filtering to trf-mod's own `-s 30`"*. There is no `-s 30` here, and
`RepeatInterval::score` is a Ruzzo–Tompa segment total, not a TRF score. So reusing the default
ships **no score gate**. Acceptable — the copy and purity floors are the real gates, and the 16/16
parity ran exactly this way — but it must be known, not inherited by accident.

### 5.1 Where the pre-filter lives

It exists only in a test file (`catalog_prefilter`,
[scanner_parity.rs:59](../../../../src/ssr/catalog/scanner_parity.rs)), and
`ssr_repeat_scanner.md` §6 says it "belongs in `catalog::run`" — written when the catalog swap was
the only consumer. There are now two. Buried inside `catalog::run`, the generator would have to copy
it, and **two copies of an admission policy is how the caller and the catalog silently diverge** —
which would invalidate the oracle this step's validation rests on.

**Recommended: a named `pub(crate)` fn in `src/ssr/catalog/`, called by both.** Still catalog policy;
just not buried in an orchestrator (§9).

### 5.2 Which repeats belong on the STR path — a starting value, not a decision

**v1 starts at period 2–6 with the catalog's copy floors, and that is scaffolding, not a finding.**
It makes v1's admitted set identical to the catalog's, which is what makes the oracle real. It is
**not evidence the grid is right, and this spec should not be read as saying so.**

**The question is two-dimensional and the code already admits it.** A 3-copy and a 40-copy period-3
tract are different problems. The frontier runs over **period × tract length** — and
`copy_number_floor` (period 2 → 5 copies, 3 → 4, rest → 3) bounded by `MIN_PERIOD`/`MAX_PERIOD`
*is* that grid, already drawn by GangSTR for GangSTR's purposes and never measured against our
question. `ng_proposal.md` §step-3 asks it: *"we should check with measurements up to which period
length it's better to go through the STR route."*

**The prior is two-sided.** The experiment that started ng found **freebayes — no STR model at all —
matching or beating HipSTR and us** on single-sample detection: evidence the generic path is
stronger in some cells than an STR specialist's policy suggests. Against that, our indel FNs
concentrate at STR and homopolymer sites. Both true; neither locates the line.

**And the answer per cell is three-way**: route to STR, leave generic, or leave generic and graft
STR-awareness onto it (DRAGstr types STR context on a period 1–8 × repeat-length 1–20 table and
modulates the indel prior instead of forking — `ng_proposal.md`'s *"cheapest way to add STR-awareness
to a general core"*).

**One mechanical fact worth keeping:** scanning 1–6 and 2–6 admit **identical loci**, because the
pre-filter drops period 1 *before* redundancy elimination — via an explicit `iv.period >= 2` clause,
**not** the copy floor (`copy_floor(1)` falls through to the catch-all `3`, which a 3 bp poly-A
clears). That ordering is load-bearing: redundancy treats every period as a multiple of 1, so a
period-1 interval reaching it would annihilate every overlapping real tract. The ranges differ only
in the *partition* — under 1–6 an over-1 kb poly-A becomes `Satellite`. v1 takes 2–6 to keep the
admitted set equal to the catalog's.

---

## 6. Shape, module, memory

**A concrete struct, no trait.** `LocusRouter`, `LocusSource`, and `LocusKind` retire, and
`RefWindow` with them — `route_locus(&RefWindow)` was its only remaining consumer, and "window"
names no concept here (production's `RefSpan` already names a sequence-carrying span if one is ever
wanted). This also closes `LocusWindow` (`read_preparation_ssr.md` §8).

A step is a trait *because* implementations compete (`module_layout.md` principle 1; "a step with no
bake-off is a file, not a folder"). Step 3 has no competitor. Active-region detection — the arch
doc's headline alternative — is **data-driven** and cannot run before any read is fetched; it belongs
inside the pileup's locus definition, and `ng_step_interfaces.md:302` **already relocated its
bake-off there**. Data-driven STR discovery is deferred by the proposal's own rule. What *does* vary
here — period range, satellite cap, admission strictness — are **config knobs, not implementations**:
swept by changing a number, not a `Box<dyn _>`. **`bench/` is deferred**, same reasoning as read
filtering: no competitors, no frontier to plot.

**Module: `src/ng/region_typing.rs`** — a file, promoted to a folder only if it grows.

**Memory: windowed, raw bytes.** Peak is *not* `window_bp`:

| held | size |
|---|---|
| the window's bases | `window_bp + 2 · max_repeat_len` (~102 kb) |
| the open bundle cluster | data-bounded; one locus in ordinary sequence |
| the open coverage run | ≤ `max_repeat_len` before the verdict is forced |
| the open generic run | **two coordinates**, however many megabases |

Two are data-bounded, not constant, and grow exactly where sequence is repeat-dense — which is where
a `Satellite` verdict is coming to release them.

**Raw bytes, not canonical — not a detail.** `RefSeq::fetch_into` serves canonical bases (ambiguity
codes folded to `N`); the catalog reads the FASTA **verbatim** and `build_loci` embeds whatever it is
handed. `Locus` compares by value, so canonical bytes would make every locus containing an IUPAC code
compare unequal to the catalog's — **silently breaking the oracle** on any assembly carrying them.

**Which accessor.** `WindowedRefSeq` is the right shape (sliding buffer + `evict_before` = this
walk's forward slide) but is **canonical-only** by design. `ref_seq.md` parked exactly this: *"Raw
from the windowed impl — YAGNI … add a `RawRefSeq` impl only if a windowed consumer ever actually
needs raw bytes."* **This walk is that consumer** — the YAGNI is spent, not violated. Leaning:
extend `WindowedRefSeq` with a raw path (§10).

**The `--regions` lesson binds in its general form**: one reader for the whole run, sliding forward,
never rebuilt per region — that is what cost 14.6 GB of peak RSS. The walk owns the slide and the
contig transition, which is `ref_seq.md` Decision 6 landing as designed ("the consumer, not the
fetcher, knows when reference bytes are no longer needed" — it named the locus-stream driver, and
this is it).

### 6.1 The scanner's region seam does not fit; its windowed *scan* does

`RegionScanner` was built for this consumer, and this consumer cannot use it: **it merges coverage
and classifies satellite before any policy can run** (§2.4's ordering), on raw permissive intervals,
and there is no way to inject the pre-filter between its detection and its merge.

**The layer below it fits exactly.** `collect_windowed`
([tandem_repeat.rs:530](../../../../src/ng/tandem_repeat.rs)) already does §2.6's hard part — core +
margin, coverage clipped to cores, intervals attributed by start — and its doc already carries the
reasoning this depends on (*"the `max_repeat_len` margins are load-bearing: Ruzzo–Tompa segmentation
is context-dependent"*). It is invariance-tested. **We should not reimplement it.** But it is
private and whole-contig-eager, so: **promote it (crate-visible) and stream it** (an iterator of
`{ core, coverage, intervals }` rather than two contig-wide `Vec`s). This is the previously-deferred
"reshape the seam", promoted to v1.

Rejected: use `RegionScanner` as-is (takes a routing verdict on noise, invisibly); reimplement the
windowed walk here (duplicates a subtle invariant that can then drift — the scanner keeps *policy*
out, and windowing is not policy).

*A seam written for an anticipated consumer did not survive contact with it, and the primitive
underneath it did.* That is the integration slice earning its keep, and an argument for §6's
scepticism about interfaces designed ahead of an implementation.

---

## 7. Cross-cutting

- **Determinism.** A pure function of the reference bases, the regions, and the config. No reads, no
  randomness, no shared state.
- **Concurrency — the walk stays single-threaded, deliberately.** Generating spans is cheap next to
  analysing them, so the parallelism that matters is **downstream, over typed regions**: one producer,
  a bounded queue, N workers. That parallelises **even on a single-contig reference**, which a
  per-contig fan-out cannot, and a single-contig cohort is not a corner case here. **A per-contig
  fan-out is rejected, not deferred**: it speeds up the fast part, caps at the contig count, and buys
  an ordering problem a single FIFO producer does not have.

  **The queue is the pipeline's, not this step's.** The walk is an `Iterator`;
  `for s in iter { tx.send(s?) }` is the whole adapter, and a crossbeam `Receiver` is itself an
  iterator of the same item — so workers see the identical type either way. A channel *inside* the
  walk would make it own a thread and a bounding policy and stop it being `collect()`-able in a
  test. Two things the pipeline must get right: the queue is **bounded** (else the cheap producer
  races ahead and puts the memory back), and the walk **owns** its inputs so it can move onto that
  thread (§3).
- **Performance.** Not a hot spot; measure before optimising. Linear per period per window, plus a
  sort of the interval set. (Honest caveat: the segment-finding pass has a worse-than-linear worst
  case on pathological sequence.)

---

## 8. Tests and the parity oracle

**The oracle checks the machinery, not the policy — and only at one configuration.** `ssr-catalog`
produces a `.cat` whose contents are a `Vec<Locus>`: same type, same reference, and — *if we
configure it so* — the same settings. Run the walk at the catalog's settings, collect every
`SsrLocus`, compare. What that proves is that the plumbing is sound: the windowed scan, the
pre-filter, the admission call, the coordinate arithmetic. **It proves nothing about whether those
settings are right**, and it is not a licence to treat the catalog's numbers as correct (§5).

At any other configuration the two *should* differ — that is what an experiment is (§5.2). So this
is a **fixed-config regression test**, and it must be pinned to the catalog's settings explicitly
rather than to whatever `Default` happens to be, or it will start failing the first time someone
moves a floor and reads it as a bug rather than a result.

`scanner_parity.rs` is the precedent (compares `Locus` sets, tolerant of boundary wobble via span
overlap, readable diff). Prefer a `.cat` built through the **scanner** path, not `trf-mod`'s — that
isolates this step's logic from the detector-swap question `ssr_repeat_scanner.md` §6 already owns.

**At those settings, the one expected divergence is the satellite cap.** Our loci are then a
**strict subset**: identical pipeline, then the cap drops loci inside satellite coverage, which the
catalog deliberately does not apply (*"a live-caller routing feature, not catalog policy"*). Assert
that shape — every catalog locus present, or absent **and inside a satellite run**. At the catalog's
settings, a locus missing for any other reason is a machinery bug. This subset property is *earned*
by §2.4's ordering: capping raw coverage would make the difference bidirectional and untestable.

**`SsrBundle` does not disturb it.** A bundle member is not a locus in either system. The
*selection* is unchanged (§2.4), so the sets stay in the subset relation. Keeping selection and
disposal separable is what preserved the oracle.

**The invariant test** (§2.3) over a synthetic **multi-contig** fixture: a clean tract; an impure and
a low-copy tract (both `Generic` — the §2.2 property); a homopolymer; a 2 kb array (`Satellite`); a
real STR *inside* a 2 kb array (swallowed — the §2.4 cost, made visible); **two tracts 10 bp apart**
(one `SsrBundle` carrying both, *not* two `Generic` regions); **three tracts chained 30 bp apart** (one
bundle of three — emergent transitivity); **a locus 60 bp from a bundle** (admitted; bundles do not
spread); a tract at position 1 (flank clamps → `Generic`); a repeat-free contig; a tract at one
contig's end abutting one at the next's start (transition arithmetic).

**Three invariance tests, each guarding a failure that is otherwise silent:**

- **window-invariance** — identical streams at two `window_bp`. Catches a bundle chain cut at an
  edge, a flank clamped against a window end (§2.6's bug), a satellite not coalesced. The fixture
  must place features **astride window boundaries**, which means a *small `window_bp`*, not a small
  fixture.
- **maximality** — a generic run longer than a window comes out as **one** span. Invariance would
  catch it too, but would report "streams disagree" rather than "every indel across a 100 kb
  boundary is callable by nobody".
- **BED-invariance** — walk whole-genome, walk with a BED; assert the loci inside are identical. The
  fixture must put a repeat cluster **astride the BED edge** (a locus just inside, its 50 bp
  neighbour just outside) — a BED whose edges land in empty sequence proves nothing. Plus: a
  straddling locus comes out whole with both flanks, **and so does a straddling bundle**; generic
  regions stop at the user's edge, not the
  grown scan edge.

---

## 9. Reuse, deferred work, and what this breaks elsewhere

| what | existing code | how |
|---|---|---|
| the windowed scan (core + margin, coverage clipped, intervals by start) | `collect_windowed` ([tandem_repeat.rs:530](../../../../src/ng/tandem_repeat.rs)) | **the primitive this is built on** — must be promoted and streamed (§6.1) |
| repeat detection | `find_tandem_repeats` ([tandem_repeat.rs:354](../../../../src/ng/tandem_repeat.rs)) | called per window by the above |
| the region tiling | `RegionScanner` ([tandem_repeat.rs:586](../../../../src/ng/tandem_repeat.rs)) | **not reused** — merges before policy (§6.1) |
| the admission policy | `build_loci` ([postprocess.rs:69](../../../../src/ssr/catalog/postprocess.rs)) | logic unchanged; windowed + returns `Admitted` (§5a) |
| the bundle selection | `drop_bundles`/`is_close` ([postprocess.rs:274](../../../../src/ssr/catalog/postprocess.rs)) | selection kept, disposal not; re-associated to stream (§2.4, §2.6) |
| the pre-filter | `catalog_prefilter` ([scanner_parity.rs:59](../../../../src/ssr/catalog/scanner_parity.rs)) | **needs a real home** (§5.1) |
| admission knobs | `CatalogParams` ([catalog/mod.rs:42](../../../../src/ssr/catalog/mod.rs)) | reuse; `bundle_threshold` collapses into `flank_bp` (§2.4) |
| the STR locus | `ssr::types::Locus` ([types.rs:136](../../../../src/ssr/types.rs)) | **used directly** — rebased to 1-based so no wrapper is needed (§4, §9) |
| what to walk | `RegionSet`/`Region`/`ContigBounds` ([regions.rs](../../../../src/regions.rs)) | wrapped by `GenomeRegions` (§2.5) |
| reference bases | `WindowedRefSeq` ([ng/ref_seq.rs](../../../../src/ng/ref_seq.rs)) | + a raw path (§6) and a `contigs()` accessor |
| the iterator seam | `ReadFilter` ([read_filtering.md](read_filtering.md) §5) | the shape to match |
| **the parity oracle** | an `ssr-catalog` `.cat` on the same reference | §8 |

**Deferred, with a home.** The **generic join** (splitting a `Generic` region into loci and gathering
reads) → the `pileup/` module and its own spec, the next slice; it also owns "what counts as a
generic locus", which is the real bake-off §6 declined to fake here. The **STR gatherer** (fetch →
prepare → tally) → its own spec. A **region-driven `RecordSource`** — the STR path needs reads *at* a
locus, but `ReadFilter` only drives off a linear scan; production's per-region path is
`AlignmentFile::get_reads_from_segment` and the per-locus precedent is `fetch_locus_reads` → the STR
gatherer's spec, and it is the piece the first integration slice will discover is missing.
**`LocusId`** →
whichever step first needs a stable cross-run id. **Shrinking the two carries** → here,
measure-first.

**What this step also changes, elsewhere.** ng is early. Where the rest of the code does not match
a decision here, **we change the code** rather than work around it — at this stage consistency is
worth more than legacy (owner, 2026-07-16). So most of this list is *work this step includes*, not
consequences to admire from a distance. The `SsrLocus` wrapper is the object lesson: it existed only
to paper over one mismatch, and it was bigger than the fix (§4).

**Rebase `ssr::types::Locus` to 1-based, and widen it to `u64`** — the one that earns the principle. 59 call sites across
`src/ssr/`, but the arithmetic **cancels**: `start`, `end`, and `ref_bytes_start` all shift together
and every slice goes through one `tract_range()`, where only the length gains a `+1`. Exactly one
site does not cancel — `build_loci`'s index into the raw `contig_seq`. The `.cat` format does **not**
move: it is a tabix-indexed BED-flavoured TSV, so convert at `catalog/io.rs`, which is the pattern
`regions.rs` already set. Do the width in the same pass — the sites are already open, and it drops
ng's last catalog-side conversion (§4). Payoff: the coordinate straddle disappears, the wrapper
disappears, `fetch_locus_reads`' inline 0→1 conversion goes away, and `build_loci` stops converting. Guards: the golden catalog fixture
(`tests/data/tandem_repeat/golden.ssr_catalog.bed.gz`), `scanner_parity.rs`,
`benchmarks/ssr_tomato1`, HipSTR concordance. **The risk is real — an off-by-one here is a silently
wrong genotype, not a crash — so it lands as its own commit, green fixture before and after, not
smuggled inside this step.**

- **`ssr_repeat_scanner.md`**, three ways: the pre-filter's home is a shared fn, not `catalog::run`
  (§5.1); the post-filter is *not* consumed "unchanged" — it is windowed and returns `Admitted`
  (§5a); `collect_windowed` must be promoted and streamed while `RegionScanner` goes unused (§6.1).
  Plus its `u32` types widen (§4).
- **`ssr_catalog.md`** — three changes to `CatalogParams`. Its coordinates and `flank_bp` **widen to
  `u64`** with `Locus` (§4), so nothing converts at the admission call. It **loses**
  `bundle_threshold` (§2.4;
  confirm no shipped catalog set it differently from `flank_bp`). And it **gains the rules that are
  hardcoded today**: `MIN_PERIOD`, `MAX_PERIOD`, and `copy_number_floor` must become parameters, and
  the pre-filter's second copy-floor table must be reconciled with them (§5, §10). Without that the
  period × length experiment — the thing this step exists to enable — cannot be expressed, let alone
  run. The catalog keeps today's values as its `Default`, so its own behaviour is unchanged.
- **`ref_seq.md`** — the parked raw-from-windowed YAGNI has fired (§6); `u32` → `u64` (§4).
- **`ng_step_interfaces.md`** — `CallerRecipe` loses its `router` field (§6). Its
  `read_preparer: Box<dyn ReadPreparer>` **cannot compile as written**: `ReadPreparer` has associated
  types, so it is not object-safe, and naming them would pin one path's types and never hold the
  other's. (The paths are selected per *region*, not per *run*, so the field may be the wrong shape
  regardless.) Also: `Position(pub u64)` is confirmed, but it was the only place ng ever wrote
  `u64` — the arch doc and the code have disagreed since the foundations landed. Its reserved
  `GenomeRegion` is **kept as-is** (§1.1), and this step is its first real use — the consolidation
  of `regions::Region` + `ContigInterval` that doc intended.
- **`read_filtering.md` §2.2** — `Bp` widens (§4).
- **`ssr_interrupted_repeat_recall.md`** — bundles are a second bucket of the same ambition
  (recovering loci our filters throw away). Different failure — purity vs anchoring — but both now
  estimate a share of the same ~35% gap, so they should know about each other.
- **The generator and the workers cannot share one accessor.** The walk runs *ahead* (it fills a
  queue) and workers take regions out of order; a sliding-window accessor is built for one monotonic
  forward reader, so the walk's `evict_before` would free bytes a worker is still behind on. Each
  needs its own. Cheap, but not automatic — and it would look like it worked until it didn't. Home:
  the pipeline's spec.

---

## 10. Open questions

- **Which repeats belong on the STR path?** — empirical, and the biggest thing this spec leaves.
  Not "is period 1 in or out": **for which (period, tract length) cells does the STR path beat the
  generic one?** v1's grid is scaffolding (§5.2). **The measurement:** GIAB HG002
  (`benchmarks/ssr_hg002/` — assembly truth + Tier BED), indel recall and genotype concordance
  **binned by period × tract length**, each cell routed generic vs STR. The binning is the point: a
  genome-wide number averages the frontier away, which is plausibly why this is still open.
  `ng_proposal.md` §step-3 asks for exactly this and it has never been run.

  **A trap:** three period-1 gates already disagree — `iv.period >= 2` (the actual gate),
  `copy_floor(1)` → catch-all `3`, `copy_number_floor(1)` → `10`, and `MIN_PERIOD` = 2 running
  *after* redundancy. Relaxing `MIN_PERIOD` alone changes nothing; relaxing the pre-filter clause
  alone re-arms the annihilation hazard (§5.2). The floors live in `postprocess.rs`, so the answer
  moves the catalog too — a reason to measure before touching, not to leave it inherited. **The
  numbers are GangSTR's; the question is ours; nobody has asked it.**

- **What should downstream do with an `SsrBundle`?** The selection and the disposal are settled
  (§2.4); the consumer does not exist. Options, each a real precedent: **treat as generic**;
  **mask** (TRTools' mark-don't-delete); **anchor further out** (freebayes' instinct — it widens its
  window while entropy stays low, commented `// a dangerous game`); or **genotype the cluster
  jointly** as one multi-block haplotype — HipSTR's design, never finished, and the only option that
  *recovers* the loci. **Leaning: none — measure first.** `ssr_bundle_bp` says how much is at stake
  and nobody has ever looked, because the answer was deleted before it could be counted. Home: the
  STR gatherer's spec.

  *The survey behind this (four vendored callers, 2026-07-15), kept because it is the evidence for
  §2.4 and expensive to reproduce — the premise is universal, the response is not:*

  | tool | what it does about a repeat with a repeat in its flank |
  |---|---|
  | **GangSTR** | drops it and the whole transitive cluster, 50 bp, at panel-build time (`2_trim.sh` step 3 → `remove_bundles.py`). **Ours came from here.** Its README markets ver13 as *"More strict removal of locus bundles"* |
  | **HipSTR** | **does not filter its panel.** Its live rule asks the *sequence*: assemble the 35 bp flank into an acyclic de Bruijn graph, abort if you cannot (`seq_stutter_genotyper.cpp:615-628`). `MIN_BLOCK_SPACING = 10` exists but is **dead code** (every `RegionGroup` holds one region) — and that machinery exists to genotype clusters **jointly** |
  | **GATK DRAGstr** | keeps the site; javadoc works the two-abutting-repeats case and tie-breaks (`DragstrReferenceAnalyzer.java:8-33`) |
  | **freebayes** | **extends into** the neighbour (`AlleleParser.cpp:1281-1290`) |
  | **TRTools** | no proximity filter — *same lab as GangSTR, later*, and did not carry the rule forward |

  No GangSTR source comment ever says *why*; the only note (`# NO NEED TO RUN DEDUP WHEN USING THIS
  SCRIPT`) suggests the four-way `is_close` doubles as dedup, so part of the rule's shape may be an
  artefact of a job we do differently (§2.4). Our own `postprocess.rs:12-15` already ties the
  threshold to `flank_bp` — §2.4 only takes that sentence seriously enough to make it the definition.

- **How much flank does the analysis actually need?** `flank_bp` = 50 is inherited and unmeasured,
  and §2.4 just made it load-bearing: it is now the primitive that *defines* an STR bundle, so
  whatever it is set to decides how much of the genome stops being STR loci. The requirement itself
  is the thing to test — the delimiter anchors on both flank junctions today, but *how much* clean
  sequence it needs to do that reliably is an empirical question, and a smaller answer would turn
  bundles back into loci. **Leaning: keep 50 for v1**; `ssr_bundle_bp` (§3.1) says what it costs.
  Note the coupling is now a feature: move `flank_bp` and the bundle definition moves with it, which
  is the whole point of collapsing the two knobs (§2.4).

- **Is 1 kb the right line between "microsatellite" and "array"?** `max_repeat_len` is the scanner's
  default and a round guess; it has never been measured against anything. §2.4's ordering makes it
  safe to reason about (it is applied to cleaned coverage, not detector noise), but that only makes
  the parameter honest, not right. **Leaning: keep 1 kb for v1**; `satellite_bp` and the loci it
  costs are counted so this can be answered rather than assumed.

- **Raw bytes from `WindowedRefSeq`** (§6) — leaning: extend it. Touches built, tested code, and
  whether to keep a second un-canonicalised buffer or serve raw from one is an implementation
  question the arch doc should settle with the code in front of it.
