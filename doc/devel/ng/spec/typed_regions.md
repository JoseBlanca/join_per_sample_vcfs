# ng step 3 — the typed-region generator

*Status: design spec (2026-07-15). **No code yet.** Third spec under
[`ng_proposal.md`](ng_proposal.md) (§1 step 3, *The locus stream*), and the first that **integrates**
rather than builds: it stands on `RefSeq` ([`ref_seq.md`](ref_seq.md)), the tandem-repeat scanner
([`ssr_repeat_scanner.md`](ssr_repeat_scanner.md)), and the STR catalog. It closes the three
questions the read-preparation specs deferred to "the router spec" — `LocusWindow`, the `SsrSegment`
shape, and `ReadPreparer::Locus` ([`read_preparation_ssr.md`](read_preparation_ssr.md) §8) — and
retires `locus_router/` from [`../arch/module_layout.md`](../arch/module_layout.md). Amends four
sibling docs; see §9. Code-facing companion: [`../arch/typed_regions.md`](../arch/typed_regions.md).
The CLI that drives this walk and writes its output to a file is
[`typed_regions_cli.md`](typed_regions_cli.md) — in `pop_var_caller_exp`, a second binary, so ng's
knobs stay out of the production command.*

*Naming: **STR** in prose, `ssr` in code. **Licensing:** implement from papers, never transliterate
— TRF-mod is AGPL-3, GangSTR GPL-3, HipSTR GPL-2. Nothing here needs their source.*

---

## Revision — 2026-07-16: ng owns its copies; production is frozen

**Owner decision, and it reverses this spec's central reuse argument.** The first draft had ng edit
`src/ssr/` — rebase `Locus`, widen `CatalogParams`, change `build_loci`'s signature — on the grounds
that "ng is early, so we change the code rather than work around it". That was wrong, and the owner
overruled it:

> *"Don't touch `src/ssr`. If you need to, copy the code to `ng/` and modify there whatever you need.
> The objective is to leave production as is and create a fresh ng caller from scratch. In the
> future, once we have done the experiments that ng is supposed to carry out, we'll decide how to
> port the functionalities to production. If you could reuse something without messing production up,
> that's fine, but if you need to change something and is not just a small tweak, let ng have its own
> version."* — owner, 2026-07-16

**Why it is right.** ng exists to *decide* things — the period × length routing frontier, the flank
size, the satellite cap (§10). Binding production to answers we do not have yet buys instability for
nothing, and it destroys the very thing that makes the answers checkable: an **independent yardstick**.
A production catalog we have not touched is a real oracle. A production catalog we rebased to suit ng
is a mirror.

**The rule.** Reuse from `src/ssr/` **only where it costs production nothing** — calling an existing
`pub(crate)` item at its existing signature. The moment reuse would need a signature change, a widened
type, or a coordinate rebase over there, **stop and copy it into `src/ng/`**. ng owning a divergent
copy is the intended outcome, not a smell. Likewise **ng must not depend on `trf-mod`**: it stays
production's catalog detector and may serve as a *comparison oracle*, never as an ng dependency. The
`ssr_repeat_scanner.md` §6–§7 production detector swap is **not** happening.

**What this rewrites below:** §4 (ng's coordinates are ng's; nothing in `src/ssr/` rebases), §5 and
§5.1 (ng **ports** the classification policy rather than reusing and reshaping it in place), §8 (parity
gains a port-fidelity check — what used to be true by construction is now a test), and §9 (the "what
this breaks elsewhere" list largely evaporates, because ng now breaks nothing elsewhere). §1–§3 and
§6–§7 — the walk itself — are unaffected: the algorithm never depended on where the code lived.

**What it costs.** Two copies of the classification policy, which §5.1 rightly warned can silently
diverge. The mitigation is that divergence is now *tested* rather than *prevented*: §8's port-fidelity
oracle pins ng's port against production's `build_loci` on shared inputs. That is a weaker guarantee
than sharing one function, and it is the price of a production tree that cannot be destabilised by an
experiment.

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
| **segment** | a physical stretch typed by what it contains — here, an STR | `SsrSegment` (ng's own, 1-based — §4) |
| **tract** | the repeat's own extent — an `SsrSegment`'s whole span, since it carries no flanks (§1.2) | `SsrSegment::start()..end()` |
| **locus** | a **genetic** object: has alleles, segregates, gets genotyped. **Not produced here** | — (a downstream concept) |
| **window** | the walk's memory unit. Never appears in output (§2.3) | `window_bp` |

**All four kinds are physical — this module produces no genetic object** (owner, 2026-07-18). It
types the genome by *what the sequence is*, and "this stretch carries an STR" is a physical statement,
exactly like "this is a satellite array" or "nothing more specific can be said". A segment becomes a
**locus** — a genetic marker — only if it turns out to be *variable*, which is decided downstream from
the data, not here from the reference:

| kind | is | why |
|---|---|---|
| `SsrSegment` | region | a stretch carrying a tandem repeat with clean flanks — a candidate STR, not yet a locus |
| `SsrBundle` | region | real repeats, none with clean flanks, so none is a usable STR segment (§2.4) |
| `Generic` | region | nothing more specific can be said from the reference alone |
| `Satellite` | region | a tandem array too long to be a microsatellite |

*(An earlier draft made `SsrSegment` "the one genetic kind, a locus". That leaked a downstream,
data-dependent concept into a reference-only typing step — corrected here.)* It also vindicates
`ng::tandem_repeat::Region`: `Repeat`/`Satellite`/`Unique` are physical classifications, so `Region`
was the right word there.

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
are ng's whole point can move it (§5, §5.2) — porting the catalog's *implementation* rather than
re-deriving it, and starting at its settings only so there is something to compare against (§8). Give
the STR path the locus it asked for — motif, borders, flanks — in a type it can already consume
(§4). Be a pure function of the reference, the regions, and the config.

**Non-goals**, deliberately: **classifying a region handed to us** (nothing is handed to this step —
it walks and produces); **finding loci in generic regions** (the pileup's job, from the data — §9);
**fetching or filtering reads** (this step never opens an alignment file); **data-driven STR
discovery** (deferred by `ng_proposal.md` §step-3's own decision rule); **a bake-off** (§6);
**parallelising the walk** (§7 — a rejection, not a deferral).

**A typed region carries coordinates, not sequence — no exceptions** (owner, 2026-07-17). This module
says what each stretch of the genome *is*; the bases are in the reference, and every consumer has it
open, because it had to in order to ask. An earlier draft made `SsrSegment` embed its tract+flanks on
the grounds that they are small (~200 bytes, capped at ~1.1 kb by the satellite cap) — which was an
argument that the payload was *affordable*, never that it was *needed*, and nothing in ng ever read
it. What it cost is in §2.6: a second reference fetch per window and a second margin to make the
embed reachable. **The flank stays as a typing criterion** (§2.4) — it just never needed the bases to
answer.

So every kind carries its region and nothing else but its own description, and an **open generic run
costs 8 bytes however many megabases it spans** — which is what makes §2.1's emission rule affordable.

---

## 2. The walk

### 2.1 The algorithm

Per window, five steps. The first three are ng's port of the catalog's *implementation* (§5), driven
by whatever settings the config carries — running them at the catalog's own settings is what makes
§8's oracle a check on the machinery.

1. **Detect** — `find_tandem_repeats(bases, periods, params)` → raw, overlapping candidate
   intervals (coordinates, period, score). No policy.
2. **Clean** — the ported pre-filter: per-period copy floor, then period-multiple redundancy
   (one tract is re-detected at every multiple of its period). **Not optional** (§5).
3. **Classify** — ng's `classify` (the ported `build_loci`) → the STR loci *and* the tracts it set aside as
   bundle members (§5).
4. **Cap** — merge cleaned intervals into coverage runs; a run over `max_repeat_len` (1 kb) is
   satellite. Drop classified loci *and bundles* inside one.
5. **Partition** — emit `SsrSegment` at each surviving tract, `SsrBundle` across each surviving
   cluster's hull, `Satellite` at each run, `Generic` across the rest.

**Emission lags decision.** A generic run is decided when a repeat ends it; a locus when the walk
has seen far enough past it to know its cluster closed and its coverage run stayed under 1 kb.
**Running out of window decides nothing** — the walk carries the open run across as many windows as
it takes. A generic region flushed at a window edge is a split run, and an indel spanning that join is
callable by neither half (§2.3). The window must never appear in the output.

**A satellite array is typed as one object, not searched for loci inside it.** Satellite runs are
excluded from classification's input, so nothing classified can land inside one and the one-label-per-base
rule needs no tie-break. That is a **typing** claim — past `max_repeat_len` the tandem structure is
an array, not a microsatellite — and `max_repeat_len` is its parameter, not a constant of nature
(§10 asks whether 1 kb is right). Whether anyone calls a satellite region is their business.

### 2.2 A rejected repeat is generic territory, not a hole

The generic path is the **default**; the other three, STR, STR bundle, and satellite are exceptions carved out of it.

**A locus's flank is generic territory, and that is not a conflict.** An STR locus *requires* ±50 bp
of clean sequence either side (§2.4), and those bases belong to the neighbouring `Generic` regions,
which own and call them; the locus's span is its tract alone. The partition is over **ownership** —
who calls variants at a base — and the flank requirement is a *test the locus had to pass*, not
territory it claims. (An earlier draft had the locus embed those bases, which is what made this look
like a conflict worth explaining. It no longer carries them — §1.2.)

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
scanner's window-count-invariance gate, extended over classification, and §2.6 earns it.

Nothing overlaps: an classified locus has clean flanks, so the next repeat either way is ≥ `flank_bp`
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
  nothing. `flank_bp` is ours: the clean sequence a locus must have either side. Both landed on 50 independently,
  and the documented `bundle_threshold >= flank_bp` invariant records that coincidence rather than
  resolving it. Collapse them — **the flank requirement is the primitive, bundle-ness is derived** —
  so when we explore what flank the delimiter actually needs (§10), the bundle definition follows
  for free. It also gives us a reason GangSTR never had: no source comment of theirs says why
  bundles are dropped at all (§10). The collapse happens in **ng's** `SsrSegmentCriteria` (§5); the
  catalog keeps both knobs, because production does not move (Revision).
- **Membership is a local test:** a repeat is a bundle member iff another repeat lies within
  `flank_bp` on either side. The cluster falls out of that — there is no separate transitive rule to
  implement. It also selects exactly the records `drop_bundles` sets aside today, so §8's oracle is
  unaffected.
- **The bundle's region is the hull of its tracts** — first start to last end. The gaps between
  members are inside it: they are shorter than a flank, so nothing can be anchored in them either.
- `build_loci` must hand the members back rather than dropping them (§5); carrying an open cluster
  across a window boundary is §2.6.

### 2.4a A satellite swallows what it touches

**Added at D3, 2026-07-16 (owner's rule), because the walk had no answer here and was quietly giving
two different wrong ones.**

**The rule: a microsatellite or a bundle too close to a satellite is swallowed by the satellite,
which expands to cover it.** "Too close" is the flank measure — fewer than `flank_bp` clean bases in
between — and it subsumes containment, so §2.1's "a locus inside a satellite is dropped" is this same
rule rather than a second one.

**Why it comes up.** A microsatellite 20 bp from a 1 kb array is not genotypeable, and the reason is
the array: the flank on that side *is* array. Two independent routes get you there, which is why the
rule is stated over regions rather than over bundling:

- the array's tract clears classification's gates, so §2.4's flank test bundles the two — and the
  cluster's hull then covers the same bases as the `Satellite`; or
- it does **not** clear them (score, compound motif), so nothing bundles — but the satellite is built
  from *coverage*, which does not care, so a lone `SsrSegment` ends up 20 bp from an array.

**Why the satellite wins.** It is already the type that says *"an array, not a microsatellite — do
not look for loci in here"* (§2.1). A repeat stuck to an array is unusable for the array's reason, so
the array's region should say so. Expanding rather than merely dropping keeps the partition honest:
those bases are accounted for, and by the region that explains them.

**Rejected: exempt the array from the flank test**, so the neighbour becomes a clean locus. Wrong for
a right-sounding reason — a tract 20 bp from a 2 kb array genuinely has no clean flank (§2.4), and
calling it a locus would hand the gatherer a locus whose reads cannot be anchored.

**What it replaced.** The walk tested a cluster's hull **start** only, so the answer depended on which
*side* of the array the microsatellite sat: on the left it emitted a bundle **overlapping** the
satellite — an invalid partition, breaking §2.3's spine; on the right the hull's start fell inside the
run, so the cluster was dropped whole and the bases silently became `Generic`. The same physical
situation, two different wrong answers, and no fixture had reached either.

**Implementation note (D3):** absorption applies to bundle **members**, before they are clustered.
That absorbs whole clusters and never part of one (`is_close` implies "within `flank_bp`", so
absorption chains along the cluster), and it is also what keeps the windowed walk correct — a window
truncates a detection at its slice edge, only an over-cap tract can be truncated, and a truncated
member would break the cluster re-derivation. Absorbed first, it never reaches it.

### 2.5 What to walk — `GenomeRegions`, and why a BED must not change the answer

A user BED is not a special case: `regions.rs` already settled it — *"'Whole genome' is not a special
case — it is the region set whose every region covers an entire contig."*

**The problem, which production does not have.** We compute
during the walk, so a BED edge can cut a decision in half: a repeat at 999 is discarded if a
neighbour sits at 1,030, but with a BED of `[1, 1000]` we never see 1,030 and **classify a locus the
whole-genome run rejects**. Same reference, different calls, because of `--regions`.

**The property: whether a base is STR / bundle / generic / satellite must not depend on the BED.**
The BED chooses what you are *shown*, never what things *are*.

**Decision (owner, 2026-07-17): scan whole contigs, emit what was asked for.** The scan set is the
**whole of every contig the request touches**; the emit set is the user's regions. Two sets still — but
the scan set is not a guess.

> **This replaces "grow each region by `max_repeat_len` and re-coalesce" (E2), which could not deliver
> the emit rule.** The margin was a *guessed reach*, and §2.6 had already said not to guess one:
> *"bundles chain — A–B–C–D each 30 bp apart runs past any margin you choose… the data tells you when
> it's over; you don't have to guess a reach."* §2.5 was the odd section out. Two measurements ended
> it:
>
> - **a 3 kb array under a BED came back `1001–2300` instead of `1001–4001`** — a `Satellite` cut to
>   1300 bp, silently. This is the **common** case, not an exotic one: a satellite is *by definition*
>   longer than `max_repeat_len`, and the margin **is** `max_repeat_len`, so no array running past it
>   could ever be returned whole.
> - **a 1.5 kb chain came back with 24 of its 30 member tracts.** This section used to wave that away —
>   *"the residual … is dense-repeat territory heading for `Satellite`; §8 tests for it directly rather
>   than trusting the argument."* The test was written and the argument was wrong: a satellite is
>   over-cap **coverage**, coverage runs merge only where they **abut** (§2.4), and every run in that
>   chain is 20 bp. It stays a bundle.
>
> **What it costs is time, not memory.** A 10 kb BED on a 90 Mb chromosome now scans 90 Mb — through
> the same window, so peak memory is unchanged (§7), and the emit set is as narrow as ever. What it
> buys is that `--regions` costs nothing in *correctness*: every finding is **the same object** a
> whole-genome run reports, asserted as equality rather than resemblance
> (`a_bed_returns_the_same_findings_the_whole_genome_run_does`).
>
> **It also deletes a rule instead of adding one.** With no scan edges but the contig's, a cluster can
> no longer be cut, so the walk needs no edge case and `bundle_clusters` gets its "a bundle has ≥ 2
> members" assert back at full strength (E2 had to weaken it).

Rejected: **scan only the user's region** (the bug above); **grow each region by a margin** (E2 — a
guessed reach, and no margin can hold a feature defined as longer than it, above).

**Which regions clip, and which come back whole (owner, 2026-07-17).** This section named only a
straddling locus or bundle, and said **nothing about `Satellite`**. The rule for all of them:

> **Every *finding* that intersects a BED edge is returned whole — microsatellite, bundle and
> satellite alike. `Generic` alone is clipped.**

The requested span **grows to hold** a straddling finding, and that grown span is the **effective
region** the invariant is stated over. It terminates in one step: no two findings overlap (§2.4), so at
most one straddles each edge.

`Generic` clips because it is the only kind that is **not a finding**: *"nothing more specific can be
said here"* is true of any stretch of a generic run, so a clipped one says exactly what it said. Each
of the other three is a **claim about its own extent**, and half of it is a different claim — a
locus's coordinates and copy number describe the whole tract; a bundle's hull has to hold its member
tracts; and a `Satellite` means *"an array too long to be a microsatellite"*, so clipping a 1.2 kb
array to a 300 bp request emits a `Satellite` region of 300 bp — **a span that contradicts the
`max_repeat_len` test that produced the label**.

*E2 first ruled the opposite, reasoning from `RegionKind`'s shape: `Generic` and `Satellite` are the
unit variants, carrying no payload, so clipping them could leave nothing misdescribed. That read the
type right and the meaning wrong — **what a region carries is not what it claims**. Recorded because
the argument is a tempting one and is available to be made again.*

*This rule is also what condemned E2's grown scan span: "whole" is a fact about the feature, and a
feature can be longer than any margin — which is where the Decision above came from.*

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
covers it — and anything longer is a satellite, handled in 3).

**This is `collect_windowed` — already written, already tested. Reuse it, don't rewrite it** (§6.1).

> **Corrected at D3, 2026-07-16 — two claims this paragraph used to make were wrong, and both were
> found by writing the test meant to confirm them.**
>
> **(a) "Keep only the repeats whose start lands in the core" — not before the policy runs.** Where
> that sentence sits, it reads as: filter to the core, then classify. Do that and the answer stops being
> window-invariant. Classification's rules are **neighbour** rules — the bundle flank test asks whether
> another repeat is within 50 bp, and the pre-filter's redundancy elimination asks whether an
> overlapping lower-period interval survived. A core tract's neighbour can sit 10 bp beyond the core
> edge, in the margin. Filter first and that neighbour is invisible, so the tract is classified as a
> clean locus instead of bundled — silently, and differently for every `window_bp`.
>
> So the margin is not only the detector's context, it is the **policy's** context, and the order is:
> detect on the slice → pre-filter on the slice → classify on the slice → **then** attribute by start.
> Attribution is the last step, not the first. (Mutation-verified: pre-filtering the core-attributed
> set alone fails `windowed_matches_the_resident_oracle_at_every_window_size`.)
>
> **(b) "It also covers classification for free … 1 kb is a lot more than 50" — true, and for a reason
> this paragraph did not have.** There is **one margin**: the scan's `max_repeat_len`. Classification runs
> over the slice the scan already read, needs each tract's own bases (motif, purity) — which are in
> that slice by construction — and answers the flank question by **arithmetic** against `contig_len`,
> which no slice could answer anyway.
>
> *D3 answered this differently and it is worth one line, because the fix that replaced it deleted
> code rather than adding it. `SsrSegment` used to embed `tract ± flank_bp` of bases, so `classify` read past
> the tract and **panicked** for tracts at the slice's own edge. D3 gave classification a second, wider
> slice — a `flank_bp` margin on top of the `max_repeat_len` one, fetched separately, every window —
> and called them "two margins, two questions". Dropping the payload (§1.2) dropped the second margin,
> the second fetch, and both panics: there was only ever one question.*

**2. Tell classification where the contig actually ends.**

The rule ng inherits from `build_loci` works out a locus's right flank as
`(new_end + flank_bp).min(contig_seq.len())` and throws the
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

Live in `src/ng/region_typing/` (§6): the walk's types in `mod.rs`, ng's ported `SsrSegment` /
`SsrSegmentCriteria` / `classify` in `segment_criteria.rs`. `GenomeRegion` is shared vocabulary and lands in
`ng::types`.

```rust
/// A genome region plus what the sequence there IS — the walk's output.
/// The span is a field, not a per-variant repeat: every typed region has one, structurally, and
/// the invariant (§2.3) reads it off directly. It is the one place ng's 1-based base is stated.
pub struct TypedRegion { pub region: GenomeRegion, pub kind: RegionKind }

/// All four kinds are physical regions (§1.1) — none is a genetic locus.
/// `Generic` and `Satellite` carry nothing because they *are* just spans.
pub enum RegionKind {
    /// ng's own `SsrSegment` — motif, borders, purity. **Coordinates, no bases** (§1.2). A physical
    /// stretch carrying an STR, *not* a locus (it becomes one only if variable, downstream). Ported
    /// from the catalog's `Locus`, born 1-based/`u64` (§4, §5) and without its `ref_bytes`;
    /// production's stays 0-based/`u32` and keeps them. No wrapper: it is 1-based like everything else
    /// in ng, and `TypedRegion` already carries the region. It is `ReadPreparer::Locus`, closing
    /// `read_preparation_ssr.md` §8 — which will need the tract's bases, and will fetch them from
    /// the reference it opens to gather reads.
    SsrSegment(SsrSegment),                               // SsrSegment = ng's, NOT ssr::types::Locus
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
**and** it is what classification must be *told* rather than infer from its slice (§2.6).

### 3.1 Config and counts

Mirrors `ReadFilterConfig`/`ReadFilterCounts` (`read_filtering.md` §2.4): `Option<T>` for absent,
defaults as named `pub const`s, `Default` = what the lab runs, no dormant knobs.

```rust
pub struct TypedRegionConfig {
    pub periods: PeriodRange,     // starts at 2..=6 — the catalog's setting, for
                                  //   comparability only. A parameter, not a policy (§5)
    pub scan: ScanParams,         // the scanner's defaults (2 / 7 / 2)
    pub max_repeat_len: Bp,       // the satellite cap AND the scan margin — one field because
                                  //   they must be the same number (§2.6). It must also be >=
                                  //   classification's flank_bp, or a window cannot see a core tract's
                                  //   bundle neighbours — asserted in the walk (D3), because §10
                                  //   sweeps both knobs and sweeps run in release
    pub window_bp: Bp,            // the memory knob; default 100 kb. Must not change the output
    pub criteria: SsrSegmentCriteria,  // ng's own (§5) — ALL of classification's rules, including the
                                  //   period scope and copy floors the catalog hardcodes. `flank_bp`
                                  //   (50) is the bundle threshold too, so there is no separate
                                  //   `bundle_threshold` knob (§2.4). Default = the catalog's
                                  //   values, for §8's comparability only
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
wrong question twice: classification trims every survivor, so a repeat that classifies one locus and sheds
200 bp contributes nothing to a per-repeat counter; and one total cannot separate a purity rejection
from a copy-floor one, which is exactly the distinction §10 needs. These are the live-caller view of
the catalog's measured ~35% STR coverage gap.

---

## 4. Coordinates

**1-based inclusive, `u64`, ids stay `u32` — everywhere in ng, with one seam.** ng decided 1-based
(`ng_step_interfaces.md` §5). ng's *own* code (the scanner, `RefSeq`, `Bp`) is 0-based/`u32` today, so
**we change it** rather than converting at every boundary forever: ng is early enough that consistency
is worth more than the legacy, and a straddle in the coordinate system is the kind of legacy that gets
more expensive every month it survives.

**Production is not ng's code, and does not move (Revision 2026-07-16).** `ssr::types::Locus` stays
0-based/`u32`; `src/regions.rs` stays `u32`. ng's `SsrSegment` is **its own type** — the same fields, born
1-based and `u64`. The rebase an earlier draft costed at "59 call sites across `src/ssr/`" is not
performed at all: there is nothing to rebase, because ng's copy is written right the first time.

That still buys what the rebase was for. An even earlier draft kept `Locus` 0-based and wrapped it in
an `ng::SsrSegment` to present a 1-based surface — a type whose *entire* purpose was to hide the
mismatch. With ng owning the type the wrapper has nothing left to do: `RegionKind::SsrSegment(SsrSegment)`
holds ng's `SsrSegment` directly, `TypedRegion` already carries the region, and `ReadPreparer::Locus` is
ng's `SsrSegment` (`read_preparation_ssr.md` §8). **The workaround was bigger than the fix** — and the fix
turned out to be cheaper still, since owning the type costs less than rebasing someone else's.

**What ng reuses from `src/ssr/` unchanged: nothing in this step, as it turns out.** `SsrSegment` is
obviously ng's — coordinates and width are all it is. `Motif`
([ssr/types.rs:36](../../../../src/ssr/types.rs)) looked like the opposite case, and this spec called
it *"the Revision's 'reuse where it costs production nothing' case exactly"* — it is `pub(crate)`,
and carries no coordinates and no width, so there is nothing to rebase.

*(Corrected 2026-07-16, A1 review — the premise was true and **the conclusion was wrong**.* `Motif`
*is `pub(crate)`; ng's `SsrSegment` is `pub` and returns one, which trips rustc's `private_interfaces`
lint. The escapes were: widen it in `src/ssr/types.rs` — **touching production, forbidden**; demote
ng's whole classification surface to `pub(crate)` — bends the ng-sibling convention and buys `dead_code`
warnings for every item until its Milestone-D consumer exists; or port 40 coordinate-free lines. So
reuse did **not** cost production nothing — it cost a visibility compromise, which is the coupling
"a fresh ng caller from scratch" exists to avoid. **ng ports `Motif` too.** The dividend: ng's
`region_typing` now names nothing from `src/ssr/` outside its `#[cfg(test)]` differential. The
general lesson for the Revision's rule: "costs production nothing" has to be checked against
**visibility**, not just against coordinates and width.)*

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

**ng's `u32` holdouts widen in this step, not later** — `RefSeq::fetch_into` and its three impls, the
scanner's `RepeatInterval`/`RegionSpan`/`SegmentOptions`, `Bp`. All of it is `src/ng/`, so all of it is
ours to change. Deferring any of it leaves ng mixed-width, which is the state the decision exists to
end. The catalog's `Locus`/`build_loci`/`CatalogParams` are **not** on this list any more — ng's copies
are born wide, and production's stay narrow.

**That leaves exactly one conversion seam**: `GenomeRegions` widening `RegionSet`'s `u32` on the way
in. One seam, one place, at the edge — which is the shape §2.5 wanted anyway.

*(Corrected 2026-07-16 at C3: this said "widening **and rebasing**". There is no rebasing.
`regions::Region` is **already 1-based inclusive** — its own doc says so and its invariant is
`1 <= start <= end` — so production and ng agree on the base already, and the seam widens only. The
spec anticipated an off-by-one at ng's busiest boundary and there is none, because the production
author had made the same call for the same reason. Pinned by `the_seam_widens_but_does_not_rebase`.)*

The line stops at `regions.rs` on purpose: it is *"a top-level peer consumed by both pipeline
stages"*, so widening it is a change to the whole production caller, not to ng. (It is also already
inconsistent with itself — `ContigBounds::length` is `u32` while `ContigEntry.length` is `u64` — so
a >4 Gb contig is unrepresentable to `RegionSet` whatever ng does.) The BAM readers are not in scope
either: this step never opens an alignment file (§1.2).

---

## 5. Classification — the policy is a parameter

**Which repeats become STR loci is the experiment, not a fixture.** ng exists to find the best way
to call SNPs, indels, and STRs; where the STR route beats the generic one is one of the things it is
here to measure (§5.2, §10). So this step does not *have* an classification policy — it takes one, and
the partition is a function of it.

**We port the catalog's implementation; we do not adopt its policy.** Those are different, and an
earlier draft of this spec collapsed them into a goal ("reuse the catalog's classification policy"),
which was backwards. `postprocess::build_loci`
([postprocess.rs:69](../../../../src/ssr/catalog/postprocess.rs)) is a working, tested
implementation of the whole rule set — period scope, score gate, compound-motif drop, bundle drop,
minimal trim, copy floor, purity floor, contig-edge drop — and **re-deriving that rule set from
scratch would be daft**. So ng takes the code. **v1 starts at its settings for one reason only:
comparability** (§8). The catalog is a yardstick, not an authority.

**One of its steps ng does not take: the flank embed** (§1.2). `build_loci` stores each tract's bases
plus a flank each side; ng's copy stores neither, because this module types regions and a payload
nobody here reads is not part of that. Every *decision* is transcribed unchanged, the flank
*requirement* included — that is what keeps §8's oracles meaningful.

**Copied into ng, not reshaped in place (Revision 2026-07-16).** An earlier draft had ng call
`build_loci` directly, which meant editing it — windowing it, widening it, rebasing it, and prising
its hardcoded rules into parameters. Production is frozen, so ng gets **its own copy**, in `src/ng/`,
carrying every change this step needs at once:

| change | why ng needs it |
|---|---|
| 1-based inclusive, `u64` | ng's coordinate system (§4); production stays 0-based/`u32` |
| takes `Vec<RepeatInterval>` | ng's detector is the scanner; `TrfRecord` is trf-mod's shape and ng must not depend on it |
| windowed — `bases_start` + `contig_len` passed in | §2.6's silent bug: a slice cannot be asked where the contig ends |
| returns `Classified { loci, bundled }` | §2.4 routes bundles instead of deleting them |
| **every rule a parameter** | §5.2's experiment (below) |

That is five changes to one function; "a small tweak" it is not, which is exactly when the owner's
rule says copy. **The logic itself is transcribed unchanged** — that is what keeps §8's oracle
meaningful, and what the port-fidelity test pins.

**Every rule must be a parameter — and in production half of them are not:**

| rule | how the catalog sets it | ng's copy |
|---|---|---|
| purity floor, score gate, flank bp, bundle radius | `CatalogParams` — real parameters | parameters |
| **period scope** (`MIN_PERIOD` = 2, `MAX_PERIOD` = 6) | **hardcoded `const`** | **parameters** |
| **copy floor per period** (`copy_number_floor`) | **hardcoded `const fn`** — and the pre-filter has a *second copy* of it (§10) | **one parameter table** (`MinCopies` — named for `min_copies`, the vocabulary ng's scanner already uses; "copy number" means CNV in this crate) |

**The two dimensions §5.2 wants measured — period and length — are precisely the two the catalog
hardcodes.** A config that cannot express the question is not a config, so ng's copy makes them
knobs. This is now a change ng makes *to its own code*, which is the whole point of the owner's rule:
the experiment no longer needs production's permission to run.

**One rule we do re-decide.** The **bundle drop** we keep as a *selection* and reject as a
*disposal* (§2.4) — the same records are set aside, but handed back and routed rather than deleted.
That is a design decision, not a parameter, and it costs no comparability because the selection is
unchanged.

**Three constraints on the port.**

**(a) It is windowed, and it returns what it set aside.** Neither is a behaviour change:

```rust
// production, untouched (src/ssr/catalog/postprocess.rs) — 0-based, u32, whole-contig
fn build_loci(recs: Vec<TrfRecord>, chrom: &str, contig_seq: &[u8], p: &CatalogParams) -> Vec<SsrSegment>

// ng's copy (src/ng/) — 1-based, u64 throughout, so nothing converts in this call.
fn classify(recs: Vec<RepeatInterval>, chrom: &str,
         bases: &[u8], bases_start: u64, contig_len: u64,
         p: &SsrSegmentCriteria) -> Classified

struct Classified {
    loci: Vec<SsrSegment>,                 // ng's SsrSegment (§4) — exactly what build_loci returns today
    bundled: Vec<RepeatInterval>,     // NEW: the cluster members build_loci silently drops today
}
```

`bases_start` + `contig_len` are the two facts `contig_seq` silently stood in for (§2.6).
**Whole-contig is the degenerate case** (`bases_start = 1`, `contig_len = bases.len()`), which is what
makes the port testable against production's whole-contig-only original at all (§8).

**Returning the bundled tracts rather than re-deriving them is what makes this safe.** The obvious
alternative — run the flank test ourselves in the pre-filter, before `build_loci` — **breaks parity
on an ordering subtlety**: compound motifs are dropped *before* bundling, so testing earlier lets a
compound record take out a real neighbour. Silently.

**(b) A pre-filter must run first** — a validated finding, not a guess. `trf-mod` hands `build_loci`
a clean set; our scanner is deliberately permissive (`min_copies = 2`), so it also emits low-copy
noise and every period-multiple of every real tract. Fed straight in, that noise trips the bundle
drop — which runs *before* the copy floor — and cascades the real loci away. With the two cleanups
(`copy floor`, period-multiple redundancy) the scanner reproduces the golden catalog at 16/16.

**The order of the three is load-bearing, and ng's copy fixes an ordering production has wrong**
(2026-07-17; §5.2's correction has the evidence). It is: **copy floor → period-multiple redundancy →
period floor.** Production gates `period >= 2` up front, which drops the period-1 interval before
redundancy elimination can use it — and period 1 is the only period that divides *every* other, so a
homopolymer's `AA` / `AAA` / `AAAAA` aliases then have nothing to eliminate them and enter the cleaned
set as real repeats. Putting the period floor last fixes it; putting the **copy** floor first is what
stops a 2 bp period-1 speck becoming an eliminator and cascading real loci away. Both floors, in that
order, or one of the two failures is back.

**(c) The port is what opens the door.** Production's `build_loci` takes `Vec<TrfRecord>`, and the
only bridge from a scanner interval is `TrfRecord::for_test`, which is `#[cfg(test)]` — so non-test ng
code **could not call it at all**. An earlier draft made this step "blocked on, or co-landing with"
the `ssr_repeat_scanner.md` §6 detector swap that deletes `TrfRecord`. That swap is cancelled
(Revision), and the block dissolves with it: **ng's copy takes `RepeatInterval` because ng wrote it
that way.** The `#[cfg(test)]` bridge survives in exactly one place it is welcome — §8's
port-fidelity test, which is itself `#[cfg(test)]`.

**A trap in the inherited defaults:** `CatalogParams::default()` sets `min_score: 0`, and the field's
comment says why — *"leaves filtering to trf-mod's own `-s 30`"*. There is no `-s 30` in ng, and
`RepeatInterval::score` is a Ruzzo–Tompa segment total, not a TRF score — the two are not on the same
scale, so a threshold carried across from TRF would not mean there what it meant here. So at
`Default` the gate **never fires**: Ruzzo–Tompa emits only positive-scoring segments (`score = r - l >
0`), so `score >= 0` cannot reject anything the scanner produces. Acceptable — the copy and purity
floors are the real gates, and the 16/16 parity ran exactly this way — but it must be known, not
inherited by accident. ng's copy states it on the field, and pins it with a test.

*(Corrected 2026-07-16, A1 review: this said "ships **no score gate**". It ships `score >= 0`, which
is a no-op **for scanner output** but is not nothing — a negative score is dropped. The test written
to prove the original claim disproved it.)*

### 5.1 Where the pre-filter lives

It exists only in a test file (`catalog_prefilter`,
[scanner_parity.rs](../../../../src/ng/scanner_parity.rs) — which B2 moved into ng, see §9), and
`ssr_repeat_scanner.md` §6 says it "belongs in `catalog::run`" — written when the catalog swap was
the only consumer.

**ng copies it, alongside the classification policy it is inseparable from** (§5b — it exists because the
bundle drop runs before the copy floor, which is `build_loci`'s ordering). Production keeps its
test-only copy; ng's lives beside ng's `classify`.

**The warning the earlier draft raised here was real, and is now accepted rather than dodged.** It
said: *"two copies of an classification policy is how the caller and the catalog silently diverge — which
would invalidate the oracle this step's validation rests on."* True. But the alternative it proposed —
one shared `pub(crate)` fn in `src/ssr/catalog/` — requires production to move whenever an experiment
does, which is the coupling the owner's rule exists to prevent. **We take the divergence risk and
test for it** (§8's port-fidelity oracle) instead of designing it away at production's expense.
Divergence caught by a test is a bug; divergence prevented by coupling is a frozen experiment.

### 5.2 Which repeats belong on the STR path — a starting value, not a decision

**v1 starts at period 2–6 with the catalog's copy floors, and that is scaffolding, not a finding.**
It makes v1's classified set identical to the catalog's, which is what makes the oracle real. It is
**not evidence the grid is right, and this spec should not be read as saying so.**

**The question is two-dimensional and the code already classifies it.** A 3-copy and a 40-copy period-3
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

**One mechanical fact worth keeping:** scanning 1–6 and 2–6 classify **identical loci**, because the
pre-filter drops period 1 *before* redundancy elimination — via an explicit `iv.period >= 2` clause,
**not** the copy floor (`copy_floor(1)` falls through to the catch-all `3`, which a 3 bp poly-A
clears). That ordering is load-bearing: redundancy treats every period as a multiple of 1, so a
period-1 interval reaching it would annihilate every overlapping real tract. The ranges differ only
in the *partition* — under 1–6 an over-1 kb poly-A becomes `Satellite`. v1 takes 2–6 to keep the
classified set equal to the catalog's.

> **The paragraph above describes a bug, which is now fixed — 2026-07-17. The sentence about
> `Satellite` was a symptom of it, and the ordering it praises was the cause.**
>
> **The bug.** A homopolymer tiles under every motif length: `AAAA…` is a perfect period-2 `AA`,
> period-3 `AAA` and period-5 `AAAAA` repeat, and the scanner emits the same span at all of them.
> Only the period-1 interval divides them all, so only it can eliminate them — and `prefilter` dropped
> period 1 **by the period floor, before redundancy elimination ever ran**. The aliases outlived their
> only eliminator. Measured on a bare 30 bp poly-A in aperiodic filler: the cleaned set came back with
> **three intervals — periods 2, 3 and 5, identical spans** (4 and 6 died as multiples of 2). They then
> fed coverage (three copies of one span), the satellite cap, and three separate `Compound` rejections.
> **`periods.min() == 2` dropped the period-1 *label*, not the homopolymer** — so the parameter did not
> mean what it said, which is the whole of the bug. That is why an over-1 kb poly-A became a
> `Satellite`: not by decision, but because its aliases were counted as tandem repeats.
>
> **The fix** ([`segment_criteria::prefilter`](../../../../src/ng/region_typing/segment_criteria.rs)) splits the one
> filter in two and puts the period floor **last**: copy floor → redundancy elimination → period floor.
> Period 1 now survives long enough to eliminate its own aliases and is dropped afterwards, so a
> homopolymer contributes **nothing** at 2..=6 and is classified as a period-1 tract under its own motif
> at 1..=6. **The range means what it says in both directions**, which is all a caller asked for.
>
> **The copy floor is what makes that safe, and it is why the ordering was as it was.** A surviving
> period-1 interval eliminates every real STR it overlaps — and aperiodic sequence is full of 2 bp
> period-1 specks (`GG`, `TT`: 45 of them per 300 bp of the test filler). `MinCopies`' period-1 entry
> is **10**, so the specks are dropped before they can eliminate anything and only a genuine
> homopolymer run becomes an eliminator. That entry was **unreachable dead weight** until this
> ordering — see the A1-review note below, which found the two copy-floor tables differ only at period
> 1 and that period 1 reached neither. It is now exactly what separates a poly-A from a speck.
>
> **What the scanner's floor is actually for**, correcting this section's *reason* rather than its
> conclusion: **an eliminator that was never detected cannot eliminate.** Scan from 2 and the aliases
> return, so the scan floor of 1 is load-bearing after all — as machinery, not policy. It is not the
> knob a caller means by "which periods to analyse", which is why
> [`typed_regions_cli.md`](typed_regions_cli.md) §2.2 exposes one period range and pins the scan floor
> at 1.
>
> **Evidence.** All 642 ng tests pass, including `.cat` parity against the trf-mod golden catalog
> (16/16), the port-fidelity differential (§8.0), window-invariance, and the partition invariant over
> the golden reference — so the fix changed no answer any oracle checks, and removed the aliases. Two
> tests failed and were rewritten *because they encoded the bug*:
> `prefilter_drops_period_one_before_it_can_eliminate_a_real_str` asserted the alias should be **kept**
> (its "real STR" was the poly-A's own period-2 alias, over the identical span), and the walk-level
> homopolymer fixture asserted the run must reach classification as `Compound`.
>
> **Consequence for §3.1's counts.** `Compound` was reachable *only* via this bug — the E1e note in §8
> below says so in as many words ("the poly-A cascade's fix is what creates this gate's only
> customer"). With the alias gone, **the walk reaches one of classification's five gates, not two**:
> `FlankClamped`. E1e's original *reasoning* about `Compound` — that a compound motif is a
> period-multiple and redundancy drops it — was right all along; it was false only of period 1, and
> only because of the ordering.

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
here — period range, satellite cap, classification strictness — are **config knobs, not implementations**:
swept by changing a number, not a `Box<dyn _>`. **`bench/` is deferred**, same reasoning as read
filtering: no competitors, no frontier to plot.

**Module: `src/ng/region_typing/`** — a folder, and the Revision is why. A step with no bake-off is a
file, not a folder (`module_layout.md` principle 1), and step 3 has no bake-off — but it now also
carries a ~500-line port of the classification policy (§5), which is a *different* concern from the walk
and has its own dense test suite. So: `region_typing/mod.rs` for the types and the walk,
`region_typing/segment_criteria.rs` for the port (ng's `SsrSegment`, `SsrSegmentCriteria`, the pre-filter, `classify`,
and §8.0's differential test against production). The folder tracks the code's size, not a bake-off —
the `tandem_repeat.rs` precedent, promoted on the same grounds.

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
handed. `SsrSegment` compares by value, so canonical bytes would make every locus containing an IUPAC code
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

### 8.0 Port fidelity — the oracle the Revision made necessary

**What used to be true by construction is now a test.** The earlier design had ng *call* production's
`build_loci`, so "ng classifies what the catalog classifies" needed no proving — it was the same function.
ng now owns a copy (§5), and a copy can drift: on transcription, and later, as experiments edit ng's
side. §5.1 named this the real cost of the Revision. This is the mitigation.

**Differential test: ng's `classify` vs production's `build_loci`, same inputs, same answers.** Both are
in this crate, so the test needs nothing new from production — `TrfRecord::for_test` is
`#[cfg(test)] pub(crate)`, and a test module is exactly where it belongs (§5c):

1. Take a set of `RepeatInterval`s (scanner output on the synthetic fixture, plus the crafted cases
   from `postprocess.rs`'s own tests).
2. Feed ng's `classify` directly; feed production's `build_loci` the same intervals bridged through
   `TrfRecord::for_test`, at the **whole-contig degenerate case** (§5a) and the catalog's settings.
3. Assert the `SsrSegment` sets are identical **modulo the coordinate base** — ng's 1-based inclusive
   `[start, end]` is production's 0-based half-open `[start, end)` shifted by exactly one on `start`,
   with `end` unchanged. Compare through one conversion helper, stated once, so the test pins the
   arithmetic rather than restating the bug.
4. **Do not compare `ref_bytes` or the flanks**: production has them, ng does not (§1.2), and that is
   a divergence rather than drift. What the differential pins is every *decision* — which tracts
   classify, at what span, motif, period and purity — which is what "faithful" has to mean once the
   payload is gone.

**What it proves and what it does not.** It proves the transcription is faithful *at the catalog's
settings, whole-contig*. It does **not** cover the windowed path (§8's window-invariance does) or any
other configuration (nothing can — that is §5.2's experiment). It is a **fixed-config regression
test**, and it is expected to be *deleted or re-pinned* the day an experiment deliberately moves ng's
classification away from the catalog's rules. Until then it is the tripwire on silent drift.

**It also outlives the port.** Production is frozen *now*; when the experiments conclude and the
port-back decision is taken (`ng_proposal.md` §3), this test is the thing that says what ng changed
and what it merely moved.

### 8.1 The `.cat` parity oracle

**The oracle checks the machinery, not the policy — and only at one configuration.** `ssr-catalog`
produces a `.cat` whose contents are a `Vec<SsrSegment>`: same type, same reference, and — *if we
configure it so* — the same settings. Run the walk at the catalog's settings, collect every
`SsrSegment`, compare. What that proves is that the plumbing is sound: the windowed scan, the
pre-filter, the classification call, the coordinate arithmetic. **It proves nothing about whether those
settings are right**, and it is not a licence to treat the catalog's numbers as correct (§5).

At any other configuration the two *should* differ — that is what an experiment is (§5.2). So this
is a **fixed-config regression test**, and it must be pinned to the catalog's settings explicitly
rather than to whatever `Default` happens to be, or it will start failing the first time someone
moves a floor and reads it as a bug rather than a result.

`scanner_parity.rs` is the precedent (compares `SsrSegment` sets, tolerant of boundary wobble via span
overlap, readable diff).

**Which `.cat`, now that the detector swap is cancelled (Revision).** An earlier draft preferred a
`.cat` built through the *scanner* path, to isolate this step from the swap question. There is no
swap: the committed golden fixture (`tests/data/tandem_repeat/golden.ssr_catalog.bed.gz`) is
**trf-mod-built**, and that is the one to use. Two consequences, both fine:

- **It is a genuinely independent yardstick** — a different detector, a different code path, nothing
  ng touched. That is worth more than isolation, and it is precisely the "trf-mod for comparison"
  use the Revision sanctioned.
- **The detector difference is already characterised**, so it is not a confound: `scanner_parity`
  measured scanner-vs-trf-mod at **16/16 recall** — 15 exact, one ±1–2 bp boundary/phase wobble, one
  genuine scanner-only locus trf-mod's significance model rejected. Read the `.cat` through
  production's `CatalogReader` (read-only; no production change), and inherit `scanner_parity`'s
  overlap tolerance rather than demanding byte-exact borders. §8.0 is what pins the *arithmetic*
  exactly; this oracle pins the *pipeline*.

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
a low-copy tract (both `Generic` — the §2.2 property, and see the correction below); a homopolymer; a 2 kb array (`Satellite`); a
real STR *inside* a 2 kb array (swallowed — the §2.4 cost, made visible); **two tracts 10 bp apart**
(one `SsrBundle` carrying both, *not* two `Generic` regions); **three tracts chained 30 bp apart** (one
bundle of three — emergent transitivity); **a locus 60 bp from a bundle** (classified; bundles do not
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

> **Two corrections to this section's fixture list, from writing it (2026-07-17).**
>
> **"An impure tract → `Generic`" is not one rule.** The *scanner* decides an impure tract's fate long
> before classification does, because Ruzzo–Tompa returns **maximal-scoring** segments. Three outcomes: a
> **small** interruption is paid for by the surrounding matches, the tract stays whole and classifies — a
> **locus**; a **large** one splits the segment, and the pure pieces are a **bundle** if they are close
> (two loci if not); only when the pieces fall under the copy floor is it **`Generic`**. The fixture
> that pins §2.2's property is the third case, and it says so.
>
> **Classification's `Purity` gate is unreachable from the walk** — a tract impure enough to fail the 0.80
> floor always contains a purer sub-segment that scores higher, so the scanner emits *that*. Measured:
> a 0.79-purity fixture comes back a **locus**, its pure core. Of classification's five gates the walk
> reaches exactly two, `Compound` (the homopolymer gate — §5b drops period 1 by the *floor*, so the
> run's period-2 multiple has no eliminator left and arrives as motif `AA`) and `FlankClamped`. See
> `segment_criteria::the_walk_reaches_only_one_of_classifications_five_gates`; §3.1's other three columns are
> structurally zero.
>
> > **Superseded 2026-07-17 — it is ONE gate, `FlankClamped`.** "The run's period-2 multiple has no
> > eliminator left" *was the bug* (§5.2's correction), and `Compound` had no other customer. With
> > `prefilter` reordered, a homopolymer never reaches classification, so nothing is turned down and nothing
> > is counted — the walk-level fixture now asserts `RejectionCounts::default()`. **The zeroes were
> > never a property of classification**: each is caused by a stage *upstream* of the gate — the scanner's
> > maximal-scoring segments for `Purity` and `NoCleanTrim`, the pre-filter for `CopyFloor` and now
> > `Compound` — so the count moves whenever that stage does, and it has now moved twice. **Four of
> > §3.1's five rejection columns are structurally zero.** The test is renamed
> > `segment_criteria::the_walk_reaches_only_one_of_classifications_five_gates`.

---

## 9. Reuse, deferred work, and what this breaks elsewhere

| what | existing code | how |
|---|---|---|
| the windowed scan (core + margin, coverage clipped, intervals by start) | `collect_windowed` ([tandem_repeat.rs:530](../../../../src/ng/tandem_repeat.rs)) | **the primitive this is built on** — must be promoted and streamed (§6.1) |
| repeat detection | `find_tandem_repeats` ([tandem_repeat.rs:354](../../../../src/ng/tandem_repeat.rs)) | called per window by the above |
| the region tiling | `RegionScanner` ([tandem_repeat.rs:586](../../../../src/ng/tandem_repeat.rs)) | **not reused** — merges before policy (§6.1) |
| the classification policy | `build_loci` ([postprocess.rs:69](../../../../src/ssr/catalog/postprocess.rs)) | **copied into ng** — every decision transcribed unchanged; windowed, `RepeatInterval`-taking, 1-based/`u64`, all-knobs, returns `Classified` (§5a), and **without the flank embed** (§1.2). Production's stays as it is |
| the bundle selection | `drop_bundles`/`is_close` ([postprocess.rs:274](../../../../src/ssr/catalog/postprocess.rs)) | **copied with it** — selection kept, disposal not; re-associated to stream (§2.4, §2.6) |
| the pre-filter | `catalog_prefilter` ([ng/scanner_parity.rs](../../../../src/ng/scanner_parity.rs)) | **copied into ng**, beside ng's `classify` (§5.1); the parity test keeps its own frozen copy, and B2 moved that test out of `src/ssr/` too (§9) |
| classification knobs | `CatalogParams` ([catalog/mod.rs:42](../../../../src/ssr/catalog/mod.rs)) | **ng's own `SsrSegmentCriteria`** — same defaults, plus the hardcoded period scope + copy floors as real knobs; no separate `bundle_threshold` (§2.4) |
| the STR locus | `ssr::types::Locus` ([types.rs:136](../../../../src/ssr/types.rs)) | **copied into ng**, born 1-based/`u64` (§4) and **without `ref_bytes`** (§1.2). Production's keeps them, 0-based/`u32` |
| the motif | `ssr::types::Motif` ([types.rs:36](../../../../src/ssr/types.rs)) | **copied into ng** — reuse looked free (no coordinates, no width) but production's is `pub(crate)` and would leak through ng's `pub` `SsrSegment` (§4, corrected at the A1 review) |
| what to walk | `RegionSet`/`Region`/`ContigBounds` ([regions.rs](../../../../src/regions.rs)) | **wrapped, read-only**, by `GenomeRegions` (§2.5); it already parses/coalesces/clamps BED. `regions.rs` does not move |
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

**What this step changes elsewhere: one test move, and nothing else (Revision 2026-07-16).**

> **The one exception, owner-approved 2026-07-16 at B2.** `src/ssr/catalog/scanner_parity.rs` moved to
> [`src/ng/scanner_parity.rs`](../../../../src/ng/scanner_parity.rs). It was **the only place in the
> whole tree where production `use`d `crate::ng`** — and that dependency pointed the wrong way: the
> file is *ng's* test (its subject is ng's scanner; production is only the yardstick), living in
> `src/ssr/` for the historical reason that the golden path was there.
>
> **It bit at B2.** Widening ng's `RepeatInterval` to `u64` (§4) could not compile without editing
> frozen production, because that file baked in its `u32`-ness — four sites, most bindingly
> `TrfRecord::for_test(start: u32, end: u32, …)`. **A freeze the experiment can break is not a
> freeze:** the entire point of leaving production alone is that ng cannot destabilise the yardstick
> it is scored against, and this made every ng type change a potential `src/ssr/` edit.
>
> So the move pays a one-time cost — delete the file, drop one `mod` declaration — to remove the
> coupling **permanently**. It is test-only: nothing shipping changed, the parity bar is the same, and
> the result is unchanged (16/16 recall, 15 exact, 1 boundary/phase wobble, 1 scanner-only). After it,
> **production depends on nothing in ng**, and this plan's "no file outside `src/ng/`" rule holds for
> everything that follows.
>
> Rejected: casting at the four sites (production keeps importing ng, so the tax recurs at every ng
> type change); and not widening at all (leaves ng mixed-width — the state §4's decision exists to
> end).

The earlier draft said the opposite — *"where the rest of the code does not match a decision here, we change the
code"* — and listed a rebase of `ssr::types::Locus` across 59 call sites, three changes to
`CatalogParams`, and a `.cat`-side conversion, all as *"work this step includes"*. **All of it is
cancelled.** Production is frozen; ng copies what it needs (§5). The list below is what remains, and
it is entirely ng's own code.

**The `SsrSegment` rebase is not performed — it is not needed.** This was the draft's riskiest single item
("an off-by-one here is a silently wrong genotype, not a crash"), scoped at 59 call sites across
`src/ssr/`, and it earns a note because *deleting it is the Revision's clearest dividend*: ng's
`SsrSegment` is born 1-based, so there is no rebase, no 59 sites, no silent-genotype risk, and no need for
the isolated-commit ritual the draft designed around it. Production's `Locus`, its `.cat` format, and
`fetch_locus_reads`' inline conversion all stay exactly as they are. The one place the arithmetic
still has to be right is ng's port — and §8.0 pins it against production directly, which is a
*stronger* check than the fixture-green-before-and-after the rebase would have had.

- **`ssr_repeat_scanner.md`**, three ways: its **§6–§7 production detector swap is cancelled** (the
  Revision) — trf-mod keeps building the catalog, `TrfRecord` and `trf.rs` stay, and ng simply never
  depends on any of it. The pre-filter's home is **ng**, not `catalog::run` (§5.1). And
  `collect_windowed` must be promoted and streamed while `RegionScanner` goes unused (§6.1). Plus its
  `u32` types widen (§4) — those are `src/ng/`, so they are ours.
- **`ssr_catalog.md`** — **no changes.** The earlier draft's three (`u64` widening, dropping
  `bundle_threshold`, hoisting `MIN_PERIOD`/`MAX_PERIOD`/`copy_number_floor` into knobs) all now land
  in ng's own `SsrSegmentCriteria` instead (§5). The catalog keeps its hardcoded rules, its `u32`, and
  its two copy-floor tables; ng's copy folds them into one **on ng's side only**. The catalog
  documents production, and production did not move.

  *(Corrected 2026-07-16, A1 review: earlier drafts of this spec called the two tables
  **disagreeing**, and A1's code inherited the word. They do not disagree — they give the same floor
  for every period either can reach. Their one numeric difference is period 1 (10 vs 3), and period 1
  reaches neither: the pre-filter gates `period >= 2` and `MIN_PERIOD` drops it again. The duplication
  is **structural, not behavioural**, so A2 deletes a copy rather than resolving a conflict — a
  smaller and better-understood job than this spec had been describing. Pinned by
  `the_two_copy_floor_tables_agree_on_every_reachable_period`.)*
- **`ref_seq.md`** — the parked raw-from-windowed YAGNI has fired (§6); `u32` → `u64` (§4). `src/ng/`,
  ours to change.
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

  *Updated 2026-07-17: **in ng this trap is gone**, and the experiment is now a one-knob change.*
  *A2 folded the two copy-floor tables into one `MinCopies`, and the pre-filter reorder (§5.2's*
  *correction) made the remaining gates consistent: the period floor is `periods.min()`, applied*
  *last, and `MinCopies::for_period(1)` = 10 is applied first. So `periods = 1..=6` is the whole of*
  *"put period 1 on the STR path" — no second clause to relax, no annihilation hazard to re-arm,*
  *because the copy floor keeps the specks out. **`copy_number_floor(1)` = 10 is no longer dead:***
  *it is what the reorder brought to life.* The trap as stated still describes **production**,
  whose `postprocess.rs` keeps both tables and the up-front gate.

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
