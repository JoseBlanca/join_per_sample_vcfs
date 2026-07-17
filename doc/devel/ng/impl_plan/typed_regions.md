# ng typed-region generator (step 3) — implementation plan

**Status:** draft, 2026-07-16 (**revised same day** — see below). Build order for **step 3, the
typed-region generator**: the `region_typing/` module and the walk that cuts the reference into
`TypedRegion`s, plus the ng-side substrate it stands on (the ported admission policy, a streamed
`collect_windowed`, ng's `u64` widening). Design is settled in
[`../spec/typed_regions.md`](../spec/typed_regions.md) (spec) and
[`../arch/typed_regions.md`](../arch/typed_regions.md) (types & interfaces). This turns that design
into build order; it is **not** a place for new design — every open item is a *parameter value*
resolved in spec §10, not an interface question.

> **Revision, 2026-07-16 — production is frozen; ng copies what it needs.** The first draft of this
> plan opened with a **Milestone A of production-side work in `src/ssr/`**: rebase `Locus` to
> 1-based/`u64` across ~59 sites, rewrite `CatalogParams`, change `build_loci`'s signature — each
> verified against the golden catalog. The owner reversed that (spec Revision): *"don't touch
> `src/ssr`… copy the code to `ng/` and modify there… leave production as is and create a fresh ng
> caller from scratch."* ng must also not depend on trf-mod, so the
> `ssr_repeat_scanner.md` §6–§7 detector swap this plan was **blocked on is cancelled**, and with it
> the blocker.
>
> **What that does to the build order.** Milestone A stops being a production rebase and becomes
> **ng's port of the admission policy** (spec §5) — same three concerns (locus type, knobs, windowed
> admission), same position in the order, but landing in `src/ng/region_typing/admission.rs` and
> verified against production *from the outside* rather than by editing it. Milestones B–E are
> **unchanged**: they were always ng's own code.
>
> **The plan gets safer, not just different.** A1's whole ceremony — its own commit, never bundled,
> golden green before *and* after, `git bisect`-able — existed because rebasing production's `Locus`
> risked *a silently wrong genotype in the shipping caller*. ng's `Locus` is born 1-based, so that
> risk does not exist: nothing shipping changes. The port's arithmetic still has to be right, and
> spec §8.0's differential test pins it against production's `build_loci` directly — a **stronger**
> check than fixture-green-before-and-after, because it compares the two implementations rather than
> one implementation to itself.

This plan touches **no file outside `src/ng/`.** Production (`src/ssr/`, `src/regions.rs`,
`Containerfile`) is read-only throughout; the only production items ng names are `RegionSet` (wrapped
read-only) and — in test code only — `build_loci` + `TrfRecord::for_test` + `CatalogReader`, as
oracles. (`ssr::types::Motif` was to be reused as-is; A1 found it is `pub(crate)`, so a `pub` ng type
cannot return it, and ported it instead — spec §4.)

---

## Scope

**In:** ng's ported admission policy (its own `Locus` 1-based/`u64`, `SsrAdmissionParams` with the
period scope + copy floors as real knobs, the pre-filter, and a windowed `admit` returning
`Admitted { loci, bundled }`); a promoted+streamed `collect_windowed`; ng's `u64` coordinate
widening + `WindowedRefSeq` raw bytes; `ng::types` `GenomeRegion`/`Position`; `region_typing/` with
`TypedRegion`/`RegionKind`/config/counts/error; `GenomeRegions`; the resident partition, bundle
detection, and the windowed walk; `TypedRegionIterator`; BED scan-wider-than-emit.

**Out (Revision — no longer this plan's work at all):** the `Locus` rebase in `src/ssr/`; any change
to `CatalogParams`, `build_loci`, `catalog/io.rs`, or `regions.rs`; the trf-mod → scanner detector
swap (`ssr_repeat_scanner.md` §6–§7), including deleting `trf.rs`/`TrfRecord` and the `Containerfile`
step. **Production does not move.**

**Out (later plans / other work):**

- **The generic join** — the pileup splitting a `Generic` region into loci from the data. The next
  slice; `pileup/`'s own spec.
- **The STR gatherer** — fetch → prepare → tally an `SsrLocus` into `LocusEvidence`; a region-driven
  `RecordSource`. Its own spec (spec §9).
- **Bundle disposal** — what a consumer does with an `SsrBundle` (spec §10). The gatherer's call; this
  plan only *emits* the region.
- **The experiments** — period × length routing frontier, flank size, satellite cap (spec §10). Config
  sweeps on the finished walk, not build steps.
- **Cleanup of `ng_step_interfaces.md` §3's stale `SsrLocus` struct / the spec's "`LocusKind` retires"
  line** — noted in the step-3 arch commit; a docs follow-up, not code.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The substrate before the walk.** Everything the walk *calls* — ng's ported `Locus` + `admit`, the
  streamed `collect_windowed`, the raw `WindowedRefSeq` — is built and verified against its own
  oracle (Milestones A–B) before a line of the ng walk exists.
- **Simplest impl first, as the next one's oracle.** The **resident** partition (D1, whole contig in
  memory) is built and proven first; the **windowed** walk (D3) is then proven *by matching it* —
  which is exactly the window-invariance the spec demands (§2.3, §2.6).
- **The algorithmic heart before the plumbing.** The partition + bundle detection (D) is decided and
  tested before the iterator, ownership, and BED narrowing (E) wire it to the surface.
- **Transcribe, then change — never both at once.** The port (A) is the plan's one silent-failure
  risk: an off-by-one in `admit` is a *wrong locus, not a crash*, and it is invisible without an
  oracle. So the port lands as a **faithful transcription first**, pinned by spec §8.0's differential
  against production's `build_loci`, and only then gains the windowing and the new knobs — each with
  the differential still green at the degenerate whole-contig case. Never transcribe and redesign in
  the same commit.
- **Verify against ground truth, from outside.** Production is frozen, so it is a real yardstick:
  **spec §8.0's differential** (ng's port vs `build_loci`, exact modulo the coordinate base) for the
  port; **`.cat` parity** against the committed trf-mod-built golden catalog (subset by the satellite
  cap) plus the partition and invariance tests for the walk. Not self-consistency.
- **Production is read-only.** No step in this plan edits a file outside `src/ng/`. If one appears to
  need to, that is a stop-and-ask, not a small deviation (spec Revision).
- **Incremental, with pauses.** One milestone, stop for review, next.
- **Container builds.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); native host build at completion.

## Preconditions (all verified in place, 2026-07-16)

- **The tandem-repeat scanner is built** — `find_tandem_repeats`, `collect_windowed` (private today),
  `RegionScanner` ([ng/tandem_repeat.rs](../../../../src/ng/tandem_repeat.rs)).
- **The port's sources are readable and the oracles are live**: `build_loci`
  ([postprocess.rs:69](../../../../src/ssr/catalog/postprocess.rs)), `catalog_prefilter`
  ([scanner_parity.rs](../../../../src/ng/scanner_parity.rs) — moved into ng at B2, see below), `Locus` + `Motif`
  ([ssr/types.rs](../../../../src/ssr/types.rs), both `pub(crate)`), `CatalogParams`, `RegionSet`
  ([regions.rs](../../../../src/regions.rs)), the **golden catalog fixture**
  (`tests/data/tandem_repeat/golden.ssr_catalog.bed.gz`), and the **`RefSeq` impls**.
- **`TrfRecord::for_test` is `#[cfg(test)] pub(crate)`** — same crate, so ng's own test modules can
  drive production's `build_loci` as spec §8.0's differential oracle **without touching it**. This is
  what lets the Revision freeze production and still pin the port.
- **No swap precondition any more.** The first draft required the trf-mod → scanner catalog swap to
  have landed (it has not, and it is now cancelled — Revision). ng's `admit` takes `RepeatInterval`
  because ng wrote it that way, so spec §5c's "`build_loci`'s door is not open" does not apply.

---

## The steps

### Milestone A — ng's admission port (`src/ng/region_typing/admission.rs`; each vs production, from outside)

*Revised 2026-07-16: was "catalog coordinates & knobs (production)". Same three concerns — locus type,
knobs, windowed admission — now landing in ng instead of `src/ssr/`. Production is read-only here;
it appears **only** in `#[cfg(test)]` code, as the oracle.*

**A1. ng's `Locus` + the faithful transcription of `build_loci`, with the differential that pins it.**  ✅
Scaffold `src/ng/region_typing/` (`mod.rs` + `admission.rs`, wired into `ng/mod.rs`). ng's `Locus`:
production's fields, **1-based inclusive and `u64`**, private fields + validating `new`. (`Motif` was
to be reused as-is; it is `pub(crate)`, so a `pub` ng `Locus` cannot return it — ported too, spec §4.)
Then transcribe `build_loci` → `admit` **at its current shape** —
whole-contig, `Vec<RepeatInterval>` in place of `Vec<TrfRecord>` (it only ever reads
`start/end/period/score`), `SsrAdmissionParams` starting as a straight copy of `CatalogParams`'s
fields — carrying `drop_bundles`/`is_close`/`minimal_trim`/purity-recompute/`ref_bytes`-embed
**unchanged in logic**, and the pre-filter (`catalog_prefilter`) beside it.

**This is the plan's one silent-failure step** — an off-by-one yields a wrong locus, not a crash — so
it lands as its own commit and its oracle is built *in the same commit*: **spec §8.0's differential**,
ng's `admit` vs production's `build_loci` (bridged via `TrfRecord::for_test`) on the scanner's
intervals over the synthetic fixture **plus** the crafted cases from `postprocess.rs`'s own tests,
asserting identical `Locus` sets **modulo the coordinate base**, through one stated conversion helper.
No behaviour change is attempted here: the only intended difference from production is the coordinate
base and the width. *Source:* spec §4, §5, §5.1, §8.0; arch §types, §recon.

**A2. `SsrAdmissionParams` — collapse the redundant knob, expose the hidden ones.**  ✅
Now that the transcription is pinned, change it — **on ng's side only**. Fold `bundle_threshold` into
`flank_bp` (spec §2.4 — one number). Hoist `MIN_PERIOD` / `MAX_PERIOD` / `copy_number_floor` from
hardcoded consts into fields (`Default` = today's values), and reconcile the pre-filter's second,
disagreeing copy-floor table with them — the reconciliation the catalog never got to make (spec §10).
State the inherited `min_score: 0` trap explicitly (spec §5c: there is no trf-mod `-s 30` in ng, so
`Default` ships no score gate — known, not accidental). *Depends:* A1. *Source:* spec §2.4, §5, §5c.
Verified: **A1's differential still green at `Default`** (that is what "same defaults" means, and it
is the whole reason A1 came first); a non-default period scope / copy floor demonstrably changes the
admitted set.

**A3. `admit` windowed, and returns what it set aside.**  ✅
Add `bases_start` + `contig_len` params — the two facts a bare slice silently stood in for (spec §2.6);
whole-contig (`bases_start = 1`, `contig_len = bases.len()`) is the degenerate case → unchanged output.
Return `Admitted { loci, bundled }`, where `bundled` is the cluster members `build_loci` silently drops
today (spec §6a). *Depends:* A1, A2. *Source:* spec §5a, §6a. Verified: **A1's differential still green
at the degenerate case** (this is what makes windowing provably behaviour-preserving); `bundled`
non-empty on a clustered fixture; a locus near a slice edge is **not** dropped when `contig_len` says
the contig continues — the §2.6 bug, tested directly rather than argued.

> **Checkpoint A — REACHED 2026-07-16.** ng owns a 1-based/`u64` admission policy whose every rule is a
> knob, windowable and bundle-reporting, **pinned to production's `build_loci` by a green differential**
> — and `src/ssr/` has not been touched. Pause for review.
>
> **Routed forward from A3's review, so they are not lost:**
> - **C1** — a `Window { bases, bases_start, contig_len }` newtype (or `Position`/`Bp` doing the same
>   job). `admit`'s three window params are co-dependent and two are same-typed `u64`s that transpose
>   silently; more importantly, the `contig_len` guard is **inherently incomplete** — a caller passing
>   the window's own end as `contig_len` is arithmetically legal, and only *provenance* (a contig table)
>   can catch it. C1's newtypes are already scheduled and need no spec amendment.
> - **D1** — assert the partition on the scanner-output fixture; today it projects `.loci` and throws
>   `bundled` away, so the only realistic-density fixture asserts nothing about the whole partition.
> - **Own step** — a `proptest` property over `assert_agrees`. `proptest` is already a dev-dependency
>   and `assert_agrees` is a ready-made property body; it is the strongest remaining lever on the port.

### Milestone B — the windowed, wide, raw substrate the walk stands on

**B1. Promote and stream `collect_windowed`.**  ✅
`private → pub(crate)`; return an **iterator** of per-window `{ core, coverage, intervals }` instead of
two whole-contig `Vec`s (spec §6.1). `RegionScanner` (which consumed it) stays as-is over the streamed
form. *Source:* spec §6.1. Verified: the existing **window-count-invariance** test still green; streamed
output == the old collected output.

**B2. Widen ng's own coordinates to `u64`.**  ✅
*(Blocked on a conflict the plan did not foresee: `src/ssr/catalog/scanner_parity.rs` imported
`crate::ng` and baked in `RepeatInterval`'s `u32`, so widening ng meant editing frozen production —
i.e. **production depended on ng**, the exact coupling the freeze exists to prevent. Owner-approved
fix: move that test into `src/ng/` (commit `d097ebf`). `rg 'use crate::ng' src/ --glob '!src/ng/**'`
now matches nothing. Also decided here, owner-agreed: **`RepeatInterval` stays 0-based half-open**
against spec §4's "1-based everywhere" — it is a detector/slice type whose coordinates index the
byte slice, so 1-basing it would add a `- 1` at every site to buy consistency with a genetic
coordinate it never expresses. The 1-based boundary belongs at `Locus`/`GenomeRegion`.)*
`RefSeq::fetch_into` + its three impls; the scanner's `RepeatInterval` / `RegionSpan` / `SegmentOptions`;
`Bp`. **Delete `ref_seq.rs:87`'s `unwrap_or(u32::MAX)` silent clamp** (a >4 Gb contig now represents,
not clamps). *Source:* spec §4, §9. Verified: the ng test suite compiles and passes.

**B3. `WindowedRefSeq` — raw bytes + contig table.**  ✅
*(Under-sized by this plan. "Add a `RawRefSeq` impl (un-canonicalised buffer)" assumed ng could reach
the raw bytes; production's `ManualEvictChromRefFetcher` **canonicalises into its buffer**, so they
are gone before ng sees them — and production has no windowed raw reader at all, with its `.fai`
machinery private. Owner chose to copy it: `src/ng/raw_chrom_reader.rs`, production's design with one
line different (`dst.push(canonicalise(b))` → `dst.push(b)`). This also settles spec §10's open
question — **raw is the buffer, canonical the derived view**, because two readers is what §6's 14.6 GB
lesson forbids and `canonicalise` only runs one way.)*
Add a `RawRefSeq` impl (un-canonicalised buffer — `ref_seq.md`'s parked YAGNI, now spent, spec §6) and a
`contigs() -> &ContigList` accessor. *Depends:* B2. *Source:* spec §6, §8.3. Verified: raw fetch returns
verbatim bytes matching `ResidentRefSeq`'s raw path; `contigs()` returns the table.

> **Checkpoint B:** the scan streams per window; ng speaks `u64` end to end (no clamp); `WindowedRefSeq`
> serves raw bytes and its contig table. Pause for review.

### Milestone C — ng vocabulary + local types (types, no logic)

**C1. `GenomeRegion` + `Position` in `ng::types`.**  ✅
*(Also lands A3's routed follow-up: `admit`'s two adjacent, same-typed, silently-transposable `u64`s
become `Position` + `Bp`, so the compiler catches always what the runtime guard caught sometimes —
the reviewer's own preferred option, since it needs no spec amendment. The **provenance** half is
not a signature problem and is not fixed here: nothing `admit` is handed proves where `contig_len`
came from, so it is a contract stated at `admit` and discharged by the **walk** (D), which is the
only party holding both the reference table (`contigs()`, B3a) and the window.)*
1-based inclusive, `u64` (`Position(u64)`; `GenomeRegion { contig, start, end }`). The consolidation
`ng_step_interfaces.md` §6 reserved — this step's first real use. *Source:* arch §types.

**C2. `region_typing/mod.rs` walk types.**  ✅
The module was scaffolded in A1; this adds the walk's own types. `TypedRegion`, `RegionKind`
(`SsrLocus(Locus)` — ng's, from `admission.rs` / `SsrBundle` / `Generic` / `Satellite`),
`TypedRegionConfig` (+`Default` = the catalog's settings, for comparability, holding A2's
`SsrAdmissionParams`), `TypedRegionCounts`, `TypedRegionError`. No logic. *Depends:*
A, B, C1. *Source:* arch §types. Verified: `Default` config test.

**C3. `GenomeRegions` wrapping `RegionSet`.**  ✅
*(Found: the seam **widens only** — `regions::Region` is already 1-based inclusive, so spec §4's
"and rebasing" was wrong and is corrected there.)*
`whole_contigs` (default), `from_bed_path`, `iter()` widening `RegionSet`'s `u32` → `u64`. Reimplements
nothing — `RegionSet` already parses/coalesces/clamps/converts-BED/drops-empty-contigs. *Depends:* C1.
*Source:* spec §2.5, arch §types. Verified: `whole_contigs` over a contig table; a BED round-trip.

> **Checkpoint C:** ng types compile; `GenomeRegions` wraps `RegionSet` and widens. Pause for review.

### Milestone D — the walk (the heart); resident first as the windowed oracle

**D1. The resident partition — the five steps, whole contig in memory.**  ✅
*(`SsrBundle` is D2's; until then the cluster members admission sets aside fall into `Generic`, which
is where a repeat rejected for any reason belongs (spec §2.2) — so the partition is already complete,
and D2 only makes it sharper. **One known-untested claim recorded in the code**: the cap applies to
the *cleaned* coverage, not the raw; mutation shows no fixture discriminates the two, because they
differ only where raw-only coverage runs ≥ `max_repeat_len` contiguously.)*
`detect → clean → admit → cap → partition` over a contig held resident, producing `Vec<TypedRegion>`
(no windowing). Uses `find_tandem_repeats` + ng's pre-filter + ng's `admit` (whole-contig form) +
the coverage merge + satellite cap — all of Milestone A. *Depends:* A, B, C. *Source:* spec §2.1, §2.2,
§8. Verified: the **partition invariant** (contiguous / non-overlapping / complete / maximal) **and
`.cat` parity** against the committed trf-mod-built golden catalog (strict subset — every catalog locus
present, or absent *and* inside a satellite run; `scanner_parity`'s overlap tolerance for detector
wobble, spec §8.1) on a small reference. **This is D3's oracle.**

**D2. Bundle detection — the flank test, and the rejected-repeat split.**  ✅
A repeat with another repeat within `flank_bp` on either side is a bundle member; the cluster (hull of
its tracts) is one `SsrBundle`; a repeat admission turns down for any *other* reason is `Generic`, not a
hole (spec §2.2). Consumes `admit`'s `Admitted.bundled` (A3). *Depends:* D1. *Source:* spec §2.2,
§2.4. Verified: two tracts 10 bp apart → one `SsrBundle` carrying both; three chained 30 bp apart → one
bundle of three; a locus 60 bp from a bundle → admitted; an impure/low-copy repeat → `Generic`.

**D3. The windowed walk.**  ✅
*(Landed 2026-07-16. Three findings worth keeping. **(1) The spec's windowing recipe was wrong in
the same paragraph, twice** — §2.6's "keep only the repeats whose start lands in the core" reads as
*filter, then admit*, which is window-**variant**: admission's rules are neighbour rules (the bundle
flank test, redundancy elimination), and a core tract's neighbour lives in the margin. Attribution
is the **last** step, not the first. And "1 kb of margin covers admission for free" is false —
admission reads `tract ± flank_bp` out of the bases it was handed and **panics** rather than quietly
drop a locus, and the tracts at the slice's own edge have no margin of their own. Fix: the **scan**
gets core ± `max_repeat_len`, the **bases handed to admission** are that slice grown by a further
`flank_bp`. Spec §2.6 corrected. **(2) `max_repeat_len >= flank_bp` is a real constraint**, now
asserted (release — §10 sweeps both). **(3) D1's known-untested claim is now tested**: the cap
applies to the *cleaned* coverage, and D3's 6 kb fixture is the first to discriminate it. Also:
`ScannedWindow` was reshaped — the fetched span plus **every** detection in it, with the two window
rules as methods over a caller-chosen interval set, so `RegionScanner` applies them to raw intervals
and the walk to pre-filtered ones, from one copy of each rule (§6.1). `ContigTable` (a trait on the
reference) is how the walk reads `contig_len`, which discharges A3's provenance obligation.)*

> **B1 and B3 contradict each other, and neither step could see it alone.**
>
> - **B1** promoted `scan_windowed` as spec §6.1 asked, keeping its source: it reads through
>   production's `ChromRefFetcher`, which **canonicalises into its buffer**.
> - **B3** established that the walk needs **raw** bytes: the STR catalog reads the FASTA verbatim
>   and `Locus` compares **by value**, so canonical bytes make every IUPAC-carrying locus compare
>   unequal to the catalog's — silently breaking the `.cat` oracle on any assembly carrying them
>   (spec §6). That is why `RawChromReader` and `WindowedRefSeq: RawRefSeq` exist.
>
> **So D3 cannot both reuse `scan_windowed` and read raw.** `WindowedRefSeq` is a `RawRefSeq`, not a
> `ChromRefFetcher`; `scan_windowed` cannot be handed it. And spec §6.1 is explicit that
> reimplementing the windowing is the thing *not* to do — *"duplicates a subtle invariant that can
> then drift"* — while spec §6 forbids holding two readers (*"one reader for the whole run … that is
> what cost 14.6 GB of peak RSS"*).
>
> **RESOLVED (owner, 2026-07-16): option 1, landed.** `scan_windowed` now takes `contig_len` and an
> `FnMut(start, len, &mut Vec<u8>) -> Result<(), E>`; both sources pass their own bytes through one
> copy of the window rules, and the error type is the caller's. All 31 B1 tests pass unchanged.
>
> **Options considered:**
> 1. **Generalise `scan_windowed` over its byte source** ✅ — a trait or a `FnMut(start, len) ->
>    Result<Vec<u8>, _>`, so ng's raw reader and production's fetcher both fit. Keeps §6.1's reuse,
>    keeps one reader, and is a small change to B1's own code. **Recommended.**
> 2. The walk windows for itself — duplicates exactly the invariant §6.1 says not to duplicate.
> 3. Two sources (canonical for detection, raw for admission) — what §6's 14.6 GB lesson forbids.
>
> *Detection itself is width- and case-agnostic (`find_tandem_repeats` is case-insensitive and
> non-ACGT never matches), so only **admission** needs raw — which is why nothing before D3 tripped
> on this.*
The three carries (open cluster, open coverage run, open generic run); the `max_repeat_len` detection
margin; the contig length **passed in** to `admit`, never inferred from the slice (spec §2.6 — the
silent window-boundary bug). Built on B1's streamed `collect_windowed`. *Depends:* D1, D2. *Source:*
spec §2.3, §2.6. Verified: **window-invariance** — output identical to D1 resident at two `window_bp`,
with features placed **astride** window boundaries (small `window_bp`, not a small fixture); the
**maximality** case (a generic run longer than a window emits as **one** region); the contig-end-clamp
guard (a locus near a window edge is not dropped).

> **Checkpoint D — REACHED 2026-07-16.** The partition is correct resident (parity + invariant),
> bundles route, and the windowed walk matches the resident oracle across window sizes and on the
> golden reference. Pause for review.
>
> **Found at D3, decided by the owner the same day, and landed — spec §2.4a:**
> - **A satellite swallows what it touches.** The walk had no answer for a microsatellite beside a
>   satellite and was giving two different wrong ones depending on which *side* it sat: on the left,
>   an `SsrBundle` **overlapping** the `Satellite` (an invalid partition); on the right, the cluster
>   dropped whole and the bases silently `Generic`. Probed rather than argued, which is what showed
>   the asymmetry. **Owner's rule:** a microsatellite *or* a bundle within `flank_bp` of a satellite
>   is swallowed, and the satellite **expands** to cover it. Both kinds, because both arise by
>   different routes (the array's tract may or may not clear admission's gates, but the satellite is
>   built from *coverage*, which does not care). Absorption subsumes §2.1's swallow — containment and
>   adjacency are one predicate — and applies to bundle **members** before clustering, which is also
>   what keeps the windowed walk correct (a truncated over-cap member would break the cluster
>   re-derivation; `bundle_clusters`' `debug_assert` caught exactly that, as D2 built it to).
> - **The stakes of the cleaned-vs-raw cap are still undemonstrated.** D3 pinned the *choice*
>   (`the_satellite_cap_applies_to_the_cleaned_coverage_not_the_raw`), but only via a boundary
>   difference. A fixture with ≥ 1 kb of contiguous low-copy noise would show the failure the spec's
>   argument is actually about — noise inventing a satellite and swallowing the loci beneath it.
> - **Still open from Checkpoint A:** a `proptest` property over `assert_agrees` (the strongest
>   remaining lever on the port); D1's scanner-output fixture asserts nothing about `bundled`.

### Milestone E — public surface, BED, integration

**E1. `TypedRegionIterator`.**  ✅
Owns `reference` + `spans` (by value — moves onto a producer thread, spec §7); `over_regions(...)`;
`Iterator` with `Item = Result<TypedRegion, TypedRegionError>`, fused, fatal-in-stream; `counts()`
running tally. Wires D3 to the surface. *Depends:* D3. *Source:* arch §interface, spec §8.2.
Verified: iterate a small reference; `counts` tally; a reference-read failure surfaces as
`Some(Err(_))` once then `None`.

*(Landed 2026-07-16. **`partition_windowed` is now this iterator collected**, so D3's whole suite —
window-invariance, `.cat` parity through the oracle, absorption — runs through the shipping code
rather than beside it. Four departures from the arch doc, each forced and each recorded:*
1. ***The constructor is fallible.*** The arch says infallible because `GenomeRegions` validated the
   spans. It validated them against *a* contig table; nothing ties that one to the reference's, so a
   span naming a contig this reference lacks is still reachable — caught up front, once.
2. ***It is generic over the reference, not concrete.*** The arch's reason for concrete was that the
   walk needs raw bytes and eviction, "impl capabilities not trait methods". They are trait methods
   now (`RawRefSeq` + `ContigTable` + the new `EvictableRefSeq`), which is what lets **the same walk**
   run over an `InMemoryRefSeq` in tests — and that is what makes window-invariance against
   `partition_resident` testable at all.
3. ***`WindowCursor` + `scan_window` were extracted from `scan_windowed`*** (no behaviour change, B1's
   tests green). The iterator **owns** its reference, so it cannot hold `scan_windowed`'s iterator,
   which borrows one — that is a self-referential struct. A cursor both can drive keeps §6.1's one
   copy of the geometry.
4. ***`rejected_by_reason` was absent, not zero*** — `admit` reported no reasons, so spec §3.1's
   breakdown could not be filled, and a zero would have read as "nothing rejected for that reason": a
   wrong answer to §10's question dressed as a measured one. **Closed by E1e** (below).
   `repeat_bp_with_no_locus` *is* filled (coverage minus loci) and pinned exactly against an
   independent computation.

**E1e. Admission reports its rejections.**  ✅
`finish_locus` returns `Result<Locus, RejectionReason>` instead of `Option<Locus>`; `admit` returns
`Admitted.rejected: Vec<(RepeatInterval, RejectionReason)>`; `TypedRegionCounts::rejected_by_reason`
is restored and filled. Same gates, same order, same outcome — pinned by A1's differential against
production's `build_loci`, still green. *Source:* spec §3.1.

**Per record, not a tally, and that is the whole design.** `admit` runs over a window's entire fetched
slice, margins included, so every window that can see a repeat rejects it again; a tally taken inside
`admit` would be counted once per window and **the number would move with `window_bp`** — a memory
knob changing a reported number (spec §2.3). Handing the records back lets the walk attribute each to
the core holding its start, exactly as it does for loci and bundle members. Verified by
`the_tally_does_not_depend_on_the_window`; mutation-verified (drop the attribution and it fails — at
`window_bp = 100` the tract at base 1 sits in ten windows' fetched slices).

> **⚠ FINDING — the pre-filter makes three of admission's five gates unreachable from the walk.**
> Found by writing the per-reason walk test and watching it come back **all zeroes**. `prefilter`
> runs first and is not optional (spec §5b):
> - **`CopyFloor`** — `prefilter` applies the **same `MinCopies` table** with the same arithmetic, so
>   nothing under the floor survives to be turned down by it. A2 folded production's two copy-floor
>   tables into one; both *call sites* remain, and this is the visible consequence.
> - **`Compound`** — a compound motif is by definition a period-multiple of a shorter one, which is
>   exactly what redundancy elimination drops.
> - **`NoCleanTrim`** — reachable in principle; a tract with no motif pair in its trim window scores
>   badly enough that the scanner segments it away first.
>
> `Purity` and `FlankClamped` are the live ones (the scanner tolerates ≈0.78, admission demands 0.80,
> so tracts land in the gap; and contig ends are nobody else's business). **Three columns of §3.1's
> breakdown are therefore structurally zero in the walk** — "no compound rejections in this genome" is
> a wrong thing to read from them. Each reason is pinned at admission's own level, where all five are
> reachable, and the *reason* the three are zero is pinned too, so that a pre-filter change fails
> where the explanation is.

*`EvictableRefSeq` is new and load-bearing: the walk releases what it has passed, without which a
`WindowedRefSeq`'s buffer grows to the whole contig and "holds one window" is a claim rather than a
fact. No existing test could catch that — the other impls have nothing to evict — so a recording
reference pins it.)*

**E2. BED scan-wider-than-emit.**  ✅
Grow each requested region by `max_repeat_len` and coalesce (the **scan** set); emit only what overlaps
the **requested** set; clip `Generic` to the requested edge; grow the effective span to hold a straddling
locus **or bundle** whole (spec §2.5). *Depends:* E1. *Source:* spec §2.5, §8. Verified:
**BED-invariance** — loci inside a BED region identical to the whole-genome run, with a cluster placed
**astride the BED edge**; a straddling `SsrLocus` *and* a straddling `SsrBundle` each emitted whole;
`Generic` clipped at the user's edge.

*(Landed 2026-07-17. Two spec corrections, both found by writing the acceptance test.*
1. ***"The walk needs no special logic, because a grown region is just a longer continuous run" is
   wrong.*** A scan span has **edges a whole contig does not**, and a cluster can straddle one: the
   member inside is attributed, its partner outside never is, and `bundle_clusters` re-derives a
   **one-member cluster**. Three of five E2 tests tripped D2's `debug_assert` on it. It needs no
   *policy*, but it does need a **rule at the edge** — a cut cluster is dropped and its bases fall to
   `Generic` (§2.2), which is safe *because* of the margin: a cluster reaching the scan edge cannot
   also reach the emit set, short of the chain-the-whole-way residual §2.5 already names. The check
   moved from `bundle_clusters` to the walk, because only the walk holds the scan span and so only the
   walk can tell an edge artefact (drop) from a real disagreement (panic — how D3's truncated member
   showed up).
2. ***§2.5 never said what `Satellite` does at an edge.*** E2 ruled that it clips, reasoning from
   `RegionKind`'s shape (`Generic` and `Satellite` carry no payload, so clipping could misdescribe
   nothing). **Owner, 2026-07-17: wrong — every *finding* comes back whole; `Generic` alone clips.**
   A `Satellite` means *"an array too long to be a microsatellite"*, so the extent **is** the claim:
   clipping the fixture's 1.2 kb array to a 300 bp request emits a `Satellite` of 300 bp, a span that
   contradicts the `max_repeat_len` test that produced the label. The E2 argument read the type right
   and the meaning wrong — *what a region carries is not what it claims* — and is recorded in §2.5
   because it is tempting and available to make again. `Generic` clips because it is the only kind
   that is not a finding: "nothing more specific can be said" survives clipping intact.

*Mutation-verified: scanning only what is emitted (margin 0) changes what things ARE — the spec's
central claim, now a tested one; clipping objects; clipping territory to the scan span instead of to
each requested span.)*

**E3. Integration anchor.**  ✅
Walk a real small **multi-contig** reference end to end; assert `.cat` parity (subset by cap), the
partition invariant, window-invariance, BED-invariance, and the edge cases (a tract at position 1 →
`Generic`; a repeat-free contig → one `Generic`; a tract at one contig's end abutting one at the next's
start). *Depends:* E1, E2. *Source:* spec §8. The port anchor.

*(Landed 2026-07-17 as `src/ng/region_typing/anchor.rs`. It drives the **shipping stack**: a real
multi-contig FASTA written to disk with its `.fai`, read through `WindowedRefSeq` — file-backed and
evicting — through `TypedRegionIterator`, and nothing else. That is what no unit test reaches:
`RawChromReader`'s windowed reads, the contig table's provenance, eviction, the iterator's ownership,
the scan/emit split and the walk, at once, on sequence nobody wrote to make a point.*

***In-crate `#[cfg(test)]`, not `tests/`***, *and both reasons are real: the golden catalog's
`CatalogReader` is `pub(crate)` and production is frozen, so an out-of-crate test cannot open the
oracle; and this plan must touch **no file outside `src/ng/`** — a `tests/` file would have broken that
on the last step. It is the `scanner_parity` shape exactly (B2, `d097ebf`), and it still consumes only
what a caller could.*

***The control is what makes the edge cases mean anything.*** *Three of them assert "not a locus", and
a tract that was never admissible satisfies all three for free — which is how D1's satellite test once
passed for the wrong reason. A fifth contig holds the **same** `tract(10)`, flanks either side: it comes
back a locus, so the absences are the contig ends doing their job.*

*Mutation-verified: passing the window's end as `contig_len` fails the invariance tests (not parity —
the golden contigs are 2.4 kb, so at the default 100 kb window `core.end == contig_len` and the mistake
is invisible there, which is itself worth knowing); emitting no loci fails parity, the counts, and the
control.)*

> **Checkpoint E — REACHED 2026-07-17.** Step 3 runs end to end; `.cat` parity + all three invariance
> tests green, through the shipping stack. **Step 3 is complete.** Pause for review.
>
> **Open, and needing an owner decision:**
> - ~~`Satellite` clips at a BED edge~~ — **decided (owner, 2026-07-17): it does not.** Every finding
>   comes back whole; `Generic` alone clips. See E2 above.
> - **Three of §3.1's five rejection columns are structurally zero in the walk** (E1e), because the
>   pre-filter gets there first. Not a bug; a thing a reader of the counts must know.
> - **Still open from Checkpoint A:** a `proptest` property over `assert_agrees` — the strongest
>   remaining lever on the port.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | **spec §8.0's differential against production's `build_loci`** — identical `Locus` sets modulo the coordinate base — green at A1 (the transcription), and **still green** after A2's knobs at `Default` and A3's windowing at the degenerate case. The oracle for the plan's one silent step, and it compares two implementations rather than one to itself |
| B | window-count-invariance still green (streamed `collect_windowed`); ng suite green (`u64`, clamp gone); raw-fetch unit + `contigs()` (`WindowedRefSeq`) |
| C | type tests — `Default` config; `GenomeRegions` `whole_contigs` + BED round-trip |
| D | **partition invariant + `.cat` parity** vs the trf-mod-built golden (resident, D1); bundle-routing cases (D2); **window-invariance** — windowed == resident (D3) |
| E | iterator + running counts + fatal-in-stream (E1); **BED-invariance** + straddle-whole (E2); **integration:** `.cat` parity + all invariants on a real multi-contig reference (E3) |

**Production untouched, and checked.** `git diff --stat` over the whole plan must show **no file
outside `src/ng/`** (plus this plan's own docs and reports). That is a mechanical check, run at every
commit — the Revision's guarantee is only worth what it is verified by.

## Out of scope (next plans)

- **The pileup / generic join** — splitting a `Generic` region into loci from the data (`pileup/`'s
  spec); it owns "what counts as a generic locus", the bake-off this step declined to fake.
- **The STR gatherer** — `SsrLocus` region → `LocusEvidence` (fetch + prepare + tally), and the
  region-driven `RecordSource` it needs. Its own spec.
- **Bundle disposal** — the `SsrBundle` consumer (spec §10); the gatherer's call.
- **The routing/threshold experiments** — period × length, flank size, satellite cap (spec §10); config
  sweeps on the finished walk, scored on the standards, not build order.
