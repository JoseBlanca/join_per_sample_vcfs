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

**B1. Promote and stream `collect_windowed`.**  ☐
`private → pub(crate)`; return an **iterator** of per-window `{ core, coverage, intervals }` instead of
two whole-contig `Vec`s (spec §6.1). `RegionScanner` (which consumed it) stays as-is over the streamed
form. *Source:* spec §6.1. Verified: the existing **window-count-invariance** test still green; streamed
output == the old collected output.

**B2. Widen ng's own coordinates to `u64`.**  ☐
`RefSeq::fetch_into` + its three impls; the scanner's `RepeatInterval` / `RegionSpan` / `SegmentOptions`;
`Bp`. **Delete `ref_seq.rs:87`'s `unwrap_or(u32::MAX)` silent clamp** (a >4 Gb contig now represents,
not clamps). *Source:* spec §4, §9. Verified: the ng test suite compiles and passes.

**B3. `WindowedRefSeq` — raw bytes + contig table.**  ☐
Add a `RawRefSeq` impl (un-canonicalised buffer — `ref_seq.md`'s parked YAGNI, now spent, spec §6) and a
`contigs() -> &ContigList` accessor. *Depends:* B2. *Source:* spec §6, §8.3. Verified: raw fetch returns
verbatim bytes matching `ResidentRefSeq`'s raw path; `contigs()` returns the table.

> **Checkpoint B:** the scan streams per window; ng speaks `u64` end to end (no clamp); `WindowedRefSeq`
> serves raw bytes and its contig table. Pause for review.

### Milestone C — ng vocabulary + local types (types, no logic)

**C1. `GenomeRegion` + `Position` in `ng::types`.**  ☐
1-based inclusive, `u64` (`Position(u64)`; `GenomeRegion { contig, start, end }`). The consolidation
`ng_step_interfaces.md` §6 reserved — this step's first real use. *Source:* arch §types.

**C2. `region_typing/mod.rs` walk types.**  ☐
The module was scaffolded in A1; this adds the walk's own types. `TypedRegion`, `RegionKind`
(`SsrLocus(Locus)` — ng's, from `admission.rs` / `SsrBundle` / `Generic` / `Satellite`),
`TypedRegionConfig` (+`Default` = the catalog's settings, for comparability, holding A2's
`SsrAdmissionParams`), `TypedRegionCounts`, `TypedRegionError`. No logic. *Depends:*
A, B, C1. *Source:* arch §types. Verified: `Default` config test.

**C3. `GenomeRegions` wrapping `RegionSet`.**  ☐
`whole_contigs` (default), `from_bed_path`, `iter()` widening `RegionSet`'s `u32` → `u64`. Reimplements
nothing — `RegionSet` already parses/coalesces/clamps/converts-BED/drops-empty-contigs. *Depends:* C1.
*Source:* spec §2.5, arch §types. Verified: `whole_contigs` over a contig table; a BED round-trip.

> **Checkpoint C:** ng types compile; `GenomeRegions` wraps `RegionSet` and widens. Pause for review.

### Milestone D — the walk (the heart); resident first as the windowed oracle

**D1. The resident partition — the five steps, whole contig in memory.**  ☐
`detect → clean → admit → cap → partition` over a contig held resident, producing `Vec<TypedRegion>`
(no windowing). Uses `find_tandem_repeats` + ng's pre-filter + ng's `admit` (whole-contig form) +
the coverage merge + satellite cap — all of Milestone A. *Depends:* A, B, C. *Source:* spec §2.1, §2.2,
§8. Verified: the **partition invariant** (contiguous / non-overlapping / complete / maximal) **and
`.cat` parity** against the committed trf-mod-built golden catalog (strict subset — every catalog locus
present, or absent *and* inside a satellite run; `scanner_parity`'s overlap tolerance for detector
wobble, spec §8.1) on a small reference. **This is D3's oracle.**

**D2. Bundle detection — the flank test, and the rejected-repeat split.**  ☐
A repeat with another repeat within `flank_bp` on either side is a bundle member; the cluster (hull of
its tracts) is one `SsrBundle`; a repeat admission turns down for any *other* reason is `Generic`, not a
hole (spec §2.2). Consumes `admit`'s `Admitted.bundled` (A3). *Depends:* D1. *Source:* spec §2.2,
§2.4. Verified: two tracts 10 bp apart → one `SsrBundle` carrying both; three chained 30 bp apart → one
bundle of three; a locus 60 bp from a bundle → admitted; an impure/low-copy repeat → `Generic`.

**D3. The windowed walk.**  ☐
The three carries (open cluster, open coverage run, open generic run); the `max_repeat_len` detection
margin; the contig length **passed in** to `admit`, never inferred from the slice (spec §2.6 — the
silent window-boundary bug). Built on B1's streamed `collect_windowed`. *Depends:* D1, D2. *Source:*
spec §2.3, §2.6. Verified: **window-invariance** — output identical to D1 resident at two `window_bp`,
with features placed **astride** window boundaries (small `window_bp`, not a small fixture); the
**maximality** case (a generic run longer than a window emits as **one** region); the contig-end-clamp
guard (a locus near a window edge is not dropped).

> **Checkpoint D:** the partition is correct resident (parity + invariant), bundles route, and the
> windowed walk matches the resident oracle across window sizes. Pause for review.

### Milestone E — public surface, BED, integration

**E1. `TypedRegionIterator`.**  ☐
Owns `reference` + `spans` (by value — moves onto a producer thread, spec §7); `over_regions(...)`
(infallible — `GenomeRegions` validated everything); `Iterator` with `Item = Result<TypedRegion,
TypedRegionError>`, fused, fatal-in-stream; `counts()` running tally. Wires D3 to the surface. *Depends:*
D3. *Source:* arch §interface, spec §8.2. Verified: iterate a small reference; `counts` tally; a
reference-read failure surfaces as `Some(Err(_))` once then `None`.

**E2. BED scan-wider-than-emit.**  ☐
Grow each requested region by `max_repeat_len` and coalesce (the **scan** set); emit only what overlaps
the **requested** set; clip `Generic` to the requested edge; grow the effective span to hold a straddling
locus **or bundle** whole (spec §2.5). *Depends:* E1. *Source:* spec §2.5, §8. Verified:
**BED-invariance** — loci inside a BED region identical to the whole-genome run, with a cluster placed
**astride the BED edge**; a straddling `SsrLocus` *and* a straddling `SsrBundle` each emitted whole;
`Generic` clipped at the user's edge.

**E3. Integration anchor.**  ☐
Walk a real small **multi-contig** reference end to end; assert `.cat` parity (subset by cap), the
partition invariant, window-invariance, BED-invariance, and the edge cases (a tract at position 1 →
`Generic`; a repeat-free contig → one `Generic`; a tract at one contig's end abutting one at the next's
start). *Depends:* E1, E2. *Source:* spec §8. The port anchor.

> **Checkpoint E:** step 3 runs end to end; `.cat` parity + all three invariance tests green. **Step 3
> is complete.** Pause for review.

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
