# ng typed-region generator (step 3) — implementation plan

**Status:** draft, 2026-07-16. Build order for **step 3, the typed-region generator**: the
`region_typing.rs` module and the walk that cuts the reference into `TypedRegion`s, plus the
production-side changes it stands on (the `Locus` rebase, a windowable `build_loci`, a streamed
`collect_windowed`, ng's `u64` widening). Design is settled in
[`../spec/typed_regions.md`](../spec/typed_regions.md) (spec) and
[`../arch/typed_regions.md`](../arch/typed_regions.md) (types & interfaces). This turns that design
into build order; it is **not** a place for new design — every open item is a *parameter value*
resolved in spec §10, not an interface question.

Much of this plan is *enabling work in existing code* — ng is early, so where the code doesn't match
a step-3 decision we change it (spec §9) rather than work around it. Those changes each verify
against an existing oracle before the ng walk is written.

---

## Scope

**In:** the `Locus` rebase (1-based + `u64`); `CatalogParams` knobs; a windowed `build_loci`
returning `Admitted { loci, bundled }`; a promoted+streamed `collect_windowed`; ng's `u64` coordinate
widening + `WindowedRefSeq` raw bytes; `ng::types` `GenomeRegion`/`Position`; `region_typing.rs` with
`TypedRegion`/`RegionKind`/config/counts/error; `GenomeRegions`; the resident partition, bundle
detection, and the windowed walk; `TypedRegionIterator`; BED scan-wider-than-emit.

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
- **The substrate before the walk.** Everything the walk *calls* — the rebased `Locus`, the windowed
  `build_loci`, the streamed `collect_windowed`, the raw `WindowedRefSeq` — is built and re-verified
  against its own oracle (Milestones A–B) before a line of the ng walk exists.
- **Simplest impl first, as the next one's oracle.** The **resident** partition (D1, whole contig in
  memory) is built and proven first; the **windowed** walk (D3) is then proven *by matching it* —
  which is exactly the window-invariance the spec demands (§2.3, §2.6).
- **The algorithmic heart before the plumbing.** The partition + bundle detection (D) is decided and
  tested before the iterator, ownership, and BED narrowing (E) wire it to the surface.
- **Isolate the silent step.** The `Locus` rebase (A1) is a coordinate change: a bug is a *wrong
  genotype, not a crash*. It lands as its own commit, golden fixture green **before and after**, never
  bundled — so a `git bisect` can find it.
- **Verify against ground truth.** North star for the catalog changes: **byte-parity with the golden
  catalog**. For the walk: **`.cat` parity** (subset by the satellite cap) plus the partition and
  invariance tests. Not self-consistency.
- **Incremental, with pauses.** One milestone, stop for review, next.
- **Container builds.** All `cargo` via `./scripts/dev.sh` (CLAUDE.md); native host build at completion.

## Preconditions (already in place, or the prior plan runs first)

- **The tandem-repeat scanner is built** — `find_tandem_repeats`, `collect_windowed` (private today),
  `RegionScanner` ([ng/tandem_repeat.rs](../../../../src/ng/tandem_repeat.rs)).
- **The trf-mod → scanner catalog swap has landed** (`ssr_repeat_scanner.md` §6–§7): `build_loci` takes
  `Vec<RepeatInterval>`, `TrfRecord` is gone, `catalog_prefilter` is a real `pub(crate)` fn in
  `src/ssr/catalog/`. If not, **that plan runs first** — step 3 cannot call `build_loci` from non-test
  code until it does (spec §5c).
- **`RegionSet`** ([regions.rs](../../../../src/regions.rs)), **`Locus` / `build_loci` / `CatalogParams`**,
  the **golden catalog fixture** (`tests/data/tandem_repeat/golden.ssr_catalog.bed.gz` + `scanner_parity.rs`)
  as the oracle, and the **`RefSeq` impls** all exist.

---

## The steps

### Milestone A — catalog coordinates & knobs (production; each vs the golden catalog)

**A1. Rebase `Locus` to 1-based inclusive, widen to `u64`.**  ☐
All ~59 `start()`/`end()`/`ref_bytes_start` sites in `src/ssr/`; the arithmetic cancels through the one
`tract_range()` (only the length gains a `+1`); the single non-cancelling site is `build_loci`'s index
into the raw `contig_seq`. The `.cat` format stays 0-based — convert at `catalog/io.rs`, the pattern
`regions.rs` set. **Silent failure: own commit, do not bundle. Golden catalog + `scanner_parity` +
`benchmarks/ssr_tomato1` green before AND after.** *Source:* spec §4, §9; arch §recon.

**A2. `CatalogParams` — collapse the redundant knob, expose the hidden ones, widen.**  ☐
Fold `bundle_threshold` into `flank_bp` (spec §2.4 — one number; confirm no shipped catalog set them
apart). Hoist `MIN_PERIOD` / `MAX_PERIOD` / `copy_number_floor` from hardcoded consts into fields
(`Default` = today's values) and reconcile the pre-filter's second copy-floor table with them. Widen
to `u64`. *Depends:* A1. *Source:* spec §2.4, §5, §9. Verified: golden **byte-identical at `Default`**.

**A3. `build_loci` windowed, and returns what it set aside.**  ☐
Add `bases_start` + `contig_len` params (whole-contig = `0, len`, the degenerate case → unchanged
output). Return `Admitted { loci, bundled }`, where `bundled` is the cluster members it silently drops
today (spec §6a). *Depends:* A1, A2. *Source:* spec §5a, §6a. Verified: golden `loci` unchanged;
`bundled` non-empty on a clustered fixture.

> **Checkpoint A:** the catalog is 1-based/`u64`, its rules are parameters, and `build_loci` is
> windowable and reports bundles — golden green throughout. Pause for review.

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

**C2. `region_typing.rs` scaffold + local types.**  ☐
`src/ng/region_typing.rs` (+ `#[cfg(test)]`), wired into `ng/mod.rs`. `TypedRegion`, `RegionKind`
(`SsrLocus(Locus)` / `SsrBundle` / `Generic` / `Satellite`), `TypedRegionConfig` (+`Default` = the
catalog's settings, for comparability), `TypedRegionCounts`, `TypedRegionError`. No logic. *Depends:*
A, B, C1. *Source:* arch §types. Verified: `Default` config test.

**C3. `GenomeRegions` wrapping `RegionSet`.**  ☐
`whole_contigs` (default), `from_bed_path`, `iter()` widening `RegionSet`'s `u32` → `u64`. Reimplements
nothing — `RegionSet` already parses/coalesces/clamps/converts-BED/drops-empty-contigs. *Depends:* C1.
*Source:* spec §2.5, arch §types. Verified: `whole_contigs` over a contig table; a BED round-trip.

> **Checkpoint C:** ng types compile; `GenomeRegions` wraps `RegionSet` and widens. Pause for review.

### Milestone D — the walk (the heart); resident first as the windowed oracle

**D1. The resident partition — the five steps, whole contig in memory.**  ☐
`detect → clean → admit → cap → partition` over a contig held resident, producing `Vec<TypedRegion>`
(no windowing). Uses `find_tandem_repeats` + `catalog_prefilter` + `build_loci` (whole-contig form) +
the coverage merge + satellite cap. *Depends:* A, B, C. *Source:* spec §2.1, §2.2, §8. Verified: the
**partition invariant** (contiguous / non-overlapping / complete / maximal) **and `.cat` parity** (strict
subset — every catalog locus present, or absent *and* inside a satellite run) on a small reference.
**This is D3's oracle.**

**D2. Bundle detection — the flank test, and the rejected-repeat split.**  ☐
A repeat with another repeat within `flank_bp` on either side is a bundle member; the cluster (hull of
its tracts) is one `SsrBundle`; a repeat admission turns down for any *other* reason is `Generic`, not a
hole (spec §2.2). Consumes `build_loci`'s `Admitted.bundled` (A3). *Depends:* D1. *Source:* spec §2.2,
§2.4. Verified: two tracts 10 bp apart → one `SsrBundle` carrying both; three chained 30 bp apart → one
bundle of three; a locus 60 bp from a bundle → admitted; an impure/low-copy repeat → `Generic`.

**D3. The windowed walk.**  ☐
The three carries (open cluster, open coverage run, open generic run); the `max_repeat_len` detection
margin; the contig length **passed in** to `build_loci`, never inferred from the slice (spec §2.6 — the
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
| A | **golden catalog byte-parity** through the `Locus` rebase, the knob changes (identical at `Default`), and `build_loci` (identical `loci`, new `bundled`) — the isolation oracle for the silent step |
| B | window-count-invariance still green (streamed `collect_windowed`); ng suite green (`u64`, clamp gone); raw-fetch unit + `contigs()` (`WindowedRefSeq`) |
| C | type tests — `Default` config; `GenomeRegions` `whole_contigs` + BED round-trip |
| D | **partition invariant + `.cat` parity** (resident, D1); bundle-routing cases (D2); **window-invariance** — windowed == resident (D3) |
| E | iterator + running counts + fatal-in-stream (E1); **BED-invariance** + straddle-whole (E2); **integration:** `.cat` parity + all invariants on a real multi-contig reference (E3) |

## Out of scope (next plans)

- **The pileup / generic join** — splitting a `Generic` region into loci from the data (`pileup/`'s
  spec); it owns "what counts as a generic locus", the bake-off this step declined to fake.
- **The STR gatherer** — `SsrLocus` region → `LocusEvidence` (fetch + prepare + tally), and the
  region-driven `RecordSource` it needs. Its own spec.
- **Bundle disposal** — the `SsrBundle` consumer (spec §10); the gatherer's call.
- **The routing/threshold experiments** — period × length, flank size, satellite cap (spec §10); config
  sweeps on the finished walk, scored on the standards, not build order.
