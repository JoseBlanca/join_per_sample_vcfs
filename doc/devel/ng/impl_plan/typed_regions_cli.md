# ng `type-regions` CLI (`pop_var_caller_exp`) — implementation plan

**Status:** draft, 2026-07-19. The build order for **the `pop_var_caller_exp` binary and its
`type-regions` subcommand**: a second binary in the same crate, its CLI surface, the `--min-copies`
parser, the output writer, and the `run_typed_regions` driver that streams step 3's walk to a file.
Design is settled in [`typed_regions_cli.md`](../spec/typed_regions_cli.md) (spec) and
[`../arch/typed_regions_cli.md`](../arch/typed_regions_cli.md) (types & interfaces), which stand on the
step-3 walk ([`typed_regions.md`](../spec/typed_regions.md)). This roadmap turns that design into build
order; it is **not** a place for new design — every fork is resolved in the spec (§3.1 format, §2.1
`--min-copies`, §9 release artefacts).

---

## Scope

**In:** the `Cargo.toml` `[[bin]]`, `src/main_exp.rs`, `src/pop_var_caller_exp/{mod,cli,typed_regions}.rs`
and `cli/parsers.rs`; the `PopVarCallerExpCommand` enum, `TypedRegionsArgs`, `TypedRegionsCliError`; the
`--min-copies` value_parser; the header + row output writer; the `run_typed_regions` driver; the
`format_error_chain` hoist (the one production edit — a pure move); the T7 guard-grep re-aim.

**Out (later plans / deferred, spec §8):**

- **ng's other exp drivers** — a `pileup` subcommand, the gatherer, the sweep runners → their own
  specs; they land in this binary (spec §2).
- **A reader for this format** and **a tabix/CSI index (+ bgzip)** → whoever first consumes the file
  from disk / needs random access (spec §8).
- **Applying the §2.3 short-read default changes to the library `DEFAULT_*` consts** → this CLI slice
  applies them (spec §2.3, kept-as-record 2026-07-19): the const edits + test re-pinning
  (delete `default_matches_the_frozen_catalog_params`, re-pin the homopolymer test, pin the `.cat`
  anchor's knobs from its header) land in **Milestone B here**, not upstream.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The algorithmic heart before the plumbing.** The novel logic is the **output writer** — the typed
  columns, the JSON `members` cell, and above all the T4 coordinate conversion (a *silent* failure:
  an off-by-one is a wrong region, not a panic). It is built and unit-tested against hand-built /
  iterator `TypedRegion`s (Milestone C) *before* the `run` driver wires the reference and file I/O
  (Milestone E).
- **Reuse over rewrite.** This is a **driver**: it calls `TypedRegionIterator::over_regions`,
  `WindowedRefSeq::new`, `GenomeRegions`, `reference_info`, and the `SsrSegment`/`RepeatInterval`
  accessors as-is (arch §7). No walk logic is re-derived; the only new *code* is the CLI shell, the
  parser, and the formatter.
- **Isolate the silent step.** The T4 coordinate conversion lands as **its own commit, its round-trip
  oracle green before and after** — never bundled — so a `git bisect` can find it if a coordinate
  moves (skill: isolate a step whose failure is silent).
- **Verify against ground truth.** The north-star test is the **round-trip**: parse the written rows
  back and assert they equal the iterator's own output, field for field (spec §7) — not
  self-consistency.
- **Incremental, with pauses.** One milestone, stop for review, then the next.
- **Ungated / container builds.** `cargo` via `./scripts/dev.sh` (CLAUDE.md), a native host build at
  completion.

## Preconditions (already in place)

- **Step 3 is complete** — `TypedRegionIterator::over_regions`, `WindowedRefSeq::new`, `GenomeRegions`,
  `TypedRegionConfig`, `TypedRegionCounts`, `TypedRegion`/`RegionKind`, `SsrSegment`, `RepeatInterval`,
  `MinCopies::new`, and the anchor fixture ([`anchor.rs`](../../../../src/ng/region_typing/anchor.rs))
  that drives the shipping stack from a real FASTA + `.fai`.
- **`reference_info` merged from `ng-contig-table`** (owner, 2026-07-19: merges to `main` before this
  is built) — `read_reference_verifying_or_creating_fai`, `ReferenceInfoCache`,
  `ReferenceInfo::contig_list()`, `VerificationHandle::join`, `ReferenceInfoError`. **The executor
  confirms this is on `main` before step D1.**
- **Templates to copy:** `ssr-catalog` for the subcommand shape
  ([`ssr_catalog.rs`](../../../../src/pop_var_caller/ssr_catalog.rs)); `pileup`'s `.tmp`+rename for the
  atomic write ([`cli.rs:621-643`](../../../../src/pop_var_caller/cli.rs)); `ssr-catalog`'s header
  ([`io.rs:56-69`](../../../../src/ssr/catalog/io.rs)); `format_error_chain`
  ([`main.rs:30-44`](../../../../src/main.rs)).

---

## The steps

### Milestone A — the binary skeleton (scaffold, no logic)

**A1. The second binary exists and runs.**  ✅
`Cargo.toml` `[[bin]] name = "pop_var_caller_exp", path = "src/main_exp.rs"`; `src/main_exp.rs` thin
(parse → dispatch → render error → exit, like `main.rs:1-4`); `pub mod pop_var_caller_exp;` in
`lib.rs`; `src/pop_var_caller_exp/{mod,cli}.rs` with the top-level `Parser` and
`PopVarCallerExpCommand { TypeRegions(TypedRegionsArgs) }` (a stub `run` for now). Prove:
`pop_var_caller_exp --help` and `type-regions --help` run. *Source:* arch §Module home, spec §2.

**A2. Hoist `format_error_chain` into the library.**  ✅ — **own commit (production edit, pure move).**
Move the fn (and its `#[cfg(test)]`) from `src/main.rs` into a library module both binaries import;
`main.rs` now calls the hoisted one. No behaviour change — the existing test stays green. *Source:*
spec T7a, arch §6. *Depends:* A1.

**A3. Re-aim the B2 dependency guard.**  ✅
`rg 'use crate::ng' src/ --glob '!src/ng/**'` gains `--glob '!src/pop_var_caller_exp/**'` (wherever
the guard lives). Confirm production still names nothing in ng. *Source:* spec T7, arch §6. *Depends:*
A1.

> **Checkpoint A:** the binary builds and `--help`s; both binaries share one error renderer; the guard
> passes with the exp path excluded. Pause for review.

### Milestone B — the input surface (Args, parser, defaults)

**B1. Apply the §2.3 short-read defaults to the library consts.**  ✅ — **own commit; oracle: the ng
test suite green after re-pinning.**
Edit `DEFAULT_MIN_PERIOD` 2→1, `DEFAULT_FLANK_BP` 50→30, `DEFAULT_MAX_STR_LEN` 1000→100,
`MinCopies::default()` `[10,5,4,3,3,3]`→`[6,4,4,3,3,3]`; delete
`default_matches_the_frozen_catalog_params`; re-pin `a_homopolymer_is_generic_and_does_not_take_its_neighbour_with_it`
to `--min-period 2` and add the "homopolymer is a period-1 locus" test; pin the `.cat` anchor's
`min_copies`/`max_str_len` from its header (`typed_regions.md` §8.1). *Source:* spec §2.3. *Depends:* —
(library-only; independent of A).

**B2. `TypedRegionsArgs` + `TypedRegionsCliError`.**  ✅
The `Args` struct — `--reference`/`--output`/`--regions` + every knob `default_value_t` the `DEFAULT_*`
const under `help_heading = "Advanced"` — and the `#[non_exhaustive]` error enum (`Reference`,
`ContigTooLong`, `Bed`, `Walk`, `Output`; `#[from]` where it fits). No logic. Test: parses through the
top-level CLI with defaults (mirror `ssr_catalog`'s `parses_through_the_top_level_cli`), asserting the
resolved defaults are the §2.3 values. *Source:* arch §1, §3; spec §2.1. *Depends:* A1, B1.

**B3. The `--min-copies` value_parser.**  ✅
In `pop_var_caller_exp/cli/parsers.rs`: parse **exactly six** comma-separated `u32`s → `MinCopies::new([…;6], 3)`
(the `for_wider_periods` supplied inert); **any other count is a hard parse error**. Wire it as the
`min_copies` field's `value_parser`. Unit tests: six values parse to the right table; five and seven
error; non-numeric errors; the error is a clap usage failure. *Source:* spec §2.1, arch §2. *Depends:*
B2.

> **Checkpoint B:** the CLI parses every knob at its (now short-read) defaults; `--min-copies` accepts
> six values and hard-errors on any other count; the ng suite is green after the default re-pinning.
> Pause for review.

### Milestone C — the output writer (the heart)

**C1. The header block.**  ☐
Write `## key: value` lines from the *resolved* `TypedRegionConfig` (every knob, fixed key order),
then the `#`-prefixed column header. **No `date`, no `reference_md5`** (determinism). Unit test:
byte-identical across two calls; every config knob appears; `window_bp` present. *Source:* spec §3.4,
§6; arch §5. *Depends:* B2.

**C2. The row formatter + the T4 coordinate conversion.**  ☐ — **own commit, do not bundle; round-trip
oracle green before and after.**
`TypedRegion` → one BED row: typed columns with `.` for absent; `motif` via `str::from_utf8` on a bound
`motif()` (T6); `copies` = `tract_len()/period()`, fractional (T5); `purity` from `purity_fraction()`;
`members` a **JSON array** `[{"start","end","period"}, …]` for a bundle (deterministic key order), `.`
otherwise. The conversion, in **one helper**: hull 1-based-inclusive → BED 0-based (`start − 1`);
member `RepeatInterval`s already 0-based half-open contig coords → **no shift** (T4). Oracle: a
round-trip unit test over hand-built `Generic`/`Satellite`/`SsrSegment`/`SsrBundle` rows (a bundle
**required**) — format, parse back, assert equal to the input. *Source:* spec §3.1, §4, T4/T5/T6; arch
§5. *Depends:* C1.

> **Checkpoint C:** the writer emits a deterministic header and one row per kind; the round-trip test
> pins the T4 conversion (both coordinate systems, a bundle's members). Pause for review.

### Milestone D — the fallible setup assembly

**D1. Reference read + contig table + region set.**  ☐
`read_reference_verifying_or_creating_fai(&cache, args.reference.clone())?` (one held
`Arc<ReferenceInfoCache>`) → `(Arc<ReferenceInfo>, Option<VerificationHandle>)`; keep the handle for
E1. `info.contig_list()` → one `ContigList`, used **twice** (T8): `WindowedRefSeq::new` and
`GenomeRegions`. The `u64 → u32` contig-length narrowing lives here, via `u32::try_from` →
`ContigTooLong` on overflow, never `as` (T2). `GenomeRegions::from_bed_path` for `--regions`, else
`whole_contigs`. Unit tests: T2 rejects a >4 Gb contig; a BED and the whole-genome path both build;
the same `ContigList` feeds both; **`.fai` absent → the sibling `.fai` is written**. *Source:* spec
T1/T2/T8, arch §4. *Depends:* B2; **precondition: `reference_info` on `main`.**

> **Checkpoint D:** the reference becomes a `ContigList` via `reference_info` (`.fai` read foreground +
> FASTA verified in background, or FASTA scanned and `.fai` written); the narrowing fails loud; regions
> build for both paths. Pause for review.

### Milestone E — the driver, end to end

**E1. `run_typed_regions`.**  ☐
Assemble `TypedRegionConfig` from the args (wrapping raw knobs in `Bp`/`Position`),
`WindowedRefSeq::new`, `TypedRegionIterator::over_regions` (fallible setup — the flank/margin pair
surfaces here, T3), then **stream** `for r in iter { writer.write_row(r?)?; }` to a `.tmp` file;
**`if let Some(h) = verify { h.join()?; }` (the FASTA-verification barrier, T1), then rename** (atomic;
no sort, no collect) — the join goes *before* the rename so a stale `.fai` aborts with nothing
published. A stderr announce line at the start and a `key: k=v` summary from `counts()` at the end —
print **all five** rejection counters, labelled, hardcoding none (T9). *Depends:* B, C, D. *Source:*
spec §6, arch §4.

**E2. Integration tests — the file, on the anchor fixture.**  ☐
On the anchor FASTA+`.fai`: (1) **round-trip** — written rows parse back field-for-field equal to the
iterator's output (bundle fixture required); (2) **partition invariant** — concatenated spans
reconstruct the requested regions; (3) **determinism** — two runs byte-identical; (4) **`--regions`
read-back** — feed the output back as `--regions`, assert the same spans; (5) **the flag-pair guard is
a CLI error, not a panic** (T3) — mutation-verify: remove the propagation and it must panic, not pass;
(6) **a stale `.fai` aborts before publishing** — corrupt the sibling `.fai`, assert the run errors and
no output file is left behind (the join barrier, T1). Home: `typed_regions.rs`'s `#[cfg(test)]`.
*Source:* spec §7, arch §9. *Depends:* E1.

> **Checkpoint E:** `type-regions` runs end-to-end on a real FASTA, writes a deterministic partition,
> and round-trips through its own `--regions`. **The subcommand is complete.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | `--help` runs for binary + subcommand; `format_error_chain` test green after the move; the guard grep passes |
| B | ng suite green after the default re-pinning; top-level-CLI parse test at the §2.3 defaults; `--min-copies` count-mismatch unit tests |
| C | **round-trip unit test** pinning the T4 conversion over all four kinds incl. a bundle; deterministic-header test |
| D | unit tests — T2 narrowing rejects >4 Gb; BED + whole-genome build; one `ContigList` feeds both; `.fai`-absent writes the sibling index |
| E | **integration on the anchor fixture:** round-trip, partition-invariant, determinism, `--regions` read-back, flag-pair-guard-is-an-error, stale-`.fai`-aborts-before-publishing (all mutation/fixture-verified) |

## Out of scope (next plans)

- **ng's other exp drivers** (`pileup`, gatherer, sweep runners) — their own specs; they land in this
  binary (spec §2).
- **A reader for this format**, **a tabix/CSI index + bgzip**, **deriving the file from a non-FASTA
  source** — deferred with homes in spec §8.
- **The first §10 sweep** that consumes the file — the experiment this artefact exists to enable, not
  a build step.
