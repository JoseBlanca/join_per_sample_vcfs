# `pop_var_caller_exp` / `type-regions` — types & interfaces

*Status: architecture draft (2026-07-19), companion to the spec
[`../spec/typed_regions_cli.md`](../spec/typed_regions_cli.md) (the design and its rationale), which
in turn drives the step-3 walk settled in [`../spec/typed_regions.md`](../spec/typed_regions.md) and
its arch companion [`typed_regions.md`](typed_regions.md). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): domain nouns for types,
verbs for functions; **STR** in prose ↔ `ssr` in code. Signatures are illustrative; the **contract**
is the deliverable. See the spec for the "why" behind every decision here — this doc carries the code
shape and points back.*

## Module home

A **second binary in the same crate** (spec §2): `pop_var_caller_exp`, so ng's experiment knobs stay
out of the production `pop_var_caller` CLI. The layout mirrors `src/pop_var_caller/`
([`mod.rs:6-31`](../../../../src/pop_var_caller/mod.rs)):

```
Cargo.toml            [[bin]] name = "pop_var_caller_exp", path = "src/main_exp.rs"   (NEW)
src/main_exp.rs       thin, like src/main.rs (main.rs:1-4): parse, dispatch, render error, exit  (NEW)
src/lib.rs            + pub mod pop_var_caller_exp;                                               (NEW)
src/pop_var_caller_exp/
  mod.rs              the module root
  cli.rs              top-level Parser + the PopVarCallerExpCommand subcommand enum
  typed_regions.rs    the `type-regions` subcommand: Args + run + error + #[cfg(test)]
  cli/parsers.rs      the --min-copies value_parser (its own; not an edit to production's)
```

Staying in-crate is what keeps it cheap: `pop_var_caller::common` is `pub(crate)`
([`mod.rs:8-9`](../../../../src/pop_var_caller/mod.rs)), so `DEFAULT_BUFFERED_IO_CAPACITY`, `basename`
and `current_command_line` are reused without a production edit (spec §2). **Do not** move `src/ng/`
into a crate of its own (spec §2).

## 1. The CLI surface

The shape copies `ssr-catalog` ([`ssr_catalog.rs:14-53`](../../../../src/pop_var_caller/ssr_catalog.rs)):
`--reference` / `--output` bare `PathBuf`s, `--regions` optional, every knob `default_value_t` a named
`pub const` under `help_heading = "Advanced"`.

```rust
/// The exp binary's subcommands. One inhabitant today (spec §9's naming question is deferred).
/// `TypeRegions` kebab-cases to the `type-regions` subcommand, as `SsrCatalog` → `ssr-catalog`.
pub enum PopVarCallerExpCommand {
    TypeRegions(TypedRegionsArgs),
}

/// `type-regions` arguments — the authoritative knob list; `run_typed_regions` translates it into a
/// `TypedRegionConfig`. Defaults are the library `pub const`s (spec §2.1); the four tagged §2.3 carry
/// the short-read values the consts do not yet hold (kept as record — the CLI slice applies them).
pub struct TypedRegionsArgs {
    pub reference: PathBuf,               // required; its .fai is read, or created if absent (spec T1)
    pub output: PathBuf,                  // required
    pub regions: Option<PathBuf>,         // BED; None = whole genome. Costs scan time, not memory (spec T10)

    // knobs — help_heading = "Advanced"
    pub min_period: u8,                   // DEFAULT_MIN_PERIOD  (spec §2.2)
    pub max_period: u8,                   // DEFAULT_MAX_PERIOD
    pub max_str_len: u64,                 // DEFAULT_MAX_STR_LEN
    pub window_bp: u64,                   // DEFAULT_WINDOW_BP   (memory only; not a comparison key)
    pub flank_bp: u64,                    // DEFAULT_FLANK_BP
    pub min_purity: f32,                  // DEFAULT_MIN_PURITY
    pub min_score: i32,                   // DEFAULT_MIN_SCORE
    pub min_copies: MinCopies,            // parsed by cli/parsers.rs (§2 below)
    pub scan_match_reward: i32,           // ScanParams defaults (tandem_repeat.rs:105-113)
    pub scan_mismatch_penalty: i32,
    pub scan_min_copies: u32,
}
```

**Contract.** `TypedRegionsArgs` is a plain clap `Args`; it holds no invariant of its own. The one
cross-knob rule — `max_str_len >= flank_bp` — is **not** checked here: the walk refuses the pair with
`TypedRegionError::MarginNarrowerThanFlank` before any work
([`mod.rs:884-889`](../../../../src/ng/region_typing/mod.rs)), so a second CLI check would be a second
place to drift (spec T3). Every knob is a flag, so no user input may panic (spec §6).

## 2. The `--min-copies` parser

`--min-copies` is a **table, not a number**, and no existing subcommand takes one (spec §2.1). Syntax
(decided, spec §2.1): **exactly six comma-separated values, one per period 1..6 — `6,4,4,3,3,3` — and
any other count is a hard parse error** (clap `value_parser` failure: usage, exit 2, before `run`).
No catch-all, no colon: `MinCopies::for_wider_periods` is unreachable here (the walk asserts the
period ceiling ≤ `MAX_MOTIF_LEN`, [`segment_criteria.rs:940-945`](../../../../src/ng/region_typing/segment_criteria.rs)),
so the parser supplies a fixed inert `3` for it. The parser lives in the exp binary's own
`cli/parsers.rs`, producing a `MinCopies` via `MinCopies::new`
([`segment_criteria.rs:419-430`](../../../../src/ng/region_typing/segment_criteria.rs)); production's
`cli/parsers.rs` is the *pattern* to copy, not a place to add to.

## 3. Errors

House pattern, inherited whole from `ssr-catalog`
([`ssr_catalog.rs:56-74`](../../../../src/pop_var_caller/ssr_catalog.rs)): a `#[non_exhaustive]`
`thiserror` enum, `run_typed_regions(&args) -> Result<(), TypedRegionsCliError>`. `main_exp.rs` walks
the source chain and exits 1, the way `main.rs` does — see the T7a decision on **sharing** that walk.

```rust
#[non_exhaustive]
pub enum TypedRegionsCliError {
    Reference(ReferenceInfoError),             // #[from]  reference_info (T1) — read / verify / write-fai
    ContigTooLong { contig: String, len: u64 }, // u64 → u32 narrowing for RegionSet (spec T2)
    Bed(BedError),                             // #[from]  (regions.rs:315-379)
    Walk(TypedRegionError),                    // #[from]  the walk's own errors, incl. the flank pair (T3)
    Output(std::io::Error),                    // #[from]  writing / atomic rename
}
```

**Contract.** Every variant states when it fires (per-variant doc, per `RefSeqError`'s shape). **Not**
a `--max-str-len`/`--flank-bp` variant — `TypedRegionError::MarginNarrowerThanFlank` already carries
both numbers (spec T3). No `anyhow` in the CLI path (spec §6).

## 4. Interfaces

```rust
/// Run the walk and write the partition. Streams — never collects (spec §6).
pub fn run_typed_regions(args: &TypedRegionsArgs) -> Result<(), TypedRegionsCliError>;
```

The reference step is **`reference_info`, reused as-is** (spec T1) — one batteries-included call
`read_reference_verifying_or_creating_fai(cache: &Arc<ReferenceInfoCache>, fasta: PathBuf) ->
(Arc<ReferenceInfo>, Option<VerificationHandle>)` (ng's `reference_info`, `ng-contig-table`,
merge-bound). `.fai` present → the index is read in the foreground and the FASTA is verified on a
background thread (`Some(handle)`); `.fai` absent → the FASTA is read once, the sibling `.fai` is
written, and `None` comes back. `.contig_list()` gives the `ContigList`. The `VerificationHandle` is
`#[must_use]` and **must be `join()`ed before the output is committed** — see step 6 and the decision
in §6.

**The `run_typed_regions` contract**, step by step (each a spec trap it discharges):

1. `let (info, verify) = read_reference_verifying_or_creating_fai(&cache, args.reference.clone())?;`
   then `info.contig_list()` → one `ContigList`, used **twice**: for `WindowedRefSeq::new` *and* for
   `GenomeRegions` (spec T8 — else the walk's `UnknownContig` error renders a `ContigId` index, not a
   name). Hold `verify: Option<VerificationHandle>` for step 6.
2. `GenomeRegions` from `--regions` (`from_bed_path`) or whole-genome (`whole_contigs`), over
   `&[ContigBounds]` derived from the same `ContigList`. The `u64 → u32` contig-length narrowing lives
   here and fails loudly (`u32::try_from`, spec T2) — never `as`.
3. `TypedRegionConfig` from the knobs; `WindowedRefSeq::new(reference.clone(), contigs)`
   ([`ref_seq.rs:524`](../../../../src/ng/ref_seq.rs)).
4. `TypedRegionIterator::over_regions(reference, spans, config)`
   ([`mod.rs:838`](../../../../src/ng/region_typing/mod.rs)) — fallible setup (contig cross-check, the
   flank/margin pair).
5. Write the header block, then stream `for r in iter { writer.write_row(r?)?; }` to a `.tmp` file.
   **Then `if let Some(h) = verify { h.join()?; }` — the FASTA-verification barrier — and only then
   rename** `.tmp` into place (spec §6 — a half-written *or* stale-`.fai` partition is silently valid
   otherwise; [`cli.rs:621-643`](../../../../src/pop_var_caller/cli.rs) is the atomic-write precedent).
   The background check overlapped the whole walk, so the join costs almost nothing. No sort (the walk emits
   in genomic order, [`mod.rs:575`](../../../../src/ng/region_typing/mod.rs)); no collect.
6. Announce line to stderr at the start, a `key: k=v` summary at the end from
   `TypedRegionIterator::counts()` — subject to spec T9 (four of five rejection counters are
   structurally zero; print all, labelled, do not hard-code which).

## 5. The output format (spec §3)

Typed columns, `.` for absent; a bundle's `members` cell carries a **JSON array** (spec §3.1). Plain
text, not bgzip (spec §3.2): read end to end, and it stays a BED our own `--regions` accepts. The row
maps from `TypedRegion { region: GenomeRegion, kind: RegionKind }`
([`mod.rs:134`](../../../../src/ng/region_typing/mod.rs)):

| kind | columns 5-9 (`motif period copies purity members`) |
|---|---|
| `Generic`, `Satellite` | all `.` |
| `SsrSegment(SsrSegment)` | `motif`=`SsrSegment::motif()` bytes, `period`=`motif.period()`, `copies`=`tract_len()/period()` (fractional, spec T5), `purity`=`purity_fraction()`, `members`=`.` ([`segment_criteria.rs:284-330`](../../../../src/ng/region_typing/segment_criteria.rs)) |
| `SsrBundle { tracts }` | `members`=JSON array `[{"start","end","period"}, …]` over `RepeatInterval` ([`tandem_repeat.rs:207-217`](../../../../src/ng/tandem_repeat.rs)); rest `.` |

**Contract — coordinates.** Output is BED (0-based half-open). The hull region is 1-based inclusive
([`types.rs:55-59`](../../../../src/ng/types.rs)) → `start - 1`; the member `RepeatInterval`s are
already 0-based half-open and contig-based ([`mod.rs:1458-1465`](../../../../src/ng/region_typing/mod.rs))
→ **no shift**. Two directions in one row: write the conversion **once, in one helper** (spec T4).
`SsrSegment::chrom()` carries the contig **name**; `Generic`/`Satellite`/`SsrBundle` do not, so name
them from the shared `ContigList` (spec §4). `Motif` has no `Display` — bind `motif()` (it is `Copy`)
then `str::from_utf8` ([`segment_criteria.rs:184-193`](../../../../src/ng/region_typing/segment_criteria.rs), spec T6).

**Contract — header + determinism.** A `## key: value` block (every resolved config value, so two
files compare and one reproduces), one `#`-prefixed column header, then rows — `ssr-catalog`'s shape
([`io.rs:56-69`](../../../../src/ssr/catalog/io.rs)). **No `date`, no `reference_md5`** (spec §3.4,
§6): the walk is a pure function of reference+regions+config, so the same three inputs must give a
byte-identical file — fixed key order, no timestamp, and the `members` JSON serialises
deterministically (fixed key order, no incidental whitespace). `window_bp` is recorded but is not a
comparison key (memory knob).

## 6. Design decisions — decided

Distilled from the spec; see it for the reasoning. Open items marked `OPEN:`.

- **Second binary, same crate — decided.** Splits the *command surface*, not the code; both binaries
  link the one library that contains ng. Keeps ng's knobs off the production CLI (spec §2).
- **Output = typed columns + JSON `members` cell — decided.** Tabular where the data is scalar
  (greppable, valid BED, pandas-direct, `ssr-catalog` precedent), JSON in the one genuinely-nested
  cell. Rejected: a pure-JSON row (forfeits terminal/`awk` + the `--regions` round-trip for one cell's
  benefit; Python is a wash) (spec §3.1).
- **Plain text, not bgzip — decided.** Read end to end, not random-access; keeps `awk`/`diff`/`head`
  and the `--regions` round-trip. Rejected: bgzip for `ssr-catalog` consistency (buys an unused index)
  (spec §3.2).
- **Stream, do not collect; write atomically — decided.** The walk holds ~102 kb and three
  coordinates; collecting a whole-genome partition puts back the memory it exists to avoid. `.tmp` +
  rename because a truncated partition is otherwise indistinguishable from a complete smaller one (spec
  §6).
- **Reference read = verify-in-background / create-`.fai`, joined at the commit — decided.**
  `read_reference_verifying_or_creating_fai` (owner, 2026-07-19): a present `.fai` is read in the
  foreground and the FASTA verified on a background thread; an absent one is created from a FASTA scan.
  The `#[must_use]` `VerificationHandle` is `join()`ed **immediately before the atomic rename**, so a
  stale `.fai` aborts the run with the output still under its temp name — never a partition built
  against a wrong contig table renamed into place. The background pass overlaps the walk, so the
  barrier is nearly free. Rejected: the bare `ReferenceSource::Fai` (skips verification) and `Fasta`
  (blocks the walk on the whole-genome read) (spec T1, §6).
- **Propagate the flank/margin guard, don't re-check it — decided.** `TypedRegionError` carries both
  numbers; a CLI check is a second place to drift (spec T3).
- **Hoist `format_error_chain` into the library — decided.** It is private to `src/main.rs`
  ([`main.rs:30-44`](../../../../src/main.rs)); `main_exp.rs` needs the same rendering. A pure move (no
  behaviour change) so both binaries share one rule. Rejected: copy the 15 lines (two copies of an
  error-rendering rule diverge the first time one is fixed) (spec T7a). **The one production file this
  work touches, and only by moving a function.**
- **Re-aim the B2 guard, don't widen it loosely — decided.** `rg 'use crate::ng' src/ --glob
  '!src/ng/**'` gains `--glob '!src/pop_var_caller_exp/**'`: the exp binary depending on ng is the safe
  direction (production still depends on nothing in ng). The excluded path is now load-bearing (spec
  T7). Mechanically, name `crate::ng::region_typing::…` in full or add the re-export (`ng/mod.rs`
  re-exports only `ref_seq`/`types`).
- **`--min-copies` = exactly six values, no catch-all — decided.** One value per period 1..6; any
  other count is a hard parse error, and the unreachable `for_wider_periods` is supplied inert by the
  parser rather than exposed. The user never types a field that decides nothing (spec §2.1).
- **`pop_var_caller_exp` ships in release artefacts — decided.** The §2 separation is of the command
  surface, not the shipped file list; the tree packages nothing today, so there is nothing to exclude
  it from (spec §9).

## 7. Reconciliation with existing code

Verify each row against the code when implementing (all cited `file:line` were read for this doc).

| CLI name | existing code | action |
|---|---|---|
| binary skeleton (`main_exp.rs`, `cli.rs`, subcommand mod) | `src/main.rs`, `src/pop_var_caller/{mod,cli,ssr_catalog}.rs` | **new**; mirror the shape |
| `TypedRegionsArgs` / `run_typed_regions` / `TypedRegionsCliError` | `SsrCatalogArgs` / `run_ssr_catalog` / `SsrCatalogCliError` ([ssr_catalog.rs:14-110](../../../../src/pop_var_caller/ssr_catalog.rs)) | **new**; same pattern |
| `--min-copies` parser | `pop_var_caller/cli/parsers.rs` | **new** parser in `pop_var_caller_exp/cli/parsers.rs`; copy the pattern |
| reference → `ContigList` + verify | `reference_info::read_reference_verifying_or_creating_fai` + `ReferenceInfo::contig_list()` + `VerificationHandle::join` (ng `reference_info`, `ng-contig-table`) | reuse as-is; join the handle before the rename (T1, §6) |
| walk driver | `TypedRegionIterator::over_regions` ([mod.rs:838](../../../../src/ng/region_typing/mod.rs)) | reuse as-is |
| reference accessor | `WindowedRefSeq::new(PathBuf, ContigList)` ([ref_seq.rs:524](../../../../src/ng/ref_seq.rs)) | reuse as-is |
| region set | `GenomeRegions::{whole_contigs, from_bed_path}` ([mod.rs:102](../../../../src/ng/region_typing/mod.rs)), over `ContigBounds` ([regions.rs:57](../../../../src/regions.rs)) | reuse; owns the `u32` narrowing (T2) |
| config / counts | `TypedRegionConfig` ([mod.rs:209](../../../../src/ng/region_typing/mod.rs)), `TypedRegionCounts` ([mod.rs:265](../../../../src/ng/region_typing/mod.rs)) | reuse; serialise config into the header |
| row source | `TypedRegion` / `RegionKind` ([mod.rs:134-180](../../../../src/ng/region_typing/mod.rs)), `SsrSegment` accessors ([segment_criteria.rs:284-330](../../../../src/ng/region_typing/segment_criteria.rs)), `RepeatInterval` ([tandem_repeat.rs:207-217](../../../../src/ng/tandem_repeat.rs)) | reuse; format to columns |
| error renderer | `format_error_chain` ([main.rs:30-44](../../../../src/main.rs)) | **hoist to lib**, share (T7a) |
| header / atomic write idioms | `ssr-catalog` header ([io.rs:56-69](../../../../src/ssr/catalog/io.rs)); `pileup` `.tmp`+rename ([cli.rs:621-643](../../../../src/pop_var_caller/cli.rs)) | copy the idiom |
| in-crate reuse | `pop_var_caller::common` `pub(crate)` ([mod.rs:8-9](../../../../src/pop_var_caller/mod.rs)) | reuse (no production edit) |

## 8. Open items

- **Impl-time confirmations (not design):** `reference_info` merges from `ng-contig-table` before this
  is built (owner, 2026-07-19), so `read_reference_verifying_or_creating_fai`, `ReferenceInfoCache`,
  `VerificationHandle::join` and `ReferenceInfoError` are as read from that branch (the CLI holds one
  `Arc<ReferenceInfoCache>`); whether `run_typed_regions` also wants `over_contig` (the whole-genome
  path already works via `whole_contigs` + `over_regions`).

## 9. Test & bench shape (spec §7)

Tests live in `typed_regions.rs`'s own `#[cfg(test)]` module (per-subcommand convention), standing on
the anchor fixture that already drives the shipping stack from a real multi-contig FASTA + `.fai`
([`anchor.rs:125-138`](../../../../src/ng/region_typing/anchor.rs)). They test **the file**, not the
walk again:

- **Round-trip** — parse the written rows back, assert field-for-field against the iterator's own
  output (pins T4's two-direction conversion; needs a fixture with a bundle).
- **Partition invariant in the file** — concatenated spans reconstruct the requested regions (catches
  a dropped generic row, a miscounted header line, a T4 off-by-one).
- **Determinism** — two runs byte-identical (§5).
- **`--regions` read-back** — feed the output back as `--regions`, assert it walks the same spans (§3.2's
  claimed property).
- **Flag-pair guard is a CLI error, not a panic** (T3) — mutation-verify: remove the check and it must
  fail with a panic, not pass.
- **Stale `.fai` aborts before publishing** (T1) — corrupt the sibling `.fai`, assert the run errors
  and leaves no output file (the verification join precedes the rename).

No `bench/`: the subcommand orchestrates a single-threaded walk; the reference verification runs on one
background thread inside `reference_info`, not the walk. No competing implementation.
