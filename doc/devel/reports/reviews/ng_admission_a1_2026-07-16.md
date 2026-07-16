# Code review — ng typed-region generator, Milestone A1 (the admission port)

**Date:** 2026-07-16 · **Scope:** working-tree diff (uncommitted) · **Reviewer:** 9 parallel
category sub-agents + synthesis · **Audit trail:** `tmp/review_2026-07-16_ng-admission-a1/`

**Scope files**

- `src/ng/region_typing/admission.rs` (new — the port)
- `src/ng/region_typing/mod.rs` (new — module scaffold)
- `src/ng/mod.rs` (modified — wires in the module)

**Design authority:** [spec](../../ng/spec/typed_regions.md) (Revision 2026-07-16, §4, §5, §5.1, §8.0),
[arch](../../ng/arch/typed_regions.md), [impl plan](../../ng/impl_plan/typed_regions.md) (Milestone A).

---

## Verdict: **approve with changes** — 1 Blocker, 9 Major, 21 Minor, ~28 Nits

The Blocker and every Major but one are **test-coverage findings, not logic defects**. That split is
the report's headline and it is worth stating plainly:

**The port is correct. The evidence that it is correct is weaker than it looks.**

Three sub-agents (`refactor_safety`, `smells`, `extras`) independently transcription-checked
`admit`/`finish_locus` and all seven helpers against `postprocess.rs:69-301`, and `prefilter` against
`scanner_parity.rs:48-82`, by mechanical diff rather than by eye. **All three found no transcription
defect.** The five declared shape divergences are the only ones; the gate order is preserved exactly
(period-1 drop before `drop_bundles`; copy floor after it — both load-bearing); the coordinate rebase
is confined to the single site the docs claim (`admission.rs:688-696`); and the `u32 → u64` widening
changes no result for any `u32`-representable value.

## The one thing to understand: why the differential is weaker than it appears

A1's whole safety argument rests on `differential_vs_production_build_loci_*` — ng's `admit` and
production's `build_loci` driven from the same inputs, asserted equal modulo the coordinate base.
`refactor_safety` confirmed the harness itself is sound: `assert_same_locus` compares ng against
production's *values* (not a re-derivation), and its `ref_tract`/`left_flank`/`right_flank` byte
comparisons cross-validate `tract_range()`, so a doubled or misplaced `+1` would fail it.

But `reliability` identified the structural limit, and it reframes every coverage finding below:

> Because ng and production share transcribed logic and the differential feeds both sides *the same
> values*, an untested input class isn't merely uncovered — it is **invisible**. Both sides would be
> wrong together and the test stays green.

A differential proves the copy matches the original. It cannot prove either is right, and **it goes
blind exactly where no input was supplied**. `extras` gave the same fact its sharpest form:
*deleting the `min_score` gate (`admission.rs:597`) leaves the entire suite green.*

So the fix set below is overwhelmingly "supply the missing inputs", not "change the code".

---

## Blocker

### B1 — `admission.rs:596,644,684` — case normalization has zero test coverage
**Categories:** reliability · **Confidence:** High

`admit`'s doc promises "any case — the tract, motif, and `ref_bytes` are upper-cased here for
case-stable identity", across three `upper()` call sites. **Every test contig in the suite is an
upper-case literal**, and the fixture was verified: the only lower-case bytes in
`tests/data/tandem_repeat/synthetic_ref.fa` are the two contig *names*. `upper()` is therefore an
identity function in 100% of the suite, including both differentials.

**Why it matters.** Real references **soft-mask repeats** — lower case is not an exotic edge case
here, it is the mainstream input class for exactly the sequence this module exists to process. The
failure mode is a lower-case `Motif` and lower-case `ref_bytes` reaching Milestone D: `Locus`
compares by value, so a soft-masked locus silently fails to match the catalog's. **Wrong locus, no
crash** — the precise failure this milestone is built to prevent.

**Fix:** a lower-case (and mixed-case) contig through `assert_agrees`, plus a direct assertion that a
soft-masked tract yields an upper-case motif and upper-case `ref_bytes`.

---

## Major

### M1 — `admission.rs:597` — the `min_score` gate is untested, and the differential is structurally blind to it
**Categories:** reliability, refactor_safety, extras (**convergent ×3**) · **Confidence:** High

Every test reaching `admit` runs `min_score: 0` against `score: 100`. The differential runs *both*
sides at `min_score = 0`, so the gate never rejects on either side. Deleting line 597 leaves the whole
suite green. Worse: A2's stated acceptance criterion ("differential still green at `Default`") is
*also* `min_score: 0`, so a mis-transcription survives A1 **and** A2 and first bites when an
experiment sets a real floor.

### M2 — `admission.rs:668` — the purity floor has no test at all
**Categories:** reliability, refactor_safety (**convergent ×2**) · **Confidence:** High

Every crafted tract is a perfect tiling (purity 1.0) against a floor of 0.8, so the gate never
rejects. Production's **one documented divergence from GangSTR** — "imperfect single-motif loci are
kept" (`postprocess.rs:20-21`) — is consequently not pinned at all.

### M3 — `admission.rs:688` — `minimal_trim` never moves the start in any crafted case
**Categories:** reliability · **Confidence:** High (crafted) / Medium (scanner)

All four surviving crafted cases were hand-traced: every one yields `st = 0, en = rep.len()`. With
`st == 0`, **deleting `+ st` from `new_start` is undetectable** — and the test named
`admit_emits_a_clean_locus_with_flanks_one_based`, which the code calls "the rebase's headline case",
is itself an `st = 0` case. The rebase's most delicate term is unpinned by the test that claims to
pin it.

**Fix:** a case with a genuine trim (`ATC` + `(AT)*6` → `st = 3`), driven through `assert_agrees`.

### M4 — `admission.rs:1419` — the differential is blind to `prefilter` by construction
**Categories:** reliability · **Confidence:** High

`differential_vs_production_build_loci_on_scanner_output` computes `cleaned = prefilter(&intervals)`
**once** and hands the same set to both sides. The pre-filter transcription is therefore never
compared against anything — a blind spot built into the harness. (Production's `catalog_prefilter` is
private inside a frozen `#[cfg(test)] mod`, so a direct differential is not available; the mitigation
is stronger direct unit tests.)

### M5 — `admission.rs:697` — `Locus::new(...).ok()` discards the port's only invariant check
**Categories:** errors, extras (**convergent ×2**) · **Confidence:** High

Faithful to production, but wrong *here*. Every `Err` arm was traced unreachable given the gates
above — which is exactly why discarding the verdict costs nothing on valid input and disarms the file
on invalid input. An impossible-coordinate locus currently exits through the same `None` as a routine
below-copy-floor rejection.

**Fix:** `debug_assert!(false, …)` + `None`. Release behaviour stays byte-identical (production not
diverged from, differential unaffected), but the differential — which runs in debug — fails loudly.

### M6 — `admission.rs:543` — `prefilter` underflows on a malformed interval
**Categories:** errors, reliability (**convergent ×2**) · **Confidence:** High

`iv.end - iv.start` on `RepeatInterval`'s **public** fields, which carry no ordering invariant. This
was harmless as a private `#[cfg(test)]` helper fed trusted input; A1 made it `pub`, changing the
threat model. Debug panic, or release wrap-to-huge — and a wrapped value sails through the copy floor
the filter exists to enforce. **Fix:** one-line predicate reorder; behaviour-preserving for
well-formed intervals.

### M7 — `admission.rs:69,457` — `MAX_MOTIF_LEN` and `MAX_PERIOD` must agree; nothing states or checks it
**Categories:** defaults, naming, smells (**convergent ×3**) · **Confidence:** High

Same value (6), same concept per the file's own docs, two unlinked names. A2 makes `MAX_PERIOD` a
knob: raise it to 8 and period-8 tracts pass the gate, `Motif::new` then rejects them, and
`finish_locus` silently drops them via `.ok()?`. **The knob appears to move and does not.**
**Fix:** `const _: () = assert!(MAX_PERIOD as usize <= MAX_MOTIF_LEN);` + the invariant in both docs.

### M8 — `admission.rs:410-419` — `SsrAdmissionParams::default()`'s stated rationale is unrealized
**Categories:** defaults, refactor_safety, extras (**convergent ×3**) · **Confidence:** High

`rg` shows **zero callers of `::default()` in the tree**. Both oracles hand-build `flank_bp: 5`; the
default is `50`. So "it exists so §8.0 compares like with like" describes something that does not
happen, and no test would fail if a field were changed to any value. **Fix:** pin it field-for-field
against `CatalogParams::default()` — production is frozen, so it can only break by a deliberate
ng-side edit, which is exactly the edit that should be loud.

### M9 — `admission.rs:410` — `min_score: 0`'s trap is documented in the one place nobody reads
**Categories:** defaults · **Confidence:** High

The prose is good — it names trf-mod's `-s 30`, ng's lack of it, and the Ruzzo–Tompa/TRF scale
mismatch. But it lives on `fn default()`, which rustdoc collapses. Someone writing
`SsrAdmissionParams { min_score: 0, .. }` by hand reads the **field's** doc, which says only "records
below are dropped" — equally true of `30` and of `0`. **Fix:** `DEFAULT_MIN_SCORE` const + the trap on
the field doc + a test pinning "at `Default`, even `i32::MIN` is admitted". Converts *documented* into
*announced*.

---

## Minor (selected; full set in the audit trail)

- **Mi1** — `admission.rs:523-527` vs `:435-437`: **the two copy-floor tables do not disagree.**
  *(defaults, extras — convergent ×2.)* They agree on every reachable period; the only numeric
  difference (period 1: 10 vs 3) is unreachable behind `prefilter`'s `iv.period >= 2` guard. They are
  **observationally identical** — the duplication is structural, not behavioural. The port is
  faithful; **the commentary overreaches**, and it inherits the overreach from spec §5/§10, which
  describe reconciling a disagreement that is not there. Fix the code doc *and* the spec.
- **Mi2** — `admission.rs:397`: `SsrAdmissionParams` documents a range (`min_purity ∈ [0,1]`) and a
  cross-field invariant (`bundle_threshold >= flank_bp`) and enforces neither. Sharp contrast with its
  sibling `Locus`, which rejects a `NaN` purity: the floor purity is *compared against* is unchecked,
  so a `NaN` `min_purity` silently disables the gate.
- **Mi3** — `admission.rs:1084-1168` vs `:1304-1384`: seven test fixtures written out twice. Both
  suites must exist (they assert different things), but the *inputs* are one thing said twice — edit
  one and the differential silently stops covering what the absolute test covers.
- **Mi4** — `admission.rs:583`: `admit` takes `Vec<RepeatInterval>` by value but only iterates it, and
  `RepeatInterval` is `Copy`. Production's `Vec<TrfRecord>` had a real reason (owned buffer, cloned in
  `drop_bundles`); `Candidate: Copy` removed it. The tests pay twice (`to_vec()`, `clone()`), and
  `prefilter` in the same file already takes `&[RepeatInterval]`.
- **Mi5** — `admission.rs:485`: `From<RepeatInterval>` reads fields instead of destructuring —
  matters because B2 widens that type and a destructure would force the error at the right place.
- **Mi6** — `is_close`'s strict `<` boundary: a `<=` slip passes the entire suite, differential
  included. *(reliability.)*
- **Mi7** — no property test, though `proptest` is already a dev-dependency and `assert_agrees` is a
  ready-made property body. *(reliability.)*
- **Mi8** — `Motif`/`MotifError`/`MAX_MOTIF_LEN` landed in the policy file **by default, not by
  decision** (spec §4 expected reuse, so no doc assigned the ported type a home). `Motif` is the one
  type in the file with zero admission content. Zero consumers today = cheapest moment to decide.
  *(module_structure.)*
- **Mi9** — the plan and spec §4 still say `Motif` is **reused**; the code ports it. The code is
  right; the authority docs need amending in the same commit. *(extras, module_structure.)*
- **Mi10** — SSR/STR prose drift (`:7`, `:68`, `:75`, `:127`, `:453`), including a user-facing error
  string. House rule: **STR** in prose, `ssr` in code. *(naming.)*
- **Mi11** — bare `2` at `:543` duplicates `MIN_PERIOD` and will not follow A2's knob. *(naming.)*
- **Mi12** — `admit`'s "must be pre-filtered" contract is carried only by prose, with a documented
  0/16 silent-cascade consequence. *(smells, reliability, naming — convergent ×3.)* A3 reshapes the
  signature anyway, so a `PrefilteredIntervals` newtype is cheap then.

## Nits

Grouped, not enumerated: the docs name a test (`differential_vs_production_build_loci`) that exists
only under two suffixed names; `Candidate`'s "widened to `u64`" is true of 2 of 4 fields; `Display for
Locus` sits ~450 lines from `impl Locus`; the scanner differential runs `admit` twice; the
frozen-production rationale is restated in full at all three module levels; `module_layout.md:38,150`
still draws step 3 as a file (doc stale, code right).

## What's good

- **The transcription itself** — three independent mechanical diffs, zero defects, on ~500 lines of
  coordinate-sensitive logic. The "transcribe, then change; never both in one commit" discipline is
  what made that checkable, and it paid.
- **The rebase is confined to one site**, exactly as spec §4 predicted, with all interior arithmetic
  left in production's 0-based space. That is why the port is diffable at all.
- **`assert_same_locus` states the conversion once** and cross-validates it through byte comparisons —
  it is not restating the same bug on both sides.
- **The boundary is clean.** After the `Motif` port, `region_typing/` names nothing from `src/ssr/`
  outside `#[cfg(test)]`. `git status` confirms production is untouched — the plan's hard constraint
  holds.
- **Params are passed explicitly** (`&SsrAdmissionParams`), so no hidden defaults, no `Option`
  fallback, all-`pub` + `Debug` = runtime-inspectable.
- **The idiomatic findings land almost entirely on ng's own lines, not the transcribed ones** — the
  signature of a port that stayed faithful.

## Out of scope observations

- `examples/ssr_psp_seqdump.rs:41` — pre-existing clippy `sort_by_key` error (`-D warnings` fails
  crate-wide because of it). Reproduces on a clean tree.
- `benches/psp_writer_perf.rs:386` — pre-existing panic under `cargo test --all-targets`. Reproduces
  on a clean tree.
- `cargo-mutants` is not installed. M1 is precisely the surviving mutant it would have found; worth
  considering for the port's later milestones.
