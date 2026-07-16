# Code review + fixes — ng admission port, Milestone A2 (the knobs)

**Date:** 2026-07-16 · **Scope:** working-tree diff vs `3e38e98` (A1) · **Reviewers:** 4 triaged
category sub-agents + synthesis · **Audit trail:** `tmp/review_2026-07-16_ng-admission-a2/`

**Category triage.** `refactor_safety` (A2 claims behaviour-preservation), `defaults` (the change *is*
configuration), `reliability`, `naming`+`idiomatic`. Skipped with reason: `module_structure` (no
module change), `unsafe_concurrency` (none present), `tooling` (no manifest change), `errors` (no
error-type change), `smells` (the duplication question was settled at A1 and is owner policy).

**Verdict: approve with changes — 0 Blockers, 6 Major, 9 Minor, ~8 Nits.** All applied.

---

## The headline

**A2 changed the shape, not the behaviour** — `refactor_safety` verified each hoist against the
*frozen originals* rather than against A2's claims: the `u16::from(u8)` gate was lossless, the folded
table reproduced both originals including the `_ => 3` arm, the `bundle_threshold` collapse is
lossless at `Default`, and `r.period as u8` was safe in release because it rested on the gate, not on
the `debug_assert`.

**But two of A2's own changes quietly weakened the evidence, and both were caught.** That is the
pattern of this review: A1's lesson — *a differential goes blind where both sides move together* —
recurred, and A2 walked into it.

## Major

### M1 · The `flank_bp` collapse silently disarmed the bundle-drop fixture
**Categories:** refactor_safety, naming+idiomatic (**convergent ×2**) · **Applied**

A1's differential ran `flank_bp: 5, bundle_threshold: 50`. Collapsing them makes the bundle radius
**5** — and that fixture's closest pair is **30 bp apart**, built to sit inside 50. At radius 5 no
`is_close` clause fires: nothing bundles, the case named *"bundle drop keeps only the isolated tract"*
admits all three, and **the test stays green because both sides moved together**. Step 3 of the policy
stopped being differentially compared, silently.

*This was reasoned about during implementation and then dropped* — the note "use a 50 bp radius for
that case" was written and never applied.

**Fix:** the bundle case now runs `assert_agrees_at` at a 50 bp radius, plus a new
`admit_bundles_at_the_flank_radius` pinning the collapse's live wire directly (30 bp apart: bundled at
50, both admitted at 5, production agreeing at both). Mutation-verified: `drop_bundles(kept, 50)`
instead of `p.flank_bp`, and deleting the bundle drop entirely, both now fail.

### M2 · `debug_assert` was the wrong guard for the period ceiling
**Categories:** refactor_safety, reliability (**convergent ×2**) · **Applied**

`PeriodRange::new` bounds neither end (verified in source), so `PeriodRange::new(2, 7)` is legal and
the path is reachable. In release, period-7 tracts clear the gate and vanish at `Motif::new` through
the same `None` as a policy rejection.

The argument that settles it is one the implementation did not have: **this knob exists to be swept by
spec §5.2's routing experiment, and sweeps run in release.** A debug-only guard lets that experiment
record *"period 7 admits nothing"* when the code never tried — **a wrong scientific result, not a
missing panic**. The "unvalidated like production's `CatalogParams`" defence does not transfer:
production has no period knob.

**Fix:** a real `assert!` (once per `admit` call, free next to the scan it guards), plus
`#[should_panic]` coverage.

### M3 · The per-period minimum was never tested *per period*
**Categories:** reliability · **Applied**

`raising_the_copy_floor_changes_what_is_admitted` uses `uniform(9)`, which by construction makes every
entry identical — so it pins *that* the table is read, never *which entry*. A `for_period` returning
`by_period[1]` regardless of its argument passed all six callers. **Non-uniform is the only shape
§5.2's frontier needs.**

**Fix:** `the_minimum_is_looked_up_by_the_tracts_own_period` — a period-2 and a period-3 tract under a
table that admits one and rejects the other, then mirrored so the surviving locus swaps, through both
`admit` and `prefilter`. Mutation-verified.

### M4 · The `Default` period pin compared a constant with itself
**Categories:** defaults · **Applied**

`default_matches_the_frozen_catalog_params` asserted `ours.periods.min() == DEFAULT_MIN_PERIOD`, but
`Default` is *defined* as `PeriodRange::new(DEFAULT_MIN_PERIOD, DEFAULT_MAX_PERIOD)` — a tautology
that cannot fail. Production's `MIN_PERIOD`/`MAX_PERIOD` are private consts ng cannot read, so
`DEFAULT_MAX_PERIOD` was the one A2 knob free to drift while the test built to catch it stayed green,
silently decalibrating both §8 oracles.

*(The implementation's mutation run did kill `ceiling 6 → 5` — but via the **differential**, which
catches it only because the fixture happens to contain period-6 tracts. Incidental coverage, not a
pin. The reviewer's finding is better than the implementer's evidence.)*

**Fix:** restate production's consts as literals, the pattern
`the_folded_copy_floor_table_reproduces_both_of_productions` already uses for the minimums.

### M5 · No `#[should_panic]` tests existed, so every new guard was free to delete
**Categories:** reliability · **Applied**

Zero in the file. **Fix:** `admit_rejects_a_period_ceiling_no_motif_can_hold` and
`admit_rejects_a_nan_purity_floor`. Noted and accepted: `finish_locus`'s `Motif::new` guard is
unreachable through `admit` in a debug build (the ceiling assert fires first), so it is a
belt-and-braces guard for direct callers, not independently testable through the public path.

### M6 · The `prefilter`-ceiling rationale was false in both halves
**Categories:** refactor_safety, reliability (**convergent ×2**) · **Applied**

The doc claimed applying the period ceiling in `prefilter` "would change which survivors eliminate
which" — **it cannot**: `floored` sorts ascending by period, so an out-of-scope interval only ever
sorts *after* in-scope ones and can only be a victim, never an eliminator; everything it could
eliminate is its own multiple, which `admit`'s ceiling drops anyway. And it claimed "the differential
would catch it otherwise" — **it cannot**, by construction, contradicting the module header two
screens up.

**The decision is right; the reasons were wrong.** Fix: the real reason (production's pre-filter gates
only `>= 2`, and no test can compare `prefilter`, so matching the original by inspection *is* the
evidence) now stands alone, with both refuted arguments recorded so they are not re-invented.

## Minor (applied)

- **`CopyNumberFloors` → `MinCopies`** *(naming)*. It coined a third word for a quantity ng already
  names `min_copies` (`tandem_repeat.rs`), whose own doc calls itself "a permissive floor" and defers
  to exactly this. Worse: **"copy number" means CNV in this crate** (`relative_copy_number`,
  `carrier_copy_numbers` in `paralog/`), so for the geneticist reader it pointed at the wrong concept.
  Free to fix — the differential pins values, not identifiers. With it: `.floor()` → `.for_period()`,
  `beyond_table` → `for_wider_periods`, the field → `min_copies`.
- **`Candidate.period: u16` → `u8`** *(idiomatic)*. A widening with no source and no sink —
  `RepeatInterval::period` is `u8`, `PeriodRange` is `u8`, the lookup takes `u8`. A1 chose `u16` to
  mirror `TrfRecord`, but `Candidate` is not a `TrfRecord`. Removes two `u16::from`s and an unchecked
  `as u8`; **port cost is negative**.
- **The dead index 0** *(defaults, naming, reliability — ×3)*. `by_period` is now indexed `period - 1`,
  so the literal is exactly GangSTR's dict values with no slot to explain, and `for_period(0)` falls to
  `for_wider_periods` — which is what **both** originals' `_ => 3` arms return for it. The total
  function now agrees with production everywhere, not just on the reachable range.
- **`Default`'s `expect` cited an assert that does not exist** *(naming+idiomatic, reliability — ×2)*.
  The const asserts prove `MAX <= MAX_MOTIF_LEN` and `MIN <= MAX`; they never proved
  `1 <= DEFAULT_MIN_PERIOD`, which is the first thing `PeriodRange::new` rejects. Harmless at `MIN = 2`
  — but §10's experiment is precisely about moving that floor toward 1. Added.
- **`for_wider_periods` was indistinguishable from the last named period** *(reliability, defaults)* —
  `Default` and `uniform` both set them equal, so a `for_period` ignoring it passed everything. Added
  `min_copies_uses_the_fallback_only_beyond_the_table`. Mutation-verified.
- **`prefilter(…, p: &SsrAdmissionParams)`** *(naming)* — `p` is the letter this module's docs use for
  *period*, one screen from `p.min_copies.for_period(iv.period)`. Renamed `params`.

## Deferred

| finding | why | home |
|---|---|---|
| `matched_params`'s `..default()` hides a new *ng* field while a new *production* field breaks the build | Real asymmetry, but `default_matches_the_frozen_catalog_params` pins every defaulted field, and A3 reshapes this helper anyway. | A3 |
| `SsrAdmissionParams` has `pub` unvalidated fields where `periods` is a checked type — an API asymmetry | A genuine design question (checked constructor vs `assert` in `admit`) and a `TypedRegionConfig`-level one, since C2 wraps this type. | C2 |
| `r.period as u8` is a B2 trap once `RepeatInterval` widens | Now moot — `Candidate.period` is `u8`; if B2 widens `RepeatInterval::period`, `From` breaks at compile time, which is the right place. | — (closed) |
| no `proptest` property over `assert_agrees` | Still the strongest remaining lever on the port. | own step |

## Validation

```
cargo fmt -- --check                clean
cargo clippy --lib --tests          clean (0 in scope)
cargo test --lib ng::region_typing  44 passed; 0 failed      (was 39)
cargo test --lib                    1810 passed; 0 failed; 4 ignored
cargo doc --no-deps                 0 errors in region_typing
```

**Mutation-verified after the fixes** — each of these previously survived:

| mutant | result |
|---|---|
| `for_period` ignores its argument | killed (3 tests) |
| `drop_bundles` no longer reads `flank_bp` | killed |
| the bundle drop deleted entirely | killed (2 tests) |
| the period-ceiling `assert` deleted | killed |
| `for_wider_periods` ignored | killed |

**Production untouched** — `git status` shows nothing outside `src/ng/` and the docs/reports.
