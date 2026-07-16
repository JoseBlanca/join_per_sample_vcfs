# Fixes applied — ng admission port, Milestone A1

**Date:** 2026-07-16 · **Review:** [ng_admission_a1_2026-07-16.md](ng_admission_a1_2026-07-16.md)
(1 Blocker, 9 Major, 21 Minor, ~28 Nits from 9 parallel category agents)

**Applied:** the Blocker, 9 of 9 Majors, 7 Minors, 4 Nits. **Deferred:** the rest, each with a reason
and a home. **No `Ask`-class findings** — nothing needed an owner decision.

---

## The shape of this fix set

The review's verdict was *"the port is correct; the evidence that it is correct is weaker than it
looks"* — three independent mechanical diffs found **no transcription defect**, while the Blocker and
7 of 9 Majors were **missing test inputs**. So this is overwhelmingly a test commit. The code changes
are four small hardening edits; everything else supplies inputs the differential was never driven
with.

**Every coverage fix was verified by mutation, not by assertion.** "The suite stays green" is the
exact failure being repaired, so claiming a test closes a gap without breaking the code to prove it
would repeat the mistake:

| mutation | before | after |
|---|---|---|
| delete the `min_score` gate | **survives** | killed |
| delete `+ st` from `new_start` (the trim term) | **survives** | killed |
| make `upper()` the identity | **survives** | killed |
| delete the purity gate | **survives** | killed |
| `is_close` clause 1 (`start`–`start`) `<` → `<=` | **survives** | killed |
| `is_close` clause 3 (gap) `<` → `<=` | **survives** | killed |
| `is_close` clause 4 (`end`–`end`) `<` → `<=` | **survives** | killed |
| `is_close` clause 2 (`start`–`end`) `<` → `<=` | survives | **survives — provably equivalent** (below) |
| double the rebase `+1` | killed | killed |

## Two claims the fixes *disproved* rather than implemented

Both were mine, and both are now corrected in the code **and** in the spec they came from.

1. **"At `Default`, ng ships no score gate at all"** (spec §5c, and the `min_score` docs). The test
   written to prove it **failed**: the gate is `score >= 0`, so `i32::MIN` *is* rejected. The true
   claim is narrower — Ruzzo–Tompa emits only positive-scoring segments
   ([`tandem_repeat.rs`](../../../../src/ng/tandem_repeat.rs): `score = r - l > 0`), so the floor is
   unreachable **for scanner output**. A no-op in practice, not an absent gate. Test renamed
   `admit_at_default_never_rejects_a_score_the_scanner_can_emit` and asserts both the reachable range
   and the literal boundary.
2. **"The two copy-floor tables disagree"** (spec §5/§10, inherited into the code docs; independently
   flagged by `defaults` and `extras`). They do not. They give the same floor for **every reachable
   period**; their one numeric difference is period 1 (10 vs 3), which neither gate can reach
   (`prefilter` requires `period >= 2`, `MIN_PERIOD` drops it again). The duplication is
   **structural, not behavioural** — so A2 deletes a copy rather than resolving a conflict, which is
   a smaller and better-understood job than the spec had been describing. Pinned by
   `the_two_copy_floor_tables_agree_on_every_reachable_period`.

## Applied — Blocker

**B1 · case normalization had zero coverage.** `upper()` was an identity function in 100% of the
suite; the fixture's only lower-case bytes are its contig *names*. Added
`admit_upper_cases_a_soft_masked_tract_and_its_flanks`: a fully soft-masked contig and a mixed-case
one (unique flanks upper, repeat lower — real soft-masking's actual shape), each asserted to produce
an upper-case motif/tract/flanks, to compare equal to the all-upper locus, and to agree with
production through `assert_agrees`. Mutation-verified.

## Applied — Major

- **M1 · `min_score` untested, differential blind** *(×3 convergent)* —
  `admit_applies_the_score_gate_when_the_floor_actually_bites` drives both sides at a floor of 50
  with scores 100 / 49 / 50, pinning above, below, and the `>=` boundary. This required
  `assert_agrees_at`, so the differential can run at a configuration where a gate actually fires (see
  "harness" below).
- **M2 · purity floor had no test** *(×2)* —
  `admit_keeps_an_imperfect_tract_above_the_floor_and_drops_one_below`: an interrupted tract
  (purity 0.9375) kept at a 0.8 floor, dropped at 0.99, both through the differential. This also
  pins production's **one documented divergence from GangSTR** ("imperfect single-motif loci are
  kept"), which nothing covered.
- **M3 · `minimal_trim` never moved the start** — every crafted case had `st = 0`, so `+ st` was
  deletable undetected, *including by the test named "the rebase's headline case"*.
  `admit_reports_the_trimmed_start_not_the_detected_one` uses an `ATC`-prefixed tract (`st = 3`) and
  asserts the trimmed start (9), not the detected one (6).
- **M4 · differential blind to `prefilter` by construction** — unavoidable for the differential
  itself (production's `build_loci` needs pre-filtered input too, and `catalog_prefilter` is private
  inside a frozen `#[cfg(test)] mod`, so no direct comparison exists). Mitigated with direct unit
  tests instead: the inverted-interval case (M6) plus the existing multiple/noise cases. **Residual
  recorded** — `prefilter`'s transcription rests on inspection + its own units, not on a differential.
- **M5 · `Locus::new(...).ok()` discarded the port's only invariant check** *(×2)* — now an explicit
  `match`: `Ok` → `Some`, `Err` → `debug_assert!(false, …)` + `None`. Release behaviour is
  byte-identical to production (so the port stays faithful and the differential is unaffected), but
  the `Err` arm is now only reachable if **the transcription is wrong** — and the differential runs
  in debug, so it fails loudly instead of quietly disagreeing.
- **M6 · `prefilter` underflowed on a malformed interval** *(×2)* — `iv.end > iv.start` now guards
  the subtraction, ordered first in the predicate. Behaviour-preserving for well-formed intervals.
  Test: `prefilter_drops_an_inverted_interval_instead_of_underflowing`. The threat model changed when
  A1 made this `pub`; the docs now say so.
- **M7 · `MAX_PERIOD`/`MAX_MOTIF_LEN` invariant unstated** *(×3)* — now
  `const _: () = assert!(MAX_PERIOD as usize <= MAX_MOTIF_LEN)` with the silent failure mode spelled
  out (a wider tract clears the gate, `Motif::new` rejects it, `finish_locus` drops it through the
  same `.ok()?` as a policy rejection — the knob appears to move and does not). Plus
  `MIN_PERIOD <= MAX_PERIOD`. Compile-time, so A2 cannot introduce it by accident.
- **M8 · `Default`'s rationale unrealized** *(×3)* — it had **zero callers**.
  `default_matches_the_frozen_catalog_params` now pins every field against `CatalogParams::default()`
  (production is frozen → can only break by a deliberate ng-side edit), plus
  `bundle_threshold == flank_bp` (spec §2.4).
- **M9 · `min_score`'s trap documented where nobody reads it** — extracted `DEFAULT_MIN_PURITY` /
  `DEFAULT_MIN_SCORE` / `DEFAULT_FLANK_BP` consts, moved the trap onto the **field** doc (the surface
  actually read when hand-building the struct), and pinned it with a test. Corrected en route — see
  "disproved" above.

## Applied — Minor / Nit

- **Mi1** — the "disagreeing tables" correction, in code **and** spec §5/§10 (see above).
- **Mi2** — `SsrAdmissionParams`' documented-but-unenforced contracts: `admit` now `debug_assert`s
  `min_purity` finite in `[0, 1]` and `bundle_threshold >= flank_bp`. The `NaN` case is the one that
  bites — every `purity < p.min_purity` is then false, so the gate silently passes everything.
  (Chose `debug_assert` in `admit` over a validating constructor: production's `CatalogParams` is
  also `pub`-field/unvalidated, and A2 reshapes this type anyway.)
- **Mi6** — `is_close`'s strict `<`, and **more than the review asked**: the first fixture pinned only
  the *gap* clause, and mutation showed 3 of 4 clauses still survived. Added an overlapping
  equal-length pair that puts three clauses on the boundary simultaneously. The fourth
  (`a.start`/`b.end`) is **provably an equivalent mutant** — its boundary requires
  `b.end == a.start + thresh`, which forces `b.start - a.start < thresh`, so clause 1 has already
  fired; callers hold `a.start <= b.start` because `admit` sorts first. Documented in the test rather
  than left as an apparent hole.
- **Mi9** — spec §4, the arch recon table + decisions, and the plan all said `Motif` is **reused**;
  the code ports it. All three corrected, with the general lesson recorded: *"costs production
  nothing" must be checked against **visibility**, not just coordinates and width.*
- **Mi10** — STR/`ssr` prose drift fixed at 8 sites, including the user-facing `MotifError` string.
- **Nits** — the module doc named a test that never existed under that name; the differential's blind
  spot is now stated in the module docs (it is not obvious, and it is why the new cases exist).

### Harness changes (enabling the above)

- **`matched_params(...)` replaces two hand-maintained literal tables.** `refactor_safety` flagged
  `params()`/`prod_params()` as duplicated by hand; that is exactly how the two sides would quietly
  stop running the same settings while the test kept passing. Production's params are now *derived*
  from ng's, so they cannot drift.
- **`assert_agrees_at(ng_p, prod_p, …)`** lets the differential run at a configuration where a gate
  fires; `assert_agrees` stays the catalog-settings default.

## Deferred (with reasons)

| finding | why deferred | home |
|---|---|---|
| Mi3 — seven test fixtures written out twice | Real, but a shared case table is a test-structure refactor best done once A3 settles the signature; the new cases are single-use and not duplicated. | A3 |
| Mi4 — `admit` takes `Vec` by value though `Candidate: Copy` | Spec §5a and arch both spell the signature `admit(recs: Vec<RepeatInterval>, …)`; **A3 reshapes it anyway** (adds `bases_start`/`contig_len`, returns `Admitted`). Changing it now diverges from the doc for one milestone. | A3 |
| Mi5 — `From<RepeatInterval>` reads fields instead of destructuring | Worth doing exactly when B2 widens `RepeatInterval` — a destructure forces the compile error at the right place then. | B2 |
| Mi7 — no property test (`proptest` is already a dev-dep; `assert_agrees` is a ready-made property body) | Genuinely valuable and the single strongest remaining lever on the port. Needs a sequence+interval generator, which is more than a fix-application. **Recommended as its own step.** | follow-up |
| Mi8 — `Motif`'s module home | Landed in the policy file *by default, not by decision*. C1/C2 place ng's shared vocabulary (`GenomeRegion`, `Position`) and are the natural moment to decide; still zero consumers. | C |
| Mi11 — bare `2` at `prefilter` duplicating `MIN_PERIOD` | **Partly applied** — now `MIN_PERIOD as u8`. The nested `copy_floor` table's own literals follow A2's knob. | A2 |
| Mi12 — `admit`'s "must be pre-filtered" contract is prose only | A `PrefilteredIntervals` newtype is cheap **when A3 reshapes the signature**; doing it now churns it twice. | A3 |
| naming — renaming inherited names (`recs`, `copy_number_floor` vs `copy_floor`, …) | The `naming` agent's own framing: inherited names defer to A2, because renaming fights the line-comparable transcription the differential depends on. Port-introduced names were reviewed and kept. | A2 |
| idiomatic — `pub` vs `pub(crate)` across ng | Tree-wide convention, not an A1 defect (the agent filed it Nit for that reason). | owner |
| idiomatic — index loops in `drop_bundles`/`count_motif`/`minimal_trim` | **Not filed by the agent, and correctly so** — verbatim production; the "idiomatic" forms are not always equivalent (`minimal_trim`'s in-loop `return None` cannot be a `find`). | — |

## Validation (container, `./scripts/dev.sh`)

```
cargo fmt -- --check                clean
cargo clippy --lib --tests          clean (0 warnings in scope)
cargo test --lib ng::region_typing  36 passed; 0 failed      (was 27)
cargo test --lib                    1802 passed; 0 failed; 4 ignored
cargo test --tests                  all suites pass
cargo doc --no-deps                 0 errors in region_typing
```

**Production untouched** — `git status` shows only `src/ng/`, the three design docs, and these
reports.

**Pre-existing, out of scope, verified on a stashed clean tree:** `examples/ssr_psp_seqdump.rs:41`
clippy `sort_by_key` (fails `-D warnings` crate-wide); `benches/psp_writer_perf.rs:386` panic under
`--all-targets`; 8 `cargo doc` unresolved-link errors elsewhere in the crate.
