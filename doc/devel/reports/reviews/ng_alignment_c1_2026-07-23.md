# Code Review: ng_alignment_c1
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step C1 of `doc/devel/ng/impl_plan/alignment_best_path.md` — band the delimiter
**Status:** Approve-with-changes (all applied — see §Author response)

---

### 1. Scope

The uncommitted diff for C1 (HEAD `c32282f`): one file,
[ssr_best_path_flat_gap.rs](../../../../src/ng/alignment/ssr_best_path_flat_gap.rs) — a `BAND_HEADROOM`
constant, a `band_width` function, the banding logic inside `delimit` (in-band cells computed,
out-of-band written `UNREACHABLE`), an updated `ViterbiScratch` doc, and new tests.

Two reviewers: one **attacking the band's correctness by experiment** (this is a silent-failure step),
one on naming/idiomatic/smells/refactor_safety.

### 2. Verdict

**Approve-with-changes.** No Blockers, no Majors. **The band survived every attack.**

### 3. The band is correct — established by experiment, not reading

The reliability reviewer narrowed each term of `band_width` in the container and confirmed the parity
oracle catches it, then tried to defeat the formula and could not:

| experiment | result |
|---|---|
| drop the `left + right` run-off terms (the spec's original model) | fails at seed `0x5eed0001` **case 1745** — a silent 4-base-short measurement `(0,38)` vs production `(0,42)` |
| drop `BAND_HEADROOM` entirely | fails at seed `0xc0ffee42` case 288 — the `+8` is load-bearing |
| halve the floor term | subtract-with-overflow panic — corner-reachability breaks (the benign failure mode) |
| **`BAND_HEADROOM = 1`** | **passes the 12k default, fails the 200k soak at case 28,307** |
| `BAND_HEADROOM = 2` | passes the soak |

So the empirical minimum headroom is **2**, and the shipped **8** carries a deliberate 4× margin —
"generous" is accurate, not hand-waving. Two structural conclusions the reviewer reached and I concur
with:

- **`left + right` is provably sufficient**, not merely enough for this fixture. The stray past the
  `|read − reference|` corridor decomposes into the net-length difference (the floor), a forced
  run-off bounded by *each flank's length*, and self-limiting net-zero bows (the constant). Because
  the run-off term scales *with* the flanks, it cannot be beaten by enlarging them — no counterexample
  exists.
- **Scratch staleness is neutralised.** Every column is written every row (out-of-band as
  `UNREACHABLE`, not skipped), so all three neighbour reads are fresh, including the low-band-edge
  `current[column-1]` and high-edge `previous[column]`. Out-of-band backpointers are stale but never
  read: the final cell `(m, n)` is always in-band and finite, and a finite cell has only finite
  (in-band) predecessors, so the traceback cannot step out of band.

### 4. Findings

#### Minor (applied)

**C1-1: the parity fixture's *default* cannot see a one-cell-too-narrow band.** A gross deficiency
(dropping the run-off terms) fails on the first seed inside 3,000 cases, but a band one cell short of
correct first diverges near case 28,307 — past the 12k default, needing the soak. Not a defect in the
band (which has a 4× margin), but a real property of the oracle that a future narrowing could trip
over silently.
**Fix applied:** documented on `cases_per_seed` in the parity harness — the default's detection floor
is coarser for banding than for the recurrence, and any change narrowing the band toward its true
minimum must be validated at soak size.

**C1-2 (from the style pass): the spec-departure deferral lived only in prose.** Flagging in code that
the band formula departs from spec §3/§9 (the run-off correction is a finding the spec did not model)
is a strength — but "recorded for the owner" with no grep-able marker risks being lost.
**Fix applied:** a `SPEC-FOLLOWUP(alignment §3/§9)` marker in the `BAND_HEADROOM` doc, cross-referenced
to the impl report and `PROJECT_STATUS.md`.

#### Nits (applied or accepted)

- **`band_width`'s four positional `usize`s** are a nominal transposition hazard, but harmless by
  construction: the formula is symmetric under both the read/reference and left/right swaps. **Fix:**
  a sentence on the function saying exactly that, so the bare primitives read as deliberate.
- One concept wore three names (`band_width` / `band` / "half-width"). **Fix:** dropped "half-width".
- The per-cell `in_band` `abs_diff` over a precomputed `[lo, hi]` range: **kept** — the band must write
  every out-of-band cell, so no cells are skipped and a range buys nothing; it also mirrors
  production's `if i.abs_diff(j) > band { continue; }`.
- The triple `[UNREACHABLE; 3]` write: **kept inlined** — a helper would obscure three genuinely
  distinct sites (row 0, column 0, inner loop).
- The long `BAND_HEADROOM` doc **earns its length** per the file's convention for subtle constants
  (both reviewers).

### 5. Missing tests added now

The two the reliability reviewer proposed:
- `a_read_off_both_ends_of_an_expanded_tract_is_measured` — a pure-tract read that deletes *both*
  flanks, stressing the band on both sides at once.
- `scratch_survives_an_extreme_size_drop_under_banding` — a huge read then a tiny one on the same
  scratch, pinning that out-of-band cells are written rather than skipped.

Plus the step's own required pair: `a_long_allele_at_the_extreme_is_measured_not_collapsed` (the
plan's long-allele extreme), `a_run_off_flank_read_is_measured_despite_the_bow` (the case-1745
finding as a fixed-input regression), `band_width_is_geometry_not_the_slip_cutoff`, and
`the_band_always_reaches_the_far_corner`.

### 6. What's good

- The band's three terms are each justified and each *proven load-bearing by a failing mutation* —
  there is no dead width.
- The doc states plainly where the formula departs from the spec, rather than silently contradicting
  the design doc — the honest handling of a spec model that turned out incomplete.
- Writing out-of-band cells `UNREACHABLE` rather than skipping them is the choice that makes banding
  safe on a reused scratch, and it is documented as such.

### 7. Commands to re-verify

```
./scripts/dev.sh cargo test --lib
./scripts/dev.sh env PVC_PARITY_CASES=50000 cargo test --release --lib ng::alignment::delimit_parity
```

### Author response

Both Minors and every actioned Nit **fixed in this step's commit**; the accepted Nits are recorded
above with their reasons. The band's correctness rests on the parity oracle, which the review
exercised directly rather than trusting.

Per-category audit trail at `tmp/review_2026-07-23_ng-alignment-c1/`.
