# ng typed-region generator — Milestone A complete (the admission port)

**Date:** 2026-07-16 · **Plan:** [typed_regions.md](../../ng/impl_plan/typed_regions.md) Milestone A
(A1 ✅ A2 ✅ A3 ✅) · **Commits:** `977b713` (design), `3e38e98` (A1), `0b2109b` (A2), A3
· **Status:** ⏸ **Checkpoint A** — the plan's hard pause.

## What Milestone A is

`src/ng/region_typing/admission.rs` — ng's **own copy** of the STR admission policy: the rules
deciding which detected tandem repeats become genotypeable STR loci. Ported from frozen production
(`ssr::catalog::postprocess::build_loci`, itself a port of GangSTR's `minimal_trim.py` /
`remove_bundles.py`) with the **logic transcribed unchanged** and five intended shape changes.

| step | what landed |
|---|---|
| **A1** | ng's `Motif` + `Locus` (1-based/`u64`), `SsrAdmissionParams`, `prefilter`, `admit` — the faithful transcription, and the differential that pins it |
| **A2** | every rule a knob: `periods`, `min_copies` (both floor tables folded into one), `bundle_threshold` collapsed into `flank_bp` |
| **A3** | windowing (`bases_start` + `contig_len`) and `Admitted { loci, bundled }` |

## Why a copy (owner, 2026-07-16)

> *"Don't touch `src/ssr`… copy the code to `ng/` and modify there… leave production as is and create
> a fresh ng caller from scratch."*

ng needs admission windowed, 1-based/`u64`, `RepeatInterval`-driven, all-knobs, and handing bundle
members back — five changes to one function, which is not "a small tweak", which is when the rule says
copy. **This milestone touched no file outside `src/ng/`**, verified at every commit.

The decision reversed the step-3 spec, which had argued the opposite at length (rebase
`ssr::types::Locus` across ~59 sites, widen `CatalogParams`, reshape `build_loci` in place). The
design docs were amended first, in `977b713`, before any code.

## How it is verified, and what that is worth

**Spec §8.0's differential.** `admit` and production's **real `build_loci`**, driven from the same
intervals, compared field-by-field modulo the coordinate base. Available without touching production
because `TrfRecord::for_test` is `#[cfg(test)] pub(crate)` in the same crate.

This is **stronger than the plan's original Milestone A**, which proved a production rebase by "the
golden fixture is still green" — one implementation against itself. This compares two implementations.
Three review agents independently mechanically-diffed the transcription and found **no defect**; a
fourth verified A3's degenerate-case algebra by hand for all inputs.

**And its limit, which shaped everything after A1.** Both sides run the same inputs through
transcribed logic, so an input class no test supplies is not merely uncovered — it is **invisible**:
both implementations are wrong together and the suite stays green. Three corollaries, each of which
bit at least once:

- `prefilter` is compared against **nothing**, by construction — `assert_agrees` computes it once and
  hands the result to both sides.
- The differential can only ever run **whole-contig**, so every windowed path is invisible to it.
- A gate whose floor is a no-op at the tested configuration is never exercised on either side.

**So every coverage claim in this milestone is mutation-verified, not asserted.** 20 mutants across
A1–A3; each was confirmed to survive the suite *before* the corresponding test existed.

## What the process actually caught

Worth recording, because it is the argument for the discipline rather than for the code:

- **A1** — `upper()` was an identity function in 100% of the suite (real references soft-mask repeats,
  so lower case is the *mainstream* input); the `+ st` trim term was deletable undetected, including
  by the test named "the rebase's headline case"; the purity floor and `min_score` gate had no test at
  all.
- **A2** — deleting `prefilter`'s period-floor gate left all 38 tests green, and that gate is
  load-bearing: **period 1 divides every period**, so a surviving homopolymer eliminates any real STR
  it overlaps (production's documented poly-A cascade). Then the review caught that **A2 itself
  re-committed A1's mistake**: collapsing `flank_bp` moved the bundle radius to 5, the fixture's 30 bp
  pair stopped bundling, and the test stayed green *because both sides moved together*.
- **A3** — clamping the left flank at `bases_start` instead of `1` survived: indistinguishable
  whenever a margin exists, differing only by silently truncating a flank instead of panicking.

**Four documented claims were disproved by the tests written to prove them:**

1. *"At `Default`, ng ships no score gate"* — it ships `score >= 0`. Ruzzo–Tompa emits only
   positive-scoring segments, so it is a no-op **for scanner output**, not an absent gate.
2. *"The two copy-floor tables disagree"* (spec §5/§10) — they agree on every reachable period. The
   duplication is **structural**, so A2 *deleted a copy* rather than resolving a conflict.
3. *"Reusing `Motif` costs production nothing"* (spec §4) — it is `pub(crate)`, so a `pub` ng `Locus`
   cannot return it. Ported.
4. *"Applying the period ceiling in `prefilter` would change which survivors eliminate which"* — it
   cannot; `floored` sorts ascending by period, so an out-of-scope interval can only ever be a victim.

All four are corrected in the spec, not only in the code.

## The error model, settled across A2–A3

One principle, arrived at the hard way: **a swept knob's guard must be a real `assert!`, because
sweeps run in release.** A debug-only guard lets an experiment record *"period 7 admits nothing"* when
the code never tried — a **wrong scientific result**, not a missing panic. `periods.max()`,
`contig_len`, `flank_bp >= 1`, and `min_purity` are all real asserts for that reason.

Everything else keeps production's model: a policy rejection is `Option::None`; `Locus::new` /
`Motif::new` failures are `debug_assert!(false, …) + None`, which stays release-identical to
production while failing loudly in the differential (which runs in debug).

**`cargo test --release` is a real check here**, and it found what no agent did: a debug-profile test
**cannot distinguish `assert!` from `debug_assert!`** — both fire. `region_typing` is now the one
module in the crate whose `should_panic` guards hold in release (crate-wide release failures 7 → 6;
the remainder are pre-existing debug-only tests, named for it).

## Deviations from the plan (recorded, not escalated)

- **`Motif` ported, not reused** (A1) — spec §4's premise was right and its conclusion wrong.
  `region_typing/` now names nothing from `src/ssr/` outside its `#[cfg(test)]` differential.
- **`region_typing/` is a folder**, where the arch specified a file — the port is a second concern with
  its own suite. Not a bake-off; `module_layout.md` principle 1 unbent.
- **`Candidate`** (internal) stands in for `TrfRecord`, so the `RepeatInterval` widening happens once
  at the entry.
- **`MinCopies`**, not the spec's implied "copy floor" vocabulary — ng's scanner already says
  `min_copies`, and **"copy number" means CNV in this crate** (`paralog/`), so the obvious name pointed
  the geneticist reader at the wrong concept.

## Carried into later milestones (from the reviews, with homes)

- **C1** — a `Window`/newtype for `admit`'s three co-dependent params. Load-bearing, not cosmetic: the
  `contig_len` guard is **inherently incomplete** (`<=` admits a caller passing the *window's own end*
  — production's exact mistake — and no arithmetic can catch that; only **provenance**, a contig table,
  can). Spec §5a spells the flat signature, so this is a design change, not a cleanup.
- **D1** — assert the partition on the scanner-output fixture; it currently projects `.loci` and throws
  `bundled` away, so the only realistic-density fixture asserts nothing about the whole partition.
- **Own step** — a `proptest` property over `assert_agrees`. `proptest` is already a dev-dependency and
  `assert_agrees` is a ready-made property body; the strongest remaining lever on the port.
- **A2/A3 successors** — `admit` takes `Vec` by value; a `PrefilteredIntervals` newtype; test-fixture
  case tables.

## Validation

```
cargo fmt -- --check                 clean
cargo clippy --lib --tests           clean (0 in scope)
cargo test --lib ng::region_typing   54 passed; 0 failed
cargo test --release --lib ng::region_typing   54 passed; 0 failed
cargo test --lib                     1820 passed; 0 failed; 4 ignored
cargo doc --no-deps                  0 errors in region_typing
git status                           nothing outside src/ng/
```

Pre-existing and out of scope (all reproduce on a stashed clean tree):
`examples/ssr_psp_seqdump.rs:41` clippy error; `benches/psp_writer_perf.rs:386` panic; 8 `cargo doc`
unresolved links; 6 release-only `should_panic` failures over `debug_assert`s.

## Next — Milestone B (after Checkpoint A)

**B1** promote and stream `collect_windowed`; **B2** widen ng's own coordinates to `u64` and delete
`ref_seq.rs`'s `unwrap_or(u32::MAX)` silent clamp; **B3** `WindowedRefSeq` raw bytes + `contigs()`.

**One question B2 must answer that the docs leave open:** spec §4 says ng is 1-based *everywhere*, and
B2's plan text says only "widen". `RepeatInterval` is 0-based today; it is both a slice-offset type
(where 0-based is natural) and the payload of `RegionKind::SsrBundle` (where genomic 1-based is what
the walk emits). A3 side-stepped this by keeping `Admitted::bundled` in the input's space and leaving
the re-basing to the walk. **B2 should decide it explicitly rather than inherit it.**
