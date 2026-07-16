# Code review + fixes — ng admission port, Milestone A3 (windowing + `Admitted`)

**Date:** 2026-07-16 · **Scope:** working-tree diff vs `0b2109b` (A2) · **Reviewers:** 3 triaged
category sub-agents (5 checklists) + synthesis · **Audit trail:**
`tmp/review_2026-07-16_ng-admission-a3/`

**Category triage.** `refactor_safety` (A3 claims degenerate-case behaviour-preservation),
`reliability`, and `errors`+`naming`+`idiomatic` (A3 introduces the module's first real panics).
Skipped with reason: `defaults` (A2's work, no new knob), `module_structure` (no module change),
`unsafe_concurrency`, `tooling`, `smells` (settled at A1; owner policy).

**Verdict: approve with changes — 0 Blockers, 7 Major, 12 Minor, ~9 Nits.** All Majors applied.

---

## The headline

**The degenerate case reduces to production exactly, verified by hand rather than by argument.**
`refactor_safety` walked the algebra for all inputs, including the cases the implementer was least
sure of:

- `new_start == flank - 1`: `saturating_sub` gives 0, `.max(1)` gives 1, production gives 0 → `prod
  + 1` ✓. `flank == 0` drops on both sides.
- `new_end == 0` is **unreachable** (`minimal_trim` guarantees `st < en`).
- Empty `bases` does **not** underflow the entry guard — `assert!(bases_start >= 1)` runs first and
  `+` binds before `-`.
- `split_bundles` is production's selection **member for member**; `bundled` is coordinate-ordered and
  correctly excludes other-reason rejects (spec §2.2's `Generic` territory).
- `From<Candidate>`'s `expect` is sound — confirmed by grep that `Candidate::from` is the only
  constructor.

**Every Major is about the error model, and they all trace to one inconsistency.** A2 established the
principle in code — `periods.max()` is a real `assert!` *because the knob is swept and sweeps run in
release, so a debug-only guard would let an experiment record a wrong result*. A3 wrote that
reasoning down and then failed to apply it to the two guards on either side of it.

## Major

### M1 · The `contig_len` guard was a `debug_assert`, and the module's own argument says it must not be
**Categories:** refactor_safety, reliability, errors (**convergent ×3**) · **Applied**

`reliability` traced the release path concretely: window `[901, 1076]`, tract `[951, 966]`,
`contig_len = 1000` → `ref_end` clamps early, the `ref_end == tract_end` test still passes, and
`admit` **emits a locus with a 34 bp right flank instead of 50, silently**. That is spec §2.6's bug
reappearing *inside the code written to kill it*, and §10's flank and routing sweeps drive this path
in release. Now an `assert!`.

**And it is inherently incomplete** (`errors`) — `<=` admits equality, so a caller passing *the
window's own end* as `contig_len` (production's exact mistake) is arithmetically legal and **no check
here can catch it**. That half needs **provenance**, not arithmetic: a `contig_len` read from the
reference's contig table rather than derived from the slice in hand. Recorded in the code and routed
to Milestone C, where the walk acquires `contigs()`.

### M2 · `min_purity`'s guard was debug-only, and its `#[should_panic]` test failed in release
**Categories:** (found by the fixes' own release run) · **Applied**

Not caught by any agent — surfaced by running `cargo test --release` to *verify M1's fix*, which is
the only way to tell an `assert!` from a `debug_assert!` (a debug-profile test cannot: both fire).
`admit_rejects_a_nan_purity_floor` asserted a panic that only existed in debug. Same argument as M1:
the purity floor is a swept knob (spec §10), and a release sweep with a `NaN` floor **silently admits
every tract** — it would report a finding about the data that is really a fact about a `NaN`. Now an
`assert!`; the module's three config guards are consistent.

### M3 · `flank_bp = 0` silently returns nothing, from any input
**Categories:** reliability · **Applied**

Every tract then fails its own flank test (`ref_start == tract_start`), so every locus is dropped and
nothing is even bundled — an empty result, no error, **from a knob spec §10 plans to sweep**. The M7
shape a third time: the knob appears to move and instead switches the step off. Now an `assert!` +
test.

### M4 · Both margin tests matched the same panic substring
**Categories:** reliability · **Applied**

Both declared `#[should_panic(expected = "windowed without a margin")]`, and that phrase was in **both**
messages — so the left-flank test (added specifically to kill the surviving `.max(bases_start)`
mutant) could be satisfied by a *right*-flank panic. It killed the mutant by luck of ordering; the
match was looser than the claim. Messages are now distinct (`flank runs LEFT/RIGHT of the window`),
both name `flank_bp` — "grow the window by *n* bp" is the actionable fix and the number is how much.

### M5 · `bundled`'s coordinate convention was pinned by nothing
**Categories:** refactor_safety, reliability (**convergent ×2**) · **Applied**

Every test reading `bundled` ran at `bases_start = 1`, where window offsets and genomic coordinates
**coincide** — so the documented "window offsets, not genomic" convention could silently flip and no
test would notice. Milestone D re-bases these to build each `SsrBundle` hull; if they were already
genomic it would **double-shift every hull**. Added
`bundled_coordinates_are_window_offsets_not_genomic` at `bases_start = 901`, asserting both the
positive and the negative. Mutation-verified.

### M6 · `ref_end - bases_start` was unchecked
**Categories:** errors · **Applied**

A wrong `contig_len` could wrap it into a plausible byte count, making the right-hand `assert!` fire
with a **garbage number and the wrong diagnosis**. Now `checked_sub`/`checked_add`.

### M7 · `assert!(bases_start >= 1)` had no test
**Categories:** reliability · **Applied**

## Minor (applied)

- **`From<Candidate> for RepeatInterval` read fields instead of destructuring** *(refactor_safety)* —
  the impl's contract is "must round-trip exactly", yet a new `Candidate` field would be silently
  dropped with no compile error, while a new `RepeatInterval` field *is* caught. Now destructured.
- **`fallible_impl_from = "warn"`** is enabled in `Cargo.toml` *(errors)* — a panicking `From` is
  against a policy this crate opted into, and clippy does not catch `expect`. Added the `PANIC-FREE:`
  marker the house style uses, with the reasoning and B2's retirement noted.
- **`ref_lo`/`ref_hi` → a `Range`** *(idiomatic)* — follows `Locus::tract_range()`'s existing idiom.

## Deferred, with reasons

| finding | why | home |
|---|---|---|
| **`Window { bases, bases_start, contig_len }`** — the three are co-dependent (the code proves it by asserting the relation on entry), two same-typed `u64`s transpose silently, and a `Window` built from the contig table is the only defence against M1's uncheckable half | **Warranted on the merits**, but spec §5a *and* the arch doc both spell the flat signature — so this is a **design change, not a cleanup**, and this skill does not make those. The reviewer's own preferred option is also the cheapest: **defer to C1's `Position`/`Bp` newtypes**, which the plan already schedules and which need no spec amendment. | C1 |
| `Admitted` is a participle whose `bundled` field is precisely what was *not* admitted | Fair, and the doc has to supply the missing noun in prose. But the name is spec §5a's. | spec, if ever |
| `finish_locus` carries four start/end pairs across two coordinate spaces with nothing in the names separating offset from genomic | Real, and in the one function whose whole job is not mixing them. C1's newtypes make the distinction type-level rather than nominal — the better fix, and scheduled. | C1 |
| The scanner-output differential throws `bundled` away (`admit_whole_contig` projects `.loci`), so the only realistic-density fixture asserts nothing about the partition | Genuine. D1's partition invariant is where this belongs — it is the milestone that owns "contiguous / non-overlapping / complete". | D1 |
| `..SsrAdmissionParams::default()` in `matched_params` | Pre-existing A2; `default_matches_the_frozen_catalog_params` pins every defaulted field. | A3-successor / C2 |
| no `proptest` property over `assert_agrees` | Still the strongest remaining lever on the port. | own step |

## Validation

```
cargo fmt -- --check                 clean
cargo clippy --lib --tests           clean (0 in scope)
cargo test --lib ng::region_typing   54 passed; 0 failed      (was 50)
cargo test --lib                     1820 passed; 0 failed; 4 ignored
cargo test --release --lib ng::region_typing   54 passed; 0 failed
cargo doc --no-deps                  0 errors in region_typing
```

**A profile note worth keeping.** `cargo test --release --lib` fails **6** tests crate-wide — all
pre-existing `should_panic` tests over `debug_assert`s, named for it
(`lgamma_panics_on_non_positive_in_debug`, `sdust_mask_debug_asserts_on_tiny_window`, …). They fail
identically on a stashed clean tree, so release is not a project gate. The counts are the point,
though: **the clean tree fails 7 and this tree fails 6** — the one that disappeared is
`admit_rejects_a_nan_purity_floor` (M2). `region_typing` is now the one module whose `should_panic`
guards hold in the profile the experiments actually run in.

**Mutation-verified after the fixes:**

| mutant | result |
|---|---|
| clamp the right flank at the window end (**the §2.6 bug**) | killed |
| clamp the left flank at `bases_start` instead of `1` | killed (needed a dedicated left-edge test) |
| ignore `bases_start` (window-relative output) | killed |
| `bundled` always empty | killed |
| `split_bundles` returns the cluster as isolated | killed |
| `bundled` leaks genomic coordinates | killed |
| `flank_bp >= 1` guard deleted | killed |
| `bases_start >= 1` guard deleted | killed |
| `contig_len` guard back to `debug_assert` | killed **in release** (debug cannot distinguish them) |

**Production untouched** — nothing outside `src/ng/` and the docs/reports.
