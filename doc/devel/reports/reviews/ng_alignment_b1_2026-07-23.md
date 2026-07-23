# Code Review: ng_alignment_b1
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step B1 of `doc/devel/ng/impl_plan/alignment_best_path.md` — algorithm 3's two-regime matrix and traceback, unbanded
**Status:** Approve-with-changes (all applied — see §Author response)

---

### 1. Scope

The uncommitted diff for B1 (A0–A3 committed through `4ed24c2`): new
[ssr_best_path_flat_gap.rs](../../../../src/ng/alignment/ssr_best_path_flat_gap.rs), plus a one-line
`mod.rs` change. Three reviewers, all confirming scope via `git status`.

**Categories dispatched (7):** a dedicated **port-fidelity** pass (line-by-line against production's
`delimit_read`), reliability, naming, idiomatic, smells, module_structure, defaults.

This step is the module's **only byte-parity oracle** and the thing the STR locus generator has been
blocked on, and its failure mode is silent — a wrong recurrence is a wrong repeat length, not a
panic. The review was scoped accordingly.

### 2. Verdict

**Approve-with-changes.** **3 Blockers, 2 Majors**, plus Minors — *none of them in the transcription*.

**The port itself is faithful.** The fidelity reviewer read production's `alignment.rs:25-328` and the
port side by side, worked all nine checkpoints explicitly, and found every one matching: the
recurrences (including that the deletion candidate correctly reads `current`, not `previous`), the
transition costs, the candidate ordering at all six sites, `best_of`'s strict `>` first-on-ties
semantics, the tract-window predicate, terminal-state selection, the traceback's M-reports-`i−1` /
D-reports-`i` asymmetry and its `right_len > 0` guard, the emission wiring, and the exact-complement
mapping of `left_anchored`/`right_anchored` onto production's `left_off`/`right_off`.

**It then went further and proved it.** Production's `delimit_read` is `pub(crate)`, so the reviewer
temporarily added an in-crate differential harness running both implementations on identical
randomized inputs — random motifs, flanks, ±4-unit length changes, arbitrary indels, substitutions,
truncations, Q0–Q40. **~800,000 cases across four seeds, debug and release: zero divergences**, about
570,000 of them production `Region` results, every one measured to the same two offsets. The harness
was removed and the file restored byte-for-byte (`diff -q` clean, independently confirmed).

**Every finding below is about what the tests could not see.**

### 3. Execution status

Container, after fixes: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -- -D
warnings` exit 0; `cargo test --lib` **2185 passed, 0 failed, 4 ignored**. The two pre-existing red
commands are unchanged.

**One anomaly, recorded rather than dismissed.** Two early runs of the fidelity reviewer's small
4,000-case sweep reported disagreements, before the large sweeps. They could not be reproduced: five
of the reported cases replay as *agreeing* when pinned as fixed inputs, and 18 subsequent runs —
including six at the identical seed and filter, and three with forced recompiles — were clean. Both
anomalous runs were incremental recompiles of changed source, so a stale incremental artifact is the
plausible explanation. **B3 builds a permanent parity harness; if this is real rather than an
artifact, that is where it will resurface**, and the reproduction path is in the per-category report.

### 4. Open questions and assumptions

1. **Is production's dropped `tract_start > tract_end` branch a real gap?** Two reviewers reached the
   same answer **independently and by different routes**: the fidelity reviewer proved it unreachable
   structurally (traceback offsets are monotone in reference position, and `left + right ≤ n` is
   guaranteed by the clamps); the reliability reviewer swept 5×5×5×6 and found zero violations. So it
   is not a gap — but the reasoning was load-bearing and unrecorded, which is itself the finding.
   **Resolved by documenting the lemma and pinning it with a sweeping test, not by adding a redundant
   guard.**

### 5. Top 3 priorities

1. **B1-1** — the M > D > I tie-break, the rule this module documents most emphatically, had no test.
2. **B1-2** — both flanks were 8 bp, so transposing them passed everything.
3. **B1-3** — the tract window's right edge was untested while its left edge was.

### 6. Findings

#### Blocker — all three demonstrated by mutating the source and re-running, not by argument

**B1-1: the M > D > I tie-break had no test.**
Reversing the entire candidate order in the match recurrence to I > D > M left **all ten** original
tests passing. This is the module's own stated "correctness rule, not a nicety" — without it,
equally-scoring line-ups let two reads of the same molecule measure different repeats from identical
input — and `mod.rs`'s `TODO(Milestone B)` names the owed test by name.
**Fix applied, and the first attempt at it did not work.** A test using an emission model where a
match and a mismatch score alike (ε = 0.75, the crossover from step A1) still failed to catch the
mutation, because the three match candidates remain strictly ordered by their **transition** costs
(`ln_match_to_match` ≈ ln(0.99994) against `ln_gap_close` = ln(0.632)) — with a strict winner the
candidate order is never consulted. The shipped test flattens **both** halves: all transitions cost
zero and all emissions are equal, so every path scores identically and only the preference order
decides. **Verified by re-running the mutation: the new test fails under it and passes without it.**
A direct unit test of `best_of`'s first-on-ties semantics was added alongside, since that is the
mechanism a `>=` slip would break.

**B1-2: both flanks were 8 bp, so the two were interchangeable.**
`frame()` was always called with `ACGTACGT` / `TTGGTTGG` — different sequences, **same length**, and
only the lengths reach the algorithm. Transposing `left_flank_len` and `right_flank_len` passed all
ten tests. This is the third time in this plan that equal-valued fixtures hid a transposition (step
A3 hid two).
**Fix applied:** the right flank is now 10 bp, and `the_two_flanks_are_not_interchangeable` asserts
the swap changes the measurement.

**B1-3: the tract window's right edge was untested.**
Widening `column <= tract_last_column` to `<` passed all ten tests, while the equivalent mutation on
the *left* edge failed two — asymmetric coverage of a two-sided window.
**Fix applied:** `the_tract_window_includes_its_last_column_and_excludes_the_flank` checks all four
columns either side of both junctions.

#### Major

**B1-4: the dropped `tract_start > tract_end` rejection was undocumented.** Not a defect (see open
question 1), but the omission was silent, and its consequence would be severe: step B2 turns these
offsets into a `RepeatSpan`, where an inverted pair becomes a `Between(..)` whose `measured_length()`
saturates to `Some(0)` — **a confident zero-unit allele**, which is far worse than a rejection.
**Fix applied:** the lemma is written out at the site, a `debug_assert!` guards it, and
`the_tract_offsets_are_never_inverted` sweeps flank sizes, tract sizes and reads off either end.

**B1-5: the junction-gap 5′ rule and other determinism properties had no coverage.** Partly addressed
by B1-1's fix; `two_reads_of_the_same_molecule_measure_the_same_repeat` added.

#### Minor (applied)

- **A broken intra-doc link** (`BestPathAligner::Scratch`) — invisible because `cargo doc` is already
  red on pre-existing links elsewhere, so clippy passing did not catch it.
- **`Bp::get() as usize` can wrap on a 32-bit target**, and the following `.min()` would turn a
  wrapped value into a *plausible* flank length. Now `try_from(..).unwrap_or(usize::MAX)`, so the
  clamp stays meaningful. The only truncating cast of the ~20 audited; the rest widen.
- **The two `.min()` flank clamps were undocumented** and have no production counterpart. Now
  documented as no-ops whenever the precondition holds, with the cost of the case where they bind
  stated plainly.
- **`best_of` dropped production's empty-slice `debug_assert`** while keeping the panic. Restored,
  with a `# Panics` section.
- **The traceback's `_ =>` arm silently accepts states 1 and 3..=255 as "insertion".** A
  `debug_assert_eq!` now keeps it honest. *(The reviewer's larger suggestion — replace the
  `const usize` state indices with a `#[repr(u8)]` enum — is deferred: it is a genuine improvement
  but it is a refactor of a port whose byte-identity was just established over 800,000 cases, and
  this plan's rule is transcribe first, change separately.)*
- **The scratch-reuse test varied only the read length**, so `stride` never changed and the
  shrinking-reference ordering was untested. Added. Also documented **why** grow-without-clear is
  safe today — the fill is exhaustive — and that **Milestone C makes it a live hazard**, since
  banding leaves cells unwritten.
- **Degenerate inputs were unpinned**: empty read, read shorter than the flanks, a frame with no
  flanks at all. All probed by the reviewer and found non-panicking with defined answers; now pinned.

#### Notable non-findings

- **No numerical defect, and the reviewer tried hard to find one** — emission models returning `-inf`
  and `+inf` (the latter genuinely producing `NaN` via `inf + -inf`), at read lengths shorter than,
  equal to and longer than the reference. No panic, no underflow, no infinite loop. The invariant
  that saves it is that **deletions carry no emission**, so row 0's deletion chain is finite whatever
  the emission model does, and no all-`UNREACHABLE` cell is ever selected.
- **The match→match inconsistency is reproduced exactly**, and the module doc's numbers all check
  out: 1.019942 ≈ 1.02, 344.8× ≈ 350×, 1.2183 ≈ 1.2×, 1.8076 ≈ 1.8×.
- **The renames carry no semantic drift** — the renamed state constants hold production's numeric
  values.

### 7. Out of scope observations

The two pre-existing red validation commands. Also flagged, not applied: an all-Q0 read shifts the
measured tract by two bases, which is correct behaviour but worth knowing when B3 chooses fixtures.

### 8. Missing tests added now

`best_of_keeps_the_first_candidate_on_a_tie`, `ties_are_broken_toward_match_then_deletion_then_insertion`,
`two_reads_of_the_same_molecule_measure_the_same_repeat`,
`the_tract_window_includes_its_last_column_and_excludes_the_flank`,
`the_tract_offsets_are_never_inverted`, `the_two_flanks_are_not_interchangeable`,
`degenerate_reads_have_defined_answers`, plus the shrinking-scratch case.

### 9. What's good

- The module doc states the match→match inconsistency, quantifies its bias, and says why it is
  reproduced rather than fixed — the reviewer verified every number.
- The traceback preserves production's subtle M/D offset asymmetry exactly, which is the single
  easiest thing to get wrong in this port and the hardest to notice.

### 10. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib
```

### Author response

All three Blockers, both Majors and every Minor **fixed in this step's commit**. One suggestion
deliberately deferred (the state enum), recorded above with its reason.

Per-category audit trail at `tmp/review_2026-07-23_ng-alignment-b1/`.
