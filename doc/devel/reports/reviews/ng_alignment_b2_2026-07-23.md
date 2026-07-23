# Code Review: ng_alignment_b2
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step B2 of `doc/devel/ng/impl_plan/alignment_best_path.md` — the `RepeatSpan` readout, the widening
**Status:** Approve-with-changes (all applied — see §Author response)

---

### 1. Scope

The uncommitted diff for B2 (B1 committed as `748df78`): `TractReadout::classify`, the
`BestPathAligner` impl, the trait's `Context` gaining a lifetime, and new tests. Two reviewers, seven
categories: reliability, errors, naming, idiomatic, smells, refactor_safety, module_structure. B1's
recurrence and traceback were explicitly out of scope — committed, and already verified byte-identical
to production over ~800,000 differential cases.

### 2. Verdict

**Approve-with-changes.** **3 Majors**, no Blockers. The most important one was found **independently
by both reviewers**, and the second demonstrated it live.

### 3. Execution status

Container, after fixes: `cargo fmt --check` exit 0; `cargo clippy --all-targets --all-features -- -D
warnings` exit 0; `cargo test --lib` **2194 passed, 0 failed, 4 ignored**. Both reviewers restored the
working tree and verified it.

### 4. The side mapping — verified by mutation, not by argument

The step's whole risk is that a mis-assigned side is a wrong observation class on a read that still
looks perfectly good. One reviewer ran four mutations, each reverted:

| mutation | caught by |
|---|---|
| swap `FromLeft`/`FromRight` in `classify` | 2 tests |
| swap the `left_anchored`/`right_anchored` initialisers in `delimit` | 3 tests |
| `(false, false)` returns a span instead of `Unanchored` | 2 tests |
| transpose the span endpoints in `classify` | 6 tests |

**The mapping is correct**: `FromLeft` ← `(left_anchored = true, right_anchored = false)`, so it does
mean "the left flank held, and the read ran off its own 3′ end". The fixtures avoid the
paired-equal-values pattern that hid transpositions in three earlier steps of this plan.

### 5. Findings

#### Major

**B2-1: a flankless geometry turned an unmeasured read into a `Between`.**
**Categories:** reliability, idiomatic — **found independently by both reviewers.** Demonstrated live:
frame `("", "CAGCAGCAGCAG", "")`, read `CAGCAG` → `Between(0..6)`, `is_measurement() == true`,
`measured_length() == Some(6)` — **a short allele nothing observed.**

The cause is mine, from B1: the convention "a flank that does not exist cannot fail to hold" was
harmless while `left_anchored`/`right_anchored` were two raw bits, but **B2 is the step that promotes
them to a measurement claim**. `RepeatGeometry`'s own field docs admit zero-length flanks at contig
edges, so this is reachable at exactly the geometry where evidence is thinnest.

**Fix applied — the convention is inverted: a flank that does not exist cannot *anchor*.** The
reasoning is that at a flankless edge the two cases that matter produce an **identical** readout — a
read that ended because the tract ended, and one that ran out mid-tract both give
`tract_end == read_len` — and nothing can tell them apart. Calling that a measurement over-claims for
one of the two. The cost is that a flankless locus can never yield a measurement, which is honest:
with no flank in the frame there is nothing to pin the tract's edge against. **This does not disturb
the parity oracle** — production's `Delimited::Region` makes no measurement claim, and B3 compares the
measured *bytes*, which are unchanged.

**B2-2: `the_classified_span_agrees_with_the_raw_readout` was a tautology.**
**Categories:** reliability. `align` *is* `delimit(..).map_or(Unanchored, classify)`, and the test
compared its result against `readout.classify()` — `classify(x)` versus `classify(x)`. **Proven: the
most dangerous mutation, swapping `FromLeft` and `FromRight`, left it passing.**
**Fix applied:** rewritten to check the part that is not definitionally true — that the wiring carries
B1's offsets through unchanged, and that the case follows from anchoring flags spelled out
independently in the test rather than by calling `classify`.

**B2-3: `classify` carried no `tract_start <= tract_end` assertion.**
**Categories:** errors. `delimit` has one, but `TractReadout` is public with public fields — B3's
parity harness needs the raw offsets — and the new tests already hand-build one. A hand-built inverted
readout gives `Between(10..4)` → `measured_length() == Some(0)`, the "confident zero-unit allele" B1's
own comment calls far worse than a rejection.
**Fix applied:** `debug_assert!` at the classification boundary, with the reasoning.

#### Minor (applied)

- **`length_lower_bound` misleads for an empty reference.** The `Unanchored` *variant* is right (arch
  §3), but the accessor's `read_len` fallback would assert a repeat at least as long as the read **at
  a locus with no bases at all**. Probed live: 20-base read, empty reference → bound of 20. Documented
  as a caller-side distinction the accessor cannot make.
- **`type Context<'a>` lacked the `Copy` bound its own justification rests on** — a generic
  one-locus-many-reads helper would not compile. Added, with the reason.
- The zero-flank behaviour had no test at all (both reviewers). Added.

#### Deferred, recorded

The stale `TODO(Milestone B)` in `mod.rs` is now partly paid: the tie-break and empty-reference tests
exist. **Two of its four remain genuinely owed** — the 5′-junction-gap rule and the debug-assert
precondition tests — verified by grep, and left in place rather than quietly trimmed.

### 6. On the GAT — assessed and upheld

Both reviewers were asked to judge the deviation from arch §3's plain `type Context` taken by
reference. One assessed all three alternatives and rejected each: **separate borrowed arguments**
deletes `Context = ()` and forces the affine aligner to accept a geometry it does not want;
**`trait BestPathAligner<'a>`** moves a per-call fact into the type and makes every generic consumer
carry an HRTB; **`RepeatContext<'static>`** contradicts arch §2.2 outright. The other verified by
compilation that `type Context<'a> = ()` works for a throwaway impl, so the gated affine aligner pays
nothing. The GAT states the truth — the impl is lifetime-free, each call picks its own — and there is
exactly one implementor, so this was the cheapest moment to settle it.

Recorded as a Nit for the next plan: `MarginalAligner` (arch §3) will hit the identical wall, and the
rationale currently lives only on `BestPathAligner`.

### 7. Notable non-findings

- `classify` as an inherent method beats both `From<TractReadout>` (`.into()` has no name at the call
  site — the wrong spelling for this module's one silent-failure risk) and a free function.
- The `match (bool, bool)` is exhaustive over a real 2×2, not the co-dependent-bool smell. One thing
  worth knowing: exhaustiveness protects the *input* shape; a fifth `RepeatSpan` variant, which
  `mod.rs` explicitly contemplates, would produce no diagnostic in the sole producer.
- The `Unanchored` arm dropping the span is right, and `length_lower_bound`'s fallback reproduces it
  exactly — but only by construction (`left_anchored == false` ⟹ `tract_start == 0`, and likewise on
  the right), which nothing pinned.

### 8. Commands to re-verify

```
./scripts/dev.sh cargo fmt --check
./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
./scripts/dev.sh cargo test --lib
```

### Author response

All three Majors and every Minor **fixed in this step's commit**. Per-category audit trail at
`tmp/review_2026-07-23_ng-alignment-b2/`.
