# Review: ng_alignment_d2 (research harness + finding)
**Date:** 2026-07-23
**Reviewer:** rust-code-review skill (single-agent validation)
**Scope:** step D2 — `examples/ssr_delimiter_comparison.rs` and the report
[ssr_delimiter_3v4_comparison_2026-07-23.md](../research/ssr_delimiter_3v4_comparison_2026-07-23.md)
**Status:** Approve-with-changes (all applied)

---

For a research deliverable the review question is not code style but **is the finding valid, and is
the comparison fair** — a harness bug could make one aligner look artificially good, and the report
draws a strong conclusion that *contradicts the spec's stated hypothesis*.

### Verdict: the finding is valid, not a harness artifact

The reviewer confirmed the comparison is **fair** — both aligners score the identical read bytes
against the identical ground truth with symmetric metrics — reproduced the table on the shipped seed,
and confirmed the qualitative result under a third independent seed. Findings #1–#3 (algorithm 4
uniformly better; algorithm 3 carries a small reference-pull bias; the §4.2 rounding-away fear does
not materialise) hold as stated.

### Two defects, neither inverting the conclusion — both applied

**Major (report accuracy): finding #4's stated *mechanism* was wrong.** The report claimed algorithm
3 loses the period-1 different-base cell because reference-pull reinterprets the +1 foreign insertion
as a same-length substitution. A hand-read sweep showed algorithm 3 measures the +1 **correctly at
every interior position** and errs **only at the left tract boundary**, collapsing by a jitter-scaled
1–3 bp — a boundary effect, not a uniform substitution, and the boundary model reproduces the −0.13
aggregate while the substitution story does not. The *finding* (algorithm 4 wins, biggest at
different-base, the split is mandatory) follows from the table and is unaffected; only the
explanation was wrong. **Fixed:** the report now states the boundary-collapse mechanism and flags the
position-0 tract-vs-flank definitional edge.

**Minor (harness fairness confound): the flanks had a tiling collision.** The comment claimed "no
motif collision" but the flanks touched a motif base at a junction. Investigating this surfaced that
the *first* corrected pick was worse — its right flank started with period 3's own `motif[0]`,
sliding the boundary a whole unit into the flank and collapsing period-3 clean accuracy to 83%. That
made the real collision condition explicit: a **tiling continuation** (`motif[0]` at the right junction,
`motif[period−1]` at the left), not a coincidental shared base. **Fixed:** flanks now genuinely
non-colliding (left ends in `C`, right starts in `G`), the comment states the tiling rule and records
the 83% cautionary tale, and the table was re-run (the finding is unchanged: period-1 different-base
92.6% → 99.5%).

Other checks passed: ground truth correct for the ordinary scenarios; sequencing error and jitter
reach both aligners as identical bytes; metric math correct; near-zero unmeasured, and the same read
for both aligners.

### Author response

The Major and the Minor **fixed before commit**; the report's numbers and mechanism now match a
re-run on non-colliding flanks, and the finding is stronger for having survived the scrutiny.
