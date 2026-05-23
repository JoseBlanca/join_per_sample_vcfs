# Fix Application Report: ref_fetcher_2026-05-23.md

**Date:** 2026-05-23
**Source review:** `doc/devel/reports/reviews/ref_fetcher_2026-05-23.md`
**Source state reviewed against:** commit `5b9e420` (post-PSP H1+H3)
**Execution mode:** interactive
**Overall status:** *in progress*

---

## 1. Executive summary

### Review totals
- Blockers: 3
- Majors: 23 (M1–M27 with M3, M6, M7, M11, M16, M17, M21, M23 noted as convergent/subsumed)
- Minors: 19 (Mi1–Mi19)
- Nits: grouped

### Open-question answers from the author

1. **Migration completion**: complete it; "we'll discuss the fix" if some functionality can't be served by the new functions.
2. **`#[non_exhaustive]` on `ChromRefFetchError`**: trust the reviewer — apply.
3. **`ChromRefFetcher` sealing**: trait is internal-to-crate, not necessary at the module level — seal it.
4. **IUPAC handling**: fold non-ACGT to N.
5. **Walker migration shape**: B — replace `RefSeqFetcher` with a typed-error trait `MultiChromRefFetcher`.

*(detailed per-finding outcomes appear in §4 below)*

## 2. Findings table

*(populated incrementally as each finding is processed)*

## 3. Questions asked and answers

1. **Migration shape (precondition for M12 + B2 + M3 + M4 + M6 + M7 + M16 + M21 + M23 + Mi8 + Mi17 + Mi18)** — see *Open-question answers* §1.5. Answer: option B (replace with typed-error trait).

## 4. Per-finding log

*(populated incrementally)*

## 5–12. (sections to be filled at end of run)
