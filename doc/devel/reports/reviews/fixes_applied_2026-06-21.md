# Fixes applied — ssr-call reading layer review (2026-06-21)

Applies [ssr_call_reading_2026-06-21.md](ssr_call_reading_2026-06-21.md). Branch
`ssr-cohort`, fix commit `c15280c`. All gates green: fmt / clippy `-D warnings` / doc
clean; **1165 lib tests, 0 failed** (+6).

## Open questions — resolved

1. **Same catalog + sorted = a caller contract, but `ssr-call` enforces it with hard
   errors.** All inputs are assumed genotyped against the one sorted catalog passed to
   `ssr-call`; any violation is a typed hard error (not a panic, not a silent drop).
   No catalog-content hash (the boundary errors are the enforcement). → drives B1, M1.
2. **`.ssr.psp` is treated as untrusted / possibly corrupt** (defensive by default —
   may not be ours, may be bit-rotted). Malformed coordinates are rejected at the
   **decode boundary**, so every reader is protected. → drives M2.

## Findings

| id | resolution |
|---|---|
| **B1** | **Fixed.** `evidence_at` `q > held` → `SsrCohortReadError::LocusNotInCatalog` (typed hard error: the sample holds a locus the sorted catalog lacks) instead of `panic!`. Test `skipping_a_held_locus_is_a_hard_error`. |
| **M1** | **Fixed.** `CohortMerger::next_locus` validates strict catalog monotonicity → `SsrMergeError::UnsortedCatalog` before any cursor is asked. The cursor rewind `debug_assert` is now belt-and-suspenders (unreachable from data). Test `errors_on_an_unsorted_catalog`. |
| **M2** | **Fixed at the decode boundary.** `SsrDecoder::next_record` rejects `start == 0` / `span == 0` (`BlockStructureInvalid`), symmetric to the writer's `validate_locus`. The cohort `-1` shift is now safe by construction (`adapt`'s `debug_assert` documents it). Tests `next_record_rejects_zero_start` / `_zero_span`. |
| **M3** | **Fixed.** `stage1_to_cohort_coordinate_round_trip` writes a `SsrLocusObs` through the real Stage-1 `to_container_record` and reads it back through the cursor, asserting original 0-based coords — pins the `+1`/`-1` shifts to one test. |
| **M4** | **Fixed.** Added `owned_records_iter_matches_borrowing_for_ssr_kind` (the `SsrKind` decoder the cohort actually runs, multi-block) + `owned_records_iter_on_empty_file_yields_none`. |
| **M5** | **Fixed.** `run_ssr_call` warns (stderr) when `--threads` / `--queue-depth` get non-default values; help text marks them `RESERVED`; `--output` help + module banner now say TSV, not VCF (covers **Mi3**). `DEFAULT_THREADS` const added. |
| **Mi1** | **Fixed.** `CohortLocus.ref_tract` → `ref_frame` (it holds tract **+ flanks**). |
| **Mi2** | **Fixed.** `SsrMergeError::Read` now carries the `sample` label (open + streaming), so a many-file cohort points at the offending input. |
| **Mi3** | **Fixed** (with M5). |
| **Mi4** | **Deferred.** The tree-wide `#![allow(dead_code)]` on `src/ssr/mod.rs` is a tracked temporary (removed once the genotyping EM consumes the cohort surface); narrowing it now would cascade warnings across the still-filling `ssr` tree. |
| **Mi5** | **Fixed.** `catalog_reference_md5` lifted to `registry_ssr::CATALOG_REFERENCE_MD5_PARAM`; merge.rs, test_support.rs, and the Stage-1 writer all import it. |
| **Mi6** | **Fixed.** Exhaustive `SsrQc` construction in the driver test (no `..default()`). |
| **Mi7** | **Fixed.** `// PANIC-FREE:` invariant comments on the two `evidence_at` `expect`s. |
| **Mi8** | **Fixed.** `// REVIEW ON UPGRADE:` marker on the `ParameterValue` catch-all. |
| Nits | **Applied:** `alleles=` → `distinct=` in the dump (it counts distinct observations); `chrom_names` dense-id invariant documented; the `reader.rs` coordinate comment corrected (with M2). **Left:** `from_utf8_lossy` on the motif (safe — `Motif` is validated ACGT); the `(u32,u32,usize)` driver-test tuple and the repeated `path.display()` (cosmetic). |

## New tests (6)

`next_record_rejects_zero_start`, `next_record_rejects_zero_span` (registry_ssr);
`stage1_to_cohort_coordinate_round_trip`, `skipping_a_held_locus_is_a_hard_error`
(cohort reader — the latter replaces the old panic test); `errors_on_an_unsorted_catalog`
(merge); `owned_records_iter_matches_borrowing_for_ssr_kind`,
`owned_records_iter_on_empty_file_yields_none` (psp reader).

## Not done

- **Mi4** (tree-wide `allow(dead_code)`) — deferred, see above.
- Cosmetic nits (motif lossy-utf8, test tuple, repeated `display()`).
