# SSR Stage-1 — remaining work (reference checklist)

**As of:** 2026-06-15, branch `ssr-architecture`.
**Purpose:** the running to-do for SSR Stage 1 now that the `.psp` **container
generalization (architecture §10) is complete**. A reference to follow, not a
design doc — each build item gets its own plan/report when it's tackled.
Cross-refs: architecture [§8 roadmap](../architecture/ssr_genotyping_architecture.md)
(item 2 = Stage 1) and [§10](../architecture/ssr_genotyping_architecture.md);
the per-step reports under [`ia/reports/implementations/`](../../../ia/reports/implementations/)
(`psp_container_generalization_step{1a,1b,2,3,4}_*`).

---

## 0. DONE: code review of the §10 container work + fixes (2026-06-15)

Reviewed (rust-code-review skill, 10 categories) →
[psp_container_generalization_2026-06-15.md](../reports/reviews/psp_container_generalization_2026-06-15.md):
**Approve-with-changes** (0 Blockers, 5 Major, 11 Minor). SNP path verified
sound; every Major was in the unshipped SSR schema. Fixes applied →
[fixes_applied_2026-06-15.md](../reports/reviews/fixes_applied_2026-06-15.md):
all 5 Majors + 9 Minors Applied (M1 inclusive `last_pos`; M2/M3 typed
structural errors + `sum(n_spanning)` check; M4 poison-on-`KindMismatch`;
M5 `column_key!` macro; Mi3 mandatory `kind`; Mi2/Mi5–Mi11). Deferred: Mi1
(`n_total_alleles`→`n_entries` rename), Mi4 (kind-taxonomy module move), and the
M2/M3/Mi9 malformed-block rejection tests (need a raw-block fixture). Gates
green; 1142 lib tests. The judgment calls are settled there — including the two
`SsrLocusRecord` types (kept as a deliberate name vs container mirror; the §2
driver adapts between them).

---

## 1. Done (foundation)

**Stage-1 compute-path units** (all unit-tested, no driver yet):
`types`/`Allele`, shared `norm_seqs`, `count_repeats` (parked fast-path),
`pair_hmm` + `score_candidates`, `candidate_generation`, `triage`,
`read_analysis` (`analyze_read`), `locus_record` (`aggregate`, all-CSR),
`fetch_reads` **reservoir only** (Algorithm R + per-locus seed).

**`.psp` container (§10)** — generic core hosts `snp` + `ssr`: `PspKind` /
`BlockAccumulator` / `BlockDecoder`, `kind` header tag + registry-by-kind,
interval region query, `registry_ssr` (`SsrKind`/`SsrBlock`/`SsrDecoder` +
`write_locus` + `records_of::<SsrKind>()`). SNP byte-identical; SSR round-trips.

---

## 2. Remaining build pieces (dependency order)

Each lands with its own implementation plan + report; each carries Bucket-1
synthetic tests.

1. **Stage 0 — `ssr-catalog`** *(prerequisite: the fetcher needs a locus
   source).* Plan: [ssr_catalog.md](ssr_catalog.md). Building incrementally:
   - **`catalog/io.rs` — DONE (2026-06-15)**: the format contract
     (`CatalogHeader`/`CatalogWriter`/`CatalogReader` + `Locus`⇄row, bgzip TSV,
     6 round-trip tests). Report:
     [ssr_catalog_io_2026-06-15.md](../reports/implementations/ssr_catalog_io_2026-06-15.md).
   - **`catalog/postprocess.rs` — DONE (2026-06-15)**: `build_loci` — the
     period≤6 → drop-compound → drop-bundle → end-trim → recompute-purity →
     embed-`ref_seq` pipeline (faithful GangSTR `minimal_trim`/`remove_bundles`
     port; consumes `TrfRecord`, emits `Locus`; 11 unit tests). Also `trf.rs`
     `TrfRecord` (the type postprocess consumes). Report:
     [ssr_catalog_postprocess_2026-06-15.md](../reports/implementations/ssr_catalog_postprocess_2026-06-15.md).
   - **`catalog/trf.rs` spawn/parse — DONE (2026-06-15)**: `locate_trf_mod`
     (override → exe-sibling → PATH), `version` (`trf-mod -v`), `run_on_contig`
     (temp-file spawn, no pipes), `parse_bed_line` (10-col BED, format pinned
     against `trf_print_bed`). 4 tests incl. a live integration run of trf-mod
     on a synthetic CAG repeat. Report:
     [ssr_catalog_trf_2026-06-15.md](../reports/implementations/ssr_catalog_trf_2026-06-15.md).
   - **`run()` orchestrator + `ssr-catalog` CLI — NEXT (last Stage-0 piece)**:
     per-contig fan-out/collect (rayon, order-preserving) reading the reference
     (noodles-fasta) → `run_on_contig` → `build_loci` → `CatalogWriter`; compute
     `reference_md5`, build the `CatalogHeader`, write the CSI index; wire the
     `ssr-catalog` subcommand. Open: `min_score`/`flank_bp`/`bundle_threshold`
     defaults (pin here).

2. **`fetch_reads` I/O driver** *(the last missing primitive).* Add to
   [`fetch_reads.rs`](../../../src/ssr/pileup/fetch_reads.rs) (reservoir already
   there): the catalog walk + per-locus
   [`AlignmentMergedReader::query`](../../../src/bam/alignment_input.rs)
   (mirror SNP `run_pileup`: load handles once, share the FASTA repo, `clear()`
   per contig transition), the coordinate-reach admission gate (reuse triage's
   footprint), the per-locus reservoir cap, the per-locus read bundle. Build
   **single-threaded first** (trivially deterministic).

3. **Stage-1 driver** — orchestrate fetch → `triage` → `analyze_read` →
   `aggregate` → container `write_locus`. New seams this needs:
   - **name → chrom_id adapter**: Stage-1 `SsrLocusRecord` (chrom name) →
     container `SsrLocusRecord` (chrom_id), via the header chromosome table —
     the SNP-symmetric boundary (gated by the review's verdict on the two-type
     question).
   - **`WriterHeader` / chromosome-table construction** for the `.ssr.psp`
     (from the reference / catalog), incl. md5s — analogue of the SNP pileup's
     header build.

4. **`ssr-pileup` CLI subcommand** — wire the driver behind the CLI
   (mirror `pileup`); plan [ssr_pileup.md](ssr_pileup.md).

5. **Parallelism** (architecture §7 / §8.4) — fetcher-thread / bounded-queue /
   worker-pool / ordered-collector topology, **plus the determinism gate:**
   byte-identical `.ssr.psp` across `--threads ∈ {1, N}`, *including* the
   reservoir subsample (per-locus seed + total read order). Built after the
   single-threaded path works.

---

## 3. Deferred (explicitly not now)

- **Off-ladder candidate generation + `offl_*` columns** — the container has no
  off-ladder columns yet (the in-memory `spanning` profiles are on-ladder only;
  `candidate_generation` keeps the verbatim-tract off-ladder key parked). Add
  the columns + wiring when off-ladder generation lands.
- **The measured fast path (`count_repeats`)** — realign-everything is the
  default; `count_repeats` is the **MUST-TRY-AND-MEASURE** shortcut the user
  requires be tried + benchmarked against realign-everything before adoption.
- **Spec §4.3 amendment** — update
  [`per_sample_pileup_format.md`](../specs/per_sample_pileup_format.md)'s SSR
  section (or a new `per_sample_ssr_format.md`) to the **all-CSR** reality:
  drop `hist_*` and `n_flank_indel`; document the `amb-lengths`/`amb-logliks`
  CSR + the QC scalars as built. Owed since the all-CSR decision.

---

## 4. Beyond Stage 1

Stage 2 (`ssr-call`: EM + stutter kernel + VCF) and the crate/test **simulator**
(two emission levels, anti-tautology boundary) — architecture §8 items 3–4,
each designed just-in-time before building.
