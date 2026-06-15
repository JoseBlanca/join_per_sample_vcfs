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

## 0. Pending: code review of the §10 container work — **do first, fresh conversation**

The container generalization landed as 5 commits (`1eae1e9..f0dfffc`). Review it
with the code-review skill in a **new conversation** (large diff, own context)
before Stage-1 builds on top. Judgment calls to scrutinise:

- **Two `SsrLocusRecord` types** — Stage-1 chrom-**name**-keyed
  ([`src/ssr/pileup/locus_record.rs`](../../../src/ssr/pileup/locus_record.rs))
  vs container chrom_**id**-keyed
  ([`src/psp/registry_ssr.rs`](../../../src/psp/registry_ssr.rs)). The most
  consequential call — it defines the §2 driver's adapter. Confirm or unify.
- **`ColumnDef.key` removed** → schema-agnostic `ColumnDef`; each schema
  dispatches via its own `from_tag` (`ColumnKey` / `SsrColumnKey`). Check the
  `tag()`-exhaustive vs `from_tag`-array duplication and the M4 guarantee.
- **Block-header invariant relaxed** (`n_total_alleles >= n_records` dropped;
  `AllelesLessThanRecords` retained but never raised). Confirm SNP integrity is
  still fully covered by write-side + read-side checks.
- **`#[allow(private_bounds)]`** on the `RecordsIter` impls (encapsulating the
  `pub(crate) BlockDecoder` bound vs widening the public API).
- **SSR column schema** — `span` vs storing `end`; `n-spanning` doubling as the
  per-record grouping count; no separate `profile-len` column (CSR offsets
  carry it); `finite_constraint: false` on `amb-logliks` (no NaN/inf sweep on
  the SSR decode path).

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
   source).* Plan drafted ([ssr_catalog.md](ssr_catalog.md)) but **verify
   whether it is built**; the work pivoted to the Stage-1 compute path, so the
   catalog may not exist yet. If not, building it (or pinning the locus-input
   format the fetcher consumes) is the first gate. Open: TRF-mod BED column
   order, post-process order, `flank_bp` default.

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
