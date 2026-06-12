# Stage 0 — `ssr-catalog` (the catalog builder)

**Status:** first draft, 2026-06-12 — a discussion starter, to work on together.
The first *stage* pass after the shared types
([ssr_shared_types.md](ssr_shared_types.md)); follows the overall
[architecture](ssr_genotyping_architecture.md) (§8 module 1) and the spec
([ssr_genotyping.md](../specs/ssr_genotyping.md) §3). The spec settles *what* the
catalog is and *why*; this designs *how the builder is structured*. **Open points
and recommendations are flagged to argue with.**

---

## 1. What Stage 0 does

`ssr-catalog` turns a **reference FASTA** into the **locus catalog** — the map of
every SSR worth genotyping, which Stages 1–2 then read (and read *only*, never
the genome — §3.2 of the spec, the embedded-reference decision). It runs **once
per reference**, offline, and its output is a self-contained, self-describing
artifact.

```
reference FASTA ──► [ detect ] ──► [ post-process ] ──► [ embed ref_seq ] ──► catalog.ssr_catalog.bed.gz (+ .tbi)
                     TRF-mod          period≤6, purity,      tract + flank          self-describing TSV
                   genome-wide        merge, split           per locus              (## header, bgzip+tabix)
```

It owns one responsibility the spec recently added: **embedding the local
reference** (`ref_seq` + `ref_seq_start`) into each row, so it is the *only*
stage that reads the FASTA for the SSR algorithm (§5 below).

---

## 2. How TRF enters the build — DECIDED: lh3/TRF-mod binary, shell-out

Two decisions (2026-06-12): **use Heng Li's [lh3/TRF-mod](https://github.com/lh3/TRF-mod)**
(the same TRF algorithm, more current and pipeline-friendly), and **shell out to
its binary** rather than reimplement it. FFI was already ruled out.

### 2.1 Why TRF-mod, and why the binary (not a port)

- **TRF-mod is the identical TRF algorithm with a pipeline-friendly interface.**
  One simple CLI — `trf-mod ref.fa > out.bed`, sensible default parameters — and
  **BED-like output on stdout**, versus upstream TRF's eight-argument invocation
  and HTML/`.dat` files. Ideal for spawning from Rust. Same detection quality
  (identical algorithm), better I/O.
- **Don't reimplement.** The project ports *small, hot-path* C tools (sdust →
  [dust_filter.rs](../../src/var_calling/dust_filter.rs), "ports `sdust_core` line
  by line"). TRF is the opposite on both axes — a large validated statistical
  algorithm run **once per genome, offline** — so a port buys no hot-path win and
  risks **silent divergence** from the detector that validates every STR catalog.
  (Same anti-reimplementation logic as the indel-norm decision, §4 of the types
  doc, but stronger: TRF *defines* what a locus is.)
- **Parallelism without a port.** A port would have given in-thread parallelism
  (as the BAQ port did). The binary recovers it at the **process level** — fan
  out one `trf-mod` per contig (rayon over contigs; GangSTR's `xargs -P`
  approach), which fits a once-per-genome batch perfectly (§8). We lose nothing
  that matters here.

### 2.2 Output-format consequence

TRF-mod emits **BED-like records on stdout**, *not* Benson's 15-column `.dat`. So
the parser (§3) targets TRF-mod's columns — it still computes coordinates, the
consensus motif, and a percent-match field (identical algorithm), which is all
the catalog needs. **Exact column order: confirm against TRF-mod itself when
building** (the one open detail).

### 2.3 Licensing — our code stays MIT

TRF-mod is **AGPL-3.0**. We want **our code MIT**, and that is fully compatible
with using *and even shipping* TRF-mod, because the copyleft boundary is the
**process/linking** boundary, not the **shipped-together** boundary:

- **Shelling out is arm's-length** — `trf-mod` is a separate process invoked over
  the CLI, not linked/FFI'd. AGPL does **not** propagate to our code. (Exactly why
  no-FFI matters: linking raises the derivative-work question, a subprocess
  doesn't.)
- **Distributing it alongside us is "mere aggregation"** — GPL/AGPL §5 states that
  including a covered work in an *aggregate* does **not** apply the license to the
  other parts. So we may ship `trf-mod` with our tool *and* keep our code MIT; the
  licenses stay separate. We honor AGPL *for that binary only* (ship its LICENSE,
  point to lh3/TRF-mod source + pinned version, keep notices).

> ⚠ **Not legal advice.** Standard, well-trodden pattern (countless MIT/BSD tools
> invoke GPL/AGPL binaries), but worth a review before a public release.

### 2.4 The user runs only `ssr-catalog` — auto-launch via layered discovery

The tool **finds and launches `trf-mod` itself**; the user never installs, names,
or invokes it. `ssr-catalog` resolves the binary in priority order and spawns it:

1. **`--trf-mod-path <path>`** if given — explicit override (and the
   bring-your-own escape hatch, e.g. AGPL-averse packaging or a custom build).
2. **A copy shipped with our tool**, located *relative to our own executable*
   (`std::env::current_exe()` → a sibling / `libexec`-style path). The default
   seamless case — invoked by absolute path, so no `PATH` clash and no need to
   rename the binary.
3. **`trf-mod` on `PATH`** — if a conda env or the user provided one.
4. Otherwise a **clear, actionable error** ("`ssr-catalog` needs `trf-mod`; pass
   `--trf-mod-path` or reinstall").

This discovery logic is **independent of how the binary is distributed** — it
works whether `trf-mod` was bundled, conda-installed, or hand-placed. Only
`ssr-catalog` triggers it; `ssr-pileup`/`ssr-call` never need `trf-mod`.

> **Deferred — distribution & container (not now).** *How* the per-platform
> `trf-mod` binary is produced and shipped (bundled, conda dependency,
> build-from-source, the dev-container install) is a **packaging-time** decision,
> out of scope for this design. The discovery logic above is what we build now;
> the binary's provenance is decided later. (The one inherent consequence already
> accepted by choosing shell-out over a port: the catalog step is no longer a
> single self-contained pure-Rust binary.)

---

## 3. TRF-mod invocation & output parsing (`catalog/trf.rs`)

*Format confirmed against the vendored source —
[`TRF-mod/src/trfrun.h`](../../TRF-mod/src/trfrun.h) `trf_print_bed()` and
[`TRF-mod/src/trf.c`](../../TRF-mod/src/trf.c). Not from the README (its
"BED-like" prose is loose and misled an earlier draft).*

- **Invocation** per contig: spawn the `trf-mod` resolved by §2.4 as
  `trf-mod <contig.fa>` and capture stdout. **BED is the default output**
  (`paramset.bedon = 1`; `-d`/`-n`/`-h` switch it off — we pass none). Default
  parameters are sensible (ULTRA-paper-based); pin the version.
- **The output is a real tab-separated BED**, one row per repeat:
  ```
  ctg   start   end   period   copyNum   fracMatch   fracGap   score   entropy   pattern
  ```
  Three things fall out *for free* (verified in the `fprintf`):
  - **`ctg` is column 1**, the contig name truncated at first whitespace
    (`>chr1 desc` → `chr1`) — each row is **self-contained**; no `@`-marker
    tracking (that's the `-n` *NGS* format, which is what the committed
    `t/*.out` fixtures use — not our path).
  - **`start = first − 1`, `end = last` → already 0-based half-open**, matching
    our catalog convention. **No coordinate conversion.**
  - **`fracMatch = matches × 0.01` → already a fraction in [0,1]** (a *sanity
    check* for purity — but `purity_fraction` is **recomputed** from the tract
    after split, §4, not taken from here).
- **Catalog mapping:** `ctg`→`chrom`; `start`/`end` as-is; `period`→the locus
  period; `pattern`→motif source (see below); `score` available as an early §4
  accept-gate; `fracMatch` is a purity sanity-check (purity is recomputed, §4).
  Ignore `copyNum`/`fracGap`/`entropy`.
- ⚠ **`period` ≠ `len(pattern)` is allowed** — the help warns "length of pattern
  may differ from period," and the source carries them as independent fields. So
  we **must not** derive period from the motif. **Resolution:** take **`period`
  from TRF's column** as authoritative, and set the catalog `motif` = the **first
  `period` bases of the reference tract** (from `ref_seq`). That is exactly our
  "verbatim, reference-strand, phase-faithful motif" decision (types §5) *and*
  guarantees `len(motif) == period`, so the `period = len(motif)` derivation
  stays valid; TRF's `pattern` is kept only as a sanity-check. (No types/spec
  change — this is how the builder *constructs* `motif`.)
- **Scope:** filter `period ≤ 6` in post-process (§4), not via TRF's `-p` knob —
  completeness over build speed (build is once-per-genome).
- **Golden fixtures:** the committed `t/*.out` are NGS format, so we generate our
  own BED goldens (run default `trf-mod` on `TRF-mod/t/*.fasta`, capture), pin
  the TRF-mod version, and assert column count/types at parse time — fail loudly
  on any upstream format change.

`catalog/trf.rs` spawns `trf-mod` per contig and parses its tab-separated BED
stdout into typed records for the post-processor.

---

## 4. Post-processing (`catalog/postprocess.rs`)

In Rust, over the parsed TRF-mod records, per contig. **Order (decided):**

1. **Period ≤ 6** — drop TRF calls outside SSR scope (cheap volume cut; acts on
   TRF's per-call period).
2. **Split compound loci** — separate overlapping *different-motif* calls into
   single-motif sub-loci (TRF does **not** do this — its redundancy elimination
   is period/score-based, §3). The **inner-flank consequence** (spec §3.1): a
   sub-locus's inner flank is the neighbour's repeat, not unique sequence.
   **Split records the overlapping sequence** (the shared inner region between
   adjacent sub-loci) so that structure is captured for Stage 1 (which then
   anchors on the *outer* flank). *(That overlap also lands inside each
   sub-locus's embedded `ref_seq` flank, §5 — the recording makes it explicit.)*
3. **Recompute purity, then filter.** **`purity_fraction` is recomputed here from
   the (sub-)tract bytes** — *not* taken from TRF. After split, TRF's `fracMatch`
   was computed on the original call and no longer describes a sub-locus's
   boundaries; and recomputing lets us use the spec's **exact** definition
   (fraction of the tract matching a perfect motif tiling, §3.2) rather than
   TRF's alignment-based `%match`. Recompute is **uniform** (split *and* unsplit)
   for one consistent definition; TRF's `fracMatch` is demoted to a sanity-check.
   The worker holds the contig resident, so the bytes are at hand. Then drop loci
   below the purity floor (an accuracy knob, §3.4). **Imperfect single-motif loci
   are kept** — the filter is a degeneracy floor, not perfect-only.

**No own merge — trust TRF (decided).** TRF-mod eliminates redundancy by default
(`IsRedundant`: period-multiples + same-period duplicates, §3); we **do not** add
our own merge. (A light same-motif boundary cleanup is added later *only if*
measurement shows a gap.)

**Why this order:** purity **after** split (a compound's combined purity is
misleadingly low → filtering first would wrongly drop real loci); period ≤ 6
first (cheap); merge delegated to TRF. *(A TRF-intrinsic `score` floor, if
wanted, fits as an early accept-gate on raw TRF calls; the recomputed `purity`
is the late, locus-level filter.)*

> **No mappability filter (decided).** An earlier draft had a
> mappability/unique-flank filter; it's **dropped**. **MAPQ already handles
> mappability** — it is paralog-aware by construction (a near-identical copy
> elsewhere gives the read a competing alignment, shrinking the best-vs-second
> gap and lowering MAPQ) and, for paired-end data, incorporates mate/insert-size
> rescue that a static single-end reference map cannot. So unmappable/paralogous
> loci self-suppress in Stage 1 (low-MAPQ reads → low depth → no-call), with no
> reference-side filter, no self-similarity search, and no segdup-track
> dependency. The catalog's locus universe is *all detected SSRs*; mappability is
> a per-read, per-experiment concern owned entirely by MAPQ. (Consistency note:
> the spec's §5.8 Stage-2 "segdup/mappability exclusion" should be revisited the
> same way when we reach the Stage-2 pass.)

---

## 5. Embedding the local reference (`ref_seq`) — Stage 0's new job

Per the spec §3.2 decision, each catalog row carries `ref_seq` (the tract + a
`flank_bp` margin, upper-cased, clamped at contig ends) and `ref_seq_start` (the
genomic coordinate of `ref_seq[0]`). Stage 0 is the natural and *only* place this
happens: it already holds the reference open (TRF ran on it), so after
post-processing settles each locus's `[start, end)`, the builder slices
`reference[start − flank_bp .. end + flank_bp]` (clamped) into `ref_seq`.

- **`flank_bp`** is a build parameter, recorded in the `##` header. It is sized
  to **Stage 1's** read-anchoring + pair-HMM band need — so its *value* is
  pinned when we design `ssr-pileup`, but the *mechanism* (store tract + margin)
  is fixed here. Lean: a generous default (tens of bp) with a CLI override.
- **md5 binding** — the header carries `reference_md5` (the spec's
  upper-cased-content convention), so any later mismatch between `ref_seq` and a
  declared reference is detectable.
- **Upper-casing** — `ref_seq` is upper-cased for consistency with the md5
  convention and so reconstruction/equality are case-stable.

---

## 6. No sdust prefilter — TRF-mod alone, genome-wide (decided)

An earlier draft kept the spec's *optional* sdust pre-search (mask the genome,
restrict TRF to masked windows) as a speed lever. **Dropped.** Reasons:

- **Speed-only, and premature** — its sole purpose is to run TRF on less
  sequence; we have no evidence genome-wide TRF is too slow, and the build is
  once-per-genome, already parallelized by contig (§8).
- **Lossy** — sdust isn't motif-aware, so restricting TRF to its windows risks
  missing/mis-bounding loci. A recall risk traded for unmeasured speed is a bad
  default for a tool meant to be trustworthy.
- **Non-standard** — TRF genome-wide is the validated path (GangSTR/HipSTR).
- **A strictly better lever already exists** — if build time ever hurts, the
  deferred **intra-contig chunking** (§8.5) is *lossless* (split the work, drop
  no sequence), unlike the prefilter. So the "what if it's slow?" case has a
  better answer than this.

Consequence: **SSR uses sdust nowhere** — the prefilter was its only touchpoint.
sdust remains the SNP path's masker; it simply drops out of the SSR reuse
surface.

---

## 7. Catalog format & writer (`catalog/format.rs`)

Mostly settled in spec §3.2; the builder's job is to write it:

- **One self-describing bgzip+tabix BED-like TSV.** `##` metadata header
  (reference path + `reference_md5`, TRF-mod version + params, post-process filters,
  `flank_bp`, tool/version, date), a `#` column header, then rows. Tabix skips
  comment lines. No sidecar.
- **Columns:** `chrom  start  end  motif  purity_fraction  ref_seq_start
  ref_seq`. Everything else (`period`, `ref_copies`, perfect/imperfect,
  canonical motif class, `locus_id`) is **derived**, not stored.
- **bgzip + tabix** so Stages 1–2 can region-query the catalog (the
  `--regions` path), keyed on the tract interval.

---

## 8. Data-flow & worker architecture

Contigs are **fully independent** — TRF runs per contig, post-processing is
per-locus within a contig — so Stage 0 is a clean fan-out/collect: one **worker
per contig**, a pool bounded by `--num-chroms-in-parallel`, and an **ordered
collector** that reassembles the reference's contig order.

```
                ┌──────────────────────── worker (one per contig) ────────────────────────┐
reference  ──►  │ read contig seq → spawn trf-mod → parse BED → post-process → slice ref_seq│ ──┐
  FASTA    ──►  │   (held resident)   (subprocess)    (§3)        (§4)         (§5, has bytes)│   │
  (per     ──►  └──────────────────────────────────────────────────────────────────────────┘   │  finished, start-sorted
  contig)        … up to --num-chroms-in-parallel workers concurrently (work-stealing) …        │  rows per contig
                                                                                                 ▼
                              ┌─────────────── ordered collector ───────────────┐
                              │ buffer by contig index; emit in reference contig │ ──►  catalog.ssr_catalog.bed.gz (+ .tbi)
                              │ order, each contig already start-sorted          │
                              └──────────────────────────────────────────────────┘
```

```
src/ssr/catalog/
├── trf.rs          # spawn `trf-mod` per contig + parse BED-like stdout → typed records
├── postprocess.rs  # period≤6, purity, merge, split compounds (no mappability — MAPQ owns it)
└── format.rs       # self-describing bgzip+tabix TSV read + write (the collector + writer)
```
(plus a small `ssr-catalog` subcommand under `pop_var_caller/cli/`.)

### 8.1 The worker

Each worker owns one contig end-to-end: read its sequence (held resident,
because it feeds `trf-mod` *and* is sliced for `ref_seq`), spawn `trf-mod` on it
(§2.4 discovery), parse the BED (§3), post-process (§4), embed `ref_seq` (§5).
Output: that contig's **finished, start-sorted** catalog rows. No shared state —
this is **process-level parallelism**, recovering what an in-process TRF port
(cf. the BAQ port) would have given us, which is the right shape for a
once-per-genome batch.

### 8.2 The ordered collector — an output *correctness* requirement

Workers finish out of order (a tiny unplaced scaffold beats chr1), so the
collector **buffers results by contig index and emits in reference-FASTA contig
order**, each contig already start-sorted. This is not cosmetic: the catalog is
**bgzip + tabix**, which *requires* coordinate-sorted input for the index to be
valid — so `(contig-order, start-order)` is a hard constraint, not a nicety. The
catalog *rows* are small, so buffering them in memory is fine; if it ever isn't,
spill to per-contig temp files (in the project-local `tmp/`, never `/tmp`) and
concatenate in order (GangSTR's pattern).

### 8.3 `--num-chroms-in-parallel` — a speed ⇄ RAM knob

The degree-of-parallelism cap on concurrent `trf-mod` workers. Two notes:

- **Work-stealing, not static partition.** Contig sizes are wildly uneven; a
  rayon iterator over contigs with a bounded pool lets a worker that finishes a
  small contig grab the next, instead of sitting idle while one chews chr1.
- **It is also a memory knob.** Each active worker holds its **whole contig
  sequence** resident (plus `trf-mod`'s process memory), so peak RAM ≈
  `num-chroms-in-parallel × (largest active contigs)`. On a big-chromosome plant
  genome that is real (a handful of large contigs at once = GBs of sequence).
  Per the project's memory thesis, document this as a **speed⇄RAM trade**, and
  don't default it blindly to core count when contigs are huge. (The collector's
  row buffer is *not* the cost — the in-flight contig sequences are.)

### 8.4 Determinism (an invariant)

The catalog is **byte-identical regardless of `--num-chroms-in-parallel`**.
Independent workers + emit-by-`(contig, start)` + a deterministic tie-break give
this for free; it is a **stated invariant** and a regression gate, matching the
project's output-determinism bar elsewhere.

### 8.5 The longest contig is a wall-time floor (deferred)

This model can't split a contig across workers, so once everything else is done
the build waits on chr1 alone. Acceptable for a once-per-genome step. If build
time ever hurts, the **deferred** escape hatch is **intra-contig chunking** — run
`trf-mod` on overlapping windows of a large contig and stitch at the cuts with a
halo rule for repeats spanning a boundary (the pattern the cohort path used for
DUST). Not built now; named so we know where to reach.

---

## 9. Accuracy harness (pointer — spec §3.4)

**Catalog accuracy** (owned by the test/validation layer, not this stage's code)
— randomized-sequence FP rate (any "SSR" in shuffled sequence is a false positive
→ calibrates the purity/score thresholds); simulation recall by
motif × copy × purity; boundary accuracy; recovery of **published tomato
capillary SSR markers** (the species anchor). These gate the parameter choices in
§3–§4; the build code just has to expose those parameters. *(The DUST-prefilter
recall study is gone with the prefilter, §6.)*

---

## 10. Open questions to work through together

*Decided: detector = lh3/TRF-mod binary, shell-out, **genome-wide, no sdust
prefilter** (§6), process-level parallelism, code-MIT / binary-AGPL (§2.3), and
**auto-launch via layered discovery** so the user runs only `ssr-catalog` (§2.4).
**Post-process order pinned** (§4): period≤6 → split (records the overlap
sequence) → recompute purity → filter; **merge delegated to TRF** (trust its
`-r`-default redundancy elimination); **no mappability**. **TRF-mod BED format
confirmed from source** — tab-separated, 0-based, `ctg` col 1; `purity_fraction`
recomputed from the tract (TRF `fracMatch` = sanity-check); motif = first-`period`
tract bases (§3). Distribution/container is **deferred** (§2.4). Still open:*

1. **`flank_bp` default** — a value now (revisited at Stage 1), or leave it
   purely a Stage-1-driven parameter with a placeholder default?
