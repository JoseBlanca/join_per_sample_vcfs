# SSR Stage 1 — `ssr-pileup` Mark-2 (architecture draft)

**Status:** draft, 2026-06-17, branch `ssr-pileup-mark2`. The Stage-1 architecture
for the **Mark-2 empirical-candidate model**
([ssr_ladder_model.md](ssr_ladder_model.md)). Mark-2 makes candidate alleles
*observed sequences* (not reference rungs) and defers all likelihood to Stage 2,
so Stage 1's job collapses to: **fetch reads → delimit each read's repeat region →
quality-gate it → tally observed sequences per locus → write**. No HMM scoring, no
on/off-ladder, no candidate ladder here.

It **reuses the Mark-1 read-I/O layer almost verbatim** (the part that pulls reads
from BAM/CRAM) and swaps out the per-read compute and the storage. Supersedes the
Mark-1 [ssr_pileup.md](ssr_pileup.md) for the compute/storage; the read-I/O
sections of that doc remain valid.

Model decisions this builds on (from [ssr_ladder_model.md](ssr_ladder_model.md)):
Q1 first-quartile quality gate, Q2 Option-A delimitation (one Viterbi+traceback vs
the reference frame; tie-break `Match > Deletion > Insertion`; junction indel →
preceding 5′ block), Q3 store per-locus `(sequence, count)` with no quality.

---

## 1. Pipeline overview

```
 one sample's            ┌──────────────── per locus (one work unit) ───────────────┐
 BAM/CRAM(s) ─┐          │ fetch_locus_reads  (REUSE: reservoir + reach gate)        │
 reference ───┼─► driver │   → delimit each read   (NEW: Viterbi+traceback, Q2)      │ ─► sample.ssr.psp
 .ssr.catalog ┘  (loci)  │   → quality-gate        (NEW: first-quartile, Q1)         │   (per-locus
                         │   → tally sequences     (NEW: HashMap<seq,count>, Q3)     │    sequence+count,
                         │   → locus record        (NEW: distinct (seq,count)+QC)    │    new schema)
                         └───────────────────────────────────────────────────────────┘
```

One locus = one self-contained unit (as in Mark-1), so the driver stays
embarrassingly parallel across loci with byte-identical output at any thread count.

---

## 2. Reuse map — lifted from `ssr_mark1` almost verbatim

The user's call: the read-fetching path copies over with little change. Concretely,
two layers:

**(a) The shared BAM/CRAM machinery — used directly, not copied** (it already lives
under `src/bam/`, not under `ssr`): `load_pileup_inputs` / `build_fasta_repository`
/ `PileupInputs`, `AlignmentFile` + `WorkerReader` (the decode-caching CRAM /
pooled BAM read source), the segment query (`fetch_mapped_reads`), `FilterCounts`,
`SegmentReadFilter`. Mark-2 calls these unchanged.

**(b) The SSR fetch wrapper + driver skeleton — copy `ssr_mark1` → `ssr`, repoint
`crate::ssr_mark1` → `crate::ssr`, swap the worker body.** What copies essentially
verbatim:

| from `ssr_mark1` | reuse in Mark-2 |
|---|---|
| `fetch_reads.rs`: `Reservoir<T>` (Algorithm R + SplitMix64), `locus_seed` (FNV-1a), `MAX_READS_PER_LOCUS` | **verbatim** — the deterministic per-locus depth cap |
| `fetch_reads.rs`: `fetch_locus_reads` + `LocusReads { reads, yielded, filtered }` | **verbatim** — segment-query the embedded window, admission-gate, reservoir-cap |
| `triage.rs`: `reaches_locus` (the cheap "could this read span the locus?" admission gate — a footprint/clip test, no alignment) + the footprint helpers it needs | **fold into `fetch_reads.rs`** — it *is* the reservoir admission gate, so it lives with fetching, not in a separate module. The footprint geometry is shared with the delimiter (§3). The rest of `triage.rs` (`triage_read`, `SpanningRead`, the rung window centre, `find_longest_stretch`) is Mark-1-specific and **dropped** — Mark-2 delimits by alignment (§3); `find_longest_stretch` may return later only as the P6 DP-bounding optimisation. |
| `driver.rs`: input loading, per-file `AlignmentFile`, the **batched `par_chunks`** loop (warm decode cache per chunk, `LOCUS_BATCH`/`CHUNKS_PER_THREAD`/`MIN_FETCH_CHUNK`), `build_ssr_writer_header`, atomic temp-write+rename, the name→chrom_id adapter *shape*, `SsrPileupConfig` | **skeleton verbatim**; the per-locus body (`process_locus`) and the record type change (§3–§5), and the header params change (drop `window`; add the quality threshold) |
| `driver.rs`: determinism approach (per-locus reservoir seed + atomic per-locus scoring + catalog-order writes) | **verbatim** — the §7 invariant |

What does **not** carry over: `pair_hmm.rs` `forward`/`score_candidates` (Mark-1
likelihood scoring — Mark-2 has no Stage-1 scoring), `candidate_generation.rs`
(rungs/off-ladder — gone), `locus_record.rs` `aggregate`/`prune_and_renormalize`
(replaced by the sequence tally), and the Mark-1 `Allele`/`NormalizedSeq` types
(no on/off-ladder; a repeat-region allele is just a byte sequence). The HMM
*emission/transition model + scratch discipline* from `pair_hmm.rs` **is** reused —
by the new delimiter (§3).

---

## 3. New: per-read delimitation (`alignment.rs`)

Replaces Mark-1 triage's region+observed_count extraction. Per read (Q2):

- **Align the whole read against the single reference frame** `left_flank +
  ref_tract + right_flank` (from `Locus`) with a pair-HMM in **Viterbi (max-path) +
  traceback** mode, unbanded. *New code* — Mark-1's pair-HMM is forward-only — but
  it reuses Mark-1's emission table (per-Q `EMISSION_LN`) and transition model
  (`HmmModel`) and the grow-and-keep scratch ethos; the addition is a backpointer
  matrix + a traceback. (Per-Q emissions are fine here: the read's quals are still
  in hand at Stage 1 — we only *store* without them.)
- **Read the boundary off the traceback**: the flank-junction haplotype columns are
  known a-priori (`left_flank.len()`, `hap.len() − right_flank.len()`); the read
  positions aligned there split the read into `left_border | repeat_region |
  right_border`. The repeat-region read bases are the observed sequence.
- **Determinism**: fixed traceback tie-break `Match > Deletion > Insertion`;
  junction indel assigned to the **preceding (5′) block** (left-junction indel →
  flank, right-junction → repeat). This is the cross-sample identity guarantee
  (same molecule → same bytes).
- **Coverage classification** (replaces Spanning/Flanking/InRepeat): a read whose
  traceback places **both** flank junctions inside the read is *usable*; if a flank
  runs off the read end (allele ≥ read length, spec §1.4) the read is **counted for
  QC, not used**.

> Emission model: **per-Q** (decided, P2) — reuse Mark-1's table; quals are in hand
> at Stage 1. DP: **unbanded for v1** (decided, P6) — the per-locus reservoir depth
> cap (reused from v1) bounds the work to `cap × read × frame`, so banding is a
> deferred, measured optimisation.

---

## 4. New: quality gate (alongside the delimiter in `alignment.rs`)

After delimitation (Q1): compute the **first quartile of the repeat-region base
qualities** and **drop the read** if it is below a threshold. Survivors are treated
as uniform-quality, which is what licenses storing sequences without quals (Q3).
Drops are tallied (a new QC scalar). **First value: Phred 15** (a starting point —
must be calibrated on real data; §8-P1).

---

## 5. New: per-locus tally + record (`locus_tally.rs`)

Replaces `locus_record.rs`. For a locus:

- Tally surviving reads' repeat-region sequences into a `HashMap<Box<[u8]>, u32>`
  (sequence → count). Order-independent (commutative counts), so it composes with
  the per-locus-atomic threading.
- Emit the **Mark-2 locus record**: the distinct `(sequence, count)` pairs **sorted
  by sequence bytes** (deterministic storage order — the byte-identity requirement,
  §7) + the QC scalars.

QC scalars (DECIDED, P3) — the **lean set, nothing derivable**: `depth`,
`mapped_reads`, `n_filtered` (carried over from the fetch pass) + the two Mark-2
dropout counts `n_low_quality` (quality-gate drops) and `n_border_off_end` (a flank
ran off the read end). **Not stored** because they're derivable: `n_usable` = sum
of the observed counts; Mark-1's `n_spanning` likewise. `n_flanking`/`n_frr` are
gone (Mark-1 classifications with no Mark-2 analogue beyond `n_border_off_end`).

```rust
// shape sketch — not final
struct SsrLocusObs {
    chrom: Box<str>, start: u32, end: u32,
    // QC scalars: depth, mapped_reads, n_filtered, n_low_quality, n_border_off_end
    observed: Vec<(Box<[u8]>, u32)>,   // distinct repeat-region seq -> count, sorted by bytes
}
```

---

## 6. New: container schema (Mark-2 `.ssr.psp`)

A **new schema**, replacing Mark-1's `registry_ssr` (`amb-lengths`/`amb-logliks`
CSR). Per locus, all-CSR over *observed sequences*:

- `obs-seq` — `Bytes { length_column: "obs-seq-len" }`: the distinct repeat-region
  sequences (verbatim bytes), reusing the container's existing `Bytes` codec (the
  SNP `allele-seq` path).
- `obs-seq-len` — per-sequence byte length (drives the `Bytes` chunking).
- `obs-count` — per-sequence observation count (varint), parallel to the sequences.
- a per-locus grouping count (how many distinct sequences this locus has — the
  `n-alleles`-style group key) + the QC scalars (`delta-start`, `span`, `depth`,
  `n-filtered`, `mapped-reads`, …).

Rides existing codecs (`Bytes` + CSR + varint/scalar) — no new wire code, mirroring
the Mark-1 schema's "a table + a record mapping" approach. **No `schema_version`
bump** (P4): pre-alpha, no existing `.ssr.psp` to stay compatible with, so the
Mark-2 columns simply replace the Mark-1 ones in place — rewrite `registry_ssr`, no
`_v2`.

zstd dedups identical sequences across loci/samples, so storing verbatim is cheap —
consistent with the Mark-2 "counts carry the signal, let zstd handle redundancy"
choice.

---

## 7. Determinism invariant (carried over verbatim)

`.ssr.psp` must be **byte-identical across `--threads ∈ {1, N}`** (the Mark-1 gate).
Mark-2 preserves every lever: per-locus reservoir seed (`locus_seed`); each locus
scored atomically on one thread; catalog-order writes; **plus** two Mark-2
additions that must also be deterministic — the delimitation tie-break (§3) and the
**sorted** storage order of distinct sequences (§5). The reservoir's
order-sensitivity caveat (offer reads in a fixed total order) carries over unchanged.

---

## 8. Module layout + open items

**Proposed `src/ssr/` (Mark-2) layout** — copy the reusable files, add the new ones:

```
src/ssr/
  types.rs          # Locus, Motif  (reuse from ssr_mark1; DROP Allele/NormalizedSeq)
  pileup/
    fetch_reads.rs  # COPY: Reservoir, locus_seed, fetch_locus_reads + the
                    #       reaches_locus admission gate & footprint helpers
                    #       (folded in from ssr_mark1 triage)
    alignment.rs    # NEW: viterbi_align + traceback + boundary read-off (§3) +
                    #      the first-quartile quality gate (§4); reuses the
                    #      footprint geometry from fetch_reads for region extraction
    locus_tally.rs  # NEW: sequence tally -> SsrLocusObs (§5)
    driver.rs       # COPY skeleton; swap process_locus body + record type
  catalog/          # COPY verbatim from ssr_mark1 (Stage 0 is model-agnostic, P5)
src/psp/registry_ssr.rs    # REWRITE in place to the Mark-2 schema, no version bump (§6, P4)
```

**Open items (Stage-1-specific; Stage-2 S1/S2/S3 are out of scope here):**

- **P1 — quality threshold value** (Q1): **first value Phred 15** (decided
  2026-06-17 as a starting point); must be calibrated on real data, and measure the
  long-allele dropout.
- **P2 — delimiter emission model. DECIDED (2026-06-17): per-Q.** The read's quals
  are in hand at delimitation, so the Viterbi uses Mark-1's per-Q emission table —
  more robust boundary placement (a low-Q base near a junction won't shift it),
  free (table already exists). The uniform-quality assumption applies only to
  *storage* + *Stage-2 scoring*, not to the Stage-1 delimiter.
- **P3 — QC scalar set. DECIDED (2026-06-17): the lean set, nothing derivable** —
  `depth`, `mapped_reads`, `n_filtered`, `n_low_quality`, `n_border_off_end`.
  `n_usable`/`n_spanning` are derivable (sum of observed counts) → not stored;
  `n_flanking`/`n_frr` retired.
- **P4 — schema rollout. DECIDED (2026-06-17): no `schema_version` bump** (pre-alpha,
  no existing files); the Mark-2 columns replace the Mark-1 ones in place — rewrite
  `registry_ssr`, no `_v2`.
- **P5 — catalog (Stage 0). DECIDED (2026-06-17): copy `ssr_mark1/catalog` verbatim**
  into `ssr/` as the first step (Stage 0 is model-agnostic — no Mark-1/Mark-2
  difference).
- **P6 — unbanded DP cost. DECIDED (2026-06-17): unbanded for v1; defer banding
  behind a measurement.** The per-individual **depth cap** — the same one-pass
  reservoir from v1 (`Reservoir` / `MAX_READS_PER_LOCUS`, deterministic per-locus
  seed), reused verbatim — already bounds the per-locus DP work to
  `cap × (read × frame)` regardless of locus depth, so the unbanded full DP has a
  hard per-locus ceiling. If profiling still shows the delimiter binds, the band
  lever (the content pre-probe, already in `ssr_mark1`) is the deferred optimisation.

**Not in Stage 1** (Stage 2, per the model doc): candidate assembly (S1), stutter
reachability between sequences (S2), the likelihood HMM (S3), EM, VCF.
