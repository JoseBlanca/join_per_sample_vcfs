# SSR Stage 1 — `ssr-pileup` Mark-2 implementation plan

**As of:** 2026-06-17, branch `ssr-pileup-mark2`. Builds the Stage-1 pipeline
designed in [architecture/ssr_pileup_mark2.md](../architecture/ssr_pileup_mark2.md)
over the model in [architecture/ssr_ladder_model.md](../architecture/ssr_ladder_model.md).
All Stage-1 design items (P1–P6, Q1–Q3) are settled; this is the build sequence.

**Working style:** one step at a time, each compiling (`./scripts/dev.sh cargo
check`) and unit-tested before the next; pause between steps. `ssr_mark1` keeps the
CLI working throughout and is deleted only at the cutover (step 7).

---

## 1. What gets copied almost verbatim (the read-I/O path)

These move from `src/ssr_mark1/` → `src/ssr/` with only `crate::ssr_mark1` →
`crate::ssr` repointing (and dropping the Mark-1-only items noted):

| source | lands as | change |
|---|---|---|
| `types.rs` (`Locus`, `Motif`) | `ssr/types.rs` | **drop** `Allele`, `NormalizedSeq` (no on/off-ladder); keep `Locus`/`Motif` + tests verbatim |
| `catalog/**` (Stage 0: TRF spawn/parse, postprocess, io, run) | `ssr/catalog/**` | **verbatim** (Stage 0 is model-agnostic, P5) |
| `pileup/fetch_reads.rs` (`Reservoir`, `SplitMix64`, `locus_seed`, `LocusReads`, `fetch_locus_reads`, `MAX_READS_PER_LOCUS`) | `ssr/pileup/fetch_reads.rs` | **verbatim** + fold in the footprint geometry (next row) |
| `pileup/triage.rs` → only `reaches_locus`, `read_footprint`, `brackets`, `ref_to_read`, `extract_region` (+ their tests) | **into** `ssr/pileup/fetch_reads.rs` | keep the **footprint geometry + admission gate**; **drop** `triage_read`, `SpanningRead`, `TriageResult`, `find_longest_stretch`, `ProbeHit` (Mark-1 rung-window machinery) |
| `pileup/driver.rs` skeleton (input loading, `AlignmentFile`/`WorkerReader`, batched `par_chunks`, `build_ssr_writer_header`, atomic temp+rename, `SsrPileupConfig`, `SsrPileupError`) | `ssr/pileup/driver.rs` | skeleton kept; **swap** `process_locus` body, the record type, `to_container_record`, `LocusScratch`, header params (§3) |
| the **shared** `src/bam/**` readers (`load_pileup_inputs`, `AlignmentFile`, `WorkerReader`, `FilterCounts`, `SegmentReadFilter`) | — | used directly, **not** copied (already shared) |
| `pileup/pair_hmm.rs` → only `HmmModel`, `EMISSION_LN`, `INS_EMIT_LN`, transition constants, `ln_*` helpers, the grow-and-keep scratch idea | **into** `ssr/pileup/alignment.rs` | reuse the **emission/transition model**; the forward (`forward`, `score_candidates`, prefix-seam) is **dropped** — replaced by Viterbi+traceback (§2) |

**Dropped entirely (Mark-1-only):** `candidate_generation.rs`, `read_analysis.rs`,
`locus_record.rs` (`aggregate`/`prune_and_renormalize`), `count_repeats.rs`,
`bench_harness.rs` (re-derive a Mark-2 bench later if needed), the Mark-1 columns in
`psp/registry_ssr.rs`.

---

## 2. New code — `src/ssr/pileup/alignment.rs` (the substantive piece)

The per-Q Viterbi+traceback delimiter (Q2/P2) + the first-quartile quality gate
(Q1/P1). Reuses `HmmModel`/`EMISSION_LN` from Mark-1's pair-HMM; adds the max-path
DP with a full backpointer matrix and the traceback.

```rust
/// First-quartile repeat-region quality cutoff (Q1/P1) — drop a read whose
/// repeat-region base-quality first quartile is below this. First value; calibrate.
pub(crate) const MIN_REGION_Q1: u8 = 15;

/// Per-worker Viterbi scratch: rolling score rows (max, not log-sum-exp) + a FULL
/// backpointer matrix for traceback (unlike the forward's two-row-only scratch).
/// Grow-and-keep, like the Mark-1 forward scratch.
pub(crate) struct ViterbiScratch {
    prev: Vec<[f64; 3]>,            // [M,I,D] best scores, previous read row
    cur:  Vec<[f64; 3]>,            //                       current read row
    back: Vec<[u8; 3]>,            // (rows+1)*(cols+1) cells; per state, the
                                    // predecessor state (M/I/D) that won the max
    cols: usize,                    // frame length + 1 (row stride)
}

/// Outcome of delimiting one read region against a locus.
pub(crate) enum Delimited {
    /// Repeat region = `region[range]` (read-coordinate, relative to the slice
    /// handed in). Both flank junctions landed inside the read.
    Region(std::ops::Range<usize>),
    /// A flank ran off the read end (allele >= read length, spec §1.4) — counted,
    /// not used.
    BorderOffEnd,
}

/// Align the read region (the locus-window slice, clips included — extracted by
/// `fetch_reads::extract_region`) against the locus reference frame
/// `left_flank + ref_tract + right_flank` with a per-Q Viterbi + traceback, and
/// return the repeat-region span between the two flank-junction columns.
///
/// Determinism (Q2): max tie-break `Match > Deletion > Insertion`; a junction
/// indel is assigned to the preceding (5') block (left junction -> flank, right
/// junction -> repeat).
pub(crate) fn delimit_read(
    region_seq:  &[u8],
    region_qual: &[u8],
    locus: &Locus,
    model: &HmmModel,
    scratch: &mut ViterbiScratch,
) -> Delimited;

/// First quartile of `quals` >= `threshold`. The Q1 gate (computed on the
/// delimited repeat region's quals).
pub(crate) fn passes_quality_gate(quals: &[u8], threshold: u8) -> bool;
```

Internal helpers (not exported): `fill_row_viterbi` (max + record backpointer,
mirrors Mark-1's `fill_row`), `traceback` (walk `back` from the argmax final cell,
emit the flank-junction read positions), and `frame(locus) -> (left_len, tract,
right_len)` to locate the junction columns (`left_flank.len()`, `hap.len() -
right_flank.len()`).

---

## 3. New code — `src/ssr/pileup/locus_tally.rs`

```rust
/// Per-read outcome the tally folds (the Mark-2 analogue of Mark-1 `ReadOutcome`).
pub(crate) enum ReadObs {
    Sequence(Box<[u8]>),  // usable, quality-passing repeat-region bytes
    LowQuality,           // dropped by the Q1 gate
    BorderOffEnd,         // a flank ran off the read end
}

/// One sample's observed evidence at one locus (in-memory, chrom-NAME-keyed,
/// 0-based; the driver's adapter shifts to the container's id-keyed/1-based form).
pub(crate) struct SsrLocusObs {
    pub(crate) chrom: Box<str>,
    pub(crate) start: u32,
    pub(crate) end: u32,
    // QC scalars (P3 lean set; nothing derivable)
    pub(crate) depth: u32,
    pub(crate) n_filtered: u32,
    pub(crate) mapped_reads: u32,
    pub(crate) n_low_quality: u32,
    pub(crate) n_border_off_end: u32,
    /// Distinct repeat-region sequences -> observation count, SORTED BY BYTES
    /// (deterministic storage order — the byte-identity invariant).
    pub(crate) observed: Vec<(Box<[u8]>, u32)>,
}

/// Fold per-read outcomes + the fetch-pass QC into the locus record: tally
/// sequences in a `HashMap`, count the two dropout classes, then emit `observed`
/// sorted by bytes. `depth`/`n_filtered`/`mapped_reads` come from `qc`.
pub(crate) fn tally(locus: &Locus, outcomes: &[ReadObs], qc: QcCounts) -> SsrLocusObs;
```

`QcCounts` is reused from `fetch_reads`/the driver (`depth`, `n_filtered`,
`mapped_reads`); `n_low_quality`/`n_border_off_end` are derived from `outcomes`.

---

## 4. New code — `src/psp/registry_ssr.rs` (rewrite in place, no version bump, P4)

Mark-2 schema: per-locus QC scalars + a per-observation CSR keyed by sequence.
Mirrors the Mark-1 `SsrKind`/`SsrBlock`/`SsrDecoder` structure; only the columns
change. Rides the existing `Bytes`+CSR+varint codecs (no new wire code).

```rust
pub enum SsrColumnKey {           // per-record [0x01..0x0F], per-obs [0x10..]
    DeltaStart = 0x01, Span = 0x02, NObs = 0x03,     // NObs = grouping count
    Depth = 0x04, NFiltered = 0x05, MappedReads = 0x06,
    NLowQuality = 0x07, NBorderOffEnd = 0x08,
    ObsCount = 0x10,              // List<Varint> CSR, per locus
    ObsSeqLen = 0x11,            // List<Varint> CSR, per locus (drives ObsSeq)
    ObsSeq = 0x12,               // Bytes { length_column: "obs-seq-len" }
}

/// Container-form record (chrom_id-keyed, 1-based) — the Mark-2 mirror of
/// SsrLocusObs. n_obs is derived from observed.len() (grouping count).
pub struct SsrLocusRecord {
    pub chrom_id: u32, pub start: u32, pub end: u32,
    pub depth: u32, pub n_filtered: u32, pub mapped_reads: u32,
    pub n_low_quality: u32, pub n_border_off_end: u32,
    pub observed: Vec<(Box<[u8]>, u32)>,   // (sequence, count), sorted by bytes
}
```

`SsrBlock` accumulates: the per-locus scalar vecs + the obs CSR slabs
(`obs_offsets: Vec<u32>` shared by count/seq-len, `obs_count_data: Vec<u64>`,
`obs_seq_len_data: Vec<u32>`, `obs_seq_bytes: Vec<u8>`). `SsrDecoder` mirrors and
reassembles `observed`. Keep the Mark-1 structural checks generalised: `obs_offsets
.len() == n_records'+1`-style, `sum(obs-seq-len) == obs-seq.len()`, each sequence
`<= MAX_ALLELE_SEQ_LEN`, and the `sum(n-obs) == n_total` grouping check.

---

## 5. Driver changes — `src/ssr/pileup/driver.rs`

```rust
struct LocusScratch {
    region_seq: Vec<u8>,          // reused extract_region output buffers
    region_qual: Vec<u8>,
    viterbi: ViterbiScratch,
    tally: HashMap<Box<[u8]>, u32>,   // cleared per locus
}

fn process_locus(readers, locus, cap, model, scratch) -> SsrLocusObs {
    let fetched = fetch_locus_reads(readers, locus, cap)?;     // REUSE
    let outcomes = fetched.reads.iter().map(|read| {
        let region = extract_region(read, locus);              // REUSE (footprint)
        match delimit_read(&read.seq[region], &read.qual[region], locus, model, &mut scratch.viterbi) {
            Delimited::Region(r) => {
                let q = &read.qual[region][r.clone()];
                if passes_quality_gate(q, MIN_REGION_Q1) { ReadObs::Sequence(read.seq[region][r].into()) }
                else { ReadObs::LowQuality }
            }
            Delimited::BorderOffEnd => ReadObs::BorderOffEnd,
        }
    }).collect();
    tally(locus, &outcomes, qc_counts(&fetched))               // REUSE qc_counts
}
```

Plus: `to_container_record(SsrLocusObs, name_to_id)` (name→id, +1 coord shift,
drop `n_obs` — derived); `build_ssr_writer_header` params drop `window`, add
`quality_q1_threshold`. The `run` batched-`par_chunks` loop is unchanged.

---

## 6. Build order (each step compiles + is tested before the next)

0. **Scaffold** `src/ssr/` + `pub mod ssr;` in `lib.rs` (empty `mod.rs`). green.
1. **types + Stage 0**: copy `types.rs` (drop `Allele`/`NormalizedSeq`) + `catalog/`
   verbatim. Tests: the copied catalog/types tests pass under `crate::ssr`. green.
2. **fetch_reads**: copy `fetch_reads.rs` + fold in the footprint geometry +
   `reaches_locus`/`extract_region`. Tests: reservoir determinism, `fetch_locus_reads`
   BAM fixture, `extract_region` (incl. clip-extended). green.
3. **alignment**: `alignment.rs` — Viterbi+traceback delimiter + quality gate.
   Tests (anti-tautology): clean read → repeat region = the tract; interior-indel
   read → indel-bearing region extracted; a read with a flank off the end →
   `BorderOffEnd`; tie-break determinism (two equal alignments → identical span);
   gate at Q1=15. green.
4. **locus_tally**: `SsrLocusObs` + `tally`. Tests: sequence histogram, sorted
   byte order, the two dropout counts, shuffle-invariance of the record. green.
5. **registry_ssr (Mark-2)**: rewrite the schema. Tests: multi-block round-trip
   (incl. a locus with zero observations, a many-sequence locus), structural
   rejections, wrong-kind poison (mirror the Mark-1 tests). green.
6. **driver**: swap `process_locus`/`to_container_record`/`LocusScratch`/header;
   keep the `run` skeleton. Tests: catalog→reference→BAM→`.ssr.psp` end-to-end
   (one clean spanning read → one observed sequence with the right bytes) +
   thread-count byte-identity. green.
7. **cutover + cleanup**: point the `ssr-catalog`/`ssr-pileup` CLI wrappers at
   `crate::ssr`; **delete `src/ssr_mark1/`**; remove any Mark-1 `registry_ssr`
   remnants. `cargo check --all-targets` + full `cargo test` green; update
   `ssr_stage1_remaining.md`.

---

## 7. Invariants to hold at every step

- **Determinism / byte-identity across `--threads`** (the v1 gate): per-locus
  reservoir seed (reused), per-locus atomic scoring, catalog-order writes, the
  Viterbi tie-break (§2), and the sorted `observed` order (§3).
- **Anti-tautology** (spec §7): delimiter/gate tests use data whose error pattern
  is defined independently of the code under test; the simulator-level end-to-end
  test (read/BAM level) lands with step 6 (reuse the Mark-1 BAM-fixture helpers).
- **Catalog binding**: the `.ssr.psp` header records `reference_md5` + `flank_bp`
  from the catalog (reused header build).

---

## 8. Not in this plan (Stage 2)

Candidate assembly from the pooled observed sequences (model S1), stutter
reachability between sequences (S2), the flat-error likelihood HMM (S3), EM, and
the VCF — all `ssr-call` (Stage 2), designed just-in-time after Stage 1 lands.
