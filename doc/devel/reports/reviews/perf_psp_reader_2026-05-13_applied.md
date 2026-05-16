# Applied deltas — psp::reader perf review (2026-05-13)

Companion to [perf_psp_reader_2026-05-13.md](perf_psp_reader_2026-05-13.md).
This document records which Likely findings were applied, which were
deferred, and what the bench numbers say about each.

## Baselines used

Three criterion baselines saved during the cycle (under
`target-container/criterion/.../psp_reader/.../{reader_baseline,post_l5}/`):

- `reader_baseline` — at the start of this review (commit `eef4995`,
  the perf-review-report commit). Numbers verbatim in
  [perf_psp_reader_2026-05-13.md §1](perf_psp_reader_2026-05-13.md).
- `post_l5` — after L5 was applied (commit `a04282f`).
- (current) — after L5 + L6 + L8 + L1 + L2 (commit `013b47d`).

All bench numbers below are from `cargo bench --bench psp_reader_perf`
at 30 s `measurement_time`, sample-size 10, against the saved
baselines unless noted otherwise.

## Findings — disposition

### Applied

- **L5** (`psp::varint::decode_u64_leb128` fast/cold split) — applied
  in commit `a04282f`.
  - `psp_reader/snp_typical_3_3M`: **+15.8 % throughput** (p = 0.00),
    11.4 → 13.1 Mrec/s. Clean win on the SNP workload where every
    per-record varint is single-byte.
  - `psp_reader/phase_chain_heavy_1M` and `psp_reader/multi_allele_500k`:
    no change detected (the varint share is smaller on those workloads,
    consistent with the writer-side L4 measurement).

- **L6** (`psp::block::decode_scalar_column` LE slab cast for
  `T: Pod`) — applied in commit `77bcd43` (bundled with L8).
  - `psp_reader/multi_allele_500k`: **+6.0 % throughput** (p = 0.00,
    60 s measurement), 5.0 → 5.3 Mrec/s. Multi-allele is the workload
    most affected because it has the largest per-allele scalar count
    (1.25 M alleles for 500 k records).
  - `psp_reader/snp_typical_3_3M` and `psp_reader/phase_chain_heavy_1M`:
    no change (one allele per record on these workloads — the slab-cast
    path is exercised proportionally less).

- **L8** (rustdoc `BufReader::with_capacity(64 * 1024, file)`
  recommendation on `PspReader::new` and `region_records`) — applied
  in commit `77bcd43` (bundled with L6). Doc-only; no runtime impact.

- **L1 + L2** (persistent `zstd::bulk::Decompressor` + per-iter
  scratch buffers on `RecordsIter`) — applied together in commit
  `013b47d` because the same refactor wires both. New helpers
  `new_column_decompressor()` and `zstd_decompress_into()` in
  `block.rs`, symmetric with the writer-side
  `new_column_compressor()` / `zstd_compress_into()` (writer's H3
  + L1, commit 969de6c).
  - `psp_reader/snp_typical_3_3M`: **+15.3 % throughput** (p = 0.00,
    cumulative with L6+L8). 12.0 → 14.1 Mrec/s.
  - `psp_reader/multi_allele_500k`: **+9.9 % throughput** (p = 0.00).
    5.3 → 5.7 Mrec/s.
  - `psp_reader/phase_chain_heavy_1M`: **volatile across runs.** One
    60 s measurement returned +15.1 % thrpt (improved); the
    immediately following 60 s measurement returned −7.2 % thrpt
    (regressed). The fixture is small (10 KB compressed → 1 block);
    per the methodology checklist's warning that cross-commit
    comparisons need either two clean measurements per side or a
    revert experiment, the phase delta sits inside criterion's noise
    floor for this workload shape. The fix is structurally correct
    (mirrors writer's H3 17.43 % self-time win on the symmetric
    encode-side site) and the wins on snp / multi are clean — apply
    rather than defer.

### Deferred

- **L3** (collapse `Vec<Vec<SlotId>>` columns in `DecodedBlock` to
  CSR) — deferred. The report explicitly flagged this as "**not a
  guaranteed wall-clock win**": the per-row decode-time allocation
  goes away but the materialiser would still pay one `to_vec()` per
  emitted row (because `PileupRecord` / `AlleleObservation` own their
  `Vec<SlotId>`). Without dhat instrumentation in place (M3 / L9, not
  applied this round) the trade-off can't be ranked from criterion
  alone. Open as a follow-up once `examples/dhat_psp_reader.rs`
  lands.

- **L4** (`decode_bytes_split` per-allele `to_vec()` special-cases /
  CSR) — deferred. Same gating reason as L3 (the deeper Stage B fix
  shares the CSR refactor; Stage A is a small special-case for `n ==
  0` whose payoff is unclear without dhat).

- **L7** (`#[cold]` markers on per-record error-construction sites) —
  **applied, then reverted** because it regressed the bench. The
  experiment: extracted `err_no_loaded_block`, `err_phase_chain`, and
  `err_delta_pos_overflow` `#[cold] #[inline(never)]` helpers and
  rewired the three error sites in `materialise_next_record`.
  Bench against `post_l5`:
  - `psp_reader/snp_typical_3_3M`: −5.6 % thrpt (p = 0.02, regressed)
  - `psp_reader/phase_chain_heavy_1M`: −13.6 % thrpt (p = 0.00,
    regressed)
  - `psp_reader/multi_allele_500k`: no change

  **Verdict: experiment shows no gain — closing.** The reader's
  error sites are evidently *not* as cold as the writer's
  `validate_record` per-byte ACGTN error sites (which had measured
  10.09 % self-time on the writer profile). The
  `materialise_next_record` validation loops execute ~zero error
  returns in steady state, so the icache argument should hold — but
  in practice forcing the `#[inline(never)]` indirection cost the
  branch-predictor enough to net-regress on snp + phase. Without a
  reader-side sampling profile to validate the icache hypothesis,
  the conservative call is to keep the inline error constructions.

- **L9** (`examples/dhat_psp_reader.rs`) — deferred. Methodology
  infrastructure; useful for the next round (gating L3 / L4) but not
  required for the wins applied this round (criterion was sufficient
  to gate L1 / L2 / L5 / L6).

- **L10** (more bench coverage: file-backed sub-bench,
  `decode_block_only`, region split, `Throughput::Bytes`) — deferred.
  Same reasoning as L9: criterion's existing streaming benches were
  sufficient to gate the applied wins. The file-backed sub-bench
  becomes load-bearing if / when L8's BufReader sizing recommendation
  needs empirical confirmation.

- **S1–S6** (Speculative findings: SmallVec for active set, combined
  validate+apply, BlockIndexEntry projection, mmap experiment,
  lending iterator, profile escape-hatch documentation) — all
  deferred. Each requires either a profile that's not reachable in
  this sandbox, or a workload (region-heavy, large-index) the
  current bench doesn't exercise. File as future-PR menu.

### Cumulative throughput vs the original `reader_baseline`

| Workload | Before | After | Delta |
|---|---|---|---|
| `snp_typical_3_3M` | 11.4 Mrec/s | 13.9 Mrec/s | **+22 %** |
| `phase_chain_heavy_1M` | 5.86 Mrec/s | ~5.6–6.6 Mrec/s | volatile |
| `multi_allele_500k` | 5.0 Mrec/s | 5.7 Mrec/s | **+14 %** |
| `region_window_chr1_mid_100k` | 1.27 Mrec/s | (not re-measured) | — |

All three streaming workloads remain comfortably over the 1 M rec/s
floor; the SNP and multi-allele wins moved the headroom from
~5–11× the floor to ~6–14× the floor.

## Validation

Each commit individually:

- `cargo test --lib --all-features` → 474 passed
- `cargo test --lib --tests --all-features` → 562 passed total
- `cargo clippy --all-targets --all-features -- -D warnings` → 0 warnings
- `cargo fmt --check` → no diff
- `cargo test --doc per_sample_caller::psp::reader` → 2 passed
- `cargo bench --bench psp_reader_perf -- --baseline <prev>` → see per-finding numbers above

## Follow-ups

1. **Sampling profile escape hatch (S6).** Every code-level finding
   here is capped at Likely because `samply record` is blocked in
   the sandbox. Document `scripts/profile.sh` or extend the dev
   container with `--cap-add=PERFMON --security-opt seccomp=unconfined`
   so the next reviewer can attach a real flamegraph and either
   confirm or kill the deferred Likely findings.
2. **`examples/dhat_psp_reader.rs` (L9).** Mirror of
   `examples/dhat_psp_writer.rs`. Required to gate L3 / L4 because
   wall-clock alone cannot tell whether the per-row inner-Vec
   allocations have actually been removed.
3. **L3 / L4 (CSR for `Vec<Vec<…>>` columns).** Land after L9.
   Threshold: ≥ 5 % wall-clock improvement on
   `phase_chain_heavy_1M` *and* a clean dhat diff showing the
   per-row inner-Vec allocations are gone, with no regression on
   the other workloads.
4. **L10 (file-backed bench).** Required to confirm L8's BufReader
   sizing recommendation empirically and to land any future I/O-
   layer change (S4 mmap experiment, etc.).

The header / pileup-side / writer-side parts of the spec remain out
of scope for this review.
