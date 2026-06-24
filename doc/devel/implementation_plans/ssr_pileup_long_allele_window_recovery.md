# ssr-pileup — long-allele window recovery + truncation safeguard

**Status:** planned (branch `ssr-cohort`)
**Owner area:** Stage 1 (`ssr-pileup`) — read delimitation
**Relates to:** [ssr_pileup_mark2.md](../architecture/ssr_pileup_mark2.md) §3.2 (region extraction), the
`alignment.rs` delimiter, `footprint.rs::extract_region`.

## 1. Problem

The per-read realignment window ([`extract_region`](../../../src/ssr/pileup/footprint.rs)) is sized
to the locus's **reference** window (`ref_bytes` = `flank_bp` + ref-tract + `flank_bp`) and mapped to
read coordinates through the read's CIGAR. A read carrying an allele **longer** than the reference
needs its extra tract bases to arrive as a soft-clip or a CIGAR insertion; an aligner that instead
represents the long allele as **all-Match-with-mismatches** (a documented STR misalignment mode —
`footprint.rs` already notes mappers "get [indel placement] wrong inside repeats") leaves the window
ref-sized. The displaced tract then pushes the far flank *out* of the window, so the pair-HMM cannot
anchor it and silently **collapses the allele toward the reference length**.

Demonstrated (synthetic, all-M `CA×10` read vs a `CA×8` reference):

```
extracted window: GGGGGG CACACACACACACACACACA TT     ← right flank truncated to 2 bp
delimit_read    : tract = CA×8                       ← collapsed (lost 2 units), silently
full flanks     : tract = CA×10                      ← recovers when the flank is present
```

The delimiter itself is correct given full flanks; the defect is the window. Short alleles are
unaffected (a ref-sized window is always big enough for a shorter tract + full flanks).

## 2. Goal

Never silently under-call a long allele. Two-part defence (the agreed **option 1 + option 4**):

- **Recover (4):** when the window looks truncated, re-delimit in a window widened by the locus's own
  flank length and keep the recovered tract.
- **Account (1):** any read still flank-truncated after widening is **dropped, never tallied**, and
  every detection is **counted** in the `.ssr.psp` so a systematically-misbehaving mapper is visible
  in QC.

Non-goals: neighbour-locus-aware widening (the catalog *does* hold neighbours, but `flank_bp` is
provably safe and far simpler — revisit only if real data shows alleles needing more headroom);
changing the delimiter's HMM; a run-aborting hard error (too brittle for a genome-scale Stage 1 — the
defensive guarantee is "no silent corruption + visible counts", achieved by drop+count).

## 3. Design

### 3.1 Flank-completeness check (the detector)

After `delimit_read` returns `Region(r)` on a region of length `m`, with `left_len` / `right_len` the
locus flank lengths:

- `left_flank_bytes  = r.start`         (region-relative tract start)
- `right_flank_bytes = m - r.end`
- a side is **window-truncated** (more read bases exist there, we just didn't include them) iff
  `extract_region` was CIGAR-bounded on that side — derivable from the returned range: left bounded ⇔
  `region.start > 0`, right bounded ⇔ `region.end < read_len`.
- **suspicious** ⇔ a window-truncated side has fewer than its full flank length of flank bytes
  (`left_flank_bytes < left_len` on a bounded left, or `right_flank_bytes < right_len` on a bounded
  right).

Properties: no false positives for faithful/short alleles (they retain full flanks); no false
negatives for this failure mode (a truncated long allele must lose flank on ≥1 bounded side). Distinct
from the existing `BorderOffEnd` (tract runs to the *actual* read end — genuinely too-short read,
already counted). Cost: a handful of integer comparisons per read — O(1) on top of the Viterbi that
already ran.

### 3.2 Widened re-delimitation (the recovery)

On a suspicious read, extend the read window by the flank length on each side (the `flank_bp`-bounded
safeguard) and re-delimit:

```
widened = (region.start - left_len).max(0) .. (region.end + right_len).min(read_len)
```

Extending in **read** coordinates avoids re-walking the unreliable CIGAR. The margin is the locus's
own flank length (`left_flank().len()` / `right_flank().len()`); the catalog's bundle guarantee
(`bundle_threshold ≥ flank_bp` clean sequence beyond each flank) makes this provably safe — the
widened window cannot reach a neighbouring tract. Re-run `delimit_read` + the quality gate on the
widened slice, then re-apply the §3.1 check:

- flanks now complete → keep the recovered tract; count `n_widened`.
- still flank-short (allele longer than the read / the `flank_bp` headroom can hold) → **drop**;
  count `n_window_truncated`.

Only suspicious reads pay the second delimit; the faithful majority are byte-for-byte unchanged and
cost nothing extra.

### 3.3 New per-locus counters (the accounting)

Two scalar columns added to the `.ssr.psp` SSR record (pre-alpha schema, no back-compat needed):

- **`n_widened`** — reads whose tract was recovered by widening (the mapper likely mis-placed a long
  allele; the evidence was *kept*).
- **`n_window_truncated`** — reads dropped because the allele exceeded even the widened window (kept
  out of the tally so a collapsed allele can never reach the VCF; surfaced for QC).

"How many reads had the problem" = `n_widened + n_window_truncated`; "how many were discarded" =
`n_window_truncated`.

## 4. Touch points

1. **`src/ssr/pileup/footprint.rs`** — keep `extract_region`; add `widen_region(range, read_len,
   left_len, right_len) -> Range` (pure, read-coord margin), plus unit tests (faithful read
   unchanged; truncated long-allele region widens to include the full far flank; clamps at read
   bounds; left/right asymmetry).
2. **`src/ssr/pileup/alignment.rs`** — no HMM change. (The detector reads `Delimited::Region`'s range;
   the flank-completeness logic lives at the call site, which knows `left_len`/`right_len`/`read_len`.)
3. **`src/ssr/pileup/locus_tally.rs`** — extend `ReadObs` with `WidenedSequence(Box<[u8]>)` (tallied
   as an observed sequence **and** counted) and `WindowTruncated` (counted, not tallied); add
   `n_widened` / `n_window_truncated` to `SsrLocusObs`; count both in `tally`. Update its tests.
4. **`src/ssr/pileup/driver.rs`** — in `process_locus`, insert the detect → widen → re-check → classify
   logic between `delimit_read` and the `ReadObs` push; thread the two counts through
   `to_container_record`.
5. **`src/psp/registry_ssr.rs`** — add the two scalar columns (`SsrColumnKey`, the block SoA encode +
   decode, `SsrLocusRecord` fields), mirroring `n_border_off_end`. Round-trip test.
6. **Run summary** — include the two new counts wherever `ssr-pileup` reports per-run QC totals
   (alongside `n_low_quality` / `n_border_off_end`).
7. **Tests** —
   - footprint: `widen_region` units (above).
   - driver `process_locus`: an all-M long-allele read is *recovered* (observed = the full long tract,
     `n_widened == 1`); an allele longer than the read/margin is *dropped* (`n_window_truncated == 1`,
     not in `observed`); a faithful read is unchanged (`n_widened == 0`, identical tally to today).
   - end-to-end ([src/ssr/end_to_end_tests.rs](../../../src/ssr/end_to_end_tests.rs)): re-add a
     `CA×10` cohort with all-M reads and assert the long allele is now *recovered* in the VCF (and the
     `n_widened` count is set), turning the original collapse into a regression guard.
   - cross-thread byte-identity: the widen + counts are per-read pure → the `.ssr.psp` stays identical
     across thread counts (extend an existing thread-invariance test to a cohort that triggers a
     widen).

## 5. Edge cases / open details

- **Short flank but already at the read end** (`region.end == read_len`, can't widen): not the
  window-truncation case — the read genuinely lacks the flank. Keep current behaviour (delimiter's
  `BorderOffEnd` or the existing accept); do **not** fold into `n_window_truncated`. Note in code.
- **Margin sufficiency:** `flank_bp` recovers alleles up to ~`flank_bp` longer than the reference
  (~`flank_bp/period` extra units; with the default `flank_bp = 50`, ≈ 25 dinucleotide units). Longer
  alleles → `n_window_truncated` (visible, not silent). Acceptable; revisit with real data only if the
  truncation count is non-trivial.
- **Determinism:** all new logic is a pure function of one read + the locus → no effect on the
  cross-thread byte-identity contract.

## 6. Order of work

1. `widen_region` + footprint tests.
2. `ReadObs`/`tally` counters + tally tests.
3. `process_locus` detect→widen→classify; driver tests.
4. `registry_ssr` columns + round-trip; run-summary wiring.
5. End-to-end `CA×10` recovery test; thread-invariance check.
6. Gate (fmt / clippy `-D warnings` / doc / full test); update PROJECT_STATUS + an impl report.
