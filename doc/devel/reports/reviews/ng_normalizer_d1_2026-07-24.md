# Code Review: ng normalizer — D1 (the differ-at-all screen + stress generator)

**Date:** 2026-07-24
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** two new example binaries — `examples/ng_normalizer_screen.rs`, `examples/ng_synthesize_stress_reads.rs`
**Status:** Approve-with-changes

---

### 1. Scope
- D1: the differ-at-all screen (plan step D1; spec §6/§10.3) + a synthetic stress-read generator (added at the owner's request to show the screen is discriminating).
- Categories dispatched (6, across two parallel agents): reliability, errors, naming, idiomatic, smells (+ the intent/logic trace). module_structure/defaults/tooling/unsafe_concurrency: N/A (example binaries).

### 2. Verdict
Approve-with-changes.

### 3. Execution status
- `cargo fmt --check` → clean · `cargo clippy --all-targets --all-features -- -D warnings` → clean · `cargo test --example ng_normalizer_screen` → 5/5 · `cargo test --example ng_synthesize_stress_reads` → 3/3 · `cargo test --lib` → 2331 passed.
- Ran on the host: GIAB HG002 10x → 63,757 indel reads, **0 disagreements**, 6 moved. Synthetic stress set → 64 reads, **28 disagreements** (1a-vs-1b, 1b-vs-1c), **0** (1a-vs-1c), **36 cap-hits** — matching a hand analysis exactly.

### 4. Open questions / assumptions
1. **The measurement logic and the generator's read construction are correct** — a reviewer traced both by hand against the normalizer sources and reproduced the 28/36 numbers; the reference window, `reference_offset = 0`, the `catch_unwind` on 1c, and the tallies are all sound. **The GIAB-0 and synthetic-28 numbers can be trusted.**

### 5. Findings

**Blocker (per rubric wording — a proof-bearing code path with no test)**
- **B1: `build_loci` had no tests** — the generator underpins the whole "screen is discriminating" proof, and a silent off-by-one would mis-calibrate the stimulus while still producing plausible numbers. **Applied**: `every_read_is_well_formed_and_placed_rightmost` (read/ref consumption balance, all-`A` middle) + `the_shift_distances_straddle_the_cap_as_reported` (asserts 36/28 from the lengths) + `the_flanks_do_not_extend_the_homopolymer`.

**Major**
- **M1: the flank-non-`A` invariant lived only in a comment.** **Applied** (the flank test above).
- **M2: the screen had no test for the D==20 cap-hit-without-disagreement boundary.** **Applied**: `an_indel_exactly_at_the_cap_hits_it_but_still_agrees`.

**Minor**
- **Mi1: `Tally` derived `Copy`** on an ~88-byte mutable accumulator only ever taken `&mut` — a silent-copy footgun. **Applied**: dropped `Copy`.
- **Mi2: error prints used bare `{error}`**, dropping the `source()` chain. **Applied**: an `error_chain` helper prints the full chain.
- **Mi3: `matched` (a base count) is a bare participle.** **Applied**: `build_loci` refactored (`append_locus`/`read_sequence` helpers de-duplicate the arms), removing the binding.

**Nits / recorded**
- `Err(_)` in the screen's per-read recovery discards *why* a read/window failed (kept as counts; 0 observed on GIAB) → recorded, not expanded (example tool; the counts sufficed). `main` length → the tested core is `screen_one`, already factored + tested; left as-is (example-appropriate). Canonical-vs-raw window fetch → canonical (uppercase ACGTN) is correct against the uppercase `MappedRead.seq`, and all three normalizers see identical input regardless, so it does not affect the counts.

### 6. What's good
- The screen's core (`screen_one`) is a small pure function, unit-tested with fixtures that pin *both* the agree case and each disagreement mechanism (1b cap, complex D/I) — so the discriminating power is a committed test, not just a run artifact.
- The generator's numbers are derived from first principles and pinned by a test that recomputes 36/28 from the length list — a self-checking stimulus.

### 7. Commands to re-verify
- `./scripts/dev.sh cargo test --example ng_normalizer_screen --example ng_synthesize_stress_reads`
- `cargo run --release --example ng_synthesize_stress_reads -- <dir>` then `cargo run --release --example ng_normalizer_screen -- <dir>/stress_ref.fa <dir>/stress_reads.bam`
