# ng reference_info — Milestones D & E (cache + entry points); module complete

*Implementation report, 2026-07-19. Plan-driven, closing the reference-info reader: D1 (the
single-flight cache), E1 and E2 (the two composing entry points). Design authority: spec
[`../../ng/spec/reference_info.md`](../../ng/spec/reference_info.md) (§3.7, §3.10, §3.11) and
arch [`../../ng/arch/reference_info.md`](../../ng/arch/reference_info.md). Code in
[`src/ng/reference_info.rs`](../../../../src/ng/reference_info.rs). This report closes the
module; earlier milestones: [A](ng_reference_info_milestone_a_2026-07-19.md),
[B](ng_reference_info_milestone_b_2026-07-19.md), C (commit `add2d2f`).*

## Changes made

- **D1** (`60d5482`) — `ReferenceInfoCache`: a caller-held, thread-safe, compute-once
  single-flight cache over `read_reference_info` (spec §3.7). Two-level lock (map lock only to
  get-or-create a key's `Slot`; slot lock held across the read = the single-flight). Key =
  `(source discriminant, per-file (path, size, mtime))`, so `Fai`/`Fasta{None}`/`Fasta{Some}`
  never cross-satisfy. Successes only cached; an un-stattable file bypasses (never fails).
  `get_or_read` delegates to a private `get_or_read_with` parameterised over the read — the
  test seam for counting/barrier-blocked reads.
- **E1** (`f6b00e3`) — `read_fai_verify_in_background` + `VerificationHandle`: the `.fai` now
  (foreground, `md5: None`), the FASTA verified on a background thread, **both through the
  cache** so the verify populates it. `#[must_use]` handle; `is_finished`; `join` returns the
  verified info by return; `Drop` without `join` warns.
- **E2** (`7d9f123`) — `read_reference_verifying_or_creating_fai`: derive `<fasta>.fai`, then
  present → background verify (`Some(handle)`, no write); absent → scan + `write_fai` +
  `(info, None)`. The one sanctioned discovery-and-write point; a write failure is fatal
  (`FaiWrite`), with the `get_or_read(Fasta{None})` escape hatch.

## Recorded deviations

- **E1 uses `std::thread::JoinHandle`, not a crossbeam channel** (arch sketched a channel). A
  `JoinHandle` *is* the result channel — the verified `Result` arrives by `join()` — so this
  meets the contract (return-by-join, `is_finished`, `#[must_use]`, a dedicated thread rather
  than a rayon worker) with no new dependency. The message-passing-not-mutation and
  not-rayon decisions stand.
- **The E2 fatal-write test is host-validated, skipped under root.** The dev container runs as
  root, which bypasses directory permission bits, so a read-only-dir write still succeeds
  there; and the write branch requires the sibling `.fai` to be *absent*, which rules out the
  container-safe "target is a directory" trick used for `write_fai`'s own atomicity test.
  The test probes whether read-only dirs are enforced and skips with a logged note when they
  are not; it runs for real on the host (non-root), where it passes.
- **`clippy::type_complexity`** on the nested cache map → a `Slot` type alias.

## Tests added (25 across D+E; 52 in the module)

D: hit-computes-once (`Arc::ptr_eq`); single-flight under 8 threads (one read, shared Arc);
different keys do not serialize (deterministic, via a channel); source-in-key (no poisoning);
error-not-cached; no-key bypass; `Send + Sync`. E1: info-before-join + join-upgrade; stale
`.fai` errors on join (mutation-verified); populates the cache; a concurrent reader
coordinates through the slot; failure does not poison; missing `.fai` errors before spawn;
abandoning warns. E2: both branches (present → verify/no-write; absent → scan + samtools-
identical `.fai`); the escape hatch; the fatal write (host).

## Validation

Every step in the container: fmt clean; clippy `--lib --tests --all-features -D warnings`
clean; full lib suite **1943 passed**. **On plan completion, verified natively on the host**
(cargo 1.95.0): fmt + clippy(lib,tests) clean, the full lib suite 1943 passed — and the
host run (non-root) exercises the E2 fatal-write path for real.

**Pre-existing, unrelated** (unchanged all session): `cargo clippy --all-targets -D warnings`
fails only on `examples/ssr_psp_seqdump.rs:41` (`unnecessary_sort_by`, a toolchain-1.95 patch
bump; on the clean tree; outside `src/ng`). Left untouched per the production freeze; awaiting
the owner's call on whether to fix that one line separately.

## Module status

`src/ng/reference_info.rs` is **complete** against the plan: the types, both reader arms, the
`.fai` writer, the single-flight cache, and the two entry points, with 52 tests. Out of scope
(next plans, per the spec): `@SQ`↔`.fai` reconciliation (read ingestion), consolidating ng's
fai-driven readers onto `ContigInfo` (an ng refactor that retires the `contig_list()` bridge),
and bgzip support / the two-phase parallel hash (deferred, error today).
