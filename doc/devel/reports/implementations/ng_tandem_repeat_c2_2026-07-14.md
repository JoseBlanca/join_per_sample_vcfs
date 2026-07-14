# ng tandem-repeat scanner — Milestone C2 (windowed streaming)

*Implementation report, 2026-07-14. Plan:
[`doc/devel/ng/impl_plan/ssr_repeat_scanner.md`](../../ng/impl_plan/ssr_repeat_scanner.md)
(Milestone C, step C2). Design: [spec](../../ng/spec/ssr_repeat_scanner.md) §3.6 +
[arch](../../ng/arch/ssr_repeat_scanner.md) §2.2, §3.*

## What landed

`RegionScanner::stream(fetcher, periods, params, opts) -> Result<Self, ScanError>` — the
memory-bounded windowed path of the region seam, over the sealed `ChromRefFetcher`
(`StreamingChromRefFetcher`, ~1 MB resident). It scans a whole contig in `window_bp` cores
without holding it in memory, and yields the same `Region` iterator as `over_slice`.

The region seam was **restructured into a two-level design** (`build_regions`) so both paths
share one region builder:

- **coverage** (span-level, repeat-covered positions) classifies `Repeat`/`Satellite`/`Unique`
  and fixes the boundaries; and
- **exact STR intervals** (start-attributed) decorate each `Repeat`.

`collect_windowed` fetches each core ± a `max_repeat_len` margin, then per found interval adds
its **core-clipped** span to coverage (so adjacent cores tile and a satellite rejoins across the
windows it spans) and, if its **start is in the core**, the **whole** interval to the STR list
(an STR ≤ `max_repeat_len` is captured whole by the margins, so segmented identically to a
whole-contig scan and attributed once).

## The bugs the invariance test caught (and the fixes)

The `window_bp`-invariance test (streamed tiling must equal the resident oracle for every
window size) earned its keep — it caught **three** issues in a row that inspection missed:

1. **Online Ruzzo–Tompa flush was under-tested, not wrong here.** The first failure looked like
   the RT flush; removing it didn't fix the failure, and a strengthened property test (len < 160,
   2000 cases) confirmed the RT is correct. The flush was dropped anyway — the pass is now plain
   offline RT (O(n) stack; an online-finalisation memory optimisation is deferred behind the
   stronger property test).
2. **Left margin too small.** With a `periods.max` (6 bp) left margin, a window starting mid-repeat
   *re-segmented* a long repeat (Ruzzo–Tompa segmentation is context-dependent) and emitted a
   spurious in-core fragment. Fix: **both margins are `max_repeat_len`**, so any STR-sized repeat
   touching a core is seen whole and segmented identically.
3. **Satellite coalescing needed the two-level split.** Start-attribution alone drops interior
   fragments of a satellite (longer than any window's reach). Fix: coverage is **clipped to the
   core** (every core contributes its slice), so the union rejoins the satellite; STR intervals
   stay start-attributed and exact.

## The invariance guarantee (and its honest limit)

`stream` is **window-count invariant for the STR (`Repeat`/`Unique`) tiling** — verified by
`stream_is_window_invariant_across_truncated_windows` (a 30 bp cap < the contig, so small windows
genuinely truncate the fetch; two isolated STRs each captured whole across windows). A
**satellite** longer than a window's reach *may* segment differently between windowed and resident
scans, because its internal interruption-bridging depends on global context a bounded window
cannot see. This is **acceptable and documented**: satellites are mask/skip regions, not
genotyping targets. `stream_detects_a_satellite_within_a_window` asserts the fits-in-window case.

## Deviations (recorded)

- **`ChromRefFetcher` is sealed** → no in-memory impl; `stream` takes `impl ChromRefFetcher` and
  tests build one from a project-local temp FASTA + hand-computed `.fai`. `RegionScanner` is a
  plain (non-generic) struct; `over_slice` operates on the slice directly (C1).
- **`stream` is eager-windowed**, not lazily interleaved: it scans all windows in the constructor
  (returning `Result` — a fetch error is fail-fast) and iterates the resulting regions. Memory is
  still bounded (one window fetch + the small interval/region lists); truly-lazy per-window yield
  is a possible later refinement.
- **`window_bp` should be ≫ `max_repeat_len`** for satellites to fit a window (production default
  100 kb ≫ 1 kb). Mega-satellites exceeding `window_bp` fall under the satellite limit above.

## Validation

`cargo fmt --check` clean; `cargo clippy --lib --tests -- -D warnings` clean on the module
(one `type_complexity` lint resolved with a `WindowedScan` type alias);
`cargo test --lib ng::tandem_repeat` → **25 pass** (incl. the strengthened RT property test and
the four streaming tests). Temp FASTAs land in the gitignored project-local `tmp/`.

## Next

**Milestone C is complete.** Next is Milestone D — validate the scanner against the `trf-mod`
golden catalog (parity), with production left on `trf-mod`.
