# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for pop_var_caller's cohort .psp -> VCF
pipeline at synthetic cohort sizes far beyond the real tomato1 set
(50, 200, 1000 vs the existing 1..26 sweep in perf_ours_joint.py).

Each pass is **idempotent**: an existing samply profile.json.gz or
dhat heap.json at the conventional output path is reused (re-parsed
for the per-bucket tables) instead of being rebuilt. To force a fresh
recording, delete the artefact before re-running. Bare runs always
re-execute because peak-RSS is the only thing we can't reconstruct
from an on-disk file. The TSVs are rewritten after every completed
per-N step so a kill mid-sweep keeps everything finished so far.

Replicates a single real .psp into N synthetic samples via the
examples/profile_cohort_e2e binary, then measures the joint step end
to end. Captures up to three passes per N:

  1. Bare wall + peak RSS via psutil polling
     -> benchmarks/tomato1/results/perf/scaling_synthetic.tsv
  2. samply CPU profile (host-side; macOS or Linux)
     -> tmp/scaling_synthetic/samply/N<n>.profile.json.gz
        + inclusive per-module CPU share -> scaling_synthetic_cpu.tsv
  3. dhat heap profile (inside the dev container via scripts/dev.sh)
     -> tmp/scaling_synthetic/dhat/N<n>.heap.json
        + per-module total/peak bytes -> scaling_synthetic_heap.tsv

Synthetic replication is structurally valid for both memory and CPU
distribution: the records carry no per-sample identity, so duplicating
them stresses the joint pipeline's per-sample state identically to a
real cohort. Per-record CPU cost is data-dependent but the *shape* of
the cost surface (which stages dominate, how they scale with N) is
not.

Source PSP defaults to the first sorted tomato1 PSP under the new
512 KiB block target (post-fcef495).

Env overrides:
  PROFILE_BIN     host-native examples/profile_cohort_e2e
                  (default: $PROJECT_ROOT/target/release/examples/profile_cohort_e2e)
  DHAT_BIN        container-built examples/dhat_var_calling
                  (default: $PROJECT_ROOT/target-container/release/examples/dhat_var_calling)
  SOURCE_PSP      single .psp to replicate (default: first sorted PSP
                  in benchmarks/tomato1/results/ours/cohort/psp/)
  REFERENCE       SL4.0 fasta (.fai sibling required)
  THREADS         rayon worker count (default: 4)
  SIZES           comma-separated cohort sizes (default: 50,200,1000)
  MAX_DHAT_N      skip dhat at any N strictly greater than this
                  (default: 200 — dhat's per-allocation backtrace
                  overhead at N=1000 needs >150 GB container heap;
                  the default container cap is 16 GiB)
  DEV_MEM_DHAT    Apple-container memory cap for dhat passes
                  (default: 48g; pass-through to scripts/dev.sh's
                  DEV_MEM env)
  SKIP_SAMPLY     non-empty => skip the samply pass
  SKIP_DHAT       non-empty => skip the dhat pass
  DEV_EXTRA_MOUNT path to bind-mount into the container (default:
                  $HOME/genomes — needed so the reference fasta is
                  visible inside the dhat run's container)

Invoke (host-side; macOS w/ Apple container or Linux dev box):
  uv run --script benchmarks/tomato1/scripts/perf_scaling_synthetic.py
"""

import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR, PROJECT_ROOT, PSP_DIR,
    Measurement, banner, check_exists, measure, write_tsv,
)

CALLER = "scaling_synthetic"
DEFAULT_SIZES = [50, 200, 1000]

SCRATCH = PROJECT_ROOT / "tmp" / "scaling_synthetic"
COHORT_DIR = SCRATCH / "cohort"
VCF_DIR = SCRATCH / "vcf"
SAMPLY_DIR = SCRATCH / "samply"
DHAT_DIR = SCRATCH / "dhat"
DHAT_RUN_BASE = SCRATCH / "dhat_runs"

PROFILE_BIN = Path(os.environ.get(
    "PROFILE_BIN",
    str(PROJECT_ROOT / "target" / "release" / "examples" / "profile_cohort_e2e"),
))
DHAT_BIN = Path(os.environ.get(
    "DHAT_BIN",
    str(PROJECT_ROOT / "target-container" / "release" / "examples" / "dhat_var_calling"),
))


def _default_source_psp() -> Path:
    override = os.environ.get("SOURCE_PSP")
    if override:
        return Path(override)
    psps = sorted(PSP_DIR.glob("*.psp"))
    if not psps:
        sys.exit(
            f"no source .psp in {PSP_DIR}; set SOURCE_PSP=<path> or "
            "rebuild tomato1 PSPs (run_ours_cohort.sh stage 1)"
        )
    return psps[0]


SOURCE_PSP = _default_source_psp()
DEV_SH = PROJECT_ROOT / "scripts" / "dev.sh"


# Order matters: more specific prefixes first so a leaf frame inside
# e.g. posterior_engine isn't swallowed by the broader var_calling /
# pop_var_caller buckets that come later in this list.
MODULE_PREFIXES: list[tuple[str, str]] = [
    ("posterior_engine",          "pop_var_caller::var_calling::posterior_engine"),
    ("per_group_merger",          "pop_var_caller::var_calling::per_group_merger"),
    ("dust_filter",               "pop_var_caller::var_calling::dust_filter"),
    ("per_position_merger",       "pop_var_caller::var_calling::per_position_merger"),
    ("variant_grouping",          "pop_var_caller::var_calling::variant_grouping"),
    ("vcf_writer",                "pop_var_caller::var_calling::vcf_writer"),
    ("contamination_estimation",  "pop_var_caller::var_calling::contamination_estimation"),
    ("var_calling_other",         "pop_var_caller::var_calling"),
    ("psp",                       "pop_var_caller::psp"),
    ("ref_fetcher",               "pop_var_caller::per_sample_pileup::ref_fetcher"),
    ("per_sample_pileup_other",   "pop_var_caller::per_sample_pileup"),
    ("bam",                       "pop_var_caller::bam"),
    ("cli_driver",                "pop_var_caller::pop_var_caller"),
    ("project_other",             "pop_var_caller::"),
]

ALLOCATOR_LIB_HINTS = ("libsystem_malloc", "libjemalloc", "libmimalloc",
                       "tcmalloc", "libc.so")
KERNEL_LIB_HINTS = ("libsystem_kernel", "libsystem_platform",
                    "linux-vdso", "vdso")
THREADING_LIB_HINTS = ("libsystem_pthread", "libpthread")


_V0_PATH_TOKEN = re.compile(r"(\d+)([A-Za-z_][A-Za-z0-9_]*)")


def _normalize_symbol(func_name: str) -> str:
    """Canonicalise a Mach-O / ELF symbol string into a form whose
    module path uses `::` separators, so a single MODULE_PREFIXES list
    matches all three mangling forms atos / nm emit:

      1. Modern demangled        — `pop_var_caller::var_calling::...`
         (atos already demangles `_ZN`-style symbols this way)
      2. Legacy trait-impl form  — `_$LT$pop_var_caller..var_calling..`
         (atos leaves `<` / `>` / `::` as `$LT$` / `$GT$` / `..` here)
      3. Rust v0 mangling        — `_RINv...11var_calling16posterior_engine`
         (atos doesn't demangle v0; we extract length-prefixed
         segments and join them with `::`)

    Forms 1 + 2 are handled by a single replace; form 3 is appended
    as an extra synthetic suffix so substring matching against
    'pop_var_caller::var_calling::posterior_engine' still hits."""
    canonical = func_name.replace("..", "::")
    if "_R" in func_name and "pop_var_caller" in func_name:
        tokens = [name for _, name in _V0_PATH_TOKEN.findall(func_name)]
        if tokens:
            canonical = canonical + " " + "::".join(tokens)
    return canonical


def bucket_for(func_name: str) -> str:
    """Map a fully-qualified function name to a project module bucket
    or to 'other'. Allocator / kernel / threading frames are
    classified at the library level (see parse_samply_profile) and
    never reach this fn — by the time we're here we've already
    confirmed the frame is in the project binary's text segment."""
    canonical = _normalize_symbol(func_name)
    for label, prefix in MODULE_PREFIXES:
        if prefix in canonical:
            return label
    return "other"


# ---- atos-based symbolication for samply Firefox-profiler JSON ---------
#
# samply on macOS records frames as `(library-name, hex-offset-string)`
# and does symbolication on demand when the user opens the profile in
# the UI. The raw JSON therefore has no function names — just
# 'profile_cohort_e2e' + '0xefaaf' style entries — and any text-based
# bucketing of those strings buckets 100% into 'other'.
#
# To bucket programmatically we resolve the offsets ourselves with
# macOS's `atos`. The binary's __TEXT segment loads at vmaddr
# 0x100000000 (confirmed via `otool -l`); add that to each
# library-relative offset, batch the resulting absolute addresses into
# one atos call per binary, parse the output into a {addr: name} map,
# then bucket-match against MODULE_PREFIXES.

TEXT_VMADDR_DEFAULT = 0x100000000  # __TEXT segment base for our Mach-O


def _text_vmaddr(binary: Path) -> int:
    """Resolve the binary's __TEXT segment vmaddr via `otool -l`.
    Falls back to TEXT_VMADDR_DEFAULT if otool isn't usable (e.g.
    cross-platform run; samply parsing won't work in that case
    anyway so the fallback is mostly cosmetic)."""
    if shutil.which("otool") is None:
        return TEXT_VMADDR_DEFAULT
    try:
        out = subprocess.check_output(["otool", "-l", str(binary)], text=True)
    except subprocess.CalledProcessError:
        return TEXT_VMADDR_DEFAULT
    in_text = False
    for line in out.splitlines():
        s = line.strip()
        if s == "segname __TEXT":
            in_text = True
        elif in_text and s.startswith("vmaddr "):
            return int(s.split()[1], 0)
    return TEXT_VMADDR_DEFAULT


def _atos_resolve(binary: Path, vmaddr: int, offsets: list[int]) -> dict[int, str]:
    """Resolve a batch of library-relative offsets to function-name
    strings via `atos`. Returns {offset: symbol_string}. Offsets that
    fail to resolve (atos echoes the address back) map to the original
    hex string so they fall into the 'other' bucket gracefully."""
    if not offsets:
        return {}
    addrs = [f"0x{vmaddr + off:x}" for off in offsets]
    # atos accepts a long argv. Mach-O Rust binaries on this project
    # have ~10 K unique frame offsets per profile, well under the
    # ARG_MAX limit.
    cmd = ["atos", "-o", str(binary), *addrs]
    try:
        out = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        out = exc.output or ""
    lines = out.splitlines()
    result: dict[int, str] = {}
    for off, addr_str, line in zip(offsets, addrs, lines):
        result[off] = line if line and line != addr_str else f"0x{off:x}"
    return result


# ---- samply Firefox-profiler JSON parser --------------------------------

def parse_samply_profile(path: Path, binary: Path) -> dict[str, float]:
    """Return inclusive % CPU share per module bucket from a samply
    Firefox-profiler .json.gz dump. Inclusive = a sample's stack
    crediting every bucket it visits. Percentages can sum to >100%
    because a single sample lives in multiple buckets simultaneously
    (a frame in posterior_engine that called into allocator credits
    both).

    Library-level buckets (allocator / kernel / threading) come from
    the per-frame resource entry. Project-binary frames are resolved
    to a Rust mangled name via atos, then matched against
    MODULE_PREFIXES."""
    with gzip.open(path, "rt") as fh:
        prof = json.load(fh)

    binary_name = binary.name  # e.g. 'profile_cohort_e2e'
    vmaddr = _text_vmaddr(binary)

    total_weight = 0.0
    inclusive: dict[str, float] = {}

    hex_pat = re.compile(r"^0x[0-9a-fA-F]+$")

    for thread in prof.get("threads", []):
        strs = thread.get("stringArray") or thread.get("stringTable", [])
        func_table = thread.get("funcTable", {})
        frame_table = thread.get("frameTable", {})
        stack_table = thread.get("stackTable", {})
        resource_table = thread.get("resourceTable", {})
        samples = thread.get("samples", {})

        # 1. Resolve every (resource_idx, func_name_str) frame into a
        #    bucket label up-front, then walk the call-stack DAG.
        res_names = resource_table.get("name", [])  # idx -> stringArray idx
        res_lib_str = {
            i: strs[res_names[i]] if i < len(res_names) else "<no-resource>"
            for i in range(len(res_names))
        }

        # Bucket library frames: allocator/kernel/threading/binary/other.
        def lib_bucket(lib_str: str) -> str | None:
            if lib_str == binary_name:
                return None  # caller resolves via atos
            if any(h in lib_str for h in ALLOCATOR_LIB_HINTS):
                return "allocator"
            if any(h in lib_str for h in KERNEL_LIB_HINTS):
                return "kernel"
            if any(h in lib_str for h in THREADING_LIB_HINTS):
                return "threading"
            return "other"

        # 2. Collect every funcTable entry that needs atos resolution
        #    (frames inside the project binary with a hex placeholder).
        func_names = func_table["name"]
        func_resources = func_table.get("resource", [None] * len(func_names))
        offsets_to_resolve: set[int] = set()
        per_func_bucket: list[str] = [None] * len(func_names)  # type: ignore[assignment]
        for fi in range(len(func_names)):
            nm = strs[func_names[fi]]
            res_idx = func_resources[fi] if fi < len(func_resources) else None
            if res_idx is not None and res_idx >= 0:
                lib = res_lib_str.get(res_idx, "<no-resource>")
                lb = lib_bucket(lib)
                if lb is not None:
                    per_func_bucket[fi] = lb
                    continue
            if hex_pat.match(nm):
                offsets_to_resolve.add(int(nm, 16))

        # 3. Batch-resolve via atos.
        resolved = _atos_resolve(binary, vmaddr, sorted(offsets_to_resolve))

        # 4. Second pass over funcTable: project-binary frames bucket
        #    via the resolved name; remaining unbucketed frames default
        #    to 'other'.
        for fi in range(len(func_names)):
            if per_func_bucket[fi] is not None:
                continue
            nm_idx = func_names[fi]
            nm = strs[nm_idx]
            if hex_pat.match(nm):
                off = int(nm, 16)
                resolved_name = resolved.get(off, nm)
                per_func_bucket[fi] = bucket_for(resolved_name)
            else:
                per_func_bucket[fi] = bucket_for(nm)

        # 5. Walk the stackTable DAG once, caching the bucket-set per
        #    stack index (the table is a prefix-shared linked list, so
        #    every visit memoises).
        prefix_col = stack_table.get("prefix", [])
        frame_col = stack_table.get("frame", [])
        frame_func_col = frame_table.get("func", [])

        stack_buckets: dict[int, frozenset] = {}

        def buckets_for_stack(s):
            if s is None:
                return frozenset()
            cached = stack_buckets.get(s)
            if cached is not None:
                return cached
            visited: set[str] = set()
            cur = s
            while cur is not None and cur >= 0:
                frame_idx = frame_col[cur]
                func_idx = frame_func_col[frame_idx]
                visited.add(per_func_bucket[func_idx])
                cur = prefix_col[cur]
            result = frozenset(visited)
            stack_buckets[s] = result
            return result

        stacks = samples.get("stack", [])
        weights = samples.get("weight") or [1] * len(stacks)
        for i, s in enumerate(stacks):
            w = weights[i] if i < len(weights) else 1
            if w is None:
                w = 1
            total_weight += w
            for bucket in buckets_for_stack(s):
                inclusive[bucket] = inclusive.get(bucket, 0.0) + w

    if total_weight == 0:
        return {}
    return {b: 100.0 * cnt / total_weight for b, cnt in inclusive.items()}


# ---- dhat heap-profile parser -------------------------------------------

def parse_dhat_heap(path: Path) -> dict[str, dict[str, int]]:
    """Return {bucket: {total_bytes, max_bytes, total_blocks}} from a
    dhat-heap.json. Each program point is credited to the leafmost
    in-project frame on its stack (i.e. closest project caller of the
    allocation); 'other' / 'allocator' buckets are skipped while
    searching, so allocator-internal frames don't shadow the user's
    code."""
    with path.open() as fh:
        prof = json.load(fh)

    ftbl = prof.get("ftbl", [])
    pps = prof.get("pps", [])

    agg: dict[str, dict[str, int]] = {}
    for pp in pps:
        fs = pp.get("fs", [])
        chosen = None
        for idx in fs:
            if idx < 0 or idx >= len(ftbl):
                continue
            frame_str = ftbl[idx]
            b = bucket_for(frame_str)
            if b not in ("other", "allocator"):
                chosen = b
                break
        if chosen is None:
            chosen = "other"
        rec = agg.setdefault(chosen, {"total_bytes": 0, "max_bytes": 0, "total_blocks": 0})
        rec["total_bytes"] += pp.get("tb", 0)
        rec["max_bytes"] += pp.get("mb", 0)
        rec["total_blocks"] += pp.get("tbk", 0)
    return agg


# ---- Per-N driver --------------------------------------------------------

def cmd_profile(n: int, vcf_path: Path, threads: int) -> list[str]:
    return [
        str(PROFILE_BIN),
        "--psp", str(SOURCE_PSP),
        "--n-samples", str(n),
        "--reference", str(DEFAULT_REFERENCE),
        "--output", str(vcf_path),
        "--threads", str(threads),
        "--cohort-dir", str(COHORT_DIR),
    ]


def run_bare(n: int, threads: int) -> Measurement:
    vcf = VCF_DIR / f"N{n:04d}.vcf"
    wall, peak_bytes, exit_code = measure(cmd_profile(n, vcf, threads))
    return Measurement(CALLER, n, wall, peak_bytes / 1024 / 1024, exit_code)


def run_samply(n: int, threads: int) -> tuple[float, int, Path, bool]:
    """Returns (wall_seconds, exit_code, output_path, reused_existing).
    If the output file already exists, returns immediately with
    wall_seconds=0 and reused_existing=True; the per-bucket table can
    still be re-derived from the existing artefact."""
    out = SAMPLY_DIR / f"N{n:04d}.profile.json.gz"
    if out.exists():
        return 0.0, 0, out, True
    vcf = VCF_DIR / f"N{n:04d}.samply.vcf"
    cmd = [
        "samply", "record", "--save-only", "--no-open",
        "-o", str(out),
        "--",
        *cmd_profile(n, vcf, threads),
    ]
    wall, _peak, exit_code = measure(cmd)
    return wall, exit_code, out, False


def run_dhat(n: int, threads: int) -> tuple[float, int, Path, bool]:
    """Run dhat_var_calling inside the dev container. dhat writes
    `dhat-heap.json` in the process's cwd, so we run it from a per-N
    scratch dir and move the file under DHAT_DIR afterwards.

    Returns (wall_seconds, exit_code, output_path, reused_existing).
    If the output file already exists, returns immediately."""
    out_target = DHAT_DIR / f"N{n:04d}.heap.json"
    if out_target.exists():
        return 0.0, 0, out_target, True
    vcf = VCF_DIR / f"N{n:04d}.dhat.vcf"
    run_dir = DHAT_RUN_BASE / f"N{n:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    # The reference fasta typically lives at $HOME/genomes/... which is
    # outside the project tree; dev.sh bind-mounts a single extra path
    # only if DEV_EXTRA_MOUNT is set. Default to $HOME/genomes here so
    # the dhat run can read the FASTA inside the container.
    if "DEV_EXTRA_MOUNT" not in env:
        if not str(DEFAULT_REFERENCE).startswith(str(PROJECT_ROOT)):
            env["DEV_EXTRA_MOUNT"] = str(Path.home() / "genomes")
    # Bump container memory cap for dhat — every allocation gets a
    # backtrace + per-pp record, so memory overhead is many × the
    # production heap. dev.sh's DEV_MEM default of 16g is too tight
    # for anything past N≈200 on this workload.
    env.setdefault("DEV_MEM", os.environ.get("DEV_MEM_DHAT", "48g"))

    inner_cmd = (
        f"cd {run_dir} && "
        f"{DHAT_BIN} "
        f"--psp-dir {COHORT_DIR} "
        f"--n-samples {n} "
        f"--reference {DEFAULT_REFERENCE} "
        f"--output {vcf} "
        f"--threads {threads}"
    )
    cmd = [str(DEV_SH), "sh", "-c", inner_cmd]
    wall, _peak, exit_code = measure(cmd, env=env)

    src_heap = run_dir / "dhat-heap.json"
    if src_heap.exists():
        shutil.move(str(src_heap), str(out_target))
    return wall, exit_code, out_target, False


def write_cpu_tsv(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    bucket_cols = sorted({
        k for r in rows for k in r if k not in ("n_samples", "wall_seconds")
    })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        header = ["caller", "n_samples", "wall_seconds"] + [f"pct_{b}" for b in bucket_cols]
        fh.write("\t".join(header) + "\n")
        for r in rows:
            vals = [CALLER, str(r["n_samples"]), f"{r['wall_seconds']:.3f}"]
            for b in bucket_cols:
                vals.append(f"{r.get(b, 0.0):.2f}")
            fh.write("\t".join(vals) + "\n")


def write_heap_tsv(rows: list[dict], path: Path) -> None:
    if not rows:
        return
    bucket_cols = sorted({
        k for r in rows for k in r if k not in ("n_samples", "wall_seconds")
    })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        header = ["caller", "n_samples", "wall_seconds"] + bucket_cols
        fh.write("\t".join(header) + "\n")
        for r in rows:
            vals = [CALLER, str(r["n_samples"]), f"{r['wall_seconds']:.3f}"]
            for b in bucket_cols:
                vals.append(str(r.get(b, 0)))
            fh.write("\t".join(vals) + "\n")


def flush_tsvs(rows: list[Measurement],
               cpu_rows: list[dict],
               heap_rows: list[dict]) -> None:
    """Rewrite all three TSVs from the in-memory accumulators, sorted
    ascending by N. Called after every per-N step so a kill mid-sweep
    preserves whatever ran. Cheap (three small writes)."""
    PERF_DIR.mkdir(parents=True, exist_ok=True)
    sorted_rows = sorted(rows, key=lambda m: m.n_samples)
    write_tsv(sorted_rows, PERF_DIR / f"{CALLER}.tsv")
    if cpu_rows:
        write_cpu_tsv(sorted(cpu_rows, key=lambda r: r["n_samples"]),
                      PERF_DIR / f"{CALLER}_cpu.tsv")
    if heap_rows:
        write_heap_tsv(sorted(heap_rows, key=lambda r: r["n_samples"]),
                       PERF_DIR / f"{CALLER}_heap.tsv")


def main() -> int:
    sizes_raw = os.environ.get("SIZES")
    sizes_set = (
        {int(s) for s in sizes_raw.split(",") if s.strip()}
        if sizes_raw else
        set(DEFAULT_SIZES)
    )
    # Iterate descending so the largest-N pass populates the cohort
    # replica dir first; smaller-N passes then see all replicas already
    # on disk and their bare wall time excludes replica-write I/O.
    sizes = sorted(sizes_set, reverse=True)
    skip_samply = bool(os.environ.get("SKIP_SAMPLY"))
    skip_dhat = bool(os.environ.get("SKIP_DHAT"))
    threads = int(os.environ.get("THREADS", str(DEFAULT_THREADS)))
    max_dhat_n = int(os.environ.get("MAX_DHAT_N", "200"))

    check_exists(PROFILE_BIN, SOURCE_PSP, DEFAULT_REFERENCE,
                 Path(str(DEFAULT_REFERENCE) + ".fai"))
    if not skip_dhat:
        check_exists(DHAT_BIN, DEV_SH)
    if not skip_samply:
        if shutil.which("samply") is None:
            print("warning: 'samply' not found on PATH; skipping the samply pass",
                  file=sys.stderr)
            skip_samply = True

    for d in (SCRATCH, COHORT_DIR, VCF_DIR, SAMPLY_DIR, DHAT_DIR, DHAT_RUN_BASE):
        d.mkdir(parents=True, exist_ok=True)

    banner(CALLER, sizes)
    print(f"source psp: {SOURCE_PSP}")
    print(f"cohort dir: {COHORT_DIR}")
    print(f"profile bin (host): {PROFILE_BIN}")
    if not skip_dhat:
        print(f"dhat bin (container): {DHAT_BIN}")
    print(f"samply pass: {'enabled' if not skip_samply else 'skipped'}")
    print(f"dhat pass:   {'enabled' if not skip_dhat else 'skipped'}")
    print()

    rows: list[Measurement] = []
    cpu_rows: list[dict] = []
    heap_rows: list[dict] = []

    for n in sizes:
        print(f"=== N={n} ===")

        # Pass 1: bare wall + peak RSS — populates COHORT_DIR for the
        # following passes as a side effect of profile_cohort_e2e.
        print(f"  [bare] T={threads} run + lazy replica synthesis")
        m = run_bare(n, threads)
        rows.append(m)
        print(f"  -> wall={m.wall_seconds:.1f}s  peak={m.peak_rss_mb:.0f} MB  exit={m.exit_code}")
        if m.exit_code != 0:
            print(f"  !! bare run failed (exit {m.exit_code}); skipping samply/dhat for N={n}")
            continue

        # Pass 2: samply CPU profile (idempotent — reuses an existing
        # profile.json.gz if one is on disk).
        if not skip_samply:
            print(f"  [samply] recording CPU profile")
            s_wall, s_exit, s_path, s_reused = run_samply(n, threads)
            if s_exit == 0 and s_path.exists():
                tag = "reused" if s_reused else f"wall={s_wall:.1f}s"
                try:
                    shares = parse_samply_profile(s_path, PROFILE_BIN)
                except Exception as exc:  # noqa: BLE001
                    print(f"  !! samply parse failed: {exc}")
                    shares = {}
                cpu_rows.append({"n_samples": n, "wall_seconds": s_wall, **shares})
                top = sorted(shares.items(), key=lambda kv: -kv[1])[:6]
                print(f"  -> samply {tag}; top inclusive: "
                      + ", ".join(f"{k}={v:.1f}%" for k, v in top))
            else:
                print(f"  !! samply failed (exit {s_exit})")

        # Pass 3: dhat heap profile (idempotent; skipped above MAX_DHAT_N
        # because dhat's per-allocation bookkeeping needs many × the
        # production heap and OOMs the container at large N).
        if not skip_dhat:
            if n > max_dhat_n:
                print(f"  [dhat] skipped: N={n} exceeds MAX_DHAT_N={max_dhat_n}")
            else:
                print(f"  [dhat] recording heap profile inside container")
                d_wall, d_exit, d_path, d_reused = run_dhat(n, threads)
                if d_exit == 0 and d_path.exists():
                    tag = "reused" if d_reused else f"wall={d_wall:.1f}s"
                    try:
                        heap = parse_dhat_heap(d_path)
                    except Exception as exc:  # noqa: BLE001
                        print(f"  !! dhat parse failed: {exc}")
                        heap = {}
                    row: dict = {"n_samples": n, "wall_seconds": d_wall}
                    for b, sums in heap.items():
                        row[f"tb_{b}"] = sums["total_bytes"]
                        row[f"mb_{b}"] = sums["max_bytes"]
                        row[f"blocks_{b}"] = sums["total_blocks"]
                    heap_rows.append(row)
                    top = sorted(heap.items(), key=lambda kv: -kv[1]["total_bytes"])[:5]
                    print(f"  -> dhat {tag}; top tb (MB): "
                          + ", ".join(f"{k}={s['total_bytes']/1e6:.0f}" for k, s in top))
                else:
                    print(f"  !! dhat failed (exit {d_exit})")

        # Incremental TSV write: rewrites every TSV after each per-N
        # step finishes, so a kill mid-sweep keeps whatever ran.
        flush_tsvs(rows, cpu_rows, heap_rows)
        sys.stdout.flush()

    # Final flush: idempotent with the per-N incremental writes above,
    # but cheap insurance in case the last step bailed before it ran.
    flush_tsvs(rows, cpu_rows, heap_rows)
    print(f"\nwrote {PERF_DIR / f'{CALLER}.tsv'}"
          + (f", {PERF_DIR / f'{CALLER}_cpu.tsv'}" if cpu_rows else "")
          + (f", {PERF_DIR / f'{CALLER}_heap.tsv'}" if heap_rows else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
