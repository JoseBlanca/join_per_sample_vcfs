# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for GATK *direct* joint calling (CRAM(s) ->
multi-sample VCF in ONE HaplotypeCaller process, no per-sample GVCF
intermediate).

This is the "direct cram -> vcf" GATK mode: a single HaplotypeCaller
invocation given every sample's CRAM via repeated `--input`, emitting
one joint multi-sample VCF directly. The pair-HMM is threaded with
`--native-pair-hmm-threads THREADS` (GATK's only built-in intra-process
parallelism), matching the THREADS budget freebayes and our caller get.

It serves two dashboard sections:
  - Section 1 (single-sample direct comparison) reads the N=1 row.
  - Section 3 (scaling cram -> vcf) reads the whole sweep.

Contrast with `perf_gatk_joint.py`, which measures the *recommended*
GVCF flow (per-sample HaplotypeCaller -ERC GVCF, then CombineGVCFs +
GenotypeGVCFs). Direct multi-sample HaplotypeCaller is the
non-incremental "call everything from scratch" path; it is expected to
be slower and heavier as N grows, which is exactly the contrast the
scaling panel shows.

Note on availability:
  GATK isn't in Homebrew. On macOS run this script inside the
  container (scripts/dev.sh) with the reference genome bind-mounted;
  on Linux either path works.

Sample sizes: 1, 2, 4, 8, 12, 16, 20, 24, 26 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/crams/*.bench.cram
              benchmarks/tomato1/regions.bed
Output:       benchmarks/tomato1/results/perf/gatk_direct.tsv
              benchmarks/tomato1/results/perf/gatk_direct/N<nn>.vcf

Env overrides:
  GATK_BIN    binary (default: /opt/gatk/gatk)
  REFERENCE   SL4.0 fasta (.fai + .dict siblings required)
  JAVA_HEAP   -Xmx for the HaplotypeCaller process (default: 8g)
  THREADS     --native-pair-hmm-threads (default: 4)
  SIZES       comma-separated sample sizes (default: 1,2,4,8,12,16,20,24,26)

Invoke (container path on macOS):
  DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh \\
      python3 benchmarks/tomato1/scripts/perf_gatk_direct.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    CRAM_DIR, DEFAULT_BED, DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR,
    Measurement, banner, check_exists, list_inputs, measure, pick_subset,
    sizes_from_env, write_tsv,
)

CALLER = "gatk_direct"
BIN = Path(os.environ.get("GATK_BIN", "/opt/gatk/gatk"))
JAVA_HEAP = os.environ.get("JAVA_HEAP", "8g")
OUT_DIR = PERF_DIR / CALLER


def main() -> int:
    sizes = sizes_from_env()
    ref_dict = DEFAULT_REFERENCE.with_suffix(".dict")
    check_exists(
        BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"),
        ref_dict, DEFAULT_BED,
    )
    crams = list_inputs(CRAM_DIR, ".bench.cram")
    banner(CALLER, sizes)
    print(f"pair-hmm threads: {DEFAULT_THREADS}")
    print()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[Measurement] = []
    for n in sizes:
        subset = pick_subset(crams, n)
        out_vcf = OUT_DIR / f"N{n:02d}.vcf"
        cmd = [
            str(BIN), "--java-options", f"-Xmx{JAVA_HEAP}", "HaplotypeCaller",
            "--reference", str(DEFAULT_REFERENCE),
            "--intervals", str(DEFAULT_BED),
            *sum((["--input", str(c)] for c in subset), []),
            "--output", str(out_vcf),
            "--native-pair-hmm-threads", str(DEFAULT_THREADS),
        ]
        print(
            f"[{CALLER}] N={n:>2}: one HaplotypeCaller over {len(subset)} "
            f"CRAMs (heap={JAVA_HEAP}, pair-hmm-threads={DEFAULT_THREADS})"
        )
        wall, peak_bytes, exit_code = measure(cmd)
        peak_mb = peak_bytes / 1024 / 1024
        rows.append(Measurement(CALLER, n, wall, peak_mb, exit_code))
        status = "ok" if exit_code == 0 else f"FAILED (exit {exit_code})"
        print(f"  -> {wall:.1f}s, peak {peak_mb:.0f} MB  [{status}]")

    tsv = PERF_DIR / f"{CALLER}.tsv"
    write_tsv(rows, tsv)
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
