# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Scaling perf experiment for GATK joint genotyping (CombineGVCFs +
GenotypeGVCFs, per-sample GVCFs -> joint VCF).

Comparable in shape to our PSP -> VCF mode (per-sample intermediate
-> joint VCF); the (perf_gatk_joint, perf_ours_joint) pair shows how
the two joint-genotyping stages scale relative to each other. The
HaplotypeCaller GVCF step is excluded — this measures the joint
stages only.

Wall time is reported as the sum of the two stages; peak RSS as the
max (they run sequentially, not concurrently). GVCF inputs must
already exist; they're built by stage 1 of
`../../lib/run_gatk.sh bench.config.sh cohort`.

Note on availability:
  GATK isn't in Homebrew. On macOS you either install it natively
  (manual download + Java runtime) or run this script inside the
  container — in which case scripts/dev.sh needs a second bind
  mount for the reference genome (which lives outside the project
  tree, default $HOME/genomes/...).

Sample sizes: 1, 2, 4, 8, 12, 16, 18 (full tomato1 cohort).
Inputs:       benchmarks/tomato1/results/gatk/cohort/gvcf/*.g.vcf.gz
              benchmarks/tomato1/regions.bed
Output:       benchmarks/tomato1/results/perf/gatk_joint.tsv
              benchmarks/tomato1/results/perf/gatk_joint/N<nn>.vcf
              benchmarks/tomato1/results/perf/gatk_joint/N<nn>.combined.g.vcf.gz

Env overrides:
  GATK_BIN    binary (default: /opt/gatk/gatk)
  REFERENCE   SL4.0 fasta (.fai + .dict siblings required)
  JAVA_HEAP   -Xmx for both stages (default: 4g)
  SIZES       comma-separated sample sizes (default: 1,2,4,8,12,16,18)

Invoke:
  uv run --script benchmarks/tomato1/scripts/perf_gatk_joint.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    DEFAULT_BED, DEFAULT_REFERENCE, GVCF_DIR, PERF_DIR,
    Measurement, banner, check_exists, list_inputs, measure, pick_subset,
    sizes_from_env, write_tsv,
)

CALLER = "gatk_joint"
BIN = Path(os.environ.get("GATK_BIN", "/opt/gatk/gatk"))
JAVA_HEAP = os.environ.get("JAVA_HEAP", "4g")
OUT_DIR = PERF_DIR / CALLER


def main() -> int:
    sizes = sizes_from_env()
    # GATK uniquely requires a .dict sibling for the FASTA.
    ref_dict = DEFAULT_REFERENCE.with_suffix(".dict")
    check_exists(
        BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"),
        ref_dict, DEFAULT_BED,
    )
    gvcfs = list_inputs(GVCF_DIR, ".g.vcf.gz")
    banner(CALLER, sizes)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[Measurement] = []
    for n in sizes:
        subset = pick_subset(gvcfs, n)
        combined = OUT_DIR / f"N{n:02d}.combined.g.vcf.gz"
        out_vcf = OUT_DIR / f"N{n:02d}.vcf"

        combine_cmd = [
            str(BIN), "--java-options", f"-Xmx{JAVA_HEAP}", "CombineGVCFs",
            "--reference", str(DEFAULT_REFERENCE),
            "--intervals", str(DEFAULT_BED),
            *sum((["--variant", str(g)] for g in subset), []),
            "--output", str(combined),
        ]
        genotype_cmd = [
            str(BIN), "--java-options", f"-Xmx{JAVA_HEAP}", "GenotypeGVCFs",
            "--reference", str(DEFAULT_REFERENCE),
            "--intervals", str(DEFAULT_BED),
            "--variant", str(combined),
            "--output", str(out_vcf),
        ]

        print(f"[{CALLER}] N={n:>2}: CombineGVCFs over {len(subset)} GVCFs")
        cb_wall, cb_peak, cb_ec = measure(combine_cmd)
        print(f"  combine  -> {cb_wall:.1f}s, peak {cb_peak / 1024 / 1024:.0f} MB")

        if cb_ec != 0:
            rows.append(Measurement(CALLER, n, cb_wall, cb_peak / 1024 / 1024, cb_ec))
            print(f"  CombineGVCFs FAILED (exit {cb_ec}); skipping GenotypeGVCFs for this N")
            continue

        print(f"[{CALLER}] N={n:>2}: GenotypeGVCFs")
        gt_wall, gt_peak, gt_ec = measure(genotype_cmd)
        print(f"  genotype -> {gt_wall:.1f}s, peak {gt_peak / 1024 / 1024:.0f} MB")

        total_wall = cb_wall + gt_wall
        total_peak_mb = max(cb_peak, gt_peak) / 1024 / 1024
        rows.append(Measurement(CALLER, n, total_wall, total_peak_mb, gt_ec))
        status = "ok" if gt_ec == 0 else f"FAILED (exit {gt_ec})"
        print(f"  TOTAL    -> {total_wall:.1f}s, peak {total_peak_mb:.0f} MB  [{status}]")

    tsv = PERF_DIR / f"{CALLER}.tsv"
    write_tsv(rows, tsv)
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
