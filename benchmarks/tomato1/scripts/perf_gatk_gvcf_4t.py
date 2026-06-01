# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "psutil",
# ]
# ///
"""Section-2 perf: build ONE sample's per-sample intermediate (a GVCF)
with the full thread budget — GATK `HaplotypeCaller -ERC GVCF
--native-pair-hmm-threads THREADS` on a single CRAM.

GATK half of the "create one per-sample intermediate file, one sample,
THREADS threads" bar-chart comparison; the our-caller half is
`perf_ours_psp_4t.py` (pileup -> .psp). Unlike
`perf_gatk_haplotype_caller.py` (which runs N single-pair-hmm-thread
HaplotypeCallers concurrently for the cohort-scaling panel), here a
single HaplotypeCaller process gets all THREADS — the realistic "build
one new sample's GVCF as fast as I can" operation.

Single process, so wall and peak RSS characterise HaplotypeCaller itself.

Note on availability: GATK isn't in Homebrew. On macOS run this script
in the container; on Linux either path works.

Inputs:  benchmarks/tomato1/crams/*.bench.cram (first sample by name)
         benchmarks/tomato1/regions.bed
Output:  benchmarks/tomato1/results/perf/gatk_gvcf_4t.tsv
         benchmarks/tomato1/results/perf/gatk_gvcf_4t/<sample>.g.vcf.gz

Env overrides:
  GATK_BIN    binary (default: /opt/gatk/gatk)
  REFERENCE   SL4.0 fasta (.fai + .dict siblings required)
  JAVA_HEAP   -Xmx for the HaplotypeCaller process (default: 4g)
  THREADS     --native-pair-hmm-threads (default: 4)

Invoke (container path on macOS):
  DEV_EXTRA_MOUNT=$HOME/genomes ./scripts/dev.sh \\
      python3 benchmarks/tomato1/scripts/perf_gatk_gvcf_4t.py
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from perf_common import (  # noqa: E402
    CRAM_DIR, DEFAULT_BED, DEFAULT_REFERENCE, DEFAULT_THREADS, PERF_DIR,
    Measurement, check_exists, list_inputs, measure, write_tsv,
)

CALLER = "gatk_gvcf_4t"
BIN = Path(os.environ.get("GATK_BIN", "/opt/gatk/gatk"))
JAVA_HEAP = os.environ.get("JAVA_HEAP", "4g")
OUT_DIR = PERF_DIR / CALLER


def main() -> int:
    ref_dict = DEFAULT_REFERENCE.with_suffix(".dict")
    check_exists(
        BIN, DEFAULT_REFERENCE, Path(str(DEFAULT_REFERENCE) + ".fai"),
        ref_dict, DEFAULT_BED,
    )
    crams = list_inputs(CRAM_DIR, ".bench.cram")
    cram = crams[0]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(
        f"=== {CALLER}: HaplotypeCaller -ERC GVCF one sample @ "
        f"--native-pair-hmm-threads {DEFAULT_THREADS} ==="
    )
    print(f"sample:    {cram.name}")
    print(f"reference: {DEFAULT_REFERENCE}")
    print()

    base = cram.name.removesuffix(".bench.cram")
    gvcf = OUT_DIR / f"{base}.g.vcf.gz"
    cmd = [
        str(BIN), "--java-options", f"-Xmx{JAVA_HEAP}", "HaplotypeCaller",
        "--reference", str(DEFAULT_REFERENCE),
        "--input", str(cram),
        "--intervals", str(DEFAULT_BED),
        "--output", str(gvcf),
        "--emit-ref-confidence", "GVCF",
        "--native-pair-hmm-threads", str(DEFAULT_THREADS),
    ]
    print(f"[{CALLER}] HaplotypeCaller {cram.name} -> {gvcf.name}")
    wall, peak_bytes, exit_code = measure(cmd)
    peak_mb = peak_bytes / 1024 / 1024
    status = "ok" if exit_code == 0 else f"FAILED (exit {exit_code})"
    print(f"  -> {wall:.1f}s, peak {peak_mb:.0f} MB  [{status}]")

    rows = [Measurement(CALLER, 1, wall, peak_mb, exit_code)]
    tsv = PERF_DIR / f"{CALLER}.tsv"
    write_tsv(rows, tsv)
    print(f"\nwrote {tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
