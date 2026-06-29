#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["polars", "numpy", "pyfaidx"]
# ///
"""Step 2 of option (c): assemble the per-sample window coverage TSVs into one
table and attach per-window GC% from the reference FASTA (no BAM access).

Output: results/window_cov.w<W>.parquet  (long form, one row per window x sample)
  chrom, win_start, sample, n_cov, sum_depth, mean_depth, breadth, gc

We only ATTACH GC here; we do not yet normalise depth by it. The point is to
look at the raw depth-vs-GC relationship first and decide the normalisation."""
import argparse
import glob
import os
import numpy as np
import polars as pl
from pyfaidx import Fasta

REF = os.path.expanduser("~/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa")


def gc_of(seq):
    s = seq.upper()
    g = s.count("G") + s.count("C")
    at = s.count("A") + s.count("T")
    n = g + at
    return g / n if n else np.nan


def main(w, indir, out):
    files = sorted(glob.glob(f"{indir}/*.w{w}.tsv"))
    if not files:
        raise SystemExit(f"no *.w{w}.tsv in {indir} — run build_window_coverage.sh first")
    print(f"loading {len(files)} sample window files (W={w})")

    frames = []
    for f in files:
        sample = os.path.basename(f).replace(f".w{w}.tsv", "")
        df = pl.read_csv(f, separator="\t", has_header=False,
                         new_columns=["chrom", "win_start", "n_cov", "sum_depth"])
        df = df.with_columns(pl.lit(sample).alias("sample"))
        frames.append(df)
    long = pl.concat(frames)
    long = long.with_columns([
        (pl.col("sum_depth") / pl.col("n_cov")).alias("mean_depth"),
        (pl.col("n_cov") / w).alias("breadth"),
    ])

    # GC per unique window from the FASTA (computed once, then joined)
    wins = long.select(["chrom", "win_start"]).unique().sort(["chrom", "win_start"])
    print(f"computing GC for {wins.height} unique windows")
    fa = Fasta(REF, sequence_always_upper=True)
    gc = []
    for chrom, start in wins.iter_rows():
        seq = fa[chrom][start - 1: start - 1 + w]
        gc.append(gc_of(str(seq)))
    wins = wins.with_columns(pl.Series("gc", gc))

    out_df = long.join(wins, on=["chrom", "win_start"], how="left")
    out_df.write_parquet(out)
    print(f"wrote {out}: {out_df.height:,} rows "
          f"({out_df['sample'].n_unique()} samples x {wins.height} windows)")
    print(out_df.select(["mean_depth", "breadth", "gc"]).describe())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--w", type=int, default=2000)
    ap.add_argument("--indir", default="benchmarks/tomato2/results/window_cov")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or f"benchmarks/tomato2/results/window_cov.w{a.w}.parquet"
    main(a.w, a.indir, out)
