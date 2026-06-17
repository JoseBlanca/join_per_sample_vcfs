//! `ssr-catalog` subcommand — build the per-genome SSR locus catalog (Stage 0
//! of the SSR caller). Detects tandem repeats in a reference with `trf-mod`,
//! post-processes them to di/tri/tetra+ microsatellite loci, and writes a
//! self-describing bgzip TSV catalog. See [`crate::ssr_mark1::catalog`].

use std::path::PathBuf;

use clap::Args;
use thiserror::Error;

use crate::ssr_mark1::catalog::{self, CatalogConfig, CatalogParams};

/// Arguments for the `ssr-catalog` subcommand. The struct is the authoritative
/// knob list; [`run_ssr_catalog`] translates it into a [`CatalogConfig`].
#[derive(Debug, Args, Clone)]
pub struct SsrCatalogArgs {
    /// Reference FASTA to scan for tandem repeats.
    #[arg(long)]
    pub reference: PathBuf,

    /// Output catalog path (a bgzip TSV, e.g. `catalog.ssr_catalog.bed.gz`).
    #[arg(long)]
    pub output: PathBuf,

    /// Explicit path to the `trf-mod` binary. If unset it is discovered beside
    /// this executable, then on `PATH`.
    #[arg(long)]
    pub trf_mod_path: Option<PathBuf>,

    /// Disk-backed scratch root for per-contig `trf-mod` temp files
    /// (CWD-relative by default; keep it off RAM-backed `/tmp`).
    #[arg(long, default_value = "ssr-catalog-tmp")]
    pub temp_dir: PathBuf,

    /// Flank margin (bp) embedded each side of every tract in `ref_seq`.
    #[arg(long, default_value_t = 50, help_heading = "Advanced")]
    pub flank_bp: u32,

    /// Drop a locus (and its whole cluster) if any other detected repeat is
    /// within this many bp. Must be `>= --flank-bp` for clean survivor flanks.
    #[arg(long, default_value_t = 50, help_heading = "Advanced")]
    pub bundle_threshold: u32,

    /// Purity floor in `[0, 1]`: a degeneracy cutoff applied to the recomputed
    /// per-tract purity (imperfect-but-above-floor loci are kept).
    #[arg(long, default_value_t = 0.8, help_heading = "Advanced")]
    pub min_purity: f32,

    /// Early accept-gate on TRF's alignment score (`trf-mod` already applies its
    /// own `-s 30`; this is an extra floor, `0` = off).
    #[arg(long, default_value_t = 0, help_heading = "Advanced")]
    pub min_score: i32,
}

/// Errors from the `ssr-catalog` subcommand.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum SsrCatalogCliError {
    /// `--bundle-threshold` is smaller than `--flank-bp`, which would let a
    /// survivor's embedded flank overlap a neighbouring repeat.
    #[error(
        "--bundle-threshold ({bundle}) must be >= --flank-bp ({flank}) so survivor flanks are clean"
    )]
    BundleThresholdTooSmall { bundle: u32, flank: u32 },

    /// The catalog build itself failed. The boxed source carries the typed
    /// `CatalogError` cause (kept boxed so this public error does not name the
    /// crate-internal type).
    #[error("catalog build failed")]
    Catalog {
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },
}

/// Build an SSR catalog from a reference FASTA.
pub fn run_ssr_catalog(args: &SsrCatalogArgs) -> Result<(), SsrCatalogCliError> {
    if args.bundle_threshold < args.flank_bp {
        return Err(SsrCatalogCliError::BundleThresholdTooSmall {
            bundle: args.bundle_threshold,
            flank: args.flank_bp,
        });
    }

    // The library reads no clock; the CLI stamps the build date (YYYY-MM-DD).
    let date = super::common::rfc3339_now()
        .split('T')
        .next()
        .unwrap_or_default()
        .to_string();

    let cfg = CatalogConfig {
        reference: args.reference.clone(),
        output: args.output.clone(),
        trf_mod_path: args.trf_mod_path.clone(),
        temp_dir: args.temp_dir.clone(),
        params: CatalogParams {
            min_purity: args.min_purity,
            min_score: args.min_score,
            flank_bp: args.flank_bp,
            bundle_threshold: args.bundle_threshold,
        },
        tool_version: env!("CARGO_PKG_VERSION").to_string(),
        date,
    };

    catalog::run(&cfg).map_err(|e| SsrCatalogCliError::Catalog {
        source: Box::new(e),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(flank_bp: u32, bundle_threshold: u32) -> SsrCatalogArgs {
        SsrCatalogArgs {
            reference: "ref.fa".into(),
            output: "out.bed.gz".into(),
            trf_mod_path: None,
            temp_dir: "tmp".into(),
            flank_bp,
            bundle_threshold,
            min_purity: 0.8,
            min_score: 0,
        }
    }

    /// The `bundle_threshold >= flank_bp` clean-flank invariant is enforced at
    /// the CLI, before any FASTA/trf-mod work.
    #[test]
    fn rejects_bundle_threshold_below_flank() {
        let err = run_ssr_catalog(&args(50, 40)).unwrap_err();
        assert!(
            matches!(
                err,
                SsrCatalogCliError::BundleThresholdTooSmall {
                    bundle: 40,
                    flank: 50
                }
            ),
            "got {err:?}"
        );
    }

    /// The subcommand parses through the top-level CLI with its defaults.
    #[test]
    fn parses_through_the_top_level_cli() {
        use crate::pop_var_caller::cli::{Cli, PopVarCallerCommand};
        use clap::Parser;
        let cli = Cli::try_parse_from([
            "pop_var_caller",
            "ssr-catalog",
            "--reference",
            "r.fa",
            "--output",
            "o.bed.gz",
        ])
        .unwrap();
        match cli.cmd {
            PopVarCallerCommand::SsrCatalog(a) => {
                assert_eq!(a.flank_bp, 50);
                assert_eq!(a.bundle_threshold, 50);
                assert_eq!(a.min_purity, 0.8);
            }
            other => panic!("expected SsrCatalog, got {other:?}"),
        }
    }
}
