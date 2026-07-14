//! Stage 0 — `ssr-catalog`: build the per-genome SSR locus catalog.
//!
//! The catalog is the SSR pipeline's first artefact: one self-describing,
//! bgzip-wrapped TSV listing every short-tandem-repeat locus in a reference,
//! each row carrying the tract coordinates, the repeat motif, a recomputed
//! purity, and the embedded local reference (`ref_seq` + `ref_seq_start`). It
//! is the *only* reference-bearing input the downstream stages need
//! (architecture [`ssr_catalog.md`](../../../doc/devel/architecture/ssr_catalog.md)).
//!
//! ```text
//! reference FASTA ─► [TRF-mod detect] ─► [post-process] ─► [embed ref_seq] ─► catalog.ssr_catalog.bed.gz (+ index)
//! ```
//!
//! **Build status (incremental).** The format I/O layer ([`io`]) lands first —
//! it is the cross-stage contract Stage 1's `fetch_reads` reader consumes, and
//! it is buildable and testable without the external `trf-mod` binary. The
//! detection/post-processing front-end follows:
//!
//! - [`io`] — **built**: [`io::CatalogHeader`], [`io::CatalogWriter`],
//!   [`io::CatalogReader`], and the `Locus` ⇄ row serialisation + round-trip.
//! - [`postprocess`] — **built**: [`postprocess::build_loci`] — the period≤6 →
//!   drop-compound → drop-bundle → end-trim → recompute-purity → embed-`ref_seq`
//!   pipeline (faithful GangSTR `minimal_trim`/`remove_bundles` port).
//! - [`trf`] — **built**: [`trf::locate_trf_mod`], [`trf::version`],
//!   [`trf::run_on_contig`] (temp-file spawn, no pipes), and
//!   [`trf::parse_bed_line`] (the 10-column BED) → [`trf::TrfRecord`].
//! - `run()` orchestrator + the `ssr-catalog` CLI subcommand — *pending* (the
//!   last Stage-0 piece: per-contig fan-out/collect + header build + CSI index).

pub mod io;
pub mod postprocess;
pub mod trf;

#[cfg(test)]
mod scanner_parity;

/// Build/accept parameters — the post-process knobs that both drive
/// [`postprocess::build_loci`] and are recorded in the catalog header
/// ([`io::CatalogHeader`]) so a reader sees exactly how it was built. The
/// per-period copy-number floors are fixed constants, not a knob.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CatalogParams {
    /// Purity floor applied after recomputation (a degeneracy cutoff in
    /// `[0, 1]`); imperfect-but-above-floor loci are kept.
    pub min_purity: f32,
    /// Early accept-gate on TRF's `score`; records below are dropped.
    pub min_score: i32,
    /// Flank margin (bp) embedded each side of the tract in `ref_seq`.
    pub flank_bp: u32,
    /// Bundle-drop radius (bp). `>= flank_bp` guarantees clean survivor flanks.
    pub bundle_threshold: u32,
}

impl Default for CatalogParams {
    /// Pinned Stage-0 defaults. `min_score = 0` leaves filtering to trf-mod's
    /// own `-s 30`; `bundle_threshold == flank_bp` satisfies the
    /// `bundle_threshold >= flank_bp` clean-flank invariant.
    fn default() -> Self {
        Self {
            min_purity: 0.8,
            min_score: 0,
            flank_bp: 50,
            bundle_threshold: 50,
        }
    }
}

/// Errors building or reading an SSR catalog (Stage 0).
///
/// `#[non_exhaustive]`: the detection / post-processing front-end will add
/// `trf-mod`-spawn and FASTA variants in later increments; only the format-I/O
/// variants exist today.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub(crate) enum CatalogError {
    /// An underlying read/write/bgzf failure. `context` names the operation.
    #[error("catalog I/O failed ({context})")]
    Io {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },

    /// The `##` metadata header is missing a required key or carries an
    /// unparsable value.
    #[error("malformed catalog header: {reason}")]
    HeaderParse { reason: String },

    /// A data row does not have the expected tab-separated column shape, or a
    /// numeric field failed to parse.
    #[error("malformed catalog row at line {line}: {reason}")]
    RowParse { line: usize, reason: String },

    /// A row's `motif` field is not a valid SSR period.
    #[error("invalid motif at line {line}")]
    InvalidMotif {
        line: usize,
        #[source]
        source: crate::ssr::types::MotifError,
    },

    /// A row's coordinates / purity violate the [`Locus`](crate::ssr::types::Locus)
    /// invariants.
    #[error("invalid locus at line {line}")]
    InvalidLocus {
        line: usize,
        #[source]
        source: crate::ssr::types::LocusError,
    },

    /// `trf-mod` could not be located: no usable override, not found beside our
    /// own executable, and not on `PATH` (architecture §2.4 layered discovery).
    #[error("trf-mod binary not found (checked: override, exe-sibling, PATH)")]
    TrfModNotFound,

    /// Spawning or waiting on `trf-mod` failed at the OS level.
    #[error("failed to run trf-mod ({context})")]
    TrfSpawn {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },

    /// `trf-mod` exited unsuccessfully while processing a contig.
    #[error("trf-mod failed on contig {contig:?} ({status})")]
    TrfRun { contig: String, status: String },

    /// `trf-mod -v` produced no recognisable version line.
    #[error("could not parse a version from `trf-mod -v`")]
    TrfVersion,

    /// A `trf-mod` BED line did not match the expected 10-column layout, a field
    /// failed to parse, or the contig-name column disagreed.
    #[error("malformed trf-mod BED at line {line}: {reason}")]
    TrfParse { line: usize, reason: String },

    /// Reading the reference FASTA failed (open or record decode).
    #[error("reading the reference FASTA failed ({context})")]
    Fasta {
        context: &'static str,
        #[source]
        source: std::io::Error,
    },
}

// ---------------------------------------------------------------------
// Orchestrator
// ---------------------------------------------------------------------

use std::path::PathBuf;

use io::{CatalogHeader, CatalogWriter};

/// Inputs for one catalog build. Plain config (the `ssr-catalog` CLI maps its
/// clap args onto this); `date` is caller-supplied so this module reads no clock.
#[derive(Debug, Clone)]
pub(crate) struct CatalogConfig {
    /// Reference FASTA to scan.
    pub reference: PathBuf,
    /// Output catalog path (a bgzip TSV).
    pub output: PathBuf,
    /// Optional explicit `trf-mod` path (else discovered, §2.4).
    pub trf_mod_path: Option<PathBuf>,
    /// Disk-backed root for per-contig trf-mod temp files (NOT tmpfs).
    pub temp_dir: PathBuf,
    /// Build/accept parameters.
    pub params: CatalogParams,
    /// Building tool version, recorded in the header.
    pub tool_version: String,
    /// Build date (`YYYY-MM-DD`), recorded in the header.
    pub date: String,
}

/// Build a catalog: detect repeats per contig with `trf-mod`, post-process to
/// loci, and write the bgzip TSV. **Single-threaded** for now — the per-contig
/// rayon fan-out (architecture §8) is a follow-up; the contig order is the
/// reference's, so output is already coordinate-sorted.
pub(crate) fn run(cfg: &CatalogConfig) -> Result<(), CatalogError> {
    use md5::{Digest, Md5};
    use std::fs::File;
    use std::io::BufReader;

    let bin = trf::locate_trf_mod(cfg.trf_mod_path.as_deref())?;
    let trf_version = trf::version(&bin)?;
    std::fs::create_dir_all(&cfg.temp_dir).map_err(|source| CatalogError::Io {
        context: "create catalog temp dir",
        source,
    })?;

    // Stream the reference: per contig, detect + post-process into loci, and
    // fold the upper-cased bases into the whole-reference md5 (the spec's
    // upper-cased-content convention — recomputable by Stages 1-2 for an
    // integrity check).
    let file = File::open(&cfg.reference).map_err(|source| CatalogError::Fasta {
        context: "open reference",
        source,
    })?;
    let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));
    let mut md5 = Md5::new();
    let mut loci = Vec::new();
    for result in reader.records() {
        let rec = result.map_err(|source| CatalogError::Fasta {
            context: "decode reference record",
            source,
        })?;
        let name = String::from_utf8_lossy(rec.name()).into_owned();
        let seq = rec.sequence().as_ref();
        for &b in seq {
            md5.update([b.to_ascii_uppercase()]);
        }
        let recs = trf::run_on_contig(&bin, &name, seq, &cfg.temp_dir)?;
        loci.extend(postprocess::build_loci(recs, &name, seq, &cfg.params));
    }
    let reference_md5 = crate::pop_var_caller::common::format_md5_hex(md5.finalize().into());

    let header = CatalogHeader {
        tool_version: cfg.tool_version.clone(),
        reference: cfg.reference.display().to_string(),
        reference_md5,
        trf_mod_version: trf_version,
        params: cfg.params.clone(),
        date: cfg.date.clone(),
    };
    let sink = File::create(&cfg.output).map_err(|source| CatalogError::Io {
        context: "create catalog output",
        source,
    })?;
    let mut writer = CatalogWriter::new(sink, &header)?;
    for locus in &loci {
        writer.write_locus(locus)?;
    }
    let file = writer.finish()?;
    file.sync_all().map_err(|source| CatalogError::Io {
        context: "sync catalog output",
        source,
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::catalog::io::CatalogReader;
    use std::io::Write;

    /// End-to-end: write a tiny reference with one clean (CAG)*40 tract, run the
    /// orchestrator (driving the installed trf-mod), and read the catalog back.
    /// Skips (does not fail) if trf-mod is unavailable — exercised in the dev
    /// container.
    #[test]
    fn run_builds_a_catalog_for_a_synthetic_reference() {
        if trf::locate_trf_mod(None).is_err() {
            eprintln!("skipping run() test: trf-mod not found on this host");
            return;
        }
        std::fs::create_dir_all("tmp").unwrap();
        let dir = tempfile::tempdir_in("tmp").unwrap();
        let ref_path = dir.path().join("ref.fa");
        let out_path = dir.path().join("catalog.bed.gz");

        // (CAG)*40 (120 bp) bounded by 100 bp poly-T flanks. trf-mod reports the
        // flanks as period-1 homopolymers, but they are dropped (MIN_PERIOD = 2)
        // before bundling, so the CAG tract survives cleanly at [100, 220) — this
        // exercises the homopolymer-drop decision end-to-end.
        let mut seq = vec![b'T'; 100];
        for _ in 0..40 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend(std::iter::repeat_n(b'T', 100));
        {
            let mut f = std::fs::File::create(&ref_path).unwrap();
            writeln!(f, ">ctg1").unwrap();
            f.write_all(&seq).unwrap();
            writeln!(f).unwrap();
        }

        let cfg = CatalogConfig {
            reference: ref_path,
            output: out_path.clone(),
            trf_mod_path: None,
            temp_dir: PathBuf::from("tmp"),
            params: CatalogParams::default(),
            tool_version: "0.0.0-test".to_string(),
            date: "2026-06-16".to_string(),
        };
        run(&cfg).expect("catalog build");

        let mut reader = CatalogReader::new(std::fs::File::open(&out_path).unwrap()).unwrap();
        assert_eq!(reader.header().reference_md5.len(), 32);
        assert!(reader.header().trf_mod_version.contains("Version"));
        let loci = reader.read_all().unwrap();
        // No period-1 homopolymer survived; the CAG tract did, at exact coords.
        assert!(
            loci.iter().all(|l| l.period() >= 2),
            "period-1 homopolymers are excluded"
        );
        let cag = loci
            .iter()
            .find(|l| l.motif().as_bytes() == b"CAG")
            .expect("the CAG locus is in the catalog");
        assert_eq!(cag.chrom(), "ctg1");
        assert_eq!(cag.start(), 100);
        assert_eq!(cag.end(), 220);
        assert_eq!(cag.purity_fraction(), 1.0);
    }
}
