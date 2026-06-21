//! The `ssr-call` driver — open the cohort, run the merge, write the output (arch
//! `doc/devel/architecture/ssr_call_reading.md` §5).
//!
//! **Phase 3: single-threaded.** The producer→queue→worker-pool→writer topology and
//! the genotyping EM are not built here. The worker would exist to overlap *expensive
//! EM work* with decode; with no EM yet (the genotyping doc owns it) and decode
//! parallelism being the separate Phase-5 prefetch pool, a thread pool around a
//! formatting stub would be unverifiable complexity. So the driver streams the merged
//! `CohortLocus`es straight to a **catalog-ordered TSV dump** of the reading layer — a
//! placeholder until the EM + VCF land, and a useful way to inspect a cohort's merged
//! evidence. `--threads` / `--queue-depth` are accepted but reserved.

use std::fs::File;
use std::io::{BufWriter, Read, Seek, Write};
use std::path::PathBuf;

use crate::ssr::cohort::merge::{CohortMerger, SsrMergeError};
use crate::ssr::cohort::types::CohortLocus;

/// Inputs for an `ssr-call` run.
pub(crate) struct SsrCallConfig {
    /// The shared `.ssr.catalog`.
    pub(crate) catalog: PathBuf,
    /// The per-sample `.ssr.psp` evidence files.
    pub(crate) psp_files: Vec<PathBuf>,
    /// Output path (currently the TSV dump; the VCF lands with the EM).
    pub(crate) output: PathBuf,
    /// Reserved — the EM worker pool is not built yet (Phase 3 is single-threaded).
    pub(crate) threads: usize,
    /// Reserved — the bounded producer→worker queue is not built yet.
    pub(crate) queue_depth: usize,
}

/// Errors from an `ssr-call` run.
#[derive(Debug, thiserror::Error)]
pub(crate) enum SsrCallError {
    /// Opening / merging the cohort failed.
    #[error(transparent)]
    Merge(#[from] SsrMergeError),
    /// Writing the output failed.
    #[error("writing ssr-call output")]
    Write(#[from] std::io::Error),
}

/// Open the cohort, merge, and write the catalog-ordered dump to `config.output`.
pub(crate) fn run(config: &SsrCallConfig) -> Result<(), SsrCallError> {
    let merger = CohortMerger::open(&config.catalog, &config.psp_files)?;
    let chrom_names = merger.chrom_names();
    let mut out = BufWriter::new(File::create(&config.output)?);
    write_dump(merger, &chrom_names, &mut out)?;
    out.flush()?;
    Ok(())
}

/// Stream every emitted `CohortLocus` to `out` as one TSV row, in catalog order.
/// Generic over the merge sources so it is testable in memory.
fn write_dump<R: Read + Seek, C: Read, W: Write>(
    merger: CohortMerger<R, C>,
    chrom_names: &[String],
    out: &mut W,
) -> Result<(), SsrCallError> {
    writeln!(out, "#chrom\tstart\tend\tmotif\tn_present\tsamples")?;
    for item in merger {
        let (_seq, cohort) = item?;
        writeln!(out, "{}", format_locus(&cohort, chrom_names))?;
    }
    Ok(())
}

/// One TSV row for a locus: coordinates + motif + the present samples, each as
/// `idx:depth=…,alleles=…`. Pure, so the formatting is unit-testable on its own.
fn format_locus(cohort: &CohortLocus, chrom_names: &[String]) -> String {
    let chrom = chrom_names
        .get(cohort.locus.chrom_id as usize)
        .map(String::as_str)
        .unwrap_or("?");
    let motif = String::from_utf8_lossy(cohort.motif.as_bytes());
    let samples = cohort
        .present
        .iter()
        .zip(&cohort.samples)
        .map(|(idx, evidence)| {
            format!(
                "{idx}:depth={},alleles={}",
                evidence.qc.depth,
                evidence.seq_counts.len()
            )
        })
        .collect::<Vec<_>>()
        .join(";");
    format!(
        "{chrom}\t{}\t{}\t{motif}\t{}\t{samples}",
        cohort.locus.start,
        cohort.locus.end,
        cohort.present_count(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::test_support::{
        REF_MD5, catalog_bytes, loc, obs, reader, rec, ssr_header, ssr_psp,
    };
    use crate::ssr::cohort::types::{LocusId, SampleEvidence, SsrQc};
    use crate::ssr::types::Motif;
    use std::io::Cursor;

    fn cohort_locus(chrom_id: u32, start: u32, present: &[(u32, u32, usize)]) -> CohortLocus {
        let mut cl = CohortLocus::new(
            LocusId {
                chrom_id,
                start,
                end: start + 6,
            },
            Motif::new(b"CA").unwrap(),
            Box::from(b"GGGGGGCACACATTTTTT".as_slice()),
        );
        for &(idx, depth, n_alleles) in present {
            cl.push(
                idx,
                SampleEvidence {
                    seq_counts: (0..n_alleles)
                        .map(|i| (vec![b'C', b'A', i as u8].into_boxed_slice(), 1))
                        .collect(),
                    qc: SsrQc {
                        depth,
                        ..SsrQc::default()
                    },
                },
            );
        }
        cl
    }

    #[test]
    fn format_locus_renders_chrom_motif_and_samples() {
        let names = vec!["chr1".to_string()];
        let cl = cohort_locus(0, 16, &[(0, 30, 1), (2, 12, 3)]);
        assert_eq!(
            format_locus(&cl, &names),
            "chr1\t16\t22\tCA\t2\t0:depth=30,alleles=1;2:depth=12,alleles=3"
        );
    }

    #[test]
    fn format_locus_falls_back_when_chrom_id_unknown() {
        let cl = cohort_locus(7, 16, &[(0, 30, 1)]);
        assert!(format_locus(&cl, &[]).starts_with("?\t16\t22\tCA\t1\t"));
    }

    #[test]
    fn write_dump_emits_header_and_one_row_per_locus_in_order() {
        let loci = [loc("chr1", 16), loc("chr1", 60), loc("chr1", 100)];
        let catalog =
            crate::ssr::catalog::io::CatalogReader::new(Cursor::new(catalog_bytes(REF_MD5, &loci)))
                .unwrap();
        // Sample 0 covers 16 + 100; sample 1 covers 60 + 100.
        let a = reader(ssr_psp(
            ssr_header(&["chr1"], REF_MD5),
            &[
                rec(0, 16, obs(&[(b"CACACA", 4)])),
                rec(0, 100, obs(&[(b"CA", 2)])),
            ],
        ));
        let b = reader(ssr_psp(
            ssr_header(&["chr1"], REF_MD5),
            &[
                rec(0, 60, obs(&[(b"CACACACA", 5)])),
                rec(0, 100, obs(&[(b"CA", 3)])),
            ],
        ));
        let merger =
            CohortMerger::from_parts(catalog, vec![("A".into(), a), ("B".into(), b)]).unwrap();
        let names = merger.chrom_names();

        let mut out = Vec::new();
        write_dump(merger, &names, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();

        assert_eq!(
            text,
            "#chrom\tstart\tend\tmotif\tn_present\tsamples\n\
             chr1\t16\t22\tCA\t1\t0:depth=30,alleles=1\n\
             chr1\t60\t66\tCA\t1\t1:depth=30,alleles=1\n\
             chr1\t100\t106\tCA\t2\t0:depth=30,alleles=1;1:depth=30,alleles=1\n"
        );
    }

    #[test]
    fn run_reads_files_and_writes_the_dump() {
        let dir = tempfile::TempDir::new().unwrap();
        let catalog_path = dir.path().join("c.ssr.catalog");
        let a_path = dir.path().join("a.ssr.psp");
        let out_path = dir.path().join("out.tsv");

        std::fs::write(&catalog_path, catalog_bytes(REF_MD5, &[loc("chr1", 16)])).unwrap();
        std::fs::write(
            &a_path,
            ssr_psp(
                ssr_header(&["chr1"], REF_MD5),
                &[rec(0, 16, obs(&[(b"CACACA", 4)]))],
            ),
        )
        .unwrap();

        let config = SsrCallConfig {
            catalog: catalog_path,
            psp_files: vec![a_path],
            output: out_path.clone(),
            threads: 4,
            queue_depth: 4,
        };
        run(&config).unwrap();

        let text = std::fs::read_to_string(&out_path).unwrap();
        assert_eq!(
            text,
            "#chrom\tstart\tend\tmotif\tn_present\tsamples\n\
             chr1\t16\t22\tCA\t1\t0:depth=30,alleles=1\n"
        );
    }

    #[test]
    fn run_errors_on_a_missing_catalog() {
        let config = SsrCallConfig {
            catalog: PathBuf::from("/no/such/catalog.ssr.catalog"),
            psp_files: vec![PathBuf::from("/no/such/a.ssr.psp")],
            output: PathBuf::from("/tmp/unused.tsv"),
            threads: 4,
            queue_depth: 4,
        };
        assert!(matches!(run(&config), Err(SsrCallError::Merge(_))));
    }
}
