//! The k-way merger — the cohort producer (arch
//! `doc/devel/architecture/ssr_call_reading.md` §4).
//!
//! Drives off the **catalog** (the authoritative, ordered master locus list) and, for
//! each catalog locus, asks **every** [`SampleEvidenceCursor`] for its evidence,
//! gathering the present samples into one [`CohortLocus`]. A locus no sample covers is
//! dropped (sparse-omit); each emitted locus carries a monotonic sequence number for
//! the writer to reorder by. One `CohortLocus` is produced at a time — no batching
//! (Phase 4 of the build can add it if worker imbalance ever shows up).
//!
//! It also owns the **same-catalog precondition**: every input must have been built
//! against the catalog the merger holds (the `catalog_reference_md5` header parameter
//! must match the catalog's `reference_md5`), and the **chromosome-id reconciliation**
//! (each file's per-file chromosome ids → cohort-global ids, taken from the shared
//! chromosome table).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek};
use std::path::{Path, PathBuf};

use crate::psp::header::ParameterValue;
use crate::psp::registry_ssr::CATALOG_REFERENCE_MD5_PARAM;
use crate::psp::{PspReadError, PspReader};
use crate::ssr::catalog::CatalogError;
use crate::ssr::catalog::io::CatalogReader;
use crate::ssr::cohort::reader::{SampleEvidenceCursor, SsrCohortReadError};
use crate::ssr::cohort::types::{CohortLocus, LocusId};

/// Errors from assembling or running the cohort merge.
#[derive(Debug, thiserror::Error)]
pub(crate) enum SsrMergeError {
    /// No input `.ssr.psp` files were given.
    #[error("no input .ssr.psp files given")]
    NoInputs,
    /// Reading the catalog failed.
    #[error("reading the catalog")]
    Catalog(#[from] CatalogError),
    /// An input is not an SSR `.psp` (wrong schema `kind`).
    #[error("input {path:?} is not an SSR .psp (kind = {kind:?})")]
    NotSsr { path: String, kind: String },
    /// An input is missing the catalog-binding header parameter.
    #[error("input {path:?} is missing the {CATALOG_REFERENCE_MD5_PARAM:?} header parameter")]
    MissingCatalogMd5 { path: String },
    /// An input was built against a different catalog.
    #[error(
        "input {path:?} was built against a different catalog \
         (reference md5 {found:?}, expected {expected:?})"
    )]
    CatalogMismatch {
        path: String,
        expected: String,
        found: String,
    },
    /// An input names a chromosome absent from the cohort chromosome table.
    #[error(
        "input {path:?} declares chromosome {chrom:?}, absent from the cohort chromosome table"
    )]
    UnknownInputChrom { path: String, chrom: String },
    /// A catalog locus sits on a chromosome no input declares.
    #[error("catalog locus on chromosome {chrom:?}, which no input declares")]
    UnknownCatalogChrom { chrom: String },
    /// The catalog is not strictly coordinate-sorted in the cohort-global chromosome
    /// order — the merge (and every cursor's monotonic contract) requires it.
    #[error("catalog is not coordinate-sorted: locus {at:?} does not follow {prev:?}")]
    UnsortedCatalog { prev: LocusId, at: LocusId },
    /// Opening an input file failed.
    #[error("opening {path:?}")]
    Open {
        path: String,
        #[source]
        source: std::io::Error,
    },
    /// Parsing a `.ssr.psp` header failed.
    #[error("opening .ssr.psp {path:?}")]
    OpenPsp {
        path: String,
        #[source]
        source: PspReadError,
    },
    /// Reading a sample's evidence failed — at open or while streaming. Names the
    /// offending input so a cohort of many files points at the right one.
    #[error("reading sample {sample:?}")]
    Read {
        sample: String,
        #[source]
        source: SsrCohortReadError,
    },
}

/// The cohort producer: a catalog-driven k-way merge over N per-sample cursors,
/// yielding `(sequence number, CohortLocus)` in catalog order.
pub(crate) struct CohortMerger<R: Read + Seek, C: Read> {
    /// The merge driver — the authoritative ordered locus list.
    catalog: CatalogReader<C>,
    /// One cursor per input sample, indexed by cohort sample index.
    cursors: Vec<SampleEvidenceCursor<R>>,
    /// Per-sample labels (paths/names), parallel to `cursors` — for error messages.
    labels: Vec<String>,
    /// Chromosome name → cohort-global id (the canonical order, taken from the shared
    /// chromosome table). Maps catalog loci into the `LocusId` frame.
    name_to_global: HashMap<String, u32>,
    /// The previous catalog locus visited — to enforce that the catalog is strictly
    /// coordinate-sorted (the cursor's monotonic precondition).
    prev_locus: Option<LocusId>,
    /// Next emitted-locus sequence number (monotonic over *emitted* loci).
    next_seq: u64,
}

impl<R: Read + Seek, C: Read> CohortMerger<R, C> {
    /// Validate the inputs against the catalog and build the merger.
    ///
    /// `inputs` are `(label, reader)` pairs — the label is a path or name used only in
    /// error messages. Validates, per input: the SSR `kind`, the `catalog_reference_md5`
    /// match, and that every declared chromosome is in the cohort table; then builds a
    /// per-file chromosome remap and opens a cursor over it.
    ///
    /// The cohort-global chromosome order is taken from the **first** input's table; the
    /// same-catalog precondition makes every input's table identical, so the remaps are
    /// effectively identity — but they are built by name so a divergent table is caught.
    pub(crate) fn from_parts(
        catalog: CatalogReader<C>,
        inputs: Vec<(String, PspReader<R>)>,
    ) -> Result<Self, SsrMergeError> {
        let Some(first) = inputs.first() else {
            return Err(SsrMergeError::NoInputs);
        };
        let expected_md5 = catalog.header().reference_md5.clone();

        // Canonical chromosome order — global id = index in the first input's table.
        let name_to_global: HashMap<String, u32> = first
            .1
            .header()
            .chromosomes
            .iter()
            .enumerate()
            .map(|(idx, chrom)| (chrom.name.clone(), idx as u32))
            .collect();

        let mut cursors = Vec::with_capacity(inputs.len());
        let mut labels = Vec::with_capacity(inputs.len());
        for (path, reader) in inputs {
            let header = reader.header();
            if header.kind != "ssr" {
                return Err(SsrMergeError::NotSsr {
                    path,
                    kind: header.kind.clone(),
                });
            }
            match header.writer.parameters.get(CATALOG_REFERENCE_MD5_PARAM) {
                Some(ParameterValue::String(found)) if *found == expected_md5 => {}
                Some(ParameterValue::String(found)) => {
                    return Err(SsrMergeError::CatalogMismatch {
                        path,
                        expected: expected_md5,
                        found: found.clone(),
                    });
                }
                // REVIEW ON UPGRADE: absent, or present-but-not-a-String (a future
                // `ParameterValue` variant), is treated as "missing the binding".
                _ => return Err(SsrMergeError::MissingCatalogMd5 { path }),
            }

            let mut remap = Vec::with_capacity(header.chromosomes.len());
            for chrom in &header.chromosomes {
                match name_to_global.get(&chrom.name) {
                    Some(&id) => remap.push(id),
                    None => {
                        return Err(SsrMergeError::UnknownInputChrom {
                            path,
                            chrom: chrom.name.clone(),
                        });
                    }
                }
            }

            let cursor =
                SampleEvidenceCursor::new(reader, remap.into_boxed_slice()).map_err(|source| {
                    SsrMergeError::Read {
                        sample: path.clone(),
                        source,
                    }
                })?;
            cursors.push(cursor);
            labels.push(path);
        }

        Ok(Self {
            catalog,
            cursors,
            labels,
            name_to_global,
            prev_locus: None,
            next_seq: 0,
        })
    }

    /// Produce the next emitted `CohortLocus`, skipping catalog loci no sample covers.
    /// `Ok(None)` at the end of the catalog.
    fn next_locus(&mut self) -> Result<Option<(u64, CohortLocus)>, SsrMergeError> {
        loop {
            let locus = match self.catalog.read_locus() {
                None => return Ok(None),
                Some(Ok(locus)) => locus,
                Some(Err(err)) => return Err(err.into()),
            };

            let chrom_id = *self.name_to_global.get(locus.chrom()).ok_or_else(|| {
                SsrMergeError::UnknownCatalogChrom {
                    chrom: locus.chrom().to_string(),
                }
            })?;
            // The catalog is already in the `LocusId` frame (0-based half-open).
            let locus_id = LocusId {
                chrom_id,
                start: locus.start(),
                end: locus.end(),
            };

            // Enforce the strictly-ascending catalog contract before asking any cursor
            // (an out-of-order query would otherwise trip the cursor's monotonic guard).
            if let Some(prev) = self.prev_locus
                && locus_id <= prev
            {
                return Err(SsrMergeError::UnsortedCatalog { prev, at: locus_id });
            }
            self.prev_locus = Some(locus_id);

            let mut cohort =
                CohortLocus::new(locus_id, locus.motif(), Box::from(locus.ref_bytes()));
            let labels = &self.labels;
            for (sample_idx, cursor) in self.cursors.iter_mut().enumerate() {
                let evidence =
                    cursor
                        .evidence_at(locus_id)
                        .map_err(|source| SsrMergeError::Read {
                            sample: labels[sample_idx].clone(),
                            source,
                        })?;
                if let Some(evidence) = evidence {
                    cohort.push(sample_idx as u32, evidence);
                }
            }

            if cohort.is_empty() {
                continue; // sparse-omit: no sample covered this locus
            }
            let seq = self.next_seq;
            self.next_seq += 1;
            return Ok(Some((seq, cohort)));
        }
    }

    /// Cohort-global chromosome names, indexed by global id (the inverse of
    /// `name_to_global`) — for labelling output.
    ///
    /// Relies on `name_to_global`'s values being a dense `0..len` permutation, which
    /// holds because `from_parts` assigns them via `enumerate()` over the chromosome
    /// table; a different id assignment would leave gaps (empty names) here.
    pub(crate) fn chrom_names(&self) -> Vec<String> {
        let mut names = vec![String::new(); self.name_to_global.len()];
        for (name, &id) in &self.name_to_global {
            names[id as usize] = name.clone();
        }
        names
    }
}

impl CohortMerger<BufReader<File>, BufReader<File>> {
    /// Open the catalog and the per-sample `.ssr.psp` files from disk and build the
    /// validated merger ([`from_parts`](Self::from_parts)).
    pub(crate) fn open(catalog: &Path, psp_files: &[PathBuf]) -> Result<Self, SsrMergeError> {
        let catalog_file = File::open(catalog).map_err(|source| SsrMergeError::Open {
            path: catalog.display().to_string(),
            source,
        })?;
        let catalog_reader = CatalogReader::new(BufReader::new(catalog_file))?;

        let mut inputs = Vec::with_capacity(psp_files.len());
        for path in psp_files {
            let file = File::open(path).map_err(|source| SsrMergeError::Open {
                path: path.display().to_string(),
                source,
            })?;
            let reader =
                PspReader::new(BufReader::new(file)).map_err(|source| SsrMergeError::OpenPsp {
                    path: path.display().to_string(),
                    source,
                })?;
            inputs.push((path.display().to_string(), reader));
        }

        Self::from_parts(catalog_reader, inputs)
    }
}

impl<R: Read + Seek, C: Read> Iterator for CohortMerger<R, C> {
    type Item = Result<(u64, CohortLocus), SsrMergeError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_locus().transpose()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ssr::cohort::test_support::{
        REF_MD5, catalog_reader, loc, obs, reader, rec, ssr_header, ssr_header_without_md5, ssr_psp,
    };
    use std::io::Cursor;

    /// Collect a merge to `(seq, LocusId, present indices, per-sample first-seq bytes)`.
    fn run(
        merger: CohortMerger<Cursor<Vec<u8>>, Cursor<Vec<u8>>>,
    ) -> Vec<(u64, LocusId, Vec<u32>)> {
        merger
            .map(|r| {
                let (seq, cl) = r.unwrap();
                (seq, cl.locus, cl.present.clone())
            })
            .collect()
    }

    #[test]
    fn merges_present_samples_in_catalog_order_with_sparse_omit() {
        // Catalog: four loci on chr1; the locus at 200 is covered by nobody.
        let loci = [
            loc("chr1", 16),
            loc("chr1", 60),
            loc("chr1", 100),
            loc("chr1", 200),
        ];
        let catalog = catalog_reader(REF_MD5, &loci);

        // Sample A covers 16 and 100; sample B covers 60 and 100.
        let a = reader(ssr_psp(
            ssr_header(&["chr1"], REF_MD5),
            &[
                rec(0, 16, obs(&[(b"CACACA", 5)])),
                rec(0, 100, obs(&[(b"CA", 2)])),
            ],
        ));
        let b = reader(ssr_psp(
            ssr_header(&["chr1"], REF_MD5),
            &[
                rec(0, 60, obs(&[(b"CACACACA", 4)])),
                rec(0, 100, obs(&[(b"CA", 3)])),
            ],
        ));

        let merger =
            CohortMerger::from_parts(catalog, vec![("A".into(), a), ("B".into(), b)]).unwrap();
        let got = run(merger);

        // Locus 200 (all-absent) is dropped; seqs are dense over emitted loci.
        assert_eq!(
            got,
            vec![
                (
                    0,
                    LocusId {
                        chrom_id: 0,
                        start: 16,
                        end: 22
                    },
                    vec![0]
                ),
                (
                    1,
                    LocusId {
                        chrom_id: 0,
                        start: 60,
                        end: 66
                    },
                    vec![1]
                ),
                (
                    2,
                    LocusId {
                        chrom_id: 0,
                        start: 100,
                        end: 106
                    },
                    vec![0, 1]
                ),
            ]
        );
    }

    #[test]
    fn carries_the_catalog_frame_and_evidence() {
        let loci = [loc("chr1", 16)];
        let catalog = catalog_reader(REF_MD5, &loci);
        let a = reader(ssr_psp(
            ssr_header(&["chr1"], REF_MD5),
            &[rec(0, 16, obs(&[(b"CACACA", 5), (b"CACACACA", 1)]))],
        ));

        let mut merger = CohortMerger::from_parts(catalog, vec![("A".into(), a)]).unwrap();
        let (seq, cl) = merger.next().unwrap().unwrap();
        assert_eq!(seq, 0);
        assert_eq!(cl.motif.as_bytes(), b"CA");
        assert_eq!(cl.ref_frame.as_ref(), b"GGGGGGCACACATTTTTT");
        assert_eq!(cl.present, vec![0]);
        assert_eq!(
            cl.samples[0].seq_counts,
            obs(&[(b"CACACA", 5), (b"CACACACA", 1)])
        );
        assert_eq!(cl.samples[0].qc.depth, 30);
        assert!(merger.next().is_none());
    }

    #[test]
    fn orders_across_chromosomes() {
        let loci = [loc("chr1", 16), loc("chr2", 30)];
        let catalog = catalog_reader(REF_MD5, &loci);
        let a = reader(ssr_psp(
            ssr_header(&["chr1", "chr2"], REF_MD5),
            &[
                rec(0, 16, obs(&[(b"CACACA", 5)])),
                rec(1, 30, obs(&[(b"CACACA", 7)])),
            ],
        ));
        let got = run(CohortMerger::from_parts(catalog, vec![("A".into(), a)]).unwrap());
        assert_eq!(
            got,
            vec![
                (
                    0,
                    LocusId {
                        chrom_id: 0,
                        start: 16,
                        end: 22
                    },
                    vec![0]
                ),
                (
                    1,
                    LocusId {
                        chrom_id: 1,
                        start: 30,
                        end: 36
                    },
                    vec![0]
                ),
            ]
        );
    }

    #[test]
    fn rejects_no_inputs() {
        let catalog = catalog_reader(REF_MD5, &[loc("chr1", 16)]);
        let result = CohortMerger::<Cursor<Vec<u8>>, _>::from_parts(catalog, vec![]);
        assert!(matches!(result, Err(SsrMergeError::NoInputs)));
    }

    #[test]
    fn rejects_a_different_catalog_md5() {
        let catalog = catalog_reader(REF_MD5, &[loc("chr1", 16)]);
        let wrong = reader(ssr_psp(
            ssr_header(&["chr1"], &"b".repeat(32)),
            &[rec(0, 16, obs(&[(b"CACACA", 5)]))],
        ));
        let result = CohortMerger::from_parts(catalog, vec![("wrong".into(), wrong)]);
        match result {
            Err(SsrMergeError::CatalogMismatch { found, .. }) => assert_eq!(found, "b".repeat(32)),
            _ => panic!("expected CatalogMismatch"),
        }
    }

    #[test]
    fn rejects_a_missing_catalog_md5_param() {
        let catalog = catalog_reader(REF_MD5, &[loc("chr1", 16)]);
        let no_param = reader(ssr_psp(
            ssr_header_without_md5(&["chr1"]),
            &[rec(0, 16, obs(&[(b"CACACA", 5)]))],
        ));
        let result = CohortMerger::from_parts(catalog, vec![("p".into(), no_param)]);
        assert!(matches!(
            result,
            Err(SsrMergeError::MissingCatalogMd5 { .. })
        ));
    }

    #[test]
    fn errors_on_a_catalog_chromosome_no_input_declares() {
        // Catalog references chr2, but the only input declares chr1 only.
        let loci = [loc("chr2", 16)];
        let catalog = catalog_reader(REF_MD5, &loci);
        let a = reader(ssr_psp(ssr_header(&["chr1"], REF_MD5), &[]));
        let mut merger = CohortMerger::from_parts(catalog, vec![("A".into(), a)]).unwrap();
        let err = merger.next().unwrap().unwrap_err();
        assert!(
            matches!(err, SsrMergeError::UnknownCatalogChrom { ref chrom } if chrom == "chr2"),
            "got {err:?}"
        );
    }

    #[test]
    fn errors_on_an_unsorted_catalog() {
        // Loci written out of coordinate order (the catalog writer does not sort).
        let loci = [loc("chr1", 60), loc("chr1", 16)];
        let catalog = catalog_reader(REF_MD5, &loci);
        let a = reader(ssr_psp(ssr_header(&["chr1"], REF_MD5), &[]));
        let mut merger = CohortMerger::from_parts(catalog, vec![("A".into(), a)]).unwrap();
        // The second locus regresses → hard error (before any cursor is asked).
        let err = merger.next().unwrap().unwrap_err();
        assert!(
            matches!(err, SsrMergeError::UnsortedCatalog { .. }),
            "got {err:?}"
        );
    }
}
