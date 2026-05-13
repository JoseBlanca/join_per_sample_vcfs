// Mi20: shared test fixtures for `psp/` test modules. Hoisted from
// the per-file duplicates in `header.rs::tests` and
// `writer.rs::tests` so a new mandatory `WriterHeader` field
// updates one place. Only compiled under `#[cfg(test)]`.
//
// Mi5 routes everything through `pub(crate)` so the consolidation
// stays inside the crate.

#![cfg(test)]

use std::collections::BTreeMap;

use super::header::{ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance};

/// One-chromosome `WriterHeader` fixture. The contig length is
/// roomy enough for all in-tree tests; callers that need a
/// different length, multiple chromosomes, or extra parameters
/// build on top of this.
pub(crate) fn writer_header(n_chroms: usize) -> WriterHeader {
    let chromosomes = (0..n_chroms)
        .map(|i| ChromosomeEntry {
            name: format!("chr{}", i + 1),
            length: 1_000_000,
            md5: "0".repeat(32),
        })
        .collect();
    let mut params = BTreeMap::new();
    params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
    WriterHeader {
        format_version: (1, 0),
        sample: "sample".to_string(),
        reference: "ref.fa".to_string(),
        created: "2026-05-13T10:00:00Z".parse().unwrap(),
        chromosomes,
        writer: WriterProvenance {
            tool: "test".to_string(),
            version: "0.0.1".to_string(),
            subcommand: "per-sample".to_string(),
            input_crams: vec!["a.cram".to_string()],
            input_fasta: "ref.fa".to_string(),
            parameters: params,
        },
    }
}

/// Realistic `WriterHeader` carrying actual sample / md5 values —
/// the `header.rs::tests` baseline.
pub(crate) fn realistic_writer_header() -> WriterHeader {
    let mut params = BTreeMap::new();
    params.insert("min-mapq".to_string(), ParameterValue::Integer(30));
    params.insert("drop-duplicate".to_string(), ParameterValue::Boolean(true));
    WriterHeader {
        format_version: (1, 0),
        sample: "NA12878".to_string(),
        reference: "GRCh38.fa".to_string(),
        created: "2026-05-13T10:00:00Z".parse().unwrap(),
        chromosomes: vec![ChromosomeEntry {
            name: "chr1".to_string(),
            length: 248956422,
            md5: "6aef897c3d6ff0c78aff06ac189178dd".to_string(),
        }],
        writer: WriterProvenance {
            tool: "join_per_sample_vcfs".to_string(),
            version: "0.3.0".to_string(),
            subcommand: "per-sample".to_string(),
            input_crams: vec!["sample.cram".to_string()],
            input_fasta: "GRCh38.fa".to_string(),
            parameters: params,
        },
    }
}
