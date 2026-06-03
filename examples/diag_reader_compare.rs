//! Diagnostic: compare the post-filter read sequence produced by the
//! streaming `AlignmentMergedReader::new` against the indexed
//! `AlignmentMergedReader::query` for a single contig.
//!
//! The whole-genome pileup byte-identity investigation established that
//! both paths admit the *same number* of reads but the walker emits
//! different records. This isolates the two halves: if the read
//! sequences printed here differ (count, order, or content), the
//! divergence is in the reader; if they're identical, it's the walker.
//!
//! Run on the host (reference lives outside the container):
//! ```text
//! cargo build --release --example diag_reader_compare
//! ./target/release/examples/diag_reader_compare \
//!     --reference <ref.fna> --contig chr1 <sample.cram>
//! ```

use std::path::PathBuf;

use clap::Parser;

use pop_var_caller::bam::alignment_input::{
    AlignmentMergedReader, AlignmentMergedReaderConfig, MappedRead, load_pileup_inputs,
};

#[derive(Parser, Debug)]
struct Cli {
    #[arg(long)]
    reference: PathBuf,
    /// Contig name to compare (e.g. `chr1`).
    #[arg(long)]
    contig: String,
    /// One alignment file (CRAM/BAM) with a sibling index.
    alignment_file: PathBuf,
}

/// A compact, order-sensitive fingerprint of one read.
fn fingerprint(r: &MappedRead) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{:?}\t{:?}",
        String::from_utf8_lossy(&r.qname),
        r.flag,
        r.pos,
        r.mapq,
        r.cigar,
        r.adaptor_boundary,
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let files = vec![cli.alignment_file.clone()];
    let cfg = AlignmentMergedReaderConfig::default();

    let inputs = load_pileup_inputs(&files, &cli.reference, false)?;
    let target_id = inputs
        .contigs
        .entries
        .iter()
        .position(|c| c.name == cli.contig)
        .ok_or("contig not in list")?;
    eprintln!("contig {} = ref_id {}", cli.contig, target_id);

    // Streaming: read until we've passed the target contig (reads are
    // sorted by (ref_id, pos), so once ref_id exceeds the target there
    // are no more target reads).
    let streaming: Vec<String> = {
        let reader = AlignmentMergedReader::new(&files, &cli.reference, cfg)?;
        let mut out = Vec::new();
        let mut seen = false;
        for r in reader {
            let r = r?;
            if r.ref_id == target_id {
                seen = true;
                out.push(fingerprint(&r));
            } else if seen {
                break;
            }
        }
        out
    };

    // Indexed: query the target contig directly.
    let indexed: Vec<String> = {
        let reader = AlignmentMergedReader::query(
            &files,
            &cli.reference,
            inputs.contigs.clone(),
            inputs.sample_name.clone(),
            &inputs.headers,
            &inputs.indexes,
            &cli.contig,
            None,
            cfg,
        )?;
        reader
            .map(|r| r.map(|x| fingerprint(&x)))
            .collect::<Result<_, _>>()?
    };

    eprintln!(
        "streaming reads: {}, indexed reads: {}",
        streaming.len(),
        indexed.len()
    );

    let mut diffs = 0usize;
    let n = streaming.len().max(indexed.len());
    for i in 0..n {
        let a = streaming.get(i);
        let b = indexed.get(i);
        if a != b {
            if diffs < 20 {
                eprintln!("DIFF @ index {i}:");
                eprintln!("  streaming: {}", a.map(String::as_str).unwrap_or("<none>"));
                eprintln!("  indexed  : {}", b.map(String::as_str).unwrap_or("<none>"));
            }
            diffs += 1;
        }
    }
    eprintln!("total positions compared: {n}, differing: {diffs}");
    if diffs == 0 {
        eprintln!("READ SEQUENCES IDENTICAL → divergence is in the walker, not the reader");
    } else {
        eprintln!("READ SEQUENCES DIFFER → divergence is in the reader");
    }
    Ok(())
}
