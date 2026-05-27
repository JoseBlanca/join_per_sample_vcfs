//! Rewrite a directory of `.psp` files with a custom
//! `TARGET_BLOCK_BYTES`. Used to sweep block-size settings and measure
//! the (compressed-size, per-reader-memory, wall-time) trade-off
//! curve.
//!
//! Reads every `*.psp` in `--input-dir`, decodes it once into memory,
//! and writes a sibling `*.psp` under `--output-dir` using
//! `PspWriter::new_with_block_target(target_bytes)` instead of the
//! default 16 MiB. Files that already exist at the output path are
//! skipped (lets reruns of the sweep at the same target be fast).
//!
//! ```text
//! ./scripts/dev.sh ./target-container/release/examples/blocksize_rewrite \
//!     --input-dir benchmarks/tomato1/results/ours/cohort/psp \
//!     --output-dir tmp/blocksize_sweep/1mb \
//!     --target-bytes 1048576
//! ```

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::time::Instant;

use clap::Parser;

use pop_var_caller::pileup_record::PileupRecord;
use pop_var_caller::psp::PspReader;
use pop_var_caller::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use pop_var_caller::psp::writer::PspWriter;

#[derive(Parser, Debug)]
#[command(about = "Rewrite a directory of .psp files at a custom block-bytes target")]
struct Cli {
    /// Source directory containing `*.psp` files to rewrite.
    #[arg(long)]
    input_dir: PathBuf,

    /// Destination directory. Created if missing. Existing files at the
    /// target path are skipped.
    #[arg(long)]
    output_dir: PathBuf,

    /// Target uncompressed bytes per block (the constant the writer
    /// flushes against). The writer's default is 16 MiB = 16777216.
    /// Use this flag to sweep.
    #[arg(long)]
    target_bytes: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    std::fs::create_dir_all(&args.output_dir)?;

    let mut psps: Vec<PathBuf> = std::fs::read_dir(&args.input_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("psp"))
        .collect();
    psps.sort();

    eprintln!(
        "rewriting {} *.psp from {} → {} at target_bytes={}",
        psps.len(),
        args.input_dir.display(),
        args.output_dir.display(),
        args.target_bytes,
    );

    let sweep_start = Instant::now();
    for src in &psps {
        let file_name = src.file_name().unwrap();
        let dst = args.output_dir.join(file_name);
        if dst.exists() {
            eprintln!("  skip (exists): {}", dst.display());
            continue;
        }
        let t0 = Instant::now();
        let (records, header) = read_psp(src)?;
        let read_elapsed = t0.elapsed();
        let t1 = Instant::now();
        write_psp(&dst, header, &records, args.target_bytes)?;
        let write_elapsed = t1.elapsed();
        let in_bytes = std::fs::metadata(src)?.len();
        let out_bytes = std::fs::metadata(&dst)?.len();
        let ratio = out_bytes as f64 / in_bytes as f64;
        eprintln!(
            "  {} {} recs  in={:>10} bytes  out={:>10} bytes  (×{:.2})  read={:.1?}  write={:.1?}",
            file_name.to_string_lossy(),
            records.len(),
            in_bytes,
            out_bytes,
            ratio,
            read_elapsed,
            write_elapsed,
        );
    }
    eprintln!("sweep total: {:.1?}", sweep_start.elapsed());
    Ok(())
}

fn read_psp(path: &Path) -> Result<(Vec<PileupRecord>, WriterHeader), Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let mut reader = PspReader::new(BufReader::with_capacity(64 * 1024, file))?;
    let parsed = reader.header().clone();
    let writer_header = WriterHeader {
        format_version: parsed.format_version,
        sample: parsed.sample.clone(),
        reference: parsed.reference.clone(),
        created: parsed.created,
        chromosomes: parsed
            .chromosomes
            .iter()
            .map(|c| ChromosomeEntry {
                name: c.name.clone(),
                length: c.length,
                md5: c.md5.clone(),
            })
            .collect(),
        writer: WriterProvenance {
            tool: parsed.writer.tool.clone(),
            version: parsed.writer.version.clone(),
            subcommand: parsed.writer.subcommand.clone(),
            input_crams: parsed.writer.input_crams.clone(),
            input_fasta: parsed.writer.input_fasta.clone(),
            parameters: clone_params(&parsed.writer.parameters),
        },
    };
    let mut records = Vec::with_capacity(1024 * 1024);
    for item in reader.records() {
        records.push(item?);
    }
    Ok((records, writer_header))
}

fn write_psp(
    path: &Path,
    header: WriterHeader,
    records: &[PileupRecord],
    target_bytes: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(path)?;
    let sink = BufWriter::with_capacity(64 * 1024, file);
    let mut writer = PspWriter::new_with_block_target(sink, header, target_bytes)?;
    for record in records {
        writer.write_record(record)?;
    }
    let buf = writer.finish()?;
    buf.into_inner()?.sync_all()?;
    Ok(())
}

fn clone_params(src: &BTreeMap<String, ParameterValue>) -> BTreeMap<String, ParameterValue> {
    src.iter().map(|(k, v)| (k.clone(), v.clone())).collect()
}
