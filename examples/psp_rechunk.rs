//! Re-chunk existing `.psp` files with a new block-cut window, without
//! re-running pileup (scaling work, 2026-06-02). Reads every record
//! from each input and rewrites it through a `PspWriter` configured
//! with `--block-window-bp = <window>` — identical records, only the
//! on-disk block boundaries change. Lets us A/B the genomic-window
//! block alignment against the existing (byte-cut, misaligned) cohort.
//!
//! Run:
//!   cargo build --release --example psp_rechunk
//!   ./target/release/examples/psp_rechunk <window_bp> <out_dir> <in1.psp> [in2.psp ...]

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};

use pop_var_caller::psp::PspReader;
use pop_var_caller::psp::header::{
    ChromosomeEntry, ParameterValue, WriterHeader, WriterProvenance,
};
use pop_var_caller::psp::writer::{MAX_BLOCK_TARGET_BYTES, PspWriter};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let window: u32 = args
        .next()
        .ok_or("usage: psp_rechunk <window_bp> <out_dir> <in.psp ...>")?
        .parse()?;
    let out_dir = PathBuf::from(args.next().ok_or("missing out_dir")?);
    let inputs: Vec<PathBuf> = args.map(PathBuf::from).collect();
    if inputs.is_empty() {
        return Err("no input .psp files".into());
    }
    std::fs::create_dir_all(&out_dir)?;

    for input in &inputs {
        let out = out_dir.join(input.file_name().unwrap());
        rechunk(input, &out, window)?;
        println!("rechunked {} -> {}", input.display(), out.display());
    }
    Ok(())
}

fn rechunk(input: &Path, output: &Path, window: u32) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::open(input)?;
    let mut reader = PspReader::new(BufReader::with_capacity(64 * 1024, file))?;
    let parsed = reader.header().clone();
    let header = WriterHeader {
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
            command_line: parsed.writer.command_line.clone(),
            parameters: parsed
                .writer
                .parameters
                .iter()
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect::<BTreeMap<String, ParameterValue>>(),
        },
    };

    let out_file = File::create(output)?;
    let sink = BufWriter::with_capacity(64 * 1024, out_file);
    // Window is the primary cut; cap generous so it never interferes.
    let mut writer =
        PspWriter::new_with_block_layout(sink, header, MAX_BLOCK_TARGET_BYTES, window)?;
    for item in reader.records() {
        writer.write_record(&item?)?;
    }
    let buf = writer.finish()?;
    buf.into_inner()?.sync_all()?;
    Ok(())
}
