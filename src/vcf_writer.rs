use std::fs::File;
use std::io::{self, BufWriter, Write};

use flate2::write::GzEncoder;
use flate2::Compression;

use crate::gvcf_parser::{VcfResult, Variant};

/// A writer that outputs merged variants in VCF format.
///
/// Supports writing to files (plain or gzipped) and to stdout.
pub struct VcfWriter {
    writer: Box<dyn Write>,
    ploidy: usize,
    broken_pipe: bool,
}

impl VcfWriter {
    /// Create a VCF writer that writes to the given file path.
    ///
    /// If the path ends with `.gz`, the output will be gzip-compressed.
    pub fn from_path(path: &str, samples: &[String]) -> io::Result<Self> {
        let writer: Box<dyn Write> = if path.ends_with(".gz") {
            let file = File::create(path)?;
            Box::new(BufWriter::new(GzEncoder::new(file, Compression::default())))
        } else {
            let file = File::create(path)?;
            Box::new(BufWriter::new(file))
        };

        let mut vcf_writer = VcfWriter { writer, ploidy: 2, broken_pipe: false };
        vcf_writer.write_header(samples)?;
        Ok(vcf_writer)
    }

    /// Create a VCF writer that writes to stdout.
    pub fn from_stdout(samples: &[String]) -> io::Result<Self> {
        let writer: Box<dyn Write> = Box::new(BufWriter::new(io::stdout().lock()));
        let mut vcf_writer = VcfWriter { writer, ploidy: 2, broken_pipe: false };
        vcf_writer.write_header(samples)?;
        Ok(vcf_writer)
    }

    /// Create a VCF writer from any `Write` implementor (useful for testing).
    pub fn from_writer(writer: Box<dyn Write>, samples: &[String]) -> io::Result<Self> {
        let mut vcf_writer = VcfWriter { writer, ploidy: 2, broken_pipe: false };
        vcf_writer.write_header(samples)?;
        Ok(vcf_writer)
    }

    /// Write to the underlying writer, silently stopping on broken pipe.
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        if self.broken_pipe {
            return Ok(());
        }
        match self.writer.write_all(buf) {
            Err(e) if e.kind() == io::ErrorKind::BrokenPipe => {
                self.broken_pipe = true;
                Ok(())
            }
            other => other,
        }
    }

    fn write_header(&mut self, samples: &[String]) -> io::Result<()> {
        self.write_all(b"##fileformat=VCFv4.2\n")?;
        self.write_all(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
        for sample in samples {
            self.write_all(format!("\t{}", sample).as_bytes())?;
        }
        self.write_all(b"\n")?;
        Ok(())
    }

    /// Write a single variant record.
    pub fn write_variant(&mut self, variant: &Variant) -> io::Result<()> {
        if self.broken_pipe {
            return Ok(());
        }

        let ref_allele = &variant.alleles[0];
        let alt_alleles = if variant.alleles.len() > 1 {
            variant.alleles[1..].join(",")
        } else {
            ".".to_string()
        };

        self.write_all(
            format!(
                "{}\t{}\t.\t{}\t{}\t.\t.\t.\tGT",
                variant.chrom, variant.pos, ref_allele, alt_alleles
            )
            .as_bytes(),
        )?;

        for sample_idx in 0..variant.n_samples {
            let gt_start = sample_idx * self.ploidy;
            let gt = &variant.genotypes[gt_start..gt_start + self.ploidy];
            let phase_sep = if variant.phase[sample_idx] { "|" } else { "/" };

            self.write_all(b"\t")?;
            for (i, &a) in gt.iter().enumerate() {
                if i > 0 {
                    self.write_all(phase_sep.as_bytes())?;
                }
                if a < 0 {
                    self.write_all(b".")?;
                } else {
                    self.write_all(format!("{}", a).as_bytes())?;
                }
            }
        }
        self.write_all(b"\n")?;
        Ok(())
    }

    /// Write all variants from an iterator.
    pub fn write_variants<I>(&mut self, variants: I) -> VcfResult<()>
    where
        I: Iterator<Item = VcfResult<Variant>>,
    {
        for result in variants {
            if self.broken_pipe {
                return Ok(());
            }
            let variant = result?;
            self.write_variant(&variant).map_err(crate::errors::VcfParseError::from)?;
        }
        Ok(())
    }

    /// Flush the underlying writer. Silently ignores broken pipe.
    pub fn flush(&mut self) -> io::Result<()> {
        if self.broken_pipe {
            return Ok(());
        }
        match self.writer.flush() {
            Err(e) if e.kind() == io::ErrorKind::BrokenPipe => {
                self.broken_pipe = true;
                Ok(())
            }
            other => other,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gvcf_parser::Variant;
    use flate2::read::MultiGzDecoder;
    use std::io::Read;
    use tempfile::NamedTempFile;

    fn make_variant(
        chrom: &str,
        pos: u32,
        alleles: &[&str],
        genotypes: Vec<i8>,
        phase: Vec<bool>,
        n_samples: usize,
    ) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            pos,
            alleles: alleles.iter().map(|s| s.to_string()).collect(),
            ref_allele_len: alleles[0].len() as u8,
            qual: f32::NAN,
            genotypes,
            phase,
            gt_format_fields: Vec::new(),
            sample_gt_fields: Vec::new(),
            n_samples,
        }
    }

    #[test]
    fn test_write_single_variant() {
        let samples = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
        let buf: Vec<u8> = Vec::new();
        let mut writer =
            VcfWriter::from_writer(Box::new(buf), &samples).unwrap();

        let variant = make_variant(
            "chr1",
            100,
            &["A", "T"],
            vec![0, 1, 1, 1],
            vec![false, true],
            2,
        );
        writer.write_variant(&variant).unwrap();
        writer.flush().unwrap();

        // Recover the buffer from the writer
        // We need a different approach - use a shared buffer
    }

    #[test]
    fn test_write_to_buffer() {
        let samples = vec!["SA".to_string(), "SB".to_string()];
        let buf = Vec::new();
        let cursor = io::Cursor::new(buf);
        let mut writer = VcfWriter {
            writer: Box::new(cursor),
            ploidy: 2,
            broken_pipe: false,
        };
        writer.write_header(&samples).unwrap();

        let variant = make_variant(
            "chr1",
            100,
            &["A", "T"],
            vec![0, 1, 1, 1],
            vec![false, true],
            2,
        );
        writer.write_variant(&variant).unwrap();
        writer.flush().unwrap();
    }

    /// Helper: write variants to an in-memory buffer and return the output as a string.
    fn write_to_string(samples: &[String], variants: &[Variant]) -> String {
        let buf: Vec<u8> = Vec::new();
        let shared_buf = std::sync::Arc::new(std::sync::Mutex::new(buf));
        let shared_clone = shared_buf.clone();

        // Use a wrapper that writes to the shared buffer
        struct SharedWriter(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl Write for SharedWriter {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let mut writer =
            VcfWriter::from_writer(Box::new(SharedWriter(shared_clone)), samples).unwrap();
        for v in variants {
            writer.write_variant(v).unwrap();
        }
        writer.flush().unwrap();
        drop(writer);

        let buf = shared_buf.lock().unwrap();
        String::from_utf8(buf.clone()).unwrap()
    }

    #[test]
    fn test_header_format() {
        let samples = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
        let output = write_to_string(&samples, &[]);
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert_eq!(
            lines[1],
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2"
        );
    }

    #[test]
    fn test_single_sample_het() {
        let samples = vec!["S1".to_string()];
        let variant = make_variant("chr1", 42, &["A", "G"], vec![0, 1], vec![false], 1);
        let output = write_to_string(&samples, &[variant]);
        let data_line = output.lines().nth(2).unwrap();
        assert_eq!(data_line, "chr1\t42\t.\tA\tG\t.\t.\t.\tGT\t0/1");
    }

    #[test]
    fn test_single_sample_phased_het() {
        let samples = vec!["S1".to_string()];
        let variant = make_variant("chr1", 42, &["A", "G"], vec![0, 1], vec![true], 1);
        let output = write_to_string(&samples, &[variant]);
        let data_line = output.lines().nth(2).unwrap();
        assert_eq!(data_line, "chr1\t42\t.\tA\tG\t.\t.\t.\tGT\t0|1");
    }

    #[test]
    fn test_hom_ref() {
        let samples = vec!["S1".to_string()];
        let variant = make_variant("chr1", 10, &["C", "T"], vec![0, 0], vec![false], 1);
        let output = write_to_string(&samples, &[variant]);
        let data_line = output.lines().nth(2).unwrap();
        assert_eq!(data_line, "chr1\t10\t.\tC\tT\t.\t.\t.\tGT\t0/0");
    }

    #[test]
    fn test_hom_alt() {
        let samples = vec!["S1".to_string()];
        let variant = make_variant("chr1", 10, &["C", "T"], vec![1, 1], vec![false], 1);
        let output = write_to_string(&samples, &[variant]);
        let data_line = output.lines().nth(2).unwrap();
        assert_eq!(data_line, "chr1\t10\t.\tC\tT\t.\t.\t.\tGT\t1/1");
    }

    #[test]
    fn test_missing_genotype() {
        let samples = vec!["S1".to_string()];
        let variant = make_variant("chr1", 10, &["A", "T"], vec![-1, -1], vec![false], 1);
        let output = write_to_string(&samples, &[variant]);
        let data_line = output.lines().nth(2).unwrap();
        assert_eq!(data_line, "chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t./.");
    }

    #[test]
    fn test_no_alt_allele() {
        let samples = vec!["S1".to_string()];
        let variant = make_variant("chr1", 10, &["A"], vec![0, 0], vec![false], 1);
        let output = write_to_string(&samples, &[variant]);
        let data_line = output.lines().nth(2).unwrap();
        assert_eq!(data_line, "chr1\t10\t.\tA\t.\t.\t.\t.\tGT\t0/0");
    }

    #[test]
    fn test_multi_sample() {
        let samples = vec!["S1".to_string(), "S2".to_string(), "S3".to_string()];
        let variant = make_variant(
            "chr2",
            500,
            &["G", "A", "C"],
            vec![0, 1, 1, 2, 0, 0],
            vec![false, true, false],
            3,
        );
        let output = write_to_string(&samples, &[variant]);
        let data_line = output.lines().nth(2).unwrap();
        assert_eq!(
            data_line,
            "chr2\t500\t.\tG\tA,C\t.\t.\t.\tGT\t0/1\t1|2\t0/0"
        );
    }

    #[test]
    fn test_multiple_variants() {
        let samples = vec!["S1".to_string()];
        let v1 = make_variant("chr1", 10, &["A", "T"], vec![0, 1], vec![false], 1);
        let v2 = make_variant("chr1", 20, &["G", "C"], vec![1, 1], vec![false], 1);
        let output = write_to_string(&samples, &[v1, v2]);
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 4); // header + column header + 2 data lines
        assert_eq!(lines[2], "chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/1");
        assert_eq!(lines[3], "chr1\t20\t.\tG\tC\t.\t.\t.\tGT\t1/1");
    }

    #[test]
    fn test_write_variants_from_iterator() {
        let samples = vec!["S1".to_string()];
        let v1 = make_variant("chr1", 10, &["A", "T"], vec![0, 1], vec![false], 1);
        let v2 = make_variant("chr1", 20, &["G", "C"], vec![1, 1], vec![false], 1);

        struct SharedWriter(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl Write for SharedWriter {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let shared_buf = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut writer =
            VcfWriter::from_writer(Box::new(SharedWriter(shared_buf.clone())), &samples).unwrap();

        let iter = vec![Ok(v1), Ok(v2)].into_iter();
        writer.write_variants(iter).unwrap();
        writer.flush().unwrap();
        drop(writer);

        let buf = shared_buf.lock().unwrap();
        let output = String::from_utf8(buf.clone()).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines.len(), 4);
    }

    #[test]
    fn test_write_to_plain_file() {
        let tmp = NamedTempFile::with_suffix(".vcf").unwrap();
        let path = tmp.path().to_str().unwrap().to_string();
        drop(tmp); // close so VcfWriter can create it

        let samples = vec!["S1".to_string()];
        {
            let mut writer = VcfWriter::from_path(&path, &samples).unwrap();
            let variant =
                make_variant("chr1", 100, &["A", "T"], vec![0, 1], vec![false], 1);
            writer.write_variant(&variant).unwrap();
            writer.flush().unwrap();
        }

        let contents = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = contents.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert_eq!(lines[1], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1");
        assert_eq!(lines[2], "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1");

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_write_to_gzipped_file() {
        let tmp = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        let path = tmp.path().to_str().unwrap().to_string();
        drop(tmp);

        let samples = vec!["S1".to_string()];
        {
            let mut writer = VcfWriter::from_path(&path, &samples).unwrap();
            let variant =
                make_variant("chr1", 200, &["G", "A"], vec![1, 1], vec![true], 1);
            writer.write_variant(&variant).unwrap();
            writer.flush().unwrap();
        }

        // Read back the gzipped file and verify
        let file = File::open(&path).unwrap();
        let mut decoder = MultiGzDecoder::new(file);
        let mut contents = String::new();
        decoder.read_to_string(&mut contents).unwrap();

        let lines: Vec<&str> = contents.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert_eq!(lines[1], "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1");
        assert_eq!(lines[2], "chr1\t200\t.\tG\tA\t.\t.\t.\tGT\t1|1");

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn test_write_variants_error_propagation() {
        let samples = vec!["S1".to_string()];

        struct SharedWriter(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl Write for SharedWriter {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let shared_buf = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut writer =
            VcfWriter::from_writer(Box::new(SharedWriter(shared_buf.clone())), &samples).unwrap();

        let iter = vec![Err(crate::errors::VcfParseError::NotVariable)].into_iter();
        let result = writer.write_variants(iter);
        assert!(result.is_err());
    }

    #[test]
    fn test_partial_missing_genotype() {
        let samples = vec!["S1".to_string()];
        let variant = make_variant("chr1", 10, &["A", "T"], vec![0, -1], vec![false], 1);
        let output = write_to_string(&samples, &[variant]);
        let data_line = output.lines().nth(2).unwrap();
        assert_eq!(data_line, "chr1\t10\t.\tA\tT\t.\t.\t.\tGT\t0/.");
    }

    /// A writer that returns BrokenPipe after a configured number of bytes.
    struct BrokenPipeWriter {
        bytes_until_break: usize,
        bytes_written: usize,
    }

    impl BrokenPipeWriter {
        fn new(bytes_until_break: usize) -> Self {
            Self {
                bytes_until_break,
                bytes_written: 0,
            }
        }
    }

    impl Write for BrokenPipeWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.bytes_written >= self.bytes_until_break {
                return Err(io::Error::new(io::ErrorKind::BrokenPipe, "pipe closed"));
            }
            let allowed = (self.bytes_until_break - self.bytes_written).min(buf.len());
            self.bytes_written += allowed;
            Ok(allowed)
        }
        fn flush(&mut self) -> io::Result<()> {
            if self.bytes_written >= self.bytes_until_break {
                return Err(io::Error::new(io::ErrorKind::BrokenPipe, "pipe closed"));
            }
            Ok(())
        }
    }

    #[test]
    fn test_broken_pipe_on_write_variant_does_not_error() {
        let samples = vec!["S1".to_string()];
        // Allow enough bytes for the header, then break on the first variant
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
        let mut writer = VcfWriter::from_writer(
            Box::new(BrokenPipeWriter::new(header.len())),
            &samples,
        )
        .unwrap();

        let v1 = make_variant("chr1", 10, &["A", "T"], vec![0, 1], vec![false], 1);
        let v2 = make_variant("chr1", 20, &["G", "C"], vec![1, 1], vec![false], 1);

        // First write triggers the broken pipe — should return Ok
        assert!(writer.write_variant(&v1).is_ok());
        // Subsequent writes are silent no-ops
        assert!(writer.write_variant(&v2).is_ok());
        // Flush is also a no-op
        assert!(writer.flush().is_ok());
    }

    #[test]
    fn test_broken_pipe_on_write_variants_stops_silently() {
        let samples = vec!["S1".to_string()];
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
        let mut writer = VcfWriter::from_writer(
            Box::new(BrokenPipeWriter::new(header.len())),
            &samples,
        )
        .unwrap();

        let v1 = make_variant("chr1", 10, &["A", "T"], vec![0, 1], vec![false], 1);
        let v2 = make_variant("chr1", 20, &["G", "C"], vec![1, 1], vec![false], 1);
        let v3 = make_variant("chr1", 30, &["C", "A"], vec![0, 0], vec![false], 1);

        let iter = vec![Ok(v1), Ok(v2), Ok(v3)].into_iter();
        // The whole iterator should complete without error
        assert!(writer.write_variants(iter).is_ok());
        assert!(writer.flush().is_ok());
    }

    #[test]
    fn test_broken_pipe_mid_variant_stops_silently() {
        let samples = vec!["S1".to_string()];
        // Break after header + a few bytes into the first variant
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
        let mut writer = VcfWriter::from_writer(
            Box::new(BrokenPipeWriter::new(header.len() + 5)),
            &samples,
        )
        .unwrap();

        let v1 = make_variant("chr1", 10, &["A", "T"], vec![0, 1], vec![false], 1);
        // Pipe breaks mid-write — should still return Ok
        assert!(writer.write_variant(&v1).is_ok());
    }

    #[test]
    fn test_broken_pipe_immediately_on_flush() {
        let samples = vec!["S1".to_string()];
        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
        let mut writer = VcfWriter::from_writer(
            Box::new(BrokenPipeWriter::new(header.len())),
            &samples,
        )
        .unwrap();

        // No variants written, but flush hits broken pipe — should be Ok
        assert!(writer.flush().is_ok());
    }
}
