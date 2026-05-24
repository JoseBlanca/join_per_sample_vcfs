//! Cross-subcommand helpers — small utilities that every
//! `pop_var_caller` subcommand needs (path basename, MD5 hex
//! formatting, current command line, RFC3339 timestamps, civil-date
//! conversion). Mi8 from the 2026-05-19 cohort CLI review
//! deduplicates these — they previously lived in 2-4 copies across
//! `cli.rs`, `var_calling.rs`, `var_calling_from_bam.rs`, and
//! `estimate_contamination.rs`.
//!
//! Also hosts the project-wide buffered-I/O capacity constant
//! ([`DEFAULT_BUFFERED_IO_CAPACITY`]) — Mi19 — used by every
//! `BufReader::with_capacity` / `BufWriter::with_capacity` call site
//! in the `pop_var_caller` and `vcf` modules.
//!
//! Everything here is `pub(crate)`; the helpers are not part of the
//! library's public surface (they would just be re-implemented by
//! any downstream consumer needing them).

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::time::{SystemTime, UNIX_EPOCH};

use md5::{Digest, Md5};
use noodles_fasta::fai;
use rayon::prelude::*;

use crate::psp::header::ParsedChromosome;

/// Window size used by [`compute_contig_md5_streaming`]. Bigger than
/// the page (4 KiB) so each `read` syscall amortises across many
/// kernel-side bytes; small enough that the buffer sits on the worker
/// stack instead of the heap.
const FASTA_MD5_BUFFER_SIZE: usize = 64 * 1024;

/// Buffered-I/O capacity used by every `BufReader::with_capacity` /
/// `BufWriter::with_capacity` in the cohort CLI's file I/O. 64 KiB
/// matches the value the per-sample `.psp` reader/writer documents
/// as the recommended minimum (see
/// [`crate::psp::reader::PspReader::new`]).
///
/// Distinct from `BLOCK_HEADER_READ_CAP` (in the psp module), which
/// sizes a different — and intentionally smaller — short-read
/// buffer for fixed-size block headers.
pub(crate) const DEFAULT_BUFFERED_IO_CAPACITY: usize = 64 * 1024;

/// Strip the directory portion from `p`, returning the file-name
/// portion as a lossy-UTF-8 `String`. Falls back to the full path
/// (also lossy) when the path has no file-name component (e.g. it
/// ends in `/` or `..`).
pub(crate) fn basename(p: &Path) -> String {
    p.file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.to_string_lossy().into_owned())
}

/// Format a 16-byte MD5 as 32 lowercase hex characters, matching
/// SAM `@SQ M5` and the `.psp` `chromosome.md5` field.
pub(crate) fn format_md5_hex(bytes: [u8; 16]) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(32);
    for b in bytes {
        // PANIC-FREE: `std::fmt::Write` on a `String` is infallible
        // — the only failure mode is OOM, which the runtime turns
        // into an abort, not an `Err`.
        write!(&mut out, "{b:02x}").expect("writing to a String never fails");
    }
    out
}

/// Reconstruct the invoking process's command line as a single
/// space-joined string suitable for `WriterProvenance.command_line`
/// / VCF `##commandline=`. Lossy-UTF-8 for non-Unicode argv entries.
pub(crate) fn current_command_line() -> String {
    std::env::args_os()
        .map(|a: OsString| a.to_string_lossy().into_owned())
        .collect::<Vec<_>>()
        .join(" ")
}

/// Format the current UTC time as a TOML-compatible RFC3339 string
/// (`YYYY-MM-DDTHH:MM:SSZ`). Done by hand so we don't pull in a date
/// crate just for this one call; uses Howard Hinnant's civil-date
/// algorithm ([`civil_from_days`]) to turn days-since-epoch into
/// (year, month, day).
pub(crate) fn rfc3339_now() -> String {
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let sod = secs % 86_400;
    let (y, m, d) = civil_from_days(days);
    let h = sod / 3600;
    let min = (sod % 3600) / 60;
    let s = sod % 60;
    format!("{y:04}-{m:02}-{d:02}T{h:02}:{min:02}:{s:02}Z")
}

/// Tracks whether [`configure_rayon_pool`] has already taken the
/// process-global rayon pool. Set once per process; never cleared.
static RAYON_POOL_CONFIGURED: OnceLock<()> = OnceLock::new();

/// Size rayon's process-global thread pool to `n` worker threads,
/// or leave rayon's default in place when `n` is `None`. Idempotent:
/// only the first call with `Some(_)` in a given process actually
/// touches the global pool; subsequent calls (with any value)
/// return `Ok(())` without touching anything. M13 from the
/// 2026-05-19 cohort CLI review — the policy is **silent no-op on
/// second call**, locked in plan-review on the same date.
///
/// Rayon's `ThreadPoolBuilder::build_global()` is itself a once-per-
/// process operation that errors on the second call. Centralising
/// the gate here means library consumers, test runners, and any
/// future multi-subcommand driver can call as many `run_*` helpers
/// as they like back-to-back without managing the first-call /
/// second-call discrimination themselves.
///
/// The first caller wins; subsequent callers' `n` is ignored. This
/// is the cheapest correct semantics — callers that care about a
/// specific thread count should set it via the binary entry-point
/// before anything else runs.
pub(crate) fn configure_rayon_pool(n: Option<usize>) -> Result<(), rayon::ThreadPoolBuildError> {
    let Some(n) = n else { return Ok(()) };
    if RAYON_POOL_CONFIGURED.get().is_some() {
        // Already configured by an earlier call in this process —
        // silent no-op rather than error.
        return Ok(());
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build_global()?;
    // `set` can race with another thread that won the
    // `build_global` race above, but our gate happened before
    // `build_global` succeeded, so at most one thread reaches this
    // line per process. Discard the result so a vanishingly-rare
    // double-set is silent.
    let _ = RAYON_POOL_CONFIGURED.set(());
    Ok(())
}

/// Why [`verify_fasta_matches_psp_chromosomes`] rejected the
/// supplied FASTA. **M5 follow-up** wraps two distinct failure
/// modes — actual MD5 mismatch (FASTA had the contig but the
/// bytes don't match) vs an I/O / lookup failure (contig name
/// missing from the FASTA index, truncated FASTA, etc.). Both
/// fail the cross-check but the user-facing remediation differs.
#[derive(Debug)]
pub(crate) enum FastaVerifyError {
    /// Contig is present in the FASTA but its uppercase-bases MD5
    /// does not match the value carried by the `.psp` header.
    Md5Mismatch {
        /// Contig name (shared between the `.psp` header and the
        /// FASTA — the upstream basename check guarantees they
        /// agree on naming).
        contig: String,
        /// Lowercase 32-hex MD5 of the contig's uppercase bases as
        /// read from the supplied FASTA.
        fasta_md5: String,
        /// 32-hex MD5 the `.psp` header carries (originally sourced
        /// from the producing CRAM's `@SQ M5`).
        psp_md5: String,
    },
    /// Reading the contig's bases from the FASTA failed — the
    /// contig name is missing from the `.fai`, the `.fai` itself is
    /// unreadable, the FASTA ended before the declared number of
    /// bases was consumed, or some other I/O error. Carries the
    /// contig name + the inner `io::Error` for context. For the
    /// `.fai` read failure, `contig` carries the literal placeholder
    /// `"<reading .fai>"`.
    FetchFailed { contig: String, source: io::Error },
}

/// Cross-check that each `.psp` header's per-contig MD5 matches the
/// MD5 of the corresponding uppercase contig bases in the supplied
/// FASTA. **M5 follow-up** from the 2026-05-19 cohort CLI review —
/// the basename cross-check (`run_var_calling` /
/// `run_estimate_contamination` step 3) is necessary but not
/// sufficient; this helper closes the "right basename, wrong
/// bytes" failure mode the project's "no silent intermediates"
/// principle exists to prevent.
///
/// Streaming + parallel (Phase A of the
/// [`reference_fasta_streaming`](../../../doc/devel/implementation_plans/reference_fasta_streaming.md)
/// plan): each contig's bases are streamed through `Md5::update` in
/// 64 KiB windows, so per-worker resident memory is one window — not
/// one contig. Per-contig work runs in parallel under rayon's global
/// pool; wall time becomes `max_chrom_size / min(threads, n_chroms)`,
/// not the serial sum. The previous "warms the `SyncRefFetcher`
/// cache as a side-effect" coupling is gone by design — this helper
/// never allocates a contig-sized buffer and never touches the
/// runtime fetcher's cache.
pub(crate) fn verify_fasta_matches_psp_chromosomes(
    fasta_path: &Path,
    chromosomes: &[ParsedChromosome],
) -> Result<(), FastaVerifyError> {
    // .fai sits at "<fasta>.fai" by convention. Read it once on the
    // caller thread; rayon workers below share it by reference.
    let mut fai_pathbuf = fasta_path.as_os_str().to_os_string();
    fai_pathbuf.push(".fai");
    let fai_path = PathBuf::from(fai_pathbuf);
    let index = fai::fs::read(&fai_path).map_err(|source| FastaVerifyError::FetchFailed {
        contig: String::from("<reading .fai>"),
        source,
    })?;

    chromosomes
        .par_iter()
        .try_for_each(|contig| -> Result<(), FastaVerifyError> {
            let record = index
                .as_ref()
                .iter()
                .find(|r| AsRef::<[u8]>::as_ref(r.name()) == contig.name.as_bytes())
                .ok_or_else(|| FastaVerifyError::FetchFailed {
                    contig: contig.name.clone(),
                    source: io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("contig {} not in FASTA index", contig.name),
                    ),
                })?;
            let fasta_md5 =
                compute_contig_md5_streaming(fasta_path, record.offset(), contig.length).map_err(
                    |source| FastaVerifyError::FetchFailed {
                        contig: contig.name.clone(),
                        source,
                    },
                )?;
            if fasta_md5 != contig.md5 {
                return Err(FastaVerifyError::Md5Mismatch {
                    contig: contig.name.clone(),
                    fasta_md5,
                    psp_md5: contig.md5.clone(),
                });
            }
            Ok(())
        })
}

/// Compute the MD5 of one contig's uppercase bases by streaming
/// fixed-size windows from disk. Matches `Md5::digest(uppercase_bytes)`
/// bit-for-bit; resident memory is one [`FASTA_MD5_BUFFER_SIZE`]
/// window regardless of contig size.
///
/// `offset` is the byte position of the contig's first base in the
/// FASTA (from the `.fai`). `expected_length` is the `.psp`-declared
/// number of bases; the function consumes exactly that many bases,
/// errors with [`io::ErrorKind::UnexpectedEof`] if the file ends
/// first, and ignores any trailing bytes if the FASTA contains more.
fn compute_contig_md5_streaming(
    fasta_path: &Path,
    offset: u64,
    expected_length: u32,
) -> io::Result<String> {
    let mut file = File::open(fasta_path)?;
    file.seek(SeekFrom::Start(offset))?;

    let mut md5 = Md5::new();
    let mut remaining: usize = expected_length as usize;
    // Read buffer holds raw FASTA bytes (sequence + newlines); upper
    // buffer holds the uppercased, newline-stripped slice fed to md5.
    // Both stack-allocated.
    let mut read_buf = [0u8; FASTA_MD5_BUFFER_SIZE];
    let mut upper_buf = [0u8; FASTA_MD5_BUFFER_SIZE];

    while remaining > 0 {
        let n = file.read(&mut read_buf)?;
        if n == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "FASTA ended with {remaining} base(s) still expected (offset {offset}, \
                     declared length {expected_length})"
                ),
            ));
        }
        let mut out_len = 0usize;
        for &b in &read_buf[..n] {
            if remaining == 0 {
                // All declared bases consumed; ignore trailing
                // bytes in this window (newlines, the next contig's
                // header, etc.).
                break;
            }
            if b == b'\n' || b == b'\r' {
                continue;
            }
            upper_buf[out_len] = b.to_ascii_uppercase();
            out_len += 1;
            remaining -= 1;
        }
        if out_len > 0 {
            md5.update(&upper_buf[..out_len]);
        }
    }

    let digest: [u8; 16] = md5.finalize().into();
    Ok(format_md5_hex(digest))
}

/// Howard Hinnant's `civil_from_days` — translate days-since-epoch
/// (1970-01-01 = 0) into a proleptic Gregorian (year, month, day).
/// Verified by external reference; see
/// <http://howardhinnant.github.io/date_algorithms.html>.
pub(crate) fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365; // [0, 399]
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32; // [1, 12]
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn basename_strips_directory() {
        assert_eq!(basename(Path::new("/tmp/a/b/sample.cram")), "sample.cram");
        assert_eq!(basename(Path::new("ref.fa")), "ref.fa");
    }

    #[test]
    fn basename_falls_back_to_full_path_on_trailing_slash() {
        // Path::file_name() returns None on a trailing slash; the
        // fallback path keeps the helper total.
        let p = PathBuf::from("/some/dir/");
        let s = basename(&p);
        // No assertion on exact value beyond "non-empty + UTF-8".
        // The contract is only that the helper does not panic.
        assert!(!s.is_empty());
    }

    #[test]
    fn format_md5_hex_is_32_lowercase() {
        let bytes = [
            0x6a, 0xef, 0x89, 0x7c, 0x3d, 0x6f, 0xf0, 0xc7, 0x8a, 0xff, 0x06, 0xac, 0x18, 0x91,
            0x78, 0xdd,
        ];
        assert_eq!(format_md5_hex(bytes), "6aef897c3d6ff0c78aff06ac189178dd");
    }

    #[test]
    fn civil_from_days_matches_known_dates() {
        // Unix epoch
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        // 1970-12-31 (day 364, 1970 is not a leap year)
        assert_eq!(civil_from_days(364), (1970, 12, 31));
        // 2000-01-01 (day 10957 since epoch, common Y2K test)
        assert_eq!(civil_from_days(10957), (2000, 1, 1));
        // 2024-02-29 — leap day
        assert_eq!(civil_from_days(19782), (2024, 2, 29));
    }

    #[test]
    fn rfc3339_now_parses_as_toml_datetime() {
        let s = rfc3339_now();
        let _: toml::value::Datetime = s.parse().expect("must parse as toml Datetime");
    }

    #[test]
    fn configure_rayon_pool_none_is_always_ok() {
        // `None` is the always-safe path: never touches rayon's
        // global pool and always returns `Ok(())`, regardless of
        // process state. M13's idempotency contract is partially
        // tested here (the `Some(_)` path is hard to exercise
        // deterministically under cargo-test's shared-process
        // model because other tests lazy-init rayon's pool; see
        // the wave's notes).
        configure_rayon_pool(None).unwrap();
        configure_rayon_pool(None).unwrap();
        configure_rayon_pool(None).unwrap();
    }

    // ---------------------------------------------------------------------
    // Phase A — streaming MD5 verify tests. The helper below writes a
    // FASTA + sibling `.fai` with caller-chosen line wrapping so the
    // newline-stripping logic in `compute_contig_md5_streaming` is
    // exercised against multi-line contigs (the production case —
    // Ensembl / Gencode FASTAs wrap at 60 chars).
    // ---------------------------------------------------------------------

    use std::fs::File;
    use std::io::Write;

    /// One contig in a synthetic FASTA fixture.
    struct ContigFixture<'a> {
        name: &'a str,
        /// Bases as a single contiguous string; the helper will line-
        /// wrap them to `line_bases` per line on disk.
        seq: &'a str,
        line_bases: usize,
    }

    /// Build a FASTA + `.fai` with the given contigs. Returns the
    /// tempdir (keep it alive for the duration of the test) and the
    /// FASTA path. Line-wraps each contig's bases per its
    /// `line_bases`. The last line may be short.
    fn build_wrapped_fasta(contigs: &[ContigFixture<'_>]) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("tempdir");
        let fasta_path = dir.path().join("ref.fa");
        let fai_path = dir.path().join("ref.fa.fai");

        let mut fa = File::create(&fasta_path).expect("fa");
        let mut fai = File::create(&fai_path).expect("fai");
        let mut offset: u64 = 0;

        for c in contigs {
            let header = format!(">{}\n", c.name);
            fa.write_all(header.as_bytes()).expect("hdr");
            offset += header.len() as u64;

            // Sequence start offset is the .fai's `offset` column.
            let seq_offset = offset;

            let bytes = c.seq.as_bytes();
            let mut written = 0usize;
            while written < bytes.len() {
                let end = (written + c.line_bases).min(bytes.len());
                fa.write_all(&bytes[written..end]).expect("seq line");
                fa.write_all(b"\n").expect("seq nl");
                offset += (end - written) as u64 + 1;
                written = end;
            }

            let length = c.seq.len();
            let line_width = c.line_bases + 1; // +1 for '\n'
            writeln!(
                fai,
                "{}\t{}\t{}\t{}\t{}",
                c.name, length, seq_offset, c.line_bases, line_width
            )
            .expect("fai entry");
        }

        (dir, fasta_path)
    }

    fn one_shot_md5(uppercase_bytes: &[u8]) -> String {
        let digest: [u8; 16] = Md5::digest(uppercase_bytes).into();
        format_md5_hex(digest)
    }

    #[test]
    fn streaming_md5_matches_one_shot_on_single_line_contig() {
        // Single-line contig — no newline-stripping needed inside the
        // contig body. Sanity check that the basic streaming pipeline
        // matches `Md5::digest`.
        let (_dir, fasta_path) = build_wrapped_fasta(&[ContigFixture {
            name: "chr0",
            seq: "ACGTACGTACGT",
            line_bases: 12, // single line
        }]);
        let index = fai::fs::read(fasta_path.with_extension("fa.fai")).expect("fai");
        let record = &index.as_ref()[0];

        let streaming = compute_contig_md5_streaming(&fasta_path, record.offset(), 12).unwrap();
        let one_shot = one_shot_md5(b"ACGTACGTACGT");
        assert_eq!(streaming, one_shot);
    }

    #[test]
    fn streaming_md5_matches_one_shot_on_line_wrapped_contig() {
        // 60-char-per-line wrap is the Ensembl/Gencode default; this
        // is the production case. The streaming function must skip
        // every embedded `\n` for the digest to match.
        let seq: String = (0..250).map(|i| b"ACGT"[i % 4] as char).collect();
        let (_dir, fasta_path) = build_wrapped_fasta(&[ContigFixture {
            name: "chr0",
            seq: &seq,
            line_bases: 60,
        }]);
        let index = fai::fs::read(fasta_path.with_extension("fa.fai")).expect("fai");
        let record = &index.as_ref()[0];

        let streaming = compute_contig_md5_streaming(&fasta_path, record.offset(), 250).unwrap();
        let one_shot = one_shot_md5(seq.as_bytes());
        assert_eq!(streaming, one_shot);
    }

    #[test]
    fn streaming_md5_uppercases_soft_masked_bases() {
        // Soft-masked FASTAs (Ensembl/Gencode default) encode repeat
        // regions as lowercase `acgtn`. The MD5 contract is "MD5 of
        // the uppercase bases" — the streaming pipeline must
        // uppercase before feeding md5.
        let (_dir, fasta_path) = build_wrapped_fasta(&[ContigFixture {
            name: "chr0",
            seq: "ACGTacgtNn",
            line_bases: 5,
        }]);
        let index = fai::fs::read(fasta_path.with_extension("fa.fai")).expect("fai");
        let record = &index.as_ref()[0];

        let streaming = compute_contig_md5_streaming(&fasta_path, record.offset(), 10).unwrap();
        let one_shot = one_shot_md5(b"ACGTACGTNN");
        assert_eq!(streaming, one_shot);
    }

    #[test]
    fn streaming_md5_truncated_contig_errors() {
        // .psp claims the contig is 100 bases but the FASTA only has
        // 10. The streaming function consumes whatever is there,
        // hits EOF with `remaining > 0`, and surfaces UnexpectedEof.
        let (_dir, fasta_path) = build_wrapped_fasta(&[ContigFixture {
            name: "chr0",
            seq: "ACGTACGTAC",
            line_bases: 10,
        }]);
        let index = fai::fs::read(fasta_path.with_extension("fa.fai")).expect("fai");
        let record = &index.as_ref()[0];

        let err = compute_contig_md5_streaming(&fasta_path, record.offset(), 100).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn verify_detects_md5_mismatch() {
        // End-to-end through `verify_fasta_matches_psp_chromosomes`:
        // build a FASTA whose MD5 we know, but tell the function to
        // expect a different MD5 (representing a "right basename,
        // wrong genome build" failure). Must surface `Md5Mismatch`
        // with the right contig name + actual vs. expected.
        let (_dir, fasta_path) = build_wrapped_fasta(&[ContigFixture {
            name: "chr0",
            seq: "ACGTACGTACGT",
            line_bases: 12,
        }]);
        let actual = one_shot_md5(b"ACGTACGTACGT");
        let chromosomes = vec![ParsedChromosome {
            name: "chr0".into(),
            length: 12,
            md5: "deadbeefdeadbeefdeadbeefdeadbeef".into(),
        }];
        let err = verify_fasta_matches_psp_chromosomes(&fasta_path, &chromosomes).unwrap_err();
        match err {
            FastaVerifyError::Md5Mismatch {
                contig,
                fasta_md5,
                psp_md5,
            } => {
                assert_eq!(contig, "chr0");
                assert_eq!(fasta_md5, actual);
                assert_eq!(psp_md5, "deadbeefdeadbeefdeadbeefdeadbeef");
            }
            other => panic!("expected Md5Mismatch, got {other:?}"),
        }
    }

    #[test]
    fn verify_missing_contig_in_fasta() {
        // .psp declares `chr_missing` but FASTA only has `chr0`. The
        // FASTA-index lookup must return `FetchFailed { source: kind
        // NotFound }`.
        let (_dir, fasta_path) = build_wrapped_fasta(&[ContigFixture {
            name: "chr0",
            seq: "ACGTACGTACGT",
            line_bases: 12,
        }]);
        let chromosomes = vec![ParsedChromosome {
            name: "chr_missing".into(),
            length: 12,
            // Value irrelevant — the lookup fails before the MD5 is
            // computed.
            md5: "deadbeefdeadbeefdeadbeefdeadbeef".into(),
        }];
        let err = verify_fasta_matches_psp_chromosomes(&fasta_path, &chromosomes).unwrap_err();
        match err {
            FastaVerifyError::FetchFailed { contig, source } => {
                assert_eq!(contig, "chr_missing");
                assert_eq!(source.kind(), io::ErrorKind::NotFound);
            }
            other => panic!("expected FetchFailed, got {other:?}"),
        }
    }

    #[test]
    fn verify_accepts_matching_md5_across_multiple_contigs() {
        // Positive path: build a FASTA with three differently-wrapped
        // contigs, compute the expected MD5 for each via the one-shot
        // path, and assert that `verify` accepts them. Exercises the
        // rayon `par_iter` happy path with > 1 worker.
        let seq_a = "ACGTACGTACGT".to_string();
        let seq_b: String = (0..150).map(|i| b"NACGT"[i % 5] as char).collect();
        let seq_c = "T".to_string();
        let (_dir, fasta_path) = build_wrapped_fasta(&[
            ContigFixture {
                name: "chr0",
                seq: &seq_a,
                line_bases: 12,
            },
            ContigFixture {
                name: "chr1",
                seq: &seq_b,
                line_bases: 60,
            },
            ContigFixture {
                name: "chr2",
                seq: &seq_c,
                line_bases: 1,
            },
        ]);
        let chromosomes = vec![
            ParsedChromosome {
                name: "chr0".into(),
                length: seq_a.len() as u32,
                md5: one_shot_md5(seq_a.as_bytes()),
            },
            ParsedChromosome {
                name: "chr1".into(),
                length: seq_b.len() as u32,
                md5: one_shot_md5(seq_b.as_bytes()),
            },
            ParsedChromosome {
                name: "chr2".into(),
                length: seq_c.len() as u32,
                md5: one_shot_md5(seq_c.as_bytes()),
            },
        ];
        verify_fasta_matches_psp_chromosomes(&fasta_path, &chromosomes).expect("verify");
    }
}
