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
//! in the `pop_var_caller` and `var_calling::vcf_writer` modules.
//!
//! Everything here is `pub(crate)`; the helpers are not part of the
//! library's public surface (they would just be re-implemented by
//! any downstream consumer needing them).

use std::ffi::OsString;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

/// Buffered-I/O capacity used by every `BufReader::with_capacity` /
/// `BufWriter::with_capacity` in the cohort CLI's file I/O. 64 KiB
/// matches the value the per-sample `.psp` reader/writer documents
/// as the recommended minimum (see
/// [`crate::per_sample_pileup::psp::reader::PspReader::new`]).
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
}
