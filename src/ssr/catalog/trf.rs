//! TRF-mod invocation + BED parsing (architecture
//! [`ssr_catalog.md`](../../../doc/devel/architecture/ssr_catalog.md) §2-3).
//!
//! Locates the `trf-mod` binary ([`locate_trf_mod`], §2.4 layered discovery),
//! reads its version ([`version`]), and runs it per contig via temp files
//! ([`run_on_contig`]) — no stdin/stdout pipes, so no deadlock and termination
//! is a clean exit-status check. Its tab-separated BED is parsed by
//! [`parse_bed_line`] into [`TrfRecord`]s for [`super::postprocess`].
//!
//! trf-mod (v4.10.0, lh3 fork) emits 10 columns
//! (`trfrun.h::trf_print_bed`): `ctg  start  end  period  copyNum  fracMatch
//! fracGap  score  entropy  pattern`, where `start = first - 1` (0-based) and
//! `end = last`, so `[start, end)` is 0-based half-open. We keep `start`/`end`,
//! the authoritative `period`, `fracMatch` (a sanity field — purity is
//! recomputed from the tract, §4), `score` (an early accept-gate), and the
//! consensus `pattern` (sanity; the catalog motif comes from the reference
//! tract, §3). `copyNum` / `fracGap` / `entropy` are parsed-and-dropped.

/// One parsed TRF-mod BED row, reduced to what the catalog needs.
/// Coordinates are **0-based half-open** (`[start, end)`), as trf-mod emits.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TrfRecord {
    /// Tract start (0-based, inclusive).
    pub(crate) start: u32,
    /// Tract end (0-based, exclusive).
    pub(crate) end: u32,
    /// Authoritative repeat period (motif length in bases). May differ from
    /// `pattern.len()` (architecture §3), so the catalog motif is taken as the
    /// first `period` bases of the tract, not from `pattern`.
    pub(crate) period: u16,
    /// TRF percent-match fraction in `[0, 1]` — a sanity field only; the
    /// catalog `purity_fraction` is recomputed from the trimmed tract (§4).
    pub(crate) frac_match: f32,
    /// TRF alignment score — used as an early accept-gate.
    pub(crate) score: i32,
    /// TRF consensus pattern (sanity only; not the catalog motif).
    pub(crate) pattern: Box<[u8]>,
}

#[cfg(test)]
impl TrfRecord {
    /// Test constructor — builds a record from the fields `postprocess`
    /// reads. (`frac_match` is unused by the pipeline, so it is fixed to 1.0.)
    pub(crate) fn for_test(start: u32, end: u32, period: u16, score: i32, pattern: &[u8]) -> Self {
        Self {
            start,
            end,
            period,
            frac_match: 1.0,
            score,
            pattern: pattern.to_vec().into_boxed_slice(),
        }
    }
}

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use super::CatalogError;

/// The binary name we look for (and the file we install in the dev image).
const TRF_MOD_BIN: &str = "trf-mod";

/// Number of tab-separated columns trf-mod's BED carries.
const BED_COLUMNS: usize = 10;

/// Locate the `trf-mod` binary (architecture §2.4): an explicit `override`
/// first, then a copy beside our own executable, then `PATH`. Returns the first
/// existing path, or [`CatalogError::TrfModNotFound`].
pub(crate) fn locate_trf_mod(override_path: Option<&Path>) -> Result<PathBuf, CatalogError> {
    // 1. Explicit override — must exist, else a hard error (the user named it).
    if let Some(p) = override_path {
        return if p.is_file() {
            Ok(p.to_path_buf())
        } else {
            Err(CatalogError::TrfModNotFound)
        };
    }

    // 2. A copy shipped beside our own executable.
    if let Ok(exe) = std::env::current_exe()
        && let Some(dir) = exe.parent()
    {
        let cand = dir.join(TRF_MOD_BIN);
        if cand.is_file() {
            return Ok(cand);
        }
    }

    // 3. `PATH`.
    if let Some(paths) = std::env::var_os("PATH") {
        for dir in std::env::split_paths(&paths) {
            let cand = dir.join(TRF_MOD_BIN);
            if cand.is_file() {
                return Ok(cand);
            }
        }
    }

    Err(CatalogError::TrfModNotFound)
}

/// Read trf-mod's version string (`trf-mod -v`) for the catalog header. The
/// version line — `"Tandem Repeats Finder, Version 4.10.0"` — may land on
/// stdout or stderr depending on the build, so both are searched.
pub(crate) fn version(bin: &Path) -> Result<String, CatalogError> {
    let out = Command::new(bin)
        .arg("-v")
        .output()
        .map_err(|source| CatalogError::TrfSpawn {
            context: "spawn trf-mod -v",
            source,
        })?;
    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);
    stdout
        .lines()
        .chain(stderr.lines())
        .map(str::trim)
        .find(|l| l.contains("Version"))
        .map(|l| l.to_string())
        .ok_or(CatalogError::TrfVersion)
}

/// Run trf-mod on one contig via temp files under `temp_root` and parse its BED.
/// Collision-free across parallel contigs (each task gets a unique temp dir).
/// No stdin/stdout pipes: the contig is written to `input.fa`, trf-mod's stdout
/// is redirected to `output.bed`, and a non-zero exit becomes a hard error.
pub(crate) fn run_on_contig(
    bin: &Path,
    name: &str,
    seq: &[u8],
    temp_root: &Path,
) -> Result<Vec<TrfRecord>, CatalogError> {
    let dir = tempfile::Builder::new()
        .prefix("ssr-catalog-")
        .tempdir_in(temp_root)
        .map_err(|source| CatalogError::TrfSpawn {
            context: "create trf-mod temp dir",
            source,
        })?;
    let input = dir.path().join("input.fa");
    let output = dir.path().join("output.bed");

    // Write the contig as a one-record FASTA.
    {
        let mut f = std::fs::File::create(&input).map_err(|source| CatalogError::TrfSpawn {
            context: "create trf-mod input.fa",
            source,
        })?;
        let write = |f: &mut std::fs::File| -> std::io::Result<()> {
            f.write_all(b">")?;
            f.write_all(name.as_bytes())?;
            f.write_all(b"\n")?;
            f.write_all(seq)?;
            f.write_all(b"\n")
        };
        write(&mut f).map_err(|source| CatalogError::TrfSpawn {
            context: "write trf-mod input.fa",
            source,
        })?;
    }

    // Spawn trf-mod, redirecting its BED stdout to output.bed and discarding its
    // progress stderr. `status()` waits — a crash surfaces as a non-zero exit.
    let out_file = std::fs::File::create(&output).map_err(|source| CatalogError::TrfSpawn {
        context: "create trf-mod output.bed",
        source,
    })?;
    let status = Command::new(bin)
        .arg(&input)
        .stdout(std::process::Stdio::from(out_file))
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(|source| CatalogError::TrfSpawn {
            context: "spawn trf-mod",
            source,
        })?;
    if !status.success() {
        return Err(CatalogError::TrfRun {
            contig: name.to_string(),
            status: status.to_string(),
        });
    }

    // Parse the BED. trf-mod writes only data rows to stdout (progress goes to
    // stderr), so every non-blank line must parse — a malformed line fails
    // loudly (catches an upstream format change), per architecture §3.
    let content = std::fs::read_to_string(&output).map_err(|source| CatalogError::TrfSpawn {
        context: "read trf-mod output.bed",
        source,
    })?;
    let mut recs = Vec::new();
    for (i, line) in content.lines().enumerate() {
        let line = line.trim_end_matches(['\n', '\r']);
        if line.is_empty() {
            continue;
        }
        recs.push(parse_bed_line(line, name, i + 1)?);
    }
    Ok(recs)
}

/// Parse one trf-mod BED line into a [`TrfRecord`], asserting the 10-column
/// layout and that the contig-name column matches `expect_ctg`. `line` is the
/// 1-based source line for diagnostics.
pub(crate) fn parse_bed_line(
    line: &str,
    expect_ctg: &str,
    line_no: usize,
) -> Result<TrfRecord, CatalogError> {
    let cols: Vec<&str> = line.split('\t').collect();
    if cols.len() != BED_COLUMNS {
        return Err(CatalogError::TrfParse {
            line: line_no,
            reason: format!("expected {BED_COLUMNS} columns, got {}", cols.len()),
        });
    }
    if cols[0] != expect_ctg {
        return Err(CatalogError::TrfParse {
            line: line_no,
            reason: format!("contig column {:?} != expected {expect_ctg:?}", cols[0]),
        });
    }

    let field = |idx: usize| cols[idx];
    let num_err = |name: &str, raw: &str| CatalogError::TrfParse {
        line: line_no,
        reason: format!("column {name} = {raw:?} did not parse"),
    };

    let start = field(1)
        .parse::<u32>()
        .map_err(|_| num_err("start", field(1)))?;
    let end = field(2)
        .parse::<u32>()
        .map_err(|_| num_err("end", field(2)))?;
    let period = field(3)
        .parse::<u16>()
        .map_err(|_| num_err("period", field(3)))?;
    let frac_match = field(5)
        .parse::<f32>()
        .map_err(|_| num_err("fracMatch", field(5)))?;
    let score = field(7)
        .parse::<i32>()
        .map_err(|_| num_err("score", field(7)))?;
    let pattern = field(9).as_bytes().to_vec().into_boxed_slice();

    if end <= start {
        return Err(CatalogError::TrfParse {
            line: line_no,
            reason: format!("end {end} <= start {start}"),
        });
    }

    Ok(TrfRecord {
        start,
        end,
        period,
        frac_match,
        score,
        pattern,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A real trf-mod BED row (from `trf-mod TRF-mod/t/small_test.fasta`):
    /// `0  0  35  7  5.00  1.00  0.00  70  1.95  TCATCGG`.
    #[test]
    fn parse_bed_line_extracts_fields() {
        let row = "0\t0\t35\t7\t5.00\t1.00\t0.00\t70\t1.95\tTCATCGG";
        let rec = parse_bed_line(row, "0", 1).unwrap();
        assert_eq!(rec.start, 0);
        assert_eq!(rec.end, 35);
        assert_eq!(rec.period, 7);
        assert_eq!(rec.frac_match, 1.0);
        assert_eq!(rec.score, 70);
        assert_eq!(&*rec.pattern, b"TCATCGG");
    }

    #[test]
    fn parse_bed_line_rejects_wrong_column_count() {
        let err = parse_bed_line("chr1\t0\t35\t7", "chr1", 4).unwrap_err();
        assert!(
            matches!(err, CatalogError::TrfParse { line: 4, .. }),
            "got {err:?}"
        );
    }

    #[test]
    fn parse_bed_line_rejects_contig_name_mismatch() {
        let row = "chrX\t0\t35\t7\t5.00\t1.00\t0.00\t70\t1.95\tTCATCGG";
        let err = parse_bed_line(row, "chr1", 1).unwrap_err();
        assert!(
            matches!(err, CatalogError::TrfParse { ref reason, .. } if reason.contains("contig")),
            "got {err:?}"
        );
    }

    /// Integration: actually run the installed `trf-mod` on a synthetic contig
    /// and parse its BED. Skips (does not fail) when trf-mod is absent — the
    /// dev container has it at `/usr/local/bin/trf-mod`, so CI exercises this.
    #[test]
    fn run_on_contig_detects_a_synthetic_repeat() {
        let bin = match locate_trf_mod(None) {
            Ok(b) => b,
            Err(_) => {
                eprintln!("skipping run_on_contig test: trf-mod not found on this host");
                return;
            }
        };
        // 60 copies of "CAG" — comfortably above trf-mod's default min score.
        let mut seq = Vec::new();
        for _ in 0..60 {
            seq.extend_from_slice(b"CAG");
        }
        std::fs::create_dir_all("tmp").unwrap();
        let recs = run_on_contig(&bin, "synth", &seq, Path::new("tmp")).unwrap();
        assert!(
            recs.iter().any(|r| r.period == 3 && r.end - r.start >= 150),
            "expected a period-3 CAG repeat, got {recs:?}"
        );
    }
}
