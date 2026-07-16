//! A windowed, caller-evictable, **raw** single-contig reference reader.
//!
//! ng's copy of `fasta::ManualEvictChromRefFetcher`, with exactly one behavioural
//! difference: **it does not canonicalise**. The bases land in the buffer verbatim
//! — soft-mask, IUPAC ambiguity codes and all.
//!
//! ## Why a copy, when the difference is one line
//!
//! Production's fetcher folds every byte to `{A,C,G,T,N}` **as it fills its
//! buffer** (`state.buf.push(canonicalise(b))`), so by the time ng can see the
//! bytes the raw ones are gone. There is no flag, and adding one would edit frozen
//! production (owner, 2026-07-16). Nor is there anything else to reach for:
//! production has no *windowed* raw reader at all (`RawContigRefCache` holds a
//! whole contig — that is `ResidentRefSeq`'s shape), and its `.fai` machinery
//! (`ContigFai`, `open_contig`) is **private**, so even the offset arithmetic
//! cannot be reused across the fence.
//!
//! **Why raw matters** (`typed_regions.md` §6): the STR catalog reads the FASTA
//! verbatim and embeds whatever it is handed, and `Locus` compares **by value** —
//! so canonical bytes would make every locus containing an IUPAC code compare
//! unequal to the catalog's. That silently breaks the parity oracle on any
//! assembly carrying them. `ref_seq.md` parked this as a YAGNI ("add a `RawRefSeq`
//! impl only if a windowed consumer ever actually needs raw bytes"); the
//! typed-region walk is that consumer, so the YAGNI is **spent, not violated**.
//!
//! ## Why the buffer is raw and canonical is the view
//!
//! Two readers were rejected on the spec's own grounds (`typed_regions.md` §6:
//! *"one reader for the whole run, sliding forward, never rebuilt per region —
//! that is what cost 14.6 GB of peak RSS"*). With one buffer the direction is
//! forced: `canonicalise` is a pure per-byte function, so **canonical is derivable
//! from raw and raw is not derivable from canonical**. Hence raw in the buffer,
//! canonical produced on the way out ([`super::ref_seq::WindowedRefSeq`]).
//!
//! ## What is *not* different
//!
//! Everything else is production's design, deliberately: the `.fai` validation
//! (which exists because `--reference` is attacker-influenced and a `line_bases`
//! of 0 divides by zero), the base→file-offset arithmetic, the bidirectional
//! extend, and `evict_before`'s drain-keeping-capacity. That arithmetic is fiddly
//! and well-tested, and this copy stays diffable against it on purpose.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use noodles_fasta::fai;

use crate::fasta::ChromRefFetchError;

/// Raw FASTA read chunk. Production's `STREAMING_REF_FILE_READ_CHUNK` is
/// `pub(crate)`-invisible here; the value is the same 64 KiB and is not tuned.
const FILE_READ_CHUNK: usize = 64 * 1024;

/// Pared-down `.fai` record — just the fields the reader needs.
#[derive(Debug, Clone, Copy)]
struct ContigFai {
    /// Byte offset of the first base in the FASTA file.
    seq_offset: u64,
    /// Total bases in the contig.
    length: u32,
    /// Bases per text line (excluding the trailing newline).
    line_bases: u32,
    /// Line width including the trailing newline.
    line_width: u32,
}

impl ContigFai {
    /// Validate a parsed `.fai` record. noodles' parser accepts values that would
    /// crash the byte-offset math: `line_bases = 0` divides by zero in
    /// [`Self::base_to_file_offset`], and `line_width < line_bases` yields wrong
    /// offsets (the trailing newline cannot have negative width).
    ///
    /// The `--reference` path is attacker-influenced, so this guards every
    /// constructor from a panic-on-construction surface.
    fn validate(&self, chrom_name: &str) -> Result<(), io::Error> {
        if self.line_bases == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "malformed .fai for contig {chrom_name}: line_bases = 0 \
                     (would divide-by-zero in offset arithmetic)"
                ),
            ));
        }
        if self.line_width < self.line_bases {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "malformed .fai for contig {chrom_name}: line_width ({}) \
                     < line_bases ({}) — line_width must include the trailing newline",
                    self.line_width, self.line_bases,
                ),
            ));
        }
        Ok(())
    }

    /// File-byte offset of the (1-based) base coordinate `b`.
    fn base_to_file_offset(&self, b: u32) -> u64 {
        let zero_based = (b - 1) as u64;
        let line_idx = zero_based / self.line_bases as u64;
        let in_line = zero_based % self.line_bases as u64;
        self.seq_offset + line_idx * self.line_width as u64 + in_line
    }
}

/// Open `chrom_name` in the FASTA at `fasta_path`, reading its sibling `.fai`
/// (or `fai_path`).
fn open_contig(
    fasta_path: &Path,
    fai_path: Option<&Path>,
    chrom_name: &str,
) -> Result<(ContigFai, File), ChromRefFetchError> {
    let mut fai_pathbuf: OsString;
    let fai_path_owned: PathBuf;
    let fai_path_resolved: &Path = match fai_path {
        Some(p) => p,
        None => {
            fai_pathbuf = fasta_path.as_os_str().to_os_string();
            fai_pathbuf.push(".fai");
            fai_path_owned = PathBuf::from(fai_pathbuf);
            &fai_path_owned
        }
    };
    let wrap_io = |source: io::Error| ChromRefFetchError::Io {
        chrom_name: chrom_name.to_string(),
        source,
    };
    let index = fai::fs::read(fai_path_resolved).map_err(wrap_io)?;
    let record = index
        .as_ref()
        .iter()
        .find(|r| AsRef::<[u8]>::as_ref(r.name()) == chrom_name.as_bytes())
        .ok_or_else(|| {
            wrap_io(io::Error::new(
                io::ErrorKind::NotFound,
                format!("contig {chrom_name} not in FASTA index"),
            ))
        })?;
    let length = u32::try_from(record.length()).map_err(|_| {
        wrap_io(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "contig {} length {} exceeds u32::MAX",
                chrom_name,
                record.length()
            ),
        ))
    })?;
    let line_bases = u32::try_from(record.line_bases()).map_err(|_| {
        wrap_io(io::Error::new(
            io::ErrorKind::InvalidData,
            ".fai line_bases exceeds u32::MAX",
        ))
    })?;
    let line_width = u32::try_from(record.line_width()).map_err(|_| {
        wrap_io(io::Error::new(
            io::ErrorKind::InvalidData,
            ".fai line_width exceeds u32::MAX",
        ))
    })?;
    let fai = ContigFai {
        seq_offset: record.offset(),
        length,
        line_bases,
        line_width,
    };
    fai.validate(chrom_name).map_err(wrap_io)?;
    let file = File::open(fasta_path).map_err(wrap_io)?;
    Ok((fai, file))
}

/// Read `n_bases` newline-stripped bases into `dst`.
///
/// **The one line that differs from production**, and the reason this file
/// exists: production's `read_uppercased_bases` ends `dst.push(canonicalise(b))`;
/// this pushes `b`. Bases arrive **verbatim** — lower-case soft-mask preserved,
/// IUPAC codes preserved — because the walk's oracle compares locus bytes against
/// a catalog built from the verbatim FASTA (`typed_regions.md` §6).
fn read_raw_bases(reader: &mut File, dst: &mut Vec<u8>, n_bases: usize) -> io::Result<()> {
    let want_total = dst.len() + n_bases;
    let mut read_buf = [0u8; FILE_READ_CHUNK];
    while dst.len() < want_total {
        let n = reader.read(&mut read_buf)?;
        if n == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "FASTA ended with {} base(s) still expected (want_total {want_total})",
                    want_total - dst.len(),
                ),
            ));
        }
        for &b in &read_buf[..n] {
            if dst.len() >= want_total {
                break;
            }
            if b == b'\n' || b == b'\r' {
                continue;
            }
            // Verbatim — NOT `canonicalise(b)`. See the fn docs.
            dst.push(b);
        }
    }
    Ok(())
}

/// A single contig's raw bases, buffered and caller-evictable.
///
/// Not `Sync`: one instance per thread; concurrent fetches on one instance are
/// forbidden. `Send`, because every field is.
pub struct RawChromReader {
    chrom_name: String,
    fai: ContigFai,
    file: File,
    /// Newline-stripped bases currently resident, **verbatim** (not uppercased,
    /// not canonicalised) — the difference from production's fetcher.
    buf: Vec<u8>,
    /// 1-based contig coordinate of `buf[0]`. Meaningful only when `!buf.is_empty()`.
    buf_start_base: u32,
}

impl RawChromReader {
    /// Open the FASTA at `fasta_path` (its sibling `.fai` is used), bind
    /// `chrom_name`, and return a reader with an empty buffer.
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::Io`] (`NotFound`) if `chrom_name` is not in the index.
    /// - [`ChromRefFetchError::Io`] (`InvalidData`) if the `.fai` entry is malformed.
    /// - [`ChromRefFetchError::Io`] propagating a `.fai`-read or `File::open` failure.
    pub fn for_contig(fasta_path: &Path, chrom_name: &str) -> Result<Self, ChromRefFetchError> {
        let (fai, file) = open_contig(fasta_path, None, chrom_name)?;
        Ok(Self {
            chrom_name: chrom_name.to_string(),
            fai,
            file,
            buf: Vec::new(),
            buf_start_base: 1,
        })
    }

    /// Total bases in the bound contig.
    pub fn length(&self) -> u32 {
        self.fai.length
    }

    /// Raw bases `[start_1based, start_1based + length)` as a borrowed slice,
    /// extending the buffer forward or backward as needed.
    ///
    /// # Errors
    ///
    /// - [`ChromRefFetchError::InvalidStart`] if `start_1based == 0`.
    /// - [`ChromRefFetchError::OutOfBounds`] if the window runs past the contig.
    /// - [`ChromRefFetchError::Io`] on a refill failure or a `u32` range overflow.
    pub fn fetch(&mut self, start_1based: u32, length: u32) -> Result<&[u8], ChromRefFetchError> {
        if start_1based == 0 {
            return Err(ChromRefFetchError::InvalidStart);
        }
        let end_exclusive =
            start_1based
                .checked_add(length)
                .ok_or_else(|| ChromRefFetchError::Io {
                    chrom_name: self.chrom_name.clone(),
                    source: io::Error::new(io::ErrorKind::InvalidInput, "fetch range overflow"),
                })?;
        if end_exclusive.saturating_sub(1) > self.fai.length {
            return Err(ChromRefFetchError::OutOfBounds {
                chrom_name: self.chrom_name.clone(),
                chrom_length: self.fai.length,
                start: start_1based,
                end: end_exclusive,
            });
        }

        if self.buf.is_empty() {
            self.read_into_buffer_at(start_1based, length as usize)?;
            return Ok(&self.buf[..length as usize]);
        }

        let buf_end_exclusive = self.buf_start_base as u64 + self.buf.len() as u64;
        let want_start = start_1based as u64;
        let want_end_exclusive = end_exclusive as u64;
        let buf_start = self.buf_start_base as u64;

        if want_end_exclusive > buf_end_exclusive {
            let extra_bases = (want_end_exclusive - buf_end_exclusive) as usize;
            self.append_forward(extra_bases)?;
        }
        if want_start < buf_start {
            let extra_bases = (buf_start - want_start) as usize;
            self.prepend_backward(start_1based, extra_bases)?;
        }

        let local_start = (start_1based - self.buf_start_base) as usize;
        let local_end = local_start + length as usize;
        Ok(&self.buf[local_start..local_end])
    }

    /// Drop buffered bases whose 1-based coordinate is `< pos`, keeping the
    /// allocation. The caller-driven memory bound. No-op when the buffer is empty
    /// or `pos <= buf_start_base`. Correctness is preserved: a later fetch of an
    /// evicted position re-reads it.
    pub fn evict_before(&mut self, pos: u32) {
        if self.buf.is_empty() || pos <= self.buf_start_base {
            return;
        }
        let buf_end_exclusive = self.buf_start_base as u64 + self.buf.len() as u64;
        if (pos as u64) >= buf_end_exclusive {
            self.buf.clear();
            self.buf_start_base = pos;
            return;
        }
        let drop_count = (pos - self.buf_start_base) as usize;
        self.buf.drain(..drop_count);
        self.buf_start_base = pos;
    }

    /// Current resident buffer length. Test/diagnostic.
    #[cfg(test)]
    pub fn buf_len(&self) -> usize {
        self.buf.len()
    }

    /// Read `n_bases` into an empty buffer starting at 1-based `start`.
    fn read_into_buffer_at(
        &mut self,
        start_1based: u32,
        n_bases: usize,
    ) -> Result<(), ChromRefFetchError> {
        debug_assert!(self.buf.is_empty());
        let offset = self.fai.base_to_file_offset(start_1based);
        self.seek_to(offset)?;
        let remaining = (self.fai.length - start_1based + 1) as usize;
        let target_len = n_bases.min(remaining);
        self.buf.reserve(target_len);
        self.buf_start_base = start_1based;
        let chrom = self.chrom_name.clone();
        read_raw_bases(&mut self.file, &mut self.buf, target_len).map_err(|source| {
            ChromRefFetchError::Io {
                chrom_name: chrom,
                source,
            }
        })
    }

    /// Append `extra_bases` past the buffer's current end.
    fn append_forward(&mut self, extra_bases: usize) -> Result<(), ChromRefFetchError> {
        let next_base = self.buf_start_base + self.buf.len() as u32;
        let offset = self.fai.base_to_file_offset(next_base);
        self.seek_to(offset)?;
        let remaining = (self.fai.length - next_base + 1) as usize;
        let take = extra_bases.min(remaining);
        self.buf.reserve(take);
        let chrom = self.chrom_name.clone();
        read_raw_bases(&mut self.file, &mut self.buf, take).map_err(|source| {
            ChromRefFetchError::Io {
                chrom_name: chrom,
                source,
            }
        })
    }

    /// Prepend `extra_bases` before the buffer's current start, splicing them in.
    fn prepend_backward(
        &mut self,
        new_start: u32,
        extra_bases: usize,
    ) -> Result<(), ChromRefFetchError> {
        debug_assert!(new_start < self.buf_start_base);
        debug_assert_eq!((self.buf_start_base - new_start) as usize, extra_bases);

        let offset = self.fai.base_to_file_offset(new_start);
        self.seek_to(offset)?;
        let mut prefix = Vec::with_capacity(extra_bases);
        let chrom = self.chrom_name.clone();
        read_raw_bases(&mut self.file, &mut prefix, extra_bases).map_err(|source| {
            ChromRefFetchError::Io {
                chrom_name: chrom,
                source,
            }
        })?;
        // `splice` does the memmove and reuses capacity where it can.
        self.buf.splice(..0, prefix);
        self.buf_start_base = new_start;
        Ok(())
    }

    fn seek_to(&mut self, offset: u64) -> Result<(), ChromRefFetchError> {
        self.file
            .seek(SeekFrom::Start(offset))
            .map(|_| ())
            .map_err(|source| ChromRefFetchError::Io {
                chrom_name: self.chrom_name.clone(),
                source,
            })
    }
}
