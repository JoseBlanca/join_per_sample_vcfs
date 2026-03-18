use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum MagicByteError {
    #[error("Insufficient bytes: got {got}, need at least {need}")]
    InsufficientBytes { got: usize, need: usize },

    #[error("There was a problem opening the file '{path}': {reason}")]
    ProblemOpeningFile { path: String, reason: String },

    #[error("There was a problem reading the buffer for the file '{path}': {reason}")]
    ProblemFillingBuffer { path: String, reason: String },
}

pub fn are_gzipped_magic_bytes(first_bytes: &[u8]) -> Result<bool, MagicByteError> {
    if first_bytes.len() < 2 {
        return Err(MagicByteError::InsufficientBytes {
            got: first_bytes.len(),
            need: 2,
        });
    }
    Ok(first_bytes[0] == 0x1f && first_bytes[1] == 0x8b)
}

fn read_first_bytes<P: AsRef<Path>>(path: &P, num_bytes: usize) -> Result<Vec<u8>, MagicByteError> {
    let file = File::open(path).map_err(|e| MagicByteError::ProblemOpeningFile {
        path: path.as_ref().to_string_lossy().to_string(),
        reason: e.to_string(),
    })?;
    let mut buf_reader = BufReader::new(file);

    let buffer = buf_reader
        .fill_buf()
        .map_err(|e| MagicByteError::ProblemFillingBuffer {
            path: path.as_ref().to_string_lossy().to_string(),
            reason: e.to_string(),
        })?;
    Ok(buffer[..num_bytes.min(buffer.len())].to_vec())
}

pub fn file_is_gzipped<P: AsRef<Path>>(path: &P) -> Result<bool, MagicByteError> {
    let first_bytes = read_first_bytes(path, 4)?;
    are_gzipped_magic_bytes(&first_bytes)
}
