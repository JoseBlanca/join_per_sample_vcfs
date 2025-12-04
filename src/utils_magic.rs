use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::Path;

pub fn are_gzipped_magic_bytes(first_bytes: &[u8]) -> bool {
    first_bytes.len() >= 2 && first_bytes[0] == 0x1f && first_bytes[1] == 0x8b
}

fn read_first_bytes(path: impl AsRef<Path>, num_bytes: usize) -> io::Result<Vec<u8>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    let mut buffer = vec![0u8; num_bytes];
    let bytes_read = reader.read(&mut buffer)?;
    buffer.truncate(bytes_read);

    Ok(buffer)
}

pub fn file_is_gzipped(path: impl AsRef<Path>) -> io::Result<bool> {
    let first_bytes = read_first_bytes(path, 4)?;
    Ok(are_gzipped_magic_bytes(&first_bytes))
}
