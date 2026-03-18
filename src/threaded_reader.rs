use std::io::{self, Read};
use std::sync::mpsc::{self, Receiver, SyncSender};
use std::thread::{self, JoinHandle};

const CHUNK_SIZE: usize = 64 * 1024;
const CHANNEL_CAPACITY: usize = 4;

/// A reader that performs I/O on a background thread, allowing the consumer
/// to process previously-read data while the next chunk is being produced.
///
/// This is useful for wrapping slow readers (e.g. gzip decompressors) so that
/// decompression and parsing can happen in parallel on separate threads.
pub struct ThreadedReader {
    receiver: Receiver<io::Result<Vec<u8>>>,
    current_chunk: Vec<u8>,
    pos: usize,
    handle: Option<JoinHandle<()>>,
}

impl ThreadedReader {
    pub fn new<R: Read + Send + 'static>(mut reader: R) -> Self {
        let (sender, receiver): (SyncSender<io::Result<Vec<u8>>>, _) =
            mpsc::sync_channel(CHANNEL_CAPACITY);

        let handle = thread::spawn(move || {
            let mut buf = vec![0u8; CHUNK_SIZE];
            loop {
                match reader.read(&mut buf) {
                    Ok(0) => break,
                    Ok(n) => {
                        if sender.send(Ok(buf[..n].to_vec())).is_err() {
                            break;
                        }
                    }
                    Err(e) => {
                        let _ = sender.send(Err(e));
                        break;
                    }
                }
            }
        });

        ThreadedReader {
            receiver,
            current_chunk: Vec::new(),
            pos: 0,
            handle: Some(handle),
        }
    }
}

impl Read for ThreadedReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos >= self.current_chunk.len() {
            match self.receiver.recv() {
                Ok(Ok(chunk)) => {
                    self.current_chunk = chunk;
                    self.pos = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Ok(0),
            }
        }

        let remaining = &self.current_chunk[self.pos..];
        let to_copy = remaining.len().min(buf.len());
        buf[..to_copy].copy_from_slice(&remaining[..to_copy]);
        self.pos += to_copy;
        Ok(to_copy)
    }
}

impl Drop for ThreadedReader {
    fn drop(&mut self) {
        // Drop the receiver first so the sender fails and the thread exits
        drop(std::mem::replace(
            &mut self.receiver,
            mpsc::sync_channel(0).1,
        ));
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}