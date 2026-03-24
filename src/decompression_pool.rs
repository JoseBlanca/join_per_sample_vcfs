use std::collections::VecDeque;
use std::io::{self, Read};
use std::sync::{Arc, Condvar, Mutex};
use std::thread::{self, JoinHandle};

const CHUNK_SIZE: usize = 64 * 1024;
const BUFFER_CAPACITY: usize = 4;

/// Shared state for a single reader's output buffer.
struct SharedBuffer {
    inner: Mutex<BufferInner>,
    /// Notifies the consumer when data becomes available.
    data_available: Condvar,
}

struct BufferInner {
    chunks: VecDeque<io::Result<Vec<u8>>>,
    done: bool,
    cancelled: bool,
}

/// A decompression job: a reader paired with its output buffer.
struct Job {
    reader: Box<dyn Read + Send>,
    buffer: Arc<SharedBuffer>,
}

/// Shared state for the pool's work queue.
struct PoolInner {
    jobs: VecDeque<Job>,
    shutdown: bool,
}

struct PoolShared {
    inner: Mutex<PoolInner>,
    /// Notifies workers when new work is available or buffers have space.
    work_available: Condvar,
}

/// A pool of threads that decompress gzipped VCF data.
///
/// Instead of spawning one thread per file, a fixed number of worker threads
/// service all readers cooperatively. Each worker reads one chunk at a time
/// from a reader, sends it to the reader's buffer, then moves on to the next
/// reader. This bounds thread count regardless of how many files are open.
pub struct DecompressionPool {
    shared: Arc<PoolShared>,
    workers: Vec<JoinHandle<()>>,
}

impl DecompressionPool {
    /// Create a new pool with the given number of worker threads.
    pub fn new(num_threads: usize) -> Self {
        let shared = Arc::new(PoolShared {
            inner: Mutex::new(PoolInner {
                jobs: VecDeque::new(),
                shutdown: false,
            }),
            work_available: Condvar::new(),
        });

        let workers = (0..num_threads)
            .map(|_| {
                let shared = shared.clone();
                thread::spawn(move || worker_loop(shared))
            })
            .collect();

        DecompressionPool { shared, workers }
    }

    /// Register a reader with the pool and return a `PooledReader` that
    /// implements `Read`. The pool's worker threads will decompress data
    /// from the reader in the background.
    pub fn register<R: Read + Send + 'static>(&self, reader: R) -> PooledReader {
        let buffer = Arc::new(SharedBuffer {
            inner: Mutex::new(BufferInner {
                chunks: VecDeque::with_capacity(BUFFER_CAPACITY),
                done: false,
                cancelled: false,
            }),
            data_available: Condvar::new(),
        });

        let job = Job {
            reader: Box::new(reader),
            buffer: buffer.clone(),
        };

        {
            let mut pool = self.shared.inner.lock().unwrap();
            pool.jobs.push_back(job);
        }
        self.shared.work_available.notify_one();

        PooledReader {
            buffer,
            pool_shared: self.shared.clone(),
            current_chunk: Vec::new(),
            pos: 0,
        }
    }
}

impl Drop for DecompressionPool {
    fn drop(&mut self) {
        {
            let mut pool = self.shared.inner.lock().unwrap();
            pool.shutdown = true;
            // Mark all queued jobs as done so consumers don't block.
            for job in pool.jobs.drain(..) {
                let mut inner = job.buffer.inner.lock().unwrap();
                inner.done = true;
                job.buffer.data_available.notify_one();
            }
        }
        self.shared.work_available.notify_all();
        for handle in self.workers.drain(..) {
            let _ = handle.join();
        }
    }
}

enum JobResult {
    /// Chunk sent, requeue for more.
    Progress,
    /// Buffer full, requeue but skip for now.
    BufferFull,
    /// Reader is done (EOF, error, or cancelled).
    Done,
}

fn process_one_chunk(job: &mut Job, buf: &mut Vec<u8>) -> JobResult {
    // Check if cancelled or buffer full before reading.
    {
        let inner = job.buffer.inner.lock().unwrap();
        if inner.cancelled {
            return JobResult::Done;
        }
        if inner.chunks.len() >= BUFFER_CAPACITY {
            return JobResult::BufferFull;
        }
    }

    // Read one chunk (no lock held — only this worker touches the reader).
    match job.reader.read(buf) {
        Ok(0) => {
            let mut inner = job.buffer.inner.lock().unwrap();
            inner.done = true;
            job.buffer.data_available.notify_one();
            JobResult::Done
        }
        Ok(n) => {
            let chunk = buf[..n].to_vec();
            let mut inner = job.buffer.inner.lock().unwrap();
            if inner.cancelled {
                return JobResult::Done;
            }
            inner.chunks.push_back(Ok(chunk));
            job.buffer.data_available.notify_one();
            JobResult::Progress
        }
        Err(e) => {
            let mut inner = job.buffer.inner.lock().unwrap();
            inner.chunks.push_back(Err(e));
            inner.done = true;
            job.buffer.data_available.notify_one();
            JobResult::Done
        }
    }
}

fn mark_job_done(job: &Job) {
    let mut inner = job.buffer.inner.lock().unwrap();
    inner.done = true;
    job.buffer.data_available.notify_one();
}

fn worker_loop(shared: Arc<PoolShared>) {
    let mut buf = vec![0u8; CHUNK_SIZE];

    loop {
        // Pop a job that can make progress.
        let mut job = {
            let mut pool = shared.inner.lock().unwrap();
            loop {
                if pool.shutdown {
                    return;
                }
                if let Some(job) = pool.jobs.pop_front() {
                    break job;
                }
                // No jobs — wait for work.
                pool = shared.work_available.wait(pool).unwrap();
            }
        };

        match process_one_chunk(&mut job, &mut buf) {
            JobResult::Progress => {
                let mut pool = shared.inner.lock().unwrap();
                if pool.shutdown {
                    mark_job_done(&job);
                } else {
                    pool.jobs.push_back(job);
                    // Wake another worker if there are more jobs.
                    shared.work_available.notify_one();
                }
            }
            JobResult::BufferFull => {
                let mut pool = shared.inner.lock().unwrap();
                if pool.shutdown {
                    mark_job_done(&job);
                    continue;
                }
                // Can't make progress — put the job back and try another.
                // To avoid busy-looping when ALL jobs are buffer-full, we check
                // if there's any other job we could service. If not, we park.
                let queue_was_empty = pool.jobs.is_empty();
                pool.jobs.push_back(job);
                if queue_was_empty {
                    // We're the only job and buffer is full. Wait until the
                    // consumer frees space (it will notify work_available).
                    pool = shared.work_available.wait(pool).unwrap();
                }
                // Otherwise there are other jobs — loop around to try one.
            }
            JobResult::Done => {
                // Don't requeue.
            }
        }
    }
}

/// A reader backed by a `DecompressionPool`. Implements `Read` by pulling
/// pre-decompressed chunks from a shared buffer that pool workers fill.
pub struct PooledReader {
    buffer: Arc<SharedBuffer>,
    pool_shared: Arc<PoolShared>,
    current_chunk: Vec<u8>,
    pos: usize,
}

impl Read for PooledReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos >= self.current_chunk.len() {
            let mut inner = self.buffer.inner.lock().unwrap();
            loop {
                if let Some(result) = inner.chunks.pop_front() {
                    drop(inner);
                    // Notify pool that buffer space is available.
                    self.pool_shared.work_available.notify_one();
                    self.current_chunk = result?;
                    self.pos = 0;
                    break;
                }
                if inner.done {
                    return Ok(0);
                }
                inner = self.buffer.data_available.wait(inner).unwrap();
            }
        }

        let remaining = &self.current_chunk[self.pos..];
        let to_copy = remaining.len().min(buf.len());
        buf[..to_copy].copy_from_slice(&remaining[..to_copy]);
        self.pos += to_copy;
        Ok(to_copy)
    }
}

impl Drop for PooledReader {
    fn drop(&mut self) {
        // Signal the pool to stop working on this reader.
        let mut inner = self.buffer.inner.lock().unwrap();
        inner.cancelled = true;
        drop(inner);
        // Wake workers so they can skip this cancelled job.
        self.pool_shared.work_available.notify_all();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_single_reader_small_data() {
        let data = b"hello world";
        let pool = DecompressionPool::new(2);
        let mut reader = pool.register(Cursor::new(data.to_vec()));

        let mut output = Vec::new();
        reader.read_to_end(&mut output).unwrap();
        assert_eq!(output, data);
    }

    #[test]
    fn test_single_reader_large_data() {
        // Data larger than CHUNK_SIZE to test multi-chunk reading.
        let data: Vec<u8> = (0..CHUNK_SIZE * 3 + 42).map(|i| (i % 256) as u8).collect();
        let pool = DecompressionPool::new(2);
        let mut reader = pool.register(Cursor::new(data.clone()));

        let mut output = Vec::new();
        reader.read_to_end(&mut output).unwrap();
        assert_eq!(output, data);
    }

    #[test]
    fn test_multiple_readers_concurrent() {
        let pool = DecompressionPool::new(2);
        let data_a: Vec<u8> = vec![0xAA; CHUNK_SIZE * 2];
        let data_b: Vec<u8> = vec![0xBB; CHUNK_SIZE * 2];
        let data_c: Vec<u8> = vec![0xCC; CHUNK_SIZE * 2];

        let mut reader_a = pool.register(Cursor::new(data_a.clone()));
        let mut reader_b = pool.register(Cursor::new(data_b.clone()));
        let mut reader_c = pool.register(Cursor::new(data_c.clone()));

        let mut out_a = Vec::new();
        let mut out_b = Vec::new();
        let mut out_c = Vec::new();

        // Read interleaved to exercise pool scheduling.
        let mut buf = [0u8; 1024];
        loop {
            let mut progress = false;
            for (reader, out) in [
                (&mut reader_a, &mut out_a),
                (&mut reader_b, &mut out_b),
                (&mut reader_c, &mut out_c),
            ] {
                let n = reader.read(&mut buf).unwrap();
                if n > 0 {
                    out.extend_from_slice(&buf[..n]);
                    progress = true;
                }
            }
            if !progress {
                break;
            }
        }

        assert_eq!(out_a, data_a);
        assert_eq!(out_b, data_b);
        assert_eq!(out_c, data_c);
    }

    #[test]
    fn test_many_readers_few_threads() {
        // More readers than pool threads — should not deadlock.
        let pool = DecompressionPool::new(2);
        let n_readers = 20;
        let data: Vec<u8> = (0..CHUNK_SIZE + 100).map(|i| (i % 256) as u8).collect();

        let mut readers: Vec<_> = (0..n_readers)
            .map(|_| pool.register(Cursor::new(data.clone())))
            .collect();

        for reader in &mut readers {
            let mut output = Vec::new();
            reader.read_to_end(&mut output).unwrap();
            assert_eq!(output, data);
        }
    }

    #[test]
    fn test_reader_drop_cancels_job() {
        let pool = DecompressionPool::new(2);
        let data = vec![0u8; CHUNK_SIZE * 10];

        {
            let mut reader = pool.register(Cursor::new(data));
            let mut buf = [0u8; 1024];
            // Read just a bit, then drop.
            reader.read(&mut buf).unwrap();
        }
        // Pool should not hang — the cancelled job is cleaned up.

        // Verify pool still works after cancellation.
        let mut reader2 = pool.register(Cursor::new(b"still works".to_vec()));
        let mut output = Vec::new();
        reader2.read_to_end(&mut output).unwrap();
        assert_eq!(output, b"still works");
    }

    #[test]
    fn test_empty_reader() {
        let pool = DecompressionPool::new(2);
        let mut reader = pool.register(Cursor::new(Vec::new()));

        let mut output = Vec::new();
        reader.read_to_end(&mut output).unwrap();
        assert!(output.is_empty());
    }

    #[test]
    fn test_error_propagation() {
        struct FailingReader;
        impl Read for FailingReader {
            fn read(&mut self, _buf: &mut [u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::Other, "read failed"))
            }
        }

        let pool = DecompressionPool::new(2);
        let mut reader = pool.register(FailingReader);

        let mut buf = [0u8; 1024];
        let result = reader.read(&mut buf);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().to_string(), "read failed");
    }

    #[test]
    fn test_pool_drop_with_active_readers() {
        let data = vec![0u8; CHUNK_SIZE * 5];
        let mut reader;
        {
            let pool = DecompressionPool::new(2);
            reader = pool.register(Cursor::new(data));
            let mut buf = [0u8; 1024];
            reader.read(&mut buf).unwrap();
            // Pool drops here while reader is still alive.
        }
        // Reader should still be able to drain what was already buffered
        // or get EOF once the pool is gone.
        let mut output = Vec::new();
        // This shouldn't hang — the pool threads are joined, the buffer
        // either has data or is marked done/cancelled.
        let _ = reader.read_to_end(&mut output);
    }

    #[test]
    fn test_gzip_decompression() {
        use flate2::write::GzEncoder;
        use flate2::read::MultiGzDecoder;
        use flate2::Compression;
        use std::io::Write;

        // Compress some data.
        let original = b"The quick brown fox jumps over the lazy dog.\n".repeat(1000);
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(&original).unwrap();
        let compressed = encoder.finish().unwrap();

        // Decompress through the pool.
        let pool = DecompressionPool::new(2);
        let gz_reader = MultiGzDecoder::new(Cursor::new(compressed));
        let mut reader = pool.register(gz_reader);

        let mut output = Vec::new();
        reader.read_to_end(&mut output).unwrap();
        assert_eq!(output, original);
    }
}
