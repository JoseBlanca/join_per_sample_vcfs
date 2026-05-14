//! Impedance-match between the upstream BAQ / merged-CRAM stream
//! (`Iterator<Item = Result<PreparedRead, CramInputError>>`) and the
//! pileup walker (which takes `Iterator<Item = PreparedRead>` — no
//! `Result`).
//!
//! The walker is shape-locked to a non-`Result` item type so its body
//! does not have to carry an upstream-error generic; making that
//! change in the walker would touch every call site. The adapter
//! localises the mismatch: it forwards `Ok` items through, and on the
//! first `Err` it stashes the error inside a shared cell and reports
//! end-of-stream. The walker drains the records it has produced so
//! far cleanly and the orchestrator surfaces the stashed error
//! *after* the seam returns.

use std::cell::RefCell;
use std::rc::Rc;

use crate::per_sample_caller::errors::CramInputError;
use crate::per_sample_caller::pileup::PreparedRead;

/// Shared handle the orchestrator keeps so it can `take()` the
/// stashed error once the walker exhausts. Cheap to clone (one
/// `Rc::clone` per handle).
#[derive(Clone, Default)]
pub struct ErrorHandle(Rc<RefCell<Option<CramInputError>>>);

impl ErrorHandle {
    /// Take the stashed error, if any. Idempotent: subsequent calls
    /// return `None`.
    pub fn take(&self) -> Option<CramInputError> {
        self.0.borrow_mut().take()
    }
}

/// Adapter wrapping a fallible `PreparedRead` source. `Ok` items
/// pass through unchanged; the first `Err` is stashed and the
/// iterator reports end-of-stream forever after.
pub struct ErrorSheddingAdapter<I>
where
    I: Iterator<Item = Result<PreparedRead, CramInputError>>,
{
    inner: I,
    handle: ErrorHandle,
    done: bool,
}

impl<I> ErrorSheddingAdapter<I>
where
    I: Iterator<Item = Result<PreparedRead, CramInputError>>,
{
    pub fn new(inner: I) -> Self {
        Self {
            inner,
            handle: ErrorHandle::default(),
            done: false,
        }
    }

    /// Cloneable handle that the orchestrator inspects after the
    /// walker exhausts. Holds a refcounted reference into the same
    /// cell the adapter writes to.
    pub fn error_handle(&self) -> ErrorHandle {
        self.handle.clone()
    }
}

impl<I> Iterator for ErrorSheddingAdapter<I>
where
    I: Iterator<Item = Result<PreparedRead, CramInputError>>,
{
    type Item = PreparedRead;

    fn next(&mut self) -> Option<PreparedRead> {
        if self.done {
            return None;
        }
        match self.inner.next() {
            Some(Ok(r)) => Some(r),
            Some(Err(e)) => {
                *self.handle.0.borrow_mut() = Some(e);
                self.done = true;
                None
            }
            None => {
                self.done = true;
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use std::sync::Arc;

    use super::*;
    use crate::per_sample_caller::cram_input::CigarOp;
    use crate::per_sample_caller::pileup::MateRole;

    fn dummy_prepared(pos: u32) -> PreparedRead {
        PreparedRead {
            chrom_id: 0,
            alignment_start: pos,
            alignment_end: pos + 4,
            cigar: vec![CigarOp::Match(5)],
            seq: vec![b'A'; 5],
            bq_baq: vec![30; 5],
            mq_log_err: -3.0,
            is_reverse_strand: false,
            qname: Arc::<str>::from("r"),
            mate_role: MateRole::Solo,
            adaptor_boundary: None,
        }
    }

    fn dummy_err() -> CramInputError {
        CramInputError::NoInputs
    }

    /// Ok items pass through unchanged; the handle is empty after a
    /// clean exhaustion.
    #[test]
    fn ok_items_pass_through() {
        let input: Vec<Result<PreparedRead, CramInputError>> =
            vec![Ok(dummy_prepared(1)), Ok(dummy_prepared(2))];
        let mut adapter = ErrorSheddingAdapter::new(input.into_iter());
        let handle = adapter.error_handle();

        let pulled: Vec<_> = adapter.by_ref().collect();
        assert_eq!(pulled.len(), 2);
        assert!(handle.take().is_none(), "handle should be empty");
        assert!(adapter.next().is_none(), "adapter is exhausted");
    }

    /// The first `Err` is stashed and the adapter reports
    /// end-of-stream from that point on.
    #[test]
    fn first_error_is_stashed_then_iterator_ends() {
        let input: Vec<Result<PreparedRead, CramInputError>> = vec![
            Ok(dummy_prepared(1)),
            Err(dummy_err()),
            // Items past the error must not be observed even if
            // upstream would still yield them.
            Ok(dummy_prepared(99)),
        ];
        let mut adapter = ErrorSheddingAdapter::new(input.into_iter());
        let handle = adapter.error_handle();

        assert!(adapter.next().is_some(), "first Ok comes through");
        assert!(adapter.next().is_none(), "error becomes end-of-stream");
        assert!(adapter.next().is_none(), "still end-of-stream");

        let stashed = handle.take().expect("error must be stashed");
        assert!(matches!(stashed, CramInputError::NoInputs));
    }

    /// `ErrorHandle::take` is idempotent — calling it after a
    /// successful take returns `None`.
    #[test]
    fn handle_take_is_idempotent() {
        let _ = PathBuf::new(); // silence unused import on cold compiles
        let input: Vec<Result<PreparedRead, CramInputError>> = vec![Err(dummy_err())];
        let mut adapter = ErrorSheddingAdapter::new(input.into_iter());
        let handle = adapter.error_handle();
        assert!(adapter.next().is_none());
        assert!(handle.take().is_some());
        assert!(handle.take().is_none());
    }
}
