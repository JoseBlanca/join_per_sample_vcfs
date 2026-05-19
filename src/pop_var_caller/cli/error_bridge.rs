//! Generic error-shedding adapter for "fallible item iterator →
//! non-`Result` item iterator" seams.
//!
//! Several downstream consumers in this crate (the pileup walker
//! taking `Iterator<Item = PreparedRead>`, the `PerPositionMerger`
//! input that the from-bam shim wraps a walker into) want a
//! non-`Result` item type but the upstream is fallible
//! (`Iterator<Item = Result<T, E>>`). [`ErrorSheddingAdapter`]
//! localises that mismatch: it forwards `Ok` items through, and on
//! the first `Err` it stashes the error inside a shared cell and
//! reports end-of-stream. The downstream drains the items it has
//! produced so far cleanly and the orchestrator surfaces the
//! stashed error *after* the seam returns.
//!
//! Parametric over both the item type and the error type so the
//! same machinery serves the Stage 1 CRAM-input seam
//! (`PreparedRead` / `CramInputError`) and the from-bam walker shim
//! (`PileupRecord` / `WalkerError`). M9 follow-up from the
//! 2026-05-19 cohort CLI review.

use std::cell::RefCell;
use std::rc::Rc;

/// Shared handle the orchestrator keeps so it can `take()` the
/// stashed error once the wrapped iterator exhausts. Cheap to clone
/// (one `Rc::clone` per handle).
///
/// Generic over the error type `E` — typically
/// [`CramInputError`](crate::per_sample_pileup::errors::CramInputError)
/// for the Stage 1 CRAM-input seam or
/// [`WalkerError`](crate::per_sample_pileup::pileup::WalkerError) for
/// the from-bam walker shim.
pub struct ErrorHandle<E>(Rc<RefCell<Option<E>>>);

impl<E> Default for ErrorHandle<E> {
    fn default() -> Self {
        Self(Rc::new(RefCell::new(None)))
    }
}

impl<E> Clone for ErrorHandle<E> {
    // Manual `Clone` because `#[derive(Clone)]` would add a `E: Clone`
    // bound that `Rc::clone` does not actually need — the smart
    // pointer is what gets cloned, not the contained `Option<E>`.
    fn clone(&self) -> Self {
        Self(Rc::clone(&self.0))
    }
}

impl<E> ErrorHandle<E> {
    /// Take the stashed error, if any. Idempotent: subsequent calls
    /// return `None`.
    pub fn take(&self) -> Option<E> {
        self.0.borrow_mut().take()
    }
}

/// Adapter wrapping a fallible item source. `Ok` items pass
/// through unchanged; the first `Err` is stashed in the
/// [`ErrorHandle`] and the iterator reports end-of-stream forever
/// after.
///
/// `T` is the item type the downstream wants
/// (e.g. [`PreparedRead`](crate::per_sample_pileup::pileup::PreparedRead),
/// [`PileupRecord`](crate::per_sample_pileup::pileup::PileupRecord)).
/// `E` is the upstream error type stored in the handle.
pub struct ErrorSheddingAdapter<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    inner: I,
    handle: ErrorHandle<E>,
    done: bool,
}

impl<I, T, E> ErrorSheddingAdapter<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    pub fn new(inner: I) -> Self {
        Self {
            inner,
            handle: ErrorHandle::default(),
            done: false,
        }
    }

    /// Cloneable handle that the orchestrator inspects after the
    /// downstream consumer exhausts. Holds a refcounted reference
    /// into the same cell the adapter writes to.
    pub fn error_handle(&self) -> ErrorHandle<E> {
        self.handle.clone()
    }
}

impl<I, T, E> Iterator for ErrorSheddingAdapter<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    type Item = T;

    fn next(&mut self) -> Option<T> {
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
    use crate::per_sample_pileup::cram_input::CigarOp;
    use crate::per_sample_pileup::errors::CramInputError;
    use crate::per_sample_pileup::pileup::{MateRole, PreparedRead};

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
    /// clean exhaustion. Concrete `T = PreparedRead, E = CramInputError`
    /// — the original Stage 1 CRAM-input use case.
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

    /// Generic over arbitrary `T` and `E`. Locks the M9 follow-up
    /// fix in place — the adapter must be reusable for the from-bam
    /// walker shim (`T = PileupRecord, E = WalkerError`) and any
    /// future seam, not coupled to `PreparedRead` / `CramInputError`.
    /// Uses primitive `T = i32, E = String` to decouple from any
    /// specific domain type.
    #[test]
    fn parametric_over_item_and_error_types() {
        let input: Vec<Result<i32, String>> = vec![Ok(1), Ok(2), Err("boom".to_string()), Ok(99)];
        let mut adapter = ErrorSheddingAdapter::new(input.into_iter());
        let handle = adapter.error_handle();

        assert_eq!(adapter.next(), Some(1));
        assert_eq!(adapter.next(), Some(2));
        assert_eq!(adapter.next(), None);
        assert_eq!(adapter.next(), None);

        let stashed = handle.take().expect("error must be stashed");
        assert_eq!(stashed, "boom");
    }
}
