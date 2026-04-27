//! A peek+buffer adapter for fallible iterators.
//!
//! See `ia/specs/buffered_peekable.md` for the full specification.

use std::collections::VecDeque;

pub const DEFAULT_BUFFER_SIZE: usize = 100;

/// Wraps an `Iterator<Item = Result<T, E>>` to add a transposed `peek`
/// returning `Result<Option<&T>, E>` and a configurable look-ahead
/// buffer that pulls from the inner iterator in batched bursts.
///
/// Invariant: `buffer` holds the prefix of the inner iterator's output
/// that has been pulled but not yet consumed by the caller, in source
/// order. Both `Ok` and `Err` items are stored verbatim.
pub struct BufferedPeekable<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    inner: I,
    buffer: VecDeque<Result<T, E>>,
    buffer_size: usize,
}

impl<I, T, E> BufferedPeekable<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    /// Wrap `inner` with the default buffer size.
    pub fn new(inner: I) -> Self {
        Self::with_buffer_size(inner, DEFAULT_BUFFER_SIZE)
    }

    /// Wrap `inner` with an explicit `buffer_size`.
    ///
    /// Panics if `buffer_size` is zero.
    pub fn with_buffer_size(inner: I, buffer_size: usize) -> Self {
        assert!(buffer_size >= 1, "buffer_size must be at least 1");
        Self {
            inner,
            buffer: VecDeque::with_capacity(buffer_size),
            buffer_size,
        }
    }

    /// Look at the next item without consuming it.
    ///
    /// Returns `Ok(Some(&t))` when the head of the stream is a record,
    /// `Ok(None)` when the stream is exhausted, or `Err(e)` when the
    /// next event is an error. Errors are taken out of the buffer when
    /// surfaced (a peek that returns `Err` advances past the error);
    /// records remain in the buffer until consumed via `next`.
    pub fn peek(&mut self) -> Result<Option<&T>, E> {
        if self.buffer.is_empty() {
            self.refill();
        }
        if matches!(self.buffer.front(), Some(Err(_))) {
            let Some(Err(e)) = self.buffer.pop_front() else {
                unreachable!("front is Err but pop_front returned otherwise")
            };
            return Err(e);
        }
        Ok(match self.buffer.front() {
            Some(Ok(t)) => Some(t),
            None => None,
            Some(Err(_)) => unreachable!("Err case handled above"),
        })
    }

    fn refill(&mut self) {
        while self.buffer.len() < self.buffer_size {
            match self.inner.next() {
                Some(item) => self.buffer.push_back(item),
                None => break,
            }
        }
    }
}

impl<I, T, E> Iterator for BufferedPeekable<I, T, E>
where
    I: Iterator<Item = Result<T, E>>,
{
    type Item = Result<T, E>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer.is_empty() {
            self.refill();
        }
        self.buffer.pop_front()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_peekable(
        items: Vec<Result<i32, String>>,
    ) -> BufferedPeekable<std::vec::IntoIter<Result<i32, String>>, i32, String> {
        BufferedPeekable::new(items.into_iter())
    }

    fn make_peekable_sized(
        items: Vec<Result<i32, String>>,
        buffer_size: usize,
    ) -> BufferedPeekable<std::vec::IntoIter<Result<i32, String>>, i32, String> {
        BufferedPeekable::with_buffer_size(items.into_iter(), buffer_size)
    }

    #[test]
    fn empty_iterator_yields_none_forever() {
        let mut p = make_peekable(vec![]);
        assert!(matches!(p.peek(), Ok(None)));
        assert!(matches!(p.peek(), Ok(None)));
        assert!(p.next().is_none());
        assert!(p.next().is_none());
        assert!(matches!(p.peek(), Ok(None)));
    }

    #[test]
    fn all_ok_in_order() {
        let mut p = make_peekable(vec![Ok(1), Ok(2), Ok(3)]);
        assert_eq!(p.next().unwrap().unwrap(), 1);
        assert_eq!(p.next().unwrap().unwrap(), 2);
        assert_eq!(p.next().unwrap().unwrap(), 3);
        assert!(p.next().is_none());
    }

    #[test]
    fn peek_does_not_consume_record() {
        let mut p = make_peekable(vec![Ok(7), Ok(8)]);
        assert_eq!(*p.peek().unwrap().unwrap(), 7);
        assert_eq!(*p.peek().unwrap().unwrap(), 7);
        assert_eq!(*p.peek().unwrap().unwrap(), 7);
        // Record still there for next.
        assert_eq!(p.next().unwrap().unwrap(), 7);
        assert_eq!(*p.peek().unwrap().unwrap(), 8);
    }

    #[test]
    fn peek_then_next_returns_same_record() {
        let mut p = make_peekable(vec![Ok(42)]);
        let peeked = *p.peek().unwrap().unwrap();
        let nexted = p.next().unwrap().unwrap();
        assert_eq!(peeked, nexted);
    }

    #[test]
    fn single_error_preserves_order() {
        let mut p = make_peekable(vec![Ok(1), Ok(2), Err("boom".into()), Ok(3)]);
        assert_eq!(p.next().unwrap().unwrap(), 1);
        assert_eq!(p.next().unwrap().unwrap(), 2);
        // Next event in the stream is the error.
        let err = p.next().unwrap().unwrap_err();
        assert_eq!(err, "boom");
        // After the error, the stream continues.
        assert_eq!(p.next().unwrap().unwrap(), 3);
        assert!(p.next().is_none());
    }

    #[test]
    fn error_first_then_recovers() {
        let mut p = make_peekable(vec![Err("first".into()), Ok(99)]);
        // Peek surfaces the error and advances past it.
        let err = p.peek().unwrap_err();
        assert_eq!(err, "first");
        // Next event is the Ok.
        assert_eq!(*p.peek().unwrap().unwrap(), 99);
        assert_eq!(p.next().unwrap().unwrap(), 99);
        assert!(p.next().is_none());
    }

    #[test]
    fn multiple_errors_each_surface_in_turn() {
        let mut p = make_peekable(vec![Err("a".into()), Err("b".into()), Ok(5)]);
        assert_eq!(p.peek().unwrap_err(), "a");
        assert_eq!(p.peek().unwrap_err(), "b");
        assert_eq!(*p.peek().unwrap().unwrap(), 5);
        assert_eq!(p.next().unwrap().unwrap(), 5);
        assert!(p.next().is_none());
    }

    #[test]
    fn peek_consumes_error_only_once() {
        let mut p = make_peekable(vec![Err("once".into()), Ok(1)]);
        // First peek surfaces the error.
        assert_eq!(p.peek().unwrap_err(), "once");
        // Second peek must not return the same error — buffer has moved on.
        assert_eq!(*p.peek().unwrap().unwrap(), 1);
    }

    #[test]
    fn next_surfaces_error_when_called_directly() {
        let mut p = make_peekable(vec![Ok(1), Err("mid".into()), Ok(2)]);
        assert_eq!(p.next().unwrap().unwrap(), 1);
        assert_eq!(p.next().unwrap().unwrap_err(), "mid");
        assert_eq!(p.next().unwrap().unwrap(), 2);
    }

    #[test]
    fn buffer_size_one_works() {
        let mut p = make_peekable_sized(vec![Ok(1), Ok(2), Ok(3)], 1);
        assert_eq!(*p.peek().unwrap().unwrap(), 1);
        assert_eq!(p.next().unwrap().unwrap(), 1);
        assert_eq!(*p.peek().unwrap().unwrap(), 2);
        assert_eq!(p.next().unwrap().unwrap(), 2);
        assert_eq!(p.next().unwrap().unwrap(), 3);
        assert!(p.next().is_none());
    }

    #[test]
    fn buffer_size_larger_than_stream() {
        let mut p = make_peekable_sized(vec![Ok(1), Ok(2)], 1000);
        assert_eq!(p.next().unwrap().unwrap(), 1);
        assert_eq!(p.next().unwrap().unwrap(), 2);
        assert!(p.next().is_none());
    }

    #[test]
    #[should_panic(expected = "buffer_size must be at least 1")]
    fn with_buffer_size_zero_panics() {
        let _ = make_peekable_sized(vec![Ok(1)], 0);
    }

    #[test]
    fn next_after_exhaustion_keeps_returning_none() {
        let mut p = make_peekable(vec![Ok(1)]);
        assert_eq!(p.next().unwrap().unwrap(), 1);
        assert!(p.next().is_none());
        assert!(p.next().is_none());
        assert!(p.next().is_none());
        assert!(matches!(p.peek(), Ok(None)));
    }
}
