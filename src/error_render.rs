//! Shared CLI error rendering. Both binaries (`pop_var_caller` and
//! `pop_var_caller_exp`) report an error the same way, so the rule lives
//! here in the library rather than being copied into each `main` — two
//! copies of an error-rendering heuristic that must agree is how they stop
//! agreeing (spec T7a).

use std::error::Error;

/// Walk the [`Error::source`] chain joining messages with `: ` so the user
/// sees the leaf cause, not just the outermost wrapper.
///
/// Several of our errors (and noodles') use the `thiserror` `: {source}`
/// pattern that already embeds the child's message in the parent's
/// Display; skip a level when its message is already a substring of the
/// previous one to avoid `invalid record: invalid record: invalid value`
/// noise.
pub fn format_error_chain(err: &(dyn Error + 'static)) -> String {
    let mut out = err.to_string();
    let mut prev = out.clone();
    let mut cur = err.source();
    while let Some(e) = cur {
        let msg = e.to_string();
        if !prev.contains(&msg) {
            out.push_str(": ");
            out.push_str(&msg);
        }
        prev = msg;
        cur = e.source();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::format_error_chain;
    use std::error::Error;
    use std::fmt;

    #[derive(Debug)]
    struct TestErr {
        msg: &'static str,
        src: Option<Box<TestErr>>,
    }

    impl fmt::Display for TestErr {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.write_str(self.msg)
        }
    }

    impl Error for TestErr {
        fn source(&self) -> Option<&(dyn Error + 'static)> {
            self.src.as_deref().map(|e| e as &(dyn Error + 'static))
        }
    }

    fn err(msg: &'static str, src: Option<TestErr>) -> TestErr {
        TestErr {
            msg,
            src: src.map(Box::new),
        }
    }

    #[test]
    fn joins_distinct_messages_with_colon() {
        let leaf = err("duplicate tag: PL", None);
        let mid = err("invalid read group", Some(leaf));
        let top = err("CRAM input", Some(mid));
        assert_eq!(
            format_error_chain(&top),
            "CRAM input: invalid read group: duplicate tag: PL"
        );
    }

    #[test]
    fn skips_level_already_embedded_in_parent() {
        // thiserror's `: {source}` pattern bakes the child's Display
        // into the parent's. Walking the chain naively would then
        // re-print the child — the dedup must drop it.
        let leaf = err("inner cause", None);
        let outer = err("outer prefix: inner cause", Some(leaf));
        assert_eq!(format_error_chain(&outer), "outer prefix: inner cause");
    }

    #[test]
    fn surfaces_deeper_message_when_intermediate_was_a_duplicate() {
        // Mirrors the real noodles chain: an outer error whose
        // Display embeds a generic middle ("invalid record") while
        // the actual useful detail lives one level deeper.
        let leaf = err("duplicate tag: PL", None);
        let mid = err("invalid record", Some(leaf));
        let top = err("failed to open CRAM 'x.cram': invalid record", Some(mid));
        assert_eq!(
            format_error_chain(&top),
            "failed to open CRAM 'x.cram': invalid record: duplicate tag: PL"
        );
    }
}
