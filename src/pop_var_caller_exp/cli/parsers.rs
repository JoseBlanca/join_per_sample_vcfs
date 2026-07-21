//! Value parsers for the `pop_var_caller_exp` CLI. Production's
//! [`crate::pop_var_caller::cli::parsers`] is the *pattern* this copies, not a
//! place to add to (spec §2.1) — ng's experiment knobs stay out of the
//! production CLI.

use crate::ng::region_typing::segment_criteria::{MAX_MOTIF_LEN, MinCopies};

/// The floor [`MinCopies`] applies to periods wider than the table, which
/// `--min-copies` deliberately does not expose.
///
/// It is **structurally unreachable** from this command: the walk asserts the
/// period ceiling is `<= MAX_MOTIF_LEN`
/// ([`segment_criteria`](crate::ng::region_typing::segment_criteria)), so no
/// period-7+ interval can ever reach the floor. Rather than make the user type a
/// value that decides nothing, the parser supplies this inert one (spec §2.1).
const INERT_WIDER_PERIOD_FLOOR: u32 = 3;

/// The arity is [`MAX_MOTIF_LEN`], but the flag's `default_value` string, its
/// help text and this module's prose all spell "six" / `1..=6` literally. If the
/// motif ceiling ever moved, those literals would be wrong *and* every run
/// without an explicit `--min-copies` would fail at parse time (clap resolves
/// `default_value` through this parser). Fail at compile time instead.
const _: () = assert!(
    MAX_MOTIF_LEN == 6,
    "--min-copies' default string and help text are written for six periods"
);

/// Parse `--min-copies` — the per-period copy-number floor table.
///
/// **Exactly [`MAX_MOTIF_LEN`] comma-separated values, one per period 1..=6**
/// (e.g. `6,4,4,3,3,3`); any other count is a hard parse error, surfaced by clap
/// as a usage failure (exit 2) before `run_typed_regions` is reached. No
/// catch-all and no colon syntax: the unreachable `for_wider_periods` slot is
/// filled with [`INERT_WIDER_PERIOD_FLOOR`] (spec §2.1, arch §2).
pub fn parse_min_copies(s: &str) -> Result<MinCopies, String> {
    let values = s
        .split(',')
        .map(|field| {
            let field = field.trim();
            field.parse::<u32>().map_err(|_| {
                format!("'{field}' is not a non-negative integer (--min-copies is a comma-separated list)")
            })
        })
        .collect::<Result<Vec<u32>, String>>()?;

    let by_period: [u32; MAX_MOTIF_LEN] = values.as_slice().try_into().map_err(|_| {
        format!(
            "--min-copies needs exactly {MAX_MOTIF_LEN} comma-separated values, one per period \
             1..={MAX_MOTIF_LEN} (e.g. 6,4,4,3,3,3) — got {}",
            values.len()
        )
    })?;

    Ok(MinCopies::new(by_period, INERT_WIDER_PERIOD_FLOOR))
}

/// Parse `--min-purity` — a **finite** fraction in `[0, 1]`.
///
/// clap can range-check integers but not floats, and the walk guards this knob
/// with a release `assert!` rather than an error (the knob is swept, and sweeps
/// run in `--release`). An assert is right for a library invariant and wrong for
/// a flag, so the bound is enforced here and a typo is a usage failure rather
/// than a backtrace (spec §6).
///
/// **`nan` is the case this exists for.** It parses happily as an `f32`, and a
/// NaN floor compares false against everything — so rather than rejecting
/// impure tracts it would pass *every* tract, silently inverting the gate.
pub fn parse_min_purity(s: &str) -> Result<f32, String> {
    let value: f32 = s
        .trim()
        .parse()
        .map_err(|_| format!("'{s}' is not a number"))?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(format!(
            "--min-purity must be a finite fraction in [0, 1] (e.g. 0.8), got {value}"
        ));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Six values become the table, read back per period.
    #[test]
    fn six_values_parse_to_the_right_table() {
        let table = parse_min_copies("6,4,4,3,3,3").expect("six values parse");
        let floors: Vec<u32> = (1..=6).map(|p| table.for_period(p)).collect();
        assert_eq!(floors, vec![6, 4, 4, 3, 3, 3]);
    }

    /// The values are read positionally, one per period — a table that is not
    /// uniform proves the order is not lost (a `uniform(n)` fixture could not).
    #[test]
    fn values_map_to_periods_in_order() {
        let table = parse_min_copies("11,12,13,14,15,16").expect("six values parse");
        for (i, expected) in [11, 12, 13, 14, 15, 16].into_iter().enumerate() {
            let period = u8::try_from(i + 1).expect("period fits");
            assert_eq!(table.for_period(period), expected, "period {period}");
        }
    }

    /// The unreachable wider-period slot is supplied inert, not exposed.
    ///
    /// **Deliberately a non-uniform table.** At `6,4,4,3,3,3` the period-6 floor
    /// is itself `3`, so a regression that *clamped* `for_period(7)` to the last
    /// table slot instead of falling through to `for_wider_periods` would still
    /// return `3` and this test would stay green. `11..16` tells the two apart.
    #[test]
    fn the_wider_period_slot_is_filled_inert() {
        let table = parse_min_copies("11,12,13,14,15,16").expect("six values parse");
        assert_eq!(
            table.for_period(7),
            INERT_WIDER_PERIOD_FLOOR,
            "period 7+ is structurally unreachable; the parser fills it rather than \
             making the user type a value that decides nothing — and it falls THROUGH \
             the table (16 would mean it clamped to the last slot)"
        );
    }

    /// **The parsed short-read default must equal the library's**
    /// [`MinCopies::default`]. The flag's `default_value` is a literal string (a
    /// `MinCopies` has no `Display` to derive `default_value_t` from), so nothing
    /// but this ties the two together: sweep the library floors (spec §10) without
    /// updating the flag and the CLI would silently apply the old table.
    #[test]
    fn the_default_string_parses_to_the_library_default() {
        assert_eq!(
            parse_min_copies("6,4,4,3,3,3").expect("the flag's default parses"),
            MinCopies::default(),
            "--min-copies' default_value string has drifted from MinCopies::default()"
        );
    }

    /// Too few values is a hard error, and the message says what was wanted.
    #[test]
    fn five_values_are_rejected() {
        let err = parse_min_copies("6,4,4,3,3").expect_err("five values must fail");
        assert!(err.contains("exactly 6"), "got {err}");
        assert!(err.contains("got 5"), "the message names the count: {err}");
    }

    /// Too many values is a hard error too — no silent truncation.
    #[test]
    fn seven_values_are_rejected() {
        let err = parse_min_copies("6,4,4,3,3,3,3").expect_err("seven values must fail");
        assert!(err.contains("exactly 6"), "got {err}");
        assert!(err.contains("got 7"), "the message names the count: {err}");
    }

    /// A non-numeric field is rejected, naming the offending field.
    #[test]
    fn a_non_numeric_value_is_rejected() {
        let err = parse_min_copies("6,4,x,3,3,3").expect_err("non-numeric must fail");
        assert!(err.contains('x'), "the message names the bad field: {err}");
    }

    /// A negative value is rejected — the floor is a count.
    #[test]
    fn a_negative_value_is_rejected() {
        let err = parse_min_copies("6,4,-1,3,3,3").expect_err("negative must fail");
        assert!(err.contains("-1"), "the message names the bad field: {err}");
    }

    /// Surrounding whitespace is tolerated, so `--min-copies "6, 4, 4, 3, 3, 3"`
    /// works as typed.
    #[test]
    fn whitespace_around_values_is_tolerated() {
        let table = parse_min_copies(" 6, 4 ,4,3,3,3 ").expect("whitespace is trimmed");
        assert_eq!(table.for_period(1), 6);
        assert_eq!(table.for_period(2), 4);
    }

    /// An empty string is rejected — as an empty *field*, not a count error
    /// (`"".split(',')` yields one empty field, which fails to parse as a number).
    /// Asserted precisely so the test cannot drift from the message it documents.
    #[test]
    fn an_empty_string_is_rejected() {
        let err = parse_min_copies("").expect_err("an empty string must fail");
        assert!(
            err.contains("not a non-negative integer"),
            "an empty string fails as an empty field, not on the count: {err}"
        );
    }

    /// A trailing comma leaves an empty field, and is rejected for that reason —
    /// never silently treated as five values or padded to six.
    #[test]
    fn a_trailing_comma_is_rejected() {
        let err = parse_min_copies("6,4,4,3,3,3,").expect_err("a trailing comma must fail");
        assert!(err.contains("not a non-negative integer"), "got {err}");
    }
}
