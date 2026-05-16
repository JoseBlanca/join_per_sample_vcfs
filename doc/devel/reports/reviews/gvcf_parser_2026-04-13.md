# Code Review: gvcf_parser.rs and gvcf_parser_test.rs
**Date:** 2026-04-13  
**Reviewer:** GitHub Copilot  
**Module:** GVCF Parser  
**Status:** Request-changes

---

## 1. Scope
- **What was reviewed:** Review of production parser code and test suite for GVCF file parsing.
- **Reviewed against:** as-provided
- **In-scope files:**
  - [src/gvcf_parser.rs](../src/gvcf_parser.rs)
  - [tests/gvcf_parser_test.rs](../tests/gvcf_parser_test.rs)
- **Deliberately out of scope:** All other files in the crate.

---

## 2. Verdict
**Request-changes**

The parser has critical correctness issues in genotype parsing and FORMAT field access that require fixes before merge. Test coverage is insufficient to catch these defects.

---

## 3. Execution Status

| Command | Exit Code | Result |
|---------|-----------|--------|
| `cargo test --test gvcf_parser_test` | 0 | 17 tests passed; 1 warning (unused_mut at line 217) |
| `rustfmt --check` | 1 | Formatting issues in src/gvcf_parser.rs (lines 306–600) |
| `cargo clippy --test gvcf_parser_test -- -D warnings` | 101 | Lint failure at src/gvcf_parser.rs:719–722; other failures out of scope |
| `cargo doc --no-deps` | not run | Deferred; depends on overall project state |
| `cargo audit` | not run | Deferred; not critical for parser review |

- **Findings marked "Needs verification":** 0

---

## 4. Open Questions & Assumptions

1. **Ploidy contract:** The parser hardcodes ploidy=2 (diploid) at [src/gvcf_parser.rs](../src/gvcf_parser.rs#L65–L66) and [src/gvcf_parser.rs](../src/gvcf_parser.rs#L537–L538). 
We support different ploidies, so we should not hardcode any particular ploidy. Ploidy=2 could be a default, but the user of the API should be able to modify it.
The test suit must cover other possible ploidies.

---

## 5. Top 3 Priorities

1. **B1 – Allele index corruption:** Indices 128–255 wrap to negative i8 and silently become indistinguishable from missing alleles. Causes data loss in multiallelically rich VCFs. Fix: reject indices ≥ 128 with a typed error.

2. **B2 – Invalid GT character acceptance:** Malformed genotype strings (e.g., `0/a`, `0//1`, `0/1x`) are partially accepted because non-ASCII bytes are silently dropped. Fix: validate each byte, reject unexpected separators and empty tokens, and fail on any unparseable character. The error returned should include the genotype string that has caused the problem.

3. **M3 – Nominal genotype tests don't verify output:** `test_genotype_parsing_*` only assert variant count, not actual genotypes, phases, or FORMAT field values. Defects in B1 and B2 can pass a fully-green test suite. Fix: add assertions on `genotypes`, `phases`, and `get_gt_field_by_index()` return values.

---

## 6. Findings

### Blocker

#### B1: [src/gvcf_parser.rs](../src/gvcf_parser.rs#L152–L205) — Allele indices 128–255 silently corrupt genotypes

- **Confidence:** High  
- **Assumptions:** None.  
- **Problem:**  
  The `parse_genotype` function accumulates allele indices in a `u8` and casts to `i8` without validation:
  ```rust
  // Line 163–169
  current = (current as u8)
      .checked_mul(10)
      .and_then(|v| v.checked_add(byte - b'0'))
      .ok_or_else(|| VcfParseError::AlleleIndexOverflow { ... })? as i8;
  ```
  Indices 0–127 cast correctly, but 128–255 wrap to negative values: `128i8 = -128`, `255i8 = -1`. Since `-1` is the sentinel for missing alleles (`MISSING_ALLELE`), a multiallelically variant with allele index 255 becomes indistinguishable from a missing genotype.

- **Why it matters:**  
  This is silent data corruption. An input VCF with sufficiently many ALT alleles (≥128) produces wrong genotypes without error, panic, or warning. Downstream analyses will inherit invalid call data.

- **Suggested fix:**  
  Before casting, validate that `current ≤ 127`:
  ```rust
  current = (current as u8)
      .checked_mul(10)
      .and_then(|v| v.checked_add(byte - b'0'))
      .ok_or_else(|| VcfParseError::AlleleIndexOverflow { ... })?;
  if current > 127 {
      return Err(VcfParseError::AlleleIndexOutOfRange {
          gt: gt_str.to_string(),
          index: current as u32,
      });
  }
  current = current as i8;
  ```

---

#### B2: [src/gvcf_parser.rs](../src/gvcf_parser.rs#L161–L197) — Malformed GT strings accepted without validation

- **Confidence:** High  
- **Assumptions:** None.  
- **Problem:**  
  The byte-by-byte parser in `parse_genotype` silently ignores unrecognized characters:
  ```rust
  // Line 161–197
  for &byte in bytes {
      match byte {
          b'0'..=b'9' => { /* accumulate */ }
          b'.' => { /* handle missing */ }
          b'|' | b'/' => { /* handle separator */ }
          _ => {} // <-- silently accepts any other byte
      }
  }
  ```
  - Input `"0/a"`: the `a` is ignored, output is `[0]` (1 allele instead of 2).
  - Input `"0//1"`: consecutive `/` with the second separator processed as "no number", output is `[0, 1]` but phases may be confused.
  - Input `"0/1x"`: trailing `x` is silently dropped, output is `[0, 1]`.

  Additionally, the `seen_dot` flag on [line 191–193] attempts to handle `"./"` but doesn't enforce proper structure, so `"./0/.0"` parses without error.

- **Why it matters:**  
  A corrupted input file can produce plausibly valid-looking genotypes that mask the underlying defect. Bioinformatic pipelines downstream will silently process invalid data, leading to wrong results that are hard to trace back to the source.

- **Suggested fix:**  
  Rewrite the parser to enforce strict grammar:
  1. Each GT string is `<allele> ( ('/' | '|') <allele> )*`.
  2. Each `<allele>` is a (possibly multi-digit) decimal number or exactly `"."`.
  3. Any other character is an error, the error message should include the malformed genotype string to inform the user.
  4. Ploidy is validated at the call site (already done in `parse_genotypes` at line 284–292).

  Example replacement:
  ```rust
  fn parse_genotype(gt_str: &str, expected_ploidy: u8) -> VcfResult<(Vec<i8>, bool)> {
      let mut alleles = Vec::with_capacity(expected_ploidy as usize);
      let mut phased = false;
      let mut current_is_missing = false;
      let mut last_sep: Option<char> = None;
      let mut i = 0;
      let bytes = gt_str.as_bytes();
      let n = bytes.len();
      
      while i < n {
          let byte = bytes[i];
          match byte {
              b'.' => {
                  if i + 1 < n && bytes[i + 1].is_ascii_digit() {
                      return Err(VcfParseError::RuntimeError {
                          message: format!("Invalid GT: '.' must be standalone, got '{}'", gt_str),
                      });
                  }
                  alleles.push(MISSING_ALLELE);
                  current_is_missing = true;
                  i += 1;
              }
              b'0'..=b'9' => {
                  let start = i;
                  while i < n && bytes[i].is_ascii_digit() {
                      i += 1;
                  }
                  let num_str = std::str::from_utf8(&bytes[start..i])
                      .map_err(|_| VcfParseError::RuntimeError {
                          message: format!("Non-UTF8 in GT: {}", gt_str),
                      })?;
                  let idx: u8 = num_str.parse()
                      .map_err(|_| VcfParseError::RuntimeError {
                          message: format!("Allele index out of range: {}", num_str),
                      })?;
                  if idx > 127 {
                      return Err(VcfParseError::AlleleIndexOutOfRange {
                          gt: gt_str.to_string(),
                          index: idx as u32,
                      });
                  }
                  alleles.push(idx as i8);
                  current_is_missing = false;
              }
              b'/' | b'|' => {
                  if alleles.is_empty() {
                      return Err(VcfParseError::RuntimeError {
                          message: format!("Invalid GT: separator at start: {}", gt_str),
                      });
                  }
                  phased = byte == b'|';
                  if last_sep.is_some() && last_sep.unwrap() != (byte as char) {
                      return Err(VcfParseError::RuntimeError {
                          message: format!(
                              "Invalid GT: mixed separators (/ and |): {}",
                              gt_str
                          ),
                      });
                  }
                  last_sep = Some(byte as char);
                  i += 1;
              }
              _ => {
                  return Err(VcfParseError::RuntimeError {
                      message: format!(
                          "Invalid GT: unexpected character '{}' (byte {}) in: {}",
                          byte as char, byte, gt_str
                      ),
                  });
              }
          }
      }
      
      Ok((alleles, phased))
  }
  ```

---

### Major

#### M1: [src/gvcf_parser.rs](../src/gvcf_parser.rs#L399–L450) — Last sample's FORMAT field contains trailing newline

- **Confidence:** High  
- **Assumptions:** None.  
- **Problem:**  
  1. `read_line` preserves the trailing newline in the line buffer.
  2. `from_line` splits the line with `splitn(10, '\t')`, which does not trim.
  3. The 10th field (sample fields) is split by tab at [line 437] but not trimmed.
  4. The last sample's last field retains the newline: if input is `…\t15\n`, the DP value is `"15\n"`, not `"15"`.
  5. `get_gt_field_by_index` at [lines 496–500] returns raw split results without trim, so `get_gt_field_by_index(2)` on the last sample returns `"15\n"`, not `"15"`.

- **Why it matters:**  
  - Breaks the public API contract of `get_gt_field_by_index`. Callers expect clean strings but get values with newlines.
  - Downstream comparisons fail: `"15\n" != "15"`.
  - Serialization / logging produces malformed output.
  - Subtle and hard to debug because tests using `==` on counts or booleans won't catch it.

- **Suggested fix:**  
  Apply `trim_end()` once to the line after `read_line` and before `splitn(10, '\t')`:
  ```rust
  // In from_line, after reading the line:
  let mut fields = line.trim_end().splitn(10, '\t');
  ```
  Alternatively, trim the last field explicitly when constructing `sample_gt_fields`:
  ```rust
  let sample_fields = sample_fields.trim_end();
  ```
  Add a regression test:
  ```rust
  #[test]
  fn test_get_gt_field_last_sample_no_newline() {
      let vcf_data = "##fileformat=VCFv4.2\n\
                       #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n\
                       chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT:DP\t0/1:15\n";
      let reader = BufReader::new(vcf_data.as_bytes());
      let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
      let variants: Vec<_> = parser.filter_map(Result::ok).collect();
      assert_eq!(variants.len(), 1);
      
      let v = &variants[0];
      if let Some(dp_idx) = v.gt_field_index("DP") {
          let dp_values = v.get_gt_field_by_index(dp_idx);
          // Last sample's DP must be exactly "15", not "15\n"
          assert_eq!(dp_values[0], "15");
      }
  }
  ```

---

#### M2: [src/gvcf_parser.rs](../src/gvcf_parser.rs#L468–L477) — `get_span` may overflow without checked arithmetic

- **Confidence:** High  
- **Assumptions:** None.  
- **Problem:**  
  `get_span` computes `self.pos + max_allele_len as u32 - 1` without overflow checking. In contrast, `get_var_end` at [lines 329–335] uses `checked_add` and returns an explicit error on overflow. The inconsistency and the unchecked arithmetic are a bug waiting to happen.

  If `pos` is large (near `u32::MAX`) and `max_allele_len - 1` is also large, the addition wraps silently in release mode or panics in debug mode.

- **Why it matters:**  
  - In debug builds: panic on variant with `pos = u32::MAX - 10` and `max_allele_len = 20`.
  - In release builds: wrap to a small number, producing an incorrect span.
  - Error handling is inconsistent with the public API (callers expect `VcfResult`, not a panic-or-wrong-answer).

- **Suggested fix:**  
  Use the same defensive pattern as `get_var_end`:
  ```rust
  pub fn get_span(&self) -> VcfResult<(u32, u32)> {
      let max_allele_len = self.alleles.iter().map(|a| a.len()).max()
          .ok_or(VcfParseError::RuntimeError {
              message: "There should be at least one allele".to_string(),
          })?;
      if max_allele_len == 0 {
          return Err(VcfParseError::RuntimeError {
              message: format!("Empty allele at {}:{}", self.chrom, self.pos),
          });
      }
      let max_len_u32: u32 = max_allele_len.try_into()
          .map_err(|_| VcfParseError::RuntimeError {
              message: format!("Allele length overflow at {}:{}", self.chrom, self.pos),
          })?;
      if max_len_u32 == 1 {
          Ok((self.pos, self.pos))
      } else {
          let end = self.pos.checked_add(max_len_u32 - 1)
              .ok_or(VcfParseError::RuntimeError {
                  message: format!(
                      "Position overflow computing span at {}:{}",
                      self.chrom, self.pos
                  ),
              })?;
          Ok((self.pos, end))
      }
  }
  ```

---

#### M3: [tests/gvcf_parser_test.rs](../tests/gvcf_parser_test.rs#L247–L315) — Genotype tests only verify variant count, not genotype correctness

- **Confidence:** High  
- **Assumptions:** None.  
- **Problem:**  
  Tests `test_genotype_parsing_basic` (line 247–258), `test_genotype_parsing_missing` (line 261–271), `test_genotype_parsing_multiallelic` (line 274–284), and `test_genotype_with_format_fields` (line 287–297) all:
  1. Create a VCF with known genotype fields.
  2. Parse it.
  3. Assert `variants.len()` or `n_samples`.
  4. Do **not** check `genotypes`, `phases`, or `get_gt_field_by_index()` values.

  This means:
  - Defect B1 (corruption of indices ≥128) is invisible in the test suite.
  - Defect B2 (silent acceptance of malformed GT) is invisible.
  - Defect M1 (newline in FORMAT fields) is invisible.
  - A test passes as long as **some** output is produced, regardless of correctness.

- **Why it matters:**  
  Test names imply coverage that does not exist. A reviewer trusting the test suite may accept regressions introduced by future changes. The current implementation can have silent data corruption while passing all tests.

- **Suggested fix:**  
  Replace each test with assertions on actual genotype values. Example rewrite:
  ```rust
  #[test]
  fn test_genotype_parsing_basic() {
      let vcf_data = "##fileformat=VCFv4.2\n\
                       #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3\n\
                       chr1\t100\t.\tA\tC\t30\tPASS\t.\tGT\t0/0\t0/1\t1/1\n\
                       chr1\t200\t.\tG\tT\t40\tPASS\t.\tGT\t0|0\t0|1\t1|1";
      
      let reader = BufReader::new(vcf_data.as_bytes());
      let parser = VarIterator::from_reader(reader).expect("Failed to create parser");
      let variants: Vec<_> = parser.filter_map(Result::ok).collect();
      
      assert_eq!(variants.len(), 2);
      
      // Variant 1: chr1:100, unphased
      let v1 = &variants[0];
      assert_eq!(v1.pos, 100);
      // Ploidy = 2, 3 samples => genotypes.len() == 6
      assert_eq!(v1.genotypes.len(), 6);
      // Sample 1: 0/0 -> [0, 0]
      assert_eq!(v1.genotypes[0], 0);
      assert_eq!(v1.genotypes[1], 0);
      assert_eq!(v1.phases[0], false);
      // Sample 2: 0/1 -> [0, 1]
      assert_eq!(v1.genotypes[2], 0);
      assert_eq!(v1.genotypes[3], 1);
      assert_eq!(v1.phases[1], false);
      // Sample 3: 1/1 -> [1, 1]
      assert_eq!(v1.genotypes[4], 1);
      assert_eq!(v1.genotypes[5], 1);
      assert_eq!(v1.phases[2], false);
      
      // Variant 2: chr1:200, phased
      let v2 = &variants[1];
      assert_eq!(v2.pos, 200);
      // Sample 1: 0|0 -> [0, 0] phased
      assert_eq!(v2.genotypes[0], 0);
      assert_eq!(v2.genotypes[1], 0);
      assert_eq!(v2.phases[0], true);
      // Sample 2: 0|1 -> [0, 1] phased
      assert_eq!(v2.genotypes[2], 0);
      assert_eq!(v2.genotypes[3], 1);
      assert_eq!(v2.phases[1], true);
      // Sample 3: 1|1 -> [1, 1] phased
      assert_eq!(v2.genotypes[4], 1);
      assert_eq!(v2.genotypes[5], 1);
      assert_eq!(v2.phases[2], true);
  }
  ```

  Similar rewrites for `test_genotype_parsing_missing`, `test_genotype_parsing_multiallelic`, and `test_genotype_with_format_fields`.

---

### Nits

- **[tests/gvcf_parser_test.rs](../tests/gvcf_parser_test.rs#L217):** Unnecessary `mut` keyword on `parser`. Use `cargo fix --test gvcf_parser_test` to remove it.
- **[src/gvcf_parser.rs](../src/gvcf_parser.rs#L719–L722):** Use `if n_items_added == 0 { return None; }` instead of `match` to simplify.
- **[src/gvcf_parser.rs](../src/gvcf_parser.rs):** Does not pass `rustfmt --check`. Run `cargo fmt` to reformat lines 306–600.

---

## 7. Out of Scope Observations

The following lints were raised by `cargo clippy --test gvcf_parser_test` but fall outside GVCF parser scope:
- [src/decompression_pool.rs](../src/decompression_pool.rs#L133): `&mut Vec` should be `&mut [u8]`
- [src/genotype_merging.rs](../src/genotype_merging.rs#L176): Needless range loop; use `enumerate()`
- [src/genotype_merging.rs](../src/genotype_merging.rs#L212): Use `repeat_n()` instead of `repeat().take()`
- [src/variant_grouping.rs](../src/variant_grouping.rs#L88, L281, L377): Collapsible if statements

---

## 8. Missing Tests to Add Now

### For `parse_genotype`:
- **`parse_rejects_invalid_characters`**: Input `"0/a"`, `"0//1"`, `"0/1x"`. Expected: error, not partial acceptance.
- **`parse_rejects_allele_index_out_of_range`**: Input `"128/129"` (or higher to exceed i8::MAX). Expected: error, not silent wrap to negative.
- **`parse_rejects_mixed_separators`**: Input `"0/1|2"`. Expected: error, not mixed phasing.
- **`parse_rejects_leading_separator`**: Input `"/0/1"`. Expected: error.

### For `Variant::get_gt_field_by_index`:
- **`get_gt_field_by_index_trims_newline_last_sample`**: Parse a VCF with `GT:DP` on the last sample. Call `get_gt_field_by_index(1)` and assert the last element is exactly `"15"`, not `"15\n"`.

### For `Variant::get_span`:
- **`get_span_returns_error_on_position_overflow`**: Create a variant with `pos = u32::MAX - 5` and `max_allele_len = 10`. Expected: `Err(VcfParseError::RuntimeError)`, not wrap or panic.

### For `parse_genotypes` (integration):
- **`parse_genotypes_with_missing_and_phased`**: Parse `"./.  |  0|1"`. Expected: `genotypes = [-1, -1, 0, 1]`, `phases = [false, true]`.
- **`parse_genotypes_validates_sample_count`**: Parse with 2 samples declared but 3 in the line. Expected: `MalformedLine` error.

---

## 9. What's Good

- **[src/gvcf_parser.rs](../src/gvcf_parser.rs#L308–L336)** `get_var_end` is a model of defensive programming: checks for empty alleles, validates allele length, and uses `checked_add` for overflow. Apply this pattern to `get_span` and `parse_genotype`.

- **[src/gvcf_parser.rs](../src/gvcf_parser.rs#L426–L434)** FORMAT field caching with a small LRU cache is a sensible optimization for a hot path and is correctly scoped.

- **[tests/gvcf_parser_test.rs](../tests/gvcf_parser_test.rs#L35–L67)** `test_peek_variant` and `test_peek_variant_exhausted` thoroughly exercise the iterator's semantics. This is the kind of precision coverage needed for the genotype tests.

---

## 10. Commands to Re-verify

After fixes, re-run:
```bash
cargo test --test gvcf_parser_test
rustfmt --check src/gvcf_parser.rs tests/gvcf_parser_test.rs
cargo clippy --test gvcf_parser_test -- -D warnings
```

Confirm all 17+ tests pass (including new tests for edge cases), formatting clean, and no clippy warnings.

---

## Author Response Template

Once the author addresses findings, respond with:

### B1 – Allele index corruption
`fixed in <commit>` / `disputed because …` / `deferred to <issue>`

### B2 – Malformed GT acceptance
`fixed in <commit>` / `disputed because …` / `deferred to <issue>`

### M1 – Newline in FORMAT field
`fixed in <commit>` / `disputed because …` / `deferred to <issue>`

### M2 – `get_span` overflow
`fixed in <commit>` / `disputed because …` / `deferred to <issue>`

### M3 – Genotype test coverage
`fixed in <commit>` / `disputed because …` / `deferred to <issue>`

### Open Question: Ploidy contract
Answer: `Only diploid supported` / `Polyploidy support planned; current docs should be updated` / other

