//! Contamination-estimates artefact (TOML).
//!
//! On-disk format produced by the `estimate-contamination` subcommand
//! and consumed by `var-calling`. Single file, normalised (no
//! redundancy between sample rows and per-batch contaminant
//! distributions). Schema and rationale are pinned in the cohort CLI
//! plan §"Contamination artefact — TOML schema".
//!
//! ```toml
//! [provenance]
//! tool       = "pop_var_caller"
//! version    = "0.1.0"
//! subcommand = "estimate-contamination"
//! created    = 2026-05-19T14:32:11Z
//!
//! [provenance.inputs]
//! reference        = "grch38.fa"
//! input_psps       = ["NA12878.psp", "NA12891.psp"]
//! batch_assignment = "lanes.tsv"        # absent when omitted
//!
//! [parameters]
//! stopping_mode = "convergence"
//! block_size    = 1000
//! # ... every effective knob
//!
//! [[batches]]
//! id                         = "lane_3"
//! contaminant_ref_prob       = 0.9989
//! contaminant_snp_alt_prob   = 0.0010
//! contaminant_indel_alt_prob = 0.0001
//!
//! [[samples]]
//! name                   = "NA12878"
//! batch                  = "lane_3"
//! contamination_fraction = 0.0123
//! ```
//!
//! ## Validation
//!
//! At load time (and re-checked at write time):
//!
//! - `batches[*].id` is unique within the file.
//! - Every `samples[*].batch` resolves into an existing
//!   `batches[*].id`.
//! - Every probability (`contamination_fraction`,
//!   `contaminant_*_prob`) is finite and in `[0.0, 1.0]`.
//! - For each batch the three `contaminant_*_prob` values sum to
//!   `1.0` within `1e-9`, **or** every field is exactly `0.0` — the
//!   explicit all-zero row is the engine's signal that the batch was
//!   floored (singleton or below `min_batch_size`) and its `q_b` is
//!   unused. Compound alleles are implicitly 0 (Assumption 2 of the
//!   contamination impl report).
//! - `samples[*].name` is unique within the file.
//!
//! Sample-name reconciliation against the cohort `.psp` inputs lives
//! in [`Self::to_estimates_for_samples`] — extras in the artefact
//! are silently ignored; missing-from-artefact `.psp` samples are a
//! hard error.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::var_calling::contamination_estimation::{
    ContaminationEstimates, ContaminationEstimationError, N_ALLELE_CLASSES,
};

/// Probability simplex check tolerance for each batch's
/// `contaminant_*_prob` row. Matches the tolerance the engine-side
/// [`ContaminationEstimates::from_user_supplied`] uses (`1e-6`)
/// loosened by an additional safety factor — round-trip via
/// `toml::Value::Float` (f64) can shave a few ULP off, so we keep some
/// slack on top of the engine bound.
pub const Q_B_SIMPLEX_TOLERANCE: f64 = 1e-9;

// ---------------------------------------------------------------------
// On-disk types — serde-derived. Public so consumers can build /
// inspect, but every load / write goes through the validator.
// ---------------------------------------------------------------------

/// Top-level artefact. Construct with the public fields, then write
/// with [`Self::write`]; load with [`Self::read`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContaminationArtefact {
    /// Tool / version / timestamp / input identification.
    pub provenance: Provenance,
    /// Every effective parameter from the producing run, recorded for
    /// reproducibility. `toml::Value` so any TOML-compatible value
    /// type round-trips.
    pub parameters: BTreeMap<String, toml::Value>,
    /// Per-batch contaminant background distributions.
    pub batches: Vec<BatchEntry>,
    /// Per-sample contamination fractions, with batch labels that
    /// resolve into [`Self::batches`].
    pub samples: Vec<SampleEntry>,
}

/// Producer-tool provenance — mirrors the .psp `WriterProvenance`
/// shape.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Provenance {
    pub tool: String,
    pub version: String,
    pub subcommand: String,
    /// RFC 3339 timestamp; round-trips as `toml::Datetime`.
    pub created: toml::value::Datetime,
    pub inputs: ProvenanceInputs,
}

/// Identifying paths from the producing run. Path values are
/// **basenames only**, matching the `.psp` provenance convention.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProvenanceInputs {
    pub reference: String,
    pub input_psps: Vec<String>,
    /// `None` when `--batch-assignment` was omitted (all-samples-one-
    /// batch default).
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub batch_assignment: Option<String>,
}

/// One row of `[[batches]]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BatchEntry {
    pub id: String,
    pub contaminant_ref_prob: f64,
    pub contaminant_snp_alt_prob: f64,
    pub contaminant_indel_alt_prob: f64,
}

/// One row of `[[samples]]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SampleEntry {
    pub name: String,
    pub batch: String,
    pub contamination_fraction: f64,
}

// ---------------------------------------------------------------------
// I/O + validation
// ---------------------------------------------------------------------

impl ContaminationArtefact {
    /// Load and validate the artefact at `path`. Returns the parsed
    /// artefact with every invariant from the module docs guaranteed.
    pub fn read(path: &Path) -> Result<Self, ContaminationArtefactError> {
        let contents =
            fs::read_to_string(path).map_err(|source| ContaminationArtefactError::Io {
                path: path.to_path_buf(),
                source,
            })?;
        let me: Self =
            toml::from_str(&contents).map_err(|source| ContaminationArtefactError::Toml {
                path: path.to_path_buf(),
                source,
            })?;
        me.validate()?;
        Ok(me)
    }

    /// Validate then atomically write the artefact to `path` via
    /// `<path>.tmp` → `fs::rename`. The atomic-rename idiom matches
    /// the rest of the pipeline's output discipline (see
    /// `run_pileup`); a killed process never leaves a half-written
    /// artefact at the final path.
    pub fn write(&self, path: &Path) -> Result<(), ContaminationArtefactError> {
        self.validate()?;
        let body = toml::to_string_pretty(self).map_err(|source| {
            ContaminationArtefactError::Serialize {
                path: path.to_path_buf(),
                source,
            }
        })?;
        let mut tmp_os = path.as_os_str().to_os_string();
        tmp_os.push(".tmp");
        let tmp = PathBuf::from(tmp_os);
        fs::write(&tmp, &body).map_err(|source| ContaminationArtefactError::Io {
            path: tmp.clone(),
            source,
        })?;
        fs::rename(&tmp, path).map_err(|source| ContaminationArtefactError::Io {
            path: path.to_path_buf(),
            source,
        })?;
        Ok(())
    }

    /// Run every validation rule documented at the module level.
    /// Called automatically by [`Self::read`] and [`Self::write`];
    /// public so callers that build an artefact in-memory can re-check
    /// after edits.
    pub fn validate(&self) -> Result<(), ContaminationArtefactError> {
        // Per-batch checks: id uniqueness + probability range + simplex.
        let mut batch_ids: BTreeMap<&str, usize> = BTreeMap::new();
        for (idx, b) in self.batches.iter().enumerate() {
            if batch_ids.insert(b.id.as_str(), idx).is_some() {
                return Err(ContaminationArtefactError::DuplicateBatchId { id: b.id.clone() });
            }
            for (name, value) in [
                ("contaminant_ref_prob", b.contaminant_ref_prob),
                ("contaminant_snp_alt_prob", b.contaminant_snp_alt_prob),
                ("contaminant_indel_alt_prob", b.contaminant_indel_alt_prob),
            ] {
                if !(value.is_finite() && (0.0..=1.0).contains(&value)) {
                    return Err(ContaminationArtefactError::ProbabilityOutOfRange {
                        batch_id: b.id.clone(),
                        field: name,
                        got: value,
                    });
                }
            }
            // An explicit all-zero row signals "batch floored, q_b
            // unused" — engine convention (mirrors
            // ContaminationEstimates::from_user_supplied). Anything
            // else must lie on the probability simplex.
            let all_zero = b.contaminant_ref_prob == 0.0
                && b.contaminant_snp_alt_prob == 0.0
                && b.contaminant_indel_alt_prob == 0.0;
            let sum =
                b.contaminant_ref_prob + b.contaminant_snp_alt_prob + b.contaminant_indel_alt_prob;
            if !all_zero && (sum - 1.0).abs() > Q_B_SIMPLEX_TOLERANCE {
                return Err(ContaminationArtefactError::BatchProbabilitiesDoNotSum {
                    batch_id: b.id.clone(),
                    sum,
                });
            }
        }

        // Per-sample checks: name uniqueness + batch resolves +
        // contamination_fraction range.
        let mut sample_names: BTreeMap<&str, ()> = BTreeMap::new();
        for s in &self.samples {
            if sample_names.insert(s.name.as_str(), ()).is_some() {
                return Err(ContaminationArtefactError::DuplicateSampleName {
                    name: s.name.clone(),
                });
            }
            if !batch_ids.contains_key(s.batch.as_str()) {
                return Err(ContaminationArtefactError::DanglingSampleBatch {
                    sample: s.name.clone(),
                    batch: s.batch.clone(),
                });
            }
            if !(s.contamination_fraction.is_finite()
                && (0.0..=1.0).contains(&s.contamination_fraction))
            {
                return Err(
                    ContaminationArtefactError::ContaminationFractionOutOfRange {
                        sample: s.name.clone(),
                        got: s.contamination_fraction,
                    },
                );
            }
        }
        Ok(())
    }

    /// Build the engine-side [`ContaminationEstimates`] for a
    /// **specific** cohort sample order. The output's
    /// `c_s_per_sample[i]`, `sample_to_batch[i]`, and
    /// `q_b_per_batch[sample_to_batch[i]]` describe the sample at
    /// position `i` of `sample_names`.
    ///
    /// Sample-name reconciliation:
    ///
    /// - Every entry in `sample_names` must appear in
    ///   [`Self::samples`]; missing samples surface
    ///   [`ContaminationArtefactError::SampleMissingFromArtefact`].
    /// - Extras in [`Self::samples`] that are not in `sample_names`
    ///   are silently dropped — users may run contamination over a
    ///   superset cohort and then call on a subset.
    ///
    /// The returned `ContaminationEstimates::source` is
    /// [`ContaminationEstimateSource::UserSupplied`]
    /// (the artefact was loaded from disk; from the engine's
    /// perspective it is externally supplied).
    pub fn to_estimates_for_samples(
        &self,
        sample_names: &[&str],
    ) -> Result<ContaminationEstimates, ContaminationArtefactError> {
        // Build batch_id -> dense_batch_idx for the subset of batches
        // referenced by the cohort. Batches that don't show up under
        // any cohort sample are also dropped (no engine-side overhead).
        let by_name: BTreeMap<&str, &SampleEntry> =
            self.samples.iter().map(|s| (s.name.as_str(), s)).collect();
        let by_batch: BTreeMap<&str, &BatchEntry> =
            self.batches.iter().map(|b| (b.id.as_str(), b)).collect();

        // First pass: enumerate the batches the cohort actually uses,
        // in stable artefact order. A batch appears in the dense table
        // only if at least one cohort sample references it.
        let mut dense_batches: Vec<&BatchEntry> = Vec::new();
        let mut batch_idx_for_id: BTreeMap<&str, usize> = BTreeMap::new();
        for cohort_name in sample_names {
            let s = by_name.get(cohort_name).ok_or_else(|| {
                ContaminationArtefactError::SampleMissingFromArtefact {
                    sample: (*cohort_name).to_string(),
                }
            })?;
            if let std::collections::btree_map::Entry::Vacant(e) =
                batch_idx_for_id.entry(s.batch.as_str())
            {
                let b = by_batch.get(s.batch.as_str()).ok_or_else(|| {
                    // Caught by `validate()` already, but defensive
                    // for the case where a caller hand-builds the
                    // artefact without calling `validate()`.
                    ContaminationArtefactError::DanglingSampleBatch {
                        sample: s.name.clone(),
                        batch: s.batch.clone(),
                    }
                })?;
                e.insert(dense_batches.len());
                dense_batches.push(b);
            }
        }

        // Second pass: build the three engine-side vectors aligned
        // with `sample_names`.
        let mut c_s_per_sample: Vec<Option<f64>> = Vec::with_capacity(sample_names.len());
        let mut sample_to_batch: Vec<usize> = Vec::with_capacity(sample_names.len());
        for cohort_name in sample_names {
            let s = by_name[cohort_name];
            c_s_per_sample.push(Some(s.contamination_fraction));
            sample_to_batch.push(batch_idx_for_id[s.batch.as_str()]);
        }
        let q_b_per_batch: Vec<[f64; N_ALLELE_CLASSES]> = dense_batches
            .iter()
            .map(|b| {
                [
                    b.contaminant_ref_prob,
                    b.contaminant_snp_alt_prob,
                    b.contaminant_indel_alt_prob,
                ]
            })
            .collect();

        ContaminationEstimates::from_user_supplied(c_s_per_sample, q_b_per_batch, sample_to_batch)
            .map_err(ContaminationArtefactError::EngineConversion)
    }
}

// ---------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------

/// Errors raised when reading, writing, validating, or converting a
/// [`ContaminationArtefact`].
#[non_exhaustive]
#[derive(Error, Debug)]
pub enum ContaminationArtefactError {
    #[error("contamination artefact {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("contamination artefact {path}: TOML parse failed — {source}")]
    Toml {
        path: PathBuf,
        #[source]
        source: toml::de::Error,
    },

    #[error("contamination artefact {path}: TOML serialise failed — {source}")]
    Serialize {
        path: PathBuf,
        #[source]
        source: toml::ser::Error,
    },

    #[error("contamination artefact: batch id `{id}` appears twice")]
    DuplicateBatchId { id: String },

    #[error("contamination artefact: sample `{name}` appears twice")]
    DuplicateSampleName { name: String },

    #[error("contamination artefact: sample `{sample}` references undefined batch `{batch}`")]
    DanglingSampleBatch { sample: String, batch: String },

    #[error(
        "contamination artefact: batch `{batch_id}` field `{field}` = {got} \
         is not a finite probability in [0.0, 1.0]"
    )]
    ProbabilityOutOfRange {
        batch_id: String,
        field: &'static str,
        got: f64,
    },

    #[error(
        "contamination artefact: batch `{batch_id}` contaminant_*_prob fields \
         sum to {sum}; expected 1.0 within {Q_B_SIMPLEX_TOLERANCE:e}"
    )]
    BatchProbabilitiesDoNotSum { batch_id: String, sum: f64 },

    #[error(
        "contamination artefact: sample `{sample}` contamination_fraction = {got} \
         is not a finite probability in [0.0, 1.0]"
    )]
    ContaminationFractionOutOfRange { sample: String, got: f64 },

    /// A cohort `.psp` sample is missing from the artefact's
    /// `[[samples]]` table. Extras are tolerated; absences are not.
    #[error(
        "contamination artefact: cohort sample `{sample}` is missing from the \
         estimates file (extras are tolerated, absences are not)"
    )]
    SampleMissingFromArtefact { sample: String },

    /// Engine-side conversion via
    /// [`ContaminationEstimates::from_user_supplied`] failed despite
    /// artefact-level validation passing. Indicates an invariant the
    /// engine enforces that the artefact validator did not catch —
    /// typically a recently added check; should be promoted to an
    /// artefact-level rule.
    #[error("contamination artefact: engine-side conversion failed — {0}")]
    EngineConversion(ContaminationEstimationError),
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn fixture() -> ContaminationArtefact {
        ContaminationArtefact {
            provenance: Provenance {
                tool: "pop_var_caller".into(),
                version: "0.1.0".into(),
                subcommand: "estimate-contamination".into(),
                created: "2026-05-19T14:32:11Z".parse().unwrap(),
                inputs: ProvenanceInputs {
                    reference: "grch38.fa".into(),
                    input_psps: vec!["NA12878.psp".into(), "NA12891.psp".into()],
                    batch_assignment: Some("lanes.tsv".into()),
                },
            },
            parameters: {
                let mut m = BTreeMap::new();
                m.insert("block_size".into(), toml::Value::Integer(1000));
                m.insert(
                    "stopping_mode".into(),
                    toml::Value::String("convergence".into()),
                );
                m
            },
            batches: vec![BatchEntry {
                id: "lane_3".into(),
                contaminant_ref_prob: 0.9989,
                contaminant_snp_alt_prob: 0.001,
                contaminant_indel_alt_prob: 0.0001,
            }],
            samples: vec![
                SampleEntry {
                    name: "NA12878".into(),
                    batch: "lane_3".into(),
                    contamination_fraction: 0.0123,
                },
                SampleEntry {
                    name: "NA12891".into(),
                    batch: "lane_3".into(),
                    contamination_fraction: 0.0031,
                },
            ],
        }
    }

    // ---- happy path round-trip -------------------------------------

    #[test]
    fn validate_accepts_fixture() {
        fixture().validate().unwrap();
    }

    #[test]
    fn round_trips_through_disk() {
        let dir = tempdir().unwrap();
        let out = dir.path().join("estimates.toml");
        let original = fixture();
        original.write(&out).unwrap();
        let loaded = ContaminationArtefact::read(&out).unwrap();
        assert_eq!(loaded.batches, original.batches);
        assert_eq!(loaded.samples, original.samples);
        assert_eq!(loaded.provenance.tool, original.provenance.tool);
        assert_eq!(loaded.provenance.subcommand, "estimate-contamination");
        assert_eq!(loaded.provenance.inputs.input_psps.len(), 2);
        assert_eq!(
            loaded.provenance.inputs.batch_assignment.as_deref(),
            Some("lanes.tsv")
        );
    }

    #[test]
    fn omits_batch_assignment_when_none() {
        let mut a = fixture();
        a.provenance.inputs.batch_assignment = None;
        let dir = tempdir().unwrap();
        let out = dir.path().join("estimates.toml");
        a.write(&out).unwrap();
        let loaded = ContaminationArtefact::read(&out).unwrap();
        assert!(loaded.provenance.inputs.batch_assignment.is_none());
    }

    // ---- validator rejections --------------------------------------

    #[test]
    fn rejects_duplicate_batch_id() {
        let mut a = fixture();
        let mut dup = a.batches[0].clone();
        dup.id = "lane_3".into();
        a.batches.push(dup);
        assert!(matches!(
            a.validate().unwrap_err(),
            ContaminationArtefactError::DuplicateBatchId { .. }
        ));
    }

    #[test]
    fn rejects_duplicate_sample_name() {
        let mut a = fixture();
        a.samples.push(SampleEntry {
            name: "NA12878".into(),
            batch: "lane_3".into(),
            contamination_fraction: 0.05,
        });
        assert!(matches!(
            a.validate().unwrap_err(),
            ContaminationArtefactError::DuplicateSampleName { .. }
        ));
    }

    #[test]
    fn rejects_dangling_sample_batch() {
        let mut a = fixture();
        a.samples[0].batch = "lane_42".into();
        assert!(matches!(
            a.validate().unwrap_err(),
            ContaminationArtefactError::DanglingSampleBatch { .. }
        ));
    }

    #[test]
    fn rejects_non_finite_prob_in_batch() {
        let mut a = fixture();
        a.batches[0].contaminant_ref_prob = f64::NAN;
        assert!(matches!(
            a.validate().unwrap_err(),
            ContaminationArtefactError::ProbabilityOutOfRange { .. }
        ));
    }

    #[test]
    fn rejects_out_of_range_prob_in_batch() {
        let mut a = fixture();
        a.batches[0].contaminant_ref_prob = 1.5;
        assert!(matches!(
            a.validate().unwrap_err(),
            ContaminationArtefactError::ProbabilityOutOfRange { .. }
        ));
    }

    #[test]
    fn rejects_batch_probs_not_summing_to_one() {
        let mut a = fixture();
        a.batches[0].contaminant_snp_alt_prob = 0.5; // sum becomes ~1.499
        assert!(matches!(
            a.validate().unwrap_err(),
            ContaminationArtefactError::BatchProbabilitiesDoNotSum { .. }
        ));
    }

    #[test]
    fn accepts_all_zero_qb_row_for_floored_batch() {
        // Engine convention: floored batches carry an all-zero q_b
        // vector. The validator must accept this even though it does
        // not lie on the probability simplex.
        let mut a = fixture();
        a.batches[0].contaminant_ref_prob = 0.0;
        a.batches[0].contaminant_snp_alt_prob = 0.0;
        a.batches[0].contaminant_indel_alt_prob = 0.0;
        a.validate().unwrap();
    }

    #[test]
    fn rejects_contamination_fraction_negative() {
        let mut a = fixture();
        a.samples[0].contamination_fraction = -0.001;
        assert!(matches!(
            a.validate().unwrap_err(),
            ContaminationArtefactError::ContaminationFractionOutOfRange { .. }
        ));
    }

    #[test]
    fn rejects_contamination_fraction_nan() {
        let mut a = fixture();
        a.samples[0].contamination_fraction = f64::NAN;
        assert!(matches!(
            a.validate().unwrap_err(),
            ContaminationArtefactError::ContaminationFractionOutOfRange { .. }
        ));
    }

    // ---- sample-name reconciliation in to_estimates_for_samples ----

    #[test]
    fn to_estimates_with_exact_cohort() {
        let a = fixture();
        let cohort = ["NA12878", "NA12891"];
        let est = a.to_estimates_for_samples(&cohort).unwrap();
        // c_s aligned with cohort order.
        assert_eq!(est.effective_c_s(0), 0.0123);
        assert_eq!(est.effective_c_s(1), 0.0031);
        // Both samples in the same batch → q_b is shared.
        let q0 = est.q_b_for_sample(0);
        let q1 = est.q_b_for_sample(1);
        assert_eq!(q0, q1);
    }

    #[test]
    fn to_estimates_rejects_missing_sample() {
        let a = fixture();
        let cohort = ["NA12878", "NA00000"];
        let err = a.to_estimates_for_samples(&cohort).unwrap_err();
        match err {
            ContaminationArtefactError::SampleMissingFromArtefact { sample } => {
                assert_eq!(sample, "NA00000");
            }
            other => panic!("expected SampleMissingFromArtefact, got {other:?}"),
        }
    }

    #[test]
    fn to_estimates_silently_drops_extras() {
        let mut a = fixture();
        // Artefact has an extra sample the cohort doesn't ask for.
        a.samples.push(SampleEntry {
            name: "NA_EXTRA".into(),
            batch: "lane_3".into(),
            contamination_fraction: 0.07,
        });
        // Cohort asks for only one sample.
        let cohort = ["NA12878"];
        let est = a.to_estimates_for_samples(&cohort).unwrap();
        assert_eq!(est.effective_c_s(0), 0.0123);
    }

    #[test]
    fn to_estimates_with_reordered_cohort_aligns_correctly() {
        let a = fixture();
        let cohort = ["NA12891", "NA12878"]; // reversed
        let est = a.to_estimates_for_samples(&cohort).unwrap();
        assert_eq!(est.effective_c_s(0), 0.0031);
        assert_eq!(est.effective_c_s(1), 0.0123);
    }

    // ---- file-level error paths ------------------------------------

    #[test]
    fn read_io_error_on_missing_file() {
        let err = ContaminationArtefact::read(Path::new("/no/such/file.toml")).unwrap_err();
        assert!(matches!(err, ContaminationArtefactError::Io { .. }));
    }

    #[test]
    fn read_parse_error_on_garbage() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("bad.toml");
        fs::write(&path, "not = valid = toml").unwrap();
        let err = ContaminationArtefact::read(&path).unwrap_err();
        assert!(matches!(err, ContaminationArtefactError::Toml { .. }));
    }

    #[test]
    fn write_then_read_validates() {
        // An artefact that fails validation at write time must surface
        // the failure rather than silently producing an unreadable
        // file.
        let mut a = fixture();
        a.batches[0].contaminant_ref_prob = 2.0; // out of range
        let dir = tempdir().unwrap();
        let path = dir.path().join("bad.toml");
        let err = a.write(&path).unwrap_err();
        assert!(matches!(
            err,
            ContaminationArtefactError::ProbabilityOutOfRange { .. }
        ));
        // And the file should not exist.
        assert!(!path.exists());
    }
}
