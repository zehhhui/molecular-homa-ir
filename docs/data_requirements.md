# Data and cohort requirements

Only synthetic data can be generated from this repository. Original participant data, fitted models, individual predictions and the private fold plan are not included. Real input/output folders must remain outside files committed to a public repository.

## Required input files

All CSV files use UTF-8 and a header. `sample_id` is an opaque participant key treated as text, not a sample accession or a name. It must identify the same person at weeks 0 and 24. Each `(sample_id, week)` pair is unique. The code uses this key for pairing and leakage checks; it must not be published for real participants.

| File | Required columns | Meaning |
| --- | --- | --- |
| `phenotype.csv` | `sample_id, week, trial, arm, sex, age, homa_ir, bmi` | Long-format visit metadata. `homa_ir` is on its original positive scale; the code takes its natural logarithm. `bmi` is in kg/m². |
| `metabolites.csv` | `sample_id, week`, then metabolite feature IDs | Processed numeric measurements, on the same assay scale at both visits. |
| `proteins.csv` | `sample_id, week`, then protein feature IDs | Processed numeric measurements, on the same assay scale at both visits. |
| `features.csv` | `feature_id, omics, label` | Unique feature IDs across both platforms; `omics` is `Metabolomics` or `Proteomics`; `label` is the analyte name. |
| Optional `folds.csv` | `sample_id, repeat, fold` | One assignment for every baseline participant in each repetition. Labels start at 1 and must match the configuration. |

`week` is 0 or 24. For the trial-specific contribution workflow, `trial` is `CMAP` or `CMAS`, and `arm` is `CG32`, `CG48`, or the trial-specific control (`PB` in CMAP; `SP` in CMAS). `sex` is 1 for male and 2 for female, following the original data encoding. Age is in years. Treatment, trial and sex must agree across visits; longitudinal treatment adjustment uses baseline age and sex.

The optional phenotype columns `ir_status` and `bmi_category` supply the original fold-balancing labels. Otherwise, measured IR is derived using HOMA-IR thresholds of 4.31 in males and 4.51 in females (`IR` versus `IS` internally), and BMI categories use cut points 18.5, 23 and 25 kg/m² (`Underweight`, `Normal`, `Overweight`, `Obese`). These categories are split-balancing variables, not additional predictors in the molecular model. Supplied fold assignments must use the configured repetition/fold counts. Relative paths in JSON configurations are resolved from the repository root, not the configuration file's directory.

Numeric molecular missing values can be left empty or represented as `NA`; they are imputed separately within each training split. Clinical cohort eligibility must be determined before analysis; clinical missingness must not be silently replaced with invented values. In particular, a completely unmeasured platform is not equivalent to a few missing analytes.

## Reproducing the manuscript population

The manuscript baseline population comprises 705 participants with eligible baseline molecular and clinical data. The actual training entry point additionally required complete positive BMI, waist circumference, triglycerides, HDL cholesterol, ALT and fasting plasma glucose, in addition to age, sex and positive measured HOMA-IR. All required metabolomics blocks and proteomics had to be available. Ancestry PC availability was not required.

The common week-24 analysis comprises 593 members of this baseline population with the required paired molecular and clinical observations: 219 in CMAP and 374 in CMAS. Cohort construction and assay-level QC are upstream of this portable package. Providing the four input CSVs does not by itself prove that these original eligibility rules have been reproduced. Supply the already eligible visits and the curated 217/837 feature panel for exact study reuse.

The same transformed assay definition must be used at weeks 0 and 24; no week-specific scaling or normalization should be introduced by the user. Do not supply unprocessed mass-spectrometry intensities or log-transform already transformed inputs again. Analytical platform preprocessing, peak annotation and assay QC are outside the two workflows released here.

## Exact versus illustrative reproduction

The synthetic example checks executable behavior, model bookkeeping and mathematical identities. Its counts, coefficients and qualifying feature numbers have no biological interpretation and will not equal the manuscript results.

Exact numerical reproduction requires authorized original processed measurements, the original eligible visit list, feature order and baseline row ordering, the private outer-fold assignments, and a compatible software environment. The repository's aggregate reference JSON is a validation target for those authorized inputs, not training data and not an expected result for the synthetic demonstration.
