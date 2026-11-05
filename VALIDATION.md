# Validation of this code package

Validated locally on 2026-11-05. No original participant-level data are distributed with the package.

## Executable synthetic workflow

The complete smoke command passed using 120 simulated participants, paired weeks 0 and 24, eight simulated metabolites and six simulated proteins. The reduced configuration used two repetitions, three outer folds, two inner folds and two integration cross-fitting folds, followed by 100 arm-stratified bootstrap replicates per trial.

- All six outer jobs completed; every participant had one held-out prediction per repetition at both visits.
- Extracted molecular contributions reconstructed the averaged score change with maximum absolute error below `1.0e-14`.
- Summed adjusted molecular contrasts reconstructed the total adjusted score contrasts with maximum absolute error below `2.8e-15`.
- The R invariant test perturbed all outer-test outcomes/assays and all week-24 outcomes/assays, refitted the entire first outer job including tuning and integration cross-fitting, and verified unchanged fitted preprocessing, selected hyperparameters, coefficients and stack weights.
- Changing week-24 outcome labels alone left frozen predictions unchanged.
- Nine Python unit tests passed. They check contrast signs, linear decomposition, the gross-lowering denominator, absence of display-subset renormalization, every-repetition eligibility, invalid/missing selection records, and a deterministic joint participant bootstrap.

The test configuration and its resulting feature counts are synthetic; they are not clinical results.

## Comparison with the archived study implementation

Using authorized original inputs outside the public package:

1. **Representative model job:** repetition 1, outer fold 1. Frozen metabolomics and proteomics predictions matched the archived outputs exactly. The metabolomics model was retuned across the full five-fold, six-alpha, 60-lambda grid; tuning MSE, selected settings and refitted coefficients matched the original implementation and archived fit exactly. OLS stack weights also matched the original implementation exactly. This check is not a complete refit of all 200 study jobs.
2. **All contribution contrasts:** all 4,216 feature–contrast rows were compared with the final no-PC analysis. Maximum absolute difference was `1.01e-16` for log-scale effects and `5.33e-15` for lowering percentages.
3. **Full bootstrap:** both trials were rerun with the original 2,000 replicates and seeds. Maximum absolute difference was `3.45e-15` for bootstrap SE and `5.33e-15` for percentile interval endpoints.
4. **Feature selection and ranking:** the 148 eligible analytes, ordered 85/80 trial-specific contributor lists and 52-analyte intersection matched exactly.
5. **Adjusted molecular changes:** all 495 trial/arm/qualifying-feature heatmap rows were checked. Maximum absolute differences were `2.20e-15` for adjusted means and `1.12e-16` for SE.
6. **Complete frozen-model extraction:** the new R extractor was checked against all 200 archived outer models for 593 paired participants and 1,054 analytes. The contribution matrix differed by at most `8.28e-25`; standardized-change matrices and every feature-by-repetition selection count matched exactly. These archived models were read and applied, not all retrained.
7. **Input adapter:** a private CSV round trip of the 705/593 input population preserved participant/feature ordering, missingness, clinical values, IR/BMI strata and inner-fold assignments. Molecular numeric round-trip error was at most `5.33e-14`; the temporary private CSVs were removed after the check.

These are implementation-equivalence and numerical conservation checks. They do not constitute an independent clinical validation or propagate molecular model-training uncertainty into the bootstrap.

## Tested runtime

| Component | Version |
| --- | --- |
| R | 4.5.3 |
| glmnet | 5.0 |
| matrixStats | 1.5.0 |
| dplyr | 1.2.1 |
| jsonlite | 2.0.0 |
| Python | 3.11.16 |
| NumPy | 2.4.6 |
| pandas | 3.0.5 |
| SciPy | 1.17.1 |

The synthetic end-to-end workflow and study-output comparisons were run in this existing environment. The Conda specification was not installed into a second clean environment during this release check. Cross-platform installation and bit-for-bit equivalence on other BLAS/R builds are not claimed.

## Publication adapters

The Code Ocean master-entry adapter passed the full synthetic workflow and scientific tests using the existing local environment. It also passed after relocation to a path containing spaces, while keeping the code tree unchanged. Synthetic input creation occurs in scratch storage; results and logs are persisted outside the code directory. The Linux post-install recipe and a remote Code Ocean reproducible run have not yet been executed.

`CITATION.cff` passed the official Citation File Format 1.2.0 JSON schema. Author, title, version and MIT license metadata agree across `CITATION.cff`, `.zenodo.json`, `VERSION` and `LICENSE`.
