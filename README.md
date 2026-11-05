# Molecular HOMA-IR modeling and treatment-effect decomposition

Code for two workflows in the ChiHOPE molecular HOMA-IR study:

1. Train metabolomics and proteomics elastic-net models for baseline log HOMA-IR, integrate their cross-fitted predictions, and apply frozen baseline models at week 24.
2. Decompose adjusted chiglitazar–control differences in the change in mol-HOMA-IR into analyte contributions, calculate lowering-direction contribution percentages, and filter/rank contributors separately in CMAP and CMAS.

This release corresponds to the **705-person baseline / 593-person paired week-24 analysis**, with **no ancestry-PC adjustment** in the treatment-contribution models. Original individual data and fitted models are not distributed. An executable synthetic demonstration is included.

[中文说明](README_zh.md) · [Input requirements](docs/data_requirements.md) · [Model methods](docs/model_methods.md) · [Contribution methods](docs/contribution_methods.md) · [Provenance](docs/provenance.md)

Version **1.0.0** · [MIT License](LICENSE) · [Software citation](CITATION.cff) · [Code Ocean setup](docs/code_ocean.md)

## Install

Use R 4.5 with `glmnet`, `matrixStats`, `dplyr` and `jsonlite`, and Python 3.11 with NumPy, pandas and SciPy. A Conda environment specification is provided:

```bash
conda env create -f environment.yml
conda activate molecular-homa-ir
```

If using an existing R/Python environment, install the packages listed in `environment.yml` and `requirements.txt`. No private dependencies, external API credentials or network calls are needed during analysis. Tested versions are recorded in `VALIDATION.md`; the environment file specifies the intended compatible environment rather than an immutable platform lockfile.

## Quick start with synthetic data

Run from the repository root:

```bash
bash tests/run_smoke.sh
```

The smoke test simulates paired data, runs a deliberately reduced nested cross-validation configuration, extracts contributions, performs trial-specific analysis and runs numerical tests. It does not reproduce clinical results. To select executables explicitly:

```bash
PYTHON=/path/to/python RSCRIPT=/path/to/Rscript bash tests/run_smoke.sh
```

The same steps can be run individually:

```bash
python examples/generate_synthetic_data.py --output-dir data/synthetic
Rscript scripts/train_models.R --config config/smoke.json
Rscript scripts/extract_contributions.R --model-run results/smoke --output results/smoke_decomposition
python scripts/treatment_contributions.py --input results/smoke_decomposition --output results/smoke_contributions --bootstrap 100 --min-folds 1 --top-n 10
Rscript tests/test_model_invariants.R results/smoke
python -m unittest discover -s tests -p 'test_*.py'
```

## Run on authorized study data

1. Prepare the eligible cohort and processed assay tables described in `docs/data_requirements.md`. Place private inputs under `data/private/`, or use paths outside the repository. Do not commit these inputs or generated individual-level results.
2. Copy/edit `config/default.json`, supplying the input paths, output directory and original outer-fold table when exact manuscript reproduction is required. The full configuration uses 20 × 10 outer cross-validation, five-fold tuning and five-fold integration cross-fitting. It is substantially more expensive than the smoke test.
3. Run the following commands, replacing the paths if your configuration uses another output folder:

```bash
Rscript scripts/train_models.R --config config/default.json
Rscript scripts/extract_contributions.R --model-run results/study --output results/study_decomposition
python scripts/treatment_contributions.py --input results/study_decomposition --output results/study_contributions --bootstrap 2000 --min-folds 1 --top-n 10
```

`--top-n` changes only the displayed subset. All analyte contrasts, the complete trial-specific ranked sets and their intersection are exported regardless. `--bootstrap 0` is available for a point-estimate audit; it does not provide the publication uncertainty estimates.

The study's original outer-fold assignments are private. A newly generated plan uses the same balancing approach but will not recreate the exact original 705-person assignments, which were inherited from a predefined 706-person plan. The reference counts below are checked outcomes, not hard-coded selection targets.

## Interpretation of the contribution percentages

For each participant, a molecular contribution is its frozen integrated coefficient multiplied by its standardized week-24-minus-baseline change, averaged over held-out repetitions. Contributions summed over **all** analytes equal the participant's mol-HOMA-IR change.

Within each trial, the model is:

```text
analyte contribution ~ treatment arm + baseline mol-HOMA-IR + sex + age
```

The contrasts are CG32 minus control and CG48 minus control. There is no BMI, BMI-change or ancestry-PC adjustment. For a contrast with analyte effects theta_j:

```text
lowering share_j (%) = 100 × max(-theta_j, 0) / sum_all_analytes max(-theta, 0)
```

Thus, the denominator is the **gross lowering contribution across all modeled analytes**, before frequency filtering or selection of display rows. It is not the net score difference or the sum restricted to displayed contributors. An analyte can contribute to a treatment–control difference because of a change in the control arm; this does not establish a within-treatment effect or causal mediation.

Eligibility requires a nonzero effective integrated coefficient in at least one outer fold in **every** repetition. Within each trial, eligible features must have a negative contribution in both dose–control contrasts. No significance threshold is used for this filter. Ranking uses the mean of the two dose-specific lowering percentages. The 2,000 participant bootstraps resample within trial and arm, refit contribution regressions and recalculate the denominator; molecular models remain frozen.

## Main outputs

| Stage | Output | Purpose |
| --- | --- | --- |
| Modeling | Fitted outer-fold models, outer assignments, per-repetition and averaged predictions, tuning/performance and audit files | Repeated held-out evaluation and frozen follow-up prediction. These outputs can contain participant keys and are private. |
| Extraction | `participant_contributions.csv`, `participant_molecular_changes.csv`, `participant_scores_covariates.csv` | Individual decomposition and paired covariates; do not publish. |
| Extraction | `feature_annotations.csv`, `selection_by_repeat.csv`, reconstruction audits | Feature identity, nonzero selection frequencies and decomposition checks. |
| Treatment analysis | `all_molecular_contributions.csv`, `total_effect_decomposition.csv` | All feature-level contrasts, lowering percentages and sum-of-contributions validation. |
| Treatment analysis | `eligible_features.csv`, `ranked_contributors.csv`, `shared_contributors.csv` | Frequency-eligible, trial-specific both-dose-lowering, and shared feature sets. |
| Treatment analysis | `displayed_features.csv`, `displayed_contributions.csv`, `heatmap_adjusted_molecular_changes.csv` | Display subsets and adjusted molecular changes; full ranked results remain available. |
| Treatment analysis | `analysis_settings.json`, `decomposition_audit.csv`, `source_sha256.json` | Model definition, sample sizes, numerical checks and input fingerprints. |

With authorized original inputs, the publication version yields 148 frequency-eligible analytes (83 metabolites, 65 proteins), 85 both-dose-lowering contributors in CMAP (51/34), 80 in CMAS (42/38), and 52 shared contributors (31/21). Aggregate reference values are supplied in `reference/manuscript_aggregate_expectations.json`.

## Scope, release and license

The package starts from processed and eligible input tables. It does not perform raw assay processing, trial recruitment/QC, all clinical association analyses, or manuscript figure assembly. It contains neither a diagnostic deployment model nor a license for clinical data.

Source code is distributed under the [MIT License](LICENSE), copyright 2026 Zehui Cao and contributors. Software citation metadata are provided in `CITATION.cff` and `.zenodo.json`; both describe the same software release. A Zenodo DOI and a Code Ocean capsule link will be added after those records are created and verified. No unissued DOI or unverified manuscript citation is supplied. The `.gitignore` excludes input data, generated results, fitted objects and local caches; check the files being uploaded because ignore rules do not remove previously tracked files.
