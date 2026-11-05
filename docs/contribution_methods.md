# Molecular contributions to the treatment contrasts

This workflow reproduces the decomposition used for the updated analyses without
ancestry principal components. It decomposes **chiglitazar minus comparator
differences in mol-HOMA-IR change**, rather than attributing within-arm changes to
the drug. A comparator-arm change can therefore contribute to a treatment
contrast even when the feature changes little in the chiglitazar arm.

## Run the two stages

After training the frozen baseline models:

```bash
Rscript scripts/extract_contributions.R \
  --model-run results/model_run \
  --output results/contribution_inputs

python scripts/treatment_contributions.py \
  --input results/contribution_inputs \
  --output results/contributions \
  --bootstrap 2000 --min-folds 1 --top-n 10
```

Use the actual training-output path in place of `results/model_run`.
`--save-bootstrap` additionally saves feature-level bootstrap arrays.
`--bootstrap 0` performs a point-estimate audit only; it does not estimate
uncertainty. `--top-n` changes the displayed subset, but never changes the
denominator or full results. The full 85/80-feature versions are available in
`ranked_contributors.csv`, `heatmap_adjusted_molecular_changes.csv` and
`all_molecular_contributions.csv` regardless of `--top-n`.

The publication uses `lambda.1se`, 20 repetitions, ten outer folds and 2,000
bootstrap replicates. Smaller synthetic examples test the implementation but do
not reproduce those publication results.

## 1. Frozen per-participant decomposition

For each baseline outer fold, the same frozen imputation, scaling, single-omics
coefficients and integration weights are applied to baseline and week-24
measurements from its held-out participants. The molecular input matrix is
already on the scale used for training; it is not transformed a second time.

Let `z_ijtr` be feature `j` for participant `i` at visit `t`, standardized using
the baseline training fold in repetition `r`. Let `b_jr` be its elastic-net
coefficient on the standardized response scale, `s_yr` the training-response
standard deviation and `w_pr` the OLS integration weight for its omics platform.
Its contribution to that participant's score change is

```text
c_ijr = w_pr × b_jr × s_yr × (z_ij,W24,r − z_ij,W0,r)
c_ij  = mean_r(c_ijr)
```

The multiplication by the response standard deviation restores log HOMA-IR
units. The integration and single-omics intercepts cancel in the within-person
change. Features with zero effective coefficients have zero contributions.
All candidate features remain in the exported matrix, including those not
selected in any model. For the original analysis these comprise 217 metabolites
and 837 proteins, or 1,054 features in total.

The extractor verifies, within each fold and for the averaged ensemble, that

```text
sum_j(c_ij) = predicted log HOMA-IR_i,W24 − predicted log HOMA-IR_i,W0
```

It also checks that held-out participants are absent from each frozen model's
training IDs. Each paired participant contributes one held-out prediction per
repetition. The extraction uses the paired W0/W24 population of the supplied
model run; the manuscript population comprises 593 participants, with 219 in
CMAP and 374 in CMAS.

For the heatmap, the extractor separately averages each feature's frozen
standardized change `z_W24 − z_W0` over repetitions. This is a molecular change
in training-fold SD units, not a score contribution weighted by coefficients.

## 2. Trial-specific treatment contrasts

Within each trial, all molecular contributions use the same OLS design:

```text
c_j ~ randomized arm + baseline mol-HOMA-IR + sex + age
```

Baseline and follow-up mol-HOMA-IR scores are on the log scale. Sex is coded
1=male and 2=female; CG32 is the reference arm. Age is standardized within the
trial using the sample standard deviation. The two contrasts are CG32 minus
the trial-specific comparator and CG48 minus the comparator: placebo in CMAP
and sitagliptin in CMAS. The analysis includes neither ancestry PCs nor BMI or
BMI change. Trial is not a covariate because the two trials are fitted
separately.

The contribution analysis uses paired molecular scores and complete baseline
sex and age. It does not require ancestry-PC availability. Exclusions for
missing required covariates are recorded in `analysis_settings.json`.

If `theta_jd` is the adjusted molecular contribution for dose contrast `d`,
linearity ensures that its sum across features is the corresponding total
mol-HOMA-IR treatment effect. The code checks this against the ANCOVA

```text
week-24 mol-HOMA-IR ~ randomized arm + baseline mol-HOMA-IR + sex + age
```

Using the follow-up score or its change as the outcome gives the same treatment
contrasts because the baseline score is included in the design. Ordinary OLS
standard errors and t-based 95% confidence intervals are exported for each
molecular contrast.

## 3. Contribution percentages and uncertainty

The percentage contribution is

```text
share_jd = 100 × max(−theta_jd, 0) / sum_k max(−theta_kd, 0)
```

The denominator includes **all features with negative adjusted contributions**
for that trial and dose, before any stability, direction or display filtering.
It is the gross lowering component, not the net treatment difference. Positive
contributions offset part of the lowering component. Consequently, percentages
sum to 100 across the lowering features, while the subset displayed in a figure
usually sums to less than 100. A contrast with no lowering component has an
undefined denominator and is reported as an error.

Participants are resampled with replacement within each trial and randomized
arm, retaining the original arm sample sizes. Every resampled participant
contributes their complete vector of molecular contributions, preserving
cross-feature dependence. Each replicate refits the contribution ANCOVAs and
recalculates the full lowering denominator. The original RNG seeds are
202610031 for CMAP and 202610032 for CMAS; draw batches contain 50 replicates.
The standard deviation of the bootstrap percentages, using `ddof=1`, gives
the bootstrap standard error. Percentile 95% intervals are also exported.

These uncertainties are conditional on the frozen molecular models. Model
training, candidate selection and ranking are not repeated within this
bootstrap; their uncertainty is not propagated. The resampling is not a
mediation analysis.

## 4. Candidate eligibility, direction filtering and ranking

A feature is eligible when its effective integrated coefficient is nonzero in
at least one of ten outer-fold models **in every one of the 20 repetitions**.
This is not an overall 10% frequency threshold pooled across 200 models.
The manuscript has 148 such features: 83 metabolites and 65 proteins.

Within each trial, an eligible feature is retained when both adjusted
dose-control contributions are negative, using a numerical tolerance of
`theta < −1e−12`. There is no P-value or confidence-interval filter. The
retained features are ranked separately in CMAP and CMAS by the arithmetic
mean of their two dose-specific lowering percentages. This gives:

| Set | Metabolites | Proteins | Total |
| --- | ---: | ---: | ---: |
| Eligible in every repetition | 83 | 65 | 148 |
| Both CMAP contrasts lower mol-HOMA-IR | 51 | 34 | 85 |
| Both CMAS contrasts lower mol-HOMA-IR | 42 | 38 | 80 |
| Shared by both trials | 31 | 21 | 52 |

These counts are reference results from the original data, not hard-coded
requirements for other datasets. The direction filter concerns **weighted
contributions to the score contrast**. It does not require the unweighted
molecular change to have the same sign in every comparison or demonstrate an
independent causal effect of the molecule.

Heatmap adjusted means use the same trial-specific design, with the averaged
standardized molecular change as the outcome. For each arm, covariates are
retained at every eligible participant's observed values while the arm is set
to that arm. Predictions are then averaged across the trial's eligible
population; ordinary OLS standard errors are exported.

## Files and reproducibility checks

The extraction stage writes the five inputs to the Python stage:

| File | Contents |
| --- | --- |
| `participant_contributions.csv` | Participant ID and all molecular contributions in log HOMA-IR units |
| `participant_molecular_changes.csv` | Participant ID and all averaged frozen-SD molecular changes |
| `participant_scores_covariates.csv` | Participant ID, trial, arm, baseline sex/age and W0/W24 mol scores |
| `feature_annotations.csv` | Feature ID, label and omics platform |
| `selection_by_repeat.csv` | Selected outer-fold count for every feature and repetition, including zeros |

The Python stage exports full contrasts and shares, eligibility and direction
audits, trial-specific rankings, shared features, heatmap means, display
subsets, aggregate decomposition checks, input-file hashes and the actual
analysis settings. Output tables do not contain participant IDs. Intermediate
extraction files and model bundles do contain participant-level data and must
remain private. Only synthetic example data belong in this public repository.

The original aggregate reconstruction errors are below `1.1e−14` at the
participant level and below `1.6e−15` for the sum of treatment contrasts.
The four total treatment contrasts, in the order CMAP CG32, CMAP CG48, CMAS
CG32 and CMAS CG48, are approximately −0.263428, −0.314353, −0.164639 and
−0.245988 log HOMA-IR units. These are implementation checks rather than
additional statistical tests.
