# Molecular HOMA-IR model

The target is the natural logarithm of measured HOMA-IR. Metabolomics and proteomics are fitted separately using elastic-net regression, and the two predictions are integrated by an outer-training-only OLS stack. This is the **direct, nonresidualized model**. BMI, age, sex and genetic ancestry PCs are not molecular model predictors. The legacy preprocessing object retains an auxiliary clinical regression, but its fitted values are not subtracted from the direct target, subtracted from the analytes, or added to the direct prediction.

## Preprocessing and tuning

Every training split estimates its own preprocessing parameters. Features need at least two finite training measurements. Missing measurements are replaced by training medians; features are centered and divided by the training sample standard deviation. Features with a standard deviation at most `1e-10` are discarded in that split. The log HOMA-IR target is also standardized using the training mean and sample standard deviation. Prediction and validation loss are returned to the original log HOMA-IR scale.

The input analytes are **already on the study's processed assay scale**. The code does not apply a second log transformation. The study candidate panel contains 217 metabolites and 837 proteins; D-glucose was excluded a priori. No univariate association screen is used to select features before elastic-net training.

Publication defaults are 20 repetitions of 10-fold outer cross-validation, five inner tuning folds, alpha values `0.1, 0.3, 0.5, 0.7, 0.9, 1`, and 60 descending lambda values from 10 to 0.001. Validation loss is the equally weighted mean of the five fold-specific MSEs. Alpha is selected at the minimum MSE across the complete grid. At that alpha, the primary model takes the largest lambda within one standard error of the minimum; the standard error is the SD of the five fold-specific MSEs divided by the square root of five. The minimum-MSE rule is retained as a secondary output. The contribution analysis uses only `lambda.1se`.

## Score integration

Within each outer training set, an additional five-fold cross-fitting loop produces held-out predictions from both single-omics models. In each cross-fitting training subset, each single-omics model is retuned using five-fold inner cross-validation. The outer-training cross-fitted predictions enter:

```text
log(HOMA-IR) = b0 + bM × met-HOMA-IR + bP × prot-HOMA-IR + error.
```

The fitted intercept and weights are applied to the two outer-test predictions. Weights are unconstrained: they need not be positive or sum to one. A nearly constant base score is treated as constant; unidentifiable coefficients are zero. This guard avoids extremely large weights on floating-point noise. Outer-test participants do not enter preprocessing, tuning, base-model fitting, or estimation of integration weights.

## Fold plan and follow-up

Split generation balances trial, measured IR category, sex and BMI category. Optional age balancing exchanges participants within the same joint stratum to improve outer-fold age means and age-quartile balance; it does not optimize molecular prediction performance. Inner splits use the original seeded allocator without the outer age-balancing swaps.

The manuscript's 705 participants retained a predefined age-balanced fold plan originally constructed in 706 participants; all molecular models were then refitted on the eligible 705. Exact reproduction therefore requires that original private plan. Generating a fresh plan directly in 705 participants is a new analysis and will not reproduce the precise published coefficients.

For week 24, each participant is evaluated with the frozen baseline model for their baseline outer-test fold. Training medians, feature and target scaling, elastic-net coefficients and integration weights are all frozen. Week-24 values are not used to train or recalibrate the model. The 20 held-out predictions are averaged on the log scale separately at each visit. Exponentiating that average gives a geometric-scale estimate, not the arithmetic average of exponentiated predictions.

## Numerical compatibility

The core fitting code is extracted from the study implementation. In particular, its historical `glmnet` invocation is retained for numerical compatibility. The tested runtime is recorded in the run output and release validation. Changes in R, `glmnet`, floating-point libraries, input row ordering or the fold plan can change numerical results; the synthetic smoke configuration is a test configuration, not the publication analysis.

The package implements repeated held-out model development and frozen follow-up evaluation. It does not provide a clinical diagnostic tool, an independently fitted deployment model, or the paper's full set of association and figure-generation analyses.
