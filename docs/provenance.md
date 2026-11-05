# Implementation provenance

This package isolates two workflows from the current 705-participant molecular HOMA-IR analysis: repeated nested model development with frozen week-24 prediction, and decomposition of adjusted chiglitazar–control score contrasts. The publication version uses 593 paired participants and does not adjust treatment contribution models for ancestry PCs.

The scientific implementation was traced to these original project scripts:

| Source | Retained role |
| --- | --- |
| `78_routine_benchmarks_705.R` | Current eligible population, baseline outer-fold refits and frozen follow-up evaluation. Clinical benchmark, ablation and augmentation branches are excluded. |
| `irscore_helpers.R` | Seeded fold allocation, training-only preprocessing, elastic-net tuning, fitting and prediction. |
| `stacking_pilot_helpers.R`, `pheno_stacking_helpers.R` | Outer-training cross-fitted single-omics scores and OLS integration. |
| `shared_panel_helpers.R` | Near-constant base-score guard and restoration of the original seeded fold allocator. |
| `08_age_balanced_sensitivity.R` | Outer-fold age balancing through within-stratum swaps. |
| `137_extract_molecular_change_705.R` | Effective integrated coefficients, paired frozen standardized molecular changes and exact contribution reconstruction. |
| `136_treatment_705_adjusted.R` and the no-PC rerun | Current 593-person trial-specific covariate design. |
| `47_treatment_molecular_contributions.py` | OLS treatment contrasts, all-feature lowering shares and arm-stratified participant bootstrap. |
| `141_molecular_contributions_80pct_trial_specific.py --min-folds 1` | At-least-one-fold-in-every-repetition eligibility, both-dose direction filtering and trial-specific ranking. |

Portable entry points replace hard-coded workstation paths and historical RDS loaders with explicit CSV inputs. Unrelated residual-score, clinical association, epigenetic aging, supplementary figure and manuscript-editing code is omitted. No original participant-level inputs or outputs are distributed.

## Version boundaries

The older broad supplementary-code archive and top-level analysis README include historical 706/594-person results. They are not the authority for this package. The older script `47_extract_frozen_molecular_changes.R` also uses that previous population and is not the current extraction workflow.

The current reference results are:

| Quantity | Current value |
| --- | --- |
| Baseline model population | 705 |
| Paired week-24 population | 593: CMAP 219; CMAS 374 |
| Candidate analytes assayed in the model | 1,054: 217 metabolites; 837 proteins |
| Eligible in at least one outer fold in every repetition | 148: 83 metabolites; 65 proteins |
| Lowering contribution for both doses in CMAP | 85: 51 metabolites; 34 proteins |
| Lowering contribution for both doses in CMAS | 80: 42 metabolites; 38 proteins |
| Common to both trial-specific sets | 52: 31 metabolites; 21 proteins |

These reference counts are outcomes of fitting and direction filtering, not numbers imposed by the code. See `reference/manuscript_aggregate_expectations.json` for aggregate effect and performance checks. No per-participant values are included in that file.

## License and citation

The author selected the MIT License and confirmed the copyright notice **2026 Zehui Cao and contributors**. The confirmed software creator is **Zehui Cao, Fudan University**. The first software version is `1.0.0`; `CITATION.cff` and `.zenodo.json` carry matching metadata. The requested GitHub destination is `zehhhui/molecular-homa-ir`. Zenodo and Code Ocean identifiers will be recorded only after they are issued and verified. This is a new software release and must not be represented as part of the older ChiHOPE-ML archive.
