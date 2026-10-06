# HCT-BE trajectory study: covariate imputation

Standalone R scripts for multiple imputation of missing covariates in the MIMIC-IV derivation cohort and the eICU-CRD validation cohort.

| Script | Cohort |
| --- | --- |
| `MIMIC_IV_CART_imputation.R` | MIMIC-IV |
| `eICU_CRD_CART_imputation.R` | eICU-CRD |

## Requirements

The analysis used R 4.5.2 and `mice` 3.19.0. Install `mice` and its dependencies before running either script. A different `mice` version may produce different results.

## Input

Each script takes a cleaned, patient-level analytic CSV for its respective cohort. Required column names and coding checks are defined in the script. Predictor-only columns, trajectory assignments, event indicators, and follow-up times must be complete. There is no fixed sample-size or missingness-percentage requirement.

The scripts do not read HCT or BE time-window values or trajectory posterior probabilities into the imputation model. Input data are not included in this repository.

## Run

Run from a terminal with `Rscript`, supplying an input CSV and a new or empty output directory:

```bash
Rscript MIMIC_IV_CART_imputation.R path/to/mimic_cleaned.csv path/to/mimic_output
Rscript eICU_CRD_CART_imputation.R path/to/eicu_cleaned.csv path/to/eicu_output
```

The two cohorts are imputed separately. Each script creates 20 completed datasets after 10 chained-equation iterations. CART uses `minbucket = 5`, `cp = 0.0001`, and `maxdepth = 30`, with seed `20260906`. No separate burn-in period is discarded.

The imputation model includes clinical covariates, trajectory group, 30-day and 60-day event indicators, and the corresponding Nelson-Aalen cumulative hazard estimates. Follow-up times are used to calculate the hazard estimates but are not direct imputation predictors. BMI and 0-48 h fluid intake per kg are recalculated in each completed dataset. In eICU-CRD, missing day-1 and day-2 fluid intake values are imputed without a separate rule based on ICU length of stay.

Missing at random is assumed conditional on the variables included in the imputation model; the scripts do not test this assumption.

## Output

Each output directory contains:

- `completed/imputed_01.csv` through `completed/imputed_20.csv`: completed analytic datasets.
- `imputation.rds`: the `mice` model object.
- `chain_traces.pdf`: chain mean and standard deviation traces.
- `convergence.csv`: numerical diagnostics for chain means.
- `logged_events.txt` or `logged_events.csv`: algorithm events, if any.
- `session_info.txt`: R session and package versions.

Inspect the diagnostic files before interpreting the imputed datasets. Creating diagnostic files does not establish convergence.

These scripts cover covariate imputation. Database extraction, cohort construction, outcome models, and pooling of model estimates are outside their scope. No patient-level input or output files are distributed here.
