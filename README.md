# MDD-only versus MDD with persistent depressive disorder

This repository contains the R analysis pipeline for a participant-level classification study comparing major depressive disorder (MDD) without persistent depressive disorder (PDD) with MDD accompanied by PDD. The models combine demographic, baseline questionnaire, passive-sensing, and ecological momentary assessment (EMA) features summarized across several observation windows.

The repository provides analysis code only. It does not contain participant data, derived participant-level predictions, or other potentially identifiable records.

## Preregistration

The study was preregistered on the Open Science Framework on November 12, 2025, after data collection but before the analyses reported in the manuscript:

> Seo H. *Optimizing temporal aggregation windows of passive sensing and baseline survey data for differentiating MDD with and without Persistent Depressive Disorder (PDD).* OSF; 2025. https://doi.org/10.17605/OSF.IO/3TJ8C

The following core elements were retained from the preregistration:

- Participant-level outcome: MDD-only versus MDD with PDD
- Observation windows: 7, 15, 30, 60, and 90 days
- Classifiers: ridge logistic regression, radial-basis-function support vector machine, and XGBoost
- Participant-level nested stratified cross-validation
- Primary performance measure: area under the receiver operating characteristic curve (AUROC)
- Secondary measures including F1 score, accuracy, and Brier score

The implemented analysis differs from the preregistered plan in several respects:

| Component | Preregistered plan | Implemented analysis |
|---|---|---|
| Predictors | Demographic, baseline survey, and specified passive-sensing features | Additional demographic, baseline, passive-sensing, and EMA PHQ-9 features |
| Window features | Mean and SD | Mean, SD, and observed-data coverage |
| Missing data | MICE with 20 imputations fitted within training folds | Training-fold median imputation after window summarization |
| Eligibility and coverage | Daily exclusions based principally on 12-hour missingness rules | At least 45 retained passive days, no more than 20% passive-cell missingness, and at least 70% window coverage for mean and SD features |
| Cross-validation | Nested stratified 5-by-5 cross-validation | Five repeats of 5-fold outer cross-validation with 5-fold inner tuning |
| Hyperparameters | Grids specified in the registration | Revised grids documented in the script |
| XGBoost preprocessing | Raw scales and native missing-value handling | The same fold-specific median imputation and standardization framework used for the other models |
| Statistical uncertainty | Paired DeLong model comparisons with Benjamini-Hochberg correction | Stratified participant-level bootstrap confidence intervals; the planned DeLong comparisons were not performed |
| Secondary analyses | Not specified | Predictor-block, calibration, probability-overlap, and ridge coefficient-stability analyses |

These deviations should be reported in the manuscript and supplement. Analyses involving added predictors, predictor blocks, calibration, and coefficient stability are nonpreregistered secondary analyses.

## Analysis file

- `analysis_peer_review.R`: complete analysis, validation, table-generation, plotting, audit, and reproducibility pipeline

The script starts with prepared daily and participant-level files. Raw smartphone and wearable processing, diagnostic ascertainment, and construction of the source CSV files are outside its scope.

## Software requirements

- R 4.1.0 or later
- R packages:
  - `dplyr`
  - `tidyr`
  - `tibble`
  - `purrr`
  - `readr`
  - `rsample`
  - `glmnet`
  - `e1071`
  - `xgboost`
  - `pROC`
  - `ggplot2`
  - `flextable`
  - `officer`

Install the required packages in R:

```r
install.packages(c(
  "dplyr", "tidyr", "tibble", "purrr", "readr", "rsample",
  "glmnet", "e1071", "xgboost", "pROC", "ggplot2",
  "flextable", "officer"
))
```

Package versions and complete R session information are recorded automatically for each run.

## Required input files

Arrange the prepared data as follows:

```text
project_directory/
├── analysis_peer_review.R
├── final_withNA.csv
├── EMA_final.csv
└── processed_csv/
    ├── demographics.csv
    └── baseline_survey.csv
```

### `final_withNA.csv`

Daily passive-sensing data containing:

- `uid`: participant identifier
- `date`: calendar date in `YYYY-MM-DD` format
- Numeric passive-sensing variables

The script excludes `tp` and `pdd` if present. All other numeric columns are treated as passive-sensing predictors, so administrative variables and numeric identifiers must be removed beforehand.

### `EMA_final.csv`

Daily EMA data containing:

- `uid`
- `date`
- Exactly nine numeric PHQ items, named `PHQ1` through `PHQ9`

An optional period, underscore, or hyphen may occur before the item number, such as `PHQ.1`, `PHQ_1`, or `PHQ-1`. Other EMA variables are excluded from the analysis.

### Participant metadata

`demographics.csv` and `baseline_survey.csv` must each contain one nonmissing, unique row per participant. Together they must supply:

- `uid`
- `pdd`: `0` for MDD-only and `1` for MDD with PDD
- Demographic and baseline variables used by the analysis

Table 2 assumes the following demographic coding:

| Variable | Coding used for display |
|---|---|
| `sex` | 0 = male; 1 = female |
| `job` | 0 = other; 1 = full-time employment |
| `living_status` | 0 = living alone; 1 = living with others |
| `race` | 0 = White; 1 = non-White |
| `education` | 1-3 = high school graduate or less; 4-8 = more than high school |
| `income` | 1-6 = less than $60,000; 7-12 = $60,000 or more |

The models currently treat numeric demographic codes as numeric predictors. Confirm that this matches the intended analysis and the source data dictionary before running the final publication analysis.

## Running the analysis

From a terminal, provide the project directory and a new output directory:

```bash
Rscript analysis_peer_review.R /path/to/project_directory /path/to/results_directory
```

The same paths can be supplied through environment variables:

```bash
export MDD_PDD_PROJECT_DIR=/path/to/project_directory
export MDD_PDD_OUTPUT_DIR=/path/to/results_directory
Rscript analysis_peer_review.R
```

If no paths are supplied, the script uses the current working directory as the project directory and writes results to `results_peer_review/` within it.

The default settings reproduce the full manuscript analysis. Set `QUICK_MODE <- TRUE` only for debugging; it reduces the number of outer repetitions, bootstrap samples, and tuning values and must not be used to generate manuscript results.

## Implemented analysis

The analysis uses the first 7, 15, 30, 60, and 90 calendar days beginning with each participant's first passive-sensing record. A complete calendar grid prevents a missing day from shifting subsequent observations into an earlier window.

Participants must have at least 45 retained passive-data days and no more than 20% missing cells across passive features on recorded days. Within each observation window, the script calculates the mean, SD, and coverage for each passive and EMA variable. Mean and SD features require at least 70% coverage, and SD additionally requires at least three observations.

All preprocessing is estimated within the applicable training fold. Predictors without training data or variance are removed, remaining missing values are replaced with training-fold medians, and predictors are standardized using training-fold means and SDs.

The outer validation uses five repeats of stratified 5-fold participant-level cross-validation. Hyperparameters are selected by AUROC using 5-fold inner cross-validation. Classification thresholds are selected within the inner cross-validation to maximize F1 score. The same outer resamples are used across classifiers and observation windows.

The primary AUROC is calculated from each participant's out-of-fold probability averaged across the five outer repetitions. Its 95% confidence interval is obtained from 2,000 stratified participant-level bootstrap samples of these averaged predictions. Models are not refitted during bootstrap resampling, so these intervals are conditional on the fitted out-of-fold predictions.

## Principal outputs

The results directory includes:

- `Table2_participant_characteristics.csv` and `.docx`
- `Table3_model_performance.csv` and `.docx`
- `Table4_predictor_block_analysis.csv` and `.docx`
- `Table5_cross_validated_coefficient_stability.csv` and `.docx`
- `Supplementary_Table_S1_secondary_metrics.csv` and `.docx`
- `main_auroc_bootstrap_ci.csv`
- `main_performance_summary.csv`
- `block_auroc_bootstrap_ci.csv`
- `paired_block_auc_differences.csv`
- `predicted_probability_overlap.png`
- `calibration_primary_ridge.png`
- Feature dictionaries, cohort attrition, missingness summaries, consistency audits, and run manifests
- `input_checksums.csv`, `package_versions.csv`, and `session_info.txt`

The script verifies that the full 90-day ridge predictions and AUROC reported in Tables 3 and 4 are identical and that each reported predictor-block AUROC change equals the arithmetic difference between the corresponding model estimates.

## Participant-level outputs and confidentiality

`EXPORT_PARTICIPANT_LEVEL` is `FALSE` by default. With this setting, tables containing `uid` and the saved cross-validation resamples are not written. This is the recommended setting for a code package shared with reviewers.

Setting `EXPORT_PARTICIPANT_LEVEL <- TRUE` additionally writes participant-level predictions, fold assignments, coefficients, and resampling objects. Review these files under the study's data-sharing and institutional requirements before distributing them.

## Reproducibility notes

- The master random seed is `2025`.
- The full analysis uses 5 outer folds, 5 outer repetitions, 5 inner folds, and 2,000 bootstrap samples.
- Input-file MD5 checksums and software versions are recorded for every run.
- The repository does not recreate raw sensor processing or the exact software environment used before this reviewer-facing revision.
- The reviewer-facing script was revised by static inspection and must be run with the study data to verify full execution in the target R environment.

## Data availability

The data are not included because they contain sensitive participant-level mental health and digital phenotyping information. Access, if available, is governed by the study consent, Dartmouth College IRB protocol STUDY00032081, and the authors' data-use procedures.

## Contact

Questions about the analysis or access to supporting materials should be directed to the corresponding author identified in the manuscript.
