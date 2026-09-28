# MDD-only versus MDD superimposed on PDD: peer-review analysis code
# ============================================================================
# PURPOSE
# Participant-level classification using demographic, baseline questionnaire,
# passive-sensing, and momentary PHQ features. Outcome: pdd = 0 (MDD-only),
# pdd = 1 (MDD + PDD). This script starts from prepared CSVs; raw sensor
# processing and diagnostic ascertainment are outside its scope.
#
# RUNNING
# Requires R >= 4.1 (native pipe syntax) and the packages listed in section 1.
# Install those packages separately in R with install.packages(c(...)).
# Command line:
#   Rscript analysis_peer_review.R /path/to/project /path/to/new/results
# Or set MDD_PDD_PROJECT_DIR and MDD_PDD_OUTPUT_DIR environment variables.
# With no arguments, the current working directory is the project directory.
# Use a new output directory for each run; matching output files are overwritten.
#
# REQUIRED INPUTS (relative to the project directory)
#   final_withNA.csv: uid, date, numeric passive-sensing variables.
#   EMA_final.csv: uid, date, nine numeric PHQ items named PHQ1 ... PHQ9
#     (one optional period, underscore, or hyphen before each digit is allowed).
#   processed_csv/demographics.csv: uid and numerically coded demographics.
#   processed_csv/baseline_survey.csv: uid and baseline questionnaire scores.
# Metadata jointly supply pdd and the predictor names in section 3. Each metadata
# file must contain one row per participant, with consistent IDs across files.
# Dates must be YYYY-MM-DD. Missing values must be blank or NA.
# Verify demographic coding against the source data dictionary: Table 2 assumes
# sex 0=male/1=female; job 0=other/1=full-time; living_status 0=alone/1=with others;
# race 0=White/1=non-White; education 1:3 versus 4:8; income 1:6 versus 7:12.
# The last two splits are labelled high school or less and below $60,000.
# Numeric demographic codes enter the models as numeric predictors, not dummy
# variables. All numeric passive columns except tp/pdd are included: input CSVs
# must not contain other numeric identifiers or administrative variables.
#
# ANALYSIS SPECIFICATION (preserved from supplied code)
# - Windows: first 7, 15, 30, 60, 90 calendar days from first passive record.
# - Eligibility: >=45 recorded passive days and <=20% missing passive cells
#   among recorded rows; wholly absent dates do not enter that eligibility rate.
# - Features: mean, SD AND coverage for each passive and EMA variable/window.
#   Mean/SD require >=70% calendar-window coverage; SD also requires >=3 values.
# - Median imputation and standardization fitted within each training fold.
# - Five repeats of five-fold outer CV; five-fold inner tuning by AUROC.
# - Threshold selected by inner-CV F1; ridge, radial SVM and XGBoost grids below.
# - Main AUROC uses each participant's OOF probability averaged over repeats.
#   Its stratified percentile bootstrap resamples participants 2,000 times;
#   models are not refitted. Intervals are conditional on the fitted predictions
#   and do not capture all model-training or tuning uncertainty.
# - Secondary metrics are mean (SD) across repeats. Repeat quantiles in auxiliary
#   CSVs are descriptive intervals, not bootstrap confidence intervals.
# - Prevalence AUROC/AP are assigned theoretical chance values by convention.
# - Average precision retains the supplied rank-based calculation; tied scores
#   can make it depend on row ordering. It is not trapezoidal PR-AUC.
# - Coefficients describe penalized multivariable fits, not independent effects.
#
# OUTPUTS AND REPRODUCIBILITY
# Tables 2-5 and S1 (CSV/DOCX), plots, aggregate performance, feature dictionaries,
# audits, run settings, input checksums and session/package versions are saved.
# Participant-level prediction CSVs and resample RDS are opt-in below; these
# contain IDs and/or outcomes and are not automatically suitable for sharing.
# No participant data are bundled in this code. Data access arrangements must
# be supplied by the authors. Original package versions were not provided;
# environment records created here describe the new run, not the original run.
# Full manuscript settings are defaults; QUICK_MODE is for debugging only.
#
# REVIEW EDITS
# Removed personal path; added portable configuration, dependency checks,
# input validation, optional participant exports and provenance records.
# Fixed the prediction-identity audit when only the primary block window runs.
# No model grids, feature definitions, seeds or primary estimators were changed.
# Execution status: edited by static inspection; not run on the study data.
# ============================================================================

# ==============================================================================
# MDD-only vs MDD superimposed on PDD
# ==============================================================================

# ==============================================================================
# 0. Configuration
# ==============================================================================

# Positional arguments take precedence over environment variables.
args <- commandArgs(trailingOnly = TRUE)
project_arg <- if (length(args) >= 1L) args[[1L]] else
  Sys.getenv("MDD_PDD_PROJECT_DIR", unset = getwd())
PROJECT_DIR <- normalizePath(project_arg, mustWork = TRUE)

PASSIVE_FILE  <- file.path(PROJECT_DIR, "final_withNA.csv")
EMA_FILE      <- file.path(PROJECT_DIR, "EMA_final.csv")
DEMO_FILE     <- file.path(PROJECT_DIR, "processed_csv", "demographics.csv")
BASELINE_FILE <- file.path(PROJECT_DIR, "processed_csv", "baseline_survey.csv")

ANALYSIS_VERSION <- "2026-08-17-v6-STABLE-TABLE5-90DAY"
OUTPUT_DIR <- if (length(args) >= 2L) args[[2L]] else
  Sys.getenv("MDD_PDD_OUTPUT_DIR", unset = file.path(PROJECT_DIR, "results_peer_review"))
EXPORT_PARTICIPANT_LEVEL <- FALSE

WINDOWS <- c(7L, 15L, 30L, 60L, 90L)
PRIMARY_WINDOW <- 90L  # Used for the block analysis and descriptive coefficients.

if (!identical(PRIMARY_WINDOW, 90L)) {
  stop(
    "This manuscript pipeline is locked to a 90-day primary analysis. ",
    "Set PRIMARY_WINDOW <- 90L before running it."
  )
}


EMA_PHQ_PATTERN <- "^PHQ[._-]?[1-9]$"

MIN_RETAINED_DAYS <- 45L
MAX_PASSIVE_CELL_MISSING <- 0.20
MIN_WINDOW_COVERAGE <- 0.70
MIN_OBS_FOR_SD <- 3L

OUTER_FOLDS <- 5L
OUTER_REPEATS <- 5L
INNER_FOLDS <- 5L
BOOTSTRAP_REPS <- 2000L
MASTER_SEED <- 2025L

# Keep a model's inner folds and stochastic fit identical wherever that model
# is run. In particular, the full 90-day ridge model in Table 3 must be the
# same fit used as the reference model in Table 4.
MODEL_SEED_OFFSET <- c(
  Prevalence = 1L,
  Ridge = 2L,
  SVM = 3L,
  XGBoost = 4L
)

RUN_NONLINEAR_MODELS <- TRUE
RUN_BLOCKS_FOR_ALL_WINDOWS <- TRUE
QUICK_MODE <- FALSE
VERBOSE <- FALSE

if (QUICK_MODE) {
  OUTER_REPEATS <- 1L
  BOOTSTRAP_REPS <- 200L
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. Packages and general helpers
# ==============================================================================

required_packages <- c(
  "dplyr", "tidyr", "tibble", "purrr", "readr", "rsample", "glmnet",
  "e1071", "xgboost", "pROC", "ggplot2", "flextable", "officer"
)


if (getRversion() < "4.1.0") stop("R 4.1.0 or later is required.")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)
)]
if (length(missing_packages)) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
readr::write_csv(
  data.frame(package = required_packages, version = vapply(
    required_packages, function(x) as.character(utils::packageVersion(x)), character(1)
  )), file.path(OUTPUT_DIR, "package_versions.csv")
)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
  library(rsample)
  library(glmnet)
  library(e1071)
  library(xgboost)
  library(pROC)
  library(ggplot2)
})

set.seed(MASTER_SEED)

assert_file <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
}

assert_columns <- function(data, columns, label) {
  absent <- setdiff(columns, names(data))
  if (length(absent) > 0L) {
    stop(label, " is missing required column(s): ", paste(absent, collapse = ", "))
  }
}

drop_csv_index <- function(data) {
  index_cols <- grep("^(X|[.][.][.]1|row.names)$", names(data), value = TRUE)
  dplyr::select(data, -dplyr::any_of(index_cols))
}

read_input_csv <- function(path) {
  assert_file(path)
  readr::read_csv(path, show_col_types = FALSE) |>
    drop_csv_index()
}

first_nonmissing <- function(x) {
  observed <- x[!is.na(x)]
  if (length(observed) == 0L) NA else observed[[1L]]
}

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

parse_date_strict <- function(x, label) {
  if (inherits(x, "Date")) return(x)
  out <- as.Date(x)
  failed <- is.na(out) & !is.na(x) & nzchar(as.character(x))
  if (any(failed)) {
    stop(
      "Could not parse ", sum(failed), " date value(s) in ", label,
      ". Convert the date column to YYYY-MM-DD before running this script."
    )
  }
  out
}

aggregate_to_one_row_per_day <- function(data, label) {
  assert_columns(data, c("uid", "date"), label)
  data$date <- parse_date_strict(data$date, label)
  if (anyNA(data$uid) || anyNA(data$date)) {
    stop(label, ": uid and date must not be missing.")
  }
  
  duplicate_n <- data |>
    count(uid, date, name = "n") |>
    filter(n > 1L) |>
    nrow()
  
  if (VERBOSE && duplicate_n > 0L) {
    message(label, ": aggregating ", duplicate_n, " duplicated participant-day(s).")
  }
  
  value_cols <- setdiff(names(data), c("uid", "date"))
  numeric_cols <- value_cols[vapply(data[value_cols], is.numeric, logical(1))]
  other_cols <- setdiff(value_cols, numeric_cols)
  
  data |>
    group_by(uid, date) |>
    summarise(
      across(all_of(numeric_cols), mean_or_na),
      across(all_of(other_cols), first_nonmissing),
      .groups = "drop"
    )
}

normalise_metadata_names <- function(data) {
  names(data) <- tolower(names(data))
  data <- drop_csv_index(data)
  
  rename_map <- c(
    bfi_ex_mean = "bfi_ex",
    bfi_ne_mean = "bfi_ne",
    bfi_con_mean = "bfi_con"
  )
  
  for (old_name in names(rename_map)) {
    new_name <- unname(rename_map[[old_name]])
    if (old_name %in% names(data) && !new_name %in% names(data)) {
      names(data)[names(data) == old_name] <- new_name
    }
  }
  data
}

combine_participant_tables <- function(demo, baseline) {
  demo <- normalise_metadata_names(demo)
  baseline <- normalise_metadata_names(baseline)
  assert_columns(demo, "uid", "Demographic data")
  assert_columns(baseline, "uid", "Baseline survey data")
  
  overlap <- setdiff(intersect(names(demo), names(baseline)), "uid")
  joined <- full_join(demo, baseline, by = "uid", suffix = c("__demo", "__base"))
  
  for (column in overlap) {
    demo_col <- paste0(column, "__demo")
    base_col <- paste0(column, "__base")
    conflict <- !is.na(joined[[demo_col]]) & !is.na(joined[[base_col]]) &
      as.character(joined[[demo_col]]) != as.character(joined[[base_col]])
    
    if (any(conflict)) {
      stop("Conflicting values for '", column, "' in demographic and baseline files.")
    }
    
    joined[[column]] <- dplyr::coalesce(joined[[demo_col]], joined[[base_col]])
    joined[[demo_col]] <- NULL
    joined[[base_col]] <- NULL
  }
  
  joined
}

write_csv_safe <- function(data, filename) {
  if (!EXPORT_PARTICIPANT_LEVEL && "uid" %in% names(data)) {
    if (VERBOSE) message("Participant-level export disabled: ", filename)
    return(invisible(NULL))
  }
  readr::write_csv(data, file.path(OUTPUT_DIR, filename), na = "")
}

# ==============================================================================
# 2. Load, validate, and align daily data
# ==============================================================================

passive_raw <- read_input_csv(PASSIVE_FILE)
ema_raw <- read_input_csv(EMA_FILE)
demo_raw <- read_input_csv(DEMO_FILE)
baseline_raw <- read_input_csv(BASELINE_FILE)

# Reject ambiguous metadata instead of silently selecting the first record.
for (metadata in list(demo_raw, baseline_raw)) {
  id_name <- names(metadata)[tolower(names(metadata)) == "uid"]
  if (length(id_name) != 1L) stop("Each metadata file must have one uid column.")
  ids <- metadata[[id_name]]
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop("Metadata uid values must be nonmissing and unique within each file.")
  }
}
input_paths <- c(PASSIVE_FILE, EMA_FILE, DEMO_FILE, BASELINE_FILE)
write_csv_safe(tibble(
  input_file = c("final_withNA.csv", "EMA_final.csv",
                 "processed_csv/demographics.csv", "processed_csv/baseline_survey.csv"),
  md5 = unname(tools::md5sum(input_paths))
), "input_checksums.csv")

names(passive_raw) <- tolower(names(passive_raw))
assert_columns(passive_raw, c("uid", "date"), "Passive daily data")

passive_raw <- passive_raw |>
  select(-any_of(c("tp", "pdd")))

passive_daily <- aggregate_to_one_row_per_day(passive_raw, "Passive daily data")

# EMA source names are preserved after the uid/date normalization and then
# prefixed, preventing confusion with lowercase baseline PHQ/GAD variables.
names(ema_raw)[tolower(names(ema_raw)) == "uid"] <- "uid"
names(ema_raw)[tolower(names(ema_raw)) == "date"] <- "date"
ema_exclude_cols <- names(ema_raw)[tolower(names(ema_raw)) %in% c("tp", "pdd")]
ema_raw <- ema_raw |>
  select(-any_of(ema_exclude_cols))
ema_daily <- aggregate_to_one_row_per_day(ema_raw, "EMA daily data")

passive_value_cols <- setdiff(names(passive_daily), c("uid", "date"))
ema_value_cols <- setdiff(names(ema_daily), c("uid", "date"))

# Explicitly limit the EMA block to momentary PHQ items 1 through 9. This both
# documents the intended model and prevents EMA GAD, Somatic, and numeric
# administrative fields from becoming predictors.
ema_phq_cols <- ema_value_cols[
  grepl(EMA_PHQ_PATTERN, ema_value_cols, ignore.case = TRUE)
]

if (length(ema_phq_cols) == 0L) {
  stop(
    "No EMA PHQ items matched EMA_PHQ_PATTERN = '",
    EMA_PHQ_PATTERN,
    "'. Check the column names in EMA_final.csv."
  )
}

if (length(ema_phq_cols) != 9L) {
  stop(
    "Expected exactly 9 EMA PHQ items but found ",
    length(ema_phq_cols),
    ": ",
    paste(ema_phq_cols, collapse = ", ")
  )
}

# Extract the terminal item digit and verify complete item coverage.
ema_phq_item_number <- as.integer(sub(".*([1-9])$", "\\1", ema_phq_cols))
if (!identical(sort(unique(ema_phq_item_number)), 1:9)) {
  stop(
    "EMA PHQ columns must contain each item from PHQ-1 through PHQ-9 exactly once: ",
    paste(ema_phq_cols, collapse = ", ")
  )
}

unexpected_ema_cols <- setdiff(ema_value_cols, ema_phq_cols)
if (VERBOSE && length(unexpected_ema_cols) > 0L) {
  message(
    "EMA columns excluded because they are not PHQ-9 items: ",
    paste(unexpected_ema_cols, collapse = ", ")
  )
}

passive_numeric <- passive_value_cols[
  vapply(passive_daily[passive_value_cols], is.numeric, logical(1))
]
ema_numeric <- ema_phq_cols[
  vapply(ema_daily[ema_phq_cols], is.numeric, logical(1))
]

if (length(passive_numeric) == 0L) stop("No numeric passive-sensing features found.")
if (length(ema_numeric) == 0L) {
  stop("The matched EMA PHQ items are not numeric.")
}

non_numeric_ema <- setdiff(ema_phq_cols, ema_numeric)
if (length(non_numeric_ema) > 0L) {
  stop(
    "These EMA PHQ items are not numeric: ",
    paste(non_numeric_ema, collapse = ", ")
  )
}

passive_daily <- passive_daily |>
  select(uid, date, all_of(passive_numeric)) |>
  rename_with(~ paste0("ps__", .x), all_of(passive_numeric))

write_csv_safe(
  tibble(
    source_variable = passive_numeric,
    predictor_block = "Passive sensing",
    used_for_participant_missingness = TRUE,
    used_for_window_features = TRUE
  ),
  "passive_variables_included.csv"
)

ema_daily <- ema_daily |>
  select(uid, date, all_of(ema_numeric)) |>
  rename_with(~ paste0("ema__", .x), all_of(ema_numeric))

write_csv_safe(
  tibble(
    source_variable = ema_numeric,
    predictor_block = "EMA",
    inclusion_rule = "Exactly one of the nine momentary PHQ items"
  ),
  "ema_variables_included.csv"
)

passive_cols_daily <- grep("^ps__", names(passive_daily), value = TRUE)
ema_cols_daily <- grep("^ema__", names(ema_daily), value = TRUE)

# Reproduce the participant-level daily-data eligibility checks before any
# modeling. If final_withNA.csv was already filtered, this should be idempotent.
passive_qc <- passive_daily |>
  group_by(uid) |>
  group_modify(function(.x, .y) {
    feature_matrix <- as.matrix(.x[, passive_cols_daily, drop = FALSE])
    tibble(
      retained_days = n_distinct(.x$date),
      passive_cell_missing = mean(is.na(feature_matrix))
    )
  }) |>
  ungroup()

eligible_sensor_uids <- passive_qc |>
  filter(
    retained_days >= MIN_RETAINED_DAYS,
    passive_cell_missing <= MAX_PASSIVE_CELL_MISSING
  ) |>
  pull(uid)

participant_data <- combine_participant_tables(demo_raw, baseline_raw) |>
  select(-any_of(c(
    "ders_awareness", "ders_clarity", "ders_goals", "ders_impulse",
    "ders_nonacceptance", "ders_strategies"
  )))

assert_columns(participant_data, c("uid", "pdd"), "Participant metadata")

participant_data <- participant_data |>
  mutate(pdd = as.integer(pdd)) |>
  filter(!is.na(pdd), pdd %in% c(0L, 1L)) |>
  distinct(uid, .keep_all = TRUE)

analytic_uids <- intersect(eligible_sensor_uids, participant_data$uid)

if (length(analytic_uids) == 0L) stop("No eligible participants remain after merging.")

passive_daily <- passive_daily |> filter(uid %in% analytic_uids)
ema_daily <- ema_daily |> filter(uid %in% analytic_uids)
participant_data <- participant_data |> filter(uid %in% analytic_uids)

# Study day is based on elapsed calendar days, not row number. A complete
# calendar grid ensures that a missing day does not shift later records into an
# earlier observation window.
study_starts <- passive_daily |>
  group_by(uid) |>
  summarise(study_start = min(date, na.rm = TRUE), .groups = "drop")

daily_grid <- study_starts[
  rep(seq_len(nrow(study_starts)), each = max(WINDOWS)),
  ,
  drop = FALSE
] |>
  mutate(study_day = rep(seq_len(max(WINDOWS)), times = nrow(study_starts))) |>
  mutate(date = study_start + study_day - 1L) |>
  select(uid, date, study_day)

daily_aligned <- daily_grid |>
  left_join(passive_daily, by = c("uid", "date")) |>
  left_join(ema_daily, by = c("uid", "date"))

attrition <- tibble(
  stage = c(
    "Participants in passive daily file",
    paste0("At least ", MIN_RETAINED_DAYS, " retained days"),
    paste0("Passive cell missingness <= ", MAX_PASSIVE_CELL_MISSING),
    "Outcome available and metadata merged"
  ),
  n = c(
    n_distinct(passive_raw$uid),
    sum(passive_qc$retained_days >= MIN_RETAINED_DAYS),
    length(eligible_sensor_uids),
    length(analytic_uids)
  )
)
write_csv_safe(attrition, "cohort_attrition.csv")

# ==============================================================================
# 3. Construct observed-day window summaries and explicit feature blocks
# ==============================================================================

coverage <- function(x) mean(!is.na(x))

covered_mean <- function(x) {
  if (coverage(x) < MIN_WINDOW_COVERAGE) NA_real_ else mean_or_na(x)
}

covered_sd <- function(x) {
  if (coverage(x) < MIN_WINDOW_COVERAGE || sum(!is.na(x)) < MIN_OBS_FOR_SD) {
    NA_real_
  } else {
    stats::sd(x, na.rm = TRUE)
  }
}

summarise_window <- function(data, window) {
  window_data <- data |> filter(study_day <= window)
  
  window_data |>
    group_by(uid) |>
    summarise(
      across(
        all_of(c(passive_cols_daily, ema_cols_daily)),
        list(mean = covered_mean, sd = covered_sd, coverage = coverage),
        .names = paste0("{.col}__{.fn}__w", window)
      ),
      .groups = "drop"
    )
}

window_summaries <- purrr::map(WINDOWS, ~ summarise_window(daily_aligned, .x)) |>
  purrr::reduce(full_join, by = "uid")

analysis_data <- participant_data |>
  inner_join(window_summaries, by = "uid")

demo_cols <- intersect(
  c("age", "sex", "job", "education", "living_status", "race", "income"),
  names(analysis_data)
)

baseline_cols <- intersect(
  c(
    "ace", "bfi_ex", "bfi_ne", "bfi_con", "iip", "phq9", "phq15",
    "gad", "audit", "ders_total"
  ),
  names(analysis_data)
)

if (length(demo_cols) == 0L) stop("No demographic predictors were found.")
if (length(baseline_cols) == 0L) stop("No baseline survey predictors were found.")

passive_window_cols <- function(window) {
  grep(paste0("^ps__.*__w", window, "$"), names(analysis_data), value = TRUE)
}

ema_window_cols <- function(window) {
  grep(paste0("^ema__.*__w", window, "$"), names(analysis_data), value = TRUE)
}

features_for_blocks <- function(window, blocks) {
  block_map <- list(
    demographics = demo_cols,
    baseline = baseline_cols,
    passive = passive_window_cols(window),
    ema = ema_window_cols(window)
  )
  unknown <- setdiff(blocks, names(block_map))
  if (length(unknown) > 0L) stop("Unknown feature block: ", paste(unknown, collapse = ", "))
  unique(unlist(block_map[blocks], use.names = FALSE))
}

all_predictor_cols <- unique(c(
  demo_cols,
  baseline_cols,
  unlist(lapply(WINDOWS, passive_window_cols), use.names = FALSE),
  unlist(lapply(WINDOWS, ema_window_cols), use.names = FALSE)
))

non_numeric_predictors <- all_predictor_cols[
  !vapply(analysis_data[all_predictor_cols], is.numeric, logical(1))
]

if (length(non_numeric_predictors) > 0L) {
  stop(
    "Predictors must use numeric coding. Convert these columns before modeling: ",
    paste(non_numeric_predictors, collapse = ", ")
  )
}

write_csv_safe(
  tibble(
    block = c(
      rep("Demographics", length(demo_cols)),
      rep("Baseline", length(baseline_cols)),
      rep("Passive example (primary window)", length(passive_window_cols(PRIMARY_WINDOW))),
      rep("EMA example (primary window)", length(ema_window_cols(PRIMARY_WINDOW)))
    ),
    predictor = c(
      demo_cols,
      baseline_cols,
      passive_window_cols(PRIMARY_WINDOW),
      ema_window_cols(PRIMARY_WINDOW)
    )
  ),
  "predictor_dictionary.csv"
)

# Exact requested feature counts before fold-specific removal of all-missing or
# zero-variance predictors. These counts should be reported in the manuscript;
# the actual usable count can vary slightly by training fold.
model_feature_counts <- purrr::map_dfr(WINDOWS, function(window) {
  tibble(
    window = window,
    demographics = length(demo_cols),
    baseline = length(baseline_cols),
    passive = length(passive_window_cols(window)),
    ema = length(ema_window_cols(window)),
    full_model = length(features_for_blocks(
      window,
      c("demographics", "baseline", "passive", "ema")
    ))
  )
})

write_csv_safe(model_feature_counts, "model_feature_counts.csv")

feature_missingness <- tibble(
  predictor = all_predictor_cols,
  missing_fraction = vapply(
    analysis_data[all_predictor_cols],
    function(x) mean(is.na(x)),
    numeric(1)
  )
) |>
  arrange(desc(missing_fraction))

write_csv_safe(feature_missingness, "participant_feature_missingness.csv")

# ==============================================================================
# 4. Training-fold-only preprocessing and metrics
# ==============================================================================

fit_preprocessor <- function(train, test, feature_cols) {
  x_train <- as.matrix(train[, feature_cols, drop = FALSE])
  x_test <- as.matrix(test[, feature_cols, drop = FALSE])
  storage.mode(x_train) <- "double"
  storage.mode(x_test) <- "double"
  
  x_train[!is.finite(x_train)] <- NA_real_
  x_test[!is.finite(x_test)] <- NA_real_
  
  has_training_data <- colSums(!is.na(x_train)) > 0L
  x_train <- x_train[, has_training_data, drop = FALSE]
  x_test <- x_test[, has_training_data, drop = FALSE]
  
  if (ncol(x_train) == 0L) stop("No usable predictors in a training fold.")
  
  medians <- apply(x_train, 2, median, na.rm = TRUE)
  for (j in seq_len(ncol(x_train))) {
    x_train[is.na(x_train[, j]), j] <- medians[[j]]
    x_test[is.na(x_test[, j]), j] <- medians[[j]]
  }
  
  sds_before_scaling <- apply(x_train, 2, stats::sd)
  keep <- is.finite(sds_before_scaling) & sds_before_scaling > 0
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- x_test[, keep, drop = FALSE]
  
  if (ncol(x_train) == 0L) stop("All predictors had zero variance in a training fold.")
  
  means <- colMeans(x_train)
  sds <- apply(x_train, 2, stats::sd)
  x_train <- sweep(sweep(x_train, 2, means, "-"), 2, sds, "/")
  x_test <- sweep(sweep(x_test, 2, means, "-"), 2, sds, "/")
  
  list(train = x_train, test = x_test, feature_names = colnames(x_train))
}

make_stratified_fold_id <- function(y, v, seed) {
  set.seed(seed)
  y <- as.integer(y)
  if (anyNA(y) || length(unique(y)) != 2L || min(table(y)) < v) {
    stop("Each outcome class must have at least v participants for inner CV.")
  }
  fold_id <- integer(length(y))
  for (class_value in sort(unique(y))) {
    indices <- sample(which(y == class_value))
    fold_id[indices] <- rep(seq_len(v), length.out = length(indices))
  }
  fold_id
}

safe_auc <- function(truth, probability) {
  valid <- is.finite(probability) & !is.na(truth)
  truth <- truth[valid]
  probability <- probability[valid]
  if (length(unique(truth)) < 2L) return(NA_real_)
  if (length(unique(probability)) < 2L) return(0.5)
  
  as.numeric(
    pROC::auc(
      pROC::roc(
        response = truth,
        predictor = probability,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      )
    )
  )
}

average_precision <- function(truth, probability) {
  order_index <- order(probability, decreasing = TRUE)
  y <- truth[order_index]
  positive_n <- sum(y == 1L)
  if (positive_n == 0L) return(NA_real_)
  precision_at_rank <- cumsum(y == 1L) / seq_along(y)
  mean(precision_at_rank[y == 1L])
}

best_f1_threshold <- function(probability, truth) {
  thresholds <- seq(0.01, 0.99, by = 0.01)
  f1 <- vapply(thresholds, function(threshold) {
    predicted <- as.integer(probability >= threshold)
    tp <- sum(predicted == 1L & truth == 1L)
    fp <- sum(predicted == 1L & truth == 0L)
    fn <- sum(predicted == 0L & truth == 1L)
    denominator <- 2 * tp + fp + fn
    if (denominator == 0L) 0 else 2 * tp / denominator
  }, numeric(1))
  
  best <- which(f1 == max(f1, na.rm = TRUE))
  # Use the largest tied threshold to avoid unnecessarily labeling everyone positive.
  thresholds[max(best)]
}

prediction_metrics <- function(truth, probability, predicted) {
  valid <- !is.na(truth) & is.finite(probability) & !is.na(predicted)
  truth <- as.integer(truth[valid])
  probability <- probability[valid]
  predicted <- as.integer(predicted[valid])
  
  tp <- sum(predicted == 1L & truth == 1L)
  tn <- sum(predicted == 0L & truth == 0L)
  fp <- sum(predicted == 1L & truth == 0L)
  fn <- sum(predicted == 0L & truth == 1L)
  
  sensitivity <- if ((tp + fn) == 0L) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0L) NA_real_ else tn / (tn + fp)
  precision <- if ((tp + fp) == 0L) NA_real_ else tp / (tp + fp)
  f1 <- if ((2 * tp + fp + fn) == 0L) 0 else 2 * tp / (2 * tp + fp + fn)
  
  clipped <- pmin(pmax(probability, 1e-6), 1 - 1e-6)
  calibration_fit <- tryCatch(
    glm(truth ~ qlogis(clipped), family = binomial()),
    error = function(e) NULL
  )
  
  tibble(
    n = length(truth),
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn,
    prevalence = mean(truth),
    auroc = safe_auc(truth, probability),
    average_precision = average_precision(truth, probability),
    accuracy = (tp + tn) / length(truth),
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1 = f1,
    brier = mean((probability - truth)^2),
    calibration_intercept = if (is.null(calibration_fit)) NA_real_ else coef(calibration_fit)[[1]],
    calibration_slope = if (is.null(calibration_fit)) NA_real_ else coef(calibration_fit)[[2]]
  )
}

# ==============================================================================
# 5. Nested model fitting
# ==============================================================================

ridge_lambda_grid <- sort(10^seq(-4, 2, length.out = if (QUICK_MODE) 15L else 40L), decreasing = TRUE)

svm_grid <- if (QUICK_MODE) {
  expand.grid(cost = c(0.5, 2), gamma = c(0.01, 0.1))
} else {
  expand.grid(cost = c(0.25, 1, 4, 16), gamma = c(0.005, 0.02, 0.08))
}

xgb_grid <- if (QUICK_MODE) {
  expand.grid(
    nrounds = c(50L, 100L), max_depth = 2L, eta = 0.10,
    min_child_weight = 1, subsample = 0.8, colsample_bytree = 0.8
  )
} else {
  expand.grid(
    nrounds = c(50L, 100L),
    max_depth = c(2L, 3L),
    eta = c(0.03, 0.10),
    min_child_weight = 1,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
}

fit_ridge_nested <- function(train, test, feature_cols, seed) {
  y <- train$pdd
  inner_id <- make_stratified_fold_id(y, INNER_FOLDS, seed)
  oof <- matrix(NA_real_, nrow(train), length(ridge_lambda_grid))
  
  for (fold in seq_len(INNER_FOLDS)) {
    inner_train <- train[inner_id != fold, , drop = FALSE]
    inner_valid <- train[inner_id == fold, , drop = FALSE]
    prep <- fit_preprocessor(inner_train, inner_valid, feature_cols)
    
    fit <- glmnet::glmnet(
      prep$train,
      inner_train$pdd,
      family = "binomial",
      alpha = 0,
      lambda = ridge_lambda_grid,
      standardize = FALSE
    )
    
    oof[inner_id == fold, ] <- as.matrix(
      predict(fit, newx = prep$test, s = ridge_lambda_grid, type = "response")
    )
  }
  
  lambda_auc <- apply(oof, 2, function(p) safe_auc(y, p))
  if (all(is.na(lambda_auc))) stop("Ridge inner CV produced no valid AUROC.")
  best_index <- which(lambda_auc == max(lambda_auc, na.rm = TRUE))[1L]
  best_lambda <- ridge_lambda_grid[[best_index]]
  threshold <- best_f1_threshold(oof[, best_index], y)
  
  prep_final <- fit_preprocessor(train, test, feature_cols)
  final_fit <- glmnet::glmnet(
    prep_final$train,
    train$pdd,
    family = "binomial",
    alpha = 0,
    lambda = best_lambda,
    standardize = FALSE
  )
  
  probability <- as.numeric(
    predict(final_fit, newx = prep_final$test, type = "response")
  )
  coefficient_matrix <- as.matrix(coef(final_fit))
  
  list(
    probability = probability,
    predicted = as.integer(probability >= threshold),
    threshold = threshold,
    best_lambda = best_lambda,
    n_features_used = length(prep_final$feature_names),
    tuning = paste0("lambda=", signif(best_lambda, 5)),
    coefficients = tibble(
      predictor = rownames(coefficient_matrix),
      coefficient = as.numeric(coefficient_matrix[, 1])
    ) |>
      filter(predictor != "(Intercept)")
  )
}

extract_svm_case_probability <- function(fit, newdata) {
  prediction <- predict(fit, newdata, probability = TRUE)
  probabilities <- attr(prediction, "probabilities")
  if (is.null(probabilities) || !"Case" %in% colnames(probabilities)) {
    stop("SVM did not return Case probabilities.")
  }
  as.numeric(probabilities[, "Case"])
}

fit_svm_nested <- function(train, test, feature_cols, seed) {
  y <- train$pdd
  inner_id <- make_stratified_fold_id(y, INNER_FOLDS, seed)
  oof <- matrix(NA_real_, nrow(train), nrow(svm_grid))
  
  for (fold in seq_len(INNER_FOLDS)) {
    inner_train <- train[inner_id != fold, , drop = FALSE]
    inner_valid <- train[inner_id == fold, , drop = FALSE]
    prep <- fit_preprocessor(inner_train, inner_valid, feature_cols)
    y_factor <- factor(inner_train$pdd, levels = c(0, 1), labels = c("Control", "Case"))
    
    for (grid_index in seq_len(nrow(svm_grid))) {
      fit <- e1071::svm(
        x = prep$train,
        y = y_factor,
        kernel = "radial",
        cost = svm_grid$cost[[grid_index]],
        gamma = svm_grid$gamma[[grid_index]],
        probability = TRUE,
        scale = FALSE
      )
      oof[inner_id == fold, grid_index] <- extract_svm_case_probability(fit, prep$test)
    }
  }
  
  grid_auc <- apply(oof, 2, function(p) safe_auc(y, p))
  if (all(is.na(grid_auc))) stop("SVM inner CV produced no valid AUROC.")
  best_index <- which(grid_auc == max(grid_auc, na.rm = TRUE))[1L]
  threshold <- best_f1_threshold(oof[, best_index], y)
  
  prep_final <- fit_preprocessor(train, test, feature_cols)
  final_fit <- e1071::svm(
    x = prep_final$train,
    y = factor(train$pdd, levels = c(0, 1), labels = c("Control", "Case")),
    kernel = "radial",
    cost = svm_grid$cost[[best_index]],
    gamma = svm_grid$gamma[[best_index]],
    probability = TRUE,
    scale = FALSE
  )
  probability <- extract_svm_case_probability(final_fit, prep_final$test)
  
  list(
    probability = probability,
    predicted = as.integer(probability >= threshold),
    threshold = threshold,
    n_features_used = length(prep_final$feature_names),
    tuning = paste0(
      "cost=", svm_grid$cost[[best_index]],
      ";gamma=", svm_grid$gamma[[best_index]]
    )
  )
}

fit_one_xgb <- function(x, y, parameters, seed) {
  xgboost::xgb.train(
    params = list(
      objective = "binary:logistic",
      eval_metric = "auc",
      max_depth = as.integer(parameters$max_depth),
      eta = parameters$eta,
      min_child_weight = parameters$min_child_weight,
      subsample = parameters$subsample,
      colsample_bytree = parameters$colsample_bytree,
      seed = as.integer(seed),
      nthread = 1
    ),
    data = xgboost::xgb.DMatrix(x, label = y),
    nrounds = as.integer(parameters$nrounds),
    verbose = 0
  )
}

fit_xgb_nested <- function(train, test, feature_cols, seed) {
  y <- train$pdd
  inner_id <- make_stratified_fold_id(y, INNER_FOLDS, seed)
  oof <- matrix(NA_real_, nrow(train), nrow(xgb_grid))
  
  for (fold in seq_len(INNER_FOLDS)) {
    inner_train <- train[inner_id != fold, , drop = FALSE]
    inner_valid <- train[inner_id == fold, , drop = FALSE]
    prep <- fit_preprocessor(inner_train, inner_valid, feature_cols)
    
    for (grid_index in seq_len(nrow(xgb_grid))) {
      fit <- fit_one_xgb(
        prep$train,
        inner_train$pdd,
        xgb_grid[grid_index, ],
        seed = seed + 1000L * fold + grid_index
      )
      oof[inner_id == fold, grid_index] <- predict(
        fit,
        xgboost::xgb.DMatrix(prep$test)
      )
    }
  }
  
  grid_auc <- apply(oof, 2, function(p) safe_auc(y, p))
  if (all(is.na(grid_auc))) stop("XGBoost inner CV produced no valid AUROC.")
  best_index <- which(grid_auc == max(grid_auc, na.rm = TRUE))[1L]
  threshold <- best_f1_threshold(oof[, best_index], y)
  
  prep_final <- fit_preprocessor(train, test, feature_cols)
  final_fit <- fit_one_xgb(
    prep_final$train,
    train$pdd,
    xgb_grid[best_index, ],
    seed = seed + 99999L
  )
  probability <- predict(final_fit, xgboost::xgb.DMatrix(prep_final$test))
  
  list(
    probability = as.numeric(probability),
    predicted = as.integer(probability >= threshold),
    threshold = threshold,
    n_features_used = length(prep_final$feature_names),
    tuning = paste(
      paste(names(xgb_grid), unlist(xgb_grid[best_index, ]), sep = "="),
      collapse = ";"
    )
  )
}

fit_prevalence_baseline <- function(train, test) {
  probability <- rep(mean(train$pdd), nrow(test))
  predicted_class <- as.integer(mean(train$pdd) >= 0.5)
  list(
    probability = probability,
    predicted = rep(predicted_class, nrow(test)),
    threshold = 0.5,
    n_features_used = 0L,
    tuning = "training-fold prevalence"
  )
}

fit_nested_model <- function(model, train, test, feature_cols, seed) {
  switch(
    model,
    Prevalence = fit_prevalence_baseline(train, test),
    Ridge = fit_ridge_nested(train, test, feature_cols, seed),
    SVM = fit_svm_nested(train, test, feature_cols, seed),
    XGBoost = fit_xgb_nested(train, test, feature_cols, seed),
    stop("Unknown model: ", model)
  )
}

# ==============================================================================
# 6. Repeated outer cross-validation
# ==============================================================================

id_data <- analysis_data |> select(uid, pdd)
if (length(unique(id_data$pdd)) != 2L || min(table(id_data$pdd)) < OUTER_FOLDS) {
  stop("Both outcome classes require at least OUTER_FOLDS participants.")
}

set.seed(MASTER_SEED)
outer_resamples <- rsample::vfold_cv(
  id_data,
  v = OUTER_FOLDS,
  repeats = OUTER_REPEATS,
  strata = pdd
)

resample_labels <- function(resamples, index) {
  if ("id2" %in% names(resamples)) {
    list(repeat_id = resamples$id[[index]], fold = resamples$id2[[index]])
  } else {
    list(repeat_id = "Repeat1", fold = resamples$id[[index]])
  }
}

run_outer_condition <- function(window, condition, blocks, models, resamples) {
  feature_cols <- features_for_blocks(window, blocks)
  prediction_output <- vector("list", nrow(resamples) * length(models))
  coefficient_output <- vector("list", nrow(resamples) * sum(models == "Ridge"))
  prediction_index <- 1L
  coefficient_index <- 1L
  
  for (resample_index in seq_len(nrow(resamples))) {
    labels <- resample_labels(resamples, resample_index)
    train_ids <- rsample::analysis(resamples$splits[[resample_index]])$uid
    test_ids <- rsample::assessment(resamples$splits[[resample_index]])$uid
    train <- analysis_data |> filter(uid %in% train_ids)
    test <- analysis_data |> filter(uid %in% test_ids)
    
    for (model in models) {
      if (!model %in% names(MODEL_SEED_OFFSET)) {
        stop("No stable seed offset has been defined for model: ", model)
      }
      seed <- MASTER_SEED +
        10000L * window +
        100L * resample_index +
        unname(MODEL_SEED_OFFSET[[model]])
      if (VERBOSE) {
        message(
          "Window ", window, ", ", condition, ", ", labels$repeat_id,
          "/", labels$fold, ", ", model
        )
      }
      
      fitted <- fit_nested_model(model, train, test, feature_cols, seed)
      prediction_output[[prediction_index]] <- tibble(
        uid = test$uid,
        truth = test$pdd,
        probability = fitted$probability,
        predicted = fitted$predicted,
        threshold = fitted$threshold,
        model = model,
        condition = condition,
        window = window,
        repeat_id = labels$repeat_id,
        fold = labels$fold,
        tuning = fitted$tuning,
        n_features_requested = length(feature_cols),
        n_features_used = fitted$n_features_used
      )
      prediction_index <- prediction_index + 1L
      
      if (model == "Ridge") {
        coefficient_output[[coefficient_index]] <- fitted$coefficients |>
          transmute(
            condition = .env$condition,
            window = .env$window,
            repeat_id = labels$repeat_id,
            fold = labels$fold,
            outer_fit_id = paste(labels$repeat_id, labels$fold, sep = "_"),
            predictor,
            coefficient,
            absolute_coefficient = abs(coefficient),
            best_lambda = fitted$best_lambda,
            n_features_requested = length(feature_cols),
            n_features_used = fitted$n_features_used
          )
        coefficient_index <- coefficient_index + 1L
      }
    }
  }
  
  list(
    predictions = bind_rows(prediction_output),
    coefficients = bind_rows(coefficient_output)
  )
}

main_models <- c("Prevalence", "Ridge")
if (RUN_NONLINEAR_MODELS) main_models <- c(main_models, "SVM", "XGBoost")

main_outer_results <- purrr::map(WINDOWS, function(window) {
  run_outer_condition(
    window = window,
    condition = "Full model including EMA",
    blocks = c("demographics", "baseline", "passive", "ema"),
    models = main_models,
    resamples = outer_resamples
  )
})

main_predictions <- purrr::map_dfr(main_outer_results, "predictions")
main_ridge_coefficients <- purrr::map_dfr(main_outer_results, "coefficients")

write_csv_safe(main_predictions, "main_outer_predictions.csv")
write_csv_safe(main_ridge_coefficients, "main_outer_ridge_coefficients.csv")

performance_by_repeat <- function(predictions) {
  predictions |>
    group_by(condition, window, model, repeat_id) |>
    group_modify(function(.x, .y) {
      metrics <- prediction_metrics(.x$truth, .x$probability, .x$predicted)
      if (.y$model == "Prevalence") {
        metrics$auroc <- 0.5
        metrics$average_precision <- metrics$prevalence
      }
      metrics
    }) |>
    ungroup()
}

summarise_performance <- function(performance) {
  metric_cols <- c(
    "auroc", "average_precision", "accuracy", "balanced_accuracy",
    "sensitivity", "specificity", "precision", "f1", "brier",
    "calibration_intercept", "calibration_slope"
  )
  
  performance |>
    pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") |>
    group_by(condition, window, model, metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      lower_95 = as.numeric(quantile(value, 0.025, na.rm = TRUE)),
      upper_95 = as.numeric(quantile(value, 0.975, na.rm = TRUE)),
      .groups = "drop"
    )
}

average_repeated_predictions <- function(predictions) {
  predictions |>
    group_by(uid, truth, condition, window, model) |>
    summarise(probability = mean(probability), .groups = "drop")
}

bootstrap_auc_ci <- function(data, reps, seed, force_chance = FALSE) {
  if (force_chance) {
    return(tibble(auroc = 0.5, lower_95 = 0.5, upper_95 = 0.5, bootstrap_reps = reps))
  }
  
  set.seed(seed)
  case_rows <- which(data$truth == 1L)
  control_rows <- which(data$truth == 0L)
  bootstrap_auc <- replicate(reps, {
    sampled <- c(
      sample(case_rows, length(case_rows), replace = TRUE),
      sample(control_rows, length(control_rows), replace = TRUE)
    )
    safe_auc(data$truth[sampled], data$probability[sampled])
  })
  
  tibble(
    auroc = safe_auc(data$truth, data$probability),
    lower_95 = as.numeric(quantile(bootstrap_auc, 0.025, na.rm = TRUE)),
    upper_95 = as.numeric(quantile(bootstrap_auc, 0.975, na.rm = TRUE)),
    bootstrap_reps = reps
  )
}

main_performance_repeat <- performance_by_repeat(main_predictions)
main_performance_summary <- summarise_performance(main_performance_repeat)

write_csv_safe(main_performance_repeat, "main_performance_by_repeat.csv")
write_csv_safe(main_performance_summary, "main_performance_summary.csv")

main_average_predictions <- average_repeated_predictions(main_predictions)
main_auc_bootstrap <- main_average_predictions |>
  group_by(condition, window, model) |>
  group_modify(function(.x, .y) {
    bootstrap_auc_ci(
      .x,
      reps = BOOTSTRAP_REPS,
      seed = MASTER_SEED + 100L * .y$window + match(.y$model, main_models),
      force_chance = .y$model == "Prevalence"
    )
  }) |>
  ungroup()

write_csv_safe(main_auc_bootstrap, "main_auroc_bootstrap_ci.csv")

# Manuscript-friendly long table including the prevalence baseline.
table_model_performance <- main_performance_summary |>
  filter(metric %in% c(
    "auroc", "average_precision", "accuracy", "balanced_accuracy",
    "sensitivity", "specificity", "f1", "brier"
  )) |>
  mutate(
    estimate = sprintf("%.3f", mean),
    repeat_interval = sprintf("%.3f to %.3f", lower_95, upper_95)
  ) |>
  select(model, metric, window, estimate, repeat_interval)

write_csv_safe(table_model_performance, "table_model_performance.csv")

# ==============================================================================
# 7. Cumulative-addition and block-removal ridge analyses, including EMA
# ==============================================================================

block_windows <- if (RUN_BLOCKS_FOR_ALL_WINDOWS) WINDOWS else PRIMARY_WINDOW

block_specs <- list(
  "Demographics only" = c("demographics"),
  "Demographics + baseline" = c("demographics", "baseline"),
  "Demographics + baseline + passive" = c("demographics", "baseline", "passive"),
  "Full: demographics + baseline + passive + EMA" = c(
    "demographics", "baseline", "passive", "ema"
  ),
  "Full minus demographics" = c("baseline", "passive", "ema"),
  "Full minus baseline" = c("demographics", "passive", "ema"),
  "Full minus passive" = c("demographics", "baseline", "ema"),
  "Full minus EMA" = c("demographics", "baseline", "passive")
)

full_condition_name <- "Full: demographics + baseline + passive + EMA"

block_predictions <- purrr::map_dfr(block_windows, function(window) {
  # Do not refit the full ridge model here. Reuse the exact held-out predictions
  # already generated for Table 3 and change only the condition label. This
  # makes the full-model AUROC in Tables 3 and 4 identical by construction.
  full_predictions_from_main <- main_predictions |>
    filter(
      .data$window == .env$window,
      model == "Ridge",
      condition == "Full model including EMA"
    ) |>
    mutate(condition = full_condition_name)
  
  reduced_and_cumulative_predictions <- purrr::imap_dfr(
    block_specs[names(block_specs) != full_condition_name],
    function(blocks, condition_name) {
      run_outer_condition(
        window = window,
        condition = condition_name,
        blocks = blocks,
        models = "Ridge",
        resamples = outer_resamples
      )$predictions
    }
  )
  
  bind_rows(reduced_and_cumulative_predictions, full_predictions_from_main)
})

write_csv_safe(block_predictions, "block_outer_predictions.csv")

# Row-level audit of the shared full-model predictions. This is stronger than
# comparing only the final AUROC: every held-out probability, threshold, tuning
# value, and fold assignment must be identical.
main_full_prediction_rows <- main_predictions |>
  filter(model == "Ridge", condition == "Full model including EMA",
         window %in% block_windows) |>
  arrange(window, repeat_id, fold, uid) |>
  select(-condition)

block_full_prediction_rows <- block_predictions |>
  filter(model == "Ridge", condition == full_condition_name) |>
  arrange(window, repeat_id, fold, uid) |>
  select(-condition)

if (!isTRUE(all.equal(
  main_full_prediction_rows,
  block_full_prediction_rows,
  tolerance = 0,
  check.attributes = TRUE
))) {
  stop(
    "The full ridge prediction rows differ between the main and block analyses."
  )
}

full_prediction_identity_audit <- tibble(
  analysis_version = ANALYSIS_VERSION,
  rows_compared = nrow(main_full_prediction_rows),
  predictions_identical = TRUE,
  maximum_absolute_probability_difference = 0
)

write_csv_safe(
  full_prediction_identity_audit,
  "Table3_Table4_full_prediction_identity_audit.csv"
)

block_performance_repeat <- performance_by_repeat(block_predictions)
block_performance_summary <- summarise_performance(block_performance_repeat)

write_csv_safe(block_performance_repeat, "block_performance_by_repeat.csv")
write_csv_safe(block_performance_summary, "block_performance_summary.csv")

paired_bootstrap_auc_difference <- function(data_a, data_b, reps, seed) {
  paired <- inner_join(
    data_a |> select(uid, truth, probability_a = probability),
    data_b |> select(uid, truth_b = truth, probability_b = probability),
    by = "uid"
  )
  
  if (any(paired$truth != paired$truth_b)) stop("Outcome mismatch in paired comparison.")
  paired <- paired |> select(-truth_b)
  
  point <- safe_auc(paired$truth, paired$probability_a) -
    safe_auc(paired$truth, paired$probability_b)
  
  set.seed(seed)
  case_rows <- which(paired$truth == 1L)
  control_rows <- which(paired$truth == 0L)
  bootstrap_delta <- replicate(reps, {
    sampled <- c(
      sample(case_rows, length(case_rows), replace = TRUE),
      sample(control_rows, length(control_rows), replace = TRUE)
    )
    safe_auc(paired$truth[sampled], paired$probability_a[sampled]) -
      safe_auc(paired$truth[sampled], paired$probability_b[sampled])
  })
  
  tibble(
    auc_difference = point,
    lower_95 = as.numeric(quantile(bootstrap_delta, 0.025, na.rm = TRUE)),
    upper_95 = as.numeric(quantile(bootstrap_delta, 0.975, na.rm = TRUE)),
    bootstrap_reps = reps
  )
}

comparison_pairs <- tribble(
  ~comparison, ~comparison_model, ~reference_model,
  "Gain from baseline surveys",
  "Demographics + baseline",
  "Demographics only",
  "Gain from passive sensing",
  "Demographics + baseline + passive",
  "Demographics + baseline",
  "Gain from EMA",
  "Full: demographics + baseline + passive + EMA",
  "Demographics + baseline + passive",
  "Change after removing demographics",
  "Full minus demographics",
  "Full: demographics + baseline + passive + EMA",
  "Change after removing baseline",
  "Full minus baseline",
  "Full: demographics + baseline + passive + EMA",
  "Change after removing passive sensing",
  "Full minus passive",
  "Full: demographics + baseline + passive + EMA",
  "Change after removing EMA",
  "Full minus EMA",
  "Full: demographics + baseline + passive + EMA"
)

block_average <- average_repeated_predictions(block_predictions)

# Absolute AUROCs for Table 4 use the same participant-level probabilities as
# the paired AUROC changes: each participant's five outer-CV probabilities are
# first averaged, and AUROC is then calculated once from those averaged OOF
# probabilities. Consequently, every reported change is the exact arithmetic
# difference between the two point estimates shown in Table 4.
block_nonfull_auc_bootstrap <- block_average |>
  filter(condition != full_condition_name) |>
  group_by(condition, window, model) |>
  group_modify(function(.x, .y) {
    bootstrap_auc_ci(
      .x,
      reps = BOOTSTRAP_REPS,
      seed = MASTER_SEED + 500000L +
        100L * .y$window +
        match(.y$condition, names(block_specs))
    )
  }) |>
  ungroup()

# Reuse the exact full-ridge AUROC and bootstrap interval already reported in
# Table 3. Thus Table 4 cannot silently substitute a mean of repeat-specific
# AUROCs or regenerate a slightly different confidence interval.
block_full_auc_bootstrap <- main_auc_bootstrap |>
  filter(model == "Ridge", window %in% block_windows) |>
  transmute(
    condition = full_condition_name,
    window,
    model,
    auroc,
    lower_95,
    upper_95,
    bootstrap_reps
  )

block_auc_bootstrap <- bind_rows(
  block_nonfull_auc_bootstrap,
  block_full_auc_bootstrap
)

write_csv_safe(block_auc_bootstrap, "block_auroc_bootstrap_ci.csv")

paired_block_auc <- purrr::map_dfr(block_windows, function(window) {
  purrr::pmap_dfr(comparison_pairs, function(comparison, comparison_model, reference_model) {
    data_a <- block_average |>
      filter(
        .data$window == .env$window,
        .data$condition == .env$comparison_model
      )
    data_b <- block_average |>
      filter(
        .data$window == .env$window,
        .data$condition == .env$reference_model
      )
    paired_bootstrap_auc_difference(
      data_a,
      data_b,
      reps = BOOTSTRAP_REPS,
      seed = MASTER_SEED + window + match(comparison, comparison_pairs$comparison)
    ) |>
      mutate(
        window = .env$window,
        comparison = .env$comparison,
        comparison_model = .env$comparison_model,
        reference_model = .env$reference_model,
        .before = 1
      )
  })
})

write_csv_safe(paired_block_auc, "paired_block_auc_differences.csv")

# Mechanical audit: every paired point change must equal the displayed AUROC
# of the comparison model minus that of its stated reference model. This guards
# against the earlier Table 3/Table 4 mismatch and against reversed removal signs.
block_change_audit <- paired_block_auc |>
  left_join(
    block_auc_bootstrap |>
      select(window, comparison_model = condition,
             comparison_model_auroc = auroc),
    by = c("window", "comparison_model")
  ) |>
  left_join(
    block_auc_bootstrap |>
      select(window, reference_model = condition,
             reference_model_auroc = auroc),
    by = c("window", "reference_model")
  ) |>
  mutate(
    arithmetic_difference = comparison_model_auroc - reference_model_auroc,
    discrepancy = auc_difference - arithmetic_difference
  )

if (
  any(!is.finite(block_change_audit$discrepancy)) ||
  max(abs(block_change_audit$discrepancy)) > 1e-12
) {
  stop(
    "A paired AUROC change does not equal the arithmetic difference between ",
    "the two displayed AUROC point estimates."
  )
}

write_csv_safe(block_change_audit, "block_change_consistency_audit.csv")

# ==============================================================================
# 8. Cross-validated ridge coefficient stability for the primary window
# ==============================================================================

predictor_block <- function(predictor) {
  case_when(
    predictor %in% demo_cols ~ "Demographics",
    predictor %in% baseline_cols ~ "Baseline",
    grepl("^ps__", predictor) ~ "Passive sensing",
    grepl("^ema__", predictor) ~ "EMA",
    TRUE ~ "Other"
  )
}

# Table 5 is based on the 25 ridge models fitted during the repeated outer CV,
# not on a separate full-sample refit. These coefficients are captured at the
# moment each Table 3 model is fitted, so no additional fitting is required and
# the coefficient summaries necessarily correspond to Table 3.
outer_ridge_coefficients <- main_ridge_coefficients |>
  filter(
    window == PRIMARY_WINDOW,
    condition == "Full model including EMA"
  ) |>
  select(-condition)

write_csv_safe(
  outer_ridge_coefficients,
  "ridge_outer_fold_coefficients_primary_window.csv"
)

total_outer_ridge_fits <- n_distinct(outer_ridge_coefficients$outer_fit_id)
expected_outer_ridge_fits <- OUTER_FOLDS * OUTER_REPEATS

if (total_outer_ridge_fits != expected_outer_ridge_fits) {
  stop(
    "Expected ", expected_outer_ridge_fits,
    " outer ridge fits but found ", total_outer_ridge_fits, "."
  )
}

ridge_coefficient_stability <- outer_ridge_coefficients |>
  group_by(predictor) |>
  summarise(
    models_available = n_distinct(outer_fit_id),
    median_coefficient = median(coefficient),
    coefficient_q1 = as.numeric(quantile(coefficient, 0.25)),
    coefficient_q3 = as.numeric(quantile(coefficient, 0.75)),
    median_absolute_coefficient = median(absolute_coefficient),
    positive_models = sum(coefficient > 0),
    negative_models = sum(coefficient < 0),
    zero_models = sum(coefficient == 0),
    .groups = "drop"
  ) |>
  mutate(
    block = predictor_block(predictor),
    median_direction = case_when(
      median_coefficient > 0 ~ "Positive",
      median_coefficient < 0 ~ "Negative",
      TRUE ~ "No consistent direction"
    ),
    sign_consistency_percent = case_when(
      median_coefficient > 0 ~ 100 * positive_models / models_available,
      median_coefficient < 0 ~ 100 * negative_models / models_available,
      TRUE ~ 100 * pmax(positive_models, negative_models) / models_available
    ),
    fit_availability_percent = 100 * models_available / total_outer_ridge_fits
  ) |>
  arrange(
    desc(median_absolute_coefficient),
    desc(sign_consistency_percent),
    predictor
  )

write_csv_safe(
  ridge_coefficient_stability,
  "ridge_coefficient_stability_primary_window.csv"
)

# ==============================================================================
# 9. Descriptive characteristics without multiplicity-heavy p-value screening
# ==============================================================================

continuous_characteristic_cols <- intersect(
  c("age", "ace", "bfi_ex", "bfi_ne", "bfi_con", "iip", "phq9", "phq15", "gad", "audit", "ders_total"),
  names(analysis_data)
)

continuous_characteristics <- purrr::map_dfr(continuous_characteristic_cols, function(variable) {
  x0 <- analysis_data[[variable]][analysis_data$pdd == 0L]
  x1 <- analysis_data[[variable]][analysis_data$pdd == 1L]
  n0 <- sum(!is.na(x0))
  n1 <- sum(!is.na(x1))
  pooled_sd <- sqrt(((n0 - 1) * var(x0, na.rm = TRUE) +
                       (n1 - 1) * var(x1, na.rm = TRUE)) / (n0 + n1 - 2))
  
  tibble(
    characteristic = variable,
    mdd_only_mean = mean(x0, na.rm = TRUE),
    mdd_only_sd = sd(x0, na.rm = TRUE),
    mdd_pdd_mean = mean(x1, na.rm = TRUE),
    mdd_pdd_sd = sd(x1, na.rm = TRUE),
    standardized_mean_difference = (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / pooled_sd
  )
})

categorical_characteristic_cols <- intersect(
  c("sex", "job", "education", "living_status", "race", "income"),
  names(analysis_data)
)

categorical_characteristics <- purrr::map_dfr(categorical_characteristic_cols, function(variable) {
  analysis_data |>
    transmute(pdd, level = as.character(.data[[variable]])) |>
    filter(!is.na(level)) |>
    count(pdd, level, name = "n") |>
    group_by(pdd) |>
    mutate(percent = 100 * n / sum(n), characteristic = variable, .before = 1) |>
    ungroup()
})

write_csv_safe(continuous_characteristics, "continuous_characteristics.csv")
write_csv_safe(categorical_characteristics, "categorical_characteristics.csv")

# ==============================================================================
# 10. Probability overlap and calibration outputs for the full ridge model
# ==============================================================================

primary_full_oof <- block_average |>
  filter(
    window == PRIMARY_WINDOW,
    condition == "Full: demographics + baseline + passive + EMA"
  ) |>
  mutate(group = factor(truth, levels = c(0, 1), labels = c("MDD-only", "MDD + PDD")))

probability_plot <- ggplot(primary_full_oof, aes(x = probability, fill = group, color = group)) +
  geom_density(alpha = 0.20, linewidth = 0.8) +
  labs(
    x = "Cross-validated predicted probability",
    y = "Density",
    fill = NULL,
    color = NULL,
    title = paste0("Predicted-probability overlap: ", PRIMARY_WINDOW, "-day full ridge model")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave(
  file.path(OUTPUT_DIR, "predicted_probability_overlap.png"),
  probability_plot,
  width = 7,
  height = 5,
  dpi = 300
)

calibration_data <- primary_full_oof |>
  mutate(bin = ntile(probability, 10L)) |>
  group_by(bin) |>
  summarise(
    mean_predicted = mean(probability),
    observed_proportion = mean(truth),
    n = n(),
    .groups = "drop"
  )

write_csv_safe(calibration_data, "calibration_deciles_primary_ridge.csv")

calibration_plot <- ggplot(calibration_data, aes(mean_predicted, observed_proportion)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(aes(size = n), color = "#1769aa") +
  geom_line(color = "#1769aa") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Mean predicted probability",
    y = "Observed proportion",
    size = "N",
    title = paste0("Calibration: ", PRIMARY_WINDOW, "-day full ridge model")
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(OUTPUT_DIR, "calibration_primary_ridge.png"),
  calibration_plot,
  width = 6,
  height = 6,
  dpi = 300
)

# ==============================================================================
# 11. Manuscript-ready Tables 2-5 (CSV and Word)
# ==============================================================================

format_p_value <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

export_manuscript_table <- function(
    data,
    title,
    filename_stem,
    note,
    landscape = FALSE
) {
  write_csv_safe(data, paste0(filename_stem, ".csv"))
  
  table_object <- flextable::flextable(data) |>
    flextable::add_header_lines(values = title) |>
    flextable::theme_booktabs() |>
    flextable::bold(part = "header") |>
    flextable::fontsize(size = 9, part = "all") |>
    flextable::padding(padding = 4, part = "all") |>
    flextable::align(align = "center", part = "all") |>
    flextable::align(j = 1, align = "left", part = "body") |>
    flextable::add_footer_lines(values = note) |>
    flextable::autofit()
  
  if ("Panel" %in% names(data)) {
    table_object <- table_object |>
      flextable::merge_v(j = "Panel") |>
      flextable::valign(j = "Panel", valign = "top")
  }
  
  output_path <- file.path(OUTPUT_DIR, paste0(filename_stem, ".docx"))
  
  if (landscape) {
    landscape_section <- officer::prop_section(
      page_size = officer::page_size(orient = "landscape"),
      page_margins = officer::page_mar(
        top = 0.5,
        bottom = 0.5,
        left = 0.5,
        right = 0.5
      )
    )
    flextable::save_as_docx(
      table_object,
      path = output_path,
      pr_section = landscape_section
    )
  } else {
    flextable::save_as_docx(table_object, path = output_path)
  }
}

# ---- Table 2: Participant characteristics ------------------------------------

table2_source <- analysis_data |>
  mutate(
    sex_display = factor(sex, levels = c(0, 1), labels = c("Male", "Female")),
    employment_display = factor(
      job,
      levels = c(0, 1),
      labels = c("Other", "Full-time job")
    ),
    education_display = factor(
      case_when(
        education %in% 1:3 ~ "High school graduate or less",
        education %in% 4:8 ~ "More than high school",
        TRUE ~ NA_character_
      ),
      levels = c("High school graduate or less", "More than high school")
    ),
    living_display = factor(
      living_status,
      levels = c(0, 1),
      labels = c("Living alone", "Living with others")
    ),
    race_display = factor(
      race,
      levels = c(0, 1),
      labels = c("White", "Non-White")
    ),
    income_display = factor(
      case_when(
        income %in% 1:6 ~ "Less than $60,000",
        income %in% 7:12 ~ "$60,000 or more",
        TRUE ~ NA_character_
      ),
      levels = c("Less than $60,000", "$60,000 or more")
    )
  )

table2_continuous_labels <- c(
  age = "Age",
  ace = "ACE",
  bfi_ex = "BFI-2 Extraversion",
  bfi_ne = "BFI-2 Negative Emotionality",
  bfi_con = "BFI-2 Conscientiousness",
  iip = "IIP-32",
  phq9 = "PHQ-9",
  phq15 = "PHQ-15",
  gad = "GAD",
  audit = "AUDIT",
  ders_total = "DERS"
)
table2_continuous_labels <- table2_continuous_labels[
  names(table2_continuous_labels) %in% names(table2_source)
]

table2_continuous <- purrr::imap_dfr(
  table2_continuous_labels,
  function(label, variable) {
    x0 <- table2_source[[variable]][table2_source$pdd == 0L]
    x1 <- table2_source[[variable]][table2_source$pdd == 1L]
    n0 <- sum(!is.na(x0))
    n1 <- sum(!is.na(x1))
    pooled_sd <- sqrt(
      ((n0 - 1) * var(x0, na.rm = TRUE) + (n1 - 1) * var(x1, na.rm = TRUE)) /
        (n0 + n1 - 2)
    )
    p_value <- tryCatch(stats::t.test(x0, x1)$p.value, error = function(e) NA_real_)
    
    tibble(
      Characteristic = label,
      `MDD-only` = sprintf("%.2f (%.2f)", mean(x0, na.rm = TRUE), sd(x0, na.rm = TRUE)),
      `MDD + PDD` = sprintf("%.2f (%.2f)", mean(x1, na.rm = TRUE), sd(x1, na.rm = TRUE)),
      `p-value` = format_p_value(p_value),
      SMD = if (is.finite(pooled_sd) && pooled_sd > 0) {
        sprintf("%.3f", (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / pooled_sd)
      } else {
        ""
      }
    )
  }
)

table2_categorical_specs <- tribble(
  ~variable, ~label,
  "sex_display", "Sex",
  "employment_display", "Employment",
  "education_display", "Education",
  "living_display", "Living status",
  "race_display", "Race",
  "income_display", "Household income"
)

make_table2_categorical_block <- function(variable, label) {
  values <- table2_source[[variable]]
  valid <- !is.na(values) & !is.na(table2_source$pdd)
  values <- droplevels(values[valid])
  outcome <- table2_source$pdd[valid]
  cross_table <- table(values, outcome)
  
  p_value <- tryCatch({
    chi <- suppressWarnings(stats::chisq.test(cross_table))
    if (any(chi$expected < 5)) stats::fisher.test(cross_table)$p.value else chi$p.value
  }, error = function(e) NA_real_)
  
  header <- tibble(
    Characteristic = label,
    `MDD-only` = "",
    `MDD + PDD` = "",
    `p-value` = format_p_value(p_value),
    SMD = ""
  )
  
  levels_rows <- purrr::map_dfr(levels(values), function(level) {
    n0 <- sum(values == level & outcome == 0L)
    n1 <- sum(values == level & outcome == 1L)
    denominator0 <- sum(outcome == 0L)
    denominator1 <- sum(outcome == 1L)
    
    tibble(
      Characteristic = paste0("    ", level),
      `MDD-only` = sprintf("%d (%.1f%%)", n0, 100 * n0 / denominator0),
      `MDD + PDD` = sprintf("%d (%.1f%%)", n1, 100 * n1 / denominator1),
      `p-value` = "",
      SMD = ""
    )
  })
  
  bind_rows(header, levels_rows)
}

table2_categorical <- purrr::pmap_dfr(
  table2_categorical_specs,
  make_table2_categorical_block
)

# Preserve the draft's ordering: demographics first, then clinical measures.
table2_data <- bind_rows(
  table2_continuous |> filter(Characteristic == "Age"),
  table2_categorical,
  table2_continuous |> filter(Characteristic != "Age")
)

n_mdd <- sum(table2_source$pdd == 0L)
n_pdd <- sum(table2_source$pdd == 1L)
names(table2_data)[names(table2_data) == "MDD-only"] <- paste0("MDD-only (N = ", n_mdd, ")")
names(table2_data)[names(table2_data) == "MDD + PDD"] <- paste0("MDD + PDD (N = ", n_pdd, ")")

export_manuscript_table(
  table2_data,
  "Table 2. Demographic and clinical characteristics of the analytic sample",
  "Table2_participant_characteristics",
  paste(
    "Note. Continuous variables are mean (SD); categorical variables are n (%).",
    "SMD is calculated as MDD + PDD minus MDD-only. P-values are descriptive;",
    "Fisher's exact test replaces the chi-squared test when expected cell counts are below 5."
  )
)

# ---- Table 3: Model performance across observation windows ------------------

model_display_labels <- c(
  Prevalence = "Prevalence",
  Ridge = "Logistic (ridge)",
  SVM = "SVM",
  XGBoost = "XGBoost"
)

# Main Table 3 is deliberately restricted to the primary metric. Its AUROCs
# are calculated from participant-level OOF probabilities averaged across the
# five outer-CV repetitions, exactly as in Table 4.
table3_data <- main_auc_bootstrap |>
  filter(model != "Prevalence") |>
  mutate(
    `Observation window` = paste(window, "days"),
    Model = unname(model_display_labels[model]),
    `AUROC (95% CI)` = sprintf(
      "%.3f (%.3f–%.3f)",
      auroc,
      lower_95,
      upper_95
    ),
    window_order = match(window, WINDOWS),
    model_order = match(model, c("Ridge", "SVM", "XGBoost"))
  ) |>
  arrange(window_order, model_order) |>
  select(`Observation window`, Model, `AUROC (95% CI)`) |>
  pivot_wider(names_from = Model, values_from = `AUROC (95% CI)`)

export_manuscript_table(
  table3_data,
  "Table 3. Cross-validated AUROC by observation window and classifier",
  "Table3_model_performance",
  paste(
    "Note. Values are AUROCs with stratified participant-level bootstrap 95%",
    "confidence intervals based on 2,000 bootstrap samples. AUROCs were calculated",
    "from participant-level out-of-fold predicted probabilities averaged across",
    "five repeated outer cross-validation runs. All models included demographic,",
    "baseline survey, passive-sensing, and EMA features. The prevalence comparator",
    "is omitted because its AUROC is 0.500 for every observation window."
  )
)

# Secondary performance measures are kept out of the main table and reported
# in Supplementary Table S1 as mean (SD) across the five outer-CV repetitions.
# AUROC is intentionally omitted because it is the sole metric in main Table 3.
secondary_metric_labels <- c(
  average_precision = "Average precision",
  f1 = "F1 score",
  accuracy = "Accuracy",
  balanced_accuracy = "Balanced accuracy",
  sensitivity = "Sensitivity",
  specificity = "Specificity",
  brier = "Brier score"
)

supplementary_s1 <- main_performance_summary |>
  filter(metric %in% c(
    "average_precision", "f1", "accuracy", "balanced_accuracy",
    "sensitivity", "specificity", "brier"
  )) |>
  mutate(
    Model = unname(model_display_labels[model]),
    Metric = unname(secondary_metric_labels[metric]),
    `Mean (SD)` = sprintf("%.3f (%.3f)", mean, sd),
    `Observation window` = paste(window, "days"),
    model_order = match(model, c("Prevalence", "Ridge", "SVM", "XGBoost")),
    metric_order = match(metric, names(secondary_metric_labels)),
    window_order = match(window, WINDOWS)
  ) |>
  arrange(model_order, metric_order, window_order) |>
  select(Model, Metric, `Observation window`, `Mean (SD)`) |>
  pivot_wider(names_from = `Observation window`, values_from = `Mean (SD)`)

export_manuscript_table(
  supplementary_s1,
  "Supplementary Table S1. Secondary cross-validated performance measures by observation window and classifier",
  "Supplementary_Table_S1_secondary_metrics",
  paste(
    "Note. Values are means (SDs) across five repeated outer cross-validation runs.",
    "For the learned classifiers, the classification threshold for F1 score, accuracy,",
    "balanced accuracy, sensitivity, and specificity was selected within each outer",
    "training set using inner cross-validation. Lower Brier scores indicate better",
    "probabilistic accuracy.",
    "AUROC is reported separately in main Table 3."
  ),
  landscape = TRUE
)

# ---- Table 4: Cumulative-addition and block-removal ridge analyses ----------

main_full_primary_auc <- main_auc_bootstrap |>
  filter(window == PRIMARY_WINDOW, model == "Ridge") |>
  select(auroc, lower_95, upper_95)

block_full_primary_auc <- block_auc_bootstrap |>
  filter(
    window == PRIMARY_WINDOW,
    model == "Ridge",
    condition == "Full: demographics + baseline + passive + EMA"
  ) |>
  select(auroc, lower_95, upper_95)

if (
  nrow(main_full_primary_auc) != 1L ||
  nrow(block_full_primary_auc) != 1L ||
  !isTRUE(all.equal(main_full_primary_auc, block_full_primary_auc, tolerance = 1e-12))
) {
  stop(
    "The full ridge AUROC or CI differs between the main and block analyses. ",
    "Check the shared resamples, model seeds, and feature definitions."
  )
}

full_auc_identity_audit <- tibble(
  analysis_version = ANALYSIS_VERSION,
  observation_window = PRIMARY_WINDOW,
  table3_auroc = main_full_primary_auc$auroc,
  table4_full_auroc = block_full_primary_auc$auroc,
  table3_lower_95 = main_full_primary_auc$lower_95,
  table4_lower_95 = block_full_primary_auc$lower_95,
  table3_upper_95 = main_full_primary_auc$upper_95,
  table4_upper_95 = block_full_primary_auc$upper_95,
  values_identical = TRUE
)

write_csv_safe(
  full_auc_identity_audit,
  "Table3_Table4_full_auroc_identity_audit.csv"
)

table4_auc <- block_auc_bootstrap |>
  filter(window == PRIMARY_WINDOW, model == "Ridge") |>
  transmute(
    model_input = condition,
    auroc,
    auc_lower_95 = lower_95,
    auc_upper_95 = upper_95
  )

table4_changes <- paired_block_auc |>
  filter(window == PRIMARY_WINDOW) |>
  select(comparison, auc_difference, lower_95, upper_95)

table4_drop_sign_check <- paired_block_auc |>
  filter(
    window == PRIMARY_WINDOW,
    grepl("^Change after removing", comparison)
  ) |>
  transmute(
    comparison,
    formula = paste0(
      "AUROC(", comparison_model, ") - AUROC(", reference_model, ")"
    ),
    signed_change = auc_difference,
    lower_95,
    upper_95
  )

write_csv_safe(table4_drop_sign_check, "Table4_drop_sign_check.csv")

table4_spec <- tribble(
  ~Panel, ~model_input, ~comparison,
  "A. Cumulative addition", "Demographics only", NA_character_,
  "A. Cumulative addition", "Demographics + baseline", "Gain from baseline surveys",
  "A. Cumulative addition", "Demographics + baseline + passive", "Gain from passive sensing",
  "A. Cumulative addition", "Full: demographics + baseline + passive + EMA", "Gain from EMA",
  "B. Block removal", "Full: demographics + baseline + passive + EMA", NA_character_,
  "B. Block removal", "Full minus demographics", "Change after removing demographics",
  "B. Block removal", "Full minus baseline", "Change after removing baseline",
  "B. Block removal", "Full minus passive", "Change after removing passive sensing",
  "B. Block removal", "Full minus EMA", "Change after removing EMA"
)

table4_data <- table4_spec |>
  left_join(table4_auc, by = "model_input") |>
  left_join(table4_changes, by = "comparison") |>
  mutate(
    `Model input` = model_input,
    `AUROC (95% CI)` = sprintf(
      "%.3f (%.3f–%.3f)",
      auroc,
      auc_lower_95,
      auc_upper_95
    ),
    Comparison = if_else(is.na(comparison), "Reference", comparison),
    `AUROC change (95% CI)` = if_else(
      is.na(auc_difference),
      "Reference",
      sprintf("%+.3f (%.3f to %.3f)", auc_difference, lower_95, upper_95)
    )
  ) |>
  select(Panel, `Model input`, `AUROC (95% CI)`, Comparison,
         `AUROC change (95% CI)`)

export_manuscript_table(
  table4_data,
  paste0(
    "Table 4. Changes in cross-validated AUROC following cumulative addition and ",
    "removal of predictor blocks in the ",
    PRIMARY_WINDOW,
    "-day ridge logistic model"
  ),
  "Table4_predictor_block_analysis",
  paste(
    "Note. Absolute AUROCs and paired changes were calculated from the same",
    "participant-level out-of-fold probabilities averaged across five repeated outer",
    "cross-validation runs. For cumulative addition, change was calculated as AUROC of",
    "the expanded model minus AUROC of the preceding model. For block removal, change",
    "was calculated as AUROC of the reduced model minus AUROC of the full model.",
    "Accordingly, positive cumulative values indicate improvement after addition,",
    "whereas negative removal values indicate poorer discrimination after removal.",
    "Confidence intervals used 2,000 stratified participant-level bootstrap samples."
  )
)

# ---- Table 5: Cross-validated ridge coefficient stability -------------------

clean_predictor_label <- function(x) {
  is_ema <- grepl("^ema__", x)
  source <- case_when(
    grepl("^ps__", x) ~ "Passive: ",
    is_ema ~ "EMA: ",
    TRUE ~ ""
  )
  label <- gsub("^(ps|ema)__", "", x)
  label <- gsub(paste0("__mean__w", PRIMARY_WINDOW, "$"), " mean", label)
  label <- gsub(paste0("__sd__w", PRIMARY_WINDOW, "$"), " SD", label)
  label <- gsub(paste0("__coverage__w", PRIMARY_WINDOW, "$"), " coverage", label)
  label <- gsub("_", " ", label)
  label[is_ema] <- gsub(
    "^PHQ[. -]*",
    "PHQ-",
    label[is_ema],
    ignore.case = TRUE
  )
  label[tolower(label) == "gad"] <- "GAD"
  label[tolower(label) == "iip"] <- "IIP"
  label[tolower(label) == "phq9"] <- "PHQ-9"
  label[tolower(label) == "phq15"] <- "PHQ-15"
  paste0(source, label)
}

invalid_ema_coefficients <- outer_ridge_coefficients |>
  filter(
    grepl("^ema__", predictor),
    !grepl(
      paste0(
        "^ema__PHQ[._-]?[1-9]__(mean|sd|coverage)__w",
        PRIMARY_WINDOW,
        "$"
      ),
      predictor,
      ignore.case = TRUE
    )
  ) |>
  distinct(predictor)

if (nrow(invalid_ema_coefficients) > 0L) {
  stop(
    "Non-PHQ predictors entered the EMA coefficient table: ",
    paste(invalid_ema_coefficients$predictor, collapse = ", ")
  )
}

table5_data <- ridge_coefficient_stability |>
  slice_head(n = 10L) |>
  mutate(
    Rank = row_number(),
    Predictor = clean_predictor_label(predictor),
    `Median beta (IQR)` = sprintf(
      "%.4f (%.4f–%.4f)",
      median_coefficient,
      coefficient_q1,
      coefficient_q3
    ),
    `Sign consistency (%)` = sprintf("%.0f", sign_consistency_percent),
    `Models (n)` = paste0(models_available, "/", total_outer_ridge_fits)
  ) |>
  select(
    Rank,
    Predictor,
    Block = block,
    `Median beta (IQR)`,
    `Sign consistency (%)`,
    `Models (n)`
  )

export_manuscript_table(
  table5_data,
  paste0(
    "Table 5. Top 10 predictors by cross-validated ridge coefficient magnitude in the ",
    PRIMARY_WINDOW,
    "-day model"
  ),
  "Table5_cross_validated_coefficient_stability",
  paste(
    "Note. Predictors were ranked by the median absolute standardized coefficient across",
    total_outer_ridge_fits,
    "outer cross-validation models. The displayed coefficient is the median signed",
    "coefficient, with the 25th and 75th percentiles in parentheses. Sign consistency is",
    "the percentage of contributing models in which the coefficient had the same direction",
    "as its median. Predictors were imputed and standardized separately within each outer",
    "training set, and the penalty was selected by inner cross-validation. These descriptive",
    "penalized coefficients should not be interpreted as independent or inferential effects."
  )
)

# Save enough metadata to reproduce and accurately describe the run.
run_manifest <- tibble(
  setting = c(
    "analysis_version", "analysis_date", "primary_window", "windows", "minimum_retained_days",
    "maximum_passive_cell_missing", "minimum_window_coverage",
    "outer_folds", "outer_repeats", "inner_folds", "bootstrap_reps",
    "daily_missing_data_strategy", "participant_level_missing_data_strategy",
    "ema_in_full_model", "ema_variables", "threshold_selection",
    "xgboost_feature_selection", "table3_table4_full_model_source",
    "table5_coefficient_summary"
  ),
  value = c(
    ANALYSIS_VERSION, as.character(Sys.Date()), as.character(PRIMARY_WINDOW),
    paste(WINDOWS, collapse = ","),
    as.character(MIN_RETAINED_DAYS), as.character(MAX_PASSIVE_CELL_MISSING),
    as.character(MIN_WINDOW_COVERAGE), as.character(OUTER_FOLDS),
    as.character(OUTER_REPEATS), as.character(INNER_FOLDS),
    as.character(BOOTSTRAP_REPS),
    "Observed-day summaries; no cross-participant daily imputation",
    "Training-fold median imputation after window summarization",
    "Yes", "Momentary PHQ-1 through PHQ-9 only",
    "Inner-CV F1 optimization", "None; tuning only",
    "One shared set of 90-day full-ridge OOF predictions",
    paste0(
      "Median signed standardized coefficient, IQR, and sign consistency across ",
      OUTER_FOLDS * OUTER_REPEATS,
      " outer ridge models"
    )
  )
)

write_csv_safe(run_manifest, "run_manifest.csv")
if (EXPORT_PARTICIPANT_LEVEL) {
  saveRDS(outer_resamples, file.path(OUTPUT_DIR, "outer_resamples.rds"))
}

if (VERBOSE) {
  message("Analysis complete. Results written to: ", OUTPUT_DIR)
}

# Additional demographic predictor audit
expected_demographics <- c(
  "age", "sex", "job", "education",
  "living_status", "race", "income"
)

demographic_audit <- tibble(
  predictor = expected_demographics,
  present_in_analysis_data = expected_demographics %in% names(analysis_data),
  included_in_demographic_block = expected_demographics %in% demo_cols,
  missing_fraction = purrr::map_dbl(
    expected_demographics,
    function(variable) {
      if (!variable %in% names(analysis_data)) return(NA_real_)
      mean(is.na(analysis_data[[variable]]))
    }
  ),
  unique_observed_values = purrr::map_int(
    expected_demographics,
    function(variable) {
      if (!variable %in% names(analysis_data)) return(NA_integer_)
      dplyr::n_distinct(analysis_data[[variable]], na.rm = TRUE)
    }
  )
) |>
  left_join(
    ridge_coefficient_stability |>
      filter(block == "Demographics") |>
      select(
        predictor,
        models_available,
        median_coefficient,
        coefficient_q1,
        coefficient_q3,
        median_absolute_coefficient,
        sign_consistency_percent
      ),
    by = "predictor"
  )

if (VERBOSE) print(demographic_audit)

readr::write_csv(
  demographic_audit,
  file.path(OUTPUT_DIR, "demographic_predictor_audit.csv")
)

# Record operational settings alongside the preserved analysis manifest.
write_csv_safe(tibble(
  setting = c("review_revision", "master_seed", "quick_mode", "nonlinear_models",
              "blocks_all_windows", "export_participant_level", "minimum_sd_observations"),
  value = as.character(c("2026-09-28", MASTER_SEED, QUICK_MODE, RUN_NONLINEAR_MODELS,
                         RUN_BLOCKS_FOR_ALL_WINDOWS, EXPORT_PARTICIPANT_LEVEL, MIN_OBS_FOR_SD))
), "review_run_settings.csv")
writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "session_info.txt"))
message("Analysis complete. Results written to: ", OUTPUT_DIR)
