# eICU-CRD covariate imputation for the HCT-BE trajectory analysis.
# Usage: Rscript eICU_CRD_CART_imputation.R cleaned.csv output_directory
# Requires mice 3.19.0 and its dependencies. Input column names are listed below.
# MAR is assumed conditional on the variables used by the imputation model.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript eICU_CRD_CART_imputation.R cleaned.csv output_directory")
}
if (!requireNamespace("mice", quietly = TRUE)) {
  stop("Install the mice package before running this script.")
}
if (as.character(utils::packageVersion("mice")) != "3.19.0") {
  warning("The study used mice 3.19.0; results may vary with another version.")
}

input_file <- args[1]
output_dir <- args[2]
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("The output directory must be empty: ", output_dir)
}

targets <- c(
  "admit_weight", "admissionheight", "max_respiratoryrate", "max_heartrate",
  "min_map", "min_spo2", "max_temperature", "bun_max", "creatinine_min",
  "inr_max", "lactate_max", "albumin_min", "alt_max_24h", "platelet_min",
  "wbc_max", "intake_day1_ml", "intake_day2_ml", "ast_max_24h"
)
clinical_predictors <- c(
  "gender", "age", "sofa_score", "unittype", "pmh_hypertension",
  "pmh_diabetes", "pmh_cad", "pmh_copd", "pmh_ckd",
  "sodium_bicarbonate_48h", "rbc_transfusion_48h",
  "vasoactive_drug_48h", "imv_48h", "rrt_48h"
)
event_columns <- c("death_2_30d", "death_2_60d")
time_columns <- c("survival_time_2_30d_days", "survival_time_2_60d_days")
hazard_columns <- c("nelson_aalen_30", "nelson_aalen_60")
group_column <- "trajectory_group_raw"
factor_columns <- c(
  "gender", "unittype", "pmh_hypertension", "pmh_diabetes", "pmh_cad",
  "pmh_copd", "pmh_ckd", "sodium_bicarbonate_48h",
  "rbc_transfusion_48h", "vasoactive_drug_48h", "imv_48h", "rrt_48h",
  event_columns, group_column
)
binary_columns <- c(
  "pmh_hypertension", "pmh_diabetes", "pmh_cad", "pmh_copd", "pmh_ckd",
  "sodium_bicarbonate_48h", "rbc_transfusion_48h", "vasoactive_drug_48h",
  "imv_48h", "rrt_48h", event_columns
)

read_input <- function(path) {
  for (encoding in c("UTF-8-BOM", "UTF-8", "GB18030")) {
    x <- tryCatch(
      utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                      strip.white = TRUE, na.strings = c("", "NA", "NaN", "NULL"),
                      fileEncoding = encoding),
      error = function(e) NULL
    )
    if (is.data.frame(x) && "patientunitstayid" %in% names(x)) return(x)
  }
  stop("Could not read a CSV containing patientunitstayid: ", path)
}

data <- read_input(input_file)
required <- unique(c("patientunitstayid", targets, clinical_predictors,
                     event_columns, time_columns, group_column))
absent <- setdiff(required, names(data))
if (length(absent)) stop("Missing input columns: ", paste(absent, collapse = ", "))
if (anyDuplicated(names(data)) || anyNA(data$patientunitstayid) ||
    anyDuplicated(data$patientunitstayid)) {
  stop("Column names and patientunitstayid must be unique; the ID cannot be missing.")
}
if (any(hazard_columns %in% names(data))) {
  stop("Use the cleaned input before Nelson-Aalen predictors were generated.")
}
fixed <- unique(c(clinical_predictors, event_columns, time_columns, group_column))
incomplete <- fixed[vapply(data[fixed], anyNA, logical(1))]
if (length(incomplete)) {
  stop("Predictor-only columns contain missing values: ",
       paste(incomplete, collapse = ", "))
}
numeric_columns <- unique(c(targets, time_columns,
                            setdiff(clinical_predictors, factor_columns)))
invalid_numeric <- numeric_columns[!vapply(data[numeric_columns], is.numeric, logical(1))]
if (length(invalid_numeric)) {
  stop("Expected numeric columns: ", paste(invalid_numeric, collapse = ", "))
}
if (any(vapply(data[binary_columns], function(x) !all(x %in% 0:1), logical(1)))) {
  stop("Binary predictors and event indicators must use 0/1 coding.")
}
if (!all(data$gender %in% c(1, 2))) stop("gender must use 1/2 coding.")
if (!all(data[[group_column]] %in% 1:4)) {
  stop("trajectory_group_raw must use values 1, 2, 3, and 4.")
}
if (any(data$admit_weight <= 0, na.rm = TRUE) ||
    any(data$admissionheight <= 0, na.rm = TRUE)) {
  stop("Observed weight and height must be positive.")
}
if (any(data[c(time_columns, "intake_day1_ml", "intake_day2_ml")] < 0,
        na.rm = TRUE)) {
  stop("Follow-up times and fluid intake cannot be negative.")
}

# Follow-up time enters only through its Nelson-Aalen estimate.
for (j in seq_along(event_columns)) {
  survival_data <- data.frame(time = data[[time_columns[j]]],
                              event = data[[event_columns[j]]])
  data[[hazard_columns[j]]] <- mice::nelsonaalen(survival_data, time, event)
}
if (anyNA(data[hazard_columns]) ||
    any(!is.finite(as.matrix(data[hazard_columns])))) {
  stop("Nelson-Aalen estimates must be finite and complete.")
}

# HCT/BE windows, posterior probabilities, identifiers and follow-up times
# are not in the imputation model. BMI and fluid intake per kg are recalculated.
model_columns <- unique(c(targets, clinical_predictors, group_column,
                          event_columns, hazard_columns))
model_data <- data[model_columns]
for (name in factor_columns) model_data[[name]] <- factor(model_data[[name]])
method <- stats::setNames(rep("", ncol(model_data)), names(model_data))
method[targets] <- "cart"
predictors <- matrix(0L, ncol(model_data), ncol(model_data),
                     dimnames = list(names(model_data), names(model_data)))
for (name in targets) predictors[name, setdiff(names(model_data), name)] <- 1L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
imp <- mice::mice(
  model_data, m = 20L, maxit = 10L, method = method,
  predictorMatrix = predictors, visitSequence = targets,
  seed = 20260906L, minbucket = 5L, cp = 1e-4, maxdepth = 30L,
  printFlag = TRUE
)
saveRDS(imp, file.path(output_dir, "imputation.rds"))

completed_dir <- file.path(output_dir, "completed")
dir.create(completed_dir)
for (i in seq_len(20L)) {
  completed <- data
  completed[targets] <- mice::complete(imp, i)[targets]
  completed$bmi <- completed$admit_weight / (completed$admissionheight / 100)^2
  completed$fluid_intake_0_48h_per_kg <-
    (completed$intake_day1_ml + completed$intake_day2_ml) /
    completed$admit_weight
  if (anyNA(completed[c(targets, "bmi", "fluid_intake_0_48h_per_kg")]) ||
      any(!is.finite(as.matrix(completed[c(targets, "bmi",
                                           "fluid_intake_0_48h_per_kg")])))) {
    stop("Incomplete or nonfinite values in completed dataset ", i)
  }
  utils::write.csv(completed, file.path(completed_dir, sprintf("imputed_%02d.csv", i)),
                   row.names = FALSE, na = "")
}

utils::write.csv(mice::convergence(imp),
                 file.path(output_dir, "convergence.csv"), row.names = FALSE)
grDevices::pdf(file.path(output_dir, "chain_traces.pdf"), width = 12, height = 9)
print(plot(imp, targets, layout = c(3, 3)))
grDevices::dev.off()
if (is.null(imp$loggedEvents)) {
  writeLines("No logged events.", file.path(output_dir, "logged_events.txt"))
} else {
  utils::write.csv(imp$loggedEvents, file.path(output_dir, "logged_events.csv"),
                   row.names = FALSE)
}
writeLines(capture.output(utils::sessionInfo()),
           file.path(output_dir, "session_info.txt"))
message("Saved 20 completed datasets and diagnostics to: ", output_dir)
