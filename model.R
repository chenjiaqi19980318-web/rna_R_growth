# ============================================================
# Seed1111 Radiomics Modeling and Final Validation Pipeline
# ============================================================
# Endpoint:
#   Label = 0: slow growth, VDT > 400 days
#   Label = 1: fast growth, VDT <= 400 days
#
# Required input files in the working directory:
#   1. Clinical_data.xlsx
#   2. intra_data.xlsx
#   3. peri_data.xlsx
#
# Cohort definition in the original files:
#   Cohort = 2: prospective temporal validation cohort, retained unchanged
#   Cohort != 2: randomly split into training/test cohorts using seed 1111
#
# Main workflow:
#   1. Seed 1111 stratified train/test split; original Cohort=2 retained as validation
#   2. Clinical, intratumor, and peritumor features processed separately
#   3. Missing values imputed using training medians and applied to train/test/validation
#   4. Pearson pre-filtering in the training cohort only
#   5. High-correlation removal in the training cohort only
#      If |r| >= 0.80, remove the feature with lower |correlation with Label|
#   6. VIF removal in the training cohort only
#      Iteratively remove features with VIF >= 20
#   7. Z-score parameters estimated in the training cohort only
#      Training mean and SD applied to train/test/validation
#   8. LASSO performed only in standardized training data
#   9. Clinical Score: multivariable logistic regression using clinical lambda.1se variables
#  10. Intra-Radscore and Peri-Radscore: LASSO linear predictors using standardized features
#  11. Seven score-based models are built and tested
#  12. Bootstrap internal validation is performed in the training cohort
#  13. Final Clinical_Intra_Peri model is evaluated in train/test/validation cohorts
#  14. ROC, calibration, DCA, clinical impact curve, and nomogram are generated
#
# Recommended citation wording:
#   The validation cohort should be described as a prospective temporal validation cohort,
#   not an external validation cohort.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "readxl", "openxlsx", "dplyr", "tidyr", "stringr", "ggplot2",
  "glmnet", "pROC", "pheatmap", "RColorBrewer", "tibble",
  "scales", "MASS", "rms"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

while (!is.null(dev.list())) {
  dev.off()
}


# ============================================================
# 1. Parameters
# ============================================================

clinical_file <- "Clinical_data.xlsx"
intra_file <- "intra_data.xlsx"
peri_file <- "peri_data.xlsx"

clinical_sheet <- 1
intra_sheet <- 1
peri_sheet <- 1

id_col <- "PatientID"
label_col <- "Label"
cohort_col <- "Cohort"

seed_split <- 1111
train_ratio <- 0.70

train_value <- 0
test_value <- 1
validation_value <- 2

pearson_r_cutoff <- 0.10
high_corr_cutoff <- 0.80
vif_cutoff <- 20
max_features_for_vif <- 100

lasso_nfolds <- 10
bootstrap_B <- 1000
calibration_bootstrap_B <- 300

final_model <- "Clinical_Intra_Peri"
final_predictors <- c("Clinical_Score", "Intra_Radscore", "Peri_Radscore")
final_pred_col <- paste0(final_model, "_Pred")

output_dir <- "Seed1111_Radiomics_Clinical_Intra_Peri_Modeling_Validation"

preprocess_dir <- file.path(output_dir, "01_Preprocessed")
lasso_dir <- file.path(output_dir, "02_LASSO")
score_dir <- file.path(output_dir, "03_Seven_Models")
bootstrap_dir <- file.path(output_dir, "04_Bootstrap_Internal_Validation")
final_validation_dir <- file.path(output_dir, "05_Final_Model_Clinical_Intra_Peri_Validation")
calibration_plot_dir <- file.path(output_dir, "06_Final_Model_Calibration_Plot")
validation_calibration_dir <- file.path(output_dir, "07_Validation_Calibration_Plot_Final_Format")

for (dd in c(
  output_dir, preprocess_dir, lasso_dir, score_dir, bootstrap_dir,
  final_validation_dir, calibration_plot_dir, validation_calibration_dir
)) {
  dir.create(dd, showWarnings = FALSE, recursive = TRUE)
}

set.seed(seed_split)


# ============================================================
# 2. General helper functions
# ============================================================

safe_var <- function(x) {
  paste0("`", gsub("`", "``", x), "`")
}

safe_formula <- function(response, predictors) {
  predictors <- unique(predictors)
  predictors <- predictors[predictors != ""]
  if (length(predictors) == 0) return(NULL)
  as.formula(paste(safe_var(response), "~", paste(safe_var(predictors), collapse = " + ")))
}

clean_id <- function(x) {
  trimws(as.character(x))
}

clip_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

logit_prob <- function(p) {
  qlogis(clip_prob(p))
}

safe_bind_rows <- function(x) {
  x <- x[!sapply(x, is.null)]
  if (length(x) == 0) return(data.frame())
  dplyr::bind_rows(x)
}

add_sheet <- function(wb, sheet_name, df) {
  sheet_name <- substr(sheet_name, 1, 31)
  if (sheet_name %in% names(wb)) {
    sheet_name <- substr(make.unique(c(names(wb), sheet_name), sep = "_")[length(names(wb)) + 1], 1, 31)
  }
  openxlsx::addWorksheet(wb, sheet_name)
  if (is.null(df) || nrow(df) == 0) {
    openxlsx::writeData(wb, sheet_name, data.frame(Message = "No result"))
  } else {
    openxlsx::writeData(wb, sheet_name, df)
  }
}

save_pdf_safe <- function(plot_obj, filename, width = 7.2, height = 6.2) {
  grDevices::pdf(filename, width = width, height = height, onefile = TRUE, useDingbats = FALSE)
  print(plot_obj)
  grDevices::dev.off()
  if (!file.exists(filename)) stop("PDF was not generated: ", filename)
  invisible(TRUE)
}

save_base_pdf_png <- function(draw_fun, pdf_file, png_file, pdf_width = 7, pdf_height = 6.5,
                              png_width = 2100, png_height = 1950, png_res = 300) {
  grDevices::pdf(pdf_file, width = pdf_width, height = pdf_height, useDingbats = FALSE)
  draw_fun()
  grDevices::dev.off()
  
  grDevices::png(png_file, width = png_width, height = png_height, res = png_res)
  draw_fun()
  grDevices::dev.off()
}

read_input_data <- function(file, sheet = 1) {
  if (!file.exists(file)) stop("Cannot find file: ", file)
  df <- readxl::read_excel(file, sheet = sheet, guess_max = 100000)
  df <- as.data.frame(df)
  colnames(df) <- trimws(colnames(df))
  
  required_cols <- c(id_col, label_col, cohort_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(file, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df[[id_col]] <- clean_id(df[[id_col]])
  df[[label_col]] <- as.numeric(df[[label_col]])
  df[[cohort_col]] <- suppressWarnings(as.numeric(df[[cohort_col]]))
  
  if (!all(df[[label_col]] %in% c(0, 1))) {
    stop(file, " Label must contain only 0 and 1.")
  }
  if (anyDuplicated(df[[id_col]]) > 0) {
    stop(file, " contains duplicated PatientID values.")
  }
  
  df
}

check_same_patients <- function(a, b, name_a, name_b) {
  aa <- a[, c(id_col, label_col), drop = FALSE] %>% dplyr::arrange(.data[[id_col]])
  bb <- b[, c(id_col, label_col), drop = FALSE] %>% dplyr::arrange(.data[[id_col]])
  if (!identical(aa[[id_col]], bb[[id_col]])) stop(name_a, " and ", name_b, " PatientID values are inconsistent.")
  if (!identical(aa[[label_col]], bb[[label_col]])) stop(name_a, " and ", name_b, " Label values are inconsistent.")
}

make_seed1111_split <- function(df) {
  set.seed(seed_split)
  split_base <- df[, c(id_col, label_col, cohort_col), drop = FALSE]
  
  if (!any(split_base[[cohort_col]] == validation_value, na.rm = TRUE)) {
    stop("No Cohort=2 samples detected. Cohort=2 should be the validation cohort.")
  }
  
  valid_df <- split_base %>% dplyr::filter(.data[[cohort_col]] == validation_value)
  non_valid_df <- split_base %>% dplyr::filter(.data[[cohort_col]] != validation_value | is.na(.data[[cohort_col]]))
  
  train_ids <- c()
  test_ids <- c()
  
  for (lab in c(0, 1)) {
    ids_lab <- non_valid_df %>%
      dplyr::filter(.data[[label_col]] == lab) %>%
      dplyr::pull(.data[[id_col]])
    n_train <- floor(length(ids_lab) * train_ratio)
    train_lab <- sample(ids_lab, size = n_train, replace = FALSE)
    test_lab <- setdiff(ids_lab, train_lab)
    train_ids <- c(train_ids, train_lab)
    test_ids <- c(test_ids, test_lab)
  }
  
  split_df <- split_base %>%
    dplyr::mutate(
      Cohort_new = dplyr::case_when(
        .data[[id_col]] %in% train_ids ~ train_value,
        .data[[id_col]] %in% test_ids ~ test_value,
        .data[[id_col]] %in% valid_df[[id_col]] ~ validation_value,
        TRUE ~ NA_real_
      )
    )
  
  if (any(is.na(split_df$Cohort_new))) stop("Some samples were not assigned to a cohort.")
  
  split_df %>%
    dplyr::select(all_of(id_col), all_of(label_col), Cohort = Cohort_new)
}

apply_split <- function(df, split_df) {
  df <- df %>%
    dplyr::select(-dplyr::any_of(cohort_col)) %>%
    dplyr::left_join(split_df[, c(id_col, cohort_col)], by = id_col)
  df[[cohort_col]] <- as.numeric(df[[cohort_col]])
  df
}

convert_clinical_to_numeric <- function(df) {
  exclude_cols <- c(id_col, label_col, cohort_col)
  
  for (v in setdiff(colnames(df), exclude_cols)) {
    if (is.character(df[[v]]) || is.factor(df[[v]])) {
      x_raw <- as.character(df[[v]])
      x_raw[is.na(x_raw) | trimws(x_raw) == ""] <- "Missing"
      x <- as.factor(x_raw)
      
      if (length(levels(x)) <= 1) {
        df[[v]] <- NULL
      } else if (length(levels(x)) == 2) {
        df[[v]] <- as.numeric(x) - 1
      } else {
        mm <- model.matrix(~ x - 1)
        mm <- as.data.frame(mm)
        colnames(mm) <- paste0(v, "_", make.names(levels(x)))
        df[[v]] <- NULL
        df <- cbind(df, mm)
      }
    } else {
      df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
    }
  }
  
  df
}


# ============================================================
# 3. Train-only preprocessing functions
# ============================================================

impute_train_median_apply_all <- function(train_df, all_df, features) {
  out <- all_df
  if (length(features) == 0) {
    return(list(data = out, parameters = data.frame(Feature = character(0), Train_Median = numeric(0))))
  }
  
  para <- list()
  for (v in features) {
    med <- median(as.numeric(train_df[[v]]), na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    out[[v]][is.na(out[[v]])] <- med
    para[[v]] <- data.frame(Feature = v, Train_Median = med, stringsAsFactors = FALSE)
  }
  
  list(data = out, parameters = dplyr::bind_rows(para))
}

label_correlation_table <- function(train_df, features) {
  y <- as.numeric(train_df[[label_col]])
  r <- sapply(features, function(v) {
    x <- as.numeric(train_df[[v]])
    if (length(unique(x[!is.na(x)])) < 2) return(NA_real_)
    suppressWarnings(cor(x, y, method = "pearson", use = "complete.obs"))
  })
  r <- as.numeric(r)
  names(r) <- features
  r[!is.finite(r)] <- NA_real_
  data.frame(Feature = features, Pearson_r_with_Label = as.numeric(r[features]), Abs_r = abs(as.numeric(r[features])), stringsAsFactors = FALSE)
}

pearson_prefilter <- function(train_df, features, cutoff = 0.10) {
  cor_tbl <- label_correlation_table(train_df, features)
  keep <- cor_tbl %>% dplyr::filter(!is.na(Abs_r), Abs_r >= cutoff) %>% dplyr::pull(Feature)
  removed <- cor_tbl %>% dplyr::filter(is.na(Abs_r) | Abs_r < cutoff) %>% dplyr::mutate(Reason = paste0("|Pearson r with Label| < ", cutoff))
  list(keep = keep, removed = removed, table = cor_tbl)
}

high_corr_filter <- function(train_df, features, label_cor_table, cutoff = 0.80) {
  if (length(features) <= 1) {
    return(list(
      keep = features,
      removed = data.frame(Feature = character(0), Reason = character(0)),
      pairs = data.frame(Feature1 = character(0), Feature2 = character(0), Correlation = numeric(0), Removed = character(0), Kept = character(0), Reason = character(0))
    ))
  }
  
  x <- train_df[, features, drop = FALSE]
  cat("  calculating correlation matrix for ", length(features), " features...\n")
  flush.console()
  
  cm <- suppressWarnings(cor(x, use = "pairwise.complete.obs", method = "pearson"))
  cm[is.na(cm)] <- 0
  diag(cm) <- 0
  
  label_abs <- label_cor_table$Abs_r
  names(label_abs) <- label_cor_table$Feature
  label_abs[is.na(label_abs)] <- 0
  
  upper_idx <- which(abs(cm) >= cutoff & upper.tri(cm), arr.ind = TRUE)
  if (nrow(upper_idx) == 0) {
    return(list(
      keep = features,
      removed = data.frame(Feature = character(0), Reason = character(0)),
      pairs = data.frame(Feature1 = character(0), Feature2 = character(0), Correlation = numeric(0), Removed = character(0), Kept = character(0), Reason = character(0))
    ))
  }
  
  pair_df <- data.frame(
    Feature1 = rownames(cm)[upper_idx[, 1]],
    Feature2 = colnames(cm)[upper_idx[, 2]],
    Correlation = cm[upper_idx],
    stringsAsFactors = FALSE
  )
  
  pair_df$LabelAbsCor_Feature1 <- label_abs[pair_df$Feature1]
  pair_df$LabelAbsCor_Feature2 <- label_abs[pair_df$Feature2]
  pair_df$Removed <- ifelse(pair_df$LabelAbsCor_Feature1 >= pair_df$LabelAbsCor_Feature2, pair_df$Feature2, pair_df$Feature1)
  pair_df$Kept <- ifelse(pair_df$Removed == pair_df$Feature1, pair_df$Feature2, pair_df$Feature1)
  pair_df$Reason <- paste0("|feature correlation| >= ", cutoff)
  
  removed <- unique(pair_df$Removed)
  keep <- setdiff(features, removed)
  removed_df <- data.frame(Feature = removed, Reason = paste0("High pairwise correlation |r| >= ", cutoff), stringsAsFactors = FALSE)
  
  list(keep = keep, removed = removed_df, pairs = pair_df)
}

vif_filter <- function(train_df, features, cutoff = 20, max_features_for_vif = 100) {
  features_now <- features
  
  if (length(features_now) == 0) {
    return(list(keep = character(0), removed = data.frame(Iteration = integer(0), Feature = character(0), VIF = numeric(0), Reason = character(0)), final_vif = data.frame(Feature = character(0), VIF = numeric(0))))
  }
  
  if (length(features_now) == 1) {
    return(list(keep = features_now, removed = data.frame(Iteration = integer(0), Feature = character(0), VIF = numeric(0), Reason = character(0)), final_vif = data.frame(Feature = features_now, VIF = 1)))
  }
  
  if (length(features_now) > max_features_for_vif) {
    warning("VIF skipped: ", length(features_now), " features remain after high-correlation filtering, more than max_features_for_vif = ", max_features_for_vif, ".")
    return(list(
      keep = features_now,
      removed = data.frame(Iteration = integer(0), Feature = character(0), VIF = numeric(0), Reason = character(0)),
      final_vif = data.frame(Feature = features_now, VIF = NA_real_, Note = paste0("VIF skipped because number of features > ", max_features_for_vif))
    ))
  }
  
  removed_list <- list()
  iter <- 0
  
  repeat {
    iter <- iter + 1
    if (length(features_now) <= 1) break
    
    x <- as.matrix(train_df[, features_now, drop = FALSE])
    x <- scale(x)
    cor_mat <- suppressWarnings(cor(x, use = "pairwise.complete.obs"))
    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- 1
    
    inv_cor <- tryCatch(solve(cor_mat), error = function(e) NULL)
    if (is.null(inv_cor)) inv_cor <- tryCatch(MASS::ginv(cor_mat), error = function(e) NULL)
    if (is.null(inv_cor)) {
      warning("VIF calculation failed. VIF step stopped.")
      break
    }
    
    vv <- diag(inv_cor)
    names(vv) <- features_now
    vv[!is.finite(vv)] <- Inf
    max_vif <- max(vv, na.rm = TRUE)
    
    if (!is.finite(max_vif) && !is.infinite(max_vif)) break
    if (max_vif < cutoff) break
    
    remove_f <- names(which.max(vv))
    removed_list[[length(removed_list) + 1]] <- data.frame(Iteration = iter, Feature = remove_f, VIF = max_vif, Reason = paste0("VIF >= ", cutoff), stringsAsFactors = FALSE)
    features_now <- setdiff(features_now, remove_f)
  }
  
  final_vif <- if (length(features_now) > 0) {
    x_final <- as.matrix(train_df[, features_now, drop = FALSE])
    x_final <- scale(x_final)
    cor_final <- suppressWarnings(cor(x_final, use = "pairwise.complete.obs"))
    cor_final[is.na(cor_final)] <- 0
    diag(cor_final) <- 1
    inv_final <- tryCatch(solve(cor_final), error = function(e) NULL)
    if (is.null(inv_final)) inv_final <- tryCatch(MASS::ginv(cor_final), error = function(e) NULL)
    if (!is.null(inv_final)) {
      data.frame(Feature = features_now, VIF = as.numeric(diag(inv_final)), stringsAsFactors = FALSE)
    } else {
      data.frame(Feature = features_now, VIF = NA_real_, stringsAsFactors = FALSE)
    }
  } else {
    data.frame(Feature = character(0), VIF = numeric(0))
  }
  
  list(keep = features_now, removed = safe_bind_rows(removed_list), final_vif = final_vif)
}

zscore_train_apply_all <- function(train_df, all_df, features) {
  out <- all_df
  if (length(features) == 0) {
    return(list(data = out, parameters = data.frame(Feature = character(0), Train_Mean = numeric(0), Train_SD = numeric(0))))
  }
  
  para <- list()
  for (v in features) {
    mu <- mean(as.numeric(train_df[[v]]), na.rm = TRUE)
    sdv <- sd(as.numeric(train_df[[v]]), na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    if (!is.finite(sdv) || sdv == 0) sdv <- 1
    out[[v]] <- (as.numeric(out[[v]]) - mu) / sdv
    para[[v]] <- data.frame(Feature = v, Train_Mean = mu, Train_SD = sdv, stringsAsFactors = FALSE)
  }
  
  list(data = out, parameters = dplyr::bind_rows(para))
}


# ============================================================
# 4. LASSO and preprocessing pipeline
# ============================================================

make_stratified_foldid <- function(y, nfolds = 10) {
  y <- as.numeric(y)
  foldid <- rep(NA_integer_, length(y))
  min_class <- min(table(y))
  k <- min(nfolds, min_class)
  if (k < 3) k <- 3
  
  for (lab in unique(y)) {
    idx <- which(y == lab)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  
  foldid
}

run_lasso <- function(train_df, features, output_prefix) {
  if (length(features) == 0) stop(output_prefix, " has no variables available for LASSO.")
  
  x <- as.matrix(train_df[, features, drop = FALSE])
  y <- as.numeric(train_df[[label_col]])
  foldid <- make_stratified_foldid(y, lasso_nfolds)
  
  cvfit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 1,
    foldid = foldid,
    standardize = FALSE,
    type.measure = "deviance"
  )
  
  coef_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
  coef_min <- as.matrix(coef(cvfit, s = "lambda.min"))
  
  coef_1se_df <- data.frame(Feature = rownames(coef_1se), Coefficient = as.numeric(coef_1se[, 1]), stringsAsFactors = FALSE)
  coef_min_df <- data.frame(Feature = rownames(coef_min), Coefficient = as.numeric(coef_min[, 1]), stringsAsFactors = FALSE)
  
  selected_1se <- coef_1se_df %>% dplyr::filter(Feature != "(Intercept)", Coefficient != 0)
  selected_min <- coef_min_df %>% dplyr::filter(Feature != "(Intercept)", Coefficient != 0)
  
  save_base_pdf_png(
    draw_fun = function() {
      plot(cvfit)
      abline(v = log(cvfit$lambda.1se), lty = 2, col = "red")
      abline(v = log(cvfit$lambda.min), lty = 2, col = "blue")
    },
    pdf_file = file.path(lasso_dir, paste0(output_prefix, "_LASSO_CV_curve_deviance.pdf")),
    png_file = file.path(lasso_dir, paste0(output_prefix, "_LASSO_CV_curve_deviance.png")),
    pdf_width = 7,
    pdf_height = 6,
    png_width = 1800,
    png_height = 1500,
    png_res = 220
  )
  
  save_base_pdf_png(
    draw_fun = function() {
      plot(cvfit$glmnet.fit, xvar = "lambda", label = FALSE)
      abline(v = log(cvfit$lambda.1se), lty = 2, col = "red")
      abline(v = log(cvfit$lambda.min), lty = 2, col = "blue")
    },
    pdf_file = file.path(lasso_dir, paste0(output_prefix, "_LASSO_coefficient_path.pdf")),
    png_file = file.path(lasso_dir, paste0(output_prefix, "_LASSO_coefficient_path.png")),
    pdf_width = 7,
    pdf_height = 6,
    png_width = 1800,
    png_height = 1500,
    png_res = 220
  )
  
  if (nrow(selected_1se) == 0 && nrow(selected_min) > 0) {
    warning(output_prefix, " lambda.1se selected no non-zero variables; lambda.min variables are used instead.")
    selected_1se <- selected_min
  }
  
  if (nrow(selected_1se) == 0) stop(output_prefix, " selected no variables at lambda.1se or lambda.min.")
  
  list(
    cvfit = cvfit,
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    coef_min = coef_min_df,
    coef_1se = coef_1se_df,
    selected_min = selected_min,
    selected_1se = selected_1se
  )
}

preprocess_and_select <- function(df, dataset_name, is_clinical = FALSE) {
  cat("\n================ ", dataset_name, " preprocessing ================\n")
  flush.console()
  
  df0 <- df
  
  if (is_clinical) {
    df0 <- convert_clinical_to_numeric(df0)
  } else {
    for (v in setdiff(colnames(df0), c(id_col, label_col, cohort_col))) {
      df0[[v]] <- suppressWarnings(as.numeric(df0[[v]]))
    }
  }
  
  features0 <- setdiff(colnames(df0), c(id_col, label_col, cohort_col))
  features0 <- features0[sapply(df0[, features0, drop = FALSE], function(x) !all(is.na(x)))]
  train_df0 <- df0 %>% dplyr::filter(.data[[cohort_col]] == train_value)
  
  cat(dataset_name, ": initial features = ", length(features0), "\n")
  impute_obj <- impute_train_median_apply_all(train_df = train_df0, all_df = df0, features = features0)
  df1 <- impute_obj$data
  train_df1 <- df1 %>% dplyr::filter(.data[[cohort_col]] == train_value)
  
  cat(dataset_name, ": Pearson pre-filtering\n")
  pearson_obj <- pearson_prefilter(train_df = train_df1, features = features0, cutoff = pearson_r_cutoff)
  features1 <- pearson_obj$keep
  if (length(features1) == 0) stop(dataset_name, " has no variables after Pearson pre-filtering.")
  
  cat(dataset_name, ": retained after Pearson = ", length(features1), "\n")
  label_cor_after_pearson <- pearson_obj$table %>% dplyr::filter(Feature %in% features1)
  highcorr_obj <- high_corr_filter(train_df = train_df1, features = features1, label_cor_table = label_cor_after_pearson, cutoff = high_corr_cutoff)
  features2 <- highcorr_obj$keep
  if (length(features2) == 0) stop(dataset_name, " has no variables after high-correlation filtering.")
  
  cat(dataset_name, ": retained after high correlation = ", length(features2), "\n")
  vif_obj <- vif_filter(train_df = train_df1, features = features2, cutoff = vif_cutoff, max_features_for_vif = max_features_for_vif)
  features3 <- vif_obj$keep
  if (length(features3) == 0) stop(dataset_name, " has no variables after VIF filtering.")
  
  cat(dataset_name, ": retained after VIF = ", length(features3), "\n")
  z_obj <- zscore_train_apply_all(train_df = train_df1, all_df = df1, features = features3)
  df2 <- z_obj$data
  train_df2 <- df2 %>% dplyr::filter(.data[[cohort_col]] == train_value)
  
  cat(dataset_name, ": LASSO\n")
  lasso_obj <- run_lasso(train_df = train_df2, features = features3, output_prefix = dataset_name)
  cat(dataset_name, " lambda.1se selected variables: ", nrow(lasso_obj$selected_1se), "\n")
  
  list(
    data = df2,
    train_data = train_df2,
    initial_features = features0,
    imputation = impute_obj$parameters,
    pearson_table = pearson_obj$table,
    pearson_removed = pearson_obj$removed,
    high_corr_removed = highcorr_obj$removed,
    high_corr_pairs = highcorr_obj$pairs,
    vif_removed = vif_obj$removed,
    final_vif = vif_obj$final_vif,
    zscore_parameters = z_obj$parameters,
    retained_after_pearson = features1,
    retained_after_highcorr = features2,
    retained_before_lasso = features3,
    lasso = lasso_obj
  )
}


# ============================================================
# 5. Score, evaluation, bootstrap, and DeLong functions
# ============================================================

fit_logistic_score <- function(train_df, full_df, selected_features, score_name) {
  selected_features <- intersect(selected_features, colnames(full_df))
  if (length(selected_features) == 0) stop(score_name, " has no available variables.")
  
  fit <- glm(safe_formula(label_col, selected_features), data = train_df, family = binomial())
  full_df[[score_name]] <- as.numeric(predict(fit, newdata = full_df, type = "link"))
  list(data = full_df, fit = fit)
}

make_radscore_from_lasso <- function(full_df, lasso_selected, lasso_coef_all, score_name) {
  selected <- intersect(lasso_selected$Feature, colnames(full_df))
  intercept <- lasso_coef_all %>% dplyr::filter(Feature == "(Intercept)") %>% dplyr::pull(Coefficient)
  if (length(intercept) == 0) intercept <- 0
  
  coef_df <- lasso_selected %>% dplyr::filter(Feature %in% selected)
  score <- rep(intercept, nrow(full_df))
  
  if (nrow(coef_df) > 0) {
    for (i in seq_len(nrow(coef_df))) {
      score <- score + full_df[[coef_df$Feature[i]]] * coef_df$Coefficient[i]
    }
  }
  
  full_df[[score_name]] <- as.numeric(score)
  full_df
}

get_auc_ci <- function(y, p) {
  idx <- !is.na(y) & !is.na(p)
  y <- y[idx]
  p <- p[idx]
  
  if (length(unique(y)) < 2 || length(unique(p)) < 2) {
    return(list(roc = NULL, auc = NA_real_, ci_l = NA_real_, ci_u = NA_real_, auc_ci = NA_character_))
  }
  
  roc_obj <- pROC::roc(response = y, predictor = p, levels = c(0, 1), direction = "<", quiet = TRUE)
  auc_value <- as.numeric(pROC::auc(roc_obj))
  ci_obj <- as.numeric(pROC::ci.auc(roc_obj))
  
  list(
    roc = roc_obj,
    auc = auc_value,
    ci_l = ci_obj[1],
    ci_u = ci_obj[3],
    auc_ci = paste0(sprintf("%.3f", auc_value), " (", sprintf("%.3f", ci_obj[1]), "-", sprintf("%.3f", ci_obj[3]), ")")
  )
}

get_youden_threshold <- function(y, p) {
  roc_info <- get_auc_ci(y, p)
  if (is.null(roc_info$roc)) return(NA_real_)
  as.numeric(pROC::coords(roc_info$roc, x = "best", best.method = "youden", ret = "threshold", transpose = FALSE)[1])
}

classification_metrics <- function(y, p, threshold) {
  idx <- !is.na(y) & !is.na(p)
  y <- y[idx]
  p <- p[idx]
  
  if (length(y) == 0 || is.na(threshold)) {
    return(data.frame(Sensitivity = NA_real_, Specificity = NA_real_, Accuracy = NA_real_, PPV = NA_real_, NPV = NA_real_, TP = NA_integer_, TN = NA_integer_, FP = NA_integer_, FN = NA_integer_))
  }
  
  pred <- ifelse(p >= threshold, 1, 0)
  TP <- sum(pred == 1 & y == 1)
  TN <- sum(pred == 0 & y == 0)
  FP <- sum(pred == 1 & y == 0)
  FN <- sum(pred == 0 & y == 1)
  
  data.frame(
    Sensitivity = ifelse(TP + FN == 0, NA, TP / (TP + FN)),
    Specificity = ifelse(TN + FP == 0, NA, TN / (TN + FP)),
    Accuracy = (TP + TN) / length(y),
    PPV = ifelse(TP + FP == 0, NA, TP / (TP + FP)),
    NPV = ifelse(TN + FN == 0, NA, TN / (TN + FN)),
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN
  )
}

evaluate_model <- function(df, pred_col, model_name, dataset_name, threshold) {
  roc_info <- get_auc_ci(df[[label_col]], df[[pred_col]])
  cls <- classification_metrics(df[[label_col]], df[[pred_col]], threshold)
  
  data.frame(
    Model = model_name,
    Dataset = dataset_name,
    N = sum(!is.na(df[[pred_col]])),
    N_Label0 = sum(df[[label_col]] == 0, na.rm = TRUE),
    N_Label1 = sum(df[[label_col]] == 1, na.rm = TRUE),
    AUC = roc_info$auc,
    AUC_CI_Lower = roc_info$ci_l,
    AUC_CI_Upper = roc_info$ci_u,
    AUC_95CI = roc_info$auc_ci,
    Threshold = threshold,
    cls
  )
}

bootstrap_auc_optimism <- function(train_df, predictors, B = 1000) {
  y_orig <- train_df[[label_col]]
  fit_app <- glm(safe_formula(label_col, predictors), data = train_df, family = binomial())
  p_app <- as.numeric(predict(fit_app, train_df, type = "response"))
  auc_app <- get_auc_ci(y_orig, p_app)$auc
  
  optimism <- c()
  boot_auc <- c()
  original_auc <- c()
  n <- nrow(train_df)
  
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    boot_df <- train_df[idx, , drop = FALSE]
    if (length(unique(boot_df[[label_col]])) < 2) next
    
    fit_b <- tryCatch(glm(safe_formula(label_col, predictors), data = boot_df, family = binomial()), error = function(e) NULL)
    if (is.null(fit_b)) next
    
    p_boot <- as.numeric(predict(fit_b, boot_df, type = "response"))
    p_orig <- as.numeric(predict(fit_b, train_df, type = "response"))
    auc_boot <- get_auc_ci(boot_df[[label_col]], p_boot)$auc
    auc_orig <- get_auc_ci(train_df[[label_col]], p_orig)$auc
    
    if (is.finite(auc_boot) && is.finite(auc_orig)) {
      boot_auc <- c(boot_auc, auc_boot)
      original_auc <- c(original_auc, auc_orig)
      optimism <- c(optimism, auc_boot - auc_orig)
    }
  }
  
  mean_optimism <- mean(optimism, na.rm = TRUE)
  
  list(
    summary = data.frame(Apparent_AUC = auc_app, Optimism = mean_optimism, Optimism_Corrected_AUC = auc_app - mean_optimism, Successful_Bootstrap = length(optimism)),
    detail = data.frame(Iteration = seq_along(optimism), Bootstrap_AUC = boot_auc, Original_AUC = original_auc, Optimism = optimism)
  )
}

delong_matrix <- function(df, model_list) {
  models <- names(model_list)
  mat <- matrix(NA_real_, nrow = length(models), ncol = length(models), dimnames = list(models, models))
  rocs <- list()
  
  for (mn in models) {
    pred_col <- paste0(mn, "_Pred")
    rocs[[mn]] <- get_auc_ci(df[[label_col]], df[[pred_col]])$roc
  }
  
  for (i in seq_along(models)) {
    for (j in seq_along(models)) {
      if (i == j) {
        mat[i, j] <- 1
      } else if (!is.null(rocs[[models[i]]]) && !is.null(rocs[[models[j]]])) {
        mat[i, j] <- tryCatch(pROC::roc.test(rocs[[models[i]]], rocs[[models[j]]], method = "delong")$p.value, error = function(e) NA_real_)
      }
    }
  }
  
  mat
}


# ============================================================
# 6. Model visualization functions
# ============================================================

plot_lasso_coef_bar <- function(coef_df, output_pdf, output_png) {
  plot_df <- coef_df %>%
    dplyr::filter(Feature != "(Intercept)", Coefficient != 0) %>%
    dplyr::mutate(
      Region = factor(Region, levels = c("Clinical", "Intra", "Peri")),
      Feature_Wrapped = stringr::str_wrap(Feature, 38)
    )
  
  if (nrow(plot_df) == 0) return(NULL)
  
  p <- ggplot(plot_df, aes(x = reorder(Feature_Wrapped, Coefficient), y = Coefficient, fill = Region)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
    geom_col(color = "black", linewidth = 0.2, width = 0.72) +
    coord_flip() +
    facet_grid(Region ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("Clinical" = "#E64B35FF", "Intra" = "#4DBBD5FF", "Peri" = "#00A087FF")) +
    labs(title = "LASSO lambda.1se selected variables", x = NULL, y = "LASSO coefficient") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.y = element_text(size = 8), legend.position = "none", strip.background = element_rect(fill = "grey95", color = "grey70"), strip.text = element_text(face = "bold"))
  
  height_use <- max(7, 0.23 * nrow(plot_df) + 3)
  ggsave(output_pdf, p, width = 10, height = height_use)
  ggsave(output_png, p, width = 10, height = height_use, dpi = 300)
  p
}

plot_three_score_violin <- function(df, output_pdf, output_png) {
  plot_df <- df %>%
    dplyr::filter(.data[[cohort_col]] %in% c(train_value, test_value)) %>%
    dplyr::mutate(
      Dataset = ifelse(.data[[cohort_col]] == train_value, "Train", "Test"),
      Label_Group = ifelse(.data[[label_col]] == 1, "Fast growth", "Slow growth")
    ) %>%
    dplyr::select(all_of(c(id_col, label_col, cohort_col)), Dataset, Label_Group, Clinical_Score, Intra_Radscore, Peri_Radscore) %>%
    tidyr::pivot_longer(cols = c(Clinical_Score, Intra_Radscore, Peri_Radscore), names_to = "Score", values_to = "Score_Value") %>%
    dplyr::mutate(
      Score = dplyr::case_when(
        Score == "Clinical_Score" ~ "Clinical Score",
        Score == "Intra_Radscore" ~ "Intra-Radscore",
        Score == "Peri_Radscore" ~ "Peri-Radscore"
      ),
      Dataset = factor(Dataset, levels = c("Train", "Test")),
      Label_Group = factor(Label_Group, levels = c("Slow growth", "Fast growth")),
      Score = factor(Score, levels = c("Clinical Score", "Intra-Radscore", "Peri-Radscore"))
    )
  
  pval_df <- plot_df %>%
    dplyr::group_by(Dataset, Score) %>%
    dplyr::summarise(
      P = tryCatch(wilcox.test(Score_Value ~ Label_Group, exact = FALSE)$p.value, error = function(e) NA_real_),
      Y = max(Score_Value, na.rm = TRUE) + 0.1 * diff(range(Score_Value, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(P_label = ifelse(is.na(P), "P = NA", ifelse(P < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", P)))))
  
  pval_df$Y[!is.finite(pval_df$Y)] <- 0
  
  p <- ggplot(plot_df, aes(x = Label_Group, y = Score_Value, fill = Label_Group)) +
    geom_violin(trim = FALSE, alpha = 0.55, color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85, color = "black", linewidth = 0.35) +
    geom_jitter(aes(color = Label_Group), width = 0.10, size = 0.85, alpha = 0.40, show.legend = FALSE) +
    geom_text(data = pval_df, aes(x = 1.5, y = Y, label = P_label), inherit.aes = FALSE, size = 3.3, fontface = "bold") +
    facet_grid(Score ~ Dataset, scales = "free_y") +
    scale_fill_manual(values = c("Slow growth" = "#4DBBD5FF", "Fast growth" = "#E64B35FF")) +
    scale_color_manual(values = c("Slow growth" = "#4DBBD5FF", "Fast growth" = "#E64B35FF")) +
    labs(title = "Distribution of the three scores by growth group", x = NULL, y = "Score value") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "bottom", strip.background = element_rect(fill = "grey95", color = "grey70"), strip.text = element_text(face = "bold"))
  
  ggsave(output_pdf, p, width = 9.5, height = 8.8)
  ggsave(output_png, p, width = 9.5, height = 8.8, dpi = 300)
  list(plot = p, pvalue = pval_df)
}

plot_bootstrap_visualization <- function(bootstrap_summary, bootstrap_detail, output_pdf, output_png, model_colors) {
  summary_df <- bootstrap_summary %>% dplyr::mutate(Model = factor(Model, levels = names(model_colors)))
  detail_df <- bootstrap_detail %>% dplyr::mutate(Model = factor(Model, levels = names(model_colors)))
  
  p1 <- ggplot(summary_df, aes(x = Model, y = Apparent_AUC, fill = Model)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.25, alpha = 0.85) +
    geom_point(aes(y = Optimism_Corrected_AUC), size = 3, color = "black") +
    geom_segment(aes(x = Model, xend = Model, y = Apparent_AUC, yend = Optimism_Corrected_AUC), linewidth = 0.7, color = "black") +
    scale_fill_manual(values = model_colors) +
    coord_cartesian(ylim = c(0.45, 1.0)) +
    labs(title = "Bootstrap internal validation", subtitle = "Bars: apparent AUC; black dots: optimism-corrected AUC", x = NULL, y = "AUC") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
  
  p2 <- ggplot(detail_df, aes(x = Model, y = Optimism, fill = Model)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, color = "black", linewidth = 0.35) +
    geom_jitter(width = 0.12, size = 0.45, alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = model_colors) +
    labs(title = "Bootstrap optimism distribution", x = NULL, y = "Optimism") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
  
  grDevices::pdf(output_pdf, width = 11, height = 9, useDingbats = FALSE)
  print(p1)
  print(p2)
  grDevices::dev.off()
  
  grDevices::png(output_png, width = 3300, height = 2700, res = 300)
  print(p1)
  print(p2)
  grDevices::dev.off()
}

plot_auc_fusion_bar <- function(metrics_df, bootstrap_summary_df, output_pdf, output_png) {
  train_auc_df <- metrics_df %>% dplyr::filter(Dataset == "Train") %>% dplyr::select(Model, AUC) %>% dplyr::mutate(AUC_Type = "Training AUC")
  boot_auc_df <- bootstrap_summary_df %>% dplyr::select(Model, Optimism_Corrected_AUC) %>% dplyr::rename(AUC = Optimism_Corrected_AUC) %>% dplyr::mutate(AUC_Type = "Bootstrap-corrected AUC")
  test_auc_df <- metrics_df %>% dplyr::filter(Dataset == "Test") %>% dplyr::select(Model, AUC) %>% dplyr::mutate(AUC_Type = "Test AUC")
  
  plot_df <- dplyr::bind_rows(train_auc_df, boot_auc_df, test_auc_df) %>%
    dplyr::mutate(
      Model = factor(Model, levels = c("Clinical", "Intratumor", "Peritumor", "Clinical_Intratumor", "Clinical_Peritumor", "Intra_Peri", "Clinical_Intra_Peri")),
      AUC_Type = factor(AUC_Type, levels = c("Training AUC", "Bootstrap-corrected AUC", "Test AUC")),
      AUC_Label = sprintf("%.3f", AUC)
    )
  
  p <- ggplot(plot_df, aes(x = Model, y = AUC, fill = AUC_Type)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.68, color = "black", linewidth = 0.25) +
    geom_text(aes(label = AUC_Label), position = position_dodge(width = 0.78), vjust = -0.35, size = 3.1, fontface = "bold") +
    coord_cartesian(ylim = c(0.45, 0.90)) +
    scale_fill_manual(values = c("Training AUC" = "#4DBBD5FF", "Bootstrap-corrected AUC" = "#00A087FF", "Test AUC" = "#E64B35FF")) +
    labs(title = "AUC comparison across training, bootstrap internal validation, and test cohorts", x = NULL, y = "AUC", fill = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 35, hjust = 1), axis.text = element_text(color = "black"), axis.title = element_text(color = "black", face = "bold"), legend.position = "bottom")
  
  ggsave(output_pdf, p, width = 12, height = 6.8)
  ggsave(output_png, p, width = 12, height = 6.8, dpi = 300)
  plot_df
}

plot_roc_7 <- function(df, dataset_name, output_pdf, output_png, model_list, model_colors) {
  draw_roc <- function() {
    first <- TRUE
    legends <- c()
    cols <- c()
    for (mn in names(model_list)) {
      pred_col <- paste0(mn, "_Pred")
      roc_info <- get_auc_ci(df[[label_col]], df[[pred_col]])
      if (is.null(roc_info$roc)) next
      if (first) {
        plot(roc_info$roc, legacy.axes = TRUE, col = model_colors[mn], lwd = 2.2, main = paste0(dataset_name, " ROC curves"))
        first <- FALSE
      } else {
        plot(roc_info$roc, add = TRUE, col = model_colors[mn], lwd = 2.2)
      }
      legends <- c(legends, paste0(mn, ": ", roc_info$auc_ci))
      cols <- c(cols, model_colors[mn])
    }
    abline(0, 1, lty = 2, col = "grey60")
    legend("bottomright", legend = legends, col = cols, lwd = 2.2, cex = 0.72, bty = "n")
  }
  
  save_base_pdf_png(draw_roc, output_pdf, output_png, pdf_width = 7.2, pdf_height = 6.8, png_width = 2200, png_height = 2000)
}

plot_delong_heatmap_both <- function(mat, output_pdf, output_png, title_text) {
  mat_plot <- mat
  mat_plot[is.na(mat_plot)] <- 1
  display_mat <- matrix(ifelse(is.na(mat), "NA", ifelse(mat < 0.001, "<0.001", sprintf("%.3f", mat))), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  
  draw_heatmap <- function() {
    pheatmap::pheatmap(
      mat_plot,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = display_mat,
      fontsize_number = 8,
      main = title_text,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
      breaks = seq(0, 1, length.out = 101),
      na_col = "grey90",
      silent = FALSE
    )
  }
  
  grDevices::pdf(output_pdf, width = 8.8, height = 8.0, useDingbats = FALSE)
  draw_heatmap()
  grDevices::dev.off()
  
  grDevices::png(output_png, width = 2600, height = 2400, res = 300)
  draw_heatmap()
  grDevices::dev.off()
}


# ============================================================
# 7. Calibration and DCA helper functions
# ============================================================

make_logistic_calibration_curve <- function(y, p, curve_name = "Logistic calibration", grid = seq(0.001, 0.999, 0.002)) {
  idx <- !is.na(y) & !is.na(p)
  y <- y[idx]
  p <- clip_prob(p[idx])
  df_cal <- data.frame(y = y, lp = logit_prob(p))
  fit <- tryCatch(glm(y ~ lp, data = df_cal, family = binomial()), error = function(e) NULL)
  if (is.null(fit)) return(data.frame(Predicted = numeric(0), Observed = numeric(0), Curve = character(0)))
  data.frame(Predicted = grid, Observed = clip_prob(predict(fit, newdata = data.frame(lp = logit_prob(grid)), type = "response")), Curve = curve_name)
}

make_nonparametric_calibration_curve <- function(y, p, curve_name = "Nonparametric", grid = seq(0.001, 0.999, 0.002)) {
  idx <- !is.na(y) & !is.na(p)
  y <- y[idx]
  p <- clip_prob(p[idx])
  
  if (length(unique(y)) < 2 || length(unique(p)) < 5) {
    return(make_logistic_calibration_curve(y, p, curve_name, grid))
  }
  
  fit <- tryCatch(
    loess(y ~ p, span = 0.75, degree = 1, family = "symmetric", control = loess.control(surface = "direct")),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(make_logistic_calibration_curve(y, p, curve_name, grid))
  
  obs <- tryCatch(predict(fit, newdata = data.frame(p = grid)), error = function(e) rep(NA_real_, length(grid)))
  data.frame(Predicted = grid, Observed = pmin(pmax(obs, 0), 1), Curve = curve_name) %>% dplyr::filter(!is.na(Observed))
}

make_bootstrap_bias_corrected_calibration <- function(train_df, predictors, pred_col, B = 300, grid = seq(0.01, 0.99, 0.005)) {
  y_train <- train_df[[label_col]]
  p_train <- clip_prob(train_df[[pred_col]])
  apparent_curve <- make_nonparametric_calibration_curve(y_train, p_train, "Apparent", grid)
  
  if (nrow(apparent_curve) == 0) return(data.frame(Predicted = numeric(0), Observed = numeric(0), Curve = character(0)))
  
  n <- nrow(train_df)
  optimism_mat <- matrix(NA_real_, nrow = length(grid), ncol = B)
  set.seed(seed_split)
  
  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    boot_df <- train_df[idx, , drop = FALSE]
    if (length(unique(boot_df[[label_col]])) < 2) next
    
    fit_b <- tryCatch(glm(safe_formula(label_col, predictors), data = boot_df, family = binomial()), error = function(e) NULL)
    if (is.null(fit_b)) next
    
    p_boot_boot <- clip_prob(predict(fit_b, newdata = boot_df, type = "response"))
    p_boot_orig <- clip_prob(predict(fit_b, newdata = train_df, type = "response"))
    
    curve_boot <- make_nonparametric_calibration_curve(boot_df[[label_col]], p_boot_boot, "boot", grid)
    curve_orig <- make_nonparametric_calibration_curve(train_df[[label_col]], p_boot_orig, "orig", grid)
    if (nrow(curve_boot) == 0 || nrow(curve_orig) == 0) next
    
    boot_obs <- approx(curve_boot$Predicted, curve_boot$Observed, xout = grid, rule = 2)$y
    orig_obs <- approx(curve_orig$Predicted, curve_orig$Observed, xout = grid, rule = 2)$y
    optimism_mat[, b] <- boot_obs - orig_obs
  }
  
  mean_optimism <- rowMeans(optimism_mat, na.rm = TRUE)
  app_obs <- approx(apparent_curve$Predicted, apparent_curve$Observed, xout = grid, rule = 2)$y
  data.frame(Predicted = grid, Observed = pmin(pmax(app_obs - mean_optimism, 0), 1), Curve = "Bias-corrected")
}

calc_calibration_index <- function(y, p) {
  idx <- !is.na(y) & !is.na(p)
  y <- as.numeric(y[idx])
  p <- clip_prob(p[idx])
  
  auc <- get_auc_ci(y, p)$auc
  dxy <- 2 * (auc - 0.5)
  brier <- mean((p - y)^2)
  
  cal_df <- data.frame(y = y, lp = logit_prob(p))
  cal_fit <- tryCatch(glm(y ~ lp, data = cal_df, family = binomial()), error = function(e) NULL)
  
  intercept <- slope <- s_z <- s_p <- NA_real_
  if (!is.null(cal_fit)) {
    intercept <- unname(coef(cal_fit)[1])
    slope <- unname(coef(cal_fit)[2])
    se <- tryCatch(summary(cal_fit)$coefficients[2, "Std. Error"], error = function(e) NA_real_)
    if (is.finite(se) && se > 0) {
      s_z <- (slope - 1) / se
      s_p <- 2 * pnorm(-abs(s_z))
    }
  }
  
  r2 <- D <- U <- Q <- NA_real_
  emax <- e90 <- eavg <- NA_real_
  
  if (requireNamespace("rms", quietly = TRUE)) {
    vp <- tryCatch(rms::val.prob(p = p, y = y, m = 10, pl = FALSE), error = function(e) NULL)
    if (!is.null(vp)) {
      vp_names <- names(vp)
      if ("R2" %in% vp_names) r2 <- as.numeric(vp["R2"])
      if ("D" %in% vp_names) D <- as.numeric(vp["D"])
      if ("U" %in% vp_names) U <- as.numeric(vp["U"])
      if ("Q" %in% vp_names) Q <- as.numeric(vp["Q"])
      if ("Brier" %in% vp_names) brier <- as.numeric(vp["Brier"])
      if ("Intercept" %in% vp_names) intercept <- as.numeric(vp["Intercept"])
      if ("Slope" %in% vp_names) slope <- as.numeric(vp["Slope"])
      if ("Emax" %in% vp_names) emax <- as.numeric(vp["Emax"])
      if ("E90" %in% vp_names) e90 <- as.numeric(vp["E90"])
      if ("Eavg" %in% vp_names) eavg <- as.numeric(vp["Eavg"])
      if ("S:z" %in% vp_names) s_z <- as.numeric(vp["S:z"])
      if ("S:p" %in% vp_names) s_p <- as.numeric(vp["S:p"])
    }
  }
  
  np <- make_nonparametric_calibration_curve(y, p, "Nonparametric", seq(0.01, 0.99, 0.01))
  if (nrow(np) > 0) {
    err <- abs(np$Observed - np$Predicted)
    emax <- max(err, na.rm = TRUE)
    e90 <- as.numeric(quantile(err, 0.90, na.rm = TRUE))
    eavg <- mean(err, na.rm = TRUE)
  }
  
  data.frame(Dxy = dxy, C_ROC = auc, R2 = r2, D = D, U = U, Q = Q, Brier = brier, Intercept = intercept, Slope = slope, Emax = emax, E90 = e90, Eavg = eavg, S_z = s_z, S_p = s_p)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x) | !is.finite(x), "NA", sprintf(paste0("%.", digits, "f"), x))
}

calculate_dca <- function(y, p, thresholds) {
  idx <- !is.na(y) & !is.na(p)
  y <- y[idx]
  p <- clip_prob(p[idx])
  n <- length(y)
  prevalence <- mean(y == 1)
  
  out <- lapply(thresholds, function(pt) {
    pred_pos <- p >= pt
    TP <- sum(pred_pos & y == 1)
    FP <- sum(pred_pos & y == 0)
    data.frame(
      Threshold = pt,
      Net_Benefit_Model = (TP / n) - (FP / n) * (pt / (1 - pt)),
      Net_Benefit_All = prevalence - (1 - prevalence) * (pt / (1 - pt)),
      Net_Benefit_None = 0
    )
  })
  dplyr::bind_rows(out)
}

calculate_clinical_impact <- function(y, p, thresholds) {
  idx <- !is.na(y) & !is.na(p)
  y <- y[idx]
  p <- clip_prob(p[idx])
  n <- length(y)
  
  out <- lapply(thresholds, function(pt) {
    high <- p >= pt
    data.frame(
      Threshold = pt,
      Number_High_Risk = sum(high),
      Number_True_Positive = sum(high & y == 1),
      Rate_High_Risk = sum(high) / n,
      Rate_True_Positive = sum(high & y == 1) / n
    )
  })
  dplyr::bind_rows(out)
}


# ============================================================
# 8. Read data, split cohorts, and preprocess features
# ============================================================

clinical_raw <- read_input_data(clinical_file, clinical_sheet)
intra_raw <- read_input_data(intra_file, intra_sheet)
peri_raw <- read_input_data(peri_file, peri_sheet)

check_same_patients(clinical_raw, intra_raw, "Clinical", "Intra")
check_same_patients(clinical_raw, peri_raw, "Clinical", "Peri")

split_df <- make_seed1111_split(clinical_raw)

clinical_raw <- apply_split(clinical_raw, split_df)
intra_raw <- apply_split(intra_raw, split_df)
peri_raw <- apply_split(peri_raw, split_df)

split_table <- split_df %>%
  dplyr::mutate(Dataset = dplyr::case_when(Cohort == train_value ~ "Training", Cohort == test_value ~ "Test", Cohort == validation_value ~ "Validation")) %>%
  dplyr::count(Dataset, .data[[label_col]])

openxlsx::write.xlsx(
  list(patient_split_seed1111 = split_df, split_table = split_table),
  file.path(preprocess_dir, "PatientID_split_seed1111.xlsx"),
  overwrite = TRUE
)

clinical_obj <- preprocess_and_select(clinical_raw, "Clinical", is_clinical = TRUE)
intra_obj <- preprocess_and_select(intra_raw, "Intra", is_clinical = FALSE)
peri_obj <- preprocess_and_select(peri_raw, "Peri", is_clinical = FALSE)

clinical_data <- clinical_obj$data
intra_data <- intra_obj$data
peri_data <- peri_obj$data

openxlsx::write.xlsx(
  list(Clinical_data_processed = clinical_data, Intra_data_processed = intra_data, Peri_data_processed = peri_data, split_table = split_table),
  file.path(preprocess_dir, "Processed_Clinical_Intra_Peri_Data.xlsx"),
  overwrite = TRUE
)

wb_lasso <- openxlsx::createWorkbook()
for (obj_name in c("clinical_obj", "intra_obj", "peri_obj")) {
  obj <- get(obj_name)
  prefix <- tools::toTitleCase(gsub("_obj", "", obj_name))
  add_sheet(wb_lasso, paste0(prefix, "_imputation"), obj$imputation)
  add_sheet(wb_lasso, paste0(prefix, "_Pearson"), obj$pearson_table)
  add_sheet(wb_lasso, paste0(prefix, "_Pearson_removed"), obj$pearson_removed)
  add_sheet(wb_lasso, paste0(prefix, "_highcorr_removed"), obj$high_corr_removed)
  add_sheet(wb_lasso, paste0(prefix, "_highcorr_pairs"), obj$high_corr_pairs)
  add_sheet(wb_lasso, paste0(prefix, "_VIF_removed"), obj$vif_removed)
  add_sheet(wb_lasso, paste0(prefix, "_final_VIF"), obj$final_vif)
  add_sheet(wb_lasso, paste0(prefix, "_Zscore_parameters"), obj$zscore_parameters)
  add_sheet(wb_lasso, paste0(prefix, "_after_pearson"), data.frame(Feature = obj$retained_after_pearson))
  add_sheet(wb_lasso, paste0(prefix, "_after_highcorr"), data.frame(Feature = obj$retained_after_highcorr))
  add_sheet(wb_lasso, paste0(prefix, "_before_lasso"), data.frame(Feature = obj$retained_before_lasso))
  add_sheet(wb_lasso, paste0(prefix, "_lasso_coef_1se"), obj$lasso$coef_1se)
  add_sheet(wb_lasso, paste0(prefix, "_lasso_selected_1se"), obj$lasso$selected_1se)
  add_sheet(wb_lasso, paste0(prefix, "_lasso_coef_min"), obj$lasso$coef_min)
  add_sheet(wb_lasso, paste0(prefix, "_lasso_selected_min"), obj$lasso$selected_min)
}
openxlsx::saveWorkbook(wb_lasso, file.path(lasso_dir, "LASSO_and_Filtering_Results.xlsx"), overwrite = TRUE)

lasso_coef_bar_df <- dplyr::bind_rows(
  clinical_obj$lasso$selected_1se %>% dplyr::mutate(Region = "Clinical"),
  intra_obj$lasso$selected_1se %>% dplyr::mutate(Region = "Intra"),
  peri_obj$lasso$selected_1se %>% dplyr::mutate(Region = "Peri")
)

plot_lasso_coef_bar(
  lasso_coef_bar_df,
  file.path(lasso_dir, "LASSO_lambda1se_Selected_Variables_Coefficient_Barplot.pdf"),
  file.path(lasso_dir, "LASSO_lambda1se_Selected_Variables_Coefficient_Barplot.png")
)


# ============================================================
# 9. Construct Clinical Score, Intra-Radscore, and Peri-Radscore
# ============================================================

add_prefix_except <- function(df, prefix) {
  vars <- setdiff(colnames(df), c(id_col, label_col, cohort_col))
  colnames(df)[match(vars, colnames(df))] <- paste0(prefix, vars)
  df
}

clinical_p <- add_prefix_except(clinical_data, "clin_")
intra_p <- add_prefix_except(intra_data, "intra_")
peri_p <- add_prefix_except(peri_data, "peri_")

base_df <- clinical_p[, c(id_col, label_col, cohort_col)]

model_df <- base_df %>%
  dplyr::left_join(clinical_p %>% dplyr::select(-all_of(c(label_col, cohort_col))), by = id_col) %>%
  dplyr::left_join(intra_p %>% dplyr::select(-all_of(c(label_col, cohort_col))), by = id_col) %>%
  dplyr::left_join(peri_p %>% dplyr::select(-all_of(c(label_col, cohort_col))), by = id_col)

clinical_selected <- intersect(paste0("clin_", clinical_obj$lasso$selected_1se$Feature), colnames(model_df))
intra_selected <- intersect(paste0("intra_", intra_obj$lasso$selected_1se$Feature), colnames(model_df))
peri_selected <- intersect(paste0("peri_", peri_obj$lasso$selected_1se$Feature), colnames(model_df))

intra_coef_1se <- intra_obj$lasso$coef_1se %>% dplyr::mutate(Feature = ifelse(Feature == "(Intercept)", Feature, paste0("intra_", Feature)))
peri_coef_1se <- peri_obj$lasso$coef_1se %>% dplyr::mutate(Feature = ifelse(Feature == "(Intercept)", Feature, paste0("peri_", Feature)))

intra_selected_df <- intra_obj$lasso$selected_1se %>% dplyr::mutate(Feature = paste0("intra_", Feature))
peri_selected_df <- peri_obj$lasso$selected_1se %>% dplyr::mutate(Feature = paste0("peri_", Feature))

train_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == train_value)

clinical_score_obj <- fit_logistic_score(
  train_df = train_df,
  full_df = model_df,
  selected_features = clinical_selected,
  score_name = "Clinical_Score"
)
model_df <- clinical_score_obj$data

model_df <- make_radscore_from_lasso(model_df, intra_selected_df, intra_coef_1se, "Intra_Radscore")
model_df <- make_radscore_from_lasso(model_df, peri_selected_df, peri_coef_1se, "Peri_Radscore")

train_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == train_value)
test_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == test_value)
valid_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == validation_value)

score_violin_obj <- plot_three_score_violin(
  model_df,
  file.path(score_dir, "ThreeScores_Train_Test_ViolinBoxplot_with_Pvalue.pdf"),
  file.path(score_dir, "ThreeScores_Train_Test_ViolinBoxplot_with_Pvalue.png")
)


# ============================================================
# 10. Build seven score-based models and bootstrap validation
# ============================================================

model_list <- list(
  Clinical = c("Clinical_Score"),
  Intratumor = c("Intra_Radscore"),
  Peritumor = c("Peri_Radscore"),
  Clinical_Intratumor = c("Clinical_Score", "Intra_Radscore"),
  Clinical_Peritumor = c("Clinical_Score", "Peri_Radscore"),
  Intra_Peri = c("Intra_Radscore", "Peri_Radscore"),
  Clinical_Intra_Peri = c("Clinical_Score", "Intra_Radscore", "Peri_Radscore")
)

model_colors <- c(
  Clinical = "#E64B35FF",
  Intratumor = "#4DBBD5FF",
  Peritumor = "#00A087FF",
  Clinical_Intratumor = "#3C5488FF",
  Clinical_Peritumor = "#F39B7FFF",
  Intra_Peri = "#8491B4FF",
  Clinical_Intra_Peri = "#91D1C2FF"
)

model_fits <- list()
metrics_list <- list()
bootstrap_summary_list <- list()
bootstrap_detail_list <- list()

for (mn in names(model_list)) {
  predictors <- model_list[[mn]]
  fit <- glm(safe_formula(label_col, predictors), data = train_df, family = binomial())
  model_fits[[mn]] <- fit
  
  pred_col <- paste0(mn, "_Pred")
  model_df[[pred_col]] <- as.numeric(predict(fit, newdata = model_df, type = "response"))
  
  train_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == train_value)
  test_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == test_value)
  valid_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == validation_value)
  
  threshold <- get_youden_threshold(train_df[[label_col]], train_df[[pred_col]])
  
  metrics_list[[mn]] <- dplyr::bind_rows(
    evaluate_model(train_df, pred_col, mn, "Train", threshold),
    evaluate_model(test_df, pred_col, mn, "Test", threshold)
  )
  
  boot <- bootstrap_auc_optimism(train_df = train_df, predictors = predictors, B = bootstrap_B)
  boot$summary$Model <- mn
  boot$detail$Model <- mn
  bootstrap_summary_list[[mn]] <- boot$summary
  bootstrap_detail_list[[mn]] <- boot$detail
}

metrics_df <- dplyr::bind_rows(metrics_list)
bootstrap_summary_df <- dplyr::bind_rows(bootstrap_summary_list) %>% dplyr::select(Model, everything())
bootstrap_detail_df <- dplyr::bind_rows(bootstrap_detail_list) %>% dplyr::select(Model, everything())

plot_bootstrap_visualization(
  bootstrap_summary_df,
  bootstrap_detail_df,
  file.path(bootstrap_dir, "Bootstrap_Internal_Validation_Visualization.pdf"),
  file.path(bootstrap_dir, "Bootstrap_Internal_Validation_Visualization.png"),
  model_colors
)

auc_fusion_df <- plot_auc_fusion_bar(
  metrics_df,
  bootstrap_summary_df,
  file.path(score_dir, "Train_BootstrapCorrected_Test_AUC_Fusion_Barplot.pdf"),
  file.path(score_dir, "Train_BootstrapCorrected_Test_AUC_Fusion_Barplot.png")
)

plot_roc_7(train_df, "Training cohort", file.path(score_dir, "Train_ROC_7models.pdf"), file.path(score_dir, "Train_ROC_7models.png"), model_list, model_colors)
plot_roc_7(test_df, "Test cohort", file.path(score_dir, "Test_ROC_7models.pdf"), file.path(score_dir, "Test_ROC_7models.png"), model_list, model_colors)

train_delong <- delong_matrix(train_df, model_list)
test_delong <- delong_matrix(test_df, model_list)

plot_delong_heatmap_both(train_delong, file.path(score_dir, "Train_DeLong_Heatmap_7models.pdf"), file.path(score_dir, "Train_DeLong_Heatmap_7models.png"), "Training cohort DeLong P values")
plot_delong_heatmap_both(test_delong, file.path(score_dir, "Test_DeLong_Heatmap_7models.pdf"), file.path(score_dir, "Test_DeLong_Heatmap_7models.png"), "Test cohort DeLong P values")

score_data <- model_df %>%
  dplyr::select(all_of(c(id_col, label_col, cohort_col)), Clinical_Score, Intra_Radscore, Peri_Radscore, dplyr::ends_with("_Pred"))

wb <- openxlsx::createWorkbook()
add_sheet(wb, "split_table", split_table)
add_sheet(wb, "patient_split_seed1111", split_df)
add_sheet(wb, "selected_features", lasso_coef_bar_df)
add_sheet(wb, "Clinical_selected_1se", clinical_obj$lasso$selected_1se)
add_sheet(wb, "Intra_selected_1se", intra_obj$lasso$selected_1se)
add_sheet(wb, "Peri_selected_1se", peri_obj$lasso$selected_1se)
add_sheet(wb, "metrics_7models_train_test", metrics_df)
add_sheet(wb, "bootstrap_summary", bootstrap_summary_df)
add_sheet(wb, "bootstrap_detail", bootstrap_detail_df)
add_sheet(wb, "AUC_fusion_barplot_data", auc_fusion_df)
add_sheet(wb, "train_delong", as.data.frame(train_delong) %>% tibble::rownames_to_column("Model"))
add_sheet(wb, "test_delong", as.data.frame(test_delong) %>% tibble::rownames_to_column("Model"))
add_sheet(wb, "score_violin_pvalues", score_violin_obj$pvalue)
add_sheet(wb, "score_prediction_data", score_data)

openxlsx::saveWorkbook(wb, file.path(output_dir, "Seed1111_SevenModels_All_Results.xlsx"), overwrite = TRUE)


# ============================================================
# 11. Final validation for Clinical_Intra_Peri model
# ============================================================

if (!final_pred_col %in% colnames(model_df)) {
  stop("model_df does not contain ", final_pred_col, ". Please check whether the seven-model analysis completed successfully.")
}

train_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == train_value)
test_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == test_value)
valid_df <- model_df %>% dplyr::filter(.data[[cohort_col]] == validation_value)

if (nrow(valid_df) == 0) stop("Validation cohort is empty. Please confirm that original Cohort=2 was retained.")

final_threshold <- get_youden_threshold(train_df[[label_col]], train_df[[final_pred_col]])

final_metrics <- dplyr::bind_rows(
  evaluate_model(train_df, final_pred_col, final_model, "Train", final_threshold),
  evaluate_model(test_df, final_pred_col, final_model, "Test", final_threshold),
  evaluate_model(valid_df, final_pred_col, final_model, "Validation", final_threshold)
)

print(final_metrics)

plot_final_roc_one <- function(df, dataset_name, output_pdf, output_png) {
  roc_info <- get_auc_ci(df[[label_col]], df[[final_pred_col]])
  draw_roc <- function() {
    plot(roc_info$roc, legacy.axes = TRUE, col = "#3C5488FF", lwd = 2.8, main = paste0(dataset_name, " ROC curve: ", final_model))
    abline(0, 1, lty = 2, col = "grey60")
    legend("bottomright", legend = paste0("AUC = ", roc_info$auc_ci), col = "#3C5488FF", lwd = 2.8, bty = "n", cex = 0.95)
  }
  save_base_pdf_png(draw_roc, output_pdf, output_png, pdf_width = 7.0, pdf_height = 6.5, png_width = 2100, png_height = 1950)
}

plot_final_roc_three <- function(output_pdf, output_png) {
  roc_train <- get_auc_ci(train_df[[label_col]], train_df[[final_pred_col]])
  roc_test <- get_auc_ci(test_df[[label_col]], test_df[[final_pred_col]])
  roc_valid <- get_auc_ci(valid_df[[label_col]], valid_df[[final_pred_col]])
  
  draw_roc <- function() {
    plot(roc_train$roc, legacy.axes = TRUE, col = "#4DBBD5FF", lwd = 2.5, main = paste0("ROC curves of ", final_model))
    plot(roc_test$roc, add = TRUE, col = "#00A087FF", lwd = 2.5)
    plot(roc_valid$roc, add = TRUE, col = "#E64B35FF", lwd = 2.5)
    abline(0, 1, lty = 2, col = "grey60")
    legend(
      "bottomright",
      legend = c(paste0("Train: AUC ", roc_train$auc_ci), paste0("Test: AUC ", roc_test$auc_ci), paste0("Validation: AUC ", roc_valid$auc_ci)),
      col = c("#4DBBD5FF", "#00A087FF", "#E64B35FF"),
      lwd = 2.5,
      bty = "n",
      cex = 0.78
    )
  }
  
  save_base_pdf_png(draw_roc, output_pdf, output_png, pdf_width = 7.4, pdf_height = 6.8, png_width = 2300, png_height = 2100)
}

plot_final_roc_one(valid_df, "Validation cohort", file.path(final_validation_dir, "Validation_ROC_Clinical_Intra_Peri.pdf"), file.path(final_validation_dir, "Validation_ROC_Clinical_Intra_Peri.png"))
plot_final_roc_three(file.path(final_validation_dir, "Train_Test_Validation_ROC_Clinical_Intra_Peri.pdf"), file.path(final_validation_dir, "Train_Test_Validation_ROC_Clinical_Intra_Peri.png"))


# ============================================================
# 12. Final model calibration plots
# ============================================================

grid_cal <- seq(0.01, 0.99, 0.005)

cal_train_app <- make_nonparametric_calibration_curve(train_df[[label_col]], train_df[[final_pred_col]], "Training apparent", grid_cal)
cal_train_boot <- make_bootstrap_bias_corrected_calibration(train_df, final_predictors, final_pred_col, calibration_bootstrap_B, grid_cal) %>% dplyr::mutate(Curve = "Bootstrap-corrected training")
cal_test <- make_nonparametric_calibration_curve(test_df[[label_col]], test_df[[final_pred_col]], "Test", grid_cal)
cal_valid_logistic <- make_logistic_calibration_curve(valid_df[[label_col]], valid_df[[final_pred_col]], "Validation logistic", grid_cal)
cal_valid_np <- make_nonparametric_calibration_curve(valid_df[[label_col]], valid_df[[final_pred_col]], "Validation nonparametric", grid_cal)

cal_all_df <- dplyr::bind_rows(cal_train_app, cal_train_boot, cal_test, cal_valid_logistic, cal_valid_np)

cal_index_df <- dplyr::bind_rows(
  calc_calibration_index(train_df[[label_col]], train_df[[final_pred_col]]) %>% dplyr::mutate(Dataset = "Train"),
  calc_calibration_index(test_df[[label_col]], test_df[[final_pred_col]]) %>% dplyr::mutate(Dataset = "Test"),
  calc_calibration_index(valid_df[[label_col]], valid_df[[final_pred_col]]) %>% dplyr::mutate(Dataset = "Validation")
) %>% dplyr::select(Dataset, everything())

p_cal_all <- ggplot(cal_all_df, aes(x = Predicted, y = Observed, color = Curve, linetype = Curve)) +
  geom_abline(aes(color = "Ideal", linetype = "Ideal"), intercept = 0, slope = 1, linewidth = 0.75) +
  geom_line(linewidth = 1.1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_color_manual(values = c("Ideal" = "black", "Training apparent" = "#4DBBD5FF", "Bootstrap-corrected training" = "#00A087FF", "Test" = "#3C5488FF", "Validation logistic" = "#E64B35FF", "Validation nonparametric" = "#F39B7FFF")) +
  scale_linetype_manual(values = c("Ideal" = "dashed", "Training apparent" = "solid", "Bootstrap-corrected training" = "solid", "Test" = "solid", "Validation logistic" = "solid", "Validation nonparametric" = "solid")) +
  labs(title = "Calibration curves of Clinical_Intra_Peri model", x = "Predicted Probability", y = "Observed Probability", color = NULL, linetype = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.line = element_line(color = "black"), axis.text = element_text(color = "black"), axis.title = element_text(color = "black", face = "bold"), legend.position = "bottom")

ggsave(file.path(final_validation_dir, "Calibration_Train_BootstrapCorrected_Test_Validation_Clinical_Intra_Peri.pdf"), p_cal_all, width = 8.2, height = 7.4)
ggsave(file.path(final_validation_dir, "Calibration_Train_BootstrapCorrected_Test_Validation_Clinical_Intra_Peri.png"), p_cal_all, width = 8.2, height = 7.4, dpi = 300)

# Publication-style calibration plot: Ideal + Apparent + Test + Bias-corrected with top rug
cal_app <- make_nonparametric_calibration_curve(train_df[[label_col]], train_df[[final_pred_col]], "Apparent", grid_cal)
cal_test_simple <- make_nonparametric_calibration_curve(test_df[[label_col]], test_df[[final_pred_col]], "test", grid_cal)
cal_bias <- make_bootstrap_bias_corrected_calibration(train_df, final_predictors, final_pred_col, calibration_bootstrap_B, grid_cal)

cal_plot_df <- dplyr::bind_rows(cal_app, cal_test_simple, cal_bias)
cal_plot_df$Curve <- factor(cal_plot_df$Curve, levels = c("Apparent", "test", "Bias-corrected"))

rug_df <- model_df %>%
  dplyr::filter(!is.na(.data[[final_pred_col]])) %>%
  dplyr::select(all_of(final_pred_col)) %>%
  dplyr::rename(Predicted = all_of(final_pred_col))

p_cal_final <- ggplot() +
  geom_abline(aes(color = "Ideal", linetype = "Ideal"), intercept = 0, slope = 1, linewidth = 0.85) +
  geom_line(data = cal_plot_df, aes(x = Predicted, y = Observed, color = Curve, linetype = Curve), linewidth = 1.25) +
  geom_rug(data = rug_df, aes(x = Predicted), inherit.aes = FALSE, sides = "t", color = "black", alpha = 0.65, linewidth = 0.35) +
  scale_color_manual(name = NULL, breaks = c("Ideal", "Apparent", "test", "Bias-corrected"), values = c("Ideal" = "black", "Apparent" = "#00A087FF", "test" = "#F39B7FFF", "Bias-corrected" = "#E64B35FF")) +
  scale_linetype_manual(name = NULL, breaks = c("Ideal", "Apparent", "test", "Bias-corrected"), values = c("Ideal" = "dashed", "Apparent" = "solid", "test" = "solid", "Bias-corrected" = "solid")) +
  scale_x_continuous(name = "Predicted Probability", limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(name = "Observed/Actual Probability", limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(title = "Calibration Plot", subtitle = "Clinical_Intra_Peri prediction model") +
  theme_classic(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA), axis.line = element_line(color = "black", linewidth = 0.65), axis.ticks = element_line(color = "black", linewidth = 0.55), axis.text = element_text(color = "black", size = 12), axis.title = element_text(color = "black", size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "black"), plot.subtitle = element_text(hjust = 0.5, size = 13, color = "black"), legend.position = c(0.74, 0.18), legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3), legend.key = element_rect(fill = "white", color = NA), legend.text = element_text(size = 11, color = "black"))

ggsave(file.path(calibration_plot_dir, "Calibration_Plot_Clinical_Intra_Peri_Ideal_Apparent_Test_BiasCorrected.pdf"), p_cal_final, width = 7.2, height = 6.6)
ggsave(file.path(calibration_plot_dir, "Calibration_Plot_Clinical_Intra_Peri_Ideal_Apparent_Test_BiasCorrected.png"), p_cal_final, width = 7.2, height = 6.6, dpi = 300)
openxlsx::write.xlsx(cal_plot_df, file.path(calibration_plot_dir, "Calibration_Plot_Clinical_Intra_Peri_Curve_Data.xlsx"), overwrite = TRUE)


# ============================================================
# 13. Validation-only calibration plot with logistic and nonparametric curves
# ============================================================

valid_df_cal <- valid_df %>% dplyr::filter(!is.na(.data[[label_col]]), !is.na(.data[[final_pred_col]]))
if (nrow(valid_df_cal) == 0) stop("Validation calibration data are empty.")
if (length(unique(valid_df_cal[[label_col]])) < 2) stop("Validation Label has only one class; calibration curve cannot be plotted.")

y_valid <- valid_df_cal[[label_col]]
p_valid <- clip_prob(valid_df_cal[[final_pred_col]])
grid_valid <- seq(0.001, 0.999, 0.002)

cal_logistic <- make_logistic_calibration_curve(y_valid, p_valid, "Logistic calibration", grid_valid)
cal_nonparametric <- make_nonparametric_calibration_curve(y_valid, p_valid, "Nonparametric", grid_valid)
cal_valid_plot_df <- dplyr::bind_rows(cal_logistic, cal_nonparametric)
cal_valid_plot_df$Curve <- factor(cal_valid_plot_df$Curve, levels = c("Logistic calibration", "Nonparametric"))

valid_cal_metric <- calc_calibration_index(y_valid, p_valid)
metric_text <- paste0(
  "Dxy: ", fmt_num(valid_cal_metric$Dxy), "\n",
  "C (ROC): ", fmt_num(valid_cal_metric$C_ROC), "\n",
  "R\u00b2: ", fmt_num(valid_cal_metric$R2), "\n",
  "D: ", fmt_num(valid_cal_metric$D), "\n",
  "U: ", fmt_num(valid_cal_metric$U), "\n",
  "Q: ", fmt_num(valid_cal_metric$Q), "\n",
  "Brier: ", fmt_num(valid_cal_metric$Brier), "\n",
  "Intercept: ", fmt_num(valid_cal_metric$Intercept), "\n",
  "Slope: ", fmt_num(valid_cal_metric$Slope), "\n",
  "Emax: ", fmt_num(valid_cal_metric$Emax), "\n",
  "E90: ", fmt_num(valid_cal_metric$E90), "\n",
  "Eavg: ", fmt_num(valid_cal_metric$Eavg), "\n",
  "S:z: ", fmt_num(valid_cal_metric$S_z), "\n",
  "S:p: ", fmt_num(valid_cal_metric$S_p)
)

p_validation_calibration <- ggplot() +
  geom_abline(aes(color = "Ideal", linetype = "Ideal"), intercept = 0, slope = 1, linewidth = 0.85) +
  geom_line(data = cal_valid_plot_df, aes(x = Predicted, y = Observed, color = Curve, linetype = Curve), linewidth = 1.25) +
  geom_rug(data = data.frame(Predicted = p_valid), aes(x = Predicted), inherit.aes = FALSE, sides = "b", color = "black", alpha = 0.70, linewidth = 0.35) +
  annotate("text", x = 0.035, y = 0.965, label = metric_text, hjust = 0, vjust = 1, size = 3.45, color = "black", lineheight = 0.95) +
  scale_color_manual(name = NULL, breaks = c("Ideal", "Logistic calibration", "Nonparametric"), values = c("Ideal" = "black", "Logistic calibration" = "#00A087FF", "Nonparametric" = "#E64B35FF")) +
  scale_linetype_manual(name = NULL, breaks = c("Ideal", "Logistic calibration", "Nonparametric"), values = c("Ideal" = "dashed", "Logistic calibration" = "solid", "Nonparametric" = "solid")) +
  scale_x_continuous(name = "Predicted Probability", limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(name = "Observed/Actual Probability", limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(title = "Validation Calibration Plot", subtitle = "Clinical_Intra_Peri prediction model") +
  theme_classic(base_size = 14) +
  theme(plot.background = element_rect(fill = "white", color = NA), panel.background = element_rect(fill = "white", color = NA), axis.line = element_line(color = "black", linewidth = 0.65), axis.ticks = element_line(color = "black", linewidth = 0.55), axis.text = element_text(color = "black", size = 12), axis.title = element_text(color = "black", size = 14, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "black"), plot.subtitle = element_text(hjust = 0.5, size = 13, color = "black"), legend.position = c(0.73, 0.16), legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3), legend.key = element_rect(fill = "white", color = NA), legend.text = element_text(size = 11, color = "black"))

ggsave(file.path(validation_calibration_dir, "Validation_Calibration_Plot_Clinical_Intra_Peri_Logistic_Nonparametric.pdf"), p_validation_calibration, width = 7.2, height = 6.8)
ggsave(file.path(validation_calibration_dir, "Validation_Calibration_Plot_Clinical_Intra_Peri_Logistic_Nonparametric.png"), p_validation_calibration, width = 7.2, height = 6.8, dpi = 300)

openxlsx::write.xlsx(
  list(
    calibration_curve_data = cal_valid_plot_df,
    calibration_metrics = valid_cal_metric,
    validation_prediction_data = valid_df_cal %>% dplyr::select(all_of(c(id_col, label_col, cohort_col)), Clinical_Score, Intra_Radscore, Peri_Radscore, all_of(final_pred_col))
  ),
  file.path(validation_calibration_dir, "Validation_Calibration_Plot_Clinical_Intra_Peri_Data.xlsx"),
  overwrite = TRUE
)


# ============================================================
# 14. Decision curve analysis and clinical impact curve
# ============================================================

thresholds_dca <- seq(0.01, 0.99, 0.01)

dca_all <- dplyr::bind_rows(
  calculate_dca(train_df[[label_col]], train_df[[final_pred_col]], thresholds_dca) %>% dplyr::mutate(Dataset = "Train"),
  calculate_dca(test_df[[label_col]], test_df[[final_pred_col]], thresholds_dca) %>% dplyr::mutate(Dataset = "Test"),
  calculate_dca(valid_df[[label_col]], valid_df[[final_pred_col]], thresholds_dca) %>% dplyr::mutate(Dataset = "Validation")
)

dca_long <- dca_all %>%
  dplyr::select(Dataset, Threshold, Model = Net_Benefit_Model, Treat_All = Net_Benefit_All, Treat_None = Net_Benefit_None) %>%
  tidyr::pivot_longer(cols = c(Model, Treat_All, Treat_None), names_to = "Strategy", values_to = "Net_Benefit")

p_dca <- ggplot(dca_long, aes(x = Threshold, y = Net_Benefit, color = Strategy, linetype = Strategy)) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.45) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ Dataset, nrow = 1) +
  scale_color_manual(values = c("Model" = "#3C5488FF", "Treat_All" = "#F39B7FFF", "Treat_None" = "grey45"), labels = c("Model" = final_model, "Treat_All" = "Treat all", "Treat_None" = "Treat none")) +
  scale_linetype_manual(values = c("Model" = "solid", "Treat_All" = "dashed", "Treat_None" = "dotted"), labels = c("Model" = final_model, "Treat_All" = "Treat all", "Treat_None" = "Treat none")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 0.5)) +
  labs(title = "Decision Curve Analysis", subtitle = final_model, x = "Threshold Probability", y = "Net Benefit", color = NULL, linetype = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5), axis.line = element_line(color = "black"), axis.text = element_text(color = "black"), axis.title = element_text(color = "black", face = "bold"), legend.position = "bottom", strip.background = element_rect(fill = "grey95", color = "grey70"), strip.text = element_text(face = "bold"))

ggsave(file.path(final_validation_dir, "DCA_Train_Test_Validation_Clinical_Intra_Peri_y_minus0.1_to_0.5.pdf"), p_dca, width = 11.5, height = 5.8)
ggsave(file.path(final_validation_dir, "DCA_Train_Test_Validation_Clinical_Intra_Peri_y_minus0.1_to_0.5.png"), p_dca, width = 11.5, height = 5.8, dpi = 300)

impact_valid <- calculate_clinical_impact(valid_df[[label_col]], valid_df[[final_pred_col]], thresholds_dca)

impact_long <- impact_valid %>%
  dplyr::select(Threshold, `High-risk patients` = Number_High_Risk, `True positives` = Number_True_Positive) %>%
  tidyr::pivot_longer(cols = c(`High-risk patients`, `True positives`), names_to = "Curve", values_to = "Number")

p_impact <- ggplot(impact_long, aes(x = Threshold, y = Number, color = Curve, linetype = Curve)) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = c("High-risk patients" = "#3C5488FF", "True positives" = "#E64B35FF")) +
  scale_linetype_manual(values = c("High-risk patients" = "solid", "True positives" = "dashed")) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Validation Clinical Impact Curve", subtitle = final_model, x = "Threshold Probability", y = "Number of patients", color = NULL, linetype = NULL) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5), axis.line = element_line(color = "black"), axis.text = element_text(color = "black"), axis.title = element_text(color = "black", face = "bold"), legend.position = "bottom")

ggsave(file.path(final_validation_dir, "Validation_Clinical_Impact_Curve_Clinical_Intra_Peri.pdf"), p_impact, width = 7.6, height = 6.2)
ggsave(file.path(final_validation_dir, "Validation_Clinical_Impact_Curve_Clinical_Intra_Peri.png"), p_impact, width = 7.6, height = 6.2, dpi = 300)


# ============================================================
# 15. Nomogram
# ============================================================

nomogram_data <- train_df[, c(label_col, final_predictors), drop = FALSE]
dd <- rms::datadist(nomogram_data)
options(datadist = "dd")

nom_fit <- rms::lrm(safe_formula(label_col, final_predictors), data = nomogram_data, x = TRUE, y = TRUE)
nom <- rms::nomogram(nom_fit, fun = plogis, funlabel = "Risk of fast growth", lp = FALSE)

save_base_pdf_png(
  draw_fun = function() plot(nom, xfrac = 0.34, cex.axis = 0.75, cex.var = 0.85),
  pdf_file = file.path(final_validation_dir, "Nomogram_Clinical_Intra_Peri.pdf"),
  png_file = file.path(final_validation_dir, "Nomogram_Clinical_Intra_Peri.png"),
  pdf_width = 10,
  pdf_height = 7,
  png_width = 2600,
  png_height = 1800,
  png_res = 300
)

nomogram_coef_df <- data.frame(Feature = names(coef(nom_fit)), Coefficient = as.numeric(coef(nom_fit)), row.names = NULL)


# ============================================================
# 16. Export final results and RData objects
# ============================================================

final_score_data <- model_df %>%
  dplyr::select(all_of(c(id_col, label_col, cohort_col)), Clinical_Score, Intra_Radscore, Peri_Radscore, all_of(final_pred_col))

wb_final <- openxlsx::createWorkbook()
add_sheet(wb_final, "final_metrics", final_metrics)
add_sheet(wb_final, "calibration_index", cal_index_df)
add_sheet(wb_final, "calibration_curve_all", cal_all_df)
add_sheet(wb_final, "calibration_plot_data", cal_plot_df)
add_sheet(wb_final, "validation_cal_curve", cal_valid_plot_df)
add_sheet(wb_final, "validation_cal_metrics", valid_cal_metric)
add_sheet(wb_final, "DCA_all", dca_all)
add_sheet(wb_final, "DCA_long", dca_long)
add_sheet(wb_final, "clinical_impact_valid", impact_valid)
add_sheet(wb_final, "final_score_prediction", final_score_data)
add_sheet(wb_final, "nomogram_coefficients", nomogram_coef_df)

openxlsx::saveWorkbook(
  wb_final,
  file.path(final_validation_dir, "Clinical_Intra_Peri_Final_Validation_Results.xlsx"),
  overwrite = TRUE
)

save(
  clinical_raw,
  intra_raw,
  peri_raw,
  split_df,
  clinical_obj,
  intra_obj,
  peri_obj,
  model_df,
  train_df,
  test_df,
  valid_df,
  model_fits,
  model_list,
  metrics_df,
  bootstrap_summary_df,
  bootstrap_detail_df,
  auc_fusion_df,
  train_delong,
  test_delong,
  final_model,
  final_predictors,
  final_pred_col,
  final_threshold,
  final_metrics,
  cal_all_df,
  cal_index_df,
  cal_plot_df,
  cal_valid_plot_df,
  valid_cal_metric,
  dca_all,
  dca_long,
  impact_valid,
  nom_fit,
  nom,
  nomogram_coef_df,
  file = file.path(output_dir, "Seed1111_Radiomics_Clinical_Intra_Peri_All_Objects.RData")
)


# ============================================================
# 17. Console summary
# ============================================================

cat("\n================ Seed1111 radiomics pipeline completed ================\n\n")
cat("Output directory:\n")
cat(normalizePath(output_dir), "\n\n")

cat("Random seed: ", seed_split, "\n")
cat("Pearson pre-filtering cutoff: |r with Label| >= ", pearson_r_cutoff, "\n")
cat("High-correlation cutoff: |r| >= ", high_corr_cutoff, "\n")
cat("VIF cutoff: ", vif_cutoff, "\n\n")

cat("Sample split:\n")
print(split_table)

cat("\nLASSO lambda.1se selected variables:\n")
cat("Clinical: ", nrow(clinical_obj$lasso$selected_1se), "\n")
cat("Intra: ", nrow(intra_obj$lasso$selected_1se), "\n")
cat("Peri: ", nrow(peri_obj$lasso$selected_1se), "\n\n")

cat("Seven-model train/test performance:\n")
print(metrics_df)

cat("\nBootstrap internal validation:\n")
print(bootstrap_summary_df)

cat("\nFinal Clinical_Intra_Peri train/test/validation performance:\n")
print(final_metrics)

cat("\nCalibration indices:\n")
print(cal_index_df)

cat("\nMain output files:\n")
cat(file.path(output_dir, "Seed1111_SevenModels_All_Results.xlsx"), "\n")
cat(file.path(final_validation_dir, "Clinical_Intra_Peri_Final_Validation_Results.xlsx"), "\n")
cat(file.path(output_dir, "Seed1111_Radiomics_Clinical_Intra_Peri_All_Objects.RData"), "\n")
cat(file.path(score_dir, "Train_ROC_7models.pdf"), "\n")
cat(file.path(score_dir, "Test_ROC_7models.pdf"), "\n")
cat(file.path(score_dir, "Train_BootstrapCorrected_Test_AUC_Fusion_Barplot.pdf"), "\n")
cat(file.path(bootstrap_dir, "Bootstrap_Internal_Validation_Visualization.pdf"), "\n")
cat(file.path(final_validation_dir, "Validation_ROC_Clinical_Intra_Peri.pdf"), "\n")
cat(file.path(final_validation_dir, "Train_Test_Validation_ROC_Clinical_Intra_Peri.pdf"), "\n")
cat(file.path(calibration_plot_dir, "Calibration_Plot_Clinical_Intra_Peri_Ideal_Apparent_Test_BiasCorrected.pdf"), "\n")
cat(file.path(validation_calibration_dir, "Validation_Calibration_Plot_Clinical_Intra_Peri_Logistic_Nonparametric.pdf"), "\n")
cat(file.path(final_validation_dir, "DCA_Train_Test_Validation_Clinical_Intra_Peri_y_minus0.1_to_0.5.pdf"), "\n")
cat(file.path(final_validation_dir, "Validation_Clinical_Impact_Curve_Clinical_Intra_Peri.pdf"), "\n")
cat(file.path(final_validation_dir, "Nomogram_Clinical_Intra_Peri.pdf"), "\n")

cat("\nWorkflow notes:\n")
cat("1. Seed 1111 was used for stratified train/test splitting among non-validation samples.\n")
cat("2. Original Cohort=2 was retained as the prospective temporal validation cohort.\n")
cat("3. Clinical, intratumor, and peritumor data were processed separately.\n")
cat("4. Imputation, Pearson filtering, high-correlation filtering, VIF filtering, and Z-score scaling were estimated only in the training cohort.\n")
cat("5. LASSO was performed only in the standardized training data.\n")
cat("6. The final Clinical_Intra_Peri model was evaluated in training, internal test, and prospective temporal validation cohorts.\n")

cat("\n================ Done ================\n")
