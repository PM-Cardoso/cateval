#' Estimate Observed vs. Predicted Treatment Benefit by Calibration Group
#'
#' This function evaluates how well predicted treatment benefits correspond to observed differences in outcomes.
#' Patients are divided into calibration groups based on predicted benefit, and observed outcome differences
#' between concordant and discordant treatments are regressed on the predicted benefit difference.
#' Concordance is defined either by receiving the top predicted treatment or any within a specified tolerance of the best.
#' Matching is used to adjust for confounding.
#'
#' @param data A data frame containing treatment assignments, outcomes, covariates, and predicted benefits.
#' @param drug_var Character. Name of the column indicating the actual treatment received.
#' @param outcome_var Character. Name of the outcome variable.
#' @param cal_groups Integer. Number of calibration groups to divide the population into (e.g., 5 for quintiles).
#' @param pred_cols Character vector. Column names of predicted outcomes (or risks) for each treatment. Should share a common prefix (e.g., "pred_GLP1", "pred_SGLT2").
#' @param conc_tolerance Optional numeric. If supplied, defines the absolute tolerance for determining concordance with the best predicted treatment. If NULL, concordance is based on the top predicted treatment only.
#' @param matching_var Character vector. Covariate names used for Mahalanobis distance matching.
#' @param match.exact Optional character vector. Variables required to match exactly across groups. The best predicted treatment is always included.
#' @param match.antiexact Optional character vector. Variables for anti-exact matching. The actual treatment variable is always included.
#' @param extract.match.cohort Optional logical. If TRUE, returns matched population used for calculating overall benefits.
#'
#' @return A data frame with one row per calibration group, containing:
#' \describe{
#'   \item{mean}{Average predicted benefit in the group (discordant minus concordant).}
#'   \item{coef}{Estimated observed difference in outcome (regression coefficient).}
#'   \item{coef_low}{Lower bound of the 95% confidence interval.}
#'   \item{coef_high}{Upper bound of the 95% confidence interval.}
#'   \item{n_groups}{Total number of calibration groups.}
#'   \item{n}{Number of matched pairs in each group.}
#' }
#'
#' @details
#' The function uses MatchIt for nearest-neighbor matching to control confounding. Concordant patients are those who received their predicted-optimal treatment, while discordant patients received a suboptimal one. Within each calibration group, observed outcomes are compared using regression.
#'
#' @importFrom dplyr mutate filter group_by ungroup distinct rename_with rowwise ungroup select ntile bind_rows across
#' @importFrom tidyr drop_na unnest_wider
#' @importFrom stringr str_detect str_split
#' @importFrom purrr map set_names
#' @importFrom MatchIt matchit get_matches
#' @importFrom stats lm predict
#' @export
compute_overall_benefit <- function(data, 
                                    drug_var, 
                                    outcome_var, 
                                    cal_groups, 
                                    pred_cols = NULL, 
                                    conc_tolerance = NULL,
                                    matching_var = NULL, 
                                    match.exact = NULL, 
                                    match.antiexact = NULL,
                                    extract.match.cohort = FALSE) {
  
  `%>%` <- dplyr::`%>%`
  
  # Input checks
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  if (!all(pred_cols %in% colnames(data))) stop("Some pred_cols not found in data")
  if (!is.null(conc_tolerance) && !is.numeric(conc_tolerance)) stop("conc_tolerance must be numeric")
  if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
  if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
  if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  if (!is.numeric(cal_groups)) stop("cal_groups must be numeric")
  if (!is.logical(extract.match.cohort)) stop("extract.match.cohort must be logical")
  
  # Rename input vars
  pre_data <- data %>%
    dplyr::rename_with(~ "dataset_drug_var", dplyr::all_of(drug_var)) %>%
    dplyr::rename_with(~ "dataset_outcome_var", dplyr::all_of(outcome_var))
  
  # Extract prediction prefix
  common_prefix <- function(strings) {
    min_len <- min(nchar(strings))
    prefix <- ""
    for (i in seq_len(min_len)) {
      current_chars <- substr(strings, i, i)
      if (length(unique(current_chars)) == 1) {
        prefix <- paste0(prefix, current_chars[1])
      } else break
    }
    prefix
  }
  prefix <- common_prefix(pred_cols)
  drug_names <- gsub(prefix, "", pred_cols)
  
  # Get rank-1 predicted drug
  interim_dataset <- get_best_drugs(
    data = pre_data,
    rank = 1,
    column_names = pred_cols,
    final_var_name = prefix
  ) %>%
    dplyr::rename("function_rank1_drug_name" := paste0(prefix, "rank1_drug_name"))
  
  # Define concordance
  if (is.null(conc_tolerance)) {
    interim_dataset <- interim_dataset %>%
      dplyr::mutate(conc_disc_label = ifelse(dataset_drug_var == function_rank1_drug_name, 1, 0))
  } else {
    n_drugs_required <- length(pred_cols)
    tolerance_vars <- paste0("tolerance_drug_", 1:n_drugs_required)
    
    interim_dataset <- get_best_drugs(
      data = interim_dataset,
      tolerance = conc_tolerance,
      column_names = pred_cols,
      final_var_name = prefix
    ) %>%
      dplyr::rename("function_tolerance_drug_name" := paste0(prefix, "within_", conc_tolerance, "_of_best_drug_name")) %>%
      dplyr::mutate(
        conc_disc_label = ifelse(
          stringr::str_detect(function_tolerance_drug_name, paste0("\\b", dataset_drug_var, "\\b")),
          1, 0
        ),
        drug_list = stringr::str_split(function_tolerance_drug_name, "\\s*[;,\\s]\\s*")
      ) %>%
      dplyr::mutate(drug_list = purrr::map(drug_list, ~ rep(.x, length.out = n_drugs_required))) %>%
      dplyr::mutate(drug_list = purrr::map(drug_list, ~ purrr::set_names(.x, paste0("tolerance_drug_", seq_along(.x))))) %>%
      tidyr::unnest_wider(drug_list) %>%
      dplyr::mutate(dplyr::across(
        dplyr::all_of(tolerance_vars),
        ~ dplyr::if_else(conc_disc_label == 0, dataset_drug_var, .x)
      ))
  }
  
  # Drop NA for matching
  interim_dataset <- tidyr::drop_na(interim_dataset, dplyr::all_of(matching_var))
  
  # Matching formula
  categorical_vars <- matching_var[sapply(interim_dataset[matching_var], function(x) is.factor(x) || is.character(x))]
  cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))
  
  matching_formula <- paste("conc_disc_label ~", paste(cont_vars, collapse = " + "))
  for (v in categorical_vars) {
    if (length(unique(interim_dataset[[v]])) > 1) {
      matching_formula <- paste(matching_formula, "+", v)
    }
  }
  
  match_model_antiexact_vars <- if (!is.null(conc_tolerance)) {
    unique(c(tolerance_vars, "dataset_drug_var", match.antiexact))
  } else {
    unique(c("dataset_drug_var", match.antiexact))
  }
  
  # Matching
  match_model <- MatchIt::matchit(
    formula = as.formula(matching_formula),
    data = interim_dataset,
    method = "nearest",
    distance = "mahalanobis",
    replace = FALSE,
    exact = unique(c("function_rank1_drug_name", match.exact)),
    antiexact = match_model_antiexact_vars
  )
  
  matched_data <- MatchIt::get_matches(match_model, data = interim_dataset)
  
  # Process matched pairs
  processed_data <- matched_data %>%
    dplyr::group_by(subclass) %>%
    dplyr::mutate(
      concordant_drugclass = dataset_drug_var[1],
      discordant_drugclass = dataset_drug_var[2],
      calibration_obs = diff(dataset_outcome_var)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(subclass, .keep_all = TRUE) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      calibration_pred = get(paste0(prefix, discordant_drugclass)) - get(paste0(prefix, concordant_drugclass))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(grouping = dplyr::ntile(calibration_pred, cal_groups)) %>%
    dplyr::select(calibration_pred, calibration_obs, grouping)
  
  # Initialize results
  coef      <- rep(NA_real_, cal_groups)
  coef_low  <- rep(NA_real_, cal_groups)
  coef_high <- rep(NA_real_, cal_groups)
  mean_vals <- rep(NA_real_, cal_groups)
  n_group   <- rep(0, cal_groups)
  
  # Fit model in each group
  for (cg in seq_len(cal_groups)) {
    group_data <- processed_data %>% dplyr::filter(grouping == cg)
    
    mean_vals[cg] <- mean(group_data$calibration_pred, na.rm = TRUE)
    n_group[cg]   <- nrow(group_data)
    
    model <- stats::lm(calibration_obs ~ calibration_pred, data = group_data)
    
    predicted_value <- stats::predict(
      model,
      newdata = as.data.frame(t(colMeans(group_data, na.rm = TRUE))),
      interval = "confidence"
    ) %>% unlist()
    
    coef[cg]      <- predicted_value[1]
    coef_low[cg]  <- predicted_value[2]
    coef_high[cg] <- predicted_value[3]
  }
  
  # Return summary
  return(data.frame(
    mean     = mean_vals,
    coef     = coef,
    coef_low = coef_low,
    coef_high= coef_high,
    n_groups = cal_groups,
    n        = n_group
  ))
}
