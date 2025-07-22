#' Estimate Observed vs Predicted Treatment Benefit Calibration
#'
#' This function evaluates how well predicted treatment benefits align with observed outcome differences. It groups patients by predicted benefit and compares outcomes between patients who received treatments concordant or discordant with the predicted optimal treatment. Concordance can be defined strictly by the top predicted treatment or flexibly within a specified tolerance of the best predicted outcome.
#'
#' To control for confounding, the function performs nearest-neighbor matching based on specified covariates before estimating observed treatment effect differences within matched pairs. The observed differences are then regressed against predicted differences to assess calibration of predicted treatment benefits.
#'
#' @param data A data frame containing individual-level data including treatment assignments, outcomes, and predicted treatment benefits.
#' @param drug_var Character. Column name indicating the treatment actually received by each patient.
#' @param outcome_var Character. Column name for the outcome variable to assess treatment effect.
#' @param cal_groups Integer. Number of calibration groups to stratify patients based on predicted benefit (currently not directly used but reserved for future enhancements).
#' @param pred_cols Character vector. Column names of predicted outcomes or risks for each treatment. These should share a common prefix (e.g., "pred_GLP1", "pred_SGLT2").
#' @param conc_tolerance Optional numeric scalar. If provided, patients are considered concordant if they received any treatment within this absolute difference of the best predicted outcome. If NULL (default), concordance requires receiving the top predicted treatment exactly.
#' @param matching_var Character vector. Names of covariates used for Mahalanobis distance matching between concordant and discordant patients.
#' @param match.exact Optional character vector. Variables for exact matching. The best predicted treatment is always included automatically.
#' @param match.antiexact Optional character vector. Variables for anti-exact matching (i.e., variables on which matches must differ). The actual treatment variable is always included automatically.
#'
#' @return 
#' A named list with two elements:
#' \describe{
#'   \item{calibration_intercept}{A one-row data frame with the intercept of the linear regression model (`calibration_obs ~ calibration_pred`), representing the average deviation from perfect calibration when the predicted difference is zero. Includes:}
#'     \describe{
#'       \item{value}{Estimated intercept.}
#'       \item{lci}{Lower 95% confidence interval for the intercept.}
#'       \item{uci}{Upper 95% confidence interval for the intercept.}
#'     }
#'
#'   \item{calibration_slope}{A one-row data frame with the slope of the regression model, representing how well predicted differences in treatment effect correspond to observed differences. Includes:}
#'     \describe{
#'       \item{value}{Estimated slope. A value near 1 indicates good calibration.}
#'       \item{lci}{Lower 95% confidence interval for the slope.}
#'       \item{uci}{Upper 95% confidence interval for the slope.}
#'     }
#' }
#'
#' Each estimate reflects model fit across matched patient pairs grouped by predicted benefit.
#'
#' @details
#' This function uses nearest-neighbor matching (via the MatchIt package) on a set of covariates to create matched groups of concordant and discordant patients. Within each calibration group defined by predicted benefit, the observed difference in outcome is regressed on the predicted difference to assess calibration.
#'
#' @importFrom dplyr mutate filter group_by ungroup distinct rowwise select rename_with across
#' @importFrom tidyr drop_na unnest_wider
#' @importFrom purrr map set_names
#' @importFrom stringr str_detect str_split
#' @importFrom MatchIt matchit get_matches
#' @importFrom stats coef confint lm
#' @export
#'
#' @examples
#' \dontrun{
#' result <- compute_overall_benefit_performance(
#'   data = data_example,
#'   drug_var = "treatment",
#'   outcome_var = "outcome",
#'   pred_cols = c("pred_GLP1", "pred_SGLT2", "pred_DPP4"),
#'   conc_tolerance = 0.05,
#'   matching_var = c("age", "sex"),
#'   match.exact = NULL,
#'   match.antiexact = NULL
#' )
#' print(result)
#' }
compute_overall_benefit_performance <- function(data, 
                                                  drug_var, 
                                                  outcome_var, 
                                                  pred_cols = NULL, 
                                                  conc_tolerance = NULL,
                                                  matching_var = NULL, 
                                                  match.exact = NULL, 
                                                  match.antiexact = NULL) {
  
  `%>%` <- dplyr::`%>%`
  
  # Input validation
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  if (!all(pred_cols %in% colnames(data))) stop("Some pred_cols not found in data")
  if (!is.null(conc_tolerance) && !is.numeric(conc_tolerance)) stop("conc_tolerance must be numeric")
  if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
  if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
  if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  
  # Rename for consistency
  pre_data <- data %>%
    dplyr::rename_with(~ "dataset_drug_var", dplyr::all_of(drug_var)) %>%
    dplyr::rename_with(~ "dataset_outcome_var", dplyr::all_of(outcome_var))
  
  # Helper to get common prefix
  common_prefix <- function(strings) {
    min_len <- min(nchar(strings))
    prefix <- ""
    for (i in seq_len(min_len)) {
      chars <- substr(strings, i, i)
      if (length(unique(chars)) == 1) {
        prefix <- paste0(prefix, chars[1])
      } else break
    }
    prefix
  }
  
  prefix <- common_prefix(pred_cols)
  drug_names <- gsub(prefix, "", pred_cols)
  
  # Get best predicted drugs
  interim_dataset <- get_best_drugs(
    data = pre_data,
    rank = 1,
    column_names = pred_cols,
    final_var_name = prefix
  ) %>%
    dplyr::rename("function_rank1_drug_name" := paste0(prefix, "rank1_drug_name"))
  
  # Determine concordance
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
        )
      ) %>%
      dplyr::mutate(
        drug_list = stringr::str_split(function_tolerance_drug_name, "\\s*[;,\\s]\\s*")
      ) %>%
      dplyr::mutate(
        drug_list = purrr::map(drug_list, ~ rep(.x, length.out = n_drugs_required))
      ) %>%
      dplyr::mutate(
        drug_list = purrr::map(drug_list, ~ purrr::set_names(.x, tolerance_vars))
      ) %>%
      tidyr::unnest_wider(drug_list) %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(tolerance_vars),
          ~ dplyr::if_else(conc_disc_label == 0, dataset_drug_var, .x)
        )
      )
  }
  
  # Drop NA in matching vars
  interim_dataset <- interim_dataset %>% tidyr::drop_na(dplyr::all_of(matching_var))
  
  categorical_vars <- matching_var[sapply(interim_dataset[matching_var], function(x) is.factor(x) || is.character(x))]
  cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))
  
  # Build matching formula
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
    dplyr::select(calibration_pred, calibration_obs)
  
  linear_model <- stats::lm(calibration_obs ~ calibration_pred, data = processed_data)
  
  result <- list(
    calibration_intercept = data.frame(
      row.names = "Intercept",
      value = stats::coef(linear_model)[1],
      lci = stats::confint(linear_model)[1, 1],
      uci = stats::confint(linear_model)[1, 2]
    ),
    calibration_slope = data.frame(
      row.names = "Slope",
      value = stats::coef(linear_model)[2],
      lci = stats::confint(linear_model)[2, 1],
      uci = stats::confint(linear_model)[2, 2]
    )
  )
  
  return(result)
  
}
