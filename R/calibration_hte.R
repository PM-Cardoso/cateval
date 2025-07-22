#' Estimate Heterogeneous Treatment Effects Across Calibration Groups
#'
#' Estimates heterogeneous treatment effects between two drugs by stratifying patients into calibration groups based on predicted benefit scores. Optionally performs covariate matching before estimating treatment effects.
#'
#' If multiple values are provided to `cal_groups`, the function will repeat the process for each value and return a combined result.
#'
#' @param data A data frame containing observed outcomes, treatment assignments, and predicted benefits.
#' @param drug_var Character string. Column name for the treatment assignment variable (e.g., drug).
#' @param drugs Character vector of length 2. Names of the two drugs to compare.
#' @param benefit_var Character string. Column name containing predicted benefit scores.
#' @param outcome_var Character string. Column name for the outcome variable (e.g., clinical measurement).
#' @param cal_groups Numeric or numeric vector. Number(s) of calibration groups (e.g., quantiles) to divide the data into based on predicted benefit.
#' @param matching Logical. Whether to perform covariate matching using the `MatchIt` package before estimating treatment effects.
#' @param adjustment_var Optional character vector. Names of covariates to include as adjustment variables in the regression model.
#' @param matching_var Optional character vector. Covariates to use for matching. Defaults to `adjustment_var` if not specified.
#' @param match.exact Optional character vector. Variables for exact matching. Matching on best predicted drug automatically added.
#' @param match.antiexact Optional character vector. Variables for anti-exact matching. `drug_var` automatically added.
#'
#' @return 
#' A data frame where each row corresponds to a calibration group within each `cal_groups` setting. 
#' Columns include:
#' 
#' \describe{
#'   \item{mean}{The mean predicted benefit score for patients in the calibration group.}
#'   \item{coef}{Estimated average treatment effect (regression coefficient) comparing the two drugs, adjusted if covariates are specified.}
#'   \item{coef_low}{Lower bound of the 95% confidence interval for the treatment effect.}
#'   \item{coef_high}{Upper bound of the 95% confidence interval for the treatment effect.}
#'   \item{n_groups}{Number of calibration groups (i.e., the value of cal_groups) used to create this stratification.}
#'   \item{drug1}{Name of the first drug in the comparison (from drugs[1]).}
#'   \item{n_drug1}{Number of patients receiving drug1 within the calibration group.}
#'   \item{drug2}{Name of the second drug in the comparison (from drugs[2]).}
#'   \item{n_drug2}{Number of patients receiving drug2 within the calibration group.}
#' }
#' 
#' @import dplyr
#' @importFrom stats glm confint.default
#' @importFrom MatchIt matchit get_matches
#'
#' @examples
#' # Basic usage without matching or adjustment
#' result <- calibration_hte(
#'   data = test_data,
#'   drug_var = "drugclass",
#'   drugs = c("SGLT2", "DPP4"),
#'   benefit_var = "benefit_score",
#'   outcome_var = "posthba1cfinal",
#'   cal_groups = 5
#' )
#' 
#' # With adjustment
#' result_adj <- calibration_hte(
#'   data = test_data,
#'   drug_var = "drugclass",
#'   drugs = c("SGLT2", "DPP4"),
#'   benefit_var = "benefit_score",
#'   outcome_var = "posthba1cfinal",
#'   cal_groups = 3,
#'   adjustment_var = c("age", "bmi", "sex")
#' )
#'
#' # With matching and exact matching
#' result_match <- calibration_hte(
#'   data = test_data,
#'   drug_var = "drugclass",
#'   drugs = c("SGLT2", "DPP4"),
#'   benefit_var = "benefit_score",
#'   outcome_var = "posthba1cfinal",
#'   cal_groups = 4,
#'   matching = TRUE,
#'   adjustment_var = c("age", "bmi", "sex"),
#'   match.exact = "sex"
#' )
#'
#' @export
calibration_hte <- function(data,
                            drug_var,
                            drugs,
                            benefit_var,
                            outcome_var,
                            cal_groups,
                            matching = FALSE,
                            adjustment_var = NULL,
                            matching_var = adjustment_var,
                            match.exact = NULL, 
                            match.antiexact = NULL) {
  
  # Load Required Libraries (granularly)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("MatchIt", quietly = TRUE)) stop("Package 'MatchIt' is required.")
  
  # Use package-qualified functions throughout
  `%>%` <- dplyr::`%>%`
  
  # Validate Inputs ----
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  if (!(benefit_var %in% colnames(data))) stop("benefit_var not found in data")
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  if (!is.null(adjustment_var) && !all(adjustment_var %in% colnames(data))) stop("Some adjustment_var not in data")
  if (isTRUE(matching)) {
    if (length(matching_var) == 0) stop("Provide at least one matching_var")
    if (!all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
    if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
    if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
  }
  if (length(drugs) != 2) stop("Exactly two drugs must be specified")
  if (!all(drugs %in% unique(data[[drug_var]]))) stop("Some specified drugs not present in drug_var column")
  if (!is.numeric(cal_groups)) stop("cal_groups must be numeric")
  if (!is.logical(matching)) stop("matching must be TRUE or FALSE")
  
  result <- NULL
  
  for (cg in cal_groups) {
    if (cg <= 0) stop("Each value in cal_groups must be a positive number")
    
    initial_dataset <- data %>%
      dplyr::rename(
        dataset_benefit = !!benefit_var,
        dataset_drug_var = !!drug_var,
        dataset_outcome_var = !!outcome_var
      ) %>%
      dplyr::filter(dataset_drug_var %in% drugs)
    
    if (isTRUE(matching)) {
      initial_dataset <- initial_dataset %>%
        dplyr::mutate(conc_disc_label = ifelse(dataset_benefit <= 0, 1, 0)) %>%
        tidyr::drop_na(dplyr::all_of(matching_var))
      
      if (length(unique(initial_dataset$conc_disc_label)) < 2) next
      
      categorical_vars <- matching_var[sapply(initial_dataset[matching_var], function(x) is.factor(x) || is.character(x))]
      cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))
      
      matching_formula <- paste("conc_disc_label ~", paste(cont_vars, collapse = " + "))
      for (v in categorical_vars) {
        if (length(unique(initial_dataset[[v]])) > 1) {
          matching_formula <- paste(matching_formula, "+", v)
        }
      }
      
      match_model <- MatchIt::matchit(
        formula = as.formula(matching_formula),
        data = initial_dataset,
        method = "nearest",
        distance = "mahalanobis",
        replace = FALSE,
        exact = match.exact,
        antiexact = match.antiexact
      )
      
      calibration_data <- MatchIt::get_matches(match_model, data = initial_dataset)
      
    } else {
      calibration_data <- initial_dataset
    }
    
    coef      <- rep(NA_real_, cg)
    coef_low  <- rep(NA_real_, cg)
    coef_high <- rep(NA_real_, cg)
    mean_vals <- rep(NA_real_, cg)
    n_drug1   <- rep(0, cg)
    n_drug2   <- rep(0, cg)
    
    calibration_data <- calibration_data %>%
      dplyr::mutate(grouping = dplyr::ntile(dataset_benefit, cg))
    
    for (g in seq_len(cg)) {
      group_data <- calibration_data %>%
        dplyr::filter(grouping == g) %>%
        dplyr::mutate(dataset_drug_var = factor(dataset_drug_var, levels = rev(drugs)))
      
      mean_vals[g] <- mean(group_data$dataset_benefit, na.rm = TRUE)
      n_drug1[g] <- nrow(dplyr::filter(group_data, dataset_drug_var == drugs[1]))
      n_drug2[g] <- nrow(dplyr::filter(group_data, dataset_drug_var == drugs[2]))
      
      if (length(unique(group_data$dataset_drug_var)) < 2) next
      
      formula_str <- "dataset_outcome_var ~ dataset_drug_var"
      
      if (!is.null(adjustment_var)) {
        cat_vars  <- adjustment_var[sapply(group_data[, adjustment_var], function(x) is.factor(x) || is.character(x))]
        cont_vars <- setdiff(adjustment_var, cat_vars)
        
        if (length(cont_vars) > 0) {
          formula_str <- paste(formula_str, paste(cont_vars, collapse = " + "), sep = " + ")
        }
        for (v in cat_vars) {
          if (length(unique(group_data[[v]])) > 1) {
            formula_str <- paste(formula_str, "+", v)
          }
        }
      }
      
      model <- stats::glm(as.formula(formula_str), data = group_data)
      ci <- suppressMessages(stats::confint.default(model))
      
      coef[g]      <- stats::coef(model)[2]
      coef_low[g]  <- ci[2, 1]
      coef_high[g] <- ci[2, 2]
    }
    
    result <- dplyr::bind_rows(
      result,
      data.frame(
        mean      = mean_vals,
        coef      = coef,
        coef_low  = coef_low,
        coef_high = coef_high,
        n_groups  = cg,
        drug1     = drugs[1],
        n_drug1   = n_drug1,
        drug2     = drugs[2],
        n_drug2   = n_drug2
      )
    )
  }
  
  return(result)
}
