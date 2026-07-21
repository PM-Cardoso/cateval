#' Perform Pairwise Heterogeneous Treatment Effect Calibration Across Multiple Drugs
#'
#' This function computes heterogeneous treatment effect calibration curves for all pairwise comparisons
#' between a specified set of drugs. For each pair, it calculates the predicted benefit (difference between
#' predicted outcomes), partitions patients into calibration groups based on predicted benefit quantiles,
#' and estimates treatment effects within each group.
#'
#' Internally, the function calls `pairwise_calibration_pair()` for each drug pair to obtain
#' subgroup treatment effect estimates and confidence intervals.
#'
#' @param data A data.frame containing observed outcomes, treatment assignments, and predicted outcomes for each drug.
#' @param drug_var Character string specifying the column name of the treatment variable (e.g., assigned drug).
#' @param drugs Character vector of two or more drug names to be compared pairwise.
#' @param prediction_vars Character vector of column names corresponding to predicted outcomes for each drug.
#'   Must be the same length and order as `drugs`.
#' @param outcome_var Character string specifying the column name of the observed outcome variable.
#' @param cal_groups Numeric scalar or vector specifying the number(s) of calibration groups (quantiles) to partition
#'   the population for subgroup analysis.
#' @param matching Logical; if `TRUE`, performs covariate matching using the `MatchIt` package before estimating effects.
#' @param adjustment_var Optional character vector of covariate column names to adjust for in regression models.
#' @param matching_var Optional character vector specifying variables to match on. Defaults to `adjustment_var` if NULL.
#' @param match.exact Optional character vector of variables for exact matching.
#' @param match.antiexact Optional character vector of variables for anti-exact matching.
#' @param caliper Optional named numeric vector passed to \code{MatchIt::matchit}
#'   (e.g. \code{c(prehba1c = 0.4)}); interpreted in standard-deviation units by
#'   default (see \code{std.caliper}).
#' @param std.caliper Logical. If \code{TRUE} (default), \code{caliper} is read
#'   in standard-deviation units; if \code{FALSE}, in the variable's own units.
#' @param seed Optional integer. Seed set before matching each drug pair, so
#'   that the matched cohorts can be reproduced exactly.
#' @param method Character. Matching method passed to \code{MatchIt::matchit}
#'   (default \code{"nearest"}).
#' @param distance Character. Distance measure passed to \code{MatchIt::matchit}
#'   (default \code{"mahalanobis"}).
#' @param replace Logical. Whether a control patient may be reused in more than
#'   one matched set (default \code{TRUE}).
#' @param ratio Numeric. Number of controls matched to each treated patient
#'   (default \code{1}).
#' @param match_args Optional named list of further arguments for
#'   \code{MatchIt::matchit}. Applied last, so it overrides the arguments above
#'   and gives access to any matching option not exposed here.
#' @param return_matched Logical. If \code{TRUE}, return a list that also
#'   contains the matched cohort and \code{matchit} object for each drug pair,
#'   so the matching can be inspected. Default \code{FALSE}.
#' @param bootstrap Logical. If \code{FALSE} (default), confidence intervals are
#'   the model's Wald intervals. If \code{TRUE}, they are percentile bootstrap
#'   intervals. See \code{\link{pairwise_calibration_pair}} for details.
#' @param n_boot Integer. Number of bootstrap replicates (used when
#'   \code{bootstrap = TRUE}).
#' @param boot_seed Optional integer. Seed set before bootstrapping for
#'   reproducibility.
#' @param conf Numeric. Confidence level for the intervals (default 0.95).
#'
#' @return By default, a data.frame with one row per calibration group and drug pair combination, including:
#' \describe{
#'   \item{cal_groups}{Number of calibration groups used.}
#'   \item{grouping}{Calibration group identifier (e.g., 1 through `cal_groups`).}
#'   \item{mean}{Mean predicted benefit in the calibration group.}
#'   \item{coef}{Estimated treatment effect (regression coefficient) within the group.}
#'   \item{coef_low}{Lower bound of the 95% confidence interval for the estimate.}
#'   \item{coef_high}{Upper bound of the 95% confidence interval.}
#'   \item{drug1}{Name of the first drug in the pair.}
#'   \item{n_drug1}{Number of patients in the group receiving the first drug.}
#'   \item{drug2}{Name of the second drug in the pair.}
#'   \item{n_drug2}{Number of patients in the group receiving the second drug.}
#' }
#'
#' If \code{return_matched = TRUE}, a list with elements \code{results} (the
#' data.frame above), \code{matched} and \code{matchit}. The last two are named
#' lists with one entry per drug pair (e.g. \code{"SGLT2_vs_GLP1"}), and are
#' \code{NULL} when \code{matching = FALSE}.
#'
#' @seealso \code{\link{pairwise_calibration_pair}} for the single-pair engine
#'   and the details of the bootstrap.
#'
#'
#' @examples
#' \dontrun{
#' # Example usage comparing 3 drugs with predicted outcomes
#' results <- pairwise_calibration(
#'   data = mydata,
#'   drug_var = "treatment",
#'   drugs = c("DrugA", "DrugB", "DrugC"),
#'   prediction_vars = c("pred_drugA", "pred_drugB", "pred_drugC"),
#'   outcome_var = "outcome",
#'   cal_groups = c(3, 5),        # test 3 and 5 calibration groups
#'   adjustment_var = c("age", "sex"),
#'   matching = TRUE
#' )
#' head(results)
#' }
#'
#' @import dplyr
#' @importFrom utils combn
#' 
#' @export
pairwise_calibration <- function(data, 
                               drug_var, 
                               drugs, 
                               prediction_vars, 
                               outcome_var, 
                               cal_groups,
                               matching = FALSE,
                               adjustment_var = NULL,
                               matching_var = adjustment_var,
                               match.exact = NULL,
                               match.antiexact = NULL,
                               caliper = NULL,
                               std.caliper = TRUE,
                               seed = NULL,
                               method = "nearest",
                               distance = "mahalanobis",
                               replace = TRUE,
                               ratio = 1,
                               match_args = list(),
                               return_matched = FALSE,
                               bootstrap = FALSE,
                               n_boot = 1000,
                               boot_seed = NULL,
                               conf = 0.95) {
  
  # Define pipe
  `%>%` <- dplyr::`%>%`
  
  # Input validation
  if (!(drug_var %in% colnames(data))) stop("`drug_var` not found in data.")
  if (!all(prediction_vars %in% colnames(data))) stop("Some `prediction_vars` not found in data.")
  if (!(outcome_var %in% colnames(data))) stop("`outcome_var` not found in data.")
  if (!is.null(adjustment_var) && !all(adjustment_var %in% colnames(data))) stop("Some `adjustment_var` columns not found in data.")
  if (isTRUE(matching)) {
    if (length(matching_var) == 0) stop("Provide at least one matching_var when matching=TRUE.")
    if (!is.null(matching_var) && !all(matching_var %in% colnames(data))) stop("Some matching_var columns not found in data.")
    if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data.")
    if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data.")
    if (!is.null(caliper) && !is.numeric(caliper)) stop("`caliper` must be a named numeric vector.")
    if (!is.logical(std.caliper)) stop("`std.caliper` must be TRUE or FALSE.")
    if (!is.character(method) || length(method) != 1) stop("`method` must be a single character string.")
    if (!is.character(distance) || length(distance) != 1) stop("`distance` must be a single character string.")
    if (!is.logical(replace)) stop("`replace` must be TRUE or FALSE.")
    if (!is.numeric(ratio) || length(ratio) != 1 || ratio < 1) stop("`ratio` must be a single number of at least 1.")
    if (!is.list(match_args)) stop("`match_args` must be a list.")
  }
  if (!is.logical(return_matched)) stop("`return_matched` must be TRUE or FALSE.")
  if (!is.logical(bootstrap)) stop("`bootstrap` must be TRUE or FALSE.")
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1) stop("`n_boot` must be a single positive number.")
  if (!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) stop("`conf` must be a single number between 0 and 1.")
  if (length(drugs) < 2) stop("At least two drugs must be specified.")
  if (!all(drugs %in% unique(data[[drug_var]]))) stop("Some `drugs` not present in the `drug_var` column.")
  if (length(drugs) != length(prediction_vars)) stop("`drugs` and `prediction_vars` must have the same length.")
  if (!is.numeric(cal_groups)) stop("`cal_groups` must be numeric.")
  if (!is.logical(matching)) stop("`matching` must be TRUE or FALSE.")
  
  # Prepare dataset
  calibration_data <- data %>%
    dplyr::rename(
      dataset_drug_var = !!drug_var,
      dataset_outcome_var = !!outcome_var
    ) %>%
    dplyr::filter(dataset_drug_var %in% drugs)
  
  # Map prediction variables
  prediction_map <- stats::setNames(prediction_vars, drugs)
  output_table <- NULL

  # Matched cohorts and matchit objects, collected per drug pair when requested.
  matched_list <- list()
  matchit_list <- list()

  # Loop over unique drug pairs
  drug_combinations <- utils::combn(drugs, 2, simplify = FALSE)

  for (pair in drug_combinations) {
    
    pair_data <- calibration_data %>%
      dplyr::filter(dataset_drug_var %in% pair) %>%
      dplyr::mutate(
        dataset_drug_var = factor(dataset_drug_var, levels = rev(pair)),
        drug_1_pred = .[[prediction_map[[pair[1]]]]],
        drug_2_pred = .[[prediction_map[[pair[2]]]]],
        benefit = drug_1_pred - drug_2_pred
      )
    
    calibration_result <- pairwise_calibration_pair(
      data = pair_data,
      drug_var = "dataset_drug_var",
      drugs = pair,
      benefit_var = "benefit",
      outcome_var = "dataset_outcome_var",
      cal_groups = cal_groups,
      matching = matching,
      adjustment_var = adjustment_var,
      matching_var = matching_var,
      match.exact = match.exact,
      match.antiexact = match.antiexact,
      caliper = caliper,
      std.caliper = std.caliper,
      seed = seed,
      method = method,
      distance = distance,
      replace = replace,
      ratio = ratio,
      match_args = match_args,
      return_matched = isTRUE(return_matched) && isTRUE(matching),
      bootstrap = bootstrap,
      n_boot = n_boot,
      boot_seed = boot_seed,
      conf = conf
    )

    if (isTRUE(return_matched) && isTRUE(matching)) {
      pair_label <- paste(pair, collapse = "_vs_")
      matched_list[[pair_label]] <- calibration_result$matched
      matchit_list[[pair_label]] <- calibration_result$matchit
      calibration_result <- calibration_result$results
    }

    output_table <- dplyr::bind_rows(output_table, calibration_result)
  }

  if (isTRUE(return_matched)) {
    if (!isTRUE(matching)) {
      warning("`return_matched = TRUE` but `matching = FALSE`; no matched cohorts to return.")
      return(list(results = output_table, matched = NULL, matchit = NULL))
    }
    return(list(results = output_table, matched = matched_list, matchit = matchit_list))
  }

  return(output_table)
}