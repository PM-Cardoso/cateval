#' Estimate Heterogeneous Treatment Effects Across Calibration Groups
#'
#' Estimates the treatment effect between two drugs within groups of predicted
#' benefit. Patients are ordered by their predicted benefit, split into equal-sized
#' calibration groups, and the observed treatment effect is estimated within each
#' group. Comparing that observed effect with the predicted benefit is what shows
#' whether the model is well calibrated. Covariate matching can optionally be
#' applied before the effects are estimated.
#'
#' The predicted benefit is expected on the scale \code{drugs[1]} minus
#' \code{drugs[2]} (as produced by \code{\link{pairwise_calibration}}), so a
#' negative value means \code{drugs[1]} is predicted to give the lower (better)
#' outcome. \code{coef} is the observed difference in outcome on
#' \code{drugs[1]} relative to \code{drugs[2]}, on the same scale.
#'
#' If several values are supplied to `cal_groups`, the process is repeated for
#' each value and the results are combined. Matching, when requested, is carried
#' out once and the same matched cohort is reused for every value of
#' `cal_groups`, so the different groupings are directly comparable.
#'
#' @param data A data frame containing observed outcomes, treatment assignments, and predicted benefits.
#' @param drug_var Character string. Column name for the treatment assignment variable (e.g., drug).
#' @param drugs Character vector of length 2. Names of the two drugs to compare.
#' @param benefit_var Character string. Column name containing predicted benefit scores.
#' @param outcome_var Character string. Column name for the outcome variable (e.g., clinical measurement).
#' @param cal_groups Numeric or numeric vector. Number(s) of calibration groups (e.g., quantiles) to divide the data into based on predicted benefit.
#' @param matching Logical. Whether to perform covariate matching using the `MatchIt` package before estimating treatment effects.
#' @param adjustment_var Optional character vector. Names of covariates to include as adjustment variables in the regression model. Must not include `drug_var` or `outcome_var`.
#' @param matching_var Optional character vector. Covariates to use for matching. Defaults to `adjustment_var` if not specified. Must not include `drug_var` or `outcome_var`. Patients with a missing value in any `matching_var`, `match.exact` or `match.antiexact` variable are excluded before matching.
#' @param match.exact Optional character vector. Variables for exact matching. Must not include `drug_var`, since matched patients are always required to have received different drugs.
#' @param match.antiexact Optional character vector. Additional variables for anti-exact matching. The received drug is always included (see Details).
#' @param caliper Optional named numeric vector passed to \code{MatchIt::matchit}
#'   (e.g. \code{c(prehba1c = 0.4)}); interpreted in standard-deviation units by
#'   default (see \code{std.caliper}).
#' @param std.caliper Logical. If \code{TRUE} (default), \code{caliper} is read
#'   in standard-deviation units; if \code{FALSE}, in the variable's own units.
#' @param seed Optional integer. Seed set before matching, so that a matched
#'   cohort can be reproduced exactly.
#' @param method Character. Matching method passed to \code{MatchIt::matchit}
#'   (default \code{"nearest"}).
#' @param distance Character. Distance measure passed to \code{MatchIt::matchit}
#'   (default \code{"mahalanobis"}).
#' @param replace Logical. Whether a control patient may be reused in more than
#'   one matched set (default \code{TRUE}). See the note on confidence intervals
#'   below.
#' @param ratio Numeric. Number of controls matched to each treated patient
#'   (default \code{1}). With \code{ratio > 1} every matched patient enters the
#'   outcome model with equal weight (MatchIt's matching weights are not used),
#'   so a treated patient counts once alongside each of its controls.
#' @param match_args Optional named list of further arguments for
#'   \code{MatchIt::matchit}. Applied after the arguments above, so it overrides
#'   them and gives access to any matching option not exposed here. The
#'   anti-exact constraint on the received drug is re-applied afterwards and
#'   cannot be removed; an \code{antiexact} entry adds to it.
#' @param return_matched Logical. If \code{TRUE}, return a list that also
#'   contains the matched cohort and the \code{matchit} object, so the matching
#'   can be inspected (for example with a balance plot). Default \code{FALSE}.
#' @param bootstrap Logical. If \code{FALSE} (default), confidence intervals are
#'   the model's Wald intervals. If \code{TRUE}, they are percentile bootstrap
#'   intervals (see Details).
#' @param n_boot Integer. Number of bootstrap replicates (used when
#'   \code{bootstrap = TRUE}).
#' @param boot_seed Optional integer. Seed set before bootstrapping for
#'   reproducibility.
#' @param conf Numeric. Confidence level for the intervals (default 0.95).
#'
#' @return By default, a data frame where each row corresponds to a calibration group within each `cal_groups` setting. Columns include:
#'
#' \describe{
#'   \item{mean}{The mean predicted benefit score for patients in the calibration group.}
#'   \item{coef}{Estimated average treatment effect (regression coefficient) comparing the two drugs, adjusted if covariates are specified. \code{NA} when the group does not contain patients on both drugs.}
#'   \item{coef_low}{Lower bound of the confidence interval for the treatment effect. \code{NA} when the effect, or its interval, cannot be estimated (see Details).}
#'   \item{coef_high}{Upper bound of the confidence interval for the treatment effect.}
#'   \item{n_groups}{Number of calibration groups (i.e., the value of cal_groups) used to create this stratification.}
#'   \item{drug1}{Name of the first drug in the comparison (from `drugs[1]`).}
#'   \item{n_drug1}{Number of patients receiving drug1 within the calibration group.}
#'   \item{drug2}{Name of the second drug in the comparison (from `drugs[2]`).}
#'   \item{n_drug2}{Number of patients receiving drug2 within the calibration group.}
#' }
#'
#' If \code{return_matched = TRUE}, a list with elements \code{results} (the data
#' frame above), \code{matched} (the matched cohort) and \code{matchit} (the
#' \code{matchit} object). The last two are \code{NULL} when
#' \code{matching = FALSE}, since no matching was performed.
#'
#' @details
#' \strong{Matching.} When \code{matching = TRUE}, patients predicted to do
#' better on \code{drugs[1]} (\code{benefit_var <= 0}) form the treated group
#' and are matched to patients predicted to do better on \code{drugs[2]}. The
#' match is always anti-exact on the received drug, so every matched set
#' contains patients who received different drugs, the same requirement as in
#' \code{\link{match_benefit_pairs}}. Arguments are layered as package defaults,
#' then the named arguments (\code{method}, \code{distance}, \code{replace},
#' \code{ratio}, \code{caliper}, \code{std.caliper}, \code{match.exact},
#' \code{match.antiexact}), then \code{match_args}; the anti-exact constraint
#' on the received drug is re-applied last. Only pair-matching methods (e.g.
#' \code{"nearest"}, \code{"optimal"}, \code{"genetic"}) can be used. If no
#' matches can be formed (for example, every patient predicted to do better on
#' one drug received the same drug as every patient predicted to do better on
#' the other), \code{NULL} is returned with a warning.
#'
#' \strong{Groups too small for an interval.} A group's Wald interval needs at
#' least one residual degree of freedom in the outcome model (more patients than
#' model coefficients). When that is not met the point estimate is still
#' returned but \code{coef_low} and \code{coef_high} are \code{NA}, with a
#' warning naming the group.
#'
#' \strong{Bootstrap intervals.} The calibration groups are fixed from the
#' original cohort before bootstrapping, so each replicate re-estimates the
#' treatment effect for the *same* groups rather than shifting the group
#' boundaries between replicates. When matching has been applied the resampling
#' unit is the matched set, which keeps matched patients together; otherwise
#' patients are resampled individually. The matching itself is not repeated
#' within each replicate, so the intervals reflect sampling variability in the
#' outcome model rather than uncertainty in the matching. A replicate in which a
#' group happens to contain only one of the two drugs contributes no estimate and
#' is ignored.
#'
#' \strong{Confidence intervals under matching with replacement.} When
#' \code{matching = TRUE} and \code{replace = TRUE} (the default), a control
#' patient can appear in several matched sets, so the matched observations are
#' not independent. Neither the Wald nor the bootstrap interval widens to reflect
#' this, so both should be read as approximate and the limitation stated when
#' they are reported. Setting \code{replace = FALSE} removes the dependence, at
#' the cost of fewer and generally poorer-quality matches.
#'
#' @import dplyr
#' @importFrom stats glm confint.default
#' @importFrom MatchIt matchit get_matches
#'
#' @examples
#' \dontrun{
#' # Basic usage without matching or adjustment
#' result <- pairwise_calibration_pair(
#'   data = test_data,
#'   drug_var = "drugclass",
#'   drugs = c("SGLT2", "DPP4"),
#'   benefit_var = "benefit_score",
#'   outcome_var = "posthba1cfinal",
#'   cal_groups = 5
#' )
#'
#' # With adjustment
#' result_adj <- pairwise_calibration_pair(
#'   data = test_data,
#'   drug_var = "drugclass",
#'   drugs = c("SGLT2", "DPP4"),
#'   benefit_var = "benefit_score",
#'   outcome_var = "posthba1cfinal",
#'   cal_groups = 3,
#'   adjustment_var = c("age", "bmi", "sex")
#' )
#'
#' # With matching, a caliper and a seed, keeping the matched cohort
#' result_match <- pairwise_calibration_pair(
#'   data = test_data,
#'   drug_var = "drugclass",
#'   drugs = c("SGLT2", "DPP4"),
#'   benefit_var = "benefit_score",
#'   outcome_var = "posthba1cfinal",
#'   cal_groups = 4,
#'   matching = TRUE,
#'   adjustment_var = c("age", "bmi", "sex"),
#'   match.exact = "sex",
#'   caliper = c(age = 0.25),
#'   seed = 42,
#'   return_matched = TRUE
#' )
#' plot(summary(result_match$matchit))
#' }
#'
#' @export
pairwise_calibration_pair <- function(data,
                            drug_var,
                            drugs,
                            benefit_var,
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

  # Use package-qualified functions throughout
  `%>%` <- dplyr::`%>%`

  # Validate Inputs ----
  if (!(drug_var %in% colnames(data))) stop("drug_var not found in data")
  if (!(benefit_var %in% colnames(data))) stop("benefit_var not found in data")
  if (!(outcome_var %in% colnames(data))) stop("outcome_var not found in data")
  if (!is.null(adjustment_var) && !all(adjustment_var %in% colnames(data))) stop("Some adjustment_var not in data")
  .check_reserved_vars(
    list(adjustment_var = adjustment_var),
    c(drug_var = drug_var, outcome_var = outcome_var)
  )
  if (isTRUE(matching)) {
    if (length(matching_var) == 0) stop("Provide at least one matching_var")
    if (!all(matching_var %in% colnames(data))) stop("Some matching_var not in data")
    if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some match.exact variables not in data")
    if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some match.antiexact variables not in data")
    if (!is.null(caliper) && !is.numeric(caliper)) stop("`caliper` must be a named numeric vector.")
    if (!is.logical(std.caliper)) stop("`std.caliper` must be TRUE or FALSE.")
    if (!is.character(method) || length(method) != 1) stop("`method` must be a single character string.")
    if (!is.character(distance) || length(distance) != 1) stop("`distance` must be a single character string.")
    if (!is.logical(replace)) stop("`replace` must be TRUE or FALSE.")
    if (!is.numeric(ratio) || length(ratio) != 1 || ratio < 1) stop("`ratio` must be a single number of at least 1.")
    if (!is.list(match_args)) stop("`match_args` must be a list.")
    .check_named_list(match_args, "match_args")
    .check_reserved_vars(
      list(matching_var = matching_var, match.exact = match.exact),
      c(drug_var = drug_var, outcome_var = outcome_var)
    )
    .check_reserved_vars(list(match.antiexact = match.antiexact), c(outcome_var = outcome_var))
    # The received drug is anti-exact internally (on a copy of drug_var), so an
    # exact match on drug_var under its own name is the same contradiction.
    if (drug_var %in% .spec_vars(match_args$exact)) {
      stop("Variable(s) ", drug_var, " cannot be matched both exactly and anti-exactly.")
    }
  }
  if (length(drugs) != 2) stop("Exactly two drugs must be specified")
  if (!all(drugs %in% unique(data[[drug_var]]))) stop("Some specified drugs not present in drug_var column")
  if (!is.numeric(cal_groups)) stop("cal_groups must be numeric")
  if (any(cal_groups <= 0)) stop("Each value in cal_groups must be a positive number")
  if (!is.logical(matching)) stop("matching must be TRUE or FALSE")
  if (!is.logical(return_matched)) stop("`return_matched` must be TRUE or FALSE.")
  if (!is.logical(bootstrap)) stop("`bootstrap` must be TRUE or FALSE.")
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1) stop("`n_boot` must be a single positive number.")
  if (!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) stop("`conf` must be a single number between 0 and 1.")

  # Prepare the two-drug cohort ----
  # The working columns are added as copies rather than renaming the originals,
  # so every column name the caller may refer to (in matching_var, match.exact,
  # match.antiexact, caliper or match_args) is still present under that name.
  internal_cols <- c(dataset_benefit = benefit_var, dataset_drug_var = drug_var,
                     dataset_outcome_var = outcome_var)
  clash <- names(internal_cols)[names(internal_cols) %in% colnames(data) &
                                  names(internal_cols) != internal_cols]
  if (length(clash) > 0) {
    stop("`data` must not contain column(s) ", paste(clash, collapse = ", "),
         ", which are used internally.")
  }

  initial_dataset <- data
  for (col in names(internal_cols)) {
    initial_dataset[[col]] <- data[[internal_cols[[col]]]]
  }
  initial_dataset <- initial_dataset[initial_dataset$dataset_drug_var %in% drugs, , drop = FALSE]

  # Matching does not depend on cal_groups, so it is done once here and the same
  # matched cohort is reused for every requested number of calibration groups.
  match_model <- NULL

  if (isTRUE(matching)) {
    initial_dataset <- initial_dataset %>%
      dplyr::mutate(conc_disc_label = ifelse(dataset_benefit <= 0, 1, 0)) %>%
      tidyr::drop_na(dplyr::all_of(unique(c(matching_var, match.exact, match.antiexact))))

    # Matching needs both groups present; if only one is, there is nothing to
    # compare and no results can be produced for this drug pair.
    if (length(unique(initial_dataset$conc_disc_label)) < 2) {
      if (isTRUE(return_matched)) {
        return(list(results = NULL, matched = NULL, matchit = NULL))
      }
      return(NULL)
    }

    matching_formula <- .build_match_formula(
      initial_dataset, matching_var, match.exact, match.antiexact,
      treat_var = "conc_disc_label"
    )

    # The received drug is always anti-exact, so a matched set never compares
    # patients on the same drug; it is re-applied after match_args.
    # When no patient can be matched (e.g. a tight caliper, or no admissible
    # control on the other drug), MatchIt stops with "No units were matched".
    # That is reported as a warning so that, inside pairwise_calibration(), one
    # unmatchable drug pair does not abort the remaining pairs.
    forced_antiexact <- "dataset_drug_var"
    match_model <- tryCatch(
      .run_match(
        initial_dataset, matching_formula,
        exact_vars = match.exact,
        antiexact_vars = unique(c(setdiff(match.antiexact, drug_var), forced_antiexact)),
        caliper = caliper, seed = seed,
        method = method, distance = distance, replace = replace, ratio = ratio,
        std.caliper = std.caliper, match_args = match_args,
        forced_antiexact = forced_antiexact
      ),
      error = function(e) {
        if (grepl("no units were matched", conditionMessage(e), ignore.case = TRUE)) return(NULL)
        stop(e)
      }
    )

    calibration_data <- if (is.null(match_model)) initial_dataset[0, , drop = FALSE] else
      .get_matched_data(match_model, initial_dataset)

    if (nrow(calibration_data) == 0) {
      warning("No matches could be formed for ", drugs[1], " vs ", drugs[2],
              "; no results are returned for this pair.")
      if (isTRUE(return_matched)) {
        return(list(results = NULL, matched = NULL, matchit = match_model))
      }
      return(NULL)
    }

  } else {
    calibration_data <- initial_dataset
  }

  result <- NULL

  for (cg in cal_groups) {

    coef      <- rep(NA_real_, cg)
    coef_low  <- rep(NA_real_, cg)
    coef_high <- rep(NA_real_, cg)
    mean_vals <- rep(NA_real_, cg)
    n_drug1   <- rep(0, cg)
    n_drug2   <- rep(0, cg)
    no_ci_groups <- integer(0)

    grouped_data <- calibration_data %>%
      dplyr::mutate(grouping = dplyr::ntile(dataset_benefit, cg))

    for (g in seq_len(cg)) {
      group_data <- grouped_data %>%
        dplyr::filter(grouping == g) %>%
        dplyr::mutate(dataset_drug_var = factor(dataset_drug_var, levels = rev(drugs)))

      mean_vals[g] <- mean(group_data$dataset_benefit, na.rm = TRUE)
      n_drug1[g] <- nrow(dplyr::filter(group_data, dataset_drug_var == drugs[1]))
      n_drug2[g] <- nrow(dplyr::filter(group_data, dataset_drug_var == drugs[2]))

      if (length(unique(group_data$dataset_drug_var)) < 2) next

      model <- stats::glm(
        stats::as.formula(.pairwise_group_formula(group_data, adjustment_var)),
        data = group_data
      )

      coef[g] <- stats::coef(model)[2]

      if (!isTRUE(bootstrap)) {
        # With no residual degrees of freedom the dispersion, and so the
        # interval, cannot be estimated.
        if (model$df.residual >= 1 && !is.na(coef[g])) {
          ci <- suppressMessages(stats::confint.default(model, level = conf))
          coef_low[g]  <- ci[2, 1]
          coef_high[g] <- ci[2, 2]
        } else if (!is.na(coef[g])) {
          no_ci_groups <- c(no_ci_groups, g)
        }
      }
    }

    if (length(no_ci_groups) > 0) {
      warning(
        drugs[1], " vs ", drugs[2], ", cal_groups = ", cg, ": group(s) ",
        paste(no_ci_groups, collapse = ", "),
        " have too few patients to estimate a confidence interval ",
        "(coef_low / coef_high set to NA)."
      )
    }

    # Percentile bootstrap intervals ----
    # The group labels above are held fixed, so every replicate re-estimates the
    # effect for the same groups. Matched sets are resampled as a unit so that
    # matched patients stay together.
    if (isTRUE(bootstrap)) {
      if (!is.null(boot_seed)) set.seed(boot_seed)

      if (isTRUE(matching)) {
        unit <- as.character(grouped_data$subclass)
      } else {
        unit <- as.character(seq_len(nrow(grouped_data)))
      }
      unit_rows <- split(seq_along(unit), unit)
      unique_units <- names(unit_rows)

      boot_coef <- matrix(NA_real_, nrow = n_boot, ncol = cg)

      for (b in seq_len(n_boot)) {
        sampled <- sample(unique_units, length(unique_units), replace = TRUE)
        rep_data <- grouped_data[unlist(unit_rows[sampled], use.names = FALSE), , drop = FALSE]

        for (g in seq_len(cg)) {
          boot_coef[b, g] <- .pairwise_group_coef(
            rep_data[rep_data$grouping == g, , drop = FALSE],
            drugs, adjustment_var
          )
        }
      }

      for (g in seq_len(cg)) {
        if (is.na(coef[g]) || all(is.na(boot_coef[, g]))) next
        ci <- .percentile_ci(boot_coef[, g], conf = conf)
        coef_low[g]  <- ci[1]
        coef_high[g] <- ci[2]
      }
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

  if (isTRUE(return_matched)) {
    if (!isTRUE(matching)) {
      warning("`return_matched = TRUE` but `matching = FALSE`; no matched cohort to return.")
      return(list(results = result, matched = NULL, matchit = NULL))
    }
    return(list(results = result, matched = calibration_data, matchit = match_model))
  }

  return(result)
}
