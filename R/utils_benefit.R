# Internal helpers for the overall-benefit calibration workflow.
#
# These functions are not exported. They are shared by match_benefit_pairs()
# and the benefit_* summary functions so that the matching logic, the
# concordant/discordant labelling and the confidence-interval helpers live in
# one place. Everything here works on the "predicted-best drug" convention:
# for each patient the recommended drug is the one with the lowest predicted
# outcome, and a patient is concordant if they received that drug.


# Longest common prefix of the prediction column names ----------------------
#
# The prediction columns share a prefix (e.g. "pred." in "pred.SGLT2"). We
# strip it to recover the drug names and rebuild individual column names when
# looking up a specific drug's prediction.
#
# @param strings Character vector of column names.
# @return A single character string with the shared leading characters.
# @noRd
.common_prefix <- function(strings) {
  min_len <- min(nchar(strings))
  prefix <- ""
  for (i in seq_len(min_len)) {
    current_chars <- substr(strings, i, i)
    if (length(unique(current_chars)) == 1) {
      prefix <- paste0(prefix, current_chars[1])
    } else {
      break
    }
  }
  prefix
}


# Label each patient as concordant or discordant ----------------------------
#
# Adds three internal columns to the data:
#   cateval_best_drug  - the predicted-best (rank 1) drug name
#   cateval_concordant - 1 if the received drug is the recommended drug
#                        (or, with a tolerance, any drug within the tolerance
#                        of the best), 0 otherwise
#   cateval_tol_*      - only when a tolerance is used: the individual drugs in
#                        the tolerance set, used later for anti-exact matching
#
# @param data Data frame containing the drug column and the prediction columns.
# @param drug_var Name of the received-drug column.
# @param pred_cols Prediction column names.
# @param prefix Shared prefix of pred_cols.
# @param conc_tolerance Optional numeric tolerance defining the best-drug set.
# @return A list with the labelled data and the extra anti-exact variables.
# @noRd
.label_concordance <- function(data, drug_var, pred_cols, prefix, conc_tolerance = NULL) {

  `%>%` <- dplyr::`%>%`

  # Predicted-best (rank 1) drug for every patient.
  labelled <- get_best_drugs(
    data = data,
    rank = 1,
    column_names = pred_cols,
    final_var_name = prefix
  ) %>%
    dplyr::rename(cateval_best_drug = paste0(prefix, "rank1_drug_name"))

  if (is.null(conc_tolerance)) {

    # Strict definition: concordant means received the single recommended drug.
    labelled <- labelled %>%
      dplyr::mutate(cateval_concordant = ifelse(.data[[drug_var]] == cateval_best_drug, 1, 0))

    antiexact_extra <- NULL

  } else {

    # Tolerance definition: concordant means received any drug whose predicted
    # outcome is within conc_tolerance of the best. The individual drugs in the
    # tolerance set are spread across cateval_tol_* columns so that a discordant
    # match can be forced to differ from every drug in the set.
    n_drugs <- length(pred_cols)
    tol_vars <- paste0("cateval_tol_", seq_len(n_drugs))

    labelled <- get_best_drugs(
      data = labelled,
      tolerance = conc_tolerance,
      column_names = pred_cols,
      final_var_name = prefix
    ) %>%
      dplyr::rename(cateval_tol_set = paste0(prefix, "within_", conc_tolerance, "_of_best_drug_name")) %>%
      dplyr::mutate(
        cateval_concordant = ifelse(
          stringr::str_detect(cateval_tol_set, paste0("\\b", .data[[drug_var]], "\\b")),
          1, 0
        ),
        cateval_tol_list = stringr::str_split(cateval_tol_set, "\\s*[;,\\s]\\s*")
      ) %>%
      dplyr::mutate(cateval_tol_list = purrr::map(cateval_tol_list, ~ rep(.x, length.out = n_drugs))) %>%
      dplyr::mutate(cateval_tol_list = purrr::map(cateval_tol_list, ~ purrr::set_names(.x, tol_vars))) %>%
      tidyr::unnest_wider(cateval_tol_list) %>%
      dplyr::mutate(dplyr::across(
        dplyr::all_of(tol_vars),
        ~ dplyr::if_else(cateval_concordant == 0, .data[[drug_var]], .x)
      ))

    antiexact_extra <- tol_vars
  }

  list(data = labelled, antiexact_extra = antiexact_extra)
}


# Build the matching formula from the matching covariates -------------------
#
# Continuous covariates enter directly (Mahalanobis distance handles them);
# categorical covariates are added only when they have more than one level in
# the data. Variables handled by exact / anti-exact matching are dropped from
# the distance formula.
#
# @return A formula string of the form "cateval_treated ~ x1 + x2 + ...".
# @noRd
.build_match_formula <- function(data, matching_var, match.exact, match.antiexact) {
  categorical_vars <- matching_var[sapply(data[matching_var], function(x) is.factor(x) || is.character(x))]
  cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))

  formula_str <- paste("cateval_treated ~", paste(cont_vars, collapse = " + "))
  for (v in categorical_vars) {
    if (length(unique(data[[v]])) > 1) {
      formula_str <- paste(formula_str, "+", v)
    }
  }
  formula_str
}


# Run one match -------------------------------------------------------------
#
# A single matching direction. The treated group is whichever group the caller
# has put in the cateval_treated column (concordant patients in one direction,
# discordant patients in the other). The matching is exact on the predicted-best
# drug and forced to differ on the received drug; the remaining behaviour is
# controlled by the caller and defaults to 1:1 nearest-neighbour Mahalanobis
# matching with replacement.
#
# match_args is applied last, so it takes precedence over the named arguments
# and can supply any other MatchIt::matchit() argument.
#
# @return A MatchIt "matchit" object.
# @noRd
.run_benefit_match <- function(data, formula_str, exact_vars, antiexact_vars,
                               caliper, seed,
                               method = "nearest",
                               distance = "mahalanobis",
                               replace = TRUE,
                               ratio = 1,
                               std.caliper = TRUE,
                               match_args = list()) {
  if (!is.null(seed)) set.seed(seed)

  args <- list(
    formula = stats::as.formula(formula_str),
    data = data,
    method = method,
    distance = distance,
    replace = replace,
    ratio = ratio,
    exact = exact_vars,
    antiexact = antiexact_vars,
    caliper = caliper,
    std.caliper = std.caliper
  )

  do.call(MatchIt::matchit, utils::modifyList(args, match_args))
}


# Turn a matched data set into one row per matched pair ---------------------
#
# Each subclass produced by MatchIt contains exactly one concordant and one
# discordant patient, and exactly one treated ("index") and one control patient.
# The index patient is whichever group was treated in this matching direction:
# the concordant patient when matching concordant-as-treated, the discordant
# patient when matching discordant-as-treated. Each pair is anchored on the
# index patient, and we compute:
#   pred.benefit - the index patient's predicted advantage of the recommended
#                  drug over the discordant (non-recommended) drug of the pair,
#                  i.e. index_pred[discordant_drug] - index_pred[recommended_drug]
#                  (positive = model expects the recommended drug to give a
#                  lower / better outcome)
#   obs.benefit  - the observed outcome difference between the discordant and
#                  concordant patient on the same scale
#   index_drug   - the drug the index (treated) patient received; used to split
#                  the calibration by drug class
#
# @param matched Output of MatchIt::get_matches(), carrying cateval_concordant,
#   cateval_treated, cateval_row_id, the drug column, the outcome column and the
#   pred columns.
# @param keep_vars Extra covariate columns to retain for inspection.
# @return A data frame with one row per matched pair.
# @noRd
.extract_benefit_pairs <- function(matched, drug_var, outcome_var, pred_cols, prefix, keep_vars = NULL) {

  `%>%` <- dplyr::`%>%`

  # Index (treated) patient of each pair. The predicted benefit is evaluated
  # from this patient's predictions, so we keep their full prediction vector,
  # their own drug and any covariates to inspect.
  index <- matched %>%
    dplyr::filter(cateval_treated == 1) %>%
    dplyr::select(
      subclass,
      index_drug = dplyr::all_of(drug_var),
      dplyr::all_of(pred_cols),
      dplyr::all_of(keep_vars)
    )

  # Concordant patient of the pair (received the recommended drug).
  concordant <- matched %>%
    dplyr::filter(cateval_concordant == 1) %>%
    dplyr::select(
      subclass,
      cateval_row_id_conc = cateval_row_id,
      concordant_drug = dplyr::all_of(drug_var),
      outcome_conc = dplyr::all_of(outcome_var)
    )

  # Discordant patient of the pair (received a different drug).
  discordant <- matched %>%
    dplyr::filter(cateval_concordant == 0) %>%
    dplyr::select(
      subclass,
      cateval_row_id_disc = cateval_row_id,
      discordant_drug = dplyr::all_of(drug_var),
      outcome_disc = dplyr::all_of(outcome_var)
    )

  pairs <- concordant %>%
    dplyr::inner_join(discordant, by = "subclass") %>%
    dplyr::inner_join(index, by = "subclass") %>%
    dplyr::filter(concordant_drug != discordant_drug)

  # Predicted benefit from the index patient's predictions (matrix indexing is
  # fast and safe for an arbitrary number of drug classes): the recommended drug
  # is the concordant patient's drug, the alternative is the discordant drug.
  pred_matrix <- as.matrix(pairs[, pred_cols])
  colnames(pred_matrix) <- sub(prefix, "", pred_cols, fixed = TRUE)
  row_index <- seq_len(nrow(pairs))
  pred_recommended <- pred_matrix[cbind(row_index, match(pairs$concordant_drug, colnames(pred_matrix)))]
  pred_alternative <- pred_matrix[cbind(row_index, match(pairs$discordant_drug, colnames(pred_matrix)))]

  pairs$pred.benefit <- pred_alternative - pred_recommended
  pairs$obs.benefit <- pairs$outcome_disc - pairs$outcome_conc

  # A direction-independent key for a pair of patients, used to drop duplicates
  # when the two matching directions are combined.
  pairs$pair_key <- paste(
    pmin(pairs$cateval_row_id_conc, pairs$cateval_row_id_disc),
    pmax(pairs$cateval_row_id_conc, pairs$cateval_row_id_disc),
    sep = "_"
  )

  pairs
}


# Outcome-model formula for one calibration group ---------------------------
#
# The treatment effect is the second coefficient of a model of the outcome on
# the received drug, optionally adjusted for covariates. Categorical adjustment
# variables are included only when they vary within the group, since a constant
# factor cannot be estimated.
#
# @return A formula string.
# @noRd
.pairwise_group_formula <- function(group_data, adjustment_var) {
  formula_str <- "dataset_outcome_var ~ dataset_drug_var"

  if (!is.null(adjustment_var)) {
    is_cat <- sapply(
      group_data[, adjustment_var, drop = FALSE],
      function(x) is.factor(x) || is.character(x)
    )
    cat_vars  <- adjustment_var[is_cat]
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

  formula_str
}


# Treatment-effect coefficient for one calibration group --------------------
#
# Used for the bootstrap replicates, where a resampled group may happen to
# contain only one of the two drugs (or produce an unestimable model). Such a
# replicate contributes NA rather than stopping the bootstrap.
#
# @return A single number, or NA_real_ when the effect cannot be estimated.
# @noRd
.pairwise_group_coef <- function(group_data, drugs, adjustment_var) {
  group_data$dataset_drug_var <- factor(group_data$dataset_drug_var, levels = rev(drugs))

  if (length(unique(stats::na.omit(group_data$dataset_drug_var))) < 2) {
    return(NA_real_)
  }

  model <- try(
    stats::glm(stats::as.formula(.pairwise_group_formula(group_data, adjustment_var)), data = group_data),
    silent = TRUE
  )
  if (inherits(model, "try-error")) return(NA_real_)

  unname(stats::coef(model)[2])
}


# Percentile confidence interval from a vector of bootstrap statistics ------
#
# @param t Numeric vector of bootstrap replicates.
# @param conf Confidence level (e.g. 0.95).
# @return A length-2 numeric vector: lower and upper bounds.
# @noRd
.percentile_ci <- function(t, conf = 0.95) {
  alpha <- (1 - conf) / 2
  as.numeric(stats::quantile(t, probs = c(alpha, 1 - alpha), na.rm = TRUE))
}
