#' Match Concordant and Discordant Patients for Overall-Benefit Calibration
#'
#' Builds the matched patient pairs used to validate a treatment-selection
#' model's predicted *overall benefit*. Each patient is labelled concordant
#' (received the predicted-best drug) or discordant (received a different drug),
#' and concordant patients are matched to clinically similar discordant patients
#' by nearest-neighbour Mahalanobis distance, exactly on the predicted-best drug
#' and forced to differ on the drug actually received.
#'
#' Matching is run in both directions by default (\code{symmetric = TRUE}):
#' once with concordant patients as the treated group and once with the roles
#' reversed. This returns three matched populations that can each be inspected
#' or summarised separately:
#' \itemize{
#'   \item \code{concordant} - concordant patients matched to discordant controls;
#'   \item \code{discordant} - discordant patients matched to concordant controls;
#'   \item \code{combined}   - the two sets pooled, with pairs that appear in both
#'         directions counted once.
#' }
#'
#' The predicted and observed benefits are computed on the scale of
#' \code{outcome_var} (i.e. absolute outcomes): for each pair, \code{pred.benefit}
#' is the concordant patient's predicted outcome on the discordant drug minus the
#' predicted outcome on their own (recommended) drug, and \code{obs.benefit} is
#' the observed outcome of the discordant patient minus that of the concordant
#' patient. Both are positive when the model's recommended drug is expected to,
#' and does, give the better outcome.
#'
#' @param data A data frame with one row per patient, containing the received
#'   drug, the observed outcome, the predicted outcome for each drug and the
#'   matching covariates.
#' @param drug_var Character. Name of the column holding the drug actually
#'   received.
#' @param outcome_var Character. Name of the observed outcome column (e.g. the
#'   12-month HbA1c). Benefits are computed on this absolute scale.
#' @param pred_cols Character vector. Names of the predicted-outcome columns, one
#'   per drug, sharing a common prefix (e.g. \code{c("pred.SGLT2", "pred.GLP1")}).
#' @param matching_var Character vector. Covariates used for Mahalanobis distance
#'   matching.
#' @param match.exact Optional character vector. Additional variables to match on
#'   exactly. The predicted-best drug is always included.
#' @param match.antiexact Optional character vector. Additional variables on which
#'   a matched pair must differ. The received-drug variable is always included.
#' @param caliper Optional named numeric vector passed to \code{MatchIt::matchit}
#'   (e.g. \code{c(prehba1c = 0.4)}); interpreted in standard-deviation units by
#'   default (see \code{std.caliper}).
#' @param std.caliper Logical. If \code{TRUE} (default), \code{caliper} is read
#'   in standard-deviation units; if \code{FALSE}, in the variable's own units.
#' @param method Character. Matching method passed to \code{MatchIt::matchit}
#'   (default \code{"nearest"}).
#' @param distance Character. Distance measure passed to
#'   \code{MatchIt::matchit} (default \code{"mahalanobis"}).
#' @param replace Logical. Whether a control patient may be reused in more than
#'   one pair (default \code{TRUE}). See the note on confidence intervals below.
#' @param ratio Numeric. Number of controls matched to each treated patient
#'   (default \code{1}).
#' @param match_args Optional named list of further arguments for
#'   \code{MatchIt::matchit}. Applied last, so it overrides the arguments above
#'   and gives access to any matching option not exposed here.
#' @param conc_tolerance Optional numeric. If supplied, a patient is concordant
#'   when they received any drug within this tolerance of the best predicted
#'   outcome, and a discordant match must differ from every drug in that set.
#' @param symmetric Logical. If \code{TRUE} (default), match in both directions
#'   and return the pooled, de-duplicated set in \code{combined}. If \code{FALSE},
#'   only the concordant-as-treated direction is run.
#' @param seed Optional integer. Seed set before each matching call so the result
#'   is reproducible.
#'
#' @return An object of class \code{cateval_benefit_match}: a list with elements
#' \describe{
#'   \item{concordant}{Matched pairs, concordant patients as treated.}
#'   \item{discordant}{Matched pairs, discordant patients as treated
#'     (\code{NULL} if \code{symmetric = FALSE}).}
#'   \item{combined}{The concordant and discordant sets pooled, with duplicate
#'     pairs removed.}
#'   \item{matchit}{A list of the underlying \code{matchit} objects
#'     (\code{$concordant} and \code{$discordant}) for balance / love plots.}
#' }
#' Each pairs data frame has one row per matched pair, including
#' \code{concordant_drug}, \code{discordant_drug}, \code{index_drug} (the drug
#' received by the index/treated patient of the pair), \code{pred.benefit},
#' \code{obs.benefit} and the retained matching covariates.
#'
#' @details
#' The matched pairs frames are the input to \code{\link{benefit_by_group}},
#' \code{\link{benefit_overall}}, \code{\link{benefit_citl_slope}} and
#' \code{\link{benefit_by_drug}}.
#'
#' \strong{Limitation (confidence intervals under matching with replacement).}
#' Matching is done with replacement by default (\code{replace = TRUE}), so a
#' single control patient can appear in
#' several matched pairs, which makes pairs that share a control not independent.
#' The confidence intervals from the summary functions -- both the analytical
#' interval and the optional pair-level bootstrap -- treat matched pairs as
#' independent, and so they do not widen to reflect this dependence. They are
#' therefore likely to be somewhat too narrow (anti-conservative) when control
#' reuse is substantial. Valid variance estimation for matching with replacement
#' is a known hard problem; the intervals here should be read as approximate, and
#' this limitation stated when they are reported. The extent of reuse can be
#' inspected directly, e.g.
#' \code{table(table(c(pairs$cateval_row_id_conc, pairs$cateval_row_id_disc)))}.
#' Setting \code{replace = FALSE} removes this dependence (each control is used
#' at most once), at the cost of fewer and generally poorer-quality matches.
#'
#' @seealso \code{\link{benefit_by_group}}, \code{\link{benefit_overall}},
#'   \code{\link{benefit_citl_slope}}, \code{\link{benefit_by_drug}}
#'
#' @examples
#' \dontrun{
#' matched <- match_benefit_pairs(
#'   data = analysis_dataset,
#'   drug_var = "drugclass",
#'   outcome_var = "posthba1cfinal",
#'   pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "SU", "TZD")),
#'   matching_var = c("agetx", "prebmi", "prehba1c", "sex"),
#'   match.exact = c("sex"),
#'   caliper = c(prehba1c = 0.4),
#'   seed = 19840503
#' )
#'
#' # inspect balance
#' plot(summary(matched$matchit$concordant))
#' # summarise the pooled set
#' benefit_by_group(matched$combined, cal_groups = 10)
#' }
#'
#' @importFrom MatchIt matchit get_matches
#' @importFrom dplyr filter select inner_join rename mutate bind_rows all_of
#' @importFrom stats as.formula
#' @export
match_benefit_pairs <- function(data,
                                drug_var,
                                outcome_var,
                                pred_cols,
                                matching_var,
                                match.exact = NULL,
                                match.antiexact = NULL,
                                caliper = NULL,
                                std.caliper = TRUE,
                                method = "nearest",
                                distance = "mahalanobis",
                                replace = TRUE,
                                ratio = 1,
                                match_args = list(),
                                conc_tolerance = NULL,
                                symmetric = TRUE,
                                seed = NULL) {

  `%>%` <- dplyr::`%>%`

  # ---- Input validation ----
  if (!(drug_var %in% colnames(data))) stop("`drug_var` not found in data.")
  if (!(outcome_var %in% colnames(data))) stop("`outcome_var` not found in data.")
  if (!all(pred_cols %in% colnames(data))) stop("Some `pred_cols` not found in data.")
  if (length(pred_cols) < 2) stop("`pred_cols` must contain at least two drugs.")
  if (!all(matching_var %in% colnames(data))) stop("Some `matching_var` not found in data.")
  if (!is.null(match.exact) && !all(match.exact %in% colnames(data))) stop("Some `match.exact` not found in data.")
  if (!is.null(match.antiexact) && !all(match.antiexact %in% colnames(data))) stop("Some `match.antiexact` not found in data.")
  if (!is.null(conc_tolerance) && !is.numeric(conc_tolerance)) stop("`conc_tolerance` must be numeric.")
  if (!is.logical(symmetric)) stop("`symmetric` must be TRUE or FALSE.")
  if (!is.null(caliper) && !is.numeric(caliper)) stop("`caliper` must be a named numeric vector.")
  if (!is.logical(std.caliper)) stop("`std.caliper` must be TRUE or FALSE.")
  if (!is.character(method) || length(method) != 1) stop("`method` must be a single character string.")
  if (!is.character(distance) || length(distance) != 1) stop("`distance` must be a single character string.")
  if (!is.logical(replace)) stop("`replace` must be TRUE or FALSE.")
  if (!is.numeric(ratio) || length(ratio) != 1 || ratio < 1) stop("`ratio` must be a single number of at least 1.")
  if (!is.list(match_args)) stop("`match_args` must be a list.")

  prefix <- .common_prefix(pred_cols)

  # ---- Label concordance and keep only complete matching covariates ----
  labelled <- .label_concordance(data, drug_var, pred_cols, prefix, conc_tolerance)
  antiexact_extra <- labelled$antiexact_extra
  labelled <- labelled$data %>%
    tidyr::drop_na(dplyr::all_of(matching_var))

  # Stable row identifier used to de-duplicate pairs across the two directions.
  labelled$cateval_row_id <- seq_len(nrow(labelled))

  # ---- Matching specification (shared by both directions) ----
  formula_str <- .build_match_formula(labelled, matching_var, match.exact, match.antiexact)
  exact_vars <- unique(c("cateval_best_drug", match.exact))
  antiexact_vars <- unique(c(drug_var, match.antiexact, antiexact_extra))

  # Covariates to carry through to the pairs output for inspection.
  keep_vars <- setdiff(matching_var, pred_cols)

  # ---- Direction 1: concordant patients as the treated group ----
  data_conc <- labelled
  data_conc$cateval_treated <- data_conc$cateval_concordant
  match_conc <- .run_benefit_match(
    data_conc, formula_str, exact_vars, antiexact_vars, caliper, seed,
    method = method, distance = distance, replace = replace, ratio = ratio,
    std.caliper = std.caliper, match_args = match_args
  )
  pairs_conc <- .extract_benefit_pairs(
    MatchIt::get_matches(match_conc, data = data_conc),
    drug_var, outcome_var, pred_cols, prefix, keep_vars
  )

  # ---- Direction 2: discordant patients as the treated group (flipped) ----
  if (isTRUE(symmetric)) {
    data_disc <- labelled
    data_disc$cateval_treated <- 1 - data_disc$cateval_concordant
    match_disc <- .run_benefit_match(
      data_disc, formula_str, exact_vars, antiexact_vars, caliper, seed,
      method = method, distance = distance, replace = replace, ratio = ratio,
      std.caliper = std.caliper, match_args = match_args
    )
    pairs_disc <- .extract_benefit_pairs(
      MatchIt::get_matches(match_disc, data = data_disc),
      drug_var, outcome_var, pred_cols, prefix, keep_vars
    )

    # Pool the two directions, keeping every concordant-as-treated pair and only
    # the discordant-as-treated pairs that are not already represented.
    combined <- dplyr::bind_rows(
      pairs_conc,
      dplyr::filter(pairs_disc, !pair_key %in% pairs_conc$pair_key)
    )
  } else {
    match_disc <- NULL
    pairs_disc <- NULL
    combined <- pairs_conc
  }

  # ---- Assemble output ----
  out <- list(
    concordant = pairs_conc,
    discordant = pairs_disc,
    combined = combined,
    matchit = list(concordant = match_conc, discordant = match_disc)
  )
  class(out) <- "cateval_benefit_match"
  out
}
