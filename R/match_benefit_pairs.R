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
#' \code{outcome_var} (i.e. absolute outcomes). Each pair has an \emph{index}
#' patient -- the treated patient of the matching direction it came from -- and
#' \code{pred.benefit} is computed from that patient's predictions: their
#' predicted outcome on the discordant patient's drug minus their predicted
#' outcome on the concordant patient's (recommended) drug. In the
#' concordant-as-treated direction the index patient is the concordant patient;
#' in the discordant-as-treated direction it is the discordant patient. Because
#' both patients of a pair share the same predicted-best drug (exact matching),
#' either way it is the predicted advantage of the recommended drug over the
#' alternative for the index patient. \code{obs.benefit} is the observed outcome
#' of the discordant patient minus that of the concordant patient. Both are
#' positive when the model's recommended drug is expected to, and does, give the
#' better (lower) outcome.
#'
#' @param data A data frame with one row per patient, containing the received
#'   drug, the observed outcome, the predicted outcome for each drug and the
#'   matching covariates.
#' @param drug_var Character. Name of the column holding the drug actually
#'   received.
#' @param outcome_var Character. Name of the observed outcome column (e.g. the
#'   12-month HbA1c). Benefits are computed on this absolute scale.
#' @param pred_cols Character vector. Names of the predicted-outcome columns, one
#'   per drug (e.g. \code{c("pred.SGLT2", "pred.GLP1")}). Lower predicted values
#'   are taken to be better. Any naming scheme can be used; see \code{drugs} and
#'   the \emph{Drug names} section for how each column is linked to a drug. The
#'   vector may also be named by drug, e.g.
#'   \code{c(SGLT2 = "yhat_sglt2", GLP1 = "yhat_glp1")}.
#' @param matching_var Character vector. Covariates used for Mahalanobis distance
#'   matching. Must not include \code{drug_var} or \code{outcome_var}. Patients
#'   with a missing value in any \code{matching_var}, \code{match.exact} or
#'   \code{match.antiexact} variable are excluded before matching.
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
#'   (default \code{1}). With \code{ratio > 1} each treated patient contributes
#'   one pair per matched control; see \emph{Matching more than one control}.
#' @param match_args Optional named list of further arguments for
#'   \code{MatchIt::matchit}, applied to both matching directions. Overrides the
#'   arguments above and gives access to any matching option not exposed here.
#' @param match_args_conc Optional named list of \code{MatchIt::matchit}
#'   arguments applied only to the concordant-as-treated direction, on top of
#'   \code{match_args}.
#' @param match_args_disc As \code{match_args_conc}, for the
#'   discordant-as-treated direction. Ignored (with a warning) when
#'   \code{symmetric = FALSE}, since that direction is not run.
#' @param conc_tolerance Optional numeric. If supplied, a patient is concordant
#'   when they received any drug within this tolerance of the best predicted
#'   outcome, and a discordant match must differ from every drug in that set.
#' @param symmetric Logical. If \code{TRUE} (default), match in both directions
#'   and return the pooled, de-duplicated set in \code{combined}. If \code{FALSE},
#'   only the concordant-as-treated direction is run.
#' @param seed Optional integer. Seed set before each matching call so the result
#'   is reproducible.
#' @param drugs Optional vector of drug names (character, factor or numeric
#'   codes), one per element of \code{pred_cols} and in the same order, matching
#'   the values in \code{drug_var}. Takes precedence over the names of \code{pred_cols}. If
#'   neither is given, the drug names are inferred (see \emph{Drug names}).
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
#'   \item{drugs}{Named character vector linking each prediction column (names)
#'     to the drug it was taken to refer to (values), so the mapping can be
#'     checked.}
#' }
#' Each pairs data frame has one row per treated-control pair, including
#' \code{concordant_drug}, \code{discordant_drug}, \code{index_drug} (the drug
#' received by the index/treated patient of the pair), \code{pred.benefit},
#' \code{obs.benefit} and the retained matching covariates.
#'
#' @details
#' The matched pairs frames are the input to \code{\link{benefit_by_group}},
#' \code{\link{benefit_overall}}, \code{\link{benefit_citl_slope}} and
#' \code{\link{benefit_by_drug}}.
#'
#' \strong{Order of precedence for matching arguments.} Each direction is built
#' from the package defaults, then the named arguments above (\code{method},
#' \code{distance}, \code{replace}, \code{ratio}, \code{caliper},
#' \code{std.caliper}), then \code{match_args}, then that direction's own
#' \code{match_args_conc} or \code{match_args_disc}. Later links win, so a
#' direction list overrides the shared list, which overrides the named
#' arguments. Because \code{utils::modifyList} treats \code{NULL} as removal,
#' \code{match_args_disc = list(caliper = NULL)} unsets the caliper for the
#' discordant direction only.
#'
#' \strong{Protected constraints.} Matching is always exact on the
#' predicted-best drug and always anti-exact on \code{drug_var} (and on the
#' tolerance-set columns when \code{conc_tolerance} is used). These are
#' re-applied after every override, so supplying \code{exact} or
#' \code{antiexact} through \code{match_args} adds to them rather than replacing
#' them. For the same reason \code{drug_var} cannot be used in
#' \code{match.exact}, and neither \code{drug_var} nor \code{outcome_var} can be
#' a matching variable. Only pair-matching methods (e.g. \code{"nearest"},
#' \code{"optimal"}, \code{"genetic"}) can be used, since each treated patient
#' must be linked to its own control(s).
#'
#' \strong{Drug names.} The package makes no assumption about what the drugs
#' are called: they can be drug classes, molecules, trial arms or numeric
#' codes, as long as \code{drug_var} and the prediction columns refer to the
#' same set. Each column in \code{pred_cols} must be linked to the drug name
#' used in \code{drug_var}. This is taken from \code{drugs} if supplied,
#' otherwise from the names of \code{pred_cols}. Failing both, it is inferred
#' from the drugs observed in \code{drug_var}:
#' \enumerate{
#'   \item each column is assigned the drug its name ends with, so any prefix
#'     works (\code{"pred.A"}, \code{"yhat_outcome_A"}); if that does not give
#'     a unique drug for every column, the drug its name starts with
#'     (\code{"A_pred"}). Where several drugs fit, the longest is taken, so
#'     \code{"pred.BA"} is \code{BA} rather than \code{A};
#'   \item failing that, the same matching ignoring upper/lower case
#'     (\code{"pred.drug_a"} for \code{"Drug_A"}), keeping the spelling in
#'     \code{drug_var};
#'   \item a drug that has a prediction column but that no patient received
#'     cannot be found in \code{drug_var}; if the other columns all follow one
#'     pattern (e.g. \code{"pred.<drug>"}), its name is read from its column
#'     with that pattern.
#' }
#' Steps 2 and 3 print the mapping they used. If none of these gives an
#' unambiguous answer, an error asks for \code{drugs}. The mapping used is
#' always returned in the \code{drugs} element. Patients who received a drug
#' with no prediction column, or whose received drug is missing, cannot be
#' labelled or paired and are excluded with a warning.
#'
#' \strong{Matching more than one control (\code{ratio > 1}).} Each treated
#' patient is paired with every control matched to it, so it contributes up to
#' \code{ratio} rows that share the same treated patient. The summary
#' functions weight every row equally; MatchIt's matching weights are not used.
#' With a fixed ratio and every treated patient receiving its full set of
#' controls, each treated patient therefore carries the same weight. When some
#' receive fewer controls -- because of a caliper, a shortage of controls when
#' \code{replace = FALSE}, or a variable ratio set through \code{match_args}
#' -- those patients carry proportionally less weight. Rows that share a
#' treated patient are also not independent; see the note on confidence
#' intervals below.
#'
#' \strong{Limitation (confidence intervals under matching with replacement).}
#' Matching is done with replacement by default (\code{replace = TRUE}), so a
#' single control patient can appear in
#' several matched pairs, which makes pairs that share a control not independent
#' (as are pairs that share a treated patient when \code{ratio > 1}, and pairs
#' that share a patient across the two directions in \code{combined}).
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
#'
#' # the same specification, but matching the (larger) discordant direction
#' # without replacement and two controls per treated patient
#' matched2 <- match_benefit_pairs(
#'   data = analysis_dataset,
#'   drug_var = "drugclass",
#'   outcome_var = "posthba1cfinal",
#'   pred_cols = paste0("pred.", c("SGLT2", "GLP1", "DPP4", "SU", "TZD")),
#'   matching_var = c("agetx", "prebmi", "prehba1c", "sex"),
#'   match_args = list(m.order = "random"),
#'   match_args_disc = list(replace = FALSE, ratio = 2),
#'   seed = 19840503
#' )
#'
#' # prediction columns with a different naming scheme: name them by drug
#' matched3 <- match_benefit_pairs(
#'   data = analysis_dataset,
#'   drug_var = "drugclass",
#'   outcome_var = "posthba1cfinal",
#'   pred_cols = c(SGLT2 = "hba1c_sglt2_m12", GLP1 = "hba1c_glp1_m12",
#'                 DPP4 = "hba1c_dpp4_m12"),
#'   matching_var = c("agetx", "prebmi", "prehba1c", "sex")
#' )
#' matched3$drugs
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
                                match_args_conc = list(),
                                match_args_disc = list(),
                                conc_tolerance = NULL,
                                symmetric = TRUE,
                                seed = NULL,
                                drugs = NULL) {

  `%>%` <- dplyr::`%>%`

  # ---- Input validation ----
  if (!(drug_var %in% colnames(data))) stop("`drug_var` not found in data.")
  if (!(outcome_var %in% colnames(data))) stop("`outcome_var` not found in data.")
  if (!is.character(pred_cols)) stop("`pred_cols` must be a character vector of column names.")
  if (!all(pred_cols %in% colnames(data))) stop("Some `pred_cols` not found in data.")
  if (length(pred_cols) < 2) stop("`pred_cols` must contain at least two drugs.")
  if (anyDuplicated(pred_cols)) stop("`pred_cols` must not contain duplicated columns.")
  if (length(matching_var) == 0) stop("Provide at least one `matching_var`.")
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
  if (!is.list(match_args_conc)) stop("`match_args_conc` must be a list.")
  if (!is.list(match_args_disc)) stop("`match_args_disc` must be a list.")
  .check_named_list(match_args, "match_args")
  .check_named_list(match_args_conc, "match_args_conc")
  .check_named_list(match_args_disc, "match_args_disc")
  .check_reserved_vars(
    list(matching_var = matching_var, match.exact = match.exact),
    c(drug_var = drug_var, outcome_var = outcome_var)
  )
  .check_reserved_vars(list(match.antiexact = match.antiexact), c(outcome_var = outcome_var))
  if (!isTRUE(symmetric) && length(match_args_disc) > 0) {
    warning("`match_args_disc` is ignored when `symmetric = FALSE`; the discordant-as-treated direction is not run.")
  }

  # ---- Link each prediction column to a drug name in drug_var ----
  drug_names <- .resolve_drug_names(pred_cols, drugs, observed = data[[drug_var]])
  pred_cols <- unname(pred_cols)

  # Patients whose received drug has no prediction column (or is missing) can
  # never be concordant and have no prediction for their own drug, so their
  # pairs would have no predicted benefit. Exclude them up front.
  no_pred <- is.na(data[[drug_var]]) | !(as.character(data[[drug_var]]) %in% drug_names)
  if (any(no_pred)) {
    unknown <- unique(as.character(data[[drug_var]][no_pred]))
    warning(
      sum(no_pred), " patient(s) excluded because their `drug_var` value has no ",
      "prediction column: ", paste(ifelse(is.na(unknown), "<NA>", unknown), collapse = ", "),
      ". Supply `drugs` if this is a naming mismatch."
    )
    data <- data[!no_pred, , drop = FALSE]
  }

  # ---- Label concordance and keep only complete matching covariates ----
  labelled <- .label_concordance(data, drug_var, pred_cols, drug_names, conc_tolerance)
  antiexact_extra <- labelled$antiexact_extra
  labelled <- labelled$data %>%
    tidyr::drop_na(dplyr::all_of(unique(c(matching_var, match.exact, match.antiexact))))

  if (length(unique(labelled$cateval_concordant)) < 2) {
    stop("Matching needs both concordant and discordant patients; all ", nrow(labelled),
         " patients are ", if (all(labelled$cateval_concordant == 1)) "concordant" else "discordant",
         ". Check that `drugs` / `pred_cols` refer to the drugs in `drug_var`.")
  }

  # Stable row identifier used to de-duplicate pairs across the two directions.
  labelled$cateval_row_id <- seq_len(nrow(labelled))

  # ---- Matching specification (shared by both directions) ----
  formula_str <- .build_match_formula(labelled, matching_var, match.exact, match.antiexact,
                                     treat_var = "cateval_treated")

  # Split into the constraints the method requires and the caller's additions,
  # so the required ones can be re-applied after any match_args override.
  forced_exact <- "cateval_best_drug"
  forced_antiexact <- unique(c(drug_var, antiexact_extra))
  exact_vars <- unique(c(forced_exact, match.exact))
  antiexact_vars <- unique(c(drug_var, match.antiexact, antiexact_extra))

  # Covariates to carry through to the pairs output for inspection.
  keep_vars <- setdiff(matching_var, pred_cols)

  # ---- Direction 1: concordant patients as the treated group ----
  data_conc <- labelled
  data_conc$cateval_treated <- data_conc$cateval_concordant
  match_conc <- .run_match(
    data_conc, formula_str, exact_vars, antiexact_vars, caliper, seed,
    method = method, distance = distance, replace = replace, ratio = ratio,
    std.caliper = std.caliper, match_args = match_args,
    match_args_dir = match_args_conc,
    forced_exact = forced_exact, forced_antiexact = forced_antiexact
  )
  pairs_conc <- .extract_benefit_pairs(
    .get_matched_data(match_conc, data_conc),
    drug_var, outcome_var, pred_cols, drug_names, keep_vars
  )

  # ---- Direction 2: discordant patients as the treated group (flipped) ----
  if (isTRUE(symmetric)) {
    data_disc <- labelled
    data_disc$cateval_treated <- 1 - data_disc$cateval_concordant
    match_disc <- .run_match(
      data_disc, formula_str, exact_vars, antiexact_vars, caliper, seed,
      method = method, distance = distance, replace = replace, ratio = ratio,
      std.caliper = std.caliper, match_args = match_args,
      match_args_dir = match_args_disc,
      forced_exact = forced_exact, forced_antiexact = forced_antiexact
    )
    pairs_disc <- .extract_benefit_pairs(
      .get_matched_data(match_disc, data_disc),
      drug_var, outcome_var, pred_cols, drug_names, keep_vars
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
    matchit = list(concordant = match_conc, discordant = match_disc),
    drugs = stats::setNames(drug_names, pred_cols)
  )
  class(out) <- "cateval_benefit_match"
  out
}
