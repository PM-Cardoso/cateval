#' Overall-Benefit Calibration Split by Drug Class
#'
#' Runs \code{\link{benefit_deciles}} separately within each drug class, where a
#' matched pair is assigned to the drug class received by its index (matched)
#' patient. This shows whether the model's benefit predictions are well
#' calibrated within the group of pairs anchored on each drug class.
#'
#' @param pairs A matched-pairs data frame from \code{\link{match_benefit_pairs}}
#'   containing \code{index_drug}, \code{pred.benefit} and \code{obs.benefit}.
#' @param drugs Optional character vector of all drug classes to report, in the
#'   desired order. Supply the full set (rather than relying on the drugs that
#'   happen to appear in \code{pairs}) so that the output and any plot show a
#'   consistent set of drug classes even when a drug class never anchors a
#'   matched pair. Drug classes with no matched pairs cannot be summarised and
#'   are dropped with a warning; the returned \code{index_drug} column is a
#'   factor whose levels follow \code{drugs}. If \code{NULL} (default), the drug
#'   classes present in \code{pairs} are used.
#' @param cal_groups Numeric scalar or vector. Number(s) of groups per drug class
#'   (e.g. \code{10}, or \code{c(3, 5, 10)}). Each value is passed in turn to
#'   \code{\link{benefit_deciles}} and the results are stacked, with an
#'   \code{n_groups} column recording which setting produced each row.
#' @param bootstrap Logical. Use pair-level bootstrap confidence intervals.
#' @param n_boot Integer. Number of bootstrap replicates.
#' @param boot_seed Optional integer. Seed for reproducibility.
#' @param conf Numeric. Confidence level (default 0.95).
#'
#' @return A data frame stacking the \code{\link{benefit_deciles}} output for each
#'   drug class and each requested \code{cal_groups} value, with an
#'   \code{index_drug} column identifying the drug class that anchors each group
#'   and an \code{n_groups} column giving the number of calibration groups used.
#'   Suitable for faceting.
#'
#' @seealso \code{\link{benefit_deciles}}, \code{\link{match_benefit_pairs}}
#'
#' @examples
#' \dontrun{
#' benefit_by_drug(matched$combined, cal_groups = 10)
#' benefit_by_drug(matched$combined, cal_groups = c(3, 5, 10))
#' }
#'
#' @importFrom dplyr filter bind_rows select
#' @export
benefit_by_drug <- function(pairs,
                            drugs = NULL,
                            cal_groups = 10,
                            bootstrap = FALSE,
                            n_boot = 1000,
                            boot_seed = NULL,
                            conf = 0.95) {

  `%>%` <- dplyr::`%>%`

  # ---- Input validation ----
  if (!("index_drug" %in% colnames(pairs))) {
    stop("`pairs` must contain `index_drug` (use match_benefit_pairs()).")
  }
  if (!is.numeric(cal_groups) || length(cal_groups) < 1 || any(cal_groups < 1)) {
    stop("`cal_groups` must be one or more positive numbers.")
  }

  # Use the supplied full drug list where given, so the set and ordering of drug
  # classes does not depend on which drugs happen to appear in this particular
  # data set. Fall back to the drug classes present in the pairs.
  if (is.null(drugs)) {
    drugs <- sort(unique(pairs$index_drug))
  }

  # A drug class that anchors no matched pairs has nothing to summarise; report
  # it and carry on with the rest.
  present_drugs <- unique(pairs$index_drug)
  drugs_with_pairs <- drugs[drugs %in% present_drugs]
  drugs_missing <- setdiff(drugs, present_drugs)
  if (length(drugs_missing) > 0) {
    warning(
      "No matched pairs for: ", paste(drugs_missing, collapse = ", "),
      ". These drug classes are omitted from the summary."
    )
  }

  # Loop over each requested number of calibration groups and, within that, each
  # drug class. The n_groups column records which cal_groups setting produced each
  # row so all requested settings can be returned in a single stacked data frame.
  result <- NULL
  for (cg in cal_groups) {
    for (drug in drugs_with_pairs) {
      drug_pairs <- dplyr::filter(pairs, index_drug == drug)

      drug_deciles <- benefit_deciles(
        pairs = drug_pairs,
        cal_groups = cg,
        bootstrap = bootstrap,
        n_boot = n_boot,
        boot_seed = boot_seed,
        conf = conf
      )
      drug_deciles$index_drug <- drug
      drug_deciles$n_groups <- cg

      result <- dplyr::bind_rows(result, drug_deciles)
    }
  }

  # Keep the full drug list as factor levels so downstream ordering / faceting is
  # consistent regardless of which drug classes were summarised.
  result$index_drug <- factor(result$index_drug, levels = drugs)

  # Lead with the identifying columns (drug class and grouping setting).
  result %>%
    dplyr::select(index_drug, n_groups, group, pred.benefit, obs.benefit, n,
                  lower.ci, upper.ci)
}
