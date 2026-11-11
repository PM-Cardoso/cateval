#' Group Matched Pairs by Predicted Benefit and Summarise Observed Benefit
#'
#' Divides matched pairs into equal-sized groups (e.g. deciles) of predicted
#' benefit and returns the mean observed benefit in each group, with a
#' confidence interval. This is the data behind the overall-benefit calibration
#' plot: the mean observed benefit in each group is compared with the mean
#' predicted benefit and the line of identity.
#'
#' @param pairs A matched-pairs data frame from \code{\link{match_benefit_pairs}}
#'   (its \code{concordant}, \code{discordant} or \code{combined} element),
#'   containing \code{pred.benefit} and \code{obs.benefit}.
#' @param cal_groups Numeric scalar or vector. Number(s) of groups to split the
#'   predicted benefit into (e.g. \code{10} for deciles, or \code{c(3, 5, 10)}).
#'   Each value is summarised in turn and the results are stacked, with an
#'   \code{n_groups} column recording which setting produced each row, as in
#'   \code{\link{benefit_by_drug}} and \code{\link{pairwise_calibration}}.
#' @param bootstrap Logical. If \code{FALSE} (default) the confidence interval is
#'   the analytical t-interval of the group mean. If \code{TRUE} the interval is a
#'   percentile bootstrap that resamples the matched pairs. The pairs are
#'   resampled independently (one row = one pair), so for a group mean this
#'   interval is usually very close to the analytical t-interval; it does not by
#'   itself widen for the dependence created when a control patient is reused
#'   across several pairs (matching with replacement). It differs from the
#'   analytical interval mainly when a group has few pairs or a skewed
#'   distribution of observed benefit.
#' @param n_boot Integer. Number of bootstrap replicates (used when
#'   \code{bootstrap = TRUE}).
#' @param boot_seed Optional integer. Seed set before bootstrapping for
#'   reproducibility.
#' @param conf Numeric. Confidence level for the interval (default 0.95).
#'
#' @return A data frame with one row per group and \code{cal_groups} setting:
#' \describe{
#'   \item{n_groups}{The \code{cal_groups} value that produced the row.}
#'   \item{group}{Group index, ordered by increasing predicted benefit.}
#'   \item{pred.benefit}{Mean predicted benefit in the group.}
#'   \item{obs.benefit}{Mean observed benefit in the group.}
#'   \item{n}{Number of matched pairs in the group with an observed benefit
#'     (pairs whose \code{obs.benefit} is missing are not counted, since they do
#'     not contribute to the mean or its interval).}
#'   \item{lower.ci, upper.ci}{Confidence-interval bounds for the mean observed
#'     benefit. \code{NA} for a group with fewer than two pairs (see Details).}
#' }
#'
#' @details
#' The group boundaries are fixed from the predicted benefit of \code{pairs}
#' before any bootstrapping, so that each bootstrap replicate re-estimates the
#' mean observed benefit for the *same* group of patients rather than shifting
#' the boundaries between replicates.
#'
#' \strong{Small groups.} A confidence interval needs at least two pairs with
#' an observed benefit. For a group with fewer, the mean is still reported but
#' \code{lower.ci} and \code{upper.ci} are \code{NA} (for both the analytical
#' and the bootstrap interval, since a bootstrap of a single pair collapses to a
#' point). Groups that receive no pairs at all, which can happen with few pairs
#' or heavily tied predictions, are omitted from the output. Both cases are
#' reported in a single warning per \code{cal_groups} value. Reducing
#' \code{cal_groups} avoids them.
#'
#' When several \code{cal_groups} values are given, each is an independent
#' grouping of the same pairs; a plot should show one setting at a time (for
#' example the largest number of groups for which every group still has enough
#' pairs), which is why \code{n_groups} is returned alongside \code{group}.
#'
#' @seealso \code{\link{match_benefit_pairs}}, \code{\link{benefit_by_drug}}
#'
#' @examples
#' matched <- match_benefit_pairs(
#'   cateval_example,
#'   drug_var     = "treatment",
#'   outcome_var  = "outcome",
#'   pred_cols    = c("pred_A", "pred_B", "pred_C"),
#'   matching_var = c("age", "sex", "baseline", "marker"),
#'   seed         = 1
#' )
#'
#' # ten groups of predicted benefit
#' benefit_by_group(matched$combined, cal_groups = 10)
#'
#' # several groupings at once, told apart by n_groups
#' benefit_by_group(matched$combined, cal_groups = c(3, 5))
#'
#' # bootstrap intervals
#' benefit_by_group(matched$combined, cal_groups = 5,
#'                  bootstrap = TRUE, n_boot = 200, boot_seed = 1)
#'
#' @importFrom dplyr group_by summarise mutate n bind_rows select
#' @importFrom stats quantile sd qt
#' @importFrom boot boot
#' @export
benefit_by_group <- function(pairs,
                             cal_groups = 10,
                             bootstrap = FALSE,
                             n_boot = 1000,
                             boot_seed = NULL,
                             conf = 0.95) {

  `%>%` <- dplyr::`%>%`

  # ---- Input validation ----
  if (!all(c("pred.benefit", "obs.benefit") %in% colnames(pairs))) {
    stop("`pairs` must contain `pred.benefit` and `obs.benefit` (use match_benefit_pairs()).")
  }
  if (!is.numeric(cal_groups) || length(cal_groups) < 1 || any(cal_groups < 1)) {
    stop("`cal_groups` must be one or more positive numbers.")
  }

  # Summarise each requested number of groups in turn and stack the results;
  # n_groups records which setting produced each row.
  result <- NULL
  for (cg in cal_groups) {
    one <- .benefit_by_group_one(
      pairs = pairs,
      cal_groups = cg,
      bootstrap = bootstrap,
      n_boot = n_boot,
      boot_seed = boot_seed,
      conf = conf
    )
    one$n_groups <- cg
    result <- dplyr::bind_rows(result, one)
  }

  result %>%
    dplyr::select(n_groups, group, pred.benefit, obs.benefit, n, lower.ci, upper.ci)
}


# Summarise the matched pairs for a single number of groups. Called once per
# cal_groups value by benefit_by_group().
.benefit_by_group_one <- function(pairs, cal_groups, bootstrap, n_boot, boot_seed, conf) {

  `%>%` <- dplyr::`%>%`

  alpha <- (1 - conf) / 2

  # ---- Fixed group boundaries from predicted benefit ----
  # When predicted benefit is heavily tied (e.g. within a small treatment
  # subset) some quantiles coincide, which would give non-unique breakpoints and
  # fewer groups than requested. Keep the distinct boundaries, derive the actual
  # number of groups from them, and warn if fewer groups were formed.
  breaks <- stats::quantile(pairs$pred.benefit, probs = seq(0, 1, length.out = cal_groups + 1), na.rm = TRUE)
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  breaks <- unique(breaks)
  n_groups <- length(breaks) - 1
  if (n_groups < cal_groups) {
    warning(
      "Requested ", cal_groups, " groups but only ", n_groups,
      " could be formed from the distinct predicted-benefit values (ties)."
    )
  }
  pairs$group <- cut(pairs$pred.benefit, breaks = breaks, labels = FALSE, include.lowest = TRUE)

  # ---- Point estimates per group ----
  # Note: `sd` is taken before `obs.benefit` is overwritten by its group mean,
  # because summarise() evaluates its arguments in order.
  # `n` counts pairs with an observed benefit, which is the sample size behind
  # the group mean and its interval.
  summary_tbl <- pairs %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      pred.benefit = mean(pred.benefit, na.rm = TRUE),
      sd = stats::sd(obs.benefit, na.rm = TRUE),
      n = sum(!is.na(obs.benefit)),
      obs.benefit = mean(obs.benefit, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(group))

  # ---- Groups too small for an interval ----
  # Fewer than two pairs: no spread to estimate (analytical), or a degenerate
  # point (bootstrap). Empty groups never appear in summary_tbl.
  too_small <- summary_tbl$n < 2
  empty_groups <- setdiff(seq_len(n_groups), summary_tbl$group)
  if (any(too_small) || length(empty_groups) > 0) {
    msg <- character(0)
    if (any(too_small)) {
      msg <- c(msg, paste0(
        "group(s) ", paste(summary_tbl$group[too_small], collapse = ", "),
        " have fewer than 2 pairs with an observed benefit, so no confidence ",
        "interval can be computed (lower.ci / upper.ci set to NA)"
      ))
    }
    if (length(empty_groups) > 0) {
      msg <- c(msg, paste0(
        "group(s) ", paste(empty_groups, collapse = ", "),
        " contain no pairs and are omitted"
      ))
    }
    warning("cal_groups = ", cal_groups, ": ", paste(msg, collapse = "; "),
            ". Consider fewer groups.", call. = FALSE)
  }

  if (isTRUE(bootstrap)) {

    # Pair-level resampling; the fixed groups keep each bootstrap column tied to
    # the same set of patients.
    if (!is.null(boot_seed)) set.seed(boot_seed)
    stat_fn <- function(d, idx) {
      s <- d[idx, ]
      as.numeric(tapply(s$obs.benefit, factor(s$group, levels = seq_len(n_groups)), mean, na.rm = TRUE))
    }
    bt <- boot::boot(pairs[, c("obs.benefit", "group")], statistic = stat_fn, R = n_boot)
    ci_mat <- apply(bt$t, 2, stats::quantile, probs = c(alpha, 1 - alpha), na.rm = TRUE)

    summary_tbl$lower.ci <- ci_mat[1, summary_tbl$group]
    summary_tbl$upper.ci <- ci_mat[2, summary_tbl$group]
    summary_tbl$lower.ci[too_small] <- NA_real_
    summary_tbl$upper.ci[too_small] <- NA_real_

  } else {

    # Analytical t-interval of the group mean.
    # Computed only where n >= 2, so qt() is never asked for zero degrees of
    # freedom.
    t_crit <- rep(NA_real_, nrow(summary_tbl))
    t_crit[!too_small] <- stats::qt(1 - alpha, summary_tbl$n[!too_small] - 1)
    summary_tbl <- summary_tbl %>%
      dplyr::mutate(
        se = sd / sqrt(n),
        lower.ci = obs.benefit - t_crit * se,
        upper.ci = obs.benefit + t_crit * se
      )
  }

  summary_tbl %>%
    dplyr::select(group, pred.benefit, obs.benefit, n, lower.ci, upper.ci)
}
