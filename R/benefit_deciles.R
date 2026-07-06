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
#' @param cal_groups Integer. Number of groups to split the predicted benefit
#'   into (e.g. 10 for deciles).
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
#' @return A data frame with one row per group:
#' \describe{
#'   \item{group}{Group index, ordered by increasing predicted benefit.}
#'   \item{pred.benefit}{Mean predicted benefit in the group.}
#'   \item{obs.benefit}{Mean observed benefit in the group.}
#'   \item{n}{Number of matched pairs in the group.}
#'   \item{lower.ci, upper.ci}{Confidence-interval bounds for the mean observed
#'     benefit.}
#' }
#'
#' @details
#' The group boundaries are fixed from the predicted benefit of \code{pairs}
#' before any bootstrapping, so that each bootstrap replicate re-estimates the
#' mean observed benefit for the *same* group of patients rather than shifting
#' the boundaries between replicates.
#'
#' @seealso \code{\link{match_benefit_pairs}}, \code{\link{plot_benefit_calibration}}
#'
#' @examples
#' \dontrun{
#' benefit_deciles(matched$combined, cal_groups = 10)
#' benefit_deciles(matched$combined, cal_groups = 10, bootstrap = TRUE, boot_seed = 42)
#' }
#'
#' @importFrom dplyr group_by summarise mutate n
#' @importFrom stats quantile sd qt
#' @importFrom boot boot
#' @export
benefit_deciles <- function(pairs,
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
  if (!is.numeric(cal_groups) || length(cal_groups) != 1 || cal_groups < 1) {
    stop("`cal_groups` must be a single positive number.")
  }

  alpha <- (1 - conf) / 2

  # ---- Fixed group boundaries from predicted benefit ----
  # When predicted benefit is heavily tied (e.g. within a small drug-class
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
  summary_tbl <- pairs %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      pred.benefit = mean(pred.benefit, na.rm = TRUE),
      sd = stats::sd(obs.benefit, na.rm = TRUE),
      n = dplyr::n(),
      obs.benefit = mean(obs.benefit, na.rm = TRUE),
      .groups = "drop"
    )

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

  } else {

    # Analytical t-interval of the group mean.
    summary_tbl <- summary_tbl %>%
      dplyr::mutate(
        se = sd / sqrt(n),
        lower.ci = obs.benefit - stats::qt(1 - alpha, n - 1) * se,
        upper.ci = obs.benefit + stats::qt(1 - alpha, n - 1) * se
      )
  }

  summary_tbl %>%
    dplyr::select(group, pred.benefit, obs.benefit, n, lower.ci, upper.ci)
}
