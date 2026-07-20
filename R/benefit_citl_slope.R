#' Calibration-in-the-Large and Calibration Slope of Predicted Benefit
#'
#' Regresses observed benefit on predicted benefit across the matched pairs and
#' returns the two calibration summaries: the calibration-in-the-large (the
#' intercept, i.e. the average observed benefit when the predicted benefit is
#' zero) and the calibration slope (how much the observed benefit changes per
#' unit of predicted benefit; 1 indicates perfect calibration).
#'
#' @param pairs A matched-pairs data frame from \code{\link{match_benefit_pairs}}
#'   containing \code{pred.benefit} and \code{obs.benefit}.
#' @param bootstrap Logical. If \code{FALSE} (default) the intervals are the
#'   analytical regression confidence intervals. If \code{TRUE} they are
#'   percentile bootstrap intervals that resample the matched pairs independently
#'   (one row = one pair). These relax the regression's homoscedastic-error
#'   assumption, but do not by themselves account for the dependence from reusing
#'   a control patient across pairs (matching with replacement).
#' @param n_boot Integer. Number of bootstrap replicates.
#' @param boot_seed Optional integer. Seed for reproducibility.
#' @param conf Numeric. Confidence level (default 0.95).
#'
#' @return A data frame with two rows (\code{calibration_in_the_large} and
#'   \code{calibration_slope}) and columns \code{estimate}, \code{lower.ci} and
#'   \code{upper.ci}.
#'
#' @seealso \code{\link{match_benefit_pairs}}, \code{\link{benefit_by_group}}
#'
#' @examples
#' \dontrun{
#' benefit_citl_slope(matched$combined)
#' benefit_citl_slope(matched$combined, bootstrap = TRUE, boot_seed = 42)
#' }
#'
#' @importFrom stats lm coef confint
#' @importFrom boot boot
#' @export
benefit_citl_slope <- function(pairs,
                                bootstrap = FALSE,
                                n_boot = 1000,
                                boot_seed = NULL,
                                conf = 0.95) {

  # ---- Input validation ----
  if (!all(c("pred.benefit", "obs.benefit") %in% colnames(pairs))) {
    stop("`pairs` must contain `pred.benefit` and `obs.benefit` (use match_benefit_pairs()).")
  }

  fit <- stats::lm(obs.benefit ~ pred.benefit, data = pairs)
  estimate <- as.numeric(stats::coef(fit))

  if (isTRUE(bootstrap)) {
    if (!is.null(boot_seed)) set.seed(boot_seed)
    stat_fn <- function(d, idx) stats::coef(stats::lm(obs.benefit ~ pred.benefit, data = d[idx, ]))
    bt <- boot::boot(pairs[, c("obs.benefit", "pred.benefit")], statistic = stat_fn, R = n_boot)
    ci_int <- .percentile_ci(bt$t[, 1], conf)
    ci_slope <- .percentile_ci(bt$t[, 2], conf)
    lower <- c(ci_int[1], ci_slope[1])
    upper <- c(ci_int[2], ci_slope[2])
  } else {
    ci <- stats::confint(fit, level = conf)
    lower <- ci[, 1]
    upper <- ci[, 2]
  }

  data.frame(
    metric = c("calibration_in_the_large", "calibration_slope"),
    estimate = estimate,
    lower.ci = as.numeric(lower),
    upper.ci = as.numeric(upper),
    row.names = NULL
  )
}
