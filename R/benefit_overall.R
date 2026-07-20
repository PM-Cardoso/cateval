#' Mean Observed Benefit of Following the Model
#'
#' Summarises a set of matched pairs into a single number: the average observed
#' benefit of concordance (receiving the predicted-best drug), with a confidence
#' interval. This is the headline "how much does following the model pay off"
#' estimate, and is typically reported for the pooled (\code{combined}) matched set.
#'
#' @param pairs A matched-pairs data frame from \code{\link{match_benefit_pairs}}
#'   containing \code{obs.benefit}.
#' @param bootstrap Logical. If \code{FALSE} (default) the interval is the
#'   analytical t-interval of the mean. If \code{TRUE} it is a percentile
#'   bootstrap that resamples the matched pairs independently (one row = one
#'   pair). For a mean over many pairs this is usually very close to the
#'   analytical interval; it does not by itself account for the dependence from
#'   reusing a control patient across pairs (matching with replacement).
#' @param n_boot Integer. Number of bootstrap replicates.
#' @param boot_seed Optional integer. Seed for reproducibility.
#' @param conf Numeric. Confidence level (default 0.95).
#'
#' @return A one-row data frame with \code{observed} (mean observed benefit),
#'   \code{lower.ci}, \code{upper.ci} and \code{n} (number of pairs).
#'
#' @seealso \code{\link{match_benefit_pairs}}, \code{\link{benefit_by_group}}
#'
#' @examples
#' \dontrun{
#' benefit_overall(matched$combined)
#' benefit_overall(matched$combined, bootstrap = TRUE, boot_seed = 42)
#' }
#'
#' @importFrom stats sd qt quantile
#' @importFrom boot boot
#' @export
benefit_overall <- function(pairs,
                            bootstrap = FALSE,
                            n_boot = 1000,
                            boot_seed = NULL,
                            conf = 0.95) {

  # ---- Input validation ----
  if (!("obs.benefit" %in% colnames(pairs))) {
    stop("`pairs` must contain `obs.benefit` (use match_benefit_pairs()).")
  }

  alpha <- (1 - conf) / 2
  obs <- pairs$obs.benefit
  n <- sum(!is.na(obs))
  observed <- mean(obs, na.rm = TRUE)

  if (isTRUE(bootstrap)) {
    if (!is.null(boot_seed)) set.seed(boot_seed)
    bt <- boot::boot(obs, statistic = function(d, idx) mean(d[idx], na.rm = TRUE), R = n_boot)
    ci <- .percentile_ci(bt$t, conf)
  } else {
    se <- stats::sd(obs, na.rm = TRUE) / sqrt(n)
    ci <- c(observed - stats::qt(1 - alpha, n - 1) * se,
            observed + stats::qt(1 - alpha, n - 1) * se)
  }

  data.frame(
    observed = observed,
    lower.ci = ci[1],
    upper.ci = ci[2],
    n = n,
    row.names = NULL
  )
}
