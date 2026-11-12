# The pairwise regression shared across calibration groups (cateval 2.6.0),
# its HC3 intervals, the spline option, and HC3 in benefit_citl_slope().

ab_data <- function() {
  d <- cateval_example[cateval_example$treatment %in% c("A", "B"), ]
  d$benefit <- d$pred_A - d$pred_B
  d
}

pair_fit <- function(d, ...) {
  args <- list(d, drug_var = "treatment", drugs = c("A", "B"), benefit_var = "benefit",
               outcome_var = "outcome", cal_groups = 5)
  do.call(pairwise_calibration_pair, utils::modifyList(args, list(...)))
}

adj <- c("age", "sex", "baseline", "marker")


test_that("the effects come from one regression with shared adjustment coefficients", {
  d <- ab_data()
  res <- pair_fit(d, adjustment_var = adj)

  d$grp <- factor(dplyr::ntile(d$benefit, 5))
  d$tA <- as.numeric(d$treatment == "A")
  f <- stats::lm(outcome ~ 0 + grp + grp:tA + age + sex + baseline + marker, data = d)
  expect_equal(res$coef, unname(stats::coef(f)[paste0("grp", 1:5, ":tA")]))
})


test_that("the intervals use HC3 standard errors", {
  skip_if_not_installed("sandwich")
  d <- ab_data()
  res <- pair_fit(d, adjustment_var = adj, conf = 0.9)

  d$grp <- factor(dplyr::ntile(d$benefit, 5))
  d$tA <- as.numeric(d$treatment == "A")
  f <- stats::lm(outcome ~ 0 + grp + grp:tA + age + sex + baseline + marker, data = d)
  nm <- paste0("grp", 1:5, ":tA")
  se <- unname(sqrt(diag(sandwich::vcovHC(f, type = "HC3")))[nm])
  expect_equal(res$coef_low,  res$coef - stats::qnorm(0.95) * se)
  expect_equal(res$coef_high, res$coef + stats::qnorm(0.95) * se)
})


test_that("one group, or no adjustment, gives the estimates of separate regressions", {
  d <- ab_data()
  one <- pair_fit(d, cal_groups = 1, adjustment_var = adj)
  d$tA <- as.numeric(d$treatment == "A")
  expect_equal(one$coef, unname(stats::coef(stats::lm(outcome ~ tA + age + sex + baseline + marker, d))["tA"]))

  unadj <- pair_fit(d)
  d$grp <- dplyr::ntile(d$benefit, 5)
  diffs <- vapply(1:5, function(g) {
    s <- d[d$grp == g, ]
    mean(s$outcome[s$treatment == "A"]) - mean(s$outcome[s$treatment == "B"])
  }, numeric(1))
  expect_equal(unadj$coef, diffs)
})


test_that("output columns and counts are unchanged", {
  d <- ab_data()
  res <- pair_fit(d, adjustment_var = adj, cal_groups = c(3, 5))
  expect_named(res, c("mean", "coef", "coef_low", "coef_high", "n_groups",
                      "drug1", "n_drug1", "drug2", "n_drug2"))
  expect_equal(sum(res$n_drug1[res$n_groups == 5]), sum(d$treatment == "A"))
  expect_equal(sum(res$n_drug2[res$n_groups == 3]), sum(d$treatment == "B"))
})


test_that("a drug with one patient in a group gives NA bounds there only", {
  d <- ab_data()
  d$grp <- dplyr::ntile(d$benefit, 5)
  # relabel all but one patient on A in group 1 as B (the groups stay the same)
  relabel <- which(d$grp == 1 & d$treatment == "A")[-1]
  d$treatment[relabel] <- "B"
  expect_warning(res <- pair_fit(d, adjustment_var = adj), "too few patients")
  g1 <- which(res$n_drug1 == 1)
  expect_length(g1, 1)
  expect_false(is.na(res$coef[g1]))
  expect_true(is.na(res$coef_low[g1]))
  others <- setdiff(seq_len(nrow(res)), g1)
  expect_true(all(is.finite(res$coef_low[others])))
})


test_that("a group with one drug only gets no estimate, without breaking the others", {
  d <- ab_data()
  d$grp <- dplyr::ntile(d$benefit, 5)
  d$treatment[d$grp == 5] <- "A"
  res <- pair_fit(d, adjustment_var = adj)
  expect_true(is.na(res$coef[5]))
  expect_true(all(is.finite(res$coef_low[1:4])))
})


test_that("unused factor levels and categorical variables with one level are handled", {
  d <- ab_data()
  d$sex <- factor(d$sex, levels = c("Female", "Male", "Other"))
  d$site <- "one site"
  res <- pair_fit(d, adjustment_var = c(adj, "site"))
  ref <- pair_fit(ab_data(), adjustment_var = adj)
  expect_equal(res$coef, ref$coef)
})


test_that("rows with a missing adjustment variable are left out of the regression only", {
  d <- ab_data()
  d$age[1:50] <- NA
  res <- pair_fit(d, adjustment_var = adj)
  expect_equal(sum(res$n_drug1) + sum(res$n_drug2), nrow(d))
  expect_true(all(is.finite(res$coef)))
})


test_that("spline_df: NULL and 1 give straight lines; other values change the fit", {
  d <- ab_data()
  straight <- pair_fit(d, adjustment_var = adj)
  expect_equal(pair_fit(d, adjustment_var = adj, spline_df = 1)$coef, straight$coef)
  curved <- pair_fit(d, adjustment_var = adj, spline_df = 3)
  expect_false(isTRUE(all.equal(curved$coef, straight$coef)))
})


test_that("spline_df: one number covers every eligible variable", {
  d <- ab_data()
  d$visits <- sample(1:4, nrow(d), replace = TRUE) # too few distinct values
  all3 <- pair_fit(d, adjustment_var = c(adj, "visits"), spline_df = 3)
  named <- pair_fit(d, adjustment_var = c(adj, "visits"),
                    spline_df = c(age = 3, baseline = 3, marker = 3))
  expect_equal(all3$coef, named$coef)
  expect_equal(all3$coef_low, named$coef_low)
})


test_that("spline_df: a named vector changes only the named variables", {
  d <- ab_data()
  one <- pair_fit(d, adjustment_var = adj, spline_df = c(baseline = 3))
  d$tA <- as.numeric(d$treatment == "A")
  d$grp <- factor(dplyr::ntile(d$benefit, 5))
  f <- stats::lm(outcome ~ 0 + grp + grp:tA + age + sex + splines::ns(baseline, df = 3) + marker, data = d)
  expect_equal(one$coef, unname(stats::coef(f)[paste0("grp", 1:5, ":tA")]))
})


test_that("spline_df: clear errors for wrong input", {
  d <- ab_data()
  d$visits <- sample(1:4, nrow(d), replace = TRUE)
  expect_error(pair_fit(d, spline_df = 3), "needs `adjustment_var`")
  expect_error(pair_fit(d, adjustment_var = adj, spline_df = c(3, 4)), "single number")
  expect_error(pair_fit(d, adjustment_var = adj, spline_df = 2.5), "whole number")
  expect_error(pair_fit(d, adjustment_var = adj, spline_df = 0), "whole number")
  expect_error(pair_fit(d, adjustment_var = adj, spline_df = c(bmi = 3)), "not in `adjustment_var`")
  expect_error(pair_fit(d, adjustment_var = adj, spline_df = c(sex = 3)), "not numeric")
  expect_error(pair_fit(d, adjustment_var = c(adj, "visits"), spline_df = c(visits = 3)),
               "4 distinct values in the A vs B data")
  expect_error(pair_fit(d, adjustment_var = adj, spline_df = c(age = 3, age = 2)), "once")
})


test_that("spline_df: a message when no variable can take a spline", {
  d <- ab_data()
  expect_message(pair_fit(d, adjustment_var = "sex", spline_df = 3), "straight lines")
})


test_that("pairwise_calibration passes spline_df to every pair", {
  d <- ab_data()
  all_pairs <- pairwise_calibration(
    cateval_example, "treatment", c("A", "B", "C"), c("pred_A", "pred_B", "pred_C"),
    "outcome", 5, adjustment_var = adj, spline_df = c(baseline = 3)
  )
  one <- pair_fit(d, adjustment_var = adj, spline_df = c(baseline = 3))
  ab <- all_pairs[all_pairs$drug1 == "A" & all_pairs$drug2 == "B", ]
  rownames(ab) <- NULL
  expect_equal(ab$coef, one$coef)
  expect_error(
    pairwise_calibration(cateval_example, "treatment", c("A", "B"), c("pred_A", "pred_B"),
                         "outcome", 5, spline_df = 3),
    "needs `adjustment_var`"
  )
})


test_that("calls written for 2.5.0, with arguments by position, still run", {
  d <- ab_data()
  res <- pairwise_calibration_pair(d, "treatment", c("A", "B"), "benefit", "outcome", 5,
                                   FALSE, adj)
  expect_equal(res$coef, pair_fit(d, adjustment_var = adj)$coef)
})


test_that("the bootstrap refits the shared regression and keeps the estimates", {
  d <- ab_data()
  boot <- pair_fit(d, adjustment_var = adj, spline_df = 3, bootstrap = TRUE,
                   n_boot = 30, boot_seed = 1)
  wald <- pair_fit(d, adjustment_var = adj, spline_df = 3)
  expect_equal(boot$coef, wald$coef)
  expect_true(all(boot$coef_low < boot$coef & boot$coef < boot$coef_high))
  expect_equal(boot, pair_fit(d, adjustment_var = adj, spline_df = 3, bootstrap = TRUE,
                              n_boot = 30, boot_seed = 1))
})


test_that("benefit_citl_slope uses HC3 standard errors", {
  skip_if_not_installed("sandwich")
  m <- match_benefit_pairs(cateval_example, "treatment", "outcome",
                           pred_cols = c("pred_A", "pred_B", "pred_C"),
                           matching_var = adj, seed = 1)
  cal <- benefit_citl_slope(m$combined, conf = 0.9)
  f <- stats::lm(obs.benefit ~ pred.benefit, data = m$combined)
  se <- unname(sqrt(diag(sandwich::vcovHC(f, type = "HC3"))))
  t_crit <- stats::qt(0.95, f$df.residual)
  expect_equal(cal$lower.ci, unname(stats::coef(f)) - t_crit * se)
  expect_equal(cal$upper.ci, unname(stats::coef(f)) + t_crit * se)
})


test_that("the HC3 helper matches sandwich, also with an aliased column", {
  skip_if_not_installed("sandwich")
  set.seed(2)
  x <- rnorm(60); z <- rnorm(60); y <- 1 + x + rnorm(60, sd = abs(x) + 0.5)
  X <- cbind(1, x, z, 2 * z) # last column aliased
  f <- stats::lm.fit(X, y)
  V <- cateval:::.vcov_hc3(X, f)
  ref <- sandwich::vcovHC(stats::lm(y ~ x + z), type = "HC3")
  expect_equal(unname(V[1:3, 1:3]), unname(ref))
  expect_true(all(is.na(V[4, ])))
})
