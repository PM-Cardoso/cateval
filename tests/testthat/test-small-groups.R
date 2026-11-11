# Groups too small for a confidence interval, and pairs with ratio > 1.

small_pairs <- function(pred, obs) data.frame(pred.benefit = pred, obs.benefit = obs)


test_that("single-pair groups get NA bounds and a clear warning", {
  # five pairs into five groups: one pair per group
  pairs <- small_pairs(pred = 1:5, obs = c(1, 2, 3, 4, 5))

  for (boot in c(FALSE, TRUE)) {
    expect_warning(
      res <- benefit_by_group(pairs, cal_groups = 5, bootstrap = boot, n_boot = 20, boot_seed = 1),
      "group\\(s\\) 1, 2, 3, 4, 5 have fewer than 2 pairs"
    )
    expect_equal(res$n, rep(1, 5))
    expect_equal(res$obs.benefit, 1:5)
    expect_true(all(is.na(res$lower.ci) & is.na(res$upper.ci)))
  }
})


test_that("groups of two or more pairs are unaffected", {
  pairs <- small_pairs(pred = 1:6, obs = c(1, 3, 2, 5, 4, 7))
  for (boot in c(FALSE, TRUE)) {
    expect_no_warning(
      res <- benefit_by_group(pairs, cal_groups = 3, bootstrap = boot, n_boot = 20, boot_seed = 1)
    )
    expect_equal(res$n, c(2, 2, 2))
    expect_true(all(is.finite(res$lower.ci) & is.finite(res$upper.ci)))
  }
})


test_that("empty groups are reported and omitted", {
  # quantile breaks at 0 and 10 leave the top group (10, Inf] with no pairs
  pairs <- small_pairs(pred = c(0, 0, 10, 10), obs = c(1, 2, 3, 4))
  expect_warning(
    res <- benefit_by_group(pairs, cal_groups = 3),
    "group\\(s\\) 3 contain no pairs and are omitted"
  )
  expect_equal(res$group, c(1, 2))
  expect_true(all(is.finite(res$lower.ci)))
})


test_that("tied predictions report both the reduced and the small groups", {
  pairs <- small_pairs(pred = c(1, 1, 1, 1, 2), obs = c(1, 2, 3, 4, 5))
  warned <- character(0)
  res <- withCallingHandlers(
    benefit_by_group(pairs, cal_groups = 4),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("only 2 could be formed", warned)))
  expect_true(any(grepl("group\\(s\\) 2 have fewer than 2 pairs", warned)))
  expect_equal(nrow(res), 2)
})


test_that("n counts only pairs with an observed benefit", {
  pairs <- small_pairs(pred = 1:6, obs = c(1, NA, 3, 4, 5, 6))
  res <- benefit_by_group(pairs, cal_groups = 2)
  expect_equal(res$n, c(2, 3))
  # the interval uses that n
  g1 <- c(1, 3)
  half <- stats::qt(0.975, 1) * stats::sd(g1) / sqrt(2)
  expect_equal(res$lower.ci[1], mean(g1) - half)
})


test_that("benefit_by_drug names the drug class in small-group warnings", {
  pairs <- data.frame(
    index_drug = c(rep("A", 6), "B"),
    pred.benefit = c(1:6, 1),
    obs.benefit = c(1:6, 1)
  )
  expect_warning(
    res <- benefit_by_drug(pairs, cal_groups = 2),
    "index_drug = B: cal_groups = 2"
  )
  expect_true(is.na(res$lower.ci[res$index_drug == "B"]))
  expect_true(all(is.finite(res$lower.ci[res$index_drug == "A"])))
})


test_that("benefit_overall and benefit_citl_slope handle too few pairs", {
  expect_warning(ov <- benefit_overall(small_pairs(1, 2)), "Fewer than 2")
  expect_equal(ov$observed, 2)
  expect_true(is.na(ov$lower.ci) && is.na(ov$upper.ci))

  expect_warning(cs <- benefit_citl_slope(small_pairs(c(1, 2), c(1, 3))), "Too few")
  expect_true(all(is.na(cs$lower.ci)))
  expect_warning(cs <- benefit_citl_slope(small_pairs(c(1, 1, 1), c(1, 3, 2))), "no variation")
})


test_that("pairwise groups too small for a Wald interval warn and give NA bounds", {
  d <- make_test_data()
  d <- d[d$drugclass %in% c("SGLT2", "GLP1"), ][1:12, ]
  expect_warning(
    res <- pairwise_calibration(
      d, "drugclass", c("SGLT2", "GLP1"), c("pred.SGLT2", "pred.GLP1"),
      "posthba1cfinal", 4, adjustment_var = c("agetx", "prebmi")
    ),
    "too few patients"
  )
  estimated_no_ci <- !is.na(res$coef) & is.na(res$coef_low)
  expect_true(any(estimated_no_ci))
})


test_that("ratio > 1 gives one row per treated-control pair", {
  d <- make_test_data()
  m <- match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
    matching_var = test_matching_var, ratio = 2, seed = 1, symmetric = FALSE
  )
  pairs <- m$concordant
  # each treated (concordant) patient appears at most `ratio` times, with
  # distinct controls
  per_treated <- table(pairs$cateval_row_id_conc)
  expect_true(all(per_treated <= 2))
  expect_true(any(per_treated == 2))
  expect_false(anyDuplicated(pairs$pair_key) > 0)
})
