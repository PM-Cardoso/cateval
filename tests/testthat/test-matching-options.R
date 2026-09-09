# Tests for the configurable matching options, the returned matched cohort and
# the pairwise bootstrap. The defaults must reproduce the previous behaviour, so
# several tests assert that explicitly passing the defaults changes nothing.

pairwise_args <- function(d, ...) {
  args <- list(
    data = d,
    drug_var = "drugclass",
    drugs = c("SGLT2", "GLP1"),
    prediction_vars = c("pred.SGLT2", "pred.GLP1"),
    outcome_var = "posthba1cfinal",
    cal_groups = 3
  )
  utils::modifyList(args, list(...))
}


test_that("explicitly passing the matching defaults reproduces the defaults", {
  d <- make_test_data()

  implicit <- do.call(pairwise_calibration, pairwise_args(d, matching = TRUE, matching_var = test_matching_var, seed = 1))
  explicit <- do.call(pairwise_calibration, pairwise_args(
    d, matching = TRUE, matching_var = test_matching_var, seed = 1,
    method = "nearest", distance = "mahalanobis", replace = TRUE, ratio = 1
  ))

  expect_equal(implicit, explicit)
})


test_that("pairwise matching is reproducible with a seed", {
  d <- make_test_data()

  a <- do.call(pairwise_calibration, pairwise_args(d, matching = TRUE, matching_var = test_matching_var, seed = 42))
  b <- do.call(pairwise_calibration, pairwise_args(d, matching = TRUE, matching_var = test_matching_var, seed = 42))

  expect_equal(a, b)
})


test_that("pairwise accepts a caliper", {
  d <- make_test_data()

  res <- do.call(pairwise_calibration, pairwise_args(
    d, matching = TRUE, matching_var = test_matching_var, seed = 1,
    caliper = c(prehba1c = 0.5)
  ))

  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0)
})


test_that("one matched cohort is shared across several cal_groups values", {
  d <- make_test_data()

  res <- do.call(pairwise_calibration, pairwise_args(
    d, cal_groups = c(2, 4), matching = TRUE, matching_var = test_matching_var,
    seed = 1, return_matched = TRUE
  ))

  # A single matched cohort per drug pair, not one per cal_groups value.
  expect_named(res, c("results", "matched", "matchit"))
  expect_length(res$matched, 1L)
  expect_setequal(unique(res$results$n_groups), c(2, 4))
})


test_that("return_matched gives back the matched cohort and matchit object", {
  d <- make_test_data()

  res <- do.call(pairwise_calibration, pairwise_args(
    d, matching = TRUE, matching_var = test_matching_var, seed = 1, return_matched = TRUE
  ))

  expect_named(res$matched, "SGLT2_vs_GLP1")
  expect_s3_class(res$matched$SGLT2_vs_GLP1, "data.frame")
  expect_s3_class(res$matchit$SGLT2_vs_GLP1, "matchit")
  # The matchit object is what allows balance to be inspected.
  expect_s3_class(summary(res$matchit$SGLT2_vs_GLP1), "summary.matchit")
})


test_that("return_matched without matching warns and returns NULL cohorts", {
  d <- make_test_data()

  expect_warning(
    res <- do.call(pairwise_calibration, pairwise_args(d, matching = FALSE, return_matched = TRUE)),
    "no matched cohorts"
  )
  expect_null(res$matched)
  expect_null(res$matchit)
  expect_s3_class(res$results, "data.frame")
})


test_that("pairwise bootstrap runs and returns usable intervals", {
  d <- make_test_data()

  res <- do.call(pairwise_calibration, pairwise_args(
    d, bootstrap = TRUE, n_boot = 25, boot_seed = 1
  ))

  estimated <- !is.na(res$coef)
  expect_true(any(estimated))
  expect_true(all(res$coef_low[estimated] <= res$coef_high[estimated]))
})


test_that("pairwise bootstrap is reproducible and differs from the Wald interval", {
  d <- make_test_data()

  boot1 <- do.call(pairwise_calibration, pairwise_args(d, bootstrap = TRUE, n_boot = 25, boot_seed = 7))
  boot2 <- do.call(pairwise_calibration, pairwise_args(d, bootstrap = TRUE, n_boot = 25, boot_seed = 7))
  wald  <- do.call(pairwise_calibration, pairwise_args(d))

  expect_equal(boot1, boot2)
  # Point estimates are unchanged; only the interval method differs.
  expect_equal(boot1$coef, wald$coef)
  expect_false(isTRUE(all.equal(boot1$coef_low, wald$coef_low)))
})


test_that("conf controls the width of the pairwise Wald interval", {
  d <- make_test_data()

  narrow <- do.call(pairwise_calibration, pairwise_args(d, conf = 0.80))
  wide   <- do.call(pairwise_calibration, pairwise_args(d, conf = 0.99))

  estimated <- !is.na(narrow$coef)
  expect_true(all((wide$coef_high - wide$coef_low)[estimated] >
                    (narrow$coef_high - narrow$coef_low)[estimated]))
})


test_that("match_benefit_pairs accepts the new matching options", {
  d <- make_test_data()

  # replace = FALSE warns when a stratum runs out of controls; that is expected.
  m <- suppressWarnings(match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
    matching_var = test_matching_var, seed = 1,
    replace = FALSE, ratio = 1
  ))

  expect_s3_class(m, "cateval_benefit_match")
  expect_true(nrow(m$combined) > 0)
})


test_that("match_args reaches matchit and overrides the named arguments", {
  d <- make_test_data()

  named <- suppressWarnings(match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
    matching_var = test_matching_var, seed = 1, replace = FALSE
  ))
  # replace = TRUE is deliberately overridden by match_args, so both runs match
  # without replacement and must agree.
  via_args <- suppressWarnings(match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
    matching_var = test_matching_var, seed = 1, replace = TRUE,
    match_args = list(replace = FALSE)
  ))

  expect_equal(named$combined, via_args$combined)
})


# Per-direction matching arguments -----------------------------------------
# match_args applies to both directions; match_args_conc / match_args_disc
# override it for one. Empty lists must reproduce the previous behaviour.

benefit_args <- function(d, ...) {
  args <- list(
    data = d,
    drug_var = "drugclass",
    outcome_var = "posthba1cfinal",
    pred_cols = paste0("pred.", test_drugs),
    matching_var = test_matching_var,
    seed = 1
  )
  utils::modifyList(args, list(...))
}


test_that("empty per-direction lists reproduce the previous behaviour", {
  d <- make_test_data()

  implicit <- do.call(match_benefit_pairs, benefit_args(d))
  explicit <- do.call(match_benefit_pairs, benefit_args(
    d, match_args_conc = list(), match_args_disc = list()
  ))

  expect_equal(implicit$concordant, explicit$concordant)
  expect_equal(implicit$discordant, explicit$discordant)
  expect_equal(implicit$combined, explicit$combined)
})


test_that("match_args_conc affects only the concordant direction", {
  d <- make_test_data()

  base <- do.call(match_benefit_pairs, benefit_args(d))
  conc <- suppressWarnings(do.call(match_benefit_pairs, benefit_args(
    d, match_args_conc = list(replace = FALSE)
  )))

  expect_false(identical(base$concordant, conc$concordant))
  expect_equal(base$discordant, conc$discordant)
})


test_that("match_args_disc affects only the discordant direction", {
  d <- make_test_data()

  base <- do.call(match_benefit_pairs, benefit_args(d))
  disc <- suppressWarnings(do.call(match_benefit_pairs, benefit_args(
    d, match_args_disc = list(replace = FALSE)
  )))

  expect_equal(base$concordant, disc$concordant)
  expect_false(identical(base$discordant, disc$discordant))
})


test_that("a direction list overrides the shared list", {
  d <- make_test_data()

  # shared replace = TRUE, overridden to FALSE for the discordant direction only
  layered <- suppressWarnings(do.call(match_benefit_pairs, benefit_args(
    d,
    match_args = list(replace = TRUE),
    match_args_disc = list(replace = FALSE)
  )))
  # the same thing said directly
  direct <- suppressWarnings(do.call(match_benefit_pairs, benefit_args(
    d, match_args_disc = list(replace = FALSE)
  )))

  expect_equal(layered$concordant, direct$concordant)
  expect_equal(layered$discordant, direct$discordant)
})


test_that("NULL in a direction list unsets that argument for one direction", {
  d <- make_test_data()

  # caliper set for both, then removed for the discordant direction; that
  # direction must therefore agree with a run that never had a caliper
  unset <- do.call(match_benefit_pairs, benefit_args(
    d,
    caliper = c(prehba1c = 0.4),
    match_args_disc = list(caliper = NULL)
  ))
  no_caliper <- do.call(match_benefit_pairs, benefit_args(d))
  with_caliper <- do.call(match_benefit_pairs, benefit_args(
    d, caliper = c(prehba1c = 0.4)
  ))

  expect_equal(unset$discordant, no_caliper$discordant)
  expect_equal(unset$concordant, with_caliper$concordant)
})


test_that("the forced exact / antiexact constraints survive a match_args override", {
  d <- make_test_data()

  m <- do.call(match_benefit_pairs, benefit_args(
    d, match_args = list(exact = "sex", antiexact = NULL)
  ))

  # anti-exact on the received drug: a pair must have taken two different drugs
  expect_true(all(m$combined$concordant_drug != m$combined$discordant_drug))
  # exact on the predicted-best drug: still honoured, and the caller's own
  # exact variable is added rather than replacing it. matchit() stores `exact`
  # as a terms formula, so the variable names come from all.vars().
  exact_used <- all.vars(m$matchit$concordant$exact)
  expect_true("cateval_best_drug" %in% exact_used)
  expect_true("sex" %in% exact_used)
})


test_that("match_args_disc warns when symmetric = FALSE", {
  d <- make_test_data()

  expect_warning(
    do.call(match_benefit_pairs, benefit_args(
      d, symmetric = FALSE, match_args_disc = list(replace = FALSE)
    )),
    "symmetric"
  )
})


test_that("unnamed matching-argument lists are rejected", {
  d <- make_test_data()

  expect_error(
    do.call(match_benefit_pairs, benefit_args(d, match_args = list("nearest"))),
    "fully named"
  )
  expect_error(
    do.call(match_benefit_pairs, benefit_args(d, match_args_conc = list(1))),
    "fully named"
  )
})
