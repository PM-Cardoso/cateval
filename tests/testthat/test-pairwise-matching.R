# Pairwise matching: treated = predicted to do better on drugs[1], always
# anti-exact on the received drug, with the user's own column names usable in
# every matching argument.

pw_matched <- function(d, ...) {
  args <- list(
    data = d, drug_var = "drugclass",
    drugs = c("SGLT2", "GLP1", "DPP4"),
    prediction_vars = c("pred.SGLT2", "pred.GLP1", "pred.DPP4"),
    outcome_var = "posthba1cfinal", cal_groups = 3,
    matching = TRUE, matching_var = test_matching_var, seed = 1,
    return_matched = TRUE
  )
  do.call(pairwise_calibration, utils::modifyList(args, list(...)))
}


test_that("every matched set contains patients on different drugs", {
  res <- pw_matched(make_test_data())

  for (pair in names(res$matched)) {
    m <- res$matched[[pair]]
    drugs_per_set <- tapply(as.character(m$drugclass), m$subclass, function(x) length(unique(x)))
    expect_true(all(drugs_per_set == 2), info = pair)
    expect_true("dataset_drug_var" %in% res$matchit[[pair]]$info$antiexact)
  }
})


test_that("the received-drug anti-exact constraint survives match_args overrides", {
  d <- make_test_data()
  base <- pw_matched(d)

  # NULL would remove it if it were not protected
  removed <- pw_matched(d, match_args = list(antiexact = NULL))
  expect_true("dataset_drug_var" %in% removed$matchit[[1]]$info$antiexact)
  expect_equal(removed$results, base$results)

  # a caller's antiexact adds to it
  added <- pw_matched(d, match_args = list(antiexact = "sex"))
  expect_setequal(added$matchit[[1]]$info$antiexact, c("sex", "dataset_drug_var"))

  # naming the drug column itself is accepted and changes nothing
  named <- pw_matched(d, match.antiexact = "drugclass")
  expect_equal(named$results, base$results)
})


test_that("matching arguments can use the caller's own column names", {
  d <- make_test_data()
  res <- pw_matched(d, match.exact = "sex", caliper = c(prehba1c = 0.5))
  expect_s3_class(res$results, "data.frame")
  expect_true("sex" %in% all.vars(res$matchit[[1]]$exact))

  # the matched cohort keeps the original columns as well as the working copies
  m <- res$matched[[1]]
  expect_true(all(c("drugclass", "posthba1cfinal", "dataset_drug_var",
                    "dataset_outcome_var", "dataset_benefit") %in% names(m)))
})


test_that("the drug and outcome columns are rejected as matching / adjustment variables", {
  d <- make_test_data()
  expect_error(pw_matched(d, match.exact = "drugclass"), "match.exact")
  expect_error(pw_matched(d, matching_var = c("agetx", "posthba1cfinal")), "matching_var")
  expect_error(pw_matched(d, adjustment_var = "drugclass", matching_var = "agetx"), "adjustment_var")
  expect_error(pw_matched(d, match_args = list(exact = "drugclass")), "both exactly and anti-exactly")
})


test_that("an unmatchable drug pair warns rather than stopping the run", {
  d <- make_test_data()
  warned <- character(0)
  res <- withCallingHandlers(
    pairwise_calibration(
      d, "drugclass", c("SGLT2", "GLP1", "DPP4"), c("pred.SGLT2", "pred.GLP1", "pred.DPP4"),
      "posthba1cfinal", 3, matching = TRUE, matching_var = "agetx",
      caliper = c(agetx = 1e-9), seed = 1
    ),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # one warning per drug pair, and an empty result rather than an error
  expect_length(grep("No matches could be formed", warned), 3)
  expect_equal(nrow(res), 0)

  # the pair engine itself returns NULL for an unmatchable pair
  expect_warning(
    one <- pairwise_calibration_pair(
      transform(d, b = pred.SGLT2 - pred.GLP1), "drugclass", c("SGLT2", "GLP1"),
      "b", "posthba1cfinal", 3, matching = TRUE, matching_var = "agetx",
      caliper = c(agetx = 1e-9), seed = 1
    ),
    "No matches"
  )
  expect_null(one)
})


test_that("non-pair matching methods give a clear error", {
  d <- make_test_data()
  expect_error(
    suppressWarnings(match_benefit_pairs(
      d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
      matching_var = c("sex"), match.exact = "sex",
      match_args = list(method = "exact", distance = NULL)
    )),
    "pair-matching method"
  )
})
