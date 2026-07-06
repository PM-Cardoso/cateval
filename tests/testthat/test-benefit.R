test_that("get_best_drugs labels the rank-1 drug", {
  d <- make_test_data()
  out <- get_best_drugs(
    d, rank = 1,
    column_names = paste0("pred.", test_drugs),
    final_var_name = "pred."
  )
  expect_true("pred.rank1_drug_name" %in% names(out))
  expect_true(all(out$pred.rank1_drug_name %in% test_drugs))
})

test_that("match_benefit_pairs returns the three matched sets and matchit objects", {
  d <- make_test_data()
  m <- match_benefit_pairs(
    data = d, drug_var = "drugclass", outcome_var = "posthba1cfinal",
    pred_cols = paste0("pred.", test_drugs),
    matching_var = test_matching_var,
    seed = 1
  )
  expect_s3_class(m, "cateval_benefit_match")
  expect_true(all(c("concordant", "discordant", "combined", "matchit") %in% names(m)))
  expect_gt(nrow(m$combined), 0)
  expect_true(all(
    c("pred.benefit", "obs.benefit", "index_drug", "concordant_drug", "discordant_drug")
    %in% names(m$combined)
  ))
  expect_s3_class(m$matchit$concordant, "matchit")
  # a matched pair must have taken two different drugs
  expect_true(all(m$combined$concordant_drug != m$combined$discordant_drug))
})

test_that("benefit summaries return the expected shapes", {
  d <- make_test_data()
  m <- match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
    matching_var = test_matching_var, seed = 1
  )

  dec <- benefit_deciles(m$combined, cal_groups = 5)
  expect_true(all(
    c("group", "pred.benefit", "obs.benefit", "n", "lower.ci", "upper.ci") %in% names(dec)
  ))
  expect_lte(nrow(dec), 5)

  ov <- benefit_overall(m$combined)
  expect_equal(nrow(ov), 1L)
  expect_true(all(c("observed", "lower.ci", "upper.ci", "n") %in% names(ov)))

  cal <- benefit_calibration(m$combined)
  expect_equal(nrow(cal), 2L)
  expect_setequal(cal$metric, c("calibration_in_the_large", "calibration_slope"))
})

test_that("benefit_by_drug groups by index_drug with the full drug list as levels", {
  d <- make_test_data()
  m <- match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
    matching_var = test_matching_var, seed = 1
  )
  byd <- suppressWarnings(benefit_by_drug(m$combined, drugs = test_drugs, cal_groups = 3))
  expect_true("index_drug" %in% names(byd))
  expect_true(is.factor(byd$index_drug))
  expect_setequal(levels(byd$index_drug), test_drugs)
})

test_that("bootstrap option runs and returns finite intervals", {
  d <- make_test_data()
  m <- match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
    matching_var = test_matching_var, seed = 1
  )
  ov <- benefit_overall(m$combined, bootstrap = TRUE, n_boot = 50, boot_seed = 1)
  expect_true(is.finite(ov$lower.ci) && is.finite(ov$upper.ci))

  cal <- benefit_calibration(m$combined, bootstrap = TRUE, n_boot = 50, boot_seed = 1)
  expect_true(all(is.finite(cal$lower.ci)) && all(is.finite(cal$upper.ci)))
})

test_that("match_benefit_pairs errors on bad input", {
  d <- make_test_data()
  expect_error(
    match_benefit_pairs(
      d, "not_a_column", "posthba1cfinal", paste0("pred.", test_drugs),
      matching_var = test_matching_var
    ),
    "drug_var"
  )
})
