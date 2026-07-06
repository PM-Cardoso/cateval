# Small synthetic cohort used by the tests: one row per patient, a received
# drug class, predicted 12-month HbA1c on each of five drugs, an observed
# outcome and a few baseline covariates for matching.
make_test_data <- function(n = 500, seed = 1) {
  set.seed(seed)
  drugs <- c("SGLT2", "GLP1", "DPP4", "SU", "TZD")
  d <- data.frame(
    drugclass = sample(drugs, n, replace = TRUE),
    agetx = rnorm(n, 60, 10),
    prebmi = rnorm(n, 30, 5),
    prehba1c = rnorm(n, 75, 10),
    sex = sample(c("Male", "Female"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  for (dr in drugs) d[[paste0("pred.", dr)]] <- rnorm(n, 60, 6)
  d$posthba1cfinal <- 55 + rnorm(n, 0, 8)
  d
}

test_drugs <- c("SGLT2", "GLP1", "DPP4", "SU", "TZD")
test_matching_var <- c("agetx", "prebmi", "prehba1c", "sex")
