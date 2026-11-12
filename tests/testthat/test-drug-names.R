# Linking prediction columns to drug names. Any column naming scheme must work,
# either inferred from the drugs in drug_var or given explicitly, and every
# scheme must give exactly the same pairs as the standard "pred.<drug>" names.

rename_preds <- function(d, new_names) {
  names(d)[match(paste0("pred.", test_drugs), names(d))] <- new_names
  d
}

standard_pairs <- function(d, ...) {
  match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
    matching_var = test_matching_var, seed = 1, ...
  )
}

expect_same_pairs <- function(a, b) {
  expect_equal(a$combined$pred.benefit, b$combined$pred.benefit)
  expect_equal(a$combined$obs.benefit, b$combined$obs.benefit)
  expect_equal(a$combined$pair_key, b$combined$pair_key)
}


test_that("drug names are inferred for any prefix or suffix", {
  d <- make_test_data()
  ref <- standard_pairs(d)
  expect_equal(unname(ref$drugs), test_drugs)
  expect_named(ref$drugs, paste0("pred.", test_drugs))

  for (cols in list(paste0("pred_hba1c_", test_drugs),   # longer prefix
                    paste0("p", test_drugs),             # no separator
                    paste0(test_drugs, "_pred"))) {      # suffix
    r <- match_benefit_pairs(
      rename_preds(d, cols), "drugclass", "posthba1cfinal", cols,
      matching_var = test_matching_var, seed = 1
    )
    expect_equal(unname(r$drugs), test_drugs)
    expect_same_pairs(r, ref)
  }
})


test_that("drugs whose names share leading characters are not truncated", {
  # "pred.SGLT2" / "pred.SU" share the characters "pred.S"; the drug names
  # must still be SGLT2 and SU.
  d <- make_test_data()
  d <- d[d$drugclass %in% c("SGLT2", "SU"), ]

  r <- match_benefit_pairs(
    d, "drugclass", "posthba1cfinal", c("pred.SGLT2", "pred.SU"),
    matching_var = c("agetx", "prehba1c"), seed = 1
  )
  expect_equal(unname(r$drugs), c("SGLT2", "SU"))
  expect_gt(nrow(r$combined), 0)
  expect_true(all(r$combined$concordant_drug %in% c("SGLT2", "SU")))
})


test_that("`drugs` or a named `pred_cols` covers names that cannot be inferred", {
  d <- make_test_data()
  ref <- standard_pairs(d)
  cols <- paste0("hba1c_", tolower(test_drugs), "_m12")
  d2 <- rename_preds(d, cols)

  expect_error(
    match_benefit_pairs(d2, "drugclass", "posthba1cfinal", cols,
                        matching_var = test_matching_var),
    "Supply `drugs`"
  )

  via_drugs <- match_benefit_pairs(
    d2, "drugclass", "posthba1cfinal", cols,
    matching_var = test_matching_var, seed = 1, drugs = test_drugs
  )
  via_names <- match_benefit_pairs(
    d2, "drugclass", "posthba1cfinal", stats::setNames(cols, test_drugs),
    matching_var = test_matching_var, seed = 1
  )
  expect_same_pairs(via_drugs, ref)
  expect_same_pairs(via_names, ref)
  expect_equal(via_drugs$drugs, stats::setNames(test_drugs, cols))
})


test_that("`drugs` is validated", {
  d <- make_test_data()
  pc <- paste0("pred.", test_drugs)
  expect_error(standard_pairs(d, drugs = test_drugs[1:4]), "one drug name per")
  expect_error(standard_pairs(d, drugs = c("SGLT2", "SGLT2", "DPP4", "SU", "TZD")), "duplicated")
  expect_error(standard_pairs(d, drugs = c("SGLT2", "GLP1", "DPP4", "SU", NA)), "missing or empty")
})


test_that("drug names with spaces, hyphens and brackets work, including with a tolerance", {
  d <- make_test_data()
  relabel <- c(SGLT2 = "SGLT2i", GLP1 = "GLP-1 RA", DPP4 = "DPP-4i", SU = "SU", TZD = "TZD (pio)")
  e <- d
  e$drugclass <- unname(relabel[e$drugclass])
  e <- rename_preds(e, paste0("pred.", relabel))

  for (tol in list(NULL, 3)) {
    ref <- standard_pairs(d, conc_tolerance = tol)
    r <- match_benefit_pairs(
      e, "drugclass", "posthba1cfinal", paste0("pred.", relabel),
      matching_var = test_matching_var, seed = 1, conc_tolerance = tol
    )
    expect_same_pairs(r, ref)
    expect_true(all(r$combined$concordant_drug %in% relabel))
  }
})


test_that("patients on a drug without predictions are excluded with a warning", {
  d <- make_test_data()
  d$drugclass[1:5] <- "INSULIN"
  d$drugclass[6] <- NA

  expect_warning(
    r <- standard_pairs(d),
    "6 patient\\(s\\) excluded.*INSULIN"
  )
  expect_false(anyNA(r$combined$pred.benefit))
})


test_that("get_best_drugs strips the prefix literally and accepts `drugs`", {
  # "." in the prefix must not act as a regular-expression wildcard.
  preds <- data.frame(predXA = c(1, 3), predXB = c(2, 1))
  out <- get_best_drugs(preds, column_names = names(preds), final_var_name = "pred.")
  expect_equal(out$pred.rank1_drug_name, c("predXA", "predXB"))

  out <- get_best_drugs(preds, column_names = names(preds), final_var_name = "best_",
                        drugs = c("A", "B"))
  expect_equal(out$best_rank1_drug_name, c("A", "B"))

  out <- get_best_drugs(preds, column_names = names(preds), final_var_name = "tol_",
                        drugs = c("A", "B"), tolerance = 1)
  expect_equal(out$tol_within_1_of_best_drug_name, c("A,B", "B"))
})


# The package is general: drug identifiers are whatever the user's data uses.
# Each scheme below relabels the same cohort, so the matched pairs (and the
# pairwise results) must be identical to the standard run.

relabel_drugs <- function(d, new_drugs, cols, as = identity) {
  map <- stats::setNames(new_drugs, test_drugs)
  d$drugclass <- as(unname(map[d$drugclass]))
  names(d)[match(paste0("pred.", test_drugs), names(d))] <- cols
  d
}

test_that("arbitrary drug identifiers give identical results", {
  d <- make_test_data()
  ref <- standard_pairs(d)
  ref_pw <- pairwise_calibration(
    d, "drugclass", test_drugs, paste0("pred.", test_drugs), "posthba1cfinal", 3
  )

  schemes <- list(
    molecules = list(c("sertraline", "fluoxetine", "citalopram", "escitalopram", "venlafaxine"),
                     function(x) paste0("yhat_", x)),
    nested    = list(c("SU", "PSU", "A", "BA", "CBA"), function(x) paste0("pred.", x)),
    numeric   = list(c(1, 2, 3, 10, 11), function(x) paste0("pred", x), as.numeric),
    factor    = list(c("Drug A", "Drug B", "Drug C", "Drug D", "Drug E"),
                     function(x) paste0("pred.", x), factor),
    dotted    = list(c("tx.1", "tx.2", "tx.3", "tx.4", "tx.5"), function(x) paste0("pred.", x)),
    unicode   = list(c("café", "naïve", "Ω-3", "β-blocker", "α"),
                     function(x) paste0("p_", x)),
    suffix    = list(c("arm1", "arm2", "arm3", "arm4", "arm5"), function(x) paste0(x, "_predicted"))
  )

  for (nm in names(schemes)) {
    sc <- schemes[[nm]]
    new_drugs <- sc[[1]]
    cols <- sc[[2]](new_drugs)
    e <- relabel_drugs(d, new_drugs, cols, if (length(sc) > 2) sc[[3]] else identity)

    r <- match_benefit_pairs(e, "drugclass", "posthba1cfinal", cols,
                             matching_var = test_matching_var, seed = 1)
    expect_equal(unname(r$drugs), as.character(new_drugs), info = nm)
    expect_same_pairs(r, ref)

    pw <- pairwise_calibration(e, "drugclass", new_drugs, cols, "posthba1cfinal", 3)
    expect_equal(pw$coef, ref_pw$coef, info = nm)
    expect_equal(as.character(pw$drug1), as.character(new_drugs)[match(ref_pw$drug1, test_drugs)], info = nm)

    byd <- suppressWarnings(benefit_by_drug(r$combined, drugs = new_drugs, cal_groups = 3))
    expect_setequal(levels(byd$index_drug), as.character(new_drugs))
  }
})


test_that("numeric drug codes can be supplied through `drugs`", {
  d <- make_test_data()
  codes <- c(1, 2, 3, 10, 11)
  cols <- paste0("score_", letters[1:5])
  e <- relabel_drugs(d, codes, cols, as.numeric)

  r <- match_benefit_pairs(e, "drugclass", "posthba1cfinal", cols,
                           matching_var = test_matching_var, seed = 1, drugs = codes)
  expect_same_pairs(r, standard_pairs(d))
  expect_equal(unname(r$drugs), as.character(codes))
})


test_that("drug names that differ only in case are matched, with a message", {
  d <- make_test_data()
  ref <- standard_pairs(d)

  # columns in lower case, drug_var in upper case
  cols <- paste0("pred.", tolower(test_drugs))
  e <- relabel_drugs(d, test_drugs, cols)
  expect_message(
    r <- match_benefit_pairs(e, "drugclass", "posthba1cfinal", cols,
                             matching_var = test_matching_var, seed = 1),
    "ignoring case"
  )
  expect_equal(unname(r$drugs), test_drugs)
  expect_same_pairs(r, ref)

  # one level spelled differently in drug_var: kept as spelled there, and no
  # patient is excluded
  e <- d
  e$drugclass[e$drugclass == "SU"] <- "su"
  expect_no_warning(suppressMessages(
    r <- match_benefit_pairs(e, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
                             matching_var = test_matching_var, seed = 1)
  ))
  expect_equal(unname(r$drugs), c("SGLT2", "GLP1", "DPP4", "su", "TZD"))

  # case-insensitive matching is not used when it would be ambiguous: drug_var
  # holds both "SGLT2" and "sglt2"
  cols <- paste0("pred.", tolower(test_drugs))
  e <- relabel_drugs(d, test_drugs, cols)
  e$drugclass[1:20] <- "sglt2"
  expect_error(
    suppressWarnings(match_benefit_pairs(e, "drugclass", "posthba1cfinal", cols,
                                         matching_var = test_matching_var)),
    "Supply `drugs`"
  )
})


test_that("a drug nobody received is read from the column pattern, with a message", {
  d <- make_test_data()
  d <- d[d$drugclass != "TZD", ]

  expect_message(
    r <- match_benefit_pairs(d, "drugclass", "posthba1cfinal", paste0("pred.", test_drugs),
                             matching_var = test_matching_var, seed = 1),
    "pred.<drug>.*TZD"
  )
  expect_equal(unname(r$drugs), test_drugs)
  # TZD can still be the predicted-best drug, so it appears among the
  # alternatives that a discordant patient's predictions are compared with
  expect_false(anyNA(r$combined$pred.benefit))

  # suffix-style columns too
  cols <- paste0(test_drugs, "_hat")
  e <- relabel_drugs(d, test_drugs, cols)
  expect_message(
    r <- match_benefit_pairs(e, "drugclass", "posthba1cfinal", cols,
                             matching_var = test_matching_var, seed = 1),
    "<drug>_hat"
  )
  expect_equal(unname(r$drugs), test_drugs)
})


test_that("plot_optimal_drugs uses neutral titles and keeps its summary tables", {
  sets <- c("A", "A,B", "B", "A,B,C", "C", "A,C")
  groups <- list("1 treatment" = 1, "2 treatments" = 2, "3 treatments" = 3)

  p <- plot_optimal_drugs(sets, groups)
  expect_s3_class(p, "patchwork")
  # long titles wrap; compare the text without line breaks
  expect_equal(gsub("\n", " ", p[[1]]$labels$title), "Number of near-best treatments")

  p2 <- plot_optimal_drugs(sets, groups, title_count = "Custom", subtitle = "Sub")
  expect_equal(p2[[1]]$labels$title, "Custom")
  expect_equal(p2[[1]]$labels$subtitle, "Sub")

  tabs <- plot_optimal_drugs(sets, groups, plot = FALSE)
  expect_named(tabs, c("overall", "breakdown"))
  expect_equal(sum(tabs$overall$Count), length(sets))
})
