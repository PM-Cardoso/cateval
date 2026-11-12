# Internal helpers for the overall-benefit calibration workflow.
#
# These functions are not exported. They are shared by match_benefit_pairs(),
# pairwise_calibration_pair() and the benefit_* summary functions so that the
# matching logic, the concordant/discordant labelling and the
# confidence-interval helpers live in one place. Everything here works on the "predicted-best drug" convention:
# for each patient the recommended drug is the one with the lowest predicted
# outcome, and a patient is concordant if they received that drug.


# Work out which drug each prediction column refers to ---------------------
#
# The overall-benefit workflow needs, for every prediction column, the drug name
# exactly as it appears in the received-drug column, so that "received the
# predicted-best drug" can be checked and a pair's predictions looked up. Drug
# names are whatever the user's data contains (classes, molecules, trial arms,
# numeric codes...), so nothing here assumes particular names.
#
# In order of preference:
#   1. `drugs`, supplied by the caller, one per column (by position);
#   2. the names of a named `pred_cols` vector, e.g. c(A = "pred_A");
#   3. inferred from the drugs observed in the received-drug column. Each column
#      is assigned the observed drug its name ends with (any prefix, e.g.
#      "pred.A", "yhat_outcome_A"), or failing that the observed drug its name
#      starts with (any suffix, e.g. "A_pred"); when several observed drugs fit,
#      the longest is taken (so "pred.BA" is BA, not A). If that fails, the same
#      is tried ignoring case, keeping the spelling used in the received-drug
#      column (only when no two observed drugs differ by case alone).
#      A drug that has predictions but that nobody received cannot be found
#      this way; when at least two columns are matched and all share one
#      prefix/suffix pattern, the remaining columns are read with that pattern.
# The two less direct routes (ignoring case, reading a pattern) report the
# mapping they used in a message. If no unambiguous assignment exists the
# caller is asked to supply `drugs` rather than a guess being made.
#
# Inference is anchored on the observed drug names rather than on a common
# prefix of the column names, because a common prefix can run into the drug
# names themselves (e.g. "pred.SGLT2" / "pred.SU" share "pred.S").
#
# @param pred_cols Prediction column names (optionally named by drug).
# @param drugs Optional vector of drug names, one per column.
# @param observed Values of the received-drug column.
# @return A character vector of drug names, one per element of pred_cols.
# @noRd
.resolve_drug_names <- function(pred_cols, drugs = NULL, observed = NULL) {

  if (is.null(drugs) && !is.null(names(pred_cols)) && all(names(pred_cols) != "")) {
    drugs <- names(pred_cols)
  }

  if (!is.null(drugs)) {
    .check_drug_names(drugs, pred_cols)
    return(as.character(drugs))
  }

  observed <- unique(as.character(observed[!is.na(observed)]))
  observed <- observed[observed != ""]
  cols <- unname(pred_cols)

  # Case-insensitive matching is only safe when the observed drugs stay
  # distinct once case is ignored.
  case_modes <- if (anyDuplicated(tolower(observed))) FALSE else c(FALSE, TRUE)

  match_columns <- function(ignore_case, at_end) {
    fold <- if (ignore_case) tolower else identity
    test <- if (at_end) endsWith else startsWith
    vapply(cols, function(col) {
      hits <- observed[test(fold(col), fold(observed))]
      if (length(hits) == 0) NA_character_ else hits[which.max(nchar(hits))]
    }, character(1), USE.NAMES = FALSE)
  }
  modes <- expand.grid(at_end = c(TRUE, FALSE), ignore_case = case_modes)

  # Pass 1: every column matched to an observed drug. All matching variants are
  # tried before any pattern is read, so a drug that *is* in drug_var (perhaps
  # in a different case) is never treated as one nobody received.
  for (i in seq_len(nrow(modes))) {
    candidate <- match_columns(modes$ignore_case[i], modes$at_end[i])
    if (!anyNA(candidate) && !anyDuplicated(candidate)) {
      if (modes$ignore_case[i]) .report_drug_mapping(cols, candidate, "matching drug names ignoring case")
      return(candidate)
    }
  }

  # Pass 2: some columns matched. If they all follow one pattern, read the
  # remaining columns (drugs that nobody in the data received) with it.
  for (i in seq_len(nrow(modes))) {
    at_end <- modes$at_end[i]
    candidate <- match_columns(modes$ignore_case[i], at_end)
    matched <- !is.na(candidate)
    if (sum(matched) < 2 || !any(!matched) || anyDuplicated(candidate[matched])) next

    n_col <- nchar(cols[matched])
    n_drug <- nchar(candidate[matched])
    prefix <- if (at_end) substr(cols[matched], 1, n_col - n_drug) else rep("", sum(matched))
    suffix <- if (at_end) rep("", sum(matched)) else substr(cols[matched], n_drug + 1, n_col)
    if (length(unique(prefix)) != 1 || length(unique(suffix)) != 1) next

    pre <- prefix[1]
    suf <- suffix[1]
    rest <- cols[!matched]
    fits <- startsWith(rest, pre) & endsWith(rest, suf) & nchar(rest) > nchar(pre) + nchar(suf)
    if (!all(fits)) next

    candidate[!matched] <- substr(rest, nchar(pre) + 1, nchar(rest) - nchar(suf))
    if (anyDuplicated(candidate)) next

    .report_drug_mapping(
      cols, candidate,
      paste0("reading the column pattern \"", pre, "<drug>", suf, "\" for drug(s) ",
             "that no patient in `drug_var` received (",
             paste(candidate[!matched], collapse = ", "), ")")
    )
    return(candidate)
  }

  stop(
    "Could not work out which drug each column in `pred_cols` refers to from ",
    "the drugs in `drug_var` (",
    paste(observed, collapse = ", "), "). Supply `drugs`, one name per ",
    "`pred_cols` column in the same order and spelled as in `drug_var`, or ",
    "name the vector by drug, e.g. pred_cols = c(<drug> = \"<column>\", ...).",
    call. = FALSE
  )
}


# Tell the user which drug was taken for each column, when it was inferred
# indirectly. @noRd
.report_drug_mapping <- function(cols, drug_names, how) {
  message(
    "Drug names for `pred_cols` were inferred by ", how, ": ",
    paste0(cols, " -> ", drug_names, collapse = ", "),
    ". Supply `drugs` to set them explicitly."
  )
}


# Label each patient as concordant or discordant ----------------------------
#
# Adds internal columns to the data:
#   cateval_best_drug  - the predicted-best (rank 1) drug name
#   cateval_concordant - 1 if the received drug is the recommended drug
#                        (or, with a tolerance, any drug within the tolerance
#                        of the best), 0 otherwise
#   cateval_tol_*      - only when a tolerance is used: the individual drugs in
#                        the tolerance set, used later for anti-exact matching
# plus the cateval_rank1_* / cateval_within_* columns written by
# get_best_drugs(), which are kept for inspection.
#
# Drug-set membership is tested literally (split on the "," that
# get_ranked_or_tolerant_drugs() uses to join the set), so drug names that
# contain spaces, hyphens or regular-expression characters (e.g. "GLP-1 RA")
# are handled correctly.
#
# @param data Data frame containing the drug column and the prediction columns.
# @param drug_var Name of the received-drug column.
# @param pred_cols Prediction column names.
# @param drug_names Drug name for each element of pred_cols.
# @param conc_tolerance Optional numeric tolerance defining the best-drug set.
# @return A list with the labelled data and the extra anti-exact variables.
# @noRd
.label_concordance <- function(data, drug_var, pred_cols, drug_names, conc_tolerance = NULL) {

  `%>%` <- dplyr::`%>%`

  label_prefix <- "cateval_"
  received <- as.character(data[[drug_var]])

  # Predicted-best (rank 1) drug for every patient.
  labelled <- get_best_drugs(
    data = data,
    rank = 1,
    column_names = pred_cols,
    final_var_name = label_prefix,
    drugs = drug_names
  ) %>%
    dplyr::rename(cateval_best_drug = paste0(label_prefix, "rank1_drug_name"))

  if (is.null(conc_tolerance)) {

    # Strict definition: concordant means received the single recommended drug.
    labelled$cateval_concordant <- as.numeric(received == labelled$cateval_best_drug)

    antiexact_extra <- NULL

  } else {

    # Tolerance definition: concordant means received any drug whose predicted
    # outcome is within conc_tolerance of the best. The individual drugs in the
    # tolerance set are spread across cateval_tol_* columns so that a discordant
    # match can be forced to differ from every drug in the set.
    n_drugs <- length(pred_cols)
    tol_vars <- paste0("cateval_tol_", seq_len(n_drugs))

    labelled <- get_best_drugs(
      data = labelled,
      tolerance = conc_tolerance,
      column_names = pred_cols,
      final_var_name = label_prefix,
      drugs = drug_names
    ) %>%
      dplyr::rename(cateval_tol_set = paste0(label_prefix, "within_", conc_tolerance, "_of_best_drug_name"))

    tol_list <- strsplit(labelled$cateval_tol_set, ",", fixed = TRUE)
    labelled$cateval_concordant <- as.numeric(mapply(function(r, set) r %in% set, received, tol_list))

    # One column per position in the set, recycled to a fixed width. For a
    # discordant patient every column holds the drug they received, so the
    # anti-exact constraint on all columns means a discordant match's drug is
    # not in the concordant patient's set.
    for (k in seq_len(n_drugs)) {
      labelled[[tol_vars[k]]] <- ifelse(
        labelled$cateval_concordant == 0,
        received,
        vapply(tol_list, function(set) rep(set, length.out = n_drugs)[k], character(1))
      )
    }

    antiexact_extra <- tol_vars
  }

  list(data = labelled, antiexact_extra = antiexact_extra)
}


# Build the matching formula from the matching covariates -------------------
#
# Shared by both workflows. Continuous covariates enter directly (Mahalanobis
# distance handles them); categorical covariates are added only when they have
# more than one level in the data. Continuous covariates already handled by
# exact / anti-exact matching are dropped from the distance formula;
# categorical ones are not, since they are constant within an exact stratum and
# so contribute nothing to the distance anyway. With no covariates left the
# right-hand side is "1", which suits exact-only methods.
#
# @param treat_var Name of the 0/1 treatment column the formula is for.
# @return A formula string of the form "treat_var ~ x1 + x2 + ...".
# @noRd
.build_match_formula <- function(data, matching_var, match.exact, match.antiexact,
                                 treat_var = "cateval_treated") {
  categorical_vars <- matching_var[vapply(data[matching_var], function(x) is.factor(x) || is.character(x), logical(1))]
  cont_vars <- setdiff(matching_var, c(categorical_vars, match.exact, match.antiexact))

  terms <- cont_vars
  for (v in categorical_vars) {
    if (length(unique(data[[v]])) > 1) terms <- c(terms, v)
  }
  if (length(terms) == 0) terms <- "1"

  paste(treat_var, "~", paste(terms, collapse = " + "))
}


# Reject variables the method itself defines ---------------------------------
#
# The received drug and the outcome are part of the comparison being validated:
# matching or adjusting on the outcome would condition on the result, and exact
# matching on the received drug contradicts the requirement that the two sides
# of a comparison took different drugs.
#
# @param vars Named list of character vectors to check, e.g.
#   list(matching_var = ..., match.exact = ...).
# @param reserved Named character vector of reserved columns, e.g.
#   c(drug_var = "drugclass", outcome_var = "hba1c").
# @noRd
.check_reserved_vars <- function(vars, reserved) {
  for (arg in names(vars)) {
    clash <- intersect(vars[[arg]], reserved)
    if (length(clash) > 0) {
      which_role <- names(reserved)[match(clash, reserved)]
      stop("`", arg, "` must not include ",
           paste0("`", which_role, "` (", clash, ")", collapse = " or "),
           ", which defines the comparison being validated.")
    }
  }
  invisible(TRUE)
}


# Require a fully named list ------------------------------------------------
#
# modifyList() needs names to match arguments on; an unnamed element otherwise
# fails deep inside matchit() with an opaque message.
#
# @noRd
.check_named_list <- function(x, arg_name) {
  if (length(x) > 0 && (is.null(names(x)) || any(names(x) == ""))) {
    stop("`", arg_name, "` must be a fully named list.")
  }
  invisible(TRUE)
}


# Union a user specification with the variables cateval requires ------------
#
# matchit() accepts `exact` / `antiexact` as either a character vector or a
# one-sided formula, so both are reduced to variable names before the union.
# The caller's variables keep their original order and the required ones are
# appended, so when nothing has been overridden the input is returned unchanged.
#
# @return A character vector of variable names.
# @noRd
.protect_vars <- function(user, forced) {
  if (length(forced) == 0) return(user)
  unique(c(.spec_vars(user), forced))
}


# Variable names in an exact / antiexact specification (vector or formula).
# @noRd
.spec_vars <- function(x) {
  if (is.null(x)) return(character(0))
  if (inherits(x, "formula")) all.vars(x) else as.character(x)
}


# Run one match -------------------------------------------------------------
#
# A single call to MatchIt::matchit(), shared by both workflows. The treated
# group is whichever group the caller has coded 1 in the formula's treatment
# column: concordant or discordant patients in the overall-benefit workflow,
# patients predicted to do better on the first drug in the pairwise workflow.
# Behaviour defaults to 1:1 nearest-neighbour Mahalanobis matching with
# replacement.
#
# Arguments are layered in order -- the named arguments, then match_args
# (shared by both directions), then match_args_dir (this direction only) -- so
# later layers win. The overlays are applied in sequence against the
# accumulated list rather than merged with each other first, because
# modifyList() treats NULL as removal: merging first would discard a NULL in
# match_args_dir instead of using it to unset the argument it names.
#
# `forced_exact` / `forced_antiexact` are then re-applied, so no override can
# drop the constraints that make a comparison interpretable. Supplying `exact`
# or `antiexact` through the overlays therefore adds to them.
#
# @return A MatchIt "matchit" object.
# @noRd
.run_match <- function(data, formula_str, exact_vars, antiexact_vars,
                       caliper, seed,
                       method = "nearest",
                       distance = "mahalanobis",
                       replace = TRUE,
                       ratio = 1,
                       std.caliper = TRUE,
                       match_args = list(),
                       match_args_dir = list(),
                       forced_exact = character(0),
                       forced_antiexact = character(0)) {
  if (!is.null(seed)) set.seed(seed)

  args <- list(
    formula = stats::as.formula(formula_str),
    data = data,
    method = method,
    distance = distance,
    replace = replace,
    ratio = ratio,
    exact = exact_vars,
    antiexact = antiexact_vars,
    caliper = caliper,
    std.caliper = std.caliper
  )

  args <- utils::modifyList(args, match_args)
  args <- utils::modifyList(args, match_args_dir)

  args$exact     <- .protect_vars(args$exact, forced_exact)
  args$antiexact <- .protect_vars(args$antiexact, forced_antiexact)

  # A variable cannot be both matched exactly and forced to differ; matchit()
  # would otherwise return no matches with an unhelpful message.
  both <- intersect(.spec_vars(args$exact), .spec_vars(args$antiexact))
  if (length(both) > 0) {
    stop("Variable(s) ", paste(both, collapse = ", "),
         " cannot be matched both exactly and anti-exactly.")
  }

  do.call(MatchIt::matchit, args)
}


# Matched data, one row per patient per matched set --------------------------
#
# A wrapper around MatchIt::get_matches(). get_matches() adds columns named
# "id", "distance", "weights" and "subclass", and stops if the data already has
# a column with one of those names. "id" and "weights" in particular are common
# column names, so when one is taken the added column gets a cateval_ prefix
# instead; the defaults are kept otherwise, so the matched cohorts look the same
# as before. "subclass" is used by cateval itself to link the patients of a
# matched set, so a user column with that name is reported clearly.
#
# Methods that do not link each treated patient to its own controls (e.g.
# "full", "subclass", "exact", "cem") have no match.matrix; that case gets a
# message in cateval's terms. Any other error is passed on unchanged.
#
# @noRd
.get_matched_data <- function(match, data) {
  if ("subclass" %in% colnames(data)) {
    stop("`data` must not contain a column named \"subclass\"; cateval uses it to ",
         "link the patients of each matched set. Please rename that column.",
         call. = FALSE)
  }
  if (is.null(match$match.matrix)) {
    stop(
      "cateval needs each treated patient linked to its own matched control(s), ",
      "which requires a pair-matching method (e.g. method = \"nearest\", ",
      "\"optimal\" or \"genetic\"). The matching method used (\"", match$info$method,
      "\") does not produce these.",
      call. = FALSE
    )
  }
  free_name <- function(x) if (x %in% colnames(data)) paste0("cateval_", x) else x
  MatchIt::get_matches(
    match, data = data,
    id = free_name("id"),
    distance = free_name("distance"),
    weights = free_name("weights")
  )
}


# Turn a matched data set into one row per matched pair ---------------------
#
# MatchIt::get_matches() returns one matched set ("subclass") per treated
# patient: that treated ("index") patient plus the `ratio` control(s) matched
# to it (fewer if a caliper or a shortage of controls left it short). Within a
# matching direction every treated patient is in the same group -- concordant
# when matching concordant-as-treated, discordant when matching
# discordant-as-treated -- and every control is in the other group, because
# the treatment indicator *is* the concordance label (or its reverse).
#
# Each treated patient is therefore paired with each of its controls, giving
# one row per treated-control pair: one row per treated patient when
# ratio = 1, and up to `ratio` rows (all sharing the same treated patient) when
# ratio > 1. For each row:
#   pred.benefit - the index patient's predicted advantage of the concordant
#                  drug over the discordant drug of the pair, i.e.
#                  index_pred[discordant_drug] - index_pred[concordant_drug]
#                  (positive = model expects the concordant drug to give a
#                  lower / better outcome)
#   obs.benefit  - the observed outcome of the discordant patient minus that of
#                  the concordant patient, on the same scale
#   index_drug   - the drug the index (treated) patient received; used to split
#                  the calibration by drug class
#
# @param matched Output of MatchIt::get_matches(), carrying cateval_concordant,
#   cateval_treated, cateval_row_id, the drug column, the outcome column and the
#   pred columns.
# @param drug_names Drug name for each element of pred_cols.
# @param keep_vars Extra covariate columns to retain for inspection.
# @return A data frame with one row per treated-control pair.
# @noRd
.extract_benefit_pairs <- function(matched, drug_var, outcome_var, pred_cols, drug_names, keep_vars = NULL) {

  `%>%` <- dplyr::`%>%`

  # Guard the structure described above; get_matches() guarantees it for the
  # pair-matching methods, and the pairing below relies on it.
  treated_per_set <- tapply(matched$cateval_treated, matched$subclass, sum)
  if (any(treated_per_set != 1)) {
    stop("Each matched set must contain exactly one treated patient; got sets with ",
         paste(sort(unique(treated_per_set[treated_per_set != 1])), collapse = ", "), ".")
  }

  # Index (treated) patient of each pair. The predicted benefit is evaluated
  # from this patient's predictions, so we keep their full prediction vector,
  # their own drug and any covariates to inspect.
  index <- matched %>%
    dplyr::filter(cateval_treated == 1) %>%
    dplyr::select(
      subclass,
      index_drug = dplyr::all_of(drug_var),
      dplyr::all_of(pred_cols),
      dplyr::all_of(keep_vars)
    )

  # Concordant patient of the pair (received the recommended drug).
  concordant <- matched %>%
    dplyr::filter(cateval_concordant == 1) %>%
    dplyr::select(
      subclass,
      cateval_row_id_conc = cateval_row_id,
      concordant_drug = dplyr::all_of(drug_var),
      outcome_conc = dplyr::all_of(outcome_var)
    )

  # Discordant patient of the pair (received a different drug).
  discordant <- matched %>%
    dplyr::filter(cateval_concordant == 0) %>%
    dplyr::select(
      subclass,
      cateval_row_id_disc = cateval_row_id,
      discordant_drug = dplyr::all_of(drug_var),
      outcome_disc = dplyr::all_of(outcome_var)
    )

  # One side of each join has a single row per subclass (the treated patient),
  # so this yields one row per treated-control pair.
  pairs <- concordant %>%
    dplyr::inner_join(discordant, by = "subclass", relationship = "many-to-many") %>%
    dplyr::inner_join(index, by = "subclass") %>%
    dplyr::filter(as.character(concordant_drug) != as.character(discordant_drug))

  # Predicted benefit from the index patient's predictions (matrix indexing is
  # fast and safe for an arbitrary number of drug classes): the recommended drug
  # is the concordant patient's drug, the alternative is the discordant drug.
  pred_matrix <- as.matrix(pairs[, pred_cols])
  colnames(pred_matrix) <- drug_names
  row_index <- seq_len(nrow(pairs))
  pred_recommended <- pred_matrix[cbind(row_index, match(as.character(pairs$concordant_drug), drug_names))]
  pred_alternative <- pred_matrix[cbind(row_index, match(as.character(pairs$discordant_drug), drug_names))]

  pairs$pred.benefit <- pred_alternative - pred_recommended
  pairs$obs.benefit <- pairs$outcome_disc - pairs$outcome_conc

  # A direction-independent key for a pair of patients, used to drop duplicates
  # when the two matching directions are combined.
  pairs$pair_key <- paste(
    pmin(pairs$cateval_row_id_conc, pairs$cateval_row_id_disc),
    pmax(pairs$cateval_row_id_conc, pairs$cateval_row_id_disc),
    sep = "_"
  )

  pairs
}


# Minimum number of distinct values for a variable to enter as a spline.
.spline_min_distinct <- 10


# Check the form of `spline_df` ----------------------------------------------
#
# NULL (no splines), one whole number (every eligible continuous adjustment
# variable) or a vector of whole numbers named by adjustment variable. The data
# are checked later, in .pairwise_adjustment_design(), once the calibration
# cohort of the pair is known.
#
# @noRd
.check_spline_df <- function(spline_df, adjustment_var) {
  if (is.null(spline_df)) return(invisible(NULL))
  if (length(adjustment_var) == 0) {
    stop("`spline_df` needs `adjustment_var`: splines apply to adjustment variables.")
  }
  if (!is.numeric(spline_df) || length(spline_df) == 0 || anyNA(spline_df) ||
      any(spline_df < 1) || any(spline_df != round(spline_df))) {
    stop("`spline_df` must be NULL, a whole number of at least 1, or a vector of them ",
         "named by adjustment variable (e.g. c(baseline = 3)).")
  }
  nms <- names(spline_df)
  if (is.null(nms)) {
    if (length(spline_df) != 1) {
      stop("`spline_df` must be a single number, or a vector named by adjustment ",
           "variable (e.g. c(baseline = 3)).")
    }
  } else {
    if (any(nms == "") || anyDuplicated(nms)) {
      stop("Every element of `spline_df` must be named, once, by an adjustment variable.")
    }
    unknown <- setdiff(nms, adjustment_var)
    if (length(unknown) > 0) {
      stop("`spline_df` names variable(s) not in `adjustment_var`: ",
           paste(unknown, collapse = ", "), ".")
    }
  }
  invisible(NULL)
}


# Adjustment variables that enter as splines, with their df -------------------
#
# @param cc The adjustment variables on the rows used in the outcome model.
# @return A numeric vector of df, named by variable (empty for no splines).
# @noRd
.spline_terms <- function(cc, adjustment_var, spline_df, pair_label) {
  if (is.null(spline_df)) return(numeric(0))
  n_distinct <- function(x) length(unique(x))

  if (is.null(names(spline_df))) {
    need <- max(.spline_min_distinct, spline_df + 2)
    eligible <- adjustment_var[vapply(adjustment_var, function(v) {
      is.numeric(cc[[v]]) && n_distinct(cc[[v]]) >= need
    }, logical(1))]
    if (length(eligible) == 0) {
      message("`spline_df`: no adjustment variable is numeric with at least ", need,
              " distinct values", pair_label, ", so all enter as straight lines.")
    }
    return(stats::setNames(rep(spline_df, length(eligible)), eligible))
  }

  for (v in names(spline_df)) {
    x <- cc[[v]]
    if (!is.numeric(x)) {
      stop("`spline_df`: `", v, "` is not numeric, so it cannot enter as a spline.")
    }
    need <- max(.spline_min_distinct, spline_df[[v]] + 2)
    if (n_distinct(x) < need) {
      stop("`spline_df`: `", v, "` has ", n_distinct(x), " distinct values", pair_label,
           "; a spline with ", spline_df[[v]], " df needs at least ", need, ".")
    }
  }
  spline_df
}


# Adjustment columns for the pairwise outcome model ---------------------------
#
# Built once for the calibration cohort of a pair, so every number of
# calibration groups and every bootstrap replicate uses the same columns, and
# the same spline knots. Continuous variables enter as they are, or as natural
# cubic splines (splines::ns, knots at quantiles) when `spline_df` says so;
# categorical ones (factor, character, logical) as indicators, reference level
# left out, with unused levels dropped. A categorical variable with one level in
# the cohort adds nothing and is left out.
#
# @return A list: `Z`, a numeric matrix with one row per row of `data` (NA on
#   rows not used), and `complete`, TRUE for the rows the outcome model can use
#   (outcome, drug and every adjustment variable present).
# @noRd
.pairwise_adjustment_design <- function(data, adjustment_var, spline_df, pair_label = "") {
  complete <- !is.na(data$dataset_outcome_var) & !is.na(data$dataset_drug_var)
  if (length(adjustment_var) == 0) {
    return(list(Z = matrix(numeric(0), nrow = nrow(data), ncol = 0), complete = complete))
  }
  complete <- complete & stats::complete.cases(data[, adjustment_var, drop = FALSE])
  cc <- data[complete, adjustment_var, drop = FALSE]

  df_by_var <- .spline_terms(cc, adjustment_var, spline_df, pair_label)

  blocks <- list()
  for (v in adjustment_var) {
    x <- cc[[v]]
    if (v %in% names(df_by_var)) {
      b <- splines::ns(x, df = df_by_var[[v]])
      blocks[[v]] <- matrix(b, nrow = nrow(b),
                            dimnames = list(NULL, paste0(v, "_ns", seq_len(ncol(b)))))
    } else if (is.factor(x) || is.character(x) || is.logical(x)) {
      f <- droplevels(as.factor(x))
      if (nlevels(f) < 2) next
      ind <- vapply(levels(f)[-1], function(l) as.numeric(f == l), numeric(length(f)))
      blocks[[v]] <- matrix(ind, nrow = length(f),
                            dimnames = list(NULL, paste0(v, levels(f)[-1])))
    } else {
      blocks[[v]] <- matrix(as.numeric(x), ncol = 1, dimnames = list(NULL, v))
    }
  }

  Z <- matrix(NA_real_, nrow = nrow(data), ncol = 0)
  if (length(blocks) > 0) {
    zc <- do.call(cbind, unname(blocks))
    Z <- matrix(NA_real_, nrow = nrow(data), ncol = ncol(zc),
                dimnames = list(NULL, colnames(zc)))
    Z[complete, ] <- zc
  }
  list(Z = Z, complete = complete)
}


# Shared outcome model across the calibration groups --------------------------
#
# One least-squares fit with an intercept and a drug effect for each group and
# the adjustment columns `Z` shared by every group:
#   outcome ~ 0 + group + group:drug1 + Z
# A group needs patients on both drugs for its effect to be estimated; groups
# without get NA. Columns that cannot be estimated (aliased) are left out, as
# lm() does; the drug effects come before `Z`, so they are kept first.
#
# @param y Outcome; @param is_drug1 TRUE for drugs[1]; @param group Integer
#   group (1..n_groups); @param Z Adjustment columns. All on the rows used.
# @param se If TRUE, also return HC3 standard errors.
# @return A list: `coef` (effect of drugs[1] minus drugs[2], per group), `n_min`
#   (fewer of the two drugs in each group), `df_residual`, and `se` when asked.
# @noRd
.pairwise_shared_fit <- function(y, is_drug1, group, n_groups, Z, se = TRUE) {
  n1 <- tabulate(group[is_drug1], n_groups)
  n2 <- tabulate(group[!is_drug1], n_groups)
  out <- list(coef = rep(NA_real_, n_groups), n_min = pmin(n1, n2),
              df_residual = NA_real_, se = rep(NA_real_, n_groups))

  effect_groups <- which(n1 > 0 & n2 > 0)
  if (length(effect_groups) == 0) return(out)

  present <- which(n1 + n2 > 0)
  X <- cbind(
    outer(group, present, `==`) * 1,
    outer(group, effect_groups, `==`) * is_drug1,
    Z
  )
  fit <- stats::lm.fit(X, y)
  eff_idx <- length(present) + seq_along(effect_groups)
  out$coef[effect_groups] <- fit$coefficients[eff_idx]
  out$df_residual <- length(y) - fit$rank

  if (isTRUE(se)) {
    V <- .vcov_hc3(X, fit)
    out$se[effect_groups] <- sqrt(diag(V)[eff_idx])
  }
  out
}


# HC3 (heteroscedasticity-consistent) covariance matrix -----------------------
#
# Sandwich estimator (X'X)^-1 X' diag(e_i^2 / (1 - h_i)^2) X (X'X)^-1, as
# sandwich::vcovHC(type = "HC3"). h_i is the leverage of row i. A row with
# leverage 1 fits exactly (residual 0) and makes HC3 undefined; it is given
# weight 0, so one such row does not turn every standard error into NaN.
# Aliased columns (not estimated) get NA.
#
# @param X The model matrix. @param fit The lm() or lm.fit() result for X.
# @return A p x p matrix in the column order of X.
# @noRd
.vcov_hc3 <- function(X, fit) {
  p <- ncol(X)
  r <- fit$rank
  V <- matrix(NA_real_, p, p)
  if (r == 0) return(V)

  keep <- fit$qr$pivot[seq_len(r)]
  xtx_inv <- chol2inv(fit$qr$qr[seq_len(r), seq_len(r), drop = FALSE])
  h <- rowSums(qr.Q(fit$qr)[, seq_len(r), drop = FALSE]^2)
  w <- ifelse(h < 1 - 1e-8, fit$residuals / (1 - h), 0)
  meat <- crossprod(X[, keep, drop = FALSE] * w)

  V[keep, keep] <- xtx_inv %*% meat %*% xtx_inv
  V
}


# Percentile confidence interval from a vector of bootstrap statistics ------
#
# @param t Numeric vector of bootstrap replicates.
# @param conf Confidence level (e.g. 0.95).
# @return A length-2 numeric vector: lower and upper bounds.
# @noRd
.percentile_ci <- function(t, conf = 0.95) {
  alpha <- (1 - conf) / 2
  as.numeric(stats::quantile(t, probs = c(alpha, 1 - alpha), na.rm = TRUE))
}
