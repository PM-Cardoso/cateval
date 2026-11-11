#' Find the Best or Near-Best Treatments for One Patient
#'
#' This function identifies the drug(s) with the best (lowest) predicted outcome from a set of model predictions. It can return the nth-best drug or all drugs within a specified tolerance of the best drug.
#'
#' @param row A named vector or one-row data frame of predictions for a single patient.
#' @param rank Integer. Rank of drug to return (e.g., 1 = best, 2 = second-best). Ignored if `tolerance` is provided.
#' @param column_names Character vector. Names of the columns containing predicted values.
#' @param prediction_column Character string. Text to remove from each name in
#'   `column_names` to recover the drug name (e.g. `"pred_"` turns `"pred_A"`
#'   into `"A"`). It is removed literally (not as a regular expression) and only
#'   once per name, so it may sit anywhere in the name, e.g. a suffix such as
#'   `"_pred"`. Ignored when `drugs` is supplied.
#' @param tolerance Optional numeric. If provided, return all drugs within this margin above the best predicted value.
#' @param drugs Optional vector of drug names (character, factor or numeric
#'   codes), one per element of `column_names` and in the same order. Use this when the drug name cannot be
#'   recovered by removing a single piece of text from the column names (for
#'   example `c("score_a_month12", "score_b_month12")`). Drug names must not
#'   contain commas, since the tolerance output is comma-separated.
#'
#' @return A character vector with two named elements:
#' \describe{
#'   \item{value}{Comma-separated string of predicted values (or single value if using `rank`).}
#'   \item{name}{Comma-separated string of corresponding drug names (or single name if using `rank`).}
#' }
#'
#' @examples
#' # one patient's predictions (lower is better)
#' row <- c(pred_A = 50.2, pred_B = 53.0, pred_C = 50.9)
#'
#' # the best treatment
#' get_ranked_or_tolerant_drugs(row, rank = 1, column_names = names(row),
#'                              prediction_column = "pred_")
#'
#' # the second best
#' get_ranked_or_tolerant_drugs(row, rank = 2, column_names = names(row),
#'                              prediction_column = "pred_")
#'
#' # every treatment within 1 point of the best
#' get_ranked_or_tolerant_drugs(row, tolerance = 1, column_names = names(row),
#'                              prediction_column = "pred_")
#'
#' @export
#'
# Function to find nth-best drug OR drugs within a tolerance
get_ranked_or_tolerant_drugs <- function(row, rank = 1, column_names = NULL, prediction_column = NULL, tolerance = NULL, drugs = NULL) {

  # Extract Relevant Predictions ----

  # Subset only the prediction values from the specified columns
  values <- row[column_names]

  # Label each value with its drug name ----
  # Either the supplied names (by position), or the column names with the
  # prediction text removed. The removal is literal so that characters such as
  # "." in "pred." are not read as regular-expression wildcards.
  names(values) <- .drug_names_from_columns(column_names, prediction_column, drugs)

  # Sort the predictions in ascending order (lower values are better)
  sorted_values <- sort(values)

  # Tolerance Mode ----

  if (!is.null(tolerance)) {

    # Get the best (lowest) prediction value
    best_val <- sorted_values[1]

    # Select all drugs within the tolerance of the best value
    close_vals <- sorted_values[sorted_values <= best_val + tolerance]
    drug_names <- names(close_vals)

    # Return comma-separated predicted values and drug names, sorted alphabetically by name
    return(c(
      value = paste(close_vals[order(drug_names)], collapse = ","),
      name  = paste(drug_names[order(drug_names)], collapse = ",")
    ))

  } else {

    # Ranking Mode ----

    # If the requested rank is beyond available predictions, return NAs
    if (rank > length(sorted_values)) {
      return(c(value = NA, name = NA))
    }

    # Get the value and drug name at the specified rank
    best_value <- sorted_values[rank]
    best_name <- names(sorted_values)[rank]

    # Return single best value and name
    return(c(value = best_value, name = best_name))
  }
}


# Drug name for each prediction column --------------------------------------
#
# Shared by get_ranked_or_tolerant_drugs() and get_best_drugs(). `drugs`, when
# given, is used as-is (one per column, by position). Otherwise the prediction
# text is removed literally, once, from each column name.
#
# @noRd
.drug_names_from_columns <- function(column_names, prediction_column = NULL, drugs = NULL) {
  if (!is.null(drugs)) {
    .check_drug_names(drugs, column_names)
    return(as.character(drugs))
  }
  if (is.null(prediction_column) || identical(prediction_column, "")) {
    return(column_names)
  }
  sub(prediction_column, "", column_names, fixed = TRUE)
}


# Validate a vector of drug names against its prediction columns ------------
#
# @noRd
.check_drug_names <- function(drugs, column_names, arg_name = "drugs") {
  # Drug identifiers can be text, factor levels or numeric codes; they are
  # compared as text.
  if (!is.atomic(drugs) || is.logical(drugs) || is.null(drugs)) {
    stop("`", arg_name, "` must be a vector of drug names (character, factor or numeric codes).")
  }
  drugs <- as.character(drugs)
  if (length(drugs) != length(column_names)) {
    stop("`", arg_name, "` must have one drug name per prediction column (",
         length(column_names), " columns, ", length(drugs), " names).")
  }
  if (anyNA(drugs) || any(drugs == "")) {
    stop("`", arg_name, "` must not contain missing or empty names.")
  }
  if (anyDuplicated(drugs)) {
    stop("`", arg_name, "` must not contain duplicated names: ",
         paste(unique(drugs[duplicated(drugs)]), collapse = ", "), ".")
  }
  if (any(grepl(",", drugs, fixed = TRUE))) {
    stop("Drug names must not contain commas (used to separate drugs in the tolerance output).")
  }
  invisible(TRUE)
}
