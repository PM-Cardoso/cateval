#' Find the Best or Near-Best Treatments for Each Patient
#'
#' Applies \code{\link{get_ranked_or_tolerant_drugs}} across all rows of a dataset to determine the best (or nearly best) drugs based on model predictions. Adds the result as new columns to the data.
#'
#' @param data A data frame containing predicted values in specified columns.
#' @param rank Integer. Rank of the prediction to select (e.g., 1 = best, 2 = second-best). Ignored if `tolerance` is set.
#' @param column_names Character vector. Names of the columns containing predicted values.
#' @param final_var_name Character string. Base name for output columns. Also
#'   the text removed from each name in `column_names` to recover the drug name
#'   (removed literally, once per name), unless `drugs` is supplied.
#' @param tolerance Optional numeric. If set, returns all drugs within this tolerance of the best.
#' @param drugs Optional vector of drug names (character, factor or numeric
#'   codes), one per element of `column_names` and in the same order. Use this when the drug name cannot be
#'   recovered by removing `final_var_name` from the column names; the output
#'   column names are still built from `final_var_name`.
#'
#' @return A modified version of the original data frame with two new columns:
#' \describe{
#'   \item{<label>_drug_value}{Predicted value(s) of the selected drug(s).}
#'   \item{<label>_drug_name}{Name(s) of the selected drug(s).}
#' }
#'
#' @examples
#' # the treatment with the lowest prediction ("pred_" is removed from the
#' # column names to give the treatment names)
#' best <- get_best_drugs(
#'   cateval_example,
#'   column_names   = c("pred_A", "pred_B", "pred_C"),
#'   final_var_name = "pred_"
#' )
#' table(best$pred_rank1_drug_name)
#'
#' # every treatment within 2 points of the best prediction
#' near_best <- get_best_drugs(
#'   cateval_example,
#'   column_names   = c("pred_A", "pred_B", "pred_C"),
#'   final_var_name = "pred_",
#'   tolerance      = 2
#' )
#' head(near_best[, c("pred_A", "pred_B", "pred_C", "pred_within_2_of_best_drug_name")])
#'
#' # treatment names given directly, for any column naming scheme
#' preds <- data.frame(score_x = c(1, 3), score_y = c(2, 1))
#' get_best_drugs(preds, column_names = names(preds), final_var_name = "best_",
#'                drugs = c("X", "Y"))
#'
#' @export
# Main wrapper function to apply drug selection logic row-wise
get_best_drugs <- function(data, rank = 1, column_names = NULL, final_var_name = "", tolerance = NULL, drugs = NULL) {
  # Input validation ----
  
  # Check that column_names are provided and exist in the data
  if (is.null(column_names) || !all(column_names %in% colnames(data))) {
    stop("Invalid or missing column_names")
  }
  if (!is.null(drugs)) .check_drug_names(drugs, column_names)
  
  # Apply selection logic row-wise ----
  
  # Apply get_ranked_or_tolerant_drugs to each row using only the specified columns
  # Results should be a two-column matrix: value (predicted score), name (drug name)
  results <- t(apply(as.data.frame(data)[, column_names, drop = FALSE], 1, get_ranked_or_tolerant_drugs, rank = rank, tolerance = tolerance, column_names = column_names, prediction_column = final_var_name, drugs = drugs))
  
  # Set column names of the results for clarity
  colnames(results) <- c("value", "name")
  
  # Create label for output columns ----
  
  # Generate a label that describes the selection logic used (rank or tolerance)
  label <- if (!is.null(tolerance)) {
    paste0(final_var_name, "within_", tolerance, "_of_best")
  } else {
    paste0(final_var_name, "rank", rank)
  }
  
  # Append new columns to dataset ----
  
  # Add predicted value(s) to the data
  # If using rank, coerce to numeric to ensure consistent column type
  if (is.null(tolerance)) {
    data[[paste0(label, "_drug_value")]] <- as.numeric(results[, "value"])
  } else {
    data[[paste0(label, "_drug_value")]] <- results[, "value"]
  }
  
  # Add selected drug name(s) to the data
  data[[paste0(label, "_drug_name")]] <- results[, "name"]
  
  # Return updated dataset ----
  
  # Return the modified dataset including the new columns
  return(data)
}