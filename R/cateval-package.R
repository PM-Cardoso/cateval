#' @keywords internal
"_PACKAGE"

# Many functions in this package operate inside dplyr / data-masking pipelines,
# where column names are referenced as bare symbols. R CMD check's static
# analysis cannot tell these apart from undefined global variables, so it is
# declared here that they are intentional. This only silences the check NOTE;
# it does not change any behaviour.
utils::globalVariables(c(
  ".", "Combinations", "Count", "Percentage",
  "cateval_best_drug", "cateval_concordant", "cateval_row_id",
  "cateval_tol_list", "cateval_tol_set", "cateval_treated",
  "concordant_drug", "dataset_benefit", "dataset_drug_var",
  "discordant_drug", "drug", "drug_1_pred", "drug_2_pred",
  "drug_count", "drugs", "group", "index_drug", "lower.ci",
  "obs.benefit", "pair_key", "pred.benefit", "se", "subclass",
  "upper.ci"
))
