# Assign Best Drugs Based on Prediction Rankings or Tolerance

Applies
[`get_ranked_or_tolerant_drugs`](https://PM-Cardoso.github.io/cateval/reference/get_ranked_or_tolerant_drugs.md)
across all rows of a dataset to determine the best (or nearly best)
drugs based on model predictions. Adds the result as new columns to the
data.

## Usage

``` r
get_best_drugs(
  data,
  rank = 1,
  column_names = NULL,
  final_var_name = "",
  tolerance = NULL
)
```

## Arguments

- data:

  A data frame containing predicted values in specified columns.

- rank:

  Integer. Rank of the prediction to select (e.g., 1 = best, 2 =
  second-best). Ignored if `tolerance` is set.

- column_names:

  Character vector. Names of the columns containing predicted values.

- final_var_name:

  Character string. Base name for output columns.

- tolerance:

  Optional numeric. If set, returns all drugs within this tolerance of
  the best.

## Value

A modified version of the original data frame with two new columns:

- \_drug_value:

  Predicted value(s) of the selected drug(s).

- \_drug_name:

  Name(s) of the selected drug(s).

## Examples

``` r
if (FALSE) { # \dontrun{
preds <- data.frame(pred_A = c(0.1, 0.3), pred_B = c(0.2, 0.2), pred_C = c(0.3, 0.1))
get_best_drugs(preds, rank = 1, column_names = names(preds), final_var_name = "best")
get_best_drugs(preds, tolerance = 0.05, column_names = names(preds), final_var_name = "tolerant")
} # }
```
