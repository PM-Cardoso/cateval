# Get Ranked or Tolerance-Based Drug Recommendation

This function identifies the drug(s) with the best (lowest) predicted
outcome from a set of model predictions. It can return the nth-best drug
or all drugs within a specified tolerance of the best drug.

## Usage

``` r
get_ranked_or_tolerant_drugs(
  row,
  rank = 1,
  column_names = NULL,
  prediction_column = NULL,
  tolerance = NULL
)
```

## Arguments

- row:

  A named vector or one-row data frame of predictions for a single
  patient.

- rank:

  Integer. Rank of drug to return (e.g., 1 = best, 2 = second-best).
  Ignored if `tolerance` is provided.

- column_names:

  Character vector. Names of the columns containing predicted values.

- prediction_column:

  Character string. Common prefix used in prediction column names (e.g.,
  "pred\_").

- tolerance:

  Optional numeric. If provided, return all drugs within this margin
  above the best predicted value.

## Value

A character vector with two named elements:

- value:

  Comma-separated string of predicted values (or single value if using
  `rank`).

- name:

  Comma-separated string of corresponding drug names (or single name if
  using `rank`).

## Examples

``` r
if (FALSE) { # \dontrun{
row <- c(pred_A = 0.2, pred_B = 0.5, pred_C = 0.21)
get_ranked_or_tolerant_drugs(row, rank = 1, column_names = names(row), prediction_column = "pred_")
get_ranked_or_tolerant_drugs(row, tolerance = 0.05, column_names = names(row), prediction_column = "pred_")
} # }
```
