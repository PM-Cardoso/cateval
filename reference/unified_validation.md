# Perform Pairwise Heterogeneous Treatment Effect Calibration Across Multiple Drugs

This function computes heterogeneous treatment effect calibration curves
for all pairwise comparisons between a specified set of drugs. For each
pair, it calculates the predicted benefit (difference between predicted
outcomes), partitions patients into calibration groups based on
predicted benefit quantiles, and estimates treatment effects within each
group.

## Usage

``` r
unified_validation(
  data,
  drug_var,
  drugs,
  prediction_vars,
  outcome_var,
  cal_groups,
  matching = FALSE,
  adjustment_var = NULL,
  matching_var = adjustment_var,
  match.exact = NULL,
  match.antiexact = NULL
)
```

## Arguments

- data:

  A data.frame containing observed outcomes, treatment assignments, and
  predicted outcomes for each drug.

- drug_var:

  Character string specifying the column name of the treatment variable
  (e.g., assigned drug).

- drugs:

  Character vector of two or more drug names to be compared pairwise.

- prediction_vars:

  Character vector of column names corresponding to predicted outcomes
  for each drug. Must be the same length and order as `drugs`.

- outcome_var:

  Character string specifying the column name of the observed outcome
  variable.

- cal_groups:

  Numeric scalar or vector specifying the number(s) of calibration
  groups (quantiles) to partition the population for subgroup analysis.

- matching:

  Logical; if `TRUE`, performs covariate matching using the `MatchIt`
  package before estimating effects.

- adjustment_var:

  Optional character vector of covariate column names to adjust for in
  regression models.

- matching_var:

  Optional character vector specifying variables to match on. Defaults
  to `adjustment_var` if NULL.

- match.exact:

  Optional character vector of variables for exact matching. The best
  predicted drug is always included.

- match.antiexact:

  Optional character vector of variables for anti-exact matching. The
  treatment variable is always included.

## Value

A data.frame with one row per calibration group and drug pair
combination, including:

- cal_groups:

  Number of calibration groups used.

- grouping:

  Calibration group identifier (e.g., 1 through `cal_groups`).

- mean:

  Mean predicted benefit in the calibration group.

- coef:

  Estimated treatment effect (regression coefficient) within the group.

- coef_low:

  Lower bound of the 95% confidence interval for the estimate.

- coef_high:

  Upper bound of the 95% confidence interval.

- drug1:

  Name of the first drug in the pair.

- n_drug1:

  Number of patients in the group receiving the first drug.

- drug2:

  Name of the second drug in the pair.

- n_drug2:

  Number of patients in the group receiving the second drug.

## Details

Internally, the function calls
[`calibration_hte()`](https://PM-Cardoso.github.io/cateval/reference/calibration_hte.md)
for each drug pair to obtain subgroup treatment effect estimates and
confidence intervals.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage comparing 3 drugs with predicted outcomes
results <- unified_validation(
  data = mydata,
  drug_var = "treatment",
  drugs = c("DrugA", "DrugB", "DrugC"),
  prediction_vars = c("pred_drugA", "pred_drugB", "pred_drugC"),
  outcome_var = "outcome",
  cal_groups = c(3, 5),        # test 3 and 5 calibration groups
  adjustment_var = c("age", "sex"),
  matching = TRUE
)
head(results)
} # }
```
