# Estimate Heterogeneous Treatment Effects Across Calibration Groups

Estimates heterogeneous treatment effects between two drugs by
stratifying patients into calibration groups based on predicted benefit
scores. Optionally performs covariate matching before estimating
treatment effects.

## Usage

``` r
calibration_hte(
  data,
  drug_var,
  drugs,
  benefit_var,
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

  A data frame containing observed outcomes, treatment assignments, and
  predicted benefits.

- drug_var:

  Character string. Column name for the treatment assignment variable
  (e.g., drug).

- drugs:

  Character vector of length 2. Names of the two drugs to compare.

- benefit_var:

  Character string. Column name containing predicted benefit scores.

- outcome_var:

  Character string. Column name for the outcome variable (e.g., clinical
  measurement).

- cal_groups:

  Numeric or numeric vector. Number(s) of calibration groups (e.g.,
  quantiles) to divide the data into based on predicted benefit.

- matching:

  Logical. Whether to perform covariate matching using the `MatchIt`
  package before estimating treatment effects.

- adjustment_var:

  Optional character vector. Names of covariates to include as
  adjustment variables in the regression model.

- matching_var:

  Optional character vector. Covariates to use for matching. Defaults to
  `adjustment_var` if not specified.

- match.exact:

  Optional character vector. Variables for exact matching. Matching on
  best predicted drug automatically added.

- match.antiexact:

  Optional character vector. Variables for anti-exact matching.
  `drug_var` automatically added.

## Value

A data frame where each row corresponds to a calibration group within
each `cal_groups` setting. Columns include:

- mean:

  The mean predicted benefit score for patients in the calibration
  group.

- coef:

  Estimated average treatment effect (regression coefficient) comparing
  the two drugs, adjusted if covariates are specified.

- coef_low:

  Lower bound of the 95% confidence interval for the treatment effect.

- coef_high:

  Upper bound of the 95% confidence interval for the treatment effect.

- n_groups:

  Number of calibration groups (i.e., the value of cal_groups) used to
  create this stratification.

- drug1:

  Name of the first drug in the comparison (from `drugs[1]`).

- n_drug1:

  Number of patients receiving drug1 within the calibration group.

- drug2:

  Name of the second drug in the comparison (from `drugs[2]`).

- n_drug2:

  Number of patients receiving drug2 within the calibration group.

## Details

If multiple values are provided to `cal_groups`, the function will
repeat the process for each value and return a combined result.

## Examples

``` r
# Basic usage without matching or adjustment
result <- calibration_hte(
  data = test_data,
  drug_var = "drugclass",
  drugs = c("SGLT2", "DPP4"),
  benefit_var = "benefit_score",
  outcome_var = "posthba1cfinal",
  cal_groups = 5
)
#> Error: object 'test_data' not found

# With adjustment
result_adj <- calibration_hte(
  data = test_data,
  drug_var = "drugclass",
  drugs = c("SGLT2", "DPP4"),
  benefit_var = "benefit_score",
  outcome_var = "posthba1cfinal",
  cal_groups = 3,
  adjustment_var = c("age", "bmi", "sex")
)
#> Error: object 'test_data' not found

# With matching and exact matching
result_match <- calibration_hte(
  data = test_data,
  drug_var = "drugclass",
  drugs = c("SGLT2", "DPP4"),
  benefit_var = "benefit_score",
  outcome_var = "posthba1cfinal",
  cal_groups = 4,
  matching = TRUE,
  adjustment_var = c("age", "bmi", "sex"),
  match.exact = "sex"
)
#> Error: object 'test_data' not found
```
