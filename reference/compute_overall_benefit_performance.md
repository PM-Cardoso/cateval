# Estimate Observed vs Predicted Treatment Benefit Calibration

This function evaluates how well predicted treatment benefits align with
observed outcome differences. It groups patients by predicted benefit
and compares outcomes between patients who received treatments
concordant or discordant with the predicted optimal treatment.
Concordance can be defined strictly by the top predicted treatment or
flexibly within a specified tolerance of the best predicted outcome.

## Usage

``` r
compute_overall_benefit_performance(
  data,
  drug_var,
  outcome_var,
  pred_cols = NULL,
  conc_tolerance = NULL,
  matching_var = NULL,
  match.exact = NULL,
  match.antiexact = NULL
)
```

## Arguments

- data:

  A data frame containing individual-level data including treatment
  assignments, outcomes, and predicted treatment benefits.

- drug_var:

  Character. Column name indicating the treatment actually received by
  each patient.

- outcome_var:

  Character. Column name for the outcome variable to assess treatment
  effect.

- pred_cols:

  Character vector. Column names of predicted outcomes or risks for each
  treatment. These should share a common prefix (e.g., "pred_GLP1",
  "pred_SGLT2").

- conc_tolerance:

  Optional numeric scalar. If provided, patients are considered
  concordant if they received any treatment within this absolute
  difference of the best predicted outcome. If NULL (default),
  concordance requires receiving the top predicted treatment exactly.

- matching_var:

  Character vector. Names of covariates used for Mahalanobis distance
  matching between concordant and discordant patients.

- match.exact:

  Optional character vector. Variables for exact matching. The best
  predicted treatment is always included automatically.

- match.antiexact:

  Optional character vector. Variables for anti-exact matching (i.e.,
  variables on which matches must differ). The actual treatment variable
  is always included automatically.

- cal_groups:

  Integer. Number of calibration groups to stratify patients based on
  predicted benefit (currently not directly used but reserved for future
  enhancements).

## Value

A named list with two elements:

- calibration_intercept:

  A one-row data frame with the intercept of the linear regression model
  (`calibration_obs ~ calibration_pred`), representing the average
  deviation from perfect calibration when the predicted difference is
  zero. Includes:

- calibration_slope:

  A one-row data frame with the slope of the regression model,
  representing how well predicted differences in treatment effect
  correspond to observed differences. Includes:

Each estimate reflects model fit across matched patient pairs grouped by
predicted benefit.

## Details

To control for confounding, the function performs nearest-neighbor
matching based on specified covariates before estimating observed
treatment effect differences within matched pairs. The observed
differences are then regressed against predicted differences to assess
calibration of predicted treatment benefits.

This function uses nearest-neighbor matching (via the MatchIt package)
on a set of covariates to create matched groups of concordant and
discordant patients. Within each calibration group defined by predicted
benefit, the observed difference in outcome is regressed on the
predicted difference to assess calibration.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- compute_overall_benefit_performance(
  data = data_example,
  drug_var = "treatment",
  outcome_var = "outcome",
  pred_cols = c("pred_GLP1", "pred_SGLT2", "pred_DPP4"),
  conc_tolerance = 0.05,
  matching_var = c("age", "sex"),
  match.exact = NULL,
  match.antiexact = NULL
)
print(result)
} # }
```
