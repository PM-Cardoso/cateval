# Calibration-in-the-Large and Calibration Slope of Predicted Benefit

Regresses observed benefit on predicted benefit across the matched pairs
and returns the two calibration summaries: the calibration-in-the-large
(the intercept, i.e. the average observed benefit when the predicted
benefit is zero) and the calibration slope (how much the observed
benefit changes per unit of predicted benefit; 1 indicates perfect
calibration).

## Usage

``` r
benefit_calibration(
  pairs,
  bootstrap = FALSE,
  n_boot = 1000,
  boot_seed = NULL,
  conf = 0.95
)
```

## Arguments

- pairs:

  A matched-pairs data frame from
  [`match_benefit_pairs`](https://PM-Cardoso.github.io/cateval/reference/match_benefit_pairs.md)
  containing `pred.benefit` and `obs.benefit`.

- bootstrap:

  Logical. If `FALSE` (default) the intervals are the analytical
  regression confidence intervals. If `TRUE` they are percentile
  bootstrap intervals that resample the matched pairs independently (one
  row = one pair). These relax the regression's homoscedastic-error
  assumption, but do not by themselves account for the dependence from
  reusing a control patient across pairs (matching with replacement).

- n_boot:

  Integer. Number of bootstrap replicates.

- boot_seed:

  Optional integer. Seed for reproducibility.

- conf:

  Numeric. Confidence level (default 0.95).

## Value

A data frame with two rows (`calibration_in_the_large` and
`calibration_slope`) and columns `estimate`, `lower.ci` and `upper.ci`.

## See also

[`match_benefit_pairs`](https://PM-Cardoso.github.io/cateval/reference/match_benefit_pairs.md),
[`benefit_deciles`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md)

## Examples

``` r
if (FALSE) { # \dontrun{
benefit_calibration(matched$combined)
benefit_calibration(matched$combined, bootstrap = TRUE, boot_seed = 42)
} # }
```
