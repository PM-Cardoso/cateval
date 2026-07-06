# Group Matched Pairs by Predicted Benefit and Summarise Observed Benefit

Divides matched pairs into equal-sized groups (e.g. deciles) of
predicted benefit and returns the mean observed benefit in each group,
with a confidence interval. This is the data behind the overall-benefit
calibration plot: the mean observed benefit in each group is compared
with the mean predicted benefit and the line of identity.

## Usage

``` r
benefit_deciles(
  pairs,
  cal_groups = 10,
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
  (its `concordant`, `discordant` or `combined` element), containing
  `pred.benefit` and `obs.benefit`.

- cal_groups:

  Integer. Number of groups to split the predicted benefit into (e.g. 10
  for deciles).

- bootstrap:

  Logical. If `FALSE` (default) the confidence interval is the
  analytical t-interval of the group mean. If `TRUE` the interval is a
  percentile bootstrap that resamples the matched pairs. The pairs are
  resampled independently (one row = one pair), so for a group mean this
  interval is usually very close to the analytical t-interval; it does
  not by itself widen for the dependence created when a control patient
  is reused across several pairs (matching with replacement). It differs
  from the analytical interval mainly when a group has few pairs or a
  skewed distribution of observed benefit.

- n_boot:

  Integer. Number of bootstrap replicates (used when
  `bootstrap = TRUE`).

- boot_seed:

  Optional integer. Seed set before bootstrapping for reproducibility.

- conf:

  Numeric. Confidence level for the interval (default 0.95).

## Value

A data frame with one row per group:

- group:

  Group index, ordered by increasing predicted benefit.

- pred.benefit:

  Mean predicted benefit in the group.

- obs.benefit:

  Mean observed benefit in the group.

- n:

  Number of matched pairs in the group.

- lower.ci, upper.ci:

  Confidence-interval bounds for the mean observed benefit.

## Details

The group boundaries are fixed from the predicted benefit of `pairs`
before any bootstrapping, so that each bootstrap replicate re-estimates
the mean observed benefit for the *same* group of patients rather than
shifting the boundaries between replicates.

## See also

[`match_benefit_pairs`](https://PM-Cardoso.github.io/cateval/reference/match_benefit_pairs.md),
[`benefit_by_drug`](https://PM-Cardoso.github.io/cateval/reference/benefit_by_drug.md)

## Examples

``` r
if (FALSE) { # \dontrun{
benefit_deciles(matched$combined, cal_groups = 10)
benefit_deciles(matched$combined, cal_groups = 10, bootstrap = TRUE, boot_seed = 42)
} # }
```
