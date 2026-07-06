# Mean Observed Benefit of Following the Model

Summarises a set of matched pairs into a single number: the average
observed benefit of concordance (receiving the predicted-best drug),
with a confidence interval. This is the headline "how much does
following the model pay off" estimate, and is typically reported for the
pooled (`combined`) matched set.

## Usage

``` r
benefit_overall(
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
  containing `obs.benefit`.

- bootstrap:

  Logical. If `FALSE` (default) the interval is the analytical
  t-interval of the mean. If `TRUE` it is a percentile bootstrap that
  resamples the matched pairs independently (one row = one pair). For a
  mean over many pairs this is usually very close to the analytical
  interval; it does not by itself account for the dependence from
  reusing a control patient across pairs (matching with replacement).

- n_boot:

  Integer. Number of bootstrap replicates.

- boot_seed:

  Optional integer. Seed for reproducibility.

- conf:

  Numeric. Confidence level (default 0.95).

## Value

A one-row data frame with `observed` (mean observed benefit),
`lower.ci`, `upper.ci` and `n` (number of pairs).

## See also

[`match_benefit_pairs`](https://PM-Cardoso.github.io/cateval/reference/match_benefit_pairs.md),
[`benefit_deciles`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md)

## Examples

``` r
if (FALSE) { # \dontrun{
benefit_overall(matched$combined)
benefit_overall(matched$combined, bootstrap = TRUE, boot_seed = 42)
} # }
```
