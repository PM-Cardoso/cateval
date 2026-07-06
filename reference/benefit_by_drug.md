# Overall-Benefit Calibration Split by Drug Class

Runs
[`benefit_deciles`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md)
separately within each drug class, where a matched pair is assigned to
the drug class received by its index (matched) patient. This shows
whether the model's benefit predictions are well calibrated within the
group of pairs anchored on each drug class.

## Usage

``` r
benefit_by_drug(
  pairs,
  drugs = NULL,
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
  containing `index_drug`, `pred.benefit` and `obs.benefit`.

- drugs:

  Optional character vector of all drug classes to report, in the
  desired order. Supply the full set (rather than relying on the drugs
  that happen to appear in `pairs`) so that the output and any plot show
  a consistent set of drug classes even when a drug class never anchors
  a matched pair. Drug classes with no matched pairs cannot be
  summarised and are dropped with a warning; the returned `index_drug`
  column is a factor whose levels follow `drugs`. If `NULL` (default),
  the drug classes present in `pairs` are used.

- cal_groups:

  Integer. Number of groups per drug class (e.g. 10). Passed to
  [`benefit_deciles`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md).

- bootstrap:

  Logical. Use pair-level bootstrap confidence intervals.

- n_boot:

  Integer. Number of bootstrap replicates.

- boot_seed:

  Optional integer. Seed for reproducibility.

- conf:

  Numeric. Confidence level (default 0.95).

## Value

A data frame stacking the
[`benefit_deciles`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md)
output for each drug class, with an added `index_drug` column
identifying the drug class that anchors each group. Suitable for
faceting.

## See also

[`benefit_deciles`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md),
[`match_benefit_pairs`](https://PM-Cardoso.github.io/cateval/reference/match_benefit_pairs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
benefit_by_drug(matched$combined, cal_groups = 10)
} # }
```
