# cateval 2.1.1

Documentation-only fixes; no change to any function's behaviour or output.

* `pairwise_calibration()`: the `@return` documentation listed `cal_groups` and
  `grouping` columns that the function does not produce. It now documents the
  columns actually returned (`mean`, `coef`, `coef_low`, `coef_high`, `n_groups`,
  `drug1`, `n_drug1`, `drug2`, `n_drug2`). Note that the output has no per-group
  identifier column; rows are ordered by group within each `n_groups` block.
* `benefit_citl_slope()`: the `@return` documentation omitted the `metric` column
  that identifies each row (`calibration_in_the_large` / `calibration_slope`); it
  is now documented alongside `estimate`, `lower.ci` and `upper.ci`.

# cateval 2.1.0

All changes below are additive: every new argument defaults to the behaviour of
the previous release, so existing code produces the same results.

## Configurable matching

* `match_benefit_pairs()`, `pairwise_calibration()` and
  `pairwise_calibration_pair()` now expose the matching options that were
  previously fixed: `method`, `distance`, `replace` and `ratio`, together with a
  `match_args` list that is applied last and can supply any other
  `MatchIt::matchit()` argument.
* `caliper` and `std.caliper` are now available in the pairwise functions, which
  previously had no caliper support.
* `seed` is now available in the pairwise functions, so a matched cohort can be
  reproduced exactly.

## Matching is performed once per drug pair

* In the pairwise functions, matching previously happened inside the
  `cal_groups` loop and was therefore repeated for each requested number of
  groups. It is now performed once per drug pair and the same matched cohort is
  reused for every value of `cal_groups`, so the different groupings are directly
  comparable. Unseeded results may differ slightly from earlier versions.

## Inspecting the matched cohort

* `return_matched = TRUE` makes the pairwise functions return a list with
  `results`, `matched` and `matchit`, the latter two holding one entry per drug
  pair. This allows the matching to be checked (for example with a balance plot),
  as was already possible for the overall-benefit workflow.

## Bootstrap confidence intervals for the pairwise analysis

* `bootstrap`, `n_boot`, `boot_seed` and `conf` are new in the pairwise
  functions. With `bootstrap = TRUE` the intervals are percentile bootstrap
  intervals instead of the model's Wald intervals; the point estimates are
  unchanged. Calibration groups are fixed before resampling, and matched sets are
  resampled as a unit.
* `conf` also controls the level of the Wald intervals (default 0.95, as before).

# cateval 2.0.0

## Breaking: functions renamed for clarity and consistency

Several exported functions were renamed so the names describe what they do and
follow a consistent scheme (`benefit_*` for the overall-benefit family,
`pairwise_*` for the drug-pair family, `plot_*` for plotting). **Update any code
that calls the old names.**

| Old name | New name | Why |
|----------|----------|-----|
| `benefit_deciles()` | `benefit_by_group()` | It splits into `cal_groups` bins (any number), not only deciles; now parallels `benefit_by_drug()`. |
| `benefit_calibration()` | `benefit_citl_slope()` | Names the two metrics it returns (calibration-in-the-large and slope) and removes the clash with the calibration-curve data. |
| `unified_validation()` | `pairwise_calibration()` | Describes what it does — pairwise benefit calibration across all drug pairs. |
| `calibration_hte()` | `pairwise_calibration_pair()` | Reads as the single-pair engine behind `pairwise_calibration()`. |
| `optimal_drug_comparison_plot()` | `plot_optimal_drugs()` | Clearer, and groups plotting functions under a `plot_` prefix. |

No arguments, return values or behaviour changed — only the function names.

# cateval 1.6.0

## `benefit_by_drug()` supports multiple calibration groupings

* `benefit_by_drug()` now accepts a vector of `cal_groups` values (e.g.
  `c(3, 5, 10)`) as well as a single value. Each requested number of groups is
  computed for every drug class and the results are stacked into one data frame.
* The output gains an `n_groups` column recording which `cal_groups` setting
  produced each row, matching the convention already used by
  `calibration_hte()` and `unified_validation()`. The identifying columns
  (`index_drug`, `n_groups`) now lead the returned data frame.
* Existing single-value calls are unchanged apart from the added `n_groups`
  column.

# cateval 1.5.0

## Overall-benefit calibration reworked

* `match_benefit_pairs()` replaces `compute_overall_benefit()` as the entry point.
  It labels each patient concordant or discordant, runs symmetric 1:1 matching
  (concordant-as-treated and discordant-as-treated), and returns the
  `concordant`, `discordant` and de-duplicated `combined` matched populations
  together with the underlying `matchit` objects.
* New summary functions take a matched-pairs frame and return the calibration
  numbers: `benefit_deciles()` (mean observed vs predicted benefit by group),
  `benefit_by_drug()` (the same, split by the index/treated patient's drug),
  `benefit_overall()` (mean observed benefit of following the model) and
  `benefit_calibration()` (calibration-in-the-large and calibration slope).
* Each summary offers analytical confidence intervals by default or an optional
  pair-level percentile bootstrap (`bootstrap = TRUE`). Note that the bootstrap
  resamples matched pairs as independent, so it does not by itself account for
  the dependence introduced by matching with replacement.
* Predicted benefit is anchored on the index (treated) patient of each pair, and
  observed benefit is computed on the absolute outcome scale.

## Removed

* `compute_overall_benefit()`, `compute_overall_benefit_performance()` and
  `plot_benefit_calibration()` have been removed. Plotting is left to the caller;
  the package returns the underlying numbers only.

## Other

* Added a test suite and continuous integration (`R-CMD-check`).
