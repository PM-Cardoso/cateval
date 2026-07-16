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
