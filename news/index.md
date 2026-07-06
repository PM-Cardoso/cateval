# Changelog

## cateval 1.5.0

### Overall-benefit calibration reworked

- [`match_benefit_pairs()`](https://PM-Cardoso.github.io/cateval/reference/match_benefit_pairs.md)
  replaces `compute_overall_benefit()` as the entry point. It labels
  each patient concordant or discordant, runs symmetric 1:1 matching
  (concordant-as-treated and discordant-as-treated), and returns the
  `concordant`, `discordant` and de-duplicated `combined` matched
  populations together with the underlying `matchit` objects.
- New summary functions take a matched-pairs frame and return the
  calibration numbers:
  [`benefit_deciles()`](https://PM-Cardoso.github.io/cateval/reference/benefit_deciles.md)
  (mean observed vs predicted benefit by group),
  [`benefit_by_drug()`](https://PM-Cardoso.github.io/cateval/reference/benefit_by_drug.md)
  (the same, split by the index/treated patient’s drug),
  [`benefit_overall()`](https://PM-Cardoso.github.io/cateval/reference/benefit_overall.md)
  (mean observed benefit of following the model) and
  [`benefit_calibration()`](https://PM-Cardoso.github.io/cateval/reference/benefit_calibration.md)
  (calibration-in-the-large and calibration slope).
- Each summary offers analytical confidence intervals by default or an
  optional pair-level percentile bootstrap (`bootstrap = TRUE`). Note
  that the bootstrap resamples matched pairs as independent, so it does
  not by itself account for the dependence introduced by matching with
  replacement.
- Predicted benefit is anchored on the index (treated) patient of each
  pair, and observed benefit is computed on the absolute outcome scale.

### Removed

- `compute_overall_benefit()`, `compute_overall_benefit_performance()`
  and `plot_benefit_calibration()` have been removed. Plotting is left
  to the caller; the package returns the underlying numbers only.

### Other

- Added a test suite and continuous integration (`R-CMD-check`).
