# cateval

**R package for validating Conditional Average Treatment Effects (CATE)**

## Overview

`cateval` provides a comprehensive suite of functions to validate and analyze conditional average treatment effects (CATE) from causal inference studies. Whether you're comparing two treatments or analyzing multiple treatment scenarios, this package offers tools for:

- **Estimating and validating** heterogeneous treatment effects
- **Calibrating** predictions for subgroup analysis
- **Identifying optimal treatments** for individual patients or populations
- **Assessing treatment benefit** across treatment arms
- **Visualizing** treatment comparisons and outcomes

The package is designed for both **pairwise comparisons** (two treatments) and **multi-arm studies** (multiple treatments), making it flexible for various study designs.

## Key Features

### Overall-benefit calibration

- **`match_benefit_pairs()`** — Match concordant patients (received the recommended drug) to discordant patients, in both directions, returning the concordant, discordant and combined matched populations for inspection
- **`benefit_by_group()`** — Mean observed benefit by group of predicted benefit (the calibration-plot data)
- **`benefit_by_drug()`** — The same, split by recommended drug class
- **`benefit_overall()`** — Mean observed benefit of following the model
- **`benefit_citl_slope()`** — Calibration-in-the-large and calibration slope of predicted benefit

Each summary reports analytical confidence intervals by default, or optional pair-level bootstrap intervals. Plotting is left to the user (see the worked example) so the calibration plot can be built with whatever detail the analysis calls for — this package only produces the underlying group and pair-level numbers.

### Pairwise drug-class calibration

- **`pairwise_calibration()`** — Predicted vs observed treatment effects across all drug pairs
- **`pairwise_calibration_pair()`** — The underlying single-pair calibration engine

Optional covariate matching (`matching = TRUE`) is performed once per drug pair and
reused across every `cal_groups` setting. Set `return_matched = TRUE` to get the
matched cohort and `matchit` object back for balance checking, and
`bootstrap = TRUE` for percentile bootstrap intervals in place of the Wald intervals.

### Customising the matching

Both workflows use `MatchIt` and default to 1:1 nearest-neighbour Mahalanobis matching
with replacement. `method`, `distance`, `replace`, `ratio`, `caliper`, `std.caliper` and
`seed` are all available, and `match_args` passes any other `MatchIt::matchit()` argument
straight through:

``` r
match_benefit_pairs(
  ...,
  replace    = FALSE,                  # each control used at most once
  ratio      = 2,                      # two controls per treated patient
  caliper    = c(prehba1c = 0.4),      # in SD units by default
  match_args = list(m.order = "random")
)
```

### Predicted-best drug utilities

- **`get_best_drugs()`** — Identify the recommended drug (or the near-best set) for each patient
- **`get_ranked_or_tolerant_drugs()`** — Rank or select treatments meeting a tolerance threshold
- **`plot_optimal_drugs()`** — Visualise the distribution of recommended drug combinations

### Built-in Support for

- Continuous outcomes
- Matching-based causal inference (via `MatchIt` integration)
- Publication-ready visualizations with `ggplot2`
- Tidy data workflows with `dplyr` and `tidyr`

## Documentation website

-   Full API reference and examples: <https://PM-Cardoso.github.io/cateval/>

## Installation

Install the latest development version from GitHub:

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("PM-Cardoso/cateval", dependencies = TRUE, build_vignettes = FALSE)
```

## Quick Start

```r
library(cateval)

drugs <- c("SGLT2", "GLP1", "DPP4", "SU", "TZD")

# Pairwise heterogeneous treatment effect calibration across all drug pairs
pairwise_results <- pairwise_calibration(
  data            = your_data,
  drug_var        = "drugclass",
  drugs           = drugs,
  prediction_vars = paste0("pred.", drugs),
  outcome_var     = "posthba1cfinal",
  cal_groups      = c(3, 5, 10),
  adjustment_var  = c("agetx", "sex", "prehba1c")
)

# Overall (concordant vs discordant) benefit calibration
# 1. match, keeping the matched populations for inspection
matched <- match_benefit_pairs(
  data         = your_data,
  drug_var     = "drugclass",
  outcome_var  = "posthba1cfinal",
  pred_cols    = paste0("pred.", drugs),
  matching_var = c("agetx", "sex", "prehba1c"),
  match.exact  = "sex",
  caliper      = c(prehba1c = 0.4),
  seed         = 19840503
)

# 2. summarise (set bootstrap = TRUE for pair-level bootstrap intervals)
by_group <- benefit_by_group(matched$combined, cal_groups = 10)
overall  <- benefit_overall(matched$combined)
slope    <- benefit_citl_slope(matched$combined)

# 3. plot the calibration yourself: group means + identity line + a LOESS
# curve fitted on the group point estimates
library(ggplot2)
ggplot(by_group, aes(x = pred.benefit, y = obs.benefit)) +
  geom_abline(intercept = 0, slope = 1, colour = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE, span = 1, colour = "blue") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2) +
  geom_point(size = 2.5)

# Identify the predicted-optimal treatment(s) for each patient
best_drugs <- get_best_drugs(
  data         = your_data,
  rank         = 1,
  column_names = paste0("pred.", drugs),
  final_var_name = "pred."
)

# Visualise the distribution of predicted-optimal drug combinations
plot_optimal_drugs(
  data   = best_drugs$pred.within_3_of_best_drug_name,
  groups = list("1-drug" = 1, "2-drug" = 2, "3+-drug" = 3:5)
)
```

## Citation

If you use `cateval` in your research, please cite:

```
@software{cardoso2024cateval,
  author = {Cardoso, Pedro},
  title = {cateval: Validation of Conditional Average Treatment Effects},
  year = {2024},
  url = {https://github.com/PM-Cardoso/cateval}
}
```

## Author

**Pedro Cardoso**  
University of Exeter  
Email: p.cardoso@exeter.ac.uk  
ORCID: [0000-0002-1014-9058](https://orcid.org/0000-0002-1014-9058)

## License

GPL (≥ 3). See [LICENSE.md](LICENSE.md) for details.

## Contributing

Found a bug or have suggestions? Please open an issue on [GitHub Issues](https://github.com/PM-Cardoso/cateval/issues).