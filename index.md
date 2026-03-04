# cateval

**R package for validating Conditional Average Treatment Effects
(CATE)**

## Overview

`cateval` provides a comprehensive suite of functions to validate and
analyze conditional average treatment effects (CATE) from causal
inference studies. Whether you’re comparing two treatments or analyzing
multiple treatment scenarios, this package offers tools for:

- **Estimating and validating** heterogeneous treatment effects
- **Calibrating** predictions for subgroup analysis
- **Identifying optimal treatments** for individual patients or
  populations
- **Assessing treatment benefit** across treatment arms
- **Visualizing** treatment comparisons and outcomes

The package is designed for both **pairwise comparisons** (two
treatments) and **multi-arm studies** (multiple treatments), making it
flexible for various study designs.

## Key Features

### Core Functionality

- **[`calibration_hte()`](https://PM-Cardoso.github.io/cateval/reference/calibration_hte.md)**
  — Assess calibration of heterogeneous treatment effect predictions
- **[`compute_overall_benefit()`](https://PM-Cardoso.github.io/cateval/reference/compute_overall_benefit.md)**
  — Calculate aggregate treatment benefit across populations
- **[`compute_overall_benefit_performance()`](https://PM-Cardoso.github.io/cateval/reference/compute_overall_benefit_performance.md)**
  — Evaluate performance metrics of benefit estimates
- **[`get_best_drugs()`](https://PM-Cardoso.github.io/cateval/reference/get_best_drugs.md)**
  — Identify optimal treatments based on predicted outcomes
- **[`get_ranked_or_tolerant_drugs()`](https://PM-Cardoso.github.io/cateval/reference/get_ranked_or_tolerant_drugs.md)**
  — Rank or select treatments meeting tolerance thresholds
- **[`optimal_drug_comparison_plot()`](https://PM-Cardoso.github.io/cateval/reference/optimal_drug_comparison_plot.md)**
  — Visualize treatment comparisons and performance
- **[`unified_validation()`](https://PM-Cardoso.github.io/cateval/reference/unified_validation.md)**
  — Comprehensive validation framework for treatment effect models

### Built-in Support for

- Continuous outcomes
- Matching-based causal inference (via `MatchIt` integration)
- Publication-ready visualizations with `ggplot2`
- Tidy data workflows with `dplyr` and `tidyr`

## Documentation website

- Full API reference and examples:
  <https://PM-Cardoso.github.io/cateval/>

## Installation

Install the latest development version from GitHub:

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("PM-Cardoso/cateval", dependencies = TRUE, build_vignettes = FALSE)
```

## Quick Start

``` r
library(cateval)

# Validate heterogeneous treatment effects calibration
calibration_results <- calibration_hte(data = your_data, 
                                        treatment = treatment_col,
                                        outcome = outcome_col)

# Identify optimal treatments for your population
best_drugs <- get_best_drugs(predictions = your_predictions,
                              treatment_names = treatment_names)

# Generate comparison visualizations
optimal_drug_comparison_plot(data = comparison_data,
                             treatment = treatment_col)
```

## Citation

If you use `cateval` in your research, please cite:

    @software{cardoso2024cateval,
      author = {Cardoso, Pedro},
      title = {cateval: Validation of Conditional Average Treatment Effects},
      year = {2024},
      url = {https://github.com/PM-Cardoso/cateval}
    }

## Author

**Pedro Cardoso**  
University of Exeter  
Email: <p.cardoso@exeter.ac.uk>  
ORCID: [0000-0002-1014-9058](https://orcid.org/0000-0002-1014-9058)

## License

GPL (≥ 3). See
[LICENSE.md](https://PM-Cardoso.github.io/cateval/LICENSE.md) for
details.

## Contributing

Found a bug or have suggestions? Please open an issue on [GitHub
Issues](https://github.com/PM-Cardoso/cateval/issues).
