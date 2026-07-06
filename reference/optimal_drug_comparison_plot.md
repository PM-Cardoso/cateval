# Plotting optimal drug combinations

This function generates a side-by-side bar plot visualizing:

1.  The proportion of individuals with 1 to 5 optimal drugs identified.

2.  The most frequent drug combinations by category (e.g., 1-drug,
    2-drug).

## Usage

``` r
optimal_drug_comparison_plot(data, groups, plot = TRUE)
```

## Arguments

- data:

  A character vector where each element is a comma-separated string of
  drugs considered optimal for an individual.

- groups:

  A named list where each name represents a category (e.g., "1-drug
  combination") and each value is a vector of drug counts (e.g., 1, 2,
  4:5).

- plot:

  Logical. If TRUE (default), returns a combined patchwork plot. If
  FALSE, returns a list of summary data frames.

## Value

A patchwork object combining two ggplot2 bar plots side-by-side.

## Examples

``` r
groups <- list(
  "1-drug combination" = 1,
  "2-drug combinations" = 2,
  "3-drug combinations" = 3,
  "3/4/5-drug combinations" = 4:5
)
optimal_drug_comparison_plot(c("SGLT2,DPP4,GLP1", "GLP1", "TZD,DPP4"), groups)

```
