# cateval

R package for validating Conditional Average Treatment Effects


### Installation

-   Install latest development version from GitHub (requires [devtools](https://github.com/hadley/devtools) package):

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("PM-Cardoso/cateval", dependencies = TRUE, build_vignettes = FALSE)
```