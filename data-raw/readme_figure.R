# Builds man/figures/README-calibration.png, the figure on the home page.
# Run from the package root after installing or loading cateval:
#   source("data-raw/readme_figure.R")

library(cateval)
library(ggplot2)

matched <- match_benefit_pairs(
  cateval_example,
  drug_var     = "treatment",
  outcome_var  = "outcome",
  pred_cols    = c("pred_A", "pred_B", "pred_C"),
  matching_var = c("age", "sex", "baseline", "marker"),
  seed         = 1
)
by_group <- benefit_by_group(matched$combined, cal_groups = 10)
lims <- range(0, by_group$pred.benefit, by_group$lower.ci, by_group$upper.ci)

p <- ggplot(by_group, aes(x = pred.benefit, y = obs.benefit)) +
  geom_abline(intercept = 0, slope = 1, colour = "#898781") +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0, colour = "#2a78d6") +
  geom_point(size = 2.5, colour = "#2a78d6") +
  coord_equal(xlim = lims, ylim = lims) +
  labs(
    x = "Predicted benefit of the recommended treatment",
    y = "Observed benefit (95% CI)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
    axis.text  = element_text(colour = "#52514e"),
    axis.title = element_text(colour = "#52514e"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("man/figures/README-calibration.png", p, width = 5, height = 4.6, dpi = 150,
       device = ragg::agg_png)
