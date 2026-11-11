#' Plot How Many Treatments Are Near-Best, and Which
#'
#' Summarises the sets of near-best treatments produced by
#' \code{\link{get_best_drugs}} with a \code{tolerance}, and draws two bar
#' charts side by side:
#' \enumerate{
#'   \item the share of patients with 1, 2, 3, ... treatments within the
#'     tolerance of their best predicted outcome;
#'   \item the most common treatment sets, split into the groups given in
#'     \code{groups} (for example single treatments vs pairs).
#' }
#'
#' @param data A character vector with one element per patient: the
#'   comma-separated set of near-best treatments, e.g. \code{"A,C"} (the
#'   \code{*_drug_name} column from \code{get_best_drugs(tolerance = ...)}).
#' @param groups A named list. Each name labels a panel of the second chart
#'   (e.g. \code{"1 treatment"}) and each value gives the set sizes it covers
#'   (e.g. \code{1}, \code{2} or \code{3:5}).
#' @param plot Logical. If \code{TRUE} (default), return the combined plot. If
#'   \code{FALSE}, return the two summary tables instead.
#' @param title_count,title_sets Titles of the two charts.
#' @param subtitle Optional subtitle for the first chart, e.g. to state the
#'   tolerance used.
#'
#' @return If \code{plot = TRUE}, a patchwork object combining the two ggplot2
#'   charts. If \code{plot = FALSE}, a list with \code{overall} (patients by
#'   number of near-best treatments) and \code{breakdown} (patients by
#'   treatment set and group).
#'
#' @importFrom dplyr mutate filter count group_by ungroup arrange rename
#' @importFrom ggplot2 ggplot aes geom_text geom_col geom_hline coord_flip scale_y_continuous labs theme_minimal theme scale_x_continuous
#' @importFrom stringr str_count
#' @importFrom tibble tibble
#' @importFrom patchwork wrap_plots
#' @importFrom purrr map_dfr
#' @importFrom scales percent
#' @importFrom stats reorder
#' @importFrom tidyr drop_na
#'
#' @examples
#' near_best <- get_best_drugs(
#'   cateval_example,
#'   column_names   = c("pred_A", "pred_B", "pred_C"),
#'   final_var_name = "pred_",
#'   tolerance      = 2
#' )
#' sets <- near_best$pred_within_2_of_best_drug_name
#' head(sets)
#'
#' plot_optimal_drugs(
#'   sets,
#'   groups   = list("1 treatment" = 1, "2 treatments" = 2, "3 treatments" = 3),
#'   subtitle = "Treatments within 2 points of each patient's best prediction"
#' )
#'
#' @export
plot_optimal_drugs <- function(data, groups, plot = TRUE,
                               title_count = "Number of near-best treatments",
                               title_sets = "Most common near-best sets",
                               subtitle = NULL) {
  
  # Use pipe locally
  `%>%` <- dplyr::`%>%`
  
  # ---- Summarize number of drugs per individual ----
  summary_number <- data.frame(drugs = data) %>%
    dplyr::mutate(drug_count = ifelse(drugs == "", 0L, stringr::str_count(drugs, ",") + 1)) %>%
    dplyr::filter(drug_count != 0) %>%
    dplyr::count(drug = drug_count, name = "Count") %>%
    dplyr::mutate(Percentage = round(Count / length(data), 4))
  
  # ---- Summarize most common drug combinations by group ----
  summary_combinations <- purrr::map_dfr(
    names(groups),
    function(label) {
      tibble::tibble(drugs = data) %>%
        dplyr::mutate(drug_count = ifelse(drugs == "", 0L, stringr::str_count(drugs, ",") + 1)) %>%
        dplyr::filter(drug_count %in% groups[[label]]) %>%
        dplyr::count(drug = drugs, name = "Count") %>%
        dplyr::mutate(
          Percentage = round(Count / length(data), 4),
          Combinations = label
        )
    }
  ) %>%
    dplyr::group_by(Combinations) %>%
    dplyr::arrange(Combinations, dplyr::desc(Percentage), .by_group = TRUE) %>%
    dplyr::mutate(Group_Percentage = round(sum(Count) / length(data), 4)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(Percentage >= 0.001)
  
  # ---- Shared look ----
  # One accent colour, hairline vertical grid, muted axis text.
  bar_colour <- "#2a78d6"
  theme_bars <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "#52514e"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "#e1e0d9", linewidth = 0.3),
      plot.title = ggplot2::element_text(face = "bold", colour = "#0b0b0b"),
      plot.subtitle = ggplot2::element_text(colour = "#52514e"),
      strip.text = ggplot2::element_text(hjust = 0, face = "bold", colour = "#0b0b0b")
    )
  pct_axis <- ggplot2::scale_x_continuous(
    labels = scales::percent,
    expand = ggplot2::expansion(mult = c(0, 0.25))
  )

  # ---- Plot number of drugs ----
  plot_a <- ggplot2::ggplot(
    summary_number,
    ggplot2::aes(x = Percentage, y = factor(drug, levels = rev(sort(unique(drug)))))
  ) +
    ggplot2::geom_col(fill = bar_colour, width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(Percentage, accuracy = 0.1)),
      hjust = -0.15, size = 3.2, colour = "#52514e"
    ) +
    pct_axis +
    ggplot2::labs(title = title_count, subtitle = subtitle) +
    theme_bars

  # ---- Prepare combination plot data ----
  combi_order <- names(groups)
  df2 <- summary_combinations %>%
    dplyr::mutate(Combinations = factor(Combinations, levels = combi_order)) %>%
    dplyr::group_by(Combinations) %>%
    dplyr::arrange(dplyr::desc(Percentage), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(drug = factor(drug, levels = rev(unique(drug)))) %>%
    dplyr::arrange(Combinations, dplyr::desc(Percentage)) %>%
    dplyr::mutate(drug = factor(drug, levels = rev(drug)))

  # ---- Plot combinations by group ----
  # One panel per group; panel heights follow the number of sets in each.
  plot_b <- ggplot2::ggplot(df2, ggplot2::aes(x = Percentage, y = drug)) +
    ggplot2::geom_col(fill = bar_colour, width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(Percentage, accuracy = 0.1)),
      hjust = -0.15, size = 3.2, colour = "#52514e"
    ) +
    ggplot2::facet_grid(Combinations ~ ., scales = "free_y", space = "free_y", switch = "y") +
    pct_axis +
    ggplot2::labs(title = title_sets) +
    theme_bars +
    ggplot2::theme(
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, vjust = 1,
                                                face = "bold", colour = "#0b0b0b")
    )

  # ---- Output ----
  if (isTRUE(plot)) {
    final_plot <- patchwork::wrap_plots(plot_a, plot_b, ncol = 2, nrow = 1)
    return(final_plot)
  } else {
    return(list(overall = summary_number, breakdown = df2))
  }
}
