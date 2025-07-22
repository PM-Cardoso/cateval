#' Plotting optimal drug combinations
#'
#' This function generates a side-by-side bar plot visualizing:
#' 1. The proportion of individuals with 1 to 5 optimal drugs identified.
#' 2. The most frequent drug combinations by category (e.g., 1-drug, 2-drug).
#'
#' @param data A character vector where each element is a comma-separated string of drugs considered optimal for an individual.
#' @param groups A named list where each name represents a category (e.g., "1-drug combination") and each value is a vector of drug counts (e.g., 1, 2, 4:5).
#' @param plot Logical. If TRUE (default), returns a combined patchwork plot. If FALSE, returns a list of summary data frames.
#'
#' @return A patchwork object combining two ggplot2 bar plots side-by-side.
#'
#' @importFrom dplyr mutate filter count group_by ungroup arrange rename
#' @importFrom ggplot2 ggplot aes geom_text geom_col geom_hline coord_flip scale_y_continuous labs theme_minimal theme scale_x_continuous
#' @importFrom stringr str_count
#' @importFrom tibble tibble
#' @importFrom patchwork wrap_plots
#' @importFrom purrr map_dfr
#' @importFrom scales percent
#' @importFrom tidyr drop_na
#'
#' @examples
#' groups <- list(
#'   "1-drug combination" = 1,
#'   "2-drug combinations" = 2,
#'   "3-drug combinations" = 3,
#'   "3/4/5-drug combinations" = 4:5
#' )
#' optimal_drug_comparison_plot(c("SGLT2,DPP4,GLP1", "GLP1", "TZD,DPP4"), groups)
#' 
#' @export
#' 
optimal_drug_comparison_plot <- function(data, groups, plot = TRUE) {
  
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
  
  # ---- Plot number of drugs ----
  plot_a <- summary_number %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(drug, rev(drug)), y = Percentage)) +
    ggplot2::geom_text(ggplot2::aes(label = drug), hjust = -0.5, size = 8) +
    ggplot2::geom_col(fill = "#076fa2") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      breaks = seq(0, max(summary_number$Percentage) + 0.05, by = 0.05)
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Number of predicted optimal drug",
      subtitle = "Proportion of individuals with optimal drugs (>3 mmol/mol benefit)"
    ) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
  
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
  plot_b <- ggplot2::ggplot(df2, ggplot2::aes(x = Percentage, y = drug)) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(Percentage, accuracy = 0.1)),
      hjust = -0.1, size = 3
    ) +
    ggplot2::geom_col(fill = "#076fa2") +
    ggplot2::labs(title = "Predicted HbA1c-optimal drug classes") +
    ggplot2::scale_x_continuous(
      labels = scales::percent,
      breaks = seq(0, max(df2$Percentage) + 0.05, by = 0.05),
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  
  # ---- Output ----
  if (isTRUE(plot)) {
    final_plot <- patchwork::wrap_plots(plot_a, plot_b, ncol = 2, nrow = 1)
    return(final_plot)
  } else {
    return(list(overall = summary_number, breakdown = df2))
  }
}
