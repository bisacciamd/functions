library(ggplot2)
library(gt)

# theme for {ggplot} objects
ehj_palette <- c(
  "#5f3393", "#06b6b3",
  "#ad1122", "#000000",
  "#005694", "#fbb800", "#ef7918"
)

theme_bisaccia <- function() {
  theme_classic() +
    theme(
      axis.title = element_text(family = "archivo-narrow"),
      title = element_text(family = "archivo-narrow")
    )
}

# theme for {gt} tables
gt_theme_bisaccia <- function(gt_tbl) {
  n_rows <- nrow(gt_tbl$`_data`)

  gt_tbl %>%
    tab_options(
      data_row.padding = px(6),
      heading.align = "left",
      column_labels.background.color = "#4a6d88",
      heading.title.font.size = px(26),
      heading.subtitle.font.size = px(14),
      table_body.hlines.width = px(0)
    ) %>%
    tab_style(
      style = cell_text(
        color = "#4a6d88",
        weight = "bold",
        font = "Archivo Narrow"
      ),
      locations = cells_title(groups = "title")
    ) %>%
    tab_style(
      style = cell_fill(color = "#c6cdd7"),
      # Use number of rows from data set
      locations = cells_body(rows = seq(1, n_rows, 2))
    ) %>%
    opt_table_font(
      google_font("Archivo Narrow")
    )
}

