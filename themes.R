# Load necessary libraries
library(sysfonts)
library(showtext)
library(ggplot2)
library(gt)

# Add Google font "Archivo Narrow" and enable it for rendering
font_add_google("Archivo Narrow", "archivo-narrow")
showtext_auto()

# Define color palette
ehj_palette <- c(
  "#5f3393", "#06b6b3",
  "#ad1122", "#000000",
  "#005694", "#fbb800", "#ef7918"
)

# Define ggplot theme with Archivo Narrow font
theme_bisaccia <- function() {
  theme_classic() %+replace%  # Use theme_classic as base
    theme(
      axis.title = element_text(family = "archivo-narrow"),
      axis.text = element_text(family = "archivo-narrow"),
      title = element_text(family = "archivo-narrow"),
      plot.title = element_text(family = "archivo-narrow", size = 14, face = "bold"),
      plot.subtitle = element_text(family = "archivo-narrow", size = 12),
      plot.caption = element_text(family = "archivo-narrow", size = 10)
    )
}

# Function for adding color scale
scale_color_ehj <- function() {
  scale_color_manual(values = ehj_palette)
}

# Function for adding fill scale
scale_fill_ehj <- function() {
  scale_fill_manual(values = ehj_palette)
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

