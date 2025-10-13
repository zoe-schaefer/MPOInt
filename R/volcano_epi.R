#' Generate a volcano plot for epiproteome data.
#'
#' `volcano_epi` generates a volcano plot given data and threshold values.
#'
#' @param df A data frame containing points to be plotted.
#' @param fc_col A column in `df` with fold change values to be plotted.
#' @param p_col A column in `df` with P-values to be plotted. These values should not be log-transformed.
#' @param fc_thr A scalar value indicating the fold change threshold. Default = 0.5.
#' @param p_thr A scalar value indicating the P-value threshold. Default = 0.05.
#' @param label_col A column in `df` with which points should be labeled. If `label_col = NA`, default labels will be applied to significant points.
#' @param ... Optional.
#'
#' @returns A ggplot2 object.
#'
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @import tidyr
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' volcano_epi(df = total_raw, fc_col = "FC", p_col = "Pval", label_col = "SYMBOL2") +
#' ggtitle("Volcano plot") + xlab("log2FC") + ylab("-log10P")
#' }
volcano_epi <- function(df, fc_col, p_col, fc_thr = 0.5, p_thr = 0.05,
                        label_col = NA, ...){

  shape_col <- NULL
  x_col <- NULL
  p_log10 <- NULL
  DiffExp <- NULL

  # Check valid inputs first
  if(!is.data.frame(df)) {
    stop(paste(df, " is not a data frame.", sep = ""))
  }
  if(!is.numeric(df[[fc_col]])) {
    stop(paste(fc_col, " is not numeric.", sep = ""))
  }
  if(!is.numeric(df[[p_col]])) {
    stop(paste(p_col, " is not numeric.", sep = ""))
  }
  if(!is.numeric(fc_thr)) {
    stop(paste(fc_thr, " is not numeric.", sep = ""))
  }
  if(!is.numeric(p_thr)) {
    stop(paste(p_thr, " is not numeric.", sep = ""))
  }
  if(!is.character(label_col) & !is.na(label_col)) {
    stop(paste(label_col, " is not a character vector.", sep = ""))
  }

  df <- as.data.frame(df)
  df$shape_col <- ""

  # Set the x column to the FC column
  df$x_col <- df[[fc_col]]

  # Set up a raw p-values column to use for thresholds and labeling
  df$p_raw <- df[[p_col]]
  df$p_log10 <- -log10(df$p_raw)

  # Setting categories for differential expression:
  # not significant or P-value and FC thresholds met
  df$DiffExp <- case_when(
    abs(df$x_col) >= fc_thr & df$p_raw <= p_thr ~ "FC and P val.",
    .default = "NS"
  )
  # Setting labels - if not NA, use the column given to label hits
  if(!is.na(label_col)) {
    df$labels <- if_else(abs(df$x_col) >= fc_thr & df$p_raw <= p_thr, df[[label_col]], "")
  }

  # Remove invalid values that generate errors
  df_plot <- df %>%
    drop_na() %>%
    filter(x_col != "-Inf" & x_col != "Inf")
  df_plot$shape_col <- "16"

  FC_min <- min(df_plot$x_col)
  FC_max <- max(df_plot$x_col)

  df_invalid <- df %>%
    drop_na() %>%
    filter(x_col == "-Inf" | x_col == "Inf")
  df_invalid$x_col[df_invalid$x_col == "-Inf"] <- FC_min*1.5
  df_invalid$x_col[df_invalid$x_col == "Inf"] <- FC_max*1.5
  df_invalid$shape_col <- if (length(df_invalid$shape_col) > 0) {
    "17"
  }

  df_plot <- full_join(df_plot, df_invalid)

  # Create base ggplot object to be returned - this can be further modified with other ggplot calls
  plot <- ggplot(df_plot, aes(x = x_col, y = p_log10)) +
    geom_point(aes(color = DiffExp, shape = shape_col), na.rm = TRUE) +
    geom_vline(xintercept = c(fc_thr, -fc_thr),
               alpha = 0.5, linetype = "dashed") +
    expand_limits(x = c(1.1*FC_min, 1.1*FC_max)) +
    geom_hline(yintercept = -log10(p_thr),
               alpha = 0.5, linetype = "dashed") +
    geom_text_repel(aes(label = labels), show.legend = FALSE, size = rel(3), na.rm = TRUE) +
    guides(shape = "none", color = "none") +
    theme(legend.position = "bottom") +
    scale_color_manual(
      values = c("FC and P val." = "steelblue2", "NS" = "gray70")
    )

  return(plot)
}
