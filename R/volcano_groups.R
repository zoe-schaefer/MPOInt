#' Generate a volcano plot with groups of interest.
#'
#' `volcano_groups` generates a volcano plot given data, threshold values, and a subset of points to label.
#'
#' @param df A data frame containing points to be plotted.
#' @param fc_col A column in `df` with fold change values to be plotted.
#' @param p_col A column in `df` with P-values to be plotted. These values should not be log-transformed.
#' @param label_col A column in `df` with values to be used if groups are not provided.
#' @param fc_thr A scalar value indicating the fold change threshold. Default = 0.5.
#' @param p_thr A scalar value indicating the P-value threshold. Default = 0.05.
#' @param compute_fc_log2 A Boolean to determine if the `fc_col` values should be additionally log2-transformed. Default = FALSE.
#' @param group_list A named list of lists where each list is a group of peptides, labeled with the group name.
#' @param ... Optional.
#'
#' @returns A ggplot2 object.
#'
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' volcano_groups(df = THP1_total, fc_col = "log2(FC)", p_col = "P value",
#'   label_col = "SYMBOL2", group_list = list("Group 1" =
#'    c("MORF4L1", "SAMSN1", "DTX3L", "MBD3",  "HSF1",
#'   "EYA3", "NCAPD2", "TAF9", "PPP4C", "ZNF638"),
#'   "Group 2" = c("CD9", "CD82", "NFKB1", "BST2", "CTSH"))) +
#'   ggtitle("Grouped volcano plot") +
#'   xlab("log2FC") + ylab("-log10P")
#' }
volcano_groups <- function(df, fc_col, p_col, fc_thr = 0.5, p_thr = 0.05,
                        label_col = NA, compute_fc_log2 = FALSE, group_list = NULL, ...){
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
  if(!is.list(group_list)) {
    stop(paste(group_list, " is not a list of lists.", sep = ""))
  }
  if(!is.logical(compute_fc_log2)) {
    stop(paste(compute_fc_log2, " is not a Boolean value.", sep = ""))
  }

  df <- as.data.frame(df)

  # If FC values should be transformed, do that. Else, set the x column to the FC column.
  if(compute_fc_log2 == TRUE){
    df$x_col <- log2(df[[fc_col]])
  } else {
    df$x_col <- df[[fc_col]]
  }

  # Set up a raw p-values column to use for thresholds and labeling
  df$p_raw <- df[[p_col]]
  df$p_log10 <- -log10(df$p_raw)
#
#   # Setting categories for differential expression:
#   # not significant or P-value and FC thresholds met
#   df$DiffExp <- case_when(
#     abs(df$x_col) >= fc_thr & df$p_raw <= p_thr ~ "FC and P val.",
#     .default = "NS"
#   )

  # # Setting alpha values and labels so the significant points are most visible
  # df$DEAlpha <- 0.5
  # df$DEAlpha[abs(df$x_col) >= fc_thr & df$p_raw <= p_thr] <- 0.8

  df$group_col <- NA

  # Check for groups
  if(!is.null(group_list)) {
    # Iterate over group names
    df$group_col <- "None"
    for (list_name in names(group_list)) {
      # Create a temporary inner list
      temp_list <- group_list[[list_name]]
      df$group_col <- dplyr::case_when(
        df[[label_col]] %in% temp_list ~ list_name,
        .default = df$group_col
      )
      df$labels <- dplyr::case_when(
        df[[label_col]] %in% temp_list ~ df[[label_col]],
        .default = df$labels
      )
    }
  } else if (!is.na(label_col)) {
    # Setting labels - if no groups and this is not NA, label hits
    df$labels <- if_else(abs(df$x_col) >= fc_thr & df$p_raw <= p_thr, df[[label_col]], "")
  }

  df$group_col <- as.factor(df$group_col)
  df$group_col <- relevel(df$group_col, "None")
  df <- dplyr::arrange(df, group_col)


  # Create base ggplot object to be returned - this can be further modified with other ggplot calls
  plot <- ggplot(df, aes(x = x_col, y = p_log10)) +
    geom_point(aes(color = group_col), na.rm = TRUE) +
    geom_vline(xintercept = c(fc_thr, -fc_thr),
               alpha = 0.5, linetype = "dashed") +
    geom_hline(yintercept = -log10(p_thr),
               alpha = 0.5, linetype = "dashed") +
    geom_text_repel(aes(label = labels), show.legend = FALSE, size = rel(3), na.rm = TRUE) +
    theme(legend.position = "bottom") +
    guides(color = guide_legend("Group"))

  return(plot)
}
