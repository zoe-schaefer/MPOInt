#' Process summary data for histone proteome.
#'
#' `histone_summary` takes a dataframe with histone mass spectrometry data and outputs a log2-transformed fold change as well as P-values and a T score.
#'
#' @param df A dataframe containing histone mass spectrometry data. This dataframe should include summary statistics (average value, standard deviation) for at least two conditions.
#' @param label_col A column containing PTM/fragment names.
#' @param base_cols A list of at least one column from which to generate a baseline.
#' @param base_std If more than one column is supplied to `base_cols`, leave this as NA. Otherwise, if only one column is supplied, provide the standard deviation for that column.
#' @param target_col A column containing average values for your condition of interest.
#' @param target_std A column containing standard deviations for your condition of interest.
#' @param ... Optional.
#'
#' @returns A dataframe with labels, log2(FC), T scores, and P-values
#'
#' @import dplyr
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' histone_summary(df = MCF7_histones, label_col = "PTM",
#' base_cols = c("293 ave (5 reps)", "HaCAT ave (5 reps)", "hESC ave (6 reps)",
#' "HFF ave (7 reps)", "Mdm13 ave (5 reps)"), base_std = NA,
#' target_col = "MCF7 ave (5 reps)", target_std = "MCF7 std")
#' }
histone_summary <- function(df, label_col, base_cols, base_std = NA, target_col, target_std, ...){
  # Check valid inputs first
  if(!is.data.frame(df)) {
    stop(paste(df, " is not a data frame.", sep = ""))
  }
  if(!is.character(df[[label_col]])) {
    stop(paste(label_col, " is not a character vector.", sep = ""))
  }
  if(!is.character(base_cols)){
    stop(paste(base_cols, " is not a list.", sep = ""))
  }
  if(!is.numeric(df[[base_std]]) & !is.na(base_std)) {
    stop(paste(base_std, " is not numeric or NA.", sep = ""))
  }
  if(!is.numeric(df[[target_col]])) {
    stop(paste(target_col, " is not numeric.", sep = ""))
  }
  if(!is.numeric(df[[target_std]])) {
    stop(paste(target_std, " is not numeric.", sep = ""))
  }

  base_avg = NULL
  `T score` = NULL

  df <- as.data.frame(df)
  result_df <- data.frame()
  base_num <- length(base_cols)

  # Creating an overall baseline
  if (is.na(base_std) & length(base_cols) > 1) {
    result_df <- dplyr::mutate(df,
                               "PTM" = df[[label_col]],
                               target_col = df[[target_col]],
                               target_std = df[[target_std]],
                               base_avg = apply(dplyr::pick(base_cols), 1, mean),
                               base_std = apply(dplyr::pick(base_cols), 1, stats::sd),
                               .keep = "used")
  } else {
    result_df <- dplyr::mutate(df,
                               "PTM" = df[[label_col]],
                               target_col = df[[target_col]],
                               target_std = df[[target_std]],
                               base_avg = apply(dplyr::pick(base_cols), 1, mean),
                               base_std = df[[base_std]],
                               .keep = "used")
  }
  result_df <- dplyr::mutate(result_df,
                             "log2(FC)" = log2(target_col/base_avg),
                             "T score" = (target_col - base_avg)/(target_std/base_num),
                             "P value" = stats::pt(q = `T score`, df = (base_num - 1)))
  return(result_df)
}
