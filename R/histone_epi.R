#' Process individual data for histone proteome.
#'
#' `histone_epi` takes a dataframe with histone mass spectrometry data and outputs a log2-transformed fold change as well as P-values and a T score.
#'
#' @param df A dataframe containing histone mass spectrometry data. This dataframe should include measurements for two conditions.
#' @param label_col A column containing PTM/fragment names.
#' @param base_cols A list of columns containing baseline values.
#' @param target_cols A list of columns containing average values for your condition of interest.
#' @param ... Optional.
#'
#' @returns A dataframe with labels, log2(FC), T scores, and P-values
#'
#' @import dplyr
#' @import tibble
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' histone_epi(df = THP1_single_raw, label_col = "PTM",
#' base_cols = c("C1", "C2", "C3", "C4", "C5"),
#' target_cols = c("I1", "I2", "I3", "I4", "I5"))
#' }
histone_epi <- function(df, label_col, base_cols, target_cols, ...){
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
  if(!is.character(target_cols)){
    stop(paste(target_cols, " is not a list.", sep = ""))
  }

  target_mean <- NULL
  base_mean <- NULL
  is_constant <- NULL

  df <- as.data.frame(df)
  result_df <- data.frame()
  base_num <- length(base_cols)
  target_num <- length(target_cols)

  df <- df %>%
    # Physically group two sets
    dplyr::relocate(dplyr::all_of(target_cols), .after = last_col()) %>%
    dplyr::relocate(label_col)

  df_clean <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(base_mean = mean(dplyr::c_across(dplyr::all_of(base_cols))),
           target_mean = mean(dplyr::c_across(dplyr::all_of(target_cols))),
           "log2(FC)" = log2(dplyr::c_across(target_mean)/dplyr::c_across(base_mean)),
           # Checking if constant to avoid p value errors
           is_constant = all(dplyr::c_across(dplyr::all_of(base_cols)) == dplyr::c_across(dplyr::all_of(target_cols))))
  df_clean <- tibble::column_to_rownames(df_clean, var = "PTM")
  df_clean$PTM <- rownames(df_clean)

  df_pvals <- df_clean %>%
    dplyr::rowwise() %>%
    filter(is_constant == FALSE) %>%
    dplyr::mutate("P value" = stats::t.test(c_across(all_of(base_cols)), c_across(all_of(target_cols)),
                              alternative = "two.sided", var.equal = TRUE)$p.value)

  df_clean <- dplyr::full_join(df_clean, df_pvals)
  df_clean$"log2(FC)"[df_clean$"log2(FC)" == "NaN"] <- NA

  return(df_clean)
}
