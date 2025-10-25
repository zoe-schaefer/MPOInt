#' Clean up names from EpiProfile.
#'
#' `clean_names` takes a dataframe with histone mass spectrometry data and outputs two non-transformed dataframes with user-friendly PTM names.
#'
#' @param df A dataframe containing histone mass spectrometry data.
#' @param label_col A column containing PTM/fragment names.
#' @param ... Optional.
#'
#' @returns Two dataframes with adjusted labels, up to 4 co-occurring PTMs per fragment.
#'
#' @import dplyr
#' @import stringr
#' @import tidyr
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clean_names(df = THP1_histones, label_col = "PTM")
#' }
clean_names <- function(df, label_col, ...){
  # Check valid inputs first
  if(!is.data.frame(df)) {
    stop(paste(df, " is not a data frame.", sep = ""))
  }
  if(!is.character(df[[label_col]])) {
    stop(paste(label_col, " is not a character vector.", sep = ""))
  }

  PTM <- NULL
  PTM1 <- NULL
  PTM2 <- NULL
  PTM3 <- NULL
  PTM4 <- NULL
  single_PTM <- NULL

  df <- as.data.frame(df)
  single_df <- data.frame()
  frag_df <- data.frame()
  single_df_raw <- data.frame()
  single_df_long <- data.frame()
  single_df_unmod <- data.frame()

  frag_df <- tidyr::drop_na(df)
  # Suppresses warnings - errors come up if regex strings are empty
  suppressWarnings(frag_df$PTM <- frag_df$PTM %>%
    stringr::str_replace_all("H3_3_8", "H3 3-8") %>%
    stringr::str_replace_all("H3_9_17", "H3 9-17") %>%
    stringr::str_replace_all("H3_18_26", "H3 18-26") %>%
    stringr::str_replace_all("H3_27_40", "H3 27-40") %>%
    stringr::str_replace_all("H33_27_40", "H3.3 27-40") %>%
    stringr::str_replace_all("H3_41_49", "H3 41-49") %>%
    stringr::str_replace_all("H3_54_63", "H3 54-63") %>%
    stringr::str_replace_all("H3_73_83", "H3 73-83") %>%
    stringr::str_replace_all("H3_117_128", "H3 117-128") %>%
    stringr::str_replace_all("H4_4_17", "H4 4-17") %>%
    stringr::str_replace_all("H4_20_23", "H4 20-23") %>%
    stringr::str_replace_all("H4_24_35", "H4 24-35") %>%
    stringr::str_replace_all("H4_40_45", "H4 40-45") %>%
    stringr::str_replace_all("H4_79_92", "H4 79-92") %>%
    stringr::str_replace_all("H1_1_35", "H1 1-35") %>%
    stringr::str_replace_all("H1_54_81", "H1 54-71") %>%
    stringr::str_replace_all("H2A1_36_42", "H2A1 36-42") %>%
    stringr::str_replace_all("H2A3_36_42", "H2A3 36-42") %>%
    stringr::str_replace_all("H2AX_36_42", "H2AX 36-42") %>%
    stringr::str_replace_all("H2A1_4_11", "H2A1 4-11") %>%
    stringr::str_replace_all("H2AJ_4_11", "H2AJ 4-11") %>%
    stringr::str_replace_all("H2AX_4_11", "H2AX 4-11") %>%
    stringr::str_replace_all("H2A1_1_11", "H2A1 1-11") %>%
    stringr::str_replace_all("H2AV_1_19", "H2AV 1-19") %>%
    stringr::str_replace_all("H2AZ_1_19", "H2AZ 1-19") %>%
    stringr::str_replace_all("H2A1_12_17", "H2A1 12-17") %>%
    stringr::str_replace_all("H2A3_12_17", "H2A3 12-17") %>%
    stringr::str_replace_all("H2A1_72_77", "H2A1 72-77") %>%
    stringr::str_replace_all("H2A_82_88", "H2A 82-88") %>%
    stringr::str_replace_all("H2A_1_88", "H2A 1-88") %>%
    stringr::str_replace_all("H2B_1_29", "H2B 1-29"))

  frag_df <- dplyr::distinct(frag_df, PTM, .keep_all = TRUE)

  single_df_raw <- frag_df
  suppressWarnings(single_df_raw <- single_df_raw[!grepl("\\.", single_df_raw$PTM),])
  suppressWarnings(single_df_unmod <- single_df_raw[grepl("unmod", single_df_raw$PTM),])
  suppressWarnings(single_df_raw <- single_df_raw[!grepl("unmod", single_df_raw$PTM),])

  suppressWarnings(single_df_raw <- single_df_raw %>%
                     dplyr::mutate("subset" = stringr::str_extract(`PTM`, "([^\ ]+$)")) %>%
                     dplyr::mutate("PTM4" = stringr::str_split_i(`subset`, "(?=K)", i = -1)) %>%
                     dplyr::mutate("PTM3" = stringr::str_split_i(`subset`, "(?=K)", i = -2)) %>%
                     dplyr::mutate("PTM2" = stringr::str_split_i(`subset`, "(?=K)", i = -3)) %>%
                     dplyr::mutate("PTM1" = stringr::str_remove(`subset`, `PTM2`)) %>%
                     dplyr::mutate("PTM1" = stringr::str_remove(`PTM1`, `PTM3`)) %>%
                     dplyr::mutate("PTM1" = stringr::str_remove(`PTM1`, `PTM4`)))

  single_df_raw <- replace(single_df_raw, single_df_raw == "", NA)

  single_df_long <- single_df_raw %>%
    tidyr::pivot_longer(cols = c("PTM4", "PTM3", "PTM2", "PTM1"), names_to = NULL,
                 values_to = .data$"single_PTM", values_drop_na = TRUE) %>%
    dplyr::mutate("PTM" = stringr::str_remove(`PTM`, `subset`)) %>%
    tidyr::unite(col = "PTM", c(`PTM`, `single_PTM`), sep = " ", remove = TRUE) %>%
    dplyr::mutate("PTM" = stringr::str_squish(PTM))

  single_df_long$PTM <- stringr::str_remove(single_df_long$PTM, " .+ ") %>%
    stringr::str_replace("me1", "me")
  single_df <- dplyr::distinct(dplyr::bind_rows(single_df_long, single_df_unmod), PTM, .keep_all = TRUE)

  return(list(single_df, frag_df))
}
