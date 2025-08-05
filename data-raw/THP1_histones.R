## Code to prepare `THP1_histones` dataset
THP1_histones <- readxl::read_excel("data-raw/histone_ratios.xls", range = cell_limits(c(4, 1), c(NA, 11)), col_names = c("PTM", "C1", "I1", "C2", "I2", "C3", "I3", "C4", "I4", "C5", "I5"))

usethis::use_data(THP1_histones, overwrite = TRUE)

