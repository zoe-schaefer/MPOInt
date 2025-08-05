## Code to prepare `THP1_phospho` dataset
THP1_phospho <- readxl::read_excel("data-raw/PXD_raw.xlsx", range = cell_cols("B:E"))

usethis::use_data(THP1_phospho, overwrite = TRUE)
