## Code to prepare `THP1_total` dataset
# Reading in total proteome data, selecting a range to reduce initial variable size
THP1_total <- readxl::read_excel("data-raw/33060_Control_vs_Infection_JMI_abundances_normalizedabundances.xlsx",
                                 range = readxl::cell_cols("D:BI")
)
# Selecting just the columns we're interested in
THP1_total <- dplyr::select(THP1_total, all_of(c(
  "Accession", "Abundance Ratio (log2): (JMI) / (Control)",
  "Abundance Ratio Adj. P-Value: (JMI) / (Control)",
  "Found in Sample: F1: Sample Control",
  "Found in Sample: F2: Sample Control",
  "Found in Sample: F3: Sample Control",
  "Found in Sample: F4: Sample JMI",
  "Found in Sample: F5: Sample JMI",
  "Found in Sample: F6: Sample JMI",
  "Found in Sample: F7: Sample JMI",
  "Found in Sample: F8: Sample JMI"
)))

usethis::use_data(THP1_total, overwrite = TRUE)

