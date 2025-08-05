## code to prepare `MCF7_histones` dataset goes here
MCF7_histones <- readxl::read_excel("data-raw/MCF_histone_frag.xlsx")

MCF7_histones <- MCF7_histones %>%
  dplyr::select(c("PTM", "MCF7 ave (5 reps)", "MCF7 std",
                  "293 ave (5 reps)", "HaCAT ave (5 reps)",
                  "hESC ave (6 reps)", "HFF ave (7 reps)",
                  "Mdm13 ave (5 reps)")) %>%
  dplyr::mutate("ctrl ave" = rowMeans(pick("293 ave (5 reps)", "HaCAT ave (5 reps)",
                                    "hESC ave (6 reps)", "HFF ave (7 reps)",
                                    "Mdm13 ave (5 reps)"))) %>%
  dplyr::mutate("log2(FC)" = log2(`MCF7 ave (5 reps)`/`ctrl ave`)) %>%
  dplyr::mutate("T score" = NA) %>%
  dplyr::mutate("P value" = NA)
MCF7_histones <- MCF7_histones %>%
  dplyr::mutate("T score" = (`MCF7 ave (5 reps)` - `ctrl ave`)/(`MCF7 std`/sqrt(5))) %>%
  dplyr::mutate("P value" = pt(q = `T score`, df = 4))

usethis::use_data(MCF7_histones, overwrite = TRUE)
