## code to prepare `MCF7_total` dataset goes here
MCF7_total <- readxl::read_excel("data-raw/MCF_total.xlsx", range = cell_cols("A:V"))

MCF7_total <- dplyr::select(MCF7_total, all_of(c("Protein accession", "M90/F Ratio", "M90/F P value")))
# Rename columns containing the log2FC (if present), Accession, and P value so they're easier to use

MCF7_total <- dplyr::rename(MCF7_total, all_of(c("Accession" = "Protein accession",
                                               "P value" = "M90/F P value")))
MCF7_total <- left_join(MCF7_total, OrganismDbi::select(
  Homo.sapiens,
  keys = keys(Homo.sapiens,keytype = "UNIPROT"),
  columns = c("ENTREZID", "SYMBOL"), keytype = "UNIPROT"), by = join_by(Accession == UNIPROT), keep = FALSE, relationship = "many-to-many")

# Drop invalid data (NA), drop any duplicate cases, rename columns for clarity
MCF7_total <- MCF7_total %>%
  tidyr::drop_na(SYMBOL) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# If needed, calculate the log2FC for mean data
MCF7_total$`log2(FC)` <- log2(MCF7_total$`M90/F Ratio`)

# Assign symbols to row names
MCF7_total <- MCF7_total %>%
  dplyr::mutate("SYMBOL2" = `SYMBOL`) %>%
  tibble::column_to_rownames(var = "SYMBOL")

usethis::use_data(MCF7_total, overwrite = TRUE)
