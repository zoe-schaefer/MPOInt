## code to prepare `MCF7_phospho` dataset goes here
MCF7_phospho <- readxl::read_excel("data-raw/MCF_phospho.xlsx", range = cell_cols("A:AA"))

# Rename columns containing the log2FC (if present), Accession, and P value so they're easier to use

MCF7_phospho <- dplyr::rename(MCF7_phospho, all_of(c("SYMBOL" = "Gene name", "Accession" = "Protein accession", "P value" = "M90/F P value")))
MCF7_phospho <- left_join(MCF7_phospho, OrganismDbi::select(
  Homo.sapiens,
  keys = keys(Homo.sapiens,keytype = "UNIPROT"),
  columns = "ENTREZID", keytype = "UNIPROT"), by = join_by(Accession == UNIPROT), keep = FALSE, relationship = "many-to-many")

# Drop invalid data (NA), drop any duplicate cases, rename columns for clarity
MCF7_phospho <- MCF7_phospho %>%
  tidyr::drop_na(SYMBOL) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

# If needed, calculate the log2FC for mean data
MCF7_phospho$`log2(FC)` <- log2(MCF7_phospho$`M90/F Ratio`)

# Assign symbols to row names
MCF7_phospho <- MCF7_phospho %>%
  dplyr::mutate("SYMBOL2" = `SYMBOL`) %>%
  tibble::column_to_rownames(var = "SYMBOL")

MCF7_phospho <- dplyr::select(MCF7_phospho, all_of(c(
  "Accession", "Position", "Amino acid", "F1", "F2", "F3", "M901", "M902", "M903", "M90/F Ratio", "P value", "ENTREZID", "log2(FC)", "SYMBOL2"
)))


usethis::use_data(MCF7_phospho, overwrite = TRUE)
