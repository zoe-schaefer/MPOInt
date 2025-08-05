## Code to prepare `genes` dataset

genes <- readr::read_delim("data-raw/genes.csv",
                           col_select = c("HGNC approved symbol", "UniProt ID (human)",
                                          "Function", "Modification", "Protein complex",
                                          "Target entity", "Product"),
                           col_types = c("c"))

usethis::use_data(genes, overwrite = TRUE)
