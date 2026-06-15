#' Generate a multi-level proteome network from total proteome, phosphoproteome, and epiproteome data.
#'
#' `multi_net` generates an R object ready to generate a network in Python given three datasets.
#'
#' @param genes_df A dataframe containing the exported information from EpiFactors.
#' @param total_df A dataframe containing total proteome data.
#' @param phospho_df A dataframe containing phosphoproteome data.
#' @param symbol_col A column in both `total_df` and `phospho_df` containing the gene symbol names for identified proteins.
#' @param fc_col A column in both `total_df` and `phospho_df` containing the log2-transformed fold change values.
#' @param epi_single_df A dataframe containing histone proteome data. If using the full pipeline, make sure this is the single PTM dataframe and not fragments.
#' @param epi_symbol_col A column in `epi_single_df` containing the gene symbol names for identified PTMs.
#' @param epi_fc_col A column in `epi_single_df` containing the log2-transformed fold change values for identified PTMs.
#' @param ... Optional.
#'
#' @returns An object with the following attributes:
#' * `grid`: A plot displaying the complexes, PTMs, and networks in an intersecting grid.
#' * `complex_list`: A list of the complexes present in the network.
#' * `member_list`: A list of the complexes, their complex members (proteins), and associated PTMs present in the network.
#' * `PTM_list`: A list of the identified PTMs and their differential expression.
#' * `gene_list`: A subset of the `genes_df` variable containing proteins identified.
#' * `net_edges`: The edges connecting network vertices, which identify member-PTM relationships.
#' * `net_vertices`: The vertices (PTMs and members) in the network.
#' * `ts`: A larger set of information on the proteins in the network. This contains information including statistical significance values, fold change values, and proteome source identification.
#'
#' @import dplyr
#' @import stringr
#' @import tidyr
#' @import ggplot2
#'
#' @importFrom rlang .data
#' @importFrom reticulate r_to_py
#' @importFrom scales breaks_width
#'
#' @export
#'
#' @examples
#' \dontrun{
#' multi_net()
#' }
multi_net <- function(genes_df, total_df, phospho_df, symbol_col, fc_col,
                      epi_single_df, epi_symbol_col, epi_fc_col, ...){
  net_obj <- list()

  # Check valid inputs first ----
  if(!is.data.frame(genes_df)) {
    stop(paste(genes_df, " is not a data frame.", sep = ""))
  }
  if(!is.data.frame(total_df)) {
    stop(paste(total_df, " is not a data frame.", sep = ""))
  }
  if(!(symbol_col %in% colnames(total_df))) {
    stop(paste(symbol_col, " is not a character vector in ", total_df, ".", sep = ""))
  }
  if(!(fc_col %in% colnames(total_df))) {
    stop(paste(fc_col, " is not a character vector in ", total_df, ".", sep = ""))
  }
  if(!is.data.frame(phospho_df)) {
    stop(paste(phospho_df, " is not a data frame.", sep = ""))
  }
  if(!(symbol_col %in% colnames(phospho_df))) {
    stop(paste(symbol_col, " is not a character vector in ", phospho_df, ".", sep = ""))
  }
  if(!(fc_col %in% colnames(phospho_df))) {
    stop(paste(fc_col, " is not a character vector in ", phospho_df, ".", sep = ""))
  }
  if(!is.data.frame(epi_single_df)) {
    stop(paste(epi_single_df, " is not a data frame.", sep = ""))
  }
  if(!(epi_symbol_col %in% colnames(epi_single_df))) {
    stop(paste(epi_symbol_col, " is not a character vector in ", epi_single_df, ".", sep = ""))
  }
  if(!(epi_fc_col %in% colnames(epi_single_df))) {
    stop(paste(epi_fc_col, " is not a character vector in ", epi_single_df, ".", sep = ""))
  }

  symbol_col <- enquo(symbol_col)
  fc_col <- enquo(fc_col)


  # Combining dataframes ----
  total_join <- data.frame()
  total_part <- data.frame()
  phospho_join <- data.frame()
  phospho_part <- data.frame()
  joined_data <- data.frame()

  total_part <- dplyr::select(total_df, (!!fc_col), (!!symbol_col)) %>%
    dplyr::rename("Symbol" = (!!symbol_col), "log2(FC)" = (!!fc_col))
  phospho_part <- dplyr::select(phospho_df, (!!fc_col), (!!symbol_col)) %>%
    dplyr::rename("Symbol" = (!!symbol_col), "log2(FC)" = (!!fc_col))

  total_join <-
    dplyr::inner_join(genes_df,
                      total_part,
                      by = join_by("HGNC approved symbol" == "Symbol"),
                      multiple = "any") %>%
    dplyr::mutate("Source" = "Total")

  phospho_join <-
    dplyr::inner_join(genes_df,
                      phospho_part,
                      by = join_by("HGNC approved symbol" == "Symbol"),
                      multiple = "any") %>%
    dplyr::mutate("Source" = "Phospho")

  joined_data <- dplyr::full_join(total_join, phospho_join)

  # Assigning differential expression identifiers
  joined_data$Diff <- dplyr::case_when(joined_data[["log2(FC)"]] > 0 ~ "Up",
                                       joined_data[["log2(FC)"]] < 0 ~ "Down")

  # Formatting histone targets ----
  targets <- NULL
  targets <- stringr::str_replace_all(
    stringr::str_flatten_comma(
      stringr::str_unique(c(joined_data$Product, joined_data$`Target entity`))
      ), ", ", "|"
    )

  # Filter by that list
  epi_filtered <- NULL
  epi_filtered <- epi_single_df %>% dplyr::filter(grepl(targets, .data$PTM))
  epi_filtered$Diff <- dplyr::case_when(epi_filtered[[epi_fc_col]] > 0 ~ "Up",
                                        epi_filtered[[epi_fc_col]] < 0 ~ "Down",
                                      .default = "0")
  epi_filtered <- epi_filtered %>% dplyr::filter(!grepl("unmod", .data$PTM))

  # Pivot and format the gene list ----
  genes_pivoted <-
    dplyr::select(genes_df, c("HGNC approved symbol","Protein complex", "Target entity", "Product")) %>%
    dplyr::rename("Symbol" = "HGNC approved symbol",
                  "Complex" = "Protein complex",
                  "Target" = "Target entity",
                  "Product" = "Product") %>%
    tidyr::separate_wider_delim(cols = c("Complex", "Target", "Product"), delim = ", ",
                         names_sep = "_", too_few = "align_start")

  genes_pivoted <- genes_pivoted %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("Complex"), names_to = NULL,
                 values_to = "Complex", values_drop_na = TRUE) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("Target"), names_to = NULL,
                 values_to = "Target", values_drop_na = TRUE) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("Product"), names_to = NULL,
                 values_to = "Product", values_drop_na = TRUE)

  genes_data_joined <- dplyr::inner_join(genes_pivoted,
                                  dplyr::select(joined_data, c("HGNC approved symbol",
                                                               "log2(FC)", "Source")),
                                  by = join_by("Symbol" == "HGNC approved symbol"),
                                  relationship = "many-to-many")

  #return(list(joined_data, epi_filtered, genes_data_joined))

  # Getting the lists into the right format for the igraph visualization ----
  genes_histone_pivoted <-
    tidyr::pivot_longer(genes_data_joined, cols = c("Target", "Product"),
                        names_to = NULL, values_to = "PTM") %>% dplyr::distinct()
  genes_histone_pivoted$Diff <-
    dplyr::case_when(genes_histone_pivoted$`log2(FC)` > 0 ~ "Up",
                     genes_histone_pivoted$`log2(FC)` < 0 ~ "Down",
                     .default = "0")

  genes_histone_joined <-
    dplyr::right_join(genes_histone_pivoted,
                      dplyr::select(epi_filtered, c("PTM", "log2(FC)", "P_unadj", "P value", "Diff")),
                      by = c("PTM" = "PTM"), suffix = c("_gene", "_hist"),
                      relationship = "many-to-many") %>% drop_na("Source") %>% drop_na("log2(FC)_hist")

  genes_histone_joined <- genes_histone_joined %>% dplyr::relocate(.data$`Complex`, .after = .data$PTM)

  # Ensuring no NA values remain, replacing "#" with a clearer identifier
  genes_histone_joined[genes_histone_joined == "#"] <- "No ID"
  genes_histone_joined$Complex <- genes_histone_joined$`Complex` %>% tidyr::replace_na("No ID")

  # Generating a dataframe for edges and a dataframe for unique vertices with metadata ----
  complexes <- unique(genes_histone_joined$`Complex`)
  n_complexes <- seq(1, length(complexes))
  complex_df <- data.frame(complexes, n_complexes)

  # Ensuring all histone PTMs have the same complex and source labels
  histone_df <- dplyr::select(epi_filtered, c("PTM", "Diff", "log2(FC)"))
  histone_df$Symbol <- histone_df$PTM
  histone_df$Complex <- "No ID"
  histone_df$Source <- "Histone"

  vertices_df <-
    dplyr::bind_rows(histone_df, genes_histone_pivoted) %>%
    dplyr::relocate("Symbol")
  vertices_df <- vertices_df %>% dplyr::filter(grepl("^H", .data$PTM))
  vertices_df$Complex <- dplyr::case_when(vertices_df$`Complex` == "#" ~ "No ID",
                                          .default = vertices_df$`Complex`)

  # Creating a dataframe connecting edges (source and target vertices) ----
  complex_vertices <- data.frame(vertices = c(unique(genes_histone_joined$`Symbol`),
                                              unique(genes_histone_joined$`PTM`)))
  edges_df <-
    full_join(genes_histone_joined, complex_df, by = join_by("Complex" == "complexes")) %>%
    relocate(c("Symbol", "PTM"))
  edges_df <-
    left_join(complex_vertices, edges_df, by = join_by("vertices" == "Symbol")) %>%
    dplyr::rename("Symbol" = "vertices")

  edges_df$Symbol <- case_when(is.na(edges_df$`Symbol`) ~ edges_df$`PTM`,
                               .default = edges_df$`Symbol`)
  edges_df$PTM <- case_when(is.na(edges_df$`PTM`) ~ edges_df$`Symbol`,
                            .default = edges_df$`PTM`)
  edges_df$Source <- case_when(is.na(edges_df$`Source`) ~ "Histone",
                               .default = edges_df$`Source`)
  edges_df$Complex <- case_when(is.na(edges_df$`Complex`) ~ "No ID",
                                .default = edges_df$`Complex`)
  edges_df$Diff <- case_when(edges_df$`Source` == "Histone" ~ edges_df$`Diff_hist`,
                             .default = edges_df$`Diff_gene`)
  edges_df$`log2(FC)` <- case_when(edges_df$Source == "Histone" ~ edges_df$`log2(FC)_hist`,
                                   .default = edges_df$`log2(FC)_gene`)
  edges_df <- edges_df %>% tidyr::drop_na(.data$Diff)

  sig_network <- filter(edges_df, edges_df$`Complex` != "No ID") %>%  group_by(.data$`Complex`)
  sig_network["weight"] <- 1
  sig_network$Complex <- factor(sig_network$`Complex`, levels = unique(as.character(sig_network$`Complex`)))

  # Creating the first plot ----
  sn_collapsed <- sig_network %>%
    group_by(.data$Symbol, .data$Complex) %>%
    summarize(PTMs = paste(unique(.data$PTM), collapse = ", "))

  grid_plot <- sn_collapsed %>% ggplot(aes(x = .data$Complex, y = .data$Symbol, fill = .data$PTMs)) +
    geom_bin_2d(binwidth = 1, center = 0, color = "black", linewidth = 0.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank()) +
    scale_fill_viridis_d() +
    scale_y_discrete(limits = rev,
                     minor_breaks = breaks_width(1, 0.5)) +
    labs(title = "Network bins") +
    scale_x_discrete(minor_breaks = breaks_width(1, 0.5))

  # Print identified complexes ----
  final_complexes <-  vertices_df %>% filter(.data$`Source` != "Histone") %>%
    filter(.data$`Complex` != "No ID") %>% drop_na() %>% distinct()
  final_complex_list <-  sort(unique(final_complexes$`Complex`))

  #comp_list <- (cat(paste("**Networked Complexes:** ", paste(final_complex_list, collapse = ", "))))

  # List members of each
  final_members <- data.frame("Complex" = "", "Members" = "",
                              "PTMs" = "")

  for (complex in final_complex_list) {
    member_df <- filter(final_complexes, final_complexes$`Complex` == complex)
    member_list <-  paste0(unique(member_df$`Symbol`), collapse = ", ")
    member_PTMs <- paste0(unique(member_df$`PTM`), collapse = ", ")
    final_members <- final_members %>%
      add_row("Complex" = complex, "Members" = member_list,
              "PTMs" = member_PTMs)
  }

  # Differential PTMs
  final_PTMs <- final_complexes %>%  distinct(dplyr::pick(.data$`PTM`, .data$`Diff`)) %>%
    dplyr::rename("Differential expression" = "Diff")

  # Related gene info
  final_proteins <- unique(final_complexes$`Symbol`)
  final_func <- genes_df %>%
    dplyr::filter(.data$`HGNC approved symbol` %in% final_proteins) %>%
    dplyr::select(!c("UniProt ID (human)"))


  # Set up network ----
  vertices_df$Shapes <-
    case_when(vertices_df$`Diff` == "Up" ~ "triangle",
              vertices_df$`Diff` == "Down" ~ "triangleDown",
              vertices_df$`Diff` == "0" ~ "ellipse")
  vertices_df$Colors <-
    case_when(vertices_df$`Source` == "Total" ~ "#21908C",
              vertices_df$`Source` == "Phospho" ~ "#FDE725",
              vertices_df$`Source` == "Histone" ~ "#440154")
  vertices_df$Border <- "black"


  sig_network_py <-
    filter(vertices_df, vertices_df$`Complex` != "No ID" | vertices_df$`Source` == "Histone") %>%
    r_to_py()

  edges_py <-
    filter(edges_df, edges_df$`Complex` != "No ID" | edges_df$`Source` == "Histone") %>%
    r_to_py()

  # Returns ----

  list("grid" = grid_plot, "complex_list" = final_complex_list,
       "member_list" = final_members, "PTM_list" = final_PTMs, "gene_list" = final_func,
       "net_edges" = edges_py, "net_vertices" = sig_network_py, "ts" = sig_network)
}
