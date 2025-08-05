#' MCF-7 total proteome dataset.
#'
#' Total proteome data from Yang et al. (2024). M90 condition: MCF-7 cells grown in 90 Pa matrix. F condition (control): MCF-7 cells grown on flat culture surface.
#'
#' @format A data frame with 8,317 rows and 6 variables:
#' \describe{
#' \item{Accession}{UniProt accession numbers.}
#' \item{M90/F Ratio}{Raw ratio for fold change calculation.}
#' \item{log2(FC)}{The log2 transformed fold change of M90/F cells.}
#' \item{P value}{Calculated unadjusted P value.}
#' \item{ENTREZID}{ENTREZ ID for pathway analysis.}
#' \item{SYMBOL2}{Common gene symbols.}
#' }
#'
#' @source \url{https://doi.org/10.1021/acs.jproteome.4c00563}
"MCF7_total"
