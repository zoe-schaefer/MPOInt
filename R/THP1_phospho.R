#' THP-1 phosphoproteome dataset.
#'
#' Phosphoproteome data from Agarwal et al. (2021). Infected condition: THP-1 macrophages exposed to M. bovis BCG for 48 hours at a multiplicity of 10.
#'
#' @format A data frame with 6,121 rows and 4 variables:
#' \describe{
#' \item{log2(FC)}{The log2 transformed fold change of infected/control cells.}
#' \item{P value}{Calculated unadjusted P value.}
#' \item{ENTREZID}{ENTREZ ID for pathway analysis.}
#' \item{SYMBOL}{Common gene symbols.}
#' }
#'
#' @source \url{https://doi.org/10.1021/acs.jproteome.9b00895}
"THP1_phospho"
