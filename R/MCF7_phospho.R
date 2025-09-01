#' MCF-7 phosphoproteome dataset.
#'
#' Phosphoproteome data from Yang et al. (2024). M90 condition: MCF-7 cells grown in 90 Pa matrix. F condition (control): MCF-7 cells grown on flat culture surface.
#'
#' @format A data frame with 4,172 rows and 13 variables:
#' \describe{
#' \item{Accession}{UniProt accession numbers.}
#' \item{Position}{Amino acid number for localization.}
#' \item{Amino acid}{Amino acid at that position.}
#' \item{Protein description}{Description of each identified protein.}
#' \item{Localization probility}{Localization probability for specific residues.}
#' \item{F1}{Raw data for F1 condition.}
#' \item{F2}{Raw data for F2 condition.}
#' \item{F3}{Raw data for F3 condition.}
#' \item{M10501}{Raw data for M1050-1 condition.}
#' \item{M10502}{Raw data for M1050-2 condition.}
#' \item{M10503}{Raw data for M1050-3 condition.}
#' \item{M4501}{Raw data for M450-1 condition.}
#' \item{M4502}{Raw data for M450-2 condition.}
#' \item{M4503}{Raw data for M450-3 condition.}
#' \item{M901}{Raw data for M90-1 condition.}
#' \item{M902}{Raw data for M90-2 condition.}
#' \item{M903}{Raw data for M90-3 condition.}
#' \item{P value Anova}{Calculated unadjusted P value for fold change of M90/F cells.}
#' \item{M1050/F Ratio}{Raw ratio for M1050/F fold change calculation.}
#' \item{M1050/F P value}{P value for M1050/F fold change calculation.}
#' \item{M450/F Ratio}{Raw ratio for M450/F fold change calculation.}
#' \item{M450/F P value}{P value for M450/F fold change calculation.}
#' \item{M90/F Ratio}{Raw ratio for M90/F fold change calculation.}
#' \item{P value}{P value for M90/F fold change calculation.}
#' \item{M1050/M450 Ratio}{Raw ratio for M1050/M450 fold change calculation.}
#' \item{M1050/M450 P value}{P value for M1050/M450 fold change calculation.}
#' \item{log2(FC)}{The log2 transformed fold change of M90/F cells.}
#' \item{ENTREZID}{ENTREZ ID for pathway analysis.}
#' \item{SYMBOL2}{Common gene symbols.}
#' }
#'
#' @source \url{https://doi.org/10.1021/acs.jproteome.4c00563}
"MCF7_phospho"
