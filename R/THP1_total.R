#' THP-1 total proteome dataset.
#'
#' Total proteome data from Schaefer et al. (2024). Infected condition: THP-1 macrophages exposed to M. bovis BCG for 24 hours at a multiplicity of 4.
#'
#' @format A data frame with 4,252 rows and 11 variables:
#' \describe{
#' \item{Accession}{UniProt accession numbers.}
#' \item{Abundance Ratio (log2): (JMI) / (Control)}{The log2 transformed fold change of infected/control cells.}
#' \item{Abundance Ratio Adj. P-Value: (JMI) / (Control)}{Calculated adjusted P value.}
#' \item{Found in Sample: F1: Sample Control}{Peak detection for each protein in sample F1.}
#' \item{Found in Sample: F2: Sample Control}{Peak detection for each protein in sample F2.}
#' \item{Found in Sample: F3: Sample Control}{Peak detection for each protein in sample F3.}
#' \item{Found in Sample: F4: Sample JMI}{Peak detection for each protein in sample F4.}
#' \item{Found in Sample: F5: Sample JMI}{Peak detection for each protein in sample F5.}
#' \item{Found in Sample: F6: Sample JMI}{Peak detection for each protein in sample F6.}
#' \item{Found in Sample: F7: Sample JMI}{Peak detection for each protein in sample F7.}
#' \item{Found in Sample: F8: Sample JMI}{Peak detection for each protein in sample F8.}
#' }
#'
#' @source \url{https://doi.org/10.1016/j.mcpro.2024.100851}
"THP1_total"
