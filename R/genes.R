#' EpiFactors gene set.
#'
#' Exported information from the EpiFactors database on epigenetic actors, complexes, and targets.
#'
#' @format A data frame with 801 rows and 7 variables:
#' \describe{
#' \item{HGNC approved symbol}{Common symbols for each member.}
#' \item{UniProt ID (human)}{Corresponding UniProt IDs.}
#' \item{Function}{Assigned function.}
#' \item{Modification}{Resulting modification.}
#' \item{Protein complex}{Complex(es) with this effector as a member.}
#' \item{Target entity}{Modification or residue targeted by each complex.}
#' \item{Product}{Resulting modification or residue.}
#' }
#'
#' @source \url{https://doi.org/10.1093/nar/gkac989}
"genes"
