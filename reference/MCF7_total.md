# MCF-7 total proteome dataset.

Total proteome data from Yang et al. (2024). M90 condition: MCF-7 cells
grown in 90 Pa matrix. F condition (control): MCF-7 cells grown on flat
culture surface.

## Usage

``` r
MCF7_total
```

## Format

A data frame with 8,317 rows and 6 variables:

- Accession:

  UniProt accession numbers.

- M90/F Ratio:

  Raw ratio for fold change calculation.

- log2(FC):

  The log2 transformed fold change of M90/F cells.

- P value:

  Calculated unadjusted P value.

- ENTREZID:

  ENTREZ ID for pathway analysis.

- SYMBOL2:

  Common gene symbols.

## Source

<https://doi.org/10.1021/acs.jproteome.4c00563>
