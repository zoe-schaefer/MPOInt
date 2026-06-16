# THP-1 phosphoproteome dataset.

Phosphoproteome data from Agarwal et al. (2021). Infected condition:
THP-1 macrophages exposed to M. bovis BCG for 48 hours at a multiplicity
of 10.

## Usage

``` r
THP1_phospho
```

## Format

A data frame with 6,121 rows and 4 variables:

- log2(FC):

  The log2 transformed fold change of infected/control cells.

- P value:

  Calculated unadjusted P value.

- ENTREZID:

  ENTREZ ID for pathway analysis.

- SYMBOL2:

  Common gene symbols.

- UNIPROT:

  UNIPROT IDs for the identified proteins.

## Source

<https://doi.org/10.1021/acs.jproteome.9b00895>
