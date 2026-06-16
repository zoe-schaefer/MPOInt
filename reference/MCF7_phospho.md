# MCF-7 phosphoproteome dataset.

Phosphoproteome data from Yang et al. (2024). M90 condition: MCF-7 cells
grown in 90 Pa matrix. F condition (control): MCF-7 cells grown on flat
culture surface.

## Usage

``` r
MCF7_phospho
```

## Format

A data frame with 4,172 rows and 13 variables:

- Accession:

  UniProt accession numbers.

- Position:

  Amino acid number for localization.

- Amino acid:

  Amino acid at that position.

- Protein description:

  Description of each identified protein.

- Localization probility:

  Localization probability for specific residues.

- F1:

  Raw data for F1 condition.

- F2:

  Raw data for F2 condition.

- F3:

  Raw data for F3 condition.

- M10501:

  Raw data for M1050-1 condition.

- M10502:

  Raw data for M1050-2 condition.

- M10503:

  Raw data for M1050-3 condition.

- M4501:

  Raw data for M450-1 condition.

- M4502:

  Raw data for M450-2 condition.

- M4503:

  Raw data for M450-3 condition.

- M901:

  Raw data for M90-1 condition.

- M902:

  Raw data for M90-2 condition.

- M903:

  Raw data for M90-3 condition.

- P value Anova:

  Calculated unadjusted P value for fold change of M90/F cells.

- M1050/F Ratio:

  Raw ratio for M1050/F fold change calculation.

- M1050/F P value:

  P value for M1050/F fold change calculation.

- M450/F Ratio:

  Raw ratio for M450/F fold change calculation.

- M450/F P value:

  P value for M450/F fold change calculation.

- M90/F Ratio:

  Raw ratio for M90/F fold change calculation.

- P value:

  P value for M90/F fold change calculation.

- M1050/M450 Ratio:

  Raw ratio for M1050/M450 fold change calculation.

- M1050/M450 P value:

  P value for M1050/M450 fold change calculation.

- log2(FC):

  The log2 transformed fold change of M90/F cells.

- ENTREZID:

  ENTREZ ID for pathway analysis.

- SYMBOL2:

  Common gene symbols.

## Source

<https://doi.org/10.1021/acs.jproteome.4c00563>
