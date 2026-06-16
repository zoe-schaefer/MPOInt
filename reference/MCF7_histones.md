# MCF-7 histone proteome dataset.

Histone proteome data from Garcia et al. (2013).

## Usage

``` r
MCF7_histones
```

## Format

A data frame with 37 rows and 12 variables:

- PTM:

  PTM descriptions.

- MCF7 std:

  Standard deviation in MCF-7 cells.

- MCF7 ave (5 reps):

  Average expression in MCF-7 cells.

- 293 ave (5 reps):

  Average expression in HEK293 cells.

- HaCAT ave (5 reps):

  Average expression in HaCAT cells.

- hESC ave (6 reps):

  Average expression in hESC cells.

- HFF ave (7 reps):

  Average expression in HFF cells.

- Mdm13 ave (5 reps):

  Average expression in Mdm13 cells.

- ctrl ave:

  Average expression in control cell types (HEK293, HaCAT, hESC, HFF,
  Mdm13).

- log2(FC):

  The log2 transformed fold change in expression between MCF-7 and
  control cells.

- T score:

  The calculated T score using standard deviation.

- P value:

  The corresponding P value.

## Source

<https://doi.org/10.1186/1756-8935-6-20>
