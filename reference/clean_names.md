# Clean up names from EpiProfile.

`clean_names` takes a dataframe with histone mass spectrometry data and
outputs two non-transformed dataframes with user-friendly PTM names.

## Usage

``` r
clean_names(df, label_col, ...)
```

## Arguments

- df:

  A dataframe containing histone mass spectrometry data.

- label_col:

  A column containing PTM/fragment names.

- ...:

  Optional.

## Value

Two dataframes with adjusted labels, up to 4 co-occurring PTMs per
fragment.

## Examples

``` r
if (FALSE) { # \dontrun{
clean_names(df = THP1_histones, label_col = "PTM")
} # }
```
