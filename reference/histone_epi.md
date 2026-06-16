# Process individual data for histone proteome.

`histone_epi` takes a dataframe with histone mass spectrometry data and
outputs a log2-transformed fold change as well as P-values and a T
score.

## Usage

``` r
histone_epi(df, label_col, base_cols, target_cols, ...)
```

## Arguments

- df:

  A dataframe containing histone mass spectrometry data. This dataframe
  should include measurements for two conditions.

- label_col:

  A column containing PTM/fragment names.

- base_cols:

  A list of columns containing baseline values.

- target_cols:

  A list of columns containing average values for your condition of
  interest.

- ...:

  Optional.

## Value

A dataframe with labels, log2(FC), T scores, and P-values

## Examples

``` r
if (FALSE) { # \dontrun{
histone_epi(df = THP1_single_raw, label_col = "PTM",
base_cols = c("C1", "C2", "C3", "C4", "C5"),
target_cols = c("I1", "I2", "I3", "I4", "I5"))
} # }
```
