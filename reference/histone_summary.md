# Process summary data for histone proteome.

`histone_summary` takes a dataframe with histone mass spectrometry data
and outputs a log2-transformed fold change as well as P-values and a T
score.

## Usage

``` r
histone_summary(
  df,
  label_col,
  base_cols,
  base_std = NA,
  target_col,
  target_std,
  ...
)
```

## Arguments

- df:

  A dataframe containing histone mass spectrometry data. This dataframe
  should include summary statistics (average value, standard deviation)
  for at least two conditions.

- label_col:

  A column containing PTM/fragment names.

- base_cols:

  A list of at least one column from which to generate a baseline.

- base_std:

  If more than one column is supplied to `base_cols`, leave this as NA.
  Otherwise, if only one column is supplied, provide the standard
  deviation for that column.

- target_col:

  A column containing average values for your condition of interest.

- target_std:

  A column containing standard deviations for your condition of
  interest.

- ...:

  Optional.

## Value

A dataframe with labels, log2(FC), T scores, and P-values

## Examples

``` r
if (FALSE) { # \dontrun{
histone_summary(df = MCF7_histones, label_col = "PTM",
base_cols = c("293 ave (5 reps)", "HaCAT ave (5 reps)", "hESC ave (6 reps)",
"HFF ave (7 reps)", "Mdm13 ave (5 reps)"), base_std = NA,
target_col = "MCF7 ave (5 reps)", target_std = "MCF7 std")
} # }
```
