# Generate a volcano plot for epiproteome data.

`volcano_epi` generates a volcano plot given data and threshold values.

## Usage

``` r
volcano_epi(df, fc_col, p_col, fc_thr = 0.5, p_thr = 0.05, label_col = NA, ...)
```

## Arguments

- df:

  A data frame containing points to be plotted.

- fc_col:

  A column in `df` with fold change values to be plotted.

- p_col:

  A column in `df` with P-values to be plotted. These values should not
  be log-transformed.

- fc_thr:

  A scalar value indicating the fold change threshold. Default = 0.5.

- p_thr:

  A scalar value indicating the P-value threshold. Default = 0.05.

- label_col:

  A column in `df` with which points should be labeled. If
  `label_col = NA`, default labels will be applied to significant
  points.

- ...:

  Optional.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
volcano_epi(df = total_raw, fc_col = "FC", p_col = "Pval", label_col = "SYMBOL2") +
ggtitle("Volcano plot") + xlab("log2FC") + ylab("-log10P")
} # }
```
