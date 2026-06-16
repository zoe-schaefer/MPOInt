# Generate a volcano plot with groups of interest.

`volcano_groups` generates a volcano plot given data, threshold values,
and a subset of points to label.

## Usage

``` r
volcano_groups(
  df,
  fc_col,
  p_col,
  fc_thr = 0.5,
  p_thr = 0.05,
  label_col = NA,
  compute_fc_log2 = FALSE,
  group_list = NULL,
  ...
)
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

  A column in `df` with values to be used if groups are not provided.

- compute_fc_log2:

  A Boolean to determine if the `fc_col` values should be additionally
  log2-transformed. Default = FALSE.

- group_list:

  A named list of lists where each list is a group of peptides, labeled
  with the group name.

- ...:

  Optional.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
volcano_groups(df = THP1_total, fc_col = "log2(FC)", p_col = "P value",
  label_col = "SYMBOL2", group_list = list("Group 1" =
   c("MORF4L1", "SAMSN1", "DTX3L", "MBD3",  "HSF1",
  "EYA3", "NCAPD2", "TAF9", "PPP4C", "ZNF638"),
  "Group 2" = c("CD9", "CD82", "NFKB1", "BST2", "CTSH"))) +
  ggtitle("Grouped volcano plot") +
  xlab("log2FC") + ylab("-log10P")
} # }
```
