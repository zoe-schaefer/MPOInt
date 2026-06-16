# Generate a multi-level proteome network from total proteome, phosphoproteome, and epiproteome data.

`multi_net` generates an R object ready to generate a network in Python
given three datasets.

## Usage

``` r
multi_net(
  genes_df,
  total_df,
  phospho_df,
  symbol_col,
  fc_col,
  epi_single_df,
  epi_symbol_col,
  epi_fc_col,
  ...
)
```

## Arguments

- genes_df:

  A dataframe containing the exported information from EpiFactors.

- total_df:

  A dataframe containing total proteome data.

- phospho_df:

  A dataframe containing phosphoproteome data.

- symbol_col:

  A column in both `total_df` and `phospho_df` containing the gene
  symbol names for identified proteins.

- fc_col:

  A column in both `total_df` and `phospho_df` containing the
  log2-transformed fold change values.

- epi_single_df:

  A dataframe containing histone proteome data. If using the full
  pipeline, make sure this is the single PTM dataframe and not
  fragments.

- epi_symbol_col:

  A column in `epi_single_df` containing the gene symbol names for
  identified PTMs.

- epi_fc_col:

  A column in `epi_single_df` containing the log2-transformed fold
  change values for identified PTMs.

- ...:

  Optional.

## Value

An object with the following attributes:

- `grid`: A plot displaying the complexes, PTMs, and networks in an
  intersecting grid.

- `complex_list`: A list of the complexes present in the network.

- `member_list`: A list of the complexes, their complex members
  (proteins), and associated PTMs present in the network.

- `PTM_list`: A list of the identified PTMs and their differential
  expression.

- `gene_list`: A subset of the `genes_df` variable containing proteins
  identified.

- `net_edges`: The edges connecting network vertices, which identify
  member-PTM relationships.

- `net_vertices`: The vertices (PTMs and members) in the network.

- `ts`: A larger set of information on the proteins in the network. This
  contains information including statistical significance values, fold
  change values, and proteome source identification.

## Examples

``` r
if (FALSE) { # \dontrun{
multi_net()
} # }
```
