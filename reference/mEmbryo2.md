# Mouse Embryo Data

This dataset contains spatial transcriptomics data from 3 mouse embryos,
with 351 genes and a total of 57536 cells. For vignettes, we subset the
data by randomly selecting 5000 cells from embryo 2, excluding cells
that were annotated as 'low quality'.

## Usage

``` r
data(mEmbryo2)
```

## Format

`me_expr` a gene expression matrix with normalised counts, where rows
indicate genes and columns indicate cells. `me_data` a data frame of
cell metadata including cell IDs, sample IDs, cell type annotations, and
x-y coordinates of cells.

## Source

Integration of spatial and single-cell transcriptomic data elucidates
mouse organogenesis, *Nature Biotechnology*, 2021. Webpage:
<https://www.nature.com/articles/s41587-021-01006-2>
