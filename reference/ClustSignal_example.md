# Example data with SpatialExperiment object

This example data was generated from the mouse embryo spatial
transcriptomics dataset of 3 mouse embryos, with 351 genes and a total
of 57536 cells. For running examples, we subset the data by selecting
1000 random cells from embryo 2, excluding any cells annotated as 'low
quality'. After subsetting, we have expression for 351 genes from 1000
cells in embryo 2.

## Usage

``` r
data(ClustSignal_example)
```

## Format

`spe` a spatialExperiment object containing gene expression matrix with
normalised counts, where rows indicate genes and columns indicate cells.
Also, contains a cell metadata including cell IDs, sample IDs, cell type
annotations, and x-y coordinates of cells. `nnCells` a matrix where each
row corresponds to a cell in spe object, and the columns correspond to
the nearest neighbors. `regXclust` a list where each element corresponds
to a cell in spe object, and contains the cluster composition
proportions.

## Source

Integration of spatial and single-cell transcriptomic data elucidates
mouse organogenesis, *Nature Biotechnology*, 2021. Webpage:
<https://www.nature.com/articles/s41587-021-01006-2>
