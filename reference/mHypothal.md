# Mouse Hypothalamus Data

This dataset contains spatial transcriptomics data from 181 mouse
hypothalamus samples, 155 genes and a total of 1,027,080 cells. For
running the vignettes, we subset the data by selecting total 6000 cells
from only 3 samples - Animal 1 Bregma -0.09 (2080 cells) and Animal 7
Bregmas 0.16 (1936 cells) and -0.09 (1984 cells), excluding cells that
were annotated as 'ambiguous', and removed 20 genes that were assessed
using a different technology.

## Usage

``` r
data(mHypothal)
```

## Format

`mh_expr` a gene expression matrix with normalised counts, where rows
indicate genes and columns indicate cells. `mh_data` a data frame of
cell metadata including cell IDs, sample IDs, cell type annotations, and
x-y coordinates of cells.

## Source

Molecular, Spatial and Functional Single-Cell Profiling of the
Hypothalamic Preoptic Region, *Science*, 2018. Webpage:
<https://www.science.org/doi/10.1126/science.aau5324>
