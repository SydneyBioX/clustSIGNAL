#' Mouse Embryo Data
#'
#' This dataset contains spatial transcriptomics data from 3 mouse embryos, with
#' 351 genes and a total of 57536 cells. For vignettes, we subset the data by
#' randomly selecting 5000 cells from embryo 2, excluding cells that were
#' annotated as 'low quality'.
#'
#' @name mEmbryo2
#' @aliases me_data me_expr
#' @docType data
#' @format
#' \code{me_expr} a gene expression matrix with normalised counts, where rows
#' indicate genes and columns indicate cells.
#' \code{me_data} a data frame of cell metadata including cell IDs, sample IDs,
#' cell type annotations, and x-y coordinates of cells.
#' @usage data(mEmbryo2)
#' @source Integration of spatial and single-cell transcriptomic data elucidates
#' mouse organogenesis, \emph{Nature Biotechnology}, 2021.
#' Webpage: \url{https://www.nature.com/articles/s41587-021-01006-2}
#' @keywords datasets
NULL


#' Example data with SpatialExperiment object
#'
#' This example data was generated from the mouse embryo spatial transcriptomics
#' dataset of 3 mouse embryos, with 351 genes and a total of 57536 cells. For
#' running examples, we subset the data by selecting 1000 random cells from
#' embryo 2, excluding any cells annotated as 'low quality'. After subsetting,
#' we have expression for 351 genes from 1000 cells in embryo 2.
#'
#' @name ClustSignal_example
#' @aliases spe nnCells regXclust
#' @docType data
#' @format
#' \code{spe} a spatialExperiment object containing gene expression matrix with
#' normalised counts, where rows indicate genes and columns indicate cells.
#' Also, contains a cell metadata including cell IDs, sample IDs, cell type
#' annotations, and x-y coordinates of cells.
#' \code{nnCells} a matrix where each row corresponds to a cell in spe object,
#' and the columns correspond to the nearest neighbors.
#' \code{regXclust} a list where each element corresponds to a cell in spe
#' object, and contains the cluster composition proportions.
#' @usage data(ClustSignal_example)
#' @source Integration of spatial and single-cell transcriptomic data elucidates
#' mouse organogenesis, \emph{Nature Biotechnology}, 2021.
#' Webpage: \url{https://www.nature.com/articles/s41587-021-01006-2}
#' @keywords datasets
NULL


#' Mouse Hypothalamus Data
#'
#' This dataset contains spatial transcriptomics data from 181 mouse
#' hypothalamus samples, 155 genes and a total of 1,027,080 cells. For running
#' the vignettes, we subset the data by selecting total 6000 cells from only 3
#' samples - Animal 1 Bregma -0.09 (2080 cells) and Animal 7 Bregmas 0.16
#' (1936 cells) and -0.09 (1984 cells), excluding cells that were annotated as
#' 'ambiguous', and removed 20 genes that were assessed using a different
#' technology.
#'
#' @name mHypothal
#' @aliases mh_data mh_expr
#' @docType data
#' @format
#' \code{mh_expr} a gene expression matrix with normalised counts, where rows
#' indicate genes and columns indicate cells.
#' \code{mh_data} a data frame of cell metadata including cell IDs, sample IDs,
#' cell type annotations, and x-y coordinates of cells.
#' @usage data(mHypothal)
#' @source Molecular, Spatial and Functional Single-Cell Profiling of the
#' Hypothalamic Preoptic Region, \emph{Science}, 2018.
#' Webpage: \url{https://www.science.org/doi/10.1126/science.aau5324}
#' @keywords datasets
NULL
