#' Mouse Embryo Data as SpatialExperiment object
#'
#' This dataset contains spatial transcriptomics data from 3 mouse embryos, with
#' 351 genes and a total of 57536 cells. For running vignettes and examples, we subset
#' the data by selecting only embryo 2 and removed all cells that were annotated
#' as 'low quality'. After subsetting, we have 14,185 cells from embryo 2 and 351
#' genes.
#'
#'
#' @name mEmbryo2
#' @aliases nnCells me_data me_expr regXclust
#' @docType data
#' @format
#' \code{me_expr} a gene expression matrix with normalised counts, where rows indicate
#' genes and columns indicate cells.
#' \code{me_data} a data frame of cell metadata including cell IDs, sample IDs,
#' cell type annotations, and x-y coordinates of cells.
#' \code{nnCells} a matrix where each row corresponds to a cell in spe object,
#' and the columns correspond to the nearest neighbors.
#' \code{regXclust} a list where each element corresponds to a cell in spe object,
#' and contains the cluster composition proportions.
#' @usage load("mEmbryo2.RData")
#' @source Integration of spatial and single-cell transcriptomic data elucidates mouse
#' organogenesis, \emph{Nature Biotechnology}, 2022.
#' Webpage: \url{https://www.nature.com/articles/s41587-021-01006-2}
#' @keywords datasets
NULL


#' Mouse Hypothalamus Data as SpatialExperiment object
#'
#' This dataset contains spatial transcriptomics data from 181 mouse hypothalamus
#' samples embryos, 155 genes and a total of 1,027,080 cells. For running the
#' vignettes, we subset the data by selecting only 3 samples - Animal 1 Bregma -0.09
#' and Animal 7 Bregmas 0.16 and -0.09, removed all cells that were annotated
#' as 'ambiguous', and removed 20 genes that were assessed using a different technology.
#' After subsetting, we have 15,848 cells from 3 mouse brain samples and 135 genes.
#'
#'
#' @name mHypothal
#' @aliases mh_data mh_expr
#' @docType data
#' @format
#' \code{mh_expr} a gene expression matrix with normalised counts, where rows indicate
#' genes and columns indicate cells.
#' \code{mh_data} a data frame of cell metadata including cell IDs, sample IDs,
#' cell type annotations, and x-y coordinates of cells.
#' @usage load("mHypothal.RData")
#' @source Molecular, Spatial and Functional Single-Cell Profiling of the
#' Hypothalamic Preoptic Region, \emph{Science}, 2018.
#' Webpage: \url{https://www.science.org/doi/10.1126/science.aau5324}
#' @keywords datasets
NULL
