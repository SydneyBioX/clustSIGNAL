#' Adaptive smoothing
#'
#' @description
#' A function to perform a weighted, adaptive smoothing of gene expression based on the heterogeneity of the cell neighbourhood. Heterogeneous neighbourhoods are smoothed less with higher weights given to cells belonging to same initial group. Homogeneous neighbourhoods are smoothed more with similar weights given to most cells.
#'
#' @param spe SpatialExperiment object with logcounts, PCA, 'putative cell type' groups, and entropy outputs included.
#' @param nnCells a character matrix of NN nearest neighbours - rows are cells and columns are their nearest neighbours ranged from closest to farthest neighbour. For sort = TRUE, the neighbours belonging to the same 'putative cell type' group as the cell are moved closer to it.
#' @param NN an integer for the number of neighbourhood cells the function should consider. The value must be greater than or equal to 1. Default value is 30.
#' @param kernel a character for type of distribution to be used. The two valid values are "G" or "E". G for Gaussian distribution, and E for exponential distribution. Default value is "G".
#' @param spread a numeric value for distribution spread, represented by standard deviation for Gaussian distribution and rate for exponential distribution. Default value is 0.05 for Gaussian distribution and 20 for exponential distribution.
#' @param cells a character vector of cell IDs of each cell. Length of vector must be equal to the number of cells in spatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param threads a numeric value for the number of CPU cores to be used for the analysis.
#'
#' @return SpatialExperiment object including smoothed gene expression values as another assay.
#'
#' @examples
#'
#' @export

#### Smoothing
adaptiveSmoothing <- function(spe, nnCells, NN, kernel, spread, cells, threads) {
    gXc <- as(logcounts(spe), "sparseMatrix")
    cl <- makeCluster(threads)
    doParallel::registerDoParallel(cl)
    # create a weights by entropy matrix
    entropies <- unique(spe$entropy)
    if (kernel == "G"){
        e <- foreach (val = c(1:length(entropies))) %dopar% {
            # normal distribution-weights
            weight <- .gauss_kernel(entropies[val], NN + 1, spread)
            colnames(weight) <- paste0("E", entropies[val])
            weight
        }
        emat <- matrix(unlist(e), ncol = length(e), dimnames = list(c(1:length(e[[1]])), sapply(e, colnames)))
    } else if (kernel == "E") {
        e <- foreach (val = c(1:length(entropies))) %dopar% {
            # exponential distribution-weights
            weight <- .exp_kernel(entropies[val], NN + 1, spread)
            colnames(weight) <- paste0("E", entropies[val])
            weight
        }
        emat <- matrix(unlist(e), ncol = length(e), dimnames = list(c(1:length(e[[1]])), sapply(e, colnames)))
    }
    # loop through each cell column
    y <- foreach (x = c(1:ncol(gXc))) %dopar% {
        # central cell name
        cell <- colnames(gXc)[x]
        # cell entropy
        ed <- paste0("E", spe$entropy[x])
        # names of central cell + NN cells
        region <- as.vector(nnCells[x, ])
        # gene expression matrix of NN cells
        inMat <- as.matrix(gXc[, region])
        smat <- .smoothedData(inMat, emat[, ed])
        colnames(smat) <- cell
        smat
    }
    stopCluster(cl)

    # QC check
    check.cells <- identical(sapply(y, colnames), cells)
    check.genes <- identical(rownames(logcounts(spe)), rownames(y[[1]]))
    check.NA <- sum(is.na(unlist(y)))
    if (check.cells == FALSE) {
        stop("ERROR: Order of cells in smoothed data does not match cell order in spatial experiment object.")
    } else if (check.genes == FALSE) {
        stop("ERROR: Order of genes in smoothed data does not match gene order in spatial experiment object.")
    } else if (check.NA != 0) {
        stop("ERROR: Missing values in smoothed data.")
    } else {
        smoothMat <- matrix(unlist(y), ncol = length(y), dimnames = list(rownames(y[[1]]), sapply(y, colnames)))
        # add data to spatial experiment
        assay(spe, i = "smoothed") <- as(smoothMat, "sparseMatrix")
        print(paste("Smoothing performed. NN =", NN, "Kernel =", kernel, "Spread =", spread, Sys.time()))
    }
    return(spe)
}
