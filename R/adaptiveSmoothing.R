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
#' @return SpatialExperiment object including smoothed gene expression values as another assay.
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assay
#' @importFrom methods show as
#' @importFrom BiocParallel bplapply
#'
#' @examples
#' data(example)
#'
#' # requires matrix containing NN nearest neighbour cell labels (nnCells),
#' # generated using the neighbourDetect() function
#' spe <- adaptiveSmoothing(spe, nnCells, NN = 30, kernel = "G", spread = 0.05)
#' spe
#'
#' @export

#### Smoothing
adaptiveSmoothing <- function(spe, nnCells, NN = 30, kernel = "G", spread = 0.05, threads = 1) {
    ed <- unique(spe$entropy)
    gXc <- as(logcounts(spe), "sparseMatrix")
    if (kernel == "G") {
        # normal distribution weights
        # for each unique entropy value
        weights <- .gauss_kernel(ed, NN + 1, spread)
    } else if (kernel == "E") {
        # exponential distribution weights
        # for each unique entropy value
        weights <- .exp_kernel(ed, NN + 1, spread)
    }

    outMats <- matrix(nrow = nrow(spe), ncol = 0)
    for (val in seq_len(length(ed))) {
        # cells with same entropy value will have same weights
        entCells <- colnames(spe[, spe$entropy == ed[val]])
        e <- paste0("E", ed[val])

        if (length(entCells) == 1) {
            region <- as.vector(nnCells[entCells, ])
            inMat <- as.matrix(gXc[, region])
            tmpMat <- .smoothedData(inMat, weights[, e])
            colnames(tmpMat) <- entCells
        } else {
            inMatList <- list()
            for (x in entCells) {
                # names of central cell + NN cells
                region <- as.vector(nnCells[x, ])
                # gene expression matrix of NN cells
                inMatList[[x]] <- as.matrix(gXc[, region])
            }
            BPPARAM <- .generateBPParam(cores = threads)
            out <- BiocParallel::bplapply(inMatList,
                                         function(imat){
                                             omat <- .smoothedData(mat = imat, weight = weights[, e])
                                             omat
                                         },
                                         BPPARAM = BPPARAM)
            # out <- lapply(inMatList, .smoothedData, weight = weights[, e])
            tmpMat <- matrix(unlist(out),
                             ncol = length(out),
                             dimnames = list(rownames(out[[1]]), sapply(out, names)))
            colnames(tmpMat) <- names(out)
        }
        outMats <- cbind(outMats, tmpMat)
    }
    # rearrange the cells in smoothed data matrix
    smoothMat <- outMats[, as.vector(colnames(spe))]

    # QC check
    check.genes <- identical(rownames(spe), rownames(smoothMat))
    check.NA <- sum(is.na(smoothMat))
    if (check.genes == FALSE) {
        stop("Order of genes in smoothed data does not match gene order in spatial experiment object.")
    } else if (check.NA != 0) {
        stop("Missing values in smoothed data.")
    } else {
        # add data to spatial experiment
        assay(spe, "smoothed") <- as(smoothMat, "sparseMatrix")
        show(paste("Smoothing performed. NN =", NN, "Kernel =", kernel, "Spread =", spread, Sys.time()))
    }
    return(spe)
}

