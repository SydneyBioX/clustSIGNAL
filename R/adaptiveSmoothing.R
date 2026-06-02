#' Adaptive smoothing
#'
#' @description
#' A function to perform a weighted, adaptive smoothing of the gene expression
#' of each cell based on the heterogeneity of its neighbourhood. Heterogeneous
#' neighbourhoods are smoothed less with higher weights given to cells belonging
#' to same initial cluster as the index cell. Homogeneous neighbourhoods are
#' smoothed more with similar weights given to most cells.
#' @param spe SpatialExperiment object containing neighbourhood entropy values
#' of each cell.
#' @param nnCells a character matrix of NN nearest neighbours - rows are index
#' cells and columns are their nearest neighbours ranging from closest to
#' farthest neighbour. For sort = TRUE, the neighbours belonging to the same
#' initial cluster as the index cell are moved closer to it.
#' @param NN an integer for the number of neighbouring cells the function
#' should consider. The value must be greater than or equal to 1. Default value
#' is 30.
#' @param kernel a character for type of distribution to be used. The two valid
#' values are "G" or "E" for Gaussian and exponential distributions,
#' respectively. Default value is "G".
#' @param spread a numeric value for distribution spread, represented by
#' standard deviation for Gaussian distribution and rate for exponential
#' distribution. Default value is 0.3 for Gaussian distribution. The recommended
#' value is 5 for exponential distribution.
#' @return SpatialExperiment object including smoothed gene expression as an
#' additional assay.
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assay
#' @importFrom methods as
#' @importFrom Matrix sparseMatrix
#' @examples
#' data(ClustSignal_example)
#'
#' # requires matrix containing NN nearest neighbour cell labels (nnCells),
#' # generated using the neighbourDetect() function
#' spe <- clustSIGNAL::adaptiveSmoothing(spe, nnCells)
#' spe
#' @export

adaptiveSmoothing <- function(spe, nnCells, NN = 30, kernel = c("G", "E"),
                              spread = 0.3) {
    # check kernel name
    kernel <- match.arg(kernel)

    # extracting normalised counts as sparse matrix
    G_mat <- as(logcounts(spe), "dgCMatrix")
    if (kernel == "G") {
        # generate normal distribution weights
        wts <-
            dnorm(sapply(spe$entropy, function(x) seq(0, x, length.out = NN+1)),
                  sd = spread)
    } else if (kernel == "E") {
        # generate exponential distribution weights
        wts <-
            dexp(sapply(spe$entropy, function(x) seq(0, x, length.out = NN+1)),
                 rate = spread)
    }
    # scale weights
    wts <- t(t(wts) / colSums(wts))

    # QC check
    if (!identical(colnames(wts), colnames(spe))) {
        stop("Cells have likely been allocated wrong weights.")
    }

    # Building a sparse weight matrix for all cells
    col_idx <- rep(seq_len(ncol(spe)), each = NN + 1)
    row_idx <- match(as.vector(t(nnCells)), colnames(spe))
    x_val <- as.vector(wts)
    W_mat <- sparseMatrix(i = row_idx,
                          j = col_idx,
                          x = x_val,
                          dims = c(ncol(spe), ncol(spe)),
                          dimnames = list(colnames(spe), colnames(spe)))
    smoothMat <- G_mat %*% W_mat
    colnames(smoothMat) <- colnames(spe)

    # QC checks
    if (!identical(rownames(spe), rownames(smoothMat))) {
        stop("Gene order in smoothed data does not match gene order
             in SpatialExperiment object.")
    }
    if (anyNA(smoothMat)) {
        stop("Smoothed data has missing values.")
    }
    SummarizedExperiment::assay(spe, "smoothed") <- smoothMat
    message(sprintf("%s Smoothing performed. NN = %d, Kernel = %s, Spread = %f",
                    format(Sys.time(), '%H:%M:%S'), NN, kernel, spread))
    return(spe)
}
