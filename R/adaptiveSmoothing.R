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
#' @importFrom methods show as
#' @importFrom BiocParallel bplapply
#' @importFrom Matrix sparseMatrix
#' @importFrom reshape2 melt
#' @examples
#' data(ClustSignal_example)
#'
#' # requires matrix containing NN nearest neighbour cell labels (nnCells),
#' # generated using the neighbourDetect() function
#' spe <- clustSIGNAL::adaptiveSmoothing(spe, nnCells)
#' spe
#' @export

adaptiveSmoothing <- function(spe, nnCells, NN = 30, kernel = "G",
                              spread = 0.3) {
    G_mat <- as(logcounts(spe), "sparseMatrix")
    if (kernel == "G") { # generate normal distribution weights
        wts <- .gauss_kernel(spe$entropy, NN + 1, spread)
    } else if (kernel == "E") { # generate exponential distribution weights
        wts <- .exp_kernel(spe$entropy, NN + 1, spread)}
    colnames(wts) <- colnames(spe)
    wts <- sweep(wts, 2, colSums(wts), FUN = "/") # scale weights
    wts_df <- cbind(reshape2::melt(t(nnCells))[, 2:3], reshape2::melt(wts)[, 3])
    colnames(wts_df) <- c("c", "ci", "wci")
    wts_df$c <- as.integer(match(wts_df$c, colnames(spe)))
    wts_df$ci <- as.integer(match(wts_df$ci, colnames(spe)))
    W_mat <- Matrix::sparseMatrix(i = wts_df$ci, j = wts_df$c,
                                  x = wts_df$wci)
    smoothMat <- G_mat %*% W_mat
    colnames(smoothMat) <- colnames(spe)
    # QC check
    check.genes <- identical(rownames(spe), rownames(smoothMat))
    check.NA <- sum(is.na(smoothMat))
    if (check.genes == FALSE) {
        stop("Gene order in smoothed data does not match gene order
             in SpatialExperiment object.")
    } else if (check.NA != 0) {
        stop("Smoothed data has missing values.")
    } else {
        SummarizedExperiment::assay(spe, "smoothed") <- smoothMat
        show(paste("Smoothing performed. NN =", NN, "Kernel =", kernel,
                   "Spread =", spread, "Time", format(Sys.time(),'%H:%M:%S')))}
    return(spe)
}
