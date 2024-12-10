#' Adaptive smoothing
#'
#' @description
#' A function to perform a weighted, adaptive smoothing of gene expression based
#' on the heterogeneity of the cell neighbourhood. Heterogeneous neighbourhoods
#' are smoothed less with higher weights given to cells belonging to same initial
#' group. Homogeneous neighbourhoods are smoothed more with similar weights given
#' to most cells.
#'
#' @param spe SpatialExperiment object with logcounts, PCA, 'putative cell type'
#' groups, and entropy outputs included.
#' @param nnCells a character matrix of NN nearest neighbours - rows are cells
#' and columns are their nearest neighbours ranged from closest to farthest
#' neighbour. For sort = TRUE, the neighbours belonging to the same 'putative
#' cell type' group as the cell are moved closer to it.
#' @param NN an integer for the number of neighbourhood cells the function should
#' consider. The value must be greater than or equal to 1. Default value is 30.
#' @param kernel a character for type of distribution to be used. The two valid
#' values are "G" or "E". G for Gaussian distribution, and E for exponential
#' distribution. Default value is "G".
#' @param spread a numeric value for distribution spread, represented by standard
#' deviation for Gaussian distribution and rate for exponential distribution.
#' Default value is 0.05 for Gaussian distribution and 20 for exponential
#' distribution.
#' @param threads a numeric value for the number of CPU cores to be used for the
#' analysis. Default value set to 1.
#' @return SpatialExperiment object including smoothed gene expression values as
#' another assay.
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assay
#' @importFrom methods show as
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr bind_rows
#' @importFrom Matrix sparseMatrix
#'
#' @examples
#' data(example)
#'
#' # requires matrix containing NN nearest neighbour cell labels (nnCells),
#' # generated using the neighbourDetect() function
#' spe <- clustSIGNAL::adaptiveSmoothing(spe, nnCells)
#' spe
#'
#' @export

#### Smoothing
adaptiveSmoothing <- function(spe, nnCells, NN = 30, kernel = "G",
                              spread = 0.05, threads = 1) {
    G_mat <- as(logcounts(spe), "sparseMatrix")
    if (kernel == "G") { # generate normal distribution weights
        wts <- .gauss_kernel(spe$entropy, NN + 1, spread)
    } else if (kernel == "E") { # generate exponential distribution weights
        wts <- .exp_kernel(spe$entropy, NN + 1, spread)}
    colnames(wts) <- colnames(spe)
    # scale weights
    wts <- sweep(wts, 2, colSums(wts), FUN = "/")
    # create a cell by cell weight matrix
    BPPARAM <- .generateBPParam(cores = threads)
    all_w <- BiocParallel::bplapply(colnames(wts), function(x) {
        neigh_c <- nnCells[x, ]
        other_c <- colnames(spe)[!colnames(spe) %in% neigh_c]
        # tmp_df <- rbind(cbind(neigh_c, rep(x, length(neigh_c)), wts[, x]),
        #                 cbind(other_c, rep(x, length(other_c)), 0)) |>
        #     as.data.frame()
        tmp_df <- rbind(cbind(neigh_c, 1, wts[, x]),
                        cbind(other_c, 1, 0)) |>
            as.data.frame()
        colnames(tmp_df) <- c("ci", "c", "wci")
        tmp_df <- tmp_df[match(colnames(spe), tmp_df$ci), ]
        tmp_df$ci <- as.integer(match(tmp_df$ci, colnames(spe)))
        # tmp_df$c <- as.integer(match(tmp_df$c, colnames(spe)))
        tmp_df$c <- as.integer(tmp_df$c)
        tmp_df$wci <- as.numeric(tmp_df$wci)
        tmp_smat <- Matrix::sparseMatrix(i = tmp_df$ci, j = tmp_df$c,
                                         x = tmp_df$wci)
        tmp_smat}, BPPARAM = BPPARAM)
    # W_mat <- dplyr::bind_cols(all_w)
    # W_mat <- rlist::list.cbind(all_w)
    W_mat <- do.call(cbind, all_w)
    # removing unwanted large variable to free memory
    rm(all_w, wts)
    invisible(gc())
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
