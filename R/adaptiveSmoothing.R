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
#' @importFrom rlist list.cbind
#'
#' @examples
#' data(example)
#'
#' # requires matrix containing NN nearest neighbour cell labels (nnCells),
#' # generated using the neighbourDetect() function
#' spe <- clustSIGNAL::adaptiveSmoothing(spe, nnCells, NN = 30, kernel = "G",
#'                                       spread = 0.05, threads = 1)
#' spe
#'
#' @export

#### Smoothing
adaptiveSmoothing <- function(spe, nnCells, NN, kernel, spread, threads) {
    ed <- unique(spe$entropy)
    gXc <- as(logcounts(spe), "sparseMatrix")
    if (kernel == "G") {
        # normal distribution weights
        weights <- .gauss_kernel(ed, NN + 1, spread)
    } else if (kernel == "E") {
        # exponential distribution weights
        weights <- .exp_kernel(ed, NN + 1, spread)}

    outMats <- matrix(nrow = nrow(spe), ncol = 0)
    BPPARAM <- .generateBPParam(cores = threads)
    out_main <- BiocParallel::bplapply(ed,
                                       function(val){
                                           entCells <- colnames(spe[, spe$entropy == val])
                                           # show(paste(val, length(entCells)))
                                           e <- paste0("E", val)
                                           if (length(entCells) == 1) {
                                               region <- as.vector(nnCells[entCells, ])
                                               inMat <- as.matrix(gXc[, region])
                                               tmpMat <- .smoothedData(inMat, weights[, e])
                                               colnames(tmpMat) <- entCells
                                               tmpMat
                                           } else {
                                               out_sub <- lapply(entCells, function(x){
                                                   # names of central cell + NN cells
                                                   region <- as.vector(nnCells[x, ])
                                                   # gene expression matrix of NN cells
                                                   inMat <- as.matrix(gXc[, region])
                                                   omat <- .smoothedData(mat = inMat,
                                                                         weight = weights[, e])
                                                   omat})
                                               tmpMat <- matrix(unlist(out_sub),
                                                                ncol = length(out_sub),
                                                                dimnames = list(rownames(out_sub[[1]]),
                                                                                entCells))
                                               tmpMat}
                                           outMats <- cbind(outMats, tmpMat)
                                           outMats
                                       },
                                       BPPARAM = BPPARAM)
    out_main_comb <- rlist::list.cbind(out_main)
    # rearrange the cells in smoothed data matrix
    smoothMat <- out_main_comb[, as.vector(colnames(spe))]

    # QC check
    check.genes <- identical(rownames(spe), rownames(smoothMat))
    check.NA <- sum(is.na(smoothMat))
    if (check.genes == FALSE) {
        stop("Gene order in smoothed data does not match gene order in SpatialExperiment object.")
    } else if (check.NA != 0) {
        stop("Smoothed data has missing values.")
    } else {
        # add data to spatial experiment
        SummarizedExperiment::assay(spe, "smoothed") <- as(smoothMat, "sparseMatrix")
        show(paste("Smoothing performed. NN =", NN, "Kernel =", kernel,
                   "Spread =", spread, "Time", format(Sys.time(),'%H:%M:%S')))}
    return(spe)
}
