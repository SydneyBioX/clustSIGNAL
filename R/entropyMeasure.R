#' Domainness measure
#'
#' @description
#' A function to measure the heterogeneity of a cell's neighbourhood in terms of
#' entropy, such that homogeneous neighbourhoods have low entropy and
#' heterogeneous neighbourhoods have high entropy.
#'
#' @param spe SpatialExperiment object with logcounts, PCA, and 'putative cell
#' type' groups included.
#' @param regXclust a list of vectors for each cell's neighbourhood composition
#' indicated by the proportion of 'putative cell type' groups it contains.
#' @param threads a numeric value for the number of CPU cores to be used for the
#' analysis.
#'
#' @return SpatialExperiment object including entropy values for each cell
#' neighbourhood.
#'
#' @importFrom BiocParallel bplapply
#' @importFrom methods show
#'
#' @examples
#' data(example)
#'
#' # requires list containing cluster proportions of each region (regXclust),
#' # generated using the neighbourDetect() function
#' spe <- clustSIGNAL::entropyMeasure(spe, regXclust)
#' head(spe$entropy)
#'
#' @export

#### Domainness measure
entropyMeasure <- function(spe, regXclust, threads = 1) {
    cell_vect <- as.vector(colnames(spe))
    BPPARAM <- .generateBPParam(cores = threads)
    regEntropy <- BiocParallel::bplapply(cell_vect, function(c){
        arr <- as.vector(unlist(regXclust[c]))
        # calculate Shannon's entropy
        y <- round(matrix(-sum(arr * log2(arr)), nrow = 1, ncol = 1), 5)
        rownames(y) <- c
        y}, BPPARAM = BPPARAM)

    # QC check
    check.cells <- identical(sapply(regEntropy, rownames), colnames(spe))
    check.NA.values <- sum(is.na(unlist(regEntropy)))
    if (check.cells == FALSE) {
        stop("Cell order in entropy data does not match cell order in Spatial Experiment object")
    } else if (check.NA.values != 0) {
        stop("Entropy values are missing for some cells.")
    } else {
        spe$entropy <- as.vector(unlist(regEntropy))
        show(paste("Region domainness calculated. Time",
                   format(Sys.time(),'%H:%M:%S')))}
    return(spe)
}
