#' Domainness measure
#'
#' @description
#' A function to measure the heterogeneity of a cell's neighbourhood in terms of entropy, such that homogeneous neighbourhoods have low entropy and heterogeneous neighbourhoods have high entropy.
#'
#' @param spe SpatialExperiment object with logcounts, PCA, and 'putative cell type' groups included.
#' @param cells a character vector of cell IDs of each cell. Length of vector must be equal to the number of cells in spatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param regXclust a list of vectors for each cell's neighbourhood composition indicated by the proportion of 'putative cell type' groups it contains.
#' @param threads a numeric value for the number of CPU cores to be used for the analysis.
#'
#' @return SpatialExperiment object including entropy values for each cell neighbourhood.
#'
#' @examples
#' data(mouseEmbryo2)
#'
#' # requires list of containing cluster proportions of each region (regXclust),
#' # generated using the neighbourDetect() function
#' spe <- entropyMeasure(spe, cells = "uniqueID", regXclust, threads = 1)
#' head(spe$entropy)
#'
#' @export

#### Domainness measure
entropyMeasure <- function(spe, cells, regXclust, threads) {
    cellsList <- as.vector(spe[[cells]])
    cl <- parallel::makeCluster(threads)
    doParallel::registerDoParallel(cl)
    regEntropy <- foreach(c = cellsList[1:length(spe[[cells]])]) %dopar% {
        arr <- as.vector(unlist(regXclust[c]))
        # calculate Shannon's entropy
        y <- matrix(-sum(arr * log2(arr)), nrow = 1, ncol = 1)
        rownames(y) <- c
        y
    }
    stopCluster(cl)

    # QC check
    check.cells <- identical(sapply(regEntropy, rownames), spe[[cells]])
    check.NA.values <- sum(is.na(unlist(regEntropy)))
    if (check.cells == FALSE) {
        stop("Order of cells in entropy data does not match cell order in spatial experiment object")
    } else if (check.NA.values != 0) {
        stop("Missing entropy values.")
    } else {
        spe$entropy <- as.matrix(unlist(regEntropy))
        print(paste("Region domainness calculated.", Sys.time()))
    }
    return(spe)
}
