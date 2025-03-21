#' Generating BPParam object
#'
#' @description
#' A utility function to generate BPPARAM object.
#' @param cores Desired number of cores for BPPARAM object.
#' @importFrom BiocParallel SerialParam SnowParam MulticoreParam bpparam
#' @return A BPPPARAM object.
#' @keywords internal

.generateBPParam <- function(cores = 1) {
    seed <- .Random.seed[1]
    if (cores == 1) {
        BPparam <- BiocParallel::SerialParam(RNGseed = seed)
    } else { # Parallel processing is desired.
        # Set BPparam RNGseed if the user ran set.seed(someNumber) themselves
        if (Sys.info()["sysname"] == "Windows") {
            # Only SnowParam suits Windows
            BPparam <- BiocParallel::SnowParam(min(
                cores,
                BiocParallel::snowWorkers("SOCK")
            ),
            RNGseed = seed
            )
        } else if (Sys.info()["sysname"] %in% c("MacOS", "Linux")) {
            BPparam <- BiocParallel::MulticoreParam(min(
                cores,
                BiocParallel::multicoreWorkers()
            ),
            RNGseed = seed
            )
            # Multicore is faster than SNOW, but it does not work on Windows
        } else {
            BPparam <- BiocParallel::bpparam()}}
    BPparam
}


#' Neighbour cell sorting
#'
#' @description
#' A function to perform neighbourhood cell sorting. Neighbourhood cells that
#' belong to the same initial cluster as the index cell are moved closer to the
#' index cell.
#' @param cell a vector of neighbourhood cell indices. The cell indices indicate
#' the row number of cells in sample metadata.
#' @param Clust a data frame of initial cluster labels of each cell in the
#' sample.
#' @return a data frame of cell IDs of sorted neighbourhood cells.
#' @keywords internal

.cellNameSort <- function (cell, Clust) {
    cellDF <- data.frame(cellName = rownames(Clust)[cell],
                         cluster = Clust[cell, 1])
    if (length(unique(Clust[cell,])) == 1) { # cells from same initial cluster
        return(rownames(Clust)[cell])
    } else if (cellDF$cluster[1] %in% cellDF$cluster[seq(2, length(cell))]) {
        centralClust <- cellDF$cluster[1]
        sameClust <- c()
        diffClust <- c()
        for (c in seq_len(nrow(cellDF))) {
            if (cellDF$cluster[c] == centralClust){
                sameClust <- append(sameClust, cellDF$cellName[c])
            } else {
                diffClust <- append(diffClust, cellDF$cellName[c])}}
        region <- c(sameClust, diffClust)
        return(region)
    } else {
        return(rownames(Clust)[cell])}
}


#' Neighbour cell naming
#'
#' @description
#' A function to fetch cell IDs.
#' @param cell a vector of neighbourhood cell indices. The cell indices indicate
#' the row number of cells in sample metadata.
#' @param Clust a data frame of initial cluster labels of each cell in the
#' sample.
#' @return a data frame of cell IDs of neighbourhood cells.
#' @keywords internal

.cellName <- function (cell, Clust) {
    return(rownames(Clust)[cell])
}


#' Neighbour cell initial subcluster label
#'
#' @description
#' A function to fetch initial subcluster label of cell.
#' @param cell a vector of neighbourhood cell indices. The cell indices indicate
#' the row number of cells in sample metadata.
#' @param subClust a data frame of initial subcluster labels of each cell in
#' the sample.
#' @return a data frame of initial subcluster labels of neighbourhood
#' cells.
#' @keywords internal

.clustNum <- function (cell, subClust) {
    return(subClust[cell,])
}


#' Cell neighbourhood composition
#'
#' @description
#' A function to calculate the cell neighbourhood composition of initial cluster
#' labels.
#' @param arr a vector of initial subcluster labels of each cell in the
#' neighbourhood.
#' @return a table of initial subcluster proportions in a neighbourhood.
#' @keywords internal

.calculateProp <- function(arr) {
    prop <- prop.table(table(arr))
    # QC check
    if (round(sum(prop)) != 1) {
        stop("Proportions are incorrect.", "row =", c, "proportion =", sum(arr))
    }
    return(prop)
}


#' Gaussian distribution weights
#'
#' @description
#' A function to generate weights from a Gaussian distribution.
#' @param ed a numeric vector of entropy values of all cell neighbourhoods.
#' @param NN an integer for the number of neighbourhood cells including the
#' index cell.
#' @param sd a numeric value for standard deviation of Gaussian distribution.
#' @importFrom stats dnorm
#' @return a matrix where the columns contain weights associated with the
#' entropy values.
#' @keywords internal

.gauss_kernel <- function(ed, NN, sd) {
    weight_mat <- dnorm(sapply(ed, function(x) seq(0, x, length.out = NN)),
                        sd = sd)
    return (weight_mat)
}


#' Exponential distribution weights
#'
#' @description
#' A function to generate weights from an exponential distribution.
#' @param ed a numeric vector of entropy values of all cell neighbourhoods.
#' @param NN an integer for the number of neighbourhood cells including the
#' index cell.
#' @param rate a numeric value for rate of exponential distribution.
#' @importFrom stats dexp
#' @return a matrix where the columns contain weights associated with the
#' entropy values.
#' @keywords internal

.exp_kernel <- function(ed, NN, rate) {
    weight_mat <- dexp(sapply(ed, function(x) seq(0, x, length.out = NN)),
                       rate = rate)
    return (weight_mat)
}
