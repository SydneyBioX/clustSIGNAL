#' clustSIGNAL
#'
#' @description
#' A clustering method for cell type classification of spatial transcriptomics data. The tool generates and uses an adaptively smoothed, spatially informed gene expression data for clustering.
#'
#' @param spe a SpatialExperiment object.
#' @param samples a character indicating name of colData(spe) column containing sample names.
#' @param cells a character indicating name of colData(spe) column containing cell IDs.
#' @param dimRed a character indicating the name of the reduced dimensions to use from the SpatialExperiment object (i.e., from reducedDimNames(spe)). Default value is 'PCA'.
#' @param NN an integer for the number of neighbourhood cells the function should consider. The value must be greater than or equal to 1. Default value is 30.
#' @param kernel a character for type of distribution to be used. The two valid values are "G" or "E". G for Gaussian distribution, and E for exponential distribution. Default value is "G".
#' @param spread a numeric value for distribution spread, represented by standard deviation for Gaussian distribution and rate for exponential distribution. Default value is 0.05 for Gaussian distribution and 20 for exponential distribution.
#' @param sort a logical parameter for whether or not to sort the neighbourhood after region description. Default value is TRUE.
#' @param threads a numeric value for the number of CPU cores to be used for the analysis. Default value set to 4 cores.
#' @param outputs a character for the type of output to return to the user. "c" for data frame of cell IDs and their respective cluster numbers (default), "n" for list of dataframe of clusters plus neighbourhood matrix, "s" for list of dataframe of clusters plus final spatialExperiment object, or "a" for list of all outputs.
#' @param ... additional parameters for TwoStepParam clustering methods. Include parameters like k for number of nearest neighbours and cluster.fun for selecting community detection method. Default values k = 5, cluster.fun = "louvain".
#'
#' @return a list of outputs
#'
#' 1. clusters: a data frame of cell names and their cluster classification.
#'
#' 2. neighbours: a matrix of cell names and the names of their NN nearest neighbour cells.
#'
#' 3. spe_final: a SpatialExperiment object with initial 'putative cell type' groups, entropy values, smoothed gene expression, post-smoothing clusters, and silhouette widths included.
#'
#' @examples
#' data(mouseEmbryo2)
#'
#' names(colData(spe))
#' # identify the column names with cell and sample labels
#' samples = "embryo"
#' cells = "uniqueID"
#' res_list <- clustSIGNAL(spe, samples, cells, dimRed = "PCA", outputs = "c")
#'
#'
#' @export

clustSIGNAL <- function (spe,
                         samples,
                         cells,
                         dimRed = "None",
                         NN = 30,
                         kernel = "G",
                         spread = 0.05,
                         sort = TRUE,
                         threads = 1,
                         outputs = "c",
                         ...) {

    # data and parameter checks
    if (is.null(spatialCoords(spe)) == TRUE){
        stop("The spe object should contain spatial coordinates as
             spatialCoords matrix in the object.")
    } else if (is.null(spe[[samples]]) == TRUE) {
        stop("The samples variable does not match any colData column name in the spe object.")
    } else if (is.null(spe[[cells]]) == TRUE) {
        stop("The cells variable does not match any colData column name in the spe object.")
    } else if (NN < 1){
        stop("Number of nearest neighbours cannot be less than 1.")
    } else if (!(kernel %in% c("G", "E"))) {
        stop("Invalid character for kernel - use 'G' or 'E' to indicate the
             weight distribution you would like to use.")
    } else if (!(outputs %in% c('c', 'n', 's', 'a'))) {
        stop("Invalid character for output type - use 'c', 'n', 's', or 'a' to
              indicate the type of output you would like.")
    }

    if (dimRed == "None") {
        print("Calculating PCA to use as reduced dimension input.")
        spe <- scater::runPCA(spe)
        dimRed = "PCA"
    } else if (!(dimRed %in% reducedDimNames(spe))){
        stop("The specified reduced dimension data does not exist.")
    }

    print(paste("clustSIGNAL run started.", Sys.time()))

    # Non-spatial clustering to identify 'putative celltype' groups
    # reclust should always be FALSE here
    spe <- nsClustering(spe, dimRed, reclust = FALSE, ...)

    # Neighborhood detection, and sorting if sort = "yes"
    outReg <- neighbourDetect(spe, samples, NN, cells, sort)

    # Calculating domainness of cell neighborhoods
    # low entropy indicates more homogeneous neighborhood
    # high entropy indicates more heterogeneous neighborhood
    spe <- entropyMeasure(spe, cells, outReg$regXclust, threads)

    # Weighted smoothing guided by neighbourhood entropy
    # Homogeneous regions are smoothed more
    # Heterogeneous regions are smoothed less
    spe <- adaptiveSmoothing(spe, outReg$nnCells, NN, kernel, spread)

    # Non-spatial clustering of adaptively smoothed expression
    # reclust should always be TRUE here
    spe <- nsClustering(spe, reclust = TRUE, ...)

    cluster_df <- data.frame("Cells" = spe[[cells]], "Clusters" = spe$reCluster)

    print(paste("clustSIGNAL run completed.", Sys.time()))

    if (outputs == "c"){
        return (list("clusters" = cluster_df))
    } else if (outputs == "n"){
        return (list("clusters" = cluster_df,
                     "neighbours" = outReg$nnCells))
    } else if (outputs == "s") {
        return (list("clusters" = cluster_df,
                     "spe_final" = spe))
    } else if (outputs == "a") {
        return (list("clusters" = cluster_df,
                     "neighbours" = outReg$nnCells,
                     "spe_final" = spe))
    }
}
