#' clustSIGNAL
#'
#' @description
#' A clustering method for cell type classification of spatial transcriptomics
#' data. The tool generates and uses an adaptively smoothed, spatially informed
#' gene expression data for clustering.
#'
#' @param spe a SpatialExperiment object.
#' @param samples a character indicating name of colData(spe) column containing
#' sample names.
#' @param cells a character indicating name of colData(spe) column containing
#' cell IDs.
#' @param dimRed a character indicating the name of the reduced dimensions to
#' use from the SpatialExperiment object (i.e., from reducedDimNames(spe)).
#' Default value is 'None'.
#' @param batch a logical parameter for whether or not to perform batch
#' correction. Default value is FALSE.
#' @param batch_by a character indicating name of colData(spe) column containing
#' the batch names.
#' @param NN an integer for the number of neighbourhood cells the function
#' should consider. The value must be greater than or equal to 1. Default value
#' is 30.
#' @param kernel a character for type of distribution to be used. The two valid
#' values are "G" or "E". G for Gaussian distribution, and E for exponential
#' distribution. Default value is "G".
#' @param spread a numeric value for distribution spread, represented by standard
#' deviation for Gaussian distribution and rate for exponential distribution.
#' Default value is 0.05 for Gaussian distribution and 20 for exponential
#' distribution.
#' @param sort a logical parameter for whether or not to sort the neighbourhood
#' after region description. Default value is TRUE.
#' @param threads a numeric value for the number of CPU cores to be used for the
#' analysis. Default value set to 1.
#' @param clustParams a list of parameters for TwoStepParam clustering methods.
#' The clustering parameters are in the order - centers (centers) for clustering
#' with KmeansParam, centers (centers) for sub-clustering clusters with
#' KmeansParam, maximum iterations (iter.max) for clustering with KmeansParam,
#' k values (k) for clustering with NNGraphParam, and community detection method
#' (cluster.fun) to use with NNGraphParam.
#' @param outputs a character for the type of output to return to the user. "c"
#' for data frame of cell IDs and their respective cluster numbers (default)
#' and "a" for list of dataframe of clusters plus final SpatialExperiment
#' object.
#'
#' @return a list of outputs
#'
#' 1. clusters: a data frame of cell names and their cluster classification.
#'
#' 2. spe_final: a SpatialExperiment object with initial 'putative cell type'
#' groups, entropy values, smoothed gene expression, post-smoothing clusters,
#' and silhouette widths included.
#'
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom methods show
#' @importFrom scater runPCA
#'
#' @examples
#' data(example)
#'
#' names(colData(spe))
#' # identify the column names with cell and sample labels
#' samples = "sample_id"
#' cells = "uniqueID"
#' res_list <- clustSIGNAL(spe, samples, cells, outputs = "c")
#'
#'
#' @export

clustSIGNAL <- function (spe, samples, cells, dimRed = "None", batch = FALSE,
                         batch_by = "None", NN = 30, kernel = "G",
                         spread = 0.05, sort = TRUE, threads = 1, outputs = "c",
                         clustParams = list(clust_c = 0, subclust_c = 0,
                                            iter.max = 30, k = 5,
                                            cluster.fun = "louvain")) {
    time_start <- Sys.time()
    # data and parameter checks
    if (length(spatialCoords(spe)) == 0){
        stop("Spatial coordinates not found.")
    } else if (is.null(logcounts(spe)) == TRUE) {
        stop("Normalised gene expression not found.")
    } else if (is.null(spe[[samples]]) == TRUE) {
        stop("samples name not found.")
    } else if (is.null(spe[[cells]]) == TRUE) {
        stop("cells name not found.")
    } else if (NN < 1){
        stop("NN cannot be less than 1.")
    } else if (!(kernel %in% c("G", "E"))) {
        stop("Invalid kernel type.")
    } else if (!(outputs %in% c('c', 'a'))) {
        stop("Invalid character for output type.")}

    if (dimRed == "None") {
        show(paste("Calculating PCA. Time", format(Sys.time(),'%H:%M:%S')))
        spe <- scater::runPCA(spe)
        dimRed <- "PCA"
    } else if (!(dimRed %in% reducedDimNames(spe))){
        stop("Specified low dimension data not found.")}
    show(paste("clustSIGNAL run started. Time", format(Sys.time(),'%H:%M:%S')))

    # Non-spatial clustering to identify initial cluster groups
    spe <- p1_clustering(spe, dimRed, batch, batch_by, clustParams)
    # Neighborhood detection, and sorting if sort = TRUE
    outReg <- neighbourDetect(spe, samples, NN, cells, sort)
    # Calculating domainness of cell neighborhoods
    spe <- entropyMeasure(spe, cells, outReg$regXclust, threads)
    # Weighted smoothing guided by neighbourhood entropy
    spe <- adaptiveSmoothing(spe, outReg$nnCells, NN, kernel, spread, threads)
    # Non-spatial clustering of adaptively smoothed expression
    spe <- p2_clustering(spe, batch, batch_by, clustParams)

    cluster_df <- data.frame("Cells" = spe[[cells]],
                             "Clusters" = spe$clustSIGNAL)
    time_end <- Sys.time()
    show(paste("clustSIGNAL run completed.", format(Sys.time(),'%H:%M:%S')))
    show(time_end - time_start)
    if (outputs == "c"){
        return (list("clusters" = cluster_df))
    } else if (outputs == "a") {
        return (list("clusters" = cluster_df,
                     "spe_final" = spe))
    }
}
