#' ClustSIGNAL
#'
#' @description
#' A clustering method for cell type classification of spatial transcriptomics
#' data. The tool generates and uses an adaptively smoothed, spatially informed
#' gene expression data for clustering.
#'
#' @param spe a SpatialExperiment object.
#' @param samples a character indicating name of colData(spe) column containing
#' sample names.
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
#' @param spread a numeric value for distribution spread, represented by
#' standard deviation for Gaussian distribution and rate for exponential
#' distribution. Default value is 0.05 for Gaussian distribution and 20 for
#' exponential distribution.
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
#' for data frame of cell IDs and their respective cluster numbers (default),
#' "n" for dataframe of clusters plus neighbourhood matrix, "s" for  dataframe
#' of clusters plus final spatialExperiment object, or "a" for all outputs.
#'
#' @return a list of outputs
#'
#' 1. clusters: a data frame of cell names and their cluster classification.
#'
#' 2. neighbours: a character matrix containing cells IDs of each cell's
#' neighbours
#'
#' 3. spe_final: a SpatialExperiment object with initial 'putative cell type'
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
#' res_list <- clustSIGNAL(spe, samples, outputs = "c")
#'
#'
#' @export

clustSIGNAL <- function (spe, samples, dimRed = "None", batch = FALSE,
                         batch_by = "None", NN = 30, kernel = "G",
                         spread = 0.05, sort = TRUE, threads = 1, outputs = "c",
                         clustParams = list(clust_c = 0, subclust_c = 0,
                                            iter.max = 30, k = 5,
                                            cluster.fun = "louvain")) {
    time_start <- Sys.time()
    # input data and parameter checks
    if (length(spatialCoords(spe)) == 0)
        stop("Spatial coordinates not found.")
    if (is.null(logcounts(spe)) == TRUE)
        stop("Normalised gene expression not found.")
    if (is.null(spe[[samples]]) == TRUE)
        stop("The column 'samples' not found.")
    if (length(unique(colnames(spe))) != ncol(spe))
        stop("The cell names are repeated. Provide spe object with unique cell
             names.")
    if (NN < 1)
        stop("NN cannot be less than 1.")
    if (!(kernel %in% c("G", "E")))
        stop("Invalid kernel type.")
    if (!(outputs %in% c('c', 'n', 's', 'a')))
        stop("Invalid character for output type.")

    if (batch == TRUE & batch_by == "None")
        stop("No group name provided for batch correction.")
    if (dimRed == "None") {
        show(paste("Calculating PCA. Time", format(Sys.time(),'%H:%M:%S')))
        spe <- scater::runPCA(spe)
        dimRed <- "PCA"
    } else if (!(dimRed %in% reducedDimNames(spe))){
        stop("Specified low dimension data not found.")}
    show(paste("ClustSIGNAL run started. Time", format(Sys.time(),'%H:%M:%S')))

    # Non-spatial clustering to identify initial cluster groups
    spe <- p1_clustering(spe, dimRed, batch, batch_by, threads, clustParams)
    # Neighborhood detection, and sorting if sort = TRUE
    outReg <- neighbourDetect(spe, samples, NN, sort, threads)
    # Calculating domainness of cell neighborhoods
    spe <- entropyMeasure(spe, outReg$regXclust, threads)
    # Weighted smoothing guided by neighbourhood entropy
    spe <- adaptiveSmoothing(spe, outReg$nnCells, NN, kernel, spread)
    # Non-spatial clustering of adaptively smoothed expression
    spe <- p2_clustering(spe, batch, batch_by, threads, clustParams)
    cluster_df <- data.frame("Cells" = colnames(spe),
                             "Clusters" = spe$ClustSIGNAL)
    time_end <- Sys.time()
    show(paste("ClustSIGNAL run completed.", format(Sys.time(),'%H:%M:%S')))
    show(time_end - time_start)
    if (outputs == "c"){
        return (list("clusters" = cluster_df))
    } else if (outputs == "n") {
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
