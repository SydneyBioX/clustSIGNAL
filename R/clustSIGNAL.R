#' ClustSIGNAL
#'
#' @description
#' A clustering method for spatially-resolved cell-type classification of
#' spatial transcriptomics data. The tool generates and uses an adaptively
#' smoothed, spatially informed gene expression data for clustering.
#' @param spe a SpatialExperiment object containing spatial coordinates in
#' 'spatialCoords' matrix and normalised gene expression in 'logcounts' assay.
#' @param samples a character indicating name of colData(spe) column containing
#' sample names.
#' @param dimRed a character indicating the name of the reduced dimensions to
#' use from the SpatialExperiment object (i.e., from reducedDimNames(spe)).
#' Default value is 'None'.
#' @param batch a logical parameter for whether to perform batch correction.
#' Default value is FALSE.
#' @param batch_by a character indicating name of colData(spe) column containing
#' the batch names. Default value is 'None'.
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
#' @param sort a logical parameter for whether to sort the neighbourhood by
#' initial clusters. Default value is TRUE.
#' @param threads a numeric value for the number of CPU cores to be used for the
#' analysis. Default value set to 1.
#' @param clustParams a list of parameters for TwoStepParam clustering methods:
#' clust_c is the number of centers to use for clustering with KmeansParam. By
#' default set to 0, in which case the method uses either 5000 centers or 1/5th
#' of the total cells in the data as the number of centers, whichever is lower.
#' subclust_c is the number of centers to use for sub-clustering the initial
#' clusters with KmeansParam. The default value is 0, in which case the method
#' uses either 1 center or half of the total cells in the initial cluster as
#' the number of centers, whichever is higher.
#' iter.max is the maximum number of iterations to perform during clustering
#' and sub-clustering with KmeansParam. Default value is 30.
#' k is a numeric value indicating the k-value used for clustering and
#' sub-clustering with NNGraphParam. Default value is 10.
#' cluster.fun is a character indicating the graph clustering method used with
#' NNGraphParam. By default, the Louvain method is used.
#' @param outputs a character for the type of output to return to the user. "c"
#' for data frame of cell IDs and their respective ClustSIGNAL cluster labels,
#' "n" for ClustSIGNAL cluster dataframe plus neighbourhood matrix, "s" for
#' ClustSIGNAL cluster dataframe plus final SpatialExperiment object, or "a" for
#' all 3 outputs.
#' @return a list of outputs depending on the type of outputs specified in the
#' main function call.
#'
#' 1. clusters: a data frame of cell names and their ClustSIGNAL cluster
#' classification.
#'
#' 2. neighbours: a character matrix containing cells IDs of each cell's
#' NN neighbours.
#'
#' 3. spe_final: a SpatialExperiment object containing the original spe object
#' data plus initial cluster and subcluster labels, entropy values, smoothed
#' gene expression, and ClustSIGNAL cluster labels.
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment reducedDimNames
#' @importFrom methods show
#' @importFrom scater runPCA
#' @examples
#' data(ClustSignal_example)
#'
#' names(colData(spe))
#' # identify the column name with sample labels
#' samples = "sample_id"
#' res_list <- clustSIGNAL(spe, samples, outputs = "c")
#' @export

clustSIGNAL <- function (spe, samples, dimRed = "None", batch = FALSE,
                         batch_by = "None", NN = 30, kernel = "G",
                         spread = 0.3, sort = TRUE, threads = 1, outputs = "c",
                         clustParams = list(clust_c = 0, subclust_c = 0,
                                            iter.max = 30, k = 10,
                                            cluster.fun = "louvain")) {
    time_start <- Sys.time()
    # input data and parameter checks
    if (length(spatialCoords(spe)) == 0) # check spatial coordinates
        stop("Spatial coordinates not found.")
    if (is.null(logcounts(spe)) == TRUE) # check normalised gene expression
        stop("Normalised gene expression not found.")
    if (is.null(spe[[samples]]) == TRUE) # check sample column
        stop("The column 'samples' not found.")
    if (length(unique(colnames(spe))) != ncol(spe)) # check cell IDs
        stop("The cell names are repeated. Provide spe object with unique cell
             names.")
    if (NN < 1) # check neighbourhood size
        stop("NN cannot be less than 1.")
    if (!(kernel %in% c("G", "E"))) # check kernel name
        stop("Invalid kernel type.")
    if (!(outputs %in% c('c', 'n', 's', 'a'))) # check output type
        stop("Invalid character for output type.")
    if (batch == TRUE & batch_by == "None") # check batch correction parameters
        stop("No group name provided for batch correction.")
    if (dimRed == "None") { # check reduced dimension data
        show(paste("Calculating PCA. Time", format(Sys.time(),'%H:%M:%S')))
        spe <- scater::runPCA(spe)
        dimRed <- "PCA"
    } else if (!(dimRed %in% reducedDimNames(spe))){
        stop("Specified low dimension data not found.")}
    show(paste("ClustSIGNAL run started. Time", format(Sys.time(),'%H:%M:%S')))

    spe <- p1_clustering(spe, dimRed, batch, batch_by, threads, clustParams)
    outReg <- neighbourDetect(spe, samples, NN, sort, threads)
    spe <- entropyMeasure(spe, outReg$regXclust, threads)
    spe <- adaptiveSmoothing(spe, outReg$nnCells, NN, kernel, spread)
    spe <- p2_clustering(spe, batch, batch_by, threads, clustParams)
    cluster_df <- data.frame("Cells" = colnames(spe),
                             "Clusters" = spe$ClustSIGNAL)
    time_end <- Sys.time()
    show(paste("ClustSIGNAL run completed.", format(Sys.time(),'%H:%M:%S')))
    show(time_end - time_start)
    if (outputs == "c"){
        return (list("clusters" = cluster_df))
    } else if (outputs == "n") {
        return (list("clusters" = cluster_df, "neighbours" = outReg$nnCells))
    } else if (outputs == "s") {
        return (list("clusters" = cluster_df, "spe_final" = spe))
    } else if (outputs == "a") {
        return (list("clusters" = cluster_df, "neighbours" = outReg$nnCells,
                     "spe_final" = spe))
    }
}
