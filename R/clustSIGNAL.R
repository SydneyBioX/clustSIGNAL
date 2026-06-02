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
#' @param dimRed_init a character indicating the name of the reduced dimensions
#' in the SpatialExperiment object (i.e., from reducedDimNames(spe)) to use for
#' initial clustering step. Default value is 'None'.
#' @param dimRed_f a character indicating the name of the reduced dimensions
#' in the SpatialExperiment object (i.e., from reducedDimNames(spe)) to use for
#' final clustering step. Two valid options are "None" (default), which triggers
#' a PCA run on smoothed expression, and "embed.smooth", which triggers a search
#' for "embed.smooth" low embedding in reducedDimNames(spe).
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
#' all 3 outputs. Default value is 'c'.
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
#' @importFrom SingleCellExperiment reducedDimNames logcounts
#' @importFrom SummarizedExperiment colData
#' @importFrom methods show
#' @examples
#' data(ClustSignal_example)
#'
#' names(colData(spe))
#' # identify the column name with sample labels
#' samples = "sample_id"
#' res_list <- clustSIGNAL(spe, samples, outputs = "c")
#' @export

clustSIGNAL <- function(spe, samples, dimRed_init = "None",
                        dimRed_f = c("None", "embed.smooth"),
                        batch = FALSE, batch_by = "None", NN = 30,
                        kernel = c("G", "E"), spread = 0.3, sort = TRUE,
                        threads = 1, outputs = c("c", "n", "s", "a"),
                        clustParams = list(clust_c = 0, subclust_c = 0,
                                           iter.max = 30, k = 10,
                                           cluster.fun = "louvain")) {
    # method run begins
    time_start <- Sys.time()

    # check spatial coordinates
    if (length(spatialCoords(spe)) == 0)
        stop("No spatial coordinates found in spatialCoords(spe).")
    # check normalised gene expression
    if (is.null(logcounts(spe)))
        stop("Normalised gene expression (logcounts) not found.")
    # check sample column
    if (is.null(colData(spe)[[samples]]))
        stop(sprintf("Column '%s' not found in colData(spe)."), samples)
    # check cell IDs
    if (anyDuplicated(colnames(spe)) > 0)
        stop("Cell IDs are repeated. Provide spe object with unique cell
             IDs.")
    # check neighbourhood size
    if (NN < 1)
        stop("NN cannot be less than 1.")
    # check batch correction parameters
    if (isTRUE(batch) & batch_by == "None")
        stop("Batch correction enabled but no 'batch_by' group provided.")
    # check kernel name
    kernel <- match.arg(kernel)
    # check output type
    outputs <- match.arg(outputs)
    message(sprintf("%s ClustSIGNAL running.", format(Sys.time(), '%H:%M:%S')))

    spe <- p1_clustering(spe, dimRed_init, batch, batch_by, threads, clustParams)
    outReg <- neighbourDetect(spe, samples, NN, sort, threads)
    spe <- entropyMeasure(spe, outReg$regXclust)
    spe <- adaptiveSmoothing(spe, outReg$nnCells, NN, kernel, spread)
    spe <- p2_clustering(spe, dimRed_f, batch, batch_by, threads, clustParams)
    cluster_df <- data.frame("Cells" = colnames(spe),
                             "Clusters" = spe$ClustSIGNAL)
    # method run ends
    time_end <- Sys.time()
    message(sprintf("%s ClustSIGNAL completed.", format(Sys.time(), '%H:%M:%S')))
    show(time_end - time_start)

    # generating outputs
    res <- list("clusters" = cluster_df)
    if (outputs %in% c("n", "a")) {
        res[["neighbours"]] <- outReg$nnCells
    }
    if (outputs %in% c("s", "a")) {
        res[["spe_final"]] <- spe
    }

    return(res)
}
