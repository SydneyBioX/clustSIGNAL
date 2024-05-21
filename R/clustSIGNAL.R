#' clustSIGNAL
#'
#' @description
#' A clustering method for cell type classification of spatial transcriptomics data. The tool generates and uses an adaptively smoothed, spatially informed gene expression data for clustering.
#'
#' @param spe SpatialExperiment object with logcounts and PCA included.
#' @param samples a character vector of sample names to which the cells belong. Length of vector must be equal to the number of cells in SpatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param cells a character vector of cell IDs of each cell. Length of vector must be equal to the number of cells in SpatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param NN an integer for the number of neighbourhood cells the function should consider. The value must be greater than or equal to 1. Default value is 30.
#' @param kernel a character for type of distribution to be used. The two valid values are "G" or "E". G for Gaussian distribution, and E for exponential distribution. Default value is "G".
#' @param spread a numeric value for distribution spread, represented by standard deviation for Gaussian distribution and rate for exponential distribution. Default value is 0.05 for Gaussian distribution and 20 for exponential distribution.
#' @param sort a logical parameter for whether or not to sort the neighbourhood after region description. Default value is TRUE.
#' @param threads a numeric value for the number of CPU cores to be used for the analysis. Default value set to 4 cores.
#' @param outputs a character for the type of output to return to the user. "c" for data frame of cell IDs and their respective cluster numbers (default), "n" for dataframe of clusters plus neighbourhood matrix, "s" for  dataframe of clusters plus final spatialExperiment object, or "a" for all outputs.
#'
#' @return a list of outputs
#'
#' 1. domainSpread: a histogram of entropy values and number of cell neighbourhoods.
#'
#' 2. domainDistribution: a spatial plot showing entropy of all cell neighbourhoods.
#'
#' 3. spatialAnnotated: a spatial plot of annotated cell types. If no celltype annotations were provides, the output will be NULL.
#'
#' 4. spatialCluster: a spatial plot of clusters generated using clustSIGNAL.
#'
#' 5. compareHeatmap: a heat map showing cell number distributions of annotated cell type and output clusters. If no celltype annotations were provides, the output will be NULL.
#'
#' 6. metrics: a table with ARI, NMI, mean silhouette width, minimum entropy, maximum entropy, and mean entropy of each sample in the dataset. If no celltype annotations were provides, the ARI and NMI values will be NULL.
#'
#' 7. neighbours: a matrix of cell names and the names of their NN nearest neighbour cells.
#'
#' 8. spe_final: a SpatialExperiment object with initial 'putative cell type' groups, entropy values, smoothed gene expression, post-smoothing clusters, and silhouette widths included.
#'
#' @examples
#'
#' # to do
#'
#' @export

clustSIGNAL <- function (spe,
                         samples,
                         cells,
                         NN = 30,
                         kernel = "G",
                         spread = 0.05,
                         sort = TRUE,
                         threads = 4,
                         outputs = "c") {

    # require(BiocNeighbors)
    # require(bluster)
    # require(DescTools)
    # require(dplyr)
    # require(foreach)
    # require(Hmisc)
    # require(parallel)
    # require(PCAtools)
    # require(scater)
    # require(scran)
    # require(scuttle)
    # require(SpatialExperiment)

    if (NN < 1){
        stop("ERROR: Number of nearest neighbours cannot be less than 1.")
    }

    print(paste("Adaptive clustering run started.", Sys.time()))

    # Non-spatial clustering
    # reclust should always be FALSE here
    spe <- nsClustering(spe, reclust = FALSE)

    # Region description and sorting, if sort = "yes"
    outReg <- neighbourDetect(spe, samples, NN, cells, sort)

    # Domainness measure and spread
    spe <- entropyMeasure(spe, cells, outReg$regXclust, threads)

    # Smoothing
    spe <- adaptiveSmoothing(spe, outReg$nnCells, NN, kernel, spread)

    # Non-spatial reclustering
    # reclust should always be TRUE here
    spe <- nsClustering(spe, reclust = TRUE)

    cluster_df <- as.data.frame(cells) %>%
        cbind(as.data.frame(spe$reCluster))

    print(paste("Adaptive clustering run completed. Outputs are ready.", Sys.time()))

    if (outputs == "c"){
        return (cluster_df)
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
    } else {
        print("ERROR: invalid character used. Use 'c', 'n', 's', or 'a' to indicate the type of output you would like." )
    }
}
