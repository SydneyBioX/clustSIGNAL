#' clustSIGNAL
#'
#' @description
#' A clustering method for cell type classification of spatial transcriptomics data. The tool generates and uses an adaptively smoothed, spatially informed gene expression data for clustering.
#'
#' @param spe SpatialExperiment object with logcounts and PCA included.
#' @param samples a character vector of sample names to which the cells belong. Length of vector must be equal to the number of cells in SpatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param cells a character vector of cell IDs of each cell. Length of vector must be equal to the number of cells in SpatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param celltypes a character vector of cell type annotation of each cell. Length of vector must be equal to the number of cells in SpatialExperiment object (i.e. the number of rows in colData(spe)).
#' @param NN an integer for the number of neighbourhood cells the function should consider. The value must be greater than or equal to 1. Default value is 30.
#' @param kernel a character for type of distribution to be used. The two valid values are "G" or "E". G for Gaussian distribution, and E for exponential distribution. Default value is "G".
#' @param spread a numeric value for distribution spread, represented by standard deviation for Gaussian distribution and rate for exponential distribution. Default value is 0.05 for Gaussian distribution and 20 for exponential distribution.
#' @param sort a logical parameter for whether or not to sort the neighbourhood after region description. Default value is TRUE.
#' @param threads a numeric value for the number of CPU cores to be used for the analysis.
#'
#' @return a list of outputs
#'
#' 1. domainSpread, a histogram of entropy values and number of cell neighbourhoods
#'
#' 2. domainDistribution, a spatial plot showing entropy of all cell neighbourhoods
#'
#' 3. spatialAnnotated, a spatial plot of annotated cell types
#'
#' 4. spatialCluster, a spatial plot of clusters generated using clustSIGNAL
#'
#' 5. compareHeatmap, a heat map showing cell number distributions of annotated cell type and output clusters
#'
#' 6. metrics, a table with ARI, NMI, mean silhouette width, minimum entropy, maximum entropy, and mean entropy of each sample in the dataset
#'
#' 7. spe_final, a SpatialExperiment object with initial 'putative cell type' groups, entropy values, smoothed gene expression, post-smoothing clusters, and silhouette widths included.
#'
#' @examples
#'
#' @export

clustSIGNAL <- function (spe,
                        samples,
                        cells,
                        celltypes,
                        NN = 30,
                        kernel = "G",
                        spread = 0.05,
                        sort = TRUE,
                        threads = 20) {

    require(BiocNeighbors)
    require(bluster)
    require(circlize)
    require(cluster)
    require(ComplexHeatmap)
    require(distances)
    require(DescTools)
    require(foreach)
    require(ggplot2)
    require(ggspavis)
    require(Hmisc)
    require(parallel)
    require(PCAtools)
    require(plotly)
    require(scater)
    require(scran)
    require(scuttle)
    require(SpatialExperiment)
    require(RColorBrewer)

    print(paste("Adaptive clustering run started.", Sys.time()))

    # Non-spatial clustering
    # reclust should always be FALSE here
    spe <- nsClustering(spe, reclust = FALSE)

    # Region description and sorting, if sort = "yes"
    outReg <- neighbourDetect(spe, samples, NN, cells, sort)

    # Domainness measure and spread
    spe <- entropyMeasure(spe, cells, outReg$regXclust, threads)
    # entropy histogram and spatial plots
    domPlotList <- plots(spe, samples, celltypes, plotType = "domain")

    # Smoothing
    spe <- adaptiveSmoothing(spe, outReg$nnCells, NN, kernel, spread, cells, threads)

    # Non-spatial reclustering
    # reclust should always be TRUE here
    spe <- nsClustering(spe, reclust = TRUE)
    # spatial plots of annotated and cluster labels
    spatialPlot <- plots(spe, samples, celltypes, plotType = "spatial")

    # Cluster metrics
    spe <- metrics(spe, samples, celltypes, cells)$spe_out
    metrics_per_sample <- metrics(spe, samples, celltypes, cells)$metrics_out
    # heat map comparison between annotated and cluster labels
    cellHeat <- plots(spe, samples, celltypes, plotType = "heatmap")

    print(paste("Adaptive clustering run completed. Outputs are ready.", Sys.time()))

    return (list("domainSpread" = domPlotList$spread,
                 "domainDistribution" = domPlotList$distribution,
                 "spatialAnnotated" = spatialPlot$cellPlot,
                 "spatialCluster" = spatialPlot$clusterPlot,
                 "compareHeatmap" = cellHeat,
                 "metrics" = metrics_per_sample,
                 "spe_final" = spe))
}
