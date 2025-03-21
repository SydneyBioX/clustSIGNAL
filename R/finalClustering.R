#' Final non-spatial clustering
#'
#' @description
#' A function to perform clustering on adaptively smoothed gene expression data
#' to generate ClustSIGNAL clusters.
#' @param spe SpatialExperiment object containing the adaptively smoothed gene
#' expression.
#' @param batch a logical parameter for whether to perform batch correction.
#' Default value is FALSE.
#' @param batch_by a character indicating name of colData(spe) column containing
#' the batch names. Default value is 'None'.
#' @param threads a numeric value for the number of CPU cores to be used for the
#' analysis. Default value set to 1.
#' @param clustParams  a list of parameters for TwoStepParam clustering methods:
#' clust_c is the number of centers to use for clustering with KmeansParam. By
#' default set to 0, in which case the method uses either 5000 centers or 1/5th
#' of the total cells in the data as the number of centers, whichever is lower.
#' subclust_c is the number of centers to use for sub-clustering the initial
#' clusters with KmeansParam. This parameter is not used in the final
#' clustering step.
#' iter.max is the maximum number of iterations to perform during clustering
#' and sub-clustering with KmeansParam. Default value is 30.
#' k is a numeric value indicating the k-value used for clustering with
#' NNGraphParam. Default value is 10.
#' cluster.fun is a character indicating the graph clustering method used with
#' NNGraphParam. By default, the Louvain method is used.
#' @return SpatialExperiment object containing clusters generated from smoothed
#' data.
#' @importFrom bluster clusterRows TwoStepParam KmeansParam NNGraphParam
#' @importFrom scater runPCA
#' @importFrom harmony RunHarmony
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom methods show
#' @importFrom stats setNames
#' @examples
#' data(ClustSignal_example)
#'
#' # For non-spatial clustering of normalised counts
#' spe <- clustSIGNAL::p2_clustering(spe)
#' spe$ClustSIGNAL |> head()
#' @export

#### Non-spatial clustering
p2_clustering <- function(spe, batch = FALSE, batch_by = "None", threads = 1,
                          clustParams = list(clust_c = 0, subclust_c = 0,
                                             iter.max = 30, k = 10,
                                             cluster.fun = "louvain")) {
    if (clustParams[[1]] == 0)
        # number of centers = 1/5th of total cells in sample
        clustVal <- min(as.integer(ncol(spe) / 5), 5000) else
            clustVal <- clustParams[[1]]
    # Clustering adaptively smoothed data
    spe <- scater::runPCA(spe, assay.type = "smoothed", name = "PCA.smooth")
    if (batch == TRUE) {
        emb <- harmony::RunHarmony(data_mat = reducedDim(spe, "PCA.smooth"),
                                   meta_data = colData(spe),
                                   vars_use = batch_by, max.iter = 20,
                                   verbose = FALSE)
        mat <- emb
    } else {
        mat <- reducedDim(spe, "PCA.smooth")}
    reClust <- bluster::clusterRows(
        mat,
        bluster::TwoStepParam(
            first = bluster::KmeansParam(centers = clustVal,
                                         iter.max = clustParams[[3]]),
            second = bluster::NNGraphParam(k = clustParams[[4]],
                                           cluster.fun = clustParams[[5]],
                                           num.threads = threads)))
    spe$ClustSIGNAL <- factor(reClust)
    show(paste("Nonspatial clustering performed on smoothed data. Clusters =",
               length(unique(reClust)), "Time",
               format(Sys.time(),'%H:%M:%S')))
    return (spe)
}
