#' Final non-spatial clustering
#'
#' @description
#' A function to perform clustering on adaptively smoothed gene expression data
#' to generate ClustSIGNAL clusters.
#' @param spe SpatialExperiment object containing the adaptively smoothed gene
#' expression.
#' @param dimRed_f a character indicating the name of the reduced dimensions
#' in the SpatialExperiment object (i.e., from reducedDimNames(spe)) to use for
#' final clustering step. Two valid options are "None" (default), which triggers
#' a PCA run on smoothed expression, and "embed.smooth", which triggers a search
#' for "embed.smooth" low embedding in reducedDimNames(spe).
#' @param batch a logical parameter for whether to perform batch correction.
#' Default value is FALSE.
#' @param batch_by a character indicating name of colData(spe) column containing
#' the batch names. Default value is 'None'.
#' @param threads a numeric value for the number of CPU cores to be used for the
#' analysis. Default value set to 1.
#' @param clustParams  a list of parameters for TwoStepParam clustering methods:
#' clust_c is the number of centers to use for clustering with KmeansParam. By
#' default set to 0, in which case the method uses either 3000 centers or 1/5th
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
p2_clustering <- function(spe, dimRed_f = c("None", "embed.smooth"),
                          batch = FALSE, batch_by = "None", threads = 1,
                          clustParams = list(clust_c = 0, subclust_c = 0,
                                             iter.max = 30, k = 10,
                                             cluster.fun = "louvain")) {

    # check kernel name
    dimRed_f <- match.arg(dimRed_f)

    # extracting parameters
    clust_c <- clustParams$clust_c
    subclust_c <- clustParams$subclust_c
    iter.max <- clustParams$iter.max
    k <- clustParams$k
    cluster.fun <- clustParams$cluster.fun

    # checking availability of low dimension data
    if (dimRed_f == "None") {
        message(sprintf("%s Calculating PCA using smoothed data.",
                        format(Sys.time(), '%H:%M:%S')))
        # Using a fast approximate PCA (Irlba) instead of exact PCA
        spe <- runPCA(spe, assay.type = "smoothed", name = "embed.smooth",
                      BSPARAM = IrlbaParam())
        dimRed_f <- "PCA"
    } else if (!(dimRed_f %in% reducedDimNames(spe))){
        stop("'embed.smooth' low dimension data not found.")}

    if (batch == TRUE) {
        emb <- RunHarmony(data_mat = reducedDim(spe, "embed.smooth"),
                          meta_data = colData(spe),
                          vars_use = batch_by, max.iter = 20,
                          verbose = FALSE)
        mat <- emb
    } else {
        mat <- reducedDim(spe, "embed.smooth")}

    # Clustering adaptively smoothed data
    # Calculating number of cluster centres
    clustVal <-
        if (clust_c == 0) min(as.integer(ncol(spe) / 5), 3000) else clust_c
    reClust <- clusterRows(
        mat,
        TwoStepParam(
            first = KmeansParam(centers = clustVal, iter.max = iter.max),
            second = NNGraphParam(k = k, cluster.fun = cluster.fun,
                                  num.threads = threads)))
    spe$ClustSIGNAL <- factor(reClust)
    message(sprintf(
        "%s Final clustering performed on smoothed data. Clusters = %d ",
        format(Sys.time(),'%H:%M:%S'), length(unique(reClust))))
    return (spe)
}
