#' Initial non-spatial clustering
#'
#' @description
#' A function to perform initial non-spatial clustering and sub-clustering of
#' normalised gene expression to generate 'initial clusters' and 'initial
#' subclusters'.
#' @param spe a SpatialExperiment object containing spatial coordinates in
#' 'spatialCoords' matrix and normalised gene expression in 'logcounts' assay.
#' @param dimRed_init a character indicating the name of the reduced dimensions
#' in the SpatialExperiment object (i.e., from reducedDimNames(spe)) to use for
#' initial clustering step. Default value is 'None'.
#' @param batch a logical parameter for whether to perform batch correction.
#' Default value is FALSE.
#' @param batch_by a character indicating name of colData(spe) column containing
#' the batch names. Default value is 'None'.
#' @param threads a numeric value for the number of CPU cores to be used for the
#' analysis. Default value set to 1.
#' @param clustParams a list of parameters for TwoStepParam clustering methods:
#' clust_c is the number of centers to use for clustering with KmeansParam. By
#' default set to 0, in which case the method uses either 3000 centers or 1/5th
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
#' @return SpatialExperiment object with initial cluster and subcluster labels
#' of each cell.
#' @importFrom bluster clusterRows TwoStepParam KmeansParam NNGraphParam
#' @importFrom scater runPCA
#' @importFrom harmony RunHarmony
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom stats setNames
#' @importFrom BiocParallel bplapply
#' @importFrom BiocSingular IrlbaParam
#' @examples
#' data(ClustSignal_example)
#'
#' spe <- clustSIGNAL::p1_clustering(spe, dimRed_init = "PCA")
#' spe$nsCluster |> head()
#' spe$initCluster |> head()
#' @export

p1_clustering <- function(spe, dimRed_init = "None", batch = FALSE,
                          batch_by = "None", threads = 1,
                          clustParams = list(clust_c = 0, subclust_c = 0,
                                             iter.max = 30, k = 10,
                                             cluster.fun = "louvain")) {

    # extracting parameters
    clust_c <- clustParams$clust_c
    subclust_c <- clustParams$subclust_c
    iter.max <- clustParams$iter.max
    k <- clustParams$k
    cluster.fun <- clustParams$cluster.fun

    # checking availability of low dimension data
    if (dimRed_init == "None") {
        message(sprintf("%s Calculating PCA.", format(Sys.time(), '%H:%M:%S')))
        spe <- runPCA(spe)
        dimRed_init <- "PCA"
    } else if (!(dimRed_init %in% reducedDimNames(spe))){
        stop("Specified low dimension data not found.")}

    # calculating number of cluster centres
    clustVal <-
        if (clust_c == 0) min(as.integer(ncol(spe) / 5), 3000) else clust_c

    # checking if batch correction needed
    if (batch == TRUE) {
        mat <- RunHarmony(data_mat = reducedDim(spe, dimRed_init),
                          meta_data = colData(spe), vars_use = batch_by,
                          max.iter = 20, verbose = FALSE)
    } else {
        mat <- reducedDim(spe, dimRed_init)
    }

    # initial clustering
    initClust <- clusterRows(
        mat,
        TwoStepParam(
            first = KmeansParam(centers = clustVal, iter.max = iter.max),
            second = NNGraphParam(k = k, cluster.fun = cluster.fun,
                                  num.threads = threads)))
    spe$initCluster <- factor(initClust)
    message(sprintf("%s Initial clustering performed. Clusters = %d",
                    format(Sys.time(), '%H:%M:%S'), length(unique(initClust))))

    # initial sub-clustering with parallel runs
    BPPARAM <- .generateBPParam(cores = threads)
    clusters <- unique(spe$initCluster)
    subclusters_list <- bplapply(clusters, function(cl) {
        speX <- spe[, spe$initCluster == cl]
        # Using a fast approximate PCA (Irlba) instead of exact PCA
        speX <- runPCA(speX, BSPARAM = IrlbaParam())
        # calculating number of sub-cluster centres
        subclustVal <-
            if (subclust_c == 0) max(as.integer(ncol(speX) / 2), 1) else subclust_c
        initClustX <- clusterRows(
            reducedDim(speX, "PCA"),
            TwoStepParam(
                first = KmeansParam(centers = subclustVal, iter.max = iter.max),
                second = NNGraphParam(k = k, cluster.fun = cluster.fun)))
        return(setNames(paste0(cl, ".", initClustX), colnames(speX)))
    }, BPPARAM = BPPARAM)
    spe$initSubcluster <- factor(unlist(subclusters_list)[colnames(spe)])
    message(
        sprintf(
            "%s Initial sub-clustering performed. Subclusters = %d ",
            format(Sys.time(), '%H:%M:%S'), length(unique(spe$initSubcluster))))
    return (spe)
}
