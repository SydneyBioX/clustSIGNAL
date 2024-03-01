#' Additional internal functions
#'
#' @description
#' .cellName: A function to fetch cell IDs. Returns a data frame of cell IDs of neighbourhood cells.
#' .cellNameSort: A function to perform neighbourhood cell sorting. Neighbourhood cells that belong to the same 'putative cell type' as the central cell are moved closer to the central cell. Returns a data frame of cell IDs of sorted neighbourhood cells.
#' .clustNum: A function to fetch 'putative cell type' assignments of cell. Returns a data frame of 'putative cell type' assignments of neighbourhood cells.
#' .calculateProp: A function to calculate the cell neighbourhood composition of 'putative' cell types. Returns a table of 'putative cell type' proportions in a neighbourhood.
#' .exp_kernel: A function to generate weights from an exponential distribution. Returns a numeric vector of weights.
#' .gauss_kernel: A function to generate weights from a Gaussian distribution. Returns a numeric vector of weights.
#' .smoothedData: A function to perform a weighted moving average of gene expression in a neighbourhood. Returns a column matrix of smoothed gene expression.
#' .getColours: Returns a vector of colours used in clustSIGNAL spatial plots by default.
#'
#' @keywords internal


# function to retrieve cell names of neighbors
.cellName <- function (cell, Clust) {
    return(rownames(Clust)[cell])
}


# function to sort neighbours and retrieving cell names of neighbors
.cellNameSort <- function (cell, Clust) {
    cellDF = data.frame(cellName = rownames(Clust)[cell], cluster = Clust[cell,1])
    if (cellDF$cluster[1] %in% cellDF$cluster[2:31]) {
        centralClust = cellDF$cluster[1]
        sameClust = c()
        diffClust = c()
        for (c in 1:nrow(cellDF)) {
            if (cellDF$cluster[c] == centralClust){
                sameClust = append(sameClust, cellDF$cellName[c])
            } else {
                diffClust = append(diffClust, cellDF$cellName[c])
            }
        }
        region = c(sameClust, diffClust)
        return(region)
    } else {
        return(rownames(Clust)[cell])
    }
}


# function to retrieve cluster numbers to which neighbors belong
.clustNum <- function (cell, subClust) {
    return(subClust[cell,])
}


# function to retrieve cluster proportions of each region
.calculateProp <- function(arr) {
    prop = prop.table(table(arr))
    if (round(sum(prop)) != 1){ # check if the proportions for each region sum up to 1
        stop(paste("ERROR: proportions are incorrect.", "row =", c, "proportion =", sum(arr)))
    }
    return(prop)
}


# function to retrieve exponential distribution weights for neighbors
.exp_kernel <- function(ed, NN, rate) {
    # distribution from 0 to entropy, with cells in smoothing radius as cut points
    cutpoints <- seq(0, ed, length.out = NN)
    return(dexp(cutpoints, rate = rate))
}


# function to retrieve normal distribution weights for neighbors
.gauss_kernel <- function(ed, NN, sd) {
    # distribution from 0 to entropy, with cells in smoothing radius as cut points
    cutpoints <- seq(0, ed, length.out = NN)
    return(dnorm(cutpoints, sd = sd))
}


# function to perform weighted moving average smoothing
.smoothedData <- function(mat, weight) {
    # column matrix of weight
    weight <- as.matrix(weight)
    # return weighted moving average
    return((mat %*% weight) / sum(weight))
}


# function to select colours for spatial plots
.getColours = function(){
    return (c("#635547", "#8EC792", "#9e6762", "#FACB12", "#3F84AA", "#0F4A9C", "#ff891c", "#EF5A9D", "#C594BF", "#DFCDE4",
              "#139992", "#65A83E", "#8DB5CE", "#005579", "#C9EBFB", "#B51D8D", "#532C8A", "#8870ad", "#cc7818", "#FBBE92",
              "#EF4E22", "#f9decf", "#c9a997", "#C72228", "#f79083", "#F397C0", "#DABE99", "#c19f70", "#354E23", "#C3C388",
              "#647a4f", "#CDE088", "#f7f79e", "#F6BFCB", "#7F6874", "#989898", "#1A1A1A", "#FFFFFF", "#e6e6e6", "#77441B",
              "#F90026", "#A10037", "#DA5921", "#E1C239", "#9DD84A"))
}
