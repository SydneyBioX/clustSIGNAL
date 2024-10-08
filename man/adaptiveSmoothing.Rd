% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptiveSmoothing.R
\name{adaptiveSmoothing}
\alias{adaptiveSmoothing}
\title{Adaptive smoothing}
\usage{
adaptiveSmoothing(spe, nnCells, NN, kernel, spread, threads)
}
\arguments{
\item{spe}{SpatialExperiment object with logcounts, PCA, 'putative cell type'
groups, and entropy outputs included.}

\item{nnCells}{a character matrix of NN nearest neighbours - rows are cells
and columns are their nearest neighbours ranged from closest to farthest
neighbour. For sort = TRUE, the neighbours belonging to the same 'putative
cell type' group as the cell are moved closer to it.}

\item{NN}{an integer for the number of neighbourhood cells the function should
consider. The value must be greater than or equal to 1. Default value is 30.}

\item{kernel}{a character for type of distribution to be used. The two valid
values are "G" or "E". G for Gaussian distribution, and E for exponential
distribution. Default value is "G".}

\item{spread}{a numeric value for distribution spread, represented by standard
deviation for Gaussian distribution and rate for exponential distribution.
Default value is 0.05 for Gaussian distribution and 20 for exponential
distribution.}

\item{threads}{a numeric value for the number of CPU cores to be used for the
analysis. Default value set to 1.}
}
\value{
SpatialExperiment object including smoothed gene expression values as
another assay.
}
\description{
A function to perform a weighted, adaptive smoothing of gene expression based
on the heterogeneity of the cell neighbourhood. Heterogeneous neighbourhoods
are smoothed less with higher weights given to cells belonging to same initial
group. Homogeneous neighbourhoods are smoothed more with similar weights given
to most cells.
}
\examples{
data(example)

# requires matrix containing NN nearest neighbour cell labels (nnCells),
# generated using the neighbourDetect() function
spe <- adaptiveSmoothing(spe, nnCells, NN = 30, kernel = "G", spread = 0.05,
                         threads = 1)
spe

}
