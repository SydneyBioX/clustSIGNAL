% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neighborDetect.R
\name{neighbourDetect}
\alias{neighbourDetect}
\title{Cell neighbourhood detection}
\usage{
neighbourDetect(spe, samples, NN = 30, sort = TRUE, threads = 1)
}
\arguments{
\item{spe}{SpatialExperiment object with initial cluster and subcluster
labels.}

\item{samples}{a character indicating name of colData(spe) column containing
sample names.}

\item{NN}{an integer for the number of neighbouring cells the function
should consider. The value must be greater than or equal to 1. Default value
is 30.}

\item{sort}{a logical parameter for whether to sort the neighbourhood by
initial clusters. Default value is TRUE.}

\item{threads}{a numeric value for the number of CPU cores to be used for the
analysis. Default value set to 1.}
}
\value{
a list containing two items:

1. nnCells, a character matrix of NN nearest neighbours - rows are index
cells and columns are their nearest neighbours ranging from closest to
farthest neighbour. For sort = TRUE, the neighbours belonging to the same
initial cluster as the index cell are moved closer to it.

2. regXclust, a list of vectors of each cell's neighbourhood composition
indicated by the proportion of initial subclusters it contains.
}
\description{
A function to identify the neighbourhood of each cell. If sort = TRUE, the
neighbourhoods are also sorted such that cells belonging to the same
'initial cluster' as the index cell are arranged closer to it.
}
\examples{
data(ClustSignal_example)

out_list <- clustSIGNAL::neighbourDetect(spe, samples = "sample_id")
out_list |> names()
}
