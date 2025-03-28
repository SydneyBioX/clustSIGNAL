% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entropyMeasure.R
\name{entropyMeasure}
\alias{entropyMeasure}
\title{Heterogeneity measure}
\usage{
entropyMeasure(spe, regXclust, threads = 1)
}
\arguments{
\item{spe}{SpatialExperiment object with initial cluster and subcluster
labels.}

\item{regXclust}{a list of vectors of each cell's neighbourhood composition
indicated by the proportion of initial subclusters it contains.}

\item{threads}{a numeric value for the number of CPU cores to be used for the
analysis. Default value set to 1.}
}
\value{
SpatialExperiment object with entropy values associated with each
cell.
}
\description{
A function to measure the heterogeneity of a cell's neighbourhood in terms of
entropy. Generally, homogeneous neighbourhoods have low entropy and
heterogeneous neighbourhoods have high entropy.
}
\examples{
data(ClustSignal_example)

# requires list containing cluster proportions of each region (regXclust),
# generated using the neighbourDetect() function
spe <- clustSIGNAL::entropyMeasure(spe, regXclust)
spe$entropy |> head()
}
