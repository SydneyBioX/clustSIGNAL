## clustSIGNAL v0.99.8 (2025-03-21)
* The neighbourhood entropy calculation now includes the index cell's initial subcluster label. This is to ensure that distinct cells existing in otherwise completely homogeneous space are taken into account when measuring neighbourhood heterogeneity.
* Following additional parameter testing on simulated data, changed the default values of spread for Gaussian distribution (from 0.05 to 0.3) and clustering k (from 5 to 10), and recommended value of spread for exponential distribution (from 20 to 5).
* Moved all helper functions to the utilities.R file.
* Minor updates to the column names in which initial cluster and subcluster labels are stored after first clustering step.
* Updated the descriptions of all functions in the package. 
* Improved the tutorial vignette content.

## clustSIGNAL v0.99.7 (2025-01-08)
* Minor updates. 
* Renamed the example dataset from example.RData to ClustSignal_example.RData
* Updated DESCRIPTION file to add bluster package version (>= 1.16.0) required to run ClustSIGNAL.

## clustSIGNAL v0.99.6 (2024-12-11)
* Modified adaptiveSmoothing() function to improve its runtime.
* Improved the tutorial vignette content.
* Added step-by-step ClustSIGNAL run guide to the tutorial vignette.
* Removed the user parameter cells. The package automatically uses the column names of the input SpatialExperiment object as the cell ID.
* Added a cell id check - throws error if duplicates found among column names of input SpatialExperiment object.

## clustSIGNAL v0.99.5 (2024-11-28)
* Fixed issue where absence of spatial coordinates was not throwing error at the beginning of the run.

## clustSIGNAL v0.99.4 (2024-11-28)
* Corrected character check for output format type.

## clustSIGNAL v0.00.3 (2024-10-31)
* Clustering parameter options added to package
* Initial and final clustering steps are now in separate functions
* Test units added to the package

## clustSIGNAL v0.00.2 (2024-10-08)
* Updated parallel runs in smoothing function for faster runs

## clustSIGNAL v0.99.0 (2024-09-19)
* Submitted to Bioconductor.
