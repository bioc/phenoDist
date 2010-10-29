\name{PDMBySvmAccuracy}
\alias{PDMBySvmAccuracy}
\title{ Compute phenotypic distance matrix by SVM classification accuracy }

\description{
  This function performs an SVM classification between two given samples, and calculates the classification accuracy via cross validation as the phenotypic distance between the two samples. For multiple samples, the function returns a pair-wise distance matrix.
}

\usage{

PDMBySvmAccuracy(x, unames, selectedCellFtrs, cross=5, verbose=FALSE, ...)
}

\arguments{
  \item{x}{An \code{imageHTS} object.}
  \item{unames}{A character vector, containing the well names from where to collect the cell features. See \code{getUnames} for details.}
  \item{selectedCellFtrs}{A character vector for cell features to be used in the calculation. If missing, all features are used.}
  \item{cross}{An interger scalar indicating how many folds of cross validation should be performed.}
  \item{verbose}{A logical scalar indicating whether progress should be reported.}
  \item{...}{Additional arguments to be passed to the \code{svm} function of the \code{e1071} package.}
}

\details{
  For every pair of wells, this function collects features of all cells from both wells, and performs a bi-class classification using Support Vector Machine (SVM). The classification accuracy is defined as the phenotypic distance for the distance matrix.
}

\value{
  A symmetric distance matrix with the number of rows equaling to the length of \code{unames}.
}

\seealso{\code{svm}}

\author{Xian Zhang}

\examples{

  library('phenoDist')

  ## load the imageHTS object
  load(system.file('kimorph', 'kimorph.rda', package='phenoDist'))
  x@localPath <- file.path(tempdir(), 'kimorph')

  ## calculate pair-wise svm distance matrix
  load(system.file('kimorph', 'selectedFtrs.rda', package='phenoDist'))
  pdm <- PDMBySvmAccuracy(x, unames=getUnames(x, plate=1, row=2:3, col=3), selectedCellFtrs=selectedCellFtrs, cross=5, verbose=FALSE, cost=1, gamma=2^-5, kernel='radial')
}
