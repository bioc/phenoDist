\name{PDMBySvmWeightVector}
\alias{PDMBySvmWeightVector}
\title{ Compute phenotypic distance matrix by SVM weight vector }

\description{
  This function performs an SVM classification between a given sample and the negative control, calculates the weight vector,  and then computes the phenotypic distance matrix based on the weight vectors.
}

\usage{

PDMBySvmWeightVector(x, unames, neg='rluc', selectedCellFtrs, distMethod=c('manhattan','euclidean', 'correlation','mahalanobis'), verbose=FALSE, kernel='linear', ...)
}

\arguments{
  \item{x}{An \code{imageHTS} object.}
  \item{unames}{A character vector, containing the well names from where to collect the cell features. See \code{getUnames} for details.}
  \item{neg}{A character string to identify the negative controls.}
  \item{selectedCellFtrs}{A character vector for cell features to be used in the calculation. If missing, all features are used.}
  \item{distMethod}{A character string indicating which distance method should be used. This must be (an abbreviation of) one of the strings 'manhattan', 'euclidean', 'correlation' or 'mahalanobis'.}
  \item{verbose}{A logical scalar indicating whether progress should be reported.}
  \item{kernel}{The kernel argument for the \code{svm} function of the \code{e1071} package.}
  \item{...}{Additional arguments to be passed to the \code{svm} function of the \code{e1071} package.}
}

\details{
  For each well, this function collects features of all cells from the well and all cells from the negative control wells, and performs a bi-class classification using Support Vector Machine (SVM). The classification weight vectors are calculated for all wells passed to \code{PDMByWellAvg} to compute the phenotypic distance matrix.
}

\value{
  A symmetric distance matrix with dimensions equaling to the length of \code{unames}.
}

\seealso{\code{svm}, \code{PDMByWellAvg}}

\author{Xian Zhang}

\examples{

  library('phenoDist')

  ## load the imageHTS object
  load(system.file('kimorph', 'kimorph.rda', package='phenoDist'))
  x@localPath <- file.path(tempdir(), 'kimorph')

  ## calculate pair-wise svm distance matrix
  load(system.file('kimorph', 'selectedFtrs.rda', package='phenoDist'))
  pdm <- PDMBySvmWeightVector(x, unames=getUnames(x,plate=1, row=2:3, col=3), neg='rluc', selectedCellFtrs=selectedCellFtrs, distMethod='euclidean', verbose=FALSE, cost=1, kernel='linear')
}
