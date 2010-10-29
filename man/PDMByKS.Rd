\name{PDMByKS}
\alias{PDMByKS}
\title{ Compute phenotypic distance matrix by Kolmogorov-Smirnov statistics }

\description{
  This function calculates the Kolmogorov-Smirnov statistic between a given sample and the negative control, for each feature, and then computes the phenotypic distance matrix based on the Kolmogorov-Smirnov statistics.
}

\usage{

PDMByKS(x, unames, neg='rluc', selectedCellFtrs, distMethod=c('manhattan', 'euclidean', 'correlation', 'mahalanobis'))
}

\arguments{
  \item{x}{An \code{imageHTS} object.}
  \item{unames}{A character vector, containing the well names from where to collect the cell features. See \code{getUnames} for details.}
  \item{neg}{A character string to identify negative controls.}
  \item{selectedCellFtrs}{A character vector for cell features to be used in the calculation. If missing, all features are used.}
  \item{distMethod}{A character string indicating which distance method should be used. This must be (an abbreviation of) one of the strings 'manhattan', 'euclidean', 'correlation' or 'mahalanobis'.}
}

\details{
  For each well, this function collects features of all cells, and performs a Kolmogorov-Smirnov test for each feature against the corresponding cell features from negative control wells, with the function \code{ks.test} of the \code{stats} package. The Kolmogorov-Smirnov statistics are collected for all wells and passed to \code{PDMByWellAvg} to calculate the phenotypic distance matrix.
}

\value{
  A symmetric distance matrix with the number of rows equaling to the length of \code{unames}.
}

\seealso{\code{ks.test}, \code{PDMByWellAvg}}

\author{Xian Zhang}

\examples{

  library('phenoDist')

  ## load the imageHTS object
  load(system.file('kimorph', 'kimorph.rda', package='phenoDist'))
  x@localPath <- file.path(tempdir(), 'kimorph')

  ## calculate pair-wise svm distance matrix
  load(system.file('kimorph', 'selectedFtrs.rda', package='phenoDist'))
  pdm <- PDMByKS(x, unames=getUnames(x,plate=1, row=2:3, col=3), neg='rluc', selectedCellFtrs=selectedCellFtrs, distMethod='euclidean')
  pdm

}
