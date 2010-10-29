\name{PDMByFactorAnalysis}
\alias{PDMByFactorAnalysis}
\title{ Calculate phenotypic distance matrix by factor analysis }

\description{
  This function transforms the cell features by factor analysis and computes the phenotypic distance matrix.
}

\usage{

PDMByFactorAnalysis(x, unames, selectedCellFtrs, distMethod=c('manhattan','euclidean', 'correlation','mahalanobis'), nFactors, scores=c('regression','Bartlett'), ...)
}

\arguments{
  \item{x}{An \code{imageHTS} object.}
  \item{unames}{A character vector, containing the well names from where to collect the cell features. See \code{getUnames} for details.}
  \item{selectedCellFtrs}{A character vector for cell features to be used in the calculation. If missing, all features are used.}
  \item{distMethod}{A character string indicating which distance method should be used. This must be (an abbreviation of) one of the strings 'manhattan', 'euclidean', 'correlation' or 'mahalanobis'.}
  \item{nFactors}{An integer scalar for the number of factors.}
  \item{scores}{A character string indicating the type of scores to be reported by factor analysis. This must be (an abbreviation of) one of the strings 'regression' or 'Bartlett'.}
  \item{...}{Additional arguments to be passed to the \code{factanal} function of the \code{stats} pacakge.}
}

\details{
  This function first collects individual cell features in all wells (which could be time and memory consuming), performs factor analysis on cell features and transforms cell features into a certain number of factors, and then the factors are averaged by well and passed to \code{PDMByWellAvg} to calculate the phenotypic distance matrix.
}

\value{
  A symmetric distance matrix with dimensions equaling to the length of \code{unames}.
}

\seealso{\code{factanal}, \code{PDMByWellAvg}}

\author{Xian Zhang}

\examples{
  library('phenoDist')

  ## load the imageHTS object
  load(system.file('kimorph', 'kimorph.rda', package='phenoDist'))
  x@localPath <- file.path(tempdir(), 'kimorph')
  
  ## segmentation and feature extraction
  unames <- setdiff(getUnames(x), getUnames(x, content='empty'))

  ## calculate pair-wise svm distance matrix
  load(system.file('kimorph', 'selectedFtrs.rda', package='phenoDist'))
  pdm <- PDMByFactorAnalysis(x, unames=getUnames(x, plate=1, row=2:3, col=3), selectedCellFtrs, distMethod='euclidean', nFactors=10, scores='regression')
  pdm
}

