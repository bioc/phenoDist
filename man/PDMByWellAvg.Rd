\name{PDMByWellAvg}
\alias{PDMByWellAvg}
\title{ Compute phenotypic distance matrix with well features }

\description{
  This function computes the phenotypic distance matrix, with cell features averaged by well.
}

\usage{
  PDMByWellAvg(profiles, selectedWellFtrs, transformMethod=c('none', 'scale', 'PCA'), distMethod=c('manhattan', 'euclidean', 'correlation', 'mahalanobis'), nPCA)
}

\arguments{
  \item{profiles}{A data frame, containing the phenotypic profiles, as returned from the \code{summarizeWells} function of the \code{imageHTS} package.}
  \item{selectedWellFtrs}{A character vector indicating well features to be used in the calculation. If missing, all features are used.}
  \item{transformMethod}{A character string indicating which transformation method should be used. This must be (an abbreviation of) one of the strings 'none', 'scale' or 'PCA'.}
  \item{distMethod}{A character string indicating which distance method should be used. This must be (an abbreviation of) one of the strings 'manhattan', 'euclidean', 'correlation' or 'mahalanobis'.}
  \item{nPCA}{An integer scalar for the number of PCA dimensions to be used in the calculation.}
}

\details{
  Pair-wise phenotypic distance measurements of the treatments in screen results in a phenotypic distance matrix. The features stored in \code{profiles} are transformed with the \code{transformMethod} and the distance matrix is calculated with the \code{distMethod}.
}

\value{
  A symmetric distance matrix with dimensions equaling to the number of rows of \code{profiles}.
}

\seealso{\code{summarizeWells}}

\author{Xian Zhang}

\examples{

  library('phenoDist')

  ## load the imageHTS object
  load(system.file('kimorph', 'kimorph.rda', package='phenoDist'))
  x@localPath <- file.path(tempdir(), 'kimorph')

  ## calculate pair-wise svm distance matrix
  load(system.file('kimorph', 'selectedFtrs.rda', package='phenoDist'))
  pdm <- PDMByWellAvg(profiles=summarizeWells(x, getUnames(x,plate=1, row=2:3, col=3), 'conf/featurepar.txt'), selectedWellFtrs, transformMethod='scale', distMethod='euclidean')
  pdm
}

