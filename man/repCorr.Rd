\name{repCorr}
\alias{repCorr}
\title{ repCorr }

\description{
  This function calculates the correlation between replicate phenotypes.
}

\usage{
repCorr(x, pheno, ...)
}

\arguments{
  \item{x}{An \code{imageHTS} object.}
  \item{pheno}{A numeric vector for distance measures to negative control, as returned from \code{distToNeg}.}
  \item{...}{Additional arguments to be passed to the \code{cor} function of the \code{stats} pacakge.}

}

\value{
  Correlation coefficient returned from the \code{cor} function of the \code{stats} pacakge. 
}

\seealso{\code{cor}, \code{distToNeg}}

\author{Xian Zhang}

\examples{

## see distToNeg.

}
