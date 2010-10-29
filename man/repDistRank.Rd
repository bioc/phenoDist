\name{repDistRank}
\alias{repDistRank}
\title{ repDistRank }

\description{
  For each treatment, this function calculates the ranking of its replicate in terms of phenotypic distance.
}

\usage{
repDistRank(x, distMatrix)
}

\arguments{
  \item{x}{An \code{imageHTS} object.}
  \item{distMatrix}{A pair-wise distance matrix or \code{dist} object.}
}

\details{
  For each treatment, this function ranks the distance measure to its replicate among the distance measurements to all other treatments.
}

\value{
  A numeric vector with length equal to the dimension of the distance matrix.
}

\seealso{\code{PDMByWellAvg}}

\author{Xian Zhang}

\examples{

## see distToNeg.

}
