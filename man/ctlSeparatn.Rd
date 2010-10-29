\name{ctlSeparatn}
\alias{ctlSeparatn}
\title{ Calculate control separation }

\description{
  This function calculates the Z'-factor between negative and positive controls, based on phenotypic distance measurements.
}

\usage{
ctlSeparatn(x, pheno, neg='rluc', pos='ubc', ...)
}

\arguments{
  \item{x}{An \code{imageHTS} object.}
  \item{pheno}{A numeric vector for distances to negative controls, as returned from \code{distToNeg}.}
  \item{neg}{A character string to identify negative controls.}
  \item{pos}{A character string to identify positive controls.}
  \item{...}{Additional arguments to be passed to the \code{Zprime} function of the \code{imageHTS} package.}
}

\details{
  Phenotypes of negative and positive controls are defined as their phenotypic distance to negative controls. The Z'-factor is a metric measuring the separation between negative controls and positive controls (Zhang et al. 1999). This function calls the \code{Zprime} function of the \code{imageHTS} package to calculate the Z'-factor. Please see \code{Zprime} for detailed description of the definition and calculation of the Z'-factor.
}

\value{
  The Z'-factor value.
}

\references{
 J. H. Zhang, T. D. Chung, K. R. Oldenburg. A Simple Statistical Parameter for Use in Evaluation and Validation of High Throughput Screening Assays. J Biomol Screening, 1999.
}

\seealso{\code{distToNeg}, \code{Zprime}}

\author{Xian Zhang}

\examples{

## see distToNeg.
}
