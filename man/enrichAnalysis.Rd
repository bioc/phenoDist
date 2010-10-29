\name{enrichAnalysis}
\alias{enrichAnalysis}
\title{ enrichAnalysis }
\description{
  This function performs enrichment analysis on genes within each cluster, with genes from all clusters as the gene universe.
}
\usage{
enrichAnalysis(x, cl, terms=c('GO','KEGG'), verbose=FALSE, ...)
}
\arguments{
  \item{x}{An \code{imageHTS} object. }
  \item{cl}{A named numeric vector, with names being gene Entrez IDs and values indicating grouping. }
  \item{terms}{A character string indicating which annotation should be used. This must be (an abbreviation of) one of the strings 'GO' or 'KEGG'.}
  \item{verbose}{A logical scalar indicating whether progress should be reported. }
  \item{...}{Additional arguments to be passed to the \code{hyperGTest} function of the \code{GOstats} package. }
}
\details{
  
  Gene enrichment analysis tests whether certain gene annotations (e.g., GO terms or KEGG ids) are enriched in genes of interest, compared with the gene universe. This function is designed to analyze gene enrichment in clusters, with genes in the tested cluster being genes of interest and genes in all clusters being the gene universe.
  
  The Hypergeometric test is performed by the \code{hyperGTest} function of the \code{GOstats} package. Please refer to \code{hyperGTest} for additional arguments.
  
}
\value{
  A list of \code{HyperGResult} instances.
}

\seealso{\code{GOstats}, \code{clusterDist}}

\author{Xian Zhang}

\examples{

## see clusterDist.

}
