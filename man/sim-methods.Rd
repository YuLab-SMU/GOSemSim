\name{sim-methods}
\docType{methods}
\alias{sim}
\alias{sim-methods}
\alias{sim,GOSet-method}
\alias{sim,GeneSet-method}
\alias{sim,GeneClusterSet-method}

\title{Methods for calculating semantic similarity}
\description{
  This generic and methods calculating semantic similarities among GO terms, gene sets and gene clusters.
}
\section{Methods}{
  \describe{
    Defined methods include:
    \item{
      \code{signature(object = "GOSet", params="Params")},
      \code{signature(object = "GeneSet", params="Params")}
	  \code{signature(object = "GeneClusterSet", params="Params")}
    }{
      These methods calculating semantic similarities about \code{GOSet},
      \code{GeneSet}, and \code{GeneClusterSet}.
    }
  }
}
\seealso{
  \code{\link{GOSet-class}}
  \code{\link{GeneSet-class}}
  \code{\link{GeneClusterSet-class}}
}
\references{
  Yu et al. (2010)
  GOSemSim: an R package for measuring semantic similarity among GO terms and gene products
  \emph{Bioinformatics} (Oxford, England), 26:7 976--978, April  2010. ISSN 1367-4803
  \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
  PMID: 20179076
}
\keyword{methods}
