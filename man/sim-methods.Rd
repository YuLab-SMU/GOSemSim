\name{sim-methods}
\docType{methods}
\alias{sim}
\alias{sim-methods}
\alias{sim,GOIdentifiers-method}
\alias{sim,GeneIdentifiers-method}
\alias{sim,GeneClusterSet-method}

\title{Methods for calculating semantic similarity}
\description{
  This generic and methods calculating semantic similarities among GO terms, gene sets and gene clusters.
}
\section{Methods}{
  \describe{
    Defined methods include:
    \item{
      \code{signature(object = "GOIdentifers", params="Params")},
      \code{signature(object = "GeneIdentifiers", params="Params")}
	  \code{signature(object = "GeneClusterSet", params="Params")}
    }{
      These methods calculating semantic similarities about \code{GOIdentifiers},
      \code{GeneIdentifiers}, and \code{GeneClusterSet}.
    }
  }
}
\seealso{
  \code{\link{GOIdentifiers-class}}
  \code{\link{GeneIdentifiers-class}}
  \code{\link{GeneClusterSet-class}}
}
\keyword{methods}
