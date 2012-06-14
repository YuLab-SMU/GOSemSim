\name{GeneSet-class}
\docType{class}
\alias{GeneSet-class}

\title{Class "GeneSet"}
\description{
  A \code{GeneSet} contains two sets of Gene identifiers.
}

\section{Slots}{
  \describe{
    \item{\code{GeneSet1}:}{containing a vector of Gene identifiers.}
	\item{\code{GeneSet2}:}{containing a vector of Gene identifiers.}
  }
}
\section{Methods}{
  GO semantic similarity measure :
  \describe{
    \item{sim}{\code{signature(object = "GeneSet", params="Params", value = "numeric")}}
  }
}
\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\linkS4class{GeneClusterSet}}
  \code{\linkS4class{Params}}
}
\examples{
	## Setting Parameters...
	params <- new("Params", ontology="MF", organism="human", method="Wang", combine="BMA")
	## Setting GeneSet...
	gs1 <- c("835", "5261","241", "994", "514", "533")
	gs2 <- c("578","582", "400", "409", "411")
	gs <- new("GeneSet", GeneSet1=gs1, GeneSet2=gs2)
	## Calculating Gene Semantic Similarities...
	#sim(gs, params)
}
\references{
  Yu et al. (2010)
  GOSemSim: an R package for measuring semantic similarity among GO terms and gene products
  \emph{Bioinformatics} (Oxford, England), 26:7 976--978, April  2010. ISSN 1367-4803
  \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
  PMID: 20179076
}
\keyword{classes}
