\name{GeneClusterSet-class}
\docType{class}
\alias{GeneClusterSet-class}

\title{Class "GeneClusterSet"}
\description{
  A \code{GeneClusterSet} contains a list of Gene clusters.

}

\section{Slots}{
  \describe{
    \item{\code{GeneClusters}:}{containing a list of Gene clusters.}
  }
}
\section{Methods}{
  Gene Clusters semantic similarity measure :
  \describe{
    \item{sim}{\code{signature(object = "GeneClusterSet", params="Params", value = "numeric")}}
  }
}
\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\linkS4class{GeneSet}}
  \code{\linkS4class{Params}}
}
\examples{
	## Setting Parameters...
	params <- new("Params", ontology="MF", organism="human", method="Wang", combine="rcmax")
	## Setting GeneClusterSet...
	cluster1 <- c("835", "5261","241", "994", "514", "517", "533")
	cluster2 <- c("578","582", "583", "400", "409", "411")
	cluster3 <- c("307", "308", "317", "321", "506", "540", "378", "388", "396")
	clusters <- list(a=cluster1, b=cluster2, c=cluster3)
	geneClusters <- new("GeneClusterSet", GeneClusters=clusters)
	## Calculating Gene Clusters Semantic Similarities...
	# sim(geneClusters, params)
}
\references{
  Yu et al. (2010)
  GOSemSim: an R package for measuring semantic similarity among GO terms and gene products
  \emph{Bioinformatics} (Oxford, England), 26:7 976--978, April  2010. ISSN 1367-4803
  \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
  PMID: 20179076
}
\keyword{classes}
