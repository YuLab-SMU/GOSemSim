\name{GOIdentifiers-class}
\docType{class}
\alias{GOIdentifiers-class}

\title{Class "GOIdentifiers"}
\description{
  A \code{GOIdentifiers} contains two sets of GO identifiers.

}

\section{Slots}{
  \describe{
    \item{\code{GOSet1}:}{containing a vector of GO identifiers.}
	\item{\code{GOSet2}:}{containing a vector of GO identifiers.}
  }
}
\section{Methods}{
  GO semantic similarity measure :
  \describe{
    \item{sim}{\code{signature(object = "GOIdentifiers", params="Params")}}
  }
}
\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\linkS4class{GeneIdentifiers}}
  \code{\linkS4class{Params}}
}
\examples{
## Setting Parameters...
params <- new("Params", ontology="MF", organism="human", method="Wang")
## Setting GOIdentifiers...
go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
go2 <- c("GO:0009055", "GO:0020037")
gos <- new("GOIdentifiers", GOSet1=go1, GOSet2=go2)
## Calculating GO Semantic Similarities...
sim(gos, params)
## Setting Combine Method
setCombineMethod(params) <- "rcmax.avg"
## Combining semantic similarity scores of multiple GO terms to one score.
sim(gos, params)
}
\keyword{classes}
