\name{GOIdentifiers}
\docType{methods}
\alias{GOIdentifiers}
\alias{GOIdentifiers-methods}
\title{Methods for GOIdentifiers class}
\description{
  Use \code{GOIdentifiers} to construct parameter class for calculating Semantic similarity.
}
\section{Slots}{
  \describe{
    \item{\code{GOSet1}:}{containing a vector of GO Terms.}
	\item{\code{GOSet2}:}{containing a vector of GO Terms.}
  }
}
\section{Methods}{
  See documentation for \code{\linkS4class{GOIdentifiers}} for examples.
}
\seealso{
  \code{\link{GOIdentifiers-class}}
}
\examples{
## Construct GOIdentifiers class
go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
go2 <- c("GO:0009055", "GO:0020037")
GOIdentifiers <- new("GOIdentifiers", GOSet1=go1, GOSet2=go2)
}
\keyword{methods}
