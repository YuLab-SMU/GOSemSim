\name{GOSet}
\docType{methods}
\alias{GOSet}
\alias{GOSet-methods}
\title{Methods for GOSet class}
\description{
  Use \code{GOSet} to construct parameter class for calculating Semantic similarity.
}
\section{Slots}{
  \describe{
    \item{\code{GOSet1}:}{containing a vector of GO Terms.}
	\item{\code{GOSet2}:}{containing a vector of GO Terms.}
  }
}
\section{Methods}{
  See documentation for \code{\linkS4class{GOSet}} for examples.
}
\seealso{
  \code{\link{GOSet-class}}
}
\examples{
## Construct GOSet class
go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
go2 <- c("GO:0009055", "GO:0020037")
GOSet <- new("GOSet", GOSet1=go1, GOSet2=go2)
}
\keyword{methods}
