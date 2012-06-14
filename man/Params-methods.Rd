\name{Params}
\docType{methods}
\alias{Params}
\alias{Params-methods}
\title{Methods for Params class}
\description{
  Use \code{Params} to construct parameter class for calculating Semantic similarity.
}

\section{Slots}{
  \describe{
    \item{\code{ontology}:}{An argument determining which ontology were to be measured.}
    \item{\code{organism}:}{Setting organism for mapping Gene IDs to GO Terms.}
    \item{\code{method}:}{Method for calculating semantic similarity.}
    \item{\code{combine}:}{Combine method for combining multiple GO semantic scores to one score.}
    \item{\code{dropCodes}:}{dropCodes for mapping Gene to GO.}
  }
}

\section{Methods}{
  See documentation for \code{\linkS4class{Params}} for examples.
}
\seealso{
  \code{\link{Params-class}}
}
\examples{
## Construct Params class
params <- new("Params", ontology="MF", organism="human", method="Wang")

## Setting params
setOntology(params) <- "BP"
setCombineMethod="rcmax"

## Accessing slots.
params["ontology"]
params["organism"]
}
\references{
  Yu et al. (2010)
  GOSemSim: an R package for measuring semantic similarity among GO terms and gene products
  \emph{Bioinformatics} (Oxford, England), 26:7 976--978, April  2010. ISSN 1367-4803
  \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
  PMID: 20179076
}
\keyword{methods}
