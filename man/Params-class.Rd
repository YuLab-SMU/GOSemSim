\name{Params-class}
\docType{class}
\alias{Params-class}
\alias{[,Params,character-method}
\alias{loadICdata,Params-method}
\alias{loadAnnoPkg,Params-method}
\alias{loadGOMap,Params-method}
\alias{computeIC,Params-method}
\alias{setOntology<-,Params-method}
\alias{setOrganism<-,Params-method}
\alias{setMethod<-,Params-method}
\alias{setCombineMethod<-,Params-method}


\alias{loadICdata}
\alias{loadAnnoPkg}
\alias{loadGOMap}
\alias{computeIC}
\alias{setOntology<-}
\alias{setOrganism<-}
\alias{setMethod<-}
\alias{setCombineMethod<-}


\title{Class "Params"}
\description{
  A \code{Params} contains parameters for calculating GO semantic similarity among GO Terms or Gene Sets.
}

\section{Slots}{
  \describe{
    \item{\code{ontology}:}{one of "MF", "BP", "CC".}
	\item{\code{organism}:}{one of supported species.}
	\item{\code{method}:}{Method for calculating GO semantic similarity, one of "Resnik", "Jiang", "Lin", "Rel", "Wang".}
	\item{\code{combine}:}{Method for combining GO semantic similarity scores, one of "avg", "max", "rcmax", "rcmax.avg"}
	\item{\code{dropCodes}:}{dropCodes for mapping Gene to GO Terms.}
  }
}
\section{Methods}{
  Slot access (e.g., \code{setOntology<-}) and retrieve
  (e.g., \code{[}) :
  \describe{
    \item{setOntology<-}{\code{signature(object = "Params")}}
    \item{setOrganism<-}{\code{signature(object = "Params")}}
    \item{setMethod<-}{\code{signature(object = "Params")}}
    \item{setCombineMethod<-}{\code{signature(object = "Params", value = "character")}}
	\item{[}{\code{signature(x="Params", i="character")}: subset the Params by index (\code{i="character"})}
  }

  Loading require data:
  \describe{
    \item{loadICdata}{\code{signature(object = "Params")}}
	\item{loadAnnoPkg}{\code{signature(object = "Params")}}
    \item{loadGOMap}{\code{signature(object = "Params")}}
  }

  Useful additional methods include:
  \describe{
    \item{computeIC}{\code{signature(object = "Params")}: compute Information Content of GO}
  }
}
\author{Guangchuang Yu <guangchuangyu@gmail.com>}
\seealso{
  \code{\linkS4class{GeneClusterSet}}
  \code{\linkS4class{Params}}
}
\examples{
## Setting Parameters...
params <- new("Params", ontology="MF", organism="human", method="Wang")
}
\keyword{classes}
