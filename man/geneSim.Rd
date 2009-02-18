\name{geneSim}
\alias{geneSim}
\title{Semantic Similarity Between two Genes}
\description{
 Given two genes, this function will calculate the semantic similarity between them,  and return their semantic similarity and the corresponding GO terms
}
\usage{
	geneSim(gene1, gene2, ont = "MF",  drop= "IEA",measure="Resnik")
}

\arguments{
  \item{gene1}{ Entrez gene id. }
  \item{gene2}{ Another entrez gene id.}
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{drop}{A set of evidence codes based on which certain annotations are
    dropped. Use NULL to keep all GO annotations.} 
  \item{measure}{One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
}
\value{
  \item{geneSim }{Semantic similarity.}
  \item{GO1 }{Corresponding GO terms for gene1.}
  \item{GO2 }{Corresponding GO terms for gene2.}
}
\seealso{
	\code{\link{goSim}}
	\code{\link{mgoSim}}
	\code{\link{mgeneSim}}
}
\examples{
	geneSim("241", "2561", measure="Wang", ont="MF")
}
\references{
Wang et al.(2007)
A new method to measure the semantic similarity of go terms.
\emph{Bioinformatics} (Oxford, England), 23:0 1274--81, May  2007. ISSN
  1460-2059. 
\url{http://www.ncbi.nlm.nih.gov/pubmed/17344234}.
PMID: 17344234.
}
\keyword{ manip }
