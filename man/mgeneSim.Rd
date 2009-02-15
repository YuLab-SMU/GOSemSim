\name{mgeneSim}
\alias{mgeneSim}
\title{Pairwise Semantic Similarity for a List of Genes}
\description{
  Given a list of genes, this function calculates pairwise semantic similarity.
}
\usage{
mgeneSim(genes, ont = "MF", drop= "IEA",  measure="Resnik")
}
\arguments{
  \item{genes}{ a list of entrez gene IDs. }
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}  
  \item{measure} {One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
}
\value{
  \item{simMatrix }{A Matrix of Pairwise Semantic Similarity}
}
\seealso{
	\code{\link{goSim}}
	\code{\link{mgoSim}}
	\code{\link{geneSim}}
}
\examples{
	mgeneSim(c("835", "5261","241", "934"), ont="MF", measure="Wang")
}
\references{
Wang et al.(2007)
A new method to measure the semantic similarity of go terms.
Bioinformatics (Oxford, England), 23:0 1274--81, May  2007. ISSN 1460-2059.
http://www.ncbi.nlm.nih.gov/pubmed/17344234.
PMID: 17344234.
}
\keyword{ manip}