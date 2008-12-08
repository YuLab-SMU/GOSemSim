\name{mgoSim}
\alias{mgoSim}
\title{Semantic Similarity Between two GO terms lists}
\description{
 Given two GO term sets, this function will calculate the semantic similarity between them,  and return their semantic similarity
}
\usage{
mgoSim(GO1, GO2, ont = "MF", weight.isa = 0.8, weight.partof = 0.6)
}

\arguments{
  \item{GO1}{ a set of go terms. }
  \item{GO2}{ another set of go terms.}
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{weight.isa} {the contribution factor for 'is-a' relation}
  \item{weight.partof} {the contribution factor for 'part-of' relation}
}
\value{
  \item{mgoSim }{Semantic similarity.}
}
\seealso{
	\code{\link{goSim}}
	\code{\link{geneSim}}
	\code{\link{mgeneSim}}
}
\examples{
	go1 <- c("GO:0004022","GO:0004024","GO:0004023")
	go2 <- c("GO:0009055","GO:0020037")
	mgoSim("GO:0003824", go2)
	mgoSim(go1, go2, ont="MF")
}
\references{
Wang et al.(2007)
A new method to measure the semantic similarity of go terms.
Bioinformatics (Oxford, England), 23:0 1274--81, May  2007. ISSN 1460-2059.
http://www.ncbi.nlm.nih.gov/pubmed/17344234.
PMID: 17344234.
}
\keyword{ manip }