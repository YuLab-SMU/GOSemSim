\name{mgoSim}
\alias{mgoSim}
\title{Semantic Similarity Between two GO terms lists}
\description{
 Given two GO term sets, this function will calculate the semantic similarity between them,  and return their semantic similarity
}
\usage{
mgoSim(GO1, GO2, ont="MF", organism="human", measure="Wang")
}

\arguments{
  \item{GO1}{A set of go terms. }
  \item{GO2}{Another set of go terms.}
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{organism}{One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish".}
  \item{measure}{One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
}
\value{
  \item{mgoSim }{Semantic similarity.}
}
\seealso{
	\code{\link{goSim}}
	\code{\link{geneSim}}
	\code{\link{mgeneSim}}
	\code{\link{clusterSim}}
	\code{\link{mclusterSim}}	
}
\examples{
	go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
	go2 <- c("GO:0009055", "GO:0020037")
	mgoSim("GO:0003824", go2, measure="Wang")
	mgoSim(go1, go2, ont="MF", organism="human", measure="Wang")
}
\references{
Wang et al.(2007)
A new method to measure the semantic similarity of go terms
\emph{Bioinformatics} (Oxford, England), 23:0 1274--81, May  2007. ISSN 1460-2059
\url{http://www.ncbi.nlm.nih.gov/pubmed/17344234}
PMID: 17344234

P W Lord et al.(2003)
Semantic similarity measures as tools for exploring the gene ontology
\url{Pacific Symposium on Biocomputing} Pacific Symposium on Biocomputing, 2003:601-12, ISSN 1793-5091
\url{http://www.ncbi.nlm.nih.gov/pubmed/12603061}
PMID: 12603061
}
\keyword{ manip }
