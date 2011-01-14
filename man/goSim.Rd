\name{goSim}
\alias{goSim}
\title{Semantic Similarity Between Two GO Terms}
\description{
Given two GO IDs, this function calculates their semantic similarity.
}
\usage{
goSim(GOID1, GOID2, ont = "MF", organism="human", measure="Wang")
}
\arguments{
  \item{GOID1}{GO ID 1.}
  \item{GOID2}{GO ID 2.}
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{organism}{One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "coelicolor", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish".}
  \item{measure}{One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
}
\value{
  \item{goSim}{Semantic similarity.}
}

\seealso{
  \code{\link{mgoSim}}
  \code{\link{geneSim}}
  \code{\link{mgeneSim}}
  \code{\link{clusterSim}}
  \code{\link{mclusterSim}}	
}
\examples{
	goSim("GO:0043121", "GO:0019838", measure="Wang")
	goSim("GO:0043121", "GO:0019838", ont="MF", organism="human", measure="Wang")
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
\keyword{manip}
