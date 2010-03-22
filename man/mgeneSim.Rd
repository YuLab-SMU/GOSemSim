\name{mgeneSim}
\alias{mgeneSim}
\title{Pairwise Semantic Similarity for a List of Genes}
\description{
  Given a list of genes, this function calculates pairwise semantic similarities.
}
\usage{
mgeneSim(genes, ont = "MF",  organism="human", measure="Wang", drop= "IEA")
}
\arguments{
  \item{genes}{A list of entrez gene IDs. }
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{measure}{One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
  \item{organism}{One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish".}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}  
}
\value{
  \item{simMatrix }{A Matrix of Pairwise Semantic Similarities}
}
\seealso{
	\code{\link{goSim}}
	\code{\link{mgoSim}}
	\code{\link{geneSim}}
	\code{\link{clusterSim}}
	\code{\link{mclusterSim}}
}
\examples{
	mgeneSim(c("835", "5261","241", "994"), ont="MF", organism="human", measure="Wang")
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
\keyword{ manip}
