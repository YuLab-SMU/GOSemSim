\name{geneSim}
\alias{geneSim}
\title{Semantic Similarity Between two Genes}
\description{
 Given two genes, this function will calculate the semantic similarity between them,  and return their semantic similarity and the corresponding GO terms
}
\usage{
	geneSim(gene1, gene2, ont = "MF", organism="human",measure="Wang",  drop= "IEA")
}

\arguments{
  \item{gene1}{ Entrez gene id. }
  \item{gene2}{ Another entrez gene id.}
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{organism}{One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish".}
  \item{measure}{One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}
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
	\code{\link{clusterSim}}
	\code{\link{mclusterSim}}	
}
\examples{
	geneSim("241", "2561", ont="MF", organism="human", measure="Wang")
	geneSim("snR18", "YPR062W", ont="MF", organism="yeast", measure="Wang")
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
