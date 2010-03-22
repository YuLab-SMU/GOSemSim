\name{mclusterSim}
\alias{mclusterSim}
\title{Pairwise Semantic Similarities for a List of Gene Clusters}
\description{
  Given a list of gene clusters, this function calculates pairwise semantic similarities.
}
\usage{
mclusterSim(clusters, ont = "MF",  organism="human", measure="Wang", drop= "IEA")
}
\arguments{
  \item{clusters}{A list of gene clusters. }
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{measure}{One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
  \item{organism}{One of "anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and "zebrafish".}
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}  
}
\value{
  \item{simmat }{A Matrix of Pairwise Semantic Similarities}
}
\seealso{
	\code{\link{goSim}}
	\code{\link{mgoSim}}
	\code{\link{geneSim}}
	\code{\link{mgeneSim}}
	\code{\link{clusterSim}}
}
\examples{
	cluster1 <- c("835", "5261","241", "994", "514", "517", "533")
	cluster2 <- c("578","582", "583", "400", "409", "411")
	cluster3 <- c("307", "308", "317", "321", "506", "540", "378", "388", "396")
	clusters <- list(a=cluster1, b=cluster2, c=cluster3)
	mclusterSim(clusters, ont="MF", organism="yeast", measure="Wang")
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
