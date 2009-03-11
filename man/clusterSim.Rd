\name{clusterSim}
\alias{clusterSim}
\title{Semantic Similarity Between Two Gene Clusters}
\description{
  Given two gene clusters, this function calculates semantic similarity between them.
}
\usage{
clusterSim(cluster1, cluster2, ont = "MF",  organism="human", measure="Wang", drop= "IEA")
}
\arguments{
  \item{cluster1}{A set of gene IDs. }
  \item{cluster2}{Another set of gene IDs. }  
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{measure}{One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
  \item{organism}{One of "human", "rat", "mouse", "fly" and "yeast".}  
  \item{drop}{A set of evidence codes based on which certain annotations are dropped. Use NULL to keep all GO annotations.}  
}
\value{
  \item{sim}{Semantic Similarity.}
}
\seealso{
	\code{\link{goSim}}
	\code{\link{mgoSim}}
	\code{\link{geneSim}}
	\code{\link{mgeneSim}}
	\code{\linke{mclusterSim}}
}
\examples{
	cluster1 <- c("snR67","snR40","snR48", "snR17a","snR8")
	cluster2 <- c("YOR251C", "YPR137C-B", "YPR010C", "YPR072W")
	clusterSim(cluster1, cluster2, ont="MF", organism="yeast", measure="Wang")
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
