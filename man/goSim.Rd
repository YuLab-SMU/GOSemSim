\name{goSim}
\alias{goSim}
\title{Semantic Similarity Between Two GO Terms}
\description{
Given two GO IDs, this function calculates their semantic similarity.
}
\usage{
goSim(GOID1, GOID2, ont = "MF", measure="Resnik")
}
\arguments{
  \item{GOID1}{GO ID 1.}
  \item{GOID2}{GO ID 2.}
  \item{ont}{One of "MF", "BP", and "CC" subontologies.}
  \item{measure} {One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.}
}
\value{
  \item{goSim}{Semantic similarity.}
}

\seealso{
  \code{\link{mgoSim}}
  \code{\link{geneSim}}
  \code{\link{mgeneSim}}
}
\examples{
	goSim("GO:0043121", "GO:0019838", measure="Resnik")
	goSim("GO:0043121", "GO:0019838", ont="MF", measure="Resnik")
}
\references{
Wang et al.(2007)
A new method to measure the semantic similarity of go terms.
Bioinformatics (Oxford, England), 23:0 1274--81, May  2007. ISSN 1460-2059.
http://www.ncbi.nlm.nih.gov/pubmed/17344234.
PMID: 17344234.
}
\keyword{manip}