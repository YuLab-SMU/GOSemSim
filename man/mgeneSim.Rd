\name{mgeneSim}
\alias{mgeneSim}
\title{Pairwise Semantic Similarity for a List of Genes}
\usage{
  mgeneSim(genes, ont = "MF", organism = "human",
    measure = "Wang", drop = "IEA", combine = "BMA",
    verbose = TRUE)
}
\arguments{
  \item{genes}{A list of entrez gene IDs.}

  \item{ont}{One of "MF", "BP", and "CC" subontologies.}

  \item{measure}{One of "Resnik", "Lin", "Rel", "Jiang" and
  "Wang" methods.}

  \item{organism}{One of "anopheles", "arabidopsis",
  "bovine", "canine", "chicken", "chimp", "coelicolor",
  "ecolik12", "ecsakai", "fly", "human", "malaria",
  "mouse", "pig", "rat", "rhesus", "worm", "xenopus",
  "yeast" and "zebrafish".}

  \item{drop}{A set of evidence codes based on which
  certain annotations are dropped. Use NULL to keep all GO
  annotations.}

  \item{combine}{One of "max", "average", "rcmax", "BMA"
  methods, for combining semantic similarity scores of
  multiple GO terms associated with protein or multiple
  proteins assiciated with protein cluster.}

  \item{verbose}{show progress bar or not.}
}
\value{
  similarity matrix
}
\description{
  Given a list of genes, this function calculates pairwise
  semantic similarities.
}
\examples{
mgeneSim(c("835", "5261","241"), ont="MF", organism="human", measure="Wang")
}
\references{
  Yu et al. (2010) GOSemSim: an R package for measuring
  semantic similarity among GO terms and gene products
  \emph{Bioinformatics} (Oxford, England), 26:7 976--978,
  April 2010. ISSN 1367-4803
  \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
  PMID: 20179076
}
\seealso{
  \code{\link{goSim}} \code{\link{mgoSim}}
  \code{\link{geneSim}} \code{\link{clusterSim}}
  \code{\link{mclusterSim}}
}
\keyword{manip}

