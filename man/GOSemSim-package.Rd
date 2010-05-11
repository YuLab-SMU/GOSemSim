\name{GOSemSim-package}
\alias{GOSemSim-package}
\alias{GOSemSim}
\docType{package}
\title{
	  Gene Ontology-based Sematic Similarity Measures
}
\description{
	Implementation of semantic similarity measures to estimate the functional similarities among Gene Ontology terms and gene products
}

\details{
	Quantitative measure of functional similarities among gene products is important for post-genomics study. and widely used in gene function prediction, cluster analysis and pathway modeling.
	This package is designed to estimate the GO terms' and genes' semantic similarities. 
	Implemented five methods proposed by Resnik, Schlicker, Jiang, Lin and Wang respectively. 
	Support many species, including Anopheles, Arabidopsis, Bovine, Canine, Chicken, Chimp, E coli strain K12 and strain Sakai, Fly, Human, Malaria, Mouse, Pig, Rhesus, Rat, Worm, Xenopus, Yeast, Zebrafish.

\tabular{ll}{
Package: \tab GOSemSim\cr
Type: \tab Package\cr
Version: \tab 1.8.0\cr
Date: \tab 5-10-2010\cr
biocViews:\tab GO, Clustering\cr
Depends:\tab GO.db\cr
Suggests:\tab Clustering, org.Hs.eg.db\cr
License: \tab GPL Version 2\cr
}
}

\author{
Guangchuang Yu

Maintainer: Guangchuang Yu <guangchuangyu@gmail.com>
}
\seealso{
	\code{\link{goSim}}
	\code{\link{mgoSim}}
	\code{\link{geneSim}}
	\code{\link{mgeneSim}}
	\code{\link{clusterSim}}
	\code{\link{mclusterSim}}	
}

\keyword{ package }
