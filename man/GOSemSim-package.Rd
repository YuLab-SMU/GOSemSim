\name{GOSemSim-package}
\alias{GOSemSim-package}
\alias{GOSemSim}
\docType{package}
\title{
	  Gene Ontology-based Sematic Similarity Measures
}
\description{
	Implementation of semantic similarity measure to estimate the functional similarity between Gene Ontology terms and gene products
}

\details{
	Quantitative measure of functional similarity between gene products is important for post-genomics study. and widely used in gene function, cluster analysis and pathway modeling.
	This package is designed to estimate the GO terms' and genes' semantic similarities. 
	This method encode the GO term's semantic into a numeric value by aggregating the semantic contributions of their ancestor terms in the GO graph, and the functional similarities of genes is based on the semantic similarities among GO terms annotating these genes. 

\tabular{ll}{
Package: \tab GOSemSim\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 03-04-2009\cr
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
