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
Version: \tab 1.9.9\cr
Date: \tab 03-18-2011\cr
biocViews:\tab GO, Clustering, Pathways, Anopheles_gambiae, Arabidopsis_thaliana, Bos_taurus, Caenorhabditis_elegans, Canis_familiaris, Danio_rerio, Drosophila_melanogaster, Escherichia_coli, Gallus_gallus, Homo_sapiens, Mus_musculus, Pan_troglodytes, Plasmodium_falciparum, Rattus_norvegicus, Saccharomyces_cerevisiae, Streptomyces_coelicolor, Sus_scrofa, Xenopus_laevis\cr
Depends:\tab \cr
Imports: \tab methods, AnnotationDbi, GO.db, org.Hs.eg.db, org.Ag.eg.db, org.At.tair.db, org.Bt.eg.db, org.Ce.eg.db, org.Cf.eg.db, org.Dm.eg.db, org.Dr.eg.db, org.EcK12.eg.db, org.EcSakai.eg.db, org.Gg.eg.db, org.Mm.eg.db, org.Mmu.eg.db, org.Pf.plasmo.db, org.Pt.eg.db, org.Rn.eg.db, org.Sc.sgd.db, org.Sco.eg.db, org.Ss.eg.db, org.Xl.eg.db\cr
Suggests:\tab Clustering\cr
License: \tab GPL Version 2\cr
}
}

\author{
Guangchuang Yu

Maintainer: Guangchuang Yu <guangchuangyu@gmail.com>
}
\seealso{
  \linkS4class{GOSet},
  \linkS4class{GeneSet}
  \linkS4class{GeneClusterSet}
}

\keyword{ package }
