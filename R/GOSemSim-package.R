##'Gene Ontology-based Sematic Similarity Measures
##'
##'Implementation of semantic similarity measures to estimate the functional
##'similarities among Gene Ontology terms and gene products
##'
##'Quantitative measure of functional similarities among gene products is
##'important for post-genomics study. and widely used in gene function
##'prediction, cluster analysis and pathway modeling.  This package is designed
##'to estimate the GO terms' and genes' semantic similarities.  Implemented five
##'methods proposed by Resnik, Schlicker, Jiang, Lin and Wang respectively.
##'Support many species, including Anopheles, Arabidopsis, Bovine, Canine,
##'Chicken, Chimp, E coli strain K12 and strain Sakai, Fly, Human, Malaria,
##'Mouse, Pig, Rhesus, Rat, Worm, Xenopus, Yeast, Zebrafish.
##'
##'\tabular{ll}{ Package: \tab GOSemSim\cr Type: \tab Package\cr Version: \tab
##'1.24.0\cr Date: \tab 09-11-2012\cr biocViews:\tab GO, Clustering, Pathways,
##'Anopheles_gambiae, Arabidopsis_thaliana, Bos_taurus, Caenorhabditis_elegans,
##'Canis_familiaris, Danio_rerio, Drosophila_melanogaster, Escherichia_coli,
##'Gallus_gallus, Homo_sapiens, Mus_musculus, Pan_troglodytes,
##'Plasmodium_falciparum, Rattus_norvegicus, Saccharomyces_cerevisiae,
##'Streptomyces_coelicolor, Sus_scrofa, Xenopus_laevis\cr Depends:\tab \cr
##'Imports: \tab methods, AnnotationDbi, GO.db, org.Hs.eg.db, org.Ag.eg.db,
##'org.At.tair.db, org.Bt.eg.db, org.Ce.eg.db, org.Cf.eg.db, org.Dm.eg.db,
##'org.Dr.eg.db, org.EcK12.eg.db, org.EcSakai.eg.db, org.Gg.eg.db, org.Mm.eg.db,
##'org.Mmu.eg.db, org.Pf.plasmo.db, org.Pt.eg.db, org.Rn.eg.db, org.Sc.sgd.db,
##'org.Sco.eg.db, org.Ss.eg.db, org.Tgondii.eg.db, org.Xl.eg.db\cr
##'Suggests:\tab clusterProfiler\cr
##'License: \tab Artistic-2.0\cr }
##'
##'@name GOSemSim-package
##'@aliases GOSemSim-package GOSemSim
##'@docType package
##'@author Guangchuang Yu
##'
##'Maintainer: Guangchuang Yu <guangchuangyu@@gmail.com>
##'@seealso \code{\link{goSim}} \code{\link{mgoSim}} \code{\link{geneSim}}
##'\code{\link{mgeneSim}} \code{\link{clusterSim}} \code{\link{mclusterSim}}
##'@references Yu et al. (2010) GOSemSim: an R package for measuring semantic
##'similarity among GO terms and gene products \emph{Bioinformatics} (Oxford,
##'England), 26:7 976--978, April 2010. ISSN 1367-4803
##'\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
##'PMID: 20179076
##'@keywords package
NULL





##'Information content of GO terms
##'
##'These datasets are the information contents of GOterms.
##'
##'
##'@name IC
##'@aliases IC Info_Contents_anopheles_MF Info_Contents_anopheles_BP
##'Info_Contents_anopheles_CC Info_Contents_arabidopsis_MF
##'Info_Contents_arabidopsis_BP Info_Contents_arabidopsis_CC
##'Info_Contents_bovine_MF Info_Contents_bovine_BP Info_Contents_bovine_CC
##'Info_Contents_canine_MF Info_Contents_canine_BP Info_Contents_canine_CC
##'Info_Contents_chicken_MF Info_Contents_chicken_BP Info_Contents_chicken_CC
##'Info_Contents_chimp_MF Info_Contents_chimp_BP Info_Contents_chimp_CC
##'Info_Contents_coelicolor_MF Info_Contents_coelicolor_BP
##'Info_Contents_coelicolor_CC Info_Contents_ecsakai_MF Info_Contents_ecsakai_BP
##'Info_Contents_ecsakai_CC Info_Contents_ecolik12_MF Info_Contents_ecolik12_BP
##'Info_Contents_ecolik12_CC Info_Contents_gondii_CC Info_Contents_gondii_MF
##'Info_Contents_gondii_BP Info_Contents_fly_MF Info_Contents_fly_BP
##'Info_Contents_fly_CC Info_Contents_human_MF Info_Contents_human_BP
##'Info_Contents_human_CC Info_Contents_malaria_MF Info_Contents_malaria_BP
##'Info_Contents_malaria_CC Info_Contents_mouse_MF Info_Contents_mouse_BP
##'Info_Contents_mouse_CC Info_Contents_worm_MF Info_Contents_worm_BP
##'Info_Contents_worm_CC Info_Contents_pig_MF Info_Contents_pig_BP
##'Info_Contents_pig_CC Info_Contents_rat_MF Info_Contents_rat_BP
##'Info_Contents_rat_CC Info_Contents_rhesus_MF Info_Contents_rhesus_BP
##'Info_Contents_rhesus_CC Info_Contents_xenopus_MF Info_Contents_xenopus_BP
##'Info_Contents_xenopus_CC Info_Contents_yeast_MF Info_Contents_yeast_BP
##'Info_Contents_yeast_CC Info_Contents_zebrafish_MF Info_Contents_zebrafish_BP
##'Info_Contents_zebrafish_CC
##'@docType data
##'@references Yu et al. (2010) GOSemSim: an R package for measuring semantic
##'similarity among GO terms and gene products \emph{Bioinformatics} (Oxford,
##'England), 26:7 976--978, April 2010. ISSN 1367-4803
##'\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
##'PMID: 20179076
##'@keywords datasets
NULL



