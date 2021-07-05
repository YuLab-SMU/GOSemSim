##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",
                "For help: https://yulab-smu.top/biomedical-knowledge-mining-book/", "\n\n")
  
  citation <- paste0("If you use ", pkgname, " in published research, please cite:\n",

                     '\033[36m', '-', '\033[39m ',
                     "Guangchuang Yu. ",
                     "Gene Ontology Semantic Similarity Analysis Using GOSemSim. ",
                     "In: Kidder B. (eds) Stem Cell Transcriptional Networks. ",
                     "Methods in Molecular Biology, 2020, 2117:207-215. ",
                     "Humana, New York, NY. ",
                     "doi:10.1007/978-1-0716-0301-7_11\n",

                     '\033[36m', '-', '\033[39m ',
                     "Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu, Shengqi Wang. ",
                     "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products ",
                     "Bioinformatics 2010, 26(7):976-978. ",
                     "doi:10.1093/bioinformatics/btq064\n\n"
                     )
  
  packageStartupMessage(paste0(msg, citation))
}

