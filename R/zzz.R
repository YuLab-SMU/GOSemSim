##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0(pkgname, " v", pkgVersion, "  ",
                "For help: https://guangchuangyu.github.io/", pkgname, "\n\n")
  
  citation <- paste0("If you use ", pkgname, " in published research, please cite:\n",
                     "Guangchuang Yu, Fei Li, Yide Qin, Xiaochen Bo, Yibo Wu, Shengqi Wang. ",
                     "GOSemSim: an R package for measuring semantic similarity among GO terms and gene products ",
                     "Bioinformatics 2010, 26(7):976-978")
  
  packageStartupMessage(paste0(msg, citation))
}
