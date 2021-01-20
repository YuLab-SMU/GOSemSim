##' Semantic Similarity Between Two GO Terms
##'
##' Given two GO IDs, this function calculates their semantic similarity.
##'
##'
##' @param GOID1 GO ID 1.
##' @param GOID2 GO ID 2.
##' @param semData GOSemSimDATA object
##' @param measure One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
##' @return similarity
##' @seealso \code{\link{mgoSim}}
##' \code{\link{geneSim}}
##' \code{\link{mgeneSim}}
##' \code{\link{clusterSim}}
##' \code{\link{mclusterSim}}
##' @references Yu et al. (2010) GOSemSim: an R package for measuring semantic
##' similarity among GO terms and gene products \emph{Bioinformatics} (Oxford,
##' England), 26:7 976--978, April 2010. ISSN 1367-4803
##' \url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
##' PMID: 20179076
##' @keywords manip
##' @export
##' @examples
##' 
##'     d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
##'	goSim("GO:0004022", "GO:0005515", semData=d, measure="Wang")
##' 
goSim <- function(GOID1, GOID2, semData, measure="Wang") {
    res <- termSim(GOID1, GOID2, semData, method=measure)
    res <- as.numeric(res)
    return(round(res,digits=3))
}

