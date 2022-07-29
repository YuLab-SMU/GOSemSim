##'Semantic Similarity Between two GO terms lists
##'
##'Given two GO term sets, this function will calculate the semantic similarity
##'between them, and return their semantic similarity
##'
##'
##'@param GO1 A set of go terms.
##'@param GO2 Another set of go terms.
##'@param semData GOSemSimDATA object
##'@param measure One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
##'@param combine One of "max", "avg", "rcmax", "BMA" methods, for combining
##'semantic similarity scores of multiple GO terms associated with protein or
##'multiple proteins assiciated with protein cluster.
##'@return similarity
##'@seealso \code{\link{goSim}} \code{\link{geneSim}} \code{\link{mgeneSim}}
##'\code{\link{clusterSim}} \code{\link{mclusterSim}}
##'@references Yu et al. (2010) GOSemSim: an R package for measuring semantic
##'similarity among GO terms and gene products \emph{Bioinformatics} (Oxford,
##'England), 26:7 976--978, April 2010. ISSN 1367-4803
##'\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
##'PMID: 20179076
##'
##'@keywords manip
##' @export
##'@examples
##'
##'     d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
##'	go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
##'	go2 <- c("GO:0009055", "GO:0020037")
##'	mgoSim("GO:0003824", go2, semData=d, measure="Wang")
##'	mgoSim(go1, go2, semData=d, measure="Wang")
##'
mgoSim <- function(GO1, GO2, semData, measure="Wang", combine="BMA"){
    scores <- termSim(GO1, GO2, semData, method=measure)
    res <- combineScores(scores, combine)
    return(round(res, digits=3))
}
