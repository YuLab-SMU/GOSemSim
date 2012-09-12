##'Semantic Similarity Between two GO terms lists
##'
##'Given two GO term sets, this function will calculate the semantic similarity
##'between them, and return their semantic similarity
##'
##'
##'@param GO1 A set of go terms.
##'@param GO2 Another set of go terms.
##'@param ont One of "MF", "BP", and "CC" subontologies.
##'@param organism One of "anopheles", "arabidopsis", "bovine", "canine",
##'"chicken", "chimp", "coelicolor", "ecolik12", "ecsakai", "fly", "human",
##'"malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast" and
##'"zebrafish".
##'@param measure One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.
##'@param combine One of "max", "average", "rcmax", "BMA" methods, for combining
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
##'	go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
##'	go2 <- c("GO:0009055", "GO:0020037")
##'	mgoSim("GO:0003824", go2, measure="Wang")
##'	mgoSim(go1, go2, ont="MF", organism="human", measure="Wang")
##'
mgoSim <- function(GO1, GO2, ont="MF", organism="human", measure="Wang", combine="BMA"){
    scores <- termSim(GO1, GO2, ont=ont, organism=organism, method=measure)
    res <- combineScores(scores, combine)
    return(round(res, digits=3))
}
