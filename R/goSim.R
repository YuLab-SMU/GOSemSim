##'Semantic Similarity Between Two GO Terms
##'
##'Given two GO IDs, this function calculates their semantic similarity.
##'
##'
##'@param GOID1 GO ID 1.
##'@param GOID2 GO ID 2.
##'@param ont One of "MF", "BP", and "CC" subontologies.
##'@param organism One of "anopheles", "arabidopsis", "bovine", "canine",
##' "chicken", "chimp", "coelicolor", "ecolik12", "ecsakai", "fly", "human",
##' "malaria", "mouse", "pig", "rat","rhesus", "worm", "xenopus", "yeast" and
##' "zebrafish".
##'@param measure One of "Resnik", "Lin", "Rel", "Jiang" and "Wang" methods.
##'@return similarity
##'@seealso \code{\link{mgoSim}}
##' \code{\link{geneSim}}
##' \code{\link{mgeneSim}}
##' \code{\link{clusterSim}}
##' \code{\link{mclusterSim}}
##'@references Yu et al. (2010) GOSemSim: an R package for measuring semantic
##'similarity among GO terms and gene products \emph{Bioinformatics} (Oxford,
##'England), 26:7 976--978, April 2010. ISSN 1367-4803
##'\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
##'PMID: 20179076
##'@keywords manip
##' @export
##'@examples
##'
##'	goSim("GO:0004022", "GO:0005515", ont="MF", measure="Wang")
##'
goSim <- function(GOID1, GOID2, ont="MF", organism="human", measure="Wang"){
    res <- termSim(GOID1, GOID2, ont=ont, organism=organism, method=measure)
    res <- as.numeric(res)
    return(round(res,digits=3))
}

