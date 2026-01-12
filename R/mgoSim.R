#' Given two sets of GO terms, calculate their semantic similarity.
#'
#'
#' @title Semantic similarity between two GO term sets
#' @param GO1 A set of GO terms
#' @param GO2 Another set of GO terms
#' @template params-measure-combine
#' @return similarity
#' @seealso [goSim()] [mgoSim()] [geneSim()] [mgeneSim()] [clusterSim()] [mclusterSim()]
#' @export
#' @examples
#' d <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)
#' go1 <- c("GO:0004022", "GO:0004024", "GO:0004023")
#' go2 <- c("GO:0009055", "GO:0020037")
#' mgoSim("GO:0003824", go2, semData = d, measure = "Wang")
#' mgoSim(go1, go2, semData = d, measure = "Wang")
#' @author Guangchuang Yu <https://yulab-smu.top>
mgoSim <- function(GO1, GO2, semData, measure="Wang", combine="BMA") {
    scores <- termSim(GO1, GO2, semData, method = measure)
    res <- combineScores(scores, combine)
    return(round(res, digits = 3))
}
