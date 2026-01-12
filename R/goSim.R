#' Given two GO IDs, calculate their semantic similarity.
#'
#'
#' @title Semantic similarity between two GO terms
#' @param GOID1 GO ID 1
#' @param GOID2 GO ID 2
#' @template params-measure
#' @return similarity
#' @seealso [goSim()] [mgoSim()] [geneSim()] [mgeneSim()] [clusterSim()] [mclusterSim()]
#' @export
#' @examples
#' d <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)
#' goSim("GO:0004022", "GO:0005515", semData = d, measure = "Wang")
#' @author Guangchuang Yu <https://yulab-smu.top>
goSim <- function(GOID1, GOID2, semData, measure="Wang") {
    res <- termSim(GOID1, GOID2, semData, method=measure)
    res <- as.numeric(res)
    return(round(res, digits = 3))
}

