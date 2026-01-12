#' Given two genes, calculate their semantic similarity and return the score with corresponding GO terms.
#'
#'
#' @title Semantic similarity between two genes
#' @param gene1 Entrez gene ID
#' @param gene2 Another Entrez gene ID
#' @template params-measure-combine
#' @param drop Evidence codes to drop; use `NULL` to keep all GO annotations
#' @return A list containing similarity value and corresponding GO terms
#' @seealso [goSim()] [mgoSim()] [mgeneSim()] [clusterSim()] [mclusterSim()]
#' @export
#' @examples
#' d <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)
#' geneSim("241", "251", semData = d, measure = "Wang")
#' @author Guangchuang Yu <https://yulab-smu.top>
geneSim <- function(gene1, gene2, semData, measure="Wang", drop="IEA", combine="BMA") {
    go1 <- gene2GO(gene1, semData, dropCodes = drop)
    go2 <- gene2GO(gene2, semData, dropCodes = drop)
    if (length(go1) == 0 || length(go2) == 0)
        return(NA)
    res <- mgoSim(go1, go2, semData = semData, measure = measure, combine = combine)
    return(list(geneSim = res, GO1 = go1, GO2 = go2))
}

