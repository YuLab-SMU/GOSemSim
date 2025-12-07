#' Semantic similarity between two gene clusters
#'
#' Given two gene clusters, calculate their semantic similarity.
#'
#' @param cluster1 A set of gene IDs
#' @param cluster2 Another set of gene IDs
#' @template params-measure-combine
#' @param drop Evidence codes to drop; use `NULL` to keep all GO annotations
#' @return similarity
#' @template seealso-similarity
#' @template references
#' @export
#' @examples
#' d <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)
#' cluster1 <- c("835", "5261", "241", "994")
#' cluster2 <- c("307", "308", "317", "321", "506", "540", "378", "388", "396")
#' clusterSim(cluster1, cluster2, semData = d, measure = "Wang")

clusterSim <- function(cluster1, cluster2, semData, measure = "Wang", drop = "IEA", combine = "BMA"){
    cgo1 <- sapply(cluster1, gene2GO, semData, dropCodes = drop)
    cgo2 <- sapply(cluster2, gene2GO, semData, dropCodes = drop)
    cgo1 <- unlist(cgo1)
    cgo2 <- unlist(cgo2)
    res <- mgoSim(cgo1, cgo2, semData, measure = measure, combine = combine)
    return(res)
}
