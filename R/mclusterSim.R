#' Calculate pairwise semantic similarities for a list of gene clusters.
#'
#'
#' @title Pairwise semantic similarities for a list of gene clusters
#' @param clusters A list of gene clusters
#' @template params-measure-combine
#' @param drop Evidence codes to drop; use `NULL` to keep all GO annotations
#' @param BPPARAM optional [BiocParallel::BiocParallelParam-class] object for
#' parallel pairwise similarity calculation. The default `NULL` uses serial
#' calculation.
#' @return similarity matrix
#' @seealso [goSim()] [mgoSim()] [geneSim()] [mgeneSim()] [clusterSim()] [mclusterSim()]
#' @export
#' @examples
#' d <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)
#' cluster1 <- c("835", "5261", "241")
#' cluster2 <- c("578", "582")
#' cluster3 <- c("307", "308", "317")
#' clusters <- list(a = cluster1, b = cluster2, c = cluster3)
#' mclusterSim(clusters, semData = d, measure = "Wang")
#' @author Guangchuang Yu <https://yulab-smu.top>
mclusterSim <- function(clusters, semData, measure="Wang", drop="IEA", combine="BMA",
                        BPPARAM = NULL) {
    n <- length(clusters)
    cluster_gos <- list()
    for (i in 1:n) {
        cluster_gos[[i]] <- sapply(clusters[[i]], gene2GO, semData, dropCodes = drop)
    }

    uniqueGO <- uniqueGOFrom(cluster_gos)
    go_matrix <- buildGoMatrix(uniqueGO, semData, measure)

    cluster_gos <- lapply(cluster_gos, function(gos) {
        gos <- unlist(gos)
        gos[!is.na(gos)]
    })

    labels <- names(clusters)
    if (is.null(labels) || all(is.na(labels)) || all(labels == "")) {
        labels <- NULL
    }
    scores <- pairwiseCombineMatrix(labels, cluster_gos, go_matrix, combine, BPPARAM = BPPARAM)

    removeRowNA <- apply(!is.na(scores), 1, sum) > 0
    removeColNA <- apply(!is.na(scores), 2, sum) > 0
    return(scores[removeRowNA, removeColNA, drop = FALSE])
}
