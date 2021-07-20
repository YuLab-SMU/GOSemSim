##'Pairwise Semantic Similarities for a List of Gene Clusters
##'
##'Given a list of gene clusters, this function calculates pairwise semantic
##'similarities.
##'
##'
##'@param clusters A list of gene clusters.
##'@param semData GOSemSimDATA object
##'@param measure One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
##'@param drop A set of evidence codes based on which certain annotations are
##'dropped. Use NULL to keep all GO annotations.
##'@param combine One of "max", "avg", "rcmax", "BMA" methods, for combining
##'semantic similarity scores of multiple GO terms associated with protein or
##'multiple proteins assiciated with protein cluster.
##'@return similarity matrix
##'@seealso \code{\link{goSim}} \code{\link{mgoSim}} \code{\link{geneSim}}
##'\code{\link{mgeneSim}} \code{\link{clusterSim}}
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
##'  d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
##'  cluster1 <- c("835", "5261","241")
##'  cluster2 <- c("578","582")
##'  cluster3 <- c("307", "308", "317")
##'  clusters <- list(a=cluster1, b=cluster2, c=cluster3)
##'  mclusterSim(clusters, semData=d, measure="Wang")
##'
mclusterSim <- function(clusters, semData, measure="Wang", drop="IEA", combine="BMA") {
    n <- length(clusters)
    cluster_gos <- list()
    for (i in 1:n) {
        cluster_gos[[i]] <- sapply(clusters[[i]], gene2GO, semData, dropCodes=drop)
    }

    uniqueGO <-  unique(unlist(cluster_gos))
    go_matrix <- mgoSim(uniqueGO, uniqueGO, semData, measure = measure, combine = NULL)

    scores <- matrix(NA, nrow=n, ncol=n)
    rownames(scores) <- names(clusters)
    colnames(scores) <- names(clusters)

    for (i in seq_along(clusters)) {
        gos1 <- unlist(cluster_gos[[i]])
        gos1 <- gos1[!is.na(gos1)]
        for (j in seq_len(i)) {
            gos2 <- unlist(cluster_gos[[j]])
            gos2 <- gos2[!is.na(gos2)]
            if (length(gos1) != 0 && length(gos2) !=0)
                scores[i,j] <- combineScores(go_matrix[gos1, gos2], combine=combine)
                scores[j,i] <- scores[i,j]
        }
    }
    removeRowNA <- apply(!is.na(scores), 1, sum)>0
    removeColNA <- apply(!is.na(scores), 2, sum)>0
    return(scores[removeRowNA, removeColNA, drop=FALSE])
}
