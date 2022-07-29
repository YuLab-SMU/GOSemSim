##'Pairwise Semantic Similarity for a List of Genes
##'
##'Given a list of genes, this function calculates pairwise semantic
##'similarities.
##'
##'
##'@param genes A list of entrez gene IDs.
##'@param semData GOSemSimDATA object
##'@param measure One of "Resnik", "Lin", "Rel", "Jiang", "TCSS" and "Wang" methods.
##'@param drop A set of evidence codes based on which certain annotations are
##'dropped. Use NULL to keep all GO annotations.
##'@param combine One of "max", "avg", "rcmax", "BMA" methods, for combining
##'semantic similarity scores of multiple GO terms associated with protein or
##'multiple proteins assiciated with protein cluster.
##' @param verbose show progress bar or not.
##'@return similarity matrix
##'@seealso \code{\link{goSim}} \code{\link{mgoSim}} \code{\link{geneSim}}
##'\code{\link{clusterSim}} \code{\link{mclusterSim}}
##'@references Yu et al. (2010) GOSemSim: an R package for measuring semantic
##'similarity among GO terms and gene products \emph{Bioinformatics} (Oxford,
##'England), 26:7 976--978, April 2010. ISSN 1367-4803
##'\url{http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/976}
##'PMID: 20179076
##'
##'@keywords manip
##' @importFrom utils setTxtProgressBar
##' @importFrom utils txtProgressBar
##'@export
##'@examples
##'
##' d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
##'	mgeneSim(c("835", "5261","241"), semData=d, measure="Wang")
##'
mgeneSim <- function (genes, semData, measure="Wang", drop="IEA", combine="BMA", verbose=TRUE) {
    genes <- unique(as.character(genes))
    n <- length(genes)
    scores <- matrix(NA, nrow=n, ncol=n)
    rownames(scores) <- genes
    colnames(scores) <- genes

    gos <- lapply(genes, gene2GO, godata=semData, dropCodes=drop)
    uniqueGO <-  unique(unlist(gos))
    go_matrix <- mgoSim(uniqueGO, uniqueGO, semData, measure = measure, combine = NULL)
    if (verbose) {
      cnt <- 1
      pb <- txtProgressBar(min=0, max=sum(1:n), style=3)
    }
    for (i in seq_along(genes)) {
        for (j in seq_len(i)){
            if (verbose) {
                setTxtProgressBar(pb, cnt)
                cnt <- cnt + 1
            }
            scores[i,j] <- combineScores(go_matrix[gos[[i]], gos[[j]]], combine = combine)
            scores[j,i] <- scores[i,j]
        }
    }
    if (verbose)
        close(pb)
    removeRowNA <- apply(!is.na(scores), 1, sum)>0
    removeColNA <- apply(!is.na(scores), 2, sum)>0
    return(scores[removeRowNA, removeColNA, drop=FALSE])
}

