#' Calculate pairwise semantic similarities for a given list of genes.
#'
#'
#' @title Pairwise semantic similarity for a list of genes
#' @param genes A list of Entrez gene IDs
#' @template params-measure-combine
#' @param drop Evidence codes to drop; use `NULL` to keep all GO annotations
#' @param verbose Whether to show a progress bar
#' @param BPPARAM optional [BiocParallel::BiocParallelParam-class] object for
#' parallel pairwise similarity calculation. The default `NULL` uses the
#' original serial implementation.
#' @details Parallel calculation is opt-in. With the default `BPPARAM = NULL`,
#' `mgeneSim()` keeps the original serial behavior and can show a progress bar
#' when `verbose = TRUE`. When a `BPPARAM` object is supplied, pairwise
#' similarities are calculated through [BiocParallel::bplapply()] and the
#' progress bar is not shown.
#' @return similarity matrix
#' @seealso [goSim()] [mgoSim()] [geneSim()] [mgeneSim()] [clusterSim()] [mclusterSim()]
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @examples
#' d <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)
#' mgeneSim(c("835", "5261", "241"), semData = d, measure = "Wang")
#' @author Guangchuang Yu <https://yulab-smu.top>
mgeneSim <- function(genes, semData, measure="Wang", drop="IEA", combine="BMA", verbose=TRUE,
                     BPPARAM = NULL) {
    genes <- unique(as.character(genes))
    n <- length(genes)
    scores <- matrix(NA, nrow=n, ncol=n)
    rownames(scores) <- genes
    colnames(scores) <- genes

    gos <- lapply(genes, gene2GO, godata = semData, dropCodes = drop)
    uniqueGO <- uniqueGOFrom(gos)
    go_matrix <- buildGoMatrix(uniqueGO, semData, measure)
    if (!is.null(BPPARAM)) {
        scores <- pairwiseCombineMatrix(genes, gos, go_matrix, combine, BPPARAM = BPPARAM)
        removeRowNA <- apply(!is.na(scores), 1, sum) > 0
        removeColNA <- apply(!is.na(scores), 2, sum) > 0
        return(scores[removeRowNA, removeColNA, drop = FALSE])
    }

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
            scores[i, j] <- subsetCombine(go_matrix, gos[[i]], gos[[j]], combine)
            scores[j, i] <- scores[i, j]
        }
    }
    if (verbose)
        close(pb)
    removeRowNA <- apply(!is.na(scores), 1, sum) > 0
    removeColNA <- apply(!is.na(scores), 2, sum) > 0
    return(scores[removeRowNA, removeColNA, drop = FALSE])
}

