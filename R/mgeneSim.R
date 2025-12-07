#' Pairwise semantic similarity for a list of genes
#'
#' Calculate pairwise semantic similarities for a given list of genes.
#'
#' @param genes A list of Entrez gene IDs
#' @template params-measure-combine
#' @param drop Evidence codes to drop; use `NULL` to keep all GO annotations
#' @param verbose Whether to show a progress bar
#' @return similarity matrix
#' @template seealso-similarity
#' @template references
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @export
#' @examples
#' d <- godata('org.Hs.eg.db', ont = "MF", computeIC = FALSE)
#' mgeneSim(c("835", "5261", "241"), semData = d, measure = "Wang")

mgeneSim <- function (genes, semData, measure = "Wang", drop = "IEA", combine = "BMA", verbose = TRUE) {
    genes <- unique(as.character(genes))
    n <- length(genes)
    scores <- matrix(NA, nrow=n, ncol=n)
    rownames(scores) <- genes
    colnames(scores) <- genes

    gos <- lapply(genes, gene2GO, godata = semData, dropCodes = drop)
    uniqueGO <- uniqueGOFrom(gos)
    go_matrix <- buildGoMatrix(uniqueGO, semData, measure)
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

