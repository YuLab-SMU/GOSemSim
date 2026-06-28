#' combining similarity matrix to similarity score
#'
#' Functions for combining similarity matrix to similarity score
#'
#'
#' @param SimScores similarity matrix
#' @param combine combine method
#' @return similarity value
#' @export
#' @author Guangchuang Yu <https://yulab-smu.top>
combineScores <- function(SimScores, combine) {

    if (length(combine) == 0) {  #if not define combine
        return(round(SimScores, digits=3))
    }

    ## if combine was defined...
    if(!sum(!is.na(SimScores))) return (NA)

    if (is.vector(SimScores) || nrow(SimScores)==1 || ncol(SimScores)==1) {
        if (combine == "avg") {
            return(round(mean(SimScores, na.rm=TRUE), digits=3))
        } else {
            return (round(max(SimScores, na.rm=TRUE), digits=3))
        }
    }


    row.na.idx <- rowSums(!is.na(SimScores)) == 0
    if (any(row.na.idx)) {
        SimScores <- SimScores[-which(row.na.idx), ]
    }

    if (! is.null(dim(SimScores)) ) {
        col.na.idx <- colSums(!is.na(SimScores)) == 0
        if (any(col.na.idx)) {
            SimScores <- SimScores[ , -which(col.na.idx)]
        }
    }
    if (is.vector(SimScores) || nrow(SimScores)==1 || ncol(SimScores)==1) {
        if (combine == "avg") {
            return(round(mean(SimScores, na.rm=TRUE), digits=3))
        } else {
            return (round(max(SimScores, na.rm=TRUE), digits=3))
        }
    }

    if (combine        == "avg") {
        result   <- mean(SimScores, na.rm=TRUE)
    } else if (combine == "max") {
        result   <- max(SimScores, na.rm=TRUE)
    } else if (combine == "rcmax") {
        rowMax <- do.call(pmax, c(as.data.frame(SimScores), na.rm = TRUE))
        colMax <- do.call(pmax, c(as.data.frame(t(SimScores)), na.rm = TRUE))
        rowScore <- mean(rowMax)
        colScore <- mean(colMax)
        result   <- max(rowScore, colScore)
    } else if (combine == "rcmax.avg" || combine == "BMA") {
        rowMax <- do.call(pmax, c(as.data.frame(SimScores), na.rm = TRUE))
        colMax <- do.call(pmax, c(as.data.frame(t(SimScores)), na.rm = TRUE))
        result   <- sum(rowMax, colMax) / sum(dim(SimScores))
    }

    return (round(result, digits=3))
}
