#' combining similarity matrix to similarity score
#'
#' Functions for combining similarity matrix to similarity score
#'
#'
#' @param SimScores similarity matrix
#' @param combine combine method, one of "avg", "max", "rcmax", "rcmax.avg"
#'   or "BMA"
#' @return similarity value
#' @export
#' @author Guangchuang Yu <https://yulab-smu.top>
combineScores <- function(SimScores, combine) {

    if (length(combine) == 0) {  #if not define combine
        return(round(SimScores, digits=3))
    }

    ## if combine was defined...
    if(!sum(!is.na(SimScores))) return (NA)

    ## a bare vector is treated as a single row, so that one set of row and
    ## column rules applies to every input shape
    if (is.null(dim(SimScores))) {
        SimScores <- matrix(SimScores, nrow = 1)
    }

    ## drop the rows and columns that carry no score at all
    row.na.idx <- rowSums(!is.na(SimScores)) == 0
    if (any(row.na.idx)) {
        SimScores <- SimScores[-which(row.na.idx), , drop = FALSE]
    }
    col.na.idx <- colSums(!is.na(SimScores)) == 0
    if (any(col.na.idx)) {
        SimScores <- SimScores[, -which(col.na.idx), drop = FALSE]
    }

    if (combine        == "avg") {
        result   <- mean(SimScores, na.rm=TRUE)
    } else if (combine == "max") {
        result   <- max(SimScores, na.rm=TRUE)
    } else if (nrow(SimScores) == 1 || ncol(SimScores) == 1) {
        ## With a single row or column the row maxima average to the max and
        ## the column maxima average to the mean, so 'rcmax' is still the max
        ## but 'BMA' / 'rcmax.avg' is not -- it reduces to (max + sum) / (n + 1)
        ## over the n cells in that one row or column. Returning the max here
        ## made a set holding one term score 1 against any set containing it.
        if (combine == "rcmax.avg" || combine == "BMA") {
            result <- (max(SimScores, na.rm=TRUE) + sum(SimScores, na.rm=TRUE)) /
                (length(SimScores) + 1)
        } else {  # "rcmax"
            result <- max(SimScores, na.rm=TRUE)
        }
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
