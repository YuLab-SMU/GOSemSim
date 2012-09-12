##'termSim
##'
##'measuring similarities between two term vectors.
##'
##'provide two term vectors, this function will calculate their similarities.
##'
##'@param t1 term vector
##'@param t2 term vector
##'@param method one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##'@param organism only "human" supported
##'@param ont ontology
##'@return score matrix
##' @export
##'@author Guangchuang Yu \url{http://ygc.name}
termSim <- function(t1,
                    t2,
                    method="Wang",
                    organism="human",
                    ont) {

    if (all(is.na(t1)) || all(is.na(t2)))
        return (NA)

    t1 <- t1[!is.na(t1)]
    t2 <- t2[!is.na(t2)]
    t1 <- unique(t1)
    t2 <- unique(t2)

    m <- length(t1)
    n <- length(t2)
    scores <- matrix(NA, nrow=m, ncol=n)
    rownames(scores) <- t1
    colnames(scores) <- t2

    ICmethods <- c("Resnik", "Jiang", "Lin", "Rel")
    ic = method %in% ICmethods

    for( i in 1:m) {
        for (j in 1:n) {
                if (ic) {
                    scores[i,j] <- infoContentMethod(t1[i],
                                                     t2[j],
                                                     ont=ont,
                                                     method=method,
                                                     organism=organism)
                }
                if (method == "Wang") {
                    scores[i,j] <- wangMethod(t1[i],
                                              t2[j],
                                              ont=ont)
                }
        }
    }
    return(scores)
}
