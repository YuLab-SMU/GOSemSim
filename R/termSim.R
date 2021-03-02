##'termSim
##'
##'measuring similarities between two term vectors.
##'
##'provide two term vectors, this function will calculate their similarities.
##'
##'@param t1 term vector
##'@param t2 term vector
##'@param method one of "Wang", "Resnik", "Rel", "Jiang", and "Lin", "TCSS".
##'@param semData GOSemSimDATA object
##'@return score matrix
##'@export
##'@author Guangchuang Yu \url{http://guangchuangyu.github.io}
termSim <- function(t1,
                    t2,
                    semData,
                    method=c("Wang","Resnik","Rel","Jiang","Lin", "TCSS")
                    ) {

    method <- match.arg(method)

    if (all(is.na(t1)) || all(is.na(t2)))
        return(NA)

    t1 <- t1[!is.na(t1)]
    t2 <- t2[!is.na(t2)]

    ## genes in cluster may share go terms,
    ## and should affect similarity of clusters in clusterSim and mclusterSim.
    ##
    ## t1 <- unique(t1)
    ## t2 <- unique(t2)

    if ( method %in% c("Resnik", "Jiang", "Lin", "Rel") ) {
        return(infoContentMethod(t1, t2, method = method, semData))
    } else if ( method == "Wang" ) {
        return(wangMethod(t1, t2, semData@ont))
    } else if ( method == "TCSS" ) {
        return(tcssMethod(t1, t2, semData))
    }
}
