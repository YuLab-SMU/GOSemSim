##'termSim
##'
##'measuring similarities between two term vectors.
##'
##'provide two term vectors, this function will calculate their similarities.
##'
##'@param t1 term vector
##'@param t2 term vector
##'@param method one of "Wang", "Resnik", "Rel", "Jiang", and "Lin".
##'@param organism about 20 species supported, please refer to the vignettes
##'@param ont ontology
##'@return score matrix
##' @export
##'@author Guangchuang Yu \url{http://ygc.name}
termSim <- function(t1,
                    t2,
                    method=c("Wang","Resnik","Rel","Jiang","Lin"),
                    organism="human",
                    ont="BP") {

    if (organism == "worm") {
        organism = "celegans"
        warning("'worm' is deprecated, please use 'celegans' instead...")
    }
    
    method <- match.arg(method)

    if (all(is.na(t1)) || all(is.na(t2)))
        return (NA)

    t1 <- t1[!is.na(t1)]
    t2 <- t2[!is.na(t2)]
    t1 <- unique(t1)
    t2 <- unique(t2)

    if ( method %in% c("Resnik", "Jiang", "Lin", "Rel") ) {
      return ( infoContentMethod( t1, t2, ont=ont, method=method, organism=organism ) )
    }
    else if ( method == "Wang" ) {
      # FIXME vectorize inside wangMethod()
      return ( matrix( mapply( wangMethod,
                                rep( t1, length(t2) ),
                                rep( t2, each=length(t1) ),
                                MoreArgs = list( ont = ont ) ),
                       dimnames = list( t1, t2 ), ncol=length(t2) ) )
    }
}
