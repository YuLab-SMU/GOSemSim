##' Method Wang for semantic similarity measuring
##'
##'
##' @title wangMethod
##' @param ID1 Ontology Term
##' @param ID2 Ontology Term
##' @param ont Ontology
##' @return semantic similarity score
##' @author Guangchuang Yu \url{http://ygc.name}
wangMethod <- function(ID1, ID2, ont="BP") {
    if (ID1 == ID2)
        return (sim=1)

    sv.a <- getSV(ID1, ont)
    sv.b <- getSV(ID2, ont)

    if(all(is.na(sv.a)) || all(is.na(sv.b)))
        return (NA)

    idx         <- intersect(names(sv.a), names(sv.b))
    inter.sva   <- sv.a[idx]
    inter.svb   <- sv.b[idx]
    if (is.null(inter.sva) ||
        is.null(inter.svb) ||
        length(inter.sva) == 0 ||
        length(inter.svb) ==0) {
        return (NA)
    } 
    
    sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
    return(sim)
}

##' @importMethodsFrom AnnotationDbi exists
##' @importMethodsFrom AnnotationDbi get
getSV <- function(ID, ont) {
    if (!exists("SemSimCache")) .initial()
    if( exists(ID, envir=SemSimCache) ) {
        sv <- get(ID, envir=SemSimCache)
    } else {
        Parents     <- .getParents(ont)
        if ( !exists(ID, Parents))
            return(NA)
        
        if (ont == "DO") {
            topNode <- "DOID:4"
        } else {
            topNode <- "all"
        }

        if (ID == topNode) {
            sv <- 1
            names(sv) <- topNode
            return (sv)
        }

        sv.name <- c(ID, getAncestors(ont)[[ID]])
        sv <- rep(NA, length(sv.name))
        names(sv) <- sv.name
        sv[ID] <- 1
        if(ont != "DO")
            sv[topNode] <- 0

        w <- c(0.8, 0.6, 0.7)
        names(w) <- c("is_a", "part_of", "other")


        pID <- ID
        while(any(is.na(sv))) {
            pp <- c()
            for (i in seq_along(pID)) {
                if (pID[i] != topNode) {
                    j <- get(pID[i], Parents)
                    idx <- which(is.na(sv[j]))
                    if (length(idx)) {
                        j <- j[idx]
                        if (is.null(names(j)) || all(is.na(names(j))) ) {
                            names(j) = "other"
                        } else {
                            names(j)[!names(j) %in% names(w)] = "other"
                        }
                        sv[j] = sv[pID[i]] * w[names(j)]
                    }
                }
                pp <- c(pp, j)
            }
            pID <- unique(pp)
            if (all(pID == topNode) || length(pID) == 0) {
                sv <- sv[!is.na(sv)]
                break
            }
        }
    }

    if( ! exists(ID, envir=SemSimCache) ) {
        assign(ID,
               sv,
               envir=SemSimCache)
    }
    return(sv)
}

