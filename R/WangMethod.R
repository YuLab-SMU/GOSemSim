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
    if (ont == "DO") {

    } else {
        if (!exists(".GOSemSimEnv")) .initial()
        rel_df <- get("godata", envir=.GOSemSimEnv)        
    }
    sv.a <- getSV(ID1, ont, rel_df)
    sv.b <- getSV(ID2, ont, rel_df)

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

getSV <- function(ID, ont, rel_df, weight=NULL) {
    if (!exists(".SemSimCache")) .initial()
    if( exists(ID, envir=.SemSimCache) ) {
        sv <- get(ID, envir=SemSimCache)
        return(sv)
    }

    if (is.null(weight)) {
        weight <- c(0.8, 0.6, 0.7)
        names(weight) <- c("is_a", "part_of", "other")
    }
    
    rel_df <- rel_df[rel_df$Ontology == ont,]
    rel_df$relationship[!rel_df$relationship %in% c("is_a", "part_of")] <- "other"
    
    id <- ID
    sv <- 1
    names(sv) <- id
    allid <- ID

    idx <- which(rel_df[,1] %in% id)
    while (length(idx) != 0) {
        p <- rel_df[idx,]
        pid <- p$parent
        allid <- c(allid, pid)
        
        sv <- c(sv, w[p$relationship]*sv[id])
        names(sv) <- allid
        idx <- which(rel_df[,1] %in% pid)
    }

    if( ! exists(ID, envir=.SemSimCache) ) {
        assign(ID,
               sv,
               envir=.SemSimCache)
    }

    return(sv)
}

