wangMethod <- function(t1, t2, ont) {
    matrix( mapply( wangMethod_internal,
                   rep( t1, length(t2) ),
                   rep( t2, each=length(t1) ),
                   MoreArgs = list( ont = ont ) ),
           dimnames = list( t1, t2 ), ncol=length(t2) ) 
}

##' Method Wang for semantic similarity measuring
##'
##'
##' @title wangMethod
##' @param ID1 Ontology Term
##' @param ID2 Ontology Term
##' @param ont Ontology
##' @return semantic similarity score
##' @author Guangchuang Yu \url{http://ygc.name}
wangMethod_internal <- function(ID1, ID2, ont="BP") {
    if (ID1 == ID2)
        return (sim=1)
    if (ont == "DO") {
        .DOSEEnv <- get(".DOSEEnv", envir=.GlobalEnv)
        rel_df <- get("dotbl", envir=.DOSEEnv)
    } else if (ont %in% c("BP", "CC", "MF")) {
        if (!exists(".GOSemSimEnv")) .initial()
        .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
        rel_df <- get("gotbl", envir=.GOSemSimEnv)
    } else {
        .meshesEnv <- get(".meshesEnv", envir=.GlobalEnv)
        rel_df <- get("meshtbl", envir=.meshesEnv)
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
    .SemSimCache <- get(".SemSimCache", envir=.GlobalEnv)
    
    if( exists(ID, envir=.SemSimCache) ) {
        sv <- get(ID, envir=.SemSimCache)
        return(sv)
    }

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
    
    if (is.null(weight)) {
        weight <- c(0.8, 0.6, 0.7)
        names(weight) <- c("is_a", "part_of", "other")
    }

    rel_df <- rel_df[rel_df$Ontology == ont,]
    if (! 'relationship' %in% colnames(rel_df))
        rel_df$relationship <- "other"
    
    rel_df$relationship[!rel_df$relationship %in% c("is_a", "part_of")] <- "other"
    

    sv <- 1
    names(sv) <- ID
    allid <- ID

    idx <- which(rel_df[,1] %in% ID)
    while (length(idx) != 0) {
        p <- rel_df[idx,]
        pid <- p$parent
        allid <- c(allid, pid)
        
        sv <- c(sv, weight[p$relationship]*sv[p[,1]])
        names(sv) <- allid
        idx <- which(rel_df[,1] %in% pid)
    }

    sv <- sv[!is.na(names(sv))]
    sv <- sv[!duplicated(names(sv))]

    if(ont != "DO")
        sv[topNode] <- 0

    if( ! exists(ID, envir=.SemSimCache) ) {
        assign(ID,
               sv,
               envir=.SemSimCache)
    }
    
    return(sv)
}

