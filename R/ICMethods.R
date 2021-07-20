##' Information Content Based Methods for semantic similarity measuring
##'
##' implemented for methods proposed by Resnik, Jiang, Lin and Schlicker.
##' @title information content based methods
##' @param ID1 Ontology Term
##' @param ID2 Ontology Term
##' @param method one of "Resnik", "Jiang", "Lin" and "Rel", "TCSS".
##' @param godata GOSemSimDATA object
##' @return semantic similarity score
##' @useDynLib GOSemSim
##' @author Guangchuang Yu \url{https://guangchuangyu.github.io}
infoContentMethod <- function(ID1,
                              ID2,
                              method,
                              godata) {
    ## IC is biased
    ## because the IC of a term is dependent of its children but not on its parents.
    ont <- godata@ont
    IC <- godata@IC

    if (length(IC) == 0) {
        stop("IC data not found, please re-generate your `semData` with `computeIC=TRUE`...")
    }

    if (ont %in% c("MF", "BP", "CC", "DO")) {
        ## .anc <- tryCatch(getAncestors(ont)[union(ID1,ID2)], error=function(e) NULL)
        ## if (is.null(.anc)) {
        ##     ## https://support.bioconductor.org/p/105822/
        ##     return(NA)
        ## }
        ## .anc <- AnnotationDbi::as.list(.anc)

        ## if some IDs are not valid, the above code will leading to return NA directly.
        .anc <- AnnotationDbi::as.list(getAncestors(ont))
        allid <- union(ID1, ID2)
        .anc <- .anc[allid]
        .anc <- .anc[!vapply(.anc, is.empty, logical(1))]

        ## invalid_ids <- c(ID1[!ID1 %in% names(.anc)], ID2[!ID2 %in% names(.anc)])
        ## if (length(invalid_ids) > 0) {
        ##     message("The following IDs are not valid and will be removed:", paste(invalid_ids, collapse=","))
        ##     # allid <- allid[!allid %in% invalid_ids]
        ## }
    } else {
        mesh_getAnc <- eval(parse(text="meshes:::getAncestors"))
        .anc <- lapply(union(ID1, ID2), mesh_getAnc)
        names(.anc) <- union(ID1, ID2)
    }
    return ( infoContentMethod_cpp( ID1, ID2,
                 .anc, IC,
                 method, ont ) )
}

is.empty <- function(x) {
    if (is.null(x)) return(TRUE)
    if (all(is.na(x))) return(TRUE)
    return(FALSE)
}


## infoContentMethod <- function(ID1,
##                               ID2,
##                               ont="DO",
##                               method,
##                               organism="human") {
##     IC <- getIC(organism, ont)

##     ## more specific term, larger IC value.
##     ## Normalized, all divide the most informative IC.
##     ## all IC values range from 0(root node) to 1(most specific node)
##     mic <- max(IC[IC!=Inf])

##     if (ont == "DO") {
##         topNode <- "DOID:4"
##     } else {
##         topNode <- "all"
##     }

##     IC[topNode] = 0

##     ic1 <- IC[ID1]/mic
##     ic2 <- IC[ID2]/mic

##     if (ic1 == 0 || ic2 == 0)
##         return (NA)

##     ancestor1 <- getAncestors(ont)[[ID1]]
##     ancestor2 <- getAncestors(ont)[[ID2]]
##     if (ID1 == ID2) {
##         commonAncestor <- ID1
##     } else if (ID1 %in% ancestor2) {
##         commonAncestor <- ID1
##     } else if (ID2 %in% ancestor1) {
##         commonAncestor <- ID2
##     } else {
##         commonAncestor <- intersect(ancestor1, ancestor2)
##     }
##     if (length(commonAncestor) == 0) return (NA)

##     ##Information Content of the most informative common ancestor (MICA)
##     mica <- max(IC[commonAncestor])/mic

##     ## IC is biased
##     ## because the IC of a term is dependent of its children but not on its parents.
##     sim <- switch(method,
##                   Resnik = mica, ## Resnik does not consider how distant the terms are from their common ancestor.
##                   ## Lin and Jiang take that distance into account.
##                   Lin = 2*mica/(ic1+ic2),
##                   Jiang = 1 - min(1, -2*mica + ic1 + ic2),
##                   Rel = 2*mica/(ic1+ic2)*(1-exp(-mica*mic))  ## mica*mic equals to the original IC value. and exp(-mica*mic) equals to the probability of the term's occurence.
##                   )
##     return (sim)
## }
