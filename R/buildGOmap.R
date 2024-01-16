##' Addding indirect GO annotation
##'
##' provided by a data.frame of GENE (column 1), GO (column 2) and ONTOLOGY (optional) that
##' describes GO direct annotation, 
##' this function will add indirect GO annotation of genes.
##' @title buildGOmap
##' @param x data.frame with two or three columns of GENE, GO and ONTOLOGY (optional)
##' @return data.frame, GO annotation with direct and indirect annotation
##' @importMethodsFrom AnnotationDbi as.list
##' @importFrom GO.db GOMFANCESTOR
##' @importFrom GO.db GOBPANCESTOR
##' @importFrom GO.db GOCCANCESTOR
##' @export
##' @author Yu Guangchuang
buildGOmap <- function(x) {
    mfanc <- as.list(GOMFANCESTOR)
    ccanc <- as.list(GOCCANCESTOR)
    bpanc <- as.list(GOBPANCESTOR)

    if (!'ONTOLOGY' %in% names(x)) {
        anc <- c(mfanc, ccanc, bpanc)
        res <- buildGOmap_internal(x, anc)
        return(res)
    }

    anc <- list(MF=mfanc, CC=ccanc, BP=bpanc)
    y <- split(x, x$ONTOLOGY)

    res <- lapply(names(y), function(i) {
        d <- buildGOmap_internal(y[[i]], anc[[i]])
        d$ONTOLOGY <- i
        return(d)
    }) |> do.call('rbind', args = _)

    return(res)
}

##' @importFrom stats setNames
##' @importFrom yulab.utils ls2df
buildGOmap_internal <- function(y, anc) {
    res <- setNames(anc[y$GO], y[,1]) |> 
        ls2df() |>
        unique()

    names(res) <- c(names(y)[1], "GO")
    res <- res[res$GO != "all", ]

    res <- rbind(y[, names(res)], res)
    return(res)
}

# old and slow version
##' @importMethodsFrom AnnotationDbi mget
##' @importFrom utils stack
buildGOmap2 <- function(gomap) {

    ## remove empty GO annotation
    gomap <- gomap[gomap[,1] != "", ]
    
    Gene2GO <- split(as.character(gomap[,1]), as.character(gomap[,2]))

    Gene2ALLGO <- lapply(Gene2GO,
                         function(i) {
                             mfans <- unlist(mget(i, GOMFANCESTOR, ifnotfound=NA))
                             bpans <- unlist(mget(i, GOBPANCESTOR, ifnotfound=NA))
                             ccans <- unlist(mget(i, GOCCANCESTOR, ifnotfound=NA))
                             ans <- c(mfans, bpans, ccans)
                             ans <- ans[ !is.na(ans) ]
                             ans <- c(i, ans)
                             ans <- unique(ans)
                             ans <- ans[ans != "all"]
                             return(ans)
                         })

    ## AMF <- as.list(GOMFANCESTOR)
    ## ACC <- as.list(GOCCANCESTOR)
    ## ABP <- as.list(GOBPANCESTOR)

    ## Gene2ALLGO <- lapply(Gene2GO, function(i) {
    ##     mfans <- AMF[i]
    ##     bpans <- ABP[i]
    ##     ccans <- ACC[i]
    ##     ans <- unlist(c(mfans,  bpans,  ccans))
    ##     ans <- ans[ !is.na(ans) ]
    ##     ans <- c(i, ans)
    ##     ans <- unique(ans)
    ##     ans <- ans[ans != "all"]
    ##     return(ans)
    ## })

    go2gene <- stack(Gene2ALLGO)
    colnames(go2gene) <- c("GO", "Gene")
    
    return(go2gene)
}
