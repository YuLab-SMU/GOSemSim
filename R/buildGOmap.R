##' Addding indirect GO annotation
##'
##' provided by a data.frame of GO TERM (column 1), GENE (column 2) and ONTOLOGY (optional) that
##' describes GO direct annotation, 
##' this function will add indirect GO annotation of genes.
##' @title buildGOmap
##' @param TERM2GENE data.frame with two or three columns of GO TERM, GENE and ONTOLOGY (optional)
##' @return data.frame, GO annotation with direct and indirect annotation
##' @importMethodsFrom AnnotationDbi as.list
##' @importFrom GO.db GOMFANCESTOR
##' @importFrom GO.db GOBPANCESTOR
##' @importFrom GO.db GOCCANCESTOR
##' @export
##' @author Yu Guangchuang
buildGOmap <- function(TERM2GENE) {
    mfanc <- as.list(GOMFANCESTOR)
    ccanc <- as.list(GOCCANCESTOR)
    bpanc <- as.list(GOBPANCESTOR)

    if (!'ONTOLOGY' %in% names(TERM2GENE)) {
        anc <- c(mfanc, ccanc, bpanc)
        res <- buildGOmap_internal(TERM2GENE, anc)
        return(res)
    }

    anc <- list(MF=mfanc, CC=ccanc, BP=bpanc)
    y <- split(TERM2GENE, TERM2GENE$ONTOLOGY)

    res <- lapply(names(y), function(i) {
        d <- buildGOmap_internal(y[[i]], anc[[i]])
        d$ONTOLOGY <- i
        return(d)
    }) |> do.call('rbind', args = _)

    return(res)
}

##' @importFrom stats setNames
##' @importFrom yulab.utils ls2df
buildGOmap_internal <- function(TERM2GENE, anc) {
    res <- setNames(anc[TERM2GENE[,1]], TERM2GENE[,2]) |> 
        ls2df() |>
        unique()

    res <- setNames(res[, c(2,1)], names(TERM2GENE)[1:2])
    res <- res[res[,1] != "all", ]
    res <- rbind(TERM2GENE[,1:2], res)
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
