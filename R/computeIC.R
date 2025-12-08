computeIC <- function(goAnno, ont) {
    ## goAnno, see godata function
    if (!exists(".GOSemSimEnv")) .initial()
    .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
    godata <- get("gotbl", envir=.GOSemSimEnv)        
    
    goids <- unique(godata[godata$Ontology == ont, "go_id"])
    ## all GO terms appearing in an given ontology ###########
    goterms=goAnno$GO
    gocount <- table(goterms)
    ## goid of specific organism and selected category.
    goname  <- names(gocount) 

    ## ensure goterms not appearing in the specific annotation have 0 frequency..
    go.diff        <- setdiff(goids, goname)
    m              <- double(length(go.diff))
    names(m)       <- go.diff
    gocount        <- as.vector(gocount)
    names(gocount) <- goname
    gocount        <- c(gocount, m)

    offspring_idx <- getOffspringIdx(ont, goids)
    gc <- gocount[goids]
    desc_sum <- vapply(goids, function(id) {
        ids <- offspring_idx[[id]]
        if (length(ids) == 0) 0 else sum(gc[ids], na.rm = TRUE)
    }, numeric(1))
    cnt <- gc + desc_sum
    names(cnt) <- goids
    
    ## the probabilities of occurrence of GO terms in a specific corpus.
    p <- cnt/sum(gocount)
    ## IC of GO terms was quantified as the negative log likelihood.
    IC <- -log(p)
    return(IC)
}


