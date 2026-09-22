computeIC <- function(goAnno, ont) {
    ## goAnno, see godata function

    ## The set of terms to score has to come from the *installed* GO.db, not from
    ## the `gotbl` table shipped in data/. That table was last rebuilt in 2021,
    ## and 524 BP terms added to GO.db since then are missing from it. They were
    ## therefore given no IC at all, and goSim()/geneSim() returned NA for them
    ## even though their descendants carry annotations and their IC is perfectly
    ## computable -- 167 of those terms are annotated in org.Hs.eg.db (issue
    ## #33). The missing terms also biased every other IC slightly: they are
    ## counted in the `sum(gocount)` denominator below but can never contribute
    ## to an ancestor's descendant sum.
    ##
    ## names(getAncestors(ont)) is exactly the set of terms of `ont` in the
    ## installed GO.db, root included, and onto_relation() caches it. Terms
    ## annotated in this organism are added on top, so that a term the local
    ## annotation package still uses keeps its IC even if GO.db has retired it.
    goids <- union(names(getAncestors(ont)), goAnno$GO)

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


