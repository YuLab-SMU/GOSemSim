computeIC <- function(goAnno, ont) {
    ## goAnno, see godata function
    get_gosemsim_env()
    gotbl_df <- yulab.utils::get_cache_element(".GOSemSimEnv", "gotbl")
    if (is.null(gotbl_df) || !is.data.frame(gotbl_df)) {
        utils::data("gotbl", package = "GOSemSim")
        gotbl_df <- get("gotbl")
        yulab.utils::update_cache_item(".GOSemSimEnv", list(gotbl = gotbl_df))
    }
    
    goids <- unique(gotbl_df[gotbl_df$Ontology == ont, "go_id"])
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


