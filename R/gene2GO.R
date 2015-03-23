##' mapping organism name to annotationDb package name
##'
##'
##' @title getDb
##' @param organism one of supported organism
##' @return annotationDb name
##' @export
##' @author Yu Guangchuang
getDb <- function(organism) {
    annoDb <- switch(organism,
                     anopheles   = "org.Ag.eg.db",
                     arabidopsis = "org.At.tair.db",
                     bovine      = "org.Bt.eg.db",
                     canine      = "org.Cf.eg.db",
                     chicken     = "org.Gg.eg.db",
                     chimp       = "org.Pt.eg.db",
                     coelicolor  = "org.Sco.eg.db", 
                     ecolik12    = "org.EcK12.eg.db",
                     ecsakai     = "org.EcSakai.eg.db",
                     fly         = "org.Dm.eg.db",
                     gondii      = "org.Tgondii.eg.db",
                     human       = "org.Hs.eg.db",
                     malaria     = "org.Pf.plasmo.db",
                     mouse       = "org.Mm.eg.db",
                     pig         = "org.Ss.eg.db",
                     rat         = "org.Rn.eg.db",
                     rhesus      = "org.Mmu.eg.db",
                     worm        = "org.Ce.eg.db",
                     xenopus     = "org.Xl.eg.db",
                     yeast       = "org.Sc.sgd.db",
                     zebrafish   = "org.Dr.eg.db",
                     )
    return(annoDb)
}

loadGOMap_internal <- function(organism){
    annoDb <- getDb(organism)

    ## loading annotation pakcage
    ## require(annoDb, character.only = TRUE)
    requireNamespace(annoDb)
    
    gomap <- switch(organism,
                    anopheles    = "org.Ag.egGO",
                    arabidopsis  = "org.At.tairGO",
                    bovine       = "org.Bt.egGO",
                    canine       = "org.Cf.egGO",
                    chicken      = "org.Gg.egGO",
                    chimp        = "org.Pt.egGO",
                    coelicolor   = "org.Sco.egGO", 
                    ecolik12     = "org.EcK12.egGO",
                    ecsakai      = "org.EcSakai.egGO",
                    fly          = "org.Dm.egGO",
                    gondii       = "org.Tgondii.egGO",
                    human        = "org.Hs.egGO",
                    malaria      = "org.Pf.plasmoGO",
                    mouse        = "org.Mm.egGO",
                    pig          = "org.Ss.egGO",
                    rat          = "org.Rn.egGO",
                    rhesus       = "org.Mmu.egGO",
                    worm         = "org.Ce.egGO",
                    xenopus      = "org.Xl.egGO",
                    yeast        = "org.Sc.sgdGO",
                    zebrafish    = "org.Dr.egGO",
                    )
    gomap <- paste0(annoDb, "::", gomap)
    gomap <- eval(parse(text=gomap))
    assign("gomap", gomap, envir=GOSemSimEnv)
    assign("gomap.flag", organism, envir=GOSemSimEnv)
}

##' loading GOMap to GOSemSimEnv
##'
##'
##' @title loadGOMap
##' @param organism one of supported organisms
##' @return envir
##' @importMethodsFrom AnnotationDbi exists
##' @importMethodsFrom AnnotationDbi get
##' @export
##' @author Yu Guangchuang
loadGOMap <- function(organism) {
    if(!exists("GOSemSimEnv")) .initial()
    if (!exists("gomap", envir=GOSemSimEnv)) {
        loadGOMap_internal(organism)
    } else {
        flag <- get("gomap.flag", envir=GOSemSimEnv)
        if (flag != organism)
            loadGOMap_internal(organism)
    }
}

##' @importMethodsFrom AnnotationDbi get
gene2GO <- function(gene, organism, ont, dropCodes) {
    gene <- as.character(gene)
    loadGOMap(organism)
    gomap <- get("gomap", envir=GOSemSimEnv)
    go <- gomap[[gene]]

    if (all(is.na(go)))
        return (NA)

    ## go.df <- ldply(go, function(i) c(GOID=i$GOID, Evidence=i$Evidence, Ontology=i$Ontology))
    ## go.df <- go.df[ !go.df$Evidence %in% dropCodes, ] ## compatible to work with NA and NULL
    ## goid <- go.df[go.df$Ontology == ont, "GOID"]
    goid <- sapply(go, function(i) i$GOID)
    evidence <- sapply(go, function(i) i$Evidence)
    ontology <- sapply(go, function(i) i$Ontology)

    idx <- ! evidence %in% dropCodes
    goid <- goid[idx] ## drop dropCodes Evidence
    ontology <- ontology[idx]
    goid <- goid[ontology == ont]

    if (length(goid) == 0)
	   return (NA)

    goid <- as.character(unique(goid))
    return (goid)
}
