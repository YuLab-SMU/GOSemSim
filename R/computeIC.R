##' @importFrom GO.db GOMFOFFSPRING
##' @importFrom GO.db GOBPOFFSPRING
##' @importFrom GO.db GOCCOFFSPRING
.getOffsprings <- function(ont="MF") {
    if(!exists("GOSemSimEnv")) .initial()
    wh_Offsprings <- switch(ont,
                            MF = "MFOffsprings",
                            BP = "BPOffsprings",
                            CC = "CCOffsprings"
                            )
    Offsprings <- switch(ont,
                         MF = AnnotationDbi::as.list(GOMFOFFSPRING) ,
                         BP = AnnotationDbi::as.list(GOBPOFFSPRING) ,
                         CC = AnnotationDbi::as.list(GOCCOFFSPRING)
                         )
    assign(eval(wh_Offsprings), Offsprings, envir=GOSemSimEnv)
}

##' @importFrom GO.db GOTERM
##' @importMethodsFrom AnnotationDbi get
##' @importMethodsFrom AnnotationDbi toTable
##' @importMethodsFrom AnnotationDbi mappedkeys
computeIC <- function(organism, ont) {
    loadGOMap(organism)
    gomap   <- get("gomap", envir=GOSemSimEnv)

    mapped_genes <- mappedkeys(gomap)
    gomap <- AnnotationDbi::as.list(gomap[mapped_genes])

    gomap <- sapply(gomap, function(x) sapply(x, function(y) y$Ontology))

    ## all GO terms appearing in an given ontology ###########
    goterms <- unlist(sapply(gomap, function(x) names(x[x == ont])), use.names=FALSE)

    ## require(GO.db)
    if ( !exists("ALLGOID", envir=GOSemSimEnv) ) {
        assign("ALLGOID", toTable(GOTERM), envir=GOSemSimEnv )
    }
    goids   <- get("ALLGOID", envir=GOSemSimEnv)
    ##goids <- toTable(GOTERM)

    ## all go terms which belong to the corresponding ontology..
    goids   <- unique(goids[goids[,"Ontology"] == ont, "go_id"])
    gocount <- table(goterms)
    goname  <- names(gocount) #goid of specific organism and selected category.

    ## ensure goterms not appearing in the specific annotation have 0 frequency..
    go.diff        <- setdiff(goids, goname)
    m              <- double(length(go.diff))
    names(m)       <- go.diff
    gocount        <- as.vector(gocount)
    names(gocount) <- goname
    gocount        <- c(gocount, m)

    Offsprings.name <- switch(ont,
                              MF = "MFOffsprings",
                              BP = "BPOffsprings",
                              CC = "CCOffsprings"
                              )
    if (!exists(Offsprings.name, envir=GOSemSimEnv)) {
        .getOffsprings(ont)
    }
    Offsprings <- get(Offsprings.name, envir=GOSemSimEnv)
    cnt        <- sapply(goids,function(x){ n=gocount[ Offsprings[[x]] ]; gocount[x]+sum(n[!is.na(n)])})
    names(cnt) <- goids
    ## the probabilities of occurrence of GO terms in a specific corpus.
    p          <- cnt/sum(gocount)
    ## IC of GO terms was quantified as the negative log likelihood.
    IC         <- -log(p)

    save(IC, file=paste(paste("Info_Contents", organism, ont, sep="_"), ".rda", sep=""), compress="xz")
}




rebuildAllICdata <- function() {
    if(!exists("GOSemSimEnv")) .initial()
    ont     <- c("BP","CC", "MF")
    species <- get("SupportedSpecies",envir=GOSemSimEnv)
    cat("------------------------------------\n")
    cat("calulating Information Content...\nSpecies:\t\tOntology\n")
    for (i in species) {
        for (j in ont) {
            cat(i)
            cat("\t\t\t")
            cat(j)
            cat("\n")
            file=paste(paste("Info_Contents", i, j, sep="_"), ".rda", sep="")
            if ( !file.exists(file) ) {
                computeIC(organism=i, ont=j)
            }
        }
    }
    cat("------------------------------------\n")
    print("done...")
}


install_dependent_db <- function() {
    ## require("BiocInstaller", character.only=TRUE)
    requireNamespace("BiocInstaller")
    biocLite <- eval(parse(text="BiocInstaller::biocLite"))

    species <- get("SupportedSpecies",envir=GOSemSimEnv)
    for (organism in species) {
        annoDb <- getDb(organism)
        ## if (!require(annoDb, character.only=TRUE)) {
        if (! requireNamespace(annoDb)) {
            biocLite(annoDb)
        }
    }
}
