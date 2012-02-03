.initial <- function() {
    assign("GOSemSimEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("ICEnv", new.env(), .GlobalEnv)
    assign("SupportedSpecies", c("anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast", "zebrafish"), envir=GOSemSimEnv)
    ## remove support of  "coelicolor",
}

.getParents <- function(ont) {
    Parents <- switch(ont,
                      MF = GOMFPARENTS,
                      BP = GOBPPARENTS,
                      CC = GOCCPARENTS
                      )
    return(Parents)
}

                                        #.getParents <- function(ont="MF") {
                                        #	if(!exists("GOSemSimEnv")) .initial()
                                        #	wh_Parents <- switch(ont,
                                        #		MF = "MFParents",
                                        #		BP = "BPParents",
                                        #		CC = "CCParents"
                                        #	)
                                        #	Parents <- switch(ont,
                                        #		MF = AnnotationDbi::as.list(GOMFPARENTS) ,
                                        #		BP = AnnotationDbi::as.list(GOBPPARENTS) ,
                                        #		CC = AnnotationDbi::as.list(GOCCPARENTS)
                                        #	)
                                        #	assign(eval(wh_Parents), Parents, envir=GOSemSimEnv)
                                        #}

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

.getAncestors <- function(ont) {
    Ancestors <- switch(ont,
                        MF = GOMFANCESTOR,
                        BP = GOBPANCESTOR,
                        CC = GOCCANCESTOR
                        )
    return(Ancestors)
}
                                        #.getAncestors <- function(ont="MF") {
                                        #	if(!exists("GOSemSimEnv")) .initial()
                                        #	wh_Ancestors <- switch(ont,
                                        #		MF = "MFAncestors",
                                        #		BP = "BPAncestors",
                                        #		CC = "CCAncestors"
                                        #	)
                                        #	Ancestors <- switch(ont,
                                        #		MF = AnnotationDbi::as.list(GOMFANCESTOR) ,
                                        #		BP = AnnotationDbi::as.list(GOBPANCESTOR) ,
                                        #		CC = AnnotationDbi::as.list(GOCCANCESTOR)
                                        #	)
                                        #	assign(eval(wh_Ancestors), Ancestors, envir=GOSemSimEnv)
                                        #}


rebuildAllICdata <- function() {
    if(!exists("GOSemSimEnv")) .initial()
    ont <- c("BP","CC", "MF")
    species <- get("SupportedSpecies",envir=GOSemSimEnv)
    cat("------------------------------------\n")
    cat("calulating Information Content...\nSpecies:\t\tOntology\n")
    params <- new("Params")
    for (i in species) {
                                        #loadAnnoPkg(params) ##load annotation pkg.
        setOrganism(params) <- i
        for (j in ont) {
            setOntology(params) <- j
            cat(i)
            cat("\t\t\t")
            cat(j)
            cat("\n")
            file=paste(paste("Info_Contents", i, j, sep="_"), ".rda", sep="")
            if ( !file.exists(file) ) {
                computeIC(params)
            }
        }
    }
    cat("------------------------------------\n")
    print("done...")
}


gene2GO <-  function(gene, params) {
    gene <- as.character(gene)
    if(!exists("GOSemSimEnv")) .initial()
    if (!exists("gomap", envir=GOSemSimEnv)) {
        loadGOMap(params)
    }
    gomap <- get("gomap", envir=GOSemSimEnv)

    qGO	<- gomap[[gene]]

    if (is.null(qGO)) {
        return (NA)
    }
    if (sum(!is.na(qGO)) == 0) {
    	return (NA)
    }

    qGO	<- qGO[qGO == params@ontology]
    if (length(qGO) == 0) {
        return (NA)
    }
    return( unique(names(qGO)) )
}
