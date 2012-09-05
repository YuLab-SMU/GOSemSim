##' Load Information Content data to DOSEEnv environment
##'
##'
##' @title Load IC data
##' @param organism "human"
##' @param ont "DO"
##' @return NULL
##' @author Guangchuang Yu \url{http://ygc.name}
loadICdata <- function(organism, ont) {
    if(!exists("ICEnv")) .initial()
    fname <- paste("Info_Contents",
                   organism,
                   ont,
                   sep="_")
    if (ont == "DO") {
        tryCatch(utils::data(list=fname,
                             package="DOSE"))
    } else {
        tryCatch(utils::data(list=fname,
                             package="GOSemSim"))
    }
    IC <- get("IC")
    org.ont.IC <- paste(organism,
                        ont,
                        "IC",
                        sep="")
    assign(eval(org.ont.IC),
           IC,
           envir=ICEnv)
    rm (IC)
}


getIC <- function(organism, ont) {
    if(!exists("ICEnv")) {
        .initial()
    }

    org.ont.IC <- paste(organism,
                        ont,
                        "IC",
                        sep="")

    if(!exists(org.ont.IC, envir=ICEnv)) {
        loadICdata(organism, ont)
    }
    IC <- get(org.ont.IC, envir=ICEnv)
    return(IC)
}

getAncestors <- function(ID, ont) {
    Ancestors <- switch(ont,
                        MF = "GOMFANCESTOR",
                        BP = "GOBPANCESTOR",
                        CC = "GOCCANCESTOR",
                        DO = "DOANCESTOR"
                        )
    if (ont == "DO") {
        db <- "DO.db"
        require(db, character.only=TRUE)
    }
    Ancestors <- eval(parse(text=Ancestors))
    ## anc <- get(ID, Ancestors)
    anc <- Ancestors[[ID]]
    return (anc)
}


##' Information Content Based Methods for semantic similarity measuring
##'
##' implemented for methods proposed by Resnik, Jiang, Lin and Schlicker.
##' @title information content based methods
##' @param ID1 Ontology Term
##' @param ID2 Ontology Term
##' @param ont Ontology
##' @param method one of "Resnik", "Jiang", "Lin" and "Rel".
##' @param organism one of supported species
##' @return semantic similarity score
##' @importFrom DO.db DOANCESTOR
##' @author Guangchuang Yu \url{http://ygc.name}
infoContentMethod <- function(ID1,
                               ID2,
                               ont="DO",
                               method,
                               organism="human") {
    IC <- getIC(organism, ont)

    ancestor1 <- getAncestors(ID1, ont)
    ancestor2 <- getAncestors(ID2, ont)

    ## IC is biased
    ## because the IC of a term is dependent of its children but not on its parents.

    sim <- .Call("infoContentMethod_cpp",
                 ID1, ID2,
                 ancestor1, ancestor2,
                 names(IC), IC,
                 method, ont,
                 package="GOSemSim"
                 )
    return (sim)
}
