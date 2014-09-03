.initial <- function() {
    assign("GOSemSimEnv", new.env(),.GlobalEnv)
    assign("SemSimCache", new.env(), .GlobalEnv)
    assign("ICEnv", new.env(), .GlobalEnv)
    assign("SupportedSpecies", c("anopheles",
                                 "arabidopsis",
                                 "bovine",
                                 "canine",
                                 "chicken",
                                 "chimp",
                                 "ecolik12",
                                 "ecsakai",
                                 "fly",
                                 "gondii",
                                 "human",
                                 "malaria",
                                 "mouse",
                                 "pig",
                                 "rat",
                                 "rhesus",
                                 "coelicolor",
                                 "worm",
                                 "xenopus",
                                 "yeast",
                                 "zebrafish"),
           envir=GOSemSimEnv)
    ## remove "coelicolor" as it is not supported by Bioconductor
}

##' get supported organisms
##'
##'
##' @title getSupported_Org
##' @return supported organisms
##' @export
##' @author Yu Guangchuang
getSupported_Org <- function() {
    if (!exists("GOSemSimEnv")) .initial()
    supported_Org <- get("SupportedSpecies", envir=GOSemSimEnv)
    return(supported_Org)
}

##' @importFrom GO.db GOMFPARENTS
##' @importFrom GO.db GOBPPARENTS
##' @importFrom GO.db GOCCPARENTS
.getParents <- function(ont) {
    Parents <- switch(ont,
                      MF = "GOMFPARENTS",
                      BP = "GOBPPARENTS",
                      CC = "GOCCPARENTS",
                      DO = "DOPARENTS"
                      )
    if (ont == "DO") {
        db <- "DO.db"
        require(db, character.only=TRUE)
    }
    Parents <- eval(parse(text=Parents))
    return(Parents)
}

