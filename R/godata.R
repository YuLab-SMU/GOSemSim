##' prepare GO DATA for measuring semantic similarity
##'
##'
##' @title godata
##' @param OrgDb OrgDb object
##' @param keytype keytype
##' @param ont one of 'BP', 'MF', 'CC'
##' @param computeIC logical, whether computer IC
##' @return GOSemSimDATA object
##' @importFrom AnnotationDbi keys
##' @importFrom AnnotationDbi select
##' @importFrom AnnotationDbi metadata
##' @importFrom methods new
##' @export
##' @author Guangchuang Yu
godata <- function(OrgDb=NULL, keytype = "ENTREZID", ont, computeIC = TRUE) {
    ont <- toupper(ont)
    ont <- match.arg(ont, c("BP", "CC", "MF"))

    if (is.null(OrgDb)) {
        return(new("GOSemSimDATA",
                   ont = ont))
    }

    OrgDb <- load_OrgDb(OrgDb)
    kk <- keys(OrgDb, keytype=keytype)
    message('preparing gene to GO mapping data...')
    goAnno <- suppressMessages(
        select(OrgDb, keys=kk, keytype=keytype,
               columns=c("GO", "ONTOLOGY")))

    goAnno <- goAnno[!is.na(goAnno$GO), ]
    goAnno <- goAnno[goAnno$ONTOLOGY == ont,]
    if (computeIC) {
        message('preparing IC data...')
        IC <- computeIC(goAnno, ont)
    }

    res <- new("GOSemSimDATA",
               keys = kk,
               ont = ont,
               geneAnno = goAnno,
               metadata = metadata(OrgDb))
    if (computeIC) {
        res@IC <- IC
    }

    return(res)
}

##' Class "GOSemSimDATA"
##' This class stores IC and gene to go mapping for semantic similarity measurement
##'
##'
##' @name GOSemSimDATA-class
##' @aliases GOSemSimDATA-class
##'   show,GOSemSimDATA-method
##'
##' @docType class
##' @slot keys gene ID
##' @slot ont ontology
##' @slot IC IC data
##' @slot geneAnno gene to GO mapping
##' @slot metadata metadata
##' @exportClass GOSemSimDATA
##' @keywords classes
##' @importFrom methods setClass
setClass("GOSemSimDATA",
         representation = representation(
             keys = "character",
             ont = "character",
             IC = "numeric",
             geneAnno = "data.frame",
             metadata = "data.frame"
         ))

##' @importFrom methods setMethod
setMethod("show", signature(object = "GOSemSimDATA"),
          function(object) {
              cat("#\n# DATA for Semantic Similarity calculation ...\n#\n")
          })

