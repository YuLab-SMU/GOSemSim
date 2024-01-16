##' prepare GO DATA for measuring semantic similarity
##'
##'
##' @title godata
##' @param OrgDb OrgDb object (will be removed in future, please use annoDb instead)
##' @param annoDb GO annotation database, 
##' can be OrgDb or a data.frame contains three columns of 'GENE', 'GO' and 'ONTOLOGY'.
##' @param keytype keytype
##' @param ont one of 'BP', 'MF', 'CC'
##' @param computeIC logical, whether computer IC
##' @param processTCSS logical, whether to process TCSS
##' @param cutoff cutoff of TCSS
##' @return GOSemSimDATA object
##' @importFrom AnnotationDbi keys
##' @importFrom AnnotationDbi select
##' @importFrom AnnotationDbi metadata
##' @importFrom methods new
##' @export
##' @author Guangchuang Yu
godata <- function(OrgDb = NULL, annoDb = NULL, keytype = "ENTREZID",
                   ont, computeIC = TRUE,
                   processTCSS = FALSE, cutoff = NULL) {
    if (processTCSS) computeIC <- TRUE

    ont <- toupper(ont)
    ont <- match.arg(ont, c("BP", "CC", "MF"))

    if (is.null(OrgDb) && is.null(annoDb)) {
        return(new("GOSemSimDATA",
        ont = ont
        ))
    }

    if (!is.null(OrgDb)) {
      warning("use 'annoDb' instead of 'OrgDb'")
      annoDb <- OrgDb
    }
    if (is.character(annoDb)) {
      annoDb <- load_OrgDb(annoDb) 
    }

    md <- data.frame()
    if (inherits(annoDb, 'OrgDb')) {
      goAnno <- parse_orgDb(annoDb, keytype)
      md <- metadata(annoDb)
    } else if (inherits(annoDb, 'gson')) {
      ## to be supported
    } else { # for data.frame
      goAnno <- check_goAnno(annoDb)
    }

    goAnno <- goAnno[goAnno$ONTOLOGY == ont, ]
    if (computeIC) {
        message("preparing IC data...")
        IC <- computeIC(goAnno, ont)
        if (processTCSS) {
            message("preparing TCSS data...")
            tcssdata <- process_tcss(ont, IC = IC, cutoff = cutoff)
        }
    }

    res <- new("GOSemSimDATA",
      keys = unique(goAnno[,1]),
      ont = ont,
      geneAnno = goAnno,
      metadata = md
    )
    if (computeIC) {
        res@IC <- IC
        if (processTCSS) {
            res@tcssdata <- tcssdata
        }
    }

    return(res)
}

check_goAnno <- function(goAnno) {
  # check whether the data frame contains neccessary columns.

  ## suppose 1st column is GENE ID and should contains GO and ONTOLOGY columns
  ## maybe we should force names(goAnno)[1] == "GENE"
  if (!all(c("GO", "ONTOLOGY") %in% names(goAnno))) {
    stop("annoDb as a data.frame should contains 'GO' and 'ONTOLOGY' columns.")
  }

  return(goAnno)
}

parse_orgDb <- function(OrgDb, keytype) {
    kk <- keys(OrgDb, keytype = keytype)
    message("preparing gene to GO mapping data...")
    goAnno <- suppressMessages(
        select(OrgDb,
        keys = kk, keytype = keytype,
        columns = c("GO", "ONTOLOGY")
        )
    )

    goAnno <- goAnno[!is.na(goAnno$GO), ]
    return(goAnno)
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
##' @slot tcssdata tcssdata
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
    tcssdata = "list",
    metadata = "data.frame"
  )
)

##' @importFrom methods setMethod
setMethod(
  "show", signature(object = "GOSemSimDATA"),
  function(object) {
    cat("#\n# DATA for Semantic Similarity calculation ...\n#\n")
  }
)
