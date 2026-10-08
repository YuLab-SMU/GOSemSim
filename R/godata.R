#' prepare GO DATA for measuring semantic similarity
#'
#'
#' @title godata
#' @param OrgDb OrgDb object (will be removed in future, please use annoDb instead)
#' @param annoDb GO annotation database,
#' can be OrgDb, a 'GSON' object of GO annotations (e.g. the output of
#' `clusterProfiler::gson_GO()`) or a data.frame contains three columns of
#' 'GENE', 'GO' and 'ONTOLOGY'.
#' @param keytype keytype
#' @param ont one of 'BP', 'MF', 'CC'
#' @param computeIC logical, whether computer IC
#' @param processTCSS logical, whether to prepare TCSS data. TCSS requires
#' `computeIC = TRUE`, which is enabled automatically when `processTCSS = TRUE`.
#' @param cutoff topology cutoff for TCSS subgraph construction. If `NULL`,
#' empirical defaults are used; users can estimate a cutoff with [tcss_cutoff()].
#' @return GOSemSimDATA object
#' @export
#' @importFrom methods new
#' @author Guangchuang Yu
godata <- function(OrgDb=NULL, annoDb=NULL, keytype = "ENTREZID",
                   ont, computeIC = TRUE, processTCSS = FALSE, cutoff = NULL) {
    if (processTCSS) computeIC <- TRUE

    ont <- toupper(ont)
    ont <- match.arg(ont, supported_GO())

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
      md <- AnnotationDbi::metadata(annoDb)
    } else if (inherits(annoDb, 'GSON')) {
      goAnno <- parse_gson(annoDb)
      md <- gson_metadata(annoDb)
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
  if (!is.data.frame(goAnno)) {
      stop("annoDb should be a data.frame, an 'OrgDb' or a 'GSON' object.")
  }

  ## suppose 1st column is GENE ID and should contains GO and ONTOLOGY columns
  ## maybe we should force names(goAnno)[1] == "GENE"
  if (!all(c("GO", "ONTOLOGY") %in% names(goAnno))) {
    stop("annoDb as a data.frame should contains 'GO' and 'ONTOLOGY' columns.")
  }

  if (ncol(goAnno) < 3) {
      warning("annoDb as a data.frame should have at least 3 columns: GENE, GO and ONTOLOGY. Assuming the first column is GENE ID.")
  }

  return(goAnno)
}

##' @importFrom methods slot
parse_gson <- function(gson) {
  gsid2gene <- slot(gson, "gsid2gene")

  if (is.null(gsid2gene) || nrow(gsid2gene) == 0) {
    stop("'gsid2gene' of the 'GSON' object is empty.")
  }

  GO <- as.character(gsid2gene[[1]])
  ## a hand-built 'GSON' may store GO IDs without the 'GO:' prefix
  bare <- grepl("^\\d{7}$", GO)
  GO[bare] <- paste0("GO:", GO[bare])
  keep <- grepl("^GO:\\d{7}$", GO)

  if (!any(keep)) {
    stop("no GO ID found in 'gsid2gene'; 'annoDb' expects a GO gene set, ",
         "e.g. the output of `clusterProfiler::gson_GO()`.")
  }

  ## 'GSON' holds no ontology information, so ONTOLOGY comes from GO.db
  ont <- AnnotationDbi::Ontology(GO.db::GOTERM)[GO[keep]]
  drop <- is.na(ont)
  if (any(drop)) {
    warning("dropped ", sum(drop),
            " annotation(s) whose GO ID is not found in 'GO.db'.")
  }

  data.frame(
    GENE = as.character(gsid2gene[[2]])[keep][!drop],
    GO = GO[keep][!drop],
    ONTOLOGY = unname(ont[!drop]),
    stringsAsFactors = FALSE
  )
}

gson_metadata <- function(gson) {
  fields <- c("species", "gsname", "version", "accessed_date",
              "keytype", "schema_version")
  values <- vapply(fields, function(f) {
    x <- slot(gson, f)
    if (is.null(x)) NA_character_ else paste(x, collapse = ", ")
  }, character(1))

  md <- data.frame(name = fields, value = values, stringsAsFactors = FALSE)
  md[!is.na(md$value), , drop = FALSE]
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

#' Class "GOSemSimDATA"
#' This class stores IC and gene to go mapping for semantic similarity measurement
#'
#'
#' @name GOSemSimDATA-class
#' @aliases GOSemSimDATA-class show,GOSemSimDATA-method
#'
#' @docType class
#' @slot keys gene ID
#' @slot ont ontology
#' @slot IC IC data
#' @slot geneAnno gene to GO mapping
#' @slot tcssdata tcssdata
#' @slot metadata metadata
#' @exportClass GOSemSimDATA
#' @keywords classes
#' @importFrom methods setClass
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

#' @importFrom methods setMethod
setMethod(
  "show", signature(object = "GOSemSimDATA"),
  function(object) {
    cat("#\n# DATA for Semantic Similarity calculation ...\n#\n")
  }
)
