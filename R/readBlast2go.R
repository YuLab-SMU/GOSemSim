##'
##' given a BLAST2GO file, this function extracts the information from it and make it use for TERM2GENE.
##' @title read.blast2go
##' @param file BLAST2GO file
##' @param add_indirect_GO whether add indirect GO annotation 
##' @importFrom rlang check_installed
##' @return a data frame with three columns: GENE, GO and ONTOLOGY
##' @export
read.blast2go <- function(file, add_indirect_GO = FALSE) {
    check_installed("readr", 'for `read.blast2go()`.')
    check_installed('tidyr', 'for `read.blast2go()`.')
    check_installed('tidyselect', 'for `read.blaset2go()`.')

    # blast2go <- utils::read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE) # has bugs
    blast2go <- yulab.utils::yread(file, readr::read_delim)
    x <- blast2go[, c("Sequence Name", "GO Accession", "GO Domains")]
    names(x) <- c("GENE", "GO", "ONTOLOGY")
    x <- x[!is.na(x[["GO"]]), ]

    y <- tidyr::separate_rows(x, tidyselect::all_of(c("GO", "ONTOLOGY")), sep = ", ")
    y$GO <- sub("^\\s+", "", y$GO)
    y$ONTOLOGY <- sub("^\\s+", "", y$ONTOLOGY)
    
    ont <- setNames(c("MF", "CC", "BP"), c("F", "C", "P"))
    y$ONTOLOGY <- ont[y$ONTOLOGY]
    y <- as.data.frame(y)

    if (!add_indirect_GO) {
        return(y)
    }

    buildGOmap(y)
}




