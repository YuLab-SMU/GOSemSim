##'
##' given a BLAST2GO file, this function extracts the information from it and make it use for TERM2GENE.
##' @title read.blast2go
##' @param file BLAST2GO file
##' @importFrom rlang check_installed
##' @return a data frame with two columns: GO and Gene
##' @export
read.blast2go <- function(file) {
    blast2go <- utils::read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
    go_annotation_data <- blast2go[, c("Blast.Top.Hit.GOs", "Sequence.Name")]
    check_installed('tidyr', 'for `read.blast2go()`.')
    check_installed('tidyselect', 'for `read.blaset2go()`.')
    go_annotation_data <- tidyr::separate_rows(go_annotation_data, tidyselect::all_of("Blast.Top.Hit.GOs"), sep = ", ")

    names(go_annotation_data) <- c("GO", "Gene")
    build.df <- buildGOmap(go_annotation_data)

    bind.info <- rbind(build.df, go_annotation_data)

    bind.info <- bind.info[order(bind.info$GO, bind.info$Gene), ]
    bind.info <- bind.info[!duplicated(bind.info), ]
    bind.info[, "Gene"] <- as.character(bind.info$Gene)
    
    # get GO banches information
    go.ALL <- AnnotationDbi::select(GO.db, keys(GO.db), columns(GO.db))
    need.anno <- go.ALL[go.ALL$GOID %in% unique(bind.info$GO), ]
     
    # get Ontology information
    if ("Ontology" %in% colnames(go.ALL)) {
        bind.info$Ontology <- need.anno$Ontology
        bind.info$Ontology <- gsub("molecular_function", "MF", bind.info$Ontology)
        bind.info$Ontology <- gsub("cellular_component", "CC", bind.info$Ontology)
        bind.info$Ontology <- gsub("biological_process", "BP", bind.info$Ontology)
    }

    return(bind.info[, c("GO", "Gene", "Ontology")])
}
