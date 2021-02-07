#' detemine the topological cutoff for TCSS method
#'
#' @param OrgDb OrgDb object
#' @param keytype keytype
#' @param ont ontology : "BP", "MF", "CC"
#' @param combine_method "max", "BMA", "avg", "rcmax", ""rcmax.avg"
#' @param IEAdrop TRUE/FALSE
#' @param ppidata A data.frame contains positive set and negative set.
#' Positive set is PPI pairs that already verified.
#' ppidata only has three columns, column 1 and 2 are character, column 3 
#' must be logical value:TRUE/FALSE.
#'
#' @return numeric, topological cutoff for given parameters
#' @export
#'
#' @examples
#' \dontrun{
#'     library(org.Hs.eg.db)
#'     library(STRINGdb)
#'
#'     string_db <- STRINGdb$new(version = "11.0", species = 9606,
#'     score_threshold = 900)
#'     string_proteins <- string_db$get_proteins()
#'
#'     #get realtionship
#'     ppi <- string_db$get_interactions(string_proteins$protein_external_id)
#'
#'     ppi$from <- vapply(ppi$from, function(e)
#'                        strsplit(e, "9606.")[[1]][2], character(1))
#'     ppi$to <- vapply(ppi$to, function(e)
#'                        strsplit(e, "9606.")[[1]][2], character(1))
#'     len <- nrow(ppi)
#'
#'     #select length
#'     s_len <- 500
#'     loca_1 <- sample(len, s_len, replace = T)
#'     #negative set
#'     loca_2 <- sample(len, s_len, replace = T)
#'     loca_3 <- sample(len, s_len, replace = T)
#'
#'     #union as ppidata
#'     ppidata <- data.frame(pro1 = c(ppi$from[loca_1], ppi$from[loca_2]),
#'      pro2 = c(ppi$to[loca_1], ppi$to[loca_3]),
#'      label = c(rep(TRUE, s_len), rep(FALSE, s_len)),
#'      stringsAsFactors = FALSE)
#'
#'     cutoff <- tcss_cutoff(OrgDb = org.Hs.eg.db, keytype = "ENSEMBLPROT",
#'     ont = "BP", combine_method = "max", IEAdrop = FALSE, ppidata)
#' }
tcss_cutoff <- function(OrgDb = NULL, keytype = "ENTREZID", ont,
                       combine_method = "max", IEAdrop = FALSE, ppidata) {

    anno_data <- godata(OrgDb, keytype = keytype, ont = ont, computeIC = TRUE,
                        processTCSS = FALSE, cutoff = NULL)
    #cutoff is in the range of ICT value
    GO <- unique(names(anno_data@IC))
    offspring <- switch(ont,
                        MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                        BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                        CC = AnnotationDbi::as.list(GOCCOFFSPRING))
    #compute ICT value for each term
    ICT <- computeICT(GO, offspring)
    #cutoffs, all possible cutoff values
    cutoffs <- seq(0.1, max(ICT), 0.1)
    #all gene/protein that has none-zero annotation
    all_pro <- unique(anno_data@geneAnno[, keytype])
    #filter the ppidata
    filtered_ppidata <- create_filtered_ppidata(all_pro, ppidata = ppidata)
    #calcualte the similarity value for filtered_ppidata
    predict_result <- lapply(cutoffs, computePre,
                             filtered_ppidata = filtered_ppidata, OrgDb = OrgDb,
                             keytype = keytype, ont = ont,
                             combine_method = combine_method, IEAdrop = IEAdrop)

    #calculate the auc valur and F1_score
    auc_F1_score <- calc_auc_F1_score(predict_result, filtered_ppidata = filtered_ppidata)
    #decide the most appropriate cutoff
    decide_cutoff(auc_F1_score, cutoffs = cutoffs)
}

#' keep the proteins with none-zero annotations
#'
#' @param all_pro all proteins that has none-zero annotation
#' @param ppidata data.frame, already verified PPI data
#'
#' @return data.frame, annotated protein pairs and their labels
#' @noRd
create_filtered_ppidata <- function(all_pro, ppidata) {
    #check data type
    if (!(is.character(ppidata[, 1]) & is.character(ppidata[, 2]) &
        is.logical(ppidata[, 3]))) {
        stop("ppidata must be a data.frame with three columns:character, character, logical")
    }

    #remove those proteins that have zero annotations
    len1 <- ppidata[, 1] %in% all_pro
    len2 <- ppidata[, 2] %in% all_pro
    ppidata_exist <- ppidata[len1 & len2, ]
    filtered_ppidata <- unique(ppidata_exist)

    if (dim(filtered_ppidata)[1] == 0) {
        stop("the length of filtered ppidata is 0, none items have GO annotation")
    }

    if (all(filtered_ppidata[, 3]) | all(!filtered_ppidata[, 3])) {
        stop("column 3 in filtered ppidata must contain TRUE and FALSE")
    }

    message(paste("positive set's length is", sum(filtered_ppidata[, 3]),
                ", negative set's length is", sum(!filtered_ppidata[, 3])))

    return(filtered_ppidata)
}

#' compute prediction value on filtered ppidata
#'
#' @param cutoff numeric, topological cutoff
#' @param filtered_ppidata data.frame, annotated protein pairs and their labels
#' @param OrgDb OrgDb object
#' @param keytype keytype
#' @param ont "BP", "MF", "CC"
#' @param combine_method "max" "BMA", "avg", "rcmax", "rcmax.avg"
#' @param IEAdrop TRUE/FALSE
#' 
#' @return list, the prediction value for the cutoff
#' @noRd
#'
computePre <- function(cutoff, filtered_ppidata, OrgDb, keytype, ont,
                       combine_method, IEAdrop) {
    #semdata is defined with this input cutoff
    suppressMessages(semdata <- godata(OrgDb = OrgDb, keytype = keytype,
                    ont = ont, computeIC = TRUE,
                    processTCSS = TRUE, cutoff = cutoff))
    #similarity value is calculated with the semdata
    mapply(function(e, f) geneSim(e, f,
                                  semData = semdata,
                                  measure = "TCSS",
                                  combine = combine_method,
                                  drop = IEAdrop),
                             filtered_ppidata[, 1], filtered_ppidata[, 2])
}


#' calculate auc and F1-score
#'
#' @param predict_result list, prediction value for all cutoffs
#' @param filtered_ppidata data.frame, annotated protein pairs and their labels
#' 
#' @return data.frame, auc and F1-score value for different cutoffs
#' @noRd
#'
calc_auc_F1_score <- function(predict_result, filtered_ppidata) {
    #geneSim returns one value and two characters in once calculation
    len <- dim(filtered_ppidata)[1]
    value_loc <- seq(from = 1, to = len * 3, by = 3)
    #just the similarity value
    pre_value <- lapply(predict_result, function(e) as.numeric(e[value_loc]))

    #auc value
    auc <- vapply(pre_value, function(e)
        ROCR::performance(ROCR::prediction(e, filtered_ppidata[, 3],
                                           label.ordering = c(FALSE, TRUE)),
                    measure = "auc")@y.values[[1]], numeric(1))

    #F1_score at different semantic similarity cutoffs
    all_F1_score <- lapply(pre_value, function(e)
        ROCR::performance(ROCR::prediction(e, filtered_ppidata[, 3],
                                           label.ordering = c(FALSE, TRUE)),
                    measure = "f")@y.values[[1]])
    #average value
    F1_score <- vapply(all_F1_score, mean, na.rm = TRUE, numeric(1))

    #save as data.frame
    return(data.frame(auc = auc,
                      F1_score = F1_score,
                      stringsAsFactors = F))
}

#' select the most appropriate cutoff
#'
#' @param auc_F1_score data.frame, auc and F1-score value for different cutoffs
#' @param cutoffs vector, all possible cutoff values
#' 
#' @return vector, topological cutoff for given parameters
#' @noRd
#'
decide_cutoff <- function(auc_F1_score, cutoffs) {
    #product value satisfies the "both maximized" requirement
    auc_mutiply_F1 <- auc_F1_score[, "auc"] * auc_F1_score[, "F1_score"]
    #get the max product value
    loca <- which(auc_mutiply_F1 == max(auc_mutiply_F1))

    if (length(loca) == 1)  return(cutoffs[loca])
    
    #if not only one pair of auc and F1-score have same product
    #take the one with larger auc
    select_auc <- auc_F1_score[loca, "auc"]

    auc_loca <- which(select_auc == max(select_auc))

    if (length(auc_loca) == 1) return(cutoffs[loca[auc_loca]])

    #if more than one pair of auc and F1-score are both same
    #take the smaller cutoff for time saving
    return(cutoffs[loca[min(auc_loca)]])
}
