#' Title detemine the topological cutoff for TCSS method
#'
#' @param OrgDb OrgDb object
#' @param keytype keytype
#' @param ont ontology : "BP", "MF", "CC"
#' @param combine_method "max", "BMA", "avg", "rcmax", ""rcmax.avg"
#' @param IEAdrop TRUE/FALSE
#' @param testdata data.frame. PPI data contains positive set and negative set.
#' Has three columns, and colnames are:"pro1", "pro2", "label".
#' Column "pro1" and "pro2" are character,
#' Column "label" must be logical value:TRUE/FALSE.
#'
#' @return cutoff
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
#'     #union as testdata
#'     testdata <- data.frame(pro1 = c(ppi$from[loca_1], ppi$from[loca_2]),
#'      pro2 = c(ppi$to[loca_1], ppi$to[loca_3]),
#'      label = c(rep(TRUE, s_len), rep(FALSE, s_len)),
#'      stringsAsFactors = FALSE)

#'     cutoff <- get_cutoff(OrgDb = org.Hs.eg.db, keytype = "ENSEMBLPROT",
#'     ont = "BP", combine_method = "max", IEAdrop = FALSE, testdata)
#' }
get_cutoff <- function(OrgDb = NULL, keytype = "ENTREZID", ont,
                       combine_method = "max", IEAdrop = FALSE, testdata) {

    anno_data <- godata(OrgDb, keytype, ont, computeIC = T, processTCSS = F,
                        cutoff = NULL)
    #cutoff is in the range of ICT value
    GO <- unique(names(anno_data@IC))
    offspring <- switch(ont,
                        MF = AnnotationDbi::as.list(GOMFOFFSPRING),
                        BP = AnnotationDbi::as.list(GOBPOFFSPRING),
                        CC = AnnotationDbi::as.list(GOCCOFFSPRING))
    #compute ICT value for each term
    ict <- computeICT(GO, offspring)
    #cutoffs
    cutoffs <- seq(0.1, max(ict), 0.1)
    #all gene/protein that has none-zero annotation
    all_pro <- unique(anno_data@geneAnno[, keytype])

    test_set <- get_test_set(all_pro, testdata = testdata)
    #calcualte the similarity value for test_set
    predict_result <- lapply(cutoffs, computePre,
                             test_set = test_set, OrgDb = OrgDb,
                             keytype = keytype, ont = ont,
                             combine_method = combine_method, IEAdrop = IEAdrop)

    #get the auc valur and F1_score
    auc_F1_score <- get_auc_F1_score(predict_result, test_set = test_set)
    #decide the most appropriate cutoff
    decide_cutoff(auc_F1_score, cutoffs = cutoffs)
}

#' Title screen the proteins with none-zero annotations
#'
#' @param all_pro all proteins have none-zero annotations
#' @param testdata data.frame
#'
#' @return data.frame
#' @noRd
get_test_set <- function(all_pro, testdata) {
    #remove those proteins that have zero annotations
    len1 <- vapply(testdata$pro1, function(e) e %in% all_pro, logical(1))
    len2 <- vapply(testdata$pro2, function(e) e %in% all_pro, logical(1))
    testdata_in <- testdata[len1 & len2, ]
    #data de-duplication
    test_set <- testdata_in[!duplicated(testdata_in), ]

    len <- dim(test_set)[1]

    if (len == 0) {
        stop("the length of test set is 0, none items have GO annotation")
    }

    freq <- table(test_set[, "label"])

    message(paste("positive set's length is", freq["TRUE"],
                ", negative set's length is", freq["FALSE"]))

    return(test_set)
}

#' Title compute prediction value
#'
#' @param test_set data.frame
#' @param OrgDb OrgDb object
#' @param keytype keytype
#' @param ont "BP", "MF", "CC"
#' @param cutoff number 
#' @param combine_method "max" "BMA", "avg", "rcmax", "rcmax.avg" 
#' @param IEAdrop TRUE/FALSE
#'
#' @return list
#' @noRd
#'
computePre <- function(test_set, OrgDb, keytype, ont, cutoff,
                       combine_method, IEAdrop) {
    #different cutoffs have different semdata
  suppressMessages(semdata <- godata(OrgDb = OrgDb, keytype = keytype,
                    ont = ont, computeIC = TRUE,
                    processTCSS = TRUE, cutoff = cutoff))
    #similarity value is calculated under this cutoff
    mapply(function(e, f) geneSim(e, f,
                                  semData = semdata,
                                  measure = "TCSS",
                                  combine = combine_method,
                                  drop = IEAdrop),
                             test_set[, "pro1"], test_set[, "pro2"])
}


#' Title calculate auc and F1-score
#'
#' @param predict_result list
#' @param test_set data.frame
#'
#' @return data.frame
#' @noRd
#'
get_auc_F1_score <- function(predict_result, test_set) {
    #geneSim returns one value and two characters in once calculation
    len <- dim(test_set)[1]
    value_loc <- seq(from = 1, to = len * 3, by = 3)
    #get the number
    pre_value <- lapply(predict_result, function(e) as.numeric(e[value_loc]))

    #auc value
    auc <- vapply(pre_value, function(e)
        ROCR::performance(ROCR::prediction(e, test_set[, "label"]),
                    measure = "auc")@y.values[[1]], numeric(1))

    #F1_score at different semantic similarity cutoffs
    all_F1_score <- lapply(pre_value, function(e)
        ROCR::performance(ROCR::prediction(e, test_set[, "label"]),
                    measure = "f")@y.values[[1]])
    #average value
    F1_score <- vapply(all_F1_score, mean, na.rm = TRUE, numeric(1))

    #save as data.frame
    return(data.frame(auc = auc,
                      F1_score = F1_score,
                      stringsAsFactors = F))
}

#' Title select the most appropriate cutoff
#'
#' @param auc_F1_score data.frame
#' @param cutoffs vector
#'
#' @return vector
#' @noRd
#'
decide_cutoff <- function(auc_F1_score, cutoffs) {
    #product value satisfies the "both maximized" requirement
    auc_mutiply_F1 <- auc_F1_score[, "auc"] * auc_F1_score[, "F1_score"]
    #get the max product value
    loca <- which(auc_mutiply_F1 == max(auc_mutiply_F1))

    if (length(loca) == 1) {
        cutoff <- cutoffs[loca]
    } else {
        #if not only one pair of auc and F1-score have same product
        #take the one with larger auc
        s_auc <- auc_F1_score[loca, "auc"]

        auc_loca <- which(s_auc == max(s_auc))

        if (length(auc_loca) == 1) {
            cutoff <- cutoffs[loca[auc_loca]]
        } else {
            #if more than one pair of auc and F1-score are both same
            #take the smaller one
            cutoff <- cutoffs[loca[min(auc_loca)]]
        }
    }
    return(cutoff)
}
