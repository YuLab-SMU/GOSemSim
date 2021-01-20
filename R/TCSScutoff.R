#' Title detemine the topological cutoff for TCSS method
#'
#' @param OrgDb OrgDb object
#' @param keytype keytype
#' @param ont ontology : "BP", "MF", "CC"
#' @param combine_method "max", "BMA", "avg", "rcmax", ""rcmax.avg"
#' @param IEAdrop TRUE/FALSE
#' @param testdata data.frame with three columns:pro1, pro2, label(TRUE/FALSE)
#'
#' @importFrom ROCR performance
#' @importFrom ROCR prediction
#' @return cutoff value
#' @export
#'
#' @examples 
#' \dontrun{
#'     library(org.Hs.eg.db)
#'     get_cutoff(OrgDb = org.Hs.eg.db, keytype = "ENSEMBLPROT",
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

#screen the proteins with none-zero annotations
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

#compute prediction value
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


#calculate auc and F1-score
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

#select the most appropriate cutoff
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
