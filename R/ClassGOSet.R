setClass("GOSet", representation(GOSet1="character", GOSet2="character"))

setMethod(
	f= "sim", 
	signature= "GOSet", 
	definition=function(object, params){
		ICmethods <- c("Resnik", "Jiang", "Lin", "Rel")
		m <- length(object@GOSet1)
		n <- length(object@GOSet2)
		scores <- matrix(nrow=m, ncol=n)
		rownames(scores) <- object@GOSet1
		colnames(scores) <- object@GOSet2
		ic = params@method %in% ICmethods
		for( i in 1:m) {
			for (j in 1:n) {
				if ( is.na(object@GOSet1[i]) || is.na(object@GOSet2[j]) ) {
					scores[i,j] <- NA
				} else {
					if (ic) {
						ONTANCESTOR <- .getAncestors(params@ontology)
						scores[i,j] <- infoContentMethod(object@GOSet1[i], object@GOSet2[j], ont=params@ontology, ONTANCESTOR, method=params@method, organism=params@organism)
					}
					if (params@method == "Wang") {
						scores[i,j] <- wangMethod(object@GOSet1[i], object@GOSet2[j], ont=params@ontology)
					}
				}
			}
		}
		result <- combineScores(scores, params@combine)
		return(result)
	}
)

