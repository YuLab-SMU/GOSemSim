setClass("GOIdentifiers", representation(GOSet1="character", GOSet2="character"))

setMethod(
	f= "sim", 
	signature= "GOIdentifiers", 
	definition=function(object, params){
		ICmethods <- c("Resnik", "Jiang", "Lin", "Rel")
		m <- length(object@GOSet1)
		n <- length(object@GOSet2)
		scores <- matrix(nrow=m, ncol=n)
		rownames(scores) <- object@GOSet1
		colnames(scores) <- object@GOSet2
		ic = params@method %in% ICmethods
		if (ic) {
			loadICdata(params)
		}
		for( i in 1:m) {
			for (j in 1:n) {
				if (ic) {
					scores[i,j] <- .infoContentMethod(object@GOSet1[i], object@GOSet2[j], ont=params@ontology, method=params@method, organism=params@organism)
				}
				if (params@method == "Wang") {
					scores[i,j] <- .wangMethod(object@GOSet1[i], object@GOSet2[j], ont=params@ontology, params@organism)
				}
			}
		}
		result <- .combineScores(scores, params@combine)
		return(result)
	}
)

