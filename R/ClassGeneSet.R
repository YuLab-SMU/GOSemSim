setClass("GeneSet", representation(GeneSet1="character", GeneSet2="character"))

setMethod(
	f= "sim", 
	signature= "GeneSet", 
	definition=function(object, params){
		if (length(params@combine)==0) {
			stop("*combine* must be setting for combining semantic similarity scores of multiple GO terms. \nUsing setCombineMethod(\"Params\") to specify which method to combine.") 
		}
		GOS1 <- lapply(object@GeneSet1, gene2GO, params)
		GOS2 <- lapply(object@GeneSet2, gene2GO, params)
		assign("GOSemSimCache", new.env(hash=TRUE),envir=.GlobalEnv)
		
		m = length(object@GeneSet1)
		n = length(object@GeneSet2)
		simScores <- matrix(NA, nrow=m, ncol=n)
		rownames(simScores) <- object@GeneSet1
		colnames(simScores) <- object@GeneSet2
		
		for (i in seq(along=object@GeneSet1)) {
			for (j in seq(along=object@GeneSet2)) {
				goids <- new("GOSet", GOSet1=GOS1[[i]], GOSet2=GOS2[[j]])
				simScores[i,j] = sim(goids, params)
			}
		}
		remove("GOSemSimCache", envir=.GlobalEnv)
		removeRowNA <- apply(!is.na(simScores), 1, sum)>0
		removeColNA <- apply(!is.na(simScores), 2, sum)>0
		return(simScores[removeRowNA, removeColNA])
	}
)