setClass("GOSet", representation(GOSet1="character", GOSet2="character"))

setMethod(
	f= "sim",
	signature= "GOSet",
          definition=function(object, params){
		ICmethods <- c("Resnik", "Jiang", "Lin", "Rel")
                go1 <- object@GOSet1
                go2 <- object@GOSet2
                go1 <- unique(go1)
                go2 <- unique(go2)

		m <- length(go1)
		n <- length(go2)
		scores <- matrix(nrow=m, ncol=n)
		rownames(scores) <- go1
		colnames(scores) <- go2
		ic = params@method %in% ICmethods
		for( i in 1:m) {
                    for (j in 1:n) {
                        if ( is.na(go1[i]) || is.na(go2[j]) ) {
                            scores[i,j] <- NA
                        } else {
                            if (ic) {
                                scores[i,j] <- infoContentMethod(go1[i],
                                                                 go2[j],
                                                                 ont=params@ontology,
                                                                 method=params@method,
                                                                 organism=params@organism)
                            }
					if (params@method == "Wang") {
                                            scores[i,j] <- wangMethod(go1[i],
                                                                      go2[j],
                                                                      ont=params@ontology)
					}
                        }
                    }
		}
		result <- combineScores(scores, params@combine)
		return(result)
            }
          )

