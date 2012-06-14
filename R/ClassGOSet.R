setClass("GOSet", representation(GOSet1="character", GOSet2="character"))

setMethod(
          f          = "sim",
          signature  = "GOSet",
          definition =function(object, params){

                go1  <- object@GOSet1
                go2  <- object@GOSet2

                ont <- params@ontology
                method <- params@method
                organism <- params@organism

                scores <- termSim(go1,go2, method, organism, ont)

		result <- combineScores(scores, params@combine)
		return(result)
            }
          )

