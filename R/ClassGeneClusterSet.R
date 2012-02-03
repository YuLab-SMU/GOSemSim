setClass("GeneClusterSet", representation(GeneClusters="list"))

setMethod(
          f= "sim",
          signature= "GeneClusterSet",
          definition=function(object, params){
              if (length(params@combine)==0) {
                  stop("Using setCombineMethod(\"Params\") to specify which method to combine.")
              }
              size <- length(object@GeneClusters)
              cluster_gos=list()
              for(i in 1:size){
                  cluster_gos[[i]]=sapply(object@GeneClusters[[i]], gene2GO, params)
              }
              assign("SemSimCache", new.env(hash=TRUE),envir=.GlobalEnv)

              simScores <- matrix(NA, nrow=size, ncol=size)
              rownames(simScores) <- names(object@GeneClusters)
              colnames(simScores) <- names(object@GeneClusters)

              for (i in seq(along=object@GeneClusters)) {
                  for (j in 1:i) {
                      gos1 <- unlist(cluster_gos[[i]])
                      gos2 <- unlist(cluster_gos[[j]])
                      gos1 <- gos1[!is.na(gos1)]
                      gos2 <- gos2[!is.na(gos2)]
                      if (length(gos1) == 0 || length(gos2)== 0) {
                          simScores[i,j] <- NA
                      } else {
                          goids <- new("GOSet", GOSet1=gos1, GOSet2=gos2)
                          simScores[i,j] = sim(goids, params)
                      }
                      if (i != j ){
                          simScores[j,i] <- simScores[i,j]
                      }
                  }
              }
              remove("SemSimCache", envir=.GlobalEnv)
              removeNA <- apply(!is.na(simScores), 1, sum)>0
              return(simScores[removeNA, removeNA])
          }
          )
