mclusterSim <- 
function(clusters, ont="MF", organism="human", measure="Wang", drop="IEA", combine="rcmax.avg") {
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	if(!exists("GOSemSimEnv")) .initial()
	wh_organism <- match.arg(organism, get("SupportedSpecies",envir=GOSemSimEnv))
	wh_combine <- match.arg(combine, c("max", "average", "rcmax", "rcmax.avg"))
	
	size <- length(clusters)
	cluster_gos=list()
	for(i in 1:size){
		cluster_gos[[i]]=sapply(clusters[[i]], function(x) ygcGetOnt(x, organism= wh_organism, ontology= wh_ont, dropCodes=drop))
	}
	simmat <- matrix(NA, nrow=size, ncol=size)
	rownames(simmat) <- names(clusters)
	colnames(simmat) <- names(clusters)
	assign("GOSemSimCache", new.env(hash=TRUE),envir=.GlobalEnv)
	for (i in 1:size) {
		for (j in 1:i) {
			gos1 <- unlist(cluster_gos[[i]])
			gos2 <- unlist(cluster_gos[[j]])
			gos1 <- gos1[!is.na(gos1)]
			gos2 <- gos2[!is.na(gos2)]
			size1 <- length(gos1)
			size2 <- length(gos2)
			if (size1 == 0 || size2 == 0) {
				simmat[i,j] <- NA
			} else {
				sim <- mgoSim(gos1,gos2, wh_ont, wh_organism, wh_measure, wh_combine)
				simmat[i,j] <- round(sim, digits=3)
			}
			if ( i != j) {
				simmat[j, i] <- simmat[i,j]
			} 
			
		}
	}	
	remove("GOSemSimCache", envir=.GlobalEnv)
	removeNA <- apply(!is.na(simmat), 1, sum) > 0
	return(simmat[removeNA, removeNA])

}
