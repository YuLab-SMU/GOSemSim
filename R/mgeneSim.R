`mgeneSim` <-
function (genes, ont="MF", drop="IEA", measure="Resnik"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	
	genes <- genes[!is.na(genes)]
	n <- length(genes)
	gos <- sapply(genes, ontology=wh_ont, dropCodes=drop, ygcGetOnt)
	gos <- gos[!is.na(gos)]
	genes <- genes[!is.na(gos)]
	
	simMatrix <- matrix(NA, nrow=n, ncol=n)
	colnames(simMatrix) <- genes
	rownames(simMatrix) <- genes
	for (i in 1:n) {
		for (j in 1:i) {
			simMatrix[i,j] <- geneSim(genes[i], genes[j], wh_ont, drop, wh_measure)$geneSim
			if (i != j)	simMatrix[j,i] <- simMatrix[i,j]
		}
	}
	removeNA <- apply(!is.na(simMatrix), 1, sum)>0

	return(simMatrix[removeNA, removeNA])
}
