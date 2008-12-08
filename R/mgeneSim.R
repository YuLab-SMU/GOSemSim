`mgeneSim` <-
function (genes, ont="MF", drop="IEA", weight.isa = 0.8, weight.partof = 0.6){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	n <- length(genes)
	gos <- sapply(genes, ontology=wh_ont, dropCodes=drop, ygcGetOnt)
	if (!exists("GOSemSimenvironment")) 
		.initial()
	Ancestor <- assign("Ancestor", GOSemSimenvironment)
	Ancestor <- ygcStoreAncestor(gos, wh_ont)
	simMatrix <- matrix(NA, nrow=n, ncol=n)
	colnames(simMatrix) <- genes
	rownames(simMatrix) <- genes
	for (i in 1:n) {
		for (j in 1:i) {
			simMatrix[i,j] <- geneSim(genes[i], genes[j], ont=wh_ont, drop="IEA", weight.isa=0.8, weight.partof=0.6)$geneSim
			if (i != j)	simMatrix[j,i] <- simMatrix[i,j]
		}
	}
	removeNA <- apply(!is.na(simMatrix), 1, sum)>0
	rm(Ancestor)
	return(simMatrix[removeNA, removeNA])
}
