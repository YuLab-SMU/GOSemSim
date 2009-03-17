`mgeneSim` <- function (genes, ont="MF", organism="human", measure="Wang", drop="IEA"){
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast"))
	
	genes <- genes[!is.na(genes)]
	n <- length(genes)
	
	simMatrix <- matrix(NA, nrow=n, ncol=n)
	colnames(simMatrix) <- genes
	rownames(simMatrix) <- genes
	for (i in 1:n) {
		for (j in 1:i) {
			simMatrix[i,j] <- geneSim(genes[i], genes[j], wh_ont, wh_organism, wh_measure, drop)$geneSim
			if (i != j)	{
				simMatrix[j,i] <- simMatrix[i,j]
			}
		}
	}
	
	removeNA <- apply(!is.na(simMatrix), 1, sum)>0

	return(simMatrix[removeNA, removeNA])
}
