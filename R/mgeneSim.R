mgeneSim <- function (genes, ont="MF", organism="human", measure="Wang", drop="IEA", combine="rcmax.avg"){
	genes <- as.character(genes)
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_measure <- match.arg(measure, c("Resnik", "Jiang", "Lin", "Rel", "Wang"))
	if(!exists("GOSemSimEnv")) .initial()
	wh_organism <- match.arg(organism, get("SupportedSpecies",envir=GOSemSimEnv))
	wh_combine <- match.arg(combine, c("max", "average", "rcmax", "rcmax.avg"))
	
	genes <- genes[!is.na(genes)]
	n <- length(genes)
	
	if (n < 2) {
		stop("gene vector must longer than one.")
	}
	
	simMatrix <- matrix(NA, nrow=n, ncol=n)
	colnames(simMatrix) <- genes
	rownames(simMatrix) <- genes

	gos <- lapply(genes, function(x) ygcGetOnt(x, organism= wh_organism, ontology= wh_ont, dropCodes=drop))

	assign("GOSemSimCache", new.env(hash=TRUE),envir=.GlobalEnv)
	for (i in 1:n) {
		for (j in 1:i) {

			#simMatrix[i,j] <- geneSim(genes[i], genes[j], wh_ont, wh_organism, wh_measure, drop)$geneSim
			
			simMatrix[i,j] <- NA
			if(any(!is.na(gos[i])) &&  any(!is.na(gos[j])))
			{
				sim <- mgoSim(gos[i],gos[j], wh_ont, wh_organism, wh_measure, wh_combine)
				sim <- round(sim, digits=3)
				simMatrix[i,j] <- sim
			}
			if (i != j)	{
				simMatrix[j,i] <- simMatrix[i,j]
			}
		}
	}
	remove("GOSemSimCache", envir=.GlobalEnv)
	
	removeNA <- apply(!is.na(simMatrix), 1, sum)>0

	return(simMatrix[removeNA, removeNA])
	
}
