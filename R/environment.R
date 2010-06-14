.initial <- function() {
	assign("GOSemSimEnv", new.env(),.GlobalEnv)
	assign("GOSemSimCache", new.env(), .GlobalEnv)
	assign("SupportedSpecies", c("anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "coelicolor", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast", "zebrafish"), envir=GOSemSimEnv)
}

ygcGetParents <- function(ont="MF") {
	if(!exists("GOSemSimEnv")) .initial()
	wh_Parents <- switch(ont,
		MF = "MFParents",
		BP = "BPParents",
		CC = "CCParents"	
	)		
	Parents <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFPARENTS) ,
		BP = AnnotationDbi::as.list(GOBPPARENTS) , 
		CC = AnnotationDbi::as.list(GOCCPARENTS)	
	)
	assign(eval(wh_Parents), Parents, envir=GOSemSimEnv)
}

ygcGetOffsprings <- function(ont="MF") {
	if(!exists("GOSemSimEnv")) .initial()
	wh_Offsprings <- switch(ont,
		MF = "MFOffsprings",
		BP = "BPOffsprings",
		CC = "CCOffsprings"	
	)		
	Offsprings <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFOFFSPRING) ,
		BP = AnnotationDbi::as.list(GOBPOFFSPRING) , 
		CC = AnnotationDbi::as.list(GOCCOFFSPRING)	
	)
	assign(eval(wh_Offsprings), Offsprings, envir=GOSemSimEnv)
}

ygcGetAncestors <- function(ont="MF") {
	if(!exists("GOSemSimEnv")) .initial()
	wh_Ancestors <- switch(ont,
		MF = "MFAncestors",
		BP = "BPAncestors",
		CC = "CCAncestors"	
	)		
	Ancestors <- switch(ont,
		MF = AnnotationDbi::as.list(GOMFANCESTOR) ,
		BP = AnnotationDbi::as.list(GOBPANCESTOR) , 
		CC = AnnotationDbi::as.list(GOCCANCESTOR)	
	)
	assign(eval(wh_Ancestors), Ancestors, envir=GOSemSimEnv)
}

ygcGetGOMap <- function(organism="human", dropCodes) {
	if(!exists("GOSemSimEnv")) .initial()
	ygcCheckAnnotationPackage(organism)
	species <- ygcConvOrgName(organism)
	gomap <- ygcGOMapName(organism)
	mapped_genes <- mappedkeys(gomap)
	gomap = AnnotationDbi::as.list(gomap[mapped_genes])
	if (!is.null(dropCodes)){
		gomap <- sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% dropCodes, y$Ontology)))
		gomap <- sapply(gomap, function(x) x[2,x[1,]=="FALSE"]) ## filt out dropCodes Evidence.
		gomap <- gomap[sapply(gomap,length) >0]		
	}else {
		gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology))
	}
	assign(eval(species), gomap, envir=GOSemSimEnv)
}

ygcLoadIC <- function(ont, organism) {
	if(!exists("GOSemSimEnv")) .initial()
	fname <- paste("Info_Contents", ont, organism, sep="_")
	tryCatch(utils::data(list=fname, package="GOSemSim"))
	IC <- get("IC")
	org.ont.IC <- paste(organism, ont, "IC", sep="")
	assign(eval(org.ont.IC), IC, envir=GOSemSimEnv)
	rm (IC)
}

