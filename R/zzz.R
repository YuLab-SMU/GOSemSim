.onLoad <- function(libname, pkgname) {
#	pkgVersion <- packageDescription(pkgname, fields="Version")
#	msg <- paste("\nWelcome to", pkgname, "version", pkgVersion, "\n")
#	packageStartupMessage(msg)
	.initial()
}

ygcCheckAnnotationPackage <- function(species){
	pkgname <- switch (species,
		human = "org.Hs.eg.db",
		fly = "org.Dm.eg.db",
		mouse = "org.Mm.eg.db",
		rat = "org.Rn.eg.db",
		yeast = "org.Sc.sgd.db",
		zebrafish = "org.Dr.eg.db",
		worm = "org.Ce.eg.db",
		arabidopsis = "org.At.tair.db",
		ecolik12 = "org.EcK12.eg.db", 
		bovine	= "org.Bt.eg.db",
		canine	= "org.Cf.eg.db", 
		anopheles	=	"org.Ag.eg.db", 
		ecsakai	=	"org.EcSakai.eg.db", 
		chicken	=	"org.Gg.eg.db", 
		chimp	=	"org.Pt.eg.db", 
		malaria	=	"org.Pf.plasmo.db", 
		rhesus	=	"org.Mmu.eg.db", 
		pig	= 	"org.Ss.eg.db", 
		coelicolor = "org.Sco.eg.db",
		xenopus	=	"org.Xl.eg.db"
	)
	p <- installed.packages()
	pn <- p[,1]
	if (sum(pn==pkgname) == 0) {
		print("The corresponding annotation package did not installed in this machine.")
		print("GOSemSim will install and load it automatically.")
		#source("http://bioconductor.org/biocLite.R")
		#biocLite(pkgname)
		install.packages(pkgname,repos="http://www.bioconductor.org/packages/release/data/annotation",type="source")
	}
	switch (species,
		human = require("org.Hs.eg.db"),
		fly = require("org.Dm.eg.db"),
		mouse = require("org.Mm.eg.db"),
		rat = require("org.Rn.eg.db"),
		yeast = require("org.Sc.sgd.db"),
		zebrafish = require("org.Dr.eg.db"),
		worm = require("org.Ce.eg.db"),
		arabidopsis = require("org.At.tair.db"),
		ecolik12 = require("org.EcK12.eg.db"),
		bovine	= require("org.Bt.eg.db"),
		canine	= require("org.Cf.eg.db"), 
		anopheles	=	require("org.Ag.eg.db"), 
		ecsakai	=	require("org.EcSakai.eg.db"), 
		chicken	=	require("org.Gg.eg.db"), 
		chimp	=	require("org.Pt.eg.db"), 
		malaria	=	require("org.Pf.plasmo.db"), 
		rhesus	=	require("org.Mmu.eg.db"), 
		pig	= require("org.Ss.eg.db"), 
		coelicolor = require("org.Sco.eg.db"),
		xenopus	=	require("org.Xl.eg.db")			
	)
}

ygcConvOrgName <- function(organism = "human") {
	species <- switch(organism,
		human = "Hs",
		fly = "Dm",
		mouse = "Mm",
		rat = "Rn",
		yeast = "Sc",
		zebrafish = "Dr",
		worm = "Ce",
		arabidopsis = "At",
		ecolik12 = "EcK12",
		bovine	= "Bt",
		canine	= "Cf", 
		anopheles	=	"Ag", 
		ecsakai	=	"EcSakai", 
		chicken	=	"Gg", 
		chimp	=	"Pt", 
		malaria	=	"Pf", 
		rhesus	=	"Mmu", 
		pig	= "Ss", 
		coelicolor = "Sco",
		xenopus	=	"Xl"
	)
	return (species)
}

ygcGOMapName <- function(organism="human") {
	gomap <- switch(organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecoli = org.EcK12.egGO,
		bovine	= org.Bt.egGO,
		canine	= org.Cf.egGO, 
		anopheles	=	org.Ag.egGO, 
		ecsakai	=	org.EcSakai.egGO, 
		chicken	=	org.Gg.egGO, 
		chimp	=	org.Pt.egGO, 
		malaria	=	org.Pf.plasmoGO, 
		rhesus	=	org.Mmu.egGO, 
		pig	= org.Ss.egGO, 
		coelicolor = org.Sco.egGO,
		xenopus	=	org.Xl.egGO		
	)
	return(gomap)
}

ygcGetOnt <-  function(gene, organism, ontology, dropCodes) {
	gene <- as.character(gene)
	if(!exists("GOSemSimEnv")) .initial()
	species <- ygcConvOrgName(organism)
	species.ont <- paste(species, ontology, sep=".")
	if (!exists(species.ont, envir=GOSemSimEnv)) {
		ygcGetGOMap(organism, dropCodes, ontology)
	}	
	gomap <- get(species.ont, envir=GOSemSimEnv)

	qGO	<- gomap[[gene]]

    if (is.null(qGO)) {
	   	return (NA)
    }
    if (sum(!is.na(qGO)) == 0) {
    	return (NA)
    }
    
	qGO	<- qGO[qGO]	
	if (length(qGO) == 0) {
		return (NA)
	}    
	return( names(qGO) )
}

ygcCompute_Information_Content <- function(dropCodes="NULL", ont, organism) {
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, get("SupportedSpecies",envir=GOSemSimEnv))
	ygcCheckAnnotationPackage(wh_organism)
	if(!exists("GOSemSimEnv")) .initial()
	species <- ygcConvOrgName(wh_organism)
	species.ont <- paste(species, ont, sep=".")
	if (!exists(species.ont, envir=GOSemSimEnv)) {
		ygcGetGOMap(organism, dropCodes, ont)
	}	
	gomap <- get(species.ont, envir=GOSemSimEnv)

	
	
	######### all GO terms appearing in an annotation ###########
	goterms<-unlist(sapply(gomap, function(x) names(x)), use.names=FALSE) 
	
	
	require(GO.db)
	goids <- toTable(GOTERM)
	# all go terms which belong to the corresponding category..
	goids <- unique(goids[goids[,"Ontology"] == wh_ont, "go_id"])  	
	gocount <- table(goterms)
	goname <- names(gocount) #goid of specific organism and selected category.
	## ensure goterms not appearing in the specific annotation have 0 frequency..
	go.diff <- setdiff(goids, goname)
	m <- double(length(go.diff)) 
	names(m) <- go.diff
	gocount <- as.vector(gocount)
	names(gocount) <- goname
	gocount <- c(gocount, m)

	Offsprings.name <- switch(wh_ont,
		MF = "MFOffsprings",
		BP = "BPOffsprings",
		CC = "CCOffsprings"	
	)	
	if (!exists(Offsprings.name, envir=GOSemSimEnv)) {
		ygcGetOffsprings(wh_ont)
	}
	Offsprings <- get(Offsprings.name, envir=GOSemSimEnv)	
	cnt <- sapply(goids,function(x){ n=gocount[unlist(Offsprings[x])]; gocount[x]+sum(n[!is.na(n)])})		
	names(cnt) <- goids	
	IC<- -log(cnt/sum(gocount))		
	save(IC, file=paste(paste("Info_Contents", wh_ont, organism, sep="_"), ".rda", sep=""))
}

rebuildICdata <- function(){
	ont <- c("MF","CC", "BP")
	species <- get("SupportedSpecies",envir=GOSemSimEnv)
	cat("------------------------------------\n")
	cat("calulating Information Content...\nSpecies:\t\tOntology\n")
	for (i in ont) {
		for (j in species) {
			cat(j)
			cat("\t\t\t")
			cat(i)
			cat("\n")
			ygcCompute_Information_Content(ont=i, organism=j)
		}
	}
	cat("------------------------------------\n")
	print("done...")
}