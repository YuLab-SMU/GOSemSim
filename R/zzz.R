.onLoad <- function(libname, pkgname) {
#	pkgVersion <- packageDescription(pkgname, fields="Version")
#	msg <- paste("\nWelcome to", pkgname, "version", pkgVersion, "\n")
#	packageStartupMessage(msg)
	.initial()
}

ygcCheckAnnotationPackage <- function(species){
	pkgname <- switch (species,
		anopheles	=	"org.Ag.eg.db",
		arabidopsis = "org.At.tair.db",
		bovine	= "org.Bt.eg.db",
		canine	= "org.Cf.eg.db", 
		chicken	=	"org.Gg.eg.db", 
		chimp	=	"org.Pt.eg.db",
		coelicolor = "org.Sco.eg.db",
		ecolik12 = "org.EcK12.eg.db", 		
		ecsakai	=	"org.EcSakai.eg.db", 
		fly = "org.Dm.eg.db",
		human = "org.Hs.eg.db",
		malaria	=	"org.Pf.plasmo.db", 
		mouse = "org.Mm.eg.db",
		pig	= 	"org.Ss.eg.db", 
		rat = "org.Rn.eg.db",
		rhesus	=	"org.Mmu.eg.db",
		worm = "org.Ce.eg.db", 
		xenopus	=	"org.Xl.eg.db",
		yeast = "org.Sc.sgd.db",
		zebrafish = "org.Dr.eg.db",
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
		anopheles	=	require("org.Ag.eg.db"), 
		arabidopsis = require("org.At.tair.db"),
		bovine	= require("org.Bt.eg.db"),
		canine	= require("org.Cf.eg.db"), 		
		chicken	=	require("org.Gg.eg.db"), 
		chimp	=	require("org.Pt.eg.db"), 
		coelicolor = require("org.Sco.eg.db"),
		ecolik12 = require("org.EcK12.eg.db"),
		ecsakai	=	require("org.EcSakai.eg.db"), 		
		fly = require("org.Dm.eg.db"),
		human = require("org.Hs.eg.db"),
		malaria	=	require("org.Pf.plasmo.db"), 
		mouse = require("org.Mm.eg.db"),
		pig	= require("org.Ss.eg.db"), 
		rat = require("org.Rn.eg.db"),
		rhesus	=	require("org.Mmu.eg.db"), 
		worm = require("org.Ce.eg.db"),
		xenopus	=	require("org.Xl.eg.db"),	
		yeast = require("org.Sc.sgd.db"),
		zebrafish = require("org.Dr.eg.db"),		
	)
}

ygcConvOrgName <- function(organism = "human") {
	species <- switch(organism,
		anopheles	=	"Ag", 
		arabidopsis = "At",
		bovine	= "Bt",
		canine	= "Cf", 
		chicken	=	"Gg", 
		chimp	=	"Pt", 
		coelicolor = "Sco",
		ecolik12 = "EcK12",
		ecsakai	=	"EcSakai", 
		fly = "Dm",
		human = "Hs",
		malaria	=	"Pf", 
		mouse = "Mm",
		pig	= "Ss", 
		rat = "Rn",
		rhesus	=	"Mmu", 
		worm = "Ce",
		xenopus	=	"Xl",
		yeast = "Sc",
		zebrafish = "Dr",			
	)
	return (species)
}

ygcGOMapName <- function(organism="human") {
	gomap <- switch(organism,
		anopheles	=	org.Ag.egGO, 
		arabidopsis = org.At.tairGO,
		bovine	= org.Bt.egGO,
		canine	= org.Cf.egGO, 
		chicken	=	org.Gg.egGO, 
		chimp	=	org.Pt.egGO, 
		coelicolor = org.Sco.egGO,
		ecolik12 = org.EcK12.egGO,
		ecsakai	=	org.EcSakai.egGO, 
		fly = org.Dm.egGO,
		human = org.Hs.egGO,
		malaria	=	org.Pf.plasmoGO,
		mouse = org.Mm.egGO,
		pig	= org.Ss.egGO, 
		rat = org.Rn.egGO,
		rhesus	=	org.Mmu.egGO, 
		worm = org.Ce.egGO,
		xenopus	=	org.Xl.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
	)
	return(gomap)
}

ygcGetOnt <-  function(gene, organism, ontology, dropCodes) {
	gene <- as.character(gene)
	if(!exists("GOSemSimEnv")) .initial()
	species <- ygcConvOrgName(organism)
	if (!exists(species, envir=GOSemSimEnv)) {
		ygcGetGOMap(organism, dropCodes)
	}	
	gomap <- get(species, envir=GOSemSimEnv)

	qGO	<- gomap[[gene]]

    if (is.null(qGO)) {
	   	return (NA)
    }
    if (sum(!is.na(qGO)) == 0) {
    	return (NA)
    }
    
	qGO	<- qGO[qGO == ontology]	
	if (length(qGO) == 0) {
		return (NA)
	}    
	return( names(qGO) )
}

ygcCompute_Information_Content <- function(dropCodes="NULL", ont, organism) {
	if(!exists("GOSemSimEnv")) .initial()
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, get("SupportedSpecies",envir=GOSemSimEnv))
	ygcCheckAnnotationPackage(wh_organism)
	species <- ygcConvOrgName(wh_organism)
	if (!exists(species, envir=GOSemSimEnv)) {
		ygcGetGOMap(wh_organism, dropCodes)
	}	
	gomap <- get(species, envir=GOSemSimEnv)
	
	
	######### all GO terms appearing in an given ontology ###########
	goterms <- unlist(sapply(gomap, function(x) names(x[x == ont])), use.names=FALSE) 
	
	
	require(GO.db)
	if ( !exists("ALLGOID", envir=GOSemSimEnv) ) {
		assign("ALLGOID", toTable(GOTERM), envir=GOSemSimEnv )
	}
	goids <- get("ALLGOID", envir=GOSemSimEnv)
	#goids <- toTable(GOTERM)
	
	# all go terms which belong to the corresponding ontology..
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
	cnt <- sapply(goids,function(x){ n=gocount[ Offsprings[[x]] ]; gocount[x]+sum(n[!is.na(n)])})		
	names(cnt) <- goids	
	# the probabilities of occurrence of GO terms in a specific corpus.
	p <- cnt/sum(gocount) 
	# IC of GO terms was quantified as the negative log likelihood. 	
	IC <- -log(p)

	save(IC, file=paste(paste("Info_Contents", wh_ont, organism, sep="_"), ".rda", sep=""))
}


rebuildICdata <- function(){
	if(!exists("GOSemSimEnv")) .initial()
	ont <- c("BP","CC", "MF")
	species <- get("SupportedSpecies",envir=GOSemSimEnv)
	cat("------------------------------------\n")
	cat("calulating Information Content...\nSpecies:\t\tOntology\n")
	for (i in ont) {
		for (j in species) {
			cat(j)
			cat("\t\t\t")
			cat(i)
			cat("\n")
			file=paste(paste("Info_Contents", i, j, sep="_"), ".rda", sep="")
			if ( !file.exists(file) ) {
				ygcCompute_Information_Content(ont=i, organism=j)
			}
		}
	}
	cat("------------------------------------\n")
	print("done...")
}
