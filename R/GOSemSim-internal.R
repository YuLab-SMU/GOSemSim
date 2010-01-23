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
		ecolik12 = "org.EcK12.eg.db"			
	)
	p <- installed.packages()
	pn <- p[,1]
	if (sum(pn==pkgname) == 0) {
		print("The corresponding annotation package did not installed in this machine.")
		print("GOSemSim will install and load it automatically.")
		source("http://bioconductor.org/biocLite.R")
		biocLite(pkgname)
	}
	switch (species,
		human = library("org.Hs.eg.db"),
		fly = library("org.Dm.eg.db"),
		mouse = library("org.Mm.eg.db"),
		rat = library("org.Rn.eg.db"),
		yeast = library("org.Sc.sgd.db"),
		zebrafish = library("org.Dr.eg.db"),
		worm = library("org.Ce.eg.db"),
		arabidopsis = library("org.At.tair.db"),
		ecolik12 = library("org.EcK12.eg.db")			
	)
}

ygcGetGOMap <- function(organism="human") {
	if(!exists("GOSemSimEnv")) .initial()
	ygcCheckAnnotationPackage(organism)
	species <- switch(organism,
		human = "Hs",
		fly = "Dm",
		mouse = "Mm",
		rat = "Rn",
		yeast = "Sc",
		zebrafish = "Dr",
		worm = "Ce",
		arabidopsis = "At",
		ecolik12 = "EcK12"
	)
	gomap <- switch(organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecoli = org.EcK12.egGO
	)
	assign(eval(species), gomap, envir=GOSemSimEnv)
}

`ygcGetOnt` <-  function(gene, organism, ontology, dropCodes) {
	if(!exists("GOSemSimEnv")) .initial()
	species <- switch(organism,
		human = "Hs",
		fly = "Dm",
		mouse = "Mm",
		rat = "Rn",
		yeast = "Sc",
		zebrafish = "Dr",
		worm = "Ce",
		arabidopsis = "At",
		ecolik12 = "EcK12"
	)
	if (!exists(species, envir=GOSemSimEnv)) {
		ygcGetGOMap(organism)
	}	
	gomap <- get(species, envir=GOSemSimEnv)
	
    allGO <- gomap[[gene]]
    if (is.null(allGO)) {
    	return (NA)
    }
    if (sum(!is.na(allGO)) == 0) {
    	return (NA)
    }
    if(!is.null(dropCodes)) { 
      evidence<-sapply(allGO, function(x) x$Evidence) 
      drop<-evidence %in% dropCodes
      allGO<-allGO[!drop]
    }
    category<-sapply(allGO, function(x) x$Ontology)
    allGO<-allGO[category %in% ontology]

    if(length(allGO)==0) return (NA)
    return (unlist(unique(names(allGO))))
}


ygcWangMethod <- function(GOID1, GOID2, ont="MF", organism="human") {
	if(!exists("GOSemSimEnv")) .initial()
	weight.isa = 0.8
	weight.partof = 0.6

	if (GOID1 == GOID2)
		return (gosim=1)
		
	Parents.name <- switch(ont,
		MF = "MFParents",
		BP = "BPParents",
		CC = "CCParents"	
	)	
	if (!exists(Parents.name, envir=GOSemSimEnv)) {
		ygcGetParents(ont)
	}
	Parents <- get(Parents.name, envir=GOSemSimEnv)			

	sv.a <- 1
	sv.b <- 1
	sw <- 1
	names(sv.a) <- GOID1
	names(sv.b) <- GOID2 
	
	sv.a <- ygcSemVal(GOID1, Parents, sv.a, sw, weight.isa, weight.partof)
	sv.b <- ygcSemVal(GOID2, Parents, sv.b, sw, weight.isa, weight.partof)
	
	sv.a <- unlist(sv.a)
	sv.b <- unlist(sv.b)
	sv.a <- uniqsv(sv.a)
	sv.b <- uniqsv(sv.b)
	
	idx <- intersect(names(sv.a), names(sv.b))
	inter.sva <- unlist(sv.a[idx])
	inter.svb <- unlist(sv.b[idx])
	sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
	return(sim)
}



uniqsv <- function(sv) {
	una <- unique(names(sv))
	sv <- unlist(sapply(una, function(x) {max(sv[names(sv)==x])}))
	return (sv)
}

ygcSemVal_internal <- function(goid, Parents, sv, w, weight.isa, weight.partof) {
	p <- Parents[goid]
	p <- unlist(p[[1]])
	if (length(p) == 0)
		return(0)
	relations <- names(p)
	old.w <- w
	for (i in 1:length(p)) {
		if (relations[i] == "isa") {
			w <- old.w * weight.isa
		} else {
			w <- old.w * weight.partof
		}
		names(w) <- p[i]
		sv <- c(sv,w)
		if (p[i] != "all") {
			sv <- ygcSemVal_internal(p[i], Parents, sv, w, weight.isa, weight.partof)
		}
	}
	return (sv)
}
ygcSemVal <- function(goid, Parents, sv, w, weight.isa, weight.partof) {
	if(!exists("GOSemSimCache")) return(ygcSemVal_internal(goid, Parents, sv, w, weight.isa, weight.partof))
	if (!exists(goid, envir=GOSemSimCache)) {
	  	value <- ygcSemVal_internal(goid, Parents, sv, w, weight.isa, weight.partof)
	  	assign(goid, value, envir=GOSemSimCache)
		#cat("recompute ", goid, value, "\n")
	}
	else{
		#cat("cache ", goid, get(goid, envir=GOSemSimCache), "\n")
	}
	return(get(goid, envir=GOSemSimCache))
}

`ygcInfoContentMethod` <- function(GOID1, GOID2, ont, measure, organism) {
	if(!exists("GOSemSimEnv")) .initial()
	fname <- paste("Info_Contents", ont, organism, sep="_")
	tryCatch(utils::data(list=fname, package="GOSemSim", envir=GOSemSimEnv))
	Info.contents <- get("IC", envir=GOSemSimEnv)

	rootCount <- max(Info.contents[Info.contents != Inf])
	
	p1 <- Info.contents[GOID1]/rootCount
	p2 <- Info.contents[GOID2]/rootCount    

	if (p1 == 0 || p2 == 0) return (NA)
	Ancestor.name <- switch(ont,
		MF = "MFAncestors",
		BP = "BPAncestors",
		CC = "CCAncestors"	
	)	
	if (!exists(Ancestor.name, envir=GOSemSimEnv)) {
		ygcGetAncestors(ont)
	}
	
	Ancestor <- get(Ancestor.name, envir=GOSemSimEnv)					
	ancestor1 <- unlist(Ancestor[GOID1])
	ancestor2 <- unlist(Ancestor[GOID2])
	if (GOID1 == GOID2) {
		commonAncestor <- GOID1
	} else if (GOID1 %in% ancestor2) {
		commonAncestor <- GOID1
	} else if (GOID2 %in% ancestor1) {
		commonAncestor <- GOID2
	} else { 
		commonAncestor <- intersect(ancestor1, ancestor2)
	}
	if (length(commonAncestor) == 0) return (NA)
	pms <- min(Info.contents[commonAncestor], na.rm=TRUE)/rootCount
	sim<-switch(measure,
   	    Resnik = pms,
   	    Lin = pms/(p1+p2),
   	    Jiang = 1 - min(1, -2*pms + p1 + p2), 
   	    Rel = 2*pms/(p1+p2)*(1-exp(-pms))
	)   	
	return (sim)
}



ygcCompute_Information_Content <- function(dropCodes="NULL", ont, organism) {
	wh_ont <- match.arg(ont, c("MF", "BP", "CC"))
	wh_organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12"))
	gomap <- switch(wh_organism,
		human = org.Hs.egGO,
		fly = org.Dm.egGO,
		mouse = org.Mm.egGO,
		rat = org.Rn.egGO,
		yeast = org.Sc.sgdGO,
		zebrafish = org.Dr.egGO,
		worm = org.Ce.egGO,
		arabidopsis = org.At.tairGO,
		ecolik12 = org.EcK12.egGO
	)
	mapped_genes <- mappedkeys(gomap)
	gomap = AnnotationDbi::as.list(gomap[mapped_genes])
	if (!is.null(dropCodes)){
		gomap<-sapply(gomap,function(x) sapply(x,function(y) c(y$Evidence %in% dropCodes, y$Ontology %in% wh_ont)))
		gomap<-sapply(gomap, function(x) x[2,x[1,]=="FALSE"])
		gomap<-gomap[sapply(gomap,length) >0]		
	}else {
		gomap <- sapply(gomap,function(x) sapply(x,function(y) y$Ontology %in% wh_ont))
	}
	
	goterms<-unlist(sapply(gomap, function(x) names(x)), use.names=FALSE) # all GO terms appearing in an annotation	
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
	cnt <- sapply(goids,function(x){ c=gocount[unlist(Offsprings[x])]; gocount[x]+sum(c[!is.na(c)])})		
	names(cnt) <- goids	
	IC<- -log(cnt/sum(gocount))		
	save(IC, file=paste(paste("Info_Contents", wh_ont, organism, sep="_"), ".rda", sep=""))
}

rebuildICdata <- function(){
	ont <- c("MF","CC", "BP")
	species <- c("human", "rat", "mouse", "fly", "yeast", "zebrafish", "arabidopsis","worm", "ecolik12") 
	cat("------------------------------------\n")
	cat("calulating Information Content...\n species: ")
	for (i in ont) {
		for (j in species) {
			cat(j)
			cat("\t ontology: ")
			cat(i)
			cat("\n")
			ygcCompute_Information_Content(ont=i, organism=j)
		}
	}
	cat("------------------------------------\n")
	print("done...")
}