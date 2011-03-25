.initial <- function() {
	assign("GOSemSimEnv", new.env(),.GlobalEnv)
	assign("GOSemSimCache", new.env(), .GlobalEnv)
	assign("SupportedSpecies", c("anopheles", "arabidopsis", "bovine", "canine", "chicken", "chimp", "coelicolor", "ecolik12", "ecsakai", "fly", "human", "malaria", "mouse", "pig", "rat", "rhesus", "worm", "xenopus", "yeast", "zebrafish"), envir=GOSemSimEnv)
}

.getParents <- function(ont="MF") {
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

.getOffsprings <- function(ont="MF") {
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

.getAncestors <- function(ont="MF") {
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


rebuildAllICdata <- function() {
	if(!exists("GOSemSimEnv")) .initial()
	ont <- c("BP","CC", "MF")
	species <- get("SupportedSpecies",envir=GOSemSimEnv)
	cat("------------------------------------\n")
	cat("calulating Information Content...\nSpecies:\t\tOntology\n")
	params <- new("Params")
	for (i in species) {
		#loadAnnoPkg(params) ##load annotation pkg.
		setOrganism(params) <- i
		for (j in ont) {
			setOntology(params) <- j
			cat(i)
			cat("\t\t\t")
			cat(j)
			cat("\n")
			file=paste(paste("Info_Contents", i, j, sep="_"), ".rda", sep="")
			if ( !file.exists(file) ) {
				computeIC(params)
			}
		}
	}
	cat("------------------------------------\n")
	print("done...")
}


gene2GO <-  function(gene, params) {
	gene <- as.character(gene)
	if(!exists("GOSemSimEnv")) .initial()
	if (!exists("gomap", envir=GOSemSimEnv)) {
		loadGOMap(params)
	}	
	gomap <- get("gomap", envir=GOSemSimEnv)

	qGO	<- gomap[[gene]]

    if (is.null(qGO)) {
	   	return (NA)
    }
    if (sum(!is.na(qGO)) == 0) {
    	return (NA)
    }
    
	qGO	<- qGO[qGO == params@ontology]	
	if (length(qGO) == 0) {
		return (NA)
	}    
	return( names(qGO) )
}

###########################################################
## Function for combine scores ###
###########################################################
.combineScores <- function(SimScores, combine) {

	if (length(combine) == 0) {  #if not define combine
		return(round(SimScores, digits=3))
	} else {
	}
	
	## combine was define...	
    if(!sum(!is.na(SimScores))) return (NA)
    if (is.vector(SimScores) || nrow(SimScores)==1 || ncol(SimScores)==1) {
        if (combine == "avg") {
            return(round(mean(SimScores), digits=3))
        } else {
            return (round(max(SimScores), digits=3)) 
        }
    }
    if (combine == "avg") {
		result <- mean(SimScores, na.rm=TRUE)
	} else if (combine == "max") {
		result <- max(SimScores, na.rm=TRUE)
	} else if (combine == "rcmax") {
		rowScore <- mean(apply(SimScores, 1, max, na.rm=TRUE))
		colScore <- mean(apply(SimScores, 2, max, na.rm=TRUE))
		result <- max(rowScore, colScore)
	} else if (combine == "rcmax.avg") {
		result <- sum( apply(SimScores, 1, max, na.rm=TRUE), apply(SimScores, 2, max, na.rm=TRUE) ) / sum(dim(SimScores))
	}
	
	return (round(result, digits=3))
}


###########################################################
## Method *Wang* for GO semantic similarity measuring ###
###########################################################
.wangMethod <- function(GOID1, GOID2, ont="MF", organism="human") {
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
		.getParents(ont)
	}
	Parents <- get(Parents.name, envir=GOSemSimEnv)
	
	sv.a <- 1
	sv.b <- 1
	sw <- 1
	names(sv.a) <- GOID1
	names(sv.b) <- GOID2 
	
	sv.a <- .SemVal(GOID1, ont, Parents, sv.a, sw, weight.isa, weight.partof)
	sv.b <- .SemVal(GOID2, ont, Parents, sv.b, sw, weight.isa, weight.partof)
	
	sv.a <- .uniqsv(sv.a)
	sv.b <- .uniqsv(sv.b)
	
	idx <- intersect(names(sv.a), names(sv.b))
	inter.sva <- unlist(sv.a[idx])
	inter.svb <- unlist(sv.b[idx])
	if (is.null(inter.sva) || is.null(inter.svb) || length(inter.sva) == 0 || length(inter.svb) ==0) {
		sim <- NA
	} else {
		sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
	}
	return(sim)
}

.uniqsv <- function(sv) {
	sv <- unlist(sv)
	una <- unique(names(sv))
	sv <- unlist(sapply(una, function(x) {max(sv[names(sv)==x])}))
	return (sv)
}

.SemVal_internal <- function(goid, ont, Parents, sv, w, weight.isa, weight.partof) {
	p <- Parents[goid]
	p <- unlist(p[[1]])
	if (length(p) == 0) {
		#warning(goid, " may not belong to Ontology ", ont)
		return(NA)
	}
	relations <- names(p)
	old.w <- w
	for (i in 1:length(p)) {
		if (relations[i] == "is_a") {
			w <- old.w * weight.isa
		} else {
			w <- old.w * weight.partof
		}
		names(w) <- p[i]
		sv <- c(sv,w)
		if (p[i] != "all") {
			sv <- .SemVal_internal(p[i], ont, Parents, sv, w, weight.isa, weight.partof)
		}
	}
	return (sv)
}

.SemVal <- function(goid, ont, Parents, sv, w, weight.isa, weight.partof) {
	if(!exists("GOSemSimCache")) return(.SemVal_internal(goid, ont, Parents, sv, w, weight.isa, weight.partof))
	goid.ont <- paste(goid, ont, sep=".")
	if (!exists(goid.ont, envir=GOSemSimCache)) {
	  	value <- .SemVal_internal(goid, ont, Parents, sv, w, weight.isa, weight.partof)
	  	assign(goid.ont, value, envir=GOSemSimCache)
		#cat("recompute ", goid, value, "\n")
	}
	else{
		#cat("cache ", goid, get(goid, envir=GOSemSimCache), "\n")
	}
	return(get(goid.ont, envir=GOSemSimCache))
}


###########################################################
## Information Content Based Methods for GO semantic similarity measuring ###
###########################################################
.infoContentMethod <- function(GOID1, GOID2, ont, method, organism) {
	if(!exists("GOSemSimEnv")) .initial()
	
	org.ont.IC <- paste(organism, ont, "IC", sep="")
	
	IC <- get(org.ont.IC, envir=GOSemSimEnv)
	
	# more specific term, larger IC value.
	# Normalized, all divide the most informative IC.
	# all IC values range from 0(root node) to 1(most specific node)	
	mic <- max(IC[IC!=Inf])	
	
	IC["all"] = 0	
	
	ic1 <- IC[GOID1]/mic
	ic2 <- IC[GOID2]/mic
	
	if (ic1 == 0 || ic2 == 0) return (NA)
	Ancestor.name <- switch(ont,
		MF = "MFAncestors",
		BP = "BPAncestors",
		CC = "CCAncestors"	
	)	
	if (!exists(Ancestor.name, envir=GOSemSimEnv)) {
		.getAncestors(ont)
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
	
	#Information Content of the most informative common ancestor (MICA)
	mica <- max(IC[commonAncestor])/mic  
	
	## IC is biased
	## because the IC of a term is dependent of its children but not on its parents.
	sim <- switch(method,
   	    Resnik = mica, ## Resnik does not consider how distant the terms are from their common ancestor.
   	    ## Lin and Jiang take that distance into account.
   	    Lin = 2*mica/(ic1+ic2),
   	    Jiang = 1 - min(1, -2*mica + ic1 + ic2), 
   	    Rel = 2*mica/(ic1+ic2)*(1-exp(-mica*mic))  ## mica*mic equals to the original IC value. and exp(-mica*mic) equals to the probability of the term's occurence. 
	)   	
	return (sim)
}