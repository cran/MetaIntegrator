


#############  ROC functions  ###############################

#Args: 
#	filterObject: An object matching the specifications for Filter
#	datasetObject: Data object of the type used as input for metaIntegrator
#Generates:
#	ROC plot
rocPlot <- function(filterObject, datasetObject, title = datasetObject$formattedName){
		#Check that our objects are the right type
		if(! checkDataObject(object = filterObject, objectType="MetaFilter")) {
			stop("filterObject that was passed to rocPlot was not appropriately formatted as a MetaFilter object")
		} 
		if(!checkDataObject(object =datasetObject, objectType =  "Dataset")){
			stop("datasetObject that was passed to rocPlot was not appropriately formatted as a Dataset object")
		}
    stopifnot(length(datasetObject$class)==dim(datasetObject$expr)[2])
    stopifnot(length(levels(factor(datasetObject$class)))==2)
    stopifnot()#create check for geneList
    
    datasetObject$score <- calculateScore(filterObject, datasetObject) 
    rocObject <- calculateROC(labels=datasetObject$class, predictions=datasetObject$score)
    cat(sprintf("For dataset %s, AUC = %s", datasetObject$formattedName, signif(rocObject$auc, 3)))
    
    .rocplot.single(rocObject, title=title)    
} 

.rocplot.single <- function(rocObject, title = "ROC Plot"){
    
    annotation <- with(rocObject, paste("AUC=",signif(auc, 2), " (95%CI ", signif(auc.CI[1], 2), " - ", signif(auc.CI[2], 2), ")", sep=""))
    
    p <- ggplot(rocObject$roc, aes_string(x = "x", y = "y")) +
        geom_line(aes(colour = "")) +
        geom_segment(x=0, y=0, xend=1, yend=1, linetype=5, colour="grey75") +
        scale_x_continuous("False Positive Rate (1-Specificity)") +
        scale_y_continuous("True Positive Rate (Sensitivity)") +
        scale_colour_manual(labels = annotation, values = "#000000") +
        ggtitle(title) +
        theme(text = element_text(size=24)) +
        theme(plot.title = element_text(size=24), 
              axis.title.x = element_text(size=24),
              axis.title.y = element_text(size=24, angle=90),
              legend.justification=c(1,0), 
              legend.position=c(1,0),
              legend.title = element_blank(),
              legend.key = element_blank(),
              legend.text=element_text(size=24),
              axis.text.x = element_text(size=24),
              axis.text.y = element_text(size=24))
    return(p)
}

#Args: 
#	labels: Vector of labels; must have exactly two unique values (ie, cases and controls).
#	predictions: Vector of predictions (for instance, test scores) to be evaluated for ability to separate the two classes.
#			Must be exactly the same length as labels. 
#Returns:
#	a list of values:
#		roc: dataframe consisting of two columns, FPR and TPR, meant for plotting
#		auc: area under the curve 
#		auc.CI: 95% confidence interval for AUC
calculateROC <- function(labels, predictions, AUConly=FALSE){
    ## Assumes that labels have been passed in as numeric
    ## If not, will create ROC of earlier-sorted against later-sorted
    ## Code borrowed from ROCR package
    stopifnot(length(unique(labels))==2)
    
    levels <- sort(unique(labels))
    labels <- ordered(labels, levels = levels)
    
    n.pos <- sum(labels == levels[2])
    n.neg <- sum(labels == levels[1])
    
    pred.order <- order(predictions, decreasing = TRUE)
    predictions.sorted <- predictions[pred.order]
    
    tp <- cumsum(labels[pred.order] == levels[2])
    fp <- cumsum(labels[pred.order] == levels[1])
    dups <- rev(duplicated(rev(predictions.sorted)))
    tp <- c(0, tp[!dups])
    fp <- c(0, fp[!dups])
    fn <- n.pos - tp
    tn <- n.neg - fp
    
    tpr <- tp / (tp + fn)
    fpr <- fp / (fp + tn)
    roc = data.frame(x = fpr, y = tpr)
    roc <- roc[order(roc$x, roc$y),]
    
    i <- 2:nrow(roc)
    auc <- (roc$x[i] - roc$x[i - 1]) %*% (roc$y[i] + roc$y[i - 1])/2
    
    if(AUConly) return(auc)
    
    ## Hanley method for SE
    q1 <- auc/(2-auc)
    q2 <- (2*auc^2)/(1+auc)
    se.auc <- sqrt( (((auc * (1 - auc)) + (n.pos -1)*(q1 - auc^2)) + ((n.neg -1)*(q2 - auc^2)))/(n.pos*n.neg) )
    ci.upper <- auc + (se.auc * 1.96)
    ci.lower <- auc - (se.auc * 1.96)
    
    return(list(roc=roc, auc=auc, auc.CI=c(ci.lower, ci.upper)))
}



#####################  Backward Search  ############################
# Args:
#		metaObject: The metaObject from the main metaIntegrator function.
#		filterObject:  An object matching the specifications for Filter
#Returns:
#		A Filter  object which has results from backward search
backwardSearch <- function(metaObject, filterObject, backThresh=0){
		#Check that our objects are the right type
		if(! checkDataObject(object = metaObject, objectType="Meta", objectStage = "Pre-Filter")) {
			stop("metaObject that was passed to backwardSearch was not appropriately formatted as a Meta object")
		} 
		if(!checkDataObject(object =filterObject, objectType =  "MetaFilter")){
			stop("filterObject that was passed to backwardSearch was not appropriately formatted as a MetaFilter object")
		}
    pos.genes <- filterObject$posGeneNames
    neg.genes <- filterObject$negGeneNames
    discovery.genes <- .convertDiscoveryListToGenes(metaObject, pos.genes, neg.genes)
    
    backSearchGenes <- .backSearchWeighted(discovery.genes, pos.genes, neg.genes, backThresh)
    
    filterObject$posGeneNames <- backSearchGenes[[1]]
    filterObject$negGeneNames <- backSearchGenes[[2]]
    
    filterObject$filterDescription <- paste(filterObject$filterDescription, ". Backward search with backThresh=", backThresh, sep="")
    
    return(filterObject)
}

.backSearchWeighted <- function(discovery.genes, pos.genes, neg.genes, backThresh){
    weights <- unlist(lapply(discovery.genes, function(GEM) length(GEM$class)))
    
    .matchAll(discovery.genes, pos.genes, neg.genes, yes.pos=NULL, yes.neg=NULL)
    
    auc.weighted <- .getWeightedAUCsGenesList(discovery.genes, pos.genes, neg.genes)
    while(TRUE){
        geneSearch <- .backsearchInner(discovery.genes, pos.genes, neg.genes, auc.weighted, weights)
        
        #Remove any genes that cause an INCREASE in avg AUC when removed
        bad<-min(geneSearch[,4])
        cat("worst AUC decrease=", bad, "\n")
        if(bad<=backThresh){
            badGene <- as.character(geneSearch[geneSearch$diff==bad, "genes"][1]) ##If multiple "worst" genes, remove the first only
            auc.weighted <- geneSearch[geneSearch$diff==bad, "search"][1]
            cat("removing ", badGene, "; new weighted AUC:", auc.weighted, "\n")
            pos.genes <- setdiff(pos.genes, badGene) 
            neg.genes <- setdiff(neg.genes, badGene) 
        } else {break}
    }
    cat("Final gene set AUCs: \n")
    .getWeightedAUCsGenesList(discovery.genes, pos.genes, neg.genes) ##Show final gene set
    print(list(pos.genes,neg.genes))
    return(list(pos.genes,neg.genes))
}


.backsearchPos <- function(discovery.genes, pos.genes, neg.genes){
    auc <- lapply(discovery.genes, function(GEM) {
        #posGenes <- GEM$genes[pos.genes, ]
        posGenes <- .getGenesData(GEM$genes, pos.genes)
        
        negScore <- 0 #in case no neg genes
        if (sum(!is.na(match(neg.genes, rownames(GEM$genes)))) > 0){
            #negGenes <- GEM$genes[neg.genes, ]
            negGenes <- .getGenesData(GEM$genes, neg.genes)
            negScore <- apply(negGenes, 2, .geomMean)
        }
        
        if(length(pos.genes) <= 1){
            totalScore = scale(-negScore)[,1]
            return(calculateROC(GEM$class, totalScore, AUConly=TRUE))
        }
        
        #Since often fewer neg genes, can make them relatively less important
        ratio <- length(neg.genes)/(length(pos.genes)-1)
        
        auc <- pos.genes #preallocate
        for(j in 1:length(pos.genes)){
            #take score without each gene (back step)
            tmp <- posGenes[rownames(posGenes)!=pos.genes[j], ]
            posScore <- apply(tmp, 2, .geomMean)
            totalScore <- scale(posScore - ratio*negScore)    
            auc[j] <- calculateROC(GEM$class, totalScore, AUConly=TRUE)
        }
        return(auc)
    }) 
    
    auc <- matrix(as.numeric(unlist(auc)), nrow=length(pos.genes), 
                  ncol=length(discovery.genes), byrow=FALSE)
    return(auc)
}

.backsearchNeg <- function(discovery.genes, pos.genes, neg.genes){
    auc <- lapply(discovery.genes, function(GEM) {
        posScore <- 0 #in case no pos genes
        if (sum(!is.na(match(pos.genes, rownames(GEM$genes)))) > 0){
            #posGenes <- GEM$genes[pos.genes, ]
            posGenes <- .getGenesData(GEM$genes, pos.genes)
            posScore <- apply(posGenes, 2, .geomMean)
        }
        
        if(length(neg.genes) <= 1){
            totalScore = scale(posScore)[,1]
            return(calculateROC(GEM$class, totalScore, AUConly=TRUE))
        }    
        #negGenes <- GEM$genes[neg.genes, ]
        negGenes <- .getGenesData(GEM$genes, neg.genes)
        
        #Since often fewer neg genes, can make them relatively less important
        # -1 since LOO below
        ratio <- (length(neg.genes)-1)/length(pos.genes)
        
        auc <- neg.genes #preallocate
        for(j in 1:length(neg.genes)){  
            #take score without each gene (LOO)
            tmp <- negGenes[rownames(negGenes)!=neg.genes[j], , drop=FALSE]
            negScore <- apply(tmp, 2, .geomMean)
            totalScore <- scale(posScore - ratio*negScore)
            
            #Get AUC for given score
            totalScore = totalScore[,1]
            auc[j] <- calculateROC(GEM$class, totalScore, AUConly=TRUE)
        }    
        return(auc)
    }) 
    auc <- matrix(as.numeric(unlist(auc)), nrow=length(neg.genes), 
                  ncol=length(discovery.genes),  byrow=FALSE)
    return(auc)
}


.backsearchInner <- function(discovery.genes, pos.genes, neg.genes, auc.weighted, weights){
    if(length(pos.genes>0)){
        geneSearchPos <- data.frame(genes=pos.genes, orig=auc.weighted, search=0)
        auc <- .backsearchPos(discovery.genes=discovery.genes, 
                             pos.genes=pos.genes, neg.genes=neg.genes)
        geneSearchPos[ ,3] <- colSums(t(auc)*weights)
    } else if(length(pos.genes) == 0){
        geneSearchPos <- data.frame(NULL)
    }
    
    if(length(neg.genes>0)){
        geneSearchNeg <- data.frame(genes=neg.genes, orig=auc.weighted, search=0)
        auc <- .backsearchNeg(discovery.genes=discovery.genes, 
                             pos.genes=pos.genes, neg.genes=neg.genes)
        geneSearchNeg[ ,3] <- colSums(t(auc)*weights)
    } else if(length(neg.genes) == 0){
        geneSearchNeg <- data.frame(NULL)
    }
    
    geneSearch <- rbind(geneSearchPos, geneSearchNeg)
    geneSearch <- data.frame(geneSearch, diff=geneSearch$orig-geneSearch$search)
    return(geneSearch)
}




#####################  Forward Search  ############################
# Args:
#		metaObject: The metaObject from the main metaIntegrator function.
#		filterObject:  An object matching the specifications for Filter
#Returns:
#		A Filter  object which has results from forward search
forwardSearch <- function(metaObject, filterObject, yes.pos=NULL, yes.neg=NULL, forwardThresh=0){
		#Check that our objects are the right type
		if(! checkDataObject(object = metaObject, objectType="Meta", objectStage = "Pre-Filter")) {
			stop("metaObject that was passed to forwardSearch was not appropriately formatted as a Meta object")
		} 
		if(!checkDataObject(object =filterObject, objectType =  "MetaFilter")){
			stop("filterObject that was passed to forwardSearch was not appropriately formatted as a MetaFilter object")
		}
		
    pos.genes <- filterObject$posGeneNames
    neg.genes <- filterObject$negGeneNames
    discovery.genes <- .convertDiscoveryListToGenes(metaObject, pos.genes, neg.genes)
    forwardSearchGenes <- .forwardSearchWeighted(discovery.genes, pos.genes, neg.genes, 
                                                yes.pos=yes.pos, yes.neg=yes.neg, forwardThresh)
    
    
    filterObject$posGeneNames <- forwardSearchGenes[[1]]
    filterObject$negGeneNames <- forwardSearchGenes[[2]]
    
    filterObject$filterDescription <- paste(filterObject$filterDescription, ". Forward search with forwardThresh=", forwardThresh, sep="")
    
    return(filterObject)
}

.forwardSearchWeighted <- function(discovery.genes, pos.genes, neg.genes, 
                                  yes.pos=NULL, yes.neg=NULL, forwardThresh){
    .matchAll(discovery.genes, pos.genes, neg.genes, yes.pos, yes.neg)
    
    weights <- unlist(lapply(discovery.genes, function(GEM) length(GEM$class)))
    
    if(is.null(yes.pos) & is.null(yes.neg)){
        auc.weighted <- 0
    } else {
        pos.genes <- setdiff(pos.genes, yes.pos)
        neg.genes <- setdiff(neg.genes, yes.neg)
        auc.weighted <- .getWeightedAUCsGenesList(discovery.genes, yes.pos, yes.neg, print=FALSE)
    }
    
    while(TRUE){
        geneSearch <- .forwardSearchInner(discovery.genes, pos.genes, neg.genes, 
                                         yes.pos, yes.neg, forwardThresh, 
                                         auc.weighted, weights)
        
        #Keep the gene that increases AUC the most, as long as it's > forwardThresh
        best <- max(geneSearch$diff, na.rm=TRUE)
        cat("next best=", best, "\n")
        if(best >= forwardThresh){
            ##If multiple "best" genes, remove the first only
            bestGene <- as.character(geneSearch[geneSearch$diff==best, 1][1]) 
            bestAUC <- geneSearch[geneSearch$diff==best, 3][1]
            cat("Adding ", as.character(bestGene), bestAUC, "\n")
            if(bestGene %in% pos.genes){
                pos.genes <- setdiff(pos.genes, bestGene)
                yes.pos <- c(yes.pos, bestGene)
            } else {
                neg.genes <- setdiff(neg.genes, bestGene)
                yes.neg <- c(yes.neg, bestGene)
            }
        } else {break} ## stop when best<forwardThresh
        auc.weighted <- bestAUC
    }
    
    cat("\nFinal discovery AUCs from forward search genes:\n")
    .getAUCsGenesList(discovery.genes, yes.pos, yes.neg, print=TRUE)
    
    return(list(yes.pos, yes.neg))
}


## Function for inner loop in forwardSearchWeighted
## Calls the helper functions forwardSearchPos and forwardSearchNeg
.forwardSearchInner <- function(discovery.genes, pos.genes, neg.genes, 
                               yes.pos, yes.neg, forwardThresh, 
                               auc.weighted, weights){
    if(length(pos.genes>0)){
        geneSearchPos <- data.frame(genes=pos.genes, orig=auc.weighted, search=0)
        auc <- .forwardSearchPos(discovery.genes=discovery.genes, 
                                pos.test=pos.genes, yes.pos=yes.pos, yes.neg=yes.neg)
        geneSearchPos[ ,3] <- colSums(t(auc)*weights)
    } else {geneSearchPos <- data.frame(genes=NULL, orig=NULL, search=NULL)}   
    
    if(length(neg.genes>0)){
        geneSearchNeg <- data.frame(genes=neg.genes, orig=auc.weighted, search=0)
        auc <- .forwardSearchNeg(discovery.genes=discovery.genes, 
                                neg.test=neg.genes, yes.pos=yes.pos, yes.neg=yes.neg) 
        geneSearchNeg[ ,3] <- colSums(t(auc)*weights)
    } else {geneSearchNeg <- data.frame(genes=NULL, orig=NULL, search=NULL)}  
    
    geneSearch <- rbind(geneSearchPos, geneSearchNeg)
    #if(dim(geneSearch)[1]==0) break;
    geneSearch <- data.frame(geneSearch, diff=geneSearch$search - geneSearch$orig)
    
    return(geneSearch)
}

##Efficient
.forwardSearchPos <- function(discovery.genes, pos.test, yes.pos, yes.neg){
    auc <- lapply(discovery.genes, function(GEM) {
        #posGenes <- GEM$genes[c(pos.test, yes.pos), ]
        posGenes <- .getGenesData(GEM$genes, c(pos.test, yes.pos))
        
        negScore <- 0 #in case no neg genes
        if (sum(!is.na(match(yes.neg, rownames(GEM$genes)))) > 0){
            #negGenes <- GEM$genes[yes.neg, ]
            negGenes <- .getGenesData(GEM$genes, yes.neg)
            negScore <- apply(negGenes, 2, .geomMean)
        }
        
        #Since often fewer neg genes, can make them relatively less important
        ratio <- length(yes.neg)/(length(yes.pos)+1)
        
        auc <- pos.test #preallocate
        for(j in 1:length(pos.test)){
            #take score without each gene (LOO)
            tmp <- posGenes[c(pos.test[j], yes.pos), , drop=FALSE]
            posScore <- apply(tmp, 2, .geomMean)
            totalScore <- scale(posScore - ratio*negScore)          
            auc[j] <- calculateROC(GEM$class, totalScore, AUConly=TRUE)
        }
        return(auc)
    }) 
    
    auc <- matrix(as.numeric(unlist(auc)), nrow=length(pos.test), 
                  ncol=length(discovery.genes), byrow=FALSE)
    auc[is.na(auc)] <- 0
    return(auc)
}


.forwardSearchNeg <- function(discovery.genes, neg.test, yes.pos, yes.neg){
    auc <- lapply(discovery.genes, function(GEM) {
        posScore <- 0 #in case no pos genes
        if (sum(!is.na(match(yes.pos, rownames(GEM$genes)))) > 0){
            #posGenes <- GEM$genes[yes.pos, ]
            posGenes <- .getGenesData(GEM$genes, yes.pos)
            posScore <- apply(posGenes, 2, .geomMean)
        }
        
        #negGenes <- GEM$genes[c(neg.test, yes.neg), ]
        negGenes <- .getGenesData(GEM$genes, c(neg.test, yes.neg))
        #Since often fewer neg genes, can make them relatively less important
        if(is.null(yes.pos)){ 
            ratio <- 1
        } else {
            ratio <- (length(yes.neg)+1)/length(yes.pos)
        }
        
        auc <- neg.test #preallocate
        for(j in 1:length(neg.test)){  
            #take score without each gene (LOO)
            tmp <- negGenes[c(neg.test[j], yes.neg), , drop=FALSE]
            negScore <- apply(tmp, 2, .geomMean)
            totalScore <- scale(posScore - ratio*negScore)
            auc[j] <- calculateROC(GEM$class, totalScore, AUConly=TRUE)
        }    
        return(auc)
    }) 
    auc <- matrix(as.numeric(unlist(auc)), nrow=length(neg.test), 
                  ncol=length(discovery.genes),  byrow=FALSE)
    auc[is.na(auc)] <- 0
    return(auc)
}

######################  Helper Functions   #############################################

.getGenesData <- function(genes.mtx, genes){
    tmp <- genes.mtx[genes, ]
    rownames(tmp) <- genes
    tmp[is.na(tmp)] <- 1
    return(tmp)
}

## Get Genes Mtx for each GEM in discAUC
## subset to only the pos.genes and neg.genes 
.convertDiscoveryListToGenes <- function(metaObject, pos.genes, neg.genes){
    cat("Converting probes:genes for gene list for all discovery GEMs")
    discovery.genes <- list()
    for (i in 1:length(metaObject$originalData)){
        #Convert to genes and subset
        GEM.genes <- getSampleLevelGeneData(metaObject$originalData[[i]], c(pos.genes, neg.genes) )
        
        #scale positive
        GEMmin <- min(GEM.genes, na.rm=TRUE)
        if (GEMmin <0) {GEM.genes <- GEM.genes + abs(GEMmin) + 1}
        
        discovery.genes[[ metaObject$originalData[[i]]$formattedName ]] <- list(genes = GEM.genes, class = metaObject$originalData[[i]]$class)
        cat(".")
    }
    rm(GEM.genes, i, GEMmin)
    cat("   done. \n")
    return(discovery.genes)
}

## Redo score function to handle genes directly
.scoreGenesMtx <- function(GeneMtx, pos.genes, neg.genes){
    ## get Geom mean of genes present
    posScore <- 0 
    if (sum(pos.genes %in% rownames(GeneMtx)) > 0){
        #posMatch <- GeneMtx[pos.genes, ]
        posMatch <- .getGenesData(GeneMtx, pos.genes)
        posScore <- apply(posMatch, 2, .geomMean)
    }
    negScore <- 0
    if (sum(neg.genes %in% rownames(GeneMtx)) > 0){
        #negMatch <- GeneMtx[neg.genes, ]
        negMatch <- .getGenesData(GeneMtx, neg.genes)
        negScore <- apply(negMatch, 2, .geomMean)
    }
    
    if(length(posScore)>1){
        totalScore <- scale(posScore - negScore)
    } else {
        totalScore <- scale(-negScore)
    }
    
    return(totalScore)
}

## AUCs for a list of GENE matrices (not probes)
.getAUCsGenesList <- function(discovery.genes, yes.pos, yes.neg, print=TRUE){
    if(length(c(yes.pos,yes.neg))==0) return(0.5)
    counter <- 1
    auc <- lapply(discovery.genes, function(GEM) {
        score <- .scoreGenesMtx(GEM$genes, yes.pos, yes.neg)
        auc <- signif(calculateROC(GEM$class, score, AUConly=TRUE), 3)
        if(print) cat(names(discovery.genes)[counter], auc, "\n")
        counter <<- counter + 1
        return(auc)
    }); 
    mean.AUC <- mean(as.numeric(unlist(auc)))
    if(print) cat("mean AUC:", mean.AUC, "\n")
    return(mean.AUC)
}

## Weighted AUCs (weight by N for each dataset)
.getWeightedAUCsGenesList <- function(genes.list, yes.pos, yes.neg, print=TRUE){
    if(length(c(yes.pos,yes.neg))==0) return(0.5)
    weights <- unlist( lapply(genes.list, function(GEM) length(GEM$class)) )
    counter <- 1
    auc <- lapply(genes.list, function(GEM) {
        score <- .scoreGenesMtx(GEM$genes, yes.pos, yes.neg)
        auc <- calculateROC(GEM$class, score, AUConly=TRUE)
        counter <- counter + 1
        return(auc)
    }); 
    auc <- as.numeric(unlist(auc))
    weightedAUC <- auc * weights
    if(print) cat("sum weighted AUC:", sum(weightedAUC), "\n")
    return(sum(weightedAUC))
}

.matchAll <- function(discovery.genes, pos.genes, neg.genes, yes.pos, yes.neg){
    discGenes <- unique(unlist(lapply(discovery.genes, function(x) rownames(x$genes))))
    if( !(all(yes.pos %in% pos.genes) )) stop("not all yes.pos %in% pos.genes")
    if( !(all(yes.neg %in% neg.genes) )) stop("not all yes.neg %in% neg.genes")
    if( !(all(neg.genes %in% discGenes) )) stop("not all neg.genes %in% discovery.genes")
    if( !(all(pos.genes %in% discGenes) )) stop("not all pos.genes %in% discovery.genes")
    if( !(all(discGenes %in% c(pos.genes, neg.genes)) )) stop("not all discovery.genes %in% c(pos.genes, neg.genes)")
}
