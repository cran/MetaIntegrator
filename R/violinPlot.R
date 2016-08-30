#This file was written by Winn Haynes
#Last updated: March 30, 2015

violinPlot <- function(filterObject, datasetObject, labelColumn= "label") {
	#Load relevant libraries
	
	#Check that our objects are the right type
	if(! checkDataObject(object = filterObject, objectType="MetaFilter")) {
		stop("filterObject that was passed to violinPlot was not appropriately formatted as a MetaFilter object")
	} 
	if(!checkDataObject(object =datasetObject, objectType =  "Dataset")){
		stop("datasetObject that was passed to violinPlot was not appropriately formatted as a Dataset object")
	}
	
	#Populate a pheno object
	datasetPheno <- datasetObject$pheno
	
	#Make sure the label column is present
	if(!labelColumn %in% colnames(datasetPheno)) {
		stop(paste("Column named \'",labelColumn, "\' not present in the datasetObject$pheno data frame", sep=""))
	}

	#Calculate the score
	datasetPheno$score <- calculateScore(filterObject, datasetObject, suppressMessages = TRUE)

	return(.createViolinPlot(datasetPheno, formattedName=datasetObject$formattedName, comparisonColumn=labelColumn))
}

#Set some visual parameters. For some reason these have to be outside the function call....

#Draw the violin plot with:
#	standard error bars
#	Wilcox rank sum comparison between groups

#segmentWidth=0.05
#jitterWidth=0.025

.createViolinPlot <- function(datasetPheno, comparisonColumn="label", scoreColumn="score", formattedName="") {
	
	#Set visualization parameters
	segmentWidth <- 0.05
	jitterWidth <- 0.025

	#Swap our column names to "x" and "y", so that we can work with them in the ggplot drawing
	datasetPheno$y <- datasetPheno[,scoreColumn]
	datasetPheno$x <- datasetPheno[,comparisonColumn]

  #Function to split names at spaces occurring every ~10 characters into new lines
  strSplitFun <- function(text) {
    prevPiece<- ""
    curPieces <- c()
    maxSplitLength <- 10
    for(textPiece in strsplit(text," ")[[1]]) {
      if(nchar(textPiece)+nchar(prevPiece) < maxSplitLength) {
        prevPiece <- paste(prevPiece, textPiece, sep=" ")
      } else {
        if(! prevPiece=="") {
          curPieces <- c(curPieces, prevPiece)
        }
        prevPiece <- textPiece
      }
    }
    curPieces <- c(curPieces, prevPiece)
    paste(curPieces, collapse="\n")
  }
  
	
	#Capitalize first letter for formatting
	simpleCap <- function(x) {
	  s <- strsplit(x, " ")[[1]]
	  paste(toupper(substring(s, 1,1)), substring(s, 2),
	        sep="", collapse=" ")
	}
	
	#Preserve the factor ordering
  if(class(datasetPheno$x) =="factor") {
  	factorLevels <- levels(datasetPheno$x)
	
		#Make sure that only levels which have contents are used
		factorLevels <- factorLevels[factorLevels %in% as.character(datasetPheno$x)]

  	factorLevels <- sapply(factorLevels, simpleCap)
  	factorLevels <- sapply(factorLevels, strSplitFun)
  	
  	#If capitalization has lead to duplicated factors, remove the duplicates
  	factorLevels <- factorLevels[!duplicated(factorLevels)]
  }
  
  #Start off with datasetPheno as a character so we can split text
  datasetPheno$x <- as.character(datasetPheno$x)
  datasetPheno$x <- sapply(datasetPheno$x, simpleCap)
  
	#Turn into a factor so that we can use group reasoning
	if(exists("factorLevels")) {
	  datasetPheno$x <- factor(sapply(datasetPheno$x, strSplitFun), levels=factorLevels)
	} else {
	  datasetPheno$x <- as.factor(sapply(datasetPheno$x, strSplitFun))
	}
	
	#Calculate pairwise group comparisons
	baseLevel <- 1
	baseFactor <- levels(datasetPheno$x)[[baseLevel]]
	for(curLevel in (2:length(levels(datasetPheno$x)))) {
		curFactor <- levels(datasetPheno$x)[[curLevel]]
		
		#Calculate a wilcoxon rank-sum test comparison of the two groups
		pVal <- wilcox.test(datasetPheno$y[which(datasetPheno$x==baseFactor)], datasetPheno$y[which(datasetPheno$x==curFactor)])$p.value
		
		#Calculate line segment locations
		yHeight <- .1/length(levels(datasetPheno$x)) * (max(datasetPheno$y) -min(datasetPheno$y)) 
		yLoc <- min(datasetPheno$y) - curLevel * yHeight
		
		#Save into data frames
		curPvalFrame <- data.frame(x=curLevel-0.5, y= yLoc-yHeight/2, label=sprintf("%.3e",pVal))
		
		curLineSegmentFrame <- data.frame(x=c(baseLevel, baseLevel, curLevel, curLevel),
																y=c(yLoc+yHeight*0.7, yLoc, yLoc, yLoc+yHeight*0.7), group=rep(curLevel,4))
		
		#Save the p-value to data frame
		if(! exists("pValFrame")) {
			pValFrame <- curPvalFrame
			lineSegmentFrame <- curLineSegmentFrame
		} else {
			pValFrame <- rbind(pValFrame, curPvalFrame)
			lineSegmentFrame <- rbind(lineSegmentFrame, curLineSegmentFrame)
		}
	}
	
	#Calculate standard errors
	standardError = Rmisc::summarySE(data.frame(x=datasetPheno$x, y=datasetPheno$y), measurevar="y", groupvars="x")[, c("x", "y", "se")]
	standardError$segmentWidth <- rep(segmentWidth, nrow(standardError))

	#Hack to supress NOTES when running CRAN checks:
	x <- NULL; y <- NULL; se <- NULL; label <- NULL; group <- NULL

	#Draw the ggplot
	plotObject <- ggplot(datasetPheno, aes(x,y)) + geom_violin(fill='grey',trim=FALSE) +
		geom_jitter(aes(color=x), size = 3, position=position_jitter(width=jitterWidth)) +
		xlab(simpleCap(comparisonColumn)) +
		ylab("Meta Score") +
		theme(text = element_text(size=24)) +
		theme(legend.title=element_blank(), legend.text=element_text(size=24)) +
		theme(axis.text.x = element_text(size=24), axis.text.y = element_text(size=24)) +
		geom_segment(data=standardError, aes(x=match(x,levels(x))-segmentWidth,
								 		xend=match(x,levels(x))+segmentWidth, y=y-se,yend=y-se), col='black') +
		geom_segment(data=standardError,  aes(x=match(x,levels(x))-segmentWidth,
								 		xend=match(x,levels(x))+segmentWidth, y=y+se,yend=y+se), col='black') +
		geom_segment(data=standardError, aes(x=match(x,levels(x)),
								 		xend=match(x,levels(x)), y=y+se,yend=y-se), col='black') +
		geom_point(data=standardError, aes(x=x, y=y), color="black", size=3) +
		geom_text(data=pValFrame, aes(x=x, y=y, label=label)) +
		geom_line(data=lineSegmentFrame, aes(x=x, y=y, group=group)) + 
		ggtitle(formattedName)
	print(plotObject)
}
