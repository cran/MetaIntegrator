#This file was written by Winn Haynes
#Last updated: March 30, 2015

#Mimics the functionality of Purvesh's forest plot script using our database backend
#Args:
#	disease: a single diseaseName
#	gene: a single geneName
#	main: (optional) a plot title
#Returns:
#	forest plot, if possible to generate
forestPlot<- function(metaObject, geneName){
	
	if(! checkDataObject(object = metaObject, objectType="Meta", objectStage="Pre-Filter")) {
		stop("metaObject that was passed to forestPlot was not appropriately formatted as a Meta object")
	} 
	
	#If there are results for our disease and gene, generate the plot
	if(geneName %in% rownames(metaObject$metaAnalysis$datasetEffectSizes)) {
		
		#Get the gene effect sizes
		studyData <- data.frame(means=metaObject$metaAnalysis$datasetEffectSizes[geneName,])
		
		#Get the standard errors
		studyData$SEs <- metaObject$metaAnalysis$datasetEffectSizeStandardErrors[geneName,]
		
		getFormattedName <- function(uglyName) {
			if(uglyName %in% names(metaObject$originalData)) {
				return(metaObject$originalData[uglyName][[1]]$formattedName)
			}
			return("")
		}
		
		#Get the study names
		studyData$names <- sapply(rownames(studyData),getFormattedName)
		
		print(class(studyData$names))
		
		#Order the results alphabetically
		studyData<- studyData[order(studyData$names),]
		
		#Pull objects out of the query results
		studyMeans <-studyData$means
		studySEs <- studyData$SEs
		studyNames <-studyData$names
		names(studyMeans) <- studyNames
		
		#Get the pooled results for our forest plot
		pooledMean<- metaObject$metaAnalysis$pooledResults[geneName,"effectSize"]
		pooledSE<- metaObject$metaAnalysis$pooledResults[geneName, "effectSizeStandardError"]
		
		metaplot( studyMeans, studySEs, labels=studyNames,
							xlab="Standardized Mean Difference (log2 scale)", ylab="",
							colors=meta.colors(box="blue", lines="lightblue",
																 zero="black", summary="orange", text="red"), 
							summn=pooledMean, sumse=pooledSE, sumnn=1/pooledSE^2, main=geneName)
	} else {
		stop(paste("No gene named \'", geneName, "\' in this metaObject. No plot generated", sep=""))
	}
}
