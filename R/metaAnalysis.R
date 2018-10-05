#######################################################################################
#metaAnalysis Functions
#2015/03/02 3:48pm @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to...
#######################################################################################

#######################################################################################
#classVectorChecker
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function performs a check on the number of cases/controls in the class vector
#######################################################################################
.classVectorChecker <- function(dataSet){

  #check if dataset has class vector of entries
  #equal to number of samples
  checkMark <- FALSE

  #
  if(length(dataSet$class) == ncol(dataSet$expr)){
    #check if there are only
    if(length(which(dataSet$class==0)) > 1 && length(which(dataSet$class==1)) > 1){
      checkMark <- TRUE
    }
  }

  #return checkMark
  return(checkMark)
}

#######################################################################################
#leaveOneOutMetaAnalysisWrapper
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function performs hypethreaded LooA by calling runMetaAnalysisCore
#######################################################################################
.leaveOneOutMetaAnalysisWrapper <- function(originalData,
                                            old = FALSE, maxCores=Inf){
  #use 10 cores at most
  max_cores <- min(length(originalData), maxCores)
  if(max_cores>10){
    max_cores <- 10
  }

  #perform parallelized lapply [expects parallel to be loaded]
  return(mclapply(1:length(originalData),
                  function(i)
                    invisible(.runMetaAnalysisCore(originalData[-i],old=old)),
                  mc.cores = max_cores))
}

#######################################################################################
#runMetaAnalysisCore
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function runs the core analysis for the metaAnalysis package
#######################################################################################
.runMetaAnalysisCore <- function(originalData,
                                 old = FALSE){

  #############################################
  #Create an annotation table for all datasets
  #############################################
  annDB <- .createAnnTable(originalData[[1]])

  if(length(originalData) > 1) {
    for(i in 2:length(originalData)) {
      tempAnnTable <- .createAnnTable(originalData[[i]])

      commonProbes <- match(annDB[,1], tempAnnTable[,1])
      commonProbes <- commonProbes[!is.na(commonProbes)]
      if(length(commonProbes) > 0) {
        cat("Found common probes in", i, "\n", sep=" ")
        tempAnnTable <- tempAnnTable[-commonProbes,]
      }
      annDB = rbind(annDB, tempAnnTable)
    }
  }
  rownames(annDB) <- annDB[,1]
  annDB           <- as.matrix(annDB[,2])
  colnames(annDB) <- c("symbol")

  ####################################################################
  #Combining effect-sizes
  ####################################################################
  cat("Computing effect sizes...")
  all.ES          <-  lapply( originalData,.effect.sizes) # computes effect size for every row (i.e., probe) in expression matrix
  output.REM      <- .combine.effect.sizes( all.ES )

  cat("\nComputing summary effect sizes...")
  pooled.ES       <- output.REM$pooled.estimates
  pooled.ES$p.fdr <- stats::p.adjust( pooled.ES$p.value, method="fdr" )
  pooled.ES       <- pooled.ES[ order(pooled.ES$p.fdr), ]

  ####################################################################
  #Combining p-values
  ####################################################################

  #declare outputFisher function
  output.Fisher <- NULL

  #old function with bug
  if(old==TRUE){
    cat("\nComputing Q-values...")
    all.Qvals     <-  lapply(originalData, .ttest.Qvalues)
    cat("\nComputing Fisher's output...")
    output.Fisher <- .sum.of.logs(all.Qvals)
  }else{
    #new function with fix
    #define new function
    .adjust.fisher <- function(output.Fisher, method="fdr"){
      F.Qval.up    <- stats::p.adjust(output.Fisher[, "F.pval.up"],   method=method)
      F.Qval.down  <- stats::p.adjust(output.Fisher[, "F.pval.down"], method=method)
      return(cbind(output.Fisher, F.Qval.up, F.Qval.down))
    }

    cat("\nComputing Fisher's output...")
    all.Pvals     <-  lapply(originalData, .ttest.Pvalues)
    output.Fisher <- .sum.of.logs(all.Pvals)
    output.Fisher <- .adjust.fisher(output.Fisher = output.Fisher)
  }

  ####################################################################
  #Generate and return output in form of a list
  ####################################################################

  #Create objects for effect sizes and their standard errors
  #note->originally matrices, convert into data.frames
  datasetEffectSizes             <- output.REM$g
  colnames(datasetEffectSizes) <- as.character(strsplit(colnames(datasetEffectSizes), "_g"))

  datasetEffectSizeStandardErrors <- output.REM$se.g
  colnames(datasetEffectSizeStandardErrors) <- as.character(strsplit(colnames(datasetEffectSizeStandardErrors), "_se.g"))
  #finally decided it would be matrices

  #Create pooledResults and merge with outputFisher into a single table
  pooledResults           <- pooled.ES[,c('summary',
                                          'se.summary',
                                          'p.value',
                                          'p.fdr',
                                          'tau2',
                                          'n.studies',
                                          'Q',
                                          'pval.het')]
  colnames(pooledResults) <- c("effectSize",
                               "effectSizeStandardError",
                               "effectSizePval",
                               "effectSizeFDR",
                               "tauSquared",
                               "numStudies",
                               "cochranesQ",
                               "heterogeneityPval")

  output.Fisher           <- output.Fisher[,c('F.stat.up',
                                              'F.pval.up',
                                              'F.Qval.up',
                                              'F.stat.down',
                                              'F.pval.down',
                                              'F.Qval.down')]
  colnames(output.Fisher) <- c("fisherStatUp",
                               "fisherPvalUp",
                               "fisherFDRUp",
                               "fisherStatDown",
                               "fisherPvalDown",
                               "fisherFDRDown")
  #merge them
  pooledResults <- merge(pooledResults,output.Fisher,by=0)#check to make sure I am merging correctly

  #assing row.names from row.names column
  row.names(pooledResults) <- pooledResults$Row.names
  pooledResults$Row.names  <- NULL

  #sort by effectSizeFDR and the effectSizePval
  pooledResults <- pooledResults[order(pooledResults$effectSizeFDR,pooledResults$effectSizePval),]

  #Create an analysis decription log
  analysisDescription <- "MetaAnalysis: Random Effects Model"

  #return output list
  return(list(datasetEffectSizes             = datasetEffectSizes,
              datasetEffectSizeStandardErrors = datasetEffectSizeStandardErrors,
              pooledResults                  = pooledResults,
              analysisDescription            = analysisDescription))
}

#######################################################################################
#originalDataNameChecker
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#check name structure of the original data
#######################################################################################
.originalDataNameChecker <- function(metaObject){
  return(all(make.names(names(metaObject$originalData)) == names(metaObject$originalData)))
}

#######################################################################################
#originalDataNameConverter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#convert name structure of the originalData
#######################################################################################
.originalDataNameConverter <- function(metaObject){

  #perform check for names in the original data
  if(.originalDataNameChecker(metaObject) == FALSE){
    #convert names and store old names
    old_names <- names(metaObject$originalData)
    new_names <- make.names(old_names)

    #store old name into data-structure
    for(i in 1:length(metaObject$originalData)){
      if(is.null(metaObject$originalData[[i]]$formattedName) | metaObject$originalData[[i]]$formattedName =="") {
        metaObject$originalData[[i]]$formattedName <- old_names[i]
      }
    }

    #update names
    names(metaObject$originalData) <- new_names
  }

  #return final object
  return(metaObject)
}

#######################################################################################
#filterMetaRun
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function filters a single metaObject meta-analysis run
#######################################################################################
.filterMetaRun <- function(metaAnalysis,
                           effectSizeThresh        = 0,
                           effectFDRSizeThresh     = 0.05,
                           fisherFDRThresh         = 0.05,
                           numberStudiesThresh     = 1,
                           heterogeneityPvalThresh = 0.05){

  #Filter by Effect Size here
  final_results <- metaAnalysis$pooledResults
  final_results <- final_results[which(abs(final_results$effectSize)   >= effectSizeThresh),       ]
  final_results <- final_results[which(final_results$effectSizeFDR     <= effectFDRSizeThresh),    ]
  final_results <- final_results[which(final_results$numStudies        >= numberStudiesThresh),    ]

  #special case-> filter for heterogeneity only if the threshold is greater than 0
  if(heterogeneityPvalThresh > 0){
    final_results <- final_results[which(final_results$heterogeneityPval >= heterogeneityPvalThresh),]
  }

  #define subsets of genes going up/down using Fisher
  posGeneNames <- row.names(final_results[intersect(which(final_results$fisherFDRUp   <= fisherFDRThresh),which(final_results$effectSize > 0)),])
  negGeneNames <- row.names(final_results[intersect(which(final_results$fisherFDRDown <= fisherFDRThresh),which(final_results$effectSize < 0)),])

  #create final list
  final_list <-  list(posGeneNames = posGeneNames,
                      negGeneNames = negGeneNames,
                      effectSizeThresh = effectSizeThresh,
                      FDRThresh  = effectFDRSizeThresh,
                      numberStudiesThresh = numberStudiesThresh,
                      heterogeneityPvalThresh = heterogeneityPvalThresh)
  #return final list
  return(final_list)
}


#######################################################################################
#plotESdistribution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function makes density plots of the effect sizes for each dataset
#######################################################################################
.plotESdistribution <- function(metaObject){
  #plot the ES histograms here
  es_plot <- ggplot(reshape2::melt(metaObject$metaAnalysis$datasetEffectSizes, varnames=c("Gene", "Study")),
                    aes_string(x      = "value",
                               colour = "Study")) +
    geom_density(size = 1.1)            +
    theme_bw()                          +
    scale_color_discrete(name = 'Dataset')
  return(es_plot)
}
