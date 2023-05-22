#' immunoStates deconvolution analysis on MetaIntegrator object(s)
#' 
#' @param metaObject a MetaIntegrator formatted Meta object. 
#' 
#' @return Results from immunoStates on the MetaIntegrator object are 
#' stored in $immunoStates of the original Meta object
#' 
#' @export
#' @import parallel data.table
#' @author Francesco Vallania
#' @examples 
#' \dontrun{
#' # Example won't work on tinyMetaObject because it requires real gene names
#' # Download the needed datasets for processing. 
#' sleData <- getGEOData(c("GSE11909","GSE50635", "GSE39088"))
#' 
#' # Run immunoStates
#' immunoStatesEstimates <- immunoStateDecov(sleData)
#' }
immunoStatesDecov <- function(metaObject){
  invisible("immunoStatesMatrix")
  
  #figure out if this is not a unix machine, in which case mc.cores must be 1
  numCores = 10
  get_os <- function() {
    if (.Platform$OS.type == "windows") { 
      "win"
    } else if (Sys.info()["sysname"] == "Darwin") {
      "mac" 
    } else if (.Platform$OS.type == "unix") { 
      "unix"
    } else {
      stop("Unknown OS")
    }
  }
  if(!get_os() %in% c("mac","unix")){numCores = 1}
  
  #compute expression matrices
  iSExpMat <- mclapply(metaObject$originalData,
                       function(i)
                         .extractDataForGenesDT(i,promiscProbes = F),
                       mc.cores = min(numCores,length(metaObject$originalData)))
  
  #run immunoStates and save output into a data.table
  outPut        <- lapply(iSExpMat[which(!sapply(iSExpMat,
                                                 function(i)
                                                   all(is.na(i))))],
                          function(j){
                            outDT <- as.data.table(iSdeconvolution(immunoStatesMatrix,j),
                                                   keep.rownames = T)
                            
                            outDT[,natural_killer_cell:=CD56bright_natural_killer_cell+CD56dim_natural_killer_cell]
                            outDT[,           monocyte:=CD14_positive_monocyte+CD16_positive_monocyte]
                            outDT[,             B_cell:=naive_B_cell+memory_B_cell+plasma_cell]
                            outDT[,             T_cell:=CD8_positive_alpha_beta_T_cell+CD4_positive_alpha_beta_T_cell+gamma_delta_T_cell]
                            outDT[,        granulocyte:=neutrophil+eosinophil+basophil]
                            
                            return(outDT)
                          })
  names(outPut) <- names(iSExpMat[which(!sapply(iSExpMat,
                                                function(i)
                                                  all(is.na(i))))])
  
  #return output into meta-object
  metaObject$immunoStates <- outPut
  return(metaObject)
}

#' Correct gene expression using cell proportions from immunoStates
#' 
#' @param metaObject a MetaIntegrator formatted Meta object. 
#' 
#' @return Results from immunoStates gene proportion correction on the MetaIntegrator 
#'    object are stored in $iScorrExp of the original Meta object
#' 
#' @export
#' @import parallel data.table
#' @author Francesco Vallania
#' copyright by Francesco Vallania
#' @examples 
#' \dontrun{
#' # Example won't work on tinyMetaObject because it requires real gene names
#' # Download the needed datasets for processing. 
#' sleData <- getGEOData(c("GSE11909","GSE50635", "GSE39088"))
#' 
#' # Run immunoStates
#' immunoStatesCorrected <- immunoStateGenePropCorr(sleData)
#' }
immunoStatesGenePropCorr <- function(metaObject){
  
  #introduce check to make sure you have run
  #immunoStates before correcting for cell proportions
  if(is.null(metaObject$immunoStates)){
    stop("Error: You are required to run immunoStates before correcting for cell proportions")
  }
  
  #get common GSEs
  gseIDs <- names(metaObject$immunoStates)
  
  #~~~~~~~~~~~~~~~~~~
  corrExp <- lapply(gseIDs,
                    function(gse){
                      #get iS data
                      #extract general cell proportion columns and make sure they add up to 1 
                      iSdataDT <- metaObject$immunoStates[[gse]][,25:29,with=F]
                      iSdataDT[,mdh:=1-(natural_killer_cell+monocyte+B_cell+T_cell+granulocyte)]
                      
                      #convert data.table into a matrix
                      iSdata <- as.matrix(iSdataDT)
                      row.names(iSdata) <- metaObject$immunoStates[[gse]]$rn
                      
                      #get expression data
                      expMat <- metaObject$originalData[[gse]]$expr
                      
                      #convert into real [non-log] space if in log space [be wary of NAs]
                      if(max(expMat,na.rm=T) < 50){
                        expMat <- 2^expMat
                      }
                      
                      #deal with NAs [assume undetectable so set to 0 for now]
                      expMat[is.na(expMat)] <- 0
                      
                      #regress out proportions for each gene from each sample
                      #and just return the residuals
                      res       <- stats::lm(t(expMat) ~ 1 + iSdata)
                      corrected <- t(sweep(stats::residuals(res),2L,stats::coef(res)[1L, ], '+'))
                      
                      #~~~~~~~~~~~~~~~
                      return(corrected)
                    })
  names(corrExp) <- gseIDs
  
  #return output into meta-object
  metaObject$iScorrExp <- corrExp
  return(metaObject)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#extractDataForGenesDT
#copyright by Francesco Vallania
#######################################################################################
#Quickly create gene matrix for a given subset of genes or all
#[note: much faster than before because through data.table]
#######################################################################################
.extractDataForGenesDT <- function(GEM,geneNames = NULL,promiscProbes = TRUE){
  
  #run stuff 
  if(!("key_comment" %in% names(GEM) && GEM$key_comment== "Annotation absent")
     && !("exp_comment" %in% names(GEM) && GEM$exp_comment == "Expression data is missing")){
    if(length(which(!is.na(GEM$keys)))>0){
      
      #update-> split keys that map to multiple genes
      #[note: use unlist2 from AnnotationDbi so that names are fine]
      #[note: make it into data.table]
      keyVector  <- NULL
      
      #choose to either include/exclude probes that map to more than one gene
      if(promiscProbes==TRUE){
        keyVector  <- AnnotationDbi::unlist2(sapply(GEM$keys,function(i)strsplit(i,",")))
      }else{
        keyVector  <- GEM$keys[grep(",",GEM$keys,invert = T)]
      }
      
      #create keyVector DT and remove NA cases
      #If duplicated probe names, keep only one copy
      keyVectorDT<- data.table(probe = names(keyVector)[!duplicated(names(keyVector))],
                               gene = keyVector[!duplicated(names(keyVector))])
      keyVectorDT<- keyVectorDT[!is.na(gene)]
      
      #expression DataTable
      #If duplicated probe names, keep only one copy
      exprDT    <- as.data.table(GEM$expr[!duplicated(GEM$expr),],keep.rownames = T)
      
      #merge them and remove probes [since they are now @ a gene leve]
      mergedDT  <- merge(exprDT,keyVectorDT,by.x='rn',by.y='probe')
      
      if(nrow(mergedDT)>0){
        
        mergedDT  <- mergedDT[,rn:=NULL]
        
        #filter genes
        if(!is.null(geneNames)){
          mergedDT  <- mergedDT[gene %in% geneNames]
        }
        
        #compute mean and cast
        meanDT    <- melt(mergedDT,id.vars = 'gene')[,mean(value),by=.(gene,variable)]
        meanDT    <- dcast.data.table(meanDT,formula = gene~variable,value.var = 'V1')
        
        #covert data.table to matrix
        tRn    <- meanDT$gene
        outMat <- as.matrix(meanDT[,gene:=NULL])
        row.names(outMat) <- tRn
        
        #return matrix
        return(outMat)
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
  }
  else{
    return(NA)
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#iSdeconvolution
#######################################################################################
#Run Linear Regression model on a Gene Expression MATRIX
#######################################################################################
iSdeconvolution <- function(basisMatrix,#default to immunoStates
                            geneExpressionMatrix){
  
  #format basis matrix and expression data
  basisMatrix          <- data.matrix(basisMatrix)
  geneExpressionMatrix <- data.matrix(geneExpressionMatrix)
  
  #order by row-names
  basisMatrix <- basisMatrix[order(rownames(basisMatrix)),]
  geneExpressionMatrix <- geneExpressionMatrix[order(rownames(geneExpressionMatrix)),]
  
  #convert into real [non-log] space if in log space [be wary of NAs]
  if(max(geneExpressionMatrix,na.rm=T) < 50){
    geneExpressionMatrix <- 2^geneExpressionMatrix
  }
  
  #run quantile normalization on gene expression matrix
  colN <- colnames(geneExpressionMatrix)
  rowN <- rownames(geneExpressionMatrix)
  geneExpressionMatrix <- preprocessCore::normalize.quantiles(geneExpressionMatrix)
  colnames(geneExpressionMatrix) <- colN
  rownames(geneExpressionMatrix) <- rowN
  
  #only pick genes in common between data and basis matrix
  bMgenes  <- row.names(basisMatrix)
  gMgenes  <- row.names(geneExpressionMatrix)
  GMintBM  <- gMgenes %in% bMgenes
  
  if(sum(GMintBM)==0) {
    stop("None of the gene names are present in the the basis matrix. Make sure your gene names are correct and standard.")
  }
  
  geneExpressionMatrix <- geneExpressionMatrix[GMintBM,]
  BMintGM              <- bMgenes %in% row.names(geneExpressionMatrix)
  basisMatrix          <- basisMatrix[BMintGM,]
  
  #standardize basis matrix [rescale globally]
  basisMatrix <- (basisMatrix-mean(basisMatrix,na.rm=T))/stats::sd(as.vector(basisMatrix),na.rm=T)
  
  #declare header for the output matrix
  header <- c('Sample',colnames(basisMatrix),"P-value","Correlation","RMSE")
  
  #declare empty output matrix
  output   <- matrix()
  pvalue   <- 9999
  
  #run for every sample
  sampleIndex <- 1
  while(sampleIndex <= ncol(geneExpressionMatrix)){
    #iterate one variable @ the time
    geneExpressionSample <- geneExpressionMatrix[,sampleIndex]
    
    #remove NAs [this is necessary to avoid issues downstream]
    basisMatrixSample    <- basisMatrix[which(!is.na(geneExpressionSample)),]
    geneExpressionSample <- geneExpressionSample[which(!is.na(geneExpressionSample))]
    
    #scale cell-mixture sample
    geneExpressionSample <- scale(geneExpressionSample)
    
    #run linear regression on a single sample
    decOut  <- DecLinearRegression(basisMatrixSample,geneExpressionSample)
    
    #create output vector
    out <- c(colnames(geneExpressionMatrix)[sampleIndex],
             decOut$props,
             pvalue,
             decOut$r,
             decOut$rmse)
    
    #update output matrix
    if(sampleIndex == 1){
      output <- out
    }else{
      output <- rbind(output, out)
    }
    
    #update sample inder
    sampleIndex <- sampleIndex + 1
  }
  
  #format matrix object containing all results
  outObj <- rbind(header,output)
  outObj <- outObj[,-1]
  outObj <- outObj[-1,]
  outObj <- matrix(as.numeric(unlist(outObj)),nrow=nrow(outObj))
  
  #object is a matrix of samples X cells [+ some output]
  rownames(outObj) <- colnames(geneExpressionMatrix)
  colnames(outObj) <- c(colnames(basisMatrix),"P-value","Correlation","RMSE")
  
  #return object
  return(outObj)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#DecLinearRegression
#######################################################################################
#Run Linear Regression model on a Gene Expression MATRIX [genes not probes]
#######################################################################################
DecLinearRegression <- function(xMat,yVector){
  
  #run linear model without intercept
  model    <- stats::lm(yVector ~ xMat -1)
  
  #go from coefficients to proportions
  coeff    <- model$coefficients
  coeff[coeff < 0] <- 0
  props    <- coeff/sum(coeff)
  
  #get RMSE and Correlation for deconvolution
  dec_r    <- stats::cor(model$fitted.values,yVector)
  dec_rmse <- Metrics::rmse(model$fitted.values,yVector)
  
  #return final list in output
  newList <- list("props" = props, "rmse" = dec_rmse, "r" = dec_r)
}


#declare global variables for variables in data.table/with notation to avoid R CMD CHECK notes
utils::globalVariables(c("gene","rn","value","variable","keys","immunoStatesMatrix","natural_killer_cell","CD56bright_natural_killer_cell",
                         "CD56dim_natural_killer_cell","monocyte","CD14_positive_monocyte","CD16_positive_monocyte","B_cell","naive_B_cell",
                         "memory_B_cell","plasma_cell","T_cell","CD8_positive_alpha_beta_T_cell","CD4_positive_alpha_beta_T_cell",
                         "gamma_delta_T_cell","granulocyte","neutrophil","eosinophil","basophil","gse","mdh"))
