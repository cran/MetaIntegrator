#######################################################################################
#metaIntegrator Functions: summarizeFilterResults
#2015/03/02 3:48pm @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to contain the summarizeFilterResults function
#######################################################################################

#######################################################################################
#summarizeFilterResults
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function returns a list of filtered genes tables based on the filter run name
#######################################################################################
summarizeFilterResults <- function(metaObject,metaFilterLabel){
  
  #check metaObject here at this step
  myObjectCheck <- checkDataObject(metaObject,'Meta','Post-Filter')
  
  #in case something did not pass just 
  if(myObjectCheck==FALSE){
    stop("Error in the input object!")
  }
  
  #in case label is missing
  if(is.null(metaObject$filterResults[[metaFilterLabel]])==TRUE){
    stop("Error! The input label does not match any filter run")
  }
  
  #extract genes from filteredObject
  pos_genes <- metaObject$filterResults[[metaFilterLabel]]$posGeneNames
  neg_genes <- metaObject$filterResults[[metaFilterLabel]]$negGeneNames
  
  #select genes
  posDF <- metaObject$metaAnalysis$pooledResults[match(pos_genes,row.names(metaObject$metaAnalysis$pooledResults)),]
  negDF <- metaObject$metaAnalysis$pooledResults[match(neg_genes,row.names(metaObject$metaAnalysis$pooledResults)),]
  
  #if empty set to NULL
  #re-sort the data.frames to rank genes by their FDRs
  if(nrow(posDF)==0){
    posDF <- NULL
  }else{
    posDF <- posDF[order(posDF$effectSizeFDR),]
  }
  if(nrow(negDF)==0){
    negDF <- NULL
  }else{
    negDF <- negDF[order(negDF$effectSizeFDR),]
  }
  
  #return data frames of positive and nagative genes
  return(list(pos=posDF,
              neg=negDF))
}
