#######################################################################################
#metaIntegrator Functions: getMostRecentFilter
#2015/03/02 3:48pm @ Stanford
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#the purpose of this code is to contain the filterGenes function
#######################################################################################

#######################################################################################
#getMostRecentFilter
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#this function returns a list of filtered genes tables 
#######################################################################################
getMostRecentFilter <- function(metaObject){
  
  #check metaObject here at this step
  myObjectCheck <- checkDataObject(metaObject,'Meta','Post-Filter')
  
  #in case something did not pass just 
  if(myObjectCheck==FALSE){
    stop("Error in the input object!")
  }
  
  #get most recent filter
  positions <- order(unlist(lapply(metaObject$filterResults,function(i) i$timestamp)),decreasing=TRUE)
  return(names(metaObject$filterResults)[positions[1]])
}
