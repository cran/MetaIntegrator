#' Get name of most recent filter 
#' @description Given a \code{metaObject} this function will look through \code{$filterResults} 
#' for the most recent filter used and return the filter name.
#' @usage
#' 	getMostRecentFilter(metaObject)
#' @param metaObject A meta object
#' @return Name of the most recent filter
#' @author Francesco Vallania
#' @examples 
#' 	getMostRecentFilter(tinyMetaObject)
#' @keywords 
#' attribute
#' utilities
#' @export

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
