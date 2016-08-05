# Erika Bongen ebongen@stanford.edu
# 4/2/2015

###############################################################################
#~ checkDataObject
###############################################################################
#~ checkDataObject looks for errors withing Meta, Dataset, MetaAnalyis, or
#~ MetaFilter objects. 
#~
#~ For MetaAnalysis and MetaFilter functions, it makes sure
#~ That each entry within the object is 1) not null and 2) the correct type.
#~
#~ For Dataset objects, it makes sure that 1) the entries are not null (expect
#~ $class, which is permitted to be NULL)  2) the entries are the correct type
#~ and 3) the sample names (within $pheno, $expr, and $class) match 4) the 
#~ probeIDs (within $expr and $keys) match. 
#~
#~ For Meta objects, it recursively checks the Dataset, MetaAnalysis, and 
#~ MetaFilter objects contained within the Meta object. The objectStage defines
#~ what entries a Meta object contains. Thus, Pre-Analysis Meta objects
#~ only contain $originalData. Pre-Filter Meta objects contain $originalData, 
#~ $metaAnalysis, and $leaveOneOutAnalysis. Post-Filter Meta objects contain
#~ $originalData, $metaAnalysis, $leaveOneOutAnalysis, and $filterResults. 
#~
#~ Inputs: 
#~    object: the object to be checked for validation
#~    objectType: one of "Meta", "Dataset", "MetaAnalysis", "MetaFilter"
#~    objectStage: if a Meta object, one of "Pre-Analysis", "Pre-Filter", or "Post-Filter". Otherwise: ""
#~
#~ Generates: Prints warning messages explaining the portion of the error checking failed
#~ Returns: True if passed error checking, false if otherwise
###############################################################################

checkDataObject <- function(object, objectType, objectStage="") {
  if (objectType == "Meta") {
    result <- .metaCheckAll(object, objectStage)
  } else if (objectType == "Dataset") {
    result <- .datasetCheckAll(object)
  } else if (objectType == "MetaAnalysis") {
    result <- (.metaAnalysisCheckNull(object) && .metaAnalysisCheckType(object))
  } else if (objectType == "MetaFilter") {
    result <- .metaFilterCheckNull(object) && .metaFilterCheckType(object)
  } else {
    result <- FALSE
    warning("Invalid object type.")
  }
  
  return(result)
}
