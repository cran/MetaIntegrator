# Erika Bongen ebongen@stanford.edu
#4/2/2015

###############################################################################
#~ .datasetCheckAll
###############################################################################
#~ .datasetCheckAll checks whether the items within a dataset object
#~ are NULL, have sample IDs listed in a consistent orders, have probe
#~ IDs listed in a consistant order, and are the correct type. 
#~
#~ Inputs:
#~    dataset: a Dataset object to be checked (e.g. that class, pheno, expr are all present)
#~ Generates: Prints warning messages explaining the portion of the error checking failed
#~ Returns: True if passed error checking, false if otherwise.
###############################################################################
.datasetCheckAll <- function(dataset) {
  
  nullOK <- .datasetCheckNull(dataset)
  rowColOK <- .datasetCheckRowColNames(dataset)
  typeOK <- .datasetCheckType(dataset)
  if(!nullOK | !rowColOK | !typeOK) {
	  return(FALSE)
  }
  keysOK <- .datasetCheckKeys(dataset)
  exprOK <- .datasetCheckExpr(dataset)
  
  return(keysOK && exprOK)
}

###############################################################################
#~ .datasetCheckKeys
###############################################################################
#~ .datasetCheckNull warns the user and returns FALSE if all keys are NA
#~
#~ Inputs: datasetObject: a dataset object
#~ Generates: Warning messages if $keys is all NA
#~ Returns: TRUE if not all NA. FALSE if all NA.
###############################################################################
.datasetCheckKeys <- function(datasetObject) {
  if(sum(!is.na(datasetObject$keys)) == 0) {
    warning(paste("datasetObject$keys are all NA for ",datasetObject$formattedName,". Check your probe mappings to gene names."))
    return(FALSE)
  } 
  return(TRUE)
}

###############################################################################
#~ .datasetCheckExpr
###############################################################################
#~ .datasetCheckNull warns the user and returns FALSE if there is an infinite 
#~ value in expr
#~
#~ Inputs: datasetObject: a dataset object
#~ Generates: Warning messages if infinite value in expr
#~ Returns: TRUE if no infinite values in expr. FALSE if infinite value in expr.
###############################################################################
.datasetCheckExpr <- function(datasetObject) {
  if(sum(is.infinite(datasetObject$expr)) > 0) {
    warning(paste("datasetObject$expr contains at least one infinite value in ",datasetObject$formattedName,"."))
    return(FALSE)
  } 
  return(TRUE)
}

###############################################################################
#~ .datasetCheckNull
###############################################################################
#~ .datasetCheckNull warns the user and returns FALSE if they have forgetten 
#~ to set expr, pheno, or keys. It does not warn the user if the class
#~ vector is null. 
#~
#~ Inputs: object: a dataset object
#~ Generates: Warning messages if $expr, $pheno, or $keys is null
#~ Returns: TRUE if nothing is null. FALSE if something is null.
###############################################################################
.datasetCheckNull <- function(datasetObject) {
  isThere <- TRUE
  if(is.null(datasetObject$expr)) {
    warning("datasetObject$expr is null. The expression data is missing.")
    isThere <- FALSE
  } 
  if(is.null(datasetObject$pheno)) {
    warning("datasetObject$pheno is null. The phenotype data is missing.")
    isThere <- FALSE
  }
  if(is.null(datasetObject$keys)) {
    warning("datasetObject$keys is null. The probe keys are missing.")
    isThere <- FALSE
  }
  if(is.null(datasetObject$formattedName)) {
    warning("datasetObject$formattedName is null. The formatted name is missing")
  }
  
  return(isThere)
}

###############################################################################
#~ .datasetCheckRowColNames
###############################################################################
#~ .datasetCheckRowColNames makes sure that expr, pheno, and class refer to the same 
#~ samples in the same order. It also ensures that expr and keys refer to the 
#~ same probes in the same order.
#~
#~ Inputs: object, a dataset object
#~ Generates: warning messages when expr, pheno, class, and keys don't match
#~ Returns: TRUE if the samples listed in expr, pheno, and class match and the
#~    probes listed in expr and keys match as well
###############################################################################
.datasetCheckRowColNames <- function(object) {
  theyMatch <- TRUE
  if (!is.null(object$expr) && !( (all(rownames(object$expr) == names(object$keys))) && (nrow(object$expr) == length(object$keys))  )) {
    warning("The probes listed in expr and the probes listed in keys do not match.")
    theyMatch <- FALSE
  }
  if(!all(colnames(object$expr) == rownames(object$pheno))) {
    warning("The samples listed in expr and pheno do not match.")
    theyMatch <- FALSE
  }
  
  if(length(names(object$class)) == 0 && !is.null(object$class)) {
    warning("Class vector is not named.")
    theyMatch <- FALSE
  } else {
    if(!all(names(object$class) == colnames(object$expr))) {
      warning("The sample names in class do not match the sample names in expr")
      theyMatch <- FALSE
    } 
  }
  return(theyMatch)
}


###############################################################################
#~ .checkClass
###############################################################################
#~ .checkClass checks to make sure that $class is either: 1) a numeric
#~ vector of 1's, 0's, -1's OR 2) NULL. 
#~
#~  1 --> Case
#~  0 --> Control
#~ -1 --> Subject to be removed from analysis (in MANATEE)
#~
#~ Inputs: dataset, a dataset object
#~
#~ Outputs: TRUE if $class is a numeric vector of 1's, 0's, -1's OR if 
#~          $class is NULL. FALSE if $class is something else. 
#~          Also gives warnings explaining how $class is the wrong type.
###############################################################################

.checkClass <- function(dataset) {
  isOK <- FALSE
  
  if(is.null(dataset$class)) {
    return(TRUE)
  } else if (is.vector(dataset$class) && is.numeric(dataset$class) ) {
    for (i in dataset$class) {
      
      if (i != 0 && abs(i) != 1) {
        warning("Values not equal to 1, 0, or -1 are present in the class vector.")
        return(FALSE)
      }
    }
    return(TRUE)
  } else {
    warning("Class vector is not a numeric vector.")
    return(FALSE)
  }
}

###############################################################################
#~ .datasetCheckType
###############################################################################
#~ .datasetCheckType makes sure that each item in dataset is the correct type
#~ pheno: data frame
#~ expr: data frame || matrix
#~ class: numeric vector containing 0's and 1's || NULL
#~ keys: vector 
#~ formattedName: string || NULL
#~ 
#~ Input: dataset: a dataset object aka a names list containing 
#~ $pheno, $expr, $class, $keys, and $formattedName. 
#~
#~ Outputs: TRUE if each item is the correct type. FALSE if any item is wrong.
#~          Also gives warnings to specify which item is wrong. 
###############################################################################
.datasetCheckType <- function(dataset) {
  phenoOK <- is.data.frame(dataset$pheno)
  exprOK  <- is.matrix(dataset$expr) && is.numeric(dataset$expr)
  classOK <- .checkClass(dataset)
  keysOK  <- is.vector(dataset$keys) 
  formattedNameOK <- !is.null(dataset$formattedName) && ((is.character(dataset$formattedName) && (length(dataset$formattedName) == 1)) || (dataset$formattedName == ""))
  
  if(!phenoOK) { warning("pheno must be a data frame.")}
  if(!exprOK) {warning("expr must be a numeric matrix.")}
  if(!classOK) {warning("class must be either: 1) a numeric vector of 1's and 0's OR 2) NULL")}
  if(!keysOK) {warning("keys must be a vector.")}
  if(!formattedNameOK) {warning("formattedName must be a string or an empty string.")}
  
  
  return(classOK && phenoOK && exprOK && keysOK && formattedNameOK)
}




###############################################################################
#~ .metaAnalysisCheckNull
###############################################################################
#~ Checks whether each $entry in a MetaAnalysis object is null
#~
#~ Inputs: a MetaAnaysis object
#~ Generates: A warning describing which $entry is null
#~ Returns: If no $entry is null, then returns TRUE
#~          If one or more $entry is null, then returns FALSE
###############################################################################
.metaAnalysisCheckNull <- function(metaAnalysisObject) {
  isNullOK <- TRUE
  
  metaAnalysisEntries <- c("datasetEffectSizes", "datasetEffectSizeStandardErrors", 
                           "pooledResults", "analysisDescription")
  pooledResultsEntries <- c("effectSize", "effectSizeStandardError", "effectSizePval", 
                            "effectSizeFDR", "tauSquared", "numStudies", "cochranesQ", 
                            "heterogeneityPval", "fisherStatUp", "fisherPvalUp", 
                            "fisherFDRUp", "fisherStatDown", "fisherPvalDown", 
                            "fisherFDRDown")
  
  
  
  for (i in metaAnalysisEntries) {
    if(is.null(metaAnalysisObject[[i]])) {
      isNullOK <- FALSE
      warning(paste(i, "is null."))
    }
  }
  
  for (i in pooledResultsEntries) {
    if(is.null(metaAnalysisObject$pooledResults[[i]])) {
      isNullOK <- FALSE
      warning(paste("$pooledResults$", i, " is null.", sep=""))
    }
  }
  
  
  return(isNullOK)
}


###############################################################################
#~ .metaAnalysisCheckType
###############################################################################
#~ .metaAnalysisCheckType makes sure that each entry within a MetaAnalysis 
#~ object is the correct type. 
#~
#~ Inputs: a MetaAnalysis object
#~ Generates: a warning if an entry within the object is the incorrect type
#~ Returns: TRUE if all entries within the object are the correct type
#~          FALSE if any entry within the object is the incorrect type
###############################################################################
.metaAnalysisCheckType <- function(metaAnalysisObject) {
  if (is.null(metaAnalysisObject)) {
    return(FALSE)
  }
  
  isTypeOK <- TRUE
  
  metaAnalysisEntries <- c("datasetEffectSizes", "datasetEffectSizeStandardErrors", 
                           "pooledResults", "analysisDescription")
  metaAnalysisTypes <- c("matrix", "matrix", "data frame", "string")
  
  
  pooledResultsEntries <- c("effectSize", "effectSizeStandardError", "effectSizePval", 
                            "effectSizeFDR", "tauSquared", "numStudies", "cochranesQ", 
                            "heterogeneityPval", "fisherStatUp", "fisherPvalUp", 
                            "fisherFDRUp", "fisherStatDown", "fisherPvalDown", 
                            "fisherFDRDown")
  pooledResultsTypes <- c("double", "double", "double", "double", "double", "integer", 
                          "double", "double", "double", "double", "double","double", 
                          "double", "double")
  pooledResultsRanges <- c("", "non-negative", "range[0,1]", "range[0,1]", 
                           "non-negative","non-negative","non-negative","range[0,1]",
                           "non-negative","range[0,1]","range[0,1]","non-negative", 
                           "range[0,1]","range[0,1]")
  
  for (i in 1:length(metaAnalysisEntries)) {
    if(!(.checkType(metaAnalysisObject[[metaAnalysisEntries[[i]]]], metaAnalysisTypes[[i]]))) {
      isTypeOK <- FALSE
      warning(metaAnalysisEntries[[i]], " is the wrong type.")
    }
  }
  
  for (i in 1:length(pooledResultsEntries)) {
    if(!.checkType(metaAnalysisObject$pooledResults[[pooledResultsEntries[[i]]]], pooledResultsTypes[[i]], pooledResultsRanges[[i]])) {
      isTypeOK <- FALSE
      warning("pooledResults$", pooledResultsEntries[[i]], " is the wrong type.")
    }
  }
  
  return(isTypeOK)
}


###############################################################################
#~ .metaFilterCheckNull
###############################################################################
#~ .metaFilterCheckNull takes a MetaFilter object and makes sure that all the 
#~ entries within it are present. It makes sure that nothing is null.
#~ Specifically, it looks for the entries listed in filterEntries within the
#~ function (aka posGeneNames, negGeneNames, FDRThresh, effectSizeThresh, etc)
#~
#~ Inputs: a MetaFilter object
#~ 
#~ Generates: warning messages describing what entries are null
#~ 
#~ Returns: TRUE, if no entry is NULL. FALSE, if anything is null. 
###############################################################################
.metaFilterCheckNull <- function(metaFilterObject) {
  isNullOK <- TRUE
  
  filterEntries <- c("posGeneNames", "negGeneNames", "FDRThresh", 
                     "effectSizeThresh", "numberStudiesThresh", "isLeaveOneOut", 
                     "heterogeneityPvalThresh","filterDescription", "timestamp")
  
  for (i in filterEntries) {
    if(is.null(metaFilterObject[[i]])) {
      isNullOK <- FALSE
      warning(paste(i, "is null."))
    }
  }
  
  return(isNullOK)
}


###############################################################################
#~ .metaFilterCheckType
###############################################################################
#~ .metaFilterCheckType, takes a MetaFilter object and checks to make sure that
#~ each entry within it are the correct type. 
#~ Specifically, it looks at the entries listed in filterEntries:
#~  (aka posGeneNames, negGeneNames, FDRThresh, effectSizeThresh, etc)
#~
#~ Inputs: a MetaFilter object
#~ Generates: warning messages describing which entries are the incorrect type
#~ Returns: TRUE, if each entry is the correct type. FALSE, if anything is the
#~          wrong type. 
###############################################################################
.metaFilterCheckType <- function(metaFilterObject) {
  isTypeOK <- TRUE
  
  filterEntries <- c("posGeneNames", "negGeneNames", "FDRThresh", 
                     "effectSizeThresh", "numberStudiesThresh", "isLeaveOneOut", 
                     "heterogeneityPvalThresh","filterDescription", "timestamp")
  
  filterTypes <- c("character vector", "character vector", "double", "double", "integer", "boolean", "double", "string", "timestamp")
  
  filterRanges <-c("", "", "range[0,1]", "non-negative", "non-negative", "integer" , "range[0,1]", "", "")
  
  
  for (j in 1:length(filterEntries)) {
    if(!.checkType(metaFilterObject[[filterEntries[[j]]]], filterTypes[[j]], filterRanges[[j]])) {
      warning(filterEntries[j], " is the wrong type.")
      isTypeOK <- FALSE
    }
  }
  
  return(isTypeOK)
}


###############################################################################
#~ .metaCheckNull
###############################################################################
#~ .metaCheckNull takes a Meta object and makes sure that each of its entries
#~ are not null. The current version does not use this function.
#~
#~ Inputs: metaObject: a Meta object
#~ Generates: warnings if any of the entries are null
#~ Returns: TRUE if none of the entries are null. FALSE if any are null
###############################################################################
.metaCheckNull <- function(metaObject) {
  isNullOK <- TRUE
  
  metaEntries <- c("originalData", "metaAnalysis", "leaveOneOutAnalysis", 
                   "filterResults")
  
  for (i in metaEntries) {
    if(is.null(metaObject[[i]])) {
      isNullOK <- FALSE
      warning(paste("metaObject$", i, " is null."))
    }
  }
  
  isNullOK <- isNullOK && 
    
    return(isNullOK)
}


###############################################################################
#~ .metaCheckType 
###############################################################################
#~ .metaCheckType takes a Meta object and makes sure that all of its entries
#~ are the correct type. Current method does not use this function. 
#~ 
#~ Entries are: $originalData, $metaAnalysis, $leaveOneOutAnalysis,
#~              $filterResults
#~
#~ Inputs: metaObject: a Meta object
#~ Generates: warnings if some entries are the wrong type
#~ Results: TRUE if all of the entries are the right type. 
#~          FALSE if any entries are the wrong type
###############################################################################
.metaCheckType <- function(metaObject) {
  oriDataOK <- TRUE
  for (aDataset in metaObject$originalData) {
    oriDataOK<- (oriDataOK && .datasetCheckType(aDataset) && .datasetCheckRowColNames(aDataset))
  }
  
  if (!oriDataOK) {
    warning("$originalData contains some datasets of the wrong type")
  }
  
  metaAnalyOK <- .metaAnalysisCheckType(metaObject$metaAnalysis)
  if(!metaAnalyOK) {
    warning("$metaAnalysis contains some entries of the wrong type")
  }
  
  looOK <- TRUE
  for(metaAnaly in metaObject$leaveOneOutAnalysis) {
    looOK <- (looOK && .metaAnalysisCheckType(metaAnaly))
  }
  if (!looOK) {
    warning("$leaveOneOutAnalysis contains a MetaAnalysis object of incorrect type.")
  }
  
  filterOK <- TRUE
  for (filterObj in metaObject$filterResults) {
    filterOK <- (filterOK && .metaFilterCheckType(filterObj))
  }
  if (!filterOK) {
    warning("$filterResults contains a MetaFilter object of the wrong type.")
  }
  
  return(oriDataOK && metaAnalyOK && looOK && filterOK)
  
}



###############################################################################
#~ .metaCheckAll 
###############################################################################
#~ .metaCheckAll takes in a Meta object and checks of each entry and sub-entry
#~ is not null and is the correct type
#~
#~ Inputs: metaObject, a Meta object
#~ Generates: warning messages denoting what entries are null or the wrong type
#~ Returns: TRUE if everything is not null and the correct type. 
#~          FALSE if anything is null or the wrong type. 
###############################################################################
.metaCheckAll <- function(metaObject, objectStage) {
  #Initialize all the checking booleans
  oriDataOK <- TRUE
  metaAnalyOK <- TRUE
  looOK <- TRUE
  filterOK <- TRUE
  
  
  # Get the list of things to check
  if (objectStage == "Pre-Analysis") {
    thingsToCheck = c("originalData")
  } else if (objectStage == "Pre-Filter") {
    thingsToCheck = c("originalData", "metaAnalysis", "leaveOneOutAnalysis")
  } else if (objectStage == "Post-Filter") {
    thingsToCheck = c("originalData", "metaAnalysis", "leaveOneOutAnalysis", "filterResults")
  } else {
    warning("objectStage is not a valid type")
    oriDataOK = FALSE
    thingsToCheck = c("")
  }
  
  
  # originalData
  if ("originalData" %in% thingsToCheck) { 
    if(is.null(metaObject$originalData)) {
      warning("$originalData is null.")
      oriDataOK <- FALSE
      
      
    } else {
      
      oriDataVector <- NULL
      
      for (datasetObj in names(metaObject$originalData)) {
        #oriDataOK<- (oriDataOK && .datasetCheckType(datasetObj) && .datasetCheckRowColNames(datasetObj) && .datasetCheckNull(datasetObj))
        #oriDataOK     <- (oriDataOK && .datasetCheckAll(metaObject$originalData[[datasetObj]]) )
        oriDataVector <- c(oriDataVector,.datasetCheckAll(metaObject$originalData[[datasetObj]])) 
      }
      
      for(i in 1:length(oriDataVector)){
        oriDataOK     <- (oriDataOK && oriDataVector[i])
      }
            
      if (!oriDataOK) {
        
        failed_datasets <- names(metaObject$originalData)[which(oriDataVector==FALSE)]
        #warning("$originalData contains some datasets that either are null, the wrong type, or have mismatched samples")
        warning(paste(c("$originalData contains some datasets that either are null, the wrong type, or have mismatched samples. These are the datasets: ",
                        paste(failed_datasets,collapse = " ")),
                      sep = " "))
      }
    }
  }
  
  
  #metaAnalysis
  if ("metaAnalysis" %in% thingsToCheck) {
    if (is.null(metaObject$metaAnalysis)) {
      warning("$metaAnalysis is null.")
      metaAnalyOK <- FALSE
    } else {
      metaAnalyOK <- (.metaAnalysisCheckType(metaObject$metaAnalysis) && .metaAnalysisCheckNull(metaObject$metaAnalysis))
      if(!metaAnalyOK) {
        warning("$metaAnalysis contains some entries that are null or of the wrong type")  
      }
    }
  }
  
  
  
  #leaveOneOutAnalysis
  if (("leaveOneOutAnalysis" %in% thingsToCheck) && !is.null(metaObject$leaveOneOutAnalysis)) {
    for(metaAnaly in metaObject$leaveOneOutAnalysis) {
      looOK <- (looOK && .metaAnalysisCheckType(metaAnaly) && .metaAnalysisCheckNull(metaAnaly))
    }
    if (!looOK) {
      warning("$leaveOneOutAnalysis contains MetaAnalysis object(s) that are null or of incorrect type.")
    }
  }
  
  
  #filterResults
  if ("filterResults" %in% thingsToCheck) {
    if (is.null(metaObject$filterResults)) {
      warning("$filterResults is null.")
      filterOK <- FALSE
    } else {
      for (filterObj in metaObject$filterResults) {
        filterOK <- (filterOK && .metaFilterCheckType(filterObj) && .metaFilterCheckNull(filterObj))
      }
      if (!filterOK) {
        warning("$filterResults contains a MetaFilter object that is null or of the wrong type.")
      }
    }
    
  }
 
  
  return(oriDataOK && metaAnalyOK && looOK && filterOK)
  
}


###############################################################################
#~ .checkType
###############################################################################
#~ .checkType checks if the type of an object is correct 
#~
#~ Inputs: object: the object to check
#~         type: string denoting the type it's supposed to be valid types are: 
#~                "integer", "double", "string", "data frame", "character vector", 
#~                 "boolean"
#~         range: (optional)  string denoting the range it's supposed to be 
#~                valid ranges are: "non-negative" and "range[0,1]"
#~
#~ Generates: a warning if the object is not the correct type or range
#~
#~ Returns: TRUE if the object is the correct type and range
#~          FALSE if the object is not the correct type or range
###############################################################################
.checkType <- function(object, type, range = "") {
  validTypes = c("integer", "double", "string", "data frame", "matrix", 
                 "dataframe or matrix" ,"character vector", "boolean", "timestamp")
  validRanges = c("non-negative", "range[0,1]")
  typeOK <- TRUE
  
  if (type == "integer") {
    if (!(is.numeric(object) && all((object %% 1 == 0)))) {
      typeOK <- FALSE 
      warning("'integer' object", deparse(substitute(object))," is not a numeric integer.")
    } 
    
  } else if (type == "double") {
    if(!is.numeric(object)) {
      typeOK <- FALSE
      warning("'double' object ",  deparse(substitute(object))," is not numeric.")
    }
    
  } else if (type == "string") {
    if (!( (is.character(object) && (length(object) == 1)) || (object == "") )) {
      typeOK <- FALSE
      warning("'string' object ", deparse(substitute(object)), " is not a character object of length 1.")
    }
    
  } else if (type == "data frame or matrix") {
    if (!(is.data.frame(object) || is.matrix(object))) {
      typeOK <- FALSE
      warning("'data frame' object ", deparse(substitute(object)), " is not a data frame or matrix.")
    }
    
  } else if (type == "character vector") {
    if (!is.character(object)) {
      typeOK <- FALSE
      warning("'character vector' object ", deparse(substitute(object)), " is not a character vector.")
    }
    
  } else if (type == "boolean") {
    if (object != TRUE && object != FALSE) {
      typeOK <- FALSE
      warning("'boolean' object ", deparse(substitute(object)), " is not a boolean.")
    }
    
  } else if (!any(type == validTypes)) {
    typeOK <- FALSE
    warning("Invalid object type.")
    
  } else if (type == "data frame") {
    if (!is.data.frame(object)) {
      typeOK <- FALSE
      warning("'data frame' object ", deparse(substitute(object)), " is not a data frame.")
    }
    
  } else if (type == "matrix") {
    if (!is.matrix(object)) {
      typeOK <- FALSE
      warning("'matrix' object ", deparse(substitute(object)), " is not a matrix.")
    }
    
  } else if (type == "timestamp") {
    if(!(all(class(object) == c("POSIXct", "POSIXt")))) {
      typeOK <- FALSE
      warning("'timestamp' object ", deparse(substitute(object)), " is not a POSIXlt or POSIXct." )
    }
  }
  
  
  
  if (range != "" && is.numeric(object)) {
    if (range == "non-negative" && any(stats::na.omit(object) < 0)) {
      typeOK <- FALSE
      warning(deparse(substitute(object)), " contains negative values.")
      
    } else if (range == "range[0,1]" && (any(stats::na.omit(object) > 1) || any(stats::na.omit(object) < 0))) {
      typeOK <- FALSE
      warning(deparse(substitute(object)), " is not between 0 and 1.")
    } else if (!any(range == validRanges)) {
      typeOK <- FALSE
      warning("Invalid range type for doubles.")
    }
  }
  
  
  return(typeOK)
}
