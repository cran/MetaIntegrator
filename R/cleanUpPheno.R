# Functions to clean up and unscramble pheno data frame
# Author: Erika Bongen
# contact: erika.bongen@gmail.com
# Date: 3/22/2018


# Background: 
#   - $pheno in a Dataset object is a data frame containing phenotypic
#      information about each sample
#   - The way pheno is downloaded using GEOquery is unsatisfying because: 
#        1) Includes a lot of useless columns, like the street address of the
#           lab that uploaded it and the date it was made public
#        2) If a sample is missing a value in the characteristics column, it 
#           causes the cells to shift over and the pheno columns become scrambled

# Purpose: 
#   1) Remove useless columns
#   2) Reformat characteristics columns to undo scrambling if it happened


# How to use: 
# Just put your dataset object through cleanUpPheno
# It will return myDatasetObj with:
#     myDatasetObj$rawPheno: your original pheno data frame
#     myDatasetObj$pheno: a reformatted pheno data frame
# myDatasetObj = cleanUpPheno(myDatasetObj)

# require(MetaIntegrator)

# .strsplit_emptyString returns an empty string if you
# try to split an empty string
#
# Inputs: 
#    myString: string, the string you want to split
#    mySplit: string, the pattern you're splitting around
#    whichPart: integer, whether you want the 1st, 2nd, 3rd, etc part around split
.strsplit_emptyString = function(myString, mySplit, whichPart) {
  # If is NA, then return NA
  if(is.na(myString)){
    return(NA)
  }
  
  # If empty string, return empty string
  if(myString == "") {
    return("")
  }
  
  # if not empty, split it!
  return(strsplit(myString, split=mySplit)[[1]][[whichPart]])
}


# strsplitVector 
# String splits all members of a vector the same way
#
# Inputs: 
#       myVector: the vector you want to split, aka column from pheno
#       mySplit: the character that you'll split each item in the vector by (e.g " ")
#       whichPart: integer, do you want the first, second, etc part of the split vector
.strsplitVector = function(myVector, mySplit, whichPart) {
  myVector = as.character(myVector)
  
  myNewVector = as.vector(sapply(myVector, function(x) .strsplit_emptyString(myString = x, mySplit = mySplit, whichPart = whichPart)))
  return(myNewVector)
}


# .removeUselessColumns
# Takes pheno, and removes columns that I've manually added to a list
# of useless columns 
#
# Inputs: 
#   - pheno, dataframe with column names
# Outputs: 
#   - pheno, with columns removed with names in my unwanted list
.removeUselessColumns = function(pheno){
  # List of columns I don't like
  unwantedColumns <- c("status", "submission_date", "last_pudate_date", "type", "channel_count", 
                       "molecule_ch1", "extract_protocol_ch1", "label_ch1", "label_protocol_ch1", "taxid_ch1", 
                       "hyb_protocol", "scan_protocol", "data_processing", "contact_name", "contact_email", 
                       "contact_institute", "contact_address", "contact_city", "contact_zip/postal_code", 
                       "contact_country", "supplementary_file", "data_row_count", "last_update_date", "contact_state", 
                       "contact_fax", "contact_phone", "channel_count", "contact_web_link", "contact_department", 
                       "molecule_ch2", "extract_protocol_ch2", "label_ch2", "label_protocol_ch2", "taxid_ch2",
                       "contact_laboratory")
  
  # Remove any unwanted columns
  pheno = pheno[,!(colnames(pheno) %in% unwantedColumns)]
  
  return(pheno)
}



# .extractOneSampleCharacteristics
# Goes through a pheno data frame, grabs the row that corresponds
# to the ith sample, and extracts the cells that contain ": ". These
# are characteristics_ch1 columns. It then creates a named vector (aka dictionary)
# where the values are the entry values for that sample, and the names are the category
# that describe those values. 
#
# For instance, "sex: female" would have a value of "female" with the corresponding name, "sex"
#
# Note: this method is robust to scrambling
#
# Input: 
#   - pheno: data frame of phenotypic information
#            contains characteristic_ch1 columns in the format "category: value"
#   - i: integer, the index of the row whose sample you want to extract
.extractOneSampleCharacteristics = function(pheno, i){
  # Grab the row for sample i
  oneSampleData = unlist(as.matrix(pheno[i,]))
  
  # Only grab the columns that have ": " 
  # aka characteristics columns
  oneSampleData = oneSampleData[grepl(pattern = ": ", x = oneSampleData)]
  
  # Grab the "category", e.g. from "group: SLE" it would take "group"
  category = .strsplitVector(oneSampleData, ": ", 1)
  
  # Grab the value for this sample in this category
  # e.g. from "group: SLE" it would grab "SLE
  oneSampleValues = .strsplitVector(oneSampleData, ": ", 2)
  
  # Make category the names of the values, to keep them connected
  names(oneSampleValues) = category
  
  return(oneSampleValues)
}


# .fleshOutCharacteristicNAs
# takes the output of .extractOneSampleCharacteristics aka
# a named vector containing information from characteristics column
# and uses a vector of all categories to 
#    1) make sure this sample has its info listed in the same order as the others
#    2) Add NAs whenever a certain category is missing
#
# Inputs: 
#   - mySampleData: named vector where values describe that one sample, and names
#                   are the corresponding category for each sample
#                   (e.g. "sex: female" sex is category, female is value)
#   - allCategories: character vector of all categories in this dataset
#                    (e.g. c("sex", "age", "group", "disease_severity"))
#
# Outputs:
#   - mySampleData_withNA: named vector just like mySample data except: 
#         1) If mySampleData was in the wrong order, it puts it in consistent order
#         2) If mySampleData was missing values, then it adds an NA where the missing
#            value should be
.fleshOutCharacteristicNAs = function(mySampleData, allCategories){
  mySampleData_withNA = mySampleData[allCategories]
  names(mySampleData_withNA) = allCategories
  return(mySampleData_withNA)
}

# cleanUpCharacteristicsCol
# Takes pheno and cleans it up by:
#    1) Reformatting characteristics columns to be more useful
#          - aka: "sex: female", becomes a column called "sex" containing
#                  values like "male" and "female"
#          - The way the reformatting is done is robust to scrambling
#            Values are assigned to columns according to the label they're found with
#            aka "female" is put in the "sex" column because it was found in the entry "sex: female"
#            even if that cell accidentally ended up in the "age" column
#    2) Converts obviously numeric columns to numeric
#    3) Removes the duplicate columns:
#          - Original untweaked characteristics_ch1 columns
#          - 
#
# Inputs: 
#   - pheno: a data frame of phenotypic information
#
# Outputs: 
#   - newPheno: same as pheno, except with better formatting
.cleanUpCharacteristicsCol = function(pheno){
  # Identify characteristics columns by columns that contain ": "
  hasColon = apply(pheno, 2, function(x) any(grepl(pattern = ": ", x = x)))
  
  # Check if there's any scrambled columns
  isScrambled = apply(pheno[,hasColon], 2, function(x) length(unique(.strsplitVector(x, ": ", 1))) >1)
  if(sum(isScrambled) >0){
    warning(paste(sum(isScrambled), "Characteristics columns may be scrambled.", 
                  "This often happens when some samples are missing values or if there's an annotation irregularity.",
                  "cleanUpPheno will unscramble the characteristics columns,", 
                  "but you might want to compare $rawPheno to the $pheno to find the annotation irregularity."))
  }
  
  # Extract the vector of characteristics column values for each sample
  allSampleValues = lapply(1:nrow(pheno), function(x) .extractOneSampleCharacteristics(pheno, x))
  names(allSampleValues) = rownames(pheno)
  
  # Grab all possible categories that were measured
  # Some samples may be missing some categories
  allCategories = c()
  for(i in 1:length(allSampleValues)){
    allCategories = c(allCategories, names(allSampleValues[[i]]))
  }
  allCategories = unique(allCategories)
  
  # Reformat each sample's vector into a data frame
  # Missing values are replaced with NAs
  newPheno = lapply(allSampleValues, function(x) .fleshOutCharacteristicNAs(x, allCategories))
  newPheno = as.data.frame(t(as.data.frame(newPheno)), stringsAsFactors = T)
  
  # Check for numeric columns and convert them to numeric
  for(i in 1:ncol(newPheno)){
    # Check if it's a numeric column
    myCol = as.character(stats::na.omit(newPheno[,i]))
    if (!any(is.na(suppressWarnings(as.numeric(myCol))))){
      # convert to numeric!
      newPheno[,i] = as.numeric(as.character(newPheno[,i]))
    }
  }
  
  # Identify columns that contain ":ch1" 
  # They are automatically parsed by GEOquery and will be duplicates to the columns I just made
  hasCh1 = grepl(pattern = ":ch1", x = colnames(pheno))
  
  
  # Combine newPheno with the non-characteristics columns in pheno
  # Exclude columns that have ": " (aka hasColon) because they're the
  #     raw version of the columns in newPheno
  # Exclude columns with names that include ":ch1" becuase they are the same columns
  #     as newPheno, just automatically parsed by GEOquery
  newPheno = cbind(newPheno, pheno[,!(hasColon | hasCh1)])
  
  # Double check if there's any duplicate columns
  if(length(unique(colnames(newPheno))) < length(colnames(newPheno))){
    warning("pheno contains duplicate column names. You might want to change that.")
  }
  
  return(newPheno)

}


# cleanUpPheno
# Takes a Dataset object and:
#   1) Saves raw version of pheno as myDataset$rawPheno
#   2) Removes useless columns
#   3) Reformats and unscrambles characteristics columns
#
# Inputs: 
#   - myDataset: Khatri lab Dataset object
# Outputs: 
#   - myDataset: Khatri lab Dataset object, but $pheno is changed and $rawPheno is added
#        - myDataset$rawPheno, original unchanged pheno dataframe
#        - myDataset$pheno: the cleaned up version of pheno

#' Automatic preprocessing of $pheno dataframe
#' 
#' @description
#' Takes a Dataset object and:
#   1) Saves raw version of pheno as myDataset$rawPheno, 
#   2) Removes useless columns, 
#   3) Reformats and unscrambles characteristics columns.
#'
#' @param myDataset a \code{datasetObject} that contains unprocessed \code{$pheno}
#'
#'
#' @return myDataset a \code{datasetObject} that contains processed \code{$pheno} and original unprocessed \code{$rawPheno}
#' 
#' @author Erika Bongen
#'
#' @examples
#' 
#' \dontrun{
#' # Download and automatically preprocess pheno
#' gse53195 = getGEOData("GSE53195")
#' gse53195 = gse53195$originalData$GSE53195
#' View(gse53195$pheno) # Original $pheno
#' gse53195 = cleanUpPheno(gse53195)
#' View(gse53195$rawPheno) # Original $pheno
#' View(gse53195$pheno) # Preprocessed $Pheno
#' }
#' 
#' @export
cleanUpPheno = function(myDataset){
  # Save original pheno as rawPheno
  myDataset$rawPheno = myDataset$pheno
  
  # Exclude unwanted columns
  myDataset$pheno = .removeUselessColumns(myDataset$pheno)
  
  # Clean up the characteristics column
  myDataset$pheno = .cleanUpCharacteristicsCol(myDataset$pheno) 
  
  # Check to make sure Datset object is still good!
  print(paste("Passes checkDataObject QC: ", checkDataObject(myDataset, "Dataset")))
  
  return(myDataset)
}

