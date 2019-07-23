#' Correct/update gene symbols in a metaObject
#' 
#' @description The gene symbols in gene expression data are sometimes outdated or incorrect, so this function goes through
#' your metaObject and updates the symbols based on the HGNChelper package, as well as correcting some other known issues.
#' @param metaObject your metaObject
#' @return A modified version of the input metaObject with updated gene symbols for each dataset in \code{metaObject$originalData}
#' @author Aditya M. Rao
#' @examples
#' tinyMetaObject = geneSymbolCorrection(tinyMetaObject)
#' @export
geneSymbolCorrection <- function(metaObject){
  metaCheck = checkDataObject(metaObject,'Meta','Pre-Analysis')
  if(!metaCheck){
    stop("The metaObject must be correctly formatted")
  }
  
  metaObject$originalData = lapply(metaObject$originalData,function(gse){
    if(is.factor(gse$keys)){
      temp.names = names(gse$keys)
      gse$keys = as.character(gse$keys)
      names(gse$keys) = temp.names
    }
    
    #Do hard coded gene name correction
    gse$keys <- correctionHelper(gse$keys,names(gse$keys))
    
    #Do HGNC gene name correction
    gse$keys <- HGNCgeneCorrection(gse$keys)
    return(gse)
  })
  
  return(metaObject)
}



#using the HGNChelper package
HGNCgeneCorrection <- function(geneNames){
  mapping = suppressWarnings(HGNChelper::checkGeneSymbols(geneNames,unmapped.as.na=FALSE))
  mapping$Suggested.Symbol[mapping$Suggested.Symbol == "NA"] = NA
  geneNames[1:length(geneNames)] = mapping$Suggested.Symbol
  return(geneNames)
}



###-###-###-###-###-###-###-##
###   correctionHelper()   ###
###-###-###-###-###-###-###-##

#DESCRIPTION
#This function was originally used to correct an issue that comes up in Excel where the
#names of some genes (e.g. MARCH1, DEC1) are automatically converted to dates
#(e.g. 1-Mar, 1-Dec). If all you have is your list of gene symbols, then it will
#fix all of the genes that could only correspond to one date. However, there are
#some genes that will result in the same date (e.g. MARCH1 and MARC1 will both
#produce 1-Mar) so if you want to fix those, you must also provide a list of
#gene IDs that corresponds to the gene symbols. The IDs that this currently
#supports are listed below.

#NOTE: This function has been updated to primarily just change gene symbols that are
#commonly found under the wrong name (e.g. converting "Septin #" genes to "SEPT#")

#PARAMETERS
#geneNames - list of gene symbols
#gene IDs - list of gene IDs
#print - if TRUE, messages will be printed to the console

#RETURN VALUE
#A modified version of the input geneNames that has been corrected to remove dates

#REQUIRED PACKAGES: None

#This is the list of identifiers that my code currently checks for:
# 1) HGNC
# 2) Entrez Gene
# 3) Ensembl
# 4) OMIM
# 5) UniProtKB
# 6) UniGene Cluster
# 7) UniGene Representative Sequence

correctionHelper <- function(geneNames, geneIDs=NULL){
  uniqueMonthNames = c("4-Feb","Septin 1","Septin 2","Septin 3","Septin 4","Septin 5","Septin 6","Septin 7","Septin 8",
                       "Septin 9","Septin 10","Septin 11","Septin 12","Septin 13","Septin 14","Selenoprotein 15","Gcom1")

  geneSymbols = c("FEB4","SEPT1","SEPT2","SEPT3","SEPT4","SEPT5","SEPT6","SEPT7","SEPT8",
                  "SEPT9","SEPT10","SEPT11","SEPT12","SEPT13","SEPT14","SEP15","GCOM1")

  for(i in 1:length(uniqueMonthNames)){
    geneNames[geneNames==uniqueMonthNames[i]] = geneSymbols[i]
  }

  if(!is.null(geneIDs)){
    if(any(geneNames=="1-Mar", na.rm = T)){
      iter = which(geneNames=="1-Mar")
      for(i in 1:length(iter)){
        if(geneIDs[iter[i]] %in% c("26077","55016","ENSG00000145416","613331","Q8TCQ1","Hs.592804","NM_001166373")){
          geneNames[iter[i]] = "MARCH1"
        } else if(geneIDs[iter[i]] %in% c("26189","64757","ENSG00000186205","614126","Q5VT66","Hs.497816","AK092439")){
          geneNames[iter[i]] = "MARC1"
        }
      }
    }
    if(any(geneNames=="2-Mar", na.rm = T)){
      iter = which(geneNames=="2-Mar")
      for(i in 1:length(iter)){
        if(geneIDs[iter[i]] %in% c("28038","51257","ENSG00000099785","613332","Q9P0N8","Hs.631861","BC032624")){
          geneNames[iter[i]] = "MARCH2"
        } else if(geneIDs[iter[i]] %in% c("26064","54996","ENSG00000117791","614127","Q969Z3","Hs.369042","AK125512")){
          geneNames[iter[i]] = "MARC2"
        }
      }
    }
  }

  return(geneNames)
}

