# Function for imputing the sex of a human sample
# Author: Erika Bongen
# Date: 4/27/2018
# Contact: erika.bongen@gmail.com


#' Imputes biological sex of each sample in a Dataset object
#'
#' @param myDataset \code{datasetObject}
#' @param femGenes vector of gene symbols of genes higher expressed in females. Defaults to NULL
#' @param malGenes vector of gene symbols of genes higher expressed in males. Defaults to NULL
#'
#' @return a vector indicating whether each sample is classified as "male" or "female" 
#'
#'@details
#' Imputes the sex of each sample in a Dataset object by performing K means
#' clustering. If genes higher expressed in females (femGenes) and genes higher
#' expressed in males (malGenes) are not supplied, then clustering
#' is performed on a default set of known X-escape genes (Tukiainen et al. 2017 Nature) 
#' and Y-chromosome genes. 
#' Genes were chosen as a subset of the immune Sex Expression Signature (iSEXS) (Bongen et al. In Prep.)
#' 
#' Known X-Escape genes: 
#' "XIST","RPS4X","CD40LG","ZRSR2","EFHC2","CA5B","ZFX","EIF1AX","CA5BP1","UBA1","SYAP1","DDX3X","FUNDC1","USP9X","SMC1A","NUP62CL","NAA10"
#' 
#' Y-Chromosome genes:
#' "KDM5D","RPS4Y1","EIF1AY","USP9Y","DDX3Y","UTY","PRKY","ZFY","TMSB4Y"
#' 
#'
#' @author Erika Bongen
#'
#' @examples
#' # Add sex labels to your dataset of choice
#' \dontrun{
#' myDatasets = getGEOData(c("GSE13485","GSE17156","GSE19442"))
#' myDatasets$originalData$GSE13485$pheno$sex = imputeSex(myDatasets$originalData$GSE13485)
#' myDatasets$originalData$GSE13485$pheno$sex
#' }
#' @export

imputeSex = function(myDataset, femGenes = NULL, malGenes = NULL){
  # If femGenes not given, use known X-escapees from iSEXS as default
  if(is.null(femGenes)){
    femGenes = c("XIST","RPS4X","CD40LG","ZRSR2",
                 "EFHC2","CA5B","ZFX","EIF1AX",
                 "CA5BP1","UBA1","SYAP1","DDX3X",
                 "FUNDC1","USP9X","SMC1A","NUP62CL","NAA10")
  }
  
  # If malGenes not give, use Y-chromosome genes from iSEXS from iSEXS as default
  if(is.null(malGenes)){
    malGenes = c("KDM5D","RPS4Y1","EIF1AY","USP9Y",
                 "DDX3Y","UTY","PRKY","ZFY","TMSB4Y")
  }
  
  # Make sure at least one of the sex annotation genes exists in this dataset
  if (!any(c(femGenes, malGenes) %in% myDataset$keys)){
    warning("Sex classification genes not present in dataset")
    return(NULL)
  }
  
  # Get sex genes in dataset
  femGenes = femGenes[which(femGenes %in% myDataset$keys)]
  malGenes = malGenes[which(malGenes %in% myDataset$keys)]
  sexGenes = c(femGenes, malGenes)
  sexExpr = getSampleLevelGeneData(myDataset, sexGenes)
  
  # Run kmeans
  kmeansResults = stats::kmeans(t(sexExpr), centers = 2)
  
  
  # Figure out which cluster is female
  ### 1 ### If both male and female genes measured: 
  if(length(femGenes) > 0 & length(malGenes)>0){
    # Get the cenetroid expression of female genes
    femCenter = kmeansResults$centers[,femGenes]
    if(!is.null(ncol(femCenter))){
      femCenter = rowMeans(femCenter)
    }
    
    # Get the centroid expression of male genes
    malCenter = kmeansResults$centers[,malGenes]
    if(!is.null(ncol(malCenter))){
      malCenter = rowMeans(malCenter)
    }
    
    # Get the centroid with high female expression and low male expression
    highFemCentroid = names(femCenter)[which(femCenter == max(femCenter))]
    lowMalCentroid = names(malCenter)[which(malCenter == min(malCenter))]
    
    # Make sure it agrees
    if(highFemCentroid != lowMalCentroid){
      warning("Strange centroid assignment in kmeans. Try different genes")
      return(NULL)
    }
    
    # Assign the femaleCentroid
    femaleCentroid = highFemCentroid
    
  }
  
  ### 2 ### if there are no male genes measured in the dataset
  if(length(femGenes) ==0 ){
    # Get the cenetroid expression of male genes
    malCenter = kmeansResults$centers[,malGenes]
    if(!is.null(ncol(malCenter))){
      malCenter = rowMeans(malCenter)
    }
    
    # Assign the female centroid
    femaleCentroid = names(malCenter)[which(malCenter == min(malCenter))]
  }
  
  ### 3 ### if there are no female genes measured in the dataset
  if(length(malGenes) ==0 ){
    # Get the cenetroid expression of female genes
    femCenter = kmeansResults$centers[,femGenes]
    if(!is.null(ncol(femCenter))){
      femCenter = rowMeans(femCenter)
    }
    
    # Assign the female centroid
    femaleCentroid = names(femCenter)[which(femCenter == max(femCenter))]
  }
  
  
  # Create imputed sex
  imputedSex = ifelse(kmeansResults$cluster == as.numeric(femaleCentroid), yes = "female", no = "male")
  
  return(imputedSex)
}
