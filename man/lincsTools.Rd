% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lincsTools.R
\name{lincsTools}
\alias{lincsTools}
\title{Run Shane's LINCS Tools on MetaIntegrator}
\usage{
lincsTools(metaObject, filterObject, report.out.folder,
  hit.number.hm = 10, hit.number.tbl = 10, resize = F,
  reportTitle = "lincsReport")
}
\arguments{
\item{metaObject}{a Meta object which must have the $originalData populated}

\item{filterObject}{a MetaFilter object containing the signature genes that will be used for calculating the score}

\item{report.out.folder}{Directory where a report with all figures and tables will be generated.}

\item{hit.number.hm}{How many hits to show in a heatmap (default:10)}

\item{hit.number.tbl}{How many hits to show in a displayed table (default:10)}

\item{resize}{Whether to resize tables in the way Purvesh prefers for figures (default: FALSE)}

\item{reportTitle}{file prefix for report outputs (default: "lincsReport")}
}
\value{
LINCS report for the data
}
\description{
Run Shane's LINCS Tools on MetaIntegrator
}
\examples{
\dontrun{
####### DATA SETUP ##########
# Example won't work on tinyMetaObject because it requires real gene names
# Download the needed datasets for processing. 
sleData <- getGEOData(c("GSE11909","GSE50635", "GSE39088"))

#Label classes in the datasets
sleData$originalData$GSE50635 <- classFunction(sleData$originalData$GSE50635, 
  column = "subject type:ch1", diseaseTerms = c("Subject RBP +", "Subject RBP -"))
sleData$originalData$GSE11909_GPL96 <- classFunction(sleData$originalData$GSE11909_GPL96, 
   column = "Illness:ch1", diseaseTerms = c("SLE"))
sleData$originalData$GSE39088 <- classFunction(sleData$originalData$GSE39088, 
   column= "disease state:ch1", diseaseTerms=c("SLE"))
 #Remove the GPL97 platform that was downloaded
sleData$originalData$GSE11909_GPL97 <- NULL

#Run Meta-Analysis
sleMetaAnalysis <- runMetaAnalysis(sleData, runLeaveOneOutAnalysis = F, maxCores = 1)

#Filter genes
sleMetaAnalysis <- filterGenes(sleMetaAnalysis, isLeaveOneOut = F, 
   effectSizeThresh = 1, FDRThresh = 0.05)
####### END DATA SETUP ##########

# Run immunoStates
 lincsTools(influenzaMeta, influenzaMeta$filterResults$FDR0.05_es0_nStudies4_looaTRUE_hetero0)
}
}
