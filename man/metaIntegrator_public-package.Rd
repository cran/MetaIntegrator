\name{MetaIntegrator}
\alias{MetaIntegrator}
\docType{package}
\title{
MetaIntegrator package for meta-analysis of gene expression data 
}
\description{
The package comprises several analysis and plot functions to perform integrated multi-cohort analysis of gene expression data (meta-analysis). 
\tabular{ll}{
Package: \tab metaIntegrator_public\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-25\cr
License: \tab LGPL \cr
}
For detailed documentation of functions and use cases read: \code{vignette(MetaIntegrator)}.
}

\details{
The advent of the gene expression microarray has allowed for a rapid increase in gene expression studies. There is now a wealth of publicly available gene expression data available for re-analysis. An obvious next step to increase statistical power in detecting changes in gene expression associated with some condition is to aggregate data from multiple studies. 

The MetaIntegrator package will perform a DerSimonian & Laird random-effects meta-analysis for each gene (not probeset) between all target studies between cases and controls; it also performs a Fischer's sum-of-logs method on the same data, and requires that a gene is significant by both methods. The resulting p-values are False discovery rate (FDR) corrected to q-values, and will evaluate the hypothesis of whether each gene is differentially expressed between cases and controls across all studies included in the analysis.

The resulting list of genes with significantly different expression between cases and controls can be used for multiple purposes, such as 
(1) a new diagnostic or prognostic test for the disease of interest, 
(2) a better understanding of the underlying biology, 
(3) identification of therapeutic targets, and multiple other applications. 

Our lab has already used these methods in a wide variety of diseases, including organ transplant reject, lung cancer, neurodegenerative disease, and sepsis (Khatri et al., J Exp Med 2013; Chen et al, Cancer Res 2014; Li et al., Acta Neur Comm 2014; Sweeney et al, Sci Trans Med 2015).
}

\author{
Winston A. Haynes, Francesco Vallania, Aurelie Tomczak, Timothy E. Sweeney, Erika Bongen, Purvesh Khatri

Maintainer: Winston A. Haynes <hayneswa@stanford.edu> 

}
\references{
%~~ package: MetaIntegrator paper ~~

Sweeney et al., Science Translational Medicine, 2015

Khatri P et al. J Exp. Med. 2013

}
\keyword{ MetaIntegrator }
\seealso{
vignette(MetaIntegrator)
}
\examples{
\dontrun{
#Run a meta analysis. 
#   maxCores is set to 1 for package guideline compliance. 
#		For personal purposes, leave parameter un-set.
runMetaAnalysis(tinyMetaObject, maxCores=1)

#### a standard meta-analysis would follow this work flow: ####

# make input metaObjects from individual GEO datasetObjects
metaObject = list()
metaObject$originalData <- tinyMetaObject$originalData
# make test datasetObject
datasetObject1 <- tinyMetaObject$originalData$Whole.Blood.Study.1

# run the meta-analysis
metaObject <- runMetaAnalysis(metaObject, maxCores=1)

# select significant genes (default parameter)
metaObject <- filterGenes(metaObject)

# print a meta-analysis result summary for selected genes
summarizeFilterResults(metaObject, getMostRecentFilter(metaObject))

# use selected genes to generate a violin plot 
violinPlot(metaObject$filterResults$FDR0.05_es0_nStudies1_looaTRUE_hetero0, datasetObject1,
					labelColumn = 'group')

# use selected genes to generate a ROC plot 
rocPlot(metaObject$filterResults$FDR0.05_es0_nStudies1_looaTRUE_hetero0, datasetObject1)

# generate a forest plot for a gene of interest with forestPlot(metaObject, geneName)
forestPlot(metaObject, "Gene27")
}
}
