#' Run Shane's LINCS Tools on MetaIntegrator
#' @param metaObject a Meta object which must have the $originalData populated
#' @param filterObject a MetaFilter object containing the signature genes that will be used for calculating the score
#' @param report.out.folder Directory where a report with all figures and tables will be generated. 
#' @param hit.number.hm How many hits to show in a heatmap (default:10)
#' @param hit.number.tbl How many hits to show in a displayed table (default:10)
#' @param resize Whether to resize tables in the way Purvesh prefers for figures (default: FALSE)
#' @param reportTitle file prefix for report outputs (default: "lincsReport")
#' 
#' @import data.table
#' @return LINCS report for the data
#' @export
#' @examples 
#' \dontrun{
#' ####### DATA SETUP ##########
#' # Example won't work on tinyMetaObject because it requires real gene names
#' # Download the needed datasets for processing. 
#' sleData <- getGEOData(c("GSE11909","GSE50635", "GSE39088"))
#' 
#' #Label classes in the datasets
#' sleData$originalData$GSE50635 <- classFunction(sleData$originalData$GSE50635, 
#'   column = "subject type:ch1", diseaseTerms = c("Subject RBP +", "Subject RBP -"))
#' sleData$originalData$GSE11909_GPL96 <- classFunction(sleData$originalData$GSE11909_GPL96, 
#'    column = "Illness:ch1", diseaseTerms = c("SLE"))
#' sleData$originalData$GSE39088 <- classFunction(sleData$originalData$GSE39088, 
#'    column= "disease state:ch1", diseaseTerms=c("SLE"))
#'  #Remove the GPL97 platform that was downloaded
#' sleData$originalData$GSE11909_GPL97 <- NULL
#' 
#' #Run Meta-Analysis
#' sleMetaAnalysis <- runMetaAnalysis(sleData, runLeaveOneOutAnalysis = F, maxCores = 1)
#' 
#' #Filter genes
#' sleMetaAnalysis <- filterGenes(sleMetaAnalysis, isLeaveOneOut = F, 
#'    effectSizeThresh = 1, FDRThresh = 0.05)
#' ####### END DATA SETUP ##########
#' 
#' # Run immunoStates
#'  lincsTools(influenzaMeta, influenzaMeta$filterResults$FDR0.05_es0_nStudies4_looaTRUE_hetero0)
#' }
lincsTools <- function(metaObject, filterObject, report.out.folder,
											 hit.number.hm = 10,hit.number.tbl = 10,resize = F, reportTitle="lincsReport") {
	dataLoad()
	
	####Set up the signature for LINCS
	signatureGenes = c(filterObject$posGeneNames, filterObject$negGeneNames)
	signatureFrame = metaObject$metaAnalysis$pooledResults[signatureGenes,]
	signature = signatureFrame[,"effectSize", drop=F]
	colnames(signature) <- c("summary")
	
	return(lincsReport(signature, hit.number.hm = hit.number.hm, hit.number.tbl = hit.number.tbl,
							report.out.folder = report.out.folder, resize = resize, reportTitle = reportTitle))
}

#' Run Shane's LINCS Correlate on MetaIntegrator
#' @param metaObject a Meta object which must have the $originalData populated
#' @param filterObject a MetaFilter object containing the signature genes that will be used for calculating the score
#' @param dataset The LINCS dataset to use. One of "CP" (drugs),"SH" (shRNA),"OE" (over-expression),
#' "LIG" (ligands),"MUT" (mutants) (default: CP)
#' @param hit.number.hm How many hits to show in a heatmap (default: 20)
#' @param direction one of "reverse", "aggravate", or "absolute" (default: "reverse") 
#' for whether you want to reverse the signature, aggravate it, or just want the top absolute hits.
#' @param cor.method method to use for correlation (pearson or spearman) (default: "pearson")
#' @param drop.string lets you include a string to drop drugs that contain a regular expression.  
#' Useful for getting rid of screening hits.  
#' One useful option is "^BRD", which gets rid of all of the Broad screening hits that aren't characterized. (default: NULL)
#' @param just_clin only consider clinically relevant results (default: FALSE)
#' @param show_clin Generate a list of clinically relevant results (default: FALSE)
#' @param gene_ann whether to annotate genes (default: FALSE)
#' 
#' @import data.table
#' @return The full list of correlations as well as the dataframe with the expression of the top hits.
#' Also generates the heatmap of the top hits.
#' @export
#' @examples 
#' \dontrun{
#' ####### DATA SETUP ##########
#' # Example won't work on tinyMetaObject because it requires real gene names
#' # Download the needed datasets for processing. 
#' sleData <- getGEOData(c("GSE11909","GSE50635", "GSE39088"))
#' 
#' #Label classes in the datasets
#' sleData$originalData$GSE50635 <- classFunction(sleData$originalData$GSE50635, 
#'   column = "subject type:ch1", diseaseTerms = c("Subject RBP +", "Subject RBP -"))
#' sleData$originalData$GSE11909_GPL96 <- classFunction(sleData$originalData$GSE11909_GPL96, 
#'    column = "Illness:ch1", diseaseTerms = c("SLE"))
#' sleData$originalData$GSE39088 <- classFunction(sleData$originalData$GSE39088, 
#'    column= "disease state:ch1", diseaseTerms=c("SLE"))
#'  #Remove the GPL97 platform that was downloaded
#' sleData$originalData$GSE11909_GPL97 <- NULL
#' 
#' #Run Meta-Analysis
#' sleMetaAnalysis <- runMetaAnalysis(sleData, runLeaveOneOutAnalysis = F, maxCores = 1)
#' 
#' #Filter genes
#' sleMetaAnalysis <- filterGenes(sleMetaAnalysis, isLeaveOneOut = F, 
#'    effectSizeThresh = 1, FDRThresh = 0.05)
#' ####### END DATA SETUP ##########
#' 
#'  lincsCorrelate( metaObject = sleMetaAnalysis, filterObject = sleMetaAnalysis$filterResults[[1]], 
#'     dataset = "CP", direction = "reverse")
#' }
lincsCorrelate <- function(metaObject, filterObject, dataset="CP", 
													 hit.number.hm = 20,
									 				 direction = "reverse",cor.method = "pearson",drop.string=NULL,
									 				 just_clin = F,show_clin = F,gene_ann =F){
	dataLoad()
	
	####Set up the signature for LINCS
	signatureGenes = c(filterObject$posGeneNames, filterObject$negGeneNames)
	signatureFrame = metaObject$metaAnalysis$pooledResults[signatureGenes,]
	signature = signatureFrame[,"effectSize", drop=F]
	colnames(signature) <- c("summary")
	
	corResults <- lincsCorrelateInternal(signature=signature,dataset=dataset,hit.number.hm = hit.number.hm,
					 direction = direction,method = cor.method,drop.string=drop.string,
					 just_clin = just_clin ,show_clin = show_clin,gene_ann =gene_ann)
	lincsHeatmap(corResults$hit.expression)
	return(corResults)
}

#' Run Shane's LINCS bait-based correlation on MetaIntegrator
#' @param metaObject a Meta object which must have the $originalData populated
#' @param filterObject a MetaFilter object containing the signature genes that will be used for calculating the score
#' @param dataset The LINCS dataset to use. One of "CP" (drugs),"SH" (shRNA),"OE" (over-expression),
#' "LIG" (ligands),"MUT" (mutants) (default: CP)
#' @param baits vector containing names of the baits being used (relevant drugs, shRNAs, etc.). See example.
#' @param just_clin only consider clinically relevant results (default: FALSE)
#' @param hit.number.hm How many hits to show in a heatmap (default: 20)
#' @param hm_baits whether or not to include the baits in the heatmap (default: FALSE)
#' @param direction one of "reverse", "aggravate", or "absolute" (default: "reverse") 
#' for whether you want to reverse the signature, aggravate it, or just want the top absolute hits.
#' @param bait_type The LINCS dataset where the baits come from. 
#' One of "CP" (drugs),"SH" (shRNA),"OE" (over-expression),
#' "LIG" (ligands),"MUT" (mutants), or NULL (don't specify) (default:NULL)
#' 
#' @description 
#' LINCS Bait Corr finds perturbagens similar to a set of interest, called baits.  
#' It searches within a defined sub space of relevant genes, usually a disease signature
#' See below for an example that recreates the work we did to find the antiviral drugs
#'  
#' @import data.table
#' @return The full list of correlations as well as the dataframe with the expression of the top hits.
#' Also generates the heatmap of the top hits.
#' @export
#' @examples 
#' \dontrun{
#' ####### DATA SETUP ##########
#' # Example won't work on tinyMetaObject because it requires real gene names
#' # Download the needed datasets for processing. 
#' sleData <- getGEOData(c("GSE11909","GSE50635", "GSE39088"))
#' 
#' #Label classes in the datasets
#' sleData$originalData$GSE50635 <- classFunction(sleData$originalData$GSE50635, 
#'   column = "subject type:ch1", diseaseTerms = c("Subject RBP +", "Subject RBP -"))
#' sleData$originalData$GSE11909_GPL96 <- classFunction(sleData$originalData$GSE11909_GPL96, 
#'    column = "Illness:ch1", diseaseTerms = c("SLE"))
#' sleData$originalData$GSE39088 <- classFunction(sleData$originalData$GSE39088, 
#'    column= "disease state:ch1", diseaseTerms=c("SLE"))
#'  #Remove the GPL97 platform that was downloaded
#' sleData$originalData$GSE11909_GPL97 <- NULL
#' 
#' #Run Meta-Analysis
#' sleMetaAnalysis <- runMetaAnalysis(sleData, runLeaveOneOutAnalysis = F, maxCores = 1)
#' 
#' #Filter genes
#' sleMetaAnalysis <- filterGenes(sleMetaAnalysis, isLeaveOneOut = F, 
#'    effectSizeThresh = 1, FDRThresh = 0.05)
#' ####### END DATA SETUP ##########
#' 
#' #Note: these are note relevant baits for SLE, just examples
#'  lincsBaitCorr(metaObject = sleMetaAnalysis, filterObject = sleMetaAnalysis$filterResults[[1]], 
#'    dataset = "CP", baits = c("NICLOSAMIDE","TYRPHOSTINA9","DISULFIRAM","SU4312","RESERPINE"))
#'  }
#' 
lincsBaitCorr <- function(metaObject, filterObject, dataset="CP", 
													baits, just_clin =F, hit.number.hm = 20,
													hm_baits = T, direction = "aggravate", bait_type=NULL){
	dataLoad()
	
	####Set up the signature for LINCS
	signatureGenes = c(filterObject$posGeneNames, filterObject$negGeneNames)
	signatureFrame = metaObject$metaAnalysis$pooledResults[signatureGenes,]
	signature = signatureFrame[,"effectSize", drop=F]
	colnames(signature) <- c("summary")
	
	corResults <- lincsBaitCorrInternal(dataset=dataset, baits=baits,sub_sig = signature, 
										just_clin =just_clin, hm_num = hit.number.hm, hm_baits = hm_baits,
										direction = direction, bait_type = bait_type)
	lincsHeatmap(corResults$heatmap.df)
	return(corResults)
}

myEnv <- new.env(parent = emptyenv())

dataLoad <- function() {
	
	# If data doesn't exist already, then load it
	if(! exists("tri.ann", envir=myEnv)) {
		tempDirPath <- tempdir()
		zipFilePath <- file.path(tempDirPath,"lincsdata.zip")
		utils::download.file(url="http://khatrilab.stanford.edu/_media/metaintegrator/lincsdata.zip", destfile = zipFilePath)
		unzippedFiles <- utils::unzip(zipFilePath, exdir = tempDirPath)
		
		all.lincs.agg.gold <- NULL
		all.lincs.qc <- NULL
		load(unzippedFiles[1])
		load(unzippedFiles[2])
		lincs.ann = utils::read.table(unzippedFiles[3], sep="\t", header = T)
		tri.ann = utils::read.csv(unzippedFiles[4],sep=",")
		gene.ann = fread(unzippedFiles[5],sep="\t")
		
		lincs.ann = lincs.ann[,c("pert_iname","pubchem_cid","pert_summary")]
		lincs.ann$merge = toupper(gsub("[[:punct:]]","",gsub(" ","",lincs.ann$pert_iname)))
		lincs.ann = lincs.ann[!duplicated(lincs.ann$merge),]
		
		tri.ann = tri.ann[,c(1,2,6,10,8,7,3:5,9,11)]
		
		setkey(gene.ann,Symbol)
		gene.ann.exp = as.data.frame(gene.ann[rownames(all.lincs.agg.gold[["CP"]]),])
		gene.ann.exp = gene.ann.exp[!duplicated(gene.ann.exp$Symbol),c(2,3,5)]
		
		gene.ann.pert = as.data.frame(gene.ann[intersect(gene.ann$Symbol
																										 ,union(colnames(all.lincs.agg.gold[["SH"]])
																										 			 ,union(colnames(all.lincs.agg.gold[["OE"]])
																										 			 			 ,colnames(all.lincs.agg.gold[["LIG"]]))
																										 )),])
		gene.ann.pert = gene.ann.pert[!duplicated(gene.ann.pert$Symbol),c(2,3,5)]
		
		myEnv$gene.ann.pert <- gene.ann.pert
		myEnv$tri.ann <- tri.ann
		myEnv$lincs.ann <- lincs.ann
		myEnv$all.lincs.agg.gold <- all.lincs.agg.gold
		myEnv$all.lincs.qc <- all.lincs.qc
	}
}

###################################################
###Function 1, lincsDataAccess
###################################################
lincsDataAccess = function(dataset,genes = "",perturbagens = "",all_genes = F,all_perturbagens = F
													 ,just_clin = F
													 ,show_clin =F){
	if(length(intersect(dataset,c("CP","SH","OE","LIG","MUT"))) == 0) return("error, incorrect dataset")
	genes = intersect(rownames(myEnv$all.lincs.agg.gold[[dataset]]),genes)
	if(all_genes) genes = rownames(myEnv$all.lincs.agg.gold[[dataset]])
	if(length(genes) == 0) return("error, no genes found")
	perturbagens = intersect(colnames(myEnv$all.lincs.agg.gold[[dataset]]),perturbagens)
	if(all_perturbagens) perturbagens = colnames(myEnv$all.lincs.agg.gold[[dataset]])
	if(just_clin) perturbagens = intersect(perturbagens,myEnv$tri.ann $merge)
	if(length(perturbagens) == 0) return("error, no pertubagens found")
	
	result = myEnv$all.lincs.agg.gold[[dataset]][genes,perturbagens,drop=F]
	result = result[order(-result[,perturbagens[1]]),,drop=F]
	if(show_clin){
		result = as.data.frame(t(result))
		result$merge = rownames(result)
		result.ann = merge(result,myEnv$tri.ann ,by = "merge",all.x = T)
		result.ann = result.ann[order(result.ann$HighestPhase),]
		result = result.ann
	}
	
	return(result)
}

###################################################
###Function 2, correlate with a signature
###################################################
lincsCorrelateInternal <- function(signature,dataset,hit.number.hm = 20
													,direction = "reverse",method = "pearson",drop.string=NULL
													,just_clin = F,show_clin = F,gene_ann =F){
	if(length(intersect(dataset,c("CP","SH","OE","LIG","MUT"))) == 0) return("error, incorrect dataset")
	dataset.used = myEnv$all.lincs.agg.gold[[dataset]]
	genes = unique(intersect(rownames(dataset.used),rownames(signature)))
	if(length(genes) == 0) return("error, no overlap between signature genes and LINCS measured genes")
	
	lincs.subset = dataset.used[genes,]
	colnames(lincs.subset) = toupper(gsub("[[:punct:]]","",colnames(lincs.subset)))
	if(dataset == "CP"){
		seps = unlist(gregexpr("^X[0-9]",colnames(lincs.subset)))
		colnames(lincs.subset)[seps != -1] = gsub("^X","",colnames(lincs.subset)[seps != -1])
	}
	
	if(!is.null(drop.string))lincs.subset = lincs.subset[,!grepl(drop.string,colnames(lincs.subset))]
	
	lincs.cor.score = as.data.frame(t(stats::cor(signature[genes,1,drop=F],lincs.subset,method = method)))
	
	###Generate p vals for the correlation tests
	lincs.cor.score$p.val = apply(lincs.subset,2,function(x){
		stats::cor.test(signature[genes,1],x,method = method)$p.value
	})
	lincs.cor.score$fdr = stats::p.adjust(lincs.cor.score$p.val,method = "fdr")
	lincs.cor.score$merge = rownames(lincs.cor.score)
	if(dataset == "CP"){
		seps = unlist(gregexpr("^X[0-9]",lincs.cor.score$merge))
		lincs.cor.score[seps != -1,"merge"] = gsub("^X","",lincs.cor.score[seps != -1,"merge"])
	}
	
	###Generate information on the number of gold and total signatures
	lincs.totsig = as.data.frame(summary(myEnv$all.lincs.qc[[dataset]]$pert_iname,maxsum = 1000000))
	lincs.totsig$merge = toupper(gsub("[[:punct:]]","",rownames(lincs.totsig)))
	lincs.totsig=  stats::aggregate(. ~ merge,lincs.totsig,sum)
	lincs.cor.score = merge(lincs.cor.score,lincs.totsig,by = "merge")
	
	lincs.goldsig = as.data.frame(myEnv$all.lincs.qc[[dataset]])
	lincs.goldsig = lincs.goldsig[lincs.goldsig$is_gold == 1,]
	lincs.goldsig = as.data.frame(summary(lincs.goldsig$pert_iname,maxsum = 1000000))
	lincs.goldsig$merge = toupper(gsub("[[:punct:]]","",rownames(lincs.goldsig)))
	lincs.goldsig=  stats::aggregate(. ~ merge,lincs.goldsig,sum)
	lincs.cor.score = merge(lincs.cor.score,lincs.goldsig,by = "merge")
	
	lincs.cor.score = lincs.cor.score[!duplicated(lincs.cor.score$merge),]
	rownames(lincs.cor.score) = lincs.cor.score$merge
	lincs.cor.score = lincs.cor.score[,-1]
	
	colnames(lincs.cor.score) = c("Sig_Corr","pval","fdr","total_sigs","gold_sigs")
	
	if(just_clin) lincs.cor.score=  lincs.cor.score[intersect(myEnv$tri.ann $merge,rownames(lincs.cor.score)),,drop=F]
	if(direction == "reverse") lincs.cor.score = lincs.cor.score[order(lincs.cor.score[,1]),,drop=F]
	if(direction == "aggravate") lincs.cor.score = lincs.cor.score[order(-lincs.cor.score[,1]),,drop=F]
	if(direction == "absolute") lincs.cor.score = lincs.cor.score[order(-abs(lincs.cor.score[,1])),,drop=F]
	
	
	lincs.subset$sig = signature[genes,1]
	lincs.subset$anti.sig = -signature[genes,1]
	
	out.list = list(hit.list = lincs.cor.score,
									hit.expression = lincs.subset[,c("sig","anti.sig",rownames(lincs.cor.score)[1:hit.number.hm])])
	
	
	if(show_clin){
		lincs.cor.ann = lincs.cor.score
		lincs.cor.ann$merge = rownames(lincs.cor.score)
		lincs.cor.ann = merge(lincs.cor.ann,myEnv$tri.ann ,by = "merge")
		lincs.cor.ann$fdr = stats::p.adjust(lincs.cor.ann$pval,method = "fdr")
		rownames(lincs.cor.ann) = lincs.cor.ann$merge
		
		if(direction == "reverse") lincs.cor.ann = lincs.cor.ann[order(lincs.cor.ann$Sig_Corr),,drop=F]
		if(direction == "aggravate") lincs.cor.ann = lincs.cor.ann[order(-lincs.cor.ann$Sig_Corr),,drop=F]
		if(direction == "absolute") lincs.cor.ann = lincs.cor.ann[order(-abs(lincs.cor.ann$Sig_Corr)),,drop=F]
		
		out.list[["hit.ann"]] = lincs.cor.ann
	}
	
	if(gene_ann){
		
		lincs.cor.ann = lincs.cor.score
		lincs.cor.ann$Symbol = rownames(lincs.cor.score)
		lincs.cor.ann = merge(lincs.cor.ann,myEnv$gene.ann.pert,by = "Symbol",all.x=  T)
		rownames(lincs.cor.ann) = lincs.cor.ann$Symbol
		lincs.cor.ann = lincs.cor.ann[rownames(lincs.cor.score),]
		out.list[["hit.ann"]] = lincs.cor.ann
	}
	
	return(out.list)
	
}

###################################################
###Function 3, lincsHeatmap
###################################################
lincsHeatmap = function(hit.expression, cutoff = NULL, pert.ordered = T,gene.ordered = F
												, transpose = T, scale = F,resize = F, hm.title = NULL, out.file = NULL
												,cluster_rows = F,cluster_cols= T){
	
	hit.expression = hit.expression[order(hit.expression$sig),]
	
	###Apply cutoff... useful for getting the data to be 0 centered
	if(!is.null(cutoff)){
		hit.expression[hit.expression > cutoff] = cutoff
		hit.expression[hit.expression < -cutoff] = -cutoff
	}
	
	
	###Scale if desired
	if(scale) hit.expression = scale(hit.expression)
	
	##Order from most to least correlated to the signature
	if(pert.ordered){
		cor.order = t(stats::cor(hit.expression[,"sig"],hit.expression))
		cor.order = cor.order[order(-cor.order[,1]),,drop=F]
		hit.expression = hit.expression[,rownames(cor.order)]
	}
	
	###Transpose if desired
	if(transpose) hit.expression = as.data.frame(t(hit.expression))
	
	hm.command = "heatmap = pheatmap::pheatmap(hit.expression
	,clustering_method = 'ward.D2',clustering_distance_cols = 'manhattan'
	,clustering_distance_rows = 'manhattan'
	,color = grDevices::colorRampPalette(c('navy', 'white', 'firebrick3'))(50)"
	if(resize)hm.command = paste(hm.command,"
															 ,cellwidth = 10, cellheight = 10",sep="")
	if(!is.null(hm.title))hm.command = paste(hm.command,"
																					 ,main = hm.title",sep="")
	if(!is.null(out.file)) hm.command = paste(hm.command,"
																						,filename = '",out.file,"'",sep="")
	
	hm.command = paste(hm.command,"
										 ,cluster_rows =cluster_rows, cluster_cols =cluster_cols)",sep="")
	
	eval(parse(text = hm.command))
	
}

###################################################
###Function 4, lincsInteractors
###################################################
lincsInteractors = function(correlate.output,signature,dataset){
	
	scores = correlate.output$hit.list
	colnames(scores)[1] = "Signature_Correlation"
	scores$pert_iname = rownames(scores)
	
	signature$Gene.Symbol = rownames(signature)
	interactor.score = merge(scores,signature,by.x = "pert_iname",by.y = "Gene.Symbol")
	interactor.score$interaction.score = interactor.score$summary * interactor.score$Signature_Correlation
	
	##Remove those hard to explain genes that seem to be moving against the signature
	if(dataset == "SH")  interactor.score = interactor.score[interactor.score$interaction.score < 0,]
	if(dataset == "OE")  interactor.score = interactor.score[interactor.score$interaction.score > 0,]
	if(dataset == "LIG")  interactor.score = interactor.score[interactor.score$interaction.score > 0,]
	interactor.score$driver.score = abs(interactor.score$interaction.score) * sign(interactor.score$summary) 
	interactor.score = interactor.score[order(-interactor.score$driver.score),]
	
	interactor.score = rbind(utils::head(interactor.score), utils::tail(interactor.score))
	interactor.score = interactor.score[!duplicated(interactor.score),]
	
	dataset.all = myEnv$all.lincs.agg.gold[[dataset]]
	genes = unique(intersect(rownames(dataset.all),rownames(signature)))
	lincs.subset = dataset.all[genes,interactor.score[,"pert_iname"]]
	lincs.subset$sig = signature[genes,1]
	lincs.subset$anti.sig = -signature[genes,1]
	
	return(list(interactor.score = interactor.score,
							interactor.expression = lincs.subset))
}

###################################################
###Function 5, lincsReport
###################################################
lincsReport = function(signature,
											 hit.number.hm = 10,hit.number.tbl = 10, method = "pearson", view.tbl = T
											 ,cutoff = 2, pert.ordered = T, gene.ordered = F,transpose = T, scale = F,resize = F  ###These are the heatmap inputs
											 ,report.out.folder, reportTitle="lincsReport"){
	
	####Run the signature through the lincsCorrelateInternal function and generate a list of hits and a heatmap
	###Generate reverse and aggravate hit lists and heatmaps for drugs, shRNA, etc.
	
	cp.conn.rev = lincsCorrelateInternal(signature,dataset = "CP", hit.number.hm = hit.number.hm, direction = "reverse",method = "pearson")
	cp.conn.rev.hits = cp.conn.rev$hit.list
	colnames(cp.conn.rev.hits)[1] = "Drug_Signature_Correlation"
	cp.conn.agg = lincsCorrelateInternal(signature,dataset = "CP", hit.number.hm = hit.number.hm, direction = "aggravate",method = "pearson")
	cp.conn.agg.hits = cp.conn.agg$hit.list
	colnames(cp.conn.agg.hits)[1] = "Drug_Signature_Correlation"
	
	sh.conn.rev = lincsCorrelateInternal(signature,dataset = "SH", hit.number.hm = hit.number.hm, direction = "reverse",method = "pearson")
	sh.conn.rev.hits = sh.conn.rev$hit.list
	colnames(sh.conn.rev.hits)[1] = "shRNA_Signature_Correlation"
	sh.conn.agg = lincsCorrelateInternal(signature,dataset = "SH", hit.number.hm = hit.number.hm, direction = "aggravate",method = "pearson")
	sh.conn.agg.hits = sh.conn.agg$hit.list
	colnames(sh.conn.agg.hits)[1] = "shRNA_Signature_Correlation"
	sh.conn.rev = lincsCorrelateInternal(signature,dataset = "SH", hit.number.hm = hit.number.hm, direction = "reverse",method = "pearson")
	#   sh.interactors.hits = lincsInteractors(sh.conn.rev,dataset = "SH",signature = signature)
	
	oe.conn.rev = lincsCorrelateInternal(signature,dataset = "OE", hit.number.hm = hit.number.hm, direction = "reverse",method = "pearson")
	oe.conn.rev.hits = oe.conn.rev$hit.list
	colnames(oe.conn.rev.hits)[1] = "OverExpression_Signature_Correlation"
	oe.conn.agg = lincsCorrelateInternal(signature,dataset = "OE", hit.number.hm = hit.number.hm, direction = "aggravate",method = "pearson")
	oe.conn.agg.hits = oe.conn.agg$hit.list
	colnames(oe.conn.agg.hits)[1] = "OverExpression_Signature_Correlation"
	oe.conn.rev = lincsCorrelateInternal(signature,dataset = "OE", hit.number.hm = hit.number.hm, direction = "reverse",method = "pearson")
	#   oe.interactors.hits = lincsInteractors(oe.conn.rev,dataset = "OE",signature = signature)
	
	lig.conn.rev = lincsCorrelateInternal(signature,dataset = "LIG", hit.number.hm = hit.number.hm, direction = "reverse",method = "pearson")
	lig.conn.rev.hits = lig.conn.rev$hit.list
	colnames(lig.conn.rev.hits)[1] = "Ligand_Signature_Correlation"
	lig.conn.agg = lincsCorrelateInternal(signature,dataset = "LIG", hit.number.hm = hit.number.hm, direction = "aggravate",method = "pearson")
	lig.conn.agg.hits = lig.conn.agg$hit.list
	colnames(lig.conn.agg.hits)[1] = "Ligand_Signature_Correlation"
	#   lig.interactors.hits = lincsInteractors(lig.conn.rev,dataset = "LIG",signature = signature)
	
	mut.conn.rev = lincsCorrelateInternal(signature,dataset = "MUT", hit.number.hm = hit.number.hm, direction = "reverse",method = "pearson")
	mut.conn.rev.hits = mut.conn.rev$hit.list
	colnames(mut.conn.rev.hits)[1] = "Mutation_Signature_Correlation"
	mut.conn.agg = lincsCorrelateInternal(signature,dataset = "MUT", hit.number.hm = hit.number.hm, direction = "aggravate",method = "pearson")
	mut.conn.agg.hits = mut.conn.agg$hit.list
	colnames(mut.conn.agg.hits)[1] = "Mutation_Signature_Correlation"
	
	###Now write a report
	cat("Writing report now.")
	system(paste("mkdir ",report.out.folder,sep=""))
	
	rdataFilePath <- file.path(report.out.folder, paste(reportTitle, "Objects.RData",sep=""))
	save(cp.conn.rev, cp.conn.rev.hits, cp.conn.agg, cp.conn.agg.hits,
			 sh.conn.rev, sh.conn.agg, sh.conn.rev.hits, sh.conn.agg.hits,
			 oe.conn.rev, oe.conn.agg, oe.conn.rev.hits, oe.conn.agg.hits,
			 lig.conn.rev, lig.conn.agg, lig.conn.rev.hits, lig.conn.agg.hits, 
			 mut.conn.rev, mut.conn.agg, mut.conn.rev.hits, mut.conn.agg.hits,
			 file=rdataFilePath)
	
	rmarkdown::render(input=system.file("lincsReport.Rmd", package="MetaIntegrator"), 
						output_format = "html_document", output_file = file.path(report.out.folder, paste(reportTitle, ".html", sep="")), 
						params = list( rdataFilePath=rdataFilePath, cutoff=cutoff, pert.ordered=pert.ordered,
													 gene.ordered=gene.ordered, transpose=transpose, scale=scale, resize=resize))
}

#######
## Function 6 - lincsQcAccess

lincsQcAccess = function(dataset,perturbagens,just_gold=F){
	
	dataset = as.data.table(myEnv$all.lincs.qc[[dataset]])
	if(just_gold) dataset = dataset[dataset$is_gold == 1,]
	
	setkey(dataset,merge)
	dataset = as.data.frame(dataset[perturbagens,])
	dataset$pert_itime = as.numeric(gsub(" h","",dataset$pert_itime))
	
	dataset = dataset[order(dataset$pert_iname
													,dataset$is_gold
													,dataset$cell_id
													,dataset$pert_idose
													,dataset$pert_itime
													,decreasing = T),]
	
	
}




# LINCS Bait Corr ---------------------------------------------------------

lincsBaitCorrInternal = function(dataset, baits, sub_sig, just_clin =F, hm_num = 20,
                                 hm_baits = T, direction = "aggravate", bait_type){
  temp.data = myEnv$all.lincs.agg.gold[[dataset]]
  temp.data = temp.data[intersect(rownames(temp.data),rownames(sub_sig)),]
  if(is.null(bait_type)){
    bait.cor = as.data.frame(t(stats::cor(temp.data[,intersect(baits,colnames(temp.data))],
                                          temp.data[,setdiff(colnames(temp.data),baits)])))
  }
  
	if(!is.null(bait_type)){
		bait.exp = myEnv$all.lincs.agg.gold[[bait_type]]
		bait.exp = bait.exp[intersect(rownames(bait.exp),rownames(sub_sig)),]
		bait.exp = bait.exp[,intersect(colnames(bait.exp),baits)]
		bait.cor = as.data.frame(stats::cor(temp.data[,setdiff(colnames(temp.data),baits)],bait.exp))
	}
	
	bait.cor$mean.cor = rowMeans(bait.cor)
	bait.cor$merge = rownames(bait.cor)
	
	if(direction == "aggravate") bait.cor = bait.cor[order(bait.cor$mean.cor,decreasing=T),]
	if(direction == "reverse") bait.cor = bait.cor[order(bait.cor$mean.cor,decreasing=F),]
	if(just_clin) bait.cor = bait.cor[intersect(union(rownames(bait.cor),baits),myEnv$tri.ann $merge),]
	
	hm.set = rownames(bait.cor)[1:hm_num]
	if(hm_baits) hm.set = intersect(union(hm.set,baits),colnames(temp.data))
	bait.hm = temp.data[,hm.set]
	bait.hm$sig = sub_sig[rownames(bait.hm),"summary"]
	bait.hm$anti.sig = -sub_sig[rownames(bait.hm),"summary"]
	
	out = list(all.hits = bait.cor, heatmap.df = bait.hm)
	if(dataset == "CP"){
		bait.clin = merge(bait.cor[,c("merge","mean.cor")],myEnv$tri.ann ,by = "merge")
		bait.clin = bait.clin[,c("pert_iname",setdiff(colnames(bait.clin),"pert_iname"))]
		bait.clin$merge = NULL
		bait.clin = bait.clin[order(bait.clin$mean.cor,decreasing=T),]
		out$hits.clin = bait.clin
		
		bait.cor = merge(bait.cor,myEnv$lincs.ann [,c("merge","pubchem_cid")])
		if(direction == "aggravate") bait.cor = bait.cor[order(bait.cor$mean.cor,decreasing=T),]
		if(direction == "reverse") bait.cor = bait.cor[order(bait.cor$mean.cor,decreasing=F),]
		out$all.hits = bait.cor
	} 
	
	return(out)
	
}


#declare global variables for variables in data.table/with notation to avoid R CMD CHECK notes
utils::globalVariables(c("Symbol"))
