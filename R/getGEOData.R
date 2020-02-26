#'GEO download/processing through GEOquery
#' @export
#' @author Francesco Vallania, Andrew Tam, Ravi Shankar, Aditya M. Rao
#'
#' @description Creates MetaIntegrator formatted objects by downloading and formatting data from GEO.
#' @param gseVector a vector of GSE ids (each a string) 
#' @param formattedNames a vector of formatted names corresponding to the GSE ids. Default: gseVector
#' @param qNorm perform quantile normalization of expression data within a dataset or not. Default: FALSE
#' @param ... will pass additional parameters to getGEO, including \code{destdir}, which specifies download location
#' 
#' @details Note: if you get the error "Error: Couldn't find driver MySQL" then just library(RMySQL) and then re-run getGEOData
#' @return a Pre-Analysis MetaObject containing the datasets loaded in $originalData
#' @export

getGEOData <- function(gseVector, formattedNames=gseVector, qNorm=FALSE, ...){
  
  #bugfix - sometimes there's a MySQL error if this isn't run
  requireNamespace(package = "RMySQL", quietly=TRUE)
  
  #Correct gses to upper case
  gseVector <- toupper(gseVector)
  
  originalData        <- lapply(gseVector,.GEOqueryCreateGEM, qNorm, ...)
  names(originalData) <- gseVector
  originalData        <- .geoquery_gems_2_singlelist(originalData, formattedNames)
  
  #add sample names on class
  originalData <- lapply(originalData,
                         function(i){
                           names(i$class) <- row.names(i$pheno)
                           return(i)
                         })
  
  #return the data in a format compatible with MetaIntegrator
  return(list(originalData = originalData))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#.GEOqueryCreateGEM
#~by Andrew Tam + Francesco Vallania [July 2014]
#~this will fetch and annotate a GEM object(s) directly from GEO using the GSE ID!!! 
#~[requires GEOquery package]
#######################################################################################
.GEOqueryCreateGEM <- function(id, qNorm, ...){
  #make sure id starts with GSE or GDS
  if(length(grep("^GSE",id,invert=T))==0 & length(grep("^GDS",id,invert=T))==0){
    id <- paste("GSE",id,sep="");
  }
  
  #fetch GEO object from GSEID
  #run function with tryCatch
  #note-> picking up only one GSE now. Fix this for multiple GPLs
  GEO <- tryCatch(GEOquery::getGEO(id, ...),
                  error   = function(e) e);
  
  if(class(GEO)[[1]]=="GDS") {
    GEO <- list(a=GEOquery::GDS2eSet(GEO))
    names(GEO)[[1]] <- id
  }
  
  #declare array object
  gem <- NULL;
  
  #if object was downloaded then go for it!
  if(class(GEO[[1]])[1] == "ExpressionSet"){
    #lapply function to GEM object
    gem <- lapply(GEO,
                  function(i){
                    .GEOqueryGEO2GEM(i,id, qNorm);
                  });
  }
  
  #return dataset
  return(gem);
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GEOqueryGEO2GEM
#~by Andrew Tam + Francesco Vallania [July 2014] + Ravi Shankar [August 2018] + Aditya M. Rao [2019]
#~this will convert the GEO object to GEM for our convenience!
#~[requires GEOquery package]
#######################################################################################
.GEOqueryGEO2GEM <- function(GEO,id,qNorm){
  
  #generate pheno matrix
  pheno = Biobase::pData(GEO);
  
  #generate expression matrix
  expr  = Biobase::exprs(GEO);
  
  if(class(expr[,1])=='character'){
    #return the new numeric matrix
    exprN = matrix(as.numeric(expr),nrow = nrow(expr))
    row.names(exprN) <- row.names(expr)
    colnames(exprN)  <-  colnames(expr)
    
    #replace original expression matrix
    expr <- exprN
  }
  
  #create gem object
  gem = list(name         = id, 
             load_comment = "Loaded with GEOquery");
  
  #parse out fData
  fData_out <- .GEO_fData_key_parser(Biobase::fData(GEO))
  
  #get platform information
  platform_info <- NA
  if (length(GEO@annotation) > 0 && grepl("^GPL", GEO@annotation)) {
    platform_info <- GEO@annotation
  }else {
    if ("platform" %in% names(GEO@experimentData@other)) {
      platform_info <- GEO@experimentData@other$platform
    }
  }
  
  #update gem object
  gem$class = rep(0, nrow(pheno))
  gem$keys = fData_out$keys
  if(length(gem$keys)>0){
    gem$key_comment = fData_out$comment
  }else{
    gem$keys = NULL
    gem$key_comment = "Annotation absent"
  }
  gem$pheno = pheno
  gem$platform = platform_info
  
  #check if there is a expression matrix
  if (dim(expr)[1] > 0) {
    gem$expr = expr
    gem$exp_comment = "Data in log scale"
    
    #check for log scale
    if (.GEM_log_check(gem) == FALSE) {
      min_value <- min(gem$expr, na.rm = T)
      
      # There is possibility of negative values in processed non-normalized 
      # data from the background subtraction algorithm. For e.g., GSE103842 (Illumina data)
      # if minimum value is negative, add small offset to all values to make them all positive
      if (min_value < 0) {
        gem$expr <- gem$expr - min_value + 1
      }
      
      gem$expr <- log2(gem$expr)
      gem$exp_comment = "Data was not in log scale originally"
    }
    
    # perform quantile normalization if desired
    if (qNorm == TRUE) {
      col_names <- colnames(gem$expr)
      row_names <- rownames(gem$expr)
      gem$expr <- preprocessCore::normalize.quantiles(gem$expr)
      colnames(gem$expr) <- col_names
      rownames(gem$expr) <- row_names
    }
  }else{
    gem$exp_comment = "Expression data is missing"
  }
  
  #return dataset
  return(gem);
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GEO_fData_key_parser
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to convert the fData input matrix into usable annotation
#~[note-> GEO_fData is fData(GEO) where GEO-> output of getGEO on a GSE object, this 
#~corresponds to its SOFT file]
#######################################################################################
.GEO_fData_key_parser <- function(GEO_fData){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Set input variables
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  symbol_pos <- NULL;
  comment    <- NULL;
  keys       <- NULL;
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Check fData length
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  check_len  <- 0;
  if(dim(GEO_fData)[2]>0){check_len <- 1};
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case1-> Simple search
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(check_len==1){
    #identify position that contains
    #gene position
    #added further term of search
    symbol_pos <- which(sapply(1:ncol(GEO_fData),
                               function(i)
                                 length(grep("^GAPDH(S)*$|^ICAM1$|^Gapdh(s)*$",GEO_fData[,i]))>0));
    #redundant fields-> pick the first one
    if(length(symbol_pos)>1){
      #if there are multiple ones pick one that matches most of the probes
      symbol_pos <- symbol_pos[which.max(sapply(symbol_pos,
                                                function(i)
                                                  length(grep("^GAPDH(S)*$|^ICAM1$|^Gapdh(s)*$",
                                                              GEO_fData[,i]))))];  
      #symbol_pos <- symbol_pos[1];
    }
    
    #generate keys vector if succesfull
    if(length(symbol_pos)==1){
      keys        <- as.character(GEO_fData[[symbol_pos]]);
      names(keys) <- as.character(GEO_fData[[1]]);
      comment     <- "Typical annotation";
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case2-> //[GENE NAME]// format
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(keys) && check_len==1){ 
    #absent field
    if(length(symbol_pos)<1){
      #identify whether it follows the //[GENE NAME]// format
      symbol_pos <- which(sapply(1:ncol(GEO_fData),
                                 function(i)
                                   length(grep("\\/\\/ GAPDH(S)* \\/\\/",GEO_fData[,i]))>0));
      
      #pick first one
      if(length(symbol_pos)> 1){symbol_pos <- symbol_pos[1]};
      if(length(symbol_pos)==1){
        
        #get row [get first one is there are multiples]
        row_pos <- grep("\\/\\/ GAPDH(S)* \\/\\/",GEO_fData[,symbol_pos]);
        if(length(row_pos)>1){row_pos <- row_pos[1]};
        
        #get exact split position [get first one is there are multiples]
        split_pos <- sapply(strsplit(as.character(GEO_fData[row_pos,symbol_pos])," /+ "),
                            function(i)
                              grep("^GAPDH(S)*$",i));
        if(length(split_pos)>1){split_pos <- split_pos[1]};
        
        #now split all lines and generate new vector
        keys <- sapply(strsplit(as.character(GEO_fData[,symbol_pos])," /+ "),
                       function(i)
                         i[split_pos]);
        names(keys) <- as.character(GEO_fData[[1]]);
        comment     <- "//[GENE NAME]// format";
      }
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case3-> genes are encoded in ENSG ID [use Human for this]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if(is.null(keys) && check_len==1){
    if(length(symbol_pos)<1){
      #identify whether it is stored as ENSG
      symbol_pos <- which(sapply(1:ncol(GEO_fData),
                                 function(i)
                                   length(grep("^ENSG00000111640$",GEO_fData[,i]))>0));
      
      if(length(symbol_pos)> 1){symbol_pos <- symbol_pos[1]};
      if(length(symbol_pos)==1){
        
        #convert genes into gene names
        ens_table   <- .ensembl_ensgID_table();
        probe2table <- match(GEO_fData[,symbol_pos],ens_table$stable_id);
        
        #store data
        keys        <- ens_table$display_label[probe2table];
        names(keys) <- as.character(GEO_fData[[1]]);
        
        comment     <- "Human ENSG format";
      }
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case4-> genes are encoded in REFSEQ format [use Human for this]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if(is.null(keys) && check_len==1){
    if(length(symbol_pos)<1){
      #identify whether it is stored as REFSEQ
      symbol_pos <- which(sapply(1:ncol(GEO_fData),
                                 function(i)#bug fix-> create pattern match for correctly mapping REFSEQ
                                   length(grep("^NM_002046.*$",GEO_fData[,i]))>0))
      
      if(length(symbol_pos)> 1){symbol_pos <- symbol_pos[1]}
      if(length(symbol_pos)==1){
        
        #convert genes into gene names by matching accession numbers
        #bug fix -> removed versioning so that it is now possible to match stuff 
        ucsc_table  <- .ucsc_refseq_table()
        
        #if REFSEQ IDs have _at at the end add it to the table
        if(length(grep("_at",GEO_fData[1,symbol_pos]))>0){
          ucsc_table$name <- gsub("$","_at",ucsc_table$name)
        }
        
        probe2table <- match(gsub("\\..$","",GEO_fData[,symbol_pos]),
                             ucsc_table$name)
        
        #store data
        keys        <- ucsc_table$name2[probe2table]
        names(keys) <- as.character(GEO_fData[[1]])
        
        comment     <- "Human REFSEQ format"
      }
    }
  }  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case5-> genes are encoded in GenBank format [use Human for this]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if(is.null(keys) && check_len==1){
    if(length(symbol_pos)<1){
      #identify whether it is stored as ENSG
      symbol_pos <- which(sapply(1:ncol(GEO_fData),
                                 function(i)
                                   length(grep("^(M61854|M29696)$",GEO_fData[,i]))>0));
      
      if(length(symbol_pos)> 1){
        symbol_pos <- symbol_pos[1]
      }
      if(length(symbol_pos)==1){
        
        #bug fix-> certain probes map to multiple genes
        gb_input_table <- GEO_fData[,symbol_pos]
        id_input_table <- GEO_fData[[1]]
        
        #
        if(length(grep(",",GEO_fData[,symbol_pos]))>0){
          #create temp lists 
          gb_input_tl <- sapply(gb_input_table,       function(i) strsplit(as.character(i),",")[[1]])
          id_input_tl <- sapply(1:length(gb_input_tl),function(i) rep(id_input_table[i],length(gb_input_tl[[i]])))
          
          #replace final values
          gb_input_table <- unlist(gb_input_tl)
          id_input_table <- unlist(id_input_tl)
        }
        
        #convert genes into gene names by matching accession numbers
        ucsc_table  <- .ucsc_genbank_table()
        probe2table <- match(gsub("\\..*$","",gb_input_table),ucsc_table$acc)
        
        #store data
        keys        <- ucsc_table$name[probe2table]
        names(keys) <- as.character(id_input_table)
        
        comment     <- "Human GenBank format"
      }
    }
  } 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case6-> genes are encoded in EntrezID format [use Human for this]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if(is.null(keys) && check_len==1){
    if(length(symbol_pos)<1){
      
      #intentify the presense of EntrezID entries
      entrez_pos <- which(sapply(1:ncol(GEO_fData),
                                 function(i)
                                   length(grep("^960$",GEO_fData[,i]))>0))
      if(length(entrez_pos)==1){
        #get table
        entrezid_table <- .ensembl_entrez_connector()
        probe2table    <- match(GEO_fData[,entrez_pos],entrezid_table$EntrezID)
        
        #store data
        keys        <- entrezid_table$GeneName[probe2table]
        names(keys) <- as.character(GEO_fData[[1]])
        
        comment     <- "Human EntrezID format"
      }
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #Case7-> MTb case handler (Rohit) using ORF column
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  if(is.null(keys) && check_len==1){
    if(length(symbol_pos)<1){
      #look for ORF column
      orf_position <- which(names(GEO_fData)=='ORF')
      
      #store data in case ORF column was found
      if(length(orf_position)>0){
        keys        <- GEO_fData[,orf_position]
        names(keys) <- as.character(GEO_fData[[1]])
        comment     <- "MTb annotation with ORF"
        #remove empty crap
        keys[which(keys=="")] <- NA
        
      }
    }
  }
  
  #Case 0a-> return absence even if there is a fData file
  if(is.null(keys) && check_len==1){
    symbol_pos  <- 1;
    comment     <- "Annotation absent";
    keys        <- as.character(GEO_fData[[symbol_pos]]);
    names(keys) <- as.character(GEO_fData[[1]]);
  }
  
  #Case 0b-> return absence of fData file
  if(is.null(keys) && check_len==0){
    comment    <- "fData file not available";
  }
  
  #replace /// into a , to make it compatible
  keys <- gsub(";| /// ",",",keys)
  
  #remove probes that do not map back to genes
  if(comment!= "MTb annotation with ORF"){
    genes      <- .ensembl_ensgID_table()
    fake_genes <- which(!(unlist(sapply(keys,function(i)strsplit(i,",")[[1]][1])) 
                          %in% 
                            genes$display_label))
    keys[fake_genes] <- NA
  }
  
  #if key vector is indeed present
  if(!is.null(keys)){
    #gene_list <- '/projects/fvallania/global_files';
    #get gene list form ensembl [replace with our own]
    #ens_g   <- ensembl_gene_connector_hs76();
    #geneset <- unique(ens_g$display_label);
    
    #set non gene keys 
    #key_map <- match(geneset,keys);
    #keys[which(is.na(key_map))] <- NA;
  }
  
  #return output
  return(list(comment=comment,
              keys   = keys));
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GEM GPL annotation functions 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#.ensembl_ensgID_table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to load an ENSGid/HUGOid table [from EnsEMBL v84]
#######################################################################################
.ensembl_ensgID_table <- function(){
  #return table from data
  return(ens_ensgID_table)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#.ensembl_entrez_table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to load an Entrez/HUGOid table [from EnsEMBL v84]
#######################################################################################
.ensembl_entrez_table <- function(){
  #return table from data
  return(ens_entrez_table)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#.ucsc_genbank_table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to load an GENBANK/HUGOid table [from UCSC]
#######################################################################################
.ucsc_genbank_table <- function(){
  #return table from data
  return(ucsc_genbank_table)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#.ucsc_refseq_table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to load an REFSEQ/HUGOid table [from UCSC]
#######################################################################################
.ucsc_refseq_table <- function(){
  #return table from data
  return(ucsc_refseq_table)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ensembl_gene_connector
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to fetch a table to convert from 
#######################################################################################
.ensembl_gene_connector_hs76 <- function(){
  
  #Establish connection to EnsEMBL
  con <- DBI::dbConnect(DBI::dbDriver("MySQL"),
                        username= 'anonymous',
                        host    = 'ensembldb.ensembl.org',
                        dbname  = 'homo_sapiens_core_76_38',
                        password= '',
                        port    = 3306);
  
  #build the parts for the query string
  query  <- "select gene.stable_id,
  xref.display_label
  from xref
  inner join gene
  where xref.xref_id = gene.display_xref_id";
  
  #run query and get table
  ens_table <- DBI::dbGetQuery(con,query);
  
  #Close connection
  DBI::dbDisconnect(con);
  
  #return table
  return(ens_table);
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ensembl_entrez_connector
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to fetch a table to convert from 
#######################################################################################
.ensembl_entrez_connector <- function(){
  
  #Establish connection to EnsEMBL
  con <- DBI::dbConnect(DBI::dbDriver("MySQL"),
                        username= 'anonymous',
                        host    = 'ensembldb.ensembl.org',
                        dbname  = 'homo_sapiens_core_76_38',
                        password= '',
                        port    = 3306)
  
  #build the parts for the query string
  query  <- "select xref.dbprimary_acc,
  xref.display_label 
  from external_db 
  join xref 
  where external_db.db_name = 'EntrezGene' AND 
  xref.external_db_id = external_db.external_db_id"
  
  #run query and get table
  ens_table           <- DBI::dbGetQuery(con,query)
  colnames(ens_table) <- c("EntrezID","GeneName")
  
  #Close connection
  DBI::dbDisconnect(con)
  
  #return table
  return(ens_table)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ucsc_genbank_connector_hg19
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to fetch a table to convert from 
#######################################################################################
.ucsc_genbank_connector_hg19 <- function(){
  
  #Establish connection to UCSC
  con <- DBI::dbConnect(DBI::dbDriver("MySQL"),
                        username= 'genome',
                        host    = 'genome-mysql.cse.ucsc.edu',
                        dbname  = 'hgFixed',
                        password= '',
                        port    = 3306);
  
  #build the parts for the query string
  query  <- "select geneName.name,
  gbCdnaInfo.acc 
  from gbCdnaInfo
  inner join geneName
  where geneName.id   = gbCdnaInfo.geneName AND
  gbCdnaInfo.organism = 448                 AND
  geneName.name      != 'n/a'               AND
  geneName.name      != 'unknown'";
  
  #run query and get table
  ucsc_tab <- DBI::dbGetQuery(con,query);
  
  #Close connection
  DBI::dbDisconnect(con);
  
  #return table
  return(ucsc_tab);
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ucsc_refseq_connector_hg19
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function is used to fetch a table to convert from 
#######################################################################################
.ucsc_refseq_connector_hg19 <- function(){
  
  #Establish connection to UCSC
  con <- DBI::dbConnect(DBI::dbDriver("MySQL"),
                        username= 'genome',
                        host    = 'genome-mysql.cse.ucsc.edu',
                        dbname  = 'hg19',
                        password= '',
                        port    = 3306);
  
  #build the parts for the query string
  query  <- "select name,name2 from refGene";
  
  #run query and get table
  ucsc_tab <- DBI::dbGetQuery(con,query);
  
  #Close connection
  DBI::dbDisconnect(con);
  
  #return table
  return(ucsc_tab);
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GEM general check functions 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GEM_log_check_v1
#~[by Francesco 2014/03/31 + Ravi Shankar August 2018]
#~this function checks if GEM is in log scale and returns TRUE/FALSE 
#######################################################################################
.GEM_log_check <- function(GEM){
  
  #remove columns that have just the same value [sd of 0] [see in GSE21126]
  goodCols <- which(sapply(1:ncol(GEM$expr),
                           function(i)
                             !stats::sd(GEM$expr[,i],na.rm=T)==0))
  
  #pick the first column that is not all NAs
  ref_p <- goodCols[which(sapply(goodCols,
                                 function(i)
                                   !all(is.na(GEM$expr[,i]))))[1]]
  
  #define input vector
  inputC <- GEM$expr[,ref_p]
  
  # Note we cannot assume that if there are negative expression values, the data is log2 normalized
  # For e.g., in the case of Illumina data, there is possibility of negative values in processed non-normalized 
  # data from the background subtraction algorithm. (GSE103842)
  # Therefore, make all values positive and apply qq_norm to check if data is log2 normalized
  min_value <- min(inputC, na.rm = T);
  
  # if minimum value is negative, add small offset to all values to make them all positive
  if (min_value < 0) {
    inputC <- inputC - min_value + 1;
  }
  
  #The data should be normally distributed
  #once in log scale
  #Get qq-norm objects
  obj_real <- stats::qqnorm(GEM$expr[,ref_p],plot.it=FALSE);
  obj_exp  <- stats::qqnorm(exp(GEM$expr[,ref_p]),plot.it=FALSE);
  obj_log  <- stats::qqnorm(log(abs(GEM$expr[,ref_p])+0.00001),plot.it=FALSE);
  
  #remove infinite
  obj_exp$x[which(is.infinite(obj_exp$x))] <- NA;
  obj_exp$y[which(is.infinite(obj_exp$y))] <- NA;  
  
  #look at correlations
  cor_real <- stats::cor(obj_real$x,obj_real$y,use='pairwise.complete');
  cor_exp  <- stats::cor(obj_exp$x, obj_exp$y, use='pairwise.complete');
  cor_log  <- stats::cor(obj_log$x ,obj_log$y ,use='pairwise.complete');
  
  #compute R^2 difference and take ratio
  log_check_score <- 0
  
  if(!all(is.na(cor_exp))){
    log_check_score <- log2(abs(cor_real**2 - cor_log**2)/abs(cor_real**2 - cor_exp**2))
    
    if(length(log_check_score) > 0 && !is.nan(log_check_score)){
      #dataset is in log scale
      if(log_check_score<=0){
        return(TRUE)
      }
      #dataset is not in log scale
      else{
        return(FALSE)
      }
    }else{
      warning("log check failed - check to make sure this dataset is log2 normalized")
    }

  }else{
    
    if(cor_real > cor_log){
      #dataset is in log scale
      return(TRUE)
    }else{
      #dataset is not in log scale
      return(FALSE)
    }
    
    
  }
  
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#geoquery_gems_2_singlelist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~this function converts a geo_query gem list into a simple list of GEM objects
#~directly compatible with the rest of the meta-analysis function
#######################################################################################
.geoquery_gems_2_singlelist <- function(geo_gem_list, formattedNames){
  print(summary(geo_gem_list))
  #make a new output list 
  gem_format <- list();
  
  #get only the gems of interest
  for(i in 1:length(names(geo_gem_list))){
    if(length(geo_gem_list[[i]])>0) {
      #get inside the object
      for(j in 1:length(geo_gem_list[[i]])){
        print(names(geo_gem_list[[i]]))
        
        #extract name
        set_name <- gsub("_series_matrix.txt.gz",
                         "",
                         names(geo_gem_list[[i]])[j]);
        #replace dashes with underscores
        set_name <- gsub("-","_",set_name);
        #save gem object in list
        gem_format[[set_name]] <- .GEM_inf_gene_remover(geo_gem_list[[i]][[j]]);
        
        gem_format[[set_name]]$formattedName <- formattedNames[[i]]
        if(length(geo_gem_list[[i]])>1) {
          gem_format[[set_name]]$formattedName <- paste(gem_format[[set_name]]$formattedName, stringr::str_match_all(set_name,"(GPL[0-9]+)")[[1]][1,2][[1]])
        }
      }
    } else {
      warning(paste("Unable to correctly download",formattedNames[[i]],". Dataset will be excluded from this object."))
    }
  }
  
  #
  return(gem_format)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GEM_inf_gene_remover
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#######################################################################################
.GEM_inf_gene_remover  <- function(GEM){ 
  #remove completely samples with -Inf, Inf, or they have some buffer overflow/underflow issues
  
  # Note: Due to a deficiency in the which() function, the commented out code fails on very large datasets
  # bad_samples           <- c(which(GEM$expr==Inf),which(GEM$expr==-Inf),which(GEM$expr >= 1.741796*10**308));
  # GEM$expr[bad_samples] <- NA;
  
  if (is.null(GEM$expr))
    return(GEM)
  
  if (ncol(GEM$expr) > 0) {
    for (i in 1:ncol(GEM$expr)) {
      GEM$expr[GEM$expr[,i]==Inf] <- NA;
      GEM$expr[GEM$expr[,i]==-Inf] <- NA;
      GEM$expr[GEM$expr[,i]>=1.741796*10**308] <- NA;
    }
  }
  
  #return GEM
  return(GEM);
}


#declare global variables for variables in data.table/with notation to avoid R CMD CHECK notes
utils::globalVariables(c("ens_ensgID_table","ens_entrez_table","ucsc_genbank_table","ucsc_refseq_table"))

