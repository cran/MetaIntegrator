#######################################################################################
#######################################################################################
#Foundation libraries
#######################################################################################
#######################################################################################

#######################################################################################
#######################################################################################
#Annotate Datasets
#######################################################################################
#######################################################################################

.createAnnTable <- function(gem.data) {
  annTable  <-cbind(rownames(gem.data$expr), gem.data$keys)
  colnames(annTable) <- c("probeid", "symbol")
  return(annTable)
}

#######################################################################################
#######################################################################################
#ES Support Functions
#######################################################################################
#######################################################################################

.check.metaGEM.study <- function(study){
  
  stopifnot( all( c("expr", "class", "keys") %in% names(study) ) )
  
  if( !all(levels(as.factor(study$class)) == c("0", "1")) )
    stop("study$class must be coded as 0 or 1")
  
  if( !is.character(study$keys) )
    stop("The keys must be stored as a character vector")
}

.cleanNA <- function(x) return( x[!is.na(x) & is.finite(x) ] )

.expand.df <- function (df, key.name = "keys", keys.sep = ",") {
  ## every row in df may have multiple identities
  ## these identities are defined in keys, with individual keys separated by keys.sep
  ## this function expands every row that has 1:n mapping to n x 1:1 rows
  
  keys           <- as.character(df[ ,key.name])
  skey           <- strsplit( keys, split = keys.sep)
  df             <- df[rep(1:nrow(df), sapply(skey, length)), ]
  df[, key.name] <- unlist(skey)
  return(df)
}

#######################################################################################
#######################################################################################
#Compute Effect Sizes
#######################################################################################
#######################################################################################

.effect.sizes <- function(study){
  # computes effect size for every row in an expression matrix (e.g., every probe)
  
  .check.metaGEM.study(study)
  
  summ <- t( apply( study$expr, 1, .getES, g = study$class ) )
  summ <- data.frame( summ, keys=study$keys )
  return(summ)
}

.getES <- function(v, g){
  ## function to calculate basic statistics for two-class comparison for a single gene
  ## v: vector of samples (one row from an expression matrix)
  ## g: vector of composed of 0,1 values for every sample in v; defines case or control
  ## 
  ## encoding: CASES should be designated with the value 1, and CONTROLS with the value 0 to preserve the 
  ## directionality of the equation (mean of cases - mean of controls)
  
  # check to make sure every sample in expression matrix is assigned case or control value
  if ( length(v) != length(g) ) {
    stop("number of samples in expression matrix is different from length of class vector -- exiting")
  } 
  
  x <- .cleanNA( v[ which(g==1) ] )
  y <- .cleanNA( v[ which(g==0) ] )
  
  n1 <- length(x); n2 <- length(y)
  if( n1 < 2 | n2 < 2 )
    return( c(n1=NA, m1=NA, sd1=NA,
              n2=NA, m2=NA, sd2=NA,
              diff=NA, pooled.sd=NA,
              g=NA, se.g=NA) )
  
  #compute mean difference
  m1 <- mean(x)
  m2 <- mean(y)
  diff <- m1 - m2
  
  #compute pooled standard deviation
  sd1  <- sd(x)
  sd2  <- sd(y)
  sp   <- sqrt( ( (n1-1)*sd1^2 + (n2-1)*sd2^2 )/( n1 + n2 - 2 ) )
  
  #correct for small sample size [from Cohen's D to Hedge's G]
  cf   <- 1 - 3/( 4*(n1 + n2) - 9 ) # small sample size correction
  g    <- cf * diff/sp
  se.g <- sqrt( (n1+n2)/(n1*n2) + 0.5*g^2 /(n1+n2-3.94) )
  
  return( c(n1=n1, m1=m1, sd1=sd1,
            n2=n2, m2=m2, sd2=sd2,
            diff=diff, pooled.sd=sp,
            g=g, se.g=se.g) ) 
}


#######################################################################################
#######################################################################################
#Compute SummaryES
#######################################################################################
#######################################################################################

.combine.effect.sizes <- function (list.of.effects,
                                   between.method = "random",
                                   within.method  = "fixed.iv",
                                   everything     =  TRUE){
  
  if( is.null(names(list.of.effects)) )
    names(list.of.effects) <- paste("data", 1:length(list.of.effects), sep="")
  
  #Compute single study effect size
  study.effects <- lapply(list.of.effects, function(effects) {
    
    effects      <- data.frame(effects)
    effects$keys <- as.character(effects$keys) # keys are genes
    
    ## remove probes that cannot be mapped or have insufficient observations to calculate effect size
    bad     <- which( is.na(effects$g) | is.na(effects$keys) | effects$keys=="NA" )
    effects <- effects[ setdiff(1:nrow(effects), bad), ]
    
    ## expand probes that maps to multiple keys
    effects <- .expand.df(effects)
    
    ## summarize multiple probes within a gene in a study
    effects <- .summ.eff.within(effects, option = within.method)
  })
  
  #
  tmp  <- .multimerge(study.effects)
  g    <- tmp[, paste(names(study.effects), "_g", sep = ""), drop=FALSE]
  se.g <- tmp[, paste(names(study.effects), "_se.g", sep = ""), drop=FALSE]
  
  #Compute PooledES
  pooled.estimates <- data.frame( .pool.inverseVar(g, se.g, method=between.method ) )
  
  if (everything) {
    return(list(g=g, se.g=se.g, pooled.estimates=pooled.estimates))
  } else {
    return(pooled.estimates)
  }
}

.summ.eff.within <- function(effects, option="fixed.iv"){
  ## summarizing multiple probes that map to the same UniGene within a study
  
  stopifnot( all( c("g", "se.g", "keys") %in% colnames(effects) ) )
  
  effects       <- effects[ , c("g", "se.g", "keys")]
  effects$keys  <- as.character(effects$keys)
  
  if( length( grep(",", effects$keys) ) > 0 ) stop("Multiple keys detected. Please expand geneID first")
  
  if(nrow(effects)==1){ ## hack when only effect size is available per study
    rownames(effects) <- effects$keys
    effects <- effects[ , c("g", "se.g")]
    return(effects)
  }
  
  ## Deal with singletons first ##
  singles.keys  <- names(which(table(effects$keys) == 1))
  singles.ind   <- which(effects$keys %in% singles.keys)
  out           <- effects[ singles.ind, ]
  rownames(out) <- out$keys
  out$keys      <- NULL
  
  ## Next, deal with multiple keys within a study ##
  multis       <- effects[ -singles.ind, ]
  multis$abs.z <- abs( multis$g/multis$se.g ) # used if extreme option is choosen
  
  if(nrow(multis) > 0){  # skip if no multiple ID found
    
    tmp <- split(multis, multis$keys)
    
    if (option == "fixed.iv") {
      out2 <- sapply(tmp, function(m) {
        unlist(meta.summaries(m$g, m$se.g, method = "fixed")[c("summary", "se.summary")])
      })
      out2 <- t(out2)
      colnames(out2) <- c("g", "se.g")
    }
    
    if (option == "extreme") {
      out2 <- lapply(tmp, function(mat) mat[which.max(mat$abs.z), ])
      out2 <- do.call(rbind, out2)
      out2 <- out2[, c("g", "se.g")]
    }
    out <- rbind(out, out2)
  }
  
  out <- out[sort(rownames(out)), ]
  return(out)
}

.pool.inverseVar <- function( g, se.g, method ){
  
  ## Modified to include various measures of heterogeneity in effect sizes
  ## 1. Q.het - Cochrane's Q
  ## 2. df.het - degrees of freedom (number of studies a gene is measured in)
  ## 3. pval.het - is heterogeneity significant?
  ## in the return value, n and df should have the same values.
  
  stopifnot( identical( rownames(g), rownames(se.g) ) )
  out <- matrix( nrow=nrow(g), ncol=8,
                 dimnames=list( rownames(g), c("n.studies", "summary", "se.summary", "tau2", "p.value", "Q", "df", "pval.het") ) )
  
  for(j in 1:nrow(g)){
    
    e  <- .cleanNA(    g[j, ] )
    se <- .cleanNA( se.g[j, ] )
    n  <- length(e)
    
    if(n==1){
      summ <- e
      se.summ <- se   
      tau2 <- NA
      Q.het = NA
      df.het = NA
      pval.het = NA
    } else {
      fit      <- meta.summaries(e, se, method = method)
      summ     <- fit$summary
      se.summ  <- fit$se.summary
      tau2     <- ifelse( method=="fixed", NA, fit$tau2 )
      Q.het    <- fit$het[1]
      df.het   <- fit$het[2]
      pval.het <- fit$het[3]
      rm(fit)
    }
    
    pval     <- 2*pnorm( abs(summ/se.summ), lower.tail=FALSE )
    
    out[j, ] <- c(n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
    rm(e, se, n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
  }
  return(out)
}


#######################################################################################
#######################################################################################
#Fisher's functions
#######################################################################################
#######################################################################################

.ttest.Qvalues <- function(study){
  
  out        <- .ttest.Pvalues(study)
  out$Q.up   <- p.adjust( out$P.up, method="fdr" )
  out$Q.down <- p.adjust( out$P.down, method="fdr" )
  
  out <- out[ , c("Q.up", "Q.down", "keys")]
  return(out)
}

.ttest.Pvalues <- function(study){
  
  .check.metaGEM.study(study)
  
  summ <- .get.ttest.P( study$expr, study$class )[ , c("P.up", "P.down")]
  summ <- data.frame( summ, keys=study$keys )
  return(summ)
}

.get.ttest.P <- function(mat, g){
  ## test statistic and DF calculated using equal variance assumption
  tstat <- mt.teststat( mat, g, test="t.equalvar" )
  df <- length(g) - 2
  
  P.both <- 2*pt( abs(tstat), df=df, lower.tail=FALSE )
  P.down <- pt( tstat, df=df, lower.tail=TRUE )
  P.up   <- pt( tstat, df=df, lower.tail=FALSE )
  
  out <- cbind(P.both, P.down, P.up)
  rownames(out) <- rownames(mat)
  return(out)
}

.sum.of.logs <- function(list.of.sigs){
  
  combsigs  <- .combine.significances( list.of.sigs )
  
  sigs.up   <- combsigs$sigs.up
  valid.up  <- rowSums( !is.na(sigs.up) )
  F.stat.up <- -2*rowSums( log(sigs.up), na.rm=TRUE )
  F.pval.up <- pchisq( F.stat.up, 2*valid.up, lower.tail=FALSE )
  
  sigs.down   <- combsigs$sigs.down
  valid.down  <- rowSums( !is.na(sigs.down) )
  F.stat.down <- -2*rowSums( log(sigs.down), na.rm=TRUE )
  F.pval.down <- pchisq( F.stat.down, 2*valid.down, lower.tail=FALSE )
  
  out <- cbind(F.stat.up, F.pval.up, F.stat.down, F.pval.down)
  return(out)
}

.combine.significances <- function(list.of.sigs){
  
  study.sigs  <- lapply(list.of.sigs, function(sigs){
    sigs      <- data.frame(sigs)
    sigs$keys <- as.character(sigs$keys)
    
    ## remove probes that cannot be mapped or have insufficient observations to calculate the p-value
    bad  <- which( is.na(sigs[ ,1]) | is.na(sigs$keys) | sigs$keys=="NA" )
    sigs <- sigs[ setdiff(1:nrow(sigs), bad), ]
    
    ## expand probes that maps to multiple keys
    sigs <- .expand.df( sigs )
    
    ## summarize multiple probes within a study
    out  <- .summ.sigs.within(sigs)
  })
  
  tmp       <- .multimerge(study.sigs)
  sigs.up   <- tmp[ , grep("\\.up$", colnames(tmp), value=TRUE), drop=FALSE]
  sigs.down <- tmp[ , grep("\\.down$", colnames(tmp), value=TRUE), drop=FALSE]
  
  return(list(sigs.up=sigs.up, sigs.down=sigs.down))
}

.summ.sigs.within <- function(sigs){
  
  pval.cols <- setdiff( 1:ncol(sigs), grep("keys", colnames(sigs)) )
  out <- sapply( sigs[ , pval.cols], function(x)
    tapply(x, sigs$keys, min) )
  return(out)
}

.multimerge <- function(mylist){
  
  ## iterative version of merge using all=T option that converts lists into a single matrix
  ## assumes rownames is unique and is the common key
  
  unames <- unique( unlist( lapply( mylist, rownames ) ) )
  n      <- length(unames)
  
  out <- lapply( mylist, function(df){
    tmp <- matrix( nrow=n, ncol=ncol(df),
                   dimnames=list( unames, colnames(df) ) )
    tmp[ rownames(df), ] <- as.matrix(df)
    return(tmp)
  })
  
  bigout <- do.call( cbind, out )
  colnames(bigout) <- paste(rep( names(mylist), sapply(mylist, ncol) ),
                            sapply(mylist, colnames), sep="_")
  return(bigout)
}
