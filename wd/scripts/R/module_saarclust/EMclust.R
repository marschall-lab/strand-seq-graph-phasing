#' EM (Expectation maximization) function
#'
#' This function performs basic steps of EM algorithm.
#'
#' @param counts.l A \code{list} of plus and minus alignments per genomic region.
#' @param theta.param A \code{list} of cell type probabilities in clusters for each single cell. 
#' \code {matrix)}(rows=clusters, cols=strand states)
#' @param pi.param A \code{vector} of estimated sizes of each cluster based on initial hard clustering.
#' @param num.iter Set number of iteration to EM algorithm.
#' @param logL.th Set the difference between objective function from the current and the previous interation for EM algorithm to converge.
#' @inheritParams countProb
#' @return A \code{list} of various exported results [pVal, pVal.logL, log.l, theta.param, pi.param].
#' @importFrom matrixStats logSumExp rowLogSumExps
#' @author David Porubsky, Maryam Ghareghani
#' @export

#saarclust <- function(tab.l, theta.l=NULL, pi.param=NULL, num.iter=100, raw.counts=NULL) {
EMclust <- function(counts.l, theta.param=NULL, pi.param=NULL, num.iter=100, alpha=0.1, 
                    logL.th=1, log.scale=FALSE, numCPU=4) {

  ptm <- startTimedMessage("making clusters for parallel processing")
  cl <- makeCluster(numCPU)
  doParallel::registerDoParallel(cl)
  stopTimedMessage(ptm)
  
  cell.names <- names(theta.param)

  if (num.iter>1) {
    message("Running EM") 
  } else {
    message("Rescaling theta parameter")
  }
  
  ## Log scale theta and pi parameters if log.scale=TRUE
  if (log.scale) {
    theta.param <- lapply(theta.param, log)
    pi.param <- log(pi.param)
  }
  
  #counts.l <- list()
  BN.probs.l <- list()
  log.like.l <- list()

  ## Run user defined number of EM iterations ##
  for (i in 1:num.iter) {
    ptm <- startTimedMessage("    Iteration num: ",i," ")  

    clusters.per.cell <- list()
    clust.gammas.colsums.l <- list() ## list of matrices per single cell: matrix rows=clusters, cols=strand states
    #clust.gammas.rowsums.l <- list()
    theta.update.l <- list() ## list of matrices per single cell: matrix rows=clusters, cols=strand states
    #num.removed.reads <- 0
    ## Loop over all Strand-seq cells
    
    parallel.results <- foreach (j = 1:length(counts.l), .packages=c('matrixStats'),
                                 .export=c('countProb','gammaFunction')) %dopar%{
      print(paste('single cell', j))
      ## Get theta estimates for a given cell
      #TODO: rename params to cell.params
      params <- theta.param[[j]]
    
      ## Get counts of W and C reads per read/genomic segment for a given cell
      counts <- counts.l[[j]]
  
      ## Calculate BN probabilities
#      if (i == 1) {
	#print(paste('i=', i, ',j=', j, ',countProb ...'))
      BN.probs <- countProb(minusCounts = counts[,1], plusCounts = counts[,2], alpha=alpha, log=log.scale)
#        BN.probs.l[[j]] <- BN.probs
	#print('countProb is done.')
#      } else {
#        BN.probs <- BN.probs.l[[j]]
#      }
      
      ## Remove BN probs and corresponding PB reads which probability is assigned NaN value or when all probabilities are zero [TODO!!! resolve this problem]
      if (log.scale) {
        mask <- which( apply(BN.probs, 1, function(x) any(is.nan(x) | logSumExp(x)==-Inf)) )
      } else {
        mask <- which( apply(BN.probs, 1, function(x) any(is.nan(x) | sum(x)==0)) )
      }
        
      #TODO: double check the potential confusion of names
      if (length(mask)>0) { #remove PB reads with extremely high StrandS read counts
        message(length(mask), " NAN numbers after binomial probability calculation!!!")
        #num.removed.reads <- num.removed.reads + length(mask)
        BN.probs <- BN.probs[-mask,]
        counts.l[[j]] <- counts[-mask,]
      }
        
      ## Each cluster in clusters is a matrix multiplication of BN.probs by theta parameter for a given cell ##
      ## dim(BN.probs)=N*t; dim(params)=k*t; dim(clusters[[i1]])=N*t 
      ## ==> clusters[[i1]]=BN.probs*t(params) = t(params*t(BN.probs)) (computation is faster in the second way)
      
      clusters <- list()
      ## clusters: represent list of matrices per cluster. Matrix: rows=reads/genomic segments, cols=strand states
      #TODO: rename clusters to ...
      # params: rows=clusters, cols=strand states
      for (i1 in 1:nrow(params)) {
        ## In the log scale change multiplication changes to summation
        if (log.scale) { # It adds params[i1,] to all columns of t(BN.probs)
          clusters[[i1]] <- t(t(BN.probs)+params[i1,])
        } else { # It multiplies params[i1,] to all columns of t(BN.probs)
          clusters[[i1]] <- t(t(BN.probs)*params[i1,])
        }
      }
      
      ## Take sums over strand states (cols) for every read/genomic segment (rows) of each cluster ##
      ## clusters.per.cell[[[j]]: a list of vectors per cluster. Each vector represent rowSums of each cluster (matrix)
      if (log.scale) {
        #clusters.per.cell[[j]] <- lapply(clusters, function(x) apply(x, 1, logSumExp))
        clusters.per.cell <- lapply(clusters, rowLogSumExps)
      } else {
        clusters.per.cell <- lapply(clusters, rowSums)
      }  
      
      ## Scale pi.param given the number of Strand-seq cells ##
      cellNum <- length(counts.l)
      if (log.scale) {
        pi.param.scaled <- ( pi.param*(1/cellNum) )
      } else {
        pi.param.scaled <- ( pi.param^(1/cellNum) )
      }  
      
      ## Multiply BN probs (multiplied previously by theta) with scaled pi.param ##
      ## clusters.scaled: list of matrices per cluster for a given cell
      ## Map multiplies each element of the list with correspoding element of the vector
      if (log.scale) {
        clusters.scaled <- Map("+", clusters, pi.param.scaled)
      } else {
        clusters.scaled <- Map("*", clusters, pi.param.scaled)
      }
      
      ## Calculate gamma function ##
      clust.gammas.l <- gammaFunction(clust.prob=clusters.scaled, pi.scaled=pi.param.scaled, cellNum=cellNum, log.scale=log.scale)
      
      ## Take sum over all reads/genomic segments ##
      if (log.scale) { 
        #clust.gammas.colsums <- lapply(clust.gammas.l, function(x) apply(x, 2, logSumExp))
        clust.gammas.colsums <- lapply( clust.gammas.l, function(x) rowLogSumExps(t(x)) )
      } else {  
        clust.gammas.colsums <- lapply(clust.gammas.l, colSums) ## clust.gammas.colsums: matrix with rows=clusters, cols=strand states
      }
      clust.gammas.colsums <- do.call(rbind, clust.gammas.colsums)
      
      ## Updating theta parameters
      if (log.scale) {
        #theta.update <- clust.gammas.colsums - apply(clust.gammas.colsums, 1, logSumExp)
        theta.update <- clust.gammas.colsums - rowLogSumExps(clust.gammas.colsums)
      } else {
        theta.update <- clust.gammas.colsums / rowSums(clust.gammas.colsums)
      }  
      
      #clust.gammas.colsums.l[[j]] <- clust.gammas.colsums
      #theta.update.l[[j]] <- theta.update
      
      list(clusters.per.cell=clusters.per.cell,
           clust.gammas.colsums.l=clust.gammas.colsums,
           theta.update.l=theta.update)
    } #end of the loop over all Strand-seq cells
    
    clusters.per.cell <- lapply(parallel.results, function(l) l$clusters.per.cell)
    clust.gammas.colsums.l <- lapply(parallel.results, function(l) l$clust.gammas.colsums.l)
    theta.update.l <- lapply(parallel.results, function(l) l$theta.update)
    
    ## Report number of removed reads
    #if (num.removed.reads > 0) {
    #  message("Warning: Removing ", num.removed.reads, " PacBio reads with all zero or NaN binom probability!!!", appendLF = FALSE)
    #}  
    
    ## Keep old pi and theta parameters
    #pi.param.old <- pi.param
    #theta.param.old <- theta.param
    
    ## Update pi parameter ##
    ## Take sum over all cell types
    if (log.scale) {
      pi.update <- rowLogSumExps(sapply(clust.gammas.colsums.l, rowLogSumExps)) #sums over all strand states per single cell. Results in matrix with rows=clusters, cols=cells
      
      ## Normalize pi parameter to 1 and update pi
      pi.norm.update <- pi.update - logSumExp(pi.update)
    } else {
      ## first rowSums = sums over all strand states per single cell. Results in matrix with rows=clusters, cols=cells
      ## second rowSums = sums over single cells
      ## pi.update: vector of lentgh of number of clusters. Nominator of pi updating formula
      pi.update <- rowSums(sapply(clust.gammas.colsums.l, rowSums))
      
      ## Normalize pi parameter to 1 and update pi
      pi.norm.update <- pi.update / sum(pi.update)
    }  
    pi.param <- pi.norm.update
    
    ## Assign upated theta param as a current theta parameter
    theta.param <- theta.update.l
    
    ## Preprocess data in order to perform soft clustering and calculate likelihood function ##
    #clust.prod <- list()
    soft.probs <- list()
    for (k in 1:length(pi.param)) {
      clust.cell <- lapply(clusters.per.cell, `[[`, k)
      if (log.scale) {
        ## accessing cluster (k) in clusters.per.cell object 
        clust.rowsums.prod <- Reduce("+", clust.cell)
        clust.soft.prob <- clust.rowsums.prod + pi.param[k] #clust.soft.prob: vector of cluster assignment likelihoods per reads/genomic segments
      } else {
        ## accessing cluster (k) in clusters.per.cell object 
        clust.rowsums.prod <- Reduce("*", clust.cell)
        clust.soft.prob <- clust.rowsums.prod * pi.param[k] #clust.soft.prob: vector of cluster assignment likelihoods per reads/genomic segments
      }
      #clust.prod[[k]] <- clust.rowsums.prod
      soft.probs[[k]] <- clust.soft.prob
    }
    clust.tab.update <- do.call(cbind, soft.probs) #rows=reads/genomic segments, cols=clusters
    
    ## Calculate the likelihood function ##
    if (log.scale) {
      clust.tab.sums <- apply(clust.tab.update, 1, logSumExp)
      log.like <- sum(clust.tab.sums)*(-1)
    } else {
      clust.tab.sums <- rowSums(clust.tab.update)
      log.like <- sum(log(clust.tab.sums))*(-1)
    }
    
    ## Check the difference between previous and last results of likelihood function
    if (length(log.like.l) > 0) {
      log.l.diff <- log.like.l[[length(log.like.l)]] - log.like
      message("[",log.l.diff,"]", appendLF = FALSE)
      if (!is.null(logL.th)) {
        if (log.l.diff <= logL.th) {
          stopTimedMessage(ptm)
          message("CONVERGENCE!!!")
          break
        }
      }
    }
    
    if (!is.infinite(log.like)) {
      log.like.l[[length(log.like.l)+1]] <- log.like
    } else {
      message("Message: LogL function numerical problem", appendLF = FALSE)
    }
    
  stopTimedMessage(ptm)
  } # end of EM iteration
  
  
  ## Perform soft clustering
  ## Scaling soft probabilities to 1
  if (log.scale) {
    soft.probs.tab.norm <- clust.tab.update - rowLogSumExps(clust.tab.update)
    soft.probs.tab.norm <- exp(soft.probs.tab.norm) #get non-log scale probabilities
    #convert theta.param and pi.param back to non-log scale probabilities
    theta.param <- lapply(theta.param, exp)
    pi.param <- exp(pi.param)
  } else {
    soft.probs.tab.norm <- clust.tab.update / rowSums(clust.tab.update)
  }
  
  names(theta.param) <- cell.names
  
  stopCluster(cl)

  message("DONE!!!")  
  return(list(soft.pVal=soft.probs.tab.norm, log.l=unlist(log.like.l), theta.param=theta.param, pi.param=pi.param))
}

