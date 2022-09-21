#' Calculate probability of PacBio reads having WW, CC or WC state.
#'
#' Import counts of directional Strand-seq reads aligned to PacBio reads and calculates binomial probabilites of PacBio reads being WW, CC or WC state.
#'
#' @param minusCounts Minus (Watson) read counts aligned to PacBio reads.
#' @param plusCounts Plus (Crick) read counts aligned to PacBio reads.
#' @param alpha Estimated level of background in Strand-seq reads.
#' @return A \code{matrix} of binomial probabilities for a given counts of plus and minus reads for a single cell. (rows=reads/genomic segments, cols=strand states)
#' @importFrom matrixStats logSumExp rowLogSumExps
#' @author David Porubsky, Maryam Ghareghani
#' @export

countProb <- function(minusCounts, plusCounts, alpha=0.1, log=FALSE) {
  
  sumCounts <- minusCounts + plusCounts
  #calculate that given PB read is WW
  prob.ww <- stats::dbinom(minusCounts, size = sumCounts, prob = 1-alpha, log = log)
    
  #calculate that given PB read is CC
  prob.cc <- stats::dbinom(minusCounts, size = sumCounts, prob = alpha, log = log)
    
  #calculate that given PB read is WC
  prob.mix <- stats::dbinom(minusCounts, size = sumCounts, prob = 0.5, log = log)
  
#  assert_that(all(!is.na(prob.ww)),
#	      all(!is.na(prob.cc)),
#	      all(!is.na(prob.mix))) %>% invisible

  prob.m <- cbind(prob.ww, prob.cc, prob.mix)
    
  return(prob.m)
}

#' Calculate gamma function of EM algorithm
#'
#' Takes binomial probabilites for every cluster and update these probabilies given the assigment of Pacbio reads to clusters.
#'
#' @param clust.prob A \code{list} of binomial probabilities (multiplied by theta & pi parameter) of each PB read in every possible cluster
#' @param pi.scaled Paramter estimate of the size of each cluster scaled by cell number
#' @param cellNum Number of Strand-seq cells used for clustering.
#' @return A \code{matrix} of unormalized probalities of cell states per cluster.
#' @author David Porubsky, Maryam Ghareghani
#' @export

gammaFunction <- function(clust.prob=NULL, pi.scaled=NULL, cellNum=NULL, log.scale=FALSE) {
  
  ## Function to take postion-wise logSumExp across all matrices (Like Reduce("+", listOfMatrices))
  reduceByLogSumExp <- function(matrix.list=NULL) {
    reducedMatrix <- matrix(nrow=nrow(matrix.list[[1]]), ncol=ncol(matrix.list[[1]]), data=0)
    for (row in 1:nrow(reducedMatrix)) {
      for (col in 1:ncol(reducedMatrix)) {
        sum <- logSumExp( sapply(matrix.list, '[[', row,col) )
        reducedMatrix[row,col] <- sum
      }
    }
    return(reducedMatrix)
  }
  
  if (log.scale) {
    #matrix.sums <- reduceByLogSumExp(clust.prob) #sum over clusters (matrices) position-wise
    matrix.sums <- apply(simplify2array(clust.prob), c(1,2), logSumExp) #sum over clusters (matrices) position-wise
    #clust.sums <- apply(matrix.sums, 1, logSumExp) #sum over rows
    clust.sums <- rowLogSumExps(matrix.sums) #sum over rows
    clust.gamma.l <- list()
    #loop over all clusters
    for (i in 1:length(clust.prob)) {
      clust.p <- clust.prob[[i]]
      clust.gamma <- clust.p - clust.sums
      clust.gamma.l[[i]] <- clust.gamma
    }
    
  } else {
    
    clust.sums <- rowSums(Reduce("+", clust.prob)) #denominator of the gamma function for a given cell, a vector over all reads/genomic segments
    #tab.colsums <- list()
    #tab.rowsums <- list()
    clust.gamma.l <- list()
    #loop over all clusters
    for (i in 1:length(clust.prob)) {
      clust.p <- clust.prob[[i]]
      #clust.gamma <- ( clust.p*( pi.scaled[i] ) ) / clust.sums #NOTE possible bug!!!
      clust.gamma <- clust.p / clust.sums #clust.gamma: rows=reads/genomic segments, cols=strand states
      #clust.gamma.colsums <- colSums(clust.gamma) # take the sum over all PB reads  #NOTE possible bug!!!
      #clust.gamma.rowsums <- rowSums(clust.gamma)
      #tab.colsums[[i]] <- clust.gamma.colsums  #NOTE possible bug!!!
      #tab.rowsums[[i]] <- clust.gamma.rowsums
      clust.gamma.l[[i]] <- clust.gamma
    }
    #gammas.colsums <- do.call(rbind, tab.colsums)  #NOTE possible bug!!!
    #gammas.rowsums <- do.call(cbind, tab.rowsums)
  }  
  
  #return(gammas.colsums) #NOTE possible bug!!!
  return(clust.gamma.l)
  #return(list(gammas.colsums=gammas.colsums, gammas.rowsums=gammas.rowsums))
} 


