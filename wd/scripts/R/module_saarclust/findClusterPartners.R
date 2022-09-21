#' Get pairs of clusters coming from the same chromosome that differs in directionality.
#'
#' This function solves problem of finding the best clustering partners using linear programing function.
#' 
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @return A \code{matrix} of pairs of clusters IDs that belong to the same chromosome.
#' @importFrom lpSolve lp.assign
#' @author David Porubsky
#' @export

findClusterPartners <- function(theta.param=NULL) {
    
	## If there is an uneven number of clusters remove the one with the most WC states
	num.clusters <- nrow(theta.param[[1]])
	rownames(theta.param[[1]]) <- 1:num.clusters

	if (num.clusters %% 2 != 0) {
		#Find cluster with WC state in majority of cells
		theta.sums <- Reduce("+", theta.param)
		remove.clust <- which.max(theta.sums[,3])
		message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
		theta.param <- lapply(theta.param, function(x) x[-remove.clust,])
	}

	## get only wc thetas
	theta.param.wc <- lapply(theta.param, function(x) x[,3])
	## cbind wc thetas for all single cells
	all.theta.param.wc <- do.call(cbind, theta.param.wc)
	## compute the pairwise distance of all clusters wc thetas
	d <- as.matrix(dist(all.theta.param.wc))
	## convert distance to a similarity measure
	d <- max(d) - d
	## set diagonal values to zero
	diag(d) <- 0
	## Find pairs of clusters with the highest similarity
	max.partners <- lpSolve::lp.assign(d, "max")
	max.partners.m <- max.partners$solution #matrix with pairs of clusters with maximal similarity
	## Extract indices of pair of clusters
	max.partners.idx <- which(max.partners.m > 0, arr.ind = TRUE)
	max.partners.idx <- max.partners.idx[max.partners.idx[,1] < max.partners.idx[,2],] #remove duplicate cluster partners
	
	get.clust.idx <- as.numeric(rownames(theta.param[[1]]))
	colnames(max.partners.idx) <- c("Cluster1", "Cluster2")
	return(matrix(get.clust.idx[max.partners.idx], ncol=2))
}

    
#' Get pairs of clusters coming from the same chromosome that differs in directionality. [DEPRECATED!!!]
#'
#' This function solves the Maximum matching problem for all possible pairs of clusters and reports pairs with highest similarity.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @return A \code{vector} of pairs of clusters IDs that belong to the same chromosome.
# @importFrom maxmatching maxmatching
# @importFrom igraph graph E
#' @author David Porubsky, Maryam Ghareghani

findClusterPartners_maxMatch <- function(theta.param=NULL) {
  
  ## Helper function
  ## Calculate euclidean distance for pair of datapoints
  euc.dist.v <- function(v) sqrt(sum((v[1] - v[2]) ^ 2))
  
  ## If there is an uneven number of clusters remove the one with the most WC states
  num.clusters <- nrow(theta.param[[1]])
  if (num.clusters %% 2 != 0) {
    #Find cluster with WC state in majority of cells
    theta.sums <- Reduce("+", theta.param)
    remove.clust <- which.max(theta.sums[,3])
    message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
    theta.param <- lapply(theta.param, function(x) x[-remove.clust,])
  }
  
  pairwise.dist <- list()
  for (i in 1:length(theta.param)) {
    cell.theta <- theta.param[[i]]
    pairs <- t(combn(nrow(cell.theta), 2))
    pairs.m <- cbind( cell.theta[pairs[,1],1], cell.theta[pairs[,2],2] )
    dist <- apply(pairs.m, 1, euc.dist.v)
    pairwise.dist[[i]] <- dist
  }
  pairwise.dist.m <- do.call(cbind, pairwise.dist)
  pairwise.simil.m <- max(pairwise.dist.m)-pairwise.dist.m
  simil.sum <- rowSums(pairwise.simil.m)
  G <- igraph::graph(c(rbind(pairs[,1], pairs[,2])), directed = FALSE)
  igraph::E(G)$weight <- simil.sum
  max.match <- maxmatching::maxmatching(G, weighted = TRUE)
  
  return(max.match$matching)
}
    

#' Get pairs of clusters likely coming from the same chromosome
#'
#' This function finds cluster likely coming from the same chromosome by looking for clusters that share the same
#' directionality across multiple Strand-seq libraries.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @return A \code{matrix} of pairs of clusters IDs that belong to the same chromosome.
#' @author David Porubsky
#' @export
#' 
findSplitedClusters <- function(theta.param=NULL, z.limit=2.58) {
  
  ## Helper function
  ## Calculate euclidean distance for pair of datapoints
  euc.dist.v <- function(v) sqrt(sum((v[1] - v[2]) ^ 2))
  
  ## If there is an uneven number of clusters remove the one with the most WC states
  num.clusters <- nrow(theta.param[[1]])
  if (num.clusters %% 2 != 0) {
    #Find cluster with WC state in majority of cells
    theta.sums <- Reduce("+", theta.param)
    remove.clust <- which.max(theta.sums[,3])
    message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
    theta.param <- lapply(theta.param, function(x) x[-remove.clust,])
  }
  
  pairwise.dist <- list()
  for (i in 1:length(theta.param)) {
    cell.theta <- theta.param[[i]]
    pairs <- t(combn(nrow(cell.theta), 2))
    pairs.w <- cbind( cell.theta[pairs[,1],1], cell.theta[pairs[,2],1] )
    dist <- apply(pairs.w, 1, euc.dist.v)
    pairwise.dist[[i]] <- dist
  }
  pairwise.dist.m <- do.call(cbind, pairwise.dist)
  pairwise.simil.m <- max(pairwise.dist.m)-pairwise.dist.m
  simil.sum <- rowSums(pairwise.simil.m)
  zscores <- (simil.sum - mean(simil.sum)) / sd(simil.sum)
  idx <- zscores > z.limit #user defined confidence interval
  simil.partners.idx <- cbind(pairs[idx,1], pairs[idx,2])
  colnames(simil.partners.idx) <- c('Cluster1', 'Cluster2')
  
  return(simil.partners.idx)
}


#' Get pairs of clusters likely coming from the same chromosome but differ in directionality
#'
#' This function finds cluster likely coming from the same chromosome by looking for clusters that have opposite
#' directionality across multiple Strand-seq libraries.
#'
#' @param theta.param A \code{list} of estimated cell types for each cluster and each cell.
#' @return A \code{matrix} of pairs of clusters IDs that belong to the same chromosome.
#' @author David Porubsky
#' @importFrom igraph graph clusters V
#' @export
#' 
findAntiparallelClusters <- function(theta.param=NULL, z.limit=2.58) {
  
  ## Helper function
  ## Calculate euclidean distance for pair of datapoints
  euc.dist.v <- function(v) sqrt(sum((v[1] - v[2]) ^ 2))
  
  ## If there is an uneven number of clusters remove the one with the most WC states
  num.clusters <- nrow(theta.param[[1]])
  if (num.clusters %% 2 != 0) {
    #Find cluster with WC state in majority of cells
    theta.sums <- Reduce("+", theta.param)
    remove.clust <- which.max(theta.sums[,3])
    message("    Removed cluster ", remove.clust, " to ensure even number of clusters!!!")
    theta.param <- lapply(theta.param, function(x) x[-remove.clust,])
  }
  
  pairwise.dist <- list()
  for (i in 1:length(theta.param)) {
    cell.theta <- theta.param[[i]]
    pairs <- t(combn(nrow(cell.theta), 2))
    pairs.w <- cbind( cell.theta[pairs[,1],1], cell.theta[pairs[,2],2] )
    dist <- apply(pairs.w, 1, euc.dist.v)
    pairwise.dist[[i]] <- dist
  }
  pairwise.dist.m <- do.call(cbind, pairwise.dist)
  pairwise.simil.m <- max(pairwise.dist.m)-pairwise.dist.m
  simil.sum <- rowSums(pairwise.simil.m)
  zscores <- (simil.sum - mean(simil.sum)) / sd(simil.sum)
  idx <- zscores > z.limit #user defined confidence interval
  
  vertices <- c(rbind(pairs[idx,1], pairs[idx,2]))
  G <- igraph::graph(vertices, directed = FALSE)
  cl <- igraph::clusters(G)
  
  cluster.l <- lapply(seq_along(cl$csize), function(x) V(G)[cl$membership %in% x])
  cluster.l <- cluster.l[lengths(cluster.l) > 1]
  
  #anti.partners.idx <- cbind(pairs[idx,1], pairs[idx,2])
  #colnames(anti.partners.idx) <- c('Cluster1', 'Cluster2')
  
  return(cluster.l)
}
