#' Calculate probability of PacBio reads having WW, CC or WC state.
#'
#' Import counts of directional Strand-seq reads aligned to PacBio reads and calculates binomial probabilites of PacBio reads being WW, CC or WC state.
#'
#' @param minusCounts Minus (Watson) read counts aligned to PacBio reads.
#' @param plusCounts Plus (Crick) read counts aligned to PacBio reads.
#' @param alpha Estimated level of background Strand-seq reads.
#' @return A \code{list} of binomial probabilities for a given counts of plus and minus reads
#' @author David Porubsky, Maryam Ghareghani

countProbOld <- function(minusCounts, plusCounts, alpha=0.1) {

  minusCounts <- as.numeric(minusCounts)
  plusCounts <- as.numeric(plusCounts)

  #calculate that given PB read is WW
  prob.w <- (1-alpha)^minusCounts
  prob.c <- alpha^plusCounts 
  prob.ww <- prob.w * prob.c

  #calculate that given PB read is CC
  prob.w <- alpha^minusCounts 
  prob.c <- (1-alpha)^plusCounts
  prob.cc <- prob.w * prob.c

  #calculate that given PB read is WC
  prob.w <- 0.5^minusCounts
  prob.c <- 0.5^plusCounts
  prob.wc <- prob.w * prob.c

  prob.mix <- prob.wc

  prob.m <- choose(n=minusCounts+plusCounts, k=plusCounts) * cbind(prob.ww, prob.cc, prob.mix)

  return(prob.m)
}


#' Import old output from the minimap
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param infile A path to the minimap file to load.
#' @return A \code{data.frame}.
#' @author David Porubsky

importOldTestData <- function(infile=NULL, removeDuplicates = TRUE) {  #TODO modify this function for input where genomic location of PB reads is unknown
  
  ptm <- startTimedMessage("Reading the data ...") 
  #data <- read.table(infile, header=F) #TODO test data.table package for faster data import
  
  #filetype = summary( file(infile) )$class #If it's gzipped, filetype will be gzfile
  if (summary( file(infile) )$class == 'gzfile') {
    data <- data.table::fread(paste0('gunzip -cq ',infile), header=F, verbose = F, showProgress = F)
  } else {
    data <- data.table::fread(infile, header=F, verbose = F, showProgress = F)
  }
  
  SSreadIDs <- as.character(data$V1)
  PBreadIDs <- as.character(data$V6)
  #SSreadIDs <- as.character(data$V6)
  #PBreadIDs <- as.character(data$V1)
  
  SSreadNames.fields <- data.table::tstrsplit(SSreadIDs, "_")
  PBreadNames.fields <- data.table::tstrsplit(PBreadIDs, "_")
  
  SSreadNames <-  SSreadNames.fields[[1]] 
  SSlibNames <- paste(SSreadNames.fields[[2]], SSreadNames.fields[[3]], SSreadNames.fields[[4]], sep="_")
  #SSlibNames <- paste(SSreadNames.fields[[2]], SSreadNames.fields[[3]], sep="_")
  
  SSflag <- SSreadNames.fields[[5]] 
  SSchrom <- SSreadNames.fields[[6]] 
  SSpos <- SSreadNames.fields[[7]] 
  #SSflag <- SSreadNames.fields[[4]] 
  #SSchrom <- SSreadNames.fields[[5]] 
  #SSpos <- SSreadNames.fields[[6]] 
  
  PBreadNames <-  paste(PBreadNames.fields[[1]], PBreadNames.fields[[2]], PBreadNames.fields[[3]], PBreadNames.fields[[4]], PBreadNames.fields[[5]], PBreadNames.fields[[6]], PBreadNames.fields[[7]], sep="_")       
  PBflag <- PBreadNames.fields[[8]] 
  PBchrom <- PBreadNames.fields[[9]] 
  PBpos <- PBreadNames.fields[[10]] 
  
  tab <- data.frame(SSreadNames=SSreadNames, 
                    SSlibNames=SSlibNames, 
                    SSflag=SSflag, 
                    SSchrom=SSchrom, 
                    SSpos=SSpos, 
                    SSreadLen=data$V2,
                    strand=factor(data$V5),
                    PBreadNames=PBreadNames,
                    PBflag=PBflag,
                    PBchrom=PBchrom,
                    PBpos=PBpos,
                    PBreadLen=data$V7,
                    MatchedBases=data$V10,
                    MatchedBasesWithGaps=data$V11,
                    stringsAsFactors = F
  )
  
  if (removeDuplicates) {
    bit.flag <- bitwAnd(1024, as.numeric(tab$SSflag))
    mask <- bit.flag == 0 	
    tab <- tab[mask,]
  }  	
  
  stopTimedMessage(ptm)
  return(tab)
}


#' Get the feature vector based on the WmiunsC ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani

WminusCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}

#' Get the feature vector based on the max{W,C} ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani

maxWandCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- counts[,2]/(counts[,2]+counts[,1])  #calculate ratio of WW reads
    ratios <- sapply(ratios, function(x) max(x, 1-x))
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}

#' Get the feature vector based on both W and C ratios
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author Maryam Ghareghani

WandCratiosFeatures <- function(counts.l=NULL) {
  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]
    
    ratios <- counts/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.nan(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }
  
  ratios.m <- do.call(cbind, ratios.l)
  #ratios.m[ratios.m<0] <- -1
  #ratios.m[ratios.m>0] <- 1
  
  ratios.m
}


#' Rescale theta values for WC cell type
#'
#' Set WC probs based on WW and CC probs distribution (if ratio of WW and CC == 1 then WC region with high prob)
#'
#' @param theta.l A \code{list} of estmated theta values for every single cell.
#' @author David Porubsky

rescaleTheta <- function(theta.l) {
  
  new.theta <- list()
  for (i in 1:length(theta.l)) {
    theta.cell <- theta.l[[i]]
    ratios <- theta.cell[,1] / theta.cell[,2]
    mask <- ratios>0.8 & ratios<1.2
    theta.cell[mask,1] <- 0.05
    theta.cell[mask,2] <- 0.05
    theta.cell[mask,3] <- 0.9
    new.theta[[i]] <- theta.cell
  }
  return(new.theta)
}


#' Merge splitted clusters
#'
#' This function merges initialy divided clusters that belongs to the same chromosome.
#'
#' @param cluster2merge ...
#' @param soft.pVal ...
#' @author David Porubsky

mergeSplitedClusters <- function(cluster2merge=NULL, soft.pVal=NULL) {
  
  merged.pVals <- list()
  for (i in 1:length(cluster2merge)) {
    merge.idx <- cluster2merge[[i]]
    merge.pVal <- apply(soft.pVal[,merge.idx], 1, max)
    merged.pVals[[i]] <- merge.pVal
  }  
  
  return(do.call(cbind,merged.pVals))
} 


#' Get pairs of clusters coming from the same chromosome with different directionalities
#'
#' This function solves the Maximum matching problem for all possible pairs of clusters and reports pairs with highest similarity.
#'
#' @param theta.param A \code{list} containing estimated cell types per cluster for each cell.
#' @return A \code{vector} of pairs of clusters IDs that belong to the same chromosome.
# @importFrom maxmatching maxmatching
# @importFrom igraph graph E
#' @author Maryam Ghareghani, David Porubsky

findClusterPartners2 <- function(theta.param=NULL) {
  
  # get only wc thetas
  theta.param.wc <- lapply(theta.param, function(x) x[,3])
  # cbid wc thetas for all single cells
  all.theta.param.wc <- do.call(cbind, theta.param.wc)
  # compute the pairwise distance of all clusters wc thetas
  d <- as.matrix(dist(all.theta.param.wc))
  # convert distance to a similarity measure
  d <- max(d) - d
  for (i in 1:nrow(d)){d[i,i] <- 0}
  G <- igraph::graph.adjacency(d, mode="undirected", weighted=TRUE)
  #TODO: replce the following code by another package or code (maxmatching doesn't work)
  max.match <- maxmatching::maxmatching(G, weighted = TRUE)
  
  return(max.match$matching)
}


#' Get pairs of clusters coming from the same chromosome with different directionalities
#'
#' This function finds the pairs by computing the pairwise distances between clusters theta_wc parameters, and pais each cluster with the cluster with min dist
#' Note that the first 23 minimum values should form a matching in the complete graph of clusters (this condition is observed to be held in practice)
#' TODO: We should replace this simple function by a proper max matching algorithm or some other solution that properly solves this problem and finds non-overlappinh pairs of clusters
#' 
#' @param theta.param A \code{list} containing estimated cell types per cluster for each cell.
#' @return A \code{vector} of pairs of clusters IDs that belong to the same chromosome.
#' @author Maryam Ghareghani
#' @export

findClusterPartners_simple <- function(theta.param=NULL, clust.to.chrom=NULL,
                                       sex = 'male', verbose=FALSE,
                                       cheating=FALSE) {
  cat('pairing clusters\n')
  # get only wc thetas
  valid.chroms <- paste0('chr', c(1:22, 'X'))
  theta.param.wc <- lapply(theta.param, function(x) x[,3])
  # cbid wc thetas for all single cells
  all.theta.param.wc <- do.call(cbind, theta.param.wc)
  # compute the pairwise distance of all wc thetas (all pairs of clusters)
  d <- as.matrix(dist(all.theta.param.wc))
  # compute the min index in each row (excluding the diogonal element)
  min.cost.clust.pair <- sapply(1:nrow(d), function(i) return(which.min(d[i,-i])))
  min.cost.clust.pair <- data.table(first_clust=1:nrow(d), second_clust=as.numeric(names(min.cost.clust.pair)))
  # put the min rank cluster in the first column
  min.cost.clust.pair[, min_id_clust:=min(first_clust, second_clust), by=1:nrow(min.cost.clust.pair)]
  min.cost.clust.pair[, max_id_clust:=max(first_clust, second_clust), by=1:nrow(min.cost.clust.pair)]
  min.cost.clust.pair[, `:=`(first_clust=min_id_clust, second_clust=max_id_clust)]
  
  min.cost.clust.pair <- unique(min.cost.clust.pair[, .(first_clust, second_clust)])
  
  if(cheating) {
    # Get assignments from alignment
    chroms <- sort(unique(clust.to.chrom$chrom))
    names(chroms) <- chroms
    num.chroms <- length(chroms)
    
    chrom_pairs <- 
      lapply(chroms, function(.x) with(clust.to.chrom, clust[chrom == .x]))
    
    chrom_pairs <-
      lapply(chrom_pairs, function(.x)
        data.frame(first_clust = min(.x), second_clust = max(.x)))
    
    min.cost.clust.pair <- as.data.table(Reduce(rbind, chrom_pairs))

  } else { # From maryam's code
    num.chroms <- 23
    if (sex=='male') {num.chroms <- 24}

    assert_that(
      nrow(min.cost.clust.pair) == num.chroms,
      msg = paste('there are not', num.chroms, 'cluster pairs')
    )
    assert_that(all(!duplicated(min.cost.clust.pair[, first_clust])) &
                  all(!duplicated(min.cost.clust.pair[, second_clust])),
                msg = 'there are duplicate clusters in clust pairing')
  }

 

  if (!is.null(clust.to.chrom)){
    min.cost.clust.pair[, clust:=first_clust]
    min.cost.clust.pair <- merge(min.cost.clust.pair, clust.to.chrom, by='clust')
    min.cost.clust.pair[, clust:=second_clust]
    min.cost.clust.pair <- merge(min.cost.clust.pair, clust.to.chrom, by='clust')
    min.cost.clust.pair[, clust:=NULL]
    # ordering rows by chromosome
    min.cost.clust.pair[, chrom:=factor(chrom.x, levels = valid.chroms)]
    setkey(min.cost.clust.pair, chrom)
    # asserting the correcness of cluster pairing
    min.cost.clust.pair[, is.pairing.correct:=all(chrom.x==chrom.y & c(0,16) %in% c(flag.x, flag.y))
                        , by=1:nrow(min.cost.clust.pair)]
    # assert_that(min.cost.clust.pair[, all(is.pairing.correct)]) %>% invisible()
    if(verbose) print(min.cost.clust.pair)
  }
  # compute the frequencies of each cluster in this pairing
#  fr <- table(c(min.cost.clust.pair$min_id_clust, min.cost.clust.pair$max_id_clust))
  
  # check whether the pairing is a matching
#  assert_that(all(sort(unique(fr))==1:3) & length(which(duplicated(min.cost.clust.pair)))==floor(nrow(min.cost.clust.pair)/2)) %>% invisible()
  
  # keep only the duplicated rows (keep the matching and kick out the garbage cluster)
#  min.cost.clust.pair <- min.cost.clust.pair[duplicated(min.cost.clust.pair)]
  
  min.cost.clust.pair[, chrom_clust:=1:nrow(min.cost.clust.pair)]
  
  return(min.cost.clust.pair[, .(first_clust, second_clust, chrom_clust)])
}
