## Load required libraries
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)

## FUNCTIONS ##
###############

strandphaser <- function(cluster1.table, cluster2.table, clusters, select.libs, output.phased.strand.states.file){
	phased.data = phase_strand_states_and_bubbles(as.data.frame(cluster1.table), as.data.frame(cluster2.table), clusters[1], clusters[2], select.libs)
	phased.strand.states = phased.data[[1]]
	phased.bubbles = phased.data[[2]]

	fwrite(phased.strand.states, file=output.phased.strand.states.file, row.names=F, sep='\t')
	#fwrite(phased.bubbles, file=output.phased.bubbles.file, row.names=F, sep='\t')
}


phase_strand_states_and_bubbles <- function(cluster1, cluster2, cluster1.name, cluster2.name, select.libs) {
	## Load bubbles into matrices
	matrices <- loadClusterBubbles(cluster1 = cluster1,
		                 cluster2 = cluster2)

	## Sort matrices
	srt.matrices <- sortClusterBubbles(matrices = matrices, select.libs = select.libs)

	## adding the haplo-phased strand states to the

	haplo.strand.states <- data.table()

	clust1.haplo.strand.states <- data.table(lib=srt.matrices$cluster1.libs, cluster=cluster1.name)
	clust1.haplo.strand.states[,`:=`(lib=lapply(lib, function(x) strsplit(x, "__")[[1]][1]), haplotype=sapply(lib, function(x) as.integer(strsplit(x, "__C")[[1]][2])-1))]
	haplo.strand.states <- rbind(haplo.strand.states, clust1.haplo.strand.states)

	clust2.haplo.strand.states = clust1.haplo.strand.states
	clust2.haplo.strand.states[, `:=`(cluster=cluster2.name, haplotype=1-haplotype)]
	haplo.strand.states <- rbind(haplo.strand.states, clust2.haplo.strand.states)

	## Get consensus
	h1.cons <- exportConsensus(data.matrix = srt.matrices$cluster1.m)
	bubble.phase <- data.table(bubbleName=h1.cons$pos, haplotype0Allele=h1.cons$hap-1)

	h2.cons <- exportConsensus(data.matrix = srt.matrices$cluster2.m)
	bubble.phase <- rbind(bubble.phase, data.table(bubbleName=h2.cons$pos, haplotype0Allele=h2.cons$hap-1))

  return(list(haplo.strand.states, bubble.phase))
}

loadClusterBubbles <- function(cluster1 = NULL, cluster2 = NULL) {

  cluster1.libs <- colnames(cluster1)[-1]
  cluster1.bubble.id <- cluster1[,1]
  cluster1 <- cluster1[,-1]
  # Switch missing values to 0 and reference alleles to 1 and alternative alleles to 2
  cluster1.m <- t(mutate_all(cluster1, recode, '-' = 0, '0' = 1, '1' = 2, '2' = 3, .default = 0))#t(apply(cluster1, 2, function(x) recode(x, '-' = 0, '0' = 1, '1' = 2, '2' = 3, .default = 0)))
  attr(cluster1.m, 'dimnames') <- NULL

  cluster2.libs <- colnames(cluster2)[-1]
  cluster2.bubble.id <- cluster2[,1]
  cluster2 <- cluster2[,-1]
  # Switch missing values to 0 and reference alleles to 1 and alternative alleles to 2
  cluster2.m <-t(mutate_all(cluster2, recode, '-' = 0, '0' = 1, '1' = 2, '2' = 3, .default = 0)) # t(apply(cluster2, 2, function(x) recode(x, '-' = 0, '0' = 1, '1' = 2, '2' = 3, .default = 0)))
  attr(cluster2.m, 'dimnames') <- NULL

  final.object <- list(cluster1.m = cluster1.m, cluster1.libs = cluster1.libs, cluster1.bubble.id = cluster1.bubble.id,
                       cluster2.m = cluster2.m, cluster2.libs = cluster2.libs, cluster2.bubble.id = cluster2.bubble.id)
  return(final.object)
}

sortClusterBubbles <- function(matrices = NULL, select.libs = NULL) {
  ## Load data
  cluster1.m <- matrices[['cluster1.m']]
  cluster1.libs <- matrices[['cluster1.libs']]
  cluster1.bubble.id <- matrices[['cluster1.bubble.id']]
  cluster2.m <- matrices[['cluster2.m']]
  cluster2.libs <- matrices[['cluster2.libs']]
  cluster2.bubble.id <- matrices[['cluster2.bubble.id']]

  ## Filter and sort shared bubbles IDs
  shared.bubble.id <- intersect(cluster1.bubble.id, cluster2.bubble.id)
  shared.bubble.id <- sort(shared.bubble.id)
  cluster1.m <- cluster1.m[,match(shared.bubble.id, cluster1.bubble.id), drop=FALSE]
  cluster1.bubble.id <- shared.bubble.id
  cluster2.m <- cluster2.m[,match(shared.bubble.id, cluster2.bubble.id), drop=FALSE]
  cluster2.bubble.id <- shared.bubble.id

  ## Filter shared cells
  shared.libs <- intersect(cluster1.libs, cluster2.libs)
  ## NOTE library names are not consistent: sometimes they include dot and sometimes don't
  if (!is.null(select.libs)) {
    shared.libs <- intersect(shared.libs, select.libs)
  }
  ## Select rows present in shared cells
  mask.rows.cluster1 <- match(shared.libs, cluster1.libs)
  mask.rows.cluster2 <- match(shared.libs, cluster2.libs)
  cluster1.m <- cluster1.m[mask.rows.cluster1, ,drop=FALSE]
  cluster2.m <- cluster2.m[mask.rows.cluster2, ,drop=FALSE]
  cluster1.libs <- cluster1.libs[mask.rows.cluster1]
  cluster2.libs <- cluster2.libs[mask.rows.cluster2]
  cluster1.libs <- paste(cluster1.libs, "C1", sep="__")
  cluster2.libs <- paste(cluster2.libs, "C2", sep="__")

  ## Filter homozygous SNV positions
  # mask.cl1 <- apply(cluster1.m, 2, function(x) any(x == 3))
  # mask.cl2 <- apply(cluster2.m, 2, function(x) any(x == 3))
  # mask <- c(which(mask.cl1 == TRUE), which(mask.cl2 == TRUE))
  # mask <- unique(mask)
  # if (length(mask) > 0){
  #   cluster1.m <- cluster1.m[,-mask, drop=FALSE]
  #   cluster2.m <- cluster2.m[,-mask, drop=FALSE]
  # }
  #


	consensus.score = data.table()
	swap.consensus.score = data.table()

  ## sort parallel matrices
  for ( i in 1:length(shared.libs)) {
    filename <- shared.libs[i]
    message("Processing ", filename, " ...")

    cluster1.pos <- which(cluster1.m[i,] > 0)
    cluster2.pos <- which(cluster2.m[i,] > 0)
    cov.pos <- union(cluster1.pos, cluster2.pos)

    ## initialize score
    cluster1.m.score <- calcMatrixScore(matrix = cluster1.m, covered.pos = cov.pos)
    cluster2.m.score <- calcMatrixScore(matrix = cluster2.m, covered.pos = cov.pos)

    ## swap rows in matrices
    cluster1.filename <- cluster1.libs[i]
    cluster2.filename <- cluster2.libs[i]
    cluster1.libs[i] <- cluster2.filename
    cluster2.libs[i] <- cluster1.filename

    cluster1.row <- cluster1.m[i,]
    cluster2.row <- cluster2.m[i,]
    cluster1.m[i,] <- cluster2.row
    cluster2.m[i,] <- cluster1.row

    ## calculate new score
    curr.cluster1.m.score <- calcMatrixScore(matrix = cluster1.m, covered.pos = cov.pos)
    curr.cluster2.m.score <-  calcMatrixScore(matrix = cluster2.m, covered.pos = cov.pos)

		## store consensus scores in data tables
		consensus.score = rbind(consensus.score, data.table(lib=filename, cons.score=cluster1.m.score+cluster2.m.score))
		swap.consensus.score = rbind(swap.consensus.score, data.table(lib=filename, swap.cons.score=curr.cluster1.m.score+curr.cluster2.m.score))

    ## compare previous matrix score with score after swapping rows
    if ( (cluster1.m.score + cluster2.m.score) < (curr.cluster1.m.score + curr.cluster2.m.score) ) {
      cluster1.libs[i] <- cluster1.filename
      cluster2.libs[i] <- cluster2.filename

      cluster1.m[i,] <- cluster1.row
      cluster2.m[i,] <- cluster2.row
    }
  }

  ## Export sorted data
  sorted.matrices <- list()
  sorted.matrices[['cluster1.m']] <- cluster1.m
  sorted.matrices[['cluster1.libs']] <- cluster1.libs
  sorted.matrices[['cluster2.m']] <- cluster2.m
  sorted.matrices[['cluster2.libs']] <- cluster2.libs
  sorted.matrices[['cluster.bubble.id']] <- shared.bubble.id
	sorted.matrices[['consensus.score']] <- consensus.score
	sorted.matrices[['swap.consensus.score']] <- swap.consensus.score
  return(sorted.matrices)
}

calcMatrixScore <- function(matrix, covered.pos) {
  sub.m <- matrix[,covered.pos,drop = FALSE]
  sub.m <- sub.m[,colSums(sub.m) > 0, drop = FALSE]
  sub.m.long <- melt(sub.m)
  base.freq <- sub.m.long %>% group_by(Var2, value) %>% summarise(count=n()) %>% ungroup()
  base.freq <- base.freq[base.freq$value > 0,]

  # print(base.freq)
  scores <- base.freq %>% group_by(Var2) %>% summarise(count = sum(count) - max(count))

  return(sum(scores$count))
}

exportConsensus <- function(data.matrix) {
  #base.freq <- apply(data.bases, 2, function(x) table(x[x>0])) ## TODO avoid warnings here
  indices <- which(data.matrix != 0,arr.ind = TRUE) # get indices of non-zero values
  values <- data.matrix[indices] # get non-zero values based on indices
  col.vals <- split(values, (indices[,2])) # split values based on column indices

  base.freq <- lapply(col.vals, table) # get base frequencies for each column

  positions <- as.numeric(names(base.freq)) # get position of each snv/column

  ## Filter SNV positions (columns) which do not pass set criteria (coverge, ambiguity)
  max.cov <- sapply(base.freq, function(x) max(x)) # get coverage of highest covered base
  score <- sapply(base.freq, function(x) sum(x)-max(x))

  if (length(positions) != 0) {

    ## Calculate entropy values for each column
    entropy <- c(rep(0, length(positions)))
    entropy <- sapply(col.vals, calcEnt)

    ## Calculate Phred score values for each column
    scores <- c(rep(0, length(positions)))

    cov = sapply(col.vals, length)

    hap <- sapply(base.freq, function(x) as.numeric(names(which.max(x))))

    assem.haps <- data.frame(pos=positions, hap=hap, cov=cov, score=scores, ent=entropy)
    rownames(assem.haps) <- NULL
    return(assem.haps)
  } else {
    return(0)
  }
}

calcEnt <- function(bases) {

  # log2 function (modif: return zero for zero values)
  getlog2 <- function(x) {
    if (x == 0 ) {
      return(0)
    } else {
      return(log2(x))
    }
  }

  counts <- c(0,0,0,0)

  # count frequency of each base in a column
  for (i in 1:4) {
    counts[i] <- length(bases[bases == i])
  }
  probs <- counts/sum(counts) # get probabilities

  ent <- -( probs[1]*getlog2(probs[1]) + probs[2]*getlog2(probs[2]) + probs[3]*getlog2(probs[3]) + probs[4]*getlog2(probs[4]) )
  return(ent)
}
