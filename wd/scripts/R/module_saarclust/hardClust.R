#' Hard clustering using k-means
#'
#' @param counts.l A \code{list}/\code{data.table} of directional read counts per PB read per library.
#' @inheritParams SaaRclust
#' @return A numeric vector including the clusters, named with object (long read/unitig) names
#' @author David Porubsky
#' @export

hard.kmeans <- function(counts.l=NULL, num.clusters=NULL, nstart=10, iter.max=10, by_chrom=F) {

  message("Hard clustering")
  ptm <- startTimedMessage("    Kmeans clustering for ",num.clusters," clusters")

  if ('data.table' %in% class(counts.l)){
    counts.l <- convert.counts.format(counts.l)
  }

  ratios.l <- list()
  for (j in 1:length(counts.l)) {
    #lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    counts <- counts.l[[j]]

    ratios <- (counts[,2]-counts[,1])/(counts[,2]+counts[,1]) #calculate ratio of WW reads
    ratios[is.na(ratios)] <- 0
    ratios.l[[j]] <- ratios
  }

  ratios.m <- do.call(cbind, ratios.l)

  if (by_chrom){
    ratios.m <- abs(ratios.m)
  }

  #hard clustering using kmeans
  km <- suppressWarnings( kmeans(ratios.m, centers = num.clusters, nstart = nstart, iter.max = iter.max) )
  ord <- km$cluster
  names(ord) <- rownames(counts.l[[1]])
  #ratios.m.ord <- ratios.m[order(ord),]
  stopTimedMessage(ptm)
  return(ord)
}

#' Getting wfrac matrix
#'
#' @param counts.dt A \code{data.table} of w/c read counts per lib/long read with the column names rname, w, c, lib.name.
#' @inheritParams SaaRclust
#' @return A \code{matrix} of w.frac values (w-c)/(w+c) with rownames=long read names and column names=single-cell library names
#' @author Maryam Ghareghani
#' @export
#'

get.wfrac.matrix <- function(counts.dt){
  w.frac <- counts.dt

  w.frac[, W.frac:=(w-c)/(w+c)]

  mat <- dcast(w.frac, rname~lib, value.var='W.frac')

  row.names <- mat[, rname]
  mat[, rname:=NULL]
  mat <- as.matrix(mat)
  rownames(mat) <- row.names
  mat[is.na(mat)] <- 0

  mat
}

#' Hard clustering using hclust
#'
#' @param counts.l A \code{list}/\code{data.table} of directional read counts per PB read per library.
#' @param num.clusters Number of output clusters
#' @param by_chrom a logical value (default FALSE) indicating whether to cluster only by chromosome
#' @param min.cluster.frequency Minimum cluster frequncy for filtering out very small clusters.
#' @param chrom.flag The ground true \code{data.table} containing the chrom/flags of long reads/unitigs
#' @param chrom.bar.plot plots the bar plot of the number of long reads/unitigs per chrom/direction. Works with the given chrom.flag input.
#' @return A numeric vector including the clusters, named with object (long read/unitig) names
#' @author Maryam Ghareghani
#' @export

hard.hclust <- function(counts.l=NULL, num.clusters=NULL,
                        by_chrom=F, min.cluster.frequency=0.002,
                        chrom.flag=NULL, chrom.bar.plot=F) {

  message("Hard clustering")
  ptm <- startTimedMessage("    Hierarchical clustering for ",num.clusters," clusters")

  if (class(counts.l)[1]=='list'){
    counts.l <- convert.counts.format(counts.l)
  }

  wfrac.matrix <- get.wfrac.matrix(counts.dt=counts.l)

  mat <- wfrac.matrix
  if (by_chrom){
    mat <- abs(wfrac.matrix)
  }

  d <- dist(mat)

  hc <- hclust(d)
  hc.clust <- cutree(hc, k=num.clusters)

  hardclust <- data.table(rname=names(hc.clust), clust=hc.clust)
  hardclust[, num.clust:=.N, by=clust]
  # computing cluster counts
  clust.num <- hardclust[, .(rname,clust,num.clust)]

  if (chrom.bar.plot){
    clust.num <- merge(clust.num, clust.to.chrom, by='clust')
    ggplot(clust.num) + geom_bar(aes(x=clust, fill=is.na(chrom)))
  }

  clust.num <- clust.num[, head(.SD,1), by=clust]

  # filtering out extremly small clusters
  min.clust.num <- nrow(hardclust)*min.cluster.frequency
  hardclust <- hardclust[num.clust > min.clust.num]

  # re-number clusters from 1 to #clusters
  hardclust[, new.clust:=frank(hardclust, clust, ties.method="dense")]
  setkey(hardclust, new.clust)

  hc.clust <- hardclust[, new.clust]
  names(hc.clust) <- hardclust[, rname]

  if (!is.null(chrom.flag))
  {
    # evaluation of hard clustering
    valid.chroms <- paste0('chr', c(1:22, "X"))
    expected.clusters <- paste0(rep(valid.chroms, each=2), rep(c('_0','_16'),23))

    if (!by_chrom){
      clust.to.chrom <- numFoundClusters(ord=hc.clust, chrom.flag=chrom.flag)
    } else {
      numFoundChromosomes(ord=hc.clust, chrom.flag=chrom.flag)
    }
  }

  stopTimedMessage(ptm)
  return(hc.clust)
}

#' Hard clustering wrapper function
#'
#' @param method hard clustering method (default=hclust). It can be "hclust" or "kmeans".
#' @inheritParams hard.kmeans
#' @inheritParams hard.hclust
#' @author Maryam Ghareghani
#' @export

hardClust <- function(counts.l=NULL, method='hclust', num.clusters=NULL, nstart=10, iter.max=10,
                      by_chrom=F, min.cluster.frequency=0.002, chrom.flag=NULL, chrom.bar.plot=F) {
  if (method=='kmeans') {
    set.seed(1000) #in order to reproduce hard clustering results
    hardClust.ord <- hard.kmeans(counts.l, num.clusters=num.clusters, nstart=nstart,
                                 iter.max=iter.max, by_chrom=by_chrom)
    return(hardClust.ord)

  } else if (method=="hclust") {
    # Hierarchical clustering
    return(hard.hclust(counts.l=counts.l, num.clusters=num.clusters, by_chrom=by_chrom,
                                 min.cluster.frequency=min.cluster.frequency, chrom.flag=chrom.flag,
                                 chrom.bar.plot=chrom.bar.plot))
  } else {
    warning("method should be kmeans or hclust")
    return(NULL)
  }
}

#' Estimate theta values based on hard clustering
#'
#' This function takes results of hard clustering and estimates majority cell types for each Strand-seq library
#'
#' @param counts.l A \code{list} of directional read counts per PB read per library.
#' @param hard.clust A \code{integer} of cluster assignments for each PacBio read.
#' @param method It can be "median" or "prob". In "median" method, thetas are estimated based on median wfrac in lib/clusters. In "prob" method, theta is estmated by binomial distribution.
#' @inheritParams SaaRclust
#' @return A \code{list} of estimated theta values for every cluster and cell.
#' @author David Porubsky
#' @export

#
estimateTheta <- function(counts=NULL, hard.clust=NULL, alpha=0.1, method='median') {
  ptm <- startTimedMessage("Estimate theta values")

  if (method=='median') {
    if (class(counts)[1]=='list'){
      counts <- convert.counts.format(counts)
    }

    if (!'W.frac' %in% colnames(counts)) {
      counts[, W.frac:=(w-c)/(w+c)]
    }
    counts[is.nan(W.frac), W.frac:=0]

    clust <- data.table(rname=names(hard.clust), clust=hard.clust)
    theta <- merge(counts[, .(rname,lib, W.frac)], clust[, .(rname, clust)], by='rname')
    theta[, median.w.frac:=median(W.frac), by=.(clust, lib)]

    theta[, `:=`(theta.ww=0, theta.cc=0, theta.wc=0)]
    theta[median.w.frac < -0.5, `:=`(theta.ww=0.05, theta.cc=0.9, theta.wc=0.05), by=.(clust,lib)]
    theta[median.w.frac >= -0.5 & median.w.frac <= 0.5, `:=`(theta.ww=0.05, theta.cc=0.05, theta.wc=0.9), by=.(clust,lib)]
    theta[median.w.frac > 0.5, `:=`(theta.ww=0.9, theta.cc=0.05, theta.wc=0.05), by=.(clust,lib)]

    theta <- theta[, head(.SD,1), by=.(clust, lib)]
    setkey(theta, clust, lib)

    theta.estim <- split(theta[, .(lib, clust, theta.ww, theta.cc, theta.wc)], by='lib') %>% invisible()

    for (i in 1:length(theta.estim))
    {
      clust.names=theta.estim[[i]][, clust]
      theta.estim[[i]][, `:=`(lib=NULL, clust=NULL)]
      theta.estim[[i]] <- as.matrix(theta.estim[[i]])
      rownames(theta.estim[[i]])=clust.names
    }
  } else { # probabilistic method

    if ('data.table' %in% class(counts)){
      counts.l <- convert.counts.format(counts)
    }

    theta.estim <- list()
    for (j in 1:length(counts.l)) {
      minus.c <- split(counts.l[[j]][,1], hard.clust)
      plus.c <- split(counts.l[[j]][,2], hard.clust)

      # adjustment for data.table count format
      if (length(minus.c)>0 & "data.table" %in% class(minus.c[[1]])){
        minus.c <- lapply(minus.c, function(x) x[,w])
        plus.c <- lapply(plus.c, function(x) x[,c])
      }

      clust.prob <- mapply(function(X,Y) { countProb(X,Y) }, X=minus.c, Y=plus.c)
      clust.prob.norm <- lapply(clust.prob, function(x) colSums(log(x)))
      estimates <- sapply(clust.prob.norm, which.max)
      # FIXME: theta_wc is always the most likely theta!!!

      #Assign cell type probs based on the majority type in each cluster
      probs <- list()
      for (i in 1:length(clust.prob)) {
        estim <- estimates[i]
        if (estim == 1) {
          theta <- c(1-alpha, alpha/2, alpha/2)
        } else if (estim == 2) {
          theta <- c(alpha/2, 1-alpha, alpha/2)
        } else {
          theta <- c(alpha/2, alpha/2, 1-alpha)
        }
        probs[[i]] <- theta
      }
      probs <- do.call(rbind, probs)
      theta.estim[[j]] <- probs
      names(theta.estim)[j] <- names(counts.l)[j]
    }
  }

  stopTimedMessage(ptm)

  return(theta.estim)
}


#' Hierarchical clustering for merging hard clusters.
#'
#' This function takes as input the hard clustering output and the initialized thetas and merges the hard clusters based on thetas
#'
#' @param theta.l A \code{list} of estimated theta values for each cluster and cell.
#' @param hard.clust The kmeans hard clustering.
#' @param k Desired number of clusters after merging.
#' @param row.clusters the cluster ids of the rows of theta.l
#' @inheritParams SaaRclust
#' @return A new hard clustering with the correct number of clusters
#' @author Maryam Ghareghani
#' @export

mergeClusters <- function(hard.clust, theta.l, k=46)
{
  ptm <- startTimedMessage("Merging clusters")

  theta.all <- do.call(cbind, theta.l)

  d <- dist(theta.all)
  hc <- hclust(d)
  merged <- cutree(hc, k=k)

  hard.clust.merged <- merged[match(hard.clust, names(merged))]
  names(hard.clust.merged) <- names(hard.clust)

  stopTimedMessage(ptm)
  return(hard.clust.merged)
}

#' Clustering SS reads based on SaaRclust clusting of long reads.
#'
#'
#' @param alignments A \code{list} of alingment data tables (cols=rname, qname, strand) per single cell
#' @param clusters \code{data.table} of clusters of long reads (cols=rname, first_clust).
#' @param clust.pairs a data table (cols=...) representing the cluster pairs.
#' @param numCPU number of threads for parallelization.
#' @return Clusters of SS reads named with the read names.
#' @author Maryam Ghareghani
#' @export

cluster.ss.reads <- function(alignments, clusters, clust.pairs, # numCPU=4,
                            ss.bam.dir=NULL, ss.bam.suffix='_haplotagged.bam',
                            clust.to.chrom=NULL){
  ##cl <- makeCluster(numCPU)
  ##doParallel::registerDoParallel(cl)

  ## parallel processing does not work!
  ##parallel.results <- foreach (i=1:length(alignments),
  ##                             .packages=c('Rsamtools', 'data.table'),
  ##                             .export='getChromFlag') %dopar%{

  message('clustering ss reads')
  ss.clust <- data.table()
  acc <- data.table()

  extend.clust.pairs <- rbind(clust.pairs, data.table(first_clust =clust.pairs$second_clust,
                                               second_clust=clust.pairs$first_clust,
                                               chrom_clust =clust.pairs$chrom_clust))

  for (i in 1:length(alignments)) {
    aln <- alignments[[i]]
    lib.name <- names(alignments)[i]
    # aln[, qname:=paste0(qname, '_', lib.name)]
    ptm <- startTimedMessage(paste("lib", lib.name))

    clust <- clusters[, .(rname, first_clust)]
    aln <- merge(aln, clust, by='rname')
    aln <- merge(aln, extend.clust.pairs, by='first_clust')

    aln[, clust:=first_clust]
    suppressWarnings(aln[strand=='-', clust:=second_clust])

    # check the uniqueness of clusters for ss reads
    aln[, ss.cluster.range:=max(clust)-min(clust), by='qname']

    # cluster only the reads that have unique clusters
    aln <- aln[ss.cluster.range==0]

    ##acc <- NULL
    if (!is.null(ss.bam.dir)){
      ss.bam.file <- file.path(ss.bam.dir, paste0(lib.name, ss.bam.suffix))
      chrom.flag <- getChromFlag(ss.bam.file)
      chrom.flag[, qname:=rname]
      clust.eval <- merge(aln, chrom.flag[, .(qname, chrom, flag)], by='qname')
      clust.eval <- merge(clust.eval, clust.to.chrom, by='clust')
      num.true <- nrow(clust.eval[chrom.x==chrom.y & flag.x==flag.y])
      num.false <- nrow(clust.eval)-num.true

      ##acc <- c(num.true=num.true, num.false=num.false)
      acc <- rbind(acc, data.table(num.true=num.true, num.false=num.false))
    }

    ## for foreach
    ##ss.clust <- aln[, .(qname, cluster)]
    ##list(ss.clust=ss.clust, acc=acc)

    ss.clust <- rbind(ss.clust, aln[, .(qname, clust, chrom_clust)])
    stopTimedMessage(ptm)
  }

  if (!is.null(ss.bam.dir)) {
    num.true <- sum(acc[, num.true])
    num.false <- sum(acc[, num.false])
    cat('number of clustered SS reads =', num.true+num.false,
        ', SS clustering accuracy =', num.true/(num.true+num.false), '\n')
  }


  ##stopCluster(cl)

  return(ss.clust)
}
