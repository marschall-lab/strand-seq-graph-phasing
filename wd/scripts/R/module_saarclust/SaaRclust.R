#' Wrapper function to run saarclust pipeline for a given number of long reads reads.
#'
#' @param minimap.file A path to the minimap file to load.
#' @param outputfolder A folder name to export the results.
#' @param num.clusters Expected number of clusters. (for 22 autosomes == 44 clusters)
#' @param minLib Minimal number of StrandS libraries being represent per long PB read
#' @param upperQ Filter out given percentage of PacBio reads with the highest number of alignments.
#' @param EM.iter Number of iteration to run EM for.
#' @param theta.constrain Recalibrate theta values to meet expected distribution of W and C strand across 
#' Strand-seq libraries.
#' @param store.counts Logical if to store raw read counts per PB read
#' @param HC.input File name where hard clustering results are stored
#' @param cellNum Specifies the number of single cells to be used in clustering
#' @inheritParams countProb
#' @inheritParams EMclust
#' @export
#' @author David Porubsky


SaaRclust <- function(minimap.file=NULL, counts.l=NULL, fileID='soft', outputfolder='SaaRclust_results', 
                      num.clusters=47, EM.iter=100, alpha=0.1, minLib=10, upperQ=0.95, theta.param=NULL, 
                      pi.param=NULL, logL.th=1, theta.constrain=FALSE, store.counts=FALSE, HC.input=NULL, 
                      cellNum=NULL, log.scale=FALSE, filter.soft.clust.input=TRUE, filter.ss.file=NULL, 
                      read.names=NULL, clust.pairs=NULL, numCPU=4) {

  print('soft clustering ...')
  #Get ID of a file to be processed
  if (is.null(fileID) & !is.null(minimap.file)){
    fileID <- basename(minimap.file)
    fileID <- strsplit(fileID, "\\.")[[1]][1]
  }

  #Prepare locations for export
  #Create a master output directory if it wasn't created before
  outputfolder <- file.path(outputfolder)
  if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
  }
  
  rawdata.store <- file.path(outputfolder, 'RawData')
  if (!file.exists(rawdata.store)) {
    dir.create(rawdata.store)
  }
  
  Clusters.store <- file.path(outputfolder, 'Clusters')
  if (!file.exists(Clusters.store)) {
    dir.create(Clusters.store)
  }
  
  plots.store <- file.path(outputfolder, 'Plots')
  if (!file.exists(plots.store)) {
    dir.create(plots.store)
  }
  
  trashbin.store <- file.path(outputfolder, 'TrashBin')
  if (!file.exists(trashbin.store)) {
    dir.create(trashbin.store)
  }
  
  ### Write README file ###
  savename <- file.path(outputfolder, 'README.txt')
  
  cat("", file=savename)
  cat("Current folder contains the following folders.\n", file=savename, append=TRUE)
  cat("==============================================\n", file=savename, append=TRUE)
  cat("- Clusters: RData files with the results of the clustering. Contains soft clustering 
      probabilities as well as estimates of theta and pi parameter.\n", file=savename, append=TRUE)
  cat("- RawData: RData files thart contains some measures of data quality. Raw alignment counts 
      are also exported here depending on option 'store.counts=FALSE'.\n", file=savename, append=TRUE)
  cat("- Plots: Some plots produced to evaluate clustering efficacy and accuracy (Left empty for now).\n", 
      file=savename, append=TRUE)
  cat("- TrashBin: This folder contains long reads that belong to the highest 5% based on the Strand-seq 
      read counts. Depends on option 'upperQ=0.95'.\n", file=savename, append=TRUE)
  
  #Load Hard clustering results and initialize parameters of EM algorithm [temporary solution for snakemake]
  #destination <- file.path(Clusters.store, HC.input)
  if (is.null(theta.param) | is.null(pi.param)) {
    if (!file.exists(HC.input)) {
      stop("Hard clustering results not available!!!")
    }
    if (class(HC.input)=="character"){ # (HC.input=filename) input as file
      hard.clust.results <- get(load(HC.input))
    } else { # HC.input=hard clust object
      hard.clust.results <- HC.input
    }
    
    #Initialize theta parameter
    theta.param <- hard.clust.results$theta.param
    #Initialize pi parameter
    pi.param <- hard.clust.results$pi.param
    
  }  
  
  if (is.null(counts.l)) {
    ### Read in minimap alignment file ###
    suppressWarnings( tab.in <- importData(infile = minimap.file) )
    #suppressWarnings( tab.in <- importTestData(infile = minimap.file, removeDuplicates = TRUE) ) 
    #SLOW because test data have to be processed differently
    
    ### get some quality measures on imported data ### [OPTIONAL]
    data.qual.measures <- getQualMeasure(tab.in)
    destination <- file.path(rawdata.store, paste0(fileID, "_dataQuals.RData"))
    save(file = destination, data.qual.measures)
    
    ### Filter imported data ###
    tab.filt <- tab.in
    if (filter.soft.clust.input)
    {
      tab.filt.l <- filterInput(inputData=tab.in, quantileSSreads = c(0, upperQ), minSSlibs = c(minLib,Inf))
      tab.filt <- tab.filt.l$tab.filt
    
      ### Store upperQ reads in a trashBin ###
      upperQ.tab <- tab.in[tab.in$PBreadNames %in% tab.filt.l$upperQ.reads,]
      #upperQ.tab$PBreadNames <- factor(upperQ.tab$PBreadNames, levels=unique(upperQ.tab$PBreadNames))
      #tab.upperQ.l <- split(upperQ.tab, upperQ.tab$SSlibNames)
      #counts.upperQ.l <- countDirectionalReads(tab.upperQ.l)
  
      ptm <- startTimedMessage("Writing upperQ reads into a file")
      destination <- file.path(trashbin.store, paste0(fileID, "_upperQreads.gz"))
      gzf = gzfile(destination, 'w')
      utils::write.table(x = upperQ.tab, file = gzf, quote = F, row.names = F)
      close(gzf)
      stopTimedMessage(ptm)
      #data.table::fwrite(upperQ.tab, destination)
      #gzip(destination)
    }
    
    if (!is.null(filter.ss.file))
    {
      filter.ss.names <- fread(filter.ss.file, header=F)[, V1]
      tab.filt <- tab.filt[!SSreadNames %in% filter.ss.names]
    }
    
    ### Sorting filtered data table by direction and chromosome ###
    ptm <- startTimedMessage("Sorting data")

    #split data by SS library
    tab.l <- split(tab.filt, tab.filt$SSlibNames)
    stopTimedMessage(ptm)
    
    ### Count directional reads ###
    counts.l <- countDirectionalReads(tab.l)
  
    # subsetting single cell libraries
    if (!is.null(cellNum)) {
      counts.l = counts.l[1:cellNum]
    }
    
    if (store.counts) {
      destination <- file.path(rawdata.store, paste0(fileID, "_counts.RData"))
      save(file = destination, counts.l)
    }
  }
  
  ### EM algorithm ###
  soft.clust.obj <- EMclust(counts.l, theta.param=theta.param, pi.param=pi.param, 
                            num.iter=EM.iter, alpha=alpha, logL.th=logL.th, 
                            log.scale=log.scale, numCPU=numCPU)
  
  #rescale theta parameter and run one more iteration to redo soft clustering [EXPERIMENTAL]
  if (theta.constrain) {
    theta.expected <- num.clusters * c(0.25,0.25,0.5)
    theta.rescaled <- thetaRescale(theta.param=soft.clust.obj$theta.param, theta.expected=theta.expected)
    soft.clust.obj <- EMclust(counts.l, theta.param=theta.rescaled, pi.param=soft.clust.obj$pi.param, 
                              num.iter=1, alpha=alpha, logL.th=logL.th, log.scale=log.scale)
  }
  
  #Get pairs of clusters coming from the same chromosome but differs in directionality of PB reads  [EXPERIMENTAL]
  #theta.sums <- Reduce("+", soft.clust.obj$theta.param)
  #remove.clust <- which.max(theta.sums[,3])
  #theta.param.filt <- lapply(soft.clust.obj$theta.param, function(x) x[-remove.clust,])
  #clust.order <- findClusterPartners(theta.param=theta.param.filt)

  ### Save final results ###
  #add known chromosome and directionality of PB reads to a final data object
#  soft.clust.obj$PBchrom <- as.character(chr.rows)
#  soft.clust.obj$PBflag <- as.character(pb.flag)
#  soft.clust.obj$pb.readLen <- tab.filt$PBreadLen[match(rownames(soft.clust.obj$soft.pVal), tab.filt$PBreadNames)]  #report PB read length
  #export data in RData object
#  destination <- file.path(Clusters.store, paste0(fileID, "_clusters.RData"))
  
  rownames(soft.clust.obj$soft.pVal) <- read.names
  ML.clust <- apply(soft.clust.obj$soft.pVal, 1, which.max)
  ML.clust <- data.table(rname=names(ML.clust), first_clust=ML.clust)
  
  if (!is.null(clust.pairs)){
    extend.clust.pairs <- rbind(clust.pairs, data.table(first_clust =clust.pairs$second_clust, 
                                                        second_clust=clust.pairs$first_clust,
                                                        chrom_clust =clust.pairs$chrom_clust))
    ML.clust <- merge(ML.clust, extend.clust.pairs, by='first_clust')
    ML.clust <- ML.clust[, .(rname, first_clust, chrom_clust)]
  }
  
  soft.clust.obj$ML.clust <- ML.clust
  
#  save(file = destination, soft.clust.obj)

  return(soft.clust.obj)
}
