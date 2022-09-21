ClustersAccuracyPerChrPerDir <- function(soft.clust.files=NULL, bam.file=NULL, thresholds=c(0.5, 0.6, 0.7, 0.8, 0.9, 0.99), minLib=NULL) {
#  Clusters2process <- list.files(file.path(inputfolder, 'Clusters'), pattern = "_clusters.RData", full.names = TRUE)
#  Quals2process <- list.files(file.path(inputfolder, 'RawData'), pattern = "dataQuals.RData", full.names = TRUE)
  #PBlen2process <- list.files("/media/daewoooo/WORK/Clustering_project/WholeGenomeAnalysis/PBreadLen/", pattern = "gz", full.names = TRUE)

  allClusters <- list()
  for (i in 1:length(soft.clust.files)) {
    data.file <- get(load(soft.clust.files[i]))
    data.file <- addChromFlag(data.file, bam.file, clust.type='soft')
#    data.qual <- get(load(Quals2process[i]))
    #PB.read.len <- data.table::fread(paste0('zcat ', PBlen2process[i]), header=T, verbose = F, showProgress = F)
    fileID <- basename(soft.clust.files[i])
    message("Processing file: ",fileID)
    
    #Sort data quals according to PB order in clusters
#    SSlib.perPB <- data.qual$SSlib.perPB
#    SSlib.perPB <- SSlib.perPB[match(rownames(data.file$soft.pVal), SSlib.perPB$PBreadNames),]
#    pb.minLib <- SSlib.perPB$counts
    
    #Select required PB read lengths
    #PB.read.len <- PB.read.len[match(rownames(data.file$soft.pVal), PB.read.len$PBreadNames),]
#    pb.readLen <- data.file$pb.readLen
    
    ##check accuracy only for autosomes and sex chrmosomes
    mask <- which(grepl('^chr[0-9X][0-9]?$', data.file$chrom))
    
    #get clusters IDs corresponding to a given chromosome
    chr.rows <- data.file$chrom[mask]
    chr.flag <- data.file$flag[mask]
    prob.tab <- data.file$soft.pVal[mask,]
#    pb.minLib <- pb.minLib[mask]
#    pb.readLen <- pb.readLen[mask]
    
    #filter out duplicates
    mask <- which(chr.flag == 16 | chr.flag == 0) 
    chr.rows <- chr.rows[mask]
    chr.flag <- chr.flag[mask]
    prob.tab <- prob.tab[mask,]
#    pb.minLib <- pb.minLib[mask]
#    pb.readLen <- pb.readLen[mask]
    
    #Find WC cluster in all cells
    theta.sums <- Reduce("+", data.file$theta.param)
    remove.clust <- which.max(theta.sums[,3])
    #Remove probabilities for always WC cluster
    prob.tab <- prob.tab[,-remove.clust]
    
#    #Remove PB reads represneted by SSlib less than minLib
#    filt <- pb.minLib >= minLib
#    prob.tab <- prob.tab[filt,]
#    chr.rows <- chr.rows[filt]
#    chr.flag <- chr.flag[filt]
#    pb.readLen <- pb.readLen[filt]
    
    #get clusters IDs corresponding to a given chromosome
    Clust.IDs <- getClusterIdentityPerChrPerDir(soft.clust=prob.tab, chr.rows=chr.rows, chr.flag=chr.flag)
    
    clust.acc.l <- list()
    for (prob.th in thresholds) {
      message("    Set threshold: ", prob.th)
      max.prob <- apply(prob.tab, 1, max)
      mask <- max.prob >= prob.th
#      pb.readLen.sub <- pb.readLen[mask]
      
      Clust.locations <- apply(prob.tab[mask,], 1, which.max) 
      
      #calculate clustering accuracy in comparison to expected values
      clust.acc <- Clust.locations == Clust.IDs[mask]
      acc.th <- table(clust.acc)
      
      clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, acc.th.match=unname(acc.th[2]), acc.th.sum=sum(acc.th), allReads=length(chr.rows)) #, seq.bases=sum(as.numeric(pb.readLen.sub))) 
      #clust.acc.l[[1+length(clust.acc.l)]] <- c(prob.th=prob.th, acc.th.match=unname(acc.th[2]), acc.th.sum=sum(acc.th), allReads=length(chr.rows))
    }
    allClusters[[fileID]] <- as.data.frame( do.call(rbind, clust.acc.l) )
  } 
  #sum all counts over all data frames (per position)
  clust.acc.df <- Reduce("+", allClusters)
  
  #calcualte accuracy percentages
  clust.acc.df$prob.th <- thresholds
  clust.acc.df$th.acc <- clust.acc.df$acc.th.match / clust.acc.df$acc.th.sum
  clust.acc.df$th.clustReads <- clust.acc.df$acc.th.sum / clust.acc.df$allReads
  
  #get genome size
  library("biovizBase")
  hg38Ideogram <- getIdeogram("hg38", cytoband = FALSE)
  hg38Ideogram <- keepSeqlevels(hg38Ideogram, paste0('chr', c(1:22,'X')), pruning.mode = 'coarse')
  genome.size <- sum(as.numeric(seqlengths(hg38Ideogram)))
#  clust.acc.df$depth <- ceiling(clust.acc.df$seq.bases/genome.size)
  
  acc.plt <- ggplot(clust.acc.df) + 
    geom_point(aes(x=th.acc, y=th.clustReads), color="deepskyblue4", size=10) + 
    geom_linerange(aes(ymin=-Inf, x=th.acc, ymax=th.clustReads),color="deepskyblue4") + 
    scale_x_continuous(limits = c(0,1)) + 
    scale_y_continuous(limits = c(0,1)) + 
    ylab("(%) evaluated PB reads") + xlab("(%) correctly assigned PB reads") + 
    geom_text(aes(x=th.acc, y=th.clustReads), label=c('all', thresholds[-1]), color="white") + 
    geom_text(aes(x=th.acc, y=th.clustReads+0.05), label=paste0(clust.acc.df$depth, "x"), color="black") +
    theme_bw()
  message("DONE!!!")
  return(list(acc.plot=acc.plt, plot.table=clust.acc.df))
}