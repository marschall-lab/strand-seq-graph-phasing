#' Import pre-processed output from the minimap
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param infile A path to the minimap file to load.
#' @return A \code{data.frame}.
#' @importFrom data.table fread
#' @author David Porubsky
#' @export

importData <- function(infile=NULL) {  #TODO modify this function for input where genomic location of PB reads is unknown

  ptm <- startTimedMessage("Reading the data")
  #data <- read.table(infile, header=F) #TODO test data.table package for faster data import

  #filetype = summary( file(infile) )$class #If it's gzipped, filetype will be gzfile
  if (summary( file(infile) )$class == 'gzfile') {
    data <- data.table::fread(paste0('zcat ',infile), header=T, verbose = F, showProgress = F)
    # select columns:
    data <- data[, .(SSreadNames, SSlibNames, SSflag, SSchrom, SSpos, strand, PBreadNames, PBflag, PBchrom, PBpos, PBreadLen, TargetCoordStart, TargetCoordend, MatchedBasesWithGaps)]
  } else {
    data <- data.table::fread(infile, header=T, verbose = F, showProgress = F)
  }

  #make sure strand info is represented as factor variable
  data$strand <- factor(data$strand)

  stopTimedMessage(ptm)
  return(data)
}


#' Import output from the minimap
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param infile A path to the minimap file to load.
#' @param removeDuplicates Set to \code{TRUE} if you want to remove duplicate reads based on the flag.
#' @return A \code{data.frame}.
#' @importFrom data.table fread tstrsplit
#' @author David Porubsky

importTestData <- function(infile=NULL, removeDuplicates = TRUE) {  #TODO modify this function for input where genomic location of PB reads is unknown

  ptm <- startTimedMessage("Reading the data")
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
  #SSlibNames <- paste(SSreadNames.fields[[2]], SSreadNames.fields[[3]], SSreadNames.fields[[4]], sep="_")
  SSlibNames <- paste(SSreadNames.fields[[2]], SSreadNames.fields[[3]], sep="_")

  #SSflag <- SSreadNames.fields[[5]]
  #SSchrom <- SSreadNames.fields[[6]]
  #SSpos <- SSreadNames.fields[[7]]
  SSflag <- SSreadNames.fields[[4]]
  SSchrom <- SSreadNames.fields[[5]]
  SSpos <- SSreadNames.fields[[6]]

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


#' Filter input data from the minimap
#'
#' This function filters loaded data from the minimap.
#'
#' @param inputData A \code{data.frame} loaded by \code{\link[SaaRclust]{importTestData}}
#' @param quantileSSreads A quantile range for number of SSreads mapped to PB read. (default: 0.4-0.9)
#' @param minSSlibs A range for the minimal and maximal number of StrandS libs being represented per PB read.
#' @return A filtered \code{data.frame}.
#' @importFrom dplyr group_by summarise
#' @author David Porubsky
#' @export


filterInput <- function(inputData=NULL, quantileSSreads=c(0,0.9), minSSlibs=c(20,25)) {

  ptm <- startTimedMessage("Filtering the data")

  #remove SS reads mapped with huge gaps (eg. summed gaps 100bp)
  gaps.perSS.mean <- round(mean(inputData$MatchedBasesWithGaps))
  inputData.filt <- inputData[inputData$MatchedBasesWithGaps <= gaps.perSS.mean,]

  print('summary(inputData.filt):')
  print(summary(inputData.filt))

  #remove SS reads mapped to multiple PB reads
  #rle.SSreads <- rle(inputData.filt$SSreadNames)
  #mask <- rle.SSreads$lengths == 1
  #single.SSreads <- unique(inputData.filt$SSreadNames)[mask]
  ##multi.SSreads <- unique(inputData.filt$SSreadNames)[!mask]
  ##single <- inputData.filt[inputData.filt$SSreadNames %in% single.SSreads,]
  ##multi <- inputData.filt[inputData.filt$SSreadNames %in% multi.SSreads,]
  #inputData.filt <- inputData.filt[inputData.filt$SSreadNames %in% single.SSreads,]

  #get number of SS reads per PB read
  if (!is.null(quantileSSreads)) {
    SSread.perPB <- sort(table(inputData.filt$PBreadNames), decreasing = T) #this won't be needed when output will be already sorted by PBreads

    #filter reads based on the mean number of SS reads aligned per each PB read
    quantile.range <- quantile(SSread.perPB, probs = quantileSSreads)
    filt <- SSread.perPB >= quantile.range[1] & SSread.perPB <= quantile.range[2]
    upperQ.reads <- names(SSread.perPB)[!filt]
    maskNames <- names(SSread.perPB)[filt]
    inputData.filt <- inputData.filt[inputData.filt$PBreadNames %in% maskNames,]
  }

  #filter reads based on the number of SS libs per PB read [FAST]
  inputData.filt %>% dplyr::group_by(PBreadNames) %>% dplyr::summarise(counts = length(unique(SSlibNames))) -> SSlib.perPB.counts
  maskNames <- SSlib.perPB.counts$PBreadNames[SSlib.perPB.counts$counts >= minSSlibs[1] & SSlib.perPB.counts$counts <= minSSlibs[2]]
  inputData.filt <- inputData.filt[inputData.filt$PBreadNames %in% maskNames,]

  #filter reads based on the number of SS libs per PB read [SLOW]
  #SSlib.perPB <- split(as.factor(inputData.filt$SSlibNames), inputData.filt$PBreadNames)
  #SSlib.perPB.counts <- sapply(SSlib.perPB, function(x) length(unique(x)))
  #maskNames <- names(SSlib.perPB.counts)[SSlib.perPB.counts >= minSSlibs]
  #inputData.filt <- inputData.filt[inputData.filt$PBreadNames %in% maskNames,]

  stopTimedMessage(ptm)
  return(list(tab.filt=inputData.filt, upperQ.reads=upperQ.reads))
}


#' Import representative PacBio alignments.
#'
#' This function selects the best representative PacBio alignements from multiple chunks
#' in order to get the best estimates of theta values using Kmeans.
#'
#' @param inputfolder  A folder name where minimap files are stored.
#' @param numAlignments  A required number of representative PacBio alignments.
#' @inheritParams filterInput
#' @return A \code{data.frame}.
#' @author David Porubsky
#' @export

getRepresentativeAlignments <- function(inputfolder=NULL, numAlignments=30000, quantileSSreads=c(0.2,0.8), minSSlibs=c(20,25)) {

  ptm <- startTimedMessage("Getting representative alignments\n")
  file.list <- list.files(path = inputfolder, pattern = "maf.gz$", full.names = TRUE)

  bestAligns <- list()
  countAligns <- 0
  for (file in file.list) {
    filename <- basename(file)
    print(paste('processing file', filename))
    ptm <- proc.time()

    suppressMessages( suppressWarnings( tab.in <- importData(infile=file) ) )
    #suppressMessages( suppressWarnings( tab.in <- importTestData(infile=file, removeDuplicates=TRUE) ) )
    print('filtering input...')
    suppressMessages( tab.filt.l <- filterInput(inputData=tab.in, quantileSSreads=quantileSSreads, minSSlibs=minSSlibs) )
    tab.filt <- tab.filt.l$tab.filt

    numAligns <- length(unique(tab.filt$PBreadNames))

    message("    Obtained ",numAligns," representative alignments", appendLF = F);
    stopTimedMessage(ptm)

    countAligns <- countAligns + numAligns
    bestAligns[[filename]] <- tab.filt
    #interupt the loop in case required amount of representative alignement was reached
    if (!is.null(numAlignments) & countAligns >= numAlignments) {
      break
    }
  }

  bestAligns.tab <- do.call(rbind, bestAligns)
  rownames(bestAligns.tab) <- NULL
  #bestAligns.tab <- bestAligns.tab[sample(nrow(bestAligns.tab)),] #shuffle rows in tab
  sample <- unique(bestAligns.tab$PBreadNames)[1:numAlignments] #get only required amount of representative alignments
  bestAligns.tab <- bestAligns.tab[bestAligns.tab$PBreadNames %in% sample,]

  return(bestAligns.tab)
}


#' Import BAM(s) and count reads
#'
#' Import aligned reads from a multiple BAM files and counts directional reads in specified genomic locations.
#' Results are stored in a \code{list} of matrices with each element of a \code{list} representing counts for single BAM file.
#'
#' @param bamfolder A folder where BAM files to be processed are stored.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param pairedEndReads Set to \code{TRUE} if you have paired-end reads in your file.
#' @param bin.length A length of a bin to count reads in.
#' @return A \code{list} of matrices (columns: minus (W) and plus (C) counts; rows: genomic regions).
#' @importFrom data.table data.table
#' @author David Porubsky
#' @export

importBams <- function(bamfolder=bamfolder, chromosomes=NULL, bin.length=1000000) {
  ## List bams present in a directory
  bamfiles <- list.files(bamfolder, pattern = '.bam$', full.names = T)

  counts.l <- list()
  for (j in 1:length(bamfiles)) {
    bam <- bamfiles[j]
    bam.name <- basename(bam)
    message("Processing ", bam.name)

    ## Load reads from BAM into GRanges object
    #fragments <- bam2GRanges(file = bam, chromosomes=chromosomes, pairedEndReads = T, min.mapq = 10, keep.duplicate.reads = F, what = 'mapq')
    suppressWarnings( fragments <- readBamFileAsGRanges(file=bam, chromosomes=chromosomes, pairedEndReads=TRUE, min.mapq=10, remove.duplicate.reads=TRUE, pair2frgm=FALSE, filtAlt=TRUE) )
    ## Sort fragments by seqlevels
    fragments <- GenomicRanges::sort(fragments, ignore.strand=TRUE)

    if (bin.length) {
      ## Get genome bins
      bins.gr <- unlist( GenomicRanges::tileGenome(seqlengths = seqlengths(fragments), tilewidth = bin.length) )
      hits <- IRanges::findOverlaps(fragments, bins.gr, select = "first") #TODO: make sure the same read can't end up in two neighbouring bins!!!

      #runLength(seqnames(fragments))[which(runValue(seqnames(fragments)) == 'Super-Scaffold_365')]

      ## Add missing values if there is no reads in any given bin
      bin.num <- 1:length(bins.gr)
      missing.bin <- which(!bin.num %in% hits)

      ## Append bin location to the contig name
      if (length(missing.bin) > 0) {
        fragments$ID <- rep(as.character(bins.gr)[-missing.bin], table(hits))
      } else {
        fragments$ID <- rep(as.character(bins.gr), table(hits))
      }

    } else {
      fragments$ID <- seqnames(fragments)
    }

    ## Transform GRanges object into a data.frame
    fragments.df <- as(fragments[,'ID'], 'data.frame')
    fragments.df$ID <- factor(fragments.df$ID, levels=as.character(bins.gr))
    fragments.df$strand <- factor(as.character(fragments.df$strand))

    #table(fragments.df[fragments.df$seqnames == 'Super-Scaffold_365',]$strand)
    #tmp <- fragments.df[fragments.df$seqnames == 'Super-Scaffold_365',]

    ## Count reads
    counts <- data.table::data.table(fragments.df)[,table(strand),by=ID]
    cov.reads <- unique(counts$ID)
    uncov.reads <- levels(cov.reads)[!levels(cov.reads) %in% cov.reads]
    counts <- rbind(matrix(counts$V1, ncol=2, byrow = T), matrix(rep(0, 2*length(uncov.reads)), ncol=2) )
    rownames(counts) <- c(as.character(unique(cov.reads)), uncov.reads)
    counts <- counts[order(match(rownames(counts),levels(fragments.df$ID))),]
    #counts.l[[j]] <- counts
    counts.l[[bam.name]] <- counts  #TEST
  }
  return(counts.l)
}

#'
#' @param bam.files A set of bam files (per single-cell) containing ss to long-read/unitig alignments
#' @return A \code{list} of \code{data.table}s.
#' @author Maryam Ghareghani
#' @export

count.wc.bam <- function(bam.files, max.unitig.cov=2, numCPU=4){
  lib.names <- sapply(bam.files, function(bam) gsub('.mdup.bam$', '', basename(bam)))

  cl <- makeCluster(numCPU)
  doParallel::registerDoParallel(cl)

  parallel.results <- foreach (bam=bam.files, .packages=c('Rsamtools', 'data.table')) %dopar%{
    cat('counting directional reads in', basename(bam), '\n')
    aln = scanBam(file = bam, param=ScanBamParam(what=c('qname','rname','strand'),
                                                 flag=scanBamFlag(isSupplementaryAlignment=F,
                                                                  isSecondaryAlignment = F,
                                                                  isDuplicate = F)))[[1]]

    aln <- data.table(rname=as.character(aln$rname), qname=aln$qname, strand=aln$strand)
    aln <- aln[strand %in% c('+','-')]
    aln[, unitig.cov:=.N, by='qname']
    aln <- aln[unitig.cov <= max.unitig.cov]

    counts <- aln
    counts[strand=='+', c:=.N, by='rname']
    counts[strand=='-', w:=.N, by='rname']

    counts[is.na(counts)] <- 0
    counts[, `:=`(c=max(c), w=max(w)), by='rname']

    counts <- counts[, head(.SD, 1), by='rname', .SDcols=c('c','w')]

    rownames(counts) <- counts[, rname]
    counts[, rname:=NULL]

    list(counts=counts, aln=aln[, .(rname, qname, strand)])
  }

  counts.l <- lapply(parallel.results, function(l) l$counts)
  alignments <- lapply(parallel.results, function(l) l$aln)

  names(counts.l) <- lib.names
  names(alignments) <- lib.names

  stopCluster(cl)

  return(list(counts.l=counts.l, alignments=alignments))
}

# converts counts format between data.table and list

convert.counts.format <- function(counts)
{
  counts.output <- NULL

  if ('data.table' %in% class(counts))
  {
    counts.output <- split(counts, by='lib')

    for (i in 1:length(counts.output)){
      rnames <- counts.output[[i]][, rname]
      counts.output[[i]] <- counts.output[[i]][, .(w,c)]
      rownames(counts.output[[i]]) <- rnames
    }
  } else if (class(counts)=='list') {
    counts.output <- data.table()

    for (i in 1:length(counts)) {
      counts.dt <- counts[[i]]
      counts.dt[, rname:=rownames(counts.dt)]
      counts.dt[, lib:=names(counts)[i]]
      counts.output <- rbind(counts.output, counts.dt)
    }
  } else {
    warning('input type should be list or data.table')
  }

  return(counts.output)
}

#'
#' @param counts.l A \code{list} of \code{data.table}s (rows=long reads/unitigs) per single cell
#' @param num.alignments the number of selected long read/unitig names with the highest ss coverage
#' @param min.ss.cov Minimum required coverage of SS reads in unitigs
#' @param chrom.bar.plot plots the bar plot of the number of long reads/unitigs per chrom/direction. Works with the given chrom.flag input.
#' @param chrom.flag The ground true \code{data.table} containing the chrom/flags of long reads/unitigs
#' @param output.type type of the output counts data. It can be \code{data.table} or \code{list} (if it is splitted by lib names)
#' @return A \code{list} of \code{data.table}s containing the selected alignments
#' @author Maryam Ghareghani
#' @export

get_representative_counts <- function(counts.l, num.alignments, min.ss.cov=30,
                                      chrom.bar.plot=FALSE, chrom.flag=NULL,
                                      output.type='data.table'){
  counts <- data.table()
  for (i in 1:length(counts.l)){
    counts.dt <- counts.l[[i]]
    lib.name=names(counts.l)[i]
    rname=rownames(counts.dt)

    counts.dt[, `:=`(rname=rname, lib=lib.name)]

    counts <- rbind(counts, counts.dt)
  }

  # expanding the counts table to have all possible unitig/libs (synchronizing the set of unitigs for all single-cells)
  rnames <- counts[, unique(rname)]
  lib.names <- counts[, unique(lib)]
  expand.counts <- data.table(expand.grid(rnames, lib.names, stringsAsFactors = F))
  colnames(expand.counts) <- c('rname', 'lib')
  counts <- merge(counts, expand.counts, by=c('rname', 'lib'), all=T)
  counts[is.na(counts)] <- 0

  counts.expand <- split(counts, by='lib')

  # subsetting the highest coverage unitigs
  counts[, ss.cov:=sum(w+c), by=rname]
  setkey(counts, ss.cov)

  counts.selected <- counts[, head(.SD, 1), by=rname]
  counts.selected <- counts.selected[ss.cov >= min.ss.cov]

  selected.rname <- counts.selected[, rname]
  if (length(selected.rname) > num.alignments){
    selected.rname <- tail(selected.rname, num.alignments)
  }

  if (chrom.bar.plot){
    counts.selected.chrom <- merge(counts.selected, chrom.flag, by='rname', all.x=T)
    ggplot(counts.selected.chrom[rname %in% selected.rname])+geom_bar(aes(x=paste(chrom, flag)))
  }

  counts.selected <- counts[rname %in% selected.rname]

  # choosing the set of reads/unitigs with a reasonable fraction of ww/cc strand states ovel all single-cells
  counts.selected[, non.zero.cov:=0][w+c>0, non.zero.cov:=1]
  counts.selected[, non.zero.cov:=sum(non.zero.cov), by=rname]
  counts.selected[, `:=`(n.ww=0, n.cc=0)]
  counts.selected[w/(w+c) < 0.2, n.cc:=1]
  counts.selected[w/(w+c) > 0.8, n.ww:=1]
  counts.selected[, n.ww.cc:=sum(n.ww+n.cc), by=rname]
  filt <- counts.selected[, head(.SD, 1), by=rname]
#  filt <- merge(filt, chrom.flag, by='rname', all.x=T)
  filt[, norm.n.ww.cc:=n.ww.cc/non.zero.cov]

  filt.rname <- filt[norm.n.ww.cc>0.25 & norm.n.ww.cc<0.75, rname]

  counts.selected <- counts.selected[rname %in% filt.rname, .(rname,w,c,lib)]
  # counts.selected <- counts.selected[, .(rname,w,c,lib)]
  #################################################

  counts.selected.l <- split(counts.selected, by='lib')

  for (i in 1:length(counts.selected.l)){
    rownames(counts.selected.l[[i]]) <- counts.selected.l[[i]][, rname]
    counts.selected.l[[i]][, `:=`(rname=NULL, lib=NULL)]

    rownames(counts.expand[[i]]) <- counts.expand[[i]][, rname]
    counts.expand[[i]][, `:=`(rname=NULL, lib=NULL)]
  }

  if (output.type=='data.table') {
    return(list(counts.expand, counts.selected))
  } else if (output.type=='list') {
    return(list(counts.expand, counts.selected.l))
  } else {
    warning('ouput.type should be data.table or list')
    return(NULL)
  }
}

get_gfa_components <- function(edge_file){
  gfa.edges <- fread(edge_file, header = F)
  edges <- as.matrix(gfa.edges[, .(V2, V4)])

  g <- graph_from_edgelist(edges, directed = F)
  #plot(g, vertex.size=1, vertex.label=NA)

  count_components(g)
  comp <- igraph::components(g)
  comp.membership <- data.table(rname=names(comp$membership), component=comp$membership)

  # barplot of graph component sizes
  # barplot(sort(comp$csize))
  comp.membership
}
