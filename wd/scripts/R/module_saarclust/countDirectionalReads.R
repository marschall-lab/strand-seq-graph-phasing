#' Count directional reads.
#'
#' This function expects output from custom minimap test dataset that contains original locations of mapped reads in the genome.
#'
#' @param tab.l A \code{list} of PB alignmetns separated per cell.
#' @return A \code{list} of matrices reporting counts of Watson('-') and Crick('+') reads aligned to each genomic segment per cell. (rows=reads or genomic segments, cols=watson & crick read counts)
#' @importFrom data.table data.table
#' @author David Porubsky
#' @export

countDirectionalReads <- function(tab.l=NULL) {
  ptm <- startTimedMessage("Counting directional reads") 
  ratios.l <- list()
  counts.l <- list()
  for (j in 1:length(tab.l)) {
    lib.name <- names(tab.l[j])
    #message("\tWorking on ",lib.name)
    lib.aligns <- tab.l[[j]]
    
    #count directional reads per PB read (SLOWEST)
    #aligns.per.read <- split(lib.aligns$strand, lib.aligns$PBreadNames)
    #counts <- t(sapply(aligns.per.read, function(x) table(x)))
    
    #count directional reads per PB read (FASTEST)  
    counts <- data.table::data.table(lib.aligns)[,table(strand),by=PBreadNames] #even faster option TEST
    cov.PBreads <- counts$PBreadNames
    uncov.PBreads <- levels(cov.PBreads)[!levels(cov.PBreads) %in% cov.PBreads]
    counts <- rbind(matrix(counts$V1, ncol=2, byrow = T), matrix(rep(0, 2*length(uncov.PBreads)), ncol=2) )
    rownames(counts) <- c(as.character(unique(cov.PBreads)), uncov.PBreads)
    counts <- counts[order(match(rownames(counts),levels(lib.aligns$PBreadNames))),]
    
    #count directional reads per PB read (MEDIUM)
    #counts <- aggregate(strand ~ PBreadNames, data=lib.aligns, FUN=table)
    #cov.PBreads <- counts$PBreadNames
    #uncov.PBreads <- levels(cov.PBreads)[!levels(cov.PBreads) %in% cov.PBreads]
    #counts <- rbind(counts[,2], matrix(rep(0, 2*length(uncov.PBreads)), ncol=2) )
    #rownames(counts) <- c(as.character(cov.PBreads), uncov.PBreads)
    #counts <- counts[order(match(rownames(counts),levels(lib.aligns$PBreadNames))),]
    
    counts.l[[j]] <- counts
    names(counts.l)[j] <- lib.name
  }
  stopTimedMessage(ptm)
  return(counts.l)
}
