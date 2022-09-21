#' Convert saarclust object into \code{\link{GRanges-class}} object.
#'
#' This function takes saarclust object and for each region reports user defined number of clusters 
#' each regions likely belongs as a \code{\link{GRanges-class}} object.
#'
#' @param saarclust.obj A path to the minimap file to load.
#' @param best.prob Set number of highest probability clusters reported per genomic region.
#' @param prob.th Only regions with this probability and higher will be reported. 
#' @return A \code{\link{GRanges-class}} object with reported cluster ID per genomic region.
#' @importFrom tidyr separate
#' @author David Porubsky
#' @export
#' 
clusters2ranges <- function(saarclust.obj = NULL, best.prob = 1, prob.th = 0.5) {
  ## TODO check if submited object is of class 'saarclust'
  
  soft.prob <- saarclust.obj$soft.pVal
  
  ## Remove segments that do not reach reguired prob.th
  row.max.prob <- apply(soft.prob, 1, max)
  mask <- row.max.prob >= prob.th
  soft.prob <- soft.prob[mask,]
  
  ## Convert rownames into GRanges object 
  df <- data.frame(name=rownames(soft.prob))
  df <- tidyr::separate(data = df, col = name, sep = ":", into = c("name", "region"))
  df <- tidyr::separate(data = df, col = region, sep = "-", into = c("start", "end"))
  regions.gr <- GenomicRanges::GRanges(seqnames=df$name, ranges=IRanges(start=as.numeric(df$start), end=as.numeric(df$end)))
  
  ## Report cluster ID based on a defined number of highest p-values
  if (best.prob == 1) {
    clust.ID <- apply(soft.prob, 1, which.max)
    regions.gr$clust.ID <- clust.ID
  } else {
    # TODO report more than one cluster per region
  }
  
  return(regions.gr)
}


#' Report W & C counts and cluster ID as \code{\link{GRanges-class}} object.
#'
#' This function takes W and C counts togther cluster probabilties and reports user defined number of cluster IDs
#' for each genomic region and every single cell as a \code{\link{GRangesList-class}} object.
#'
#' @param saarclust.obj A path to the minimap file to load.
#' @param best.prob Set number of highest probability clusters reported per genomic region.
#' @param prob.th Only regions with this probability and higher will be reported. 
#' @inheritParams EMclust
#' @return A \code{\link{GRangesList-class}} object with reported cluster ID per genomic region.
#' @importFrom tidyr separate
#' @author David Porubsky
#' @export 
#' 
counts2ranges <- function(counts.l=NULL, saarclust.obj=NULL, best.prob=1, prob.th=0.5) {
  ## TODO check if submited object is of class 'saarclust'
  
  soft.prob <- saarclust.obj$soft.pVal
  
  ## Remove segments that do not reach reguired prob.th
  row.max.prob <- apply(soft.prob, 1, max)
  mask <- row.max.prob >= prob.th
  soft.prob <- soft.prob[mask,]
  
  ## Convert rownames into GRanges object 
  df <- data.frame(name=rownames(soft.prob))
  df <- tidyr::separate(data = df, col = name, sep = ":", into = c("name", "region"))
  df <- tidyr::separate(data = df, col = region, sep = "-", into = c("start", "end"))
  regions.gr <- GenomicRanges::GRanges(seqnames=df$name, ranges=IRanges(start=as.numeric(df$start), end=as.numeric(df$end)))
  
  ## Report cluster ID based on a defined number of highest p-values
  if (best.prob == 1) {
    clust.ID <- apply(soft.prob, 1, which.max)
    regions.gr$clust.ID <- clust.ID
  } else {
    # TODO report more than one cluster per region
  }
  
  ## Go over W and C counts for every cell
  counts.grl <- GenomicRanges::GRangesList()
  for (i in seq_along(counts.l)) {
    counts.cell <- counts.l[[i]]
    ## Remove segments that do not reach reguired prob.th
    counts.cell <- counts.cell[mask,]
    counts.cell.gr <- regions.gr
    counts.cell.gr$C <- counts.cell[,1]
    counts.cell.gr$W <- counts.cell[,2]
    counts.grl[[i]] <- counts.cell.gr 
  }
  return(counts.grl)
}