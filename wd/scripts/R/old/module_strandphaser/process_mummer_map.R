library(data.table)
library(reshape2)
library(seqinr)

output_mummer_map_table <- function(mummer.map.file, bubbles.fasta.file)
{
	bubbles.fasta <- read.fasta(bubbles.fasta.file)
	bubbles <- data.table(name=names(bubbles.fasta))
	bubbles[, bubbleName:=sapply(name, function(x) strsplit(x, '_')[[1]][2])]
	bubbles[, bubbleAllele:=sapply(name, function(x) as.character(as.integer(strsplit(x, '_')[[1]][4])-1))]
	bubbles[, unitig_name:=sapply(name, function(x) strsplit(x, '_')[[1]][8])]
	bubbles[, name:=NULL]

	map <- fread(mummer.map.file, fill=TRUE)
	# checke if next raw startswith '>'
	# copy the next row in new columns
	map[, `:=`(next_V1=data.table::shift(V1, -1), next_next_V1=data.table::shift(V1, -2), next_V2=data.table::shift(V2, -1), next_V3=data.table::shift(V3, -1), next_V4=data.table::shift(V4, -1))]
	# filter out the ss names with no match (next row also an ss read)
	map <- map[V1!=">" | next_V1!=">"]
	# keep only the ss reads with unique match (two rows after that is also ss read)
	map <- map[V1==">" & next_next_V1==">"]
	map[, `:=`(V1=NULL, V4=NULL, next_next_V1=NULL)]
	colnames(map) <- c('SSname_lib', 'isReverseMapped', 'unitig_name', 'bubbleStart', 'SSstart', 'alnLen')
	map[isReverseMapped!="", isReverseMapped:="True"]
	map[isReverseMapped=="", isReverseMapped:="False"]
	map[, SSname:=strsplit(SSname_lib, '_')[[1]][1], by=SSname_lib]
	map[, SSlib:=strsplit(SSname_lib, '_')[[1]][2], by=SSname_lib]
	map[, SSname_lib:=NULL]
	map <- merge(map, bubbles, by="unitig_name", all=T)
	map[is.na(bubbleName), `:=`(bubbleName="None", bubbleAllele="None")]
	map[, ss_rep:=.N, by=SSname]
	map <- map[ss_rep==1]
	map[, ss_rep:=NULL]
	# sort columns
	map <- map[, .(SSname, SSlib, unitig_name, bubbleName, bubbleAllele, isReverseMapped, bubbleStart, SSstart, alnLen)]

	return(map)
}

output_bwa_map_table <- function(bwa.map.file, bubbles.fasta.file)
{
	bubbles.fasta <- read.fasta(bubbles.fasta.file)
	bubbles <- data.table(name=names(bubbles.fasta))
	bubbles[, bubbleName:=sapply(name, function(x) strsplit(x, '_')[[1]][2])]
	bubbles[, bubbleAllele:=sapply(name, function(x) as.character(as.integer(strsplit(x, '_')[[1]][4])-1))]
	bubbles[, unitig_name:=sapply(name, function(x) strsplit(x, '_')[[1]][8])]
	bubbles[, name:=NULL]

	map <- fread(bwa.map.file, fill=TRUE, blank.lines.skip=T)
	map <- map[V1 %in% c("SQ", "EM")]



	# checke if next raw startswith '>'
	# copy the next row in new columns
	map[, `:=`(next_V1=data.table::shift(V1, -1), next_next_V1=data.table::shift(V1, -2), next_V2=data.table::shift(V2, -1), next_V3=data.table::shift(V3, -1), next_V4=data.table::shift(V4, -1))]
	# filter out the ss names with no match (next row also an ss read)
	map <- map[V1!=">" | next_V1!=">"]
	# keep only the ss reads with unique match (two rows after that is also ss read)
	map <- map[V1==">" & next_next_V1==">"]
	map[, `:=`(V1=NULL, V4=NULL, next_next_V1=NULL)]
	colnames(map) <- c('SSname_lib', 'isReverseMapped', 'unitig_name', 'bubbleStart', 'SSstart', 'alnLen')
	map[isReverseMapped!="", isReverseMapped:="True"]
	map[isReverseMapped=="", isReverseMapped:="False"]
	map[, SSname:=strsplit(SSname_lib, '_')[[1]][1], by=SSname_lib]
	map[, SSlib:=strsplit(SSname_lib, '_')[[1]][2], by=SSname_lib]
	map[, SSname_lib:=NULL]
	map <- merge(map, bubbles, by="unitig_name", all=T)
	map[is.na(bubbleName), `:=`(bubbleName="None", bubbleAllele="None")]
	map[, ss_rep:=.N, by=SSname]
	map <- map[ss_rep==1]
	map[, ss_rep:=NULL]
	# sort columns
	map <- map[, .(SSname, SSlib, unitig_name, bubbleName, bubbleAllele, isReverseMapped, bubbleStart, SSstart, alnLen)]

	return(map)
}
# for outputting:
# map <- map[!is.na(bubbleAllele) & !is.na(SSstart)]

output_bubble_allele_coverage_matrix <- function(clusters, wc.cell.clust, ss.clust, map)
{
	colnames(ss.clust) <- c("SSname", "SSclust", "clust.backward")
	ss.clust[, SSname:=strsplit(SSname, '_')[[1]][1], by=SSname]
	clust1 <- clusters[1]

	map <- merge(map, ss.clust, by="SSname")
	
	map[bubbleAllele=='0', num.bubble.al0:=.N, .(SSlib, SSclust, bubbleName)][is.na(num.bubble.al0), num.bubble.al0:=0]
	map[bubbleAllele=='1', num.bubble.al1:=.N, .(SSlib, SSclust, bubbleName)][is.na(num.bubble.al1), num.bubble.al1:=0]
	map[, num.bubble.al0:=max(num.bubble.al0), .(SSlib, SSclust, bubbleName)]
	map[, num.bubble.al1:=max(num.bubble.al1), .(SSlib, SSclust, bubbleName)]


	# have only one row per bubble/SSclust/SSlib
	map <- map[, head(.SD, 1), .(SSlib, SSclust, bubbleName)]
	map[, al0_ratio:=num.bubble.al0/(num.bubble.al0+num.bubble.al1), .(SSlib, SSclust, bubbleName)]

	# set bubbleAllele equal to 2 if the ratio of allele0 coverage is not significantly low or high (no clear haplotype distinction in the bubble/lib)
	map[al0_ratio>0.25 & al0_ratio<0.75, bubbleAllele:='2']

	# keep a subset of columns
	map <- map[, .(SSlib, SSclust, bubbleName, bubbleAllele, clust.backward)]

	# convert long to wide data table (put different ss libs in columns)
	map <- data.table::dcast(map, SSclust+bubbleName+clust.backward~SSlib, value.var="bubbleAllele")

	map[is.na(map)] <- "-"

	# split map by SSclust
	map.sp <- split(map, map$SSclust)

	lapply(map.sp, function(x) x[, `:=`(SSclust=NULL, clust.backward=NULL)])

	if (length(map.sp) != 2){
		print(paste("warning: the size of the cluster pair is", length(map.sp)))
	}

	# make both clusters have the same set of bubbleNames
	map.sp[[1]] <- merge(map.sp[[1]], map.sp[[2]][, .SD, .SDcols="bubbleName"], by="bubbleName", all=T)
	map.sp[[2]] <- merge(map.sp[[2]], map.sp[[1]][, .SD, .SDcols="bubbleName"], by="bubbleName", all=T)

	map.sp[[1]][is.na(map.sp[[1]])] <- "-"
	map.sp[[2]][is.na(map.sp[[2]])] <- "-"

	lapply(map.sp, function(x) setkey(x, bubbleName)) %>% invisible()

	# select a subset of cells that are wc in this cluster pair
	wc.cells <- wc.cell.clust[clust.forward==clust1, lib]
	selected.col.names <- c('bubbleName', wc.cells)

	# This doesn' work because not all wc are mapping to  bubbles?
	map.sp[[1]] <- map.sp[[1]][, intersect(selected.col.names, names(map.sp[[1]])), with=FALSE]
	map.sp[[2]] <- map.sp[[2]][, intersect(selected.col.names, names(map.sp[[1]])), with=FALSE]

	return(map.sp)
}
