
# Args --------------------------------------------------------------------


args <- commandArgs(trailingOnly = FALSE)


# Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--bam',
    '--gfa',
    ## Output
    '--output-prefix',
    '--log',
    
    ## Params
    '--threads',
    '--segment-length-threshold'
  )

# Have to handle multiple inputs per tag
arg_idx <- sort(which(args %in% expected_args))
arg_idx <- c(arg_idx, length(args) + 1) # edge case of last tag

get_values <- function(arg, singular=TRUE){
  idx <- which(args == arg)
  stopifnot(length(idx) == 1)
  
  next_idx <- arg_idx[which.max(arg_idx > idx)]
  values <- args[(idx + 1):(next_idx - 1)]
  
  # More than one value? return list. One value ~ remove from list structure. It
  # is probably bad practice to return two different types like that, but most
  # arguments have single values and it makes the code less annoying to look at.
  if(singular) {
    stopifnot(length(values)==1)
    values <- values[[1]]
  } else {
    # stopifnot(length(values)>1)
  }
  
  return(values)
}


get_script_dir <- function() {
  needle <- '--file='
  script_path <- gsub(needle, '',  args[grepl(needle, args)])
  return(dirname(script_path))
}

# Log ---------------------------------------------------------------------

log_path <- get_values('--log', singular=TRUE)
log <- file(log_path, open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

print(args)

# Library -----------------------------------------------------------------

library(dplyr)
library(data.table)
library(Rsamtools)
library(doParallel)
library(matrixStats)
library(contiBAIT)

source(file.path(get_script_dir(), "module_utils/utils.R"))
# Parsing -----------------------------------------------------------------

## Input

input.alignment.files <- get_values("--bam", singular=FALSE)
gfa <- get_values("--gfa", singular=TRUE)


## Parameters
numCPU <- as.numeric(get_values("--threads"))
segment_length_threshold <- as.numeric(get_values('--segment-length-threshold'))

## Output
outputfolder <- get_values('--output-prefix')

parsed_args <-
  list(
    input.alignment.files = input.alignment.files,
    gfa=gfa,
    numCPU = numCPU,
    segment_length_threshold=segment_length_threshold,
    outputfolder = outputfolder
  )

print("Parsed Args:")
print(parsed_args)

# Functions ---------------------------------------------------------------


# get alignment direction fractions from aggregated counts
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

# contiBAIT expects rownames to have a particular format which (I believe) is
# supposed to reflect an alignment to a reference genome. These spoofing
# functions simply format names so that contiBAIT doesn't throw an error.
spoof <- function(x){
  paste0(x , ':0-0')
}

unspoof <- function(x) {
  gsub(':0-0$', '', x)
}

spoof_rownames <- function(m) {
  rownames(m) <- spoof(rownames(m))
  return(m)
}

unspoof_rownames <- function(m) {
  rownames(m) <- unspoof(rownames(m))
  return(m)
}



# Useful Variables --------------------------------------------------------


rlengths  <- read_segment_lengths_from_gfa(gfa)


rnames <- names(rlengths)

rlengths_dt <-
  data.table(rname = rnames, length = rlengths)



# Load Alignments ---------------------------------------------------------


print('counting w/c reads...\n')
lib.names <- sapply(input.alignment.files, function(bam) gsub('.mdup.bam$', '', basename(bam)))

cl <- makeCluster(numCPU)
registerDoParallel(cl)

alignments <- foreach(bam=input.alignment.files, .packages=c('Rsamtools', 'data.table')) %dopar%{
  cat('counting directional reads in', basename(bam), '\n')
  aln = scanBam(file = bam,
                param = ScanBamParam(
                  what = c('qname', 'rname', 'strand', 'pos', 'qwidth', 'mrnm', 'isize', 'mapq'),
                  flag = scanBamFlag(
                    isSupplementaryAlignment = FALSE,
                    isSecondaryAlignment = FALSE,
                    isDuplicate = FALSE,
                    # for the purpose of determining qname/direction mapping, having both mates simply provides redundant information?
                    isFirstMateRead = TRUE,
                    isProperPair = TRUE
                  )
                ))
  
  aln <- as.data.table(aln[[1]])
  # Keep only reads that successfully aligned
  aln <- aln[strand %in% c('+','-')]
  
  # filter any ss_reads that map to too many unitigs
  max.unitig.cov <- 1
  aln[, unitig.cov:=.N, by='qname']
  aln <- aln[unitig.cov <= max.unitig.cov]
  
  # to simplify, only keep alignments where both mates land on the same rname
  aln <- aln[rname == mrnm]
  
  return(aln[, .(rname, qname, strand, mapq)])
}
stopCluster(cl)

names(alignments) <- lib.names



# Count Alignments by Unitig ----------------------------------------------

cl <- makeCluster(numCPU)
doParallel::registerDoParallel(cl)
counts.l.all <-  foreach (aln=alignments, .packages=c('data.table')) %dopar%{
  counts <- aln[, .(c=sum(strand=='+'), w=sum(strand=='-')), by='rname']
  return(counts)
}
stopCluster(cl)

names(counts.l.all) <- lib.names


counts <-
  counts.l.all %>%
  rbindlist(idcol = 'lib') %>%
  tidyr::complete(lib, rname, fill=list(c=0, w=0)) %>%
  as.data.table()



# Filter Segment Length Threshold -----------------------------------------

long_rnames <-
  rnames[rlengths >= segment_length_threshold]

counts <- counts[rname %in% long_rnames] 

wfrac.matrix <-
  get.wfrac.matrix(counts) %>% 
  spoof_rownames()

# ContiBAIT QC ------------------------------------------------------------


strand.freq <- StrandFreqMatrix(wfrac.matrix)

# debugonce(preprocessStrandTable, signature = 'StrandFreqMatrix')
strand.states <- preprocessStrandTable(strand.freq)

# Contibait Chromosome Clustering -----------------------------------------

# TODO incoporporate  mapping quality as well?
# weight the unitigs that have more alignments to be more likely to be
# selected earlier by the contiBAIT clustering algorithm.
weights_dt <-
  counts[, .(coverage = mean(w + c)), by='rname'] %>%
  mutate(rname = spoof(as.character(rname)))

weights <- with(weights_dt, setNames(coverage, rname))

# arrange weights by order in wfrac matrix
weights <- weights[rownames(strand.states$strandMatrix)]

# debugonce(clusterContigs, signature = 'StrandStateMatrix')
cl <-
  clusterContigs(strand.states$strandMatrix,
                 recluster = 100,
                 randomWeight = weights,
                 clusterBy = 'hetero')


# Haploid/Sex Chromosomes ------------------------------------------------

# TODO what if multiple groups of sex detected chromosomes?

# debugonce(findSexGroups, signature = c('LinkageGroupList', 'StrandStateMatrix'))
cl <- findSexGroups(cl, strand.states$strandMatrix)

is_haploid_sex_cluster <- grepl('^sex', names(cl))
if(sum(is_haploid_sex_cluster) > 1) {
  warning('More than 1 cluster has been identified as a haploid chromosome cluster')
}

# If no haploid, need an empty vector when looking at excluded unitigs
haploid_sex_rnames <- c()
if(any(is_haploid_sex_cluster)) {
  
  haploid_sex_rnames <-
    cl[is_haploid_sex_cluster] %>%
    Reduce(c, .) %>%
    unspoof()
  
  # remove haploid clusters from rest of process.
  cl[is_haploid_sex_cluster] <- NULL
}

# Excluded Unitigs --------------------------------------------------------

# At this point, all untigs that will be filtered out, should already be
# filtered out. Can now record reason why each unitig was excluded
short_rnames <-
  setdiff(rnames, long_rnames)

high_wc_rnames <-as.character(strand.states$AWCcontigs@seqnames@values)

low_wc_rnames <- haploid_sex_rnames

excluded_unitigs <-
  list(
    data.table(unitig_name = short_rnames, exclusion_reason = paste('Length less than threshold:', as.character(segment_length_threshold))),
    data.table(unitig_name = high_wc_rnames, exclusion_reason = 'Too many WC SSlib'),
    data.table(unitig_name = low_wc_rnames, exclusion_reason = 'Too few WC SSlib')
  ) %>%
  rbindlist(fill=TRUE) %>% # for when there are some categories without any unitigs
  filter(!is.na(unitig_name))
stopifnot(all(excluded_unitigs[, .N, by='unitig_name']$N == 1))


# Unitig Orientation Detection with Inversion -----------------------------

# Add inverted version of every unitig to dataset. Guarentees that there will
# unitigs in both orientations when clustering
clusters_dt <-
  setNames(cl@.Data, cl@names) %>%
  lapply(function(x) data.table(rname = unspoof(x))) %>%
  rbindlist(idcol='clust') %>%
  inner_join(counts) %>%
  mutate(invert=0)

clusters_dt <-
  clusters_dt %>%
  mutate(n = w+c) %>%
  mutate(root_rname = rname)

# concat with inverted
inverted_clusters_dt <-
  clusters_dt %>%
  mutate(invert=1,
         rname = paste0(root_rname, '_inverted'),
         c = n-c,
         w = n-w)


clusters_dt <-
  rbind(clusters_dt, inverted_clusters_dt)

# Detect orientation within choromomsome clusters in principal components space
prcomps <-
  clusters_dt %>%
  split(by='clust') %>%
  lapply(get.wfrac.matrix) %>% 
  lapply(prcomp, center=FALSE, scale=FALSE)

# Assume each contiBAIT cluster comes from one chromosome
strand_orientation_clusters <-
  lapply(prcomps, `[[`, 'x') %>%
  # lapply(abs) %>%
  lapply(dist) %>%
  lapply(hclust) %>%
  lapply(cutree, 2) %>%
  lapply(dt_from_named_vector, 'strand_clust', 'rname') %>%
  rbindlist(idcol='clust')



# Call WC Libraries -------------------------------------------------------

# No longer using SaaRclust, now just a simple vote

strand_state_dt <-
  strand.states$strandMatrix@.Data %>% 
  as.data.table(x, keep.rownames = TRUE) %>% 
  dplyr::rename(rname = rn) %>% 
  mutate(rname  = unspoof(rname)) 

strand_state_dt <-
  strand_state_dt %>% 
  tidyr::pivot_longer(-rname, names_to = 'library', values_to = 'state')

# to get cluster labels
strand_state_dt <-
  strand_state_dt %>% 
  left_join(clusters_dt, by='rname')

# 2 ~ heterozygous call
het_fracs <-
  strand_state_dt %>% 
  group_by(clust, library) %>% 
  summarise(het_frac = mean(state==2)) %>%  
  ungroup()

# Chosen based on looking at a plot
wc_threshold <- 0.9

wc_libraries <-
  het_fracs %>% 
  dplyr::filter(het_frac >= wc_threshold) 




# Exclusion Check ---------------------------------------------------------
clustered_rnames <- 
  clusters_dt[invert==0, rname] %>% 
  unique()

excluded_rnames <- excluded_unitigs$unitig_name
stopifnot(setequal(rnames, c(clustered_rnames, excluded_rnames)))


# Export Preparation ------------------------------------------------------
# Matching old output 

old_clust <- unique(clusters_dt$clust)

clusters_dict <-
  seq_along(old_clust) %>% 
  setNames(old_clust) 

## Unitig Clustering
# names: rname, first_clust, second_clust
unitig_clusters <- 
  clusters_dt %>%
  filter(invert==0) %>% 
  distinct(rname, clust) %>% 
  left_join(strand_orientation_clusters, by = c("clust", "rname"))

## Cluster ID table
# names: first_clust, second_clust, chrom_clust
new_cluster_ids <-
  unitig_clusters %>% 
  distinct(clust) %>%
  mutate(strand_clust = list(c(1,2))) %>% 
  tidyr::unnest(cols = c(strand_clust)) %>% 
  mutate(first_clust = 1:n()) %>% 
  mutate(second_clust = ifelse(first_clust %% 2 == 0, first_clust-1, first_clust+1)) %>% 
  mutate(chrom_clust = as.integer(as.factor(clust)))

unitig_clusters <-
  unitig_clusters %>% 
  left_join(new_cluster_ids, by = c("clust", "strand_clust")) %>% 
  select(rname, first_clust, second_clust, chrom_clust)

## List of WC libraries for each cluster
##### names: lib, thetawc, clust.forward, clust.backward
temp <-
  new_cluster_ids %>% 
  distinct(clust, chrom_clust, first_clust,second_clust) %>% 
  group_by(clust) %>% 
  dplyr::slice(1) %>% 
  ungroup()

wc_libraries <-
  wc_libraries %>% 
  left_join(temp, by = "clust")  %>% 
  select(library, het_frac, first_clust, second_clust)


## Cluster each SS read aligns to
#### Names: #qname, clust, chrom_clust
clustered_ss_reads <-
  alignments %>% 
  rbindlist(idcol='library') %>% 
  select(rname, qname, strand) %>% 
  left_join(select(unitig_clusters, -chrom_clust), by = "rname") %>% 
  left_join(new_cluster_ids, by = c("first_clust", "second_clust")) 

clustered_ss_reads <-
  clustered_ss_reads %>% 
  mutate(clust = ifelse(strand == '+', first_clust, second_clust)) %>% 
  select(qname, clust, chrom_clust)

clustered_ss_reads <-
  clustered_ss_reads %>% 
  filter(!is.na(clust)) # reads mapping to ecluded untigis 

# Export ------------------------------------------------------------------

#Create a master output directory
outputfolder <- file.path(outputfolder)
dir_create_if_does_not_exist(outputfolder)

# TODO to refactor these output, it appear that one would have to learn/ remove
# the whatshap_split step, as changing the variable names/order in these tables
# can lead to unclear errors, and the scripts is very confusing. Consider
# eliminating that step all together lol.

## Unitig Clustering
# names: rname, first_clust, second_clust, chrom_clust
unitig_clusters %>% 
  dplyr::select(rname, first_clust, chrom_clust) %>% 
  dplyr::rename(`#rname` = rname) %>% 
  fwrite(file=file.path(outputfolder, 'unitig_clusters.tsv'), 
         sep='\t',
         quote = F,
         row.names = F)

## Cluster ID table
# names: first_clust, second_clust, chrom_clust
new_cluster_ids  %>% 
  select(first_clust, second_clust, chrom_clust) %>% 
  filter(first_clust < second_clust) %>% 
  fwrite( 
    file=file.path(outputfolder, 'clust_partners.tsv'), 
    sep='\t',
    quote = F,
    row.names = F)

## List of WC libraries for each cluster
##### names: lib, thetawc, clust.forward, clust.backward
fwrite(wc_libraries, 
       file=file.path(outputfolder, 'wc_libraries.tsv'), 
       sep='\t',
       quote = F,
       row.names = F)


## Cluster each SS read aligns to
#### Names: #qname, clust, chrom_clust
clustered_ss_reads %>% 
  dplyr::rename(`#qname` = qname) %>% 
  fwrite( 
    file=file.path(outputfolder, 'ss_clusters.tsv'), 
    sep='\t',
    quote = F,
    row.names = F)



# excluded unitigs:
excluded_unitigs %>% 
  fwrite(
    file=file.path(outputfolder, 'excluded_unitigs.tsv'), 
    sep='\t',
    quote = F,
    row.names = F
  )




