# Functions ---------------------------------------------------------------

strip_range <-
  function(x) {
    stringr::str_remove(x, ':.*$')
  }

make_wc_matrix <- function(watson, crick, lib, unitig) {
  stopifnot(length(unique(lengths(list(watson, crick, lib, unitig))))==1)
  
  # factors can cause issue ~ is the integer index or name used? When making
  # dimnames, it appears the name, when filling the matrix, it appears the
  # integer index is used.
  stopifnot(is.character(lib), is.character(unitig)) 
  stopifnot(is.integer(watson), is.integer(crick)) 
  stopifnot(!is_duplicate_pair(lib, unitig))
  
  w_frac <- (watson - crick) / (watson + crick)
  
  dimnames <-
    list(sort(unique(unitig)), sort(unique(lib)))
  
  mat <-
    matrix(
      nrow = length(dimnames[[1]]),
      ncol = length(dimnames[[2]]),
      dimnames = dimnames
    )
  
  for (i in seq_along(w_frac)) {
    mat[unitig[i], lib[i]] <- w_frac[i]
  }
  
  mat[is.na(mat)] <- 0
  
  return(mat)
  
}

extract_exact_matches <- function(fastmap_file) {

  # SQ tags designate a particular read and the read length. 
  # SQ - Read Name - Read Length
  
  # EM tags designate an exact match for a segment of the read
  # EM - Start pos on read - End pos on Read - Number of matches - Matched unitig(s).
  
  # The number of matched unitigs that are printed in the output is controlled
  # by the -w parameter. A '*' in the matched unitigs position indicates that
  # there are more matches than the -w value used in the fastmap call.
  # Otherwise, up to -w matches are printed.
  
  # A read may have multiple exact matches to different unitigs, if different
  # parts of the read each uniquely match to a different unitig. EG: the lines:
  
  # EM
  # SQ ending in *
  # SQ not ending in asterisk
  
  #indicates two different sections of the read aligned, the first section
  #aligned to more than one  location, and therefore ended in *, while a
  #different section of the read had a unique alignment.
  
  lines <-
    readr::read_lines(
      fastmap_file
    )
  
  # There is a weird thing where all SQ tags following the first are preceded
  # by a '//'. No idea why. The last line is also often just '//'
  is_header_line <- grepl('^(//)?SQ', lines)
  is_match_line <-  grepl('^EM', lines)
  
  
  lines <- lines[is_header_line | is_match_line]
  
  # Filter Exact Matches
  groupings <-
    cumsum(grepl('^(//)?SQ', lines)) 
  
  lines <- 
    split(lines, groupings)
  
  lines <-
    lines %>% 
    keep(function(x) length(x) == 2) # one header and one match line
  
  # SQ tags designate a particular read and the read length. 
  # SQ - Read Name - Read Length
  
  # EM tags designate an exact match for a segment of the read
  # EM - Start pos on read - End pos on Read - Number of matches - Matched unitig(s).
  
  # The number of matched unitigs that are printed in the output is controlled
  # by the -w parameter. A '*' in the matched unitigs position indicates that
  # there are more matches than the -w value used in the fastmap call.
  # Otherwise, up to -w matches are printed.
  
  # Dataframes
  
  header <-
    lines %>%
    map_chr(1)
  
  header <-
    readr::read_tsv(I(header),
                    col_names = c('tag', 'qname', 'qlen'),
                    col_types = '_c_')
  
  header <-
    header %>%
    mutate(group = names(lines))
  
  matches <-
    lines %>%
    map_chr(2)
  
  matches <-
    readr::read_tsv(
      I(matches),
      col_names = c('tag', 'lpos', 'rpos', 'n_matches', 'map_info'),
      col_types = '___ic'
    )
  
  matches <-
    matches %>%
    mutate(group = names(lines)) %>% 
    filter(n_matches == 1)
  
  # Tidy map info
  matches <-
    matches %>%
    tidyr::separate(
      col = map_info,
      into = c('unitig', 'more_info'),
      sep = ':'
    ) %>%
    tidyr::separate(
      more_info,
      into = c(
        'strand',
        "I don't know this value's meaning fastmap is poorly documented"
      ),
      sep = 1
    ) %>%
    select(unitig, strand, group)
  
  out <-
    inner_join(header,  matches, by='group') %>% 
    select(-group)
  
  return(out)
  
}

calc_concensus_margin <- function(x, ...) {
  if(all(is.na(x))) {
    return(NA)
  }
  # lower is better
  counts <- table(x, ...)
  return(sum(counts) - max(counts)) # if no names(counts) %in% values, then warning and function returns Inf
}

swap_bubbles <- function(phaser_array, ix) {
  tmp <- phaser_array[ix, , 'watson']
  phaser_array[ix, , 'watson'] <- phaser_array[ix, , 'crick']
  phaser_array[ix, , 'crick'] <- tmp
  
  return(phaser_array)
}

strandphaser_sort_array <- function(phaser_array) {
  # x[lib, bubble, watson/crick] <- unitig
  n_libs <-
    dim(phaser_array)[1]
  
  lib_swapped <- 
    logical(length = n_libs) %>% 
    set_names(dimnames(phaser_array)[[1]])
  
  for(i in seq_len(n_libs)){
    
    concensus_score <-
      phaser_array %>% 
      apply(c(2,3), calc_concensus_margin) %>% 
      sum(na.rm = TRUE)
    
    swapped_concensus_score <-
      swap_bubbles(phaser_array, i) %>% 
      apply(c(2,3), calc_concensus_margin) %>% 
      sum(na.rm = TRUE)
    
    if(swapped_concensus_score < concensus_score) {
      # print(swapped_concensus_score)
      phaser_array <- swap_bubbles(phaser_array, i)
      lib_swapped[i] <- TRUE
    } else{
      # print(concensus_score)
    }
    
  }
  
  return(lib_swapped)
  
  
}


make_bandage_colors <- function(color_hex, counts_1, counts_2) {
  
  stopifnot(length(color_hex) == 1)
  stopifnot(length(counts_1) == 1)
  stopifnot(length(counts_2) == 1)
  
  pal_f <- colorRamp(c(color_hex, 'grey', invert_hex(color_hex))) 
  
  if(is.na(counts_1) | is.na(counts_2)) {
    return('black')
  }
  
  if(counts_1 == 0 & counts_2 == 0) {
    rgb_color <- pal_f(0.5)
    alpha <- round(0.1 * 255)
  } else {
    ratio <- counts_1/(counts_1+counts_2)
    rgb_color <- pal_f(ratio)
    alpha <- round(ratio * 255)
  }
  
  color <- rgb(rgb_color[1], rgb_color[2], rgb_color[3], maxColorValue=255)
  
  return(color)
}


# Command Line ------------------------------------------------------------


## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--mem-alignment-bams',
    '--fastmap-alignments',
    '--homology',
    '--connected-components',
    ## Output
    '--output',
    
    ## Params
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


print(args)

# Library -----------------------------------------------------------------

library(dplyr)
library(purrr)

library(Rsamtools)
library(contiBAIT)

source(file.path(get_script_dir(), "module_utils/utils.R"))

# Parsing -----------------------------------------------------------------

## Input

mem_alignment_files <- get_values("--mem-alignment-bams", singular=FALSE)
fastmap_alignment_files <- get_values("--fastmap-alignments", singular=FALSE)
mashmap <- get_values("--homology", singular=TRUE)
connected_components <- get_values('--connected-components', singular=TRUE)

## Parameters
segment_length_threshold <- as.numeric(get_values('--segment-length-threshold'))

## Output
output <- get_values('--output')

# Import ------------------------------------------------------------------

## Exploded Graph ---------------------------------------------------------

components_df <- readr::read_tsv(connected_components)

## Load Headers ------------------------------------------------------------

# All bam headers should be the same right? Only need one
unitig_lengths <- scanBamHeader(mem_alignment_files[[1]], what='targets')[[1]]$targets
unitig_lengths_df <- tibble::enframe(unitig_lengths, name='unitig', value='length')

long_unitigs_df <- 
  unitig_lengths_df %>% 
  filter(length >= segment_length_threshold) %>% 
  select(unitig)

# Concatenate unitig with its length. Necessary for contiBAIT functions to not break.
contibait_names_df <- 
  unitig_lengths_df %>% 
  mutate(unitig_range = paste0(unitig, ':1-', length)) %>% 
  select(-length)


## Load Mashmap ------------------------------------------------------------
homology_df <- readr::read_delim(mashmap, delim='\t', col_names = c('unitig_1', 'unitig_2'))

# filter out identical pairs, just with unitig_1 and unitig_2 switched
homology_df <-
  homology_df %>%
  filter(!is_duplicate_pair(unitig_1, unitig_2))

# tidy
homology_df <-
  homology_df %>% 
  mutate(bubble = paste0('bubble_', 1:n())) %>% 
  tidyr::pivot_longer(cols = c('unitig_1', 'unitig_2'), names_to = 'bubble_arm', values_to ='unitig') %>% 
  mutate(bubble_arm = stringr::str_replace(bubble_arm, '^unitig_', 'arm_'))

# Remove unititgs that are homologous to more than one other unitig.
homology_df <-
  homology_df %>% 
  group_by(unitig) %>% 
  filter(n() == 1) %>% 
  ungroup()

homology_df <-
  homology_df %>% 
  group_by(bubble) %>% 
  filter(n() == 2) %>% 
  ungroup()

## Count Alignments ---------------------------------------------------------

# library(furrr)
# plan('multisession', workers=8)
# plan('sequential')
### Count mem Alignments ------------------------------------------------------
lib_names <- map_chr(mem_alignment_files, function(x) gsub('.mdup.bam$', '', basename(x)))


counts_df <- map(mem_alignment_files, function(bam){
  cat(paste('counting bwa-mem alignments in', basename(bam), '\n'))
  aln <- scanBam(file = bam,
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
  
  aln <- as_tibble(aln[[1]])
  
  # filter out alignments to short unitigs
  aln <- 
    aln %>%  
    mutate(rname = as.character(rname)) %>% # default is as.factor import?
    filter(rname %in% long_unitigs_df$unitig)
  
  
  # Keep only reads that successfully aligned
  aln <- 
    aln %>% 
    filter(strand %in% c('+','-'))
  
  # filter any ss_reads that map to too many unitigs
  # dplyr::filter with lots of groups can be very slow -> duplicated is faster
  duplicated_qnames <-
    with(aln, qname[duplicated(qname)])
  
  if(any(duplicated_qnames)) {
    aln <-
      aln %>%
      filter(!(qname %in% duplicated_qnames))
  }
  
  # to simplify, only keep alignments where both mates land on the same rname
  aln <- 
    aln %>% 
    filter(rname == mrnm)
  
  # Counting
  out <-
    aln %>%
    # filter(mapq > 0) %>% 
    group_by(rname) %>%
    summarise(c = sum(strand == '+'), w = sum(strand == '-')) %>%
    ungroup()
  
  return(out)
  
  
})


# TODO experiment with filtering on the total/average/something number of
# alignments?
counts_df <-
  counts_df %>%
  set_names(lib_names) %>% 
  bind_rows(.id = 'lib') %>%
  dplyr::rename(unitig = rname) 

counts_df <-
  counts_df %>%
  tidyr::complete(lib, unitig, fill=list(c=0, w=0))

counts_df <-
  counts_df %>% 
  mutate(n = c+w)
 
### Count fastmap Alignments ------------------------------------------------

lib_names <-
  map_chr(fastmap_alignment_files, function(x)
    gsub('_maximal_unique_exact_match.tsv$', '', basename(x)))

# There is a common warning stating "incomplete final line found on". I
# think it is connected to the source of the extra // that appear in the files
exact_match_counts_df <- map(fastmap_alignment_files, function(x) {
  cat(paste('counting bwa-fastmap alignments in', basename(x), '\n'))
  out <- extract_exact_matches(x)
  
  out <-
    out %>%
    filter(unitig %in% long_unitigs_df$unitig) %>% 
    group_by(unitig) %>%
    summarise(c = sum(strand == '+'), w = sum(strand == '-')) %>%
    ungroup()
  
  return(out)
})

exact_match_counts_df <-
  exact_match_counts_df %>%
  set_names(lib_names) %>% 
  bind_rows(.id = 'lib') %>%
  tidyr::complete(lib, unitig, fill=list(c=0, w=0))

exact_match_counts_df <-
  exact_match_counts_df %>% 
  mutate(n = c+w)

# 
# # Component Clustering ----------------------------------------------------
# 
# # As an observation, the chromosome clusters generated purely via contibait
# # clustering w/ strand-seq alignments are not-uncommonly less accurate than
# # those generated by simply assigning each non-acrocentric component to its own
# # cluster. This step essentailly puts a lot of faith in verkko to not
# # missassemble the non-acrocentric chromosomes.
# cluster_df <-
#   components_df %>%
#   semi_join(long_unitigs_df, by = 'unitig') %>%
#   mutate(cluster = ifelse(member_largest_component, NA, paste0('cluster_', component))) %>% 
#   select(unitig, cluster)
# 
# # readr::write_csv(cluster_df, 'cluster_df.csv')
# 

# contiBait Linking ----------------------------------------------------

## W Fraction Matrix -------------------------------------------------------


wfrac.matrix <-
  counts_df %>% 
  left_join(contibait_names_df, by='unitig') %>% 
  with(make_wc_matrix(w, c, lib, unitig_range))


## ContiBAIT QC ------------------------------------------------------------

strand.freq <- StrandFreqMatrix(wfrac.matrix)


# getMethod(plotWCdistribution, "StrandFreqMatrix")
# debugonce(plotWCdistribution, signature = "StrandFreqMatrix")
plotWCdistribution(strand.freq)
# 


# debugonce(preprocessStrandTable, signature = 'StrandFreqMatrix')
strand.states <- preprocessStrandTable(strand.freq, filterThreshold = 0.7)

# filterThreshold = 0.7 is a stricter values than the default (filterThreshold =
# 0.8) because, with larger, higher quality initial assemblies, I expected that
# more unitigs will be more well behaved

# TODO experiment with setting filterThreshold even more strictly (maybe ~ .65,
# .6) and splitting large unitigs that fail QC to see if segments will pass WC

# TODO consider setting filter threshold according to the data. EG, given the
# distribution of aligned counts per unitig, what is the expected variance in
# unitig fraction when the strand state is WC? Set the parameter according to
# that variance.

# TODO The counts for watson and crick following binomial distribution for when
# the state is WC appears to be underdispered. Experiment with overdispersed
# binomials (beta binomials?) to account for things like degeneracy, strand
# state changes, and incorrect assembles?


## Contibait Chromosome Clustering -----------------------------------------

# TODO incorporate  mapping quality as well? Maybe not, filtering based on
# read quality seems to produce weird results with contiBAIT, so maybe
# clustering with read quality is also not a good idea

# weight the unitigs that have more alignments to be more likely to be
# selected earlier by the contiBAIT clustering algorithm.
mean_coverage <- 
  counts_df %>% 
  group_by(unitig) %>% 
  summarise(coverage = mean(w+c))

weights <-
  mean_coverage %>% 
  left_join(contibait_names_df, by='unitig') %>% 
  with(setNames(coverage, unitig_range))

# arrange weights by order in wfrac matrix
weights <- weights[rownames(strand.states$strandMatrix)]

# getMethod(clusterContigs, "StrandStateMatrix")
# debugonce(clusterContigs, signature = 'StrandStateMatrix')
clust <-
  clusterContigs(strand.states$strandMatrix,
                 recluster = 100,
                 randomWeight = weights,
                 clusterBy = 'hetero',
                 verbose = FALSE)
# 
# ## Detect Haploid Chromosomes ----------------------------------------------
# 
# debugonce(findSexGroups, signature = c('LinkageGroupList', 'StrandStateMatrix'))
clust <- findSexGroups(clust, strand.states$strandMatrix)

if(length(grep('^sex', names(clust), value = TRUE)) > 1) {
  # TODO what if multiple groups of haploid detected chromosomes?
  warning('More than 1 cluster has been identified as a haploid chromosome cluster')
}


## Enframe -----------------------------------------------------------------

contibait_cluster_df <-
  set_names(clust@.Data, clust@names) %>% 
  map(~tibble(unitig = strip_range(.x))) %>% 
  bind_rows(.id = 'cluster')

contibait_cluster_df <-
  contibait_cluster_df %>% 
  right_join(long_unitigs_df, by='unitig')

## Label Propagation -------------------------------------------------------

# TODO add a check on the proportion of a component that is clustered? EG. A
# comopnent with only one small node clustered. Can that one node cluster the
# whole component? Or should a component be "mostly clustered" in order too add
# to the cluster in this way?

# TODO, on a related note, if there is a very small cluster attached to a very
# large one, should the largest one take over the smaller ones? Is this better
# handled with some sort of counts threshold, to better control rogue clusters
# in the telomeres and centromeres? ---  I have taken a look and unfortunately a
# countsfilter  is unlikely to work, as some cases will have > 200 alignments
# for many libraries and still cluster separately. 

# TODO What to do about components that are both large and have unassigned
# unitigs? Should I still attempt to assign them somehow? Probably yes

### Remove Micro Clusters ---------------------------------------------------


# clusters that are contained on only one component
# one_component_clusters <-
#   contibait_cluster_df %>% 
#   left_join(components_df, by='unitig') %>% 
#   filter(!is.na(cluster)) %>% 
#   distinct(component, cluster) %>% 
#   count(cluster) %>% 
#   filter(n == 1) %>% 
#   pull(cluster) 

cluster_component_fractions <-
  components_df %>% 
  left_join(contibait_cluster_df, by='unitig') %>% 
  left_join(unitig_lengths_df, by = 'unitig')

cluster_component_fractions <-
  cluster_component_fractions %>% 
  group_by(component, cluster) %>% 
  summarise(length = sum(length)) %>% 
  group_by(component) %>% 
  mutate(perc_length = length/sum(length)) %>% 
  ungroup()

# Arbitrary threshold
component_fraction_threshold <- 0.01

# TODO filter to only components with multiple clusters? What about a
# component consisting of only micro clusters? Worry about this stuff later if
# things start breaking lol. One case maybe to worry about: the PAR taking
# over both X/Y chromosomes, if no haploid chromosomes detected?

micro_component_cluster_df <-
  cluster_component_fractions %>%
  filter(!is.na(cluster)) %>%
  # filter(!(component %in% haploid_components))  %>% # Don't want to accidentally remove the PAR
  # filter(cluster %in% one_component_clusters) %>% 
  filter(perc_length <= component_fraction_threshold)# %>% 
  # pull(cluster)

micro_component_cluster_unitigs <-
  left_join(components_df, contibait_cluster_df, by='unitig') %>% 
  semi_join(micro_component_cluster_df, by = c('component', 'cluster')) %>% 
  pull(unitig)

contibait_cluster_df <-
  contibait_cluster_df  %>% 
  mutate(cluster = ifelse(unitig %in% micro_component_cluster_unitigs, NA, cluster))

### Propagate One-Cluster-Component Clusters --------------------------------

one_cluster_components_df <-
  contibait_cluster_df %>%
  left_join(components_df, by = 'unitig') %>%
  # filter(!(component %in% haploid_components)) %>%
  filter(!is.na(cluster))

one_cluster_components_df <-
  one_cluster_components_df %>%
  distinct(component, cluster) %>%
  group_by(component) %>%
  filter(n() == 1) %>%
  ungroup()

one_cluster_component_unitigs <-
  components_df %>%
  inner_join(one_cluster_components_df, by='component')

# lookup vector
one_cluster_component_unitigs <-
  with(one_cluster_component_unitigs, set_names(cluster, unitig))

contibait_cluster_df <-
  contibait_cluster_df %>%
  mutate(cluster = ifelse(
    unitig %in% names(one_cluster_component_unitigs),
    one_cluster_component_unitigs[unitig],
    cluster
  ))



## Cluster-Component Linking -------------------------------------------------



largest_component_components <- 
  contibait_cluster_df %>% 
  left_join(components_df, by = 'unitig') %>% 
  filter(member_largest_component) %>% 
  # filter(!is.na(cluster)) %>% 
  pull(component) %>% 
  unique()

acrocentric_contibait_clusters <-
  contibait_cluster_df %>% 
  left_join(components_df, by = 'unitig') %>% 
  filter(component %in% largest_component_components) %>% 
  filter(!is.na(cluster)) %>% 
  pull(cluster) %>% 
  unique()

acrocentric_contibait_cluster_components <- 
  components_df %>% 
  left_join(contibait_cluster_df, by = 'unitig') %>% 
  filter(!is.na(cluster)) %>% 
  filter(cluster %in% acrocentric_contibait_clusters) %>% 
  pull(component) %>% 
  unique()

cluster_df <-
  contibait_cluster_df %>% 
  # left_join(
  #   dplyr::rename(contibait_cluster_df, contibait_cluster = cluster),
  #   by = 'unitig'
  # ) %>% 
  left_join(components_df, by='unitig') %>% 
  dplyr::rename(contibait_cluster = cluster) 

cluster_df <-
  cluster_df %>%
  mutate(cluster = ifelse(
    component %in% acrocentric_contibait_cluster_components,
    NA,
    paste0('cluster_', component)
  ))


### Link Non-Acrocentrics ---------------------------------------------------

# one_cluster_components <-
#   one_cluster_components_df %>% 
#   # left_join(components_df, by='component') %>% 
#   # filter(!member_largest_component) %>% 
#   pull(component) %>% 
#   sort()

# one_cluster_components_df <-
#   cluster_df %>%
#   filter(!is.na(contibait_cluster)) %>%
#   distinct(component, contibait_cluster) %>%
#   group_by(component) %>%
#   filter(n() == 1) %>%
#   ungroup()

one_cluster_components <-
  one_cluster_components_df %>% 
  pull(component) %>% 
  setdiff(acrocentric_contibait_cluster_components) %>% 
  sort()

# I fusking hate this linking code it sucks right now
none_merged <- FALSE
while(!none_merged) {
  
  none_merged <- TRUE
  

  
  for(cmp in one_cluster_components) {
    
    # if(cmp == 313) {
    #   break
    # }
    contibait_cluster_to_link <-
      cluster_df %>% 
      filter(component == cmp) %>% 
      distinct(contibait_cluster) %>% 
      filter(!is.na(contibait_cluster)) %>% 
      pull()
    
    stopifnot(length(contibait_cluster_to_link) == 1)# only 1
    
    old_cluster_label <-
      cluster_df %>%
      filter(contibait_cluster == contibait_cluster_to_link) %>% 
      filter(component == cmp) %>% 
      distinct(cluster) %>% 
      pull()
    
    stopifnot(length(old_cluster_label) == 1)
    
    new_cluster_label <-
      cluster_df %>%
      filter(contibait_cluster == contibait_cluster_to_link) %>% 
      filter(component != cmp) %>% 
      filter(cluster != old_cluster_label)
    
    if(nrow(new_cluster_label) > 0) {

      
      new_cluster_label <- 
        new_cluster_label %>% 
        arrange(desc(component)) %>% 
        slice_head(n=1) %>% # arbitrary selection
        pull(cluster)
      

      if(new_cluster_label == old_cluster_label) next
      
      cat('joining cluster', old_cluster_label, 'to cluster', new_cluster_label, '\n' )
      
      cluster_df <-
        cluster_df %>% 
        mutate(cluster = ifelse(cluster == old_cluster_label, new_cluster_label, cluster))
      
      none_merged <- FALSE
    } else {
      # cat('cluster', old_cluster_label, 'stands alone', '\n' )
    }
    
    
  }
}

### Coalesce Acrocentrics------------------------------------------------------

cluster_df <-
  cluster_df %>% 
  mutate(cluster = coalesce(cluster, contibait_cluster)) 


### Clean up ----------------------------------------------------------------

cluster_df <-
  cluster_df %>% 
  select(unitig, cluster)
# 
# 
# # PAR Detection -----------------------------------------------------------
# 
# haploid_components <-
#   cluster_df %>%
#   left_join(components_df, by='unitig') %>% 
#   group_by(component) %>% 
#   filter(any(grepl('^sex', cluster))) %>% 
#   pull(component) %>% 
#   unique()
# 
# haploid_component_unititgs <-
#   components_df %>% 
#   filter(component %in% haploid_components) %>% 
#   pull(unitig)
# 
# par_clusters <- 
#   cluster_df %>% 
#   filter(unitig %in% haploid_component_unititgs) %>% 
#   filter(!grepl('^sex', cluster)) %>% 
#   filter(!is.na(cluster)) %>%
#   pull(cluster) %>% 
#   unique()
# 
# if(length(par_clusters) > 1) {
#   warning("more than 1 PAR cluster, haven't thought about what will happen in this case")
# }
# 
# par_unitigs <-
#   cluster_df %>% 
#   filter(cluster %in% par_clusters) %>% 
#   pull(unitig)

# Homology Linking --------------------------------------------------------

bubble_components <-
  components_df %>%
  left_join(homology_df, by='unitig') %>%
  filter(!is.na(bubble)) %>%
  pull(component)

solo_components <-
  components_df %>%
  count(component) %>%
  filter(n == 1) %>%
  pull(component)

solo_bubble_components <-
  intersect(bubble_components, solo_components)

for(cmp in solo_bubble_components) {
  solo_bubble_unitig <-
    components_df %>%
    filter(component == cmp) %>%
    pull(unitig)

  other_bubble_unitig <-
    homology_df %>%
    group_by(bubble) %>%
    filter(any(unitig == solo_bubble_unitig)) %>%
    filter(unitig != solo_bubble_unitig) %>%
    pull(unitig)

  target_cluster <-
    cluster_df %>%
    filter(unitig == other_bubble_unitig) %>%
    pull(cluster)

  cluster_df <-
    cluster_df %>%
    mutate(cluster = ifelse(unitig == solo_bubble_unitig, target_cluster, cluster))

}


# TODO what about non-isolated bubble unitigs?
# bubbles_across_clusters <-
#   cluster_df %>%
#   left_join(homology_df, by='unitig') %>%
#   group_by(bubble) %>%
#   filter(length(unique(cluster)) == length(cluster)) %>%
#   ungroup() %>%
#   pull(bubble) %>%
#   unique()
#
# for(bub in bubbles_across_clusters) {
#
# }



# Call WC Libraries -------------------------------------------------------

# TODO incorporate the low-quality libraries from the contiBAIT QC into this. Or
# is that done already by using the strand.states$strandMatrix? I think that is
# already handled that's nice lol.

## Unusual Unitigs ---------------------------------------------------------

high_wc_unitigs <-
  strand.states$AWCcontigs@seqnames %>% 
  as.character()

low_wc_unitigs <-
  contibait_cluster_df %>% 
  filter(grepl('^sex', cluster)) %>% 
  pull(unitig)

## Calls -------------------------------------------------------------------

strand_state_df <-
  strand.states$strandMatrix@.Data %>% 
  as_tibble(rownames = 'unitig_range')

strand_state_df <-
  strand_state_df %>% 
  tidyr::pivot_longer(-unitig_range, names_to = 'lib', values_to = 'state')

strand_state_df <-
  strand_state_df %>% 
  mutate(unitig = strip_range(unitig_range)) %>% 
  left_join(cluster_df, by='unitig')

# 2 ~ heterozygous call
het_fracs <-
  strand_state_df %>% 
  group_by(cluster, lib, .drop=FALSE) %>% 
  filter(!(unitig %in% c(low_wc_unitigs, high_wc_unitigs))) %>% 
  summarise(het_frac = mean(state==2)) %>%  
  ungroup()

# TODO, check that all clusters have at least one unitig to call strand states
# with.

# Threshold chosen based on looking at a plot
wc_threshold <- 0.9

wc_libraries_df <-
  het_fracs %>% 
  dplyr::filter(het_frac >= wc_threshold) 

ww_libraries_df <-
  het_fracs %>% 
  anti_join(wc_libraries_df, by=c('cluster', 'lib'))



# ## Propagate PAR -----------------------------------------------------------
# 
# # This step allows the hapoid XY to be phased like a diploid chromosome.
# xy_cluster_label <- 'LGXY'
# 
# if (length(par_clusters) > 0) {
#   cluster_df <-
#     cluster_df %>%
#     mutate(cluster = ifelse(
#       unitig %in% c(haploid_component_unititgs, par_unitigs), xy_cluster_label, cluster))
#   
#   wc_libraries_df <-
#     wc_libraries_df %>% 
#     mutate(cluster = ifelse(cluster %in% par_clusters, xy_cluster_label, cluster))
#   
#   ww_libraries_df <-
#     ww_libraries_df %>% 
#     mutate(cluster = ifelse(cluster %in% par_clusters, xy_cluster_label, cluster))
# }
# 
# 
# 

# 
# 
# # # Assign Rest? ------------------------------------------------------------
# 
# 
# 
# 
# 
# 
# 
# # chr15
# cluster_df %>%
#   filter(unitig == 'utig4-2934')
# 
# cluster_df %>%
#   filter(unitig == 'utig4-2459')
# 
# cluster_df %>%
#   filter(unitig == 'utig4-3364')
# 
# 
# cluster_df %>%
#   filter(unitig == 'utig4-7136')
# 
# 
# chr13 <-
#   cluster_df %>%
#   filter(cluster =='LG16 (57)' | unitig %in% c('utig4-7285', 'utig4-7304', 'utig4-6650'))
# 
# chr14 <-
#   cluster_df %>%
#   filter(cluster =='LG14 (64)')
# 
# chr15 <-
#   cluster_df %>%
#   filter(cluster == 'LG15 (62)' | unitig %in% c('utig4-2021', 'utig4-15', 'utig4-1412', 'utig4-1807'))
# 
# chr21 <-
#   cluster_df %>%
#   filter(cluster =='LG21 (28)'| unitig %in% c('utig4-6651', 'utig4-7047','utig4-2456', 'utig4-2457', 'utig4-2566'))
# 
# chr22 <-
#   cluster_df %>%
#   filter(cluster =='LG22 (14)'  | unitig %in% c('utig4-1814', 'utig4-7447', 'utig4-4170', 'utig4-4715'))
# 
# 
# # using inverted exact matches appears to be promising?
# 
# # wc_df <- 
# #   counts_df %>% 
# #   mutate(wfrac = ifelse(n > 0, (w-c)/(w+c), 0)) 
# 
# 
# wc_vec_df <-
#   counts_df %>% 
#   mutate(wfrac = ifelse(n > 0, (w-c)/(w+c), 0)) 
# 
# # %>% 
# #   group_by(lib, unitig) %>%
# #   arrange(unitig_dir) %>%
# #   summarise(wfrac_vec = wfrac[1] - wfrac[2]) %>% 
# #   ungroup()
# 
# 
# scale_vector <- function(x) {x / sqrt(sum(x^2))}
# 
# wc_vec_df <-
#   wc_vec_df %>% 
#   group_by(unitig_dir) %>% 
#   mutate(wfrac = scale_vector(wfrac)) %>% 
#   ungroup()
# 
# 
# wc_vec_mat <- 
#   wc_vec_df %>% 
#   left_join(components_df) %>%
#   filter(component == 32) %>%
#   select(lib, unitig_dir, wfrac) %>% 
#   tidyr::pivot_wider(names_from = 'lib', values_from = 'wfrac') 
# 
# rnms <- wc_vec_mat$unitig_dir
# 
# wc_vec_mat <-
#   as.matrix(wc_vec_mat[, -1])
# 
# dimnames(wc_vec_mat) <- list(unitig_dir = rnms, lib = colnames(wc_vec_mat))
# 
# 
# cosine_similarity <- function(mat) {
#   sim <- mat / sqrt(rowSums(mat * mat))
#   sim <- sim %*% t(sim)
#   return(sim)
# }
# 
# cosine_distance <- function(mat) {
#   cos_sim <- cosine_similarity(mat)
#   cos_dist <- as.dist(1 - cos_sim)
#   return(cos_dist)
# }
# 
# cos_sim <- cosine_similarity(wc_vec_mat)
# 
# which_unitig <-
#   rownames(wc_vec_mat)
# 
# which_unitig <- which_unitig[!grepl('inverted', which_unitig)]
#   
# cos_sim[which_unitig, which_unitig] %>% 
#   abs() %>% 
#   {as.dist(1-.)} %>% 
# hclust() %>% plot
# pca <-
#   prcomp(abs(cos_sim))
# 
# pca$x[which_unitig, 1:2] %>% 
#   dist() %>% 
#   hclust(method='single') %>% 
#   plot()
# 
# # # 
# # pca <-
# #   counts_df %>%
# #   #counts_df %>%
# #   left_join(components_df) %>%
# #   left_join(cluster_df) %>%
# #   # filter(component == 32 | is.na(cluster) | unitig %in% chr15$unitig) %>%
# #   # filter(lib %in% these_libs) %>%
# #   # filter(!grepl('_inverted', unitig_dir)) %>%
# #   #semi_join(bind_rows(chr21, chr15)) %>%
# #   with(make_wc_matrix(w,c,lib,unitig_dir)) %>%
# #   dist() %>%
# #   prcomp()
# 
# get_prcomp_plotdata(pca, 'unitig_dir') %>%
#   mutate(unitig = stringr::str_remove(unitig_dir, '_inverted')) %>%
#   left_join(cluster_df) %>%
#   left_join(components_df) %>% 
#   # filter(component == 32) %>% 
#   #filter(!(!is.na(cluster) & grepl('inverted', unitig_dir))) %>%
#   filter(!grepl('inverted', unitig_dir)) %>%
#   ggplot(aes(x=PC1, y=PC2, color=cluster)) +
#   geom_point() +
#   geom_text(aes(label = ifelse(is.na(cluster), unitig_dir, NA)), size=3)

## Split and recluster if any low quality? -------------------------------_

# TODO, EG ala makeChrTable(splitBy = 1e6), though this would be tiresome to
# implement...


# Orientation Detection w/ Inverted Unitigs -------------------------------

# Add inverted version of every unitig to dataset. Guarantees that there will
# unitigs in both orientations when clustering 

# TODO need to double check /rerun the strand orientation clustering. If The
# haploid and PAR clusters are oriented separately ~ chance that they could
# become unsynchronized. Is there a way to fix this? Do it at the haploid clustering step?

## Invert Unitigs ----------------------------------------------------------

counts_df <-
  counts_df %>% 
  mutate(unitig_dir = unitig)

counts_df <-
  counts_df %>% 
  bind_rows(
    counts_df %>% 
      mutate(unitig_dir = paste0(unitig, '_inverted'),
             c = n-c,
             w = n-w)
  )

## Orientation Detection -------------------------------------------------

# Use WW and CC libraries only for this step? Or does it not really matter? I
# guess the first principal component is picking out all the variation from the
# WW/CC libraries, and filtering to only the WW/CC libraries will therefore only
# have a minimal effect? Do it anyways, it is just a semi_join? This could be a
# vital step for clustering the haploid segments? It appears that for the
# haploid clusters it works well to improve the explained variance of the first
# PC, and makes it comparable to a diploid cluster.

# Is working in a continuous space (eg, principal components of W-fraction)
# instead of a discretized space (eg, k-modes clustering on strand states)
# better for tasks when there are strong prior assumptions?

prcomps <-
  counts_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  semi_join(ww_libraries_df, by=c('lib', 'cluster')) %>% 
  split(.$cluster) %>% 
  map(function(x) with(x, make_wc_matrix(w,c,lib,unitig_dir))) %>% 
  map(prcomp)


# It appears that I need to do clustering only on the first PC, which
# corresponds to orientation, as otherwise, other PCs can influence the results
# and lead to grouping of a unitig and its inversion together.

strand_orientation_clusters_df <-
  map(prcomps, 'x') %>%
  map(function(x)
    tibble(unitig_dir = rownames(x), strand_cluster = sign(x[, 'PC1']))) %>%
  bind_rows(.id = 'cluster')

strand_orientation_clusters_df <-
  strand_orientation_clusters_df %>% 
  mutate(unitig = stringr::str_remove(unitig_dir, '_inverted$'))

# Warning that checks that every unitig and its invert are in opposite clusters
 bad <-
  strand_orientation_clusters_df %>% 
  group_by(unitig, strand_cluster) %>% 
  filter(n() != 1)

 # TODO, gather warnings in a log somewhere or something?
if(any(bad)) {
  warning('A warning about inversions clustered together or something')
}



# Phase Chromosomes -------------------------------------------------------


## Orient Fastmap Counts --------------------------------------------------
# concatenate with inverted, then filter based on strand orientation

exact_match_counts_df <-
  exact_match_counts_df %>% 
  mutate(unitig_dir = unitig)

exact_match_counts_df <-
  exact_match_counts_df %>% 
  bind_rows(
    exact_match_counts_df %>% 
      mutate(unitig_dir = paste0(unitig, '_inverted'),
             c = n-c,
             w = n-w)
  )

exact_match_counts_df <-
  exact_match_counts_df %>%
  semi_join(
    filter(strand_orientation_clusters_df, strand_cluster == 1),
    by = c('unitig', 'unitig_dir')
  ) 

exact_match_counts_df <-
  exact_match_counts_df %>%
  select(-unitig_dir)


## Filter to WC libraries for each cluster ---------------------------------

bubble_coverage <-
  exact_match_counts_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  semi_join(wc_libraries_df, by=c('lib', 'cluster')) %>% 
  # TODO not inner_join here? Need to left join and check that all clusters have
  # homology? Maybe make cluster a factor to handle this or something?
  inner_join(homology_df, by='unitig')

bubble_coverage <-
  bubble_coverage %>% 
  filter(stringr::str_detect(cluster, '^sex', negate=TRUE))

## Count Exact Alignments to Bubbles ---------------------------------------

coverage_minimum <- 5

bubble_coverage <-
  bubble_coverage %>%
  mutate(crick_coverage_ratio = ifelse(n >= coverage_minimum, c / n, NA)) 

# coverage_ratio_threshold <- 0.75

whats_covered <- function(crick_coverage_ratio, coverage_ratio_threshold=0.75) {
  
  if(is.na(crick_coverage_ratio)) {
    return(NA)
  }
  
  if(crick_coverage_ratio >= coverage_ratio_threshold) {
    return('crick')
  }
  
  if(crick_coverage_ratio <= 1-coverage_ratio_threshold) {
    return('watson')
  }
  
  return(NA)
  
}

bubble_coverage <-
  bubble_coverage %>%
  mutate(covered_by = map_chr(crick_coverage_ratio, whats_covered))

# This step will also filter out libraries that are evaluated as WC but that
# cover no bubbles
bubble_coverage <-
  bubble_coverage %>%
  filter(!is.na(covered_by))


## Create StrandphaseR Arrays --------------------------------------------

strandphaser_arrays <-
  bubble_coverage %>% 
  split(.$cluster) %>% 
  map(function(x) {
    
    dimensions <-
      with(x,
           list(
             lib = sort(unique(lib)),
             bubble = sort(unique(bubble)),
             orientation = c('watson', 'crick')
           ))
    
    bubble_coverage_array <-
      array(dim = map_int(dimensions, length), dimnames = dimensions)
    
    x %>% 
      distinct(lib, bubble, covered_by, unitig) %>% 
      pwalk(function(lib, bubble, covered_by, unitig){
        bubble_coverage_array[lib, bubble, covered_by] <<- unitig
      })
    
    return(bubble_coverage_array)
  })


## Sorting StrandphaseR Arrays ---------------------------------------------

lib_swaps <- 
  strandphaser_arrays %>% 
  map(strandphaser_sort_array)

## Haplotype Marker Counts -------------------------------------------------


lib_swaps_df <-
  lib_swaps %>%
  map_dfr(function(x)
    tibble::enframe(x, name = 'lib', value = 'swapped'),
    .id = 'cluster')

marker_counts <-
  exact_match_counts_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  semi_join(wc_libraries_df, by=c('lib', 'cluster')) %>% 
  left_join(lib_swaps_df, by=c('lib', 'cluster')) 

# NA values indicate WC libraries that mapped to no bubbles, or clusters with no
# bubbles. 

#  treat NA as swapped=FALSE for now?
marker_counts <-
  marker_counts %>%
  filter(!is.na(swapped)) %>% 
  mutate(c = ifelse(swapped & !is.na(swapped), n - c, c), 
         w = ifelse(swapped & !is.na(swapped), n - w, w))

marker_counts <-
  marker_counts %>% 
  group_by(unitig) %>% 
  summarise(c = sum(c), w=sum(w)) %>% 
  ungroup() 


# Export ------------------------------------------------------------------

## Join in Excluded Unitigs ------------------------------------------------
all_unitigs <- unitig_lengths_df$unitig

short_unitigs_df <- 
  unitig_lengths_df %>% 
  anti_join(long_unitigs_df, by='unitig') %>% 
  pull(unitig)

short_unitigs_df <-
  tibble(
    unitig = short_unitigs_df,
    exclusion = paste('Length less than threshold:', segment_length_threshold)
  )

failed_qc_unitigs_df <-
  tibble(unitig = high_wc_unitigs,
         exclusion = 'Too many WC Libraries')

accounted_unitigs <- 
  c(marker_counts$unitig,
    short_unitigs_df$unitig,
    failed_qc_unitigs_df$unitig)

unaccounted_unitigs_df <-
  tibble(unitig = setdiff(all_unitigs, accounted_unitigs),
         exclusion = 'other')

exclusions_df <-
  bind_rows(short_unitigs_df, failed_qc_unitigs_df, unaccounted_unitigs_df)

# Make sure no unitig is double listed
stopifnot(length(unique(exclusions_df$unitig)) == length(exclusions_df$unitig))

# Make sure all unitigs are included
stopifnot(setequal(all_unitigs, c(exclusions_df$unitig, marker_counts$unitig)))

marker_counts <-
  full_join(marker_counts, exclusions_df, by='unitig')


## Additional Information --------------------------------------------------

marker_counts <-
  marker_counts %>% 
  left_join(cluster_df, by='unitig') %>% 
  left_join(dplyr::rename(contibait_cluster_df, contibait_cluster = cluster), by='unitig')

marker_counts <-
  marker_counts %>% 
  left_join(homology_df, by='unitig') %>% 
  select(-bubble_arm)

marker_counts <- 
  marker_counts %>% 
  left_join(filter(strand_orientation_clusters_df, strand_cluster==1), by=c('unitig', 'cluster')) %>% 
  select(-strand_cluster) %>% 
  mutate(unitig_dir = ifelse(grepl('inverted$', unitig_dir), -1, 1)) %>% 
  dplyr::rename(unitig_orientation = unitig_dir)

# first three columns: name, counts_1, counts_2 for rukki
marker_counts <-
  marker_counts %>% 
  dplyr::rename(hap_1_counts = c, hap_2_counts = w) %>% 
  select(unitig, hap_1_counts, hap_2_counts, everything(), exclusion) %>% 
  arrange(cluster, bubble)

## Bandage Cluster Colors -----------------------------------------------------

cluster_palette <-
  marker_counts %>%
  distinct(cluster) %>%
  mutate(Color = rainbow(n(), start=0.7, end=0.1))

marker_counts <-
  marker_counts %>% 
  left_join(cluster_palette, by='cluster') %>% 
  mutate(Color = pmap_chr(list(Color, hap_1_counts, hap_2_counts), make_bandage_colors))


## CSV ---------------------------------------------------------------------

readr::write_csv(marker_counts, output)
