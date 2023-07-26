# Functions ---------------------------------------------------------------

# Spoofing for contibait functions to work
spoof_range <- function(x) {
  return(paste0(x, ':0-0'))
}

spoof_rownames <- function(mat) {
  rownames(mat) <- spoof_range(rownames(mat))
  return(mat)
}

strip_range <- function(x) {
  return(stringr::str_remove(x, ':.*$'))
}

bind_with_inverted_unitigs <- function(counts_df) {
  counts_df <-
    counts_df %>% 
    # binned_counts_df %>% 
    # marginalize_wc_counts %>% 
    mutate(n=c+w, unitig_dir = unitig)
  
  counts_df <-
    counts_df %>% 
    bind_rows(
      counts_df %>% 
        mutate(unitig_dir = paste0(unitig, '_inverted'),
               c = n-c,
               w = n-w)
    )
  
  return(counts_df)
}

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
    # alpha <- round(0.1 * 255)
  } else {
    ratio <- counts_1/(counts_1+counts_2)
    rgb_color <- pal_f(ratio)
    # alpha <- round(ratio * 255)
  }
  
  color <- rgb(rgb_color[1], rgb_color[2], rgb_color[3], maxColorValue=255)
  
  return(color)
}


# Command Line ------------------------------------------------------------


## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

# source('/Users/henglinm/Documents/GitHub/strand-seq-graph-phasing/scripts/R/module_utils/phasing_test_args.R')

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--mem-counts',
    '--fastmap-counts',
    '--homology',
    '--connected-components',
    ## Output
    '--output',
    
    ## Params
    '--segment-length-threshold',
    '--expect-XY-separate',
    '--threads'
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

library(contiBAIT)

source(file.path(get_script_dir(), "module_utils/utils.R"))
source(file.path(get_script_dir(), "module_utils/phasing_utils.R"))
# Parsing -----------------------------------------------------------------

## Input

mem_counts <- get_values("--mem-counts", singular=FALSE)
fastmap_counts <- get_values("--fastmap-counts", singular=FALSE)
homology <- get_values("--homology", singular=TRUE)
connected_components <- get_values('--connected-components', singular=TRUE)

## Parameters
segment_length_threshold <- as.numeric(get_values('--segment-length-threshold'))
expect_XY_separate <- as.logical(get_values('--expect-XY-separate', singular=TRUE))
n_threads <- as.numeric(get_values('--threads'))
stopifnot(n_threads >= 1)

## Output
output <- get_values('--output')

# Import ------------------------------------------------------------------

## Connected Components in Exploded Graph ----------------------------------

cat('Reading connected components\n')

components_df <- 
  readr::read_tsv(connected_components)

## Homology ----------------------------------------------------------------

cat('Reading Homology\n')
homology_df <- readr::read_tsv(homology)

cat('No. Bubbles: ', nrow(homology_df)/2, '\n')
## Alignment Counts  ---------------------------------------------------------

counts_df <-
  readr::read_csv(mem_counts, col_types = 'cciiii')

exact_match_counts_df <-
  readr::read_csv(fastmap_counts, col_types = 'cciiii')

stopifnot(setequal(pull_distinct(counts_df, unitig), pull_distinct(exact_match_counts_df, unitig)))

unitig_lengths_df <- 
  counts_df %>% 
  distinct(unitig, length)

long_unitigs_df <-
  unitig_lengths_df %>% 
  filter(length >= segment_length_threshold) %>% 
  distinct(unitig)

counts_df <-
  counts_df %>%
  filter(length >= segment_length_threshold) %>% 
  select(-length) 

exact_match_counts_df <-
  exact_match_counts_df %>% 
  filter(length >= segment_length_threshold) %>% 
  select(-length)


# contiBAIT ---------------------------------------------------------------

## ContiBAIT QC ------------------------------------------------------------

strand.freq <- 
  with(counts_df, make_wc_matrix(w, c, lib, unitig)) %>% 
  spoof_rownames() %>% 
  StrandFreqMatrix()


# getMethod(plotWCdistribution, "StrandFreqMatrix")
# debugonce(plotWCdistribution, signature = "StrandFreqMatrix")
# plotWCdistribution(strand.freq)

# It seems that more libraries can be removed than is explicitly stated during
# contiBAIT preprocessing. EG, it will say 5 libs are removed for insufficent
# read counts, but far fewer than 5 are removed in the strand state matrix.

# UPDATE: It does appear that the additional libraries removed are indeed
# strangely or improperly behaved. EG, a normal library will show three wfrac
# peaks, at -1, 0, 1, while a library that was removed unannounced will have a
# wfrac distribution looking like a normal centered at 0. So I guess there is a
# good reason some additional libraries are excluded, but the fact that this
# isn't properly announced is so frustrating.

cat('Running contiBAIT QC\n')
# debugonce(preprocessStrandTable, signature = 'StrandFreqMatrix')
strand.states <-
  preprocessStrandTable(
    strand.freq,
    filterThreshold = 0.7,
    lowQualThreshold = 0.8,
    minLib = 20
  )

lib_names <- 
  counts_df %>% 
  pull_distinct(lib)

included_libraries <- colnames(strand.states$strandMatrix)
excluded_libraries <- setdiff(lib_names, included_libraries)
WARNINGS <- character()
WARNINGS <- c(WARNINGS, paste0(length(included_libraries)/length(lib_names), '% of Strand-seq libraries pass QC'))

counts_df <-
  counts_df %>%
  filter(lib %in% included_libraries)

exact_match_counts_df <-
  exact_match_counts_df %>%
  filter(lib %in% included_libraries)


## Contibait Chromosome Clustering -----------------------------------------

# TODO incorporate mapping quality as well? Maybe not, filtering based on
# read quality seems to produce weird results with contiBAIT, so maybe
# clustering with read quality is also not a good idea

# weight the unitigs that have more alignments to be more likely to be
# selected earlier by the contiBAIT clustering algorithm.
mean_coverage <- 
  counts_df %>% 
  tidyr::complete(lib, unitig, fill=list(c=0, w=0)) %>% 
  group_by(unitig) %>% 
  summarise(coverage = mean(w+c), .groups="drop") 

weights <-
  mean_coverage %>% 
  mutate(unitig_range = spoof_range(unitig)) %>% 
  with(set_names(coverage, unitig_range))

# arrange weights by order in wfrac matrix
weights <- weights[rownames(strand.states$strandMatrix)]

cat('Running contiBAIT clustering\n')
# getMethod(clusterContigs, "StrandStateMatrix")
# debugonce(clusterContigs, signature = 'StrandStateMatrix')
clust <-
  clusterContigs(strand.states$strandMatrix,
                 recluster = 100,
                 randomWeight = weights,
                 clusterBy = 'hetero',
                 verbose = FALSE)
 
## Detect Haploid Clusters ----------------------------------------------

# debugonce(findSexGroups, signature = c('LinkageGroupList', 'StrandStateMatrix'))
clust <- findSexGroups(clust, strand.states$strandMatrix, callThreshold = 0.333)



if(length(grep('^sex', names(clust), value = TRUE)) > 1) {
  # TODO what if multiple groups of haploid detected chromosomes?
  WARNINGS <- c(WARNINGS, 'More than 1 cluster has been identified as a haploid chromosome cluster')
}


## Enframe -----------------------------------------------------------------

cluster_df <-
  set_names(clust@.Data, clust@names) %>% 
  map(~tibble(unitig = strip_range(.x))) %>% 
  bind_rows(.id = 'cluster')

cluster_df <-
  cluster_df %>% 
  distinct(cluster, unitig)

cluster_df <-
  cluster_df %>% 
  right_join(long_unitigs_df, by='unitig')

# Cluster Refinement ------------------------------------------------------

cat('Refining Clusters\n')




## Remove Micro Clusters ---------------------------------------------------

cat('Removing micro clusters: \n')

cluster_component_fractions <-
  components_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  left_join(unitig_lengths_df, by = 'unitig')

cluster_component_fractions <-
  cluster_component_fractions %>% 
  group_by(component, cluster) %>% 
  summarise(length = sum(length), .groups="drop") %>% 
  group_by(component) %>% 
  mutate(perc_length = length/sum(length)) %>% 
  ungroup()

# Arbitrary threshold
component_fraction_threshold <- 0.02

micro_component_cluster_df <-
  cluster_component_fractions %>%
  filter(!is.na(cluster)) %>%
  # filter(!(cluster %in% par_clusters)) %>% 
  # filter(cluster %in% one_component_clusters) %>% 
  filter(perc_length <= component_fraction_threshold)

micro_component_clusters <-
  micro_component_cluster_df %>% 
  pull_distinct(cluster)

cat('No. micro clusters: ', length(micro_component_clusters), '\n')

micro_component_cluster_unitigs <-
  left_join(components_df, cluster_df, by='unitig') %>% 
  semi_join(micro_component_cluster_df, by = c('component', 'cluster')) %>% 
  pull(unitig)

cluster_df <-
  cluster_df  %>% 
  mutate(cluster = ifelse(unitig %in% micro_component_cluster_unitigs, NA, cluster))



## Cosine Cluster Merging ----------------------------------------------------

cat('Cosine cluster merging\n')

cluster_df <- merge_similar_clusters_on_components(counts_df, cluster_df, components_df, similarity_threshold = 0.5)

cluster_df <- merge_similar_clusters(counts_df, cluster_df, similarity_threshold = 0.66)

### Cosine Unassigned -----------------------------------------------------

cat('Assigning "high-noise" unitigs using cosine similarity \n')



#TODO make this parameter more visible. Explain why only 5
min_n <- 5

#TODO NA handling of values? What if all NA
similarities <- 
  counts_df %>% 
  with(make_wc_matrix(w, c, lib, unitig, min_n=min_n)) %>% 
  cosine_similarity()

# TODO
new_cluster_ix <- 0
any_assigned <- TRUE
while(any_assigned) {
  any_assigned <- FALSE
  
  unclustered_unitigs <-
    cluster_df %>% 
    filter(is.na(cluster)) %>% 
    pull(unitig)
  
  if(length(unclustered_unitigs) == 0) {
    break
  }
  
  candidate_cluster_unitigs <-
    cluster_df %>%
    left_join(components_df, by='unitig') %>%
    with(split(unitig, cluster))
  
  scores <-
    map(candidate_cluster_unitigs, function(x) {
      unclustered_unitigs <- unclustered_unitigs[!(unclustered_unitigs %in% x)]
      apply(abs(similarities[unclustered_unitigs, x, drop=FALSE]), 1, mean, na.rm=TRUE)
    })
  
  # convert to df
  scores_df <-
    scores %>% 
    map(tibble::enframe, name='unitig', value='score') %>% 
    map(~mutate(.x, unitig = as.character(unitig))) %>% 
    bind_rows(.id='cluster')
  
  
  if(all(is.na(scores_df$score))) {
    cat('No valid similarity scores\n')
    break
  }
  
  # TODO what if all NA
  max_score <-
    scores_df %>% 
    filter(score == max(score, na.rm = TRUE)) %>% 
    slice_head(n=1) # break ties
    
  score_thresh <- 0.5
  if(max_score$score > score_thresh) {
    any_assigned <- TRUE
    
    cat(
      'Assigning unitig:',
      max_score$unitig,
      'to cluster:',
      max_score$cluster,
      'similarity score:',
      max_score$score,
      '\n'
    )
    
    cluster_df <-
      cluster_df %>%
      mutate(cluster = ifelse(unitig == max_score$unitig, max_score$cluster, cluster))
  } else {
    cat(
      'No unitig similarity greater than threshold:',
      score_thresh,
      '\n'
    )
    
    # TODO justify
    length_factor <- 3
    
    long_unclustered_unitigs_df <- 
      unitig_lengths_df %>% 
      filter(unitig %in% unclustered_unitigs) %>% 
      filter(length >= length_factor * segment_length_threshold)
    
    
    if(nrow(long_unclustered_unitigs_df) > 0) {
      any_assigned <- TRUE
      new_cluster_ix <- new_cluster_ix + 1
      
      new_cluster_unitig <- slice_max(long_unclustered_unitigs_df, length, n=1)
      cat('creating new cluster with unitig:', new_cluster_unitig$unitig,
          'with length:', new_cluster_unitig$length, 
          'longer than segment length threshold:', segment_length_threshold,
          '\n')
      
      cluster_df <-
        cluster_df %>%
        mutate(cluster = ifelse(unitig == new_cluster_unitig$unitig, paste0('LGcosLong', new_cluster_ix), cluster))
    }
  }
  
  
}





## POCC-----------------------------------------------------------

cat('Propagating one cluster components\n')
cluster_df <- propagate_one_cluster_components(cluster_df, components_df)

## Bare Components ---------------------------------------------------------

cat('Assigning clusters to bare components with homology\n')

# Assign bare components w/ at least one bubble to its own cluster
bare_components <-
  components_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  group_by(component) %>% 
  filter(all(is.na(cluster))) %>% 
  ungroup() %>% 
  pull_distinct(component) 

components_with_whole_bubbles <-
  homology_df %>% 
  left_join(components_df, by='unitig') %>% 
  group_by(bubble) %>% 
  filter(!anyNA(component) & all_are_identical(component)) %>% 
  ungroup() %>% 
  pull_distinct(component) 

# acrocentric_components <-
#   components_df %>% 
#   filter(member_largest_component) %>% 
#   pull_distinct(component)

# TODO Why? There is a recurring pattern where nearly-perfectly phased unitigs
# (straight out of verkko) tend to be unclustered by contibait because the
# degenrate regions add noise to the untiig. Therefore if it is unclustered ->
# likely highly phased out of verkko -> safe to cluster together? Or maybe its
# fine to stay more conservateive if its failing it situations where additional
# phasing isnt even needed?
bare_bubble_components <- 
  intersect(bare_components, components_with_whole_bubbles) # %>%
# setdiff(acrocentric_components)

i <- 0
for(bbc in bare_bubble_components) {
  i <- i+1
  bbc_unitigs <- 
    components_df %>% 
    filter(component == bbc) %>% 
    pull(unitig)
  
  label <- paste0('LG Bare ', i)
  
  cluster_df <-
    cluster_df %>% 
    mutate(cluster = ifelse(unitig %in% bbc_unitigs, label, cluster))
}



## Cosine Cluster Merging ----------------------------------------------------

cat('Cosine cluster merging\n')

# Sometimes, some of the newly created clusters will should be merged into other
# components on cluster (centromere troubles especially)
cluster_df <- merge_similar_clusters_on_components(counts_df, cluster_df, components_df, similarity_threshold = 0.5)
## Link Homology -----------------------------------------------------------

#### Remove multi cluster bubbles --------------------------------------------

multi_cluster_bubbles <-
  homology_df %>%
  left_join(cluster_df, by='unitig') %>% 
  group_by(bubble) %>% 
  filter(all_are_unique(cluster) & !anyNA(cluster)) %>% 
  ungroup() %>% 
  pull_distinct(bubble)

for(bub in multi_cluster_bubbles) {
  bubble_clusters_df <-
    homology_df %>% 
    left_join(cluster_df, by = 'unitig') %>% 
    filter(bubble %in% bub) 
  
  bubble_clusters <-
    bubble_clusters_df %>% 
    pull(cluster) 
  
  warning(
    'clusters: ',
    bubble_clusters[1],
    ' and ',
    bubble_clusters[2],
    ' share homology via unitigs: ',
    paste(bubble_clusters_df$unitig, collapse = ', '),
    '\n'
  )
  
}

homology_df <-
  homology_df %>% 
  filter(!(bubble %in% multi_cluster_bubbles))

#### Link --------------------------------------------------------------------
cat('Linking nodes sharing homology\n')
cluster_df <- link_homology(cluster_df, homology_df, components_df)

## NA Cluster --------------------------------------------------------------
# TODO why do these form?

cluster_df <-
  cluster_df %>%
  mutate(cluster = ifelse(is.na(cluster), paste0('LGNA_', unitig), cluster))


## Unusual Unitigs ---------------------------------------------------------

high_wc_unitigs <-
  strand.states$AWCcontigs@seqnames %>% 
  as.character()

low_wc_unitigs <-
  cluster_df %>% 
  filter(grepl('^sex', cluster)) %>% 
  pull_distinct(unitig) %>% 
  setdiff(high_wc_unitigs)

## PAR Detection -----------------------------------------------------------
cat('Detecting PAR\n')

# TODO add some sort of size based threshold for the PAR (hpc and non-hpc)

# TODO, add a percentage check to declare a PAR cluster. EG the sex cluster must
# take up at least X% of the cluster on a component or something like that.

haploid_component_fractions_df <-
  cluster_df %>%
  left_join(components_df, by='unitig') %>% 
  left_join(unitig_lengths_df, b='unitig') %>% 
  group_by(component) %>% 
  filter(any(grepl('^sex', cluster))) %>% 
  ungroup()

haploid_component_fractions_df <-
  haploid_component_fractions_df %>% 
  group_by(component, cluster) %>% 
  summarise(length = sum(length), .groups="drop") %>% 
  group_by(component) %>% 
  mutate(perc = length/sum(length)) %>% 
  ungroup()

# TODO handle warning when nrow(haploid_component_fractions_df) == 0

# TODO
# Only a haploid component if x% of the             unitigs on a component are haploid clustered?
# Only a haploid component if x% of the *clustered* unitigs on a component are haploid clustered?

perc_threshold <- 0.90

# TODO experiment with filtering criteria that doesn't demand a flag.

# Haploid unititgs sometimes appear in the degenerate regions of other
# chromosomes (i've noticed them in the big circular tangle on Chr1 a couple
# times) so there needs to be a little logic before calling a cluster a PAR
# cluster. Here a genuine haploid unitig is only called if it is the largest
# cluster among clustered unitigs, and occupies more than `perc_threshold`
# percent of a component
haploid_components <-
  haploid_component_fractions_df %>% 
  filter(!is.na(cluster)) %>% 
  group_by(component) %>% 
  filter(perc == max(perc) & grepl('^sex', cluster) & perc >= perc_threshold) %>% 
  ungroup() %>% 
  pull_distinct(component) 

haploid_component_unititgs <-
  components_df %>% 
  filter(component %in% haploid_components) %>% 
  semi_join(long_unitigs_df, by='unitig') %>% 
  pull(unitig)

# expect_XY_separate ~ only perform PAR detection for assemblies where the XY
# component is likely only conntected to the PAR. Hifiasm graphs are more
# tangled than that right now and therefore PAR detection for hifiasm likely
# shoould not be attempted. The condition that the haploid cluster be the
# largest likely would prevent anything from going wrong, but this flag is here
# to be extra safe.
par_clusters <- 
  cluster_df %>% 
  left_join(unitig_lengths_df) %>%
  filter(length < 2.8e6) %>% # size check
  filter(unitig %in% haploid_component_unititgs) %>% 
  filter(!grepl('^sex', cluster)) %>% 
  filter(!is.na(cluster)) %>%
  filter(expect_XY_separate) %>% 
  pull_distinct(cluster) 



cat('No. PAR clusters detected: ', length(par_clusters), '\n')

if(length(par_clusters) > 1) {
  WARNINGS <-
    c(WARNINGS,
      "more than 1 PAR cluster, haven't thought about what will happen in this case")
  
}

### Haploid Propagation------------------------------------------------------

par_unitigs <-
  cluster_df %>% 
  filter(cluster %in% par_clusters) %>% 
  pull(unitig)

cluster_df <-
  cluster_df %>% 
  mutate(cluster = ifelse(cluster %in% par_clusters | grepl('^sex', cluster), 'LGXY', cluster))

# cluster_df <-
#   cluster_df %>%
#   mutate(cluster = ifelse(grepl('^sex', cluster), 'LGXY', cluster))

# Bad XY Bubbles ----------------------------------------------------------

# XY bubbles should only be PAR bubbles, therefore small
bad_xy_bubbles <- 
  homology_df  %>% 
  left_join(cluster_df, by='unitig') %>% 
  left_join(unitig_lengths_df, by='unitig') %>% 
  filter(cluster == 'LGXY') %>% 
  filter(length >= 3e6)
  # filter(length >= 0) 

homology_df <-
  homology_df %>% 
  anti_join(bad_xy_bubbles, by='bubble')

# Orientation Detection w/ Inverted Unitigs -------------------------------

cat('Detecting unitig orientation\n')
# Add inverted version of every unitig to dataset. Guarantees that there will
# unitigs in both orientations when clustering 

# TODO need to double check /rerun the strand orientation clustering. If The
# haploid and PAR clusters are oriented separately ~ chance that they could
# become unsynchronized. Is there a way to fix this? Do it at the haploid clustering step?

## Invert Unitigs ----------------------------------------------------------

## Orientation Detection -------------------------------------------------

# Use WW and CC libraries only for this step? Or does it not really matter? I
# guess the first principal component is picking out all the variation from the
# WW/CC libraries, and filtering to only the WW/CC libraries will therefore only
# have a minimal effect? Do it anyways, it is just a semi_join? This could be a
# vital step for clustering the haploid segments? It appears that for the
# haploid clusters it works well to improve the explained variance of the first
# PC, and makes it comparable to a diploid cluster. However, it seems to still
# work correctly in general regardless of library selection, and it simplifies
# things to not have to filter to any libraries for this step. EG, what if a
# cluster consists entirely of misbehaving unitigs that make library calling
# difficult, and having to account for that special case.

# Is working in a continuous space (eg, principal components of W-fraction)
# instead of a discretized space (eg, k-modes clustering on strand states)
# better for tasks when there are strong prior assumptions?


prcomps <-
  counts_df %>% 
  bind_with_inverted_unitigs() %>%
  left_join(cluster_df, by='unitig') %>% 
  # semi_join(ww_libraries_df, by=c('lib', 'cluster')) %>% 
  split(.$cluster) %>% 
  map(function(x) with(x, make_wc_matrix(w,c,lib,unitig_dir)))  %>% 
  # filling with 0s doesn't seem to affect first PC too much, compared to
  # probabilistic or Bayesian PCA (from pcaMethods bioconductor package)
  map(function(x) {
    x[is.na(x)] <- 0 
    return(x)
  }) %>% 
  map(prcomp)

# It appears that I need to do clustering only on the first PC, which
# corresponds to orientation, as otherwise, other dimensions can influence the
# results and lead to grouping of a unitig and its inversion together.

strand_orientation_clusters_df <-
  map(prcomps, 'x') %>%
  map(function(x)
    tibble(unitig_dir = rownames(x), strand_cluster = sign(x[, 'PC1']))) %>%
  bind_rows(.id = 'cluster')

strand_orientation_clusters_df <-
  strand_orientation_clusters_df %>% 
  mutate(unitig = stringr::str_remove(unitig_dir, '_inverted$'))

# I have discovered that some points are located truly at the origin (for
# all PCs), and therefore do not separate on the first PC. What lol. 

# TODO Figure out what leads to this. (HG03456)

strand_orientation_clusters_df <- 
  strand_orientation_clusters_df %>% 
  mutate(strand_cluster = ifelse(strand_cluster != 0, strand_cluster, ifelse(grepl('_inverted', unitig_dir), -1, 1)))

# Warning that checks that every unitig and its invert are in opposite clusters
bad <-
  strand_orientation_clusters_df %>% 
  group_by(unitig, strand_cluster) %>% 
  filter(n() != 1)

if(nrow(bad) > 0) {
  stop('Unitigs have been clustered with their inversions ~ something is wrong with unitig orientation detection')
  # WARNINGS <- c(WARNINGS, 'Unitigs have been clustered with their inversions ~ something is wrong with unitig orientation detection')
}


# Call WC Libraries -------------------------------------------------------

cat('Calling library strand states\n')
# TODO incorporate the low-quality libraries from the contiBAIT QC into this. Or
# is that done already by using the strand.states$strandMatrix? I think that is
# already handled that's nice lol.




## Call Library State contiBAIT  -------------------------------------------------------

# FIXME This will only call WC states for clusters that contiBAIT created!
strand_state_df <-
  strand.states$strandMatrix@.Data %>% 
  as_tibble(rownames = 'unitig_range')

strand_state_df <-
  strand_state_df %>% 
  tidyr::pivot_longer(-unitig_range, names_to = 'lib', values_to = 'state')

strand_state_df <-
  strand_state_df %>% 
  mutate(unitig = strip_range(unitig_range)) %>% 
  left_join(cluster_df, by='unitig') %>% # Full join to pull in missing clusters?
  left_join(counts_df, by=c('unitig', 'lib'))

# 2 ~ heterozygous call
het_fracs <-
  strand_state_df %>% 
  group_by(cluster, lib, .drop=FALSE) %>% 
  filter(!(unitig %in% c(low_wc_unitigs, high_wc_unitigs))) %>%
  summarise(het_frac = weighted.mean(state==2, n, na.rm=TRUE), .groups="drop") 

# TODO, check that all clusters have at least one unitig to call strand states
# with. Current fall through: If a cluster consists of only unusual unitigs (no
# unitigs in het_fracs) then all libraries will be called WC

# Threshold chosen based on looking at a plot
wc_threshold <- 0.9

ww_libraries_df <-
  het_fracs %>% 
  dplyr::filter(het_frac < wc_threshold) 

wc_libraries_df <-
  het_fracs %>% 
  anti_join(ww_libraries_df, by=c('cluster', 'lib'))


## Call Library State no ContiBait -----------------------------------------

# TODO If a cluster contains no unitigs that passed contibait QC, then the WC
# libraries need to be inferred somehow else.

# Use the prcomps!
missing_clusters <-
  cluster_df %>% 
  pull_distinct(cluster) %>% 
  setdiff(strand_state_df$cluster)

missing_wc <-
  prcomps[missing_clusters] %>% 
  imap_dfr(function(x, nm) {
    # Use loadings from first PC. High loadings ~ WW libraries Low Loadings WC Libraries
    loading_percentages <-
      x$rotation[, 1] ^2
    
    n_libs <- length(loading_percentages)
    
    loading_range <- 
      range(loading_percentages)
    
    threshold <-
      loading_range[1] + 0.2 * loading_range[2]
    
    wc_libs <-
      names(loading_percentages)[loading_percentages < threshold]
    
    # if(length(wc_libs) > n_libs %/% 2) {
    #   wc_libs <-
    #     sort(loading_percentages, decreasing=FALSE)[1:(n_libs %/% 2)] %>% 
    #     names()
    # }
    
    out <-
      tibble(cluster=nm, lib=wc_libs)
    
    return(out)
  })

wc_libraries_df <-
  wc_libraries_df %>% 
  bind_rows(missing_wc)


# Phase Chromosomes -------------------------------------------------------

## Orient Fastmap Counts --------------------------------------------------
# exact_match_counts_df <- raw_exact_match_counts_df

exact_match_counts_df <-
  orient_counts(exact_match_counts_df, strand_orientation_clusters_df)

## Clusters covered by bubbles ---------------------------------------------------

n_bubbles <- 1
libs_per_bubble <- 3
coverage_minimum <- 5


# Mark clusters without bubbles, or without sufficient bubble alignments
clusters_covered_with_bubbles <-
  exact_match_counts_df %>%
  inner_join(homology_df, by='unitig')  %>% 
  left_join(cluster_df, by='unitig') %>% 
  filter(n >= coverage_minimum) %>% 
  group_by(cluster, bubble, lib) %>% 
  filter(n() == 2) %>% 
  ungroup()

clusters_covered_with_bubbles <-
  clusters_covered_with_bubbles %>% 
  group_by(cluster, bubble) %>% 
  filter(length(unique(lib)) >= libs_per_bubble) %>% 
  ungroup()


clusters_covered_with_bubbles <-
  clusters_covered_with_bubbles %>% 
  group_by(cluster) %>% 
  filter(length(unique(bubble)) >= n_bubbles)%>% 
  ungroup()

clusters_covered_with_bubbles <-
  clusters_covered_with_bubbles %>% 
  pull_distinct(cluster)

clusters_covered_with_bubbles <-
  clusters_covered_with_bubbles %>% 
  # setdiff(bad_xy_bubbles) %>% 
  set_names()

clusters_not_covered_with_bubbles <-
  cluster_df %>% 
  pull_distinct(cluster) %>% 
  setdiff(clusters_covered_with_bubbles) %>% 
  set_names()


## Phase Clusters Without Homology -----------------------------------------

### No Homology Counts  ------------------------------------------------

# It appears that for XY clusters, regular counts perform better for separation
# that exact match counts.
non_sex_counts <-
  exact_match_counts_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  filter(cluster %in% clusters_not_covered_with_bubbles) %>% 
  filter(cluster != 'LGXY')

sex_counts <-
  counts_df %>% 
  orient_counts(strand_orientation_clusters_df) %>% 
  left_join(cluster_df, by='unitig') %>% 
  filter(cluster %in% clusters_not_covered_with_bubbles) %>% 
  filter(cluster =='LGXY')

no_homology_counts_df <- 
  bind_rows(non_sex_counts, sex_counts) 

### Sort Unitigs ------------------------------------------------------------

big_unitigs <-
  unitig_lengths_df %>%
  filter(length >= 1e6) %>%
  pull_distinct(unitig)

very_big_unitigs <-
  unitig_lengths_df %>%
  filter(length >= 1e7) %>%
  pull_distinct(unitig)

one_cluster_components <-
  components_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  filter(!is.na(cluster)) %>% 
  distinct(component, cluster) %>% 
  count(component) %>% 
  filter(n == 1) %>% 
  pull(component)

one_cluster_component_unitigs <-
  components_df %>% 
  filter(component %in% one_cluster_components) %>% 
  pull(unitig)

no_bubble_unitig_sorts <- 
  map(clusters_not_covered_with_bubbles, function(x) {
    cat('Sorting no-bubble cluster:', x, '\n')
  
  sort_counts <-
    no_homology_counts_df %>%
    filter(cluster == x) %>%
    semi_join(wc_libraries_df, by=c('lib', 'cluster'))

    # TODO do this more cleverly
    if(nrow(sort_counts) == 0){
      sort_counts <-
        no_homology_counts_df %>%
        filter(cluster == x)
    }
  
    unitigs <-
      sort_counts %>%
      pull_distinct(unitig)
    
    if(length(unitigs) == 1 & all(unitigs %in% one_cluster_component_unitigs)) {
      cat('One unitig cluster on a one-cluster component:', x, '\n')
      out <-
        tibble(sort = 1, unitig = unitigs)
      return(out)
    }
    
    # TODO maybe use && to short circuit? instead of all()
    if(length(unitigs) == 1 & !all(unitigs %in% one_cluster_component_unitigs)) {
      cat('One unitig cluster on a multi-cluster component; not sorting cluster:', x, '\n')
      out <-
        tibble(sort = NA, unitig = sort_counts %>% pull_distinct(unitig))
      return(out)
    }
  
  
    # Trying to sort on good, big unitigs that are often present in verkko graphs.
    unitigs <-
      sort_counts %>%
      pull_distinct(unitig) %>% 
      setdiff(c(low_wc_unitigs, high_wc_unitigs)) %>% 
      intersect(very_big_unitigs)
    
    if(length(unitigs) < 1) {
      unitigs <-
        sort_counts %>%
        pull_distinct(unitig) %>%
        intersect(big_unitigs)
      
      if(length(unitigs) < 1) {
        cat('No adequately long unitigs for homology-free sorting of cluster:', x, '\n')
        out <-
          tibble(sort = NA, unitig = sort_counts %>% pull_distinct(unitig))
        return(out)
      }
      
      # TODO collect and justify parameters
      
      ## TODO Minimum vector length needs to be scaled by num libraries
      wc <-
        sort_counts %>%
        with(make_wc_matrix(w,c,lib,unitig, min_n=1))

      unitig_vector_lengths <-
        wc %>%
        apply( 1, function(x) (sqrt(sum(x^2, na.rm=T))))

      print(unitig_vector_lengths)

      threshold <- 
        max(unitig_vector_lengths) - 0.2 *  max(unitig_vector_lengths)
      
      unitigs <-
        names(unitig_vector_lengths)[unitig_vector_lengths > 4 & unitig_vector_lengths > threshold]
      
      if(length(unitigs) < 1) {
        cat('No adequately low-noise unitigs for homology-free sorting of cluster:', x, '\n')
        out <-
          tibble(sort = NA, unitig = sort_counts %>% pull_distinct(unitig))
        return(out)
      }
    }
    # 
    # 
    # sort_counts <- 
    #   no_homology_counts_df %>%
    #   filter(cluster == x) 
    
   
    
    fake_cluster_df <-
      tibble(unitig = unitigs) %>% 
      mutate(cluster = unitig)

    # TODO justify these parameter values
    out <- 
      merge_similar_clusters2_(sort_counts, fake_cluster_df, similarity_threshold = 0.5, min_n=1, min_overlaps =2) %>% 
      mutate(sort = cluster) %>% 
      mutate(sort = as.integer(as.factor(sort))) 
    
    num_loops <-
      out %>% 
      pull_distinct(num_loops)
    
    if(num_loops == 0) {
      cat('No valid similarity scores for cluster:', x, '\n')
      out <-
        out %>% 
        mutate(sort = NA)
    }
    
    num_clusters <- 
      out %>% 
      pull_distinct(sort) %>% 
      length()
    
    
    if(!(num_clusters %in% c(1,2))) {
      # browser()
      # cat('Hierarchically dividing cluster:', x, '\n')
      # # take two most dissimilar clusters
      # wc <-
      #   sort_counts %>%
      #   with(make_wc_matrix(w,c,lib,unitig, min_n=1))
      # 
      # cs <- 
      #   wc[unitigs, ,drop=FALSE] %>% 
      #   cosine_similarity(min_overlaps=2)
      # 
      # cs[is.na(cs)] <- 0
      # # pc <- prcomp(cs)
      # # get_prcomp_plotdata(pc, 'unitig') %>% ggplot() + geom_label(aes(x=PC1, y=PC2, label=unitig))
      # 
      # clusters <-
      #   # dist(pc$x) %>% 
      #   as.dist(1-cs) %>%
      #   hclust() %>% 
      #   cutree(2)
      cat('Taking largest clusters:', x, '\n')
      
      sorts_to_keep <-
        out %>% 
        left_join(unitig_lengths_df, by='unitig') %>% 
        count(sort, wt=length) %>% 
        arrange(desc(n)) %>% 
        dplyr::slice(1:2)
      
      out <-
        out %>% 
        semi_join(sorts_to_keep, by='sort') %>% 
        mutate(sort = as.integer(as.factor(sort)))
        
    }
    

    # sort_counts %>%
    #   distinct(unitig) %>%
    #   left_join(
    #     out %>%
    #       select(sort, unitig)) %>%
    #   return()
    
    out %>%
      select(sort, unitig) %>%
      return()
  })

### Library Swapping ------------------------------------------------------


# All NA
all_na_sorts <- 
  no_bubble_unitig_sorts %>% 
  keep(function(x) all(is.na(x$sort)))

# Only one Sort:
one_cluster_sorts <-
  no_bubble_unitig_sorts %>% 
  keep(function(x) {
    sorts <- x$sort
    sorts <- sorts[!is.na(sorts)]
    return(length(unique(sorts)) == 1)
    })

two_cluster_sorts <-
  no_bubble_unitig_sorts %>% 
  keep(function(x) {
    sorts <- x$sort
    sorts <- sorts[!is.na(sorts)]
    return(length(unique(sorts)) == 2 & all(sorts %in% c(1,2)))
  })

stopifnot(setequal(c(
  names(all_na_sorts),
  names(one_cluster_sorts),
  names(two_cluster_sorts)
), names(no_bubble_unitig_sorts)))


all_na_swaps <-
  list()
if(length(all_na_sorts) >= 1) {
  
  libs <-
    counts_df %>% 
    pull_distinct(lib) %>% 
    sort()
  
    
  
  all_na_swaps <-
    map(all_na_sorts, function(x) {
      swaps <-
        tibble(lib=libs, swap=NA)%>% 
        tibble::deframe()
      
      return(swaps)
    })
}


# One swaps
one_cluster_swaps <- list()
if(length(one_cluster_sorts) >= 1) {
  
  # sort by length:
  one_cluster_sort_lengths <-
    one_cluster_sorts %>% 
    imap_dbl(function(x, nm){
      
      cluster_size <- 
        cluster_df %>% 
        filter(cluster == nm) %>% 
        left_join(unitig_lengths_df, by='unitig') %>% 
        pull_distinct(length) %>% 
        sum(na.rm=T)
      
      return(cluster_size)
    })
  
  one_cluster_sorts <- one_cluster_sorts[names(sort(one_cluster_sort_lengths))]
  
  one_cluster_swaps <-
    pmap(list(one_cluster_sorts, names(one_cluster_sorts), seq_along(one_cluster_sorts)), function(x, nm, n) {
      
      # TODO add library count check
      
      sort_counts <- 
        no_homology_counts_df %>%
        filter(cluster == nm) %>% 
        inner_join(x, by='unitig') %>% 
        semi_join(wc_libraries_df, by=c('lib', 'cluster'))
      
      # TODO do this more cleverly
      if(nrow(sort_counts) == 0){
        sort_counts <- 
          no_homology_counts_df %>%
          filter(cluster == nm) %>%
          inner_join(x, by='unitig')
      }
      
      # Component check
      # total_length <-
      #   sort_counts %>% 
      #   distinct(unitig) %>% 
      #   left_join(unitig_lengths_df, by='unitig') %>% 
      #   pull(length) %>%
      #   sum()
      
      # cluster_components <- 
      #   sort_counts %>% 
      #   left_join(components_df, by='unitig') %>% 
      #   pull_distinct(component)
      # 
      # clusters_on_cluster_component <-
      #   cluster_df %>% 
      #   left_join(components_df, by='unitig') %>% 
      #   filter(component %in% cluster_components) 
      # 
      # cluster_component_fractions <-
      #   cluster_df %>% 
      #   left_join(components_df, by='unitig') %>% 
      #   left_join(unitig_lengths_df, by='unitig') %>% 
      #   filter(component %in% cluster_components) %>% 
      #   group_by(component, cluster) %>%
      #   summarise(length = sum(length), .groups = 'drop') %>% 
      #   group_by(component) %>% 
      #   mutate(perc = length/sum(length)) %>% 
      #   ungroup()
      #   
      # n_components <- 
      #   cluster_component_fractions %>% 
      #   pull_distinct(component) %>% 
      #   length()
      # 
      # n_majority_components <-
      #   cluster_component_fractions %>% 
      #   filter(cluster == nm) %>% 
      #   filter(perc >= 0.5) %>% 
      #   nrow()
      # 
      # if(n_majority_components < n_components) {
      #   libs <- pull_distinct(sort_counts, lib)
      #   swaps <- 
      #     tibble(lib=libs, swap=NA) %>% 
      #     tibble::deframe()
      #   
      #   return(swaps)
      # }
      
      sort_counts <-
        sort_counts %>% 
        group_by(lib) %>% 
        summarise(c = sum(c), w=sum(w))
      
      if ((n %% 2) == 1) {
        swaps <- set_names(sort_counts$c > sort_counts$w, sort_counts$lib)
      } else {
        swaps <- set_names(sort_counts$c < sort_counts$w, sort_counts$lib)
      }
      
      return(swaps)
    })
}

# Two sort
two_cluster_swaps <- list()
if(length(two_cluster_sorts) >= 1) {
  two_cluster_swaps <-
    imap(two_cluster_sorts, function(x, nm) {
      # browser()
      # TODO can you make this work without wc_libraries df?
      sort_counts <- 
        no_homology_counts_df %>%
        filter(cluster == nm) %>% 
        inner_join(x, by='unitig') %>% 
        semi_join(wc_libraries_df, by=c('lib', 'cluster'))
      
      # TODO do this more cleverly
      if(nrow(sort_counts) == 0){
        sort_counts <- 
          no_homology_counts_df %>%
          filter(cluster == nm) %>%
          inner_join(x, by='unitig')
      }
      
      sort_counts <-
        sort_counts %>% 
        group_by(lib, sort) %>% 
        summarise(c = sum(c), w=sum(w), n=sum(n), .groups = 'drop')
      
      # Sort libraries so that sort 1  has more crick reads and sort 2 has more
      # watson reads
      sort_1_votes <- 
        sort_counts %>% 
        filter(sort == 1) %>% 
        mutate(vote = w-c) %>% 
        with(set_names(vote, lib))
      
      if(all_are_identical(sort_counts$sort)) {
        swaps <- 
          sign(sort_1_votes) == 1
        
        return(swaps)
      }
      
      sort_2_votes <-
        sort_counts %>% 
        filter(sort == 2) %>% 
        mutate(vote = c-w) %>% 
        with(set_names(vote, lib))
      
      swaps <- 
        sign(sort_1_votes + sort_2_votes) == 1
      
      
      return(swaps)
      
    })
}

no_bubble_lib_swaps <-
  c(
    all_na_swaps,
    one_cluster_swaps,
    two_cluster_swaps
  )



### Haplotype Marker Counts -------------------------------------------------

# no_bubble_lib_swaps <- c(no_bubble_lib_swaps, one_unitig_cluster_swaps)

cat('Counting haplotype markers\n')

# Binding an empty table ensures that the output has columns if it turns out
# there were no clusters without homology.
no_bubble_lib_swaps_df <-
  no_bubble_lib_swaps %>%
  map_dfr(function(x)
    tibble::enframe(x, name = 'lib', value = 'swapped'),
    .id = 'cluster') %>% 
  bind_rows(tibble(cluster = character(), lib=character(), swapped=logical()))

no_bubble_marker_counts <-
  no_homology_counts_df %>% 
  inner_join(no_bubble_lib_swaps_df, by=c('lib', 'cluster')) 

no_bubble_marker_counts <-
  no_bubble_marker_counts %>%
  mutate(
    c = ifelse(is.na(swapped), 0, ifelse(swapped, n-c, c)), 
    w = ifelse(is.na(swapped), 0, ifelse(swapped, w-c, w))
  )

no_bubble_marker_counts <-
  no_bubble_marker_counts %>% 
  group_by(unitig) %>% 
  summarise(c = sum(c), w=sum(w), .groups="drop")

## Phase Clusters With Homology -----------------------------------------

### Filter to WC libraries for each cluster ---------------------------------

#TODO I think there is potential here for there to be an issue if there is a
#cluster with covered bubbles but no wc libraries?
bubble_coverage <-
  exact_match_counts_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  filter(cluster %in% clusters_covered_with_bubbles) %>% 
  semi_join(wc_libraries_df, by=c('lib', 'cluster')) %>% 
  # TODO not inner_join here? Need to left join and check that all clusters have
  # homology? Maybe make cluster a factor to handle this or something?
  inner_join(homology_df, by='unitig')

### Count Exact Alignments to Bubbles ---------------------------------------

cat('Counting fastmap alignments to bubbles\n')


bubble_coverage <-
  bubble_coverage %>%
  mutate(crick_coverage_ratio = ifelse(n >= coverage_minimum, c / n, NA)) 

bubble_coverage <-
  bubble_coverage %>%
  mutate(covered_by = map_chr(crick_coverage_ratio, whats_covered))

# This step will also filter out libraries that are evaluated as WC but that
# cover no bubbles
bubble_coverage <-
  bubble_coverage %>%
  filter(!is.na(covered_by))


### Create StrandphaseR Arrays --------------------------------------------
cat('Creating StrandphaseR arrays\n')

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


### Sorting StrandphaseR Arrays ---------------------------------------------
cat('Sorting StrandphaseR arrays\n')

lib_swaps <- 
  strandphaser_arrays %>% 
  map(strandphaser_sort_array)


### Haplotype Marker Counts -------------------------------------------------

cat('Counting haplotype markers\n')

lib_swaps_df <-
  lib_swaps %>%
  map_dfr(function(x)
    tibble::enframe(x, name = 'lib', value = 'swapped'),
    .id = 'cluster')

bubble_marker_counts <-
  exact_match_counts_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  semi_join(wc_libraries_df, by=c('lib', 'cluster')) %>% 
  inner_join(lib_swaps_df, by=c('lib', 'cluster')) 

# NA values indicate WC libraries that mapped to no bubbles, or clusters with no
# bubbles. 

#  treat NA as swapped=FALSE for now?
bubble_marker_counts <-
  bubble_marker_counts %>%
  filter(!is.na(swapped)) %>% 
  mutate(c = ifelse(swapped & !is.na(swapped), n - c, c), 
         w = ifelse(swapped & !is.na(swapped), n - w, w))

bubble_marker_counts <-
  bubble_marker_counts %>% 
  group_by(unitig) %>% 
  summarise(c = sum(c), w=sum(w), .groups="drop")


## Combine Marker Counts ---------------------------------------------------

marker_counts <-
  bubble_marker_counts %>% 
  bind_rows(no_bubble_marker_counts)

# Export ------------------------------------------------------------------

cat('Exporting\n')

## Join in Excluded Unitigs ------------------------------------------------
all_unitigs <- unitig_lengths_df$unitig

short_unitigs <- 
  unitig_lengths_df %>% 
  anti_join(long_unitigs_df, by='unitig') %>% 
  pull(unitig)

short_unitigs_df <-
  tibble(
    unitig = short_unitigs,
    exclusion = paste('Length less than threshold:', segment_length_threshold)
  )

failed_qc_unitigs_df <-
  tibble(unitig = high_wc_unitigs,
         exclusion = 'Too many WC Libraries')

haploid_unitigs_df <-
  tibble(unitig = low_wc_unitigs,
         exclusion = 'Too few WC Libraries')

accounted_unitigs <- 
  c(marker_counts$unitig,
    short_unitigs_df$unitig,
    failed_qc_unitigs_df$unitig,
    haploid_unitigs_df$unitig)

unaccounted_unitigs_df <-
  tibble(unitig = setdiff(all_unitigs, accounted_unitigs),
         exclusion = 'other')

exclusions_df <-
  bind_rows(short_unitigs_df, failed_qc_unitigs_df, haploid_unitigs_df, unaccounted_unitigs_df)

# Make sure no unitig is double listed. This will also catch if a unitig is
# listed in multiple bubbles (which may be desired in some future iteration)
stopifnot(all_are_unique(exclusions_df$unitig))

# Make sure all unitigs are included
stopifnot(setequal(all_unitigs, c(exclusions_df$unitig, marker_counts$unitig)))

marker_counts <-
  full_join(marker_counts, exclusions_df, by='unitig')


## Additional Information --------------------------------------------------

marker_counts <-
  marker_counts %>% 
  left_join(cluster_df, by='unitig') 

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
  mutate(Color = rainbow(n()))

marker_counts <-
  marker_counts %>% 
  left_join(cluster_palette, by='cluster') %>% 
  mutate(Color = pmap_chr(list(Color, hap_1_counts, hap_2_counts), make_bandage_colors))

# fill NA with 0 for rukki?
marker_counts <-
  marker_counts %>%
  mutate(attempted = !is.na(hap_1_counts) & !is.na(hap_2_counts)) %>% 
  mutate(hap_1_counts = ifelse(is.na(hap_1_counts), 0, hap_1_counts),
         hap_2_counts = ifelse(is.na(hap_2_counts), 0, hap_2_counts))

stopifnot(nrow(marker_counts) == length(all_unitigs))
stopifnot(setequal(all_unitigs, marker_counts$unitig))



## Unitig-Confidence Scaling -----------------------------------------------

# This emerged from discussions with Sergey Koren with regards to rukki. It was
# discovered that very large unitigs can have very low unitig counts while still
# appearing to show a plausible phasing. EG utig4-234[01], utig4-2336,
# utig4-1216 from HG02106 show plausible marker counts despite a very low count
# of markers per unitig (<= 3). This could cause problems with rukki where,
# unless certain filtering thresholds were essentially entirely eliminated,
# those large unitigs would not be considered during path finding. After Sergey
# K. consulted with Sergey N., it was suggested that one simple trick would be
# to multiply the marker counts of large unitigs by some scaling factor. This
# would essentially allow large unitigs to be included during path finding
# without removing all filtering. This allows large unitigs to be retained even
# when smaller unitigs are discarded for low marker counts.

marker_counts <-
  marker_counts %>%
  left_join(unitig_lengths_df, by = 'unitig') %>%
  mutate(
    hap_1_counts = ifelse(length >= 1e6, 10 * hap_1_counts, hap_1_counts),
    hap_2_counts = ifelse(length >= 1e6, 10 * hap_2_counts, hap_2_counts)
  )

## CSV ---------------------------------------------------------------------

readr::write_csv(marker_counts, output)

# Warnings ----------------------------------------------------------------


for(w in WARNINGS) {
  warning(w)
}

