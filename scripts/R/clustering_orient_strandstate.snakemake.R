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

unproject_prcomp <- function(x, prc, name = 'lib', value = 'wc_ssf', n_components=2) {
  stopifnot(n_components <= length(x))
  x <- x[seq_len(n_components)]
  x <- uv(x)
  unprojected_x <- prc$rotation[, seq_len(n_components), drop=FALSE] %*% x
  
  colnames(unprojected_x) <- value
  out <-
    as_tibble(unprojected_x, rownames = name)
  
  return(out)
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
    '--connected-components',
    ## Output
    '--output-marker-counts',
    '--output-lib',

    ## Params
    '--segment-length-threshold',
    '--cluster-PAR-with-haploid',
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
connected_components <- get_values('--connected-components', singular=TRUE)

## Parameters
segment_length_threshold <- as.numeric(get_values('--segment-length-threshold'))
cluster_PAR_with_haploid <- as.logical(get_values('--cluster-PAR-with-haploid', singular=TRUE))
n_threads <- as.numeric(get_values('--threads'))
stopifnot(n_threads >= 1)

## Output
output_counts <- get_values('--output-marker-counts')
output_lib <- get_values('--output-lib')

# Import ------------------------------------------------------------------

## Connected Components in Exploded Graph ----------------------------------

cat('Reading connected components\n')

components_df <-
  readr::read_tsv(connected_components)

## Alignment Counts  ---------------------------------------------------------

# TODO, a separate file with unitig lengths
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


# SSF Matrix --------------------------------------------------------------

#TODO make this parameter more visible. Explain why only 10
min_n <- 10

ssf_mat <-
  counts_df %>%
  with(make_wc_matrix(w, c, lib, unitig, min_n=min_n))

# contiBAIT ---------------------------------------------------------------

## ContiBAIT QC ------------------------------------------------------------

strand.freq <-
  spoof_rownames(ssf_mat) %>%
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
    # filterThreshold = 0.7,
    # lowQualThreshold = 0.8,
    minLib = 10
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

ssf_mat <- ssf_mat[, colnames(ssf_mat) %in% included_libraries, drop=FALSE]

## Contibait Chromosome Clustering -----------------------------------------

# TODO incorporate mapping quality as well? Maybe not, filtering based on
# read quality seems to produce weird results with contiBAIT, so maybe
# clustering with read quality is also not a good idea

# weight the unitigs that have more alignments to be more likely to be
# selected earlier by the contiBAIT clustering algorithm.

mean_coverage <-
  counts_df %>%
  tidyr::complete(lib, unitig, fill=list(c=0, w=0)) %>%
  mutate(unitig_range = spoof_range(unitig)) %>% 
  group_by(unitig_range) %>%
  summarise(mean_coverage = mean(w+c), .groups="drop") 

weights <-
  mean_coverage %>%
  with(set_names(log2(mean_coverage), unitig_range))

# arrange weights by order in wfrac matrix
weights <- weights[rownames(strand.states$strandMatrix)]

cat('Running contiBAIT clustering\n')
# getMethod(clusterContigs, "StrandStateMatrix")
# debugonce(clusterContigs, signature = 'StrandStateMatrix')
clust <-
  clusterContigs(strand.states$strandMatrix,
                 # similarityCutoff = 0.8,
                 recluster = 1000,
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
component_fraction_threshold <- 0.15

micro_component_cluster_df <-
  cluster_component_fractions %>%
  filter(!is.na(cluster)) %>%
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



# Cosine Similarity -------------------------------------------------------


# TODO The mean cosine similarities should be weighted by the number of
# shared libraries.

# Why 5? 
min_overlaps <- 5
cosine_similarity_mat <-
  pairwise_complete_cosine_similarity(ssf_mat, min_overlaps = min_overlaps)

abs_cosine_similarity_mat <-
  pairwise_complete_cosine_similarity(abs(ssf_mat), min_overlaps = min_overlaps)

# abs_dp <- 
#   pairwise_complete_dp(abs(ssf_mat), t(abs(ssf_mat)), min_overlaps = min_overlaps)


## Seed Haploid Chromosomes ------------------------------------------------

# Because haploid chromosomes are not parallel with one another in mem counts,
# there is a possibility,eg, that each could be assembled in one contig, not
# assigned a cluster by contiBAIT, and then not clustered in the cluster_unitigs
# step due to the fact that there isn't a guarenteed pair. Accordingly, seed the
# unitigs as clusters with the highest mean absolute SSF

haploid_seed_threshold <- 0.75
# non_na_lib_theshold <- 10

mean_haploid_scores <-
  ssf_mat %>% 
  apply(1, mean_abs, na.rm=TRUE) 

non_na_lib_counts <-
  ssf_mat %>% 
  apply(1, function(x) sum(!is.na(x)))

unitigs_to_seed <-
  names(mean_haploid_scores)[(mean_haploid_scores >= haploid_seed_threshold)  & (non_na_lib_counts >= min_overlaps)]

cluster_df <-
  cluster_df %>% 
  mutate(cluster = ifelse(is.na(cluster) & (unitig %in% unitigs_to_seed), paste0('LGhap_', unitig), cluster))

## Cosine Cluster Merging ----------------------------------------------------

cat('Cosine cluster merging\n')
cluster_df <- merge_similar_clusters_on_components(cosine_similarity_mat, cluster_df, components_df, similarity_threshold = 0.5)
cluster_df <- merge_similar_clusters(cosine_similarity_mat, cluster_df, similarity_threshold = 0.66)

## High Sim Pairing --------------------------------------------------------
#TODO some sort of single linage clustering? I have noticed that often, a noisy
#unitig will be really smiilar to about a third of the unitigs in a cluster, but
#surprisingly dissmilar from the other two thirds of the unitigs in the cluster.
#When using average linkage, the two thirds unitigs will make the unitig appear
#dissimilar from the cluster. However, if you look at like 7 nearest neighbors,
#all those neighbors are way more similar than every other unitig, and come from
#the same cluster. Attempt to harness that somehow?



## Cosine Unassigned -----------------------------------------------------


cat('Assigning "high-noise" unitigs using cosine similarity \n')


cluster_df <-
  cluster_unitigs(
    cosine_similarity_mat,
    cluster_df,
    new_cluster_id = 'LGcos',
    # optimistic_unitigs = optimistic_unitigs
  )

## Max-Based Assignment ----------------------------------------------------

# TODO Assign a unitig to a cluster if it has a very high similarity with any of the individual members?

## POCC-----------------------------------------------------------

cat('Propagating one cluster components\n')
cluster_df <- propagate_one_cluster_components(cluster_df, components_df)

## Cosine Cluster Merging ----------------------------------------------------

optimistic_unitigs <-
  unitig_lengths_df %>% 
  filter(length >= 10e6) %>% 
  pull(unitig)

cluster_df <-
  cluster_df %>% 
  mutate(cluster = ifelse(is.na(cluster) & (unitig %in% optimistic_unitigs), paste0('LGoptim_', unitig), cluster))

cat('Cosine cluster merging\n')
# Sometimes, some of the newly created clusters will should be merged
# into other components on cluster (centromere troubles especially)

cluster_df <- merge_similar_clusters_on_components(cosine_similarity_mat, cluster_df, components_df, similarity_threshold = 0.5)
cluster_df <- merge_similar_clusters(cosine_similarity_mat, cluster_df, similarity_threshold = 0.66)

cluster_df <- merge_similar_clusters_on_components(abs_cosine_similarity_mat, cluster_df, components_df, similarity_threshold = 0.9)
cluster_df <- merge_similar_clusters(abs_cosine_similarity_mat, cluster_df, similarity_threshold = 0.9)

cluster_df <-
  cluster_unitigs(
    abs_cosine_similarity_mat,
    cluster_df,
    cluster_unitig_similarity_threshold = 0.9,
    unitig_unitig_similarity_threshold = 0.9,
    new_cluster_id = 'LGabscos'
  )

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

# cluster_PAR_with_haploid ~ only perform PAR detection for assemblies where the XY
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
  filter(cluster_PAR_with_haploid) %>%
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

# PAR and haploid need to be oriented together.
cluster_df <-
  cluster_df %>%
  mutate(cluster = ifelse(cluster %in% par_clusters | grepl('^sex', cluster), 'LGXY', cluster))

## Small Cluster Removal --------------------------------------------

threshold <- 1e7

cluster_sizes <-
  cluster_df %>%
  left_join(unitig_lengths_df) %>%
  group_by(cluster) %>%
  summarise(length = sum(length), .groups = 'drop')

small_clusters <-
  cluster_sizes %>%
  filter(length < threshold) %>%
  pull(cluster)

cluster_df <-
  cluster_df %>%
  mutate(cluster = ifelse(cluster %in% small_clusters, NA, cluster))

# Orientation Detection w/ Inverted Unitigs -------------------------------

cat('Detecting unitig orientation\n')
# Add inverted version of every unitig to dataset. Guarantees that there will
# unitigs in both orientations when clustering

cluster_cosine_similarities <-
  with(cluster_df, split(unitig, cluster)) %>% 
  map(function(x) cosine_similarity_mat[x, x, drop=FALSE])

# if whole row is NA ~ then it had 0 for all wfracs.
cluster_cosine_similarities <-
  map(cluster_cosine_similarities, function(x) {
    row_ix <-
      apply(x, 1, function(xx) all(is.na(xx)))
    
    col_ix <-
      apply(x, 2, function(xx) all(is.na(xx)))


    x[!row_ix, !col_ix, drop=FALSE]
  })

# concatenate w/inverse
cluster_cosine_similarities <-
  cluster_cosine_similarities %>% 
  map(function(ul) {
    
    lr <- ul
    rownames(lr) <- paste0(rownames(lr), '_inverted')
    colnames(lr) <- paste0(colnames(lr), '_inverted')
    
    ll <- -1 * ul
    rownames(ll) <- paste0(rownames(ll), '_inverted')
    
    ur <- t(ll)
    
    out <- 
      rbind(
        cbind(ul, ur),
        cbind(ll, lr)
      )
    
    return(out)
    
  })

strand_orientation_clusters_df<-
  map(cluster_cosine_similarities, pairwise_complete_hclust_n, n=2, agg_f=mean, na.rm=TRUE) 

strand_orientation_clusters_df <-
  strand_orientation_clusters_df %>% 
  map_dfr(tibble::enframe, 'unitig_dir', 'strand_cluster') %>% 
  mutate(strand_cluster = ifelse(strand_cluster==1, -1, 1))

# Warning that checks that every unitig and its invert are in opposite clusters
bad <-
  strand_orientation_clusters_df %>%
  mutate(unitig = gsub('_inverted', '', unitig_dir)) %>% 
  group_by(unitig, strand_cluster) %>%
  filter(n() != 1) %>% 
  filter(!all(is.na(strand_cluster)))

if(nrow(bad) > 0) {
  stop('Unitigs have been clustered with their inversions ~ something is wrong with unitig orientation detection')
  # WARNINGS <- c(WARNINGS, 'Unitigs have been clustered with their inversions ~ something is wrong with unitig orientation detection')
}

# Add back in excluded unitigs
strand_orientation_clusters_df <-
  left_join(long_unitigs_df, strand_orientation_clusters_df,
             by = c('unitig' = 'unitig_dir')
  )



# Phase Chromosomes -------------------------------------------------------


## Orient Counts -----------------------------------------------------------

exact_match_counts_df <-
  orient_counts(exact_match_counts_df, strand_orientation_clusters_df)

## Principal Components ----------------------------------------------------

# does this work with, eg, a 1 unitig cluster? Or does the dimension projection
# lost a lot of library information?
em_counts <-
  exact_match_counts_df %>%
  left_join(cluster_df, by='unitig') %>% 
  arrange(cluster, unitig, lib) %>%
  split(.$cluster)

em_ssf_mats <- 
  map(em_counts, function(x) {
    with(x, make_wc_matrix(w, c, lib, unitig, min_n=4)) 
  })

# Concatenate with inverted unitigs
em_ssf_mats_inverted <-
  map(em_ssf_mats, `*`, -1) %>% 
  map(function(x) {
    rownames(x) <- paste0(rownames(x), '_inverted')
    return(x)
  })

em_ssf_mats <- map2(em_ssf_mats, em_ssf_mats_inverted, rbind)

# 
# # Projected
# em_ssf_mats_proj <-
#   map2(em_ssf_mats, ww_vectors, function(x, ww_df) {
#     ww_df <-
#       ww_df %>%
#       arrange(lib)
#     
#     libs <-
#       ww_df %>%
#       pull(lib)
#     
#     stopifnot(all(colnames(x) == libs))
#     
#     ssfs <-
#       ww_df %>%
#       pull(ww_ssf)
#     
#     # apply to each row, because of need for pairwise complete case management.
#     # Additionally transpose at end, because, even though the function is
#     # applied by row, it fills the results by column...
#     x_proj <-
#       apply(x, 1, project_through, y=ssfs) %>%
#       t()
#     
#     rownames(x_proj) <- paste0(rownames(x_proj), '_projected')
#     
#     return(x_proj)
#   })



em_ssf_mats <-
  em_ssf_mats %>% 
  # filling with 0s doesn't seem to affect first PC too much, compared to
  # probabilistic or Bayesian PCA (from pcaMethods bioconductor package)
  map(function(x) {
    x[is.na(x)] <- 0 
    return(x)
  }) 

prcomps <-
  em_ssf_mats %>% 
  map(prcomp)

## Calculate WC and WW Basis Vectors ----------------------------------------


### Count and Rotate --------------------------------------------------------
oriented_counts_df <-
  orient_counts(counts_df, strand_orientation_clusters_df)

oriented_counts_df <-
  oriented_counts_df %>%  
  left_join(cluster_df) %>% 
  arrange(cluster, unitig, lib) 

oriented_counts_df <-
  oriented_counts_df %>%
  group_by(cluster, lib) %>%
  summarise(across(c(c, w, n), sum), .groups = 'drop')

ww_vectors <-
  oriented_counts_df %>%
  mutate(ww_ssf = ifelse(n < 1, NA, (w-c)/n))  %>% 
  select(cluster, lib, ww_ssf) %>% 
  split(.$cluster)

ww_vectors <- ww_vectors[names(em_counts)]

# Project WW vectors into prcomp space, to make orthognalization, which is
# equivalent to calculation of WC vector, easy.
projected_ww_vectors <- 
  map2(prcomps, ww_vectors, function(x, y){
    ssf_mat <- matrix(y$ww_ssf, nrow=1)
    colnames(ssf_mat) <- y$lib
    ssf_mat[is.na(ssf_mat)] <- 0 
    
    predict(x, ssf_mat) 
  })


projected_wc_vectors <- 
  map(projected_ww_vectors, function(x){
    wc_vec <- uv(c(x[, 2], -x[, 1]))
    wc_vec <- set_names(wc_vec, c('PC1', 'PC2'))
    return(wc_vec)
  })

# To go from PC space back to library space for the WC vector, project the first
# two PCS onto the rotated WW vector.
unprojected_wc_vectors <-
  map2(projected_wc_vectors, prcomps, unproject_prcomp, n_components=2, value = 'wc_ssf')

# Should the WW vector be recalculated by unprojecting the estimate based on the
# first two principal components? As a form of "noise correction?"
unprojected_ww_vectors <-
  map2(projected_ww_vectors, prcomps, unproject_prcomp, n_components=2, value='ww_ssf')

#### Range Balanced Estimate -------------------------------------------------

# Range-Balancing Correction for haplotype clusters of different size: With the
# haploid chromosomes, the larger chromosome will over-attract the WW vector
# towards it, leading to a biased estimate. Assume that the outermost unitigs in
# each haplotype are the same distance from the origin. Use the initial guess to
# partition unitigs into two groups, find the outermost in each group, and
# rotate to balance those two points

projected_wc_vectors_rbc <-
  map2(prcomps, projected_wc_vectors, function(x,y) {
    projections <- x$x[!grepl('inverted', rownames(x$x)), 1:2] %*% y 
    
    sgns <- sign(projections)
    if(!(1 %in% sgns & -1 %in% sgns)) {
      return(y)
    }
    
    max_pos_unitig <- rownames(projections)[which.max(projections)]
    max_neg_unitig <- rownames(projections)[which.min(projections)]
    # break ties
    
    max_pos_angle <- uv(x$x[max_pos_unitig, 1:2]) %*% uv(y)
    max_neg_angle <- uv(x$x[max_neg_unitig, 1:2]) %*% uv(y)
    rotation_amount <-  as.vector(acos(max_pos_angle) - acos(abs(max_neg_angle)))/2
    
    out <- 
      (twod_rotation_mat(rotation_amount) %*% y) %>% 
      as.vector() %>% 
      set_names(c('PC1', 'PC2'))
    
    return(uv(out))
  })

projected_ww_vectors_rbc <-
  map(projected_wc_vectors_rbc, function(x) {
    out <- c(x[2], -x[1]) %>% 
      set_names(c('PC1', 'PC2'))
    
    return(out)
    
  })

unprojected_wc_vectors_rbc <-
  map2(projected_wc_vectors_rbc, prcomps, unproject_prcomp, value = 'wc_ssf_rbc', n_components = 2)


unprojected_ww_vectors_rbc <-
  map2(projected_ww_vectors_rbc, prcomps, unproject_prcomp, value = 'ww_ssf_rbc', n_components = 2)


### GLMs --------------------------------------------------------------------

model_input <-
  prcomps %>% 
  map(get_prcomp_plotdata, 'unitig_dir')

model_input <-
  map(model_input, function(x) {
    x <-
      x %>%
      # get_prcomp_plotdata(., 'unitig_dir') %>%
      mutate(unitig = gsub('_inverted', '', unitig_dir)) %>%
      left_join(unitig_lengths_df, by = 'unitig')
    
    x <-
      x %>%
      mutate(y = grepl('inverted', unitig_dir)) 
    
    return(x)
  })

glms <-
  map(model_input, function(d) {
    w <- d$length
    d <-
      d %>% 
      select(-unitig, -unitig_dir, -length)
    
    mod <-
      glm(
        y ~ -1 + .,
        data = d,
        # weights = w,
        family = binomial(link = 'logit')
      )
    return(mod)
  })

coefs <-
  map(glms, coef) 

unprojected_ww_vectors_glm <-
  map2(coefs, prcomps, unproject_prcomp, n_components=2, value="ww_ssf_glm")

unprojected_wc_vectors_glm <-
  coefs %>% 
  map(function(x) c(x[2], -1*x[1])) %>% 
  map2(prcomps, unproject_prcomp, n_components=2, value="wc_ssf_glm")


## Bind --------------------------------------------------------------------

library_weights <-
  list(
    unprojected_ww_vectors,
    unprojected_wc_vectors,
    unprojected_ww_vectors_rbc,
    unprojected_wc_vectors_rbc,
    unprojected_ww_vectors_glm,
    unprojected_wc_vectors_glm
  ) %>% 
  map(bind_rows, .id='cluster') %>% 
  purrr::reduce(full_join, by=c('cluster', 'lib'))

### Plots -------------------------------------------------------------------

# 
# coefs <-
#   coefs %>%
#   map_dfr(function(x) as_tibble(as.list(x)), .id='cluster')
# 
# 
# 
# model_input <-
#   bind_rows(model_input, .id = 'cluster') %>%
#   select(cluster, unitig, unitig_dir, length, y, everything())
# 
# projected_wc_vectors_plot <-
#   projected_wc_vectors_rbc %>%
#   map(as.data.frame.list) %>%
#   bind_rows(.id = 'cluster')
# 
# projected_ww_vectors_plot <-
#   projected_ww_vectors %>%
#   map(as_tibble) %>%
#   bind_rows(.id = 'cluster')
# 
# # Test plots
# label_data <-
#   model_input %>%
#   group_by(cluster) %>%
#   filter(length == max(length)) %>%
#   ungroup()
# 
# model_input %>%
#   ggplot() +
#   geom_point(aes(
#     x = PC1,
#     y = PC2,
#     color = y,
#     size = length
#   )) +
#   facet_wrap(~ cluster) +
#   geom_text(aes(x = PC1, y = PC2, label=unitig_dir), size=2, color='black', data=filter(label_data, !grepl('projected', unitig_dir))) +
#   # geom_abline(aes(slope=PC2/PC1, intercept=0), data=projected_ww_vectors) +
#   geom_abline(aes(slope = PC2 / PC1, intercept = 0), data = projected_wc_vectors_plot) +
#   # geom_abline(aes(slope = -PC1 / PC2, intercept = 0), data = projected_wc_vectors_plot) +
#   geom_abline(aes(slope = -PC1 / PC2, intercept = 0), data = projected_ww_vectors_plot, linetype='dashed', color='green') +
#   # geom_abline(aes(slope = PC2 / PC1, intercept = 0), data = projected_ww_vectors_plot, linetype='dashed', color='green') +
#   geom_abline(aes(slope = -PC1 / PC2, intercept = 0), data = coefs, linetype = 'dotted', color='red') +
#   # geom_abline(aes(slope = PC2 / PC1, intercept = 0), data = coefs, linetype = 'dotted', color='red') +
# 
#   # geom_point(aes(x = PC2, y = -PC1), data = projected_ww_vectors) +
#   # facet_wrap( ~cluster + grepl('projected', unitig_dir)) +
#   scale_size_area() +
#   coord_equal()

# 
# map2(norm_ww_vectors, unprojected_ww_vectors, function(x, y) {
#   cosine_similarity_(x$ww_ssf, y$ww_ssf)
# })
## Count Haplotype Markers -------------------------------------------------

# TODO this weight thing isn't quite right.
hmc <-
  imap(em_counts, function(counts, nm){
    # browser()

    # Project weights onto the cube, to use indicator variable interpretation?
    weights <-
      unprojected_wc_vectors_glm[[nm]] %>% 
      mutate(wc_ssf = wc_ssf_glm/max(abs(wc_ssf_glm), na.rm=TRUE))
    
    out <-
      counts %>%
      right_join(weights, by='lib')

    out <-
      out %>%
      mutate( # projected_counts
        p_c = ifelse(sign(wc_ssf) == 1, wc_ssf * c, abs(wc_ssf) * w),
        p_w = ifelse(sign(wc_ssf) == 1, wc_ssf * w, abs(wc_ssf) * c)
      ) %>%
      group_by(unitig) %>%
      summarise(c = sum(p_c), w = sum(p_w)) %>% 
      ungroup() %>% 
      mutate(n=ceiling(c + w))

    return(out)
  })

marker_counts <-
  bind_rows(hmc) %>%
  bind_rows(tibble(unitig = character(), n=integer(), c=double(), w=double())) 


# One Unitig Clusters -----------------------------------------------------

#FIXME Design decision: There are no one unitig clusters!
## NA-Out One Small One-Unitig Clusters -----------------------------------

# Want to filter out one unitig clusters on a components with other clusters
one_unitig_cluster_unitigs <-
  cluster_df %>%
  group_by(cluster) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  pull_distinct(unitig)

# TODO maybe make this a percentage of component threshold, like a one unitig cluster makes
# up > 90% of a component.
multi_cluster_components <-
  components_df %>%
  left_join(cluster_df, by='unitig') %>%
  filter(!is.na(cluster)) %>%
  distinct(component, cluster) %>%
  count(component) %>%
  filter(n > 1) %>%
  pull(component)

# Often better to not attempt to assign solo unitigs, which are often
# misclustered. Often, the proper phasing can be inferred from the unitigs
# surrounding the misclustered unitig. However, sometimes half of a chromosome
# is all by itself in a cluster, and we still want to include those.
length_threshold <- 10e6

unitigs_to_zero_out <-
  cluster_df %>% 
  left_join(unitig_lengths_df, by='unitig') %>% 
  left_join(components_df, by='unitig') %>% 
  filter(unitig %in% one_unitig_cluster_unitigs) %>% 
  filter(component %in% multi_cluster_components) %>% 
  filter(length < length_threshold) %>% 
  pull_distinct(unitig)


marker_counts <-
  marker_counts %>%
  mutate(across(c(c, w), function(x) {
    ifelse(unitig %in% unitigs_to_zero_out, NA, x)
  }))

# Export ------------------------------------------------------------------

cat('Exporting\n')

## Marker Counts -----------------------------------------------------------


### Join in Excluded Unitigs ------------------------------------------------
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


### Additional Information --------------------------------------------------

marker_counts <-
  marker_counts %>%
  left_join(cluster_df, by='unitig') %>%
  left_join(unitig_lengths_df, by='unitig') %>%
  left_join(components_df, by='unitig')

marker_counts <-
  marker_counts %>%
  left_join(strand_orientation_clusters_df, by=c('unitig')) %>%
  dplyr::rename(unitig_orientation = strand_cluster)

marker_counts <-
  marker_counts %>%
  mutate(ssf = round((w-c)/(w+c), 3))
# first three columns: name, counts_1, counts_2 for rukki
marker_counts <-
  marker_counts %>%
  dplyr::rename(hap_1_counts = c, hap_2_counts = w) %>%
  select(unitig, hap_1_counts, hap_2_counts, n, ssf, everything(), exclusion) %>%
  arrange(cluster)

### Bandage Cluster Colors -----------------------------------------------------

cluster_palette <-
  marker_counts %>%
  distinct(cluster) %>%
  mutate(Color = rainbow(n()))

marker_counts <-
  marker_counts %>%
  left_join(cluster_palette, by='cluster') %>%
  mutate(Color = pmap_chr(list(Color, hap_1_counts, hap_2_counts), make_bandage_colors))

# fill NA with 0 for rukki?
# Round marker counts for rukki
marker_counts <-
  marker_counts %>%
  mutate(attempted = !is.na(hap_1_counts) & !is.na(hap_2_counts)) %>%
  mutate(hap_1_counts = ifelse(is.na(hap_1_counts), 0, ceiling(hap_1_counts)),
         hap_2_counts = ifelse(is.na(hap_2_counts), 0, ceiling(hap_2_counts)))

stopifnot(nrow(marker_counts) == length(all_unitigs))
stopifnot(setequal(all_unitigs, marker_counts$unitig))

stopifnot(!anyNA(marker_counts$hap_1_counts))
stopifnot(!anyNA(marker_counts$hap_2_counts))

## Haplotype Size Evening --------------------------------------------------

# call obvious hets? Or just any? how will this affect with HOM?
ratio_threshold <- 2

hap_calls <-
  marker_counts %>%
  filter(!is.na(cluster)) %>%
  mutate(hap_1_counts = hap_1_counts + 1,
         hap_2_counts = hap_2_counts + 1) %>%
  mutate(call = ifelse(hap_1_counts/hap_2_counts >= ratio_threshold, 1, NA)) %>%
  mutate(call = ifelse(hap_2_counts/hap_1_counts >= ratio_threshold, 0, call)) %>%
  filter(!is.na(call)) %>%
  mutate(original_call = call)

hap_size_ratio_score <-
  hap_calls %>%
  count(call, wt = length) %>%
  with(n[1]/n[2]) %>%
  log() %>%
  abs()

n_iter <- 100
clusters <-
  hap_calls %>%
  pull_distinct(cluster) %>%
  set_names()

for(unused in 1:n_iter) {
  swap_scores <-
    map_dbl(clusters, function(x) {
      tmp_hap_calls <-
        hap_calls %>%
        mutate(call = ifelse(cluster == x, abs(call-1), call))
      
      out <-
        tmp_hap_calls %>%
        count(call, wt = length) %>%
        with(n[1]/n[2]) %>%
        log() %>%
        abs()
      
      return(out)
    })
  
  min_score <-
    swap_scores[which.min(swap_scores)]
  
  if(min_score < hap_size_ratio_score) {
    cat('swapping cluster:', names(min_score), '\n')
    hap_size_ratio_score <- min_score
    hap_calls <-
      hap_calls %>%
      mutate(call = ifelse(cluster == names(min_score), abs(call-1), call))
  }
  
}

hap_calls <-
  hap_calls %>%
  mutate(swap = !(call == original_call)) %>%
  distinct(cluster, swap)

bad <-
  hap_calls %>%
  count(cluster) %>%
  filter(n> 1)

if(nrow(bad) > 0) {
  stop('hap swapping bad')
}

marker_counts <-
  marker_counts %>%
  left_join(hap_calls, by=c('cluster')) %>%
  mutate(swap = ifelse(is.na(swap), FALSE, swap)) %>%  #unitigs without clusters
  mutate(h1_temp = hap_1_counts, h2_temp = hap_2_counts) %>%
  mutate(
    hap_1_counts = ifelse(swap, h2_temp, h1_temp),
    hap_2_counts = ifelse(swap, h1_temp, h2_temp)
  ) %>%
  select(-h1_temp, -h2_temp, -swap)

stopifnot(!anyNA(marker_counts$hap_1_counts))
stopifnot(!anyNA(marker_counts$hap_2_counts))

stopifnot(all(as.integer(marker_counts$hap_1_counts) == marker_counts$hap_1_counts))
stopifnot(all(as.integer(marker_counts$hap_2_counts) == marker_counts$hap_2_counts))
### CSV ---------------------------------------------------------------------

marker_counts <-
  marker_counts %>%
  select(unitig, hap_1_counts, hap_2_counts, n, everything())

readr::write_csv(marker_counts, output_counts)

## Library Weights ---------------------------------------------------------
# add missing libraries back in
library_weights <-
  library_weights %>% 
  mutate(lib = factor(lib, levels = lib_names)) %>% 
  tidyr::complete(cluster, lib) #%>% 
# mutate(across(contains('ssf'), function(x) ifelse(is.na(x), 0, x)))

# TODO schema checks


### CSV ---------------------------------------------------------------------

readr::write_csv(library_weights, output_lib)


# Warnings ----------------------------------------------------------------


for(w in WARNINGS) {
  warning(w)
}

