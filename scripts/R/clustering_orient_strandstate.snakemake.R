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
    # filterThreshold = 0.7,
    # lowQualThreshold = 0.8,
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

mean_coverage <-
  mean_coverage %>%
  mutate(unitig_range = spoof_range(unitig))

vector_lengths <-
  strand.freq@.Data %>%
  apply(1, function(x) sqrt(sum(x * x, na.rm=TRUE))) %>%
  tibble::enframe(name = 'unitig_range', value = 'vec_length')

weights <-
  components_df %>%
  left_join(unitig_lengths_df) %>%
  group_by(component) %>%
  mutate(component_size = sum(length)) %>%
  ungroup()

weights <-
  weights %>%
  mutate(unitig_range = spoof_range(unitig)) %>%
  left_join(vector_lengths) %>%
  left_join(mean_coverage) %>%
  mutate(vec_length = round(vec_length))

weights <-
  weights %>%
  arrange(component_size, vec_length, coverage, length) %>%
  mutate(weight = 1:n())

weights <-
  weights %>%
  with(set_names(weight, unitig_range))

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
component_fraction_threshold <- 0.15

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
# Use coverage instead of length?
length_quantile_threshold <- 1 # 0.1
score_thresh <- 0.5

# TODO fix this while condition so that a break is not needed
while(any_assigned || (length_quantile_threshold <= 1 || score_thresh >= 0.5)) {


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
      sims <- abs(similarities[unclustered_unitigs, x, drop=FALSE])
      apply(sims, 1, mean, na.rm=TRUE)
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

    cat(
      'Unitig length quantile threshold:',
      length_quantile_threshold,
      '\n'
    )

    unitigs_to_consider <-
      unitig_lengths_df %>%
      slice_max(order_by=length, prop=length_quantile_threshold) %>%
      pull(unitig) %>%
      intersect(unclustered_unitigs)

    unclustered_sims <-
      similarities[unitigs_to_consider, unitigs_to_consider, drop=FALSE] %>%
      abs()

    unclustered_sims <-
      unclustered_sims * upper.tri(unclustered_sims)

    max_pair <-
      which(unclustered_sims == max(unclustered_sims, na.rm=TRUE), arr.ind = TRUE)

    max_sim <-
      unclustered_sims[max_pair[1], max_pair[2]]

    if(!is.na(max_sim) && max_sim > score_thresh) {
      any_assigned <- TRUE
      new_cluster_ix <- new_cluster_ix + 1

      new_cluster_unitigs <- c(rownames(unclustered_sims)[max_pair[1]], colnames(unclustered_sims)[max_pair[2]])
      cat('creating new cluster with unitigs:', new_cluster_unitigs,
          'with similarity:', max_sim,
          '\n')

      cluster_df <-
        cluster_df %>%
        mutate(cluster = ifelse(unitig %in% new_cluster_unitigs, paste0('LGcosLong', new_cluster_ix), cluster))
    } else {
      if(length_quantile_threshold == 1L) {
        break
      }
      length_quantile_threshold <- length_quantile_threshold + 0.1
      if(length_quantile_threshold > 1) length_quantile_threshold <- 1L
    }
  }


}


## Max-Based Assignment ----------------------------------------------------

# TODO Assign a unitig to a cluster if it has a very high similarity with any of the individual members?

## POCC-----------------------------------------------------------

cat('Propagating one cluster components\n')
cluster_df <- propagate_one_cluster_components(cluster_df, components_df)

## Cosine Cluster Merging ----------------------------------------------------

cat('Cosine cluster merging\n')
# Sometimes, some of the newly created clusters will should be merged
# into other components on cluster (centromere troubles especially)
cluster_df <- merge_similar_clusters_on_components(counts_df, cluster_df, components_df, similarity_threshold = 0.5)
cluster_df <- merge_similar_clusters(counts_df, cluster_df, similarity_threshold = 0.66)


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

# cluster_df <-
#   cluster_df %>%
#   mutate(cluster = ifelse(grepl('^sex', cluster), 'LGXY', cluster))


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
  mutate(cluster = ifelse(cluster %in% small_clusters & grepl('LGcos', cluster), NA, cluster))

# Assign NA cluster -------------------------------------------------------


cluster_df <-
  cluster_df %>%
  mutate(cluster = ifelse(is.na(cluster), paste0('LGNA2_', unitig), cluster))

# Orientation Detection w/ Inverted Unitigs -------------------------------

cat('Detecting unitig orientation\n')
# Add inverted version of every unitig to dataset. Guarantees that there will
# unitigs in both orientations when clustering

# TODO need to double check /rerun the strand orientation clustering. If The
# haploid and PAR clusters are oriented separately ~ chance that they could
# become unsynchronized. Is there a way to fix this? Do it at the haploid clustering step?

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

# TODO (for fun) an algorithm to derive the proper WC library weights vector from mem
# data. A cheap solution could be just to filter the data by mapq to "reveal"
# the WC signal. However, I am curious to see if there is a way to derive the WC
# library weights without that.

# It appears that I need to do clustering only on the first PC, which
# corresponds to orientation, as otherwise, other dimensions can influence the
# results and lead to grouping of a unitig and its inversion together. *NOTE*
# This isn't true for the case of the X and Y chromosomes. While the separation
# between autosomal chromosomes is not visible in unfiltered mem data the X and
# Y chromosomes separate (the same way that autosomal chromosomes do with
# fastmap data.) This can lead to issues with just taking the sign of the first
# PC, and has on one occasion (HGSVC 1.4.1 HG00512) lead to a misorientaiton of
# a large unitig. To be more robust to this possibility, now we take the sign of
# the projection onto the line perpendicular to the boundary plane. This works
# becuase the boundary plane rotates with the data, and should be able to better
# account for the increased spread with X and Y chromosomes. It's not foolproof
# but it should be better, in part because it accounts for the length of each
# unitig.

# TODO the supervision with PC1 method is a little brittle for X and Y. Better
# would be a way that works directly with repulsive pairs to clusters the two
# orientations. Maybe something perpendicular to the mean vector between each
# pair of inverses or something like that, though this would require some
# thought with regards to the signs of the vectors.

# FIXME The above has become a fixme, because there are now occasions where the
# PAR is mis-clustered due to improper orientation. Mean Perependicular vector?
# Thats also an imperfect solution ~ some sort of MEC style formulation?


cluster_counts <-
  counts_df %>%
  bind_with_inverted_unitigs() %>% # self-supervision
  left_join(cluster_df, by='unitig') %>%
  split(.$cluster)

# Handle one unitig clusters separately.This is a little bit of a hack to handle
# some clusters that only have one unitigs and are particularly poorly behaved.
# It probably is resulting from a lack of consistency across thresholds. EG, a
# unitig that fails all thresholds at this step somehow passed a threshold
# earlier.
n_unitigs_per_cluster <- 
  cluster_counts %>% 
  map(function(x) n_distinct(x$unitig_dir))

one_unitig_clusters <-
  n_unitigs_per_cluster %>% 
  keep(function(x) x==2) %>% 
  names()
  
wfracs <-
  cluster_counts[!names(cluster_counts) %in% one_unitig_clusters] %>%
  map(function(x) with(x, make_wc_matrix(w,c,lib,unitig_dir,min_n=min_n)))  %>%
  # filling with 0s doesn't seem to affect first PC too much, compared to
  # probabilistic or Bayesian PCA (from pcaMethods bioconductor package)
  map(function(x) {
    # TODO fill 0 or not? If not fill 0, how to handle possible NA similarities?
    # How to handle unitigs w/ a PC=0?
    x[is.na(x)] <- 0
    return(x)
  })

sims <-
  map(wfracs, cosine_similarity)

#if whole row is NA ~ then it had 0 for all wfracs.
sims <-
  map(sims, function(x) {
    row_ix <-
      apply(x, 1, function(xx) all(is.na(xx)))

    col_ix <-
      apply(x, 2, function(xx) all(is.na(xx)))


    x[!row_ix, !col_ix, drop=FALSE]
  })

dists <-
  map(sims, function(x) as.dist(1-x))


clusts <-
  map(dists, function(x) {
    hclust(x, method = 'complete') %>%
      cutree(k=2)
  })

strand_orientation_clusters_df <-
  clusts %>%
  map_dfr(tibble::enframe, 'unitig_dir', 'strand_cluster', .id='cluster') %>%
  mutate(unitig = gsub('_inverted$', '', unitig_dir)) %>%
  mutate(strand_cluster = ifelse(strand_cluster==1, -1, 1))

# add back in any unitigs removed at ealier step
strand_orientation_clusters_df <-
  strand_orientation_clusters_df %>%
  right_join(
    cluster_counts %>%
      bind_rows() %>%
      distinct(unitig, unitig_dir, cluster)
  )

strand_orientation_clusters_df <-
  strand_orientation_clusters_df %>%
  mutate(strand_cluster = ifelse(is.na(strand_cluster), 0, strand_cluster))

# TODO, a comparison to PC1 clustering, with a possible message if they differ
strand_orientation_clusters_df <-
  strand_orientation_clusters_df %>%
  select(cluster, unitig, unitig_dir, strand_cluster)

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
  filter(n() != 1 | 0 %in% strand_cluster)

if(nrow(bad) > 0) {
  stop('Unitigs have been clustered with their inversions ~ something is wrong with unitig orientation detection')
  # WARNINGS <- c(WARNINGS, 'Unitigs have been clustered with their inversions ~ something is wrong with unitig orientation detection')
}

# Phase Chromosomes -------------------------------------------------------


## Orient Counts -----------------------------------------------------------

# exact_match_counts_df <- raw_exact_match_counts_df

exact_match_counts_df <-
  orient_counts(exact_match_counts_df, strand_orientation_clusters_df)

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

## Calculate WC and WW Basis Vectors ----------------------------------------

# does this work with, eg, a 1 unitig cluster? Or does the dimension projection
# lost a lot of library information?
em_counts <-
  exact_match_counts_df %>%
  left_join(cluster_df, by='unitig') %>% 
  arrange(cluster, unitig, lib) %>%
  split(.$cluster)

ww_vectors <-
  oriented_counts_df %>%
  mutate(ww_ssf = ifelse(n < 1, NA, (w-c)/n))  %>% 
  select(cluster, lib, ww_ssf) %>% 
  split(.$cluster)

ww_vectors <- ww_vectors[names(em_counts)]

em_ssf_mats <- 
  map(em_counts, function(x) {
    with(x, make_wc_matrix(w, c, lib, unitig, min_n=1)) 
  })

em_ssf_mats_proj <-
  map2(em_ssf_mats, ww_vectors, function(x, ww_df) {
    ww_df <- 
      ww_df %>% 
      arrange(lib) 
    
    libs <- 
      ww_df %>% 
      pull(lib)
    
    stopifnot(all(colnames(x) == libs))
    
    ssfs <-
      ww_df %>% 
      pull(ww_ssf)
    
    # apply to each row, because of need for pairwise complete case management.
    # Additionally transpose at end, because, even though the function is
    # applied by row, it fills the results by column...
    x_proj <-
      apply(x, 1, project_through, y=ssfs) %>%
      t()
    
    rownames(x_proj) <- paste0(rownames(x_proj), '_projected')
    
    return(x_proj)
  })

em_ssf_mats <-
  map2(em_ssf_mats,em_ssf_mats_proj, rbind) 

# Concatenante with inverted unitigs
em_ssf_mats_inverted <-
  map(em_ssf_mats, `*`, -1) %>% 
  map(function(x) {
    rownames(x) <- paste0(rownames(x), '_inverted')
    return(x)
  })

em_ssf_mats <- map2(em_ssf_mats, em_ssf_mats_inverted, rbind)


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

# Project WW vectors into prcomp space, to make orthognalization, which is
# equivalent to calculation of WC vector, easy.
projected_ww_vectors <- 
  map2(prcomps, ww_vectors, function(x, y){
    ssf_mat <- matrix(y$ww_ssf, nrow=1)
    colnames(ssf_mat) <- y$lib
    ssf_mat[is.na(ssf_mat)] <- 0 
    
    predict(x, ssf_mat) 
  })

# Rotate and Project back to libray space to get WC Vector, using only first twp
# PCs. Should this also be done for the WW vector? as like a form of noise
# correction? Seems unnessecary, for basically all "normal" cluster, the cosine similarity is > 99%
# unprojected_ww_vectors <-
#   map2(prcomps, projected_ww_vectors, function(x, y) {
#     ww_vec <- y[, 1:2]
#     ww_uv <- ww_vec/sum(ww_vec^2, na.rm=TRUE)
#     unprojected_ww <- x$rotation[, 1:2] %*% ww_uv
#     colnames(unprojected_ww) <- 'ww_ssf'
#     out <- 
#       as_tibble(unprojected_ww, rownames='lib') %>% 
#       mutate(ww_ssf = ww_ssf/(sqrt(sum(ww_ssf^2))))
#     return(out)
#   })


# To go from PC space back to library space for the WC vector, project the first
# two PCS onto the rotated WW vector.
unprojected_wc_vectors <-
  map2(prcomps, projected_ww_vectors, function(x, y) {
    wc_vec <- c(y[, 2], -y[, 1])
    wc_uv <- wc_vec / sum(wc_vec ^ 2, na.rm = TRUE)
    unprojected_wc <- x$rotation[, 1:2] %*% wc_uv
    colnames(unprojected_wc) <- 'wc_ssf'
    out <-
      as_tibble(unprojected_wc, rownames = 'lib') %>%
      mutate(wc_ssf = wc_ssf / (sqrt(sum(wc_ssf ^ 2))))
    return(out)
  })

norm_ww_vectors <-
  ww_vectors %>%
  map(function(x) {
    x %>%
      group_by(cluster) %>%
      mutate(ww_ssf = ww_ssf / sqrt(sum(ww_ssf ^ 2, na.rm = TRUE))) %>%
      ungroup()
  })

## Count Haplotype Markers -------------------------------------------------

# TODO this weight thing isn't quite right.
hmc <-
  imap(em_counts, function(counts, nm){
    # browser()

    weights <-
      unprojected_wc_vectors[[nm]] 

    out <-
      counts %>%
      right_join(weights, by='lib')

    # TODO this weights thing isn't quite right.
    out <-
      out %>%
      mutate( # pseudo_counts
        p_c = ifelse(sign(wc_ssf) == 1, wc_ssf^2 * c, wc_ssf^2 * w),
        p_w = ifelse(sign(wc_ssf) == 1, wc_ssf^2 * w, wc_ssf^2 * c)
      ) %>%
      group_by(unitig) %>%
      summarise(p_c = sum(p_c), p_w = sum(p_w), n=sum(n)) %>%
      ungroup()

    # scale marker ratio to original N
    out <-
      out %>%
      mutate(
        c = ifelse(p_c + p_w == 0, 0, p_c/(p_c + p_w) * n),
        w = ifelse(p_c + p_w == 0, 0, p_w/(p_c + p_w) * n)
      )%>%
      select(-p_c, -p_w)
    # Previously, this was select(-p_c, -p_w, n), which appeared to act exactly
    # the same as select(-p_c, -p_w), despite the apparent bug of the n at
    # the end of the select!

    return(out)
  })

marker_counts <-
  bind_rows(hmc) %>%
  bind_rows(tibble(unitig = character(), n=integer(), c=double(), w=double())) 

## NA-Out One Unitig Clusters -----------------------------------------------

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
  left_join(filter(strand_orientation_clusters_df, strand_cluster==1), by=c('unitig', 'cluster')) %>%
  select(-strand_cluster) %>%
  mutate(unitig_dir = ifelse(grepl('inverted$', unitig_dir), -1, 1)) %>%
  dplyr::rename(unitig_orientation = unitig_dir)

# first three columns: name, counts_1, counts_2 for rukki
marker_counts <-
  marker_counts %>%
  dplyr::rename(hap_1_counts = c, hap_2_counts = w) %>%
  select(unitig, hap_1_counts, hap_2_counts, everything(), exclusion) %>%
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
  select(unitig, hap_1_counts, hap_2_counts, everything())
  
readr::write_csv(marker_counts, output_counts)


## Library Weights ---------------------------------------------------------
# TODO add missing libraries back in?
library_weights <-
  full_join(
    bind_rows(norm_ww_vectors),
    bind_rows(unprojected_wc_vectors, .id='cluster'),
    by=c('cluster', 'lib')
  )

# TODO schema checks


### CSV ---------------------------------------------------------------------

readr::write_csv(library_weights, output_lib)


# Warnings ----------------------------------------------------------------


for(w in WARNINGS) {
  warning(w)
}

