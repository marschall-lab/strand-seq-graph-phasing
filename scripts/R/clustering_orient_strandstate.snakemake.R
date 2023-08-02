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
    '--output-marker-counts',
    '--output-lib',
    
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
# homology <- get_values("--homology", singular=TRUE)
connected_components <- get_values('--connected-components', singular=TRUE)

## Parameters
segment_length_threshold <- as.numeric(get_values('--segment-length-threshold'))
expect_XY_separate <- as.logical(get_values('--expect-XY-separate', singular=TRUE))
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


### POCC-----------------------------------------------------------

cat('Propagating one cluster components\n')
cluster_df <- propagate_one_cluster_components(cluster_df, components_df)

### Cosine Cluster Merging ----------------------------------------------------

cat('Cosine cluster merging\n')

# Sometimes, some of the newly created clusters will should be merged into other
# components on cluster (centromere troubles especially)
cluster_df <- merge_similar_clusters_on_components(counts_df, cluster_df, components_df, similarity_threshold = 0.5)

### NA Cluster --------------------------------------------------------------
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

# PAR and haploid need to be oriented together.
cluster_df <-
  cluster_df %>% 
  mutate(cluster = ifelse(cluster %in% par_clusters | grepl('^sex', cluster), 'LGXY', cluster))

# cluster_df <-
#   cluster_df %>%
#   mutate(cluster = ifelse(grepl('^sex', cluster), 'LGXY', cluster))


### Small Cluster Removal --------------------------------------------

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
  mutate(cluster = ifelse(cluster %in% small_clusters & grepl('LGcos', cluster), NA, cluster)) %>%
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

mem_data_plane <- 
  get_chrom_cluster_data_planes(counts_df, cluster_df, unitig_lengths_df, supervision = 'PC1')
    
# It appears that I need to do clustering only on the first PC, which
# corresponds to orientation, as otherwise, other dimensions can influence the
# results and lead to grouping of a unitig and its inversion together.

strand_orientation_clusters_df <-
  mem_data_plane$model_input %>% 
  select(cluster, unitig, unitig_dir, PC1) %>% 
  mutate(strand_cluster = sign(PC1)) %>% 
  select(-PC1)

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

# Phase Chromosomes -------------------------------------------------------

## Orient Fastmap Counts --------------------------------------------------
# exact_match_counts_df <- raw_exact_match_counts_df

exact_match_counts_df <-
  orient_counts(exact_match_counts_df, strand_orientation_clusters_df)

fastmap_data_plane <-
  get_chrom_cluster_data_planes(exact_match_counts_df, cluster_df, unitig_lengths_df, supervision = 'inverse')

## Counts ------------------------------------------------------------------

em_counts <-
  exact_match_counts_df %>%
  bind_with_inverted_unitigs() %>%
  left_join(cluster_df, by='unitig') %>%
  split(.$cluster) 

## One Haplotype Clusters --------------------------------------------------

### Detect one Haplotype Clusters --------------------------------------------------

# NOTE. This isn't detecting just clusters with one haplotype: It will also
# detect two haplotype clusters where only one unitig is "good enough" to use
# for sorting. However, sorting those one-good-unitig two-haplotype clusters
# happens in the same way as sorting a one haplotype cluster.

is_one_haplotype_cluster <-
  imap(em_counts, function(x, nm) {
    # TODO collect and justify parameters
    cat('Detecting number of clusters in:', nm, '\n')

    sort_counts <-
      x %>%
      filter(!grepl('inverted', unitig_dir))
    
    unitigs <-
      sort_counts %>% 
      pull_distinct(unitig)
    
    if(length(unitigs) == 1) {
      cat('One unitig cluster:', nm, '\n')
      return(unitigs)
    }
    
    
    sort_counts <- 
      sort_counts %>% 
      left_join(fastmap_data_plane$weights, by=c('cluster', 'lib')) %>% 
      mutate(p_weight = ifelse(is.na(p_weight), 0, p_weight))
    
    # TODO this weight thing isn't quite right.
    sort_counts <-
      sort_counts %>% 
      mutate( # pseudo_counts
        p_c = ifelse(sign(p_weight) == 1, c/p_weight, w/abs(p_weight)),
        p_w = ifelse(sign(p_weight) == 1, w/p_weight, c/abs(p_weight))
      ) 
    
    # scale marker ratio to original N
    sort_counts <-
      sort_counts %>% 
      mutate(
        c = ifelse(p_c + p_w == 0, 0, p_c/(p_c + p_w) * n),
        w = ifelse(p_c + p_w == 0, 0, p_w/(p_c + p_w) * n)
      )%>% 
      select(-p_c, -p_w) %>% 
      mutate(c = as.integer(ceiling(c)), w=as.integer(ceiling(w)))



    # Trying to sort on good, big unitigs that are often present in verkko graphs.
    unitigs <-
      sort_counts %>%
      left_join(unitig_lengths_df, by='unitig') %>% 
      filter(length >= 1e6) %>% 
      pull_distinct(unitig) 

      if(length(unitigs) < 1) {
        cat('No adequately long unitigs for number of clusters detection:', nm, '\n')
        return(NA)
      }

      wc <-
        sort_counts %>%
        with(make_wc_matrix(w,c,lib,unitig, min_n=1))

      unitig_vector_lengths <-
        wc %>%
        apply(1, function(x) (sqrt(sum(x^2, na.rm=T))))

      print(unitig_vector_lengths)

      threshold <-
        max(unitig_vector_lengths) - 0.2 *  max(unitig_vector_lengths)
      
      #FIXME A genuine one-haplotype cluster may be flagged by this incorrectly.
      #If something fails this check, it needs to be looked at in unweighted
      #data space as well. Similarly, a one haplotype cluster may look very WC
      #after down-weighting WW libraries. OR maybe not I have run tests with
      #manually created one-haplotype clusters and it seems fine...
      
      # TODO Minimum vector length needs to be scaled by
      #num libraries
      min_vector_length <- 4
      unitigs <-
        names(unitig_vector_lengths)[unitig_vector_lengths > min_vector_length & unitig_vector_lengths > threshold]

      if(length(unitigs) < 1) {
        cat('Too few adequately low-noise unitigs for number of clusters detection:', nm, '\n')
        return(NA)
      }
    
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
      cat('No valid similarity scores for cluster:', nm, '\n')
      return(NA)
    }

    num_clusters <-
      out %>%
      pull_distinct(sort) %>%
      length()

    if(num_clusters > 1) {
      return(character())
    } else if(num_clusters == 1){
      return(out %>% pull_distinct(unitig))
    } else {
      stop("You shouldn't be here.")
    }

  })

one_haplotype_cluster_unitigs <-
  is_one_haplotype_cluster %>% 
  discard(function(x) length(x) == 0) %>% 
  discard(function(x) all(is.na(x)))

two_haplotype_clusters <- 
  cluster_df %>% 
  pull_distinct(cluster) %>% 
  setdiff(names(one_haplotype_cluster_unitigs))

### Filter One-haplotype Clusters -------------------------------------------

# Want to filter out one unitig cluster on a components with other clusters
one_unitig_cluster_unitigs <-
  cluster_df %>% 
  group_by(cluster) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  pull_distinct(unitig)

# TODO maybe make this a percentage of component threshold, like a one unitig cluster makes
# up > 90% of a component.
one_cluster_components <-
  components_df %>%
  left_join(cluster_df, by='unitig') %>%
  filter(!is.na(cluster)) %>%
  distinct(component, cluster) %>%
  count(component) %>%
  filter(n == 1) %>%
  pull(component)

one_cluster_component_unitigs <-
  cluster_df %>% 
  left_join(components_df, by='unitig') %>%
  filter(component %in% one_cluster_components) %>% 
  pull_distinct(unitig)

one_haplotype_cluster_unitigs <-
  one_haplotype_cluster_unitigs %>%
  map(function(x) {
    if (length(x) == 1 && x %in% one_unitig_cluster_unitigs) {
      if (!(x %in% one_cluster_component_unitigs)) {
        return(NA)
      }
    }
    
    return(x)
    
  })

one_haplotype_cluster_unitigs <-
  one_haplotype_cluster_unitigs %>% 
  discard(function(x) all(is.na(x))) 

### Count one-haplotype clusters --------------------------------------------


# arrange by size
one_haplotype_clusters <-
  cluster_df %>% 
  filter(cluster %in% names(one_haplotype_cluster_unitigs)) %>% 
  left_join(unitig_lengths_df, by='unitig') %>% 
  group_by(cluster) %>% 
  summarise(length=sum(length), .groups='drop') %>% 
  arrange(length) %>% 
  pull(cluster)

# Assign to haplotype alternating  by size.
one_haplotype_cluster_counts <-
  imap(one_haplotype_clusters, function(clust, n) {
  
  counts <-
    em_counts[[clust]] %>%
    filter(!grepl('inverted', unitig_dir))
  
  swapping_unitigs <- 
    one_haplotype_cluster_unitigs[[clust]]
  
  # stopifnot(length(swapping_unitig) == 1)

  swaps <- 
    counts %>% 
    filter(unitig %in% swapping_unitigs) %>% 
    group_by(lib) %>% 
    summarise(c=sum(c), w=sum(w), .groups='drop') %>% 
    mutate(swap = c > w) %>% 
    distinct(lib, swap)
  
  counts <-
    counts %>% 
    left_join(swaps, by='lib') %>% 
    mutate(tmp_c = ifelse(swap, c, w),
           tmp_w = ifelse(swap, w, c)
           )
  if ((n %% 2) == 1) {
    counts <-
      counts %>% 
      mutate(c = tmp_c, w = tmp_w)
  } else {
    counts <-
      counts %>% 
      mutate(c = tmp_w, w = tmp_c)
  }
  
  counts %>% 
    group_by(unitig) %>% 
    summarise(n = sum(n), c=sum(c), w=sum(w)) %>% 
    select(unitig, n, c, w)
})

# In case there are 0 uni clusters.
one_haplotype_cluster_counts <-
  one_haplotype_cluster_counts %>% 
  bind_rows(tibble(unitig = character(), n=integer(), c=double(), w=double()))
          

## Two Haplotype Clusters --------------------------------------------------

### Count -----------------------------------------------------------------

# TODO this weight thing isn't quite right.
hmc <-
  imap(em_counts[two_haplotype_clusters], function(counts, nm){
    # browser()
    
    weights <- 
      fastmap_data_plane$weights %>% 
      filter(cluster == nm)
    
    out <- 
      counts %>% 
      filter(!grepl('inverted', unitig_dir)) %>% 
      right_join(weights, by='lib') 

    # TODO this weights thing isn't quite right.
    out <-
      out %>% 
      mutate( # pseudo_counts
        p_c = ifelse(sign(b_weight) == 1, b_weight^2 * c, b_weight^2 * w),
        p_w = ifelse(sign(b_weight) == 1, b_weight^2 * w, b_weight^2 * c)
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
  bind_rows(tibble(unitig = character(), n=integer(), c=double(), w=double())) %>% 
  bind_rows(one_haplotype_cluster_counts)

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
  left_join(unitig_lengths_df, by='unitig')

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


### Haplotype Size Evening --------------------------------------------------

# TODO, flip haplotypes labels that make the overall haplotype sizes more even.


### CSV ---------------------------------------------------------------------

marker_counts <-
  marker_counts %>% 
  select(unitig, hap_1_counts, hap_2_counts, everything())
  
readr::write_csv(marker_counts, output_counts)


## Library Weights ---------------------------------------------------------

mem_weights <-
  mem_data_plane$weights %>% 
  select(-b_weight) %>% 
  dplyr::rename(ww_weight_mem = p_weight)

fastmap_weights <-
  fastmap_data_plane$weights %>% 
  dplyr::rename(
    ww_weight_fastmap = p_weight,
    wc_weight_fastmap = b_weight
    )

library_weights <-
  full_join(mem_weights, fastmap_weights,by=c('cluster', 'lib')) 

# TODO schema checks


### CSV ---------------------------------------------------------------------

readr::write_csv(library_weights, output_lib)


# Warnings ----------------------------------------------------------------


for(w in WARNINGS) {
  warning(w)
}

