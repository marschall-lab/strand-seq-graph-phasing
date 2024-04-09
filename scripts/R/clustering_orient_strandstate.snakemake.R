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

unproject_prcomp <- function(x, pca, name = 'lib', value = 'wc_ssf') {
  # From: http://factominer.free.fr/question/FAQ.html
  loadings <- sweep(pca$var$coord, 2,sqrt(pca$eig[1:ncol(pca$var$coord),1]),FUN="/")
  unprojected_x <- loadings[, 1:2] %*% uv(x)
  
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

# Coverage weighted mean
cwm_ <- function(unitig_coverage) {
  unitig_coverage <- unitig_coverage
  cwm<- function(mat, weight_f = sum, mat_f = identity, na.rm=FALSE) {
    stopifnot(!is.null(rownames(mat)))
    stopifnot(!is.null(colnames(mat)))
    
    weights <- matrix(nrow = nrow(mat), ncol = ncol(mat))
    dimnames(weights) <- dimnames(mat)
    for(unitig_1 in rownames(weights)) {
      for(unitig_2 in colnames(weights)) {
        weight_1 <- unitig_coverage[unitig_1]
        weight_2 <- unitig_coverage[unitig_2]
        weights[unitig_1, unitig_2] <- weight_f(weight_1, weight_2)
      }
    }
    
    mat <- mat_f(mat)
    out <- weighted.mean(mat, weights, na.rm=na.rm)
    return(out)
  }
  
  return(cwm)
}

# make alpha-numeric factor
manf <- function(x) {
  factor(x, levels = stringr::str_sort(unique(x), numeric=TRUE))
}

# Command Line ------------------------------------------------------------


## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

# source('/Users/henglinm/Documents/GitHub/strand-seq-graph-phasing/scripts/R/module_utils/phasing_test_args.R')

## Parsing Helper ----------------------------------------------------------

# Have to handle multiple inputs per tag
arg_idx <- sort(which(grepl('^--', args)))
arg_idx <- c(arg_idx, length(args) + 1) # edge case of last tag

get_values <- function(arg, singular=TRUE, null_ok=FALSE){
  idx <- which(args == arg)
  if(length(idx) == 0) {
    if(null_ok) {
      return(character())
    } else {
      stop('arg: ', arg, ' not found')
    } 
  } 
  stopifnot(length(idx) == 1)

  next_idx <- arg_idx[which.max(arg_idx > idx)]
  if(next_idx == idx + 1) {
    if(null_ok) {
      return(character())
    } else {
      stop('arg: ', arg, ' not found')
    } 
  }
  values <- args[(idx + 1):(next_idx - 1)]

  # More than one value? return list. One value ~ remove from list structure. It
  # is probably bad practice to return two different types like that, but most
  # arguments have single values and it makes the code less annoying to look at.
  if(singular) {
    stopifnot(length(values)==1)
    values <- values[[1]]
  } else {
    stopifnot(length(values)>=1)
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

# library(contiBAIT)

source(file.path(get_script_dir(), "module_utils/utils.R"))
source(file.path(get_script_dir(), "module_utils/phasing_utils.R"))
# Parsing -----------------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--mem-counts',
    '--fastmap-counts',
    '--connected-components',
    ## Output
    '--intermediate-output-dir',
    '--output-marker-counts',
    '--output-lib',
    
    ## Params
    '--segment-length-threshold',
    '--cluster-PAR-with-haploid',
    '--threads'
  )

stopifnot(all(expected_args %in% args))
## Input

mem_counts <- get_values("--mem-counts", singular=FALSE)
fastmap_counts <- get_values("--fastmap-counts", singular=FALSE)
connected_components <- get_values('--connected-components', singular=TRUE)

## Parameters
segment_length_threshold <- as.numeric(get_values('--segment-length-threshold'), singular=TRUE)
cluster_PAR_with_haploid <- as.logical(get_values('--cluster-PAR-with-haploid', singular=TRUE))
n_threads <- as.numeric(get_values('--threads'))
stopifnot(n_threads >= 1)

## Output
intermediate_output_dir <- get_values('--intermediate-output-dir')

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

# Coverage Weighted Mean --------------------------------------------------

unitig_coverage_df <-
  counts_df %>% 
  group_by(unitig) %>% 
  summarise(coverage = mean(n, na.rm=TRUE))

unitig_coverage_df <-
  bind_rows(
    unitig_coverage_df,
    mutate(unitig_coverage_df, unitig = paste0(unitig, '_inverted'))
  )

unitig_coverage <-
  with(unitig_coverage_df, set_names(coverage, unitig))

coverage_weighted_mean <- cwm_(unitig_coverage)

# SSF Matrix --------------------------------------------------------------

#TODO make this parameter more visible. Explain why only 10
min_n <- 10

ssf_mat <-
  counts_df %>%
  with(make_wc_matrix(w, c, lib, unitig, min_n=min_n))

# Library QC --------------------------------------------------------------


## ContiBAIT QC ------------------------------------------------------------

strand.freq <-
  ssf_mat

cat('Running contiBAIT QC\n')

strand.states <-
  preprocessStrandTable(
    strand.freq,
    lowQualThreshold = 0.7,
    minLib = 10
  )

lib_names <-
  counts_df %>%
  pull_distinct(lib)

temp <- get_values('--included-libraries', null_ok=TRUE)
if(length(temp) != 0) {
  cat('Importing included libraries.\n')
  included_libraries <-
    readr::read_tsv(temp, col_names = FALSE)[[1]]
} else {
  included_libraries <- colnames(strand.states$strandMatrix)
  readr::write_tsv(tibble(included_libraries), file.path(intermediate_output_dir, 'included_libraries.tsv') ,col_names = FALSE)
}

excluded_libraries <- setdiff(lib_names, included_libraries)
WARNINGS <- character()
WARNINGS <- c(WARNINGS, paste0(100*length(included_libraries)/length(lib_names), '% of Strand-seq libraries pass QC'))

counts_df <-
  counts_df %>%
  filter(lib %in% included_libraries)

exact_match_counts_df <-
  exact_match_counts_df %>%
  filter(lib %in% included_libraries)

ssf_mat <- ssf_mat[, colnames(ssf_mat) %in% included_libraries, drop=FALSE]


# Initial Clustering ------------------------------------------------------
temp <- get_values('--initial-clusters', null_ok=TRUE)
if(length(temp) != 0) {
  # TODO unitig should be first cluster in everything
  cat('Importing initial clustering.\n')
  cluster_df <- readr::read_tsv(temp, col_names = c('unitig', 'cluster'), col_types = 'cc')
} else {
  #TODO doublecheck orering with initial clusters
  cluster_df <- mutate(long_unitigs_df, cluster=NA)
  readr::write_tsv(cluster_df, file.path(intermediate_output_dir, 'initial_clusters.tsv'), col_names = FALSE)
}

# Cluster Refinement ------------------------------------------------------

# Cosine Similarity Matrix ------------------------------------------------

# TODO only calculate the cosine similarity mat if it is used (eg, if final
# clusters and strand orientation is provided)

# Why 5? 
min_overlaps <- 5
cosine_similarity_mat <-
  pairwise_complete_cosine_similarity(ssf_mat, min_overlaps = min_overlaps)

# TODO addd back micro cluster on component removal?

temp <- get_values('--final-clusters', null_ok=TRUE)
if(length(temp) != 0) {
  cat('Importing final clustering.\n')
  cluster_df <-
    readr::read_tsv(temp, col_names = c('cluster', 'unitig'), col_types = 'cc')
} else {
  
  for(prop in seq(0.1, 1, 0.1)) {
    cluster_df2 <-
      cluster_df %>% 
      left_join(unitig_coverage_df) %>% 
      left_join(tibble::enframe(strand.states$contigQA, 'unitig', 'qual')) %>% 
      slice_max(tibble(qual, coverage), prop = prop) %>% 
      select(cluster, unitig)
    
    cluster_df2 <-
      cluster_unitigs(
        cosine_similarity_mat,
        cluster_df2,
        new_cluster_id = paste0('[C', prop, ']'),
        cluster_unitig_similarity_threshold = 0.5,
        unitig_unitig_similarity_threshold = 0.5,
        agg_f = coverage_weighted_mean, 
        mat_f=abs,
        na.rm=TRUE
      )

    
    cluster_df <-
      anti_join(cluster_df, cluster_df2, by='unitig') %>% 
      bind_rows(cluster_df2)
    
    
    cluster_df <- merge_similar_clusters_on_components(cosine_similarity_mat, cluster_df, components_df, similarity_threshold = 0.5, agg_f = coverage_weighted_mean, mat_f=abs, na.rm=TRUE)
    cluster_df <- merge_similar_clusters(cosine_similarity_mat, cluster_df, similarity_threshold = 0.5, agg_f = coverage_weighted_mean, mat_f=abs, na.rm=TRUE)
    
  }

  
  ## Cosine Unassigned -----------------------------------------------------
  
  cat('Assigning "high-noise" unitigs using cosine similarity \n')
  
  
  cluster_df <- remove_small_clusters(cluster_df, unitig_lengths_df, threshold = 1e7)
  cluster_df <- propagate_one_cluster_components(cluster_df, components_df)
  
  # cat('Cosine cluster merging\n')
  cluster_df <- merge_similar_clusters_on_components(cosine_similarity_mat, cluster_df, components_df, similarity_threshold = 0.5, agg_f = coverage_weighted_mean, mat_f=abs, na.rm=TRUE)
  cluster_df <- merge_similar_clusters(cosine_similarity_mat, cluster_df, similarity_threshold = 0.5, agg_f = coverage_weighted_mean, mat_f=abs, na.rm=TRUE)

  cluster_df <-
    cluster_unitigs(
      cosine_similarity_mat,
      cluster_df,
      new_cluster_id = 'C',
      agg_f = coverage_weighted_mean,
      mat_f=abs,
      na.rm=TRUE
    )

  
  ## POCC-----------------------------------------------------------
  
  cat('Propagating one cluster components\n')
  cluster_df <- propagate_one_cluster_components(cluster_df, components_df)
  
  ## Cosine Cluster Merging ----------------------------------------------------
  
  cat('Cosine cluster merging\n')
  # Sometimes, some of the newly created clusters will should be merged
  # into other components on cluster (centromere troubles especially)

  cluster_df <- merge_similar_clusters_on_components(cosine_similarity_mat, cluster_df, components_df, similarity_threshold = 0.5, agg_f = coverage_weighted_mean, mat_f=abs, na.rm=TRUE)
  cluster_df <- merge_similar_clusters(cosine_similarity_mat, cluster_df, similarity_threshold = 0.5, agg_f = coverage_weighted_mean, mat_f=abs, na.rm=TRUE)

  readr::write_tsv(
    cluster_df,
    file.path(intermediate_output_dir, 'intermediate_clusters.tsv'),
    col_names = FALSE
  )
}



# Haploid Chromosomes -----------------------------------------------------

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
  map(function(x){
    row_weights <- with(unitig_lengths_df, set_names(length, unitig))
    col_weights <- with(strand.states$qualList, set_names(quality, library))
    row_weights <- row_weights[gsub('_inverted$', '', rownames(x))]
    col_weights <- col_weights[colnames(x)]
    out <- FactoMineR::PCA(x, ncp = 2, row.w = row_weights, col.w = col_weights,  graph = FALSE)
    return(out)
  })

component_percent_variances <- 
  prcomps %>% 
  map(function(x) x$eig[1:2, 2])

hap_threshold <- 0.20

hap_clusters <-
  component_percent_variances %>% 
  keep(function(x) x['comp 2'] <= hap_threshold * 100) %>% 
  names()

if(!(length(hap_clusters) %in% c(0, 2))){
  WARNINGS <- c(WARNINGS, paste('An unexpected number of hap clusters has been detected', length(hap_clusters)))
}

hap_name <- paste(hap_clusters, collapse='_')
cluster_df <-
  cluster_df %>%
  mutate(cluster = ifelse(cluster %in% hap_clusters, hap_name, cluster))

## Small Cluster Removal --------------------------------------------

cluster_df <- remove_small_clusters(cluster_df, unitig_lengths_df, threshold = 1e7)

# Orientation Detection ---------------------------------------------------

temp <- get_values('--unitig-orientation', null_ok=TRUE)
if(length(temp) != 0) {
  cat('Importing unitig orientation.\n')
  strand_orientation_clusters_df <-
    readr::read_tsv(temp, col_names = c('unitig', 'strand_cluster'), col_types = 'cn')
} else {
  
  cat('Detecting unitig orientation\n')
  # Add inverted version of every unitig to dataset. Guarantees that there will
  # unitigs in both orientations when clustering
  
  cluster_cosine_similarities <-
    with(cluster_df, split(unitig, cluster)) %>% 
    map(function(x) cosine_similarity_mat[x, x, drop=FALSE])
  
  
  ## Concatenate with inverted unitigs ---------------------------------------
  
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
  
  ## Hclust ------------------------------------------------------------------
  
  
  strand_orientation_clusters_df<-
    imap(cluster_cosine_similarities, function(x, nm){
      cat('Detecting Orientations in Cluster: ', nm, '\n')
      return(pairwise_complete_hclust_n(x, n=2, agg_f=coverage_weighted_mean, na.rm=TRUE))
    }) 

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
  
  readr::write_tsv(
    strand_orientation_clusters_df,
    file.path(intermediate_output_dir, 'unitig_orientation.tsv'),
    col_names = FALSE
  )
}

# Phase Chromosomes -------------------------------------------------------
cat('Computing marker counts\n')

## Orient Counts -----------------------------------------------------------

exact_match_counts_df <-
  orient_counts(exact_match_counts_df, strand_orientation_clusters_df)
#
# ## Principal Components ----------------------------------------------------
# TODO modify this so that all PCAs are note recomputed and only new clusters
# are recomputed.

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
  map(function(x){
    row_weights <- with(unitig_lengths_df, set_names(length, unitig))
    col_weights <- with(strand.states$qualList, set_names(quality, library))
    row_weights <- row_weights[gsub('_inverted$', '', rownames(x))]
    col_weights <- col_weights[colnames(x)]
    out <- FactoMineR::PCA(x, ncp = length(row_weights), row.w = row_weights, col.w = col_weights,  graph = FALSE)
    return(out)
  })

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
  summarise(across(c(c, w, n), function(x) sum(x, na.rm=TRUE)), .groups = 'drop')

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
    ww_mat <- matrix(y$ww_ssf, nrow=1)
    colnames(ww_mat) <- y$lib
    ww_mat[is.na(ww_mat)] <- 0 
    # there seems to be a bug if there is only one row?
    ww_mat <- rbind(ww_mat, ww_mat)
    predict(x, ww_mat)$coord[1, 1:2] %>% 
      set_names('PC1', 'PC2')
  })


projected_wc_vectors <- 
  map(projected_ww_vectors, function(x){
    wc_vec <- uv(c(x[2], -x[1]))
    wc_vec <- set_names(wc_vec, c('PC1', 'PC2'))
    return(wc_vec)
  })

# To go from PC space back to library space for the WC vector, project the first
# two PCS onto the rotated WW vector.
unprojected_wc_vectors <-
  map2(projected_wc_vectors, prcomps, unproject_prcomp, value = 'wc_ssf')

# Should the WW vector be recalculated by unprojecting the estimate based on the
# first two principal components? As a form of "noise correction?"
unprojected_ww_vectors <-
  map2(projected_ww_vectors, prcomps, unproject_prcomp, value='ww_ssf')

#### Range Balanced Estimate -------------------------------------------------

# Range-Balancing Correction for haplotype clusters of different size: With the
# haploid chromosomes, the larger chromosome will over-attract the WW vector
# towards it, leading to a biased estimate. Assume that the outermost unitigs in
# each haplotype are the same distance from the origin. Use the initial guess to
# partition unitigs into two groups, find the outermost in each group, and
# rotate to balance those two points

projected_wc_vectors_rbc <-
  pmap(list(prcomps, projected_wc_vectors, projected_ww_vectors), function(x, y, z) {
    wc_projections <- x$ind$coord[!grepl('inverted', rownames(x$ind$coord)), 1:2] %*% y 
    ww_projections <- x$ind$coord[!grepl('inverted', rownames(x$ind$coord)), 1:2] %*% z 
    projections <- wc_projections * ww_projections
    wc_signs <- sign(wc_projections)
    # ww_signs <- sign(ww_projections)
    if(!(1 %in% wc_signs & -1 %in% wc_signs)) {
      return(y)
    }
    
    # Sometimes, noise can lead to a unitig to appear almost improperly
    # oriented.
    # allowed_unitigs <- rownames(ww_signs)[ww_signs == 1]
    # wc_projections <- wc_projections[allowed_unitigs, ,drop=FALSE]
    
    max_pos_unitig <- rownames(projections)[which.max(projections)]
    max_neg_unitig <- rownames(projections)[which.min(projections)]
    # break ties
    
    max_pos_angle <- uv(x$ind$coord[max_pos_unitig, 1:2]) %*% uv(y)
    max_neg_angle <- uv(x$ind$coord[max_neg_unitig, 1:2]) %*% uv(y)
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
  map2(projected_wc_vectors_rbc, prcomps, unproject_prcomp, value = 'wc_ssf_rbc')


unprojected_ww_vectors_rbc <-
  map2(projected_ww_vectors_rbc, prcomps, unproject_prcomp, value = 'ww_ssf_rbc')


### GLMs --------------------------------------------------------------------
# 
model_input <-
  prcomps %>%
  map(function(x) x$ind$coord) %>%
  map(as_tibble, rownames='unitig_dir')

model_input <-
  map(model_input, function(x) {
    x <-
      x %>%
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
  map(glms, coef) %>% 
  map(`[`, 1:2)

unprojected_ww_vectors_glm <-
  map2(coefs, prcomps, unproject_prcomp, value="ww_ssf_glm")

unprojected_wc_vectors_glm <-
  coefs %>%
  map(function(x) c(x[2], -1*x[1])) %>%
  map2(prcomps, unproject_prcomp, value="wc_ssf_glm")


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


## Count Haplotype Markers -------------------------------------------------

count_haplotype_markers <- function(counts, lib_weights) {

  lib_weights <-
    lib_weights %>% 
    set_names(function(x)  gsub('wc_ssf.*', 'wc_ssf', x)) %>% 
    mutate(wc_ssf = wc_ssf/max(abs(wc_ssf), na.rm=TRUE))
  
  out <-
    counts %>%
    right_join(lib_weights, by='lib')
  
  # Project weights onto the cube, to use indicator variable interpretation?
  out <-
    out %>%
    mutate(
      tmp_c = ifelse(sign(wc_ssf) == 1, wc_ssf * c, abs(wc_ssf) * w),
      tmp_w = ifelse(sign(wc_ssf) == 1, wc_ssf * w, abs(wc_ssf) * c)
    ) %>%
    group_by(unitig) %>%
    summarise(c = sum(tmp_c, na.rm=TRUE), w = sum(tmp_w, na.rm=TRUE)) %>% 
    ungroup() %>% 
    mutate(n = ceiling(c + w)) %>% 
    mutate(c = round(c), w=round(w))
  
  return(out)
}

marker_counts <- 
  map2(em_counts, unprojected_wc_vectors, count_haplotype_markers)
marker_counts_rbc <- 
  map2(em_counts, unprojected_wc_vectors_rbc, count_haplotype_markers)
marker_counts_glm <-
  map2(em_counts, unprojected_wc_vectors_glm, count_haplotype_markers)


marker_counts <-
  list(cr = marker_counts, 
       rbc=marker_counts_rbc, 
       glm=marker_counts_glm
       ) %>% 
  map(bind_rows) %>% 
  bind_rows(.id='counting_method')

marker_counts <-
  marker_counts %>%
  mutate(ssf = round((w-c)/(w+c), 3)) %>% 
  mutate(wfrac = round(w/(w+c), 2))

marker_counts <-
  marker_counts %>% 
  tidyr::pivot_wider(names_from = counting_method, values_from = c(c, w, n, ssf, wfrac)) %>% 
  relocate(unitig, contains('cr'), contains('rbc'), contains('glm')) 

## Pick Counting Method ----------------------------------------------------

temp <- get_values('--counting-methods', null_ok=TRUE)
if(length(temp) != 0) {
  cat('Importing marker counting methods.\n')
  counting_methods_df <-
    readr::read_tsv(temp, col_names = c('cluster', 'method'), col_types = 'cc')
} else {
  
  counting_methods_df <-
    distinct(cluster_df, cluster) %>% 
    mutate(method = ifelse(grepl('LGXY', cluster) | grepl('_', cluster), 'rbc', 'rbc')) #TODO RBC for all vs cr?
  
  readr::write_tsv(
    counting_methods_df,
    file.path(intermediate_output_dir, 'counting_methods.tsv'),
    col_names = FALSE
  )
  
}


## Assign Haplotype Markers ------------------------------------------------

# first three columns: name, counts_1, counts_2 for rukki
marker_counts <-
  marker_counts %>%
  left_join(cluster_df, by='unitig') %>% 
  left_join(counting_methods_df, by='cluster') %>% 
  mutate(
    hap_1_counts = case_when(
      method == 'cr' ~  c_cr,
      method == 'rbc' ~ c_rbc,
      method == 'glm' ~ c_glm,
      TRUE~NA_real_
    ),
    hap_2_counts = case_when(
      method == 'cr' ~  w_cr,
      method == 'rbc' ~ w_rbc,
      method == 'glm' ~ w_glm,
      TRUE~NA_real_
    ),
    ssf = case_when(
      method == 'cr' ~  ssf_cr,
      method == 'rbc' ~ ssf_rbc,
      method == 'glm' ~ ssf_glm,
      TRUE~NA_real_
    ),
    wfrac = case_when(
      method == 'cr' ~  wfrac_cr,
      method == 'rbc' ~ wfrac_rbc,
      method == 'glm' ~ wfrac_glm,
      TRUE~NA_real_
    )
  ) %>% 
  mutate(n = hap_1_counts + hap_2_counts) %>% 
  relocate(unitig, hap_1_counts, hap_2_counts, n, ssf, wfrac, cluster, method) 

# Export ------------------------------------------------------------------

cat('Exporting\n')

## Marker Counts -----------------------------------------------------------


### Join in Excluded Unitigs ------------------------------------------------

marker_counts <-
  full_join(marker_counts, unitig_lengths_df, by='unitig')

### Additional Information --------------------------------------------------

marker_counts <-
  marker_counts %>%
  left_join(components_df, by='unitig') %>%
  left_join(strand_orientation_clusters_df, by=c('unitig')) %>%
  dplyr::rename(unitig_orientation = strand_cluster)

pca_perc_var_df <-
  prcomps %>% 
  map('eig') %>% 
  map(function(x) x[1:2, 2]) %>% 
  bind_rows(.id = 'cluster') %>%
  mutate(cluster = manf(cluster))

pca_perc_var_df <-
  pca_perc_var_df %>% 
  dplyr::rename(PC1_pvar = `comp 1`, PC2_pvar = `comp 2`)

marker_counts <-
  marker_counts %>%
  left_join(pca_perc_var_df, by='cluster')

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

all_unitigs <- unitig_lengths_df$unitig
stopifnot(nrow(marker_counts) == length(all_unitigs))
stopifnot(setequal(all_unitigs, marker_counts$unitig))

stopifnot(!anyNA(marker_counts$hap_1_counts))
stopifnot(!anyNA(marker_counts$hap_2_counts))

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
  tidyr::complete(cluster, lib) 

# TODO schema checks

### CSV ---------------------------------------------------------------------

readr::write_csv(library_weights, output_lib)

# Plots -------------------------------------------------------------------

library(ggplot2)

plots <- list()
cluster_plots_ssf <- list()
cluster_plots_sim <- list()
cluster_plots_ndp <- list()

## SSF Barcodes ------------------------------------------------------------

ssf_counts_df <-
  orient_counts(counts_df, strand_orientation_clusters_df) %>%
  filter(!is.na(w) & !is.na(c))

ssf_counts_df <-
  ssf_counts_df %>%
  bind_rows(anti_join(counts_df, ssf_counts_df, by = 'unitig'))

ssf_counts_df <-
  ssf_counts_df %>%
  mutate(ssf = (w - c) / n)

plot_data <-
  ssf_counts_df %>%
  left_join(cluster_df, by = 'unitig') %>%
  mutate(cluster = coalesce(cluster, unitig)) %>% 
  mutate(cluster = manf(cluster))

plot_data <-
  full_join(plot_data, long_unitigs_df) %>% 
  mutate(unitig = manf(unitig))

p <-
ggplot(plot_data) +
  geom_raster(aes(x = lib, y = unitig, fill = ssf)) +
  scale_fill_viridis_c(limits = c(-1, 1),  option = 'E', na.value = 'darkred') +
  scale_x_discrete(expand = expansion()) +
  scale_y_discrete(expand = expansion()) +
  ggtitle('SSF Values') +
  xlab('Library') +
  ylab('Cluster/Unitig') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.y = unit(2, 'points'),
    strip.text.y.left = element_text(angle = 0, size = 7)
  )

p1 <-
  p + 
  facet_grid(rows = vars(cluster),
             switch = 'y',
             scales = 'free_y',
             space = 'free') 

p2 <-
  p + 
  facet_grid(rows = vars(cluster),
             switch = 'y',
             scales = 'free_y') 

plots[['ssfb1']] <- p1
plots[['ssfb2']] <- p2


## SSF Barcode by Cluster --------------------------------------------------

# TODO height of each bar ~ size of unitig
for(clust in stringr::str_sort(pull_distinct(plot_data, cluster), numeric = TRUE)) {
  p <-
    plot_data %>%
    filter(cluster == clust) %>%
    ggplot() +
    geom_raster(aes(x = lib, y = unitig, fill = ssf)) +
    scale_fill_viridis_c(limits = c(-1, 1),  option = 'E', na.value = 'darkred') +
    facet_grid(cols = vars(cluster),
               switch = 'y',
               scales = 'free_y') +
    scale_x_discrete(expand = expansion()) +
    scale_y_discrete(expand = expansion()) +
    ggtitle('SSF Values') +
    xlab('Library') +
    ylab('Unitig') +
    theme_classic() +
    theme(
      panel.spacing.y = unit(2, 'points'),
      strip.text.y.left = element_text(angle = 0, size = 7),
      axis.text.x = element_text(angle = 360-90, hjust=0, vjust=0.5)
    )

  cluster_plots_ssf[[paste0('ssfb_', clust)]] <- p
}


## SSF Histograms ----------------------------------------------------------

## Marker Count SSF Histogram ----------------------------------------------

p <-
  marker_counts %>% 
  filter(!is.na(cluster)) %>% 
  mutate(cluster = manf(cluster)) %>% 
  ggplot() +
  geom_histogram(aes(ssf)) +
  facet_wrap(~cluster, scales = 'free_y') +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  theme_linedraw() +
  theme(
    strip.background = element_rect(fill='white'),
    strip.text = element_text(color='black')
  )

plots[['ssfh']] <- p

## Cosine Similarity Heatmaps ----------------------------------------------

make_cosine_sim_plot_data <- function(cosine_similarity_mat) {
  plot_data <-
    cosine_similarity_mat %>% 
    as_tibble(rownames = 'unitig_1') %>%
    tidyr::pivot_longer(-unitig_1, names_to = 'unitig_2', values_to = 'sim')
  
  plot_data <-
    plot_data %>% 
    full_join(long_unitigs_df, by=c('unitig_1' = 'unitig')) %>% 
    full_join(long_unitigs_df, by=c('unitig_2' = 'unitig')) 
  
  plot_data <-
    plot_data %>% 
    left_join(cluster_df, by=c('unitig_1' = 'unitig')) %>% 
    dplyr::rename(cluster_1 = cluster) %>% 
    left_join(cluster_df, by=c('unitig_2' = 'unitig')) %>% 
    dplyr::rename(cluster_2 = cluster)
  
  plot_data <-
    plot_data %>% 
    mutate(cluster_1 = coalesce(cluster_1, unitig_1)) %>% 
    mutate(cluster_2 = coalesce(cluster_2, unitig_2)) 
  
  # Floating point nonsense
  plot_data <-
    plot_data %>%
    mutate(sim = ifelse(sim > 1, 1, sim)) %>% 
    mutate(sim = ifelse(sim < -1, -1, sim))
}

make_cosine_sim_heatmap_plots <-
  function(cosine_similarity_mat,
           f_agg = function(x, ...) coverage_weighted_mean(x, mat_f = abs, ...),
           f_ind = abs,
           title = '',
           colorbar_title = '',
           limits = c(NA, NA)) {
    
  out <- list()
 
  plot_data <- make_cosine_sim_plot_data(cosine_similarity_mat)
  ### Aggregated Plots --------------------------------------------------------
  
  plot_data_tmp <-
    plot_data %>%
    group_by(cluster_1, cluster_2) %>%
    summarise(sim = f_agg(sim, na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(cluster_1), desc(cluster_2)) %>%
    mutate(
      cluster_1 = manf(cluster_1),
      cluster_2 = manf(cluster_2)
    ) %>%
    mutate(cluster_2 = factor(cluster_2, levels = rev(levels(cluster_2))))
  
  p <-
    plot_data_tmp %>%
    ggplot() +
    geom_raster(aes(x = cluster_1, y = cluster_2, fill = sim)) +
    scale_fill_viridis_c(name=colorbar_title, limits = limits, option = 'E', na.value = 'darkred') +
    ggtitle(paste('Aggregated', title)) +
    xlab('Cluster/Unitig') +
    ylab('Cluster/Unitig') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 360 - 45, hjust = 0)) 
  
  out[['acs']] <- p
  
  ### Unitig Plots ------------------------------------------------------------
  
  plot_data_tmp <-
    plot_data %>%   
    arrange(cluster_1, cluster_2) %>%
    mutate(
      cluster_1 = manf(cluster_1),
      cluster_2 = manf(cluster_2)
    ) %>% 
    mutate(
      unitig_1 = manf(unitig_1),
      unitig_2 = manf(unitig_2)
    ) %>% 
    mutate(unitig_2 = factor(unitig_2, levels = rev(levels(unitig_2))))

  p <-
    plot_data_tmp %>% 
    ggplot() +
    geom_raster(aes(x = unitig_1, y = unitig_2, fill = f_ind(sim))) +
    scale_fill_viridis_c(name=colorbar_title, limits = limits, option = 'E', na.value = 'darkred') +
    scale_x_discrete(expand = expansion()) +
    scale_y_discrete(expand = expansion()) +
    theme(axis.text.x = element_text(angle = 360 - 45, hjust = 0)) +
    ggtitle(title) +
    xlab('Unitig within Cluster') +
    ylab('Unitig within Cluster') +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.spacing = unit(2, 'points'),
      strip.text.y.left = element_text(angle = 0),
      strip.text.x = element_text(angle = 90),
      strip.text = element_text(size = 7)
    )
  
  p1 <-
    p +
    facet_grid(
      cols = vars(cluster_1),
      rows = vars(cluster_2),
      scales = 'free',
      space = 'free',
      switch = 'y'
    )
  
  p2 <-
    p +
    facet_grid(
      cols = vars(cluster_1),
      rows = vars(cluster_2),
      scales = 'free',
      # space = 'free',
      switch = 'y'
    ) 
  
  out[['cs1']] <- p1
  out[['cs2']] <- p2
  

  ## Within Cluster ----------------------------------------------------------

  plot_data_tmp <-
    plot_data_tmp %>% 
    filter(cluster_1 == cluster_2) %>% 
    filter(!(as.character(unitig_1) == as.character(cluster_1))) %>% 
    filter(!(as.character(unitig_2) == as.character(cluster_2))) 
  
  p <-
    plot_data_tmp %>% 
    ggplot() +
    geom_raster(aes(x = unitig_1, y=unitig_2, fill=f_ind(sim))) +
    facet_wrap(~cluster_1, scales = 'free') +
    scale_fill_viridis_c(name=colorbar_title, limits = limits, option = 'E', na.value = 'darkred') +
    theme_classic() +
    theme(
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 360 - 90, hjust = 0, vjust=0.5)
      ) +
    ggtitle(title) +
    xlab('Unitig') +
    ylab('Unitig') 
  
  out[['cswc']] <- p
  
  return(out)
  }

# plots_sim <-
#   make_cosine_sim_heatmap_plots(
#     cosine_similarity_mat,
#     f_agg = mean,
#     f_ind = identity,
#     title = 'Cosine Similarities',
#     colorbar_title = 'sim',
#     limits = c(-1, 1)
#   )

plots_abs_sim <-
  make_cosine_sim_heatmap_plots(
    cosine_similarity_mat,
    f_agg = mean_abs,
    f_ind = abs,
    title = 'Absolute Cosine Similarities',
    colorbar_title = 'abs(sim)',
    limits = c(0, 1)
  )


plots <- c(plots, 
           # plots_sim, 
           plots_abs_sim
           )



# Close Up Cosine Similarity Plots ----------------------------------------

make_cosine_sim_cluster_heatmap_plots <-
  function(cosine_similarity_mat,
           f_agg = function(x, ...) coverage_weighted_mean(x, mat_f = abs, ...) ,
           f_ind = abs,
           title = '',
           colorbar_title = '',
           limits = c(NA, NA)) {
    
    out <- list()
    
    plot_data <- make_cosine_sim_plot_data(cosine_similarity_mat)

    plot_data_tmp <-
      plot_data %>%   
      arrange(cluster_1, cluster_2) %>%
      mutate(
        cluster_1 = manf(cluster_1),
        cluster_2 = manf(cluster_2)
      ) %>% 
      mutate(
        unitig_1 = manf(unitig_1),
        unitig_2 = manf(unitig_2)
      ) %>% 
      mutate(unitig_2 = factor(unitig_2, levels = rev(levels(unitig_2))))
    
    plot_data_tmp <-
      plot_data_tmp %>% 
      filter(cluster_1 == cluster_2) %>% 
      group_by(cluster_1, cluster_2) %>% 
      filter(n() > 1)
    
    for(cluster in stringr::str_sort(pull_distinct(plot_data_tmp, cluster_1), numeric=TRUE)) {
      p <-
        plot_data_tmp %>% 
        filter(cluster_1 == cluster) %>% 
        ggplot() +
        geom_raster(aes(x = unitig_1, y=unitig_2, fill=f_ind(sim))) +
        facet_wrap(~cluster_1, scales = 'free') +
        scale_fill_viridis_c(name=colorbar_title, limits = limits, option = 'E', na.value = 'darkred') +
        theme_classic() +
        theme(
          axis.text = element_text(size = 5),
          axis.text.x = element_text(angle = 360 - 90, hjust = 0, vjust=0.5)
        ) +
        ggtitle(title) +
        xlab('Unitig') +
        ylab('Unitig') 
      
      out[[paste0('cs_', cluster)]] <- p
    }

    
    return(out)
  }

cluster_plots_sim <-
  make_cosine_sim_cluster_heatmap_plots(
    cosine_similarity_mat,
    f_agg = mean_abs,
    f_ind = abs,
    title = 'Absolute Cosine Similarities',
    colorbar_title = 'abs(sim)',
    limits = c(0, 1)
  )


## Phasing Planes ----------------------------------------------------------

facet_labeller <- function(x) {
  x %>% 
    left_join(pca_perc_var_df, by='cluster') %>% 
    mutate(PC1_pvar = round(PC1_pvar), PC2_pvar = round(PC2_pvar), tot=round(PC1_pvar+PC2_pvar)) %>%
    mutate(annotation=paste0('[',PC1_pvar, ':', PC2_pvar, '] [', tot, ']')) %>% 
    mutate(cluster = paste(cluster, annotation)) %>% 
    select(cluster)
}

pca_df <-
  prcomps %>% 
  map('ind') %>% 
  map('coord') %>% 
  map(function(x) {
    as_tibble(x[, 1:2], rownames='unitig') %>% 
      set_names(c('unitig_dir', 'PC1', 'PC2'))
  }) %>% 
  bind_rows(.id='cluster') %>% 
  mutate(
    y=grepl('inverted', unitig_dir)
  ) %>% 
  mutate(
    unitig = gsub('_inverted$', '', unitig_dir),
    cluster = manf(cluster)
    ) %>% 
  left_join(unitig_lengths_df) %>% 
  select(cluster, unitig, unitig_dir, length, y, everything())


projected_ww_vectors_glm_plot <-
  coefs %>%
  map_dfr(function(x) as_tibble(as.list(x)), .id='cluster') %>% 
  set_names(c('cluster', 'PC1', 'PC2'))


projected_ww_vectors_rbc_plot <-
  projected_ww_vectors_rbc %>%
  map(as.data.frame.list) %>%
  bind_rows(.id = 'cluster') %>% 
  mutate(cluster = manf(cluster))


projected_ww_vectors_plot <-
  projected_ww_vectors %>%
  map(as.list) %>% 
  map(as_tibble) %>%
  bind_rows(.id = 'cluster')%>% 
  mutate(cluster = manf(cluster))


projected_ww_vectors_plot <-
  list(`Count and Rotate` = projected_ww_vectors_plot,
       `Range Balanced Correction` = projected_ww_vectors_rbc_plot,
       `Logistic Regression` = projected_ww_vectors_glm_plot
       ) %>% 
  map(function(x) select(x, cluster, PC1, PC2)) %>% 
  bind_rows(.id='method') %>% 
  mutate(cluster = manf(cluster))


p <-
  pca_df %>%
  filter(!y) %>%
  dplyr::rename(inverted = y) %>%
  ggplot() +
  # Guidelines
  geom_hline(yintercept=0, alpha = 0.2) +
  geom_vline(xintercept=0, alpha = 0.2) +
  # WW vectors
  geom_abline(
    aes(
    slope = PC2 / PC1,
    intercept = 0,
    color = method,
    linetype = method
  ), 
  data = projected_ww_vectors_plot) +
  # WC vectors
#   geom_abline(
#     aes(
#       slope = -PC1 / PC2,
#       intercept = 0,
#       color = method,
#       linetype = method
#     ),
#     data = projected_ww_vectors_plot) +
  facet_wrap( ~ cluster, labeller = facet_labeller) +
  scale_shape_manual(guide = 'none', values = c('TRUE' = NA, 'FALSE' = 21)) +
  scale_size_area(name='Unitig Size') +
  scale_colour_brewer(name='WW Calculation Method', palette = 'Dark2') +
  scale_linetype_discrete(name='WW Calculation Method') +
  ggtitle('Fastmap Principal Components Plot') +
  theme_classic() +
  theme(panel.border = element_rect(size = 1, fill=NA))


plots[['pcp']] <- p + 
  geom_point(
    aes(
      x = PC1,
      y = PC2,
      size = length,
      shape = inverted,
    ),
    fill = 'darkgrey',
    alpha = 0.5
  )


plots[['pcpl']] <- p +
  geom_text(
    aes(
      x = PC1,
      y = PC2,
      label=unitig
    ),
    size=1,
    alpha = 0.5
  )


## Uncounted Markers -------------------------------------------------------

plot_data <-
  marker_counts %>%
  filter(length >= segment_length_threshold) %>%
  filter(hap_1_counts == 0 & hap_2_counts == 0) %>%
  mutate(clustered = !is.na(cluster))%>%
  arrange(desc(length)) %>%
  mutate(unitig = factor(unitig, levels = unique(unitig)))

p <-
  ggplot(plot_data) +
  geom_histogram(aes(length, fill = clustered), color = 'black') +
  scale_fill_discrete(name = 'Has Cluster') +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  scale_x_log10(limits = c(segment_length_threshold, NA)) +
  xlab('Log10 Length') +
  ylab('Count') +
  ggtitle(paste0('0 Marker Unitigs >= ', segment_length_threshold/1e6, 'Mbp')) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color='black', size = 1))

plots[['0_marker_unitigs_hist']] <- p


p <-
  ggplot(plot_data) +
  geom_point(aes(y = unitig, x = length, fill= clustered), shape=21) +
  scale_fill_discrete(name = 'Has Cluster') +
  scale_x_log10(limits = c(segment_length_threshold, NA)) +
  xlab('Log10 Length') +
  ylab('Unitig') +
  ggtitle(paste0('0 Marker Unitigs >= ', segment_length_threshold/1e6, 'Mbp')) +
  theme_linedraw()

plots[['0_marker_unitigs_point']] <- p


## Export ------------------------------------------------------------------

pdf(file.path(intermediate_output_dir, 'plots.pdf'), width = 15, height = 15)
iwalk(plots, function(x, nm) {
  cat('Plotting', nm, '\n')
  print(x)
})
dev.off()

pdf(file.path(intermediate_output_dir, 'plots_cluster_ssf.pdf'), width = 15, height = 15)
for(p in cluster_plots_ssf) {
  print(p)
}
dev.off()

pdf(file.path(intermediate_output_dir, 'plots_cluster_sim.pdf'), width = 15, height = 15)
for(p in cluster_plots_sim) {
  print(p)
}
dev.off()


# Warnings ----------------------------------------------------------------


for(w in WARNINGS) {
  warning(w)
}

