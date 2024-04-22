
# Import ------------------------------------------------------------------


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
  
  if(length(lines) ==0) {
    # Maybe this should be an error even?
    warning(paste('File': fastmap_file, 'has no alignments to the assembly. Is there a mismatch between the reads and assembly?'))
    out <- tibble(
      qname=character(),
      qlen=integer(),
      unitig=character(),
      strand=character(),
      lpos=integer(),
      rpos=integer(),
      n_matches=integer(),
      tstart=character() # should change it so this function returns integer here but this is what it seems to be for now.
    )
    return(out)
  }
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
                    col_types = '_ci')
  
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
      col_types = '_iiic'
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
        'tstart'
      ),
      sep = 1
    ) %>%
    select(unitig, strand, group, everything())
  
  out <-
    inner_join(header,  matches, by='group') %>% 
    select(-group)
  
  return(out)
  
}

# Counting ----------------------------------------------------------------


make_wc_matrix <- function(watson, crick, lib, unitig, min_n=10) {
  stopifnot(all_have_same_length(watson, crick, lib, unitig))
  
  # factors can cause issue ~ is the integer index or name used? When making
  # dimnames, it appears the name, when filling the matrix, it appears the
  # integer index is used.
  stopifnot(is.character(lib), is.character(unitig))
  stopifnot(is.integer(watson), is.integer(crick))
  stopifnot(!is_duplicate_pair(lib, unitig))
  
  
  w_frac <-
    ifelse((watson + crick) < min_n, NA, (watson - crick) / (watson + crick))
  
  
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
  
  # mat[is.na(mat)] <- 0
  
  return(mat)
  
}

marginalize_wc_counts <-
  function(df, do_not_marginalize = character()) {
    
    if(length(do_not_marginalize) > 1) {
      stopifnot(all(do_not_marginalize %in% df$unitig))
    }
    
    marginalized_df <-
      df %>%
      group_by(lib, unitig) %>%
      summarise(c = sum(c), w = sum(w), .groups="drop") %>%
      filter(!(unitig %in% do_not_marginalize))
    
    df <-
      df %>%
      filter(unitig %in% do_not_marginalize) %>%
      select(lib, unitig_bin, c, w) %>%
      dplyr::rename(unitig = unitig_bin)
    
    
    out <- bind_rows(marginalized_df, df)
    
    return(out)
    
  }



# Similarities ------------------------------------------------------------

make_ssf_mat <- function(counts_df, min_n) {
  
  ssf_mat <-
    counts_df %>%
    mutate(ssf = ifelse(n < min_n, NA, ifelse(n==0, NA, (w-c)/n))) %>% 
    select(lib, unitig, ssf) %>% 
    tidyr::pivot_wider(names_from = lib, values_from = ssf)
  
  rnames <- ssf_mat$unitig
  
  ssf_mat <- 
    ssf_mat %>% 
    select(-unitig) %>% 
    as.matrix()
  
  rownames(ssf_mat) <- rnames
  
  return(ssf_mat)
}


pwc_cos_sim <- function(ssf_mat, min_overlaps) {
  cosine_mat <- coop::cosine(t(ssf_mat), use = "pairwise.complete.obs")
  
  is_not_na_mat <- !is.na(ssf_mat)
  keep_sim_mat <- is_not_na_mat %*% t(is_not_na_mat)
  keep_sim_mat[keep_sim_mat < min_overlaps] <- NA
  keep_sim_mat[keep_sim_mat >= min_overlaps] <- 1
  diag(keep_sim_mat) <- 1
  cosine_mat <- cosine_mat  * keep_sim_mat
  
  return(cosine_mat)
} 



# Clustering --------------------------------------------------------------




make_contrasts <- function(x, y=NULL) {
  if(is.null(y)) {
    out <- combn(x, 2, simplify = FALSE)
  } else {
    out <-
      expand.grid(x, y, stringsAsFactors=FALSE) %>% 
      pmap(c) %>% 
      map(unname)
  }
  
  
  return(out)
}

initilize_cluster_sims_ <- function(cluster_df, sim_mat, weights_mat, sim=TRUE) {
  clustered_unitigs <- 
    cluster_df %>% 
    filter(!is.na(cluster)) %>% 
    pull_distinct(unitig)
  
  if(length(clustered_unitigs) == 0) {
    out <-
      tibble(
        cluster = character(),
        unitig = character(),
        den = double(),
        num = double(),
        sim = double()
      )
    
    return(out)
  }
  
  unclustered_unitigs <- 
    cluster_df %>% 
    filter(is.na(cluster)) %>% 
    pull_distinct(unitig)
  
  unitig_clusters <- 
    cluster_df %>% 
    filter(!is.na(cluster)) %>% 
    with(split(unitig, cluster))
  
  rws <- rownames(sim_mat) %in% unclustered_unitigs
  
  out <-
    map(unitig_clusters, function(x) {
      sims <- sim_mat[rws, x, drop = FALSE]
      wts <- weights_mat[rws, x, drop = FALSE]
      
      den <- apply(wts,        1, sum)
      num <- apply(wts * sims, 1, sum)
      
      return(
        bind_rows(
          # Need null tibble, as unitig column is not created properly if
          # length(den) == 0)
          tibble(unitig = names(den), den=den, num=num),
          tibble(unitig = character(), den = double(), num=double())
        )
      )
    })
  
  out <-
    out %>% 
    bind_rows(.id='cluster') 
  
  if(sim) {
    out <-
      out %>% 
      mutate(sim = ifelse(den > 0, num/den, NA))
  }
  
  return(out)
}

cluster_unitigs_2 <-
  function(cluster_df,
           similarity_mat,
           weights,
           sim_threshold_cu = 0.5,
           sim_threshold_uu = 0.5,
           new_cluster_id = 'C',
           new_cluster_ix = 0) {
    
    included_unitigs <- pull_distinct(cluster_df, unitig)
    sim_mat <- similarity_mat[included_unitigs, included_unitigs, drop=FALSE]
    # is this right?
    not_na_mat <- !is.na(sim_mat)
    
    # Weights is a neamed vector, where each name is a unitig, one for each unitig
    # under consideration in the whole script (all the long unitigs)
    weights <- weights[rownames(not_na_mat)]
    # vector by matrix multiplication is by column
    weights_mat <- t(t(not_na_mat) * weights)
    sim_mat[!not_na_mat] <- 0
    
    stopifnot(all(rownames(weights_mat) == rownames(sim_mat)))
    stopifnot(all(colnames(weights_mat) == colnames(sim_mat)))
    
    cu_sim_df <- initilize_cluster_sims_(cluster_df,  sim_mat, weights_mat)
    
    # Having to filter to unclustered unitigs feels clunky
    unclustered_unitigs <- 
      cluster_df %>% 
      filter(is.na(cluster)) %>% 
      pull_distinct(unitig)
    sim_mat <- sim_mat[unclustered_unitigs, unclustered_unitigs, drop=FALSE]
    weights_mat <- weights_mat[unclustered_unitigs, unclustered_unitigs, drop=FALSE]
    
    while(TRUE) {
      # if(ncol(sim_mat) != sum(is.na(cluster_df$cluster))) browser()
      stopifnot(nrow(sim_mat) == ncol(sim_mat), ncol(sim_mat) == sum(is.na(cluster_df$cluster)))
      stopifnot(nrow(weights_mat) == ncol(weights_mat), ncol(weights_mat) == sum(is.na(cluster_df$cluster)))
      
      # Unitig -> cluster
      max_sim <-
        cu_sim_df %>% 
        filter(sim >= sim_threshold_cu) %>% 
        slice_max(sim, n=1, with_ties = FALSE, na_rm = TRUE)
      
      if(nrow(max_sim) > 0) {
        
        cur_utg <- max_sim$unitig
        cur_clust <- max_sim$cluster
        
        message(
          'Adding: ',
          cur_utg,
          ' to cluster: ',
          cur_clust,
          ', sim: ',
          max_sim$sim
        )
        
        cluster_df <-
          mutate(cluster_df, cluster = ifelse(unitig == cur_utg, cur_clust, cluster))
        
        cu_sim_df <-
          filter(
            cu_sim_df, 
            !(unitig == cur_utg)
          )
        
        sims_to_add <-
          initilize_cluster_sims_(
            filter(cluster_df, is.na(cluster) | unitig == cur_utg),
            sim_mat,
            weights_mat,
            sim=FALSE
          ) 
        
        # this relies on perfect order.
        ix <- which(pull(cu_sim_df, cluster) == cur_clust)
        # if(nrow(sims_to_add) ==0) browser()
        stopifnot(all(pull(cu_sim_df, unitig)[ix] == pull(sims_to_add, unitig)))
        cu_sim_df$den[ix] <- cu_sim_df$den[ix] + sims_to_add$den
        cu_sim_df$num[ix] <- cu_sim_df$num[ix] + sims_to_add$num
        
        
        ix_to_remove <- which(rownames(sim_mat) == cur_utg)
        sim_mat <- sim_mat[-ix_to_remove, -ix_to_remove, drop=FALSE]
        weights_mat <- weights_mat[-ix_to_remove, -ix_to_remove, drop=FALSE]
        
        next
      }
      
      # Unitig -> Unitig
      
      # From https://www.tutorialspoint.com/how-to-find-the-row-and-column-index-for-upper-triangular-matrix-elements-in-r
      
      upper_tri_ix_mat <- which(upper.tri(sim_mat,diag=FALSE),arr.ind=TRUE)
      upper_tri_sims <- sim_mat[upper_tri_ix_mat]
      
      which_max <- which.max(upper_tri_sims)
      max_sim <- upper_tri_sims[which_max]
      
      
      if(length(max_sim) > 0 && max_sim >= sim_threshold_uu) {
        
        max_ix <- upper_tri_ix_mat[which_max, , drop=FALSE]
        cluster_members <- rownames(sim_mat)[max_ix]
        new_cluster_name <- paste0(new_cluster_id, new_cluster_ix)
        new_cluster_ix <- new_cluster_ix + 1
        
        message(
          'Creating a new cluster: ',
          new_cluster_name,
          ' with members: ',
          paste(cluster_members, collapse=', '),
          ', sim: ',
          max_sim
        )
        
        cluster_df <-
          mutate(cluster_df, cluster = ifelse(unitig %in% cluster_members, new_cluster_name, cluster))
        
        sims_to_add <-
          initilize_cluster_sims_(
            filter(cluster_df, is.na(cluster) | cluster == new_cluster_name),
            sim_mat,
            weights_mat 
          )
        
        cu_sim_df <-
          bind_rows(cu_sim_df, sims_to_add)
        
        cu_sim_df <-  filter(cu_sim_df, !(unitig %in% cluster_members))
        sim_mat <- sim_mat[-as.vector(max_ix), -as.vector(max_ix), drop=FALSE]
        weights_mat <- weights_mat[-as.vector(max_ix), -as.vector(max_ix), drop=FALSE]
        
        next
      } 
      
      
      break 
    }
    
    return(cluster_df)
    
  }



merge_similar_clusters_2 <- function(cluster_df,
                                     similarity_mat,
                                     weights,
                                     sim_threshold = 0.5) {
  
  included_unitigs <- pull_distinct(cluster_df, unitig)
  sim_mat <- similarity_mat[included_unitigs, included_unitigs, drop=FALSE]
  not_na_mat <- !is.na(sim_mat)
  # Weights is a nameed vector, where each name is a unitig, one for each unitig
  # under consideration in the whole script (all the long unitigs)
  weights <- weights[rownames(not_na_mat)]
  
  weights_mat <- not_na_mat * outer(weights, weights)
  sim_mat[!not_na_mat] <- 0
  
  while(TRUE) {
    
    sets <-
      cluster_df %>% 
      filter(!is.na(cluster)) %>% 
      with(split(unitig, cluster))
    
    contrasts <- 
      cluster_df %>% 
      filter(!is.na(cluster)) %>% 
      pull_distinct(cluster) %>% 
      make_contrasts()
    
    cluster_similarities <-
      map_dfr(contrasts, function(pair) {
        id_1 <- pair[1]
        id_2 <- pair[2]
        nms_1 <- sets[[id_1]]
        nms_2 <- sets[[id_2]]
        
        sim <-
          weighted.mean(
            sim_mat[nms_1, nms_2, drop=FALSE],
            weights_mat[nms_1, nms_2, drop=FALSE]
          )
        
        return(tibble(clust_1 = id_1, clust_2 = id_2, sim=sim))
        
      })
    
    if(all(is.na(cluster_similarities$sim))) {
      cat('No valid similarity scores\n')
      break
    }
    
    max_contrast <-
      cluster_similarities %>% 
      slice_max(sim, n = 1, with_ties = FALSE)
    
    max_sim <- max_contrast$sim
    # similarity_threshold chosen by reviewing ~ 30 HGSVC assemblies. A value of
    # 0.4 appears to also work, as there seems to be a sharp gulf where unitigs
    # from different clusters tend to have similarity no more than ~ 0.21,
    # while those from the same have similarity >0.5 there was one case where an
    # X and  chromosome had similarity ~ 0.447. However, will stay at 0.5 as
    # buffer for precision?
    if(max_sim > sim_threshold) {
      # sort to ensure that sex cluster (sex > LG) is not overwritten
      clusters_to_merge <- sort(c(max_contrast$clust_1, max_contrast$clust_2))
      clust_1 <- clusters_to_merge[1]
      clust_2 <- clusters_to_merge[2]
      
      cat('Merging ', clust_1, ' into ', clust_2, ', cosine similarity: ', max_sim,
          '\n')
      
      cluster_df <-
        cluster_df %>% 
        mutate(cluster = ifelse(cluster == clust_1, clust_2, cluster))
      
      next
    } 
    
    break
  }
  
  return(cluster_df)
  
}

merge_similar_clusters_on_components_2 <- 
  function(cluster_df,
           similarity_mat,
           weights,
           components_df,
           sim_threshold = 0.4) {
    
    included_unitigs <- pull_distinct(cluster_df, unitig)
    sim_mat <- similarity_mat[included_unitigs, included_unitigs, drop=FALSE]
    not_na_mat <- !is.na(sim_mat)
    # Weights is a nameed vector, where each name is a unitig, one for each unitig
    # under consideration in the whole script (all the long unitigs)
    weights <- weights[rownames(not_na_mat)]
    
    weights_mat <- not_na_mat * outer(weights, weights)
    sim_mat[!not_na_mat] <- 0
    
    while(TRUE) {
      
      components_with_multiple_clusters <-
        components_df %>% 
        inner_join(cluster_df, by='unitig') %>%
        filter(!is.na(cluster)) 
      
      components_with_multiple_clusters <-
        components_with_multiple_clusters %>% 
        group_by(component) %>% 
        filter(length(unique(cluster)) > 1) %>% 
        distinct(cluster) %>% 
        ungroup() %>% 
        arrange(component, cluster)
      
      if(nrow(components_with_multiple_clusters) < 1) {
        cat('No components with multiple clusters')
        break
      }
      
      sets <- 
        cluster_df %>% 
        semi_join(components_with_multiple_clusters, by='cluster') %>% 
        with(split(unitig, cluster))
      
      contrasts <-
        components_with_multiple_clusters %>% 
        group_split(component) %>% 
        map(function(x) make_contrasts(x$cluster)) %>% 
        flatten()
      
      cluster_similarities <-
        map_dfr(contrasts, function(pair) {
          id_1 <- pair[1]
          id_2 <- pair[2]
          nms_1 <- sets[[id_1]]
          nms_2 <- sets[[id_2]]
          
          sim <-
            weighted.mean(
              sim_mat[nms_1, nms_2, drop=FALSE],
              weights_mat[nms_1, nms_2, drop=FALSE]
            )
          
          return(tibble(clust_1 = id_1, clust_2 = id_2, sim=sim))
          
        })
      
      
      if(all(is.na(cluster_similarities$sim))) {
        cat('No valid similarity scores\n')
        break
      }
      
      max_contrast <-
        cluster_similarities %>% 
        slice_max(sim, n = 1, with_ties = FALSE)
      
      max_sim <- max_contrast$sim
      # similarity_threshold chosen by reviewing ~ 30 HGSVC assemblies. A value of
      # 0.4 appears to also work, as there seems to be a sharp gulf where unitigs
      # from different clusters tend to have similarity no more than ~ 0.21,
      # while those from the same have similarity >0.5 there was one case where an
      # X and  chromosome had similarity ~ 0.447. However, will stay at 0.5 as
      # buffer for precision?
      if(max_sim > sim_threshold) {
        # sort to ensure that sex cluster (sex > LG) is not overwritten
        clusters_to_merge <- sort(c(max_contrast$clust_1, max_contrast$clust_2))
        clust_1 <- clusters_to_merge[1]
        clust_2 <- clusters_to_merge[2]
        
        cat('Merging ', clust_1, ' into ', clust_2, ', cosine similarity: ', max_sim,
            '\n')
        
        cluster_df <-
          cluster_df %>% 
          mutate(cluster = ifelse(cluster == clust_1, clust_2, cluster))
        next
      } 
      
      break
    }
    
    return(cluster_df)
  }


get_strand_orientation_clusters <- function(similarity_mat, weights) {
  
  # included_unitigs <- pull_distinct(cluster_df, unitig)
  # sim_mat <- similarity_mat[included_unitigs, included_unitigs, drop=FALSE]
  sim_mat <- similarity_mat
  
  cluster_df <- 
    tibble(unitig = rownames(sim_mat), cluster = rownames(sim_mat))
  
  not_na_mat <- !is.na(sim_mat)
  
  # Weights is a neamed vector, where each name is a unitig, one for each unitig
  # under consideration in the whole script (all the long unitigs)
  weights <- weights[rownames(not_na_mat)]
  weights_mat <- not_na_mat * outer(weights, weights)
  diag(weights_mat) <- 0 # setting diagonal to 0 is important for the stacking
  
  sim_mat[!not_na_mat] <- 0
  
  stopifnot(all(rownames(weights_mat) == rownames(sim_mat)))
  stopifnot(all(colnames(weights_mat) == colnames(sim_mat)))
  
  # keep track of two matrices: numerator and denominator, together make the sim matrix.
  
  den_mat <- weights_mat
  num_mat <- weights_mat * sim_mat
  while(TRUE) {
    # if(ncol(sim_mat) != sum(is.na(cluster_df$cluster))) browser()
    # stopifnot(
    #   all(rownames(den_mat) == rownames(num_mat)), 
    #   all(rownames(den_mat) == rownames(sim_mat))
    # )
    # stopifnot(
    #   all(colnames(den_mat) == colnames(num_mat)), 
    #   all(colnames(den_mat) == colnames(sim_mat))
    # )
    # stopifnot(
    #   all(colnames(den_mat) == rownames(den_mat)),
    #   all(colnames(num_mat) == rownames(num_mat)),
    #   all(colnames(sim_mat) == rownames(sim_mat))
    # )
    
    n_clusters <- n_distinct(cluster_df$cluster)
    if(n_clusters <= 2) break
    
    # From https://www.tutorialspoint.com/how-to-find-the-row-and-column-index-for-upper-triangular-matrix-elements-in-r
    
    upper_tri_ix_mat <- which(upper.tri(sim_mat,diag=FALSE),arr.ind=TRUE)
    upper_tri_sims <- sim_mat[upper_tri_ix_mat]
    
    which_max <- which.max(upper_tri_sims)
    max_sim <- upper_tri_sims[which_max]
    
    if(length(max_sim) > 0) {
      # browser()
      max_ix <- upper_tri_ix_mat[which_max, , drop=FALSE]
      cluster_members <- rownames(sim_mat)[max_ix]
      
      message(
        'Merging clusters: ',
        paste(cluster_members, collapse=', '),
        ', sim: ',
        max_sim
      )
      
      remove_ix <- max(max_ix)
      keep_ix <- min(max_ix)
      remove_utg <- rownames(sim_mat)[remove_ix]
      keep_utg <- rownames(sim_mat)[keep_ix]
      
      cluster_df <-
        mutate(cluster_df, cluster = ifelse(cluster == remove_utg, keep_utg, cluster))
      
      remove_ix <- max(max_ix)
      keep_ix <- min(max_ix)
      
      den_mat[keep_ix, ] <- den_mat[keep_ix, ] + den_mat[remove_ix, ]
      den_mat[, keep_ix] <- den_mat[, keep_ix] + den_mat[, remove_ix]
      
      num_mat[keep_ix, ] <- num_mat[keep_ix, ] + num_mat[remove_ix, ]
      num_mat[, keep_ix] <- num_mat[, keep_ix] + num_mat[, remove_ix]
      
      sim_mat[keep_ix, ] <- num_mat[keep_ix, ] / den_mat[keep_ix, ]
      sim_mat[, keep_ix] <- num_mat[, keep_ix] / den_mat[, keep_ix]
      
      # Maybe faster to not remove the rows?
      den_mat <- den_mat[-remove_ix, -remove_ix, drop=FALSE]
      num_mat <- num_mat[-remove_ix, -remove_ix, drop=FALSE]
      sim_mat <- sim_mat[-remove_ix, -remove_ix, drop=FALSE]
      
      # den_mat[remove_ix, ] <- NA
      # den_mat[, remove_ix] <- NA
      # num_mat[remove_ix, ] <- NA
      # num_mat[, remove_ix] <- NA
      # sim_mat[remove_ix, ] <- NA
      # sim_mat[, remove_ix] <- NA
      next
    } 
    
    
    break 
  }
  
  # if(n_clusters != 2) browser()
  stopifnot(n_clusters == 2)
  
  out <-
    cluster_df %>% 
    mutate(cluster = as.numeric(as.factor(cluster))) %>% 
    with(set_names(cluster, unitig))
  
  return(out)
}





propagate_one_cluster_components <- function(cluster_df, components_df) {
  
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
  
  one_cluster_components_df <-
    cluster_df %>%
    left_join(components_df, by = 'unitig') %>%
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
  
  cluster_df <-
    cluster_df %>%
    mutate(cluster = ifelse(
      unitig %in% names(one_cluster_component_unitigs),
      one_cluster_component_unitigs[unitig],
      cluster
    ))
  
  return(cluster_df)
  
}

remove_small_clusters <- function(cluster_df, unitig_lengths_df, threshold=1e7) {
  
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
  
  return(cluster_df)
}



remove_small_clusters_on_components <-
  function(cluster_df,
           unitig_lengths_df,
           components_df,
           threshold = 0.02) {
    
    component_sizes <-
      components_df %>%
      left_join(unitig_lengths_df, by = 'unitig') %>%
      count(component, wt = length, name = 'component_size')
    
    component_cluster_fractions <-
      cluster_df %>%
      filter(!is.na(cluster)) %>%
      left_join(components_df, by = 'unitig') %>%
      left_join(unitig_lengths_df, by = 'unitig') %>%
      group_by(component, cluster) %>%
      summarise(length = sum(length), .groups = 'drop')
    
    component_cluster_fractions <-
      component_cluster_fractions %>%
      left_join(component_sizes, by = 'component') %>%
      mutate(frac = length / component_size)
    
    component_clusters_to_remove <-
      component_cluster_fractions %>%
      filter(frac <= threshold) %>%
      distinct(component, cluster)
    
    unitigs_to_decluster <-
      cluster_df %>%
      left_join(components_df, by = 'unitig') %>%
      semi_join(component_clusters_to_remove, by = c('component', 'cluster')) %>%
      pull_distinct(unitig)
    
    cluster_df <-
      cluster_df %>%
      mutate(cluster = ifelse(unitig %in% unitigs_to_decluster, NA, cluster))
    
    return(cluster_df)
  }
# Phasing -----------------------------------------------------------------
orient_counts <- function(counts_df, strand_orientation_clusters_df){
  out <-
    counts_df %>%
    left_join(strand_orientation_clusters_df, by='unitig')
  
  out <-
    out %>% 
    mutate(
      w_temp = ifelse(strand_cluster == 1, w, c),
      c_temp = ifelse(strand_cluster == 1, c, w),
    ) %>% 
    mutate(c = c_temp, w = w_temp) %>% 
    select(-c_temp, -w_temp, -strand_cluster)
  
  return(out)
  
}

# Misc --------------------------------------------------------------------

mean_abs <- function(x, ...) {
  mean(abs(x), ...)
}

# contiBAIT ---------------------------------------------------------------

preprocessStrandTable <-
  function (strandTable, strandTableThreshold = 0.8, filterThreshold = 0.8, 
            orderMethod = "libsAndConc", lowQualThreshold = 0.9, verbose = TRUE, 
            minLib = 10) 
  {
    
    
    # strandTable	
    # data.frame containing the strand table to use as input
    # 
    # strandTableThreshold	
    # threshold at which to call a contig WW or CC rather than WC
    # 
    # filterThreshold	
    # maximum number of libraries a contig can be NA or WC in
    # 
    # orderMethod	
    # the method to oder contigs. currently libsAndConc only option. Set to FALSE to not order contigs based on library quality
    # 
    # lowQualThreshold	
    # background threshold at which to toss an entire library. If NULL, function will not make an overall assessment of library quality. Very chimeric assemblies can appear low quality across all libraries.
    # 
    # verbose	
    # messages written to terminal
    # 
    # minLib	
    # minimum number of libraries a contig must be present in to be included in the output
    # 
    strandTableLength <- nrow(strandTable)
    lowQualList <- data.frame(library = vector(), quality = vector())
    qualList <- lowQualList
    if (!(is.null(lowQualThreshold))) {
      if (verbose) {
        message("-> Checking for high quality libraries")
      }
      for (col in seq_len(ncol(strandTable))) {
        libName <- colnames(strandTable, do.NULL = FALSE)[col]
        backGroundC <- abs(mean(strandTable[, col][which(strandTable[, 
                                                                     col] < -0.6)], na.rm = TRUE))
        backGroundW <- abs(mean(strandTable[, col][which(strandTable[, 
                                                                     col] > 0.6)], na.rm = TRUE))
        libraryQual <- round((backGroundC + backGroundW)/2, 
                             digits = 3)
        colQual <- data.frame(library = libName, quality = libraryQual)
        if (libraryQual < lowQualThreshold || libraryQual == 
            "NaN") {
          if (libraryQual == "NaN" & verbose) {
            message(paste("    -> ", libName, " has insufficient reads. Removing", 
                          sep = ""))
          }
          else {
            if (verbose) {
              message(paste("    -> ", libName, " has high background (", 
                            (1 - libraryQual) * 100, " %). Removing", 
                            sep = ""))
            }
          }
          lowQualList <- rbind(lowQualList, colQual)
        }
        else {
          qualList <- rbind(qualList, colQual)
        }
      }
      if (nrow(lowQualList) == ncol(strandTable)) {
        warning("-> WARNING! ALL LIBRARIES ARE OF LOW QUALITY!! UNABLE TO REMOVE HIGH BACKGROUND LIBRARIES!")
      }
      else if (nrow(lowQualList) == 0 & verbose) {
        message("-> All libraries of good quality")
      }
      else {
        if (verbose) {
          stCol <- ncol(strandTable)
          lqRow <- nrow(lowQualList)
          message(paste("-> Removed ", lqRow, " libraries from a total of ", 
                        stCol, ". ", stCol - lqRow, " remaining (", 
                        round((stCol - lqRow)/stCol * 100, digits = 1), 
                        "%)", sep = ""))
        }
        strandTable <- strandTable[, !(colnames(strandTable) %in% 
                                         lowQualList[, 1])]
      }
    }
    rawTable <- strandTable
    strandTable[strandTable >= strandTableThreshold] <- 1
    strandTable[strandTable <= -strandTableThreshold] <- 3
    strandTable[strandTable < strandTableThreshold & strandTable > 
                  -strandTableThreshold] <- 2
    preFilterData <- function(strandTable, filterThreshold, 
                              minLib, onlyWC = FALSE) {
      strandTable <- strandTable[which(apply(strandTable, 
                                             1, function(x) {
                                               length(which(!is.na(x)))
                                             } >= minLib)), ]
      strandTable <- strandTable[, which(apply(strandTable, 
                                               2, function(x) {
                                                 length(which(is.na(x)))
                                               } <= length(x) * filterThreshold))]
      WCvaluesCon <- apply(strandTable, 1, function(row) {
        sum(row == 2, na.rm = TRUE)/sum(is.element(row, 
                                                   c(1, 2, 3)), na.rm = TRUE)
      })
      WCvaluesLib <- apply(strandTable, 2, function(col) {
        sum(col == 2, na.rm = TRUE)/sum(is.element(col, 
                                                   c(1, 2, 3)), na.rm = TRUE)
      })
      if (onlyWC) {
        strandTable <- strandTable[which(WCvaluesCon > filterThreshold), 
        ]
      }
      else {
        strandTable <- strandTable[which(WCvaluesCon <= 
                                           filterThreshold), ]
        strandTable <- strandTable[, which(WCvaluesLib <= 
                                             filterThreshold)]
      }
      return(strandTable)
    }
    strandTableAWC <- preFilterData(strandTable, filterThreshold, 
                                    minLib, onlyWC = TRUE)
    AWCrange <- rownames(strandTableAWC)
    AWClengths <- sub(".*:", "", AWCrange)
    # AWCrange <- GRanges(seqnames = sub(":.*", "", AWCrange), 
    #                     IRanges(start = as.numeric(sub("-.*", "", AWClengths)), 
    #                             end = as.integer(sub(".*-", "", AWClengths))))
    strandTable <- preFilterData(strandTable, filterThreshold, 
                                 minLib)
    rawTable <- rawTable[rownames(strandTable), ]
    rawTable <- rawTable[, colnames(strandTable)]
    if (orderMethod == "libsAndConc") {
      if (verbose) {
        message("-> Computing QA measures for contigs and sorting by best quality first")
      }
      contigNAs <- apply(strandTable, 1, function(x) {
        length(which(!is.na(x)))
      })/ncol(strandTable)
      computeOneAgreement <- function(contigName) {
        contigCall <- strandTable[contigName, ]
        contigRaw <- rawTable[contigName, ]
        mean(ifelse(contigCall == 2, (strandTableThreshold - 
                                        abs(contigRaw))/strandTableThreshold, (abs(contigRaw) - 
                                                                                 strandTableThreshold)/(1 - strandTableThreshold))[!is.na(contigCall)])
      }
      contigAgreement <- sapply(rownames(strandTable), computeOneAgreement)
      contigQA <- contigAgreement * contigNAs
      strandTable <- strandTable[names(sort(contigQA, decreasing = TRUE)), 
      ]
    }
    strandTable[is.na(strandTable)] <- NA
    # strandTable <- StrandStateMatrix(strandTable)
    return(list(strandMatrix = strandTable, qualList = qualList, 
                lowQualList = lowQualList, AWCcontigs = AWCrange, contigQA=contigQA))
  }


findSexGroups <- function(linkageGroupList, strandStateMatrix, callThreshold=0.2)
{
  
  createTable <- function(group, allStrands)												
  {	
    groupStrands <- allStrands[group,, drop=FALSE]
    
    tables <- sapply(1:ncol(groupStrands), 
                     function(y) sapply(c(paste(c(1,2,3), collapse="|"), 2), 
                                        function(x) length(grep(x, groupStrands[,y]))))
    tabTots <- apply(tables, 1, sum)
    tabTots <- tabTots[2]/tabTots[1]
    
  }
  
  LGconsensus <- data.frame(do.call(rbind, lapply(linkageGroupList, createTable, strandStateMatrix))) 
  LGconsensus <- LGconsensus[order(LGconsensus),1, drop=FALSE]
  LGconsensus <- rownames(LGconsensus[which(LGconsensus[,1] <= callThreshold), ,drop=FALSE])
  
  if(length(LGconsensus) == 0)
  {
    warning('NO SEX LINKAGE GROUPS FOUND.')
    return(linkageGroupList)
  }else{
    LGNames <- replace(names(linkageGroupList), names(linkageGroupList) == LGconsensus, paste("sex", LGconsensus, sep='_')) 
    names(linkageGroupList) <- LGNames
    return(linkageGroupList)
  }
}