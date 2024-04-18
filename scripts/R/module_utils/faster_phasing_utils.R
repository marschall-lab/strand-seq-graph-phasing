
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



pwc_hclust_n <- function(similarity_mat, weights, n=22) {
  
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
    if(n_clusters <= n) break
    
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
  
  # if(n_clusters != n) browser()
  # stopifnot(n_clusters == n)
  
  out <-
    cluster_df %>% 
    mutate(cluster = as.numeric(as.factor(cluster))) %>% 
    with(set_names(cluster, unitig))
  
  return(out)
}

