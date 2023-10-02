
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

# Clustering --------------------------------------------------------------

merge_similar_clusters_on_components <- function(counts_df, cluster_df, components_df, similarity_threshold =0.40 ) {
  # TODO
  wfrac_matrix <- with(counts_df, make_wc_matrix(w, c, lib, unitig))
  
  any_merged <- TRUE
  while(any_merged) {
    any_merged <- FALSE
    
    components_with_multiple_clusters <-
      components_df %>% 
      left_join(cluster_df, by='unitig') %>%
      filter(!is.na(cluster)) %>% 
      group_by(component) %>% 
      filter(length(unique(cluster)) > 1) %>% 
      pull_distinct(component)
    
    if(length(components_with_multiple_clusters) < 1) {
      cat('No components with multiple clusters')
    }
    
    for(cmp in components_with_multiple_clusters) {
      cmp_clusters <-
        components_df %>% 
        left_join(cluster_df, by='unitig') %>% 
        filter(component == cmp) %>% 
        filter(!is.na(cluster)) %>% 
        pull_distinct(cluster)
      
      cluster_unitigs <-
        cmp_clusters %>%
        set_names() %>%
        map(function(x) {
          cluster_df %>%
            filter(cluster == x) %>%
            pull(unitig)
        })
      
      # TODO weighted by vector length?
      component_similarities <-
        wfrac_matrix[flatten_chr(cluster_unitigs), ] %>% 
        cosine_similarity() %>% 
        abs()
      
      
      # Average within and between clusters
      n_clusters <- length(cmp_clusters)
      
      clust_sim <- matrix(nrow=n_clusters, ncol=n_clusters)
      dimnames(clust_sim) <- list(cmp_clusters, cmp_clusters)
      
      for(i in cmp_clusters) {
        for(j in cmp_clusters) {
          i_unitigs <- cluster_unitigs[[i]]
          j_unitigs <- cluster_unitigs[[j]]
          #TODO NA handling of values? What if all NA?
          val <- mean(component_similarities[i_unitigs, j_unitigs], na.rm = TRUE)
          clust_sim[i,j] <- clust_sim[j,i] <- val
        }
      }
      
      if(all(is.na(clust_sim[upper.tri(clust_sim)]))) {
        cat('No valid similarity scores\n')
        break
      }
      #TODO NA handling of values? What if all NA
      max_sim <- max(clust_sim[upper.tri(clust_sim)], na.rm=TRUE)
      
      max_ix <-
        which(clust_sim == max_sim, arr.ind = TRUE)
      
      # sort to ensure that sex cluster (sex > LG) is not overwritten
      merge_clusters <- c(
        cmp_clusters[max_ix[1,1]],
        cmp_clusters[max_ix[1,2]]
      ) %>% 
        sort()
      
      clust_1 <- merge_clusters[1]
      clust_2 <- merge_clusters[2]
      
      # similarity_threshold chosen by reviewing ~ 30 HGSVC assemblies. A value of
      # 0.4 appears to also work, as there seems to be a sharp gulf where unitigs
      # from different clusters tend to have similarity no more than ~ 0.21,
      # while those from the same have similarity >0.5 there was one case where an
      # X and  chromosome had similarity ~ 0.447. However, will stay at 0.5 as
      # buffer for precision?
      if(max_sim > similarity_threshold) {
        any_merged <- TRUE
        
        cat('Merging ',
            clust_1,
            ', ',
            clust_2,
            ', cosine similarity: ',
            max_sim,
            '\n')
        
        cluster_df <-
          cluster_df %>% 
          mutate(cluster = ifelse(cluster == clust_1, clust_2, cluster))
        
        # a cheap fix for a bug that can occur if a cluster exists on two
        # components, and after merging due to one component, the second component
        # no longer has multiple clusters
        break
        
      } else {
        cat(
          'Not merging',
          paste(cmp_clusters, collapse=','),
          ' greatest similarity between',
          clust_1,
          ', ',
          clust_2,
          ', cosine similarity: ',
          max_sim,
          '\n'
        )
      }
    }
  }
  
  return(cluster_df)
}


merge_similar_clusters <- function(counts_df, cluster_df, similarity_threshold =0.50, min_n=10, min_overlaps=5 ) {
  cat('Creating wfrac matrix \n')
  wfrac_matrix <- 
    counts_df %>% 
    semi_join(cluster_df, by='unitig') %>% 
    with(make_wc_matrix(w, c, lib, unitig, min_n=min_n))
  
  cat('Calculating similarities \n')
  similarities <-
    wfrac_matrix %>%
    cosine_similarity(min_overlaps=min_overlaps) %>% 
    abs()
  
  any_merged <- TRUE
  cat('Merging clusters \n')
  while(any_merged) {
    any_merged <- FALSE

      clusters <-
        cluster_df %>% 
        filter(!is.na(cluster)) %>% 
        pull_distinct(cluster) %>% 
        set_names()
      
      cluster_unitigs <-
        clusters %>%
        map(function(x) {
          cluster_df %>%
            filter(cluster == x) %>%
            pull(unitig)
        })
      

      # Average within and between clusters
      n_clusters <- length(clusters)
      
      clust_sim <- matrix(nrow=n_clusters, ncol=n_clusters)
      dimnames(clust_sim) <- list(clusters, clusters)
      
      for(i in clusters) {
        for(j in clusters) {
          i_unitigs <- cluster_unitigs[[i]]
          j_unitigs <- cluster_unitigs[[j]]
          #TODO NA handling of values? What if all NA?
          val <- mean(similarities[i_unitigs, j_unitigs], na.rm = TRUE)
          clust_sim[i,j] <- clust_sim[j,i] <- val
        }
      }
      
      if(all(is.na(clust_sim[upper.tri(clust_sim)]))) {
        cat('No valid similarity scores\n')
        break
      }
      #TODO NA handling of values? What if all NA
      max_sim <- max(clust_sim[upper.tri(clust_sim)], na.rm=TRUE)
      
      max_ix <-
        which(clust_sim == max_sim, arr.ind = TRUE)
      
      clust_1 <- clusters[max_ix[1,1]]
      clust_2 <- clusters[max_ix[1,2]]
      
      if(max_sim > similarity_threshold) {
        any_merged <- TRUE
        
        cat('Merging ',
            clust_1,
            ', ',
            clust_2,
            ', cosine similarity: ',
            max_sim,
            '\n')
        
        cluster_df <-
          cluster_df %>% 
          mutate(cluster = ifelse(cluster == clust_1, clust_2, cluster))
        
      } else {
        cat(
          'Not merging',
          paste(clusters, collapse=','),
          ' greatest similarity between',
          clust_1,
          ', ',
          clust_2,
          ', cosine similarity: ',
          max_sim,
          '\n'
        )
      }
    }
  
  return(cluster_df)
}


merge_similar_clusters2_ <- function(counts_df, cluster_df, similarity_threshold =0.50, min_n=10, min_overlaps=5 ) {
  
  if(nrow(cluster_df) == 1) {
    out <-
      cluster_df %>% 
      mutate(num_loops=1) 
    
    return(out)
  }
  wfrac_matrix <- with(counts_df, make_wc_matrix(w, c, lib, unitig, min_n=min_n))
  any_merged <- TRUE
  num_loops <- 0
  while(any_merged) {
    any_merged <- FALSE
    
    clusters <-
      cluster_df %>% 
      filter(!is.na(cluster)) %>% 
      pull_distinct(cluster) %>% 
      set_names()
    
    cluster_unitigs <-
      clusters %>%
      map(function(x) {
        cluster_df %>%
          filter(cluster == x) %>%
          pull(unitig)
      })
    
    # TODO Should this abs be removed? Why does it work even when abs is there?
    # That is very confusing. I think it works because dissimilar clusters ended
    # up orthogonal no one another, where the "error" (I think it is an error)
    # of applying abs() would have little effect.
    similarities <-
      wfrac_matrix[flatten_chr(cluster_unitigs), ,drop=FALSE] %>%
      cosine_similarity(min_overlaps=min_overlaps) # %>% 
      # abs()
    
    # Average within and between clusters
    n_clusters <- length(clusters)
    
    clust_sim <- matrix(nrow=n_clusters, ncol=n_clusters)
    dimnames(clust_sim) <- list(clusters, clusters)
    
    for(i in clusters) {
      for(j in clusters) {
        i_unitigs <- cluster_unitigs[[i]]
        j_unitigs <- cluster_unitigs[[j]]
        #TODO NA handling of values? What if all NA?
        val <- mean(similarities[i_unitigs, j_unitigs], na.rm = TRUE)
        clust_sim[i,j] <- clust_sim[j,i] <- val
      }
    }
    
    if(all(is.na(clust_sim[upper.tri(clust_sim)]))) {
      cat('No valid similarity scores\n')
      break
    }
    #TODO NA handling of values? What if all NA
    max_sim <- max(clust_sim[upper.tri(clust_sim)], na.rm=TRUE)
    
    # TODO something to handle comparison of floats here, getting some 1.000 !=
    # 1 results
    max_ix <-
      which(clust_sim == max_sim, arr.ind = TRUE)
    
    # If similarity is 1, values on the diagonal can be selected. Remove diagonals
    max_ix <- max_ix[max_ix[, 'row'] != max_ix[, 'col'] ,]
    
    clust_1 <- clusters[max_ix[1,1]]
    clust_2 <- clusters[max_ix[1,2]]
    
    num_loops <- num_loops + 1
    if(max_sim > similarity_threshold) {
      
      any_merged <- TRUE
      
      cat('Merging ',
          clust_1,
          ', ',
          clust_2,
          ', cosine similarity: ',
          max_sim,
          '\n')
      
      cluster_df <-
        cluster_df %>% 
        mutate(cluster = ifelse(cluster == clust_1, clust_2, cluster))
      
    } else {
      cat(
        'Not merging',
        paste(clusters, collapse=','),
        ' greatest similarity between',
        clust_1,
        ', ',
        clust_2,
        ', cosine similarity: ',
        max_sim,
        '\n'
      )
    }
  }
  
  cluster_df %>% 
    mutate(num_loops = num_loops) %>% 
    return()
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

link_homology <- function(cluster_df, homology_df, components_df) {
  
  multi_cluster_bubbles <-
    homology_df %>%
    left_join(cluster_df, by='unitig') %>% 
    group_by(bubble) %>% 
    filter(all_are_unique(cluster)) %>% 
    ungroup() %>% 
    pull_distinct(bubble)
  
  for(bub in multi_cluster_bubbles) {
    
    # This is a fairly aggressive strategy, that puts a lot of faith in the
    # homology detection by mashmap.
    
    bubble_clusters_df <-
      homology_df %>% 
      left_join(cluster_df, by = 'unitig') %>% 
      filter(bubble %in% bub) 
    
    bubble_clusters <-
      bubble_clusters_df %>% 
      pull(cluster) 
    
    if (all(is.na(bubble_clusters))) {
      # This step shouldn't ever be reached, as both NA would not trigger a multi
      # cluster bubble?
      next
    }
    
    if(sum(is.na(bubble_clusters)) == 1) {
      # TODO expand to check for entire NA components? Or just run
      # propagate_one_cluster_components again?
      na_unitig <- 
        bubble_clusters_df %>% 
        filter(is.na(cluster)) %>% 
        pull(unitig)
      
      partner_unitig <- 
        bubble_clusters_df %>% 
        filter(!is.na(cluster)) %>% 
        pull(unitig)
      
      
      target_cluster <- 
        bubble_clusters_df %>% 
        filter(!is.na(cluster)) %>% 
        pull(cluster)
      
      cat(
        'joining NA unitig: ',
        na_unitig,
        ' to cluster:',
        target_cluster,
        'via partner: ',
        partner_unitig, 
        '\n'
      )
      cluster_df <-
        cluster_df %>% 
        mutate(cluster = ifelse(unitig == na_unitig, target_cluster, cluster))
    } else {
      
      next 
      
      # This rule is currently really brittle, and not-uncommonly to sometimes
      # entirely separate clusters to be merged, because two small telomeric
      # unitigs that appear homologous are clustered separately. Currently
      # refraining from merging clusters, relying on cosine similarity instead.
      
      # cat(
      #   'joining clusters: ',
      #   bubble_clusters[1],
      #   ' and ',
      #   bubble_clusters[2],
      #   ' via unitigs: ',
      #   bubble_clusters_df$unitig,
      #   '\n'
      # )
      # 
      # cluster_df <-
      #   cluster_df %>% 
      #   mutate(cluster = ifelse(cluster == bubble_clusters[1], bubble_clusters[2], cluster))
      
      
    }
    
  }
  
  return(cluster_df)
}


# Phasing -----------------------------------------------------------------
orient_counts <- function(counts_df, strand_orientation_clusters_df){
  out <-
    counts_df %>%
    bind_with_inverted_unitigs() %>% 
    semi_join(
      filter(strand_orientation_clusters_df, strand_cluster == 1),
      by = c('unitig', 'unitig_dir')
    )
  
  out <- 
    out %>% 
    select(-unitig_dir)
  
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
  
  
  n_bubbles <-
    dim(phaser_array)[2]
  
  # This special case can sometimes lead to a weird behavior where no swaps
  # occur, even when it should? See verkko 1.3.1 NA24385 X/Y Chromosome. It
  # appears that the NAs and bubble arms were just so distributed such that the
  # concensus score was the same before and after swapping for every row, and
  # therefore no swaps were performed, even though there was an obviously better
  # solution that involved swaps. Now I just force an arbitrary ordering in
  # these cases.
  if(n_bubbles == 1) {
    min_value = min(phaser_array, na.rm = TRUE)
    max_value = max(phaser_array, na.rm = TRUE)
    
    min_swap <- phaser_array[,,'watson'] == min_value
    max_swap <- phaser_array[,,'crick'] == max_value
    
    min_swap[is.na(min_swap)] <- TRUE
    max_swap[is.na(max_swap)] <- TRUE
    
    lib_swapped <- min_swap & max_swap
    # for if there is only one library and one bubble:
    if(is.null(names(lib_swapped))) {
      names(lib_swapped) <- dimnames(phaser_array)$lib
    }
    return(lib_swapped)
  }
  
  
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


get_chrom_cluster_data_planes <- function(counts_df, cluster_df, unitig_lengths_df, supervision = c('PC1' , 'inverse'), min_n=5) {
  
  supervision <- match.arg(supervision)
  
  cluster_counts <-
    counts_df %>% 
    bind_with_inverted_unitigs() %>% # self-supervision
    left_join(cluster_df, by='unitig') %>% 
    split(.$cluster)
  
  wfracs <-
    cluster_counts %>% 
    map(function(x) with(x, make_wc_matrix(w,c,lib,unitig_dir,min_n=min_n)))  %>% 
    # filling with 0s doesn't seem to affect first PC too much, compared to
    # probabilistic or Bayesian PCA (from pcaMethods bioconductor package)
    map(function(x) {
      x[is.na(x)] <- 0 
      return(x)
    }) 
  
  prcomps <-
    wfracs %>% 
    map(prcomp)

  model_input <-
    map(prcomps, function(x) {
      x <-
        x %>%
        get_prcomp_plotdata(., 'unitig_dir') %>%
        mutate(unitig = gsub('_inverted', '', unitig_dir)) %>%
        left_join(unitig_lengths_df, by='unitig')
      
      if(supervision == 'PC1') {
        x <-
          x %>%
          mutate(y = sign(PC1)==1) 
      } else if(supervision == 'inverse') {
        x <-
          x %>%
          mutate(y = grepl('inverted', unitig_dir)) 
      } else {
        stop('nothing yet')
      }
      
      return(x)

    })
  
  # TODO weighted PCA instead of weighted GLM after PCA?
  glms <-
    map(model_input, function(x){
      w <- x$length
      
      x <-
        x %>%
        select(-length,-unitig_dir, -unitig)
      
      glm(y ~ -1 + ., family=binomial(), weights = w, data=x)
    })
  
  plane_vectors <-
    map(glms, function(x) {
      list(
        decision_boundary = -coef(x)[1] / coef(x)[2],
        decision_boundary_perpendicular = coef(x)[2] / coef(x)[1]
      ) %>% map(function(slp) {
        vec <- c(1, slp)
        unit_vector <- vec / sqrt(sum(vec * vec))
        return(unit_vector)
      })
    })

  
  lib_weights <-
    map2(prcomps, plane_vectors, function(x, vecs) {
      projections <-
        map(vecs, function(vec) {
          x$rotation[, 1:2] %*% vec
        })
      
      # What is this check for/ how does it work?
      stopifnot(all(rownames(projections[[1]]) == rownames(projections[[2]])))
      
      out <-
        tibble(
          lib = rownames(projections[[1]]),
          b_weight = projections$decision_boundary[, 1],
          p_weight = projections$decision_boundary_perpendicular[, 1]
        )
      
      return(out)
    })
  
  lib_weights_df <-
    bind_rows(lib_weights, .id='cluster')
  
  model_input <-
    bind_rows(model_input, .id='cluster') %>% 
    select(cluster, unitig, unitig_dir, y, everything())
  
  return(list(weights = lib_weights_df, model_input = model_input, plane_vectors = plane_vectors))
  
}
# Misc --------------------------------------------------------------------


