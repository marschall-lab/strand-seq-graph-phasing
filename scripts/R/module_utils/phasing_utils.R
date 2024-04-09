
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

aggregate_over_subsets <- function(mat, sets, contrasts, col_names=c('id_1', 'id_2', 'value'), agg_f = mean, ...) {
  # mat ~ row and column names.
  # sets ~ named list, names are cluster ids, values are unitigs
  # contrasts ~ a list, with pairs of set names to be compared to one another.
  
  stopifnot(length(col_names) == 3)
  stopifnot(all(map_int(contrasts, length) == 2))
  
  out <-
    map_dfr(contrasts, function(pair) {
      id_1 <- pair[1]
      id_2 <- pair[2]
      e_1 <- sets[[id_1]]
      e_2 <- sets[[id_2]]
      value <- agg_f(mat[e_1, e_2, drop=FALSE], ...)
      out <- tibble(id_1 = id_1, id_2=id_2, value=value)
      return(out)
    })
  
 out <-
   bind_rows(
     out,
     tibble(id_1=character(), id_2=character(), value=double()) 
   )
 
  out <-
    set_names(out, col_names)
  
  
  return(out)
}
merge_similar_clusters_on_components <- function(cosine_similarity_mat, cluster_df, components_df, similarity_threshold =0.40, agg_f=mean_abs, ...) {
  any_merged <- TRUE
  while(any_merged) {
    any_merged <- FALSE
    
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
      aggregate_over_subsets(
        cosine_similarity_mat,
        sets,
        contrasts,
        col_names = c('clust_1', 'clust_2', 'sim'),
        agg_f = agg_f,
        ...
      )
    
    
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
    if(max_sim > similarity_threshold) {
      any_merged <- TRUE
      
      # sort to ensure that sex cluster (sex > LG) is not overwritten
      clusters_to_merge <- sort(c(max_contrast$clust_1, max_contrast$clust_2))
      clust_1 <- clusters_to_merge[1]
      clust_2 <- clusters_to_merge[2]
      
      cat('Merging ', clust_1, ' into ', clust_2, ', cosine similarity: ', max_sim,
          '\n')
      
      cluster_df <-
        cluster_df %>% 
        mutate(cluster = ifelse(cluster == clust_1, clust_2, cluster))
      
    } 
  }
  
  return(cluster_df)
}


merge_similar_clusters <- function(cosine_similarity_mat, cluster_df, similarity_threshold =0.50, agg_f=mean_abs, ...) {
  any_merged <- TRUE
  while(any_merged) {
    any_merged <- FALSE
    
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
      aggregate_over_subsets(
        cosine_similarity_mat,
        sets,
        contrasts,
        col_names = c('clust_1', 'clust_2', 'sim'),
        agg_f = agg_f,
        ...
      )
    
    
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
    if(max_sim > similarity_threshold) {
      any_merged <- TRUE
      
      # sort to ensure that sex cluster (sex > LG) is not overwritten
      clusters_to_merge <- sort(c(max_contrast$clust_1, max_contrast$clust_2))
      clust_1 <- clusters_to_merge[1]
      clust_2 <- clusters_to_merge[2]
      
      cat('Merging ', clust_1, ' into ', clust_2, ', cosine similarity: ', max_sim,
          '\n')
      
      cluster_df <-
        cluster_df %>% 
        mutate(cluster = ifelse(cluster == clust_1, clust_2, cluster))
      
    } 
  }
  
  return(cluster_df)
}

increment_cluster_ <- function(cluster_df, max_contrast, new_cluster_str=NULL) {
  max_sim <- max_contrast$sim
  if(is.null(new_cluster_str)) {

    max_cluster <- max_contrast$cluster
    max_unitig <- max_contrast$unitig
    
    cat(
      'Assigning unitig:', max_unitig,
      'to cluster:', max_cluster,
      'similarity score:', max_sim,
      '\n'
    )
    
    cluster_df <-
      cluster_df %>%
      mutate(cluster = ifelse(unitig == max_unitig, max_cluster, cluster))
  } else {
    new_cluster_unitigs <- c(max_contrast$unitig_1, max_contrast$unitig_2)
    
    cat('creating new cluster with unitigs:', new_cluster_unitigs,
        'with similarity:', max_sim,
        '\n')
    
    cluster_df <-
      cluster_df %>%
      mutate(cluster = ifelse(unitig %in% new_cluster_unitigs, new_cluster_str, cluster))
  }
  
  return(cluster_df)
}


cluster_unitigs <-
  function(cosine_similarity_mat,
           cluster_df,
           cluster_unitig_similarity_threshold = 0.50,
           unitig_unitig_similarity_threshold = 0.50,
           new_cluster_id = 'LGcos',
           new_cluster_start_ix = 0,
           agg_f=mean_abs, 
           ...) {
    
    
    unassigned_ids <-
      cluster_df %>% 
      filter(is.na(cluster)) %>% 
      pull_distinct(unitig)
    
    if(length(unassigned_ids) == 0) {
      return(cluster_df)
    }
    
    sets <- split(unassigned_ids, unassigned_ids)
    
    contrasts <-
      make_contrasts(unassigned_ids, unassigned_ids)
    
    contrasts <- unique(map(contrasts, sort))
    
    clustered_similarities_uu <-
      aggregate_over_subsets(
        cosine_similarity_mat,
        sets,
        contrasts,
        col_names = c('unitig_1', 'unitig_2', 'sim'),
        agg_f = agg_f,
        ...
      )
    
    clustered_similarities_uu <-
      clustered_similarities_uu %>%  
      filter(unitig_1 != unitig_2) %>% 
      filter(!is.na(sim)) 
    
    sets <- list()
    clustered_similarities_cu <- tibble(cluster=character(), unitig=character(), sim=double())
    
    cur_clusters <- 
      cluster_df %>% 
      filter(!is.na(cluster)) %>% 
      pull_distinct(cluster)
    
    new_cluster_ix <- new_cluster_start_ix
    any_assigned <- TRUE
    while(any_assigned) {
      any_assigned <- FALSE

      unassigned_ids <-
        cluster_df %>% 
        filter(is.na(cluster)) %>% 
        pull_distinct(unitig)
      
      if(length(unassigned_ids) == 0) {
        break
      }
      
      # clustered_similarities_uu <-
      #   clustered_similarities_uu %>% 
      #   filter(unitig_1 %in% unassigned_ids | unitig_2 %in% unassigned_ids)
      
      new_sets <-
        cluster_df %>% 
        filter(cluster %in% cur_clusters | unitig %in% unassigned_ids) %>% 
        mutate(cluster = coalesce(cluster, unitig)) %>% 
        with(split(unitig, cluster))
      
      sets <- c(new_sets, sets)
      
      contrasts <-
        make_contrasts(cur_clusters, unassigned_ids)
      
      new_clustered_similarities_cu <-
        aggregate_over_subsets(
          cosine_similarity_mat,
          sets,
          contrasts,
          col_names = c('cluster', 'unitig', 'sim'),
          agg_f = agg_f,
          ...
        )
      # browser()
      clustered_similarities_cu <- 
        bind_rows(new_clustered_similarities_cu, clustered_similarities_cu)

      max_contrast <-
        clustered_similarities_cu %>% 
        filter(!is.na(sim)) %>% 
        slice_max(sim, n = 1, with_ties = FALSE)
      
      max_sim <- max_contrast$sim
      if((length(max_sim) > 0) && (max_sim > cluster_unitig_similarity_threshold)) {
        any_assigned <- TRUE
        cluster_df <- increment_cluster_(cluster_df, max_contrast)
        cur_clusters <- max_contrast$cluster
        
        clustered_similarities_cu <- 
          clustered_similarities_cu %>% 
          filter(cluster != cur_clusters) %>% 
          filter(unitig != max_contrast$unitig)
        
        clustered_similarities_uu <- 
          clustered_similarities_uu %>% 
          filter(unitig_1 != max_contrast$unitig) %>% 
          filter(unitig_2 != max_contrast$unitig) 
        
        next
      } 
      
      cat(
        'No unitig-cluster similarity greater than threshold:', cluster_unitig_similarity_threshold, '\n'
      )

      max_contrast <-
        clustered_similarities_uu %>% 
        filter(unitig_1 %in% unassigned_ids) %>% 
        filter(unitig_2 %in% unassigned_ids) %>% 
        slice_max(sim, n = 1, with_ties = FALSE)

      max_sim <- max_contrast$sim
      if(length(max_sim > 0) && max_sim > unitig_unitig_similarity_threshold) {
        any_assigned <- TRUE
        new_cluster_ix <- new_cluster_ix + 1
        cur_clusters <- paste0(new_cluster_id, new_cluster_ix)
        cluster_df <- increment_cluster_(cluster_df, max_contrast, new_cluster_str = cur_clusters)
        
        clustered_similarities_cu <- 
          clustered_similarities_cu %>% 
          filter(unitig != max_contrast$unitig_1) %>% 
          filter(unitig != max_contrast$unitig_2) 
        
        clustered_similarities_uu <- 
          clustered_similarities_uu %>% 
          filter(unitig_1 != max_contrast$unitig_1) %>% 
          filter(unitig_2 != max_contrast$unitig_1) %>% 
          filter(unitig_1 != max_contrast$unitig_2) %>% 
          filter(unitig_2 != max_contrast$unitig_2)
      } 
      
    }
    return(cluster_df)
  }

pairwise_complete_hclust_n <- function(sim_mat, n = 2, agg_f=mean, ...) {
  stopifnot(setequal(rownames(sim_mat), colnames(sim_mat)))
  cluster_df <- 
    tibble(unitig = rownames(sim_mat), cluster = rownames(sim_mat))
  
  n_clusters <- n_distinct(cluster_df$cluster)
  
  cur_clusters <-  
    cluster_df %>% 
    filter(!is.na(cluster)) %>% 
    pull_distinct(cluster)
  
  sets <- list()
  cluster_similarities <- tibble(clust_1=character(), clust_2=character(), sim=double())
  while(n_clusters > n) {
    
    new_sets <-
      cluster_df %>% 
      filter(!is.na(cluster)) %>% 
      filter(cluster %in% cur_clusters) %>% 
      with(split(unitig, cluster))
    
    sets <- c(new_sets, sets)
    
    all_clusters <- 
      cluster_df %>%
      filter(!is.na(cluster)) %>%
      pull_distinct(cluster)

    contrasts <-
      make_contrasts(cur_clusters, all_clusters) %>% 
      map(sort) %>% 
      unique()
    
    new_cluster_similarities <-
      aggregate_over_subsets(
        sim_mat,
        sets,
        contrasts,
        col_names = c('clust_1', 'clust_2', 'sim'),
        agg_f = agg_f,
        ...
      )
    
    new_cluster_similarities <- 
      new_cluster_similarities %>% 
      filter(clust_1 != clust_2) 
    
    cluster_similarities <-
      bind_rows(cluster_similarities, new_cluster_similarities)
    
    
    if(all(is.na(cluster_similarities$sim))) {
      cat('No valid similarity scores\n')
      break
    }
    
    max_contrast <-
      cluster_similarities %>% 
      slice_max(sim, n = 1, with_ties = FALSE)
    
    max_sim <- max_contrast$sim
    
    # sort to ensure that sex cluster (sex > LG) is not overwritten
    clusters_to_merge <- sort(c(max_contrast$clust_1, max_contrast$clust_2))
    clust_1 <- clusters_to_merge[1]
    clust_2 <- clusters_to_merge[2]
    
    cat('Merging ', clust_1, ' into ', clust_2, ', cosine similarity: ', max_sim,
        '\n')
    
    cluster_df <-
      cluster_df %>% 
      mutate(cluster = ifelse(cluster == clust_1, clust_2, cluster))
    
    sets[clust_1] <- NULL
    sets[clust_2] <- NULL
    
    cluster_similarities <-
      cluster_similarities %>% 
      filter(clust_1 != max_contrast$clust_1) %>% 
      filter(clust_1 != max_contrast$clust_2) %>%
      filter(clust_2 != max_contrast$clust_1) %>%
      filter(clust_2 != max_contrast$clust_2) 
    
    cur_clusters <- clust_2
    n_clusters <- n_distinct(cluster_df$cluster)
  }
  
  n_clusters <- n_distinct(cluster_df$cluster)
  stopifnot(n_clusters == n)
  
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


# Misc --------------------------------------------------------------------

mean_abs <- function(x, ...) {
  mean(abs(x), ...)
}

project_onto <- function(x, y) {
  stopifnot(length(y) == length(x))
  complete_ix <- !is.na(x) & !is.na(y)
  x[!complete_ix] <- NA
  y[!complete_ix] <- NA
  y_norm <- y/sqrt(sum(y^2, na.rm=TRUE))
  x_onto_y <-  sum(x * y_norm, na.rm=TRUE) * y_norm
  
  return(x_onto_y)
}

gram_schmidt_orthognalize <- function(x, y) {
  stopifnot(length(y) == length(x))
  x_onto_y <- project_onto(x, y)
  x_ortho <- x - x_onto_y
  return(x_ortho)
}

project_through <- function(x, y) {
  return(x - 2 * gram_schmidt_orthognalize(x, y))
}


# contiBAIT ---------------------------------------------------------------

preprocessStrandTable <-
  function (strandTable, strandTableThreshold = 0.8, filterThreshold = 0.8, 
            orderMethod = "libsAndConc", lowQualThreshold = 0.9, verbose = TRUE, 
            minLib = 10) 
  {
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
