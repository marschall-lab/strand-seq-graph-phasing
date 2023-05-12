
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
      
      cat(
        'joining clusters: ',
        bubble_clusters[1],
        ' and ',
        bubble_clusters[2],
        ' via unitigs: ',
        bubble_clusters_df$unitig,
        '\n'
      )
      
      cluster_df <-
        cluster_df %>% 
        mutate(cluster = ifelse(cluster == bubble_clusters[1], bubble_clusters[2], cluster))
    }
    
  }
  
  return(cluster_df)
}


# Phasing -----------------------------------------------------------------


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
# Misc --------------------------------------------------------------------

