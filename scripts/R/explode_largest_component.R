# Library -----------------------------------------------------------------

library(dplyr)
library(purrr)
library(igraph)

# Import ------------------------------------------------------------------

extract_links_from_gfa <- function(filepath) {
  con <- file(filepath, "r")
  out <- character()

  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }

    if(stringr::str_starts(line, 'L')){
      out <- c(out, line)
    }
  }
  
  close(con)
  return(out)
}


extract_segment_sizes_from_gfa <- function(filepath) {
  con <- file(filepath, "r")
  out <- integer()
  
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    if(stringr::str_starts(line, 'S')){
      line <- strsplit(line, '\t')[[1]]
      segment_name <- line[2]
      segment_length <- nchar(line[3])
      out <- c(out, set_names(segment_length, segment_name))
    }
  }
  
  close(con)
  return(out)
}

gfa <- 'NA24385_good_frac-50_rep-1_test/assembly.gfa'
links_df <- extract_links_from_gfa(gfa)

links_df <-
  links_df %>%
  strsplit("\t") %>%
  map(as.list) %>%
  map(set_names,
      c('RecordType', 'From', 'FromOrient', 'To', 'ToOrient', 'Overlap')) %>%
  map(as_tibble) %>%
  bind_rows()

# Segment Lengths ---------------------------------------------------------

segment_sizes <-
  extract_segment_sizes_from_gfa(gfa)

# Graph Formatting --------------------------------------------------------


graph <-
  graph_from_data_frame(
    select(links_df, From, To, everything()),
    directed = FALSE,
    vertices = tibble::enframe(segment_sizes, value =
                                 'size')
  )

# Remove Small Isolated ---------------------------------------------------

small_nodes <-
  names(segment_sizes)[segment_sizes < segment_length_threshold]

isolated_small_nodes <-
  degree(graph, v = unname(small_nodes)) %>%
  keep( ~ .x == 0) %>%
  names()

graph <-
  graph %>%
  delete_vertices(isolated_small_nodes)


# Explode Largest Component -----------------------------------------------

get_nodes_in_largest_component <- function(graph) {
  ccs <- components(graph)
  
  largest_id <-
    which.max(ccs$csize)
  
  out <-
    with(ccs, names(membership)[membership == largest_id])
  
  return(out)
}

acrocentric_nodes <- get_nodes_in_largest_component(graph)

acrocentric_small_nodes <-
  intersect(small_nodes, acrocentric_nodes)

largest_tangle_nodes <-
  get_nodes_in_largest_component(subgraph(graph, acrocentric_small_nodes))

graph <-
  graph %>%
  delete_vertices(largest_tangle_nodes)


# Export ------------------------------------------------------------------

deleted_nodes <- setdiff(names(segment_sizes), names(V(graph)))

remove_nodes_from_gfa <- function(inpath, outpath, nodes_to_remove) {
  in_gfa <- file(inpath, "r")
  out_gfa <- file(outpath, "w")
  
  
  while ( TRUE ) {
    line <- readLines(in_gfa, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    split_line <- strsplit(line, '\t')[[1]]
    if(split_line[1] == 'S'){
      if(split_line[2] %in% nodes_to_remove) {
        next
      }
    } else if(split_line[1] == 'L') {
      if(split_line[2] %in% nodes_to_remove | split_line[4] %in% nodes_to_remove) {
        next
      }
    } else {
      browser()
      stop('error lol')
    }
    
    writeLines(line, out_gfa)
  }
  
  close(in_gfa)
  close(out_gfa)
  return(TRUE)
}

remove_nodes_from_gfa(gfa, 'NA24385_good_frac-50_rep-1_test/exploded.gfa', deleted_nodes)

# Connected Components in Exploded Graph ----------------------------------

components_df <-
  with(components(graph),
       tibble(unitig = names(membership), component = membership)) %>% 
  mutate(member_largest_component = unitig %in% acrocentric_nodes)



# Scratch -----------------------------------------------------------------

# 
# 
# cluster_df <-
#   set_names(clust@.Data, clust@names) %>% 
#   map(~tibble(unitig = .x)) %>% 
#   bind_rows(.id = 'cluster')
# 
# cluster_df <- 
#   cluster_df %>% 
#   mutate(unitig = strip_range(unitig))
# 
# # components_df <-
# #   components_df %>% 
# #   left_join(cluster_df, by = c('unitig')) %>% 
# #   left_join(rlengths_df, by=c('unitig' = 'rname')) %>% 
# #   left_join(reference, by=c('unitig' = 'qname'))
# 
# 
# library(ggfortify)
# 
# 
# tmp_wf <-
#   wfrac.matrix
# 
# rownames(tmp_wf) <- strip_range(rownames(tmp_wf))
# 
# 
# info_df <-
#   counts %>% 
#   ungroup() %>% 
#   mutate(rname = strip_range(rname)) %>% 
#   left_join(reference,  by = c('rname' = 'qname')) %>% 
#   left_join(components_df, by=c('rname' = 'unitig')) %>% 
#   left_join(cluster_df, by=c('rname' = 'unitig'))
# 
# 
# tmp_ulist <-info_df %>% filter(component == 32) %>% filter(qlen >= segment_length_threshold) %>%  pull(rname) 
# # highwc <- as.character(strand.states$AWCcontigs@seqnames)
# # tmp_ulist <- setdiff(tmp_ulist, highwc)
# 
# # dbs <- dbscan::hdbscan(tmp_wf[rownames(tmp_wf) %in% tmp_ulist, ], minPts = 10)
# # 
# # dbs_df 
# 
# autoplot(prcomp(abs(tmp_wf[rownames(tmp_wf) %in% tmp_ulist, ])),
#          data = info_df %>% filter(rname %in% tmp_ulist),
#          colour = 'cluster',
#          shape='tname',
#          x=1, y=2)   
# 
# 
# 
# bandage_colors <-
#   info_df %>% 
#   select(-lib, -c, -w) %>% 
#   distinct()
#   
# palette <-
#   bandage_colors %>% 
#   distinct(cluster) %>% 
#   mutate(Color = rainbow((n()))) %>% 
#   mutate(Color = ifelse(is.na(cluster), 'black', Color))
# 
# 
# bandage_colors <-
#   bandage_colors %>% 
#   left_join(palette)
# 
# readr::write_csv(bandage_colors, 'NA24385_good_frac-50_rep-1_test/colors.csv')
# 
# 
