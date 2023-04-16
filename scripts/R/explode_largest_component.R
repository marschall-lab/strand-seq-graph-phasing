# Library -----------------------------------------------------------------

library(dplyr)
library(purrr)
library(igraph)

# Import ------------------------------------------------------------------

# TODO make this import process line by line.
gfa <-
  readLines('NA24385_good_frac-50_rep-1_test/assembly.gfa')

segments <-
  gfa %>%
  keep(stringr::str_starts, 'S') %>%
  strsplit("\t")

## Remove sequence to shrink graph
# segments <-
#   map(segments, function(x) {
#     x[3] <- '*'
#     return(x)
#   })

links <-
  gfa %>%
  keep(stringr::str_starts, 'L') %>%
  strsplit("\t")

links_df <-
  links %>%
  map(as.list) %>%
  map(set_names,
      c('RecordType', 'From', 'FromOrient', 'To', 'ToOrient', 'Overlap')) %>%
  map(as_tibble) %>%
  bind_rows()
#
#   readr::read_tsv(I(links),
#                   col_names = c('RecordType', 'From', 'FromOrient', 'To', 'ToOrient', 'Overlap'))
#
# segments <-
#   readr::read_tsv(I(segments))


# Segment Lengths ---------------------------------------------------------

segment_names <-
  segments %>%
  map_chr(2)

segment_lengths <-
  segments %>%
  map_chr(3) %>%
  map_int(nchar)

vertex_information <-
  tibble(unitig = segment_names, size = segment_lengths)


# Graph Formatting --------------------------------------------------------

# nodes <-
#   with(links_df, unique(c(From, To))) %>%
#   sort()
#
# links_df <-
#   links_df %>%
#   mutate(From = paste(From, ifelse(FromOrient == '+', 'S', 'E'), sep='_'),
#          To = paste(To, ifelse(ToOrient == '+', 'S', 'E'), sep='_')) %>%
#   select(From, To)
#
# end_to_end_links_df <-
#   tibble(From = paste(nodes, 'S', sep='_'),
#          To = paste(nodes, 'E', sep='_'))
#
# links_df <-
#   links_df %>%
#   bind_rows(end_to_end_links_df)
#
graph <-
  graph_from_data_frame(select(links_df, From, To, everything()),
                        directed = FALSE,
                        vertices = vertex_information)

# Remove Small Isolated ---------------------------------------------------

small_nodes <-
  vertex_information %>%
  filter(size < segment_length_threshold) %>%
  pull(unitig)

isolated_small_nodes <-
  degree(graph, v = small_nodes) %>%
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

sg <- subgraph(graph, acrocentric_small_nodes)

largest_tangle <-
  get_nodes_in_largest_component(sg)

graph <-
  graph %>%
  delete_vertices(largest_tangle)


# Export ------------------------------------------------------------------

deleted_nodes <- setdiff(segment_names, names(V(graph)))

segments <-
  segments %>%
  discard( ~ .x[2] %in% deleted_nodes) %>%
  map_chr(paste, collapse = '\t')

links <-
  links %>%
  discard( ~ (.x[[2]] %in% deleted_nodes) | (.x[[4]] %in% deleted_nodes)) %>%
  map_chr(paste, collapse = '\t')

readr::write_lines(c(segments, links), file = 'NA24385_good_frac-50_rep-1_test/exploded.gfa')



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
