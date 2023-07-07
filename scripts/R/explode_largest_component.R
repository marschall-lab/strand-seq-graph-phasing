
# Functions ---------------------------------------------------------------

get_nodes_in_largest_component <- function(graph) {
  ccs <- components(graph)
  
  largest_id <-
    which.max(ccs$csize)
  
  out <-
    with(ccs, names(membership)[membership == largest_id])
  
  return(out)
}

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
    } else if(!(split_line[1] %in% c('S', 'L'))){
      next
    }
    
    writeLines(line, out_gfa)
  }
  
  close(in_gfa)
  close(out_gfa)
  return(TRUE)
}

# Args --------------------------------------------------------------------


args <- commandArgs(trailingOnly = FALSE)

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--gfa',
    ## Output
    '--output-gfa',
    '--output-ccs',
    ## params
    '--segment-length-threshold'
    
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
library(igraph)

source(file.path(get_script_dir(), "module_utils/utils.R"))
# Parsing -----------------------------------------------------------------

gfa <- get_values("--gfa", singular=TRUE)
output_gfa <- get_values("--output-gfa", singular=TRUE) 
output_components <- get_values("--output-ccs", singular=TRUE) 
segment_length_threshold <- as.integer(get_values("--segment-length-threshold", singular=TRUE))
# Import ------------------------------------------------------------------
# gfa <- 'hifiasm_draft/v0.19.5/HG00733/HG00733.hifiasm.bp.p_utg.noseq.gfa'
# cat('Importing gfa:', gfa, '\n')


links_df <-
  read_links_from_gfa(gfa) %>%
  strsplit("\t") %>%
  map(as.list) %>%
  map(function(x) x[1:6]) %>% # remove extra link information, found in hifiasm graphs
  map(set_names,
      c('RecordType', 'From', 'FromOrient', 'To', 'ToOrient', 'Overlap')) %>%
  map(as_tibble) %>%
  bind_rows()

# Segment Lengths ---------------------------------------------------------

segment_sizes <-
  read_segment_sizes_from_gfa(gfa)

# Graph Formatting --------------------------------------------------------


vertex_information <-
  tibble(unitig = names(segment_sizes), size = segment_sizes)

graph <-
  graph_from_data_frame(
    select(links_df, From, To, everything()),
    directed = FALSE,
    vertices = vertex_information
  )


# Remove Small Isolated ---------------------------------------------------

small_nodes <-
  names(segment_sizes)[segment_sizes < segment_length_threshold]


isolated_nodes <-
  degree(graph) %>%
  keep( ~ .x == 0) %>%
  names()
#v = unname(small_nodes)

isolated_small_nodes <-
  intersect(small_nodes, isolated_nodes)

graph <-
  graph %>%
  delete_vertices(isolated_small_nodes)


# Remove Largest Component -----------------------------------------------


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

remove_nodes_from_gfa(gfa, output_gfa, deleted_nodes)

# Connected Components in Exploded Graph ----------------------------------

components_df <-
  with(components(graph),
       tibble(unitig = names(membership), component = membership)) %>% 
  mutate(member_largest_component = unitig %in% acrocentric_nodes)

readr::write_tsv(components_df, output_components)
