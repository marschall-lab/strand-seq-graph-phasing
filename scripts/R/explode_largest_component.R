
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

# Library -----------------------------------------------------------------
library(argparse)
library(dplyr)
library(purrr)
library(igraph)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  needle <- '--file='
  script_path <- gsub(needle, '',  args[grepl(needle, args)])
  return(dirname(script_path))
}

source(file.path(get_script_dir(), "module_utils/utils.R"))

# Parsing -----------------------------------------------------------------

print(commandArgs(trailingOnly = FALSE))

parser <- ArgumentParser(description = 'Remove Tangle in Largest Component')
parser$add_argument(
  '--gfa',
  help = 'input gfa',
  required = TRUE,
  nargs = 1
  )
parser$add_argument(
  '--output-gfa', 
  required = TRUE, 
  nargs = 1
  )
parser$add_argument(
  '--output-ccs', 
  required = TRUE, 
  nargs = 1
  )
parser$add_argument(
  '--segment-length-threshold',
  help = 'Length below which a node is considered "short" and to be included in tangle detection and removal',
  required = TRUE,
  nargs = 1
)

args <- parser$parse_args()
print(args)


# Import ------------------------------------------------------------------
# gfa <- 'hifiasm_draft/v0.19.5/HG00733/HG00733.hifiasm.bp.p_utg.noseq.gfa'
# cat('Importing gfa:', gfa, '\n')


links_df <-
  read_links_from_gfa(args$gfa) %>%
  strsplit("\t") %>%
  map(as.list) %>%
  map(function(x) x[1:6]) %>% # remove extra link information, found in hifiasm graphs
  map(set_names,
      c('RecordType', 'From', 'FromOrient', 'To', 'ToOrient', 'Overlap')) %>%
  map(as_tibble) %>%
  bind_rows()

# Segment Lengths ---------------------------------------------------------

segment_sizes <-
  read_segment_sizes_from_gfa(args$gfa)

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
  names(segment_sizes)[segment_sizes < args$segment_length_threshold]


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

remove_nodes_from_gfa(args$gfa, args$output_gfa, deleted_nodes)

# Connected Components in Exploded Graph ----------------------------------

components_df <-
  with(components(graph),
       tibble(unitig = names(membership), component = membership)) %>% 
  mutate(member_largest_component = unitig %in% acrocentric_nodes)

readr::write_tsv(components_df, args$output_ccs)
