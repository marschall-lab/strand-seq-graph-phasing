args <- commandArgs(trailingOnly = FALSE)

# Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--components-gfa',
    '--ss-libs',
    '--alignment-link',
    ## Params
    '--threads',
    ## Output
    '--output-prefix',
    '--log',

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
    stopifnot(length(values)>1)
  }
  
  return(values)
}


get_script_dir <- function() {
  needle <- '--file='
  script_path <- gsub(needle, '',  args[grepl(needle, args)])
  return(dirname(script_path))
}



# Log ---------------------------------------------------------------------
## Log
log_path <- get_values('--log')
log <- file(log_path, open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

print(args)
# Library -----------------------------------------------------------------

library(magrittr)
library(dplyr)
library(purrr)
library(furrr)
library(Biostrings)

# Parsing -----------------------------------------------------------------

stopifnot(all(expected_args %in% args))

## Input
components_gfa <- get_values("--components-gfa", singular=FALSE)
ss_fasta_files <- get_values("--ss-fasta", singular=FALSE)
link <- get_values("--alignment-link")

n_threads <- as.numeric(get_values("--threads"))
## Output
output_prefix <- get_values('--output-prefix')



# Functions ---------------------------------------------------------------

read_rnames_from_gfa <- function(x) {
  grep('^S',
       readr::read_lines(x),
       value = TRUE) %>%
    strsplit("\t") %>%
    map_chr(2) # grab the second field in every line
}


basename_no_ext <- function(x) {
  tools::file_path_sans_ext(basename(x))
}
# Parallel Computing ------------------------------------------------------

if(n_threads > 1){
  if(supportsMulticore()) {
    plan(multicore,    workers=n_threads)
  } else {
    plan(multisession, workers=n_threads)
  }
}

# Main --------------------------------------------------------------------

### Associate Rnames to components
link <- readr::read_tsv(link)

component_ids <- 
  basename_no_ext(components_gfa) %>% 
  set_names()

components_gfa %<>% 
  set_names(component_ids)

component_rnames <-
  future_map(components_gfa, read_rnames_from_gfa) %>% 
  map_dfr(~tibble(rname=.x), .id='component')

link %<>% inner_join(component_rnames, by='rname')

### Associate qnames to libs

ss_fasta_files %<>% set_names(basename_no_ext(ss_fasta_files))

# DNAStringSet is basically a slightly fancy list
all_libs <- 
  future_map(ss_fasta_files, readDNAStringSet) 

# split link by ss lib
lib_links <-
  map_dfr(all_libs, ~tibble(qname=names(.x)), .id='lib') %>%  # qnames_by_lib
  right_join(link, by='qname') %>% 
  split(.$lib) %>% 
  map(~select(.x, -lib))

# split lib by component
all_libs <- 
  future_map2(all_libs, lib_links, function(lib, lib_link){
    # select and order
    lib <- lib[lib_link$qname]
    
    return(split(lib, lib_link$component))
    }) 

# Gather ss_reads by component
out <- map(component_ids, function(comp) {
  map(all_libs, ~.x[[comp]]) %>% 
    purrr::reduce(c)
})

# Handle empty components?
# export
iwalk(out, function(x, nm) {
  out_name <- paste0(output_prefix, nm, '.fasta')
  if(is.null(x)) {
    file.create(out_name)
  } else { 
    writeXStringSet(x, out_name)
  }

})


# Close Parallel Workers --------------------------------------------------

plan(sequential)