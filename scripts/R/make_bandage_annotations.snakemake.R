
# Command Line ------------------------------------------------------------


## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    '--haplotype-marker-counts',
    '--ref-alignment',
    '--output'
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
library(pafr)

# Parsing -----------------------------------------------------------------

haplotype_marker_counts <- get_values('--haplotype-marker-counts', singular=TRUE)

ref_aln <- get_values('--ref-alignment', singular=TRUE)

out_path <- get_values('--output', singular=TRUE)
# Import Reference Alignments ---------------------------------------------

ref_aln_df <-
   pafr::read_paf(ref_aln, tibble = TRUE, include_tags = FALSE)

# assign qnames to tnames based on majority alignment? Based on first listing?
ref_aln_df <-
  ref_aln_df %>%
  group_by(qname) %>%
  slice_head(n=1) %>%
  dplyr::rename(unitig = qname) 

# Import Marker Counts ----------------------------------------------------

haplotype_marker_df <- 
  readr::read_csv(haplotype_marker_counts) 

# Reference Alignment & Marker Counts -------------------------------------------------

# For looking at w/ bandage graphs

out <- 
  left_join(haplotype_marker_df, ref_aln_df, by = 'unitig') %>% 
  select(unitig, everything())

out %>% 
  readr::write_csv(out_path)
