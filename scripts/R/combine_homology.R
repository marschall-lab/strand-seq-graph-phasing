
# Functions ---------------------------------------------------------------
tidy_bubble_ <- function(l) {
  unitigs <- flatten_chr(l$inside)
  tibble(bubble_type=l$type,
         bubble = paste0('bubble_', l$id),
         bubble_arm = c('arm_1', 'arm_2'),
         unitig = c(min(unitigs), max(unitigs)))
}

tidy_bubble_chain_ <- function(l) {
  out <- map_dfr(l$bubbles, tidy_bubble_)
  out$chain_id <- l$chain_id
  
  return(out)
}

# Command Line ------------------------------------------------------------
## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--mashmap',
    '--bubblegun',
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
library(purrr)

source(file.path(get_script_dir(), "module_utils/utils.R"))
# Parsing -----------------------------------------------------------------

mashmap <- get_values("--mashmap", singular=TRUE)
bubblegun <- get_values('--bubblegun', singular=TRUE)
output <- get_values('--output')

# Mashmap -----------------------------------------------------------------

# mashmap homolgy can be very useful, as it can be used to connect disconnected
# unitigs to the rest of their cluster. However, mashmap also has a bad habit of
# connecting unrelated degenerate unitigs across chromosomes. For this reason,
# bubblegun is prioritized in homology.

# TODO experiment with what will happen if a unitig can have more than 1 partner
# bubble unitig. this may happen in the case of, eg a bubble broken in the middle

# Experimenting has occurred. it also leads to homology being mapped across
# chromosomes, especially with the more degenerate acrocentric chromosomes, but
# generally with the degenerate chromosome regions. Unfortunately, this means
# some genuine homology where a single unitig makes multiple appearence, eg like
# broken bubbles, will be thrown out to prevent consequences of such errors.
# There is likely some potential at rescuing at least some information from
# these instances, but that is currently only done by picking the one match with
# the highest ID. There is probably more than could be done by, eg, using
# component information.

# TODO check if the above is still true or just due to the filtering bug when
# mashmap was updated to 3.0 (paf format).

# mashmap <- 'homology/NA18989/NA18989_mashmap.paf'

mashmap_df <-
  pafr::read_paf(mashmap, tibble=TRUE, include_tags = TRUE)

# Filtering criteria from verkko HiC pipeline
mashmap_df <-
  mashmap_df %>% 
  filter(id > 99) %>% 
  filter(qend - qstart >= 5e5) %>% 
  filter(qname != tname)

# If a unitig is in multiple bubbles, pick the highest id?
mashmap_df <-
  mashmap_df %>% 
  arrange(id, qend-qstart) %>% 
  group_by(qname) %>% 
  slice_tail(n=1) %>% 
  ungroup()

# filter out identical pairs, just with unitig_1 and unitig_2 switched
mashmap_df <-
  mashmap_df %>%
  filter(!is_duplicate_pair(qname, tname))

# sort bubble arms
mashmap_df <-
  mashmap_df %>% 
  transmute(unitig_1 = qname, unitig_2 = tname) %>% 
  mutate(unitig_1 = map2_chr(unitig_1, unitig_2, min),
         unitig_2 = map2_chr(unitig_1, unitig_2, max))

# tidy
mashmap_df <-
  mashmap_df %>% 
  mutate(bubble = paste0('bubble_', 1:n())) %>% 
  tidyr::pivot_longer(cols = c('unitig_1', 'unitig_2'), names_to = 'bubble_arm', values_to ='unitig') %>% 
  mutate(bubble_arm = stringr::str_replace(bubble_arm, '^unitig_', 'arm_'))


# filter out bubbles that contain unitigs in multiple bubbles

# only unitigs that appear once
mashmap_df <-
  mashmap_df %>%
  group_by(unitig) %>%
  filter(n() == 1) %>%
  ungroup()

#only bubbles with both arms remaining
mashmap_df <-
  mashmap_df %>% 
  group_by(bubble) %>% 
  filter(n() == 2) %>% 
  ungroup()


# Bubblegun ---------------------------------------------------------------

# bubblegun <- 'homology/NA18989/NA18989_exploded_simplified_bubblegun.json'

bubblegun_chains_json <- jsonlite::read_json(bubblegun)

bubblegun_df <- map_dfr(bubblegun_chains_json, tidy_bubble_chain_)

bubblegun_df <-
  bubblegun_df %>% 
  filter(bubble_type == 'simple') %>% 
  select(bubble,bubble_arm,unitig)


# Deduplicate Bubbles -----------------------------------------------------

mashmap_df <-
  mashmap_df %>% 
  anti_join(bubblegun_df, by='unitig') %>%  
  group_by(bubble) %>% 
  filter(n() == 2) %>% 
  ungroup()

# renumber 
mashmap_df <-
  mashmap_df %>% 
  mutate(bubble = paste0('bubble_', as.integer(as.factor(bubble)))) %>% 
  arrange(bubble, bubble_arm)

# Rename Bubbles ----------------------------------------------------------

mashmap_df <-
  mashmap_df %>%
  mutate(bubble = gsub('bubble_', 'mashmap_', bubble))

bubblegun_df <-
  bubblegun_df %>%
  mutate(bubble = gsub('bubble_', 'bubblegun_', bubble))


# Export ------------------------------------------------------------------

homology_df <-
  bind_rows(mashmap_df, bubblegun_df)

# output <- 'homology/NA18989/NA18989_combined_homology.tsv'
readr::write_tsv(homology_df, output)
