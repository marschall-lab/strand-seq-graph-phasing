
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

tidy_multibubble_ <- function(bubble_name, arm_1_unitigs, arm_2_unitigs) {
  arm_1_unitigs <- strsplit(arm_1_unitigs, '_')[[1]]
  arm_2_unitigs <- strsplit(arm_2_unitigs, '_')[[1]]
  
  # for now only one bubble from the multibubble is kept, as I don't have a way
  # of guaranteeing that unitigs from the same merged bubble arm won't end up in
  # different haplotypes. If more than 1 bubble is desired, remove the slice_head
  
  # Sort by size
  arm_1_unitigs <-
    unitig_lengths_df %>%
    filter(unitig %in% arm_1_unitigs) %>%
    arrange(desc(length)) %>%
    slice_head(n=1) %>% 
    pull(unitig)
  
  arm_2_unitigs <-
    unitig_lengths_df %>%
    filter(unitig %in% arm_2_unitigs) %>%
    arrange(desc(length)) %>%
    slice_head(n=1) %>% 
    pull(unitig)
  
  n_pairs <- min(length(arm_1_unitigs), length(arm_2_unitigs))
  
  out_df <-
    tibble(bubble = character(),
           bubble_arm = character(),
           unitig = character())
  ix <- 0
  for (i in seq_len(n_pairs)) {
    ix <- ix + 1
    bub_df <- tibble(
      bubble = paste0(bubble_name, '_merged_', ix),
      bubble_arm = c('arm_1', 'arm_2'),
      unitig = c(arm_1_unitigs[i], arm_2_unitigs[i])
    )
    
    out_df <- bind_rows(out_df, bub_df)
  }
  
  return(out_df)
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
    '--gfa',
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
gfa <- get_values('--gfa', singular=TRUE)
output <- get_values('--output', singular=TRUE)

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

# mashmap <- 'homology/NA19317/NA19317_mashmap.paf'

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
  mutate(unitig_1 = map2_chr(qname, tname, min),
         unitig_2 = map2_chr(qname, tname, max)) %>% 
  select(unitig_1, unitig_2)

# tidy
mashmap_df <-
  mashmap_df %>% 
  mutate(bubble = paste0('bubble_', 1:n())) %>% 
  tidyr::pivot_longer(cols = c('unitig_1', 'unitig_2'), names_to = 'bubble_arm', values_to ='unitig') %>% 
  mutate(bubble_arm = stringr::str_replace(bubble_arm, '^unitig_', 'arm_'))

# filter out bubbles that contain unitigs in multiple bubbles

# only unitigs that appear once. Is this step still needed given earlier
# filtering?
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

## Read Unitig Lengths -----------------------------------------------------

# gfa <- '../wd/gfa/gfa/NA19317_exploded.gfa'
unitig_lengths_df <-
  tibble::enframe(
    read_segment_sizes_from_gfa(gfa),
    name = 'unitig',
    value = 'length'
  )


## Read Bubbles ------------------------------------------------------------
# bubblegun <- 'homology/NA19317/NA19317_exploded_simplified_bubblegun.json'

bubblegun_chains_json <- jsonlite::read_json(bubblegun)

bubblegun_df <- map_dfr(bubblegun_chains_json, tidy_bubble_chain_)

bubblegun_df <-
  bubblegun_df %>% 
  filter(bubble_type == 'simple') %>% 
  select(bubble,bubble_arm,unitig)


## Merged Bubble Handling --------------------------------------------------

# The graph simplification process now merges linear paths! That should recover
# some more homology! Need to decide how to handle many-to-many unitig
# homologies as (for now) each unitig is only appearing in one bubble.
# Currently, pair biggest unitig with biggest possible partner.

# FIXME this is an issue where locally, merged bubbles work, but on the cluster,
# they appear to be removed, eg utig4-835 sample HG04036

bubbles_with_merged_unitigs <- 
  bubblegun_df %>% 
  filter(grepl('_', unitig)) %>% 
  pull_distinct(bubble)

merged_bubble_df <- 
  bubblegun_df %>% 
  filter(bubble %in% bubbles_with_merged_unitigs) %>% 
  tidyr::pivot_wider(names_from = bubble_arm, values_from = unitig) %>% 
  pmap_dfr(tidy_multibubble_)

bubblegun_df <-
  bubblegun_df %>% 
  filter(!(bubble %in% bubbles_with_merged_unitigs)) %>% 
  bind_rows(merged_bubble_df) 
# Deduplicate Bubbles -----------------------------------------------------

mashmap_wide <-
  tidyr::pivot_wider(mashmap_df, names_from = 'bubble_arm', values_from = 'unitig') %>% 
  arrange(arm_1) %>% 
  select(-bubble)
  
bubblegun_wide <-
  tidyr::pivot_wider(bubblegun_df, names_from = 'bubble_arm', values_from = 'unitig')%>% 
  arrange(arm_1) %>% 
  select(-bubble)



shared_bubbles <-
  inner_join(mashmap_wide, bubblegun_wide, by=c('arm_1', 'arm_2')) %>% 
  mutate(bubble = paste0('mashgun_', seq_len(n())))

mashmap_unique <-
  mashmap_wide %>% 
  anti_join(shared_bubbles, by=c('arm_1', 'arm_2'))%>% 
  mutate(bubble = paste0('mashmap_',seq_len(n())))

bubblegun_unique <-
  bubblegun_wide %>% 
  anti_join(shared_bubbles, by=c('arm_1', 'arm_2'))%>% 
  mutate(bubble = paste0('bubblegun_', seq_len(n())))

# 
# mashmap_unique_df <-
#   mashmap_unique %>% 
#   tidyr::pivot_longer(cols = c(arm_1, arm_2), names_to = 'bubble_arm', values_to = 'unitig')
# 
# bubblegun_unique_df <-
#   bubblegun_unique %>% 
#   tidyr::pivot_longer(cols = c(arm_1, arm_2), names_to = 'bubble_arm', values_to = 'unitig')
# 
# shared_df <-
#   shared_bubbles %>% 
#   tidyr::pivot_longer(cols = c(arm_1, arm_2), names_to = 'bubble_arm', values_to = 'unitig')


homology_df <-
  bind_rows(
    mashmap_unique,
    bubblegun_unique,
    shared_bubbles
  ) 



# Second Check Cleanup ----------------------------------------------------

# filter out bubbles that contain unitigs in multiple bubbles
homology_df <-
  homology_df %>% 
  filter(!is_duplicate_pair(arm_1, arm_2)) %>%
  tidyr::pivot_longer(
    cols = c(arm_1, arm_2),
    names_to = 'bubble_arm',
    values_to = 'unitig'
  ) 

bubbles_to_remove <- 
  homology_df %>%
  group_by(unitig) %>%
  filter(n() != 1) %>% 
  arrange(bubble) %>%  # alphabetically, favor bubblegun over mashmap
  slice_head(n=1) %>% 
  ungroup() %>% 
  pull(bubble)

homology_df <-
  homology_df %>% 
  filter(!(bubble %in% bubbles_to_remove))


# Export ------------------------------------------------------------------


#bubble double  check
unitigs_per_bubble <-
  homology_df %>% 
  count(bubble) %>% 
  pull()

bubbles_per_unitig <-
  homology_df %>% 
  count(unitig) %>% 
  pull()

stopifnot(all(unitigs_per_bubble == 2))
stopifnot(all(bubbles_per_unitig == 1))

# output <- 'homology/NA18989/NA18989_combined_homology.tsv'
readr::write_tsv(homology_df, output)
