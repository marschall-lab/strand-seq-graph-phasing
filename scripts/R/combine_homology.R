
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
# mashmap <- 'homology/NA18989/NA18989_mashmap_unitig-matches.tsv'
mashmap_df <- readr::read_delim(mashmap, delim='\t', col_names = c('unitig_1', 'unitig_2'))

# filter out identical pairs, just with unitig_1 and unitig_2 switched
mashmap_df <-
  mashmap_df %>%
  filter(!is_duplicate_pair(unitig_1, unitig_2))

# sort bubble arms
mashmap_df <-
  mashmap_df %>% 
  mutate(unitig_1 = map2_chr(unitig_1, unitig_2, min),
         unitig_2 = map2_chr(unitig_1, unitig_2, max))

# tidy
mashmap_df <-
  mashmap_df %>% 
  mutate(bubble = paste0('bubble_', 1:n())) %>% 
  tidyr::pivot_longer(cols = c('unitig_1', 'unitig_2'), names_to = 'bubble_arm', values_to ='unitig') %>% 
  mutate(bubble_arm = stringr::str_replace(bubble_arm, '^unitig_', 'arm_'))

# TODO experiment with what will happen if a unitig can have more than 1 partner
# bubble unitig. this may happen in the case of, eg a bubble broken in the middle

# Experimenting has occurred. it also leads to homology being mapped across
# chromosomes, especially with the more degenerate acrocentric chromosomes, but
# generally with the degenerate chromosome regions. Unfortunately, this means
# some genuine homology, eg like broken bubbles, will be thrown out to prevent
# consequences of such errors. There is likely some potenential at rescuing at
# least some information from these instances, but that is not currently done.


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

# This step works because the bubble arms are sorted, such that arm_1 < arm_2
# for both dataframes
wide_bubblegun_df <-
  bubblegun_df %>% 
  tidyr::pivot_wider(names_from = 'bubble_arm', values_from = 'unitig')

wide_mashmap_df <-
  mashmap_df %>% 
  tidyr::pivot_wider(names_from = 'bubble_arm', values_from = 'unitig')

bubblegun_df <-
  anti_join(wide_bubblegun_df, wide_mashmap_df, by=c('arm_1', 'arm_2')) %>% 
  tidyr::pivot_longer(cols = c(arm_1, arm_2), names_to = 'bubble_arm', values_to = 'unitig')

# renumber 
bubblegun_df <-
  bubblegun_df %>% 
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
