args <- commandArgs(trailingOnly = FALSE)

# Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--bams',
    
    ## Output
    '--output',
    '--log',
    
    ## Params
    '--threads',
    '--max-unitig-cov',
    '--min-ss-cov'
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
library(tidyr)
library(Rsamtools)

# Parsing -----------------------------------------------------------------

stopifnot(all(expected_args %in% args))

## Input
bam_files <- get_values("--bams", singular=FALSE)

## Parameters
n_threads <-      as.numeric(get_values('--threads'))
max_unitig_cov <- as.numeric(get_values('--max-unitig-cov'))
min_ss_cov <-     as.numeric(get_values("--min-ss-cov"))

## Output
output <- get_values('--output')


# Parallel Computing ------------------------------------------------------

if(n_threads > 1){
  if(supportsMulticore()) {
    plan(multicore,    workers=n_threads)
  } else {
    plan(multisession, workers=n_threads)
  }
}
# Filtering ---------------------------------------------------------------

ss_lib_names <-
  map_chr(bam_files, function(x) gsub('.mdup.bam$', '', basename(x)))

bam_files %<>%
  set_names(ss_lib_names)

### count.wc.bam
alignments <- future_map(bam_files, function(bam) {
  aln <- scanBam(file = bam,
                 param = ScanBamParam(
                   what = c('qname', 'rname', 'strand'), #, 'flag', 'pos', 'mapq'),
                   flag = scanBamFlag(
                     isSupplementaryAlignment = F,
                     isSecondaryAlignment = F,
                     isDuplicate = F
                   )
                 ))[[1]] %>% as_tibble()
  
  # Keep only successfully mapped reads
  aln %<>%
    filter(strand %in% c('+', '-'))
  
  # Remove ss_reads that map to too many unitigs. (Slow step)
  aln %<>%
    group_by(qname) %>%
    filter(n() <= max_unitig_cov) %>%
    ungroup()
  
  return(aln)
})

alignments %<>% 
  bind_rows(.id='lib_name')

counts_df <-
  alignments %>%
  group_by(lib_name, rname) %>%
  summarise(w = sum(strand == '+'), c = sum(strand == '-')) %>%
  ungroup()

### get_representative_counts
# counts_df <- bind_rows(counts, .id = 'lib_name')

# # expanding the counts table to have all possible unitig/libs (synchronizing the
# # set of unitigs for all single-cells)
# counts <-
#   counts_df %>%
#   complete(lib_name, rname, fill = list(w = 0, c = 0)) %>%
#   group_split(lib_name)

# Selecting unitigs with at least `min_ss_cov` ss reads mapped to them
high_coverage_rnames <-
  counts_df %>%
  group_by(rname) %>%
  summarise(ss_cov = sum(w + c)) %>%
  ungroup() %>%
  filter(ss_cov >= min_ss_cov)

counts_df %<>%
  semi_join(high_coverage_rnames, by = 'rname')

# choosing the set of reads/unitigs with a reasonable fraction of ww/cc strand
# states over all single-cells. In general, it is expected that a chromosome
# should inherit ww and cc 25% of the time, and wc 50% of the time

# TODO statistical test instead of between(, 0.25 ,75)?

# TODO many libraries have very few reads that map to a unitig, which can make
# for unreliable calculations of is_probably_ww and is_probably_cc
state_counts <-
  counts_df %>%
  mutate(is_probably_ww = w / (w + c) > 0.8,
         is_probably_cc = w / (w + c) < 0.2) %>%
  group_by(rname) %>%
  summarise(
    non_zero_cov = sum((w + c) > 0), # number of libraries that cover a unitig
    probably_ww = sum(is_probably_ww),
    probably_cc = sum(is_probably_cc)
  ) %>%
  ungroup() %>%
  mutate(probably_not_wc = probably_ww + probably_cc) %>%
  mutate(probably_not_wc_frac = probably_not_wc / non_zero_cov)

state_counts %<>%
  filter(between(probably_not_wc_frac, 0.25, 0.75))

counts_df %<>%
  semi_join(state_counts, by = 'rname')

# TODO Are we going to check libraries for W/C mappings? Maybe with bubbles
# somehow? Maybe as a different rule after the bwa fastmap step?


### Export
out <- alignments %>% 
  semi_join(counts_df, by=c('lib_name', 'rname')) %>% 
  distinct(qname, rname)

readr::write_tsv(out, output)

# Close Parallel Workers --------------------------------------------------

plan(sequential)
