# Command Line ------------------------------------------------------------


## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

# source('/Users/henglinm/Documents/GitHub/strand-seq-graph-phasing/scripts/R/module_utils/phasing_test_args.R')

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    '--mem-alignment-bams',
    '--fastmap-alignments',

    '--output-mem',
    '--output-fastmap',
    '--output-mem-raw',
    '--output-fastmap-raw',

    '--threads'
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

library(Rsamtools)

source(file.path(get_script_dir(), "module_utils/utils.R"))
source(file.path(get_script_dir(), "module_utils/phasing_utils.R"))
# Parsing -----------------------------------------------------------------

## Input

mem_alignment_files <- get_values("--mem-alignment-bams", singular=FALSE)
fastmap_alignment_files <- get_values("--fastmap-alignments", singular=FALSE)

n_threads <- as.numeric(get_values('--threads'))
stopifnot(n_threads >= 1)

## Output
output_mem <- get_values('--output-mem')
output_fastmap <- get_values('--output-fastmap')
output_mem_raw <- get_values('--output-mem-raw')
output_fastmap_raw <- get_values('--output-fastmap-raw')

# Load Headers ------------------------------------------------------------

cat('Scanning Bam Header:\n')
# All bam headers should be the same right? Only need one
unitig_lengths <- scanBamHeader(mem_alignment_files[[1]], what='targets')[[1]]$targets
unitig_lengths_df <- tibble::enframe(unitig_lengths, name='unitig', value='length')


# Count Alignments ---------------------------------------------------------

if(n_threads > 1) {
  library(furrr)
  plan(multisession, workers=n_threads)
  import_mapper <- furrr::future_map  
} else {
  import_mapper <- purrr::map
}

### Count mem Alignments ------------------------------------------------------
lib_names <- map_chr(mem_alignment_files, function(x) gsub('.mdup.bam$', '', basename(x)))

counts_df <- import_mapper(mem_alignment_files, function(bam){
  cat(paste('counting bwa-mem alignments in', basename(bam), '\n'))
  aln <- scanBam(file = bam,
                 param = ScanBamParam(
                   what = c('qname', 'rname', 'strand', 'pos', 'qwidth', 'mrnm', 'isize', 'mapq'),
                   flag = scanBamFlag(
                     isSupplementaryAlignment = FALSE,
                     isSecondaryAlignment = FALSE,
                     isDuplicate = FALSE,
                     # for the purpose of determining qname/direction mapping, having both mates simply provides redundant information?
                     isFirstMateRead = TRUE,
                     isProperPair = TRUE
                   )
                 ))
  
  aln <- as_tibble(aln[[1]])
  
  #TODO add check for nrow() > 0. nrow() = 0 can happen with meme alignments if
  #one of the Strand-seq files is malformed.
  
  # filter out alignments to short unitigs
  aln <- 
    aln %>%  
    mutate(rname = as.character(rname)) # default is as.factor import?
  
  # Keep only reads that successfully aligned
  aln <- 
    aln %>% 
    filter(strand %in% c('+','-'))
  
  # filter any ss_reads that map to too many unitigs
  # dplyr::filter with lots of groups can be very slow -> duplicated is faster
  duplicated_qnames <-
    with(aln, qname[duplicated(qname)])
  
  if(any(duplicated_qnames)) {
    cat('Duplicated qnames found in:', basename(bam), '\n')
    aln <-
      aln %>%
      filter(!(qname %in% duplicated_qnames))
  }
  
  # to simplify, only keep alignments where both mates land on the same rname
  aln <- 
    aln %>% 
    filter(rname == mrnm)
  

  
  return(aln)
  
  
})

counts_df <-
  counts_df %>%
  set_names(lib_names) %>% 
  bind_rows(.id = 'lib') %>%
  dplyr::rename(unitig = rname) 

# Counting
aggregated_counts_df <-
  counts_df %>%
  # filter(mapq > 0) %>% 
  group_by(lib, unitig) %>%
  summarise(c = sum(strand == '+'), w = sum(strand == '-'), .groups="drop")

aggregated_counts_df <-
  aggregated_counts_df %>%
  right_join(unitig_lengths_df, by = 'unitig') %>% 
  tidyr::complete(lib, tidyr::nesting(unitig, length), fill=list(c=0, w=0)) %>% 
  mutate(n = c+w) %>% 
  filter(!is.na(lib))


# raw_counts_df <- counts_df
### Count fastmap Alignments ------------------------------------------------
# fastmap_alignment_files <- list('bwa_alignments/fastmap/HG00171/000000_maximal_unique_exact_match.tsv', 'bwa_alignments/fastmap/HG00171/000001_maximal_unique_exact_match.tsv')
lib_names <-
  map_chr(fastmap_alignment_files, function(x)
    gsub('_maximal_unique_exact_match.tsv$', '', basename(x)))

# There is a common warning stating "incomplete final line found on". I
# think it is connected to the source of the extra // that appear in the files
exact_match_counts_df <- import_mapper(fastmap_alignment_files, function(x) {
  cat(paste('counting bwa-fastmap alignments in', basename(x), '\n'))
  out <- extract_exact_matches(x)
  
  return(out)
})

exact_match_counts_df <-
  exact_match_counts_df %>%
  set_names(lib_names) %>% 
  bind_rows(.id = 'lib') 

aggregated_exact_match_counts_df <-
  exact_match_counts_df %>%
  group_by(lib, unitig) %>%
  summarise(c = sum(strand == '+'), w = sum(strand == '-'), .groups="drop") 

aggregated_exact_match_counts_df <-
  aggregated_exact_match_counts_df %>%
  right_join(unitig_lengths_df, by = 'unitig') %>% 
  tidyr::complete(lib, tidyr::nesting(unitig, length), fill=list(c=0, w=0)) %>% 
  mutate(n = c+w) %>% 
  filter(!is.na(lib))


# raw_aggregated_exact_match_counts_df <- aggregated_exact_match_counts_df

if(n_threads > 1) {
  # close workers
  plan(sequential)
}


# Export ------------------------------------------------------------------

readr::write_csv(aggregated_counts_df, output_mem)
readr::write_csv(aggregated_exact_match_counts_df, output_fastmap)
readr::write_csv(counts_df, output_mem_raw)
readr::write_csv(exact_match_counts_df, output_fastmap_raw)
