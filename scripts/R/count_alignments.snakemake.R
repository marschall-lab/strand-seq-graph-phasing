
# Helpers -----------------------------------------------------------------

read_targets_from_sam_lines <- function(sam_lines) {
  # sam_lines ~ from readr::read_lines of a sam file.

  targets <-
    readr::read_tsv(
      I(stringr::str_subset(sam_lines, '^\\@SQ')),
      col_names = c('linetag', 'rname', 'rlen')
    )


  targets <-
    targets %>%
    tidyr::separate_wider_delim(c(rname, rlen), delim = ':', names_sep = '_') %>%
    select(rname_2, rlen_2) %>%
    rename(unitig = rname_2, length = rlen_2)

  return(targets)

}


read_alns_from_sam_lines <- function(sam_lines) {
  sam_alns <-
    readr::read_tsv(
      I(stringr::str_subset(sam_lines, '^\\@', negate=TRUE)),
      col_names = c(
        'qname',
        'flag',
        'rname',
        'pos',
        'mapq',
        'cigar',
        'rnext',
        'pnext',
        'tlen',
        'seq',
        'qual'
      )
    )

  return(sam_alns)
}


filter_alignments <- function(aln) {

  aln <-
    aln %>%
    filter(sign(tlen) %in% c(-1, 1))

  is_duplicated_qname <- duplicated(aln$qname)

  if(any(is_duplicated_qname)) {
    cat('Duplicated qnames found in:', basename(bam), '\n')
    aln <-
      aln %>%
      filter(!is_duplicated_qname)
  }

  aln <-
    aln %>%
    filter(rnext == '=')

  return(aln)
}


# Command Line ------------------------------------------------------------


## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    '--mem-alignment-bams',
    '--fastmap-alignments',

    '--aggregate-alignments',
    '--output-mem',
    '--output-fastmap',

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

source(file.path(get_script_dir(), "module_utils/utils.R"))
source(file.path(get_script_dir(), "module_utils/phasing_utils.R"))
# Parsing -----------------------------------------------------------------

## Input

mem_alignment_files <- get_values("--mem-alignment-bams", singular=FALSE)
fastmap_alignment_files <- get_values("--fastmap-alignments", singular=FALSE)

n_threads <- as.numeric(get_values('--threads'))
stopifnot(n_threads >= 1)

aggregate_alignments <- as.logical(get_values('--aggregate-alignments'))

## Output
output_mem <- get_values('--output-mem')
output_fastmap <- get_values('--output-fastmap')




# Load Headers ------------------------------------------------------------

cat('Scanning Header:\n')
# All headers should be the same right? Only need one

sam_lines <- readr::read_lines(mem_alignment_files[[1]])
unitig_lengths_df <- read_targets_from_sam_lines(sam_lines)
rm(sam_lines)

# Count Alignments ---------------------------------------------------------

if(n_threads > 1) {
  library(furrr)
  plan(multisession, workers=n_threads)
  import_mapper <- furrr::future_map  
} else {
  import_mapper <- purrr::map
}

### Count mem Alignments ------------------------------------------------------
lib_names <- map_chr(mem_alignment_files, function(x) gsub('.mdup.filt.sam$', '', basename(x)))

counts_df <- import_mapper(mem_alignment_files, function(sam){
  cat(paste('counting bwa-mem alignments in', basename(sam), '\n'))

  aln_lines <- readr::read_lines(sam)
  aln <- read_alns_from_sam_lines(aln_lines)
  rm(aln_lines)

  # TODO add check for nrow() > 0. nrow() = 0 can happen with mem alignments if
  # one of the Strand-seq files is malformed.


  aln <- filter_alignments(aln)
  
  # match necessary columns from old naming scheme
  aln <-
    aln %>%
    mutate(strand = case_when(sign(tlen) == 1 ~ '+', sign(tlen) == -1 ~ '-', TRUE ~ NA)) %>%
    rename(unitig = rname) %>%
    select(unitig, strand, pos, mapq)
  
  
  if(aggregate_alignments) {

      aln <-
        aln %>%
        group_by(unitig) %>%
        summarise(c = sum(strand == '+'), w = sum(strand == '-'), .groups="drop")

      # explicit 0s
      aln <-
        aln %>%
        right_join(unitig_lengths_df) %>%
        mutate(c = coalesce(c, 0), w = coalesce(w, 0)) %>%
        mutate(n = c+w)
  }
  
  return(aln)
  
  
})

counts_df <-
  counts_df %>%
  set_names(lib_names) %>%
  bind_rows(.id = 'lib')  

### Count fastmap Alignments ------------------------------------------------
lib_names <-
  map_chr(fastmap_alignment_files, function(x)
    gsub('_maximal_unique_exact_match.tsv$', '', basename(x)))

# There is a common warning stating "incomplete final line found on". I
# think it is connected to the source of the extra // that appear in the files
exact_match_counts_df <- import_mapper(fastmap_alignment_files, function(x) {
  cat(paste('counting bwa-fastmap alignments in', basename(x), '\n'))
  aln <- extract_exact_matches(x)
  
  aln <-
    aln %>% 
    select(-qname)

  if(aggregate_alignments) {

      aln <-
        aln %>%
        group_by(unitig) %>%
        summarise(c = sum(strand == '+'), w = sum(strand == '-'), .groups="drop")

      aln <-
        aln %>%
        right_join(unitig_lengths_df) %>%
        mutate(c = coalesce(c, 0), w = coalesce(w, 0)) %>%
        mutate(n = c+w)
  }
  
  return(aln)
})

exact_match_counts_df <-
  exact_match_counts_df %>%
  set_names(lib_names) %>% 
  bind_rows(.id = 'lib') 

if(n_threads > 1) {
  # close workers
  plan(sequential)
}


# Export ------------------------------------------------------------------

readr::write_csv(counts_df, output_mem)
readr::write_csv(exact_match_counts_df, output_fastmap)
