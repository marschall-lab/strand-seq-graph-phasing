
# Helpers -----------------------------------------------------------------


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

# Library -----------------------------------------------------------------

library(dplyr)
library(purrr)
library(argparse)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  needle <- '--file='
  script_path <- gsub(needle, '',  args[grepl(needle, args)])
  return(dirname(script_path))
}

source(file.path(get_script_dir(), "module_utils/utils.R"))
source(file.path(get_script_dir(), "module_utils/phasing_utils.R"))
# Parsing -----------------------------------------------------------------

print(commandArgs(trailingOnly = FALSE))

parser <- ArgumentParser(description = 'Count Alignments')
parser$add_argument('--input', help = 'bwa mem sam alignment files', required=TRUE, nargs = '+')
parser$add_argument('--lengths', required=TRUE, nargs = 1)
parser$add_argument('--output', required=TRUE, nargs = 1)
parser$add_argument('--aggregate-alignments', help = 'should alignment counts be aggregated by unitig?', action = 'store_const', const=TRUE, default=FALSE)
parser$add_argument('--threads', nargs=1, default=1L, type='integer')

args <- parser$parse_args()
print(args)

stopifnot(args$threads >= 1)

# Load Headers ------------------------------------------------------------

unitig_lengths_df <- readr::read_tsv(args$lengths, col_names=c('unitig', 'length'))

# Mapper ---------------------------------------------------------

if(args$threads > 1) {
  library(furrr)
  plan(multisession, workers=args$threads)
  import_mapper <- furrr::future_map  
} else {
  import_mapper <- purrr::map
}

### Count Alignments ------------------------------------------------------
lib_names <- map_chr(args$input, function(x) gsub('.mdup.filt.sam$', '', basename(x)))

counts_df <- import_mapper(args$input, function(sam){
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
  
  
  if(args$aggregate_alignments) {

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

if(args$threads > 1) {
  # close workers
  plan(sequential)
}

counts_df <-
  counts_df %>%
  set_names(lib_names) %>%
  bind_rows(.id = 'lib')  



# Export ------------------------------------------------------------------

readr::write_csv(counts_df, args$output)

