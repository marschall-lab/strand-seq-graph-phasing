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
parser$add_argument('--input', help = 'bwa fastmap files', required=TRUE, nargs = '+')
parser$add_argument('--lengths', required=TRUE, nargs = 1)
parser$add_argument('--output', nargs = 1)
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

# Count Alignments ---------------------------------------------------------

lib_names <-
  map_chr(args$input, function(x)
    gsub('_maximal_unique_exact_match.tsv$', '', basename(x)))

# There is a common warning stating "incomplete final line found on". I
# think it is connected to the source of the extra // that appear in the files
exact_match_counts_df <- import_mapper(args$input, function(x) {
  cat(paste('counting bwa-fastmap alignments in', basename(x), '\n'))
  aln <- extract_exact_matches(x)
  
  aln <-
    aln %>% 
    select(-qname)

  if(args$aggregate_alignments) {

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

if(args$threads > 1) {
  # close workers
  plan(sequential)
}


# Export ------------------------------------------------------------------

readr::write_csv(exact_match_counts_df, args$output)
