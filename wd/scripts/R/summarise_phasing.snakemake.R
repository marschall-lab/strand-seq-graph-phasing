args <- commandArgs(trailingOnly = FALSE)

# args <-c(
#  "/gpfs/project/projects/medbioinf/projects/mihen108/wd/.snakemake/conda/43452f50fb7d157c7f3e7e0cf48c0de7/lib/R/bin/exec/R",
#  "--slave"                                                                                                                 ,
#  "--no-restore"                                                                                                            ,
#  "--vanilla"                                                                                                               ,
#  "--file=scripts/R/unmerged_SaaRclust_by_gfa.snakemake.R"                                                            ,
#  "--args"                                                                                                                  ,
#  "--phase-files"                                                                                                                   ,
#  'HG002/phasing/phased_unitigs.csv',
#  "--gfa"                                                                                                                  ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component1.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component2.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component3.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component4.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component5.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component6.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component7.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component8.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component9.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component10.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component11.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component12.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component13.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component14.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component15.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component16.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component17.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component18.gfa"                                                   ,
#  "HG002/assembly.homopolymer-compressed_graph_components/component19.gfa"                                                   ,
#  "--output"                                                                                                         ,
#  "HG002/phasing/phasing_summary.csv"                                                                                                        ,
#  "--log"                                                                                                                   ,
#  "log/SaaRclust_by_gfa_HG002_initial_clusters.log")
# Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--phase-files',
    '--gfa',
    ## Output
    '--output',
    '--log'
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



# Log ---------------------------------------------------------------------
## Log
log_path <- get_values('--log')
log <- file(log_path, open='wt')
# sink(file=log, type='message')
# sink(file=log, type='output')

print(args)

# Library -----------------------------------------------------------------

library(dplyr)
library(data.table)
# Parsing -----------------------------------------------------------------




stopifnot(all(expected_args %in% args))

## Input
gfas <- get_values("--gfa", singular=FALSE)
phase_files <- get_values("--phase-files", singular=FALSE)

## Output
output_path <- get_values('--output')

# Functions ---------------------------------------------------------------

read_segment_lengths_from_gfa <- function(x) {
  split_lines <-
    grep('^S',
       readLines(x),
       value = TRUE) %>%
    strsplit("\t") 
  
  tig_lengths <- 
    split_lines %>%
    vapply(function(x) nchar(x[[3]]), integer(1)) 
  
  tig_names <- 
    split_lines %>%
    vapply(function(x) x[[2]], character(1))
  
  names(tig_lengths) <- tig_names
  
  return(tig_lengths)
}

basename_no_ext <- function(x) {
  tools::file_path_sans_ext(basename(x))
}


dt_from_named_vector <- function(x, value='value', name='name') {
  # Probably a less dumb looking way to do this
  out <-
    data.table(value = x, name = names(x)) %>% 
    setnames(old = c('value', 'name'), new = c(value, name))

    return(out)
}

invert_sign <- function(x) {
  stopifnot(x %in% c('+', '-'))
  ifelse(x == '+', '-', '+')
}

dir_create_if_does_not_exist <- function(x) {
  if (!file.exists(x)) {
    dir.create(x)
    return(invisible(TRUE))
  } else {
    return(invisible(FALSE))
  }
}

# Segment Info ------------------------------------------------------------



rlengths_by_gfa <-
  lapply(gfas, read_segment_lengths_from_gfa) 

names(rlengths_by_gfa) <- basename_no_ext(gfas)

# filter(gfa != 'component14') %>% 
segment_info <- 
  lapply(rlengths_by_gfa, function(x) data.table(unitig_name=names(x), unitig_length=x)) %>% 
  rbindlist(idcol='gfa') # %>% 
  # filter(gfa != 'component14') 

segment_info <-
  rbindlist(list(segment_info, mutate(segment_info, gfa = 'all')))

# segment_info <-
#   segment_info %>% 
#   mutate(is_long_unitig = unitig_length >= 1e6)

# Import Phasings ---------------------------------------------------------


names(phase_files) <- basename_no_ext(phase_files)

phasings <- 
  lapply(phase_files, fread) 

# Summary Stats -----------------------------------------------------------


calc_summary_stats <- function(phase_file, segment_info){
  sucessfully_phased_values <- c('H1', 'H2')
  
  segment_info <-
    segment_info %>% 
    mutate(attempted_to_phase = unitig_name %in% phase_file$unitig_name) %>% # I don't want to rely on a NA values in any particular column to determine this
    full_join(phase_file, by = "unitig_name")
  
  segment_info <- 
    segment_info %>% 
    mutate(sucessfully_phased = haplotype %in% sucessfully_phased_values)
  
  summary_stats <-
    segment_info %>%
    group_by(gfa) %>%
    summarise(
      total_bp = sum(unitig_length),
      total_bp_attempted = sum(unitig_length[attempted_to_phase]),
      total_bp_not_attempted = sum(unitig_length[!attempted_to_phase]),
      total_bp_phased = sum(unitig_length[sucessfully_phased])
    ) 
  
  summary_stats <-
    summary_stats %>%
    mutate(
      percent_attempted = total_bp_attempted / total_bp,
      percent_not_attempted = total_bp_not_attempted / total_bp,
      percent_phased = total_bp_phased / total_bp,
      percent_phased_of_attempted = ifelse(total_bp_attempted > 0, total_bp_phased / total_bp_attempted, 0)
    ) %>% 
    mutate(across(starts_with('percent'), ~round(.x*100, 1)))
  
  summary_stats <-
    summary_stats %>%
    select(-total_bp_attempted, -total_bp_not_attempted, -total_bp_phased)
  return(summary_stats)
}



# Output ------------------------------------------------------------------


out <-
  lapply(phasings, calc_summary_stats, segment_info) %>%
  rbindlist(idcol = 'phase_file') 

dir_create_if_does_not_exist(dirname(output_path))
fwrite(out, output_path)


