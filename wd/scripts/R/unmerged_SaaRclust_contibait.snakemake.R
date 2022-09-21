args <- commandArgs(trailingOnly = FALSE)
#
# args <-c(
#  "/gpfs/project/projects/medbioinf/projects/mihen108/wd/.snakemake/conda/43452f50fb7d157c7f3e7e0cf48c0de7/lib/R/bin/exec/R",
#  "--slave"                                                                                                                 ,
#  "--no-restore"                                                                                                            ,
#  "--vanilla"                                                                                                               ,
#  "--file=scripts/R/unmerged_SaaRclust_by_gfa.snakemake.R"                                                            ,
#  "--args"                                                                                                                  ,
#  "--bam"                                                                                                                   ,
#  "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20376.mdup.bam"                          ,
#  "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20368.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20484.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20341.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20337.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20387.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20364.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20458.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20306.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20466.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20417.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20303.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20473.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20357.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20328.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20350.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20405.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20414.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20492.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20339.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20347.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20379.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20469.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20327.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20464.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20391.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20440.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20483.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20392.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20419.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20345.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20352.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20486.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20359.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20307.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20343.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20319.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20452.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20435.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20308.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20362.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20355.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20421.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20407.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20450.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20403.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20487.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20305.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20485.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20425.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20313.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20325.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20381.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20443.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20389.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20454.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20393.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20431.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20342.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20358.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20488.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20432.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20351.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20336.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20479.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20340.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20367.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20335.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20301.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20434.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20353.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20356.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20329.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20424.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20422.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20445.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20491.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20331.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20430.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20377.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20374.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20315.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20439.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20490.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20361.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20467.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20428.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20334.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20363.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20494.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20378.mdup.bam"                          ,
# "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20388.mdup.bam"                          ,
#  "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20318.mdup.bam"                          ,
#  "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_01PE20332.mdup.bam"                          ,
#  "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20456.mdup.bam"                          ,
#  "HG002/unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_hpg_ilnxs-80pe_02PE20416.mdup.bam"                          ,
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
#  "--output-prefix"                                                                                                         ,
#  "HG002/SaaRclust/"                                                                                                        ,
#  "--em-iter"                                                                                                               ,
#  "100"                                                                                                                     ,
#  "--threads"                                                                                                               ,
#  "7"                                                                                                                       ,
#  "--log"                                                                                                                   ,
#  "log/SaaRclust_by_gfa_HG002_initial_clusters.log")

print(args)
# Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--bam',
    '--gfa',
    ## Output
    '--output-prefix',
    '--log',

    ## Params
    '--em-iter',
    '--threads',
    '--segment-length-threshold'
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

log_path <- get_values('--log', singular=FALSE)
log <- file(log_path, open='wt')
sink(file=log, type='message')
sink(file=log, type='output')

print(args)

# Library -----------------------------------------------------------------

# .libPaths(c(.libPaths(), 'utils/R-packages/'))
library(dplyr)
#library(SaaRclust)
library(data.table)
library(Rsamtools)
library(doParallel)
library(matrixStats)
# library(assertthat)
#library(biovizBase)
library(contiBAIT)
# Parsing -----------------------------------------------------------------




stopifnot(all(expected_args %in% args))



## Input

input.alignment.files <- get_values("--bam", singular=FALSE)
gfas <- get_values("--gfa", singular=FALSE)

inputfolder <- dirname(input.alignment.files[[1]]) # bwa_ss_unitigs/

## Parameters
EM.iter <-       as.numeric(get_values("--em-iter"))
numCPU <-        as.numeric(get_values("--threads")) # snakemake@threads[[1]]
segment_length_threshold <- as.numeric(get_values('--segment-length-threshold'))

## Output
# outputfolder <- dirname(dirname(snakemake@output[["hard_clust"]])) # SaarClust/
outputfolder <- get_values('--output-prefix')

MLclust.file <-     paste0(outputfolder, 'Clusters/MLclust.data') # snakemake@output[["ML_clust"]]
ss.clust.file <-    paste0(outputfolder, 'Clusters/ss_clusters.data') # snakemake@output[["ss_clust"]]
clust.pairs.file <- paste0(outputfolder, 'Clusters/clust_partners.txt') # snakemake@output[["clust_pairs"]]
wc.cells.file <-    paste0(outputfolder, 'Clusters/wc_cells_clusters.data') # snakemake@output[["wc_cells_clusters"]]
exclusion_file <-    paste0(outputfolder, 'Clusters/SaaRclust_unclustered_rnames.tsv')

parsed_args <-
  list(
    inputfolder = inputfolder,
    outputfolder = outputfolder,
    gfas=gfas,
    input.alignment.files = input.alignment.files,
    EM.iter = EM.iter,
    MLclust.file = MLclust.file,
    ss.clust.file = ss.clust.file,
    clust.pairs.file = clust.pairs.file,
    wc.cells.file = wc.cells.file,
    numCPU = numCPU
  )
print("Parsed Args:")
print(parsed_args)


# Sourcing ----------------------------------------------------------------



all.sources <-
  c(
    "module_saarclust/calcProbs.R",
    "module_saarclust/countDirectionalReads.R",
    "module_saarclust/dumped_functions.R",
    "module_saarclust/EMclust.R",
    "module_saarclust/export.R",
    "module_saarclust/findClusterPartners.R",
    "module_saarclust/hardClust.R",
    "module_saarclust/helperFuctions.R",
    "module_saarclust/import.R",
    "module_saarclust/importReads.R",
    "module_saarclust/SaaRclust_evaluation_plots.R",
    "module_saarclust/SaaRclust.R",
    "module_saarclust/timedMessage.R",
    "module_saarclust/utils.R",
    "module_saarclust/wrapper_parallel.R",
    "module_saarclust/wrapper.R"
  )

# Path Handling

all.sources <- paste0(get_script_dir(), '/', all.sources)

invisible(sapply(all.sources, source))

# Functions ---------------------------------------------------------------


read_rnames_from_gfa <- function(x) {
  grep('^S',
       readLines(x),
       value = TRUE) %>%
    strsplit("\t") %>%
    vapply(function(x) x[[2]], character(1)) # grab the second field in every line
}

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

dir_create_if_does_not_exist <- function(x) {
  if (!file.exists(x)) {
    dir.create(x)
    return(invisible(TRUE))
  } else {
    return(invisible(FALSE))
  }
}

subset_rownames <- function(m, x) {

  rn <- rownames(m)
  idx <- rn %in% x
  out <- m[idx, ,drop=FALSE]
  rownames(out) <- rn[idx]
  return(out)
}

convert_dt_col_to_rownames <- function(dt, col) {
  rn <- dt[[col]]
  cols_to_keep <- names(dt) != col
  dt <- subset(dt, select=cols_to_keep)
  rownames(dt) <- rn
  return(dt)
}

spoof <- function(x){
  paste0(x , ':0-0')
}

unspoof <- function(x) {
  gsub(':0-0$', '', x)
}

spoof_rownames <- function(m) {
  rownames(m) <- spoof(rownames(m))
  return(m)
}

unspoof_rownames <- function(m) {
  rownames(m) <- unspoof(rownames(m))
  return(m)
}

spoof_names <- function(x) {
  names(x) <- spoof(names(x))
  return(x)
}

unspoof_names <- function(x) {
  names(x) <- unspoof(names(x))
  return(x)
}

get_read_counts_matrix <- function(counts.dt){
  read_counts <- counts.dt

  read_counts[, read_count := (w+c)]

  mat <- dcast(read_counts, rname~lib, value.var='read_count')

  row.names <- mat[, rname]
  mat[, rname:=NULL]
  mat <- as.matrix(mat)
  rownames(mat) <- row.names
  mat[is.na(mat)] <- 0

  mat
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

# Useful Variables --------------------------------------------------------

cl <- makeCluster(numCPU)
doParallel::registerDoParallel(cl)
rlengths_by_gfa <-
  foreach (gfa = gfas, .packages = 'dplyr') %dopar% {
    read_segment_lengths_from_gfa(gfa)
  }
stopCluster(cl)
names(rlengths_by_gfa) <- basename_no_ext(gfas)

long_rnames_by_gfa <-
  lapply(rlengths_by_gfa, function(x) names(x[x >= segment_length_threshold]))

long_rnames <-
  long_rnames_by_gfa %>%
  Reduce(c, .)

rnames_by_gfa <-
  lapply(rlengths_by_gfa, names)

rnames_by_gfa_dt <-
  lapply(rnames_by_gfa, function(x) data.table(rname=x)) %>%
  rbindlist(idcol='gfa')

all_rnames <-
  rnames_by_gfa %>%
  Reduce(c, .)

# ref_aln <- fread('ref_aln.csv')
# ref_aln <- setnames(ref_aln, old=c('qname', 'rname'), new=c('rname', 'chrom'))
# ref_aln <- ref_aln[]

# Directories -------------------------------------------------------------


#Create a master output directory
outputfolder.destination <- file.path(outputfolder)
dir_create_if_does_not_exist(outputfolder.destination)

#Directory to store raw read counts and best alignments
rawdata.store <- file.path(outputfolder.destination, 'RawData')
dir_create_if_does_not_exist(rawdata.store)

# Directory to store processed/clustered data
Clusters.store <- file.path(outputfolder.destination, 'Clusters')
dir_create_if_does_not_exist(Clusters.store)

# # Directory to store split_ss ids/clustered data
# ss_Clusters.store <- file.path(outputfolder.destination, 'Clusters', 'ss_clusters')
# dir_create_if_does_not_exist(ss_Clusters.store)

#Directory to store plots
plots.store <- file.path(outputfolder.destination, 'Plots')
dir_create_if_does_not_exist(plots.store)

#Directory to store 'difficult' PacBio reads for later processing [TODO]
trashbin.store <- file.path(outputfolder.destination, 'TrashBin')
dir_create_if_does_not_exist(trashbin.store)

# getting Strand-seq read counts
#reuse existing data if they were already created and saved in a given location
destination1 <- file.path(rawdata.store, "alignments.RData")
destination2 <- file.path(rawdata.store, "read_counts.RData")
destination3 <- file.path(rawdata.store, "read_selected_counts.RData")


# Alignments --------------------------------------------------------------

cat('counting w/c reads...\n')
bam.files <-  input.alignment.files
lib.names <- sapply(bam.files, function(bam) gsub('.mdup.bam$', '', basename(bam)))

cl <- makeCluster(numCPU)
doParallel::registerDoParallel(cl)
alignments <- foreach (bam=bam.files, .packages=c('Rsamtools', 'data.table')) %dopar%{
  cat('counting directional reads in', basename(bam), '\n')
  aln = scanBam(file = bam,
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

  aln <- as.data.table(aln[[1]])
  # Keep only reads that successfully aligned
  aln <- aln[strand %in% c('+','-')]

  # filter any ss_reads that map to too many unitigs
  max.unitig.cov <- 1
  aln[, unitig.cov:=.N, by='qname']
  aln <- aln[unitig.cov <= max.unitig.cov]

  # to simplify, only keep alignments where both mates land on the same rname
  aln <- aln[rname == mrnm]

  return(aln[, .(rname, qname, strand, mapq)])
}
stopCluster(cl)

names(alignments) <- lib.names



# Counts ------------------------------------------------------------------




cl <- makeCluster(numCPU)
doParallel::registerDoParallel(cl)
counts.l.all <-  foreach (aln=alignments, .packages=c('data.table')) %dopar%{
  counts <- aln[, .(c=sum(strand=='+'), w=sum(strand=='-')), by='rname']
  return(counts)
}
stopCluster(cl)

names(counts.l.all) <- lib.names


counts <-
  counts.l.all %>%
  rbindlist(idcol = 'lib') %>%
  tidyr::complete(lib, rname, fill=list(c=0, w=0)) %>%
  as.data.table()

# # Quality Control ---------------------------------------------------------


# subsetting the highest coverage unitigs
# Keep unititgs with minimum coverage
# min.ss.cov <- 5 * length(counts.l.all)
# counts[, ss.cov:=sum(w+c), by=rname]
# counts <- counts[ss.cov >= min.ss.cov]

# rn <- long_rnames
# wfrac.matrix <-
#   counts[rname %in% rn, ] %>%
#   get.wfrac.matrix()
#
#
# # TRUE if wc
# cell_state_matrix <-
#   abs(wfrac.matrix) < 0.8
#
# high_wc_rnames <-
#   cell_state_matrix %>%
#   apply(1, mean) %>%
#   `[`(. >= 0.8) %>%
#   names()
#
# low_wc_rnames <-  # ~ monosome  chromosomes?
#   (!cell_state_matrix) %>%
#   apply(1, mean) %>%
#   `[`(. >= 0.8) %>%
#   names()
#
# counts <- counts[!is.element(rname, c(high_wc_rnames, low_wc_rnames))]
#

# Contibait ---------------------------------------------------------------

rn <- long_rnames
wfrac.matrix <-
  counts[rname %in% rn, ] %>%
  get.wfrac.matrix()

# TODO incoporporate  mapping quality as well?
# weight the unitigs that have more alignments ~ they will be more likely to be
# selected earlier by the contibatit clustering algorithm?
weights_dt <-
  counts[rname %in% rn, .(cov = mean(w + c)), by='rname'] %>%
  mutate(rname = as.character(rname))

weights <- with(weights_dt, setNames(cov, rname))

# order and spoof
wfrac.matrix <-
  wfrac.matrix[weights_dt$rname,] %>%
  spoof_rownames()

weights <-  spoof_names(weights)


strand.freq <- StrandFreqMatrix(wfrac.matrix)

# debugonce(preprocessStrandTable, signature = 'StrandFreqMatrix')
strand.states <- preprocessStrandTable(strand.freq)

# debugonce(clusterContigs, signature = 'StrandStateMatrix')
cl <-
  clusterContigs(strand.states$strandMatrix,
                 recluster = 100,
                 randomWeight = weights[rownames(strand.states$strandMatrix)],
                 clusterBy = 'hetero')

# Monosome Sex Chromosomes ------------------------------------------------

# TODO what if multiple groups of sex detected chromosomes?

# debugonce(findSexGroups, signature = c('LinkageGroupList', 'StrandStateMatrix'))
cl <- findSexGroups(cl, strand.states$strandMatrix)

if(length(grep('^sex', names(cl))) > 1) {
  warning('More than 1 cluster has been identifies as a monosome chromosome cluster?')
}

#### Check Plotting ####
# clusters <-
#   cl %>%
#   lapply(function(x) tibble(rname=x)) %>%
#   bind_rows(.id='clust') %>%
#   mutate(rname = unspoof(rname))
#
# library(ggplot2)
# text_data <-
#   clusters %>%
#   right_join(ref_aln, by='rname') %>%
#   group_by(chrom, clust) %>%
#   summarise(weight = sum(qwidth), .groups = 'drop')
#
# ggplot(text_data) +
#   geom_col(aes(x = chrom, y=weight, fill = clust), color='black') +
#   geom_text(aes(x=chrom, y=weight, label=clust), angle=90)



# Bootstrap with Inversions -----------------------------------------------

# TODO bootstrap or just simple single inversion adequate?
bootstrap_dt <-
  setNames(cl@.Data, cl@names) %>%
  lapply(function(x) data.table(rname = unspoof(x))) %>%
  rbindlist(idcol='clust') %>%
  # left_join(ref_aln) %>%
  inner_join(counts) %>%
  mutate(invert=0)

bootstrap_dt <-
  bootstrap_dt %>%
  mutate(n = w+c) %>%
  mutate(root_rname = rname)


#### BEGIN bootstrap block ####


n_boot <- 0
bootstraps <-
  lapply(seq_len(n_boot), function(i) {
    boot_dt <-
      bootstrap_dt %>%
      mutate(rname = paste0(root_rname, '(', i, ')'))

    boot_dt <-
      boot_dt %>%
      mutate(w = rbinom(n(), size = n, prob = ifelse(n > 0, w/n, 0))) %>%
      mutate(c = n - w)

    boot_dt <-
      boot_dt %>%
      group_by(rname) %>%
      mutate(invert =  rbinom(1, 1, 0.5)) %>%
      ungroup() %>%
      mutate(w = ifelse(invert == 1, n - w, w),
             c = ifelse(invert == 1, n - c, c))

    return(boot_dt)

  })

bootstraps <-
  rbindlist(bootstraps, idcol='bootstrap') %>%
  rbind(bootstrap_dt, fill=TRUE)

#### END bootstrap block ####

#### BEGIN Simple Inversion Block ####
# concat with inverted
inverted_bootstrap_dt <-
  bootstrap_dt %>%
  mutate(invert=1,
         rname = paste0(root_rname, '(', 0, ')'),
         c = n-c,
         w = n-w)


bootstraps <-
  rbind(bootstraps, inverted_bootstrap_dt)
#### END Simple Inversion Block ####




# Monosome Sex Clustering -------------------------------------------------

mono_sex_rnames <- c()
mono_sex <- grepl('^sex', names(cl))
if(any(mono_sex)) {
  # TODO currently  assumes that all monosome clusters come from XY chromosomes,
  # and merges them together. Finda  way to handle each cluster separately in
  # case of other monosome chromsomes? And related, need a way to determine how
  # many clusters are in a given monosome cluster. Maybe other monosomes would
  # be grouped with the XY? Would monosomes be assigned a separate sex cluster
  # instead?

  # TODO the sex clusters are currently not output to a file, and probably should be moved
  # into its own function
  mono_sex_rnames <-
    cl[mono_sex] %>%
    Reduce(c, .) %>%
    unspoof()
  # cutree(4) breaks unexpectedly Error in cutree(., 4) : elements of 'k' must be between 1 and 2
  # mono_sex_counts <-
  #   bootstraps %>%
  #   filter(root_rname %in% mono_sex_rnames)
  #
  # wfm <-
  #   mono_sex_counts %>%
  #   get.wfrac.matrix()
  #
  # # working in PCA space makes things like averaging points easier to reason about?
  # pca <- prcomp(wfm, center=FALSE, scale=FALSE)
  #
  # hc <-
  #   pca$x %>%
  #   dist() %>%
  #   hclust()
  #
  # clusts <-
  #   hc %>%
  #   cutree(4) # 2 chrom x 2 unitig directions
  #
  # # cluster most similar pairs of clusters
  # centroids <-
  #   setNames(nm = unique(clusts)) %>%
  #   lapply(function(clust) {
  #     out <- pca$x[names(clusts)[clusts == clust], , drop = FALSE] %>%
  #       apply(2, mean) %>%
  #       as.matrix()
  #
  #     colnames(out) <- clust
  #
  #     return(out)
  #   }) %>%
  #   Reduce(cbind, .) %>%
  #   t()
  #
  # d <-
  #   dist(centroids) %>%
  #   as.matrix()
  #
  # first_pair <- arrayInd(which.max(d), dim(d), useNames = TRUE)[1, ]
  # second_pair <- setdiff(rownames(centroids), first_pair) %>% as.numeric()
  #
  # sex_clusters <-
  #   list(
  #     '1'=data.table(rname = names(clusts)[clusts %in% first_pair]),
  #     '2'=data.table(rname = names(clusts)[clusts %in% second_pair])
  #   ) %>%
  #   rbindlist(idcol = 'sex_clust')
  #
  # sex_clusters <-
  #   sex_clusters %>%
  #   semi_join(filter(bootstraps, root_rname == rname), by='rname')

  # remove sex clusters from rest of process.
  cl[mono_sex] <- NULL
  bootstraps <-
    bootstraps %>%
    filter(!grepl('^sex', clust))
}


# Excluded Unitigs --------------------------------------------------------

short_rnames <-
  setdiff(all_rnames, long_rnames)

high_wc_rnames <-as.character(strand.states$AWCcontigs@seqnames@values)

low_wc_rnames <- mono_sex_rnames

excluded_unitigs <-
  list(
    data.table(unitig_name = short_rnames, exclusion_reason = paste('Length less than threshold:', as.character(segment_length_threshold))),
    data.table(unitig_name = high_wc_rnames, exclusion_reason = 'Too many WC SSlib'),
    data.table(unitig_name = low_wc_rnames, exclusion_reason = 'Too few WC SSlib')
  ) %>%
  rbindlist(fill=TRUE) %>% # for when there are some categories without any unitigs
  filter(!is.na(unitig_name))
stopifnot(all(excluded_unitigs[, .N, by='unitig_name']$N == 1))
# Unitig Orientation with Bootstrapping -----------------------------------

wfm <-
  bootstraps %>%
  split(by='clust') %>%
  lapply(get.wfrac.matrix)
  # get.wfrac.matrix()
# wfrac.matrix[rn,] #

prcomps <-
  lapply(wfm, prcomp, center=FALSE, scale=FALSE)

# Assume each cluster comes from one chromosome
strand_orientation_clusters <-
  lapply(prcomps, `[[`, 'x') %>%
  # lapply(abs) %>%
  lapply(dist) %>%
  lapply(hclust) %>%
  lapply(cutree, 2) %>%
  lapply(dt_from_named_vector, 'strand_clust', 'rname') %>%
  rbindlist(idcol='clust')
#
#
#
#
#
# #### Check Plotting ####
# pca <-
#   lapply(prcomps, function(pr) {
#     data.table(x = pr$x[, 1], y = pr$x[, 2], rname = rownames(pr$x))
#   }) %>%
#   rbindlist(idcol = 'clust')
#
# ggplot(pca) +
#   geom_point(aes(x=x, y=y)) +
#   facet_wrap(~clust)
#
# roots <-
#   distinct(bootstraps, rname, root_rname)
#
# inversions <-
#   distinct(bootstraps, rname, invert)
#
# plotdata <-
#   strand_orientation_clusters %>%
#   left_join(roots) %>%
#   left_join(inversions) %>%
#   # mutate(rname = unspoof(rname)) %>%
#   left_join(pca, by = c("rname", 'clust')) %>%
#   left_join(ref_aln, by = c('root_rname' = 'rname'))
#
#
# plotdata <-
#   plotdata %>%
#   mutate(strand = ifelse(invert == 0 | is.na(invert), strand, invert_sign(strand)))
#
# plotdata %>%
#   filter(rname == root_rname) %>%
#   # filter(clust == 'sex_LG21 (4)') %>%
#   ggplot(aes(
#     x = x,
#     y = y,
#     fill = as.factor(strand_clust),
#     shape = strand
#   )) +
#   geom_point(size = 5) +
#   facet_wrap( ~ clust) +
#   scale_shape_manual(values = c('+' = 21, '-' = 24))

# # Reorient? ---------------------------------------------------------------
# # debugonce(reorientAndMergeLGs,signature=c('LinkageGroupList','StrandStateMatrix'))
# # getMethod(reorientAndMergeLGs,signature=c('LinkageGroupList','StrandStateMatrix'))
#
# cl2 <-
#   reorientAndMergeLGs(cl,
#                     strand.states$strandMatrix,
#                     cluster = 100,
#                     similarityCutoff = .99 # 1 ~ do not merge clusters?
#                     )
# # fails on component8, component11, component12, component15, component18
#
# cl_dt <-
#   setNames(cl@.Data, cl@names) %>%
#   lapply(function(x) data.table(rname = unspoof(x))) %>%
#   rbindlist(idcol='clust') %>%
#   left_join(ref_aln)
#
# single_orientation_chom <-
#   cl_dt %>%
#   dplyr::count(chrom,strand) %>%
#   dplyr::count(chrom) %>%
#   filter(n==1)
#
#
# # rn <-
# #    intersect(rownames(wfrac.matrix), spoof(rnames_by_gfa$component15))# spoof(rnames_by_gfa$component15))
# wfm <-
#   counts[rname %in% rnames_by_clust$component15] %>%
#   get.wfrac.matrix()
#   # wfrac.matrix[rn,] #
#
# prcmp <-
#   prcomp(wfm, center = FALSE, scale=FALSE)
#
# pca <-
#   data.table(x = prcmp$x[, 1], y = prcmp$x[, 2], rname = unspoof(rownames(prcmp$x)))
#
# hc <-
#   wfm %>%
#   dist() %>%
#   hclust() %>%
#   cutree(2) %>%
#   dt_from_named_vector('clust', 'rname')
#
# contibait_clusters <-
#   setNames(cl2[[3]]@.Data, cl2[[3]]@names) %>%
#   lapply(function(x) data.table(rname=unspoof(x))) %>%
#   rbindlist(idcol='contibait_clust')
#
# contibait_strand <-
#   data.table(
#     rname = cl2[[2]]@seqnames@values,
#     contibait_strand = rep(cl2[[2]]@strand@values, times = cl2[[2]]@strand@lengths)
#   ) %>%
#   left_join(contibait_clusters)
#
# plotdata <-
#   hc %>%
#   mutate(rname = unspoof(rname)) %>%
#   left_join(pca, by = "rname") %>%
#   left_join(ref_aln) %>%
#   left_join(contibait_strand)
#
# with(plotdata, table(contibait_strand, strand, chrom))
#
#
# ggplot(plotdata) +
#   geom_point(aes(x=x, y=y, shape=as.factor(clust), color=contibait_strand))




# Split Counts by Clust ---------------------------------------------------

# Split for SaaRclust
counts.l.all <-
  counts %>%
  split(by = 'lib', keep.by=FALSE) %>%
  lapply(convert_dt_col_to_rownames, 'rname')

counts.l_by_clust <-
  split(bootstraps, by='clust')

counts.l.all_by_clust <-
  lapply(counts.l_by_clust, function(x) {
    x[, .(lib, rname, w, c)] %>%
      split(by = 'lib', keep.by=FALSE) %>%
      lapply(convert_dt_col_to_rownames, 'rname')
  })

# Hclust ------------------------------------------------------------------


## Hclust
hardClust.ord_by_clust <-
  strand_orientation_clusters %>%
  split(by = 'clust') %>%
  lapply(function(x)
    with(x, setNames(strand_clust, rname)))

#
# hardClust.ord_by_clust <-
#   lapply(counts.l_by_clust, function(counts.l) {
#
#     # rnms <- rownames(counts.l)
#     # which_long <- rnms %in% long_rnames
#
#    # counts.l <- counts.l[rname %in% long_rnames, ]
#
#     hardClust(counts.l=counts.l, method='hclust',
#               num.clusters=2,
#               min.cluster.frequency=0,
#               by_chrom=FALSE,
#               chrom.flag=NULL)
#   })


# SaaRclust Input ---------------------------------------------------------

theta.estim_by_clust <-
  lapply(1:length(hardClust.ord_by_clust), function(i) {
    counts.l <- counts.l_by_clust[[i]]
    hardClust.ord <- hardClust.ord_by_clust[[i]]
    estimateTheta(counts.l,
                  hard.clust = hardClust.ord,
                  alpha = alpha,
                  method = 'median')
  }) %>%
  setNames(names(hardClust.ord_by_clust))

readsPerClusts_by_clust <- lapply(hardClust.ord_by_clust, table)
pi.param_by_clust <- lapply(readsPerClusts_by_clust, function(x) x/sum(x))


clust.pairs_by_clust <-
  lapply(theta.estim_by_clust, function(theta.param) {
    data.table(first_clust=1, second_clust=2, chrom_clust=1)
  })


# Soft Clustering ---------------------------------------------------------


soft.clust_by_clust <- lapply(1:length(counts.l.all_by_clust), function(i){
  counts.l.all <- counts.l.all_by_clust[[i]]
  theta.param <- theta.estim_by_clust[[i]]
  pi.param <- pi.param_by_clust[[i]]

  read.names <- rownames(counts.l.all[[1]])
  # clust.pairs <- clust.pairs_by_clust[[i]]

  alpha=0.01
  minLin=NULL
  upperQ=1
  logL.th=1
  theta.constrain=FALSE
  log.scale=TRUE

  SaaRclust(
    counts.l = counts.l.all,
    outputfolder = outputfolder.destination,
    num.clusters = length(pi.param),
    EM.iter = EM.iter,
    alpha = alpha,
    minLib = minLib,
    upperQ = upperQ,
    theta.param = theta.param,
    pi.param = pi.param,
    logL.th = logL.th,
    theta.constrain = theta.constrain,
    log.scale = log.scale,
    read.names = read.names,
    # clust.pairs = clust.pairs,
    numCPU = numCPU
  )

})

names(soft.clust_by_clust) <- names(counts.l.all_by_clust)


# Spoofing, because ML.clust does not appear to be accurate, even though calling WC strand state is?
ML.clust_all <-
  hardClust.ord_by_clust %>%
  lapply(function(x) {
    data.table(first_clust=x, rname=names(x), chrom_clust=1)
  }) %>%
  rbindlist(idcol='clust')

# ML.clust_all <-
#   lapply(soft.clust_by_clust, `[[`, 'ML.clust')  %>%
#   rbindlist(idcol = 'component')


#
# # Extend Clusters that Own Single Components ------------------------------
#
# # TODO add a percentage of single component threshold that has to be met in
# # order to call the whole component for a cluster?
#
# # TODO is this even a good idea? the components that were excluded had
# # properties that made them a typical and possibly low quality?
#
# rnames_by_gfa_df <-
#   rnames_by_gfa %>%
#   lapply(function(x) data.table(rname=x)) %>%
#   rbindlist(idcol='component')
#
# clusters <-
#   right_join(clusters, rnames_by_gfa_df, by='rname')
#
# clust_by_gfa <-
#   clusters %>%
#   filter(!is.na(clust)) %>%
#   distinct(component, clust)
#
# one_cluster_components <-
#   clust_by_gfa %>%
#   group_by(component) %>%
#   filter(n() == 1) %>%
#   ungroup()
#
# # Dictinonary-style vector
# one_cluster_components <-
#   setNames(one_cluster_components$clust,
#            one_cluster_components$component)
#
# clusters <-
#   clusters %>%
#   mutate(clust = ifelse(
#     component %in% names(one_cluster_components),
#     one_cluster_components[component],
#     clust
#   ))
#
# text_data <-
#   clusters %>%
#   right_join(ref_aln) %>%
#   filter(!is.na(component)) %>%
#   group_by(chrom, clust) %>%
#   summarise(weight = sum(qwidth), .groups = 'drop')
#
# ggplot(text_data) +
#   geom_col(aes(x = chrom, y=weight, fill = clust), color='black') +
#   geom_text(aes(x=chrom, y=weight, label=clust), angle=90)
# #
# #
# # # Convert to contibait format
# # wfrac.matrix <-
# #   counts[rname %in% unique(clusters$rname), ] %>%
# #   get.wfrac.matrix() %>%
# #   spoof_rownames()
# #
# # strand.freq <- StrandFreqMatrix(wfrac.matrix)
# #
# #
# # # debugonce(preprocessStrandTable, signature = 'StrandFreqMatrix')
# # strand.states <- preprocessStrandTable(strand.freq,
# #                                        filterThreshold = 1, # 1 ~ don't remove any unitigs
# #                                        orderMethod=FALSE
# # )
# #
# # cl2 <-
# #   clusters %>%
# #   mutate(rname=spoof(rname)) %>%
# #   select(clust, rname) %>%
# #   split(.$clust) %>%
# #   lapply(`[[`, 'rname') %>%
# #   LinkageGroupList(names=names(.))


# Chrom Cluster -----------------------------------------------------------



# renumber to imitate SaarClust Output
strand_clust_to_clust <-
  data.table(clust = names(soft.clust_by_clust))
strand_clust_to_clust[, chrom_clust:= 1:.N]


ML.clust_all <- merge(ML.clust_all[, !"chrom_clust"], strand_clust_to_clust[, .(clust, chrom_clust)])
ML.clust_all[, first_clust:=ifelse(first_clust == 1, 2*chrom_clust-1, 2*chrom_clust)]

# original column order
ML.clust_all <- ML.clust_all[, .(rname, first_clust, chrom_clust)]
# Remove bootstrapped samples
ML.clust_all <-
  semi_join(ML.clust_all, filter(bootstraps, rname==root_rname), by='rname')



# Exclusion Check ---------------------------------------------------------
stopifnot(setequal(all_rnames, c(ML.clust_all$rname, excluded_unitigs$unitig_name)))

# Export ------------------------------------------------------------------
fwrite(excluded_unitigs, exclusion_file, sep='\t', quote=F, row.names=F)

dt.destination <- MLclust.file
orig.colnames <- colnames(ML.clust_all)
colnames(ML.clust_all)[1] <- paste0('#', colnames(ML.clust_all)[1])
fwrite(ML.clust_all, dt.destination, sep='\t', quote=F, row.names=F)
colnames(ML.clust_all) <- orig.colnames

# renumber to imitate SaarClust Output
destination <- clust.pairs.file
clust.pairs <-
  ML.clust_all[, .(first_clust = 2*chrom_clust-1,
                   second_clust = 2*chrom_clust),
               by = chrom_clust]
clust.pairs <-
  clust.pairs[, .(first_clust, second_clust, chrom_clust)]
fwrite(clust.pairs, file=destination, quote = F, row.names = F)


# WC Cell States ----------------------------------------------------------

wc.cells.clusters_all <-
  lapply(1:length(clust.pairs_by_clust), function(i){
    get.wc.cells.clusters(soft.clust_by_clust[[i]], clust.pairs_by_clust[[i]])
  }) %>%
  setNames(names(clust.pairs_by_clust)) %>%
  rbindlist(idcol='clust')

wc.cells.clusters_all <-
  merge(wc.cells.clusters_all, strand_clust_to_clust)

wc.cells.clusters_all[,`:=`(
  clust.forward = ifelse(clust.forward == 1, chrom_clust * 2 - 1, chrom_clust * 2),
  clust.backward = ifelse(clust.backward == 1, chrom_clust * 2 - 1, chrom_clust * 2)
)]
wc.cells.clusters_all = wc.cells.clusters_all[, .(lib, thetawc, clust.forward, clust.backward)]

destination <- wc.cells.file
fwrite(wc.cells.clusters_all, file=destination, quote = F, row.names = F)


# SS Clusters -------------------------------------------------------------

destination <- ss.clust.file

ss.clust <- cluster.ss.reads(alignments, ML.clust_all, clust.pairs)
orig.colnames <- colnames(ss.clust)
colnames(ss.clust)[1] <- paste0('#', colnames(ss.clust)[1])
fwrite(ss.clust, file=destination, sep='\t', quote = F, row.names = F)
dir.name <- dirname(ss.clust.file)

#
# ss.clust.sp <- split(ss.clust, by="chrom_clust")
# lapply(names(ss.clust.sp), function(clust) {
#   component <- with(strand_clust_to_clust, component[chrom_clust==clust])
#   fwrite(
#     ss.clust.sp[[clust]],
#     file = paste0(dir.name, "/ss_clusters/ss_clusters_", component, ".data"),
#     sep = '\t',
#     quote = F,
#     row.names = F
#   )
# }) %>% invisible()
