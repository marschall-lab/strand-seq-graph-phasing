args <- commandArgs(trailingOnly = FALSE)

# Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--bam',
    '--unitigs-bam',

    ## Output
    '--output-prefix',
    '--log',

    ## Params
    '--input-type',
    '--sex',
    '--num-clusters',
    '--num-clusters-soft',
    '--num-alignments',
    '--em-iter',
    '--threads',
    '--cheating',
    '--minSScov',
    '--min-cluster-frequency'
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
  print(values)
  print(idx)
  print(next_idx)
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

# .libPaths(c(.libPaths(), 'utils/R-packages/'))
library(dplyr)
#library(SaaRclust)
library(data.table)
library(Rsamtools)
library(doParallel)
library(matrixStats)
library(assertthat)
library(biovizBase)







# Parsing -----------------------------------------------------------------




stopifnot(all(expected_args %in% args))



## Input

input.alignment.files <- get_values("--bam", singular=FALSE)
ref.aln.bam <- get_values("--unitigs-bam")

inputfolder <- dirname(input.alignment.files[[1]]) # bwa_ss_unitigs/

## Parameters
input_type <- get_values("--input-type")
sex <-        get_values("--sex")
num.clusters <-  as.numeric(get_values("--num-clusters"))
numAlignments <- as.numeric(get_values("--num-alignments"))
EM.iter <-       as.numeric(get_values("--em-iter"))
numCPU <-        as.numeric(get_values("--threads")) # snakemake@threads[[1]]
cheating <-        get_values("--cheating")
minSScov <-        as.numeric(get_values("--minSScov"))
min.cluster.frequency <-        as.numeric(get_values("--min-cluster-frequency"))

## Output
# outputfolder <- dirname(dirname(snakemake@output[["hard_clust"]])) # SaarClust/
outputfolder <- get_values('--output-prefix')

hardclust.file <-   paste0(outputfolder, 'Clusters/hard_clusters.RData') # snakemake@output[["hard_clust"]]
softclust.file <-   paste0(outputfolder, 'Clusters/soft_clusters.RData') # snakemake@output[["soft_clust"]]
MLclust.file <-     paste0(outputfolder, 'Clusters/MLclust.data') # snakemake@output[["ML_clust"]]
ss.clust.file <-    paste0(outputfolder, 'Clusters/ss_clusters.data') # snakemake@output[["ss_clust"]]
clust.pairs.file <- paste0(outputfolder, 'Clusters/clust_partners.txt') # snakemake@output[["clust_pairs"]]
wc.cells.file <-    paste0(outputfolder, 'Clusters/wc_cells_clusters.data') # snakemake@output[["wc_cells_clusters"]]

parsed_args <-
  list(
    inputfolder = inputfolder,
    outputfolder = outputfolder,
    input_type = input_type,
    input.alignment.files = input.alignment.files,
    sex = sex,
    num.clusters = num.clusters,
    EM.iter = EM.iter,
    numAlignments = numAlignments,
    hardclust.file = hardclust.file,
    softclust.file = softclust.file,
    MLclust.file = MLclust.file,
    ss.clust.file = ss.clust.file,
    clust.pairs.file = clust.pairs.file,
    wc.cells.file = wc.cells.file,
    ref.aln.bam = ref.aln.bam,
    numCPU = numCPU,
    cheating=cheating,
    minSScov = minSScov,
    min.cluster.frequency=min.cluster.frequency
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

# Running -----------------------------------------------------------------


# stop("Congragulations, you made it here!")


clust <-
  runSaaRclust(
    inputfolder = inputfolder,
    outputfolder = outputfolder,
    input_type = input_type,
    input.alignment.files = input.alignment.files,
    sex = sex,
    num.clusters = num.clusters,
    EM.iter = EM.iter,
    numAlignments = numAlignments,
    hardclust.file = hardclust.file,
    softclust.file = softclust.file,
    MLclust.file = MLclust.file,
    ss.clust.file = ss.clust.file,
    clust.pairs.file = clust.pairs.file,
    wc.cells.file = wc.cells.file,
    ref.aln.bam = ref.aln.bam,
    numCPU = numCPU,
    cheating=cheating,
    minSScov = minSScov,
    min.cluster.frequency=min.cluster.frequency
  )
