args <- commandArgs(trailingOnly = FALSE)

# Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--bam',
    '--gfas',
    ## Output
    '--output-prefix',
    '--log',
    
    ## Params
    '--em-iter',
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
gfas <- get_values("--gfas", singular=FALSE)

inputfolder <- dirname(input.alignment.files[[1]]) # bwa_ss_unitigs/

## Parameters
EM.iter <-       as.numeric(get_values("--em-iter"))
numCPU <-        as.numeric(get_values("--threads")) # snakemake@threads[[1]]


## Output
# outputfolder <- dirname(dirname(snakemake@output[["hard_clust"]])) # SaarClust/
outputfolder <- get_values('--output-prefix')

MLclust.file <-     paste0(outputfolder, 'Clusters/MLclust.data') # snakemake@output[["ML_clust"]]
ss.clust.file <-    paste0(outputfolder, 'Clusters/ss_clusters.data') # snakemake@output[["ss_clust"]]
clust.pairs.file <- paste0(outputfolder, 'Clusters/clust_partners.txt') # snakemake@output[["clust_pairs"]]
wc.cells.file <-    paste0(outputfolder, 'Clusters/wc_cells_clusters.data') # snakemake@output[["wc_cells_clusters"]]

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

# Running -----------------------------------------------------------------

#=========================#
### Create directiories ###
#=========================#


read_rnames_from_gfa <- function(x) {
  grep('^S',
       readLines(x),
       value = TRUE) %>%
    strsplit("\t") %>%
    vapply(function(x) x[[2]], character(1)) # grab the second field in every line
}

read_component_lengths_from_gfa <- function(x) {
  grep('^S',
       readLines(x),
       value = TRUE) %>%
    strsplit("\t") %>%
    vapply(function(x) nchar(x[[3]]), integer(1)) # grab the second field in every line
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


#Create a master output directory
outputfolder.destination <- file.path(outputfolder)
dir_create_if_does_not_exist(outputfolder.destination)

#Directory to store raw read counts and best alignments
rawdata.store <- file.path(outputfolder.destination, 'RawData')
dir_create_if_does_not_exist(rawdata.store)

# Directory to store processed/clustered data
Clusters.store <- file.path(outputfolder.destination, 'Clusters')
dir_create_if_does_not_exist(Clusters.store)

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

#=========================#
### SS Count Features ###
#=========================#


cat('counting w/c reads...\n')
bam.files <- input.alignment.files
counts <- count.wc.bam(bam.files, numCPU=numCPU)
counts.l.all <- counts[[1]]
alignments <- counts[[2]]
# Don't want to actually filter out any counts at this step, just using it for  table and data formatting
cat('getting representative alignments...\n')
counts <- get_representative_counts(counts.l=counts.l.all, 
                                    num.alignments = Inf, 
                                    min.ss.cov=0)
counts.l.all <- counts[[1]]
counts.l <- counts[[2]]


# Split by component ------------------------------------------------------

cl <- makeCluster(numCPU)
doParallel::registerDoParallel(cl)
rnames_by_component <-
  foreach (gfa = gfas, .packages = 'dplyr') %dopar% {
    read_rnames_from_gfa(gfa)
  }
stopCluster(cl)

names(rnames_by_component) <- basename_no_ext(gfas)

counts.l.all_by_component <-
  lapply(rnames_by_component, function(rnames) {
    lapply(counts.l.all, function(x) {
      rn <- rownames(x)
      idx <- rn %in% rnames
      out <- x[idx]
      rownames(out) <- rn[idx]
      return(out)
    })
    
  })

counts.l_by_component <-
  lapply(rnames_by_component, function(rnames) {
    counts.l[rname %in% rnames]
  })


## Hclust

hardClust.ord_by_component <- 
  lapply(counts.l_by_component, function(counts.l) {
    hardClust(counts.l=counts.l, method='hclust', 
              num.clusters=2, 
              min.cluster.frequency=0, 
              chrom.flag=NULL)
  })


theta.estim_by_component <-
  lapply(1:length(hardClust.ord_by_component), function(i) {
    counts.l <- counts.l_by_component[[i]]
    hardClust.ord <- hardClust.ord_by_component[[i]]
    estimateTheta(counts.l,
                  hard.clust = hardClust.ord,
                  alpha = alpha,
                  method = 'median')
  }) %>% 
  setNames(names(hardClust.ord_by_component))

readsPerClusts_by_component <- lapply(hardClust.ord_by_component, table)
pi.param_by_component <- lapply(readsPerClusts_by_component, function(x) x/sum(x))


clust.pairs_by_component <-
  lapply(theta.estim_by_component, function(theta.estim) {
    data.table(first_clust=1, second_clust=2, chrom_clust=1)
  })

read.names_by_component <- rnames_by_component#rownames(counts.l.all[[1]])
counts.l.all_by_component <- lapply(counts.l.all_by_component, lapply, as.data.frame)

soft.clust_by_component <- lapply(1:length(counts.l.all_by_component), function(i){
  counts.l.all <- counts.l.all_by_component[[i]]
  theta.param <- theta.estim_by_component[[i]]
  pi.param <- pi.param_by_component[[i]]
  read.names <- read.names_by_component[[i]]
  clust.pairs <- clust.pairs_by_component[[i]]
  
  alpha=0.01
  minLin=NULL
  upperQ=1
  logL.th=1
  theta.constrain=FALSE
  log.scale=TRUE
  
  SaaRclust(counts.l=counts.l.all, outputfolder=outputfolder.destination, 
            num.clusters=length(pi.param), EM.iter=EM.iter, alpha=alpha, 
            minLib=minLib, upperQ=upperQ, theta.param=theta.param, pi.param=pi.param, 
            logL.th=logL.th, theta.constrain=theta.constrain, log.scale=log.scale, 
            read.names=read.names, clust.pairs=clust.pairs,
            numCPU=numCPU)
  
})

names(soft.clust_by_component) <- names(counts.l.all_by_component)


ML.clust_all <-
  lapply(soft.clust_by_component, `[[`, 'ML.clust')  %>% 
  rbindlist(idcol = 'component')

# renumber to imitate SaarClust Output
component_to_chrom <- 
  data.table(component = names(soft.clust_by_component))
component_to_chrom[, chrom_clust:= 1:.N]
# component_to_chrom[,`:=`(first_clust=seq(1, 2*.N, by=2))]
# component_to_chrom[, second_clust:=first_clust+1]


ML.clust_all <- merge(ML.clust_all[, !"chrom_clust"], component_to_chrom[, .(component, chrom_clust)])
ML.clust_all[, first_clust:=ifelse(first_clust == 1, 2*chrom_clust-1, 2*chrom_clust)]

# original column order
ML.clust_all <- ML.clust_all[, .(rname, first_clust, chrom_clust)]


dt.destination <- MLclust.file
orig.colnames <- colnames(ML.clust_all)
colnames(ML.clust_all)[1] <- paste0('#', colnames(ML.clust_all)[1])
fwrite(ML.clust_all, dt.destination, sep='\t', quote=F, row.names=F)
colnames(ML.clust_all) <- orig.colnames

# renumber to imitate SaarClust Output
destination <- clust.pairs.file
clust.pairs <- ML.clust_all[, .(first_clust=min(first_clust), second_clust=max(first_clust)), by=chrom_clust]
clust.pairs <- clust.pairs[, .(first_clust, second_clust, chrom_clust)]
fwrite(clust.pairs, file=destination, quote = F, row.names = F)


wc.cells.clusters_all <- 
  lapply(1:length(clust.pairs_by_component), function(i){
    get.wc.cells.clusters(soft.clust_by_component[[i]], clust.pairs_by_component[[i]])
  }) %>% 
  setNames(names(clust.pairs_by_component)) %>% 
  rbindlist(idcol='component')

wc.cells.clusters_all <-
  merge(wc.cells.clusters_all, component_to_chrom)
wc.cells.clusters_all[,`:=`(
  clust.forward = ifelse(clust.forward == 1, chrom_clust * 2 - 1, chrom_clust * 2),
  clust.backward = ifelse(clust.backward == 1, chrom_clust * 2 - 1, chrom_clust * 2)
)]
wc.cells.clusters_all = wc.cells.clusters_all[, .(lib, thetawc, clust.forward, clust.backward)]

destination <- wc.cells.file
fwrite(wc.cells.clusters_all, file=destination, quote = F, row.names = F)


destination <- ss.clust.file

ss.clust <- cluster.ss.reads(alignments, ML.clust_all, clust.pairs)
orig.colnames <- colnames(ss.clust)
colnames(ss.clust)[1] <- paste0('#', colnames(ss.clust)[1])
ss.clust.sp <- split(ss.clust, by="chrom_clust")
fwrite(ss.clust, file=destination, sep='\t', quote = F, row.names = F)
dir.name <- dirname(ss.clust.file)
lapply(names(ss.clust.sp), function(clust) {
  component <- with(component_to_chrom, component[chrom_clust==clust])
  fwrite(
    ss.clust.sp[[clust]],
    file = paste0(dir.name, "/ss_clusters_", component, ".data"),
    sep = '\t',
    quote = F,
    row.names = F
  )
}) %>% invisible()

