# Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

# Parsing Helper ----------------------------------------------------------

## Expected Args
# All just single strings?
expected_args <-
  c(
    ## Input
    '--clust-pairs',
    '--wc-cell-clust',
    '--ss-clust',
    '--map',
    '--bubbles',

    ## Output
    '--phased-strand-states',
    '--log',


    ## Wildcards
    '--sample',
    '--clust'
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

# Sourcing ----------------------------------------------------------------

all.sources <-
  c(
    "module_strandphaser/bubble_phasing_lts.R",
    "module_strandphaser/process_mummer_map.R"
  )

# Path Handling
all.sources <- paste0(get_script_dir(), '/', all.sources)

invisible(sapply(all.sources, source))


# Parsing -----------------------------------------------------------------

# sample=snakemake@wildcards[['sample']]
sample <- get_values('--sample')
print(paste('sample=', sample))

# snakemake@wildcards[["clust"]]
clust <- get_values('--clust')
clust.pairs <- fread(get_values('--clust-pairs'))
clust.pairs <- clust.pairs[chrom_clust==clust] # TODO chrom_clust does not exist? data.table syntax maybe?
clusters <- c(clust.pairs$first_clust, clust.pairs$second_clust)
clusters <- as.character(clusters)
wc.cell.clust <- fread(get_values("--wc-cell-clust"))

print('got clusters')

# wc.cell.clust <- fread(snakemake@input[["wc_cell_clust"]])
ss.clust <- fread(get_values('--ss-clust'), header=F)

print(ss.clust)
map_path <- get_values("--map")
print(paste('map:', map_path))
print(paste('bubbles:', get_values("--bubbles")))


map <- fread(map_path)
map <- map[bubbleAllele!="None"]
map.sp <- output_bubble_allele_coverage_matrix(clusters, wc.cell.clust, ss.clust, map)



print('splitted map')

## Get selected library names
select.libs <- wc.cell.clust[clust.forward %in% clusters, unique(lib)]

print('select.libs')
print(select.libs)

#stop('congrats, you've made it this far')
strandphaser(map.sp[[clusters[1]]], map.sp[[clusters[2]]], clusters, select.libs, get_values("--phased-strand-states"))

print ('done')
