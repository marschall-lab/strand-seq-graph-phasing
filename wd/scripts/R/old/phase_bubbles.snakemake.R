# Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
print(args)

# args <- c(
#   "/gpfs/project/projects/medbioinf/projects/mihen108/wd/.snakemake/conda/43452f50fb7d157c7f3e7e0cf48c0de7/lib/R/bin/exec/R",
#   "--slave"                                                                                                                 ,
#   "--no-restore"                                                                                                            ,
#   "--vanilla"                                                                                                               ,
#   "--file=scripts/R/strandphaser.snakemake.R"                                                                               ,
#   "--args"                                                                                                                  ,
#   "--clust-pairs"                                                                                                           ,
#   "HG002/SaaRclust/Clusters/clust_partners.txt"                                                                             ,
#   "--wc-cell-clust"                                                                                                         ,
#   "HG002/SaaRclust/Clusters/wc_cells_clusters.data"                                                                         ,
#   "--ss-clust"                                                                                                              ,
#   "HG002/SaaRclust/Clusters/ss_clusters_component5.data"                                                                    ,
#   "--map"                                                                                                                   ,
#   "HG002/exact_match/valid_5_maximal_unique_exact_match.data"                                                               ,
#   "--bubbles"                                                                                                               ,
#   "HG002/assembly.homopolymer-compressed_graph_components/simplified_component5.gfa.json"                                   ,
#   "--sample"                                                                                                                ,
#   "HG002"                                                                                                                   ,
#   "--clust"                                                                                                                 ,
#   "4"                                                                                                                       ,
#   "--output"                                                                                                  ,
#   "HG002/phased_strand_states/haplo_strand_states_5.data"                                                                   ,
#   "--log"                                                                                                                   ,
#   "log/strandphaser_HG002_5.log"
# )

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
    '--output',
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


# Sourcing ----------------------------------------------------------------
#
# all.sources <-
#   c(
#     "module_strandphaser/bubble_phasing_lts.R",
#     "module_strandphaser/process_mummer_map.R"
#   )
#
# # Path Handling
# # all.sources <- paste0(get_script_dir(), '/', all.sources)
#
# invisible(sapply(all.sources, source))



# Lib ---------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(data.table)
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
ss.clust <- fread(get_values('--ss-clust'), header=T)
stopifnot(length(unique(ss.clust$chrom_clust)) == 1)

print(ss.clust)
map_path <- get_values("--map")
print(paste('map:', map_path))

print(paste('bubbles:', get_values("--bubbles")))


map <- fread(map_path)
map <- map[bubbleAllele!="None"]



# Bubble Coverage Matrices ------------------------------------------------


colnames(ss.clust) <- c("SSname", "SSclust", "chrom.clust")
# ssname and sslib are pasted together, separated by an underscore
ss.clust[, SSname:=strsplit(SSname, '_')[[1]][1], by=SSname]
clust1 <- clusters[1]

map <- merge(map, ss.clust, by="SSname")

coverage <-
  map[, .(
    num.bubble.al0 = sum(bubbleAllele == '0'),
    num.bubble.al1 = sum(bubbleAllele == '1')
  ), .(SSlib, SSclust, bubbleName)]

coverage[, al0_ratio:=num.bubble.al0/(num.bubble.al0+num.bubble.al1)]

# allele supported by coverage ratio.  set bubbleAllele equal to 2 if the ratio
# of allele0 coverage is not significantly low or high (no clear haplotype
# distinction in the bubble/lib)
coverage[, bubbleAllele:=ifelse(al0_ratio >= 0.75, '0', ifelse(al0_ratio < 0.25, '1', '2'))]

# make both SSclust have the same set of lib and bubbleNames
coverage <- tidyr::complete(coverage, SSclust, SSlib, bubbleName) %>% as.data.table()

# select a subset of cells that are wc in this cluster pair
wc.cells <- wc.cell.clust[clust.forward==clust1, lib]
coverage <- coverage[SSlib %in% wc.cells]


map.sp <- data.table::dcast(coverage, SSclust+bubbleName~SSlib, value.var="bubbleAllele")
# SSclust/bubblName combinations without any mappings from a SSlib
map.sp[is.na(map.sp)] <- "-"
# split map by SSclust
map.sp <- split(map.sp, by='SSclust', keep.by=FALSE)




# map.sp <- output_bubble_allele_coverage_matrix(clusters, wc.cell.clust, ss.clust, map)
#stop('congrats, you've made it this far')


# Strandphaser Prep -------------------------------------------------------

print('splitted map')

## Get selected library names
libs_in_wc_state <- wc.cell.clust[clust.forward %in% clusters, unique(lib)]

# strandphaser(map.sp[[clusters[1]]], map.sp[[clusters[2]]], clusters, libs_in_wc_state, get_values("--output"))

print('libs_in_wc_state')
print(libs_in_wc_state)

## recode and convert to matrix class
coverage_matrices <-
  lapply(map.sp, function(x) {
    bubbleNames <- x$bubbleName
    x <- x[, !"bubbleName"]
    x <- mutate_all(x, recode, '-' = 0, '0' = 1, '1' = 2, '2' = 3, .default = 0)
    x <- as.matrix(x)
    rownames(x) <- bubbleNames
    return(t(x))
  })

## Filter and sort shared bubbles ID and cells

bubble_ids <-
  lapply(coverage_matrices, colnames)
shared_bubble_ids <-
  Reduce(intersect, bubble_ids) %>%
  sort()

ss_libs <-
  lapply(coverage_matrices, rownames)
shared_ss_libs <-
  Reduce(intersect, ss_libs)

if (!is.null(libs_in_wc_state)) { # This step is stupid what
  shared_ss_libs <- intersect(shared_ss_libs, libs_in_wc_state)
}

shared_ss_libs <- sort(shared_ss_libs)

coverage_matrices <- lapply(coverage_matrices, function(x) x[shared_ss_libs, shared_bubble_ids, drop=FALSE])

rownames(coverage_matrices[[1]]) <-
  paste(rownames(coverage_matrices[[1]]),  "C1", sep = "__")

rownames(coverage_matrices[[2]]) <-
  paste(rownames(coverage_matrices[[2]]),  "C2", sep = "__")

stopifnot(all(Reduce(`==`, lapply(coverage_matrices, dim))))

# Concensus Row Swapping --------------------------------------------------

consensus.score = data.table()
swap.consensus.score = data.table()

# smaller is better
calc_concensus_margin <- function(x, values = c('1', '2')) {
  values <- as.character(values)
  counts <- table(x)
  counts <- counts[names(counts) %in% values]
  return(sum(counts) - max(counts)) # if no names(counts) %in% values, then warning and function returns Inf
}

calc_matrix_concensus_score <- function(m, positions=1:ncol(m), values =  c('1', '2')) {
  stopifnot(inherits(m, 'matrix'))

  m[, positions, drop=FALSE] %>%
  apply(2, calc_concensus_margin, values=values) %>%
    sum()
}

swap_matrix_rows <- function(l, i) {
  m1 <- l[[1]]
  m2 <- l[[2]]
  stopifnot(inherits(m1, 'matrix'))
  stopifnot(inherits(m2, 'matrix'))

  m1_row_name <- rownames(m1)[i]
  m2_row_name <- rownames(m2)[i]

  m1_row <- m1[i, ]
  m2_row <- m2[i, ]

  m1[i, ] <- m2_row
  rownames(m1)[i] <- m2_row_name

  m2[i, ] <- m1_row
  rownames(m2)[i] <- m1_row_name

  out <- list(m1, m2)
  names(out) <- names(l)

  return(out)
}

for(i in 1:nrow(coverage_matrices[[1]])) {
  # filename <- shared.libs[i]
  # message("Processing ", filename, " ...")
  print(i)
  covered_positions <-
    lapply(coverage_matrices, function(x) {
      which(x[i,] %in% c(0, 1))
    }) %>%
    Reduce(union, .)

    unswapped_concensus_score <-
      vapply(
        coverage_matrices,
        calc_matrix_concensus_score,
        positions = covered_positions,
        FUN.VALUE = double(1) # double instead of integer to allow for Inf
      )

    swapped_concensus_score <-
      vapply(
        swap_matrix_rows(coverage_matrices, i),
        calc_matrix_concensus_score,
        positions = covered_positions,
        FUN.VALUE = double(1)
      )

    if(sum(swapped_concensus_score) < sum(unswapped_concensus_score)) {
      coverage_matrices <- swap_matrix_rows(coverage_matrices, i)
    }

}



# Format Output -----------------------------------------------------------


haplo_strand_states <-
  coverage_matrices %>%
  lapply(function(x) data.table(lib_hap=rownames(x))) %>%
  rbindlist(idcol ='cluster')

haplo_strand_states <-
  tidyr::separate(
    haplo_strand_states,
    col = 'lib_hap',
    into = c('lib', 'haplotype'),
    sep = '__C'
  )

haplo_strand_states[, `:=`(haplotype=as.integer(haplotype)-1, cluster=as.integer(cluster))]


# # Export ------------------------------------------------------------------
# Skip a separate script go right to phasing
# fwrite(haplo_strand_states, file=output.phased.strand.states.file, row.names=F, sep='\t')
#



# Phase Unitigs -----------------------------------------------------------



map <-
  fread(map_path) %>%
  left_join(ss.clust, by = "SSname")

map <- map[SSlib %in% libs_in_wc_state]
map <- map[!is.na(SSclust)]

map <- left_join(map, haplo_strand_states, by=c('SSlib'='lib', 'SSclust'='cluster'))

# check every unititg only maps to one bubble/allele
stopifnot(
  map %>%
    distinct(unitig_name, bubbleName, bubbleAllele) %>%
    dplyr::count(unitig_name) %>%
    with(all(n==1))
  )


phase_counts <-
  map[, .(n=.N), by=c('unitig_name', 'bubbleName', 'bubbleAllele', 'SSclust', 'haplotype')] %>%
  # dplyr::count(unitig_name, bubbleName, bubbleAllele, SSclust, haplotype) %>%
  tidyr::complete(
    nesting(unitig_name, bubbleName, bubbleAllele),
    SSclust,
    haplotype,
    fill = list(n = 0)
  ) %>%
  as.data.table()

phase_counts <- phase_counts[, .(n=sum(n)), by=c('unitig_name', 'bubbleName', 'bubbleAllele', 'haplotype')]

phase_counts <-
  phase_counts %>%
  dcast(unitig_name + bubbleName + bubbleAllele ~ haplotype, value.var='n') %>%
  rename(haplotype_1_support = '0', haplotype_2_support='1')

phase_counts[, support_percentage := haplotype_1_support/(haplotype_1_support+haplotype_2_support)]


# Bubble Phasing Calling -------------------------------------------------
phase_counts <- phase_counts[order(bubbleName, bubbleAllele)]
bubble_phase_counts <- phase_counts[bubbleName != 'None']
bubble_phase_counts[, can_confidently_phase :=
                      ((support_percentage[1] < 0.25 & support_percentage[2] > 0.75) |
                      (support_percentage[2] < 0.25 & support_percentage[1] > 0.75)) & haplotype_1_support+haplotype_2_support > 30,
                    by = c('bubbleName')]
# Non-Bubble Phasing Calling ----------------------------------------------

non_bubble_phase_counts <- phase_counts[bubbleName == 'None']
non_bubble_phase_counts[, can_confidently_phase :=
                          (support_percentage < 0.25 | support_percentage > 0.75) & haplotype_1_support+haplotype_2_support > 30,
                        by = c('unitig_name')]


# Calling -----------------------------------------------------------------


out <- rbindlist(list(bubble_phase_counts,non_bubble_phase_counts))
out[, haplotype:=ifelse(!can_confidently_phase, NA, ifelse(support_percentage > 0.75, 'H1', 'H2'))]

# Export ------------------------------------------------------------------

outpath <- get_values('--output')
if(!dir.exists(dirname(outpath))) {
  dir.create(dirname(outpath), recursive = TRUE)
}
print(outpath)
fwrite(out, outpath, na='NA')

# bubble.allele0.haplo_coverage
# [502, 9]
#
# bubble.allele1.haplo_coverage
# [14, 506]
