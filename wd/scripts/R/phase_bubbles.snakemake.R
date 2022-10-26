# Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
print(args)

#TODO rewrite this to work without ss-clust, only unitig -clust? This would
#likely involve the skipping the output_valid_maps_step

# args <- c(
# "/opt/conda/lib/R/bin/exec/R"                                                        ,
# "--no-echo"                                                                          ,
# "--no-restore"                                                                       ,
# "--vanilla"                                                                          ,
# "--file=scripts/R/phase_bubbles.snakemake.R"                                         ,
# "--args"                                                                             ,
# "--clust-pairs"                                                                      ,
# "HG002v11nl96slt1e5/SaaRclust/Clusters/clust_partners.txt"                                        ,
# "--wc-cell-clust"                                                                    ,
#  "HG002v11nl96slt1e5/SaaRclust/Clusters/wc_cells_clusters.data"                                    ,
#  "--ss-clust"                                                                         ,
#  "HG002v11nl96slt1e5/SaaRclust/Clusters/ss_clusters.data"                                          ,
# "--unitig-clust"                                                                         ,
# "HG002v11nl96slt1e5/SaaRclust/Clusters/MLclust.data"                                          ,
#  "--map"                                                                              ,
#  "HG002v11nl96slt1e5/valid_exact_match"                                                            ,
#  "--bubbles"                                                                          ,
#  "HG002v11nl96slt1e5/assembly.homopolymer-compressed_graph_components/simplified_assembly.gfa.json",
#  "--sample"                                                                           ,
#  "HG002v11nl96slt1e5"                                                                              ,
#  "--output"                                                                           ,
#  "HG002v11nl96slt1e5/phasing/phased_unitigs.csv"                                                  ,
#  "--log"                                                                              ,
#  "log/phase_unitigs_HG002v11nl96slt1e5.log"
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
    '--unitig-clust',
    '--map',
    '--bubbles',

    ## Output
    '--output',
    '--log',

    ## Wildcards
    '--sample'
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



# Library -----------------------------------------------------------------


library(dplyr)
library(tidyr)
library(data.table)


# Functions ---------------------------------------------------------------

basename_no_ext <- function(x) {
  tools::file_path_sans_ext(basename(x))
}
# Parsing -----------------------------------------------------------------

unitig_clust <- fread(get_values('--unitig-clust'))

# sample=snakemake@wildcards[['sample']]
sample <- get_values('--sample')
print(paste('sample=', sample))

# snakemake@wildcards[["clust"]]
clust.pairs <- fread(get_values('--clust-pairs'))
# clust.pairs <- clust.pairs[chrom_clust==clust] # TODO chrom_clust does not exist? data.table syntax maybe?
# clusters <- c(clust.pairs$first_clust, clust.pairs$second_clust)
# clusters <- as.character(clusters)
wc.cell.clust <- fread(get_values("--wc-cell-clust"))
wc.cell.clust <-
  merge(
    wc.cell.clust,
    clust.pairs,
    by.x = c('clust.forward', 'clust.backward'),
    by.y = c('first_clust', 'second_clust')
  )
print('got clusters')


# wc.cell.clust <- fread(snakemake@input[["wc_cell_clust"]])
ss.clust <- fread(get_values('--ss-clust'), header=T)
colnames(ss.clust) <- c("SSname", "SSclust", "chrom.clust")

#TODO this step is by far the slowest in the whole script. It might take longer
#than the rest of the script combined. It might also appear to be unneeded? ~ any(grepl('_', ss.clust$SSname))
#
print(ss.clust$SSname[1:5])

ss.clust[, SSname:=strsplit(SSname, '_')[[1]][1], by=SSname] # ssname and sslib are pasted together, separated by an underscore?


# Import Map --------------------------------------------------------------


map_path <- get_values("--map")
print(paste('map:', map_path))

print(paste('bubbles:', get_values("--bubbles")))

map <-
  list.files(map_path, full.names = TRUE) %>%
  setNames(basename_no_ext(.)) %>%
  lapply(fread) %>%
  rbindlist()



# Bubble Allele Check -----------------------------------------------------


# map <- fread(map_path)
# map <- map[bubbleAllele!="None"]

# Join --------------------------------------------------------------------
# clust1 <- clusters[1]

map <- merge(map, ss.clust, by="SSname", all.x=TRUE)

# check every unititg only maps to one bubble/allele and chrom.clust
stopifnot(
  map %>%
    distinct(unitig_name, chrom.clust, bubbleName, bubbleAllele) %>%
    dplyr::count(unitig_name) %>%
    with(all(n==1))
)



# Filter to WC libs
map <-
  semi_join(map,
            wc.cell.clust,
            by = c('SSlib' = 'lib', 'chrom.clust' = 'chrom_clust'))


# Filter ------------------------------------------------------------------

coverage <-
  map[, .(
    num.bubble.al0 = sum(bubbleAllele == '0'),
    num.bubble.al1 = sum(bubbleAllele == '1'),
    num_reads = .N
  ), .(chrom.clust, SSlib, SSclust, bubbleName)]

# within each chrom clust, make both SSclust have the same set of lib and bubbleNames

coverage <-
  coverage %>%
  group_by(chrom.clust) %>%
  tidyr::complete(SSclust, SSlib, bubbleName, fill=list(num.bubble.al0=0,num.bubble.al1=0, num_reads=0)) %>%
  ungroup() %>%
  as.data.table()

# number of reads mapped to bubbles by a lib to allow it to be used for phasing
# 0 to perform no filtering for now.
lib_coverage_threshold <- 0 # 10

coverage[, bubble_arm_coverage := num.bubble.al0 + num.bubble.al1]
coverage[, al0_ratio := ifelse(bubble_arm_coverage >= lib_coverage_threshold, num.bubble.al0 / bubble_arm_coverage, NA)]

# allele supported by coverage ratio.  set bubbleAllele equal to 2 if the ratio
# of allele0 coverage is not significantly low or high (no clear haplotype
# distinction in the bubble/lib)
coverage[, bubbleAllele := ifelse(al0_ratio >= 0.75, '0', ifelse(al0_ratio < 0.25, '1', NA))]

valid_libs_per_bubble <-
  coverage[, .(
    valid_libs = sum(!is.na(bubbleAllele))), by=c('chrom.clust', 'bubbleName')]

# need at least 2 Libs to make comparisons. 2 is a logical minimum to vote
# conconcensus scores, again, 2 is essentially not filtering anything but the
# option is availible.
valid_lib_bubble_thresh <- 2

coverage <-
  coverage %>%
  semi_join(
    valid_libs_per_bubble[valid_libs >= valid_lib_bubble_thresh],
    by=c('chrom.clust', 'bubbleName')
  )

# TODO some sort of more specific reports about what and how much each is filtered?

# Coverage Matrices -------------------------------------------------------


make_coverage_matrices <- function(coverage) {
  coverage_matrix <-
    data.table::dcast(coverage, SSclust + bubbleName ~ SSlib, value.var = "bubbleAllele")
  # SSclust/bubblName combinations without any mappings from a SSlib
  coverage_matrix[is.na(coverage_matrix)] <- "-"
  # split map by SSclust
  return(split(coverage_matrix, by = 'SSclust', keep.by = FALSE))
}

coverage_matrices <-
  coverage %>%
  split(by = 'chrom.clust', keep.by = FALSE) %>%
  lapply(make_coverage_matrices)


# Bubble Coverage Matrices ------------------------------------------------
#
#
# make_coverage_matrices <- function(map) {
#   coverage <-
#     map[, .(
#       num.bubble.al0 = sum(bubbleAllele == '0'),
#       num.bubble.al1 = sum(bubbleAllele == '1')
#     ), .(chrom.clust, SSlib, SSclust, bubbleName)]
#
#   coverage[, al0_ratio := num.bubble.al0 / (num.bubble.al0 + num.bubble.al1)]
#
#   # allele supported by coverage ratio.  set bubbleAllele equal to 2 if the ratio
#   # of allele0 coverage is not significantly low or high (no clear haplotype
#   # distinction in the bubble/lib)
#   coverage[, bubbleAllele := ifelse(al0_ratio >= 0.75, '0', ifelse(al0_ratio < 0.25, '1', NA))]
#
#   # make both SSclust have the same set of lib and bubbleNames
#   coverage <-
#     tidyr::complete(coverage, chrom.clust, SSclust, SSlib, bubbleName) %>%
#     as.data.table()
#
#   # select a subset of cells that are wc in this cluster pair
#   # ss_clust_id <- unique(coverage$SSclust)
#   # stopifnot(length(ss_clust_id) == 2)
#   #
#   # wc.cells <-
#   #   wc.cell.clust[clust.forward %in% ss_clust_id, lib] %>%
#   #   unique()
#   #
#   # coverage <- coverage[SSlib %in% wc.cells]
#   # coverage <-
#   #   semi_join(coverage,
#   #             wc.cell.clust,
#   #             by = c('SSlib' = 'lib', 'chrom.clust' = 'chrom_clust'))
#
#   coverage_matrix <-
#     data.table::dcast(coverage, SSclust + bubbleName ~ SSlib, value.var = "bubbleAllele")
#   # SSclust/bubblName combinations without any mappings from a SSlib
#   coverage_matrix[is.na(coverage_matrix)] <- "-"
#   # split map by SSclust
#   return(split(coverage_matrix, by = 'SSclust', keep.by = FALSE))
# }
#
# coverage_matrices <-
#   map %>%
#   split(by='chrom.clust', keep.by=TRUE) %>%
#   lapply(make_coverage_matrices)

#
# coverage <-
#   map[, .(
#     num.bubble.al0 = sum(bubbleAllele == '0'),
#     num.bubble.al1 = sum(bubbleAllele == '1')
#   ), .(SSlib, SSclust, bubbleName)]
#
# coverage[, al0_ratio:=num.bubble.al0/(num.bubble.al0+num.bubble.al1)]
#
# # allele supported by coverage ratio.  set bubbleAllele equal to 2 if the ratio
# # of allele0 coverage is not significantly low or high (no clear haplotype
# # distinction in the bubble/lib)
# coverage[, bubbleAllele:=ifelse(al0_ratio >= 0.75, '0', ifelse(al0_ratio < 0.25, '1', '2'))]
#
# # make both SSclust have the same set of lib and bubbleNames
# coverage <- tidyr::complete(coverage, SSclust, SSlib, bubbleName) %>% as.data.table()
#
# # select a subset of cells that are wc in this cluster pair
# wc.cells <- wc.cell.clust[clust.forward==clust1, lib]
# coverage <- coverage[SSlib %in% wc.cells]
#
#
# map.sp <- data.table::dcast(coverage, SSclust+bubbleName~SSlib, value.var="bubbleAllele")
# # SSclust/bubblName combinations without any mappings from a SSlib
# map.sp[is.na(map.sp)] <- "-"
# # split map by SSclust
# map.sp <- split(map.sp, by='SSclust', keep.by=FALSE)




# map.sp <- output_bubble_allele_coverage_matrix(clusters, wc.cell.clust, ss.clust, map)
#stop('congrats, you've made it this far')


# Strandphaser Prep -------------------------------------------------------


code_matrix_for_strandphasing <- function(x) {
  bubbleNames <- x$bubbleName
  x <- x[, !"bubbleName"]
  x <- mutate_all(x, recode, '-' = 0, '0' = 1, '1' = 2, .default = 0)
  x <- as.matrix(x)
  rownames(x) <- bubbleNames
  return(t(x))
}

strandphaser_prep <- function(coverage_pair) {
  ## recode and convert to matrix class
  coverage_pair <-
    lapply(coverage_pair, code_matrix_for_strandphasing)

  ## Filter and sort shared bubbles ID and cells so that the matrices have matching dimensions
  bubble_ids <-
    lapply(coverage_pair, colnames)

  shared_bubble_ids <-
    Reduce(intersect, bubble_ids) %>%
    sort()

  ss_libs <-
    lapply(coverage_pair, rownames)
  shared_ss_libs <-
    Reduce(intersect, ss_libs) %>%
    sort()
  # Ordering matrices
  coverage_pair <- lapply(coverage_pair, function(x) x[shared_ss_libs, shared_bubble_ids, drop=FALSE])

  rownames(coverage_pair[[1]]) <-
    paste(rownames(coverage_pair[[1]]),  "C1", sep = "__")

  rownames(coverage_pair[[2]]) <-
    paste(rownames(coverage_pair[[2]]),  "C2", sep = "__")

  stopifnot(all(Reduce(`==`, lapply(coverage_pair, dim))))

  return(coverage_pair)
}

coverage_matrices <-
  lapply(coverage_matrices, strandphaser_prep)
#
# print('splitted map')
#
# ## Get selected library names
# libs_in_wc_state <- wc.cell.clust[clust.forward %in% clusters, unique(lib)]
#
# # strandphaser(map.sp[[clusters[1]]], map.sp[[clusters[2]]], clusters, libs_in_wc_state, get_values("--output"))
#
# print('libs_in_wc_state')
# print(libs_in_wc_state)
#
# ## recode and convert to matrix class
# coverage_matrices <-
#   lapply(map.sp, function(x) {
#     bubbleNames <- x$bubbleName
#     x <- x[, !"bubbleName"]
#     x <- mutate_all(x, recode, '-' = 0, '0' = 1, '1' = 2, '2' = 3, .default = 0)
#     x <- as.matrix(x)
#     rownames(x) <- bubbleNames
#     return(t(x))
#   })
#
# ## Filter and sort shared bubbles ID and cells
#
# bubble_ids <-
#   lapply(coverage_matrices, colnames)
# shared_bubble_ids <-
#   Reduce(intersect, bubble_ids) %>%
#   sort()
#
# ss_libs <-
#   lapply(coverage_matrices, rownames)
# shared_ss_libs <-
#   Reduce(intersect, ss_libs)
#
# if (!is.null(libs_in_wc_state)) { # This step is stupid what
#   shared_ss_libs <- intersect(shared_ss_libs, libs_in_wc_state)
# }
#
# shared_ss_libs <- sort(shared_ss_libs)
#
# coverage_matrices <- lapply(coverage_matrices, function(x) x[shared_ss_libs, shared_bubble_ids, drop=FALSE])
#
# rownames(coverage_matrices[[1]]) <-
#   paste(rownames(coverage_matrices[[1]]),  "C1", sep = "__")
#
# rownames(coverage_matrices[[2]]) <-
#   paste(rownames(coverage_matrices[[2]]),  "C2", sep = "__")
#
# stopifnot(all(Reduce(`==`, lapply(coverage_matrices, dim))))

# Concensus Row Swapping --------------------------------------------------


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

sort_matrices <- function(coverage_pair){
  # consensus.score = data.table()
  # swap.consensus.score = data.table()
  for(i in 1:nrow(coverage_pair[[1]])) {
    # filename <- shared.libs[i]
    # message("Processing ", filename, " ...")
    print(i)
    covered_positions <-
      lapply(coverage_pair, function(x) {
        which(x[i,] %in% c(1, 2))
      }) %>%
      Reduce(union, .)

    unswapped_concensus_score <-
      vapply(
        coverage_pair,
        calc_matrix_concensus_score,
        positions = covered_positions,
        FUN.VALUE = double(1) # double instead of integer to allow for Inf
      )

    swapped_concensus_score <-
      vapply(
        swap_matrix_rows(coverage_pair, i),
        calc_matrix_concensus_score,
        positions = covered_positions,
        FUN.VALUE = double(1)
      )

    if(sum(swapped_concensus_score) < sum(unswapped_concensus_score)) {
      coverage_pair <- swap_matrix_rows(coverage_pair, i)
    }

  }

  return(coverage_pair)
}


coverage_matrices <-
  lapply(coverage_matrices, sort_matrices)

# consensus.score = data.table()
# swap.consensus.score = data.table()
#
# # smaller is better
# calc_concensus_margin <- function(x, values = c('1', '2')) {
#   values <- as.character(values)
#   counts <- table(x)
#   counts <- counts[names(counts) %in% values]
#   return(sum(counts) - max(counts)) # if no names(counts) %in% values, then warning and function returns Inf
# }
#
# calc_matrix_concensus_score <- function(m, positions=1:ncol(m), values =  c('1', '2')) {
#   stopifnot(inherits(m, 'matrix'))
#
#   m[, positions, drop=FALSE] %>%
#     apply(2, calc_concensus_margin, values=values) %>%
#     sum()
# }
#
# swap_matrix_rows <- function(l, i) {
#   m1 <- l[[1]]
#   m2 <- l[[2]]
#   stopifnot(inherits(m1, 'matrix'))
#   stopifnot(inherits(m2, 'matrix'))
#
#   m1_row_name <- rownames(m1)[i]
#   m2_row_name <- rownames(m2)[i]
#
#   m1_row <- m1[i, ]
#   m2_row <- m2[i, ]
#
#   m1[i, ] <- m2_row
#   rownames(m1)[i] <- m2_row_name
#
#   m2[i, ] <- m1_row
#   rownames(m2)[i] <- m1_row_name
#
#   out <- list(m1, m2)
#   names(out) <- names(l)
#
#   return(out)
# }
#
# for(i in 1:nrow(coverage_matrices[[1]])) {
#   # filename <- shared.libs[i]
#   # message("Processing ", filename, " ...")
#   print(i)
#   covered_positions <-
#     lapply(coverage_matrices, function(x) {
#       which(x[i,] %in% c(0, 1))
#     }) %>%
#     Reduce(union, .)
#
#   unswapped_concensus_score <-
#     vapply(
#       coverage_matrices,
#       calc_matrix_concensus_score,
#       positions = covered_positions,
#       FUN.VALUE = double(1) # double instead of integer to allow for Inf
#     )
#
#   swapped_concensus_score <-
#     vapply(
#       swap_matrix_rows(coverage_matrices, i),
#       calc_matrix_concensus_score,
#       positions = covered_positions,
#       FUN.VALUE = double(1)
#     )
#
#   if(sum(swapped_concensus_score) < sum(unswapped_concensus_score)) {
#     coverage_matrices <- swap_matrix_rows(coverage_matrices, i)
#   }
#
# }



# Format Output -----------------------------------------------------------


get_strand_states <- function(coverage_pair) {
  haplo_strand_states <-
    coverage_pair %>%
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

  return(haplo_strand_states)
}


haplo_strand_states <-
  lapply(coverage_matrices, get_strand_states) %>%
  rbindlist(idcol='chrom_clust') %>%
  mutate(chrom_clust = as.integer(chrom_clust))  # for joining

#
# haplo_strand_states <-
#   coverage_matrices %>%
#   lapply(function(x) data.table(lib_hap=rownames(x))) %>%
#   rbindlist(idcol ='cluster')
#
# haplo_strand_states <-
#   tidyr::separate(
#     haplo_strand_states,
#     col = 'lib_hap',
#     into = c('lib', 'haplotype'),
#     sep = '__C'
#   )
#
# haplo_strand_states[, `:=`(haplotype=as.integer(haplotype)-1, cluster=as.integer(cluster))]


# # Export ------------------------------------------------------------------
# Skip a separate script go right to phasing
# fwrite(haplo_strand_states, file=output.phased.strand.states.file, row.names=F, sep='\t')
#





# Import ------------------------------------------------------------------

# unitig_clust <- fread(get_values('--unitig-clust'))

# map <-
#   list.files(map_path, full.names = TRUE) %>%
#   setNames(basename_no_ext(.)) %>%
#   lapply(fread) %>%
#   # rbindlist(idcol='file') %>%
#   rbindlist() %>%
#   left_join(ss.clust, by = "SSname")
#
# # map <-
# #   fread(map_path) %>%
# #   left_join(ss.clust, by = "SSname")
#
# # map <- map[SSlib %in% libs_in_wc_state]
# # map <- map[!is.na(SSclust)]
#
# map <-
#   semi_join(map,
#             wc.cell.clust,
#             by = c('SSlib' = 'lib', 'chrom.clust' = 'chrom_clust'))
# map <- map[SSlib %in% libs_in_wc_state]


# Exclusion ---------------------------------------------------------------

unitigs_receiving_reads <-
  map[, unique(unitig_name)]
  
bubble_coverage_summary <-
  valid_libs_per_bubble[, .(
    bubbles_recieving_any_ss_reads_from_wc_sslibs = sum(bubbleName != 'None'),
    bubbles_with_adequate_wc_sslib_coverage = sum(valid_libs >= valid_lib_bubble_thresh)
  ), by = 'chrom.clust']

exclusion_summary <-
  unitig_clust %>%
  rename(unitig_name = `#rname`) %>%
  full_join(bubble_coverage_summary, by = c('chrom_clust' = 'chrom.clust')) %>%
  mutate(
    exclusion = case_when(
      bubbles_with_adequate_wc_sslib_coverage > 0 & unitig_name %in% unitigs_receiving_reads ~ NA_character_,
      bubbles_with_adequate_wc_sslib_coverage > 0 & !(unitig_name %in% unitigs_receiving_reads)~ 'No reads mapped to unitig',
      bubbles_recieving_any_ss_reads_from_wc_sslibs == 0 ~ 'No bubbles in clust or no reads mapped to bubbles in clust',
      bubbles_with_adequate_wc_sslib_coverage==0 ~ 'Inadequately covered bubbles in clust',
      is.na(bubbles_recieving_any_ss_reads_from_wc_sslibs) ~ 'No reads mapped to any unitig in clust',
      TRUE ~ 'Failure'
    )
  ) %>%
  select(unitig_name, chrom_clust, exclusion) %>%
  rename(clust=chrom_clust)

if(any(exclusion_summary$exclusion == 'Failure', na.rm = TRUE)) {
  stop('Exclusion Failure')
}

# unitigs that are in clusters with bubbles, but recieve no reads at all.

# Phasing -----------------------------------------------------------------

map <-
  inner_join(
    map,
    haplo_strand_states,
    by = c(
      'SSlib' = 'lib',
      'chrom.clust' = 'chrom_clust',
      'SSclust' = 'cluster'
    )
  )


phase_counts <-
  map[, .(n=.N), by=c('unitig_name', 'bubbleName', 'bubbleAllele', 'chrom.clust', 'SSclust', 'haplotype')] %>%
  tidyr::complete(
    nesting(chrom.clust, unitig_name, bubbleName, bubbleAllele),
    SSclust,
    haplotype,
    fill = list(n = 0)
  ) %>%
  as.data.table()

phase_counts <- phase_counts[, .(n=sum(n)), by=c( 'chrom.clust', 'unitig_name', 'bubbleName', 'bubbleAllele', 'haplotype')]



#
# phase_counts <-
#   map[, .(n=.N), by=c('unitig_name', 'bubbleName', 'bubbleAllele', 'SSclust', 'haplotype')] %>%
#   # dplyr::count(unitig_name, bubbleName, bubbleAllele, SSclust, haplotype) %>%
#   tidyr::complete(
#     nesting(unitig_name, bubbleName, bubbleAllele),
#     SSclust,
#     haplotype,
#     fill = list(n = 0)
#   ) %>%
#   as.data.table()
#
# phase_counts <- phase_counts[, .(n=sum(n)), by=c('unitig_name', 'bubbleName', 'bubbleAllele', 'haplotype')]

phase_counts <-
  phase_counts %>%
  dcast(chrom.clust + unitig_name + bubbleName + bubbleAllele ~ haplotype, value.var='n') %>%
  rename(haplotype_1_support = '0', haplotype_2_support='1')

phase_counts[, support_percentage := haplotype_1_support/(haplotype_1_support+haplotype_2_support)]


# Bubble Confidence Calling -------------------------------------------------
min_reads_for_confidence <- 20
phase_counts <- phase_counts[order(bubbleName, bubbleAllele)]
bubble_phase_counts <- phase_counts[bubbleName != 'None']
bubble_phase_counts[, can_confidently_phase :=
                      ((support_percentage[1] < 0.25 & support_percentage[2] > 0.75) |
                         (support_percentage[2] < 0.25 & support_percentage[1] > 0.75)) & haplotype_1_support+haplotype_2_support >= min_reads_for_confidence,
                    by = c('bubbleName')]
# Non-Bubble Confidence Calling ----------------------------------------------

non_bubble_phase_counts <- phase_counts[bubbleName == 'None']
non_bubble_phase_counts[, can_confidently_phase :=
                          (support_percentage < 0.25 | support_percentage > 0.75) & haplotype_1_support+haplotype_2_support >= min_reads_for_confidence,
                        by = c('unitig_name')]

# Haplotype Calling -------------------------------------------------------



out <- rbindlist(list(bubble_phase_counts,non_bubble_phase_counts))
out[, haplotype:=ifelse(!can_confidently_phase, NA, ifelse(support_percentage > 0.75, 'H1', 'H2'))]

# Export ------------------------------------------------------------------

out <-
  out %>%
  arrange(chrom.clust, bubbleName, bubbleAllele) # %>%
  #select(-chrom.clust) # chrom clust from exclusion_summary

out <-
full_join(out, exclusion_summary,  by = "unitig_name")

out <-
  out %>%
  select(unitig_name, clust, everything()) %>%
  mutate(support_percentage = round(support_percentage * 100, 1))

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
