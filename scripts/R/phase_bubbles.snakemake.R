# Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
print(args)

#TODO rewrite this to work without ss-clust, only unitig -clust? This would
#likely involve the skipping the output_valid_maps_step and incorporating it
#into this script?

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

log_path <- get_values('--log')
log <- file(log_path, open='wt')
sink(file=log, type='message')
sink(file=log, type='output')


# Library -----------------------------------------------------------------


library(dplyr)
library(tidyr)
library(data.table)

source(file.path(get_script_dir(), "module_utils/utils.R"))

# Parsing -----------------------------------------------------------------



sample <- get_values('--sample')

# Cluster -> Chromosome map
clust.pairs <- fread(get_values('--clust-pairs'))
temp <- clust.pairs %>% 
  mutate(temp = first_clust) %>% 
  mutate(first_clust = second_clust,
         second_clust=temp) %>% 
  select(-temp)
clust.pairs <- rbind(clust.pairs, temp) %>%  arrange(first_clust)


# unitigs clustered by chromosome
unitig_clust <- fread(get_values('--unitig-clust'))
unitig_clust <- merge(unitig_clust, clust.pairs,    
                      by.x = c('first_clust', 'chrom_clust'),
                      by.y = c('first_clust', 'chrom_clust'),
                      all.x = TRUE)

# WC libraries for each chromosome cluster
wc.cell.clust <- fread(get_values("--wc-cell-clust"))
wc.cell.clust <-
  merge(
    wc.cell.clust,
    clust.pairs,
    by.x = c('first_clust', 'second_clust'),
    by.y = c('first_clust', 'second_clust')
  )
print('got clusters')


# Import Exact Matches ----------------------------------------------------

# 
map_files <- get_values("--map", singular=FALSE)

map <-
  # list.files(map_path, full.names = TRUE) 
  map_files %>%
  # setNames(gsub('_valid_maximal_unique_exact_match', '', basename_no_ext(.))) %>%
  lapply(fread) %>%
  rbindlist()



# Join --------------------------------------------------------------------


# join cluster ids to clusters
map <- merge(map, unitig_clust, by.x="unitig_name", by.y='#rname', all.x=TRUE)

# check every unititg only maps to one bubble/allele and chrom.clust
stopifnot(
  map %>%
    distinct(unitig_name, chrom_clust, bubbleName, bubbleAllele) %>%
    dplyr::count(unitig_name) %>%
    with(all(n==1))
)



# Filter to WC libs
map <-
  semi_join(map,
            wc.cell.clust,
            by = c('SSlib' = 'library', 'chrom_clust' = 'chrom_clust'))

# Alignment Direction -----------------------------------------------------

map <-
  map %>% 
  mutate(alignment_direction_cluster = ifelse(isReverseMapped, second_clust, first_clust))
# Coverage ----------------------------------------------------------------


# Count alignments to bubble alleles within each alignment direction cluster
coverage <-
  map[, .(
    num.bubble.al0 = sum(bubbleAllele == '0'),
    num.bubble.al1 = sum(bubbleAllele == '1'),
    num_reads = .N
  ), .(chrom_clust, SSlib, alignment_direction_cluster, bubbleName)]

# within each chrom clust, make both alignment_direction_cluster have the same set of lib and bubbleNames
coverage <-
  coverage %>%
  group_by(chrom_clust) %>%
  tidyr::complete(alignment_direction_cluster, SSlib, bubbleName, fill=list(num.bubble.al0=0,num.bubble.al1=0, num_reads=0)) %>%
  ungroup() %>%
  as.data.table()


# Filtering ---------------------------------------------------------------

# TODO infrastructure to perform filtering is here, but currently no filtering
# is actually performed, given how the threshold values are set.

# number of reads mapped to bubbles by a lib to allow the lib to be used for
# phasing 0 to perform no filtering for now.
lib_coverage_threshold <- 0 # 10

coverage[, bubble_arm_coverage := num.bubble.al0 + num.bubble.al1]
coverage[, al0_ratio := ifelse(bubble_arm_coverage >= lib_coverage_threshold, num.bubble.al0 / bubble_arm_coverage, NA)]

# allele supported by coverage ratio.  set bubbleAllele equal to NA if the ratio
# of allele0 coverage is not significantly low or high (no clear haplotype
# distinction in the bubble/lib)
coverage[, bubbleAllele := ifelse(al0_ratio >= 0.75, '0', ifelse(al0_ratio < 0.25, '1', NA))]


# Filter which bubbles to be used for phasing based on how many libraries map to
# the bubble.
valid_libs_per_bubble <-
  coverage[, .(
    valid_libs = sum(!is.na(bubbleAllele))), by=c('chrom_clust', 'bubbleName')]

# need at least 2 Libs to make comparisons. 2 is a logical minimum to vote on
# concensus scores. Terefore, 2 is essentially not filtering anything but the
# option is availible.
valid_lib_bubble_thresh <- 2

coverage <-
  coverage %>%
  semi_join(
    valid_libs_per_bubble[valid_libs >= valid_lib_bubble_thresh],
    by=c('chrom_clust', 'bubbleName')
  )

# TODO some sort of more specific report about what and how much each is
# filtered?

# All clusters should have reads in both watson and crick orientation. Filter
# out those that don't (small clusters). Sorting breaks if not filtered out.
# TODO see how this affects exclusions?
no_both_ori <-
  coverage %>%
  distinct(chrom_clust, alignment_direction_cluster) %>%
  count(chrom_clust) %>%
  filter(n < 2)

coverage <-
  anti_join(coverage, no_both_ori, by = 'chrom_clust') %>%
  as.data.table()

# Coverage Matrices -------------------------------------------------------
# Phasing is calculated on each chromosome cluster infividually -> Lots of
# function mapping coming up

make_coverage_matrices <- function(coverage) {
  coverage_matrix <-
    data.table::dcast(coverage, alignment_direction_cluster + bubbleName ~ SSlib, value.var = "bubbleAllele")
  
  # alignment_direction_cluster/bubblName combinations without any mappings from a SSlib
  coverage_matrix[is.na(coverage_matrix)] <- "-"
  
  # split map by alignment_direction_cluster
  return(split(coverage_matrix, by = 'alignment_direction_cluster', keep.by = FALSE))
}

coverage_matrices <-
  coverage %>%
  split(by = 'chrom_clust', keep.by = FALSE) %>%
  lapply(make_coverage_matrices)


# StrandphaseR Prep -------------------------------------------------------


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
  
  # Ordering matrices to have the same rows and columns
  coverage_pair <- lapply(coverage_pair, function(x) x[shared_ss_libs, shared_bubble_ids, drop=FALSE])
  
  # Original Matrix ID
  rownames(coverage_pair[[1]]) <-
    paste(rownames(coverage_pair[[1]]),  "C1", sep = "__")
  
  rownames(coverage_pair[[2]]) <-
    paste(rownames(coverage_pair[[2]]),  "C2", sep = "__")
  
  stopifnot(all(Reduce(`==`, lapply(coverage_pair, dim))))
  
  return(coverage_pair)
}

coverage_matrices <-
  lapply(coverage_matrices, strandphaser_prep)

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


# Strand States -----------------------------------------------------------


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
    ) %>% 
    as.data.table()
  
  haplo_strand_states[, `:=`(haplotype=as.integer(haplotype)-1, cluster=as.integer(cluster))]
  
  return(haplo_strand_states)
}


haplo_strand_states <-
  lapply(coverage_matrices, get_strand_states) %>%
  rbindlist(idcol='chrom_clust') %>%
  mutate(chrom_clust = as.integer(chrom_clust))  # for joining


# Exclusion ---------------------------------------------------------------

unitigs_receiving_reads <-
  map[, unique(unitig_name)]

bubble_coverage_summary <-
  valid_libs_per_bubble[, .(
    bubbles_recieving_any_ss_reads_from_wc_sslibs = sum(bubbleName != 'None'),
    bubbles_with_adequate_wc_sslib_coverage = sum(valid_libs >= valid_lib_bubble_thresh)
  ), by = 'chrom_clust']

exclusion_summary <-
  unitig_clust %>%
  dplyr::rename(unitig_name = `#rname`) %>%
  full_join(bubble_coverage_summary, by = c('chrom_clust' = 'chrom_clust')) %>%
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
  dplyr::rename(clust=chrom_clust)

if(any(exclusion_summary$exclusion == 'Failure', na.rm = TRUE)) {
  stop('Exclusion Failure')
}

# Phasing -----------------------------------------------------------------

map <-
  inner_join(
    map,
    haplo_strand_states,
    by = c(
      'SSlib' = 'lib',
      'chrom_clust' = 'chrom_clust',
      'alignment_direction_cluster' = 'cluster'
    )
  )

# Count alignments to each bubble allele from each sorted cluster.
phase_counts <-
  map[, .(n=.N), by=c('unitig_name', 'bubbleName', 'bubbleAllele', 'chrom_clust', 'alignment_direction_cluster', 'haplotype')] %>%
  tidyr::complete(
    nesting(chrom_clust, unitig_name, bubbleName, bubbleAllele),
    alignment_direction_cluster,
    haplotype,
    fill = list(n = 0)
  ) %>%
  as.data.table()

phase_counts <- phase_counts[, .(n=sum(n)), by=c( 'chrom_clust', 'unitig_name', 'bubbleName', 'bubbleAllele', 'haplotype')]


phase_counts <-
  phase_counts %>%
  dcast(chrom_clust + unitig_name + bubbleName + bubbleAllele ~ haplotype, value.var='n') %>%
  dplyr::rename(haplotype_1_support = '0', haplotype_2_support='1')

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



# Export Formatting -------------------------------------------------------
out <-
  out %>%
  arrange(chrom_clust, bubbleName, bubbleAllele) 

out <-
  full_join(out, exclusion_summary,  by = "unitig_name")

out <-
  out %>%
  select(unitig_name, clust, everything()) %>%
  mutate(support_percentage = round(support_percentage * 100, 1))


# Formatting for Rukki? From documentation for input market counts: TSV file,
# where first three columns of every line are interpreted as
# node_name\tmaternal\tpaternal, where 'maternal'/'paternal' are
# parental-specific marker counts. All columns after the third in TSV are
# ignored.
out <-
  out %>% 
  select(unitig_name, haplotype_1_support, haplotype_2_support, everything())

# Export ------------------------------------------------------------------



outpath <- get_values('--output')
if(!dir.exists(dirname(outpath))) {
  dir.create(dirname(outpath), recursive = TRUE)
}
print(outpath)

fwrite(out, outpath, sep= "\t", na='NA')
