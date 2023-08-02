
# Command Line ------------------------------------------------------------

## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

# source('/Users/henglinm/Documents/GitHub/strand-seq-graph-phasing/scripts/R/module_utils/phasing_test_args.R')

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--hmc',
    '--homology',
    ## Output
    '--output'
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


print(args)

# Library -----------------------------------------------------------------

library(dplyr)
library(purrr)

source(file.path(get_script_dir(), "module_utils/utils.R"))
source(file.path(get_script_dir(), "module_utils/phasing_utils.R"))
# Parsing -----------------------------------------------------------------

## Input
homology <- get_values("--homology", singular=TRUE)
marker_counts <- get_values('--hmc', singular=TRUE)

## Output
output <- get_values('--output')

# Import ------------------------------------------------------------------

homology_df <- readr::read_tsv(homology)

marker_counts <- readr::read_csv(marker_counts)


# Homology Filtering ------------------------------------------------------

## Remove multi cluster bubbles --------------------------------------------

multi_cluster_bubbles <-
  homology_df %>%
  left_join(marker_counts, by='unitig') %>% 
  group_by(bubble) %>% 
  filter(all_are_unique(cluster) & !anyNA(cluster)) %>% 
  ungroup() %>% 
  pull_distinct(bubble)

homology_df <-
  homology_df %>% 
  filter(!(bubble %in% multi_cluster_bubbles))

## Bad XY Bubbles ----------------------------------------------------------

# XY bubbles should only be PAR bubbles, therefore small
bad_xy_bubbles <- 
  homology_df  %>% 
  left_join(marker_counts, by='unitig') %>% 
  filter(cluster == 'LGXY') %>% 
  filter(any(length >= 3e6))


homology_df <-
  homology_df %>% 
  anti_join(bad_xy_bubbles, by='bubble')


# Marker Fudging ----------------------------------------------------------


## Count Inflation ---------------------------------------------------------

# This emerged from discussions with Sergey Koren with regards to rukki. It was
# discovered that very large unitigs can have very low unitig counts while still
# appearing to show a plausible phasing. EG utig4-234[01], utig4-2336,
# utig4-1216 from HG02106 show plausible marker counts despite a very low count
# of markers per unitig (<= 3). This could cause problems with rukki where,
# unless certain filtering thresholds were essentially entirely eliminated,
# those large unitigs would not be considered during path finding. After Sergey
# K. consulted with Sergey N., it was suggested that one simple trick would be
# to multiply the marker counts of large unitigs by some scaling factor. This
# would essentially allow large unitigs to be included during path finding
# without removing all filtering. This allows large unitigs to be retained even
# when smaller unitigs are discarded for low marker counts.

marker_counts <-
  marker_counts %>%
  mutate(
    fudge_hap_1_counts = ifelse(length >= 1e6, 10 * hap_1_counts, hap_1_counts),
    fudge_hap_2_counts = ifelse(length >= 1e6, 10 * hap_2_counts, hap_2_counts)
  )




## Length Ratio Inflation ---------------------------------------------------------

## Some very long acrocentric unitigs have relatively low marker ratios, despite
## separating bubble arms in a logical way. This assumes that no homozygous
## unitig will ever be longer than 1e7 Mbp, and amplifies the marker ratio so
## that rukki will assign the unitig to a haplotype.
# marker_counts <-
#   marker_counts %>% 
#   mutate(
#     fudge_hap_1_counts = ifelse(length < 1e7, fudge_hap_1_counts, ifelse(hap_1_counts > hap_2_counts, 10 * fudge_hap_1_counts, fudge_hap_1_counts)),
#     fudge_hap_2_counts = ifelse(length < 1e7, fudge_hap_2_counts, ifelse(hap_2_counts > hap_1_counts, 10 * fudge_hap_2_counts, fudge_hap_2_counts)),
#   )

long_unitigs_to_inflate <- character(0)
# long_unitigs_to_inflate <-
#   marker_counts %>%
#   filter(length >= 4e7) %>% # HG03683 utig4-166[5,6]
#   mutate(fudge_hap_1_counts = fudge_hap_1_counts + 1,fudge_hap_2_counts = fudge_hap_2_counts + 1 ) %>%
#   # filter((fudge_hap_1_counts / fudge_hap_2_counts > 1.25) | (fudge_hap_2_counts/fudge_hap_1_counts > 1.25)) %>%
#   pull_distinct(unitig)

## Bubble Ratio Inflation --------------------------------------------------

## This uses homology to "rescue" low marker ratio unitigs. If both the marker
## ratios and homology agree, then the marker ratio is inflated for rukki. The
## marginal gain expected from this step is quite small, as the other steps
## currently appear to be working quite well. This may mean that, for the sake
## of maintenance and development, that this step will be abandoned in the
## future, but for now, this is how homology is being used.


# bubble_counts <-
#   marker_counts %>%
#   left_join(homology_df, by='unitig') %>% 
#   filter(!is.na(bubble)) %>% 
#   mutate(wfrac = (hap_1_counts - hap_2_counts)/(hap_1_counts + hap_2_counts))
# 
# count_threshold <- 5
# bubble_counts <-
#   bubble_counts %>% 
#   filter((hap_1_counts + hap_2_counts) > count_threshold) %>% 
#   group_by(bubble) %>% 
#   filter(n() == 2) %>% 
#   ungroup()
# 
# bubble_counts <-
#   bubble_counts %>% 
#   group_by(bubble) %>% 
#   filter((1 %in% sign(wfrac)) & (-1 %in% sign(wfrac))) %>% 
#   ungroup()

bubble_unitigs_to_inflate <- character(0)
# bubble_unitigs_to_inflate <-
#   bubble_counts %>% 
#   pull_distinct(unitig)



## Inflate Ratios ----------------------------------------------------------

unitigs_to_inflate <-
  unique(c(long_unitigs_to_inflate, bubble_unitigs_to_inflate))

marker_counts <-
  marker_counts %>% 
  mutate(
    fudge_hap_1_counts = ifelse(!(unitig %in% unitigs_to_inflate), fudge_hap_1_counts, ifelse(hap_1_counts > hap_2_counts, 10 * fudge_hap_1_counts, fudge_hap_1_counts)),
    fudge_hap_2_counts = ifelse(!(unitig %in% unitigs_to_inflate), fudge_hap_2_counts, ifelse(hap_2_counts > hap_1_counts, 10 * fudge_hap_2_counts, fudge_hap_2_counts)),
  )


# Wfrac -------------------------------------------------------------------
marker_counts <-
  marker_counts %>% 
  mutate(
    fudge_wfrac = (fudge_hap_1_counts - fudge_hap_2_counts)/(fudge_hap_1_counts + fudge_hap_2_counts),
    wfrac = (hap_1_counts - hap_2_counts)/(hap_1_counts + hap_2_counts)
  ) %>% 
  mutate(
    fudge_wfrac = round(fudge_wfrac,2),
    wfrac = round(wfrac,2)
  )


# Export ------------------------------------------------------------------

marker_counts <-
  marker_counts %>% 
  select(unitig, fudge_hap_1_counts, fudge_hap_2_counts, everything())


readr::write_csv(marker_counts, output)
