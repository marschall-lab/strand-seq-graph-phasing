
# Library -----------------------------------------------------------------

library(argparse)
library(dplyr)
library(purrr)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  needle <- '--file='
  script_path <- gsub(needle, '',  args[grepl(needle, args)])
  return(dirname(script_path))
}

source(file.path(get_script_dir(), "module_utils/utils.R"))
source(file.path(get_script_dir(), "module_utils/phasing_utils.R"))
# Parsing -----------------------------------------------------------------

print_se(commandArgs(trailingOnly = FALSE))

parser <- ArgumentParser(description = 'Fudge Marker Counts for Rukki')
parser$add_argument('--hmc', help = 'haplotype marker count file', required=TRUE, nargs = 1)
parser$add_argument('--output', required=FALSE, nargs = 1)


args <- parser$parse_args()
print_se(args)


# Args --------------------------------------------------------------------


if(is.null(args$output)) {
  out <- stdout()
} else {
  out <- args$output
}

# Import ------------------------------------------------------------------

marker_counts <- readr::read_csv(args$hmc)

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


readr::write_csv(marker_counts, out)
