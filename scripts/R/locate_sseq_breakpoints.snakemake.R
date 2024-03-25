
# Functions ---------------------------------------------------------------

pooled_binomial_statistic_wc <- function(wc_signal, half_window_width=100, breakshift= c('leading', 'lagging')) {
  
  breakshift <- match.arg(breakshift)
  # For WC signals specifically. There are "fractional" attempts and successes
  successes <- wc_signal * (wc_signal > 0)
  ns <- abs(wc_signal)
  
  # The current value is included in the leading sum, but not the lagging sum
  lead_sum = slider::slide_dbl(successes, sum, .after = half_window_width - 1) 
  lead_n  = slider::slide_dbl(ns, sum, .after = half_window_width -1)
  lag_sum = slider::slide_dbl(successes, sum, .before = half_window_width, .after = -1)
  lag_n  = slider::slide_dbl(ns, sum, .before = half_window_width, .after = -1)
  
  lead_p = lead_sum/lead_n
  lag_p = lag_sum/lag_n
  lead_lag_p = (lag_n*lag_p + lead_n * lead_p)/(lead_n + lag_n)
  
  # pooled binomial test
  Z = (lag_p - lead_p)/sqrt(lead_lag_p * (1-lead_lag_p) * (1/lag_n + 1/lead_n))
  return(Z)
}

bin_loglik_wc <- 
  function(wc, p) {
    stopifnot(p >= 0, p <= 1)
    n1 <- sum(abs(wc)[wc >0])
    n2 <- sum(abs(wc)[wc <0])
    log_lik <- log(p)*n1  + log(1-p)*n2
    
    return(log_lik)
  }

neighborhood_binomial_likelihood_wc <- function(wc_signal, p, half_window_width=100) {
  slider::slide_dbl(wc_signal, bin_loglik_wc,  p=p, .before = half_window_width, .after = half_window_width )
}


# Command Line ------------------------------------------------------------


## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

# source('/Users/henglinm/Documents/GitHub/strand-seq-graph-phasing/scripts/R/module_utils/phasing_test_args.R')

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    # Input
    '--lib-weights',
    '--haplotype-marker-counts',
    '--sseq-alignments',
    
    # Params
    # '--marker-ratio-threshold',
    # '--peak-pvalue-threshold',
    # '--half-window-width',
    
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
library(magrittr)
library(dplyr)
library(purrr)
library(ggplot2)

# Parsing -----------------------------------------------------------------

# Input
# lib_weights <- "library_weights/NA21487_library_weights.csv"
# lib_weights <- "library_weights/NA19434_library_weights.csv"
# lib_weights <- "library_weights/HG01114_library_weights.csv"
lib_weights <- get_values('--lib-weights')

# raw_counts_df <-  "sseq_alignment_counts/NA21487_sseq_mem_raw.csv" 
# raw_counts_df <-  "sseq_alignment_counts/NA19434_sseq_mem_raw.csv"
# raw_counts_df <-  "sseq_alignment_counts/HG01114_sseq_mem_raw.csv"
raw_counts_df <- get_values('--sseq-alignments')

# haplotype_marker_counts <- "haplotype_marker_counts/NA21487_fudged_haplotype_marker_counts.csv"
# haplotype_marker_counts <- "haplotype_marker_counts/NA19434_fudged_haplotype_marker_counts.csv"
# haplotype_marker_counts <- "haplotype_marker_counts/HG01114_fudged_haplotype_marker_counts.csv"
haplotype_marker_counts <- get_values('--haplotype-marker-counts')

# Params
# ratio_threshold <- as.numeric(get_values('--marker-ratio-threshold'))
ratio_threshold <- 8
# pvalue_threshold <- as.numeric(get_values('--peak-pvalue-threshold'))
pvalue_threshold <- 1e-8
# half_window_width <- as.integer(get_values('--half-window-width'))
half_window_width <- 150

## Output
# output_table <- 'test_output.csv'
output_table <- get_values('--output')


# FIXME A better way of getting sample into things.Mostly used for plotting code to label facets.?
# samp <- strsplit(basename(lib_weights), '_')[[1]][1]


# Import ------------------------------------------------------------------

lib_weights <-
  readr::read_csv(lib_weights)

haplotype_marker_counts <- 
  readr::read_csv(haplotype_marker_counts)

# Import and Filter Counts ------------------------------------------------

# Attempt to filter to HOM-looking unitigs, which may also be those that contain
# suspicious haplotype switches. As long as the hap switches are balanced.
# Filter out unitigs with 0 alignments, as then there is no evidence that they
# could be suspicious?
haplotype_marker_counts <-
  haplotype_marker_counts %>% 
  filter(!(hap_1_counts + hap_2_counts == 0)) %>%
  mutate(
    hap_1_counts = hap_1_counts + 1,
    hap_2_counts = hap_2_counts + 1
  )

haplotype_marker_counts <-
  haplotype_marker_counts %>% 
  mutate(rat = hap_1_counts/hap_2_counts) %>% 
  mutate(pos_rat = exp(abs(log(rat)))) # ~ max(h1, h2)/min(h1, h2)

to_look_at <- 
  haplotype_marker_counts %>% 
  filter(pos_rat <= ratio_threshold) %>%
  select(unitig)

# Lazy loading is faster when filtering most of the data out:
mapq_threshold <- 10
raw_counts_df <-
  readr::read_csv(raw_counts_df, lazy=TRUE) %>% 
  filter(mapq >= mapq_threshold) %>% 
  semi_join(to_look_at, by=c('unitig'))
  

raw_counts_df <-
  raw_counts_df %>%
  dplyr::rename(tstart = pos) %>% 
  mutate(tstart = as.double(tstart))  



# Joining -----------------------------------------------------------------

unitig_info <-
  haplotype_marker_counts %>% 
  distinct(unitig, cluster, length, unitig_orientation)

raw_counts_df <- 
  raw_counts_df %>% 
  left_join(unitig_info, by = c('unitig')) %>% 
  filter(!is.na(cluster))

raw_counts_df <-
  raw_counts_df %>% 
  left_join(lib_weights, by=c('lib', 'cluster')) %>% 
  filter(!is.na(ww_ssf_glm))


# Strand-state specific signal --------------------------------------------

raw_counts_df <- 
  raw_counts_df %>% 
  mutate(
    y = ifelse(strand == '+', 1, -1)
  )

# TODO combine signals from reads with same tstart.
  

# raw_counts_df <-
#   raw_counts_df %>% 
#   group_by(sample, unitig) %>%
#   mutate(x = rank(tstart, ties.method = 'first')) %>% 
#   ungroup()

raw_counts_df <-
  raw_counts_df %>% 
  mutate(
    wc_signal = y * sign(wc_ssf_glm) * wc_ssf_glm^2 * unitig_orientation,
    ww_signal = y * sign(ww_ssf_glm) * ww_ssf_glm^2 * unitig_orientation
  )    

# scale to original counts
raw_counts_df <-
  raw_counts_df %>% 
  group_by(unitig) %>% 
  mutate(
    wc_signal = wc_signal * sum(abs(y)) / sum(abs(wc_signal)),
    ww_signal = ww_signal * sum(abs(y)) / sum(abs(ww_signal))
  ) %>% 
  ungroup()


# Pooled Binomial Convolution ---------------------------------------------

# TODO parallelize  
# Test statistic for binomial difference of proportions

bp_df <-
  raw_counts_df %>% 
  group_by(unitig) %>% 
  arrange(unitig, tstart) %>% # TODO if multiple reads start at the same place, how to handle?
  mutate(
    Z = pooled_binomial_statistic_wc(wc_signal, half_window_width),
    p0.05 = neighborhood_binomial_likelihood_wc(wc_signal, p=0.05, half_window_width),
    p0.50 = neighborhood_binomial_likelihood_wc(wc_signal, p=0.5, half_window_width),
    p0.95 = neighborhood_binomial_likelihood_wc(wc_signal, p=0.95, half_window_width)
  ) %>% 
  mutate(
    pv = pnorm(abs(Z), lower.tail = FALSE),
  ) %>% 
  ungroup()


# Adjust Pvalues ----------------------------------------------------------

# Adjust across whole sample at once.
bp_df <-
  bp_df %>%
  mutate(p_adj = p.adjust(pv, method = 'bonferroni'))

# Handle Edge NAs ---------------------------------------------------------

shift_values_dbl <- function(x, n=0) {
  slider::slide_dbl(x, ~.x, .before=-n, .after=n, .complete = TRUE)
}
bp_df <-
  bp_df %>%
  mutate(pv = coalesce(pv, shift_values_dbl(pv, 1)), 
         p_adj = coalesce(p_adj, shift_values_dbl(p_adj, 1)))


# breakseekR --------------------------------------------------------------

## Peak Suspicious -----------------------------------------------------------


peaks_df <-
  bp_df %>% 
  group_by(unitig) %>% 
  mutate(is_peak_suspicious = p_adj <= pvalue_threshold) %>% 
  ungroup()



## Smooth Peak Suspects ----------------------------------------------------

# TODO some better peak merging/smoothing strategy
# There will often be several small peak windows near the threshold as noise
# lifts the trend above and below the theshold rapdily, while the trend slowly
# crosses the peak.

peaks_df <-
  peaks_df %>%
  group_by(unitig) %>% 
  mutate(is_peak_break = c(FALSE, diff(is_peak_suspicious) != 0)) %>%  
  mutate(break_window = as.factor(cumsum(is_peak_break))) %>% 
  ungroup()


big_peak_width_threshold <- half_window_width/4

peaks_df <-
  peaks_df %>% 
  group_by(unitig, break_window) %>% 
  mutate(peak_size = n()) %>%
  mutate(is_big_peak = is_peak_suspicious & (peak_size >= big_peak_width_threshold)) %>% 
  mutate(is_big_peak = ifelse(is_big_peak, 'big_peak', 'not_big_peak')) %>% 
  ungroup()


# Format Output -----------------------------------------------------------

## Peak Windows ------------------------------------------------------------

# summarise then filter avoids bug when trying to select min or max of 0 length
# variable, which can happen when there are no peaks
peak_windows <-
  peaks_df %>% 
  group_by(unitig, break_window, is_peak_suspicious, is_big_peak) %>% 
  summarise(window_min = min(tstart), window_max = max(tstart), .groups='drop')  %>% 
  filter(is_peak_suspicious)


## Peak Points -------------------------------------------------------------
# 
# 
# peak_points <-
#   peaks_df %>% 
#   filter(is_peak) %>% 
#   group_by(sample, unitig, break_window) %>% 
#   filter(pv == min(pv)) %>% 
#   # select(sample, unitig, break_window, tstart) %>% 
#   ungroup()
# 
# peak_points <-
#   peak_points %>% 
#   group_by(sample, unitig, break_window) %>% 
#   mutate(trank = rank(tstart, ties.method = 'first')) %>% 
#   filter(trank == min(trank) | trank == max(trank)) %>% 
#   ungroup()
# 
# peak_points <-
#   peak_points %>% 
#   group_by(sample, unitig,break_window) %>% 
#   summarise(peak_min = min(tstart), peak_max = max(tstart), .groups='drop')
# 


# Filter w/ Binomial Likelihood -------------------------------------------
unitigs_with_big_state_switches <- 
  peaks_df %>%
  semi_join(peak_windows, by = c('unitig')) %>%
  filter(!is.na(p0.05)) %>% 
  mutate(state = pmap_chr( list(p0.05, p0.50, p0.95), function(x,y,z) c('WW', 'WC', 'CC')[which.max(c(x,y,z))])) %>% 
  count(unitig, state)  

unitigs_with_big_state_switches <-
  unitigs_with_big_state_switches %>% 
  group_by(unitig) %>% 
  mutate(frac = n/sum(n)) %>% 
  filter(frac >= 0.1) %>% 
  filter(all(c('CC', 'WW') %in% state)) %>% 
  distinct(unitig) %>% 
  pull()

peak_windows <-
  peak_windows %>%
  mutate(big_state_switch = ifelse(unitig %in% unitigs_with_big_state_switches, 'hom_state_switch', 'no_hom_state_switch')) %>% 
  select(unitig, break_window, window_min, window_max, everything())
# Output Table ------------------------------------------------------------

peak_windows %>% 
  # select(unitig, window_min, window_max) %>% 
  readr::write_csv(output_table)


# TODO This giant nested block is hideous
if(nrow(peak_windows) > 0) {
  # Plotting ----------------------------------------------------------------
  
  
  ## Orienting Tstart --------------------------------------------------------
  
  # TODO Orienting Tstart can be important if, eg, two arms of a bubble with a
  # switch error are in opposite orientations ~ In this case, then the plots will
  # show the signal switch in the same place after orientation. However, the
  # coordinates in the unitig from a mapping would be broken.
  
  plot_data <-
    peaks_df %>%
    semi_join(peak_windows, by = join_by(unitig)) %>% 
    left_join(distinct(peak_windows, unitig, big_state_switch))
  
  plot_data <-
    plot_data %>%
    group_by(unitig) %>%
    mutate(tstart = ifelse(unitig_orientation == 1, tstart, length -tstart)) %>%
    ungroup()
  
  stopifnot(all(plot_data$tstart >= 0))
  ## Plot --------------------------------------------------------------------
  p <-
    ggplot()

  
  bw <- 50e3
  
  hist_data <-
    plot_data %>%
    filter(wc_signal != 0)
  
  
  # Stupid fucking conda doesn't have ggokabeito availible on the cluster.
  okabeito_palette <-
    c("#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#999999",
      "#000000")
  
  # WC Signal
  p <- 
    p +
    geom_histogram(aes(
      tstart,
      weight = wc_signal,
      group = sign(wc_signal),
      fill = as.factor(sign(wc_signal))
    ),
    binwidth = bw,
    data=hist_data) +
    facet_wrap(~ big_state_switch + unitig, scale = 'free') +
    theme_classic(base_size = 16) +
    scale_fill_manual(values=okabeito_palette[c(1,2,3,7)])
  

  # Have to recalculate to account for inverted Tstarts
  plot_peak_windows <-
    peak_windows %>% 
    left_join(unitig_info) %>% 
    mutate(
      window_min = ifelse(unitig_orientation == 1, window_min, length-window_min),
      window_max = ifelse(unitig_orientation == 1, window_max, length-window_max)
    )
  
  # Peak Windows
  p <-
    p +
    geom_vline(aes(xintercept = window_min),
               alpha = .25,
               data = plot_peak_windows) +
    geom_vline(aes(xintercept = window_max),
               alpha = .25,
               data = plot_peak_windows) +
    geom_rect(
      aes(
        xmin = window_min,
        xmax = window_max,
        ymin = -Inf,
        ymax = Inf,
        fill = is_big_peak
      ),
      alpha = 0.5,
      data = plot_peak_windows
    )
  
  ## Output Plot -------------------------------------------------------------
  
  
  pdf(paste0(output_table, '.pdf'), width = 20)
  print(p)
  dev.off()
  
}

