setwd("~/Documents/wd")

# Library -----------------------------------------------------------------

library(magrittr)
library(dplyr)
library(purrr)
library(pafr)
library(ggplot2)
library(ggokabeito)
library(progressr)
# Import ------------------------------------------------------------------

# sample <- 'HG02769'
# sample <- 'NA19434'
# sample <- 'HG03683'
sample <- c('NA19434')


# library_weights ---------------------------------------------------------


lib_weights <-
  list.files('library_weights/',recursive = TRUE,  full.names = TRUE) %>% 
  set_names(function(x) hutils::trim_common_affixes(basename(x)))

lib_weights <-
  # lib_weights[sample] %>% 
  lib_weights %>% 
  map_dfr(readr::read_csv, .id='sample')


## Haplotype Marker Counts -------------------------------------------------


haplotype_marker_counts <- list.files('haplotype_marker_counts/', full.names = TRUE, pattern = 'fudged_haplotype_marker_counts.csv') %>% 
  set_names(hutils::trim_common_affixes(.))


haplotype_marker_counts <-
  # haplotype_marker_counts[sample] %>% 
  haplotype_marker_counts %>% 
  map(readr::read_csv) %>% 
  bind_rows(.id='sample')


# Size-Ratio Plot ---------------------------------------------------------
plot_data <-
  haplotype_marker_counts %>% 
  filter(!(hap_1_counts + hap_2_counts == 0)) # %>%
# mutate(
#   hap_1_counts = hap_1_counts + 1,
#   hap_2_counts = hap_2_counts + 1
# )


plot_data <-
  plot_data %>% 
  mutate(rat = hap_1_counts/hap_2_counts) %>% 
  mutate(log_rat = log(hap_1_counts/hap_2_counts)) %>% 
  mutate(primary_rat = exp(abs(log_rat))) %>% 
  mutate(abs_log_rat = abs(log_rat),
         log10_length =log10(length))


plot_data %>%
  # filter(log10(length) >= 5) %>%
  ggplot(aes(x = log10(length), abs(log_rat))) +
  geom_point() +
  geom_smooth(method = 'lm')


# To look at: -------------------------------------------------------------

#1.5 ~ 93%
#1 ~ 98%

ratio_threshold <- 6
to_look_at <- 
  plot_data %>% 
  filter(abs(log_rat) <= log(ratio_threshold)) %>%
  select(sample,unitig)



# Raw Counts --------------------------------------------------------------

raw_counts <-
  list.files('sseq_alignment_counts/', full.names = TRUE, pattern='mem_raw.csv') %>%
  # list.files('sseq_alignment_counts/', full.names = TRUE, pattern='fastmap_raw.csv') %>%
  set_names(function(x) hutils::trim_common_affixes(basename(x)))  %>% 
  `[`(sample)

# raw_counts <- raw_counts[sample]
# raw_counts <- raw_counts[1]

n_threads <- 1
if(n_threads > 1) {
  library(furrr)
  plan(multisession, workers=n_threads)
  import_mapper <- furrr::future_map  
} else {
  import_mapper <- purrr::map
}

with_progress({
  p <- progressor(steps = length(raw_counts))

mapq_threshold <- 10
raw_counts_df <-
  raw_counts %>% 
  import_mapper(function(x) {
    out <- readr::read_csv(x, lazy = TRUE)
    out %<>% 
      filter(mapq >= mapq_threshold) %>% 
      semi_join(to_look_at, by=c('sample', 'unitig'))
    p()
    return(out)
    }) %>% 
  bind_rows(.id = 'sample')
})

if(n_threads > 1) {
  # close workers
  plan(sequential)
}

# TODO  for mem
raw_counts_df %<>%
  dplyr::rename(tstart = pos)


raw_counts_df <- 
  raw_counts_df %>% 
  mutate(
    y = ifelse(strand == '+', 1, -1),
    tstart = as.double(tstart)
  )  %>% 
  select(-qname)

raw_counts_df <-
  raw_counts_df %>% 
  group_by(sample, unitig) %>%
  mutate(x = rank(tstart, ties.method = 'first')) %>% 
  ungroup()

SAVE <- raw_counts_df
# 
# # Size-Ratio Plot ---------------------------------------------------------
# plot_data <-
#   haplotype_marker_counts %>% 
#   filter(!(hap_1_counts + hap_2_counts == 0)) # %>%
#   # mutate(
#   #   hap_1_counts = hap_1_counts + 1,
#   #   hap_2_counts = hap_2_counts + 1
#   # )
# 
# 
# plot_data <-
#   plot_data %>% 
#   mutate(rat = hap_1_counts/hap_2_counts) %>% 
#   mutate(log_rat = log(hap_1_counts/hap_2_counts)) %>% 
#   mutate(primary_rat = exp(abs(log_rat))) %>% 
#   mutate(abs_log_rat = abs(log_rat),
#          log10_length =log10(length))
# 
# 
# plot_data %>%
#   # filter(log10(length) >= 5) %>%
#   ggplot(aes(x = log10(length), abs(log_rat))) +
#   geom_point() +
#   geom_smooth(method = 'lm')
# 
# 
# # To look at: -------------------------------------------------------------
# 
# #1.5 ~ 93%
# #1 ~ 98%
# 
# ratio_threshold <- 6
# to_look_at <- 
#   plot_data %>% 
#   filter(abs(log_rat) <= log(ratio_threshold)) %>%
#   select(sample,unitig)
# 
# raw_counts_df <- 
#   raw_counts_df %>% 
#   semi_join(to_look_at, by=c('sample', 'unitig'))
# Joining -----------------------------------------------------------------



unitig_info <-
  haplotype_marker_counts %>% 
  distinct(sample, unitig, cluster, length, unitig_orientation)

raw_counts_df <- 
  raw_counts_df %>% 
  left_join(unitig_info, by = c('sample', 'unitig')) %>% 
  filter(!is.na(cluster))

raw_counts_df <-
  raw_counts_df %>% 
  left_join(lib_weights, by=c('sample', 'lib', 'cluster')) %>% 
  filter(!is.na(ww_weight_mem))



# Weighting Counts --------------------------------------------------------

# FIXME, not max tstart, but actual tlen
raw_counts_df <-
  raw_counts_df %>% 
  group_by(sample, unitig) %>% 
  mutate(tstart = ifelse(unitig_orientation == 1, tstart, max(tstart) -tstart)) %>% 
  ungroup()

raw_counts_df <- 
  raw_counts_df %>% 
  mutate(
    wc_signal = y * sign(wc_weight_fastmap) * wc_weight_fastmap^2 * unitig_orientation,
    # ww_signal = y * sign(ww_weight_fastmap) * ww_weight_fastmap^2 * unitig_orientation
    ww_signal = y * sign(ww_weight_mem) * ww_weight_mem^2 * unitig_orientation
    )    %>% 
  select(-wc_weight_fastmap,
         -ww_weight_fastmap,
         -ww_weight_mem, 
         -unitig_orientation,
         -strand,
         -cluster)

# scale to original counts
raw_counts_df <- 
  raw_counts_df %>% 
  group_by(sample, unitig) %>% 
  mutate(
    wc_signal = wc_signal * sum(abs(y))/ sum(abs(wc_signal)),
    ww_signal = ww_signal * sum(abs(y))/sum(abs(ww_signal))
  ) %>% 
  ungroup()

# raw_counts_df <-
#   raw_counts_df %>% 
#   mutate(
#     wc_signal = round(wc_signal),
#     ww_signal = round(ww_signal)
#   )
# Haplotype Signal Plots --------------------------------------------------

 

bw <- 100e3

raw_counts_df %>%
  filter(unitig == 'utig4-1192') %>% 
  # filter(length >= 1e6) %>%
  filter(mapq >= 1) %>%
  ggplot() +
  geom_histogram(aes(
    tstart,
    weight = (y),
    group = sign(y),
    fill = as.factor(sign(y))
  ),
  binwidth = bw) +
  facet_wrap( ~ sample + unitig + lib, scale = 'free') +
  theme_grey(base_size = 9)


raw_counts_df %>%
  # filter(unitig == 'utig4-1192') %>% 
  # filter(length >= 1e6) %>%
  filter(mapq >= 1) %>%
  ggplot() +
  geom_histogram(aes(
    tstart,
    weight = (wc_signal),
    group = sign(wc_signal),
    # group = mapq > 0,
    fill = as.factor(sign(wc_signal))
    # fill = mapq > 0
  ),
  binwidth = bw,
  position = 'stack') +
  geom_freqpoly(aes(
    tstart,
    weight = (wc_signal),
  ),
  binwidth = bw,
  alpha=1/3) +
  facet_wrap( ~ sample + unitig, scales='free') +
  theme_grey(base_size = 9)

raw_counts_df %>%
  filter(unitig %in% c('utig4-1192', 'utig4-1424')) %>%
  # filter(length >= 1e6) %>%
  filter(mapq >= 1) %>%
  ggplot() +
  geom_histogram(aes(
    tstart,
    weight = ww_signal ,
    group = sign(ww_signal),
    fill = as.factor(sign(ww_signal))
  ),
  binwidth = bw) +
  geom_freqpoly(aes(
    tstart,
    weight = (ww_signal),
  ),
  binwidth = bw,
  alpha=0.33) +
  facet_wrap( ~ lib + unitig, scale = 'free') +
  theme_grey(base_size = 9)


# for mem
raw_counts_df %<>%
  filter(mapq >= 1) 

# breakpointR test --------------------------------------------------------
# 
# 
# test_df <- 
#   tibble(wc_signal = c(1,1,1,1,1,-1,-1,-1,-1,-1))
# 
# test_df %>% 
#   mutate(
#     positive_cumsum = cumsum(wc_signal * (wc_signal > 0)),
#     negative_cumsum = cumsum(-wc_signal * (wc_signal < 0))
#   ) %>% 
#   mutate(
#     positive_window=c(rep(NA, 2), diff(positive_cumsum, lag = 2)),
#     #negative_window=c(na_fill, diff(negative_cumsum, lag = half_window_size)),
#   ) %>% 
#   mutate(
#     delta_wc = c(diff(positive_window, lag = 2), rep(NA, 2))
#   )
# 
# test_df %>%
#   mutate(pv =wc_signal * (wc_signal > 0) ) %>% 
#   mutate(
#     positive_cumsum = cumsum(pv)
#   ) %>% 
#   mutate(p_lag = positive_cumsum - lag(positive_cumsum, 2),
#          p_lead = lead(positive_cumsum,2) - positive_cumsum) %>% 
#   mutate(edge = p_lead - p_lag)

# breakpointR -------------------------------------------------------------
# 
# x <- 1:100
# stats::filter(x, rep(1, 3))
# stats::filter(x, rep(1, 3), sides = 1)
# stats::filter(x, rep(1, 3), sides = 1, method='recursive')
# stats::filter(x, rep(1, 3), sides = 1, circular = TRUE)
# 
# stats::filter(presidents, rep(1, 3))


step_filter_convolution <- function(x, half_window_width = 50, circular=FALSE) {
  step_kernel <- c(rep(-1, half_window_width), 0,  rep(1, half_window_width))
  out <- stats::filter(x, step_kernel, method='convolution', sides=2, circular=circular)
  return(as.vector(out))
}

tringle_filter_convolution <- function(x, half_window_width = 50, circular=FALSE) {
  kernel <-
    c(0:(half_window_width-1), half_window_width, (half_window_width-1):0)
  kernel <- kernel / sum(kernel)
  out <- stats::filter(x, kernel, method='convolution', sides=2, circular=circular)
}

moving_average_convolution <- function(x, half_window_width = 50, circular=FALSE) {
  frac <- 1/(2*half_window_width)
  step_kernel <- c(rep(frac, half_window_width), 0,  rep(frac, half_window_width))
  out <- stats::filter(x, step_kernel, method='convolution', sides=2, circular=circular)
  return(as.vector(out))
}

moving_sum_convolution <- function(x, half_window_width = 50, circular=FALSE) {
  step_kernel <- c(rep(1, half_window_width), 0,  rep(1, half_window_width))
  out <- stats::filter(x, step_kernel, method='convolution', sides=2, circular=circular)
  return(as.vector(out))
}

moving_median <- function(x, half_window_width = 50, circular=FALSE) {
  zoo::rollapplyr(x, width = half_window_width * 2 + 1L, FUN = median, fill=NA)
}

leading_sum <- function(x, lag_size = 50) {
  zoo::rollapplyr(x, width = lag_size, FUN = sum, align='left', fill=NA)
}

lagging_sum <- function(x, lag_size = 50) {
  zoo::rollapplyr(x, width = lag_size, FUN = sum, align='right', fill=NA)
}
  
pooled_binomial_statistic_wc <- function(wc_signal, half_window_width=100) {
  # For WC signals specifically
  lead_sum = leading_sum(wc_signal * (wc_signal > 0),half_window_width)
  lag_sum = lagging_sum(wc_signal * (wc_signal > 0),half_window_width)
  lead_n  = leading_sum(abs(wc_signal), half_window_width)
  lag_n  = lagging_sum(abs(wc_signal), half_window_width)
  
  lead_p = lead_sum/lead_n
  lag_p = lag_sum/lag_n
  lead_lag_p = (lag_n*lag_p + lead_n * lead_p)/(lead_n + lag_n)
  
  # pooled binomial test
  Z = (lag_p - lead_p)/sqrt(lead_lag_p * (1-lead_lag_p) * (1/lag_n + 1/lead_n))
  return(Z)
}

# Test statistic for difference of proportions
half_window_width <- 100
bp_df <-
  raw_counts_df %>% 
  group_by(sample, unitig) %>% 
  filter(n() > 5 * half_window_width) %>% 
  arrange(sample, unitig, x) %>% 
  mutate(
    ma_wc = moving_average_convolution(wc_signal, half_window_width)
  ) %>% 
  mutate(
    delta_wc = step_filter_convolution(wc_signal, half_window_width),
    delta_ma_wc = step_filter_convolution(ma_wc, half_window_width),
    Z =pooled_binomial_statistic_wc(wc_signal, half_window_width)
  ) %>% 
  mutate(
    pv = pnorm(abs(Z), lower.tail = FALSE),
  ) %>% 
  mutate(Z_tri = tringle_filter_convolution(Z, half_window_width)) %>% 
  ungroup()



# Adjust Pvalues ----------------------------------------------------------

# Adjust across whole sample at once.
bp_df %<>%
  mutate(p_adj = p.adjust(pv, method = 'bonferroni'))
# Carry over NAs ----------------------------------------------------------


bp_df %<>%
  group_by(sample, unitig) %>%
  mutate(p_adj = zoo::na.locf(p_adj, na.rm = FALSE)) %>%
  mutate(p_adj = zoo::na.locf(p_adj, fromLast = TRUE)) %>%
  ungroup()

# Plotting ----------------------------------------------------------------

bp_df %>%
  ggplot() +
  geom_line(aes(x = x, y=delta_ma_wc)) +
  facet_wrap(~unitig, scales = 'free_x')

bp_df %>%
  ggplot() +
  geom_line(aes(x = x, y=delta_wc)) +
  facet_wrap(~unitig, scales = 'free_x')

bp_df %>%
  ggplot() +
  geom_line(aes(x = x, y=Z)) +
  facet_wrap(~unitig, scales = 'free_x')
 
bp_df %>%
  ggplot() +
  geom_line(aes(x = x, y=-log10(p_adj))) +
  facet_wrap(~unitig, scales = 'free_x')


bp_df %>%
  group_by(sample, unitig) %>% 
  mutate(l10pv = -log10(p_adj)) %>% 
  mutate(p_rel = l10pv/max(l10pv, na.rm=TRUE)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_histogram(aes(
    tstart,
    y=after_stat(ncount),
    weight = (wc_signal),
    group = sign(wc_signal),
    fill = as.factor(sign(wc_signal))
  ), bins=50) +
  geom_line(aes(x = tstart, y=p_rel, group = unitig, color = p_adj <= 1e-6)) +
  facet_wrap(~sample+unitig, scales = 'free') +
  scale_color_okabe_ito()



## breakseekR --------------------------------------------------------------

## Peak Suspicious -----------------------------------------------------------

pvalue_threshold <- 1e-6
peaks_df <-
  bp_df %>% 
  group_by(sample, unitig) %>% 
  mutate(is_peak_suspicious = p_adj <= pvalue_threshold) %>% 
  ungroup()

peaks_df %>%
  ggplot() +
  geom_line(aes(x = x, y=-log10(p_adj), group=unitig, color=as.factor(is_peak_suspicious))) +
  facet_wrap(~sample+unitig, scales = 'free_x') +
  scale_color_okabe_ito()


## Filter Peak Suspects ----------------------------------------------------

# There will often be several small peak windows near the threshold as noise
# lifts the trend above and below the theshold rapdily, while the trend slowly
# crosses the peak.

# TODO some better peak merging/smoothing strategy

peaks_df <-
  peaks_df %>%
  group_by(sample, unitig) %>% 
  mutate(is_peak_break = c(FALSE, diff(is_peak_suspicious) != 0)) %>%  
  mutate(break_window = as.factor(cumsum(is_peak_break))) %>% 
  ungroup()

peak_width_threshold <- half_window_width/4

peaks_df <-
  peaks_df %>% 
  group_by(sample, unitig, break_window) %>% 
  mutate(peak_size = n()) %>%
  mutate(is_peak = is_peak_suspicious & (peak_size >= peak_width_threshold)) %>% 
  ungroup()


peaks_df %>%
  ggplot() +
  geom_line(aes(x = x, y=-log10(p_adj), group=unitig, color=is_peak)) +
  facet_wrap(~sample+unitig, scales = 'free_x') +
  scale_color_okabe_ito()



# Refine Edge Peaks -------------------------------------------------------

edge_peaks <-
  peaks_df %>% 
  group_by(sample, unitig) %>% 
  filter(any(is_peak_suspicious & (x==max(x) | x==min(x))))







# Changepoint Haplotype Signal Plot ---------------------------------------

peak_points <-
  peaks_df %>% 
  filter(is_peak) %>% 
  group_by(sample, unitig, break_window) %>% 
  filter(pn == max(pn)) %>% 
  select(sample, unitig, break_window, tstart) %>% 
  ungroup()

peak_points <-
  peak_points %>% 
  group_by(sample, unitig, break_window) %>% 
  mutate(trank = rank(tstart, ties.method = 'first')) %>% 
  filter(trank == min(trank) | trank == max(trank)) %>% 
  ungroup()

peak_points <-
  peak_points %>% 
  group_by(sample, unitig,break_window) %>% 
  summarise(peak_min = min(tstart), peak_max = max(tstart), .groups='drop')

peak_windows <-
  peaks_df %>% 
  filter(is_peak) %>% 
  group_by(sample, unitig, break_window) %>% 
  mutate(trank = rank(tstart, ties.method = 'first')) %>% 
  filter(trank == min(trank) | trank == max(trank))

peak_windows <-
  peak_windows %>% 
  group_by(sample, unitig,break_window) %>% 
  summarise(window_min = min(tstart), window_max = max(tstart), .groups='drop')
  

to_plot <-
  peaks_df %>%
  # filter(length >= 5e5) %>%
  semi_join(peak_windows, by=c('sample', 'unitig')) %>% 
  distinct(sample, unitig) %>% 
  sample_n(40)


bw <- 100e3
peaks_df %>%
  semi_join(to_plot, by=c('sample', 'unitig')) %>%
  # filter(length >= 5e5) %>%
  semi_join(peak_windows, by=c('sample', 'unitig')) %>%
  ggplot() +
  geom_rect(aes(
    xmin = window_min,
    xmax = window_max,
    ymin = -Inf,
    ymax = Inf
  ),
  alpha = 0.5,
  data = semi_join(peak_windows, to_plot)) +
  geom_vline(aes(xintercept = peak_min), linetype = 'dashed', alpha = 0.75, data=semi_join(peak_points, to_plot)) +
  geom_histogram(aes(
    tstart,
    weight = ww_signal,
    group = sign(ww_signal),
    fill = as.factor(sign(ww_signal))
  ),
  binwidth = bw) +
  facet_wrap( ~ sample + unitig, scale = 'free') +
  theme_classic(base_size = 16)


to_plot <-
  peaks_df %>%
  # filter(length >= 5e5) %>%
  anti_join(peak_windows, by=c('sample', 'unitig')) %>% 
  distinct(sample, unitig) %>% 
  sample_n(40)

peaks_df %>%
  # filter(length >= 5e5) %>%
  semi_join(to_plot, by=c('sample', 'unitig')) %>%
  ggplot() +
  geom_histogram(aes(
    tstart,
    weight = wc_signal,
    group = sign(wc_signal),
    fill = as.factor(sign(wc_signal))
  ),
  binwidth = bw) +
  facet_wrap( ~ sample + unitig, scale = 'free') +
  theme_classic(base_size = 5)

# Pick Peak Points --------------------------------------------------------




  

