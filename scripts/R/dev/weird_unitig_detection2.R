setwd("~/Documents/wd")

# Library -----------------------------------------------------------------


library(dplyr)
library(purrr)
library(pafr)
library(ggplot2)
library(ggokabeito)
library(progressr)
# Import ------------------------------------------------------------------

sample <- 'HG02769'

## Haplotype Marker Counts -------------------------------------------------


haplotype_marker_counts <- list.files('haplotype_marker_counts/', full.names = TRUE, pattern = 'fudged_haplotype_marker_counts.csv') %>% 
  set_names(hutils::trim_common_affixes(.))


haplotype_marker_counts <-
  # haplotype_marker_counts[sample] %>% 
  haplotype_marker_counts %>% 
  map(readr::read_csv) %>% 
  bind_rows(.id='sample')

# library_weights ---------------------------------------------------------


lib_weights <-
  list.files('library_weights/',recursive = TRUE,  full.names = TRUE) %>% 
  set_names(function(x) hutils::trim_common_affixes(basename(x)))

lib_weights <-
  # lib_weights[sample] %>% 
  lib_weights %>% 
  map_dfr(readr::read_csv, .id='sample')


# Raw Counts --------------------------------------------------------------


raw_counts <-
  list.files('sseq_alignment_counts/', full.names = TRUE, pattern='fastmap_raw.csv') %>% 
  set_names(function(x) hutils::trim_common_affixes(basename(x))) # %>% 
  # `[`(sample)

# raw_counts <- raw_counts[sample]
# raw_counts <- raw_counts[1]

n_threads <- 10
if(n_threads > 1) {
  library(furrr)
  plan(multisession, workers=n_threads)
  import_mapper <- furrr::future_map  
} else {
  import_mapper <- purrr::map
}

with_progress({
  p <- progressor(steps = length(raw_counts))

raw_counts_df <-
  raw_counts %>% 
  import_mapper(function(x) {
    out <- readr::read_csv(x)
    p()
    return(out)
    }) %>% 
  bind_rows(.id = 'sample')
})

if(n_threads > 1) {
  # close workers
  plan(sequential)
}


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
# Size-Ratio Plot ---------------------------------------------------------
plot_data <-
  haplotype_marker_counts %>% 
  filter(!(hap_1_counts + hap_2_counts == 0)) %>% 
  mutate(
    hap_1_counts = hap_1_counts + 1,
    hap_2_counts = hap_2_counts + 1
  )


plot_data <-
  plot_data %>% 
  mutate(rat = hap_1_counts/hap_2_counts) %>% 
  mutate(log_rat = log(hap_1_counts/hap_2_counts)) %>% 
  mutate(primary_rat = exp(abs(log_rat))) %>% 
  mutate(abs_log_rat = abs(log_rat),
         log10_length =log10(length))


plot_data %>%
  # filter(unitig %in% c('utig4-417', 'utig4-418', 'utig4-2736', 'utig4-2737')) %>%
  # filter(sample == 'HG02769') %>%
  filter(log10(length) >= 6) %>%
  ggplot(aes(x = log10(length), abs(log_rat))) +
  geom_point() +
  geom_smooth(method = 'lm')


# To look at: -------------------------------------------------------------

#1.5 ~ 93%
#1 ~ 98%

to_look_at <- 
  plot_data %>% 
  filter(abs(log_rat) <= 2) %>% 
  select(sample,unitig)

raw_counts_df <- 
  raw_counts_df %>% 
  semi_join(to_look_at, by=c('sample', 'unitig'))
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

raw_counts_df <-
  raw_counts_df %>% 
  group_by(sample, unitig) %>% 
  mutate(tstart = ifelse(unitig_orientation == 1, tstart, max(tstart) -tstart)) %>% 
  ungroup()

raw_counts_df <- 
  raw_counts_df %>% 
  mutate(
    wc_signal = y * sign(wc_weight_fastmap) * wc_weight_fastmap^2 * unitig_orientation,
    ww_signal = y * sign(ww_weight_fastmap) * ww_weight_fastmap^2 * unitig_orientation
    )  %>% 
  select(-wc_weight_fastmap, -ww_weight_fastmap, -ww_weight_mem, -unitig_orientation, -strand, -cluster)

# scale to original counts
raw_counts_df <- 
  raw_counts_df %>% 
  group_by(sample, unitig) %>% 
  mutate(
    wc_signal = wc_signal * sum(abs(y))/ sum(abs(wc_signal)),
    ww_signal = ww_signal * sum(abs(y))/sum(abs(ww_signal))
  ) %>% 
  ungroup()

# Haplotype Signal Plots --------------------------------------------------


bw <- 100e3
raw_counts_df %>%
  filter(sample == 'HG00733red1') %>% 
  filter(unitig %in% c('utig4-35', 'utig4-373', 'utig4-374', 'utig4-377', 'utig4-386', 'utig4-53' )) %>% 
  filter(length >= 1e6) %>%
  # semi_join(unitigs_to_plot) %>%
  ggplot() +
  geom_histogram(aes(
    tstart,
    weight = ww_signal,
    group = sign(ww_signal),
    fill = as.factor(sign(ww_signal))
  ),
  binwidth = bw) +
  facet_wrap( ~ sample + unitig, scale = 'free') +
  theme_grey(base_size = 14)



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


# Read Aggregation --------------------------------------------------------
# 
# aggregated_counts_df <-
#   raw_counts_df %>% 
#   group_by(sample, unitig, tstart) %>% 
#   summarise(wc_signal = sum(wc_signal), 
#             ww_signal=sum(ww_signal), 
#             .groups = 'drop') 
# 
# aggregated_counts_df <-
#   aggregated_counts_df %>% 
#   group_by(sample, unitig) %>% 
#   mutate(x = dense_rank(tstart)) %>%
#   ungroup()
# 
# reads_per_window <- 1
# aggregated_counts_df<-
#   aggregated_counts_df %>% 
#   arrange(sample, unitig, x) %>% 
#   mutate(cum_count = cumsum(abs(wc_signal))) %>% 
#   mutate(count_mod_diff = c(0, diff(cum_count %% reads_per_window))) %>% 
#   mutate(count_block = cumsum(count_mod_diff < 0))
# 
# 
# aggregated_counts_df <-
#   aggregated_counts_df %>%
#   group_by(sample, unitig, count_block) %>%
#   summarise(
#     wc_signal = sum(wc_signal),
#     ww_signal = sum(ww_signal),
#     tmin = min(tstart),
#     tmax = max(tstart),
#     tstart = mean(tstart),
#     n = n(),
#     .groups = 'drop'
#   )
# 
# aggregated_counts_df <-
#   aggregated_counts_df %>% 
#   group_by(sample, unitig) %>% 
#   arrange(sample, unitig, tstart) %>% 
#   mutate(x = rank(tstart, ties.method='first')) %>% 
#   ungroup()
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
  step_kernel <- c(rep(1, half_window_width), 0,  rep(-1, half_window_width))
  out <- stats::filter(x, step_kernel, method='convolution', sides=2, circular=circular)
  return(as.vector(out))
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


half_window_width <- 100
bp_df <-
  raw_counts_df %>% 
  group_by(sample, unitig) %>% 
  filter(n() > 3 * half_window_width) %>% 
  arrange(sample, unitig, x) %>% 
  mutate(
    ma_wc = moving_average_convolution(wc_signal, 10)
    ) %>% 
  mutate(    
    delta_wc = step_filter_convolution(wc_signal, half_window_width),
    delta_ma_wc = step_filter_convolution(ma_wc, half_window_width)
    )

temp <- 
  bp_df %>% 
  filter(!is.na(delta_wc))%>%
  mutate(fft_delta_wc = fft(delta_wc)) #%>% 
  # filter(!is.na(delta_wc))# %>% 


temp %>%
  filter(sample == 'HG02059') %>% 
  ggplot() +
  geom_line(aes(x = x, y=Mod(fft_delta_wc))) +
  facet_wrap(~unitig, scales = 'free_x')


# 
# half_window_size <- 200
# na_fill <- rep(NA, half_window_width)
# 
# bp_df <-
#   raw_counts_df %>%
#   # aggregated_counts_df %>%
#   group_by(sample, unitig) %>%
#   filter(n() > 3 * half_window_width) %>%
#   arrange(sample, unitig, x) %>%
#   mutate(
#     positive_cumsum = cumsum(wc_signal),
#     # a_cummean = cummean(wc_signal),
#     # negative_cumsum = cumsum(-ww_signal * (ww_signal < 0))
#     ) %>%
#   mutate(
#     positive_window=c(na_fill, diff(positive_cumsum, lag = half_window_width)),
#     # negative_window=c(na_fill, diff(negative_cumsum, lag = half_window_size)),
#   ) %>%
#   mutate(
#     # p_window = positive_window / (positive_window + negative_window)
#   ) %>%
#   mutate(
#     delta_wc = c(diff(positive_window, lag = half_window_width), na_fill),
#     # delta_p =  c(diff(p_window, lag = half_window_size), na_fill)
#   ) %>%
#   ungroup()
# 
# bp_df <-
#   bp_df %>% 
#   select(-positive_cumsum, -negative_cumsum,  -negative_window) %>% 
#   # select(-positive_window) %>% 
#   filter(!is.na(delta_wc))

bp_df %>%
  filter(sample == 'HG02059') %>% 
  ggplot() +
  geom_line(aes(x = x, y=delta_wc/(2*half_window_width))) +
  facet_wrap(~unitig, scales = 'free_x')

 bp_df %>%
  filter(sample == 'HG02059') %>% 
  ggplot() +
  geom_line(aes(x = x, y=delta_ma_wc/(2*half_window_width))) +
  facet_wrap(~unitig, scales = 'free_x')
 
 bp_df %>%
   filter(sample == 'HG02059') %>% 
   ggplot() +
   geom_line(aes(x = x, y=Mod(fft_delta_wc))) +
   facet_wrap(~unitig, scales = 'free_x')
 
 

bp_df %>%
  filter(sample == 'HG02059') %>% 
  ggplot() +
  geom_line(aes(x = x, y=ma_wc)) +
  facet_wrap(~unitig, scales = 'free_x')


## breakseekR --------------------------------------------------------------


# trim per unitiug, calculate SD across whole sample? Or SD per unitig?


## SD Filtering ------------------------------------------------------------

half_quantile <- 0.025
bp_df <-
  bp_df %>%
  group_by(sample, unitig) %>%
  mutate(
    is_extreme = delta_p <= quantile(delta_p, half_quantile, na.rm = TRUE) |
                 delta_p >= quantile(delta_p, (1 - half_quantile), na.rm =
                                   TRUE)
  ) %>%
  mutate(sd = sd(delta_p[!is_extreme], na.rm = TRUE)) %>%
  ungroup()

bp_df %>%
  filter(sample == 'HG02769') %>% 
  ggplot() +
  geom_line(aes(x = x, y=delta_p, group=unitig, color=is_extreme)) +
  facet_wrap(~unitig, scales = 'free_x')


## Peak Suspicious -----------------------------------------------------------

p_threshold <- 1-1e-3
peaks_df <-
  bp_df %>% 
  group_by(sample, unitig) %>% 
  mutate(pn = pnorm(abs(delta_p), mean = 0, sd = sd, log.p = TRUE, lower.tail=TRUE)) %>% 
  mutate(is_peak_suspicious = pn >= log(p_threshold)) %>% 
  ungroup()

peaks_df %>%
  filter(sample == 'HG02769') %>% 
  ggplot() +
  geom_line(aes(x = x, y=delta_p, group=unitig, color=as.factor(is_peak_suspicious))) +
  facet_wrap(~unitig, scales = 'free_x') +
  scale_color_okabe_ito()


## Filter Peak Suspects ----------------------------------------------------

# There will often be several small peak windows near the threshold as noise
# lifts the trend above and below the theshold rapdily, while the trend slowl
# crosses the peak.

peaks_df <-
  peaks_df %>%
  group_by(sample, unitig) %>% 
  mutate(is_peak_break = c(FALSE, diff(is_peak_suspicious) != 0)) %>%  
  mutate(break_window = cumsum(is_peak_break)) %>% 
  ungroup()

peak_width_threshold <- half_window_size/4
peaks_df <-
  peaks_df %>% 
  group_by(sample, unitig, break_window) %>% 
  mutate(peak_size = n()) %>%
  # need the all() due to off by one errors that I am not accounting for.
  mutate(is_peak = is_peak_suspicious & (peak_size >= peak_width_threshold)) %>% 
  ungroup()

peaks_df %>%
  filter(sample == 'HG02769') %>% 
  ggplot() +
  geom_line(aes(x = x, y=delta_p, group=unitig, color=as.factor(is_peak))) +
  facet_wrap(~unitig, scales = 'free_x') +
  scale_color_okabe_ito()


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




  

