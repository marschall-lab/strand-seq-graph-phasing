setwd("~/Documents/wd")
# Functions ---------------------------------------------------------------

calculate_n50 <- function(lengths) {
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  total_length <- sum(sorted_lengths)
  half_length <- total_length / 2
  cumulative_length <- cumsum(sorted_lengths)
  n50_index <- min(which(cumulative_length >= half_length))
  n50 <- sorted_lengths[n50_index]
  return(n50)
}
# Library -----------------------------------------------------------------


library(dplyr)
library(purrr)
library(furrr)
library(progressr)
library(pafr)
library(ggplot2)
library(ggokabeito)

# Import ------------------------------------------------------------------


## Reference Alignments ----------------------------------------------------

reference_alignments <- list.files('reference_alignments/T2Tv11_hg002Yv2_chm13/', full.names = TRUE)



reference_alignments <-
  list.files('reference_alignments/T2Tv11_hg002Yv2_chm13/', full.names = TRUE) %>%
  # reference_alignments %>%
  set_names(hutils::trim_common_affixes(.))


n_threads <- 10
if(n_threads > 1) {
  library(furrr)
  plan(multisession, workers=n_threads)
  import_mapper <- furrr::future_map
} else {
  import_mapper <- purrr::map
}

with_progress({
  p <- progressor(steps = length(reference_alignments))

  reference_alignments <-
    reference_alignments %>%
    import_mapper(
      function(x) {
        p()
        out <- pafr::read_paf(x, tibble = TRUE, include_tags = FALSE)

        # take the first alignment listed each time. I think there is an order to them
        # as output by minimap so this seems okay?
        out <-
          out %>%
          group_by(qname) %>%
          slice_head(n=1) %>%
          dplyr::rename(unitig = qname)

        out <-
          out %>%
          mutate(tname = factor(tname, levels = c(paste0('chr', 1:22), 'chrX', 'chrY', 'chrM')))

        return(out)
      }) %>%
    bind_rows(.id='sample')
})

if(n_threads > 1) {
  # close workers
  plan(sequential)
}

## Haplotype Marker Counts -------------------------------------------------


haplotype_marker_counts <- list.files('haplotype_marker_counts/', full.names = TRUE, pattern = 'fudged_haplotype_marker_counts.csv$') %>%
  set_names(hutils::trim_common_affixes(.))


haplotype_marker_counts <-
  haplotype_marker_counts %>%
  map(readr::read_csv) %>%
  bind_rows(.id='sample')

haplotype_marker_counts <-
  haplotype_marker_counts %>%
  filter(!grepl('good-frac', sample) | grepl('frac-100', sample))
## Connected Components ----------------------------------------------------

ccs <- list.files('gfa/ccs', full.names = TRUE) %>%
  set_names(hutils::trim_common_affixes(.))

ccs <-
  ccs %>%
  map(readr::read_tsv) %>%
  bind_rows(.id='sample')

ccs <-
  ccs %>% mutate(component = ifelse(member_largest_component, 0, component))


## Rukki Paths -------------------------------------------------------------------


rukki_paths <-
  list.files('rukki',recursive = TRUE,  full.names = TRUE, pattern = '*rukki_paths.tsv') %>%
  purrr::set_names(function(x) hutils::trim_common_affixes(basename(x)))

rukki_paths <-
  rukki_paths %>%
  map_dfr(readr::read_tsv, .id='sample')

rukki_paths <-
  rukki_paths %>%
  tidyr::separate_longer_delim(path, delim = ',') %>%
  group_by(sample, name) %>%
  mutate(order = 1:n()) %>%
  ungroup() %>%
  dplyr::rename(unitig = path) %>%
  mutate(unitig = gsub("[+-]+$", '', unitig))


rukki_paths <-
  rukki_paths %>%
  mutate(is_gap = grepl('\\[', unitig)) %>%
  group_by(sample, name) %>%
  mutate(run = cumsum(is_gap)) %>%
  ungroup()


rukki_paths <-
  rukki_paths %>%
  mutate(gap_size = ifelse(is_gap, stringr::str_extract(unitig, '[0-9]+'), NA)) %>%
  mutate(gap_size = as.integer(gap_size)) %>%
  left_join(reference_alignments, by = c('sample', 'unitig')) %>%
  mutate(node_len = coalesce(qlen, gap_size)) %>%
  select(sample, name, unitig, assignment, order, is_gap, run, node_len)



rukki_paths <-
  rukki_paths %>%
  mutate(is_not_gap = !is_gap)

## Rukki Annotations -------------------------------------------------------

rukki_annotations <-
  list.files('rukki',recursive = TRUE,  full.names = TRUE, pattern = '*out_final_ann.csv') %>%
  set_names(function(x) hutils::trim_common_affixes(basename(x)))

rukki_annotations <-
  rukki_annotations %>%
  map_dfr(readr::read_tsv, .id='sample')

rukki_annotations <-
  rukki_annotations %>%
  select(sample, node, assignment, info) %>%
  dplyr::rename(unitig = node, annotation = assignment)
# ccn50 -------------------------------------------------------------------
ccn502 <-
  reference_alignments %>%
  left_join(ccs, by=c('sample', 'unitig')) %>%
  group_by(sample, component) %>%
  summarise(qlen = sum(qlen), .groups='drop') %>%
  group_by(sample) %>%
  summarise(ccn502 = calculate_n50(qlen/2)) %>%
  ungroup()

xxx <-
  ccn502 %>%
  filter(sample == 'NA24385') %>%
  pull(ccn502)

ccn502 <-
  ccn502 %>%
  mutate(ccn502 = ifelse(grepl('NA24385_', sample), xxx,ccn502))

# Join --------------------------------------------------------------------

left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
  # filter(!grepl('^Length less than', exclusion)) %>%
  split(.$sample) %>%
  iwalk(function(x, nm) {
    x %>%
      select(unitig, everything()) %>%
      readr::write_csv(glue::glue('rmc/{nm}_ref_marker_counts.csv'))
  })


# Regexes -----------------------------------------------------------------

pilot_regex <- 'HG00733(?!r)|HG02953|HG00171|NA21487|HG02666|HG00358|NA19983'

new_samples <-
  c(
    'HG00096',
    'HG00512',
    'HG00733R1',
    'HG00733R2',
    'HG01596',
    'HG01890',
    'HG02011',
    'HG02492',
    'HG03371',
    'HG03456',
    'HG03732',
    'NA18534',
    'NA19650',
    'NA19705',
    'NA20509'
  )

# Qlen Historgrams --------------------------------------------------------

reference_alignments %>%
  semi_join(haplotype_marker_counts, by='sample') %>%
  ungroup() %>%
  filter(!grepl('^NA24385-', sample)) %>%
  # filter(grepl(pilot_regex, sample, perl = TRUE)) %>%
  ggplot() +
  geom_vline(xintercept = log10(250000), linetype='dashed', alpha=0.5) +
  geom_histogram(aes(log10(qlen))) +
  facet_wrap(~sample, ncol = 5, scale='free_y') +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  xlab('Log10 Unitig Size (verkko ~ hpc, hifiasm ~ not hpc)') +
  theme_bw(base_size = 9) +
  theme(axis.title.y = element_blank())

# # Contig Histogram
# reference_alignments %>%
#   filter(!grepl('^NA24385_', sample)) %>%
#   left_join(ccs, by=c('sample', 'unitig')) %>%
#   group_by(sample, component) %>%
#   summarise(qlen = sum(qlen)) %>%
#   ggplot() +
#   geom_vline(xintercept = log10(250000), linetype='dashed', alpha=0.7) +
#   geom_histogram(aes(log10(qlen)), bins=40) +
#   facet_wrap(~sample, scales='free_y') +
#   scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
#   xlab('Log10 Connected Component Size (hpc)') +
#   theme_bw()


# Marker Ratio Histograms -------------------------------------------------
#
# # Titration improvent plot
# plot_data <-
#   left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
#   mutate(marker_ratio = log2((1+hap_1_counts)/(1+hap_2_counts)))
#
# plot_data %>%
#   filter(grepl('^NA24385_', sample)) %>%
#   filter(!(hap_1_counts==0 & hap_2_counts == 0)) %>%
#   mutate(fraction = case_when(
#     grepl('frac-0', sample) ~ 0,
#     grepl('frac-25', sample) ~ 25,
#     grepl('frac-50', sample) ~ 50,
#     grepl('frac-75', sample) ~ 75,
#     grepl('frac-100', sample) ~ 100
#   )) %>%
#   mutate(rep = case_when(
#     grepl('rep-0', sample) ~ 1,
#     grepl('rep-1', sample) ~ 1,
#     grepl('rep-2', sample) ~ 2,
#     grepl('rep-3', sample) ~ 3,
#     grepl('rep-4', sample) ~ 4
#   )) %>%
#   filter(rep == 1) %>%
#   ggplot() +
#   geom_histogram(aes(marker_ratio)) +
#   scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
#   facet_grid( cols=vars(fraction)) +
#   xlab('Log2 Marker Ratio') +
#   theme_bw(base_size=18) +
#   theme(axis.title.y = element_blank())
# 2D chromosome cluster plot ----------------------------------------------

# plot_data <-
#   full_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
#   group_by(sample, cluster, tname, tlen, attempted) %>%
#   summarise(qlen = sum(qlen)) %>%
#   ungroup()
#

# Simple Cluster-Chrom Size -----------------------------------------------
plot_data <-
  full_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
  group_by(sample, cluster, tname) %>%
  summarise(qlen = sum(qlen), n=n(),  .groups='drop')

plot_data <-
  plot_data %>%
  group_by(sample, tname) %>%
  summarise(`Total` = sum(qlen),
            `Largest Cluster` = max(qlen),
            `Below Threshold` = sum(qlen[is.na(cluster)]), .groups='drop')

plot_data <-
  plot_data %>%
  tidyr::pivot_longer(-(sample:tname), names_to = 'measure', values_to = 'bp') %>%
  mutate(Mbp = bp / 1e6)

plot_data <-  plot_data %>%
  filter(tname != 'chrM') %>%
  filter(!is.na(tname)) %>%
  mutate(tname = gsub('chr', '', tname)) %>%
  mutate(tname = factor(tname, levels = c(1:22, 'X', 'Y')))

plot_data %>%
  filter(!grepl('^NA24385_', sample)) %>%
  # filter(grepl(pilot_regex, sample, perl = TRUE)) %>%
  ggplot() +
  # geom_point(aes(x = tname, y = Mbp, color=measure), size=5.5, shape='-') +
  geom_col(aes(x = tname, y=Mbp, fill=measure), position = 'identity') +
  facet_wrap(~sample, ncol=5) +
  scale_fill_okabe_ito() +
  ylab('Mbp (hpc)') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


# # Chromosome Correctness --------------------------------------------------
# plot_data <-
#   left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
#   group_by(sample, cluster, tname, tlen) %>%
#   summarise(qlen = sum(qlen), n=n(),  .groups='drop')
#
# plot_data <-
#   plot_data %>%
#   # filter(grepl('^NA24385_', sample))
#   filter(grepl(pilot_regex, sample, perl=TRUE))
#
#
# colors <-
#   plot_data %>%
#   filter(!is.na(cluster)) %>%
#   group_by(sample, tname) %>%
#   arrange(desc(qlen)) %>%
#   mutate(color_cluster = LETTERS[1:n()]) %>%
#   ungroup()
# #
# # # some groups of size 1 are just an entire haplotype...
# # colors <-
# #   colors %>%
# #   group_by(sample, cluster) %>%
# #   mutate(color_cluster = ifelse(n == 1, NA, color_cluster)) %>%
# #   ungroup()
#
# plot_data <-
#   plot_data %>%
#   left_join(colors) %>%
#   # mutate(attempted = as.factor(attempted)) %>%
#   mutate(color_cluster = as.factor(color_cluster)) %>%
#   filter(tname != 'chrM') %>%
#   filter(!is.na(tname)) %>%
#   mutate(tname = gsub('chr', '', tname)) %>%
#   mutate(tname = factor(tname, levels = c(1:22, 'X', 'Y')))
#
# plot_data %>%
#   ggplot() +
#   geom_col(
#     aes(
#       x = tname,
#       y = qlen/1e6,
#       fill = color_cluster
#     ),
#     color = 'black',
#     linewidth = 0.1,
#     show.legend = FALSE
#   ) +
#   facet_wrap(~ sample, ncol=5) +
#   scale_fill_viridis_d(
#     option = 'D',
#     direction = -1,
#     na.value = "grey50",
#     begin = 0.1,
#     end = 0.7
#   ) +
#   ylab('Mbp (hpc)') +
#   theme_bw() +
#   theme(axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),)

#
# plot_data %>%
#   group_by(sample, tname) %>%
#   mutate(qlen = qlen / sum(qlen)) %>%
#   ungroup() %>%
#   ggplot() +
#   geom_col(
#     aes(
#       x = tname,
#       y = qlen,
#       fill = color_cluster
#     ),
#     color = 'black',
#     linewidth = 0.1,
#     show.legend = FALSE
#   ) +
#   facet_wrap(~ sample, ncol=5) +
#   scale_fill_viridis_d(
#     option = 'D',
#     direction = -1,
#     na.value = "grey50",
#     begin = 0.1,
#     end = 0.7
#   ) +
#   ylab('Mbp (hpc)') +
#   theme_bw() +
#   theme(axis.title.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),)

# Perc Chrom Correct ------------------------------------------------------
#
# plot_data <-
#   left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
#   group_by(sample, cluster, tname, tlen, attempted) %>%
#   summarise(qlen = sum(qlen),  .groups='drop') %>%
#   group_by(sample, tname) %>%
#   mutate(perc = qlen/sum(qlen))
#
# # plot_data <-
# #   plot_data %>%
# #   filter(!grepl('^NA24385_', sample))
#
#
# colors <-
#   plot_data %>%
#   filter(!is.na(cluster)) %>%
#   group_by(sample, tname) %>%
#   arrange(desc(qlen)) %>%
#   mutate(color_cluster = LETTERS[1:n()])
#
# plot_data %>%
#   left_join(colors) %>%
#   mutate(attempted = as.factor(attempted)) %>%
#   ggplot() +
#   geom_col(
#     aes(x = tname, y = perc, group = cluster),
#     color = 'black',
#     linewidth = 0.1,
#     show.legend = FALSE
#   ) +
#   facet_wrap( ~ sample) +
#   scale_fill_viridis_d(option = 'B', na.value = "grey50") +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))


# Cluster Correctness -----------------------------------------------------

plot_data <-
  full_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
  # distinct(sample, cluster, tname) %>%
  filter(!is.na(cluster)) %>%
  count(sample, cluster, tname, wt = qlen)

plot_data <-
  plot_data %>%
  # filter(!grepl('^NA24385_', sample))
  filter(!grepl(pilot_regex, sample, perl=TRUE))

# plot_data %>%
#   # filter(sample == 'NA19320') %>%
#   ggplot() +
#   geom_col(aes(x = cluster, y = n, fill = tname), show.legend = FALSE) +
#   # geom_histogram(aes(n)) +
#   facet_wrap( ~ sample, scales = 'free_x')



plot_data %>%
  group_by(sample, cluster) %>%
  mutate(perc = n / sum(n)) %>%
  filter(perc != 1) %>%
  ggplot() +
  geom_col(aes(x = cluster, y = n, fill = tname),
           show.legend = FALSE,
           color = 'black',
           linewidth=0.25) +
  # geom_histogram(aes(n)) +
  facet_wrap(~ sample, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


# Error Fractions ---------------------------------------------------------

# plot_data <-
#   left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
#   # distinct(sample, cluster, tname) %>%
#   # filter(!is.na(cluster)) %>%
#   count(sample, cluster, tname, wt = qlen)
#
#
# plot_data <-
#   plot_data %>%
#   mutate(tname = as.character(tname)) %>%
#   mutate(tname = ifelse(tname %in% c('chrX', 'chrY'), 'chrXY', tname)) %>%
#   mutate(tname = factor(tname, levels = c(paste0('chr', 1:22), 'chrXY', 'chrM'))) %>%
#   group_by(sample, tname) %>%
#   mutate(perc = n / sum(n)) %>%
#   # filter(perc != 1) %>%
#   ungroup()
#
# plot_data <-
#   plot_data %>%
#   filter(!is.na(perc)) %>%
#   group_by(sample, cluster) %>%
#   mutate(major_class=tname[which.max(perc)]) %>%
#   filter(major_class != tname) %>%
#   ungroup()
#
# plot_data %>%
#   # filter(!grepl('^NA24385_', sample)) %>%
#   filter(!is.na(cluster)) %>%
#   ggplot() +
#   geom_col(aes(x=tname, y=perc, fill=major_class)) +
#   facet_wrap(~sample) +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1))


# # Unclustered Error Fractions ---------------------------------------------
#
#
#
# plot_data <-
#   left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
#   mutate(below_threshold = exclusion == 'Length less than threshold: 250000') %>%
#   mutate(below_threshold = ifelse(is.na(below_threshold), FALSE, below_threshold)) %>%
#   mutate(unclustered = is.na(cluster) & !below_threshold) %>%
#   mutate(cluster_exclusion = ifelse(below_threshold, 'threshold', ifelse(unclustered, 'unclustered', 'clustered'))) %>%
#   # distinct(sample, cluster, tname) %>%
#   # filter(!is.na(cluster)) %>%
#   count(sample, cluster, tname, cluster_exclusion, wt = qlen)
#
#
# plot_data <-
#   plot_data %>%
#   mutate(tname = as.character(tname)) %>%
#   mutate(tname = ifelse(tname %in% c('chrX', 'chrY'), 'chrXY', tname)) %>%
#   mutate(tname = factor(tname, levels = c(paste0('chr', 1:22), 'chrXY', 'chrM'))) %>%
#   group_by(sample,cluster, tname, cluster_exclusion) %>%
#   summarise(n = sum(n), .groups = 'drop') %>%
#   group_by(sample, tname) %>%
#   mutate(perc = n / sum(n)) %>%
#   # filter(perc != 1) %>%
#   filter(!is.na(perc)) %>%
#   ungroup()
#
# plot_data <-
#   bind_rows(
#     plot_data %>%
#       filter(cluster_exclusion == 'clustered') %>%
#       group_by(sample, cluster) %>%
#       mutate(major_class=tname[which.max(perc)]) %>%
#       filter(major_class != tname) %>%
#       ungroup(),
#
#     plot_data %>%
#       filter(cluster_exclusion != 'clustered')
#   )
#
# plot_data <-
#   plot_data %>%
#   filter(!grepl('^NA24385_', sample))
#
# plot_data <-
#   plot_data %>%
#   filter(cluster_exclusion != 'threshold') %>%
#   # filter(!grepl('^NA24385_', sample))%>%
#   filter(!grepl('hifiasm', sample))%>%
#   # filter(grepl(pilot_regex, sample, perl=TRUE)) %>%
#   filter(!is.na(tname)) %>%
#   filter(tname != 'chrM')
#
# # plot_data <-
# #   plot_data %>%
# #   group_by(sample, tname) %>%
# #   summarise(perc=sum(perc)) %>%
# #   ungroup()
#
# label_data <-
#   plot_data%>%
#   group_by(sample, tname) %>%
#   summarise(perc=sum(perc)) %>%
#   group_by(sample) %>%
#   slice_max(order_by = perc) %>%
#   ungroup()
#
# #filter(!is.na(major_class)) %>%
# ggplot(mapping = aes(x=tname, y=perc)) +
#   geom_col(data=plot_data) +
#   geom_text(aes(label = paste0(round(100*perc, 1), '%')), nudge_y = 0.03, data=label_data) +
#   facet_wrap(~sample, ncol=5) +
#   # scale_y_continuous(limits = c(0, 0.3),  expand = expansion(mult=c(0,0.05))) +
#   ylab('Percent of Chromosome') +
#   xlab('Chromosome') +
#   theme_bw(base_size = 18) +
#   theme(
#     # strip.background = element_rect(color='black'),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     # panel.grid.minor.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust=1),
#   )
#
#
#
# plot_data %>%
#   filter(cluster_exclusion != 'threshold') %>%
#   filter(grepl('^NA24385_', sample)) %>%
#   mutate(fraction = case_when(
#     grepl('frac-0', sample) ~ 0,
#     grepl('frac-25', sample) ~ 25,
#     grepl('frac-50', sample) ~ 50,
#     grepl('frac-75', sample) ~ 75,
#     grepl('frac-100', sample) ~ 100
#   )) %>%
#   mutate(rep = case_when(
#     grepl('rep-0', sample) ~ 1,
#     grepl('rep-1', sample) ~ 1,
#     grepl('rep-2', sample) ~ 2,
#     grepl('rep-3', sample) ~ 3,
#     grepl('rep-4', sample) ~ 4
#   )) %>%
#   filter(!is.na(tname)) %>%
#   filter(tname != 'chrM') %>%
#   #filter(!is.na(major_class)) %>%
#   ggplot() +
#   geom_col(aes(x=tname, y=perc)) +
#   facet_grid(rows = vars(rep), cols=vars(fraction)) +
#   scale_y_continuous(limits = c(0, 1), expand = expansion(mult=c(0,0.05))) +
#   ylab('Percent of Chromosome') +
#   xlab('Chromosome') +
#   theme_bw(base_size=18) +
#   theme(
#     # strip.background = element_rect(color='black'),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     # panel.grid.minor.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust=1),
#   )
#
#

# Chrom - Error ---------------------------------------------

plot_data <-
  left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
  # filter(sample %in% new_samples) %>%
  mutate(color = case_when(qlen < 250000 ~ 'threshold',
                           hap_1_counts == 0 & hap_2_counts == 0 ~ 'uncounted',
                           TRUE ~ 'misclustered'))

plot_data <-
  plot_data %>%
  mutate(tname = as.character(tname)) %>%
  mutate(tname = ifelse(tname %in% c('chrX', 'chrY'), 'chrXY', tname)) %>%
  mutate(tname = factor(tname, levels = c(paste0('chr', 1:22), 'chrXY', 'chrM')))

plot_data <-
  plot_data %>%
  group_by(sample, cluster, tname, color) %>%
  summarise(n = sum(qlen), .groups = 'drop') %>%
  group_by(sample, tname) %>%
  mutate(perc = n / sum(n)) %>%
  # filter(perc != 1) %>%
  ungroup() %>%
  filter(!is.na(perc))

plot_data <-
  bind_rows(
    plot_data %>%
      filter(color == 'misclustered') %>%
      group_by(sample, cluster) %>%
      mutate(major_class=tname[which.max(perc)]) %>%
      filter(major_class != tname) %>%
      ungroup(),

    plot_data %>%
      filter(color != 'misclustered')
  )

# plot_data <-
#   plot_data %>%
#   filter(!grepl('^NA24385_', sample))

plot_data <-
  plot_data %>%
  # filter(color != 'threshold') %>%
  filter(!grepl('^NA24385-', sample))%>%
  # filter(!grepl('hifiasm', sample))%>%
  # filter(grepl(pilot_regex, sample, perl=TRUE)) %>%
  filter(!is.na(tname)) %>%
  filter(tname != 'chrM')

plot_data <-
  plot_data %>%
  group_by(sample, tname, color) %>%
  summarise(perc=sum(perc)) %>%
  ungroup()

label_data <-
  plot_data%>%
  group_by(sample, tname) %>%
  summarise(perc=sum(perc)) %>%
  # group_by(sample, color) %>%
  group_by(sample) %>%
  slice_max(order_by = perc) %>%
  ungroup()

#filter(!is.na(major_class)) %>%
ggplot(mapping = aes(x = tname, y = perc)) +
  geom_vline(
    xintercept = c('chr13', 'chr14', 'chr15', 'chr21', 'chr22'),
    alpha = 0.33,
    linetype = 'dashed'
  ) +
  geom_col(aes(fill = color), color='black', linewidth=0.2, data = plot_data) +
  geom_text(aes(label = paste0(round(100 * perc, 1), '%')), size=2, nudge_y = 0.1, data =
              label_data) +
  facet_wrap( ~ sample, ncol = 5) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  ylab('Percent of Chromosome') +
  xlab('Chromosome') +
  scale_fill_okabe_ito() +
  # theme_gray() +
  theme_bw(base_size = 5) +
  # theme_bw() +
  theme(
    # strip.background = element_rect(color='black'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )




# XY opposite phasing -----------------------------------------------------

plot_data <-
  left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
  # filter(sample %in% new_samples) %>%
  filter(qlen >= 250000) %>%
  filter(!grepl('hifiasm', sample))

y_chr_sums <-
  plot_data %>%
  group_by(sample) %>%
  filter(tname == 'chrY') %>%
  summarise(length = sum(qlen))

y_chr_sums <-
  y_chr_sums %>%
  filter(length > 3e7)

plot_data <-
  plot_data %>%
  # group_by(sample) %>%
  semi_join(y_chr_sums, by='sample') %>%
  filter(tname %in% c('chrY', 'chrX'))

# plot_data <-
#   # plot_data >%
#   # # filter(grepl(pilot_regex, sample, perl=TRUE)) %>%
#   # filter(!is.na(tname))

plot_data <-
  plot_data %>%
  mutate(hap_1_counts = hap_1_counts + 1,
         hap_2_counts = hap_2_counts + 1) %>%
  mutate(rat = (hap_1_counts - hap_2_counts)/(hap_1_counts+hap_2_counts))

plot_data %>%
  ggplot() +
  geom_point(aes(
    x = rat,
    y = tname,
    size = qlen,
    fill = tname,
    group = sample
  ),
  shape = 21) +
  facet_wrap( ~ sample) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  scale_size_area()  +
  theme_bw(base_size = 18) +
  # theme_bw() +
  theme(
    # strip.background = element_rect(color='black'),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
  )



# Marker Ratio Rukki Plots ------------------------------------------------------

plot_data <-
  left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>%
  left_join(rukki_paths, by=c('sample', 'unitig')) %>%
  filter(qlen >= 250000) %>%
  filter(!grepl('NA24385-', sample))

plot_data <-
  plot_data %>%
  # mutate(hap_1_counts = hap_1_counts + 1,hap_2_counts = hap_2_counts + 1,)
  mutate(rat = (hap_1_counts - hap_2_counts)/(hap_1_counts+hap_2_counts))

ggplot(plot_data) +
  geom_histogram(aes(rat, weight=qlen, fill=assignment)) +
  facet_wrap(~sample) +
  geom_vline(xintercept = 0, linetype = 'dotted')

plot_data %>%
  filter(is.na(assignment)) %>%
  ggplot() +
  geom_histogram(aes(rat, weight=qlen)) +
  facet_wrap(~sample) +
  geom_vline(xintercept = 0, linetype = 'dotted')

# Chrom-Cluster Entropies -------------------------------------------------
#
# plot_data <-
#   left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig'))
#
#
# plot_data <-
#   plot_data %>%
#   filter(!is.na(cluster))
#
# plot_data <-
#   plot_data %>%
#   mutate(tname = as.character(tname)) %>%
#   mutate(tname = ifelse(tname %in% c('chrX', 'chrY'), 'chrXY', tname)) %>%
#   mutate(tname = factor(tname, levels = c(paste0('chr', 1:22), 'chrXY', 'chrM'))) %>%
#   group_by(sample,cluster, tname) %>%
#   summarise(qlen = sum(qlen), .groups = 'drop') %>%
#   group_by(sample, tname) %>%
#   mutate(perc = qlen / sum(qlen)) %>%
#   # filter(perc != 1) %>%
#   ungroup()
#
# plot_data <-
#   plot_data %>%
#   mutate(entropy2 = -1*perc*log2(perc)) %>%
#   group_by(sample, tname) %>%
#   summarise(entropy2 = sum(entropy2, na.rm=TRUE), .groups='drop')
#
# plot_data %>%
#   filter(tname != 'chrM') %>%
#   filter(!is.na(tname)) %>%
#   ggplot() +
#   geom_col(aes(x = tname, y=entropy2)) +
#   facet_wrap(~sample) +
#   theme_bw() +
#   theme(
#     # strip.background = element_rect(color='black'),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     # panel.grid.minor.y = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust=1),
#   )

# Rukki Plots -------------------------------------------------------------


n50s <-
  rukki_paths %>%
  group_by(sample) %>%
  summarise(n50 = calculate_n50(node_len * is_not_gap))

pn50s <-
  rukki_paths %>%
  group_by(sample, name, run) %>%
  summarise(node_len = sum(node_len * is_not_gap), .groups = 'drop') %>%
  group_by(sample) %>%
  summarise(pn50 = calculate_n50(node_len))

sn50s <-
  rukki_paths %>%
  group_by(sample, name) %>%
  summarise(node_len = sum(node_len * is_not_gap), .groups = 'drop') %>%
  group_by(sample) %>%
  summarise(sn50 = calculate_n50(node_len))

length_stats <-
  n50s %>%
  full_join(pn50s, by='sample') %>%
  full_join(sn50s, by='sample')
#
# stopifnot(all(length_stats$n50 <= length_stats$pn50))
# stopifnot(all(length_stats$pn50 <= length_stats$sn50))

length_stats <-
  length_stats %>%
  tidyr::pivot_longer(cols = -sample, names_to = 'statistic', values_to = 'value') %>%
  left_join(ccn502, by='sample')

length_stats %>%
  mutate(
    statistic = forcats::fct_recode(
      statistic,
      N50 = 'n50',
      `Phased\nN50` = 'pn50',
      `Scaffold\nN50` = 'sn50'
    )
  ) %>%
  filter(!grepl('^NA24385_', sample)) %>%
  filter(!grepl('hifiasm', sample)) %>%
  ggplot() +
  geom_col(aes(x = statistic, y = value / 1e6)) +
  # geom_hline(aes(yintercept = ccn502 / 1e6),
  #            linetype = 'dashed',
  #            alpha = 0.8) +
  facet_wrap(facets = vars(sample),
             ncol = 5,
             # scales = 'free_y'
             ) +
  ylab('Mbp (hpc)') +
  theme_bw(base_size = 16) +
  theme(axis.title.x = element_blank(),
        # axis.text.x = element_text(angle = 45, hjust=1)
        )


# Rukki Histograms --------------------------------------------------------
rukki_paths %>%
  filter(!is_gap) %>%
  left_join(reference_alignments) %>%
  ggplot() +
  geom_histogram(aes(qlen, fill=assignment)) +
  geom_vline(aes(xintercept = 250000)) +
  facet_wrap(~sample, scales = 'free_y') +
  scale_x_log10()


# Rukki Assignment Sizes --------------------------------------------------

rukki_paths %>%
  filter(!is_gap) %>%
  left_join(reference_alignments) %>%
  ggplot() +
  geom_col(aes(x=assignment,y=qlen)) +
  facet_wrap(~ sample)


rukki_paths %>%
  left_join(reference_alignments) %>%
  filter(!is_gap) %>%
  count(sample, assignment, wt=qlen)  %>% filter(is.na(assignment)) %>%
  arrange(desc(n))
# Rukki Path-Marker Concordance -----------------------------------------------

# TODO what about large homozygous nodes, as found in hifiasm graphs?

marker_calls <-
  haplotype_marker_counts %>%
  mutate(
    marker_assignment =
      case_when(hap_1_counts > hap_2_counts ~ 'HAPLOTYPE1',
                hap_2_counts > hap_1_counts ~ 'HAPLOTYPE2',
                hap_1_counts == hap_2_counts ~'HOM',
                TRUE ~ NA)
    )# %>%

#select(-hap_1_counts, -hap_2_counts)



plot_data <-
  rukki_paths %>%
  mutate(assignment = ifelse(is.na(assignment), 'NA', assignment)) %>%
  mutate(
    swapped_assignment =
      case_when(
        assignment == 'HAPLOTYPE1'  ~ 'HAPLOTYPE2',
        assignment == 'HAPLOTYPE2'  ~ 'HAPLOTYPE1',
        TRUE ~ 'NA'
      )
  ) %>%
  left_join(marker_calls, by = c('sample', 'unitig')) %>%
  left_join(reference_alignments, by = c('sample', 'unitig')) %>%
  filter(is_not_gap) %>%
  # filter(!is.na(marker_assignment)) %>%
  # filter(!grepl('na_unused', name)) %>%
  group_by(sample, name) %>%
  mutate(
    disagreement = ifelse(assignment != marker_assignment, node_len, 0),
    swapped_disagreement = ifelse(swapped_assignment != marker_assignment, node_len, 0)
  ) %>%
  ungroup() # %>%
  # mutate(
  #   final_disagreement = ifelse(
  #     disagreement < swapped_disagreement,
  #     disagreement,
  #     swapped_disagreement
  #   )
  # )
  #


plot_data %>%
  left_join(rukki_annotations) %>%
  filter(!grepl('^NA24385_', sample)) %>%
  filter(!grepl('hifiasm', sample)) %>%
  # filter(annotation != 'HOM') %>%
  # filter(annotation == 'ISSUE') %>%
  mutate(sample = as.factor(sample)) %>%
  # filter(final_disagreement != 0) %>%
  filter(disagreement != 0) %>%
  filter(hap_1_counts + hap_2_counts >= 7) %>%
  ggplot() +
  geom_point(aes(
    x = sample,
    y = log10(disagreement),
    fill = (hap_1_counts - hap_2_counts) / n
  ), shape = 21) +
  # geom_hline(yintercept=log10(250000), linetype='solid') +
  scale_x_discrete(drop = FALSE) +
  scale_y_log10(limits = c(log10(250000), 8), expand = expansion(mult =c(0, 0.05))) +
  ylab('Log10 Unitig Size (hpc)') +
  theme_bw(base_size = 18) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_distiller(type = 'div')


plot_data %>%
  left_join(rukki_annotations) %>%
  # filter(annotation == 'HOM')  %>%
  # filter(annotation == 'ISSUE') %>%
  # filter(disagreement != 0) %>%
  filter(!grepl('^NA24385_', sample)) %>%
  filter(!grepl('hifiasm', sample)) %>%
  filter(hap_1_counts + hap_2_counts >= 5) %>%
  # filter(node_len >= 500000) %>%
  arrange(desc(disagreement)) %>%
  select(sample, name, unitig, assignment, node_len, hap_1_counts, hap_2_counts, marker_assignment,disagreement) %>%
  filter(sample == 'HG02587')



#
# plot_data %>%
#   filter(grepl('^NA24385_', sample)) %>%
#   # filter(final_disagreement != 0) %>%
#   filter(disagreement != 0) %>%
#   ggplot() +
#   geom_point(aes(x = sample, y = log10(disagreement))) +
#   geom_hline(yintercept = log10(250000), linetype = 'solid') +
#   scale_y_continuous(limits = c(5, 8)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Rough Rukki -----------------------------------------------------------

marker_calls <-
  haplotype_marker_counts %>%
  mutate(
    marker_assignment =
      case_when(hap_1_counts > hap_2_counts ~ 'HAPLOTYPE1',
                hap_2_counts > hap_1_counts ~ 'HAPLOTYPE2',
                hap_1_counts == hap_2_counts ~'HOM', # Given the rukki lables as HAP1, HAP2, or NA, this will always produce disagreement
                TRUE ~ NA)
  )# %>%

hom_nodes <-
  rukki_paths %>%
  filter(is_not_gap) %>%
  group_by(sample, unitig) %>%
  mutate(n = n()) %>%
  filter('HAPLOTYPE1' %in% assignment & 'HAPLOTYPE2' %in% assignment) %>%
  ungroup() %>%
  filter(n == 2) %>%
  distinct(sample, unitig) %>%
  mutate(is_hom = TRUE)

issue_nodes <-
  rukki_paths %>%
  filter(is.na(assignment)) %>%
  left_join(select(haplotype_marker_counts, sample, unitig, length)) %>%
  distinct(sample, unitig) %>%
  mutate(is_issue = TRUE)

stopifnot(nrow(semi_join(hom_nodes, issue_nodes)) == 0)

plot_data <-
  rukki_paths %>%
  filter(grepl('red', sample)) %>% 
  # filter(sample %in% new_samples) %>%
  # mutate(assignment = ifelse(is.na(assignment), 'NA', assignment)) %>%
  left_join(marker_calls, by = c('sample', 'unitig')) %>%
  left_join(reference_alignments, by = c('sample', 'unitig')) %>%
  filter(is_not_gap) %>%
  group_by(sample, name) %>%
  mutate(
    disagreement = ifelse(assignment != marker_assignment, node_len, 0)
  ) %>%
  ungroup()
# Hom Fractions
plot_data <-
  plot_data %>%
  filter(!is.na(tname)) %>%
  filter(tname != 'chrM')


plot_data <-
  plot_data %>%
  left_join(hom_nodes) %>%
  left_join(issue_nodes) %>%
  mutate(color = case_when(
    is_hom ~ 'RUKKI_HOM',
    is_issue ~ 'RUKKI_NA',
    disagreement > 0 ~ 'DISAGREEMENT',
    disagreement == 0 ~ 'AGREEMENT',
    TRUE ~ 'UNCLASSIFIED'
  ))


plot_data %>%
  filter(!grepl('^NA24385-', sample)) %>%
  group_by(sample, tname, color) %>%
  summarise(length = sum(length, na.rm = TRUE), .groups = 'drop') %>%
  ggplot() +
  geom_col(aes(x=tname, y=length, fill=color), color = 'black', linewidth=0.2) +
  facet_wrap(~sample, scales = 'free_x') +
  theme_bw(base_size = 11) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_okabe_ito(order = c(3,1,2,4,5,6,7,8))



# Autosomal Haplotype Balance
plot_data %>%
  filter(as.character(tname) %in% paste0('chr', 1:22)) %>%
  group_by(sample, tname, marker_assignment)  %>%
  summarise(length = sum(length, na.rm = TRUE), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from='marker_assignment', values_from='length') %>%
  mutate(HAPLOTYPE1 = HAPLOTYPE1 + HOM, HAPLOTYPE2 = HAPLOTYPE2 + HOM) %>%
  mutate(HAP_RAT = log2(HAPLOTYPE1 / HAPLOTYPE2)) %>%
  ggplot() +
  geom_vline(xintercept = c('chr13','chr14','chr15','chr21','chr22'), linetype='dashed') +
  geom_col(aes(x=tname, y=HAP_RAT), color = 'black', linewidth=0.2) +
  facet_wrap(~sample, scales = 'free_x') +
  theme_bw(base_size = 11) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_okabe_ito(order = c(3,1,2,4,5,6,7,8))



# Rukki NA24385 Sensitivity -----------------------------------------------

# Need to think harder, as nodes can appear in multiple paths
#
#
# bulk <-
#   rukki_paths %>%
#   mutate(
#     swapped_assignment =
#       case_when(
#         assignment == 'HAPLOTYPE1'  ~ 'HAPLOTYPE2',
#         assignment == 'HAPLOTYPE2'  ~ 'HAPLOTYPE1',
#         TRUE ~ NA
#       )
#   ) %>%
#   filter(grepl('^NA24385_', sample))  %>%
#   filter(is_not_gap)
#
# ref <-
#   bulk %>%
#   filter(grepl('frac-100', sample)) %>%
#   select(unitig, assignment)
#
# swapped_ref <-
#   bulk %>%
#   filter(grepl('frac-100', sample)) %>%
#   select(unitig, swapped_assignment) %>%
#   rename(assignment = swapped_assignment)
#   # mutate(swapped_assignment = assignment) %>%
#   # select(-assignment) %>%
#   # rename(assignment = swapped_assignment)
#
# # need to think about switched names for haplotype assignments
# ref_bulk <-
#   bulk %>%
#   anti_join(ref, by=c('unitig', 'assignment'))
#
# swapped_ref_bulk <-
#   bulk %>%
#   anti_join(swapped_ref, by=c('unitig', 'assignment'))
#
# rbc <-
#   ref_bulk %>%
#   left_join(haplotype_marker_counts, by=c('sample', 'unitig')) %>%
#   left_join(reference_alignments, by=c('sample', 'unitig')) %>%
#   group_by(sample, cluster, tname) %>%
#   summarise(qlen_ref=sum(qlen), .groups='drop')
#
# srbc <-
#   swapped_ref_bulk %>%
#   left_join(haplotype_marker_counts, by=c('sample', 'unitig')) %>%
#   left_join(reference_alignments, by=c('sample', 'unitig')) %>%
#   group_by(sample, cluster, tname) %>%
#   summarise(swapped_qlen_ref=sum(qlen), .groups='drop')
#
# plot_data <-
#   full_join(rbc, srbc, by=c('sample', 'cluster', 'tname')) %>%
#   mutate(
#     qlen_ref = tidyr::replace_na(qlen_ref, 0),
#     swapped_qlen_ref = tidyr::replace_na(swapped_qlen_ref, 0)
#     ) %>%
#   mutate(qlen = ifelse(qlen_ref < swapped_qlen_ref, qlen_ref, swapped_qlen_ref))
#
#
# plot_data <-
#   plot_data %>%
#   group_by(sample, tname) %>%
#   summarise(qlen = sum(qlen))
#
# total_assembly_size <-
#   reference_alignments %>%
#   filter(sample == 'NA24385') %>%
#   group_by(tname) %>%
#   summarise(tot_size = sum(qlen))
#
# plot_data <-
#   plot_data %>%
#   left_join(total_assembly_size) %>%
#   mutate(perc= qlen / tot_size)
#
#
# plot_data %>%
#   mutate(fraction = case_when(
#     grepl('frac-0', sample) ~ 0,
#     grepl('frac-25', sample) ~ 25,
#     grepl('frac-50', sample) ~ 50,
#     grepl('frac-75', sample) ~ 75,
#     grepl('frac-100', sample) ~ 100
#   )) %>%
#   mutate(rep = case_when(
#     grepl('rep-0', sample) ~ 1,
#     grepl('rep-1', sample) ~ 1,
#     grepl('rep-2', sample) ~ 2,
#     grepl('rep-3', sample) ~ 3,
#     grepl('rep-4', sample) ~ 4
#   )) %>%
#   filter(!grepl('good_frac', sample) |
#            grepl('frac-100', sample)) %>%
#   filter(fraction != 100) %>%
#   ggplot() +
#   geom_point(aes(x = tname, y = perc)) +
#   facet_grid(cols = vars(fraction)) +
#   scale_y_continuous(expand = expansion(mult = c(0.00, 0.1))) +
#   xlab('Fraction High-Quality Strand-Seq Libraries (96 Total)') +
#   ylab('Percent of Total Assembly Size') +
#   theme_bw(base_size = 18) +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.title.x = element_blank()
#   )
# #
# bulk <-
#   left_join(bulk, ref, by='unitig')

# Error Stats -------------------------------------------------------------


plot_data <-
  readr::read_delim('phased_assembly_stats.txt')


plot_data %>%
  # mutate(`Haplotype Group` = gsub('Trio ', 'Trio\n', `Haplotype Group`)) %>%
  tidyr::pivot_longer(
    cols = -c(Haplotype, `Haplotype Group`, `Haplotype ID`) ,
    names_to = 'statistic',
    values_to = 'value'
  ) %>%
  mutate(`Haplotype Group` = factor(
    `Haplotype Group`,
    levels = c('Trio (T2T)', 'Trio', 'HiC', 'Strand-seq')
  )) %>%
  mutate(statistic = factor(
    statistic,
    levels = c(
      "# Sequences",
      "GC(%)",
      "Haplotype Size (Mbp)",
      "Longest Sequence (Mbp)",
      "N50 (Mbp)",
      "Switch Error (%)",
      "Hamming Error (%)"
    )
  )) %>%
  ggplot(aes(x = `Haplotype Group`, y = value, shape = as.factor(`Haplotype ID`))) +
  geom_point(size=3) +
  facet_wrap( ~ statistic, scales = 'free') +
  guides(shape='none') +
  theme_bw(base_size = 20) +
  theme(axis.title = element_blank())

